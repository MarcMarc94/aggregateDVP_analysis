"""
Preprocessing of the nanoPhos data (NanophosAnalysis notebook, section 2).

Input: the site- and peptide-level outputs of collapse_precursors
(PeptideCollapse_v4.py; rows = features, columns = metadata + one log2 column
per sample) and the proteome report. The notebook calls each step with its
parameters:

    read_fasta_entries                  accession -> entry name, gene symbol
    group_gene                          gene symbol of a protein group
    contaminant_groups                  contaminant flags of the other notebooks
    apply_condition_aware_classI_mask   Class I localisation, decided per condition
    filter_per_condition                completeness filter on sites / peptides
    filter_proteins                     completeness filter on the proteome
    pivot_wide                          features x samples -> samples x features
    hybrid_impute                       kNN where a condition is mostly measured,
                                        down-shifted Gaussian where it is absent
    normalize_phospho_by_condition_mean phosphosite minus the condition mean of
                                        its parent protein
    site_var, peptide_var,
    phospho_adata, proteome_adata       AnnData objects for the analysis cells
"""

from __future__ import annotations

import re
from typing import Iterable, Mapping, Tuple, Union

import anndata as ad
import alphapepttools as at
import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

from preprocessing import flag_contaminants

# Metadata columns of a collapse_precursors output; every other column is a
# sample. A name listed here that a DataFrame does not have is ignored.
META_COLS = {
    'PTM_Collapse_key',
    'UPD_seq',
    'PTM_localization',
    'Protein_group',
    'Gene_group',
    'PG.Genes',
    'PG.ProteinGroups',
}
DEFAULT_META_COLS: Tuple[str, ...] = tuple(sorted(META_COLS))


# ------------------------------------------------------------ contaminants

def read_fasta_entries(path):
    """UniProt FASTA -> {accession: (entry name, gene symbol)}. The gene symbol
    is the GN= field of the header, '' where the header has none."""
    entries = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                _db, accession, rest = line[1:].split("|", 2)
                gene = re.search(r"\bGN=(\S+)", rest)
                entries[accession] = (rest.split()[0], gene.group(1) if gene else "")
    return entries


def group_gene(group, fasta_entries):
    """Gene symbol of a protein group: the GN= of the first accession of the
    group that has one in `fasta_entries`, '' if none has. Accessions of the
    contaminant FASTA ('Cont_') are not in the search FASTA, so a group such as
    'Cont_P02533;P02533' takes the symbol of P02533 (KRT14)."""
    for accession in str(group).split(";"):
        gene = fasta_entries.get(accession, ("", ""))[1]
        if gene:
            return gene
    return ""


def contaminant_groups(protein_groups, fasta_entries):
    """Contaminant flags of Spectronaut protein groups, by flag_contaminants
    (src/preprocessing.py), the rule of the other notebooks.

    flag_contaminants reads DIA-NN columns: uniprot_ids (accessions), proteins
    (UniProt entry names) and single_gene_name (group_gene).
    They are built here from the accessions of each group and `fasta_entries`.
    Spectronaut marks the accessions of the contaminant FASTA 'Cont_', where
    DIA-NN writes 'cRAP-'; the prefix is translated. Returns one row per
    distinct group with the five flag columns."""
    def ids(group, entry_name):
        out = []
        for acc in group.split(";"):
            if acc.startswith("Cont_"):
                out.append("cRAP-" + acc[len("Cont_"):])
            else:
                out.append(fasta_entries.get(acc, (acc, ""))[0] if entry_name else acc)
        return ";".join(out)

    groups = pd.Index(sorted(set(map(str, protein_groups))), name="PG.ProteinGroups")
    var = pd.DataFrame({
        "uniprot_ids":      [ids(g, entry_name=False) for g in groups],
        "proteins":         [ids(g, entry_name=True) for g in groups],
        "single_gene_name": [group_gene(g, fasta_entries) for g in groups],
    }, index=groups)
    return flag_contaminants(var)


# ---------------------------------------------------------------- filters

def _split_meta_quant(
    df: pd.DataFrame, meta_cols: Iterable[str]
) -> Tuple[list, list]:
    meta_present = [c for c in meta_cols if c in df.columns]
    quant_cols = [c for c in df.columns if c not in meta_present]
    return meta_present, quant_cols


def _coerce_sample_to_condition(
    sample_to_condition: Union[Mapping[str, str], pd.Series, pd.DataFrame],
    sample_col: str = "sample",
    condition_col: str = "condition",
) -> pd.Series:
    if isinstance(sample_to_condition, pd.Series):
        return sample_to_condition.astype(str)
    if isinstance(sample_to_condition, pd.DataFrame):
        if (sample_col not in sample_to_condition.columns
                or condition_col not in sample_to_condition.columns):
            raise ValueError(
                f"DataFrame must have columns {sample_col!r} and "
                f"{condition_col!r}; got {list(sample_to_condition.columns)}"
            )
        return (sample_to_condition
                .set_index(sample_col)[condition_col]
                .astype(str))
    return pd.Series(dict(sample_to_condition), dtype=str)


def filter_per_condition(
    df: pd.DataFrame,
    sample_to_condition: Union[Mapping[str, str], pd.Series, pd.DataFrame],
    threshold: float = 0.50,
    mode: str = "any",
    meta_cols: Iterable[str] = DEFAULT_META_COLS,
    sample_col: str = "sample",
    condition_col: str = "condition",
    return_mask: bool = False,
):
    """
    Apply per-condition completeness filter to a site-level matrix.

    Parameters
    ----------
    df : DataFrame
        Site-level matrix from PeptideCollapse v4.
        Rows are sites, columns are metadata + sample (numeric) columns.
    sample_to_condition : dict, Series, or DataFrame
        Maps sample name (matching df sample columns) -> condition label.
        If a DataFrame is passed, it must contain ``sample_col`` and
        ``condition_col`` columns. Samples not in the mapping are ignored
        for the filter (they are kept in the returned DataFrame).
    threshold : float
        Required fraction of present samples within a condition,
        in (0, 1]. Default 0.50.
    mode : {'any', 'all'}
        'any' (default): site kept if threshold is met in at least one
        condition (preserves condition-specific phosphorylation).
        'all': site kept only if threshold is met in every condition.
    meta_cols : iterable of str
        Metadata column names to exclude from the sample/quant block.
    return_mask : bool
        If True, also return the boolean keep-mask aligned to ``df.index``.

    Returns
    -------
    DataFrame, or (DataFrame, Series) if ``return_mask=True``.

    Notes
    -----
    For small conditions the threshold rounds in a way that can collapse
    distinct nominal thresholds onto the same effective count.
    E.g. for n=5 samples per condition:
        50% -> need >=3, 60% -> need >=3, 70% -> need >=4.
    """
    if mode not in ("any", "all"):
        raise ValueError(f"mode must be 'any' or 'all', got {mode!r}")
    if not 0 < threshold <= 1:
        raise ValueError(f"threshold must be in (0, 1], got {threshold}")

    s2c = _coerce_sample_to_condition(
        sample_to_condition, sample_col, condition_col
    )
    _, quant_cols = _split_meta_quant(df, meta_cols)

    sample_in_df = [c for c in quant_cols if c in s2c.index]
    if not sample_in_df:
        raise ValueError(
            "No overlap between sample_to_condition keys and df sample "
            f"columns. First df samples: {quant_cols[:5]}; "
            f"first mapping keys: {list(s2c.index[:5])}"
        )

    quant = df[sample_in_df]
    cond = s2c.loc[sample_in_df]

    # Per-(site, condition) fraction of present samples.
    frac_by_cond = quant.notna().T.groupby(cond).mean().T

    if mode == "any":
        keep = (frac_by_cond >= threshold).any(axis=1)
    else:
        keep = (frac_by_cond >= threshold).all(axis=1)

    filtered = df.loc[keep].copy()

    if return_mask:
        return filtered, keep
    return filtered


def filter_proteins(dataset, meta_cols, how='all', cutoff=0.7, condition_col=None):
    """Keep the protein columns quantified in >= cutoff of the samples, either
    over all samples (how='all') or in at least one condition
    (how='condition'). Rows are samples; `meta_cols` are carried along."""
    df = dataset.copy()

    if (how == 'condition') and (condition_col == None):
        raise ValueError('Please provide the name of condition if you are using condition-based filtering')

    df_num = df.drop(columns=meta_cols)
    df_meta = df[meta_cols]
    if how == 'all':
        df_tmp = df_num.loc[:, (1 - (df_num.isna().sum() / len(df_num))) >= cutoff]
        df_filtered = pd.concat([df_tmp, df_meta], axis=1)
        return df_filtered

    elif how == 'condition':
        if not any(condition_col in str(col) for col in df.columns):
            raise ValueError('Condition column is not present in the dataframe.')

        presence_dict = {}
        for condition, group in df.groupby(condition_col):
            presence = group[df_num.columns].notna().sum() / len(group)
            presence_dict[condition] = presence

        presence_df = pd.DataFrame(presence_dict)
        sites_to_keep = presence_df.max(axis=1) >= cutoff
        cols_to_keep = sites_to_keep[sites_to_keep].index.tolist()

        df_filtered = df[cols_to_keep + meta_cols]

        return df_filtered

    else:
        raise ValueError('Please select valid filtering convention: "all" or "condition".')


# --------------------------------------------------------- Class I mask

def apply_condition_aware_classI_mask(
    df_sites: pd.DataFrame,
    loc_per_run: pd.DataFrame,
    sample_to_condition: Union[Mapping[str, str], pd.Series, pd.DataFrame],
    classI_cutoff: float = 0.75,
    condition_threshold: float = 0.50,
    meta_cols: Iterable[str] = DEFAULT_META_COLS,
    drop_all_nan: bool = True,
    return_decision_table: bool = False,
):
    """
    Class I masking with the condition as the unit of decision.

    For each (site, condition), with N replicates of which N_I have
    loc_prob >= classI_cutoff:

        N_I / N >= condition_threshold: keep every value of the site in this
            condition, including replicates below the cutoff (the site is
            localised in the majority of the condition, so the other
            replicates are taken to be the same site).
        otherwise: set the values with loc_prob < classI_cutoff to NaN.

    Parameters
    ----------
    df_sites : DataFrame
        UNMASKED site-level matrix (the sites output of collapse_precursors
        in PeptideCollapse_v4.py). Rows are sites, columns are metadata +
        sample-quant columns. ``PTM_Collapse_key`` must be one of the
        metadata columns.
    loc_per_run : DataFrame
        Per-(site, run) localization probability matrix (the loc_per_run
        output of collapse_precursors). Index must align with the
        ``PTM_Collapse_key`` values in ``df_sites``; columns must align with
        the sample columns in ``df_sites``.
    sample_to_condition : dict, Series, or DataFrame
        Maps sample column name -> condition label. Samples not in the
        mapping (e.g. blank wells) are left untouched in the matrix but
        do not participate in the per-condition fraction calculation.
    classI_cutoff : float
        Localization-probability cutoff for Class I. Default 0.75
        (Olsen/Mann).
    condition_threshold : float
        Fraction of replicates in a condition that must be Class I for
        the "keep all quants" branch. Default 0.50 (majority rule).
    meta_cols : iterable of str
        Metadata columns to skip when identifying sample columns.
    drop_all_nan : bool
        Drop sites whose entire quant row becomes NaN after masking.
        Default True.
    return_decision_table : bool
        If True, also return a (site x condition) DataFrame of the
        per-(site, condition) Class I fractions, useful for diagnostics.

    Returns
    -------
    DataFrame, or (DataFrame, decision_table) if ``return_decision_table``.
    """
    if not 0 < condition_threshold <= 1:
        raise ValueError(
            f"condition_threshold must be in (0, 1], got {condition_threshold}"
        )
    if not 0 <= classI_cutoff <= 1:
        raise ValueError(
            f"classI_cutoff must be in [0, 1], got {classI_cutoff}"
        )

    s2c = _coerce_sample_to_condition(sample_to_condition)

    meta_present = [c for c in meta_cols if c in df_sites.columns]
    sample_cols = [c for c in df_sites.columns if c not in meta_present]
    if "PTM_Collapse_key" not in df_sites.columns:
        raise ValueError(
            "df_sites must contain a 'PTM_Collapse_key' column for alignment."
        )

    # Work in (site x sample) form indexed by PTM_Collapse_key for alignment.
    quant = (df_sites
             .set_index("PTM_Collapse_key")[sample_cols]
             .astype(float)
             .copy())

    # A site without a localisation row would be masked in every sample below.
    no_loc = quant.index.difference(loc_per_run.index)
    if len(no_loc):
        raise ValueError(f"{len(no_loc)} sites of df_sites have no row in loc_per_run, "
                         f"e.g. {no_loc[:3].tolist()}; the two tables must share "
                         f"PTM_Collapse_key")
    loc = loc_per_run.reindex(index=quant.index, columns=quant.columns)

    # Per-cell Class I status: True iff loc >= cutoff. NaN loc => False.
    is_classI = loc.ge(classI_cutoff).fillna(False)

    keep_mask = pd.DataFrame(True, index=quant.index, columns=quant.columns)

    decision_rows = {}
    for cond in s2c.unique():
        cond_samples = [s for s in s2c.index[s2c == cond] if s in quant.columns]
        if not cond_samples:
            continue
        N = len(cond_samples)
        N_I = is_classI[cond_samples].sum(axis=1)
        fraction = N_I / N
        decision_rows[cond] = fraction

        # Sites with fraction < threshold: in this condition only, drop
        # the non-Class-I cells (where loc < cutoff or loc is NaN).
        low_frac_sites = fraction < condition_threshold
        if low_frac_sites.any():
            block = is_classI.loc[low_frac_sites, cond_samples]
            keep_block = keep_mask.loc[low_frac_sites, cond_samples] & block
            keep_mask.loc[low_frac_sites, cond_samples] = keep_block

    masked_quant = quant.where(keep_mask, np.nan)

    out = df_sites.copy().set_index("PTM_Collapse_key")
    out[sample_cols] = masked_quant
    out = out.reset_index()

    if drop_all_nan:
        all_nan = out[sample_cols].isna().all(axis=1)
        out = out.loc[~all_nan].reset_index(drop=True)

    if return_decision_table:
        decision = pd.DataFrame(decision_rows)
        decision.index.name = "PTM_Collapse_key"
        return out, decision

    return out


# ---------------------------------------------------- reshape and impute

def pivot_wide(df_long, key_col, meta_set, well_pattern):
    """Collapse output (rows = features, columns = samples + meta) -> wide
    (rows = samples, columns = features + sample_id + well). The well is the
    first group of `well_pattern` in the sample name."""
    sample_cols = [c for c in df_long.columns if c not in meta_set]
    wide = df_long.set_index(key_col)[sample_cols].T
    wide.index.name = 'sample_id'
    wide = wide.reset_index()
    wide['well'] = [well_pattern.search(s).group(1) for s in wide['sample_id']]
    return wide


def hybrid_impute(
    df_wide,
    sample_to_condition,
    meta_cols=None,
    sample_id_col=None,
    quant_substring=None,
    knn_min_frac=0.5,
    mec_max_count=0,
    downshift_factor=1.8,
    width_factor=0.3,
    n_neighbors=None,
    random_state=42,
    verbose=True,
):
    """Two-step imputation that follows the missingness of each condition.

    Per (feature, condition): measured in >= knn_min_frac of the replicates ->
    the missing replicates are filled by kNN; measured in <= mec_max_count
    replicates -> filled from a Gaussian down-shifted by downshift_factor SD
    and of width width_factor SD of that sample's observed values. Features
    between the two stay missing.

    Either pass `meta_cols` (everything else becomes quant) OR `quant_substring`
    (cols containing this substring are quant; rest are meta).

    Sample IDs are taken from `sample_id_col` if given, otherwise from the
    DataFrame index.

    Returns the imputed DataFrame (same columns and order as the input) and
    the per-cell decision table ('knn', 'downshift', 'none').
    """
    df = df_wide.copy()

    if sample_id_col is not None:
        df = df.set_index(sample_id_col)

    if meta_cols is not None:
        meta_in_df = [c for c in meta_cols if c in df.columns]
        quant_cols = [c for c in df.columns if c not in meta_in_df]
    elif quant_substring is not None:
        quant_cols = [c for c in df.columns if quant_substring in str(c)]
        meta_in_df = [c for c in df.columns if c not in quant_cols]
    else:
        raise ValueError("Pass either meta_cols= or quant_substring=")

    mat = df[quant_cols].astype(float)
    cond_per_sample = mat.index.to_series().map(sample_to_condition)

    decision = pd.DataFrame('none', index=mat.index, columns=mat.columns)
    nan_mask = mat.isna()

    for cond_name in cond_per_sample.dropna().unique():
        rows = cond_per_sample[cond_per_sample == cond_name].index
        N = len(rows)
        n_obs = mat.loc[rows].notna().sum(axis=0)
        frac = n_obs / N

        knn_sites = frac.index[frac >= knn_min_frac]
        mec_sites = n_obs.index[n_obs <= mec_max_count]

        if len(knn_sites):
            block = nan_mask.loc[rows, knn_sites]
            decision.loc[rows, knn_sites] = (
                decision.loc[rows, knn_sites].mask(block, 'knn')
            )
        if len(mec_sites):
            block = nan_mask.loc[rows, mec_sites]
            decision.loc[rows, mec_sites] = (
                decision.loc[rows, mec_sites].mask(block, 'downshift')
            )

    if n_neighbors is None:
        n_neighbors = max(3, int(np.sqrt(len(mat))))
    knn_imputer = KNNImputer(n_neighbors=n_neighbors)
    knn_filled = pd.DataFrame(
        knn_imputer.fit_transform(mat.T).T,
        index=mat.index, columns=mat.columns,
    )
    mat = mat.where(decision != 'knn', knn_filled)

    rng = np.random.default_rng(random_state)
    for sample_id in mat.index:
        observed = mat.loc[sample_id].dropna()
        if len(observed) < 2:
            continue
        mu, sigma = observed.mean(), observed.std()
        mu_dn = mu - downshift_factor * sigma
        sigma_dn = width_factor * sigma
        dn_cells = decision.loc[sample_id] == 'downshift'
        n_dn = int(dn_cells.sum())
        if n_dn:
            mat.loc[sample_id, dn_cells] = rng.normal(mu_dn, sigma_dn, size=n_dn)

    df_imputed = mat.copy()
    for col in meta_in_df:
        df_imputed[col] = df[col]
    df_imputed = df_imputed[list(df.columns)]
    if sample_id_col is not None:
        df_imputed = df_imputed.reset_index()

    if verbose:
        n_knn  = int((decision == 'knn').sum().sum())
        n_dn   = int((decision == 'downshift').sum().sum())
        n_left = int(((decision == 'none') & nan_mask).sum().sum())
        n_obs0 = int((~nan_mask).sum().sum())
        total  = mat.size
        print("hybrid_impute summary:")
        print(f"  total cells: {total:,}")
        print(f"  originally observed: {n_obs0:,}  ({100*n_obs0/total:.1f}%)")
        print(f"  KNN-imputed (>= {knn_min_frac:.0%} obs): {n_knn:,}")
        print(f"  downshift-imputed (<= {mec_max_count} obs): {n_dn:,}")
        print(f"  left NaN (uncertain middle): {n_left:,}")

    return df_imputed, decision


# ------------------------------------------------------- normalisation

def normalize_phospho_by_condition_mean(
    phospho_df,
    protein_df,
    sample_col='well',
    condition_col='condition',
):
    """Subtract from each phospho feature the mean log2 intensity of its parent
    protein over the wells of the same condition:

        normalized[well, site] = phospho[well, site] - mean(protein[wells of the condition, parent])

    The condition mean, and not the well's own protein value, because the
    proteome comes from a separate injection of the same lysate.

    Parameters
    ----------
    phospho_df : DataFrame
        One row per well; feature columns contain '~' (PTM_Collapse_key or
        peptide_collapse_key), plus `sample_col`.
    protein_df : DataFrame
        One row per well; numeric protein columns plus `sample_col` and
        `condition_col`.

    Returns
    -------
    normalized : DataFrame
        `phospho_df` with every feature that has a parent protein
        normalised; features without a parent protein are dropped.
    phospho_to_protein : dict
        feature column -> protein column.

    The parent protein is looked up by the accession, the first accession,
    the gene and the first gene of the feature key, in that order.
    """
    for c in (sample_col,):
        if c not in phospho_df.columns:
            raise ValueError(f"'{c}' not in phospho_df.columns")
    for c in (sample_col, condition_col):
        if c not in protein_df.columns:
            raise ValueError(f"'{c}' not in protein_df.columns")

    phospho_site_cols = [c for c in phospho_df.columns if '~' in str(c)]

    protein_quant_cols = list(protein_df.select_dtypes(include='number').columns)
    _exclude = {sample_col, condition_col}
    protein_quant_cols = [c for c in protein_quant_cols if c not in _exclude]

    # every protein column under its full group name and under each accession
    prot_id_to_col = {}
    for col in protein_quant_cols:
        col_str = str(col).strip()
        if col_str and col_str not in prot_id_to_col:
            prot_id_to_col[col_str] = col
        for sub_id in col_str.split(';'):
            sub_id = sub_id.strip()
            if sub_id and sub_id not in prot_id_to_col:
                prot_id_to_col[sub_id] = col

    def _candidate_ids(site_name):
        # key format: 'Accession~Gene_AA<pos>_M<mult>'
        s = str(site_name)
        if '~' in s:
            left, right = s.split('~', 1)
            gene = right.split('_')[0]
            return [left, left.split(';')[0], gene, gene.split(';')[0]]
        return [s.split('_')[0]]

    phospho_to_protein = {}
    for site in phospho_site_cols:
        for key in _candidate_ids(site):
            if key in prot_id_to_col:
                phospho_to_protein[site] = prot_id_to_col[key]
                break

    normalizable = list(phospho_to_protein.keys())
    if not normalizable:
        raise ValueError("No phosphosites matched a parent protein. "
                         "Check protein_df columns vs phospho key format.")

    cond_mean = protein_df.groupby(condition_col)[protein_quant_cols].mean()

    well_to_cond = protein_df.set_index(sample_col)[condition_col].to_dict()
    phospho_wells = phospho_df[sample_col].astype(str).values
    phospho_conds = [well_to_cond.get(w) for w in phospho_wells]
    if not any(c is not None for c in phospho_conds):
        raise ValueError(
            f"No phospho '{sample_col}' values matched any proteome "
            f"'{sample_col}'/{condition_col} lookup."
        )

    # rows = phospho wells, columns = proteins: the condition mean of each well
    virtual_proteome = cond_mean.reindex(phospho_conds).reset_index(drop=True)

    normalized = phospho_df.copy().reset_index(drop=True)
    phospho_vals = normalized[normalizable].values
    prot_vals = np.column_stack([
        virtual_proteome[phospho_to_protein[site]].values
        for site in normalizable
    ])
    normalized[normalizable] = phospho_vals - prot_vals

    drop_cols = [c for c in phospho_site_cols if c not in phospho_to_protein]
    normalized = normalized.drop(columns=drop_cols)

    return normalized, phospho_to_protein


# ------------------------------------------------------------- AnnData

def _parse_site_key(key):
    left, right = key.split('~', 1)
    gene, site, mult = right.rsplit('_', 2)
    return {
        'Protein.Group':       left,
        'Protein.Group_first': left.split(';')[0],
        'Genes':               gene,
        'Genes_first':         gene.split(';')[0],
        'aa':                  site[0],
        'position':            int(site[1:]),
        'mult':                int(mult.lstrip('M')),
    }


def site_var(keys, mapped_keys):
    """var table of a site matrix, parsed from the PTM_Collapse_key
    ('Accession~Gene_AA<pos>_M<mult>'); has_proteome_match marks the sites
    with a parent protein."""
    var = pd.DataFrame([_parse_site_key(k) for k in keys], index=keys)
    var['has_proteome_match'] = var.index.isin(mapped_keys)
    return var


def peptide_var(keys, df_peptides, mapped_keys):
    """var table of a peptide matrix, from the metadata columns of the
    PeptideCollapse peptide output `df_peptides`."""
    info = df_peptides.set_index('peptide_collapse_key')[
        ['PG.ProteinGroups', 'PG.Genes', 'clean_sequence',
         'all_modifications', 'phospho_count', 'EG.PTMAssayProbability']
    ].reindex(keys)
    var = pd.DataFrame(index=keys)
    var['Protein.Group']       = info['PG.ProteinGroups'].values
    var['Protein.Group_first'] = info['PG.ProteinGroups'].astype(str).str.split(';').str[0].values
    var['Genes']               = info['PG.Genes'].astype(str).str.replace('#', '_').values
    var['Genes_first']         = var['Genes'].astype(str).str.split(';').str[0].values
    var['clean_sequence']      = info['clean_sequence'].values
    var['all_modifications']   = info['all_modifications'].astype(str).values
    var['phospho_count']       = info['phospho_count'].values
    var['localization_prob']   = info['EG.PTMAssayProbability'].values
    var['has_proteome_match']  = var.index.isin(mapped_keys)
    return var


def phospho_adata(wide_df, var_full, wells, obs_metadata, layer_label):
    """AnnData of a wide phospho matrix: rows in the order of `wells`, the
    feature columns ('~' in the name) as var, obs from `obs_metadata`."""
    keys = [c for c in wide_df.columns if '~' in str(c)]
    var = var_full.loc[keys]
    X = (wide_df.set_index('well')[keys]
                .astype(np.float32).reindex(wells).values)
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=wells), var=var)
    adata = at.pp.add_metadata(adata=adata, incoming_metadata=obs_metadata,
                               axis=0, keep_data_shape=False)
    adata.uns['layer_meaning'] = layer_label
    return adata


def proteome_adata(wide_df, wells, obs_metadata, label):
    """AnnData of the wide proteome matrix (numeric columns other than well and
    condition), rows in the order of `wells`."""
    quant = [c for c in wide_df.select_dtypes(include='number').columns
             if c not in ('well', 'condition')]
    var = pd.DataFrame(index=quant)
    var.index.name = None
    var['Protein.Group']       = var.index
    var['Protein.Group_first'] = [str(c).split(';')[0] for c in var.index]
    X = (wide_df.set_index('well')[quant]
                .astype(np.float32).reindex(wells).values)
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=wells), var=var)
    adata = at.pp.add_metadata(adata=adata, incoming_metadata=obs_metadata,
                               axis=0, keep_data_shape=False)
    adata.uns['layer_meaning'] = label
    return adata
