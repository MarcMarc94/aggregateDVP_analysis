"""Collapse of a Spectronaut phosphopeptide report (one row per precursor and
run) to the two matrices of section 2 of the Nanophos notebook:

    modified peptides  one row per protein group, sequence, start position and
                       set of modifications; the precursors of such a peptide
                       (charge states, and phosphate positions within the same
                       sequence) are summed per run
    phosphosites       one row per site and multiplicity (protein group, residue,
                       position, M1 / M2 / M3 = one, two, three or more phosphates
                       on the precursor); per run, the median over the
                       precursors that carry the site

Both matrices are log2 of EG.TotalQuantity (Settings), with the log2 values 0
and 1 (quantities 1 and 2) set to missing, the noise floor. No localisation
filter is applied here: collapse_precursors also returns the localisation
probability of every site in every run, from which the notebook decides
Class I per condition.

Reduced from PeptideCollapse v4, written by Denys Oliinyk, to the code path the
notebook runs (protein-group level, median over precursors for the sites).
"""
import hashlib
import logging
import re
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger("PeptideCollapse")

QUANT = "EG.TotalQuantity (Settings)"
PHOSPHO = "[Phospho (STY)]"
REQUIRED_COLUMNS = [
    "R.FileName", "EG.PrecursorId", QUANT, "PEP.PeptidePosition",
    "EG.PTMAssayProbability", "EG.PTMLocalizationProbabilities",
    "PG.Genes", "PG.ProteinGroups",
]


def collapse_precursors(report: pd.DataFrame):
    """Modified-peptide and phosphosite matrices of a Spectronaut precursor report.

    Returns
    -------
    peptides : DataFrame
        One row per modified peptide: peptide_collapse_key
        ('Accession~Gene_SEQUENCE_start_modifications'), one log2 column per
        run, then PG.Genes, PG.ProteinGroups, all_modifications,
        clean_sequence, phospho_count and EG.PTMAssayProbability (highest over
        the runs).
    sites : DataFrame
        One row per site and multiplicity: one log2 column per run, then
        UPD_seq (peptide sequence, site in lower case followed by '*'),
        PTM_localization (highest over the runs), Protein_group, Gene_group and
        PTM_Collapse_key ('Accession~Gene_AA<position>_M<multiplicity>').
    loc_per_run : DataFrame
        PTM_Collapse_key x run: the highest localisation probability of the
        site among the precursors that carry it in that run.
    """
    _setup_logging()
    missing = [c for c in REQUIRED_COLUMNS if c not in report.columns]
    if missing:
        raise ValueError(f"the report lacks the column(s) {missing}")
    logger.info("Loaded %d rows", len(report))
    logger.info("Unique samples (R.FileName): %d", report["R.FileName"].nunique())
    _check_duplicate_raw_files(report)

    df = _preprocess(report)
    peptides = _collapse_peptides(df)
    sites, loc_per_run = _collapse_sites(df)
    return peptides, sites, loc_per_run


# ------------------------------------------------------------------ set-up

def _setup_logging() -> None:
    """The counts of each step (INFO and above) to the notebook output."""
    logger.handlers.clear()
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter("[%(levelname)s] %(name)s: %(message)s"))
    logger.addHandler(console_handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False      # a root handler would print every line a second time


def _check_duplicate_raw_files(report: pd.DataFrame) -> None:
    """Warn when two runs have the same first 100 precursors (sorted), a sign
    that one raw file was searched twice under two names."""
    file_hashes = {}
    for fname in report["R.FileName"].unique():
        ids = report.loc[report["R.FileName"] == fname, "EG.PrecursorId"]
        first = sorted(ids.dropna().astype(str).tolist())[:100]
        file_hashes.setdefault(hashlib.md5("||".join(first).encode()).hexdigest(), []).append(fname)
    for fnames in file_hashes.values():
        if len(fnames) > 1:
            logger.warning("Potential file duplication detected: %s share identical "
                           "first-100 precursor patterns", fnames)


# ----------------------------------------------------------- preprocessing

def _preprocess(report: pd.DataFrame) -> pd.DataFrame:
    """Phosphorylated precursors with a peptide start position, and the
    sequence, phosphate positions and modifications parsed from EG.PrecursorId.
    Underscores in gene names become '#' while the keys are built, because the
    keys use '_' as separator; the outputs have them back."""
    df = report.copy()
    n_underscore = df["PG.Genes"].str.contains("_", na=False).sum()
    if n_underscore > 0:
        df["PG.Genes"] = df["PG.Genes"].str.replace("_", "#", regex=False)
        logger.warning("%d gene names contained underscores and have been replaced with '#'.",
                       n_underscore)

    parsed = df["EG.PrecursorId"].apply(_parse_precursor)
    df["clean_sequence"] = [p["clean_sequence"] for p in parsed]
    df["phospho_positions"] = [p["phospho_positions"] for p in parsed]
    df["phospho_count"] = [p["phospho_count"] for p in parsed]
    df["all_modifications"] = [p["all_modifications"] for p in parsed]

    n_rows = len(df)
    df = df[df["phospho_count"] > 0].copy()
    logger.info("Preprocessing: %d rows remaining after phospho filter, %d non-phospho removed",
                len(df), n_rows - len(df))

    df["peptide_start_position"] = df["PEP.PeptidePosition"].astype(str).apply(_first_position)
    n_rows = len(df)
    df = df.dropna(subset=["peptide_start_position"]).copy()
    df["peptide_start_position"] = df["peptide_start_position"].astype(int)
    if n_rows > len(df):
        logger.warning("%d rows removed due to invalid peptide start position", n_rows - len(df))
    logger.info("Preprocessing complete: %d phospho rows retained", len(df))
    return df


def _parse_precursor(precursor_id: str) -> Dict[str, Union[List[int], List[str], str, int]]:
    """Sequence without modifications, 1-based phosphate positions, number of
    phosphates and the names of all modifications of a Spectronaut precursor
    ID such as '_AS[Phospho (STY)]PEK_.2'."""
    base = precursor_id
    if len(precursor_id) >= 4:
        if base.startswith("*.") or base.startswith("_."):
            base = base[2:]
        elif base.startswith("*") or base.startswith("_"):
            base = base[1:]
        if ".*" in base:
            base = base.split(".*")[0]
        elif "._" in base:
            base = base.split("._")[0]
        elif base.endswith("*") or base.endswith("_"):
            base = base[:-1]
        parts = base.split(".")
        if len(parts) > 1 and parts[-1].isdigit():          # charge state
            base = ".".join(parts[:-1])

    phospho_only = re.sub(r"\[(?!Phospho \(STY\))[^\]]*\]", "", base)
    clean = re.sub(r"\[[^\]]*\]", "", base).replace("_", "")
    phospho_count = len(precursor_id.split(PHOSPHO)) - 1

    positions = []
    if phospho_count > 0:
        segments = phospho_only.split(PHOSPHO)
        pos = 0
        for seg in segments[:-1]:
            pos += len(seg)
            positions.append(pos)

    return {
        "clean_sequence": clean,
        "phospho_positions": positions,
        "phospho_count": phospho_count,
        "all_modifications": re.findall(r"\[([^\]]+)\]", base),
    }


def _first_position(position_str: str) -> Optional[int]:
    """First start position in PEP.PeptidePosition ('12;40' or '12,40'), or None."""
    if pd.isna(position_str) or position_str in ("None", ""):
        return None
    try:
        first = str(position_str).split(";")[0].split(",")[0].strip()
        if first and first != "None":
            return int(first)
    except (ValueError, AttributeError):
        pass
    return None


def _parse_localization_probabilities(loc_string: str) -> Dict[int, float]:
    """EG.PTMLocalizationProbabilities -> {peptide position: probability}, e.g.
    '_PVS[Phospho (STY): 92.3%]PS[Phospho (STY): 7.6%]...' -> {3: 0.923, 5: 0.076}."""
    if pd.isna(loc_string) or not isinstance(loc_string, str):
        return {}
    s = loc_string.strip("_.*  ")
    result, pos, i = {}, 0, 0
    while i < len(s):
        if s[i] == "[":
            end = s.index("]", i)
            content = s[i + 1:end]
            if content.startswith("Phospho (STY)"):
                match = re.search(r":\s*([\d.]+)%", content)
                if match:
                    result[pos] = float(match.group(1)) / 100.0
            i = end + 1
        elif s[i].isalpha():
            pos += 1
            i += 1
        else:
            i += 1
    return result


# -------------------------------------------------------- modified peptides

def _peptide_key(protein_groups, genes, clean_sequence, start, modifications) -> str:
    """'Accession~Gene_SEQUENCE_start_modifications': first accession and first
    gene of the group, modifications sorted, carbamidomethylation left out."""
    protein = "Unknown" if pd.isna(protein_groups) else str(protein_groups).split(";")[0]
    gene = "Unknown" if pd.isna(genes) else str(genes).replace("#", "_").split(";")[0]
    mods = sorted(str(m) for m in modifications if "Carbamidomethyl" not in str(m))
    return f"{protein}~{gene}_{clean_sequence}_{int(start)}_{'_'.join(mods) or 'NoMods'}"


def _collapse_peptides(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["peptide_collapse_key"] = [
        _peptide_key(*row) for row in zip(df["PG.ProteinGroups"], df["PG.Genes"],
                                          df["clean_sequence"], df["peptide_start_position"],
                                          df["all_modifications"])]

    # precursors of the same key summed per run
    quant = df.pivot_table(index="peptide_collapse_key", columns="R.FileName",
                           values=QUANT, aggfunc="sum").replace(0, np.nan)
    logger.info("Peptide pivot: %d keys x %d samples, NaN count: %d",
                quant.shape[0], quant.shape[1], int(quant.isna().sum().sum()))

    loc = df.pivot_table(index="peptide_collapse_key", columns="R.FileName",
                         values="EG.PTMAssayProbability", aggfunc="first").fillna(-1)
    loc = pd.DataFrame({"EG.PTMAssayProbability": loc.apply(max, axis=1)})

    meta_cols = ["PG.Genes", "PG.ProteinGroups", "clean_sequence", "all_modifications", "phospho_count"]
    meta = df.pivot_table(index="peptide_collapse_key", values=meta_cols, aggfunc="first")

    quant = np.log2(quant)
    quant = _noise_floor(quant)

    result = pd.concat([quant, meta, loc], axis=1).reset_index()
    result["EG.PTMAssayProbability"] = result["EG.PTMAssayProbability"].replace(-1, np.nan)
    result["PG.Genes"] = result["PG.Genes"].astype(str).str.replace("#", "_", regex=False)
    logger.info("Peptide final output: %d peptides, %d samples", len(result), quant.shape[1])
    return result


# ------------------------------------------------------------- phosphosites

def _collapse_sites(df: pd.DataFrame):
    runs = list(df["R.FileName"].unique())

    # one row per (precursor, phosphate position)
    df = df.copy()
    df["_loc_dict"] = df["EG.PTMLocalizationProbabilities"].apply(_parse_localization_probabilities).values
    df = df.explode("phospho_positions")
    df["UPD_seq"] = [_marked_sequence(s, p) for s, p in zip(df["clean_sequence"], df["phospho_positions"])]
    df["PTM_0_aa"] = [_residue(s, p) for s, p in zip(df["clean_sequence"], df["phospho_positions"])]

    # localisation probability of the site in this precursor; the joint
    # probability where the per-site string has no value for the position
    df["PTM_localization"] = [d.get(int(p), np.nan) if d else np.nan
                              for d, p in zip(df["_loc_dict"], df["phospho_positions"])]
    df = df.drop(columns=["_loc_dict"])
    fallback = df["PTM_localization"].isna()
    df.loc[fallback, "PTM_localization"] = df.loc[fallback, "EG.PTMAssayProbability"].astype(np.float64)
    logger.info("Per-site localization: %d rows got per-site prob, %d fell back to joint prob",
                int((~fallback).sum()), int(fallback.sum()))

    index = ["EG.PrecursorId", "phospho_positions"]
    quant = pd.pivot_table(df, index=index, columns=["R.FileName"], values=[QUANT],
                           aggfunc="sum").replace(0, np.nan)
    logger.info("Site pivot: %d keys x %d samples, NaN count: %d",
                quant.shape[0], quant.shape[1], int(quant.isna().sum().sum()))
    loc_max = pd.pivot_table(df, index=index, columns=["R.FileName"], values=["PTM_localization"],
                             aggfunc="first").fillna(-1)
    loc_max = pd.DataFrame({"PTM_localization": loc_max.apply(max, axis=1)})
    info = pd.pivot_table(df, index=index, aggfunc="first",
                          values=["PEP.PeptidePosition", "PG.ProteinGroups", "PG.Genes",
                                  "phospho_count", "PTM_0_aa", "UPD_seq"])
    per_precursor = pd.concat([quant, info, loc_max], axis=1).reset_index()

    # site key from the absolute position: peptide start + position in the peptide - 1
    start = (per_precursor["PEP.PeptidePosition"].astype(str)
             .apply(lambda s: list(filter(None, s.split(";")))[0])
             .apply(lambda s: list(filter(None, s.split(",")))[0]))
    per_precursor = per_precursor[start != "None"].copy()
    start = start[start != "None"].astype(np.int64)
    multiplicity = per_precursor["phospho_count"].astype(np.int64).clip(upper=3)
    per_precursor["PTM_Collapse_key"] = [
        f"{pg}~{genes}_{aa}{int(pos + st - 1)}_M{int(m)}"
        for pg, genes, aa, pos, st, m in zip(
            per_precursor["PG.ProteinGroups"], per_precursor["PG.Genes"].astype(str),
            per_precursor["PTM_0_aa"], per_precursor["phospho_positions"], start, multiplicity)]

    # localisation probability per (site, run): the highest over the precursors
    key_map = per_precursor[["EG.PrecursorId", "phospho_positions", "PTM_Collapse_key"]].drop_duplicates()
    loc_per_run = (df.reset_index(drop=False).merge(key_map, on=index, how="inner")
                     .pivot_table(index="PTM_Collapse_key", columns="R.FileName",
                                  values="PTM_localization", aggfunc="max")
                     .reindex(columns=runs))
    logger.info("Per-(site, run) localization matrix built: %d sites x %d runs (%.1f%% non-NaN)",
                loc_per_run.shape[0], loc_per_run.shape[1],
                100.0 * loc_per_run.notna().sum().sum() / loc_per_run.size if loc_per_run.size else 0)

    # intensity per (site, run): the median over the precursors that carry the site
    by_key = per_precursor.set_index("PTM_Collapse_key")
    quant_cols = [c for c in by_key.columns if isinstance(c, tuple) and c[0] == QUANT]
    quant = by_key[quant_cols].groupby(level=0).median()
    logger.info("Site aggregation: %d collapse keys, NaN count in aggregated matrix: %d",
                len(quant), int(quant.isna().sum().sum()))
    quant = _noise_floor(np.log2(quant))

    loc_site = by_key[["PTM_localization"]].groupby(level=0).max()
    meta = by_key[["PG.Genes", "PG.ProteinGroups", "UPD_seq"]].groupby(level=0).first()
    result = pd.concat([quant, meta, loc_site], axis=1).reset_index()
    result["PTM_localization"] = result["PTM_localization"].replace(-1, np.nan)
    result.columns = [c[1] if isinstance(c, tuple) else c for c in result.columns]

    full_keys = result["PTM_Collapse_key"].copy()
    result = _protein_group_keys(result)
    loc_per_run = _rekey_localization(loc_per_run, full_keys, result["PTM_Collapse_key"])
    logger.info("Site final output: %d sites, %d samples", len(result), quant.shape[1])
    return result, loc_per_run


def _rekey_localization(loc_per_run: pd.DataFrame, full_keys: pd.Series,
                        site_keys: pd.Series) -> pd.DataFrame:
    """Give loc_per_run the keys of the site table. _protein_group_keys shortens
    each site key to the first accession and the first gene of its protein group,
    and the Class I mask looks every site up by the shortened key; a site of a
    group with several accessions would otherwise find no localisation."""
    key_map = pd.Series(site_keys.to_numpy(), index=full_keys.to_numpy())
    shared = key_map[key_map.duplicated(keep=False)]
    if len(shared):
        raise ValueError(f"{shared.nunique()} site keys are produced by more than one "
                         f"protein group after shortening, e.g. "
                         f"{shared.index[:4].tolist()}")
    no_site = loc_per_run.index.difference(key_map.index)
    if len(no_site):
        raise ValueError(f"{len(no_site)} localisation rows have no site, e.g. "
                         f"{no_site[:3].tolist()}")
    logger.info("Localisation matrix re-keyed to the site keys: %d of %d keys changed",
                int((key_map.index != key_map.to_numpy()).sum()), len(key_map))
    return loc_per_run.rename(index=key_map.to_dict())


def _marked_sequence(sequence: str, position: int) -> str:
    """Sequence with the residue at `position` (1-based) in lower case followed by '*'."""
    if position < 1 or position > len(sequence):
        return sequence
    i = position - 1
    return sequence[:i] + sequence[i].lower() + "*" + sequence[i + 1:]


def _residue(sequence: str, position: int) -> str:
    """Residue at `position` (1-based), 'X' outside the sequence."""
    if position < 1 or position > len(sequence):
        return "X"
    return sequence[position - 1]


def _protein_group_keys(data: pd.DataFrame) -> pd.DataFrame:
    """Rebuild PTM_Collapse_key from the first accession and the first gene of
    the group, and add them as Protein_group and Gene_group."""
    accession = data["PTM_Collapse_key"].str.split("~").str[0]
    data["Protein_group"] = accession.str.split(";").str[0]
    gene, site, mult = (data["PTM_Collapse_key"].str.split("~").str[1]
                        .str.split("_", expand=True).T.values)
    data["Gene_group"] = pd.Series(gene, index=data.index).str.split(";").str[0]
    data = data.drop(["PTM_Collapse_key", "PG.Genes", "PG.ProteinGroups"], axis=1)
    data["PTM_Collapse_key"] = (data["Protein_group"] + "~" + data["Gene_group"] + "_"
                                + pd.Series(site, index=data.index) + "_"
                                + pd.Series(mult, index=data.index))
    data["PTM_Collapse_key"] = data["PTM_Collapse_key"].str.replace("#", "_", regex=False)
    data["Gene_group"] = data["Gene_group"].str.replace("#", "_", regex=False)
    return data


def _noise_floor(log2_matrix: pd.DataFrame) -> pd.DataFrame:
    """log2 values 0 and 1 (quantities 1 and 2) set to missing."""
    before = int(log2_matrix.notna().sum().sum())
    out = log2_matrix.replace(0, np.nan).replace(1, np.nan)
    logger.info("Noise floor filter removed %d values", before - int(out.notna().sum().sum()))
    return out
