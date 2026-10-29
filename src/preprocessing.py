"""Contaminant flags for DIA-NN protein group tables, shared by the aggregate DVP notebooks.

flag_contaminants() adds five boolean columns to a protein table (adata.var):

    is_contaminant_keratin  first gene symbol in CONTAMINANT_KERATINS
    is_crap                 every ID of the protein group is cRAP-prefixed
    is_immunoglobulin       first gene symbol starts with IGK, IGL or IGH, except IGLON
    is_extra_contaminant    any ID of the protein group is in the extra list
    is_contaminant          any of the four

The extra list differs by sample type; EXTRA_CONTAMINANTS is the default. The table
needs the column single_gene_name (the first symbol of `genes`). The ID columns
proteins and uniprot_ids are checked when present.
"""
import pandas as pd

CONTAMINANT_KERATINS = [
    "KRT1", "KRT2", "KRT5", "KRT6A", "KRT6B",
    "KRT9", "KRT10", "KRT14", "KRT16", "KRT17",
]

# UniProt entry names
EXTRA_CONTAMINANTS = [
    # Keratins
    "K1C13_HUMAN", "K1C18_HUMAN", "K1C27_HUMAN",
    "K2C1B_HUMAN", "K2C3_HUMAN", "K2C4_HUMAN",
    "K2C72_HUMAN", "K2C78_HUMAN", "K2C79_HUMAN", "K2C8_HUMAN",
    # Skin barrier / epidermal
    "HORN_HUMAN", "FILA_HUMAN", "FILA2_HUMAN",
    "DCD_HUMAN", "CDSN_HUMAN", "DMKN_HUMAN", "SBSN_HUMAN", "PROL4_HUMAN",
    "DSG1_HUMAN", "DSC1_HUMAN", "DSC3_HUMAN", "DESP_HUMAN", "PLAK_HUMAN",
    "TGM1_HUMAN", "TGM3_HUMAN",
    # Skin / secretory and S100 family
    "SPB12_HUMAN", "SG1D2_HUMAN", "LCN1_HUMAN", "PIP_HUMAN",
    "ZA2G_HUMAN", "CYTA_HUMAN", "CYTM_HUMAN", "BPIA1_HUMAN",
    "S10A7_HUMAN", "S10A8_HUMAN", "S10A9_HUMAN", "S10AB_HUMAN",
    # Serum / blood
    "ALBU_HUMAN", "HBA_HUMAN", "HBB_HUMAN", "HBD_HUMAN",
    "APOA1_HUMAN", "APOB_HUMAN",
    "TRFE_HUMAN", "TRFL_HUMAN", "LYSC_HUMAN", "ANT3_HUMAN",
    # Immunoglobulins (the IG[KLH] gene-symbol rule catches the families)
    "IGHA1_HUMAN", "IGHG1_HUMAN",
    # Laboratory / non-human
    "cRAP-TRYP_PIG", "cRAP-K2M3_SHEEP", "INS_HUMAN", "TRY3_HUMAN",
    # Matrigel coating
    "CILP2_HUMAN",
]


def flag_contaminants(var: pd.DataFrame, extra=EXTRA_CONTAMINANTS) -> pd.DataFrame:
    """Return a copy of `var` with the five contaminant columns added.

    `extra` lists UniProt entry names or gene symbols of the contaminants specific to
    the sample type."""
    var = var.copy()
    id_cols = [c for c in ("single_gene_name", "proteins", "uniprot_ids") if c in var.columns]
    extra = set(extra)

    def only_crap(row):
        """True only if every part across all ID columns is cRAP-prefixed. A
        non-cRAP part anywhere (e.g. 'P10636' alongside 'cRAP-P10636', tau)
        keeps the protein."""
        has_any_part = False
        for col in id_cols:
            val = row[col]
            if pd.isna(val):
                continue
            parts = [p.strip() for p in str(val).split(";") if p.strip()]
            for p in parts:
                has_any_part = True
                if not p.startswith("cRAP-"):
                    return False
        return has_any_part

    def has_extra(row):
        """True if any part of any ID column is in `extra`."""
        for col in id_cols:
            val = row[col]
            if pd.isna(val):
                continue
            if {p.strip() for p in str(val).split(";")} & extra:
                return True
        return False

    var["is_contaminant_keratin"] = var["single_gene_name"].isin(CONTAMINANT_KERATINS)
    var["is_crap"]                = var[id_cols].apply(only_crap, axis=1)
    # IGLON5 is a neuronal adhesion protein of the IgLON family, not an immunoglobulin.
    var["is_immunoglobulin"]      = var["single_gene_name"].astype(str).str.match(r"^IG[KLH](?!ON)")
    var["is_extra_contaminant"]   = var[id_cols].apply(has_extra, axis=1)
    var["is_contaminant"] = (
        var["is_contaminant_keratin"]
        | var["is_crap"]
        | var["is_immunoglobulin"]
        | var["is_extra_contaminant"]
    )
    return var
