"""
metabomatch
============
Match Thermo Compound Discoverer (or any generic CSV/XLSX) metabolite
names against HMDB and annotate with HMDB accession IDs, names,
formulas, and InChIKeys — for downstream multi-omics network integration.

Quick start
-----------
    from metabomatch import MatchConfig, run

    cfg = MatchConfig(
        input_path="serum_metabolomics_e101.xlsx",  # .xlsx / .xls / .csv / .tsv
        hmdb_path="serum_metabolites.zip",          # .zip or .xml
        output_path="hmdb_annotated.xlsx",          # .xlsx or .csv
        score_cutoff=80,
        top_n=3,
    )
    run(cfg)
"""

from .config import MatchConfig
from .pipeline import run

__all__ = ["MatchConfig", "run"]
__version__ = "0.2.0"
