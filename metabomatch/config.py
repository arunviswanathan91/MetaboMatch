from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union


@dataclass
class MatchConfig:
    # ── Required ────────────────────────────────────────────────────────────
    input_path: Union[str, Path]
    """Path to the metabolomics table. .xlsx / .xls / .csv / .tsv are all accepted."""

    hmdb_path: Union[str, Path]
    """Path to the HMDB reference. .zip (containing one .xml) or a raw .xml."""

    # ── Output ──────────────────────────────────────────────────────────────
    output_path: Union[str, Path] = "hmdb_annotated.xlsx"
    """Where to write the annotated table. .xlsx or .csv (a .tsv sibling is
    always written alongside .xlsx outputs for convenience)."""

    # ── Matching behaviour ──────────────────────────────────────────────────
    score_cutoff: int = 80
    """Minimum RapidFuzz token_set_ratio (0-100) for a fuzzy hit to be kept."""

    top_n: int = 3
    """Number of HMDB candidates to keep per compound."""

    # ── Column selection (generic input support) ──────────────────────────
    name_column: Optional[Union[str, int]] = None
    """
    Which input column holds the compound name.

    - None (default): auto-detect. If the header matches the Thermo
      Compound Discoverer layout, column B ("Name") is used. Otherwise
      the first column whose header looks like a name/compound/metabolite
      column is used.
    - str: an exact column header to use.
    - int: a 0-based column index to use.
    """

    # ── Semantic search (optional 4th pass) ────────────────────────────────
    use_semantic: bool = False
    """
    If True, queries that get no raw/normalised/fuzzy hit are additionally
    searched with a sentence-embedding model against HMDB canonical names.

    Requires the optional dependency: pip install metabomatch[semantic]
    (installs sentence-transformers + torch).

    Semantic hits NEVER override an exact/fuzzy hit — they only fill in
    rows that would otherwise be "no_match", and are always tagged
    Match_Method = "semantic" with status "semantic_candidate" so they
    stand out for manual review.
    """

    semantic_model: str = "sentence-transformers/all-MiniLM-L6-v2"
    """HuggingFace model name used for semantic embeddings."""

    semantic_score_cutoff: float = 0.55
    """Minimum cosine similarity (0-1) for a semantic hit to be kept."""

    semantic_top_n: int = 1
    """Number of semantic candidates to keep per unmatched compound."""

    # ── Misc ────────────────────────────────────────────────────────────────
    verbose: bool = True
    """Print progress / summary information."""

    def __post_init__(self):
        self.input_path = Path(self.input_path)
        self.hmdb_path = Path(self.hmdb_path)
        self.output_path = Path(self.output_path)

        if not 0 <= self.score_cutoff <= 100:
            raise ValueError("score_cutoff must be between 0 and 100")
        if not 0.0 <= self.semantic_score_cutoff <= 1.0:
            raise ValueError("semantic_score_cutoff must be between 0.0 and 1.0")
        if self.top_n < 1:
            raise ValueError("top_n must be >= 1")
