from __future__ import annotations

import argparse

from .config import MatchConfig
from .pipeline import run


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="metabomatch",
        description="Match metabolomics compound names against HMDB.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--input", required=True,
                    help="Input table: .xlsx, .xls, .csv, or .tsv")
    ap.add_argument("--hmdb", required=True,
                    help="HMDB reference: .zip (containing one .xml) or .xml")
    ap.add_argument("--output", default="hmdb_annotated.xlsx",
                    help="Output path: .xlsx (+ .tsv sibling) or .csv "
                         "(default: hmdb_annotated.xlsx)")
    ap.add_argument("--score", type=int, default=80,
                    help="Minimum fuzzy score 0-100 (default: 80)")
    ap.add_argument("--top", type=int, default=3,
                    help="Top N HMDB candidates per compound (default: 3)")
    ap.add_argument("--name-column", default=None,
                    help="Name of (or 0-based index into) the compound-name "
                         "column. Default: auto-detect.")

    sem = ap.add_argument_group("semantic search (optional 4th pass)")
    sem.add_argument("--semantic", action="store_true",
                     help="Enable semantic (embedding-based) fallback "
                          "matching for otherwise-unmatched compounds. "
                          "Requires: pip install metabomatch[semantic]")
    sem.add_argument("--semantic-model",
                     default="sentence-transformers/all-MiniLM-L6-v2",
                     help="HuggingFace model name for embeddings")
    sem.add_argument("--semantic-cutoff", type=float, default=0.55,
                     help="Minimum cosine similarity 0-1 (default: 0.55)")
    sem.add_argument("--semantic-top", type=int, default=1,
                     help="Top N semantic candidates (default: 1)")

    ap.add_argument("--quiet", action="store_true", help="Suppress progress output")
    return ap


def main(argv=None) -> None:
    args = build_parser().parse_args(argv)

    name_column = args.name_column
    if name_column is not None and name_column.isdigit():
        name_column = int(name_column)

    cfg = MatchConfig(
        input_path=args.input,
        hmdb_path=args.hmdb,
        output_path=args.output,
        score_cutoff=args.score,
        top_n=args.top,
        name_column=name_column,
        use_semantic=args.semantic,
        semantic_model=args.semantic_model,
        semantic_score_cutoff=args.semantic_cutoff,
        semantic_top_n=args.semantic_top,
        verbose=not args.quiet,
    )
    run(cfg)


if __name__ == "__main__":
    main()
