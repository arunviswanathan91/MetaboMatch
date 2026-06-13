from __future__ import annotations

from collections import defaultdict

from tqdm import tqdm

from .config import MatchConfig
from .io_hmdb import parse_hmdb
from .io_input import read_input_table
from .io_output import write_output
from .matching import match_name, row_status
from .normalize import is_catalog_code


def run(config: MatchConfig) -> None:
    v = config.verbose

    # ── Read input table ────────────────────────────────────────────────
    if v:
        print(f"\n[INPUT] Reading {config.input_path} …")
    rows, columns, name_col = read_input_table(config.input_path, config.name_column)

    n_named = sum(
        1 for r in rows
        if r.get(name_col) and not is_catalog_code(str(r[name_col]))
    )
    n_catalog = sum(
        1 for r in rows
        if r.get(name_col) and is_catalog_code(str(r[name_col]))
    )
    n_unnamed = sum(1 for r in rows if not r.get(name_col))

    if v:
        print(f"  {len(rows):,} total rows, name column = '{name_col}'")
        print(f"  {n_named:,} named rows        ← matched against HMDB")
        print(f"  {n_catalog:,} catalog codes     ← skipped")
        print(f"  {n_unnamed:,} unnamed rows       ← passed through")

    # ── Build HMDB indexes ──────────────────────────────────────────────
    raw_index, norm_index, records = parse_hmdb(config.hmdb_path, verbose=v)
    norm_keys = list(norm_index.keys())

    # ── Optional semantic index ─────────────────────────────────────────
    semantic_index = None
    if config.use_semantic:
        from .semantic import SemanticIndex
        semantic_index = SemanticIndex(model_name=config.semantic_model)
        semantic_index.build(records, show_progress=v)

    # ── Match ────────────────────────────────────────────────────────────
    if v:
        print(f"\n[MATCH] Processing {len(rows):,} rows …")
    all_hits: list[list[dict]] = []
    counts: dict[str, int] = defaultdict(int)

    iterator = tqdm(rows, desc="Matching", unit="row") if v else rows
    for row in iterator:
        raw = str(row.get(name_col) or "").strip()

        if not raw or is_catalog_code(raw):
            all_hits.append([])
            counts[row_status(raw, [])] += 1
            continue

        hits = match_name(
            raw, raw_index, norm_index, norm_keys,
            config.score_cutoff, config.top_n,
            semantic_index=semantic_index,
            semantic_cutoff=config.semantic_score_cutoff,
            semantic_top_n=config.semantic_top_n,
        )
        all_hits.append(hits)
        counts[row_status(raw, hits)] += 1

    if v:
        total = len(rows)

        def pct(key):
            return 100 * counts[key] / total if total else 0.0

        print(f"\n  exact               : {counts['exact']:5,}  ({pct('exact'):.1f}%)")
        print(f"  high_confidence ≥90 : {counts['high_confidence']:5,}  ({pct('high_confidence'):.1f}%)")
        print(f"  probable        ≥80 : {counts['probable']:5,}  ({pct('probable'):.1f}%)")
        print(f"  low_confidence  <80 : {counts['low_confidence']:5,}  ({pct('low_confidence'):.1f}%)")
        if config.use_semantic:
            print(f"  semantic_candidate  : {counts['semantic_candidate']:5,}  ({pct('semantic_candidate'):.1f}%)")
        print(f"  no_match            : {counts['no_match']:5,}  ({pct('no_match'):.1f}%)")
        print(f"  catalog_code        : {counts['catalog_code']:5,}  ({pct('catalog_code'):.1f}%)")
        print(f"  unnamed             : {counts['unnamed']:5,}  ({pct('unnamed'):.1f}%)")

    # ── Write output ─────────────────────────────────────────────────────
    write_output(rows, columns, name_col, all_hits,
                 config.output_path, config.top_n, verbose=v)
    if v:
        print("\nDone ✓")
