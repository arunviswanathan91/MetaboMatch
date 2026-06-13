from __future__ import annotations

import csv as csv_mod
from pathlib import Path
from typing import Union

import openpyxl

from .normalize import normalise
from .matching import row_status


def _annotation_columns(top_n: int) -> list[str]:
    extra = ["Normalised_Name", "Match_Status", "Match_Method"]
    hit_cols = []
    for i in range(1, top_n + 1):
        hit_cols += [
            f"HMDB_Accession_{i}", f"HMDB_Name_{i}",
            f"HMDB_Formula_{i}", f"HMDB_InChIKey_{i}",
            f"Match_Score_{i}", f"Match_Method_{i}",
        ]
    return extra + hit_cols


def _build_rows(
    rows: list[dict],
    columns: list[str],
    name_col: str,
    all_hits: list[list[dict]],
    top_n: int,
) -> tuple[list[str], list[list]]:
    header = list(columns) + _annotation_columns(top_n)
    out_rows = []

    for row, hits in zip(rows, all_hits):
        raw_name = str(row.get(name_col) or "")
        norm = normalise(raw_name) if raw_name else ""
        status = row_status(raw_name, hits)
        method = hits[0]["match_method"] if hits else ""

        data = [row.get(c) for c in columns]
        data += [norm, status, method]

        for i in range(top_n):
            if i < len(hits):
                h = hits[i]
                data += [h["hmdb_accession"], h["hmdb_name"],
                         h["hmdb_formula"], h["hmdb_inchikey"],
                         h["match_score"], h["match_method"]]
            else:
                data += ["", "", "", "", "", ""]
        out_rows.append(data)

    return header, out_rows


def write_output(
    rows: list[dict],
    columns: list[str],
    name_col: str,
    all_hits: list[list[dict]],
    output_path: Union[str, Path],
    top_n: int,
    verbose: bool = True,
) -> None:
    output_path = Path(output_path)
    header, data_rows = _build_rows(rows, columns, name_col, all_hits, top_n)

    if output_path.suffix.lower() == ".csv":
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv_mod.writer(f)
            writer.writerow(header)
            writer.writerows(data_rows)
        if verbose:
            print(f"[OUT] CSV saved   → {output_path}")
        return

    # Default: .xlsx (+ .tsv sibling)
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "HMDB_Annotated"
    ws.append(header)
    for r in data_rows:
        ws.append(r)

    for col in ws.columns:
        max_len = max((len(str(c.value or "")) for c in col), default=10)
        ws.column_dimensions[col[0].column_letter].width = min(max_len + 2, 60)

    wb.save(str(output_path))
    if verbose:
        print(f"[OUT] Excel saved → {output_path}")

    tsv_path = output_path.with_suffix(".tsv")
    with open(tsv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv_mod.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(data_rows)
    if verbose:
        print(f"[OUT] TSV saved   → {tsv_path}")
