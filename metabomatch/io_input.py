from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Union

import openpyxl


# Thermo Compound Discoverer column layout (24 columns, A-X) -> friendly names.
# Order matters: this is also the output column order when this layout is detected.
THERMO_COLUMN_MAP: dict[str, str] = {
    "Checked":                              "Checked",
    "Name":                                 "Name",
    "Formula":                              "Formula",
    "Annot. Source: Predicted Compositions": "Annot_Predicted",
    "Annot. Source: mzCloud Search":        "Annot_mzCloud",
    "Annot. Source: Metabolika Search":     "Annot_Metabolika",
    "Annot. Source: ChemSpider Search":     "Annot_ChemSpider",
    "Annot. DeltaMass [ppm]":               "Delta_ppm",
    "Calc. MW":                             "Calc_MW",
    "m/z":                                  "mz",
    "RT [min]":                             "RT_min",
    "# ChemSpider Results":                 "n_ChemSpider",
    "# mzCloud Results":                    "n_mzCloud",
    "mzCloud Best Match Confidence":        "mzCloud_Confidence",
    "MS2":                                  "MS2",
    "Reference Ion":                        "Reference_Ion",
}

# A handful of headers that, if present, strongly indicate a Thermo CD export.
_THERMO_SIGNATURES = {"Calc. MW", "RT [min]", "Reference Ion", "m/z"}

# Sample/solvent area columns are matched by prefix since filenames vary
# ("Area: Sample_01.raw (F3)", "Area: SolA_04.raw (F9)", ...).
_AREA_RE = re.compile(r"^Area:\s*(.+?)\.raw", re.IGNORECASE)

# Candidate header names for auto-detecting the compound-name column
# in a generic (non-Thermo) table.
_NAME_HEADER_HINTS = (
    "name", "compound", "metabolite", "compound name",
    "metabolite name", "feature", "annotation",
)


def _is_thermo_layout(header: list[str]) -> bool:
    return len(_THERMO_SIGNATURES & set(header)) >= 2


def _rename_thermo_header(header: list[str]) -> list[str]:
    out = []
    for h in header:
        if h in THERMO_COLUMN_MAP:
            out.append(THERMO_COLUMN_MAP[h])
            continue
        m = _AREA_RE.match(h or "")
        if m:
            out.append(f"Area_{m.group(1)}")
            continue
        out.append(h)
    return out


def _detect_name_column(header: list[str],
                        name_column: Optional[Union[str, int]]) -> int:
    """Return the 0-based index of the compound-name column."""
    if isinstance(name_column, int):
        if not (0 <= name_column < len(header)):
            raise ValueError(
                f"name_column index {name_column} out of range "
                f"for {len(header)} columns"
            )
        return name_column

    if isinstance(name_column, str):
        try:
            return header.index(name_column)
        except ValueError:
            raise ValueError(
                f"name_column={name_column!r} not found in header: {header}"
            )

    # Auto-detect
    if _is_thermo_layout(header):
        return header.index("Name") if "Name" in header else 1

    lower = [str(h).strip().lower() for h in header]
    for hint in _NAME_HEADER_HINTS:
        for i, h in enumerate(lower):
            if h == hint:
                return i
    for hint in _NAME_HEADER_HINTS:
        for i, h in enumerate(lower):
            if hint in h:
                return i

    raise ValueError(
        "Could not auto-detect the compound-name column. "
        f"Header was: {header}\n"
        "Set MatchConfig.name_column to the column name or index."
    )


def _read_xlsx(path: Path) -> tuple[list[list], list[str]]:
    wb = openpyxl.load_workbook(str(path), read_only=True)
    ws = wb.active
    rows_iter = ws.iter_rows(values_only=True)
    header = [str(h) if h is not None else "" for h in next(rows_iter)]
    data = [list(r) for r in rows_iter]
    return data, header


def _nan_to_none(data: list[list]) -> list[list]:
    out = []
    for row in data:
        out.append([None if (isinstance(v, float) and v != v) else v for v in row])
    return out


def _read_xls(path: Path) -> tuple[list[list], list[str]]:
    import pandas as pd
    df = pd.read_excel(str(path), engine="xlrd", header=0)
    header = [str(c) for c in df.columns]
    data = _nan_to_none(df.values.tolist())
    return data, header


def _read_csv(path: Path) -> tuple[list[list], list[str]]:
    import pandas as pd
    sep = "\t" if path.suffix.lower() == ".tsv" else None  # None -> pandas sniffs
    df = pd.read_csv(str(path), sep=sep, engine="python")
    header = [str(c) for c in df.columns]
    data = _nan_to_none(df.values.tolist())
    return data, header


_READERS = {
    ".xlsx": _read_xlsx,
    ".xlsm": _read_xlsx,
    ".xls":  _read_xls,
    ".csv":  _read_csv,
    ".tsv":  _read_csv,
}


def read_input_table(
    path: Union[str, Path],
    name_column: Optional[Union[str, int]] = None,
) -> tuple[list[dict], list[str], str]:
    """
    Read a metabolomics table.

    Returns
    -------
    rows     : list of dicts, one per data row, keyed by output column name
    columns  : output column order
    name_col : the key in each row dict that holds the compound name
    """
    path = Path(path)
    ext = path.suffix.lower()
    if ext not in _READERS:
        raise ValueError(
            f"Unsupported input file type '{ext}'. "
            f"Supported: {', '.join(sorted(_READERS))}"
        )

    raw_data, header = _READERS[ext](path)

    if _is_thermo_layout(header):
        columns = _rename_thermo_header(header)
    else:
        columns = [str(h) for h in header]

    name_idx = _detect_name_column(columns if _is_thermo_layout(header) else header,
                                    name_column)
    name_col = columns[name_idx]

    n_cols = len(columns)
    rows = []
    for r in raw_data:
        r = list(r) + [None] * (n_cols - len(r))  # pad short rows
        rows.append({col: r[i] for i, col in enumerate(columns[:n_cols])})

    return rows, columns, name_col
