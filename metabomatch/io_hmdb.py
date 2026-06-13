from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path
from typing import Union

from lxml import etree

from .normalize import norm_variants

HMDB_NS = "http://www.hmdb.ca"
_TAG_RE = re.compile(r"\{[^}]*\}")


def _local(tag: str) -> str:
    return _TAG_RE.sub("", tag)


def _text(el) -> str:
    return (el.text or "").strip()


def _hmdb_stream(hmdb_path: Union[str, Path]):
    """Return (zip_handle_or_None, file_handle) for the HMDB XML."""
    p = Path(hmdb_path)
    if p.suffix.lower() == ".zip":
        zf = zipfile.ZipFile(str(p))
        xml_names = [n for n in zf.namelist() if n.lower().endswith(".xml")]
        if not xml_names:
            sys.exit(f"No .xml found inside {p}")
        print(f"[HMDB] Streaming {xml_names[0]} from {p.name} "
              f"({p.stat().st_size / 1e9:.2f} GB zip)")
        return zf, zf.open(xml_names[0])
    elif p.suffix.lower() in (".xml",):
        print(f"[HMDB] Reading {p.name} ({p.stat().st_size / 1e9:.2f} GB)")
        return None, open(str(p), "rb")
    else:
        raise ValueError(
            f"Unsupported HMDB file type '{p.suffix}'. Expected .zip or .xml"
        )


def parse_hmdb(hmdb_path: Union[str, Path], verbose: bool = True):
    """
    Stream-parse HMDB XML and build raw + normalised alias indexes.

    Returns
    -------
    raw_index  : dict[str, dict]   alias.lower() -> record
    norm_index : dict[str, dict]   normalised alias variant -> record
    records    : dict[str, dict]   accession -> record
    """
    raw_index: dict[str, dict] = {}
    norm_index: dict[str, dict] = {}
    records: dict[str, dict] = {}

    zf, fh = _hmdb_stream(hmdb_path)
    try:
        context = etree.iterparse(
            fh,
            events=("end",),
            tag=("{%s}metabolite" % HMDB_NS, "metabolite"),
            recover=True,
        )
        n = 0
        for _, elem in context:
            rec = {
                "accession": "", "name": "", "synonyms": [],
                "iupac_name": "", "traditional_iupac": "",
                "formula": "", "inchikey": "", "smiles": "",
            }
            for child in elem:
                tag = _local(child.tag)
                t = _text(child)
                if tag == "accession" and not rec["accession"]:
                    rec["accession"] = t
                elif tag == "name":
                    rec["name"] = t
                elif tag == "iupac_name":
                    rec["iupac_name"] = t
                elif tag == "traditional_iupac":
                    rec["traditional_iupac"] = t
                elif tag == "chemical_formula":
                    rec["formula"] = t
                elif tag == "inchikey":
                    rec["inchikey"] = t
                elif tag == "smiles":
                    rec["smiles"] = t
                elif tag == "synonyms":
                    rec["synonyms"] = [_text(c) for c in child if _text(c)]

            elem.clear()

            if not rec["accession"]:
                continue

            all_aliases = (
                [rec["name"]]
                + rec["synonyms"]
                + [rec["iupac_name"], rec["traditional_iupac"]]
            )

            for alias in all_aliases:
                if not alias:
                    continue
                rk = alias.strip().lower()
                if rk and rk not in raw_index:
                    raw_index[rk] = rec

                for v in norm_variants(alias):
                    if v not in norm_index:
                        norm_index[v] = rec

            records[rec["accession"]] = rec
            n += 1
            if verbose and n % 20_000 == 0:
                print(f"  … indexed {n:,} metabolites", flush=True)

    finally:
        fh.close()
        if zf:
            zf.close()

    if verbose:
        print(
            f"[HMDB] Done — {n:,} metabolites | "
            f"{len(raw_index):,} raw aliases | {len(norm_index):,} norm variants"
        )
    return raw_index, norm_index, records
