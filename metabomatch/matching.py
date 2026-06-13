from __future__ import annotations

from rapidfuzz import fuzz, process

from .normalize import extract_core, norm_variants, is_catalog_code


def match_name(
    raw_name: str,
    raw_index: dict,
    norm_index: dict,
    norm_keys: list[str],
    score_cutoff: int,
    top_n: int,
    semantic_index=None,
    semantic_cutoff: float = 0.55,
    semantic_top_n: int = 1,
) -> list[dict]:
    """
    Match a single raw compound name against the HMDB indexes.

    Returns up to top_n hit dicts, sorted by score desc, each with keys:
      hmdb_accession, hmdb_name, hmdb_formula, hmdb_inchikey,
      match_score, match_method
    """
    hits: dict[str, tuple[float, str, dict]] = {}

    core = extract_core(raw_name)

    def _add(acc, score, method, rec):
        if acc not in hits or score > hits[acc][0]:
            hits[acc] = (score, method, rec)

    # ── Pass 1: raw exact ────────────────────────────────────────────────
    rk = core.lower()
    if rk in raw_index:
        rec = raw_index[rk]
        _add(rec["accession"], 100, "raw_exact", rec)

    # ── Pass 2: norm exact ───────────────────────────────────────────────
    variants = norm_variants(core)
    for v in variants:
        if v in norm_index:
            rec = norm_index[v]
            _add(rec["accession"], 100, "norm_exact", rec)

    # ── Pass 3: norm fuzzy (only if no perfect hit yet) ─────────────────
    best_so_far = max((s for s, _, _ in hits.values()), default=0)
    if best_so_far < 100:
        for v in variants:
            if len(v) < 4:
                continue
            results = process.extract(
                v,
                norm_keys,
                scorer=fuzz.token_set_ratio,
                score_cutoff=score_cutoff,
                limit=top_n * 3,
            )
            for matched_key, score, _ in results:
                rec = norm_index[matched_key]
                _add(rec["accession"], score, "norm_fuzzy", rec)

    # ── Pass 4: semantic (only if NOTHING at all so far, and only if enabled) ──
    if not hits and semantic_index is not None:
        sem_hits = semantic_index.search(core, top_k=semantic_top_n)
        for key, sim, rec in sem_hits:
            if sim >= semantic_cutoff:
                _add(rec["accession"], sim, "semantic", rec)

    top = sorted(hits.items(), key=lambda x: -x[1][0])[:top_n]
    return [
        {
            "hmdb_accession": acc,
            "hmdb_name": rec["name"],
            "hmdb_formula": rec["formula"],
            "hmdb_inchikey": rec["inchikey"],
            "match_score": score,
            "match_method": method,
        }
        for acc, (score, method, rec) in top
    ]


def row_status(raw_name: str, hits: list[dict]) -> str:
    """
    Granular status for a row, distinguishing WHY no HMDB hit exists:

      - "unnamed"            : the name field was empty
      - "catalog_code"       : a vendor catalog code (e.g. MFCD#####) - skipped
      - "no_match"           : a real compound name was searched but
                                nothing cleared the score cutoff (incl. semantic)
      - "semantic_candidate" : ONLY a semantic-similarity hit was found -
                                flagged for manual review, not a confirmed match
      - "exact" / "high_confidence" / "probable" / "low_confidence" :
                                based on the top hit's match_score (raw/norm/fuzzy)
    """
    raw = str(raw_name or "").strip()
    if not raw:
        return "unnamed"
    if is_catalog_code(raw):
        return "catalog_code"
    if not hits:
        return "no_match"

    top = hits[0]
    if top["match_method"] == "semantic":
        return "semantic_candidate"

    s = top["match_score"]
    if s == 100:
        return "exact"
    if s >= 90:
        return "high_confidence"
    if s >= 80:
        return "probable"
    return "low_confidence"
