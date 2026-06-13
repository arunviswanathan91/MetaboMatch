from __future__ import annotations

import re
import unicodedata

__all__ = [
    "GREEK_MAP",
    "greek_to_ascii",
    "extract_core",
    "normalise",
    "norm_variants",
    "is_catalog_code",
]

GREEK_MAP = {
    "α": "alpha",   "β": "beta",    "γ": "gamma",   "δ": "delta",
    "ε": "epsilon", "ζ": "zeta",    "η": "eta",     "θ": "theta",
    "ι": "iota",    "κ": "kappa",   "λ": "lambda",  "μ": "mu",
    "ν": "nu",      "ξ": "xi",      "ο": "omicron", "π": "pi",
    "ρ": "rho",     "σ": "sigma",   "τ": "tau",     "υ": "upsilon",
    "φ": "phi",     "χ": "chi",     "ψ": "psi",     "ω": "omega",
    "Α": "Alpha",   "Β": "Beta",    "Γ": "Gamma",   "Δ": "Delta",
    "Ε": "Epsilon", "Ζ": "Zeta",    "Η": "Eta",     "Θ": "Theta",
    "Ι": "Iota",    "Κ": "Kappa",   "Λ": "Lambda",  "Μ": "Mu",
    "Ν": "Nu",      "Ξ": "Xi",      "Ο": "Omicron", "Π": "Pi",
    "Ρ": "Rho",     "Σ": "Sigma",   "Τ": "Tau",     "Υ": "Upsilon",
    "Φ": "Phi",     "Χ": "Chi",     "Ψ": "Psi",     "Ω": "Omega",
}

_STRIP_PREFIXES = re.compile(
    r"^(L-|D-|DL-|R-|S-|RS-|sn-|N-|O-|meso-|cis-|trans-|syn-|anti-)+",
    re.IGNORECASE,
)
_SIMILAR_TO = re.compile(
    r"\[Similar\s+to:\s*(.+?)(?:;\s*[ΔD]Mass:.*)?]",
    re.IGNORECASE | re.DOTALL,
)
# Consumes stereo descriptors including the adjacent hyphen: -(2R)-, -(9Z,12E)-
_STEREO = re.compile(r"-?\(\s*(?:[0-9]+[RSEZ],?\s*)+\)-?", re.IGNORECASE)

# Thermo Compound Discoverer "no name" catalog codes, e.g. MFCD00036904
CATALOG_CODE_PATTERNS = (
    re.compile(r"^MFCD\d+$", re.IGNORECASE),
)


def is_catalog_code(raw: str) -> bool:
    """True if `raw` looks like a vendor catalog code rather than a chemical name."""
    raw = (raw or "").strip()
    return any(p.match(raw) for p in CATALOG_CODE_PATTERNS)


def extract_core(raw: str) -> str:
    """Pull the real compound name out of '[Similar to: X; ΔMass: ...]' wrappers."""
    m = _SIMILAR_TO.search(raw or "")
    return m.group(1).strip() if m else (raw or "").strip()


def greek_to_ascii(text: str) -> str:
    for ch, name in GREEK_MAP.items():
        text = text.replace(ch, name)
    return text


def normalise(raw: str) -> str:
    """Full normalisation pipeline -> Greek to ASCII, stereo stripped, lowercased."""
    if not raw:
        return ""
    s = extract_core(str(raw))
    s = greek_to_ascii(s)
    s = _STEREO.sub("", s)
    s = unicodedata.normalize("NFKD", s)
    s = s.encode("ascii", "ignore").decode("ascii")
    s = s.lower()
    s = re.sub(r"[_\s]+", " ", s)
    s = re.sub(r"[^\w\s\-\+\(\),:]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def norm_variants(raw: str) -> list[str]:
    """Generate all normalised variants of a name for index/lookup."""
    base = normalise(raw)
    vs = {base}

    stripped = _STRIP_PREFIXES.sub("", base).strip()
    if stripped and stripped != base:
        vs.add(stripped)

    short = re.sub(r"\balpha\b", "a", base)
    short = re.sub(r"\bbeta\b",  "b", short)
    short = re.sub(r"\bgamma\b", "g", short)
    if short != base:
        vs.add(short)

    vs.add(base.replace(" ", ""))
    return [v for v in vs if v]
