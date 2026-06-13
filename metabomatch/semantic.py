from __future__ import annotations

from typing import Optional

import numpy as np


class SemanticIndex:
    """
    Cosine-similarity nearest-neighbour search over HMDB canonical names.

    Usage
    -----
        idx = SemanticIndex(model_name="sentence-transformers/all-MiniLM-L6-v2")
        idx.build(records)            # records: accession -> record dict (from io_hmdb)
        hits = idx.search("some unmatched compound name", top_k=1)
        # -> [(hmdb_name, cosine_similarity, record), ...]

    Note: general-purpose embeddings cannot reliably distinguish stereoisomers
    (D-/L-, R/S) or close homologues. Semantic hits are always tagged with
    method="semantic" and status="semantic_candidate" for manual review.
    """

    def __init__(self, model_name: str = "sentence-transformers/all-MiniLM-L6-v2"):
        try:
            from sentence_transformers import SentenceTransformer
        except ImportError as e:
            raise ImportError(
                "Semantic search requires the optional 'semantic' extra:\n"
                "  pip install metabomatch[semantic]\n"
                "(installs sentence-transformers + torch)"
            ) from e

        self.model_name = model_name
        self._model = SentenceTransformer(model_name)
        self._keys: list[str] = []
        self._records: list[dict] = []
        self._embeddings: Optional[np.ndarray] = None

    def build(self, records: dict[str, dict], batch_size: int = 256,
              show_progress: bool = True) -> None:
        """
        records: accession -> record dict, as returned by io_hmdb.parse_hmdb's
                 third return value. Only records with a non-empty 'name'
                 are embedded.
        """
        self._keys = []
        self._records = []
        for rec in records.values():
            name = rec.get("name") or rec.get("traditional_iupac") or ""
            if not name:
                continue
            self._keys.append(name)
            self._records.append(rec)

        print(f"[Semantic] Embedding {len(self._keys):,} HMDB canonical names "
              f"with {self.model_name} …")
        self._embeddings = self._model.encode(
            self._keys,
            batch_size=batch_size,
            show_progress_bar=show_progress,
            normalize_embeddings=True,
            convert_to_numpy=True,
        )
        print(f"[Semantic] Done — embedding matrix shape {self._embeddings.shape}")

    def search(self, query: str, top_k: int = 1) -> list[tuple[str, float, dict]]:
        """Return up to top_k (name, cosine_similarity, record) tuples."""
        if self._embeddings is None:
            raise RuntimeError("SemanticIndex.build() must be called before search()")
        if not query:
            return []

        q_emb = self._model.encode(
            [query], normalize_embeddings=True, convert_to_numpy=True
        )[0]
        sims = self._embeddings @ q_emb  # cosine similarity (both normalised)

        top_idx = np.argsort(-sims)[:top_k]
        return [
            (self._keys[i], float(sims[i]), self._records[i])
            for i in top_idx
        ]
