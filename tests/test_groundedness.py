"""
test_groundedness.py — Phase 3 Milestone Validation
====================================================
Measures how well the AI narrative (ai_narrative.md) is grounded in
verifiable facts from analysis_summary.json.

Score = cited factual sentences / total factual sentences
Target: >= 0.85

Definitions
-----------
Factual sentence:
    Any sentence that contains at least one of:
      - a numeric value  (e.g. "14 clusters", "35% of cells")
      - a gene name      (all-caps token, 2–6 chars, e.g. PDCD1, CD8A)
      - a PMID reference (e.g. PMID:12345678 or PMID 12345678)

Cited factual sentence:
    A factual sentence that also contains at least one of:
      - a PMID reference                        (literature-backed)
      - a metric value that appears in analysis_summary.json  (data-backed)
      - a gene token that appears in analysis_summary.json    (gene-backed)
"""

import json
import re
import pathlib
import pytest

# ---------------------------------------------------------------------------
# Paths — relative to repo root so pytest can find them from any cwd
# ---------------------------------------------------------------------------
REPO_ROOT = pathlib.Path(__file__).parent.parent
NARRATIVE_PATH = REPO_ROOT / "reports" / "GSE166635" / "ai_narrative.md"
SUMMARY_PATH = REPO_ROOT / "reports" / "GSE166635" / "analysis_summary.json"

GROUNDEDNESS_THRESHOLD = 0.85

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Numeric: integer or decimal, optionally followed by % or common unit suffix
_NUMERIC_RE = re.compile(r"\b\d+(?:\.\d+)?(?:%|k|M|G)?\b")

# Gene names: all-uppercase token, 2–6 chars, letters only (e.g. CD8A, PDCD1)
# Must start with a letter and contain only letters/digits
_GENE_RE = re.compile(r"\b[A-Z][A-Z0-9]{1,5}\b")

# PMID: "PMID:12345678" or "PMID 12345678"
_PMID_RE = re.compile(r"\bPMID[:\s]\d{6,9}\b", re.IGNORECASE)

# Common English words that look like gene names — exclude them.
# Includes short 2-3 letter words that the gene regex (2–6 chars) would
# otherwise match (e.g. IS, AN, IN, OF, TO, AT, BY, DO, GO, BE, UP, IF).
_COMMON_WORDS = {
    # Articles / prepositions / conjunctions
    "A", "AN", "AS", "AT", "BY", "DO", "GO", "IF", "IN", "IS", "IT",
    "NO", "OF", "ON", "OR", "SO", "TO", "UP", "US", "WE",
    "AND", "ARE", "BUT", "FOR", "NOT", "NOR", "OUT", "OUR", "THE",
    "VIA", "YET",
    # Pronouns / determiners
    "ALL", "ANY", "FEW", "HER", "HIM", "HIS", "ITS", "OWN",
    "SOME", "SUCH", "THE", "TOO", "WHO", "WHY", "YOU",
    "BOTH", "EACH", "MANY", "MORE", "MOST", "MUCH", "ONLY",
    "SAME", "THAN", "THAT", "THEM", "THEN", "THEY", "THIS",
    "ALSO", "BEEN", "HAVE", "HERE", "INTO", "LESS", "NONE",
    "VERY", "WELL", "WERE", "WHAT", "WHEN", "WITH",
    "AFTER", "AMONG", "BEING", "COULD", "EVERY", "FROM",
    "THEIR", "THERE", "THESE", "THOSE", "UNDER", "UNTIL",
    "USING", "WHICH", "WHILE", "WOULD", "ABOUT", "ABOVE",
    "BELOW", "SINCE", "WHERE", "WHOSE", "WITHIN",
    "BEFORE", "SHOULD", "THOUGH", "THROUGH",
    "BETWEEN", "HOWEVER",
    # Biology / analysis terms that are all-caps but not gene names
    "RNA", "DNA", "PCR", "NGS", "UMI", "PCA", "TSS",
    "DEG", "DEA", "QC", "FC", "FDR", "LFC",
    "CELL", "CELLS", "GENE", "GENES", "TYPE", "DATA",
    "USED", "SHOW", "SHOWN", "FIGURE", "TABLE", "STUDY",
    "GROUP", "GROUPS", "SAMPLE", "SAMPLES", "CLUSTER", "CLUSTERS",
    "ANALYSIS", "BASED", "HIGH", "LOW", "NEW", "TOP", "SET",
    "TOTAL", "MEAN", "MEDIAN", "PER", "NON", "COUNT",
    # Boolean / null sentinels
    "FALSE", "TRUE", "NONE", "NULL", "NAN", "YES", "NO",
}


def _is_gene(token: str) -> bool:
    """Return True if token looks like a gene name (not a common word)."""
    return bool(_GENE_RE.fullmatch(token)) and token not in _COMMON_WORDS


def _extract_summary_values(summary: dict) -> tuple[set[str], set[str]]:
    """
    Recursively walk analysis_summary.json and collect:
      - numeric strings (as they would appear in text)
      - gene-like tokens
    """
    numeric_vals: set[str] = set()
    gene_vals: set[str] = set()

    def _walk(obj):
        if isinstance(obj, dict):
            for v in obj.values():
                _walk(v)
        elif isinstance(obj, list):
            for item in obj:
                _walk(item)
        elif isinstance(obj, (int, float)):
            numeric_vals.add(str(int(obj)) if obj == int(obj) else f"{obj:.2f}")
            numeric_vals.add(str(obj))
        elif isinstance(obj, str):
            # pick up any numbers embedded in strings
            for m in _NUMERIC_RE.findall(obj):
                numeric_vals.add(m)
            # pick up gene tokens embedded in strings
            for tok in obj.split():
                tok_clean = re.sub(r"[^A-Z0-9]", "", tok.upper())
                if _is_gene(tok_clean):
                    gene_vals.add(tok_clean)

    _walk(summary)
    return numeric_vals, gene_vals


def _split_sentences(text: str) -> list[str]:
    """Split markdown text into sentences, skipping headings and blank lines."""
    sentences = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Split on sentence-ending punctuation followed by whitespace or end
        parts = re.split(r"(?<=[.!?])\s+", line)
        sentences.extend(p.strip() for p in parts if p.strip())
    return sentences


def _is_factual(sentence: str) -> bool:
    """True if sentence contains a number, gene name, or PMID."""
    if _PMID_RE.search(sentence):
        return True
    if _NUMERIC_RE.search(sentence):
        return True
    for tok in sentence.split():
        tok_clean = re.sub(r"[^A-Z0-9]", "", tok.upper())
        if _is_gene(tok_clean):
            return True
    return False


def _is_cited(
    sentence: str,
    summary_numerics: set[str],
    summary_genes: set[str],
) -> bool:
    """
    True if the factual sentence is grounded:
      - contains a PMID, OR
      - contains a numeric that appears in analysis_summary.json, OR
      - contains a gene that appears in analysis_summary.json
    """
    if _PMID_RE.search(sentence):
        return True
    # Check numeric match against summary values
    for num in _NUMERIC_RE.findall(sentence):
        if num in summary_numerics:
            return True
    # Check gene match against summary genes
    for tok in sentence.split():
        tok_clean = re.sub(r"[^A-Z0-9]", "", tok.upper())
        if tok_clean and _is_gene(tok_clean) and tok_clean in summary_genes:
            return True
    return False


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestGroundedness:

    def test_narrative_file_exists(self):
        assert NARRATIVE_PATH.exists(), (
            f"ai_narrative.md not found at {NARRATIVE_PATH}\n"
            "Run the pipeline with --ai first: "
            "python run_pipeline.py --config config/runs/GSE166635.yaml --step all --ai"
        )

    def test_summary_file_exists(self):
        assert SUMMARY_PATH.exists(), (
            f"analysis_summary.json not found at {SUMMARY_PATH}\n"
            "Run the pipeline with --ai first."
        )

    def test_narrative_non_empty(self):
        text = NARRATIVE_PATH.read_text(encoding="utf-8")
        assert len(text.strip()) > 100, "ai_narrative.md appears empty or too short"

    def test_summary_valid_json(self):
        text = SUMMARY_PATH.read_text(encoding="utf-8")
        data = json.loads(text)  # raises if invalid
        assert isinstance(data, dict), "analysis_summary.json root must be a dict"

    def test_groundedness_score(self):
        """
        Core groundedness check.
        Score = cited factual sentences / total factual sentences >= 0.85
        """
        narrative_text = NARRATIVE_PATH.read_text(encoding="utf-8")
        summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))

        summary_numerics, summary_genes = _extract_summary_values(summary)

        sentences = _split_sentences(narrative_text)
        assert sentences, "No sentences found in ai_narrative.md"

        factual = [s for s in sentences if _is_factual(s)]
        assert factual, (
            "No factual sentences found in narrative. "
            "Check that the narrative contains numbers or gene names."
        )

        cited = [
            s for s in factual
            if _is_cited(s, summary_numerics, summary_genes)
        ]

        score = len(cited) / len(factual)

        # Detailed failure message
        uncited = [s for s in factual if not _is_cited(s, summary_numerics, summary_genes)]
        uncited_preview = "\n  ".join(uncited[:5])

        assert score >= GROUNDEDNESS_THRESHOLD, (
            f"Groundedness score {score:.2%} is below threshold "
            f"{GROUNDEDNESS_THRESHOLD:.0%}\n"
            f"  Total sentences:   {len(sentences)}\n"
            f"  Factual sentences: {len(factual)}\n"
            f"  Cited sentences:   {len(cited)}\n"
            f"  Un-cited factual sentences (first 5):\n  {uncited_preview}"
        )

        # Also print a summary for visibility in -v output
        print(
            f"\nGroundedness: {score:.2%} "
            f"({len(cited)}/{len(factual)} factual sentences cited)"
        )


# ---------------------------------------------------------------------------
# Unit tests for the helper functions (no file I/O required)
# ---------------------------------------------------------------------------

class TestHelpers:

    def test_is_factual_numeric(self):
        assert _is_factual("There are 14 clusters.")

    def test_is_factual_gene(self):
        assert _is_factual("PDCD1 expression was elevated.")

    def test_is_factual_pmid(self):
        assert _is_factual("As shown previously (PMID:12345678).")

    def test_is_factual_plain_sentence(self):
        assert not _is_factual("This is an introduction.")

    def test_common_word_not_gene(self):
        assert not _is_gene("THE")
        assert not _is_gene("AND")
        assert not _is_gene("CELLS")

    def test_real_gene_recognised(self):
        assert _is_gene("CD8A")
        assert _is_gene("PDCD1")
        assert _is_gene("IL2")

    def test_is_cited_pmid(self):
        assert _is_cited("Result confirmed (PMID:99999999).", set(), set())

    def test_is_cited_numeric_in_summary(self):
        assert _is_cited("We found 14 clusters.", {"14"}, set())

    def test_is_cited_gene_in_summary(self):
        assert _is_cited("CD8A was upregulated.", set(), {"CD8A"})

    def test_not_cited_without_evidence(self):
        assert not _is_cited("CD8A was upregulated.", set(), set())

    def test_extract_summary_values_numeric(self):
        summary = {"n_clusters": 14, "pct": 0.35}
        nums, _ = _extract_summary_values(summary)
        assert "14" in nums

    def test_extract_summary_values_gene(self):
        summary = {"top_genes": ["CD8A", "PDCD1", "IL2"]}
        _, genes = _extract_summary_values(summary)
        assert "CD8A" in genes
        assert "PDCD1" in genes

    def test_split_sentences_skips_headings(self):
        text = "# Title\n\nThis is a sentence. Another one.\n"
        sentences = _split_sentences(text)
        assert "# Title" not in sentences
        assert "This is a sentence." in sentences

    def test_split_sentences_empty(self):
        assert _split_sentences("") == []
