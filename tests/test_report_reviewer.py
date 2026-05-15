"""
tests/test_report_reviewer.py
------------------------------
All tests for ai/report_reviewer.py — no real API key required.
Mock pattern: patch("ai.report_reviewer.call_llm")
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

VALID_CONFIG = {
    "ai": {
        "features": True,
        "provider": "ollama",
        "model": "llama3",
        "modules": {},
    }
}

DISABLED_CONFIG = {
    "ai": {
        "features": False,
        "provider": "ollama",
        "model": "llama3",
    }
}

MODULE_DISABLED_CONFIG = {
    "ai": {
        "features": True,
        "provider": "ollama",
        "model": "llama3",
        "modules": {"report_reviewer": False},
    }
}

STUDY_CONTEXT = {
    "dataset": {"tissue": "liver", "species": "human"},
    "disease": {"context": "hepatocellular carcinoma"},
    "biological_question": "Characterise the immune microenvironment of HCC.",
}

MOCK_FLAGS = [
    {
        "category": "narrative",
        "severity": "warning",
        "description": "Cell type proportions are stated but not shown in a figure.",
        "suggestion": "Add a stacked bar chart of cell type proportions per sample.",
    },
    {
        "category": "methods",
        "severity": "info",
        "description": "Leiden resolution parameter not reported.",
        "suggestion": "State the resolution value used in the Methods section.",
    },
    {
        "category": "conclusions",
        "severity": "critical",
        "description": "Main conclusion not supported by presented DEG data.",
        "suggestion": "Add volcano plot or heatmap for the claimed top DEGs.",
    },
]

MOCK_LLM_RESPONSE = json.dumps({
    "report_flags": MOCK_FLAGS,
    "overall_report_quality": "The report is well-structured but lacks key figures.",
    "reasoning": "Flags raised based on missing figures and unsupported conclusions.",
})

MINIMAL_HTML = """
<html>
<head><title>OmicSage Report</title></head>
<body>
  <h1>Analysis Report</h1>
  <p>Tissue: liver. Disease: hepatocellular carcinoma.</p>
  <p>We identified 14 clusters using Leiden clustering.</p>
  <p>CD8+ T cells were the most abundant immune population.</p>
</body>
</html>
"""


def _make_html_file(tmp_path: Path, content: str = MINIMAL_HTML) -> Path:
    p = tmp_path / "report.html"
    p.write_text(content, encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_returns_none_when_ai_features_false(tmp_path):
    """Returns None when ai_features=False."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    result = run(
        str(html_path),
        DISABLED_CONFIG,
        STUDY_CONTEXT,
        report_dir=str(tmp_path / "reports"),
    )
    assert result is None


def test_returns_none_when_module_disabled(tmp_path):
    """Returns None when report_reviewer module is explicitly disabled."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    result = run(
        str(html_path),
        MODULE_DISABLED_CONFIG,
        STUDY_CONTEXT,
        report_dir=str(tmp_path / "reports"),
    )
    assert result is None


def test_returns_none_when_runtime_ai_false(tmp_path):
    """Returns None when runtime_ai=False."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    result = run(
        str(html_path),
        VALID_CONFIG,
        STUDY_CONTEXT,
        report_dir=str(tmp_path / "reports"),
        runtime_ai=False,
    )
    assert result is None


def test_returns_none_when_html_missing(tmp_path):
    """Returns None when html_report_path does not exist."""
    from ai.report_reviewer import run

    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        result = run(
            str(tmp_path / "nonexistent.html"),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(tmp_path / "reports"),
        )
    assert result is None


def test_returns_result_on_valid_mock_llm(tmp_path):
    """Returns ReportReviewerResult when mock LLM returns valid JSON."""
    from ai.report_reviewer import ReportReviewerResult, run

    html_path = _make_html_file(tmp_path)
    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        result = run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(tmp_path / "reports"),
        )
    assert isinstance(result, ReportReviewerResult)


def test_report_flags_parsed_correctly(tmp_path):
    """3 flags from mock JSON are parsed into ReportFlag objects."""
    from ai.report_reviewer import ReportFlag, run

    html_path = _make_html_file(tmp_path)
    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        result = run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(tmp_path / "reports"),
        )
    assert result is not None
    assert len(result.report_flags) == 3
    assert all(isinstance(f, ReportFlag) for f in result.report_flags)
    categories = {f.category for f in result.report_flags}
    assert "narrative" in categories
    assert "methods" in categories
    assert "conclusions" in categories


def test_overall_report_quality_populated(tmp_path):
    """overall_report_quality is populated from mock JSON."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        result = run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(tmp_path / "reports"),
        )
    assert result is not None
    assert len(result.overall_report_quality) > 0


def test_review_md_written_to_report_dir(tmp_path):
    """report_review.md is written to report_dir."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    report_dir = tmp_path / "reports"
    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(report_dir),
        )
    assert (report_dir / "report_review.md").exists()


def test_review_md_contains_overall_quality_section(tmp_path):
    """report_review.md contains the Overall Quality section."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    report_dir = tmp_path / "reports"
    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(report_dir),
        )
    md_text = (report_dir / "report_review.md").read_text()
    assert "## Overall Quality" in md_text


def test_review_md_contains_flag_section(tmp_path):
    """report_review.md contains at least one flag section."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    report_dir = tmp_path / "reports"
    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(report_dir),
        )
    md_text = (report_dir / "report_review.md").read_text()
    assert "**Issue:**" in md_text
    assert "**Suggestion:**" in md_text


def test_review_path_set_in_result(tmp_path):
    """review_path is set in the result and points to the written file."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    report_dir = tmp_path / "reports"
    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        result = run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(report_dir),
        )
    assert result is not None
    assert result.review_path is not None
    assert Path(result.review_path).exists()


def test_airesult_base_fields_populated(tmp_path):
    """AiResult base fields (timestamp, model, provider, skill_name, skill_version) are set."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        result = run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(tmp_path / "reports"),
        )
    assert result is not None
    assert result.timestamp
    assert result.model == "llama3"
    assert result.provider == "ollama"
    assert result.skill_name == "report_reviewer"
    assert result.skill_version == "1.0"


def test_one_call_llm_call_made(tmp_path):
    """Exactly one call_llm call is made per run."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE) as mock_llm:
        run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(tmp_path / "reports"),
        )
    assert mock_llm.call_count == 1


def test_html_tags_stripped_before_llm(tmp_path):
    """HTML tags are stripped before passing text to call_llm."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    captured_inputs: list[dict] = []

    def capture_call(skill_name, inputs, config, **kwargs):
        captured_inputs.append(inputs)
        return MOCK_LLM_RESPONSE

    with patch("ai.report_reviewer.call_llm", side_effect=capture_call):
        run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(tmp_path / "reports"),
        )

    assert len(captured_inputs) == 1
    report_text = captured_inputs[0]["report_text"]
    assert "<" not in report_text, "HTML tags must be stripped before passing to LLM"


def test_text_truncated_to_max_chars(tmp_path):
    """Plain text is truncated to 24000 chars before the LLM call."""
    from ai import report_reviewer
    from ai.report_reviewer import _MAX_CHARS, run

    # Build an HTML file whose plain text will exceed _MAX_CHARS
    long_content = "<html><body>" + ("A" * (_MAX_CHARS + 5000)) + "</body></html>"
    html_path = tmp_path / "long_report.html"
    html_path.write_text(long_content, encoding="utf-8")

    captured_inputs: list[dict] = []

    def capture_call(skill_name, inputs, config, **kwargs):
        captured_inputs.append(inputs)
        return MOCK_LLM_RESPONSE

    with patch("ai.report_reviewer.call_llm", side_effect=capture_call):
        run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(tmp_path / "reports"),
        )

    assert len(captured_inputs) == 1
    assert len(captured_inputs[0]["report_text"]) <= _MAX_CHARS


def test_invalid_json_handled_gracefully(tmp_path):
    """Invalid JSON from LLM: result is not None, report_flags is empty."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    with patch("ai.report_reviewer.call_llm", return_value="this is not json at all"):
        result = run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(tmp_path / "reports"),
        )
    assert result is not None
    assert result.report_flags == []
    # overall_report_quality gets the raw response truncated to 500 chars
    assert len(result.overall_report_quality) <= 500


def test_report_dir_created_if_missing(tmp_path):
    """report_dir is created automatically if it does not exist."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    nested_dir = tmp_path / "deep" / "nested" / "reports"
    assert not nested_dir.exists()

    with patch("ai.report_reviewer.call_llm", return_value=MOCK_LLM_RESPONSE):
        run(
            str(html_path),
            VALID_CONFIG,
            STUDY_CONTEXT,
            report_dir=str(nested_dir),
        )
    assert nested_dir.exists()
    assert (nested_dir / "report_review.md").exists()


def test_bs4_absent_falls_back_to_regex(tmp_path):
    """When bs4 is not available, regex fallback strips HTML tags correctly."""
    from ai.report_reviewer import run

    html_path = _make_html_file(tmp_path)
    captured_inputs: list[dict] = []

    def capture_call(skill_name, inputs, config, **kwargs):
        captured_inputs.append(inputs)
        return MOCK_LLM_RESPONSE

    # Simulate bs4 being absent by patching the import inside _strip_html
    import builtins
    real_import = builtins.__import__

    def mock_import(name, *args, **kwargs):
        if name == "bs4":
            raise ImportError("bs4 not available")
        return real_import(name, *args, **kwargs)

    with patch("builtins.__import__", side_effect=mock_import):
        with patch("ai.report_reviewer.call_llm", side_effect=capture_call):
            run(
                str(html_path),
                VALID_CONFIG,
                STUDY_CONTEXT,
                report_dir=str(tmp_path / "reports"),
            )

    assert len(captured_inputs) == 1
    report_text = captured_inputs[0]["report_text"]
    assert "<" not in report_text, "Regex fallback must also strip HTML tags"
    assert "Analysis Report" in report_text or "liver" in report_text
