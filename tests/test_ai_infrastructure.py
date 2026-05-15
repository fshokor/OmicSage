"""
tests/test_ai_infrastructure.py
--------------------------------
Tests for the four Phase 3 AI infrastructure modules:
    ai/_base.py
    ai/_config_gate.py
    ai/_audit_log.py
    ai/_llm_client.py

All tests run without a real API key.
LLM calls are intercepted by monkeypatching _build_conversation.

Run:
    cd ~/OmicSage
    conda activate omicsage
    python -m pytest tests/test_ai_infrastructure.py -v
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_config(
    features: bool = True,
    provider: str = "claude",
    model: str = "claude-sonnet-4-20250514",
    modules: dict | None = None,
) -> dict:
    """Return a minimal valid pipeline config dict."""
    cfg: dict = {
        "ai": {
            "features": features,
            "provider": provider,
            "model": model,
        }
    }
    if modules is not None:
        cfg["ai"]["modules"] = modules
    return cfg


def _mock_conversation(response: str = '{"cell_type": "T cell"}') -> MagicMock:
    """Return a mock BioChatter conversation that yields *response*."""
    conv = MagicMock()
    conv.query.return_value = (
        response,
        {"prompt_tokens": 10, "completion_tokens": 5},
        "",
    )
    return conv


# ---------------------------------------------------------------------------
# 1. AiResult base dataclass
# ---------------------------------------------------------------------------

class TestAiResult:
    def test_instantiates_with_defaults(self):
        from ai._base import AiResult

        result = AiResult()
        assert result.model == ""
        assert result.provider == ""
        assert result.skill_name == ""
        assert result.skill_version == ""
        assert result.reasoning == ""
        # timestamp auto-populated
        assert result.timestamp.endswith("Z")
        assert "T" in result.timestamp

    def test_fields_settable(self):
        from ai._base import AiResult

        result = AiResult(
            model="claude-sonnet-4-20250514",
            provider="claude",
            skill_name="cluster_annotator",
            skill_version="1.0",
            reasoning="Because PDCD1 is high.",
        )
        assert result.model == "claude-sonnet-4-20250514"
        assert result.provider == "claude"
        assert result.skill_name == "cluster_annotator"
        assert result.skill_version == "1.0"
        assert result.reasoning == "Because PDCD1 is high."

    def test_subclass_inherits_base_fields(self):
        from dataclasses import dataclass

        from ai._base import AiResult

        @dataclass
        class MyResult(AiResult):
            cell_type: str = ""

        r = MyResult(model="gpt-4o", cell_type="NK cell")
        assert r.model == "gpt-4o"
        assert r.cell_type == "NK cell"
        assert r.timestamp.endswith("Z")


# ---------------------------------------------------------------------------
# 2. Config gate
# ---------------------------------------------------------------------------

class TestConfigGate:
    def test_passes_when_fully_enabled(self):
        from ai._config_gate import check_ai_enabled

        cfg = _make_config(features=True)
        # Should return None silently
        assert check_ai_enabled(cfg, module="cluster_annotator") is None

    def test_raises_when_features_false(self):
        from ai._config_gate import AiDisabledError, check_ai_enabled

        cfg = _make_config(features=False)
        with pytest.raises(AiDisabledError, match="globally"):
            check_ai_enabled(cfg, module="cluster_annotator")

    def test_raises_when_module_explicitly_false(self):
        from ai._config_gate import AiDisabledError, check_ai_enabled

        cfg = _make_config(modules={"cluster_annotator": False})
        with pytest.raises(AiDisabledError, match="cluster_annotator"):
            check_ai_enabled(cfg, module="cluster_annotator")

    def test_passes_when_module_key_missing(self):
        """Missing key defaults to enabled — not an error."""
        from ai._config_gate import check_ai_enabled

        # modules dict exists but cluster_annotator not listed
        cfg = _make_config(modules={"pipeline_advisor": True})
        assert check_ai_enabled(cfg, module="cluster_annotator") is None

    def test_passes_when_modules_section_absent(self):
        """No modules dict at all → all modules default to enabled."""
        from ai._config_gate import check_ai_enabled

        cfg = _make_config()  # no modules key
        assert check_ai_enabled(cfg, module="cluster_annotator") is None

    def test_raises_when_runtime_ai_false(self):
        """runtime_ai=False disables AI regardless of config."""
        from ai._config_gate import AiDisabledError, check_ai_enabled

        cfg = _make_config(features=True)
        with pytest.raises(AiDisabledError, match="runtime"):
            check_ai_enabled(cfg, module="cluster_annotator", runtime_ai=False)

    def test_raises_when_no_ai_section(self):
        from ai._config_gate import AiDisabledError, check_ai_enabled

        cfg: dict = {}  # no "ai" key at all
        with pytest.raises(AiDisabledError, match="No 'ai' section"):
            check_ai_enabled(cfg, module="cluster_annotator")


# ---------------------------------------------------------------------------
# 3. Audit log
# ---------------------------------------------------------------------------

class TestAuditLog:
    def test_creates_file_and_writes_valid_jsonl(self, tmp_path):
        from ai._audit_log import write_audit_record

        write_audit_record(
            log_dir=tmp_path,
            module="cluster_annotator",
            skill_version="1.0",
            model="claude-sonnet-4-20250514",
            provider="claude",
            input_summary={"tissue": "liver", "cluster_id": "3"},
            token_usage={"prompt_tokens": 100, "completion_tokens": 50},
            raw_response='{"cell_type": "T cell"}',
            parsed_output={"cell_type": "T cell"},
            parse_success=True,
        )

        log_file = tmp_path / "cluster_annotator.jsonl"
        assert log_file.exists()

        lines = log_file.read_text().strip().splitlines()
        assert len(lines) == 1

        record = json.loads(lines[0])
        assert record["module"] == "cluster_annotator"
        assert record["skill_version"] == "1.0"
        assert record["model"] == "claude-sonnet-4-20250514"
        assert record["provider"] == "claude"
        assert record["parse_success"] is True
        assert record["prompt_tokens"] == 100
        assert record["completion_tokens"] == 50
        assert record["timestamp"].endswith("Z")

    def test_appends_on_second_call(self, tmp_path):
        from ai._audit_log import write_audit_record

        kwargs = dict(
            log_dir=tmp_path,
            module="pipeline_advisor",
            skill_version="1.0",
            model="llama3",
            provider="ollama",
            input_summary={"tissue": "liver"},
            token_usage=None,
            raw_response="hello",
            parsed_output=None,
            parse_success=False,
        )
        write_audit_record(**kwargs)
        write_audit_record(**kwargs)

        log_file = tmp_path / "pipeline_advisor.jsonl"
        lines = log_file.read_text().strip().splitlines()
        assert len(lines) == 2

    def test_creates_log_dir_if_missing(self, tmp_path):
        from ai._audit_log import write_audit_record

        nested = tmp_path / "logs" / "llm"
        assert not nested.exists()

        write_audit_record(
            log_dir=nested,
            module="test_module",
            skill_version="1.0",
            model="gpt-4o",
            provider="openai",
            input_summary={},
            token_usage=None,
            raw_response="test",
            parsed_output=None,
            parse_success=False,
        )

        assert nested.exists()
        assert (nested / "test_module.jsonl").exists()

    def test_never_raises_on_bad_log_dir(self, tmp_path):
        """write_audit_record must not raise even if the write fails."""
        from ai._audit_log import write_audit_record

        # Pass a file path as log_dir — mkdir will fail on it if it's a file
        bad_path = tmp_path / "a_file.txt"
        bad_path.write_text("I am a file, not a dir")

        # Should not raise
        write_audit_record(
            log_dir=bad_path / "subdir",  # can't mkdir under a file
            module="test",
            skill_version="1.0",
            model="x",
            provider="ollama",
            input_summary={},
            token_usage=None,
            raw_response="x",
            parsed_output=None,
            parse_success=False,
        )


# ---------------------------------------------------------------------------
# 4. LLM client — routing and mock integration
# ---------------------------------------------------------------------------

class TestLlmClientRouting:
    """
    All routing tests monkeypatch _build_conversation so no real API call
    is made and no API key is required.
    """

    def _mock_skill_loader(self, monkeypatch):
        """Patch load_skill to return dummy prompts."""
        monkeypatch.setattr(
            "ai._llm_client.load_skill",
            lambda name, **kwargs: ("You are a scientist.", "Annotate this cluster."),
        )

    def _mock_skill_version(self, monkeypatch):
        monkeypatch.setattr(
            "ai._llm_client._get_skill_version",
            lambda name: "1.0",
        )

    def test_routes_to_anthropic_for_claude(self, monkeypatch, tmp_path):
        from ai._llm_client import _build_conversation

        self._mock_skill_loader(monkeypatch)
        self._mock_skill_version(monkeypatch)

        mock_conv = _mock_conversation()
        captured_provider = {}

        def fake_build(provider, model, config):
            captured_provider["value"] = provider
            return mock_conv

        monkeypatch.setattr("ai._llm_client._build_conversation", fake_build)

        from ai._llm_client import call_llm

        cfg = _make_config(provider="claude")
        call_llm("cluster_annotator", {"tissue": "liver"}, cfg, log_dir=tmp_path)

        assert captured_provider["value"] == "claude"

    def test_routes_to_ollama(self, monkeypatch, tmp_path):
        from ai._llm_client import _build_conversation

        self._mock_skill_loader(monkeypatch)
        self._mock_skill_version(monkeypatch)

        mock_conv = _mock_conversation()
        captured_provider = {}

        def fake_build(provider, model, config):
            captured_provider["value"] = provider
            return mock_conv

        monkeypatch.setattr("ai._llm_client._build_conversation", fake_build)

        from ai._llm_client import call_llm

        cfg = _make_config(provider="ollama", model="llama3")
        call_llm("cluster_annotator", {"tissue": "liver"}, cfg, log_dir=tmp_path)

        assert captured_provider["value"] == "ollama"

    def test_routes_to_openai(self, monkeypatch, tmp_path):
        self._mock_skill_loader(monkeypatch)
        self._mock_skill_version(monkeypatch)

        mock_conv = _mock_conversation()
        captured_provider = {}

        def fake_build(provider, model, config):
            captured_provider["value"] = provider
            return mock_conv

        monkeypatch.setattr("ai._llm_client._build_conversation", fake_build)

        from ai._llm_client import call_llm

        cfg = _make_config(provider="openai", model="gpt-4o")
        call_llm("cluster_annotator", {"tissue": "liver"}, cfg, log_dir=tmp_path)

        assert captured_provider["value"] == "openai"

    def test_unknown_provider_raises_value_error(self):
        from ai._llm_client import _build_conversation

        with pytest.raises(ValueError, match="claude"):
            _build_conversation(
                provider="cohere",
                model="command-r",
                config={},
            )

    def test_unknown_provider_error_lists_valid_options(self):
        from ai._llm_client import _build_conversation

        with pytest.raises(ValueError) as exc_info:
            _build_conversation(provider="magic_llm", model="x", config={})

        msg = str(exc_info.value)
        assert "claude" in msg
        assert "ollama" in msg
        assert "openai" in msg

    def test_call_llm_writes_audit_record(self, monkeypatch, tmp_path):
        """call_llm must write an audit record even for mock calls."""
        self._mock_skill_loader(monkeypatch)
        self._mock_skill_version(monkeypatch)

        mock_conv = _mock_conversation()
        monkeypatch.setattr(
            "ai._llm_client._build_conversation",
            lambda provider, model, config: mock_conv,
        )

        from ai._llm_client import call_llm

        cfg = _make_config(provider="ollama", model="llama3")
        call_llm(
            "cluster_annotator",
            {"tissue": "liver", "cluster_id": "3"},
            cfg,
            log_dir=tmp_path,
            module="cluster_annotator",
        )

        log_file = tmp_path / "cluster_annotator.jsonl"
        assert log_file.exists()
        record = json.loads(log_file.read_text().strip())
        assert record["module"] == "cluster_annotator"
        assert record["provider"] == "ollama"

    def test_call_llm_returns_raw_string(self, monkeypatch, tmp_path):
        self._mock_skill_loader(monkeypatch)
        self._mock_skill_version(monkeypatch)

        expected = '{"cell_type": "NK cell", "confidence": "high"}'
        mock_conv = _mock_conversation(response=expected)
        monkeypatch.setattr(
            "ai._llm_client._build_conversation",
            lambda provider, model, config: mock_conv,
        )

        from ai._llm_client import call_llm

        cfg = _make_config(provider="ollama", model="llama3")
        result = call_llm("cluster_annotator", {"tissue": "liver"}, cfg,
                          log_dir=tmp_path)

        assert result == expected
