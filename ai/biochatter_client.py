"""
OmicSage — BioChatter AI Client
Phase 3 implementation. Currently a typed stub.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class LLMConfig:
    provider: str = "claude"
    model: str = "claude-opus-4-5"
    api_key_env: str = "ANTHROPIC_API_KEY"
    ollama_host: str = "http://localhost:11434"
    ollama_model: str = "llama3.2"
    log_calls: bool = True
    log_dir: str = "logs/llm/"


class OmicSageAIClient:
    """
    Thin wrapper around BioChatter for OmicSage-specific prompting.
    Phase 3 implementation target.
    """

    def __init__(self, config: LLMConfig):
        self.config = config
        self._client = None
        self._log_dir = Path(config.log_dir)
        self._log_dir.mkdir(parents=True, exist_ok=True)

    def suggest_qc_thresholds(self, qc_metrics: dict) -> dict:
        logger.warning("suggest_qc_thresholds: Phase 3 — not yet implemented")
        return {}

    def interpret_cluster(self, cluster_id: int, marker_genes: list[str], n_cells: int) -> dict:
        logger.warning("interpret_cluster: Phase 3 — not yet implemented")
        return {}

    def generate_narrative(self, analysis_summary: dict) -> str:
        logger.warning("generate_narrative: Phase 3 — not yet implemented")
        return ""

    def _log_call(self, method: str, prompt: str, response: str, model: str):
        if not self.config.log_calls:
            return
        entry = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "method": method,
            "model": model,
            "prompt_chars": len(prompt),
            "response_chars": len(response),
            "prompt": prompt,
            "response": response,
        }
        log_file = self._log_dir / f"llm_calls_{datetime.now().strftime('%Y%m%d')}.jsonl"
        with open(log_file, "a") as f:
            f.write(json.dumps(entry) + "\n")
