# AI/LLM tests require a running Ollama instance locally.
# Excluded from the standard test suite until the AI layer is reactivated (Phase 3).
collect_ignore = ["tests/test_groundedness.py"]
