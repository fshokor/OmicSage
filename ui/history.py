"""
OmicSage UI — Project History
==============================
Stores and retrieves run history in .omicsage_history.json at repo root.
Each entry: dataset_name, modality, config_path, reports_dir, timestamp, status.
Max 50 entries kept (oldest dropped).
"""
import json
import os
from datetime import datetime
from pathlib import Path

HISTORY_FILE = Path(".omicsage_history.json")
MAX_ENTRIES  = 50


def _load() -> list:
    if HISTORY_FILE.exists():
        try:
            return json.loads(HISTORY_FILE.read_text(encoding="utf-8"))
        except Exception:
            return []
    return []


def _save(entries: list) -> None:
    HISTORY_FILE.write_text(
        json.dumps(entries, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def add_entry(
    dataset_name: str,
    dataset_id: str,
    modality: str,
    config_path: str,
    reports_dir: str,
    status: str = "running",
) -> str:
    """
    Add a new history entry. Returns the entry ID (timestamp string).
    """
    entries = _load()
    entry_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    entry = {
        "id":           entry_id,
        "timestamp":    datetime.now().isoformat(timespec="seconds"),
        "dataset_name": dataset_name,
        "dataset_id":   dataset_id,
        "modality":     modality,
        "config_path":  config_path,
        "reports_dir":  reports_dir,
        "status":       status,
    }
    entries.insert(0, entry)          # newest first
    entries = entries[:MAX_ENTRIES]   # trim
    _save(entries)
    return entry_id


def update_status(entry_id: str, status: str) -> None:
    """Update the status of an existing entry by ID."""
    entries = _load()
    for e in entries:
        if e.get("id") == entry_id:
            e["status"] = status
            break
    _save(entries)


def get_recent(n: int = 10) -> list:
    """Return the n most recent history entries."""
    return _load()[:n]


def load_config_from_history(entry: dict) -> dict | None:
    """
    Try to load and return the YAML config dict for a history entry.
    Returns None if the file no longer exists.
    """
    import yaml
    config_path = entry.get("config_path", "")
    if config_path and Path(config_path).exists():
        try:
            return yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
        except Exception:
            return None
    return None


STATUS_ICON = {
    "running": "⚙️",
    "done":    "✅",
    "error":   "❌",
    "idle":    "⏸",
}


def status_icon(status: str) -> str:
    return STATUS_ICON.get(status, "⏸")
