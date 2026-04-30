"""OmicSage — Phase 0 Smoke Tests"""
import yaml
import pytest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent


class TestRepoStructure:
    REQUIRED_DIRS = [
        ".dev_memory", "config",
        "pipeline/modules/qc", "pipeline/modules/processing",
        "pipeline/modules/clustering", "pipeline/modules/annotation",
        "pipeline/modules/downstream", "pipeline/workflows",
        "ai", "reports/templates", "reports/slides",
        "ui", "cli", "data/benchmark", "data/references",
        "tests", "docs", "docker", ".github/workflows",
    ]
    REQUIRED_FILES = [
        "config/schema.yaml", "pipeline/main.nf",
        "pipeline/workflows/scrna.nf", "pipeline/workflows/scatac.nf",
        "pipeline/workflows/spatial.nf", "pipeline/workflows/integration.nf",
        "nextflow.config", "docker/Dockerfile.python", "docker/Dockerfile.r",
        "docker-compose.yml", "environment.yml",
        ".github/workflows/ci.yml", ".github/workflows/docker-publish.yml",
        "ai/biochatter_client.py", "cli/omicsage.py", "ui/app.py",
        "README.md", ".gitignore",
    ]

    @pytest.mark.parametrize("dir_path", REQUIRED_DIRS)
    def test_directory_exists(self, dir_path):
        assert (REPO_ROOT / dir_path).is_dir(), f"Missing directory: {dir_path}"

    @pytest.mark.parametrize("file_path", REQUIRED_FILES)
    def test_file_exists(self, file_path):
        assert (REPO_ROOT / file_path).is_file(), f"Missing file: {file_path}"


class TestConfigSchema:
    @pytest.fixture
    def schema(self):
        with open(REPO_ROOT / "config" / "schema.yaml") as f:
            return yaml.safe_load(f)

    def test_schema_loads(self, schema):
        assert schema is not None

    def test_required_top_level_keys(self, schema):
        for key in ["project", "modality", "input", "qc", "processing",
                    "clustering", "annotation", "downstream", "reports", "ai", "compute"]:
            assert key in schema, f"Missing key: {key}"

    def test_ai_enabled_is_bool(self, schema):
        assert isinstance(schema["ai"]["enabled"], bool)

    def test_modality_flags_are_bool(self, schema):
        for flag in ["scrna", "scatac", "spatial", "multiome"]:
            assert isinstance(schema["modality"][flag], bool)


class TestDevMemory:
    MEMORY_FILES = [
        ".dev_memory/NEXT_SESSION.md", ".dev_memory/CURRENT_STATUS.md",
        ".dev_memory/DECISIONS.md", ".dev_memory/PROGRESS.md",
        ".dev_memory/BLOCKERS.md",
    ]

    @pytest.mark.parametrize("file_path", MEMORY_FILES)
    def test_memory_file_exists(self, file_path):
        assert (REPO_ROOT / file_path).is_file(), f"Missing: {file_path}"

    @pytest.mark.parametrize("file_path", MEMORY_FILES)
    def test_memory_file_not_empty(self, file_path):
        path = REPO_ROOT / file_path
        if path.exists():
            assert path.stat().st_size > 50, f"{file_path} is empty"


class TestCLI:
    def test_ai_client_imports(self):
        import sys
        sys.path.insert(0, str(REPO_ROOT))
        from ai.biochatter_client import LLMConfig, OmicSageAIClient
        cfg = LLMConfig()
        assert cfg.provider == "claude"
        client = OmicSageAIClient(cfg)
        assert client is not None
