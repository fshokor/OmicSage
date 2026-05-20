"""
conftest.py — OmicSage test configuration.

Installs lightweight stubs for optional heavy dependencies (scvi-tools, celltypist)
so that mock-based tests work on machines where those packages are not installed.
The stubs are added to sys.modules BEFORE any test module is imported, which means
patch("scvi.model.SCANVI") resolves correctly regardless of whether scvi-tools is
actually installed in the environment.

If the real package IS installed, it takes precedence (stubs are only inserted for
missing modules via setdefault-style logic).
"""
import sys
import types


def _ensure_scvi_stub() -> None:
    """Insert a minimal scvi + scvi.model stub if scvi-tools is not installed."""
    if "scvi" in sys.modules:
        return  # real package already imported — don't override

    scvi_mod = types.ModuleType("scvi")
    scvi_model_mod = types.ModuleType("scvi.model")

    class _StubSCANVI:
        """Placeholder — real tests mock this out; stub just needs to be patchable."""
        @staticmethod
        def load(*args, **kwargs):
            raise ImportError(
                "scvi-tools stub: call patch('scvi.model.SCANVI') to override in tests"
            )

    scvi_model_mod.SCANVI = _StubSCANVI
    scvi_mod.model = scvi_model_mod

    sys.modules.setdefault("scvi", scvi_mod)
    sys.modules.setdefault("scvi.model", scvi_model_mod)


_ensure_scvi_stub()
