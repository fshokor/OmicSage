"""
OmicSage UI — Page 3: Run
==========================
Supports two run engines:
  - Python  : calls run_<modality>_pipeline.py directly (no Docker needed)
  - Nextflow: calls `nextflow run main.nf` (Docker + reproducibility)

Live log streaming via background reader thread + queue.
"""
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

import streamlit as st

from ui.state import *
from ui import history as _history
from ui.defaults import MODALITY_STEPS, STEP_LABELS, MODALITY_RUNNER, MODALITY_NF_NAME, NF_IMPLEMENTED

# Queue key in session state — holds lines read by the background thread
_QUEUE_KEY = "_log_queue"

# Session state key for selected engine
_ENGINE_KEY = "_run_engine"


def render():
    init_state()
    modality     = st.session_state[KEY_MODALITY]
    dataset_name = st.session_state[KEY_DATASET_NAME]
    config_path  = st.session_state.get(KEY_CONFIG_PATH)
    selected     = st.session_state.get(KEY_SELECTED_STEPS, [])
    status       = st.session_state[KEY_RUN_STATUS]

    st.markdown(f"## Run  ·  *{dataset_name}*  ·  `{modality}`")

    if not config_path:
        st.warning("No config saved — complete the Configure step first.")
        if st.button("← Back to Configure"):
            st.session_state[KEY_PAGE] = 1
            st.rerun()
        return

    all_steps     = MODALITY_STEPS[modality]
    enabled_steps = [s for s in all_steps if s in selected]

    # ── Idle / post-run controls ──────────────────────────────────────────────
    if status in ("idle", "done", "error"):
        st.divider()
        st.markdown("### Run options")

        # ── Engine selector ───────────────────────────────────────────────────
        nf_available = modality in NF_IMPLEMENTED
        engine_options = ["Python (direct)", "Nextflow (Docker)"]

        if not nf_available:
            engine_help = f"Nextflow not yet implemented for {modality} — using Python."
            engine_index = 0
            engine_disabled = True
        else:
            engine_help = (
                "**Python**: runs locally, no Docker needed, faster for development.\n\n"
                "**Nextflow**: runs in Docker, fully reproducible, supports -resume, "
                "HPC/cloud ready."
            )
            engine_disabled = False
            engine_index = engine_options.index(
                st.session_state.get(_ENGINE_KEY, "Python (direct)")
            )

        selected_engine = st.radio(
            "Run engine",
            engine_options,
            index=engine_index,
            horizontal=True,
            disabled=engine_disabled,
            help=engine_help,
            key="engine_radio",
        )
        st.session_state[_ENGINE_KEY] = selected_engine
        use_nextflow = selected_engine == "Nextflow (Docker)" and nf_available

        if use_nextflow:
            st.info(
                "⚙️ Nextflow mode — pipeline runs inside `omicsage:latest` Docker container. "
                "Use **-resume** to skip completed steps automatically.",
                icon="🐋",
            )

        st.divider()

        col1, col2 = st.columns([2, 1])
        with col1:
            if use_nextflow:
                # Nextflow handles step selection differently — always runs full pipeline
                # but -resume skips completed steps automatically
                run_mode = st.radio(
                    "What to run",
                    ["Full pipeline (-resume)", "Full pipeline (force all steps)"],
                    horizontal=True,
                    key="run_mode_radio",
                )
                step_single = step_from = step_to = None
            else:
                run_mode = st.radio(
                    "What to run",
                    ["All configured steps", "Single step", "Step range (from → to)"],
                    horizontal=True,
                    key="run_mode_radio",
                )
                step_single = step_from = step_to = None

                if run_mode == "Single step":
                    step_single = st.selectbox(
                        "Step to run", options=enabled_steps,
                        format_func=lambda s: f"{s}  —  {STEP_LABELS.get(s,'')}",
                        key="single_step_sel",
                    )
                elif run_mode == "Step range (from → to)":
                    ca, cb = st.columns(2)
                    with ca:
                        step_from = st.selectbox(
                            "From step (inclusive)", options=all_steps, index=0,
                            format_func=lambda s: f"{s}  —  {STEP_LABELS.get(s,'')}",
                            key="from_step_sel",
                        )
                    with cb:
                        from_idx  = all_steps.index(step_from) if step_from else 0
                        remaining = all_steps[from_idx:]
                        step_to   = st.selectbox(
                            "To step (inclusive)", options=remaining,
                            index=len(remaining) - 1,
                            format_func=lambda s: f"{s}  —  {STEP_LABELS.get(s,'')}",
                            key="to_step_sel",
                        )

        with col2:
            if not use_nextflow:
                force = st.toggle(
                    "⚡ Force re-run",
                    value=False,
                    help="Re-runs steps even if checkpoint .h5ad already exists.",
                    key="force_toggle",
                )
            else:
                force = False   # Nextflow resume handles this via run_mode

        # Preview
        if use_nextflow:
            preview = enabled_steps   # show what will run; Nextflow skips cached
            if run_mode == "Full pipeline (-resume)":
                st.markdown("**Will run:**  " + "  →  ".join(f"`{s}`" for s in preview))
                st.caption("Completed steps will be skipped automatically via Nextflow -resume.")
            else:
                st.markdown("**Will run (force all):**  " + "  →  ".join(f"`{s}`" for s in preview))
                st.warning("⚡ Force mode — all steps will re-run even if checkpoints exist.")
        else:
            if run_mode == "All configured steps":
                preview = enabled_steps
            elif run_mode == "Single step":
                preview = [step_single] if step_single else []
            else:
                if step_from and step_to:
                    fi = all_steps.index(step_from)
                    ti = all_steps.index(step_to)
                    preview = all_steps[fi:ti+1]
                else:
                    preview = []

            if preview:
                st.markdown("**Will run:**  " + "  →  ".join(f"`{s}`" for s in preview))
            if force:
                st.warning("⚡ Force mode on — existing checkpoints will be overwritten.")

        st.divider()
        col_back, col_run = st.columns([1, 3])
        with col_back:
            if st.button("← Back to Configure", use_container_width=True):
                st.session_state[KEY_PAGE] = 1
                st.rerun()
        with col_run:
            label = "▶ Run Pipeline" if status == "idle" else "▶ Run Again"
            if st.button(label, type="primary", use_container_width=True,
                         disabled=not preview):
                _launch(
                    modality, config_path, run_mode,
                    step_single, step_from, step_to, force,
                    use_nextflow=use_nextflow,
                )
                st.rerun()

    # ── Running ───────────────────────────────────────────────────────────────
    if status == "running":
        st.info("⚙️  Pipeline running…")

        if st.button("⏹ Stop", type="secondary", key="stop_btn"):
            proc = st.session_state.get(KEY_RUN_PROCESS)
            if proc:
                proc.terminate()
            st.session_state[KEY_RUN_STATUS]  = "error"
            st.session_state[KEY_RUN_PROCESS] = None
            st.session_state[KEY_RUN_LOG].append("— Stopped by user —")
            st.rerun()

        st.markdown("**Live output**")
        log_box = st.empty()

        q    = st.session_state.get(_QUEUE_KEY)
        proc = st.session_state.get(KEY_RUN_PROCESS)

        if q is not None:
            try:
                while True:
                    line = q.get_nowait()
                    st.session_state[KEY_RUN_LOG].append(line)
            except queue.Empty:
                pass

        _render_log(log_box)

        if proc is not None and proc.poll() is not None:
            if q:
                try:
                    while True:
                        line = q.get_nowait()
                        st.session_state[KEY_RUN_LOG].append(line)
                except queue.Empty:
                    pass
            _render_log(log_box)
            st.session_state[KEY_RUN_STATUS]  = "done" if proc.returncode == 0 else "error"
            st.session_state[KEY_RUN_PROCESS] = None
            st.session_state[_QUEUE_KEY]      = None
            st.rerun()
        else:
            time.sleep(0.5)
            st.rerun()

    # ── Done ──────────────────────────────────────────────────────────────────
    if status == "done":
        st.success("✅  Pipeline complete!")
        _show_done_actions()
        st.divider()
        st.markdown("**Full log**")
        _render_log(st.empty())

    # ── Error ─────────────────────────────────────────────────────────────────
    if status == "error":
        st.error("❌  Pipeline exited with an error — check the log below.")
        if st.button("↩ Try again", key="retry_btn"):
            reset_run()
            st.rerun()
        st.divider()
        _render_log(st.empty())


# ── Launch ────────────────────────────────────────────────────────────────────

def _launch(modality, config_path, run_mode, step_single, step_from, step_to,
            force, use_nextflow=False):

    if use_nextflow:
        cmd = _build_nextflow_cmd(modality, config_path, run_mode)
    else:
        cmd = _build_python_cmd(modality, config_path, run_mode,
                                step_single, step_from, step_to, force)

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        universal_newlines=True,
        encoding="utf-8",
        errors="replace",
    )

    log_queue: queue.Queue = queue.Queue()

    def _reader(p, q):
        try:
            for line in iter(p.stdout.readline, ""):
                q.put(line.rstrip("\n"))
        except Exception:
            pass
        finally:
            try:
                p.stdout.close()
            except Exception:
                pass

    t = threading.Thread(target=_reader, args=(proc, log_queue), daemon=True)
    t.start()

    st.session_state[KEY_RUN_PROCESS] = proc
    st.session_state[KEY_RUN_STATUS]  = "running"
    st.session_state[KEY_RUN_LOG]     = [f"$ {' '.join(cmd)}"]
    st.session_state[_QUEUE_KEY]      = log_queue

    entry_id = _history.add_entry(
        dataset_name = st.session_state.get(KEY_DATASET_NAME, ""),
        dataset_id   = st.session_state.get(KEY_DATASET_ID, ""),
        modality     = st.session_state.get(KEY_MODALITY, ""),
        config_path  = config_path,
        reports_dir  = st.session_state.get(KEY_REPORTS_DIR, ""),
        status       = "running",
    )
    st.session_state["_history_entry_id"] = entry_id


def _build_python_cmd(modality, config_path, run_mode,
                      step_single, step_from, step_to, force):
    """Build command for direct Python runner."""
    runner = MODALITY_RUNNER[modality]
    cmd    = [sys.executable, runner, "--config", config_path]

    if run_mode == "Single step" and step_single:
        cmd += ["--step", step_single]
    elif run_mode == "Step range (from → to)":
        if step_from:
            cmd += ["--from-step", step_from]
        if step_to:
            cmd += ["--to-step", step_to]

    if force:
        cmd.append("--force")

    return cmd


def _build_nextflow_cmd(modality, config_path, run_mode):
    """Build command for Nextflow engine."""
    nf_modality = MODALITY_NF_NAME[modality]
    cmd = [
        "nextflow", "run", "main.nf",
        "--config",   config_path,
        "--modality", nf_modality,
        "-profile",   "local",
        "-ansi-log",  "false",   # plain text output (no ANSI colours in log box)
    ]

    # -resume skips steps whose checkpoints already exist
    if run_mode == "Full pipeline (-resume)":
        cmd.append("-resume")

    return cmd


# ── Log renderer ──────────────────────────────────────────────────────────────

def _render_log(placeholder):
    lines = st.session_state.get(KEY_RUN_LOG, [])
    if not lines:
        return
    display = []
    step_prefixes = [
        "qc","normalize","reduce","cluster","annotate","deg","gsea",
        "harmony","pseudobulk","atac","multiome","ingest","deconv",
        "downstream","impute","ai","data_report",
    ]
    # Nextflow-specific line markers
    nf_markers = ["executor >", "process >", "SCRNA_", "N E X T F L O W", "["]

    for line in lines[-400:]:
        if any(line.startswith(f"[{s}]") for s in step_prefixes):
            display.append(f"  {line}")
        elif any(w in line for w in ("ERROR","Traceback","Error","KeyError","ValueError")):
            display.append(f"✗ {line}")
        elif any(w in line for w in ("complete","Pipeline complete","✓")):
            display.append(f"✓ {line}")
        elif line.startswith("$"):
            display.append(line)
        elif any(m in line for m in nf_markers):
            display.append(f"  {line}")
        else:
            display.append(line)

    with placeholder:
        st.code("\n".join(display), language=None)


# ── Done actions ──────────────────────────────────────────────────────────────

def _show_done_actions():
    reports_dir  = st.session_state.get(KEY_REPORTS_DIR, "")
    combined     = Path(reports_dir) / "00_combined_report.html"
    dataset_name = st.session_state.get(KEY_DATASET_NAME, "report")

    col1, col2, col3 = st.columns(3)
    with col1:
        if st.button("📊 View Report →", type="primary", use_container_width=True,
                     key="view_report_btn"):
            st.session_state[KEY_PAGE] = 3
            st.rerun()
    with col2:
        if st.button("↩ Re-run (change steps)", use_container_width=True,
                     key="rerun_btn"):
            reset_run()
            st.rerun()
    with col3:
        if combined.exists():
            with open(combined, "rb") as f:
                st.download_button(
                    "⬇ Download HTML report",
                    data=f.read(),
                    file_name=f"{dataset_name.replace(' ','_')}_report.html",
                    mime="text/html",
                    use_container_width=True,
                    key="dl_report_btn",
                )
