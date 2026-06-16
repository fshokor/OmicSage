"""
OmicSage UI — Page 3: Run
==========================
Fixes:
  - Live log streaming via background reader thread + queue
    (replaces select() which is unreliable on WSL2 subprocess pipes)
  - Page no longer flickers — st.rerun() only called when new lines arrive
    or process finishes
"""
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

import streamlit as st

from ui.state import *
from ui.defaults import MODALITY_STEPS, STEP_LABELS, MODALITY_RUNNER

# Queue key in session state — holds lines read by the background thread
_QUEUE_KEY = "_log_queue"


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

        col1, col2 = st.columns([2, 1])
        with col1:
            run_mode = st.radio(
                "What to run",
                ["All configured steps", "Single step", "Step range (from → to)"],
                horizontal=True,
                key="run_mode_radio",
            )
        with col2:
            force = st.toggle(
                "⚡ Force re-run",
                value=False,
                help="Re-runs steps even if checkpoint .h5ad already exists.",
                key="force_toggle",
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

        # Preview
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
                _launch(modality, config_path, run_mode,
                        step_single, step_from, step_to, force)
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

        # Drain the queue populated by the background reader thread
        q   = st.session_state.get(_QUEUE_KEY)
        proc = st.session_state.get(KEY_RUN_PROCESS)
        got_new = False

        if q is not None:
            try:
                while True:
                    line = q.get_nowait()
                    st.session_state[KEY_RUN_LOG].append(line)
                    got_new = True
            except queue.Empty:
                pass

        _render_log(log_box)

        # Check if process finished
        if proc is not None and proc.poll() is not None:
            # Drain any remaining items in the queue
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
            # Still running — poll every 0.5s
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


# ── Launch + background reader ────────────────────────────────────────────────

def _launch(modality, config_path, run_mode, step_single, step_from, step_to, force):
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

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,                 # line-buffered
        universal_newlines=True,
        encoding="utf-8",
        errors="replace",
    )

    # Background thread reads lines and puts them in a queue.
    # This avoids blocking the Streamlit event loop and makes
    # lines available immediately as they arrive (no select() needed).
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
    for line in lines[-400:]:
        if any(line.startswith(f"[{s}]") for s in step_prefixes):
            display.append(f"  {line}")
        elif any(w in line for w in ("ERROR","Traceback","Error","KeyError","ValueError")):
            display.append(f"✗ {line}")
        elif any(w in line for w in ("complete","Pipeline complete","✓")):
            display.append(f"✓ {line}")
        elif line.startswith("$"):
            display.append(line)
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
