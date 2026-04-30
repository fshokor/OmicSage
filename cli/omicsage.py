#!/usr/bin/env python3
"""OmicSage CLI — Project management for bioinformaticians."""
import click
import json
import shutil
from pathlib import Path
from datetime import datetime

OMICSAGE_VERSION = "0.1.0-dev"


@click.group()
@click.version_option(version=OMICSAGE_VERSION, prog_name="OmicSage")
def cli():
    """OmicSage — AI-Assisted Single-Cell Multi-Omics Platform."""
    pass


@cli.command()
@click.argument("project_name")
@click.option("--modality", "-m",
              type=click.Choice(["scrna", "scatac", "spatial", "multiome"]),
              default="scrna", show_default=True)
@click.option("--outdir", "-o", default=None)
def create_project(project_name: str, modality: str, outdir: str):
    """Create a new OmicSage analysis project."""
    project_dir = Path(outdir or project_name)
    if project_dir.exists():
        click.echo(f"Error: directory already exists: {project_dir}", err=True)
        raise click.Abort()

    for subdir in ["data/raw", "data/processed", "results", "logs/llm", "reports"]:
        (project_dir / subdir).mkdir(parents=True)

    config_src = Path(__file__).parent.parent / "config" / "schema.yaml"
    config_dst = project_dir / "config.yaml"
    if config_src.exists():
        shutil.copy(config_src, config_dst)

    meta = {
        "project_name": project_name,
        "modality": modality,
        "created": datetime.now().isoformat(),
        "omicsage_version": OMICSAGE_VERSION,
        "status": "created"
    }
    (project_dir / "project.json").write_text(json.dumps(meta, indent=2))

    click.echo(f"Project created: {project_dir}/")
    click.echo(f"  Modality : {modality}")
    click.echo(f"  Config   : {config_dst}")
    click.echo(f"\nNext: edit {config_dst}, then: omicsage run {project_dir}/")


@cli.command()
@click.argument("project_dir", type=click.Path(exists=True))
@click.option("--profile", "-p",
              type=click.Choice(["local", "docker", "singularity", "slurm"]),
              default="local", show_default=True)
@click.option("--dry-run", is_flag=True)
def run(project_dir: str, profile: str, dry_run: bool):
    """Run the OmicSage pipeline on a project directory."""
    import subprocess
    project_path = Path(project_dir)
    meta_file = project_path / "project.json"
    meta = json.loads(meta_file.read_text()) if meta_file.exists() else {}
    modality = meta.get("modality", "scrna")
    nf_main = Path(__file__).parent.parent / "pipeline" / "main.nf"
    cmd = [
        "nextflow", "run", str(nf_main),
        "-profile", profile,
        "--config", str(project_path / "config.yaml"),
        "--modality", modality,
        "--outdir", str(project_path / "results"),
    ]
    if dry_run:
        click.echo("DRY RUN: " + " ".join(cmd))
        return
    click.echo(f"Running OmicSage [{modality}] on {project_dir}")
    subprocess.run(cmd, check=True)


@cli.command("list")
@click.argument("search_dir", type=click.Path(exists=True), default=".")
def list_projects(search_dir: str):
    """List all OmicSage projects in a directory."""
    projects = list(Path(search_dir).rglob("project.json"))
    if not projects:
        click.echo("No OmicSage projects found.")
        return
    click.echo(f"{'Project':<30} {'Modality':<12} {'Status':<12} {'Created':<20}")
    click.echo("-" * 74)
    for p in sorted(projects):
        meta = json.loads(p.read_text())
        click.echo(
            f"{meta.get('project_name','?'):<30} "
            f"{meta.get('modality','?'):<12} "
            f"{meta.get('status','?'):<12} "
            f"{meta.get('created','?')[:19]:<20}"
        )


@cli.command()
@click.argument("project_dirs", nargs=-1, type=click.Path(exists=True), required=True)
@click.option("--output", "-o", default="comparison_report.html")
def compare(project_dirs: tuple, output: str):
    """Compare results across multiple OmicSage projects."""
    click.echo(f"Comparing {len(project_dirs)} projects -> {output}")
    click.echo("compare: Phase 2 feature — not yet implemented")


if __name__ == "__main__":
    cli()
