"""Summarize the completed Kuppe P1 execution pilot without loading expression."""
from pathlib import Path
import json
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def summarize():
    root = Path(__file__).resolve().parents[1]
    directory = root / "reports/kuppe_heart_P1_pilot"
    checkpoint = root / "data/processed/kuppe_heart_P1_pilot/05_deconvolved.h5ad"
    with h5py.File(checkpoint) as data:
        meta = data["uns/omicsage_spatial_deconvolve/outputs"]
        names = meta["cell_type_names"].asstr()[()].tolist()
        abundance = data["obsm/cell_type_abundances"][()]
        images = list(data["uns/spatial"])
        genes = int(meta["n_shared_genes"][()])
        finite = bool(np.isfinite(abundance).all())
        assert images == ["control_P1"], images
        assert abundance.shape == (4279, 11), abundance.shape
        assert finite and (abundance >= 0).all()
        histories = {}
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        for ax, stage in zip(axes, ("reference", "spatial")):
            group = data[f"uns/cell2location_{stage}_history"]
            for key, dataset in group.items():
                values = dataset[()]
                if len(values) and np.isfinite(values).all():
                    histories[f"{stage}/{key}"] = {
                        "n_epochs": len(values), "first": float(values[0]), "last": float(values[-1])}
                    if "train" in key:
                        ax.plot(np.arange(1, len(values)+1), values, label=key)
            ax.set(title=f"{stage.capitalize()} pilot loss", xlabel="Epoch", ylabel="Loss")
            ax.legend()
        fig.suptitle("control_P1 CPU pilot — convergence not established")
        fig.tight_layout()
        fig.savefig(directory / "pilot_training_history.png", dpi=150)
        plt.close(fig)
        summary = {"sample": "control_P1", "n_spots": abundance.shape[0],
                   "image_libraries": images, "n_cell_types": len(names),
                   "n_shared_genes": genes, "finite_nonnegative_abundance": finite,
                   "mean_q05_abundance": dict(zip(names, abundance.mean(0).astype(float))),
                   "zero_weight_types": [name for i, name in enumerate(names) if not np.any(abundance[:, i] > 0)],
                   "training_history": histories, "fit_validated": False}
    resources = root / "logs/kuppe_heart_P1_pilot_resources.txt"
    if resources.exists():
        summary["resource_measurement"] = resources.read_text()
    (directory / "pilot_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    summarize()
