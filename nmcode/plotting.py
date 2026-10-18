from pathlib import Path
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


PALETTE = ["#4C78A8", "#72B7B2", "#E45756", "#F2CF5B", "#B279A2", "#FF9D4D"]


def set_style():
    sns.set_theme(context="paper", style="ticks", font="Arial", palette=PALETTE)
    mpl.rcParams.update({
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "font.family": "Arial",
        "font.size": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.transparent": False,
    })


def finish_figure(fig, path, width=7.2, height=4.8):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(width, height)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*tight_layout.*")
        warnings.filterwarnings("ignore", message="Tight layout not applied.*")
        fig.tight_layout()
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def heatmap(ax, matrix, xlabels=None, ylabels=None, cmap="RdBu_r", center=None):
    values = np.asarray(matrix, dtype=float)
    limit = np.nanmax(np.abs(values)) if center == 0 else None
    image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=-limit if limit else None, vmax=limit)
    if xlabels is not None:
        ax.set_xticks(range(len(xlabels)), xlabels, rotation=45, ha="right")
    if ylabels is not None:
        ax.set_yticks(range(len(ylabels)), ylabels)
    return image
