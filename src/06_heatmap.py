#!/usr/bin/env python3
"""Heatmap de identidad de pares (delecion de pares) para el set9.

Usa resultados/pairwise_identidad_clustalo_set9.csv (calculado en 03_pairwise.py).
Salida: resultados/fig_heatmap_identidad.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "resultados")
FIG = RES  # las figuras van dentro de resultados/

ORDER = ["humano", "bartmosca", "mosca", "cebada", "altramuz", "abedul", "vigna",
         "ectocarpus", "rhizobium"]


def main():
    df = pd.read_csv(os.path.join(RES, "pairwise_identidad_clustalo_set9.csv"), index_col=0)
    df = df.reindex(index=ORDER, columns=ORDER)
    M = df.to_numpy(dtype=float, na_value=1.0)
    M = np.where(np.isnan(M), 1.0, M)

    fig, ax = plt.subplots(figsize=(8.2, 7.0), dpi=150)
    im = ax.imshow(M, cmap="YlGnBu", vmin=0.45, vmax=1.0)
    ax.set_xticks(range(len(ORDER)), ORDER, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(ORDER)), ORDER, fontsize=10)
    for i in range(len(ORDER)):
        for j in range(len(ORDER)):
            v = M[i, j]
            ax.text(j, i, f"{v*100:.1f}", ha="center", va="center", fontsize=8,
                    color="white" if v > 0.85 else "black")
    ax.set_title("Identidad de pares de citocromo c (set9, %)\n"
                 "delecion de pares, alineamiento Clustal Omega, 131 columnas", fontsize=11)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("identidad (%)", fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    out = os.path.join(FIG, "fig_heatmap_identidad.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("figura:", out)


if __name__ == "__main__":
    main()
