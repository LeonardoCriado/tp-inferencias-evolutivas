#!/usr/bin/env python3
"""Analisis de mosaico de la secuencia 'bartmosca'.

Ventanas deslizantes de 15 aa con paso 5 sobre el alineamiento Clustal Omega del
set9 (9 secuencias). Para cada ventana se calcula la identidad de bartmosca
contra humano y contra mosca con delecion de pares (solo columnas donde ambos
tienen residuo). Ademas, identidad por regiones: N-terminal (primeros 25
residuos de bartmosca) vs resto de la proteina.

Salidas:
  resultados/bartmosca_ventanas.csv
  resultados/bartmosca_regiones.csv
  resultados/fig_ventanas_bartmosca.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import AlignIO

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "resultados")
FIG = RES  # las figuras van dentro de resultados/
os.makedirs(FIG, exist_ok=True)

ALN = os.path.join(RES, "aln_clustalo_set9.fasta")
WIN, STEP, NTERM = 15, 5, 25


def ident(a, b):
    both = match = 0
    for x, y in zip(a, b):
        if x != "-" and y != "-":
            both += 1
            if x == y:
                match += 1
    return (match / both if both else float("nan")), both


def main():
    aln = AlignIO.read(ALN, "fasta")
    seqs = {r.id: str(r.seq).upper() for r in aln}
    L = aln.get_alignment_length()
    bart, hum, mos = seqs["bartmosca"], seqs["humano"], seqs["mosca"]

    # posicion en el alineamiento del residuo 25 de bartmosca (regiones)
    seen = 0
    cut = L
    for j, c in enumerate(bart, start=1):
        if c != "-":
            seen += 1
            if seen == NTERM:
                cut = j
                break

    rows = []
    for start in range(0, max(1, L - WIN + 1), STEP):
        end = start + WIN
        bh, nb = ident(bart[start:end], hum[start:end])
        bm, nm = ident(bart[start:end], mos[start:end])
        rows.append(dict(ventana=f"{start+1}-{end}", inicio=start + 1, fin=end,
                         pos_centro=(start + 1 + end) / 2,
                         n_posiciones_ambas_hum=nb, n_posiciones_ambas_mos=nm,
                         ident_humano=round(bh, 4), ident_mosca=round(bm, 4),
                         delta_humano_menos_mosca=round(bh - bm, 4),
                         mas_parecida_a=("humano" if bh > bm else "mosca")))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "bartmosca_ventanas.csv"), index=False)
    print(df.to_string(index=False))

    # regiones
    regions = []
    for name, (s, e) in (("N-terminal (1-25)", (0, cut)), ("resto (26-fin)", (cut, L))):
        bh, nb = ident(bart[s:e], hum[s:e])
        bm, nm = ident(bart[s:e], mos[s:e])
        regions.append(dict(region=name, columnas_alineamiento=f"{s+1}-{e}",
                            n_posiciones=nb,
                            ident_humano=round(bh, 4), ident_mosca=round(bm, 4),
                            delta_humano_menos_mosca=round(bh - bm, 4)))
    rdf = pd.DataFrame(regions)
    rdf.to_csv(os.path.join(RES, "bartmosca_regiones.csv"), index=False)
    print()
    print(rdf.to_string(index=False))
    print(f"\ncorte N-term en columna {cut} del alineamiento (25 residuos de bartmosca)")

    # ---- figura ----
    fig, ax = plt.subplots(figsize=(11, 5.2), dpi=150)
    ax.axvspan(1, cut, color="#d9e8f5", alpha=0.9, zorder=0,
               label=f"Region N-terminal de bartmosca (residuos 1-{NTERM})")
    ax.plot(df["pos_centro"], df["ident_humano"], "o-", color="#1f4e79", lw=2,
            ms=4, label="bartmosca vs humano")
    ax.plot(df["pos_centro"], df["ident_mosca"], "s-", color="#c0392b", lw=2,
            ms=4, label="bartmosca vs mosca")
    ax.fill_between(df["pos_centro"], df["ident_mosca"], df["ident_humano"],
                    where=df["ident_humano"] >= df["ident_mosca"],
                    color="#1f4e79", alpha=0.12, interpolate=True)
    ax.axhline(0.5, color="grey", ls=":", lw=1)
    ax.set_xlabel("posicion en el alineamiento (aa)", fontsize=11)
    ax.set_ylabel("identidad por ventana (delecion de pares)", fontsize=11)
    ax.set_title("Analisis de mosaico de bartmosca: ventanas de 15 aa (paso 5)\n"
                 "identidad frente a humano y mosca (alineamiento Clustal Omega, set9)",
                 fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(1, L)
    ax.grid(alpha=0.25)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.95)
    # anotacion: region donde gana humano
    ax.text(0.02, 0.03, "azul: bartmosca mas parecida al humano\n"
                        "rojo: bartmosca mas parecida a la mosca\n"
                        "las 2 primeras ventanas no tienen posiciones comparables (gaps)",
            transform=ax.transAxes, fontsize=8, va="bottom",
            bbox=dict(boxstyle="round", fc="white", alpha=0.85, ec="grey"))
    fig.tight_layout()
    out = os.path.join(FIG, "fig_ventanas_bartmosca.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("figura:", out)


if __name__ == "__main__":
    main()
