"""Mosaico de bartmosca a nivel de residuo (Desafio I).

Compara bartmosca contra humano y mosca sobre el alineamiento de Clustal Omega
del set3 y clasifica cada residuo de bartmosca en:

  - identico_a_ambos     : el mismo residuo en humano y en mosca
  - identico_a_humano    : solo coincide con humano
  - identico_a_mosca     : solo coincide con mosca
  - distinto_de_ambos    : no coincide con ninguno de los dos
  - gap_en_referencia    : humano o mosca tienen gap en esa columna

Entrada : resultados/aln_clustalo_set3.fasta
Salidas : resultados/mosaico_residuos_bartmosca.csv
          resultados/fig_mosaico_bartmosca.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import AlignIO
from matplotlib.patches import Rectangle

RES = "resultados"
aln = AlignIO.read(os.path.join(RES, "aln_clustalo_set3.fasta"), "fasta")
S = {r.id: str(r.seq) for r in aln}
hum, mos, bar = S["humano"], S["mosca"], S["bartmosca"]
GAPS = set("-.")

COLORES = {
    "identico_a_ambos": "#d9d9d9",
    "identico_a_humano": "#a8c8e8",
    "identico_a_mosca": "#f0b27a",
    "distinto_de_ambos": "#c0392b",
    "gap_en_referencia": "#7f8c8d",
}
ETIQUETAS = {
    "identico_a_ambos": "igual en humano y mosca",
    "identico_a_humano": "solo coincide con humano",
    "identico_a_mosca": "solo coincide con mosca",
    "distinto_de_ambos": "distinto de ambos",
    "gap_en_referencia": "gap en alguna referencia",
}

filas = []
pos = 0
for col, (h, m, b) in enumerate(zip(hum, mos, bar), start=1):
    if b in GAPS:
        continue
    pos += 1
    if h in GAPS or m in GAPS:
        clase = "gap_en_referencia"
    elif b == h and b == m:
        clase = "identico_a_ambos"
    elif b == h:
        clase = "identico_a_humano"
    elif b == m:
        clase = "identico_a_mosca"
    else:
        clase = "distinto_de_ambos"
    filas.append(dict(pos_bart=pos, col_alineamiento=col, bart=b, humano=h, mosca=m, clase=clase))

df = pd.DataFrame(filas)
df.to_csv(os.path.join(RES, "mosaico_residuos_bartmosca.csv"), index=False)

conteo = df["clase"].value_counts().reindex(list(COLORES), fill_value=0)
total = int(conteo.sum())
print(df["clase"].value_counts().to_string())
print("total residuos bartmosca:", total)
print("\nResiduos que no coinciden con el humano:")
print(df[~df["clase"].isin(["identico_a_ambos", "identico_a_humano"])].to_string(index=False))

n_term_col = df.loc[df["pos_bart"] <= 25, "col_alineamiento"].max()
n_term = df["pos_bart"] <= 25
ident_hum_nt = (df.loc[n_term, "clase"].isin(["identico_a_ambos", "identico_a_humano"])).mean()
ident_mos_nt = (df.loc[n_term, "clase"].isin(["identico_a_ambos", "identico_a_mosca"])).mean()
ident_hum_ct = (df.loc[~n_term, "clase"].isin(["identico_a_ambos", "identico_a_humano"])).mean()
ident_mos_ct = (df.loc[~n_term, "clase"].isin(["identico_a_ambos", "identico_a_mosca"])).mean()
print(f"\nN-terminal (1-25): ident. vs humano {ident_hum_nt:.3f} | vs mosca {ident_mos_nt:.3f}")
print(f"resto ({int(n_term.sum())+1}-{total}): ident. vs humano {ident_hum_ct:.3f} | vs mosca {ident_mos_ct:.3f}")

# ---------------------------------------------------------------- figura
fig = plt.figure(figsize=(16, 9.2))
gs = fig.add_gridspec(3, 1, height_ratios=[0.9, 1.5, 0.9], hspace=0.42)

# Panel A: tira completa de bartmosca
axA = fig.add_subplot(gs[0])
for _, r in df.iterrows():
    axA.add_patch(Rectangle((r.pos_bart - 1, 0), 1, 1,
                            facecolor=COLORES[r.clase], edgecolor="white", linewidth=0.4))
axA.set_xlim(0, total)
axA.set_ylim(-0.15, 1.35)
axA.set_yticks([])
axA.set_xlabel("posición en bartmosca (residuo)")
axA.set_title("A) Clasificación residuo por residuo de bartmosca (107 aa) frente a humano y mosca", loc="left")
for _, r in df.iterrows():
    if r.clase in ("identico_a_mosca", "distinto_de_ambos"):
        axA.annotate(f"{r.pos_bart}{r.bart}", (r.pos_bart - 0.5, 1.03), rotation=90,
                     ha="center", va="bottom", fontsize=8,
                     color=COLORES[r.clase])
axA.text(0.5, 1.30, "N-terminal (1-25)", transform=axA.get_xaxis_transform(),
         fontsize=8, color="#555555")
axA.axvline(25.5, color="#555555", lw=0.8, ls=":")

# Panel B: zoom del N-terminal, tres filas alineadas
axB = fig.add_subplot(gs[1])
zoom_ini, zoom_fin = 1, 40
sub = df[df["col_alineamiento"] <= df.loc[df.pos_bart <= zoom_fin, "col_alineamiento"].max()]
ident_cols = set(sub["col_alineamiento"])
cols = sorted(ident_cols)
x0 = 0
for c in cols:
    fila = df[df["col_alineamiento"] == c].iloc[0] if (df["col_alineamiento"] == c).any() else None
    clase = fila["clase"] if fila is not None else "gap_en_referencia"
    color = COLORES[clase]
    for y, (nombre, seq) in enumerate([("bartmosca", bar), ("humano", hum), ("mosca", mos)]):
        ch = seq[c - 1]
        axB.add_patch(Rectangle((x0, 2 - y), 1, 1, facecolor=color if ch not in GAPS else "#ffffff",
                                edgecolor="#cccccc", linewidth=0.3))
        axB.text(x0 + 0.5, 2 - y + 0.5, ch if ch not in GAPS else "–", ha="center", va="center",
                 fontsize=9, family="monospace", color="#222222")
    x0 += 1
axB.set_xlim(0, x0)
axB.set_ylim(-0.05, 3.05)
axB.set_yticks([0.5, 1.5, 2.5])
axB.set_yticklabels(["mosca", "humano", "bartmosca"])
axB.set_xticks([i + 0.5 for i in range(0, x0, 5)])
axB.set_xticklabels([str(list(sub["col_alineamiento"])[i]) for i in range(0, x0, 5)])
axB.set_xlabel("columna del alineamiento (Clustal Omega, set3)")
axB.set_title("B) N-terminal alineado: los residuos de bartmosca que no están en el humano se concentran acá", loc="left")
axB.tick_params(axis="y", length=0)

# Panel C: conteo por clase
axC = fig.add_subplot(gs[2])
y = range(len(conteo))
axC.barh(list(y), conteo.values, color=[COLORES[k] for k in conteo.index])
axC.set_yticks(list(y))
axC.set_yticklabels([f"{ETIQUETAS[k]} ({int(v)})" for k, v in conteo.items()], fontsize=9)
axC.invert_yaxis()
axC.set_xlabel("cantidad de residuos de bartmosca")
axC.set_title("C) Casi la mitad del N-terminal difiere; del resto, bartmosca es prácticamente humano", loc="left")
for i, v in enumerate(conteo.values):
    axC.text(v + 0.6, i, str(int(v)), va="center", fontsize=9)
axC.set_xlim(0, max(conteo.values) * 1.15)

fig.suptitle("El 'niño mosca' no es mitad y mitad: bartmosca contra citocromo c humano y de mosca",
             fontsize=13, y=0.98)
fig.text(0.5, 0.005,
         "Clasificación sobre el alineamiento de Clustal Omega del set3 (109 columnas). "
         "En las 24 posiciones donde humano y mosca difieren, bartmosca sigue al humano 19 veces y a la mosca 3; "
         "además tiene 3 cambios propios en posiciones conservadas y 2 residuos extra en el N-terminal.",
         ha="center", fontsize=8.5, color="#444444")
fig.savefig(os.path.join(RES, "fig_mosaico_bartmosca.png"), dpi=200, bbox_inches="tight")
print("figura:", os.path.join(RES, "fig_mosaico_bartmosca.png"))
