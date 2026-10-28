#!/usr/bin/env python3
"""Arboles y figuras del TP.

Entrada:
  resultados/iqtree_set9.treefile          (ML set9, -B 1000 -alrt 1000 -o rhizobium)
  resultados/iqtree_set9_boot1000.treefile (ML set9, bootstrap estandar -b 1000)
  resultados/iqtree_set3.treefile          (ML set3, sin outgroup)

Salida:
  resultados/tree_set9_rooted.nwk    (enraizado por outgroup rhizobium)
  resultados/tree_set9_unrooted.nwk  (copia exacta del .treefile, sin enraizar)
  resultados/tree_set3.nwk           (sin enraizar: con 3 taxones no hay outgroup)
  resultados/fig_arbol_set9.png
  resultados/fig_arbol_set3.png
  resultados/soportes_nodos.csv      (tabla de soportes por nodo)

Notas:
  - Se enraiza en la rama que lleva a rhizobium (outgroup). La rama del outgroup
    se parte a la mitad para colgar la raiz: la longitud total de la rama
    (0.8734 subs/sitio en el .treefile) se conserva.
  - El arbol de bootstrap estandar se usa solo para leer los soportes de los
    clados; la topologia dibujada es la del arbol ML con UFBoot/SH-aLRT.
"""
import os
import sys
import json
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "resultados")

OUTGROUP = "rhizobium"
SHOW_IDS = ["humano", "mosca", "bartmosca", "cebada", "altramuz", "abedul",
            "vigna", "ectocarpus", "rhizobium"]


# ----------------------------------------------------------------------------
# utilidades de arbol
# ----------------------------------------------------------------------------
def clade_key(clade):
    return frozenset(c.name for c in clade.get_terminals())


def split_support(label):
    """'78/72' -> (78.0, 72.0); '72' -> (None, 72.0)."""
    if label is None:
        return None, None
    parts = str(label).split("/")
    try:
        if len(parts) == 2:
            return float(parts[0]), float(parts[1])
        return None, float(parts[0])
    except ValueError:
        return None, None


def read_support_map(path):
    """Devuelve {frozenset(hojas): bootstrap%} a partir de un .treefile de -b."""
    if not os.path.exists(path):
        return {}
    t = Phylo.read(path, "newick")
    out = {}
    for c in t.get_nonterminals():
        if c is t.root:
            continue
        lab = c.name if c.name is not None else c.confidence
        _, bs = split_support(lab)
        if bs is not None:
            out[clade_key(c)] = bs
    return out


def root_on_outgroup(tree, outgroup=OUTGROUP):
    """Enraiza el arbol (que viene con el outgroup como hijo de la raiz
    trifurcada) partiendo por la mitad la rama del outgroup."""
    root = tree.root
    out_child = None
    for c in root.clades:
        if c.is_terminal() and c.name == outgroup:
            out_child = c
    if out_child is None:
        raise SystemExit(f"ERROR: no encontre '{outgroup}' como hijo de la raiz")
    L = out_child.branch_length or 0.0
    ingroup = Clade(clades=[c for c in root.clades if c is not out_child],
                    branch_length=L / 2.0)
    rhizo = Clade(name=outgroup, branch_length=L / 2.0)
    new_root = Clade(clades=[ingroup, rhizo])
    return Tree(root=new_root, rooted=True)


def depth_map(tree):
    return {id(c): d for c, d in tree.depths().items()} if hasattr(tree, "depths") else {}


# ----------------------------------------------------------------------------
# dibujo rectangular propio (matplotlib)
# ----------------------------------------------------------------------------
def draw_rect_tree(ax, tree, support_map, bs_map, colors=None, xpad=0.004,
                   fontsize=9, support_fontsize=8):
    depths = tree.depths()
    leaves = tree.get_terminals()
    ypos = {id(l): i for i, l in enumerate(leaves)}

    def y_of(clade):
        if id(clade) in ypos:
            return ypos[id(clade)]
        ys = [y_of(c) for c in clade.clades]
        return (min(ys) + max(ys)) / 2.0

    def walk(clade, parent=None):
        x = depths[clade]
        y = y_of(clade)
        if parent is not None:
            px = depths[parent]
            ax.plot([px, x], [y, y], color="black", lw=1.1, solid_capstyle="butt", zorder=2)
        if not clade.is_terminal():
            ys = [y_of(c) for c in clade.clades]
            ax.plot([x, x], [min(ys), max(ys)], color="black", lw=1.1, zorder=2)
            # soportes
            key = clade_key(clade)
            vals = []
            if key in support_map:
                sh, uf = support_map[key]
                if sh is not None:
                    vals.append(f"{sh:g}/{uf:g}")
                else:
                    vals.append(f"{uf:g}")
            if key in bs_map:
                vals.append(f"bs {bs_map[key]:g}")
            if vals and clade is not tree.root:
                # a la izquierda del nodo, por encima de la rama que entra
                ax.text(x - xpad, y - 0.10, "\n".join(vals),
                        fontsize=support_fontsize, ha="right", va="bottom",
                        color="#1a1a1a", zorder=6, linespacing=1.05,
                        bbox=dict(boxstyle="round,pad=0.08", fc="white",
                                  ec="none", alpha=0.85))
        for c in clade.clades:
            walk(c, clade)

    walk(tree.root)
    for l in leaves:
        x = depths[l]
        ax.plot([x, x], [ypos[id(l)], ypos[id(l)]], marker="o", ms=0)
        if colors and l.name in colors:
            ax.text(x + xpad, ypos[id(l)], l.name, fontsize=fontsize, va="center",
                    ha="left", color=colors[l.name], fontweight="bold", zorder=7)
        else:
            ax.text(x + xpad, ypos[id(l)], l.name, fontsize=fontsize, va="center",
                    ha="left")
    ax.set_ylim(len(leaves) + 1.3, -0.9)   # eje y invertido: arriba = primer taxon
    xmax = max(depths.values())
    ax.set_xlim(0, xmax * 1.25)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_yticks([])
    ax.set_xlabel("sustituciones por sitio", fontsize=10)
    return depths, leaves


def scale_bar(ax, x0=None, y=None, length=0.1, label="0.1"):
    """Barra de escala en el hueco inferior del arbol (coordenadas de datos)."""
    xmin, xmax = ax.get_xlim()
    if x0 is None:
        x0 = xmin + 0.58 * (xmax - xmin)
    if y is None:
        y = ax.get_ylim()[0] - 0.55   # ylim esta invertido: restar sube
    ax.plot([x0, x0 + length], [y, y], color="black", lw=2.6, zorder=10)
    ax.text(x0 + length / 2, y - 0.22, label, ha="center", va="bottom",
            fontsize=9, zorder=10)
    ax.text(x0 + length / 2, y + 0.22, "sustituciones/sitio", ha="center",
            va="top", fontsize=8, color="#333333", zorder=10)


# ----------------------------------------------------------------------------
# set9
# ----------------------------------------------------------------------------
def arbol_set9():
    tf = os.path.join(RES, "iqtree_set9.treefile")
    boot_tf = os.path.join(RES, "iqtree_set9_boot100.treefile")
    tree = Phylo.read(tf, "newick")

    # supports del arbol principal: SH-aLRT / UFBoot
    sup = {}
    for c in tree.get_nonterminals():
        if c is tree.root:
            continue
        if c.name:
            sup[clade_key(c)] = split_support(c.name)
    bs = read_support_map(boot_tf)

    # --- newick sin enraizar (copia del treefile) ---
    unrooted = open(tf).read().strip()
    with open(os.path.join(RES, "tree_set9_unrooted.nwk"), "w") as fh:
        fh.write(unrooted + "\n")

    # --- newick enraizado por outgroup ---
    rooted = root_on_outgroup(tree)
    rooted.ladderize()  # no cambia topologia, ordena hojas
    Phylo.write(rooted, os.path.join(RES, "tree_set9_rooted.nwk"), "newick")

    # --- tabla de soportes ---
    artificial = [c for c in rooted.root.clades if not c.is_terminal()]
    rows = []
    for c in rooted.get_nonterminals():
        if c is rooted.root or c in artificial:
            continue
        key = clade_key(c)
        sh, uf = sup.get(key, (None, None))
        names = sorted(key)
        if len(names) == 1:
            continue
        rows.append(dict(
            clado="+".join(names), n_taxones=len(names),
            sh_alrt=sh, ufboot=uf,
            bootstrap_estandar=bs.get(key),
            largo_rama_sust_por_sitio=round(c.branch_length or 0.0, 6)))
    df = pd.DataFrame(rows).sort_values("n_taxones", ascending=False)
    df.to_csv(os.path.join(RES, "soportes_nodos.csv"), index=False)
    print(df.to_string(index=False))
    if bs and len(bs) != len(sup):
        print(f"AVISO: el arbol de bootstrap tiene {len(bs)} nodos internos con soporte "
              f"y el ML {len(sup)}", file=sys.stderr)

    # --- figura ---
    colors = {"bartmosca": "#c0392b"}
    fig, ax = plt.subplots(figsize=(11.5, 6.2), dpi=150)
    draw_rect_tree(ax, rooted, sup, bs, colors=colors)
    nboot = 100
    ax.set_title("Arbol filogenetico de citocromo c (set9, 9 secuencias, 131 sitios)\n"
                 "IQ-TREE 2.3.6, ModelFinder (MFP) -> Dayhoff+G4; nodos: SH-aLRT/UFBoot "
                 f"y bootstrap estandar {nboot} replicas (bs); enraizado por outgroup rhizobium",
                 fontsize=11)
    scale_bar(ax, length=0.1, label="0.1")
    fig.text(0.01, 0.015,
             "Longitudes de rama en sustituciones por sitio. La rama del outgroup se "
             "partio a la mitad al enraizar (rama ML completa = 0.8734). "
             "Soportes: SH-aLRT/UFBoot del arbol ML y bootstrap estandar de 100 replicas.",
             fontsize=7.5, color="#444444")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out = os.path.join(RES, "fig_arbol_set9.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("figura:", out)
    return rooted, df


# ----------------------------------------------------------------------------
# set3
# ----------------------------------------------------------------------------
def arbol_set3():
    tf = os.path.join(RES, "iqtree_set3.treefile")
    if not os.path.exists(tf):
        print(f"FALTA {tf}", file=sys.stderr)
        return
    tree = Phylo.read(tf, "newick")
    newick = open(tf).read().strip()
    with open(os.path.join(RES, "tree_set3.nwk"), "w") as fh:
        fh.write(newick + "\n")

    # figura: trifurcacion no enraizada + longitudes de rama
    lengths = {c.name: (c.branch_length or 0.0) for c in tree.get_terminals()}
    print("set3 longitudes:", lengths)

    # modelo y logL del .iqtree (modelo fijado a Dayhoff+G4)
    import re as _re
    txt = open(os.path.join(RES, "iqtree_set3.iqtree")).read()
    m = _re.search(r"Model of substitution:\s*(\S+)", txt)
    modelo3 = m.group(1) if m else "?"
    m = _re.search(r"Log-likelihood of the tree:\s*(-?[\d.]+)", txt)
    logl3 = float(m.group(1)) if m else float("nan")

    import math
    maxL = max(lengths.values())
    ang = {"humano": 105, "bartmosca": 15, "mosca": 225}
    colors = {"bartmosca": "#c0392b"}

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(11.8, 5.6), dpi=150,
                                  gridspec_kw={"width_ratios": [1.15, 1.0]})

    # --- panel A: estrella no enraizada ---
    scale = 1.0 / maxL
    for name, L in lengths.items():
        a = math.radians(ang[name])
        r = L * scale
        x, y = r * math.cos(a), r * math.sin(a)
        col = colors.get(name, "black")
        ax.plot([0, x], [0, y], color=col, lw=1.9, zorder=2)
        ax.plot([x], [y], "o", color=col, ms=5, zorder=3)
        dx, dy = math.cos(a), math.sin(a)
        ha = "left" if dx > 0.2 else ("right" if dx < -0.2 else "center")
        va = "bottom" if dy > 0.2 else ("top" if dy < -0.2 else "center")
        ax.text(x + dx * 0.07, y + dy * 0.07,
                f"{name}\n{L:.4f}", ha=ha, va=va, fontsize=10,
                color=col, fontweight=("bold" if name in colors else "normal"))
    ax.plot([0], [0], "o", color="black", ms=6, zorder=4)
    ax.set_xlim(-1.02, 0.62)
    ax.set_ylim(-0.95, 0.62)
    ax.set_aspect("equal")
    ax.axis("off")
    # barra de escala (coordenadas de ejes)
    frac = 0.1 / (1.64)
    ax.plot([0.66, 0.66 + frac], [0.06, 0.06], transform=ax.transAxes,
            color="black", lw=2.6, clip_on=False)
    ax.text(0.66 + frac / 2, 0.025, "0.1 sust./sitio", transform=ax.transAxes,
            ha="center", va="top", fontsize=8.5)
    ax.set_title("A) Arbol no enraizado (trifurcacion unica)", fontsize=10.5)

    # --- panel B: longitudes de rama ---
    order = ["bartmosca", "humano", "mosca"]
    vals = [lengths[k] for k in order]
    bars = axb.barh(order, vals, color=[colors.get(k, "#4a6fa5") for k in order],
                    height=0.55)
    for b, v in zip(bars, vals):
        axb.text(v + 0.004, b.get_y() + b.get_height() / 2, f"{v:.4f}",
                 va="center", fontsize=10)
    axb.set_xlim(0, max(vals) * 1.30)
    axb.set_xlabel("longitud de rama (sustituciones por sitio)", fontsize=10)
    axb.set_title("B) Longitudes de rama terminales", fontsize=10.5)
    axb.grid(axis="x", alpha=0.3)
    axb.set_axisbelow(True)
    for sp in ("top", "right"):
        axb.spines[sp].set_visible(False)

    fig.suptitle("Set3 (humano, mosca, bartmosca): arbol NO enraizado\n"
                 f"IQ-TREE 2.3.6, modelo fijo {modelo3}, logL = {logl3:.4f}",
                 fontsize=12)
    fig.text(0.5, 0.015,
             "Con 3 taxones la topologia no enraizada es unica (trifurcacion): no hay nodo "
             "interno que resolver y no se puede enraizar por outgroup.\n"
             "Las longitudes de rama son la unica senal: bartmosca y humano tienen ramas "
             f"cortas ({lengths['bartmosca']:.4f} y {lengths['humano']:.4f}) y mosca es la "
             f"mas divergente ({lengths['mosca']:.4f}).",
             ha="center", fontsize=8, color="#333333")
    fig.tight_layout(rect=(0, 0.075, 1, 0.94))
    out = os.path.join(RES, "fig_arbol_set3.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("figura:", out)


if __name__ == "__main__":
    arbol_set9()
    arbol_set3()
