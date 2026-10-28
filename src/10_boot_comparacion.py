#!/usr/bin/env python3
"""Efecto del numero de replicas: UFBoot 100 vs UFBoot 1000.

IQ-TREE 2.3.6 NO permite `-B 100` (exige >= 1000 replicas; error
"#replicates must be >= 1000"), asi que la comparacion se hace sobre las mismas
1000 replicas de UFBoot guardadas con `--boot-trees`:

  resultados/iqtree_set9_B1000.ufboot   (1000 arboles de bootstrap)

Para cada nodo se calcula la frecuencia del split usando las **primeras 100**
replicas y las **1000** replicas (la primera es un subconjunto de la segunda, que
es exactamente el efecto de reducir el numero de replicas). Se valida contra los
soportes que IQ-TREE escribio en el .treefile.

Salida: resultados/boot_comparacion.csv
"""
import os
from Bio import Phylo
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "resultados")

NODOS = {
    "humano+bartmosca": {"humano", "bartmosca"},
    "clado animal (humano+mosca+bartmosca)": {"humano", "mosca", "bartmosca"},
    "clado de plantas verdes (cebada+altramuz+abedul+vigna)": {"cebada", "altramuz",
                                                              "abedul", "vigna"},
    "plantas verdes + ectocarpus": {"cebada", "altramuz", "abedul", "vigna",
                                    "ectocarpus"},
    "todos menos rhizobium (raiz del ingroup)": {"humano", "mosca", "bartmosca",
                                                 "cebada", "altramuz", "abedul",
                                                 "vigna", "ectocarpus"},
}


def bipartitions(tree, all_taxa):
    """Conjunto de bipartitions (lado menor, canonico) de un arbol sin enraizar."""
    n = len(all_taxa)
    splits = set()
    for c in tree.find_clades():
        s = frozenset(x.name for x in c.get_terminals())
        if not s or s == all_taxa:
            continue
        if len(s) > n / 2:
            s = all_taxa - s
        splits.add(frozenset(s))
    return splits


def target_key(taxa, all_taxa):
    s = frozenset(taxa)
    if len(s) > len(all_taxa) / 2:
        s = all_taxa - s
    return frozenset(s)


def main():
    ufboot = os.path.join(RES, "iqtree_set9_B1000.ufboot")
    trees = list(Phylo.parse(ufboot, "newick"))
    print(f"replicas UFBoot leidas: {len(trees)}")
    all_taxa = frozenset(x.name for x in trees[0].get_terminals())
    print("taxones:", sorted(all_taxa))

    splits_100 = [bipartitions(t, all_taxa) for t in trees[:100]]
    splits_1000 = [bipartitions(t, all_taxa) for t in trees]

    # validacion: soportes del .treefile
    tf = Phylo.read(os.path.join(RES, "iqtree_set9_B1000.treefile"), "newick")
    tf_sup = {}
    for c in tf.get_nonterminals():
        if c is tf.root:
            continue
        lab = c.name if c.name is not None else c.confidence
        try:
            key = target_key([x.name for x in c.get_terminals()], all_taxa)
            tf_sup[key] = float(str(lab).split("/")[0])
        except (TypeError, ValueError):
            pass

    # bootstrap estandar 100 (si ya termino)
    bs100 = {}
    p = os.path.join(RES, "iqtree_set9_boot100.treefile")
    if os.path.exists(p):
        bt = Phylo.read(p, "newick")
        for c in bt.get_nonterminals():
            if c is bt.root:
                continue
            lab = c.name if c.name is not None else c.confidence
            try:
                bs100[target_key([x.name for x in c.get_terminals()], all_taxa)] = \
                    float(str(lab).split("/")[0])
            except (TypeError, ValueError):
                pass

    rows = []
    for nombre, taxa in NODOS.items():
        key = target_key(taxa, all_taxa)
        n100 = sum(1 for s in splits_100 if key in s)
        n1000 = sum(1 for s in splits_1000 if key in s)
        rows.append({
            "nodo": nombre,
            "n_taxones": len(taxa),
            "ufboot_100_replicas": round(100 * n100 / 100, 1),
            "ufboot_1000_replicas": round(100 * n1000 / 1000, 1),
            "diferencia_pp": round(100 * n100 / 100 - 100 * n1000 / 1000, 1),
            "ufboot_1000_segun_treefile": tf_sup.get(key),
            "bootstrap_estandar_100": bs100.get(key),
        })
    df = pd.DataFrame(rows)
    out = os.path.join(RES, "boot_comparacion.csv")
    df.to_csv(out, index=False)
    print()
    print(df.to_string(index=False))
    print("\nescrito", out)

    # verificacion: los soportes calculados con 1000 replicas deben coincidir
    # con los que IQ-TREE reporta en el .treefile (los nodos del arbol ML)
    print("\nvalidacion (calculado sobre las 1000 replicas vs .treefile, que redondea):")
    for _, r in df.iterrows():
        if r["ufboot_1000_segun_treefile"] is not None:
            ok = abs(r["ufboot_1000_replicas"] - r["ufboot_1000_segun_treefile"]) < 0.6
            print(f"  {r['nodo']}: {r['ufboot_1000_replicas']:.1f} vs "
                  f"{r['ufboot_1000_segun_treefile']:.1f} -> {'OK (redondeo)' if ok else 'DIFERENCIA'}")
        else:
            print(f"  {r['nodo']}: {r['ufboot_1000_replicas']:.1f} (nodo no presente "
                  "en el .treefile: es la raiz del ingroup)")


if __name__ == "__main__":
    main()
