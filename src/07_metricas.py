#!/usr/bin/env python3
"""Genera resultados/metricas.json con los numeros clave del TP.

Lee los artefactos ya generados (alineamientos, CSV de pares, CSV del mosaico,
logs de IQ-TREE, tabla de soportes). No recalcula nada: solo agrega.

Salida: resultados/metricas.json
"""
import json
import os
import platform
import re

import pandas as pd
from Bio import AlignIO

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "resultados")


def read_iqtree_summary(prefix):
    """Extrae modelo, BIC y logL del .iqtree."""
    path = os.path.join(RES, f"{prefix}.iqtree")
    txt = open(path).read()
    out = {}
    m = re.search(r"Best-fit model according to BIC:\s*(\S+)", txt)
    out["modelo"] = m.group(1) if m else None
    if out["modelo"] is None:
        m = re.search(r"Model of substitution:\s*(\S+)", txt)
        out["modelo"] = m.group(1) if m else None
        out["modelo_fijado"] = True
    else:
        out["modelo_fijado"] = False
    m = re.search(r"Gamma shape alpha:\s*([\d.]+)", txt)
    out["alpha_gamma"] = float(m.group(1)) if m else None
    m = re.search(r"Number of parsimony informative sites:\s*(\d+)", txt)
    out["sitios_informativos"] = int(m.group(1)) if m else None
    m = re.search(r"Number of constant sites:\s*(\d+)", txt)
    out["sitios_constantes"] = int(m.group(1)) if m else None
    m = re.search(r"Total CPU time used:.*", txt)
    m = re.search(r"Log-likelihood of the tree:\s*(-?[\d.]+)", txt)
    out["logL"] = float(m.group(1)) if m else None
    m = re.search(r"Bayesian information criterion \(BIC\) score:\s*(-?[\d.]+)", txt)
    out["bic"] = float(m.group(1)) if m else None
    m = re.search(r"Total tree length \(sum of branch lengths\):\s*(-?[\d.]+)", txt)
    out["longitud_total_arbol"] = float(m.group(1)) if m else None
    m = re.search(r"Random seed number:\s*(\d+)", txt)
    out["semilla"] = int(m.group(1)) if m else None
    m = re.search(r"Input data: (\d+) sequences with (\d+) amino-acid sites", txt)
    if m:
        out["n_secuencias"] = int(m.group(1))
        out["n_sitios"] = int(m.group(2))
    return out


def read_model_table(prefix, top=6):
    """Devuelve las filas de la tabla de ModelFinder ordenadas por BIC."""
    path = os.path.join(RES, f"{prefix}.iqtree")
    lines = open(path).read().splitlines()
    rows = []
    start = None
    for i, ln in enumerate(lines):
        if ln.startswith("Model ") and "BIC" in ln:
            start = i + 1
            break
    if start is None:
        return rows
    pat = re.compile(
        r"^(\S+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*([+-])\s+([\d.eE+-]+)\s+"
        r"(-?[\d.]+)\s*([+-])\s+([\d.eE+-]+)\s+(-?[\d.]+)\s*([+-])\s+([\d.eE+-]+)")
    for ln in lines[start:]:
        m = pat.match(ln.strip())
        if not m:
            if ln.strip() == "" and rows:
                break
            continue
        rows.append(dict(modelo=m.group(1), logL=float(m.group(2)),
                         aic=float(m.group(3)), bic=float(m.group(9)),
                         w_bic=float(m.group(11))))
        if len(rows) >= top:
            break
    return rows


def pairwise_lookup(tag):
    df = pd.read_csv(os.path.join(RES, f"pairwise_{tag}.csv"))
    out = {}
    for _, r in df.iterrows():
        a, b = str(r["par"]).split("-")
        out[(a, b)] = r
        out[(b, a)] = r      # el orden del par no importa
    return out


def main():
    import Bio
    import numpy as np
    import matplotlib
    import sys

    aln_lens = {}
    for tag in ["clustalo_set3", "mafft_set3", "clustalo_set9", "mafft_set9"]:
        p = os.path.join(RES, f"aln_{tag}.fasta")
        aln = AlignIO.read(p, "fasta")
        aln_lens[tag] = {"n_secuencias": len(aln),
                         "columnas": aln.get_alignment_length()}

    pw_cl = pairwise_lookup("clustalo_set9")
    pw_ma = pairwise_lookup("mafft_set9")
    pw3_cl = pairwise_lookup("clustalo_set3")
    pw3_ma = pairwise_lookup("mafft_set3")

    def pair_block(par, pw_c, pw_m):
        return {
            "clustalo": {
                "identidad": round(float(pw_c[par]["ident_pairdel"]), 4),
                "p_distancia": round(float(pw_c[par]["p_dist"]), 4),
                "identidad_colapsada": round(float(pw_c[par]["ident_colapsada"]), 4),
                "posiciones_comparadas": int(pw_c[par]["n_posiciones_ambas"]),
                "identidad_gaps_como_dif": round(float(pw_c[par]["ident_full"]), 4),
            },
            "mafft": {
                "identidad": round(float(pw_m[par]["ident_pairdel"]), 4),
                "p_distancia": round(float(pw_m[par]["p_dist"]), 4),
                "identidad_colapsada": round(float(pw_m[par]["ident_colapsada"]), 4),
                "posiciones_comparadas": int(pw_m[par]["n_posiciones_ambas"]),
                "identidad_gaps_como_dif": round(float(pw_m[par]["ident_full"]), 4),
            },
        }

    identidades = {
        "bartmosca-humano": pair_block(("bartmosca", "humano"), pw_cl, pw_ma),
        "bartmosca-mosca": pair_block(("bartmosca", "mosca"), pw_cl, pw_ma),
        "humano-mosca": pair_block(("humano", "mosca"), pw_cl, pw_ma),
        "cebada-humano": pair_block(("cebada", "humano"), pw_cl, pw_ma),
        "altramuz-abedul": pair_block(("altramuz", "abedul"), pw_cl, pw_ma),
        "vigna-abedul": pair_block(("vigna", "abedul"), pw_cl, pw_ma),
        "cebada-ectocarpus": pair_block(("cebada", "ectocarpus"), pw_cl, pw_ma),
    }
    identidades_set3 = {
        "bartmosca-humano": pair_block(("bartmosca", "humano"), pw3_cl, pw3_ma),
        "bartmosca-mosca": pair_block(("bartmosca", "mosca"), pw3_cl, pw3_ma),
        "humano-mosca": pair_block(("humano", "mosca"), pw3_cl, pw3_ma),
    }

    # --- mosaico ---
    vent = pd.read_csv(os.path.join(RES, "bartmosca_ventanas.csv"))
    reg = pd.read_csv(os.path.join(RES, "bartmosca_regiones.csv"))
    ventanas_validas = vent.dropna(subset=["ident_humano", "ident_mosca"])
    mosaico = {
        "ventana_aa": 15, "paso_aa": 5,
        "n_ventanas": int(len(vent)),
        "n_ventanas_con_posiciones_comparables": int(len(ventanas_validas)),
        "ventanas_mas_parecidas_a_humano": int((ventanas_validas["ident_humano"] >=
                                                ventanas_validas["ident_mosca"]).sum()),
        "ventanas_mas_parecidas_a_mosca": int((ventanas_validas["ident_humano"] <
                                               ventanas_validas["ident_mosca"]).sum()),
        "ventana_max_delta_humano": {
            "ventana": str(ventanas_validas.loc[ventanas_validas["delta_humano_menos_mosca"].idxmax(), "ventana"]),
            "ident_humano": float(ventanas_validas["ident_humano"].max()),
            "delta": round(float(ventanas_validas["delta_humano_menos_mosca"].max()), 4)},
        "regiones": {r["region"]: {"columnas_alineamiento": r["columnas_alineamiento"],
                                   "ident_humano": float(r["ident_humano"]),
                                   "ident_mosca": float(r["ident_mosca"]),
                                   "delta": float(r["delta_humano_menos_mosca"])}
                     for _, r in reg.iterrows()},
        "corte_columna_alineamiento": int(str(reg.iloc[0]["columnas_alineamiento"]).split("-")[1]),
    }

    # --- IQ-TREE ---
    s9 = read_iqtree_summary("iqtree_set9")
    s9["tabla_bic_top"] = read_model_table("iqtree_set9", top=6)
    s3 = read_iqtree_summary("iqtree_set3")
    s3["tabla_bic_top"] = read_model_table("iqtree_set3_mfp", top=12)
    s3["referencia_model_finder"] = read_iqtree_summary("iqtree_set3_mfp")

    # sitios variables/informativos del set3 con el criterio util para 3 taxones
    aln3 = AlignIO.read(os.path.join(RES, "aln_clustalo_set3.fasta"), "fasta")
    S3 = {r.id: str(r.seq) for r in aln3}
    var3 = inf3 = 0
    for j in range(aln3.get_alignment_length()):
        states = [S3[k][j] for k in S3 if S3[k][j] != "-"]
        if len(set(states)) > 1:
            var3 += 1
            cnt = pd.Series(states).value_counts()
            if cnt.max() >= 2:
                inf3 += 1
    s3["sitios_variables_criterio_3_taxones"] = var3
    s3["sitios_informativos_criterio_3_taxones"] = inf3

    # comparacion UFBoot 100 vs 1000 (submuestreo de las mismas 1000 replicas)
    comp_path = os.path.join(RES, "boot_comparacion.csv")
    if os.path.exists(comp_path):
        comp = pd.read_csv(comp_path)
        s9["ufboot_100_vs_1000"] = {
            r["nodo"]: {"n_taxones": int(r["n_taxones"]),
                        "ufboot_100": None if pd.isna(r["ufboot_100_replicas"]) else float(r["ufboot_100_replicas"]),
                        "ufboot_1000": None if pd.isna(r["ufboot_1000_replicas"]) else float(r["ufboot_1000_replicas"]),
                        "diferencia_pp": None if pd.isna(r["diferencia_pp"]) else float(r["diferencia_pp"]),
                        "bootstrap_estandar_100": None if pd.isna(r["bootstrap_estandar_100"]) else float(r["bootstrap_estandar_100"])}
            for _, r in comp.iterrows()}
    # bootstrap estandar -b 100 (modelo fijo)
    bs_summary = read_iqtree_summary("iqtree_set9_boot100")
    s9["bootstrap_estandar_100"] = bs_summary
    s9["ufboot_100_no_permitido"] = ("IQ-TREE 2.3.6 rechaza -B 100 con el error "
                                     "'#replicates must be >= 1000'; la comparacion "
                                     "100 vs 1000 se hizo submuestreando las 1000 "
                                     "replicas guardadas con --boot-trees")
    s3["nota_enraizamiento"] = ("con 3 taxones el arbol no enraizado es una trifurcacion "
                                "unica; no hay outgroup posible dentro del set")

    sop = pd.read_csv(os.path.join(RES, "soportes_nodos.csv"))
    soportes = {}
    for _, r in sop.iterrows():
        key = str(r["clado"])
        soportes[key] = {
            "n_taxones": int(r["n_taxones"]),
            "sh_alrt": None if pd.isna(r["sh_alrt"]) else float(r["sh_alrt"]),
            "ufboot": None if pd.isna(r["ufboot"]) else float(r["ufboot"]),
            "bootstrap_estandar": None if pd.isna(r["bootstrap_estandar"]) else float(r["bootstrap_estandar"]),
            "largo_rama": round(float(r["largo_rama_sust_por_sitio"]), 6),
        }

    # longitudes de rama terminales (del treefile enraizado)
    from Bio import Phylo
    rooted = Phylo.read(os.path.join(RES, "tree_set9_rooted.nwk"), "newick")
    ramas_term = {}
    for c in rooted.get_terminals():
        ramas_term[c.name] = round(c.branch_length, 6) if c.branch_length else 0.0

    # rama completa del outgroup en el treefile (sin partir)
    unrooted = Phylo.read(os.path.join(RES, "tree_set9_unrooted.nwk"), "newick")
    out_branch = [c.branch_length for c in unrooted.root.clades
                  if c.is_terminal() and c.name == "rhizobium"]
    unrooted_newick = open(os.path.join(RES, "tree_set9_unrooted.nwk")).read().strip()
    rooted_newick = open(os.path.join(RES, "tree_set9_rooted.nwk")).read().strip()

    ebi = {}
    for ln in open(os.path.join(RES, "ebi_jobids.txt")):
        parts = ln.rstrip("\n").split("\t")
        if len(parts) >= 2:
            ebi[parts[0]] = parts[1]

    metricas = {
        "fecha": "2026-10-27",
        "zona_horaria": "America/Buenos_Aires",
        "versiones": {
            "iqtree": "2.3.6 (Linux x86 64-bit, 2024-08-01)",
            "mafft": "v7.526 (2024/Apr/26)",
            "clustalo": "Clustal Omega via REST de EBI (https://www.ebi.ac.uk/Tools/services/rest/clustalo)",
            "biopython": Bio.__version__,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "ebi_job_ids": ebi,
        "longitudes_alineamiento": aln_lens,
        "identidades_pares_set9": identidades,
        "identidades_pares_set3": identidades_set3,
        "mosaico_bartmosca": mosaico,
        "iqtree_set9": {
            **s9,
            "topologia_newick": unrooted_newick,
            "rama_outgroup_completa": (round(out_branch[0], 6) if out_branch else None),
            "soportes_por_nodo": soportes,
            "longitudes_rama_terminales": ramas_term,
            "newick_enraizado": rooted_newick,
        },
        "iqtree_set3": {
            **s3,
            "newick_sin_enraizar": open(os.path.join(RES, "tree_set3.nwk")).read().strip(),
            "longitudes_rama_terminales": {
                c.name: round(c.branch_length, 6)
                for c in Phylo.read(os.path.join(RES, "tree_set3.nwk"), "newick").get_terminals()},
        },
        "figuras": {name: f"resultados/{name}" for name in [
            "fig_arbol_set9.png", "fig_arbol_set3.png",
            "fig_ventanas_bartmosca.png", "fig_heatmap_identidad.png"]},
        "archivos_clave": {name: f"resultados/{name}" for name in [
            "aln_clustalo_set9.fasta", "aln_mafft_set9.fasta",
            "aln_clustalo_set3.fasta", "aln_mafft_set3.fasta",
            "pairwise_clustalo_set9.csv", "pairwise_mafft_set9.csv",
            "bartmosca_ventanas.csv", "bartmosca_regiones.csv",
            "soportes_nodos.csv",
            "tree_set9_rooted.nwk", "tree_set9_unrooted.nwk", "tree_set3.nwk",
            "iqtree_set9.treefile", "iqtree_set9.iqtree",
            "iqtree_set9_boot100.treefile", "iqtree_set9_boot100.iqtree",
            "iqtree_set9_B1000.treefile", "iqtree_set9_B1000.ufboot",
            "boot_comparacion.csv", "iqtree_set3.treefile", "iqtree_set3.iqtree",
            "iqtree_set3_mfp.iqtree"]},
    }

    out = os.path.join(RES, "metricas.json")
    with open(out, "w") as fh:
        json.dump(metricas, fh, indent=2, ensure_ascii=False)
    print("escrito", out)
    return metricas


if __name__ == "__main__":
    main()
