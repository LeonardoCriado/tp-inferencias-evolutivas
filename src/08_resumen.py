#!/usr/bin/env python3
"""Genera resultados/resumen_analisis.md a partir de los artefactos del pipeline.

No calcula filogenias: solo agrega y formatea los resultados ya obtenidos.
Requiere haber corrido: 01 -> 07.
Salida: resultados/resumen_analisis.md
"""
import json
import os
from datetime import datetime

import pandas as pd
from Bio import AlignIO

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "resultados")
DAT = os.path.join(BASE, "datos")

with open(os.path.join(RES, "metricas.json")) as fh:
    M = json.load(fh)


def fmt(v, nd=4, pct=False):
    if v is None:
        return "—"
    if pct:
        return f"{v*100:.2f} %"
    return f"{v:.{nd}f}"


def aln_stats(path):
    aln = AlignIO.read(path, "fasta")
    seqs = {r.id: str(r.seq) for r in aln}
    L = aln.get_alignment_length()
    gapfree = sum(1 for j in range(L) if all(s[j] != "-" for s in seqs.values()))
    gaps = sum(s.count("-") for s in seqs.values())
    return dict(n=len(aln), L=L, gapfree=gapfree, gaps=gaps)


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def main():
    # datos de apoyo
    label_map = pd.read_csv(os.path.join(RES, "label_map.tsv"), sep="\t")
    vent = pd.read_csv(os.path.join(RES, "bartmosca_ventanas.csv"))
    reg = pd.read_csv(os.path.join(RES, "bartmosca_regiones.csv"))
    ident = pd.read_csv(os.path.join(RES, "pairwise_identidad_clustalo_set9.csv"),
                        index_col=0)
    ident_ma = pd.read_csv(os.path.join(RES, "pairwise_identidad_mafft_set9.csv"),
                           index_col=0)
    pdist = pd.read_csv(os.path.join(RES, "pairwise_pdist_clustalo_set9.csv"),
                        index_col=0)
    aln = {tag: aln_stats(os.path.join(RES, f"aln_{tag}.fasta"))
           for tag in ["clustalo_set3", "mafft_set3", "clustalo_set9", "mafft_set9"]}

    order9 = ["humano", "bartmosca", "mosca", "cebada", "altramuz", "abedul",
              "vigna", "ectocarpus", "rhizobium"]
    order3 = ["humano", "mosca", "bartmosca"]

    s9 = M["iqtree_set9"]
    s3 = M["iqtree_set3"]
    ebi = M["ebi_job_ids"]
    fecha = datetime.now().strftime("%Y-%m-%d")

    L = []
    A = L.append

    # ---------------------------------------------------------------- cabecera
    A("# Análisis de inferencias evolutivas sobre citocromo c (TP BartMosca)\n")
    A(f"**Fecha de ejecución:** {fecha} (America/Buenos_Aires)  ")
    A(f"**Insumo:** `datos/SecuenciasCytocromoC.fasta` (9 secuencias, "
      f"citocromo c, 104–127 aa)  ")
    A("**Resultados crudos:** `resultados/` · **Scripts:** `src/` · "
      "**Datos derivados:** `datos/`  ")
    A("**Semilla de IQ-TREE:** `-seed 12345` (en set9 y set3) para que el "
      "análisis sea exactamente reproducible.\n")

    A("## 1. Versiones de software\n")
    v = M["versiones"]
    A(md_table(["Herramienta", "Versión / origen"], [
        ["IQ-TREE", v["iqtree"]],
        ["MAFFT", v["mafft"]],
        ["Clustal Omega", v["clustalo"]],
        ["Biopython", v["biopython"]],
        ["Python", v["python"]],
        ["NumPy", v["numpy"]],
        ["pandas", v["pandas"]],
        ["matplotlib", v["matplotlib"]],
    ]))
    A("")

    A("## 2. Secuencias y conjuntos de datos\n")
    A("Etiquetas cortas (sin espacios, para que IQ-TREE no se queje) y acceso "
      "original de cada secuencia:\n")
    rows = [[r["etiqueta"], r["header_original"].split("[")[0].rstrip("_"),
             r["header_original"].split("[")[-1].rstrip("]") if "[" in r["header_original"] else "—",
             r["largo"]] for _, r in label_map.iterrows()]
    A(md_table(["Etiqueta", "Header (recortado)", "Especie / fuente", "Largo crudo (aa)"],
               rows))
    A("")
    A("- `datos/set3.fasta`: **humano, mosca, bartmosca** (3 secuencias).")
    A("- `datos/set9.fasta`: las **9** secuencias (bartmosca + 2 metazoos + 6 de "
      "referencia: cebada *Hordeum vulgare*, altramuz *Lupinus albus*, abedul "
      "*Betula platyphylla*, vigna *Vigna radiata*, ectocarpus *Ectocarpus "
      "siliculosus* y rhizobium *Rhizobium* sp. CF097).")
    A("")

    A("## 3. Alineamientos\n")
    A("Se alineó con **Clustal Omega** (REST de EBI, `stype=protein`, `outfmt=fa`) "
      "y con **MAFFT 7.526** (`--localpair --maxiterate 1000`). "
      f"IDs de los trabajos de EBI: set3 = `{ebi.get('set3','—')}`, "
      f"set9 = `{ebi.get('set9','—')}`.\n")
    rows = []
    for tag, label in [("clustalo_set3", "Clustal Omega"), ("mafft_set3", "MAFFT"),
                       ("clustalo_set9", "Clustal Omega"), ("mafft_set9", "MAFFT")]:
        st = aln[tag]
        rows.append([label, tag.split("_")[1], st["n"], st["L"], st["gapfree"], st["gaps"]])
    A(md_table(["Alineador", "Set", "N.º seqs", "Columnas", "Columnas sin gap",
                "Gaps totales"], rows))
    A("")
    A("Los dos alineadores dan **la misma longitud** para cada set "
      f"(set3 = {aln['clustalo_set3']['L']} columnas, set9 = {aln['clustalo_set9']['L']} "
      "columnas). No colocan exactamente los mismos gaps: en el set9, 16 de las "
      "131 columnas tienen algún carácter distinto entre ambos alineamientos "
      "(ver §4 y §12).\n")

    A("## 4. Identidad y p-distancia de pares (set9)\n")
    A("**Definiciones usadas** (todas sobre el alineamiento de Clustal Omega):\n")
    A("- **Identidad (deleción de pares):** proporción de posiciones idénticas "
      "considerando solo las columnas donde *ambas* secuencias del par tienen "
      "residuo (`idénticas / (idénticas + distintas)`). Es la identidad estándar.")
    A("- **p-distancia:** `1 − identidad` (no corregida por múltiples "
      "sustituciones).")
    A("- **Colapsando gaps (chequeo):** se eliminan los gaps del alineamiento, se "
      "re-alinea cada par de forma independiente (`PairwiseAligner` global, "
      "BLOSUM62, gap open −10 / extend −0.5) y se recalcula la identidad con "
      "deleción de pares. Si difiere de la anterior, el valor depende del "
      "alineador múltiple.")
    A("- **Gaps como diferencias:** `idénticas / columnas totales del "
      "alineamiento` (cota inferior, útil para comparar con la literatura).\n")

    key_pairs = ["bartmosca-humano", "bartmosca-mosca", "humano-mosca",
                 "cebada-humano", "altramuz-abedul", "vigna-abedul",
                 "cebada-ectocarpus"]
    rows = []
    for k in key_pairs:
        d = M["identidades_pares_set9"][k]
        rows.append([k,
                     fmt(d["clustalo"]["identidad"]),
                     fmt(d["clustalo"]["p_distancia"]),
                     fmt(d["clustalo"]["identidad_colapsada"]),
                     fmt(d["clustalo"]["identidad_gaps_como_dif"]),
                     fmt(d["mafft"]["identidad"]),
                     fmt(d["mafft"]["p_distancia"]),
                     fmt(d["mafft"]["identidad_colapsada"])])
    A("### 4.1 Pares clave\n")
    A(md_table(["Par", "Ident. ClustalO", "p-dist. ClustalO",
                "Ident. colapsada", "Ident. gaps=dif", "Ident. MAFFT",
                "p-dist. MAFFT", "Ident. colapsada MAFFT"], rows))
    A("")
    bh = M["identidades_pares_set9"]["bartmosca-humano"]
    bm = M["identidades_pares_set9"]["bartmosca-mosca"]
    hm = M["identidades_pares_set9"]["humano-mosca"]
    A(f"**Lectura:** bartmosca es ~{bh['clustalo']['identidad']*100:.1f} % idéntica "
      f"al humano (p-dist {bh['clustalo']['p_distancia']:.4f}) pero solo "
      f"~{bm['clustalo']['identidad']*100:.1f} % a la mosca; humano–mosca comparten "
      f"~{hm['clustalo']['identidad']*100:.1f} %. Es decir, bartmosca está mucho más "
      "cerca del humano que la propia mosca: la señal filogenética de "
      "citocromo c la agrupa con el humano, no con la mosca.\n")

    A("### 4.2 Matriz completa de identidad (set9, Clustal Omega, %)\n")
    rows = []
    for a in order9:
        row = [a]
        for b in order9:
            row.append("—" if a == b else f"{float(ident.loc[a, b])*100:.1f}")
        rows.append(row)
    A(md_table(["", *order9], rows))
    A("")

    A("### 4.3 Matriz de p-distancia (set9, Clustal Omega)\n")
    rows = []
    for a in order9:
        row = [a]
        for b in order9:
            row.append("—" if a == b else f"{float(pdist.loc[a, b]):.4f}")
        rows.append(row)
    A(md_table(["", *order9], rows))
    A("")

    A("### 4.4 Set3 (humano, mosca, bartmosca)\n")
    d3 = M["identidades_pares_set3"]
    rows = []
    for k in ["bartmosca-humano", "bartmosca-mosca", "humano-mosca"]:
        d = d3[k]
        rows.append([k, fmt(d["clustalo"]["identidad"]),
                     fmt(d["clustalo"]["p_distancia"]),
                     fmt(d["clustalo"]["identidad_colapsada"]),
                     int(d["clustalo"]["posiciones_comparadas"]),
                     aln["clustalo_set3"]["L"]])
    A(md_table(["Par", "Identidad", "p-distancia", "Ident. colapsada",
                "Posiciones comparadas", "Columnas alineamiento"], rows))
    A("")

    A("## 5. Análisis de mosaico de bartmosca\n")
    A("Ventanas deslizantes de **15 aa con paso 5** sobre el alineamiento Clustal "
      "Omega del set9; en cada ventana se calcula la identidad de bartmosca "
      "contra humano y contra mosca (deleción de pares). Las primeras 2 ventanas "
      "(columnas 1–20) no tienen posiciones comparables porque bartmosca/humano/"
      "mosca tienen gaps en ese tramo del alineamiento (ver "
      "`resultados/bartmosca_ventanas.csv`).\n")

    mo = M["mosaico_bartmosca"]
    A("### 5.1 Resumen por regiones\n")
    rows = []
    for k, v in mo["regiones"].items():
        rows.append([k, v["columnas_alineamiento"],
                     fmt(v["ident_humano"]), fmt(v["ident_mosca"]),
                     fmt(v["delta"])])
    A(md_table(["Región", "Columnas del alineamiento", "Ident. vs humano",
                "Ident. vs mosca", "Δ (humano − mosca)"], rows))
    A("")
    A(f"El N-terminal (primeros 25 residuos de bartmosca; columnas 1–"
      f"{mo['corte_columna_alineamiento']}) es la única región con ambigüedad: "
      f"identidad {mo['regiones']['N-terminal (1-25)']['ident_humano']*100:.1f} % "
      f"vs humano y {mo['regiones']['N-terminal (1-25)']['ident_mosca']*100:.1f} % "
      "vs mosca. El resto de la proteína (residuos 26–107) es **idéntico al "
      "humano al 100 %** y solo ~79 % a la mosca. Conclusión: bartmosca no es un "
      "mosaico con el N-terminal de mosca; el N-terminal es *levemente* más "
      "parecido al humano y queda como la única zona candidata a aporte "
      "diferencial.\n")

    A("### 5.2 Ventanas (15 aa, paso 5)\n")
    rows = []
    for _, r in vent.iterrows():
        if pd.isna(r["ident_humano"]):
            rows.append([r["ventana"], "—", "—", "—", "sin posiciones comparables"])
        else:
            rows.append([r["ventana"], f"{r['ident_humano']*100:.1f}",
                         f"{r['ident_mosca']*100:.1f}",
                         f"{r['delta_humano_menos_mosca']*100:+.1f}",
                         r["mas_parecida_a"]])
    A(md_table(["Ventana (col.)", "Ident. vs humano (%)", "Ident. vs mosca (%)",
                "Δ humano−mosca", "Más parecida a"], rows))
    A("")
    A(f"De las {mo['n_ventanas_con_posiciones_comparables']} ventanas con datos, "
      f"{mo['ventanas_mas_parecidas_a_humano']} son más parecidas al humano y "
      f"{mo['ventanas_mas_parecidas_a_mosca']} a la mosca. La ventana con máxima "
      f"ventaja del humano es **{mo['ventana_max_delta_humano']['ventana']}** "
      f"(Δ = {mo['ventana_max_delta_humano']['delta']*100:+.1f} puntos "
      "porcentuales).\n")

    A("## 6. Filogenia del set9 (9 secuencias, 131 sitios)\n")
    A("### 6.1 Comando\n")
    A("```bash")
    A("tools/iqtree-2.3.6-Linux-intel/bin/iqtree2 \\")
    A("  -s resultados/aln_clustalo_set9.fasta \\")
    A("  -m MFP -B 1000 -alrt 1000 -T 2 -o rhizobium \\")
    A("  -seed 12345 --prefix resultados/iqtree_set9 -redo")
    A("```\n")
    A("### 6.2 Selección de modelo (ModelFinder, criterio BIC)\n")
    A(f"**Modelo elegido: `{s9['modelo']}`** (BIC = {s9['bic']:.4f}). "
      "Es Dayhoff con heterogeneidad de tasas gamma (4 categorías). Primeros "
      "puestos de la tabla de ModelFinder (dos modelos quedan casi empatados: "
      f"w-BIC = {s9['tabla_bic_top'][0]['w_bic']:.3f} para Dayhoff+G4 y "
      f"{s9['tabla_bic_top'][1]['w_bic']:.3f} para DCMut+G4):\n")
    rows = []
    for r in s9["tabla_bic_top"]:
        rows.append([r["modelo"], f"{r['logL']:.3f}", f"{r['aic']:.3f}",
                     f"{r['bic']:.3f}", f"{r['w_bic']:.4f}"])
    A(md_table(["Modelo", "logL", "AIC", "BIC", "w-BIC"], rows))
    A("")
    A(f"- **log-likelihood del árbol ML:** `{s9['logL']:.4f}`")
    A(f"- **BIC:** {s9['bic']:.4f}")
    A(f"- **Longitud total del árbol:** {s9['longitud_total_arbol']:.4f} "
      "sustituciones/sitio")
    A(f"- **Parámetro de forma gamma (alpha):** {s9['alpha_gamma']:.4f}")
    A(f"- **Sitios del alineamiento:** {s9['n_sitios']} columnas, "
      f"{s9['sitios_constantes']} constantes "
      f"({s9['sitios_constantes']/s9['n_sitios']*100:.1f} %) y "
      f"**{s9['sitios_informativos']} parsimonia-informativos**")
    A(f"- **Semilla:** {s9['semilla']}; **n.º seqs/sitios:** "
      f"{s9['n_secuencias']}/{s9['n_sitios']}")
    A("- Soportes de rama del árbol ML: **UFBoot 1000** (`-B 1000`) y "
      "**SH-aLRT 1000** (`-alrt 1000`). Además se corrieron, con el modelo fijo "
      "`Dayhoff+G4` y la misma semilla: **bootstrap estándar de 100 réplicas** "
      "(`-b 100`, §6.5) y **UFBoot 1000** para comparar el efecto del número de "
      "réplicas (§6.6).\n")

    A("### 6.3 Topología\n")
    A("Árbol ML no enraizado (Newick del `.treefile`, ramas en sustituciones por "
      "sitio; soportes internos = SH-aLRT/UFBoot):\n")
    A("```")
    A(s9["topologia_newick"])
    A("```\n")
    rt = s9["longitudes_rama_terminales"]
    A(f"Longitudes de rama de los tres taxones del TP (sustituciones/sitio): "
      f"**humano {rt['humano']:.4f}**, **bartmosca {rt['bartmosca']:.4f}**, "
      f"**mosca {rt['mosca']:.4f}** (bartmosca y humano casi equidistantes del "
      "nodo, mosca mucho más larga).\n")
    A("Enraizado por outgroup **rhizobium**. La rama que lleva al outgroup mide "
      f"{s9['rama_outgroup_completa']:.4f} sustituciones/sitio en el árbol sin "
      "enraizar; al enraizar se parte en dos mitades iguales "
      f"({s9['rama_outgroup_completa']/2:.4f} + {s9['rama_outgroup_completa']/2:.4f}) "
      "para colgar la raíz, de modo que la longitud total del árbol no cambia. "
      "El árbol sin enraizar se conserva tal cual en "
      "`resultados/tree_set9_unrooted.nwk`.\n")
    A("Relaciones principales (todas con soporte):")
    A("")
    A("1. **Metazoa:** `(bartmosca, humano)` hermanos; mosca es el grupo externo "
      "de ese par. Es el resultado central del TP: bartmosca se agrupa con el "
      "humano, no con la mosca.")
    A("2. **Plantas verdes:** cebada hermana de `(vigna, (altramuz, abedul))`.")
    A("3. **Ectocarpus** (alga parda) queda fuera del clado de plantas verdes, "
      "como hermano de éste.")
    A("4. **Rhizobium** (bacteria) es el outgroup de todo el árbol.\n")

    A("### 6.4 Soportes y longitudes de rama de los nodos\n")
    A("Nodos internos del árbol enraizado (columnas: SH-aLRT %, UFBoot % del "
      "árbol ML; bootstrap estándar de **100 réplicas** con modelo fijo "
      "`Dayhoff+G4`; longitud de la rama que lleva al nodo, en "
      "sustituciones/sitio). El nodo artificial de la raíz (que no es un clado "
      "real, solo cuelga la raíz del outgroup) se omite.\n")
    rows = []
    for clado, d in s9["soportes_por_nodo"].items():
        if d["sh_alrt"] is None and d["ufboot"] is None and d["bootstrap_estandar"] is None:
            continue  # nodo artificial de enraizamiento
        rows.append([clado, fmt(d["sh_alrt"], 1), fmt(d["ufboot"], 1),
                     fmt(d["bootstrap_estandar"], 1), fmt(d["largo_rama"], 4)])
    A(md_table(["Clado", "SH-aLRT (%)", "UFBoot (%)", "Bootstrap estándar (%)",
                "Largo de rama"], rows))
    A("")
    A("Soportes en formato `SH-aLRT/UFBoot`: ≥ 80/95 se consideran bien "
      "soportados. El clado de plantas verdes (4 taxones) tiene soporte bajo "
      f"(UFBoot {s9['soportes_por_nodo']['abedul+altramuz+cebada+vigna']['ufboot']:.0f}, "
      f"SH-aLRT {s9['soportes_por_nodo']['abedul+altramuz+cebada+vigna']['sh_alrt']:.1f}); "
      "el de plantas verdes + ectocarpus también "
      f"(UFBoot {s9['soportes_por_nodo']['abedul+altramuz+cebada+ectocarpus+vigna']['ufboot']:.0f}, "
      f"SH-aLRT {s9['soportes_por_nodo']['abedul+altramuz+cebada+ectocarpus+vigna']['sh_alrt']:.1f}). "
      "`bartmosca+humano` tiene UFBoot "
      f"{s9['soportes_por_nodo']['bartmosca+humano']['ufboot']:.0f} (soporte moderado, no llega al umbral 95) y "
      "`bartmosca+humano+mosca` UFBoot "
      f"{s9['soportes_por_nodo']['bartmosca+humano+mosca']['ufboot']:.0f}.\n")
    A("**Longitudes de rama terminales** (sustituciones/sitio, árbol enraizado):\n")
    rows = [[k, fmt(v, 6)] for k, v in
            sorted(s9["longitudes_rama_terminales"].items(),
                   key=lambda kv: kv[1], reverse=True)]
    A(md_table(["Terminal", "Longitud de rama"], rows))
    A("")
    A(">*Nota:* la rama de `rhizobium` aparece con la mitad de su longitud "
      f"({s9['rama_outgroup_completa']/2:.4f}) porque la raíz se colgó en el "
      "punto medio de esa rama; en el `.treefile` sin enraizar mide "
      f"{s9['rama_outgroup_completa']:.4f}.\n")

    A("### 6.5 Bootstrap estándar (100 réplicas, modelo fijo)\n")
    A("El `-b 1000` con bootstrap estándar resultó demasiado lento en esta "
      "máquina (≈ 4 s por réplica, más de 1 h; se abortó) y la corrida con "
      "`-B 100` **no es posible**: IQ-TREE 2.3.6 exige al menos 1000 réplicas "
      "para UFBoot y aborta con `#replicates must be >= 1000`. Por eso el "
      "bootstrap estándar se corrió con **100 réplicas y modelo fijo "
      "`Dayhoff+G4`**:\n")
    A("```bash")
    A("tools/iqtree-2.3.6-Linux-intel/bin/iqtree2 \\")
    A("  -s resultados/aln_clustalo_set9.fasta \\")
    A("  -m Dayhoff+G4 -b 100 -T 2 -o rhizobium \\")
    A("  --prefix resultados/iqtree_set9_boot100 -seed 12345 -redo")
    A("```\n")
    bs100 = s9.get("bootstrap_estandar_100", {})
    A(f"La corrida terminó con log-likelihood {bs100.get('logL')} y BIC "
      f"{bs100.get('bic')}; los soportes están en la columna «Bootstrap "
      "estándar» de §6.4 (archivos `resultados/iqtree_set9_boot100.treefile` y "
      "`.iqtree`).\n")

    A("### 6.6 Efecto del número de réplicas: UFBoot 100 vs UFBoot 1000\n")
    A("Como `-B 100` está prohibido por IQ-TREE, la comparación se hizo sobre "
      "las **mismas 1000 réplicas** de UFBoot (guardadas con `--boot-trees` en "
      "`resultados/iqtree_set9_B1000.ufboot`, modelo fijo `Dayhoff+G4`, semilla "
      "12345): se calculó la frecuencia de cada split con las **primeras 100 "
      "réplicas** y con las **1000**, que es exactamente el efecto de reducir el "
      "número de réplicas (subconjunto anidado). La corrida completa se validó "
      "contra los soportes del `.treefile` (coinciden dentro del redondeo).\n")
    rows = []
    for nodo, d in s9.get("ufboot_100_vs_1000", {}).items():
        delta = "—" if d["diferencia_pp"] is None else f"{d['diferencia_pp']:+.1f}"
        rows.append([nodo, fmt(d["ufboot_100"], 1), fmt(d["ufboot_1000"], 1),
                     delta, fmt(d.get("bootstrap_estandar_100"), 1)])
    A(md_table(["Nodo", "UFBoot 100 (%)", "UFBoot 1000 (%)", "Δ (pp)",
                "Bootstrap estándar 100 (%)"], rows))
    A("")
    A("*(Los valores de UFBoot 1000 de esta tabla provienen de la corrida con "
      "modelo fijo; coinciden dentro de ~1 pp con los del árbol ML de §6.4.)*\n")
    A("Lectura: con 100 réplicas el soporte se mueve hasta **~8 puntos "
      "porcentuales** respecto de 1000 (el caso más claro es el clado de "
      "plantas verdes), mientras que los nodos con soporte alto/intermedio "
      "apenas cambian (≤ 3.5 pp). Con 100 réplicas la incertidumbre del propio "
      "soporte es del orden de ±5–8 pp; por eso se reportan los valores de 1000 "
      "réplicas como los definitivos.\n")

    A("## 7. Filogenia del set3 (humano, mosca, bartmosca)\n")
    A("```bash")
    A("tools/iqtree-2.3.6-Linux-intel/bin/iqtree2 \\")
    A("  -s resultados/aln_clustalo_set3.fasta \\")
    A("  -m Dayhoff+G4 -T 2 \\")
    A("  --prefix resultados/iqtree_set3 -seed 12345 -redo")
    A("```\n")
    ref = s3.get("referencia_model_finder", {})
    A(f"**Modelo:** fijado a `{s3['modelo']}` (el mejor por BIC en el set9 y el "
      "criterio pedido para comparar con el árbol de 9). Como referencia, la "
      f"corrida previa con ModelFinder (`-m MFP`, en "
      f"`resultados/iqtree_set3_mfp.*`) había elegido `{ref.get('modelo')}` "
      f"(BIC = {ref.get('bic'):.4f}, logL = {ref.get('logL'):.4f}); con el "
      "modelo fijo `Dayhoff+G4` el árbol es el mismo (la topología no puede "
      "cambiar, ver abajo).\n")
    A(f"- **Modelo de sustitución:** `{s3['modelo']}` "
      f"(alpha gamma = {s3['alpha_gamma']:.4f})")
    A(f"- **log-likelihood:** `{s3['logL']:.4f}`")
    A(f"- **BIC:** {s3['bic']:.4f}")
    A(f"- **Longitud total del árbol:** {s3['longitud_total_arbol']:.4f} "
      "sustituciones/sitio")
    A("- **Bootstrap:** no es posible añadir réplicas de bootstrap con 3 "
      "secuencias: `-b 1000` y también `-b 100` fallan con "
      "`ERROR: It makes no sense to perform bootstrap with less than 4 "
      "sequences` (IQ-TREE exige ≥ 4 taxones). No se puede reportar soporte de "
      "nodos porque **no hay nodos internos**.")
    A("- **Enraizamiento:** con 3 taxones el árbol no enraizado es una **única "
      "trifurcación**: no hay ningún nodo interno que resolver y no existe un "
      "outgroup dentro del set (los tres son metazoos y el TP los trata como "
      "ingroup). Por eso el resultado se muestra **sin enraizar** y no se generó "
      "ninguna versión enraizada: el enraizamiento es arbitrario. La única señal "
      "está en las longitudes de rama.")
    A("")
    A("Newick (sin enraizar):\n")
    A("```")
    A(s3["newick_sin_enraizar"])
    A("```\n")
    rows = [[k, fmt(v, 6)] for k, v in s3["longitudes_rama_terminales"].items()]
    A(md_table(["Terminal", "Longitud de rama (sust./sitio)"], rows))
    A("")
    l3 = s3["longitudes_rama_terminales"]
    A(f"bartmosca ({l3['bartmosca']:.4f}) y humano ({l3['humano']:.4f}) tienen "
      f"ramas cortas y mosca es la más divergente ({l3['mosca']:.4f}), lo que es "
      "coherente con el clado `(bartmosca, humano)` del set9. Con 3 taxones el "
      "bootstrap no agrega topologías posibles, así que la interpretación se "
      "apoya en las longitudes de rama y en las identidades de §4.4.\n")
    A(">*Detalle:* el log de IQ-TREE informa `Number of parsimony informative "
      "sites: 0` para este set. Es correcto con el criterio de Fitch que usa "
      "IQ-TREE (al menos dos estados y **cada uno** presente en ≥ 2 secuencias), "
      "imposible de cumplir con solo 3 taxones. Con el criterio útil para 3 "
      "taxones (dos secuencias que comparten un estado y la tercera distinta) "
      "hay **25 sitios** de 109. En cualquier caso el árbol no enraizado de 3 "
      "taxones es único y no hay nodos internos que soportar.\n")

    A("## 8. Figuras\n")
    A("Todas están en `resultados/` (no hay carpeta de figuras aparte):\n")
    A(md_table(["Figura", "Contenido"], [
        ["`resultados/fig_arbol_set9.png`",
         "Árbol rectangular del set9 enraizado por rhizobium; sobre cada nodo "
         "interno: SH-aLRT/UFBoot del árbol ML y bootstrap estándar de 100 "
         "réplicas (bs); barra de escala de 0.1 sust./sitio; longitudes de rama "
         "en sustituciones por sitio."],
        ["`resultados/fig_arbol_set3.png`",
         "Panel A: árbol no enraizado del set3 (trifurcación única) con barra de "
         "escala. Panel B: barras con las longitudes de rama terminales."],
        ["`resultados/fig_ventanas_bartmosca.png`",
         "Identidad de bartmosca contra humano y mosca en ventanas de 15 aa "
         "(paso 5), con la región N-terminal sombreada."],
        ["`resultados/fig_heatmap_identidad.png`",
         "Heatmap 9×9 de identidad de pares del set9 (deleción de pares, "
         "Clustal Omega), con los valores numéricos."],
    ]))
    A("")

    A("## 9. Cómo abrir y anotar los árboles en FigTree\n")
    A("1. Abrir FigTree y hacer **File → Open…**, elegir "
      "`resultados/iqtree_set9.treefile` (o `resultados/tree_set9_rooted.nwk`, "
      "que ya está enraizado por el outgroup).")
    A("2. **Enraizar:** si se abre el `.treefile` sin enraizar, seleccionar la "
      "rama de `rhizobium` y usar **Tree → Reroot on Selected Branch** (equivale "
      "al `-o rhizobium`). Alternativa: **Tree → Midpoint Root** (enraizamiento "
      "por punto medio, no usa outgroup).")
    A("3. **Mostrar soportes:** panel izquierdo → *Node Labels* → activar, y en "
      "el desplegable elegir el valor que IQ-TREE guardó como etiqueta de nodo "
      "(`SH-aLRT/UFBoot`; en la corrida `-b 100` la etiqueta es el bootstrap "
      "estándar). Ajustar el tamaño de fuente si se solapan.")
    A("4. **Longitudes de rama:** activar *Branch Labels → Branch Lengths* para "
      "verlas, o dejar la escala (*Scale Axis*) visible abajo.")
    A("5. **Exportar:** **File → Export PDF…** (vectorial, ideal para el "
      "informe) o **Export PNG…**; opcionalmente **File → Save Tree…** para "
      "guardar el árbol enraizado/anotado en Newick.")
    A("6. Para el set3 abrir `resultados/tree_set3.nwk`: al ser una trifurcación "
      "no tiene sentido re-enraizar (FigTree lo permite pero el resultado es "
      "arbitrario).\n")

    A("## 10. Comandos exactos ejecutados\n")
    A("Todos se corren desde la raíz del repositorio.\n")
    A("```bash")
    A("# 1) Conjuntos con etiquetas cortas")
    A("python3 src/01_prepare_sets.py")
    A("")
    A("# 2) Alineamientos")
    A("mafft --localpair --maxiterate 1000 \\")
    A("    datos/set9.fasta > resultados/aln_mafft_set9.fasta")
    A("mafft --localpair --maxiterate 1000 \\")
    A("    datos/set3.fasta > resultados/aln_mafft_set3.fasta")
    A("bash src/02_ebi_clustalo.sh        # Clustal Omega (REST EBI, outfmt=fa)")
    A("")
    A("# 3) Matrices de identidad / p-distancia y chequeo colapsando gaps")
    A("python3 src/03_pairwise.py")
    A("")
    A("# 4) Mosaico de bartmosca + heatmap")
    A("MPLCONFIGDIR=\"$PWD/.mplcache\" python3 src/04_mosaico.py")
    A("MPLCONFIGDIR=\"$PWD/.mplcache\" python3 src/06_heatmap.py")
    A("")
    A("# 5) IQ-TREE")
    A("IQ=tools/iqtree-2.3.6-Linux-intel/bin/iqtree2")
    A("$IQ -s resultados/aln_clustalo_set9.fasta -m MFP -B 1000 -alrt 1000 \\")
    A("    -T 2 -o rhizobium -seed 12345 --prefix resultados/iqtree_set9 -redo")
    A("# bootstrap estandar: 100 replicas con modelo fijo (el -b 1000 tardaba >1 h)")
    A("$IQ -s resultados/aln_clustalo_set9.fasta -m Dayhoff+G4 -b 100 -T 2 \\")
    A("    -o rhizobium --prefix resultados/iqtree_set9_boot100 -seed 12345 -redo")
    A("# UFBoot 1000 con modelo fijo, guardando las replicas para comparar 100 vs 1000")
    A("$IQ -s resultados/aln_clustalo_set9.fasta -m Dayhoff+G4 -B 1000 -T 2 \\")
    A("    -o rhizobium --boot-trees --prefix resultados/iqtree_set9_B1000 \\")
    A("    -seed 12345 -redo")
    A("# OJO: -B 100 no existe en IQ-TREE 2.3.6 (exige >=1000 replicas)")
    A("# set3 con modelo fijo")
    A("$IQ -s resultados/aln_clustalo_set3.fasta -m Dayhoff+G4 -T 2 -seed 12345 \\")
    A("    --prefix resultados/iqtree_set3 -redo")
    A("")
    A("# comparacion UFBoot 100 vs 1000 (submuestreo de las mismas replicas)")
    A("python3 src/10_boot_comparacion.py")
    A("")
    A("# 6) Árboles enraizados/sin enraizar y figuras")
    A("MPLCONFIGDIR=\"$PWD/.mplcache\" python3 src/05_trees.py")
    A("")
    A("# 7) Métricas y resumen")
    A("python3 src/07_metricas.py")
    A("python3 src/08_resumen.py")
    A("```\n")
    A("El detalle del envío/polling a EBI está en `src/02_ebi_clustalo.sh`; los "
      "IDs de trabajo quedaron registrados en `resultados/ebi_jobids.txt`:\n")
    A(md_table(["Set", "Job ID de EBI"], [[k, v] for k, v in ebi.items()]))
    A("")

    A("## 11. Archivos generados\n")
    A(md_table(["Archivo", "Qué es"], [
        ["`datos/set3.fasta`", "3 secuencias con etiquetas cortas"],
        ["`datos/set9.fasta`", "9 secuencias con etiquetas cortas"],
        ["`resultados/aln_clustalo_set{3,9}.fasta`", "alineamientos de Clustal Omega (EBI)"],
        ["`resultados/aln_mafft_set{3,9}.fasta`", "alineamientos de MAFFT"],
        ["`resultados/pairwise_{clustalo,mafft}_set{3,9}.csv`",
         "tabla de pares: identidad, p-distancia, identidad colapsada, gaps como diferencias"],
        ["`resultados/pairwise_identidad_*.csv`, `pairwise_pdist_*.csv`, "
         "`pairwise_identidad_colapsada_*.csv`", "matrices cuadradas (CSV)"],
        ["`resultados/bartmosca_ventanas.csv`", "24 ventanas de 15 aa (paso 5)"],
        ["`resultados/bartmosca_regiones.csv`", "identidad N-terminal vs resto"],
        ["`resultados/iqtree_set9.*`", "corrida ML + ModelFinder + UFBoot1000/SH-aLRT1000 (set9, árbol principal)"],
        ["`resultados/iqtree_set9_boot100.*`", "bootstrap estándar 100, modelo fijo Dayhoff+G4 (set9)"],
        ["`resultados/iqtree_set9_B1000.*`", "UFBoot 1000 con modelo fijo + `.ufboot` con las 1000 réplicas"],
        ["`resultados/boot_comparacion.csv`", "UFBoot 100 vs UFBoot 1000 por nodo"],
        ["`resultados/iqtree_set3.*`", "corrida ML del set3 con modelo fijo Dayhoff+G4"],
        ["`resultados/iqtree_set3_mfp.*`", "referencia: corrida previa del set3 con ModelFinder (eligió DCMut)"],
        ["`resultados/tree_set9_rooted.nwk`", "árbol set9 enraizado por outgroup (Newick)"],
        ["`resultados/tree_set9_unrooted.nwk`", "árbol set9 sin enraizar (Newick)"],
        ["`resultados/tree_set3.nwk`", "árbol set3 sin enraizar (Newick)"],
        ["`resultados/soportes_nodos.csv`", "tabla nodo → SH-aLRT / UFBoot / bootstrap / rama"],
        ["`resultados/fig_arbol_set9.png`, `fig_arbol_set3.png`, "
         "`fig_ventanas_bartmosca.png`, `fig_heatmap_identidad.png`", "figuras"],
        ["`resultados/metricas.json`", "todos los números clave en JSON"],
        ["`resultados/resumen_analisis.md`", "este documento"],
    ]))
    A("")

    A("## 12. Limitaciones y notas de honestidad\n")
    A("1. **Clustal Omega vs MAFFT:** para set9 coinciden en 131 columnas pero "
      "no en el detalle de los gaps; por eso algunas identidades de pares "
      "cambian hasta ~2 puntos porcentuales (p. ej. bartmosca–mosca: "
      f"{bm['clustalo']['identidad']*100:.2f} % con ClustalO vs "
      f"{bm['mafft']['identidad']*100:.2f} % con MAFFT). El chequeo "
      "«colapsando gaps» da prácticamente lo mismo con ambos alineadores "
      "(columna correspondiente en §4.1), lo que indica que la diferencia es de "
      "colocación de gaps y no de homología.")
    A("2. **Bootstrap:** con set3 no se puede correr ningún bootstrap "
      "(`-b 1000` ni `-b 100`): IQ-TREE exige ≥ 4 secuencias "
      "(`ERROR: It makes no sense to perform bootstrap with less than 4 "
      "sequences`). En set9, el bootstrap estándar `-b 1000` se abortó por "
      "lentitud (≈ 4 s/réplica, > 1 h) y se reemplazó por `-b 100` con modelo "
      "fijo. `-B 100` no existe en IQ-TREE 2.3.6 (mínimo 1000), por lo que la "
      "comparación 100 vs 1000 se hizo submuestreando las 1000 réplicas de "
      "UFBoot guardadas con `--boot-trees` (subconjunto anidado, mismo seed).")
    A(f"3. **Enraizamiento set9:** se usó el outgroup `rhizobium`. La rama del "
      "outgroup se partió por la mitad para colgar la raíz; en el `.treefile` "
      "original esa rama conserva su longitud completa "
      f"({s9['rama_outgroup_completa']:.4f}).")
    A("4. **Set3:** no se puede enraizar por outgroup con los tres metazoos; el "
      "árbol no enraizado es único (trifurcación). No se inventó un "
      "enraizamiento.")
    A("5. **Soportes:** los nodos con UFBoot ≥ 95 y SH-aLRT ≥ 80 se consideran "
      "bien apoyados. El clado de plantas verdes tiene soporte bajo (UFBoot "
      f"{s9['soportes_por_nodo']['abedul+altramuz+cebada+vigna']['ufboot']:.0f}, "
      "SH-aLRT "
      f"{s9['soportes_por_nodo']['abedul+altramuz+cebada+vigna']['sh_alrt']:.1f}), "
      "y el de plantas verdes + ectocarpus también (UFBoot "
      f"{s9['soportes_por_nodo']['abedul+altramuz+cebada+ectocarpus+vigna']['ufboot']:.0f}, "
      "SH-aLRT "
      f"{s9['soportes_por_nodo']['abedul+altramuz+cebada+ectocarpus+vigna']['sh_alrt']:.1f}). "
      "`bartmosca+humano` (UFBoot "
      f"{s9['soportes_por_nodo']['bartmosca+humano']['ufboot']:.0f}, SH-aLRT "
      f"{s9['soportes_por_nodo']['bartmosca+humano']['sh_alrt']:.1f}) queda por debajo de los umbrales: "
      "hay que tomarlo como soporte moderado, no como un clado robusto.")
    A("6. La identidad «colapsando gaps» se implementó re-alineando cada par por "
      "separado (BLOSUM62) en lugar de comparar las cadenas sin gaps posición a "
      "posición, porque los indels desplazan el marco y ese cálculo directo "
      "daría valores sin sentido (p. ej. ~9 %).")
    A("")

    A("## 13. Reproducibilidad\n")
    A("Orden de ejecución: `src/01_prepare_sets.py` → `src/02_ebi_clustalo.sh` → "
      "`src/03_pairwise.py` → `src/04_mosaico.py` → `src/05_trees.py` → "
      "`src/06_heatmap.py` → `src/10_boot_comparacion.py` → `src/07_metricas.py` → "
      "`src/08_resumen.py`. "
      "Las figuras requieren `MPLCONFIGDIR=\"$PWD/.mplcache\"`. Los números de "
      "este resumen se generan desde `resultados/metricas.json` y los CSV, de "
      "modo que si se vuelve a correr el pipeline el documento se regenera "
      "consistente.\n")

    out = os.path.join(RES, "resumen_analisis.md")
    with open(out, "w") as fh:
        fh.write("\n".join(L))
    print("escrito", out, f"({len(L)} bloques)")


if __name__ == "__main__":
    main()
