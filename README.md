# TP de Inferencias Evolutivas

En este trabajo práctico analicé la secuencia problema del caso del "moscahumano" y la comparé con secuencias de citocromo c de distintos organismos. Resolví los desafíos I, II y III: calculé distancias, revisé si la secuencia podía ser un mosaico, comparé alineamientos y construí árboles filogenéticos con IQ-TREE.

La secuencia problema tiene 92,4% de identidad con el citocromo c humano. Sus últimos 82 residuos son idénticos al humano; las coincidencias con la mosca aparecen en tres residuos del extremo N-terminal (V14, Q15 y A18). En el árbol queda agrupada con el humano, aunque el soporte de ese par es moderado (UFBoot 75, SH-aLRT 77 y bootstrap estándar 59). El análisis y sus limitaciones se describen en los informes.

## Archivos principales

- `datos/`: secuencias originales y subconjuntos usados (`set3.fasta`, `set9.fasta`).
- `src/`: scripts para preparar los conjuntos, alinear las secuencias, calcular distancias, analizar el mosaico y generar árboles y figuras.
- `resultados/`: alineamientos, tablas, árboles, figuras y el resumen técnico de la corrida.
- `informe_desafio_I.md`: comparación de secuencias, distancias, análisis de mosaico y efecto de incluir el resto de las secuencias.
- `informe_desafio_II.md`: alineamientos de Clustal Omega y MAFFT, controles de calidad y revisión del extremo N-terminal.
- `informe_desafio_III.md`: árboles con IQ-TREE, modelos de sustitución, bootstrap, enraizamiento e interpretación.
- `resultados/resumen_analisis.md`: comandos, parámetros, versiones, soportes por nodo y limitaciones de la corrida.

## Cómo reproducir

Se requiere Python 3 con Biopython, pandas, numpy y matplotlib, además de [MAFFT](https://mafft.cbrc.jp/alignment/software/) e [IQ-TREE 2](http://www.iqtree.org/). El alineamiento con Clustal Omega se hizo mediante el servicio REST de EBI y ese paso requiere conexión a internet.

```bash
pip install -r requirements.txt

# IQ-TREE se descarga aparte porque el binario no está versionado
mkdir -p tools && curl -sL https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz \
  | tar xz -C tools

python3 src/01_prepare_sets.py
bash src/02_ebi_clustalo.sh
mafft --localpair --maxiterate 1000 datos/set9.fasta > resultados/aln_mafft_set9.fasta
mafft --localpair --maxiterate 1000 datos/set3.fasta > resultados/aln_mafft_set3.fasta
python3 src/03_pairwise.py
MPLCONFIGDIR="$PWD/.mplcache" python3 src/04_mosaico.py
MPLCONFIGDIR="$PWD/.mplcache" python3 src/09_mosaico_residuos.py

IQ=tools/iqtree-2.3.6-Linux-intel/bin/iqtree2
$IQ -s resultados/aln_clustalo_set9.fasta -m MFP -B 1000 -alrt 1000 -T 2 -o rhizobium \
    -seed 12345 --prefix resultados/iqtree_set9 -redo
$IQ -s resultados/aln_clustalo_set9.fasta -m Dayhoff+G4 -b 100 -T 2 -o rhizobium \
    -seed 12345 --prefix resultados/iqtree_set9_boot100 -redo
$IQ -s resultados/aln_clustalo_set9.fasta -m DCMut+G4 -T 2 -o rhizobium \
    -seed 12345 --prefix resultados/iqtree_set9_dcmut -redo
$IQ -s resultados/aln_clustalo_set9.fasta -m Dayhoff+G4 -B 1000 -T 2 -o rhizobium --boot-trees \
    -seed 12345 --prefix resultados/iqtree_set9_B1000 -redo
$IQ -s resultados/aln_clustalo_set3.fasta -m Dayhoff+G4 -T 2 -seed 12345 \
    --prefix resultados/iqtree_set3 -redo

MPLCONFIGDIR="$PWD/.mplcache" python3 src/05_trees.py
MPLCONFIGDIR="$PWD/.mplcache" python3 src/06_heatmap.py
python3 src/10_boot_comparacion.py
python3 src/07_metricas.py
python3 src/08_resumen.py
```

Los árboles se pueden abrir en [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) desde `resultados/iqtree_set9.treefile`. Para mostrarlo enraizado, se selecciona la rama de `rhizobium` y se activa *Node Labels*. Las figuras del repositorio se generaron con Biopython y matplotlib a partir de los mismos archivos.

## Enunciado y datos

El enunciado está en el repositorio del curso: [Inferencias evolutivas](https://github.com/AJVelezRueda/Introduccion_a_la_Bioinformatica/blob/master/Teorico_Practicos/Inferencias_evolutivas/inferencias_evolutivas.md). El archivo original de secuencias es [`SecuenciasCytocromoC.fasta`](https://raw.githubusercontent.com/AJVelezRueda/Introduccion_a_la_Bioinformatica/master/Teorico_Practicos/Inferencias_evolutivas/SecuenciasCytocromoC.fasta).

El enunciado menciona "Cyt P450", pero el archivo y la secuencia problema corresponden a citocromo c. Dejé la aclaración en el [informe del Desafío I](informe_desafio_I.md); el análisis se hizo con las secuencias que acompañan el TP.

## Fuente y licencia

El trabajo práctico se basa en el material de Ana Julia Velez Rueda y Guillermo I. Benítez, de Introducción a la Bioinformática. El enunciado y las secuencias originales están bajo licencia [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/); esta resolución se distribuye bajo la misma licencia.
