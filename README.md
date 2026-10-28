# TP de Inferencias Evolutivas — el caso del bartmosca

Resolución de los desafíos I, II y III del TP de *Inferencias Evolutivas* (Introducción a la Bioinformática), con todo el flujo: alineamientos (Clustal Omega y MAFFT), distancias, análisis de mosaico, filogenia con IQ-TREE, bootstrap y enraizamiento.

**El resultado en dos líneas:** la secuencia del "moscahumano" es 92,4% idéntica al citocromo c humano y sus últimos 82 residuos son idénticos al humano; el árbol la agrupa con el humano y no con la mosca (el par bartmosca–humano tiene soporte moderado: UFBoot 75, SH-aLRT 77, bootstrap 59; el clado de metazoos tiene UFBoot 86). La parte "mosca" son 3 residuos (V14, Q15, A18) del extremo N-terminal, en una región que ni siquiera coincide exactamente con la mosca de referencia.

## Contenido

- [`informe_desafio_I.md`](informe_desafio_I.md): criterios de comparación, distancias, análisis de mosaico y respuesta a la pregunta de qué pasaría al usar el resto de las secuencias.
- [`informe_desafio_II.md`](informe_desafio_II.md): alineamiento con Clustal Omega (REST de EBI) y control con MAFFT, controles de calidad del alineamiento y el problema del N-terminal.
- [`informe_desafio_III.md`](informe_desafio_III.md): árboles con IQ-TREE, modelo de sustitución, bootstrap, enraizamiento e interpretación con FigTree.
- `datos/`: fasta original del TP y los subconjuntos usados (`set3.fasta`, `set9.fasta`).
- `src/`: scripts que reproducen todo (alineamientos, matrices de distancia, mosaico, árboles, figuras).
- `resultados/`: alineamientos, tablas, árboles (`.treefile`, `.nwk`), figuras y el resumen técnico de la corrida.
- `resultados/resumen_analisis.md`: detalle técnico completo (comandos, parámetros, versiones, soportes por nodo, limitaciones).

## Cómo reproducir

Requiere Python 3 con Biopython, pandas, numpy y matplotlib, [MAFFT](https://mafft.cbrc.jp/alignment/software/) y [IQ-TREE 2](http://www.iqtree.org/). El alineamiento con Clustal Omega se hizo contra el servicio REST de EBI, así que ese paso necesita conexión.

```bash
# dependencias de Python
pip install -r requirements.txt

# IQ-TREE (binario, no se versiona por tamaño)
mkdir -p tools && curl -sL https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz \
  | tar xz -C tools

# pipeline
python3 src/01_prepare_sets.py          # datos/set3.fasta y datos/set9.fasta
bash src/02_ebi_clustalo.sh             # Clustal Omega vía EBI
mafft --localpair --maxiterate 1000 datos/set9.fasta > resultados/aln_mafft_set9.fasta
mafft --localpair --maxiterate 1000 datos/set3.fasta > resultados/aln_mafft_set3.fasta
python3 src/03_pairwise.py              # identidad y p-distancia por pares
MPLCONFIGDIR="$PWD/.mplcache" python3 src/04_mosaico.py    # ventanas deslizantes
MPLCONFIGDIR="$PWD/.mplcache" python3 src/09_mosaico_residuos.py  # mosaico residuo por residuo
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
MPLCONFIGDIR="$PWD/.mplcache" python3 src/05_trees.py      # árboles, soportes y figuras
MPLCONFIGDIR="$PWD/.mplcache" python3 src/06_heatmap.py
python3 src/10_boot_comparacion.py      # comparación UFBoot 100 vs 1000
python3 src/07_metricas.py              # resultados/metricas.json
python3 src/08_resumen.py               # resultados/resumen_analisis.md
```

Los árboles se pueden abrir y anotar en [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) a partir de `resultados/iqtree_set9.treefile` (enraizar en la rama de `rhizobium`, activar *Node Labels*, exportar). Las figuras del repo se generaron desde esos mismos archivos con Biopython y matplotlib.

## Enunciado

El enunciado original está en el repositorio del curso: [Inferencias evolutivas](https://github.com/AJVelezRueda/Introduccion_a_la_Bioinformatica/blob/master/Teorico_Practicos/Inferencias_evolutivas/inferencias_evolutivas.md). Las secuencias de citocromo c son las que acompaña el TP: [`SecuenciasCytocromoC.fasta`](https://raw.githubusercontent.com/AJVelezRueda/Introduccion_a_la_Bioinformatica/master/Teorico_Practicos/Inferencias_evolutivas/SecuenciasCytocromoC.fasta).

## Notas

- El enunciado menciona "Cyt P450", pero el archivo de secuencias y la secuencia problema son de **citocromo c**. Queda aclarado en el [informe del Desafío I](informe_desafio_I.md); el análisis se hizo sobre el marcador que efectivamente traen los datos.
- Las figuras de los árboles se generaron con Biopython y matplotlib a partir de los mismos `.treefile`; los archivos quedan listos para abrirlos y anotarlos en FigTree.
- `tools/` no se versiona (binario de IQ-TREE); se descarga con el comando de arriba.

## Fuente y licencia

Trabajo práctico basado en el material de **Ana Julia Velez Rueda y Guillermo I. Benítez** (Introducción a la Bioinformática). El enunciado y las secuencias originales están bajo licencia [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/); esta resolución se distribuye bajo la misma licencia.
