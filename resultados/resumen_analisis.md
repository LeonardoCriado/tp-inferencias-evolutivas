# Análisis de inferencias evolutivas sobre citocromo c (TP BartMosca)

**Fecha de ejecución:** 2026-10-27 (America/Buenos_Aires)  
**Insumo:** `datos/SecuenciasCytocromoC.fasta` (9 secuencias, citocromo c, 104–127 aa)  
**Resultados crudos:** `resultados/` · **Scripts:** `src/` · **Datos derivados:** `datos/`  
**Semilla de IQ-TREE:** `-seed 12345` (en set9 y set3) para que el análisis sea exactamente reproducible.

## 1. Versiones de software

| Herramienta | Versión / origen |
|---|---|
| IQ-TREE | 2.3.6 (Linux x86 64-bit, 2024-08-01) |
| MAFFT | v7.526 (2024/Apr/26) |
| Clustal Omega | Clustal Omega via REST de EBI (https://www.ebi.ac.uk/Tools/services/rest/clustalo) |
| Biopython | 1.88 |
| Python | 3.14.8 |
| NumPy | 2.5.3 |
| pandas | 3.0.6 |
| matplotlib | 3.11.2 |

## 2. Secuencias y conjuntos de datos

Etiquetas cortas (sin espacios, para que IQ-TREE no se queje) y acceso original de cada secuencia:

| Etiqueta | Header (recortado) | Especie / fuente | Largo crudo (aa) |
|---|---|---|---|
| humano | NP_061820.1_cytochrome_c | Homo_sapiens | 105 |
| mosca | NP_001285984.1_cytochrome_c_proximal,_isoform_B | Drosophila_melanogaster | 108 |
| bartmosca | bartmosca | — | 107 |
| cebada | KAE8767936.1_Cytochrome_c | Hordeum_vulgare | 112 |
| altramuz | KAE9612850.1_Cytochrome_c | Lupinus_albus | 113 |
| abedul | QEE59979.1_cytochrome_c | Betula_platyphylla | 113 |
| vigna | NP_001363126.1_cytochrome_c | Vigna_radiata | 112 |
| ectocarpus | CBJ31344.1_Cytochrome_c | Ectocarpus_siliculosus | 104 |
| rhizobium | WP_037127040.1_cytochrome_c_family_protein | Rhizobium_sp._CF097 | 127 |

- `datos/set3.fasta`: **humano, mosca, bartmosca** (3 secuencias).
- `datos/set9.fasta`: las **9** secuencias (bartmosca + 2 metazoos + 6 de referencia: cebada *Hordeum vulgare*, altramuz *Lupinus albus*, abedul *Betula platyphylla*, vigna *Vigna radiata*, ectocarpus *Ectocarpus siliculosus* y rhizobium *Rhizobium* sp. CF097).

## 3. Alineamientos

Se alineó con **Clustal Omega** (REST de EBI, `stype=protein`, `outfmt=fa`) y con **MAFFT 7.526** (`--localpair --maxiterate 1000`). IDs de los trabajos de EBI: set3 = `clustalo-R20261028-001614-0565-18522630-p1m`, set9 = `clustalo-R20261028-001625-0783-85075839-p1m`.

| Alineador | Set | N.º seqs | Columnas | Columnas sin gap | Gaps totales |
|---|---|---|---|---|---|
| Clustal Omega | set3 | 3 | 109 | 104 | 7 |
| MAFFT | set3 | 3 | 109 | 104 | 7 |
| Clustal Omega | set9 | 9 | 131 | 100 | 178 |
| MAFFT | set9 | 9 | 131 | 100 | 178 |

Los dos alineadores dan **la misma longitud** para cada set (set3 = 109 columnas, set9 = 131 columnas). No colocan exactamente los mismos gaps: en el set9, 16 de las 131 columnas tienen algún carácter distinto entre ambos alineamientos (ver §4 y §12).

## 4. Identidad y p-distancia de pares (set9)

**Definiciones usadas** (todas sobre el alineamiento de Clustal Omega):

- **Identidad (deleción de pares):** proporción de posiciones idénticas considerando solo las columnas donde *ambas* secuencias del par tienen residuo (`idénticas / (idénticas + distintas)`). Es la identidad estándar.
- **p-distancia:** `1 − identidad` (no corregida por múltiples sustituciones).
- **Colapsando gaps (chequeo):** se eliminan los gaps del alineamiento, se re-alinea cada par de forma independiente (`PairwiseAligner` global, BLOSUM62, gap open −10 / extend −0.5) y se recalcula la identidad con deleción de pares. Si difiere de la anterior, el valor depende del alineador múltiple.
- **Gaps como diferencias:** `idénticas / columnas totales del alineamiento` (cota inferior, útil para comparar con la literatura).

### 4.1 Pares clave

| Par | Ident. ClustalO | p-dist. ClustalO | Ident. colapsada | Ident. gaps=dif | Ident. MAFFT | p-dist. MAFFT | Ident. colapsada MAFFT |
|---|---|---|---|---|---|---|---|
| bartmosca-humano | 0.9238 | 0.0762 | 0.9238 | 0.7405 | 0.9333 | 0.0667 | 0.9238 |
| bartmosca-mosca | 0.7547 | 0.2453 | 0.7736 | 0.6107 | 0.7736 | 0.2264 | 0.7736 |
| humano-mosca | 0.7692 | 0.2308 | 0.7692 | 0.6107 | 0.7788 | 0.2212 | 0.7692 |
| cebada-humano | 0.6154 | 0.3846 | 0.6154 | 0.4885 | 0.6250 | 0.3750 | 0.6154 |
| altramuz-abedul | 0.9381 | 0.0619 | 0.9381 | 0.8092 | 0.9381 | 0.0619 | 0.9381 |
| vigna-abedul | 0.9464 | 0.0536 | 0.9464 | 0.8092 | 0.9464 | 0.0536 | 0.9464 |
| cebada-ectocarpus | 0.6796 | 0.3204 | 0.6796 | 0.5344 | 0.6893 | 0.3107 | 0.6796 |

**Lectura:** bartmosca es ~92.4 % idéntica al humano (p-dist 0.0762) pero solo ~75.5 % a la mosca; humano–mosca comparten ~76.9 %. Es decir, bartmosca está mucho más cerca del humano que la propia mosca: la señal filogenética de citocromo c la agrupa con el humano, no con la mosca.

### 4.2 Matriz completa de identidad (set9, Clustal Omega, %)

|  | humano | bartmosca | mosca | cebada | altramuz | abedul | vigna | ectocarpus | rhizobium |
|---|---|---|---|---|---|---|---|---|---|
| humano | — | 92.4 | 76.9 | 61.5 | 62.9 | 61.9 | 63.5 | 60.6 | 53.5 |
| bartmosca | 92.4 | — | 75.5 | 61.3 | 61.7 | 59.8 | 61.3 | 57.7 | 53.4 |
| mosca | 76.9 | 75.5 | — | 61.1 | 58.3 | 57.4 | 57.4 | 56.3 | 49.5 |
| cebada | 61.5 | 61.3 | 61.1 | — | 83.9 | 83.9 | 84.8 | 68.0 | 51.4 |
| altramuz | 62.9 | 61.7 | 58.3 | 83.9 | — | 93.8 | 95.5 | 71.2 | 48.6 |
| abedul | 61.9 | 59.8 | 57.4 | 83.9 | 93.8 | — | 94.6 | 70.2 | 47.7 |
| vigna | 63.5 | 61.3 | 57.4 | 84.8 | 95.5 | 94.6 | — | 71.8 | 46.8 |
| ectocarpus | 60.6 | 57.7 | 56.3 | 68.0 | 71.2 | 70.2 | 71.8 | — | 55.0 |
| rhizobium | 53.5 | 53.4 | 49.5 | 51.4 | 48.6 | 47.7 | 46.8 | 55.0 | — |

### 4.3 Matriz de p-distancia (set9, Clustal Omega)

|  | humano | bartmosca | mosca | cebada | altramuz | abedul | vigna | ectocarpus | rhizobium |
|---|---|---|---|---|---|---|---|---|---|
| humano | — | 0.0762 | 0.2308 | 0.3846 | 0.3714 | 0.3810 | 0.3654 | 0.3942 | 0.4653 |
| bartmosca | 0.0762 | — | 0.2453 | 0.3868 | 0.3832 | 0.4019 | 0.3868 | 0.4231 | 0.4660 |
| mosca | 0.2308 | 0.2453 | — | 0.3889 | 0.4167 | 0.4259 | 0.4259 | 0.4369 | 0.5048 |
| cebada | 0.3846 | 0.3868 | 0.3889 | — | 0.1607 | 0.1607 | 0.1518 | 0.3204 | 0.4862 |
| altramuz | 0.3714 | 0.3832 | 0.4167 | 0.1607 | — | 0.0619 | 0.0446 | 0.2885 | 0.5138 |
| abedul | 0.3810 | 0.4019 | 0.4259 | 0.1607 | 0.0619 | — | 0.0536 | 0.2981 | 0.5229 |
| vigna | 0.3654 | 0.3868 | 0.4259 | 0.1518 | 0.0446 | 0.0536 | — | 0.2816 | 0.5321 |
| ectocarpus | 0.3942 | 0.4231 | 0.4369 | 0.3204 | 0.2885 | 0.2981 | 0.2816 | — | 0.4500 |
| rhizobium | 0.4653 | 0.4660 | 0.5048 | 0.4862 | 0.5138 | 0.5229 | 0.5321 | 0.4500 | — |

### 4.4 Set3 (humano, mosca, bartmosca)

| Par | Identidad | p-distancia | Ident. colapsada | Posiciones comparadas | Columnas alineamiento |
|---|---|---|---|---|---|
| bartmosca-humano | 0.9238 | 0.0762 | 0.9238 | 105 | 109 |
| bartmosca-mosca | 0.7547 | 0.2453 | 0.7736 | 106 | 109 |
| humano-mosca | 0.7692 | 0.2308 | 0.7692 | 104 | 109 |

## 5. Análisis de mosaico de bartmosca

Ventanas deslizantes de **15 aa con paso 5** sobre el alineamiento Clustal Omega del set9; en cada ventana se calcula la identidad de bartmosca contra humano y contra mosca (deleción de pares). Las primeras 2 ventanas (columnas 1–20) no tienen posiciones comparables porque bartmosca/humano/mosca tienen gaps en ese tramo del alineamiento (ver `resultados/bartmosca_ventanas.csv`).

### 5.1 Resumen por regiones

| Región | Columnas del alineamiento | Ident. vs humano | Ident. vs mosca | Δ (humano − mosca) |
|---|---|---|---|---|
| N-terminal (1-25) | 1-46 | 0.6522 | 0.6400 | 0.0122 |
| resto (26-fin) | 47-131 | 1.0000 | 0.7901 | 0.2099 |

El N-terminal (primeros 25 residuos de bartmosca; columnas 1–46) es la única región con ambigüedad: identidad 65.2 % vs humano y 64.0 % vs mosca. El resto de la proteína (residuos 26–107) es **idéntico al humano al 100 %** y solo ~79 % a la mosca. Conclusión: bartmosca no es un mosaico con el N-terminal de mosca; el N-terminal es *levemente* más parecido al humano y queda como la única zona candidata a aporte diferencial.

### 5.2 Ventanas (15 aa, paso 5)

| Ventana (col.) | Ident. vs humano (%) | Ident. vs mosca (%) | Δ humano−mosca | Más parecida a |
|---|---|---|---|---|
| 1-15 | — | — | — | sin posiciones comparables |
| 6-20 | — | — | — | sin posiciones comparables |
| 11-25 | 50.0 | 25.0 | +25.0 | humano |
| 16-30 | 57.1 | 44.4 | +12.7 | humano |
| 21-35 | 66.7 | 57.1 | +9.5 | humano |
| 26-40 | 66.7 | 73.3 | -6.7 | mosca |
| 31-45 | 73.3 | 80.0 | -6.7 | mosca |
| 36-50 | 73.3 | 80.0 | -6.7 | mosca |
| 41-55 | 86.7 | 80.0 | +6.7 | humano |
| 46-60 | 93.3 | 80.0 | +13.3 | humano |
| 51-65 | 100.0 | 86.7 | +13.3 | humano |
| 56-70 | 100.0 | 80.0 | +20.0 | humano |
| 61-75 | 100.0 | 73.3 | +26.7 | humano |
| 66-80 | 100.0 | 61.5 | +38.5 | humano |
| 71-85 | 100.0 | 66.7 | +33.3 | humano |
| 76-90 | 100.0 | 75.0 | +25.0 | humano |
| 81-95 | 100.0 | 78.6 | +21.4 | humano |
| 86-100 | 100.0 | 86.7 | +13.3 | humano |
| 91-105 | 100.0 | 93.3 | +6.7 | humano |
| 96-110 | 100.0 | 93.3 | +6.7 | humano |
| 101-115 | 100.0 | 80.0 | +20.0 | humano |
| 106-120 | 100.0 | 66.7 | +33.3 | humano |
| 111-125 | 100.0 | 73.3 | +26.7 | humano |
| 116-130 | 100.0 | 73.3 | +26.7 | humano |

De las 22 ventanas con datos, 19 son más parecidas al humano y 3 a la mosca. La ventana con máxima ventaja del humano es **66-80** (Δ = +38.5 puntos porcentuales).

## 6. Filogenia del set9 (9 secuencias, 131 sitios)

### 6.1 Comando

```bash
tools/iqtree-2.3.6-Linux-intel/bin/iqtree2 \
  -s resultados/aln_clustalo_set9.fasta \
  -m MFP -B 1000 -alrt 1000 -T 2 -o rhizobium \
  -seed 12345 --prefix resultados/iqtree_set9 -redo
```

### 6.2 Selección de modelo (ModelFinder, criterio BIC)

**Modelo elegido: `Dayhoff+G4`** (BIC = 2375.9237). Es Dayhoff con heterogeneidad de tasas gamma (4 categorías). Primeros puestos de la tabla de ModelFinder (dos modelos quedan casi empatados: w-BIC = 0.411 para Dayhoff+G4 y 0.385 para DCMut+G4):

| Modelo | logL | AIC | BIC | w-BIC |
|---|---|---|---|---|
| Dayhoff+G4 | -1148.960 | 2329.921 | 2375.924 | 0.4110 |
| DCMut+G4 | -1149.027 | 2330.053 | 2376.056 | 0.3850 |
| Dayhoff+R2 | -1148.398 | 2330.796 | 2379.675 | 0.0630 |
| DCMut+R2 | -1148.450 | 2330.901 | 2379.779 | 0.0598 |
| Dayhoff+I+G4 | -1149.144 | 2332.288 | 2381.167 | 0.0299 |
| DCMut+I+G4 | -1149.256 | 2332.512 | 2381.390 | 0.0267 |

- **log-likelihood del árbol ML:** `-1148.9603`
- **BIC:** 2375.9237
- **Longitud total del árbol:** 2.4260 sustituciones/sitio
- **Parámetro de forma gamma (alpha):** 1.1530
- **Sitios del alineamiento:** 131 columnas, 55 constantes (42.0 %) y **44 parsimonia-informativos**
- **Semilla:** 12345; **n.º seqs/sitios:** 9/131
- Soportes de rama del árbol ML: **UFBoot 1000** (`-B 1000`) y **SH-aLRT 1000** (`-alrt 1000`). Además se corrieron, con el modelo fijo `Dayhoff+G4` y la misma semilla: **bootstrap estándar de 100 réplicas** (`-b 100`, §6.5) y **UFBoot 1000** para comparar el efecto del número de réplicas (§6.6).

### 6.3 Topología

Árbol ML no enraizado (Newick del `.treefile`, ramas en sustituciones por sitio; soportes internos = SH-aLRT/UFBoot):

```
(((humano:0.0443950531,bartmosca:0.0496997480)77.2/75:0.0846090960,mosca:0.2257317127)88.2/86:0.3134132475,((cebada:0.1016766776,((altramuz:0.0222977052,abedul:0.0419460216)75.8/71:0.0208774470,vigna:0.0068399353)85.9/76:0.0578497744)28.1/45:0.0940594585,ectocarpus:0.3013481749)35.8/56:0.1878695972,rhizobium:0.8733520988);
```

Longitudes de rama de los tres taxones del TP (sustituciones/sitio): **humano 0.0444**, **bartmosca 0.0497**, **mosca 0.2257** (bartmosca y humano casi equidistantes del nodo, mosca mucho más larga).

Enraizado por outgroup **rhizobium**. La rama que lleva al outgroup mide 0.8734 sustituciones/sitio en el árbol sin enraizar; al enraizar se parte en dos mitades iguales (0.4367 + 0.4367) para colgar la raíz, de modo que la longitud total del árbol no cambia. El árbol sin enraizar se conserva tal cual en `resultados/tree_set9_unrooted.nwk`.

Relaciones principales (todas con soporte):

1. **Metazoa:** `(bartmosca, humano)` hermanos; mosca es el grupo externo de ese par. Es el resultado central del TP: bartmosca se agrupa con el humano, no con la mosca.
2. **Plantas verdes:** cebada hermana de `(vigna, (altramuz, abedul))`.
3. **Ectocarpus** (alga parda) queda fuera del clado de plantas verdes, como hermano de éste.
4. **Rhizobium** (bacteria) es el outgroup de todo el árbol.

### 6.4 Soportes y longitudes de rama de los nodos

Nodos internos del árbol enraizado (columnas: SH-aLRT %, UFBoot % del árbol ML; bootstrap estándar de **100 réplicas** con modelo fijo `Dayhoff+G4`; longitud de la rama que lleva al nodo, en sustituciones/sitio). El nodo artificial de la raíz (que no es un clado real, solo cuelga la raíz del outgroup) se omite.

| Clado | SH-aLRT (%) | UFBoot (%) | Bootstrap estándar (%) | Largo de rama |
|---|---|---|---|---|
| abedul+altramuz+cebada+ectocarpus+vigna | 35.8 | 56.0 | 62.0 | 0.1879 |
| abedul+altramuz+cebada+vigna | 28.1 | 45.0 | 29.0 | 0.0941 |
| abedul+altramuz+vigna | 85.9 | 76.0 | 84.0 | 0.0578 |
| bartmosca+humano+mosca | 88.2 | 86.0 | 81.0 | 0.3134 |
| bartmosca+humano | 77.2 | 75.0 | 59.0 | 0.0846 |
| abedul+altramuz | 75.8 | 71.0 | 70.0 | 0.0209 |

Soportes en formato `SH-aLRT/UFBoot`: ≥ 80/95 se consideran bien soportados. El clado de plantas verdes (4 taxones) tiene soporte bajo (UFBoot 45, SH-aLRT 28.1); el de plantas verdes + ectocarpus también (UFBoot 56, SH-aLRT 35.8). `bartmosca+humano` tiene UFBoot 75 (soporte moderado, no llega al umbral 95) y `bartmosca+humano+mosca` UFBoot 86.

**Longitudes de rama terminales** (sustituciones/sitio, árbol enraizado):

| Terminal | Longitud de rama |
|---|---|
| rhizobium | 0.436676 |
| ectocarpus | 0.301348 |
| mosca | 0.225732 |
| cebada | 0.101677 |
| bartmosca | 0.049700 |
| humano | 0.044395 |
| abedul | 0.041946 |
| altramuz | 0.022298 |
| vigna | 0.006840 |

>*Nota:* la rama de `rhizobium` aparece con la mitad de su longitud (0.4367) porque la raíz se colgó en el punto medio de esa rama; en el `.treefile` sin enraizar mide 0.8734.

### 6.5 Bootstrap estándar (100 réplicas, modelo fijo)

El `-b 1000` con bootstrap estándar resultó demasiado lento en esta máquina (≈ 4 s por réplica, más de 1 h; se abortó) y la corrida con `-B 100` **no es posible**: IQ-TREE 2.3.6 exige al menos 1000 réplicas para UFBoot y aborta con `#replicates must be >= 1000`. Por eso el bootstrap estándar se corrió con **100 réplicas y modelo fijo `Dayhoff+G4`**:

```bash
tools/iqtree-2.3.6-Linux-intel/bin/iqtree2 \
  -s resultados/aln_clustalo_set9.fasta \
  -m Dayhoff+G4 -b 100 -T 2 -o rhizobium \
  --prefix resultados/iqtree_set9_boot100 -seed 12345 -redo
```

La corrida terminó con log-likelihood -1148.9603 y BIC 2375.9238; los soportes están en la columna «Bootstrap estándar» de §6.4 (archivos `resultados/iqtree_set9_boot100.treefile` y `.iqtree`).

### 6.6 Efecto del número de réplicas: UFBoot 100 vs UFBoot 1000

Como `-B 100` está prohibido por IQ-TREE, la comparación se hizo sobre las **mismas 1000 réplicas** de UFBoot (guardadas con `--boot-trees` en `resultados/iqtree_set9_B1000.ufboot`, modelo fijo `Dayhoff+G4`, semilla 12345): se calculó la frecuencia de cada split con las **primeras 100 réplicas** y con las **1000**, que es exactamente el efecto de reducir el número de réplicas (subconjunto anidado). La corrida completa se validó contra los soportes del `.treefile` (coinciden dentro del redondeo).

| Nodo | UFBoot 100 (%) | UFBoot 1000 (%) | Δ (pp) | Bootstrap estándar 100 (%) |
|---|---|---|---|---|
| humano+bartmosca | 78.0 | 74.5 | +3.5 | 59.0 |
| clado animal (humano+mosca+bartmosca) | 84.0 | 85.1 | -1.1 | 81.0 |
| clado de plantas verdes (cebada+altramuz+abedul+vigna) | 36.0 | 44.4 | -8.4 | 29.0 |
| plantas verdes + ectocarpus | 53.0 | 56.3 | -3.3 | 62.0 |
| todos menos rhizobium (raiz del ingroup) | 100.0 | 100.0 | +0.0 | — |

*(Los valores de UFBoot 1000 de esta tabla provienen de la corrida con modelo fijo; coinciden dentro de ~1 pp con los del árbol ML de §6.4.)*

Lectura: con 100 réplicas el soporte se mueve hasta **~8 puntos porcentuales** respecto de 1000 (el caso más claro es el clado de plantas verdes), mientras que los nodos con soporte alto/intermedio apenas cambian (≤ 3.5 pp). Con 100 réplicas la incertidumbre del propio soporte es del orden de ±5–8 pp; por eso se reportan los valores de 1000 réplicas como los definitivos.

## 7. Filogenia del set3 (humano, mosca, bartmosca)

```bash
tools/iqtree-2.3.6-Linux-intel/bin/iqtree2 \
  -s resultados/aln_clustalo_set3.fasta \
  -m Dayhoff+G4 -T 2 \
  --prefix resultados/iqtree_set3 -seed 12345 -redo
```

**Modelo:** fijado a `Dayhoff+G4` (el mejor por BIC en el set9 y el criterio pedido para comparar con el árbol de 9). Como referencia, la corrida previa con ModelFinder (`-m MFP`, en `resultados/iqtree_set3_mfp.*`) había elegido `DCMut` (BIC = 967.3606, logL = -476.6433); con el modelo fijo `Dayhoff+G4` el árbol es el mismo (la topología no puede cambiar, ver abajo).

- **Modelo de sustitución:** `Dayhoff+G4` (alpha gamma = 3.0901)
- **log-likelihood:** `-476.6251`
- **BIC:** 972.0155
- **Longitud total del árbol:** 0.3650 sustituciones/sitio
- **Bootstrap:** no es posible añadir réplicas de bootstrap con 3 secuencias: `-b 1000` y también `-b 100` fallan con `ERROR: It makes no sense to perform bootstrap with less than 4 sequences` (IQ-TREE exige ≥ 4 taxones). No se puede reportar soporte de nodos porque **no hay nodos internos**.
- **Enraizamiento:** con 3 taxones el árbol no enraizado es una **única trifurcación**: no hay ningún nodo interno que resolver y no existe un outgroup dentro del set (los tres son metazoos y el TP los trata como ingroup). Por eso el resultado se muestra **sin enraizar** y no se generó ninguna versión enraizada: el enraizamiento es arbitrario. La única señal está en las longitudes de rama.

Newick (sin enraizar):

```
(humano:0.0488689776,mosca:0.2752227133,bartmosca:0.0409441262);
```

| Terminal | Longitud de rama (sust./sitio) |
|---|---|
| humano | 0.048869 |
| mosca | 0.275223 |
| bartmosca | 0.040944 |

bartmosca (0.0409) y humano (0.0489) tienen ramas cortas y mosca es la más divergente (0.2752), lo que es coherente con el clado `(bartmosca, humano)` del set9. Con 3 taxones el bootstrap no agrega topologías posibles, así que la interpretación se apoya en las longitudes de rama y en las identidades de §4.4.

>*Detalle:* el log de IQ-TREE informa `Number of parsimony informative sites: 0` para este set. Es correcto con el criterio de Fitch que usa IQ-TREE (al menos dos estados y **cada uno** presente en ≥ 2 secuencias), imposible de cumplir con solo 3 taxones. Con el criterio útil para 3 taxones (dos secuencias que comparten un estado y la tercera distinta) hay **25 sitios** de 109. En cualquier caso el árbol no enraizado de 3 taxones es único y no hay nodos internos que soportar.

## 8. Figuras

Todas están en `resultados/` (no hay carpeta de figuras aparte):

| Figura | Contenido |
|---|---|
| `resultados/fig_arbol_set9.png` | Árbol rectangular del set9 enraizado por rhizobium; sobre cada nodo interno: SH-aLRT/UFBoot del árbol ML y bootstrap estándar de 100 réplicas (bs); barra de escala de 0.1 sust./sitio; longitudes de rama en sustituciones por sitio. |
| `resultados/fig_arbol_set3.png` | Panel A: árbol no enraizado del set3 (trifurcación única) con barra de escala. Panel B: barras con las longitudes de rama terminales. |
| `resultados/fig_ventanas_bartmosca.png` | Identidad de bartmosca contra humano y mosca en ventanas de 15 aa (paso 5), con la región N-terminal sombreada. |
| `resultados/fig_heatmap_identidad.png` | Heatmap 9×9 de identidad de pares del set9 (deleción de pares, Clustal Omega), con los valores numéricos. |

## 9. Cómo abrir y anotar los árboles en FigTree

1. Abrir FigTree y hacer **File → Open…**, elegir `resultados/iqtree_set9.treefile` (o `resultados/tree_set9_rooted.nwk`, que ya está enraizado por el outgroup).
2. **Enraizar:** si se abre el `.treefile` sin enraizar, seleccionar la rama de `rhizobium` y usar **Tree → Reroot on Selected Branch** (equivale al `-o rhizobium`). Alternativa: **Tree → Midpoint Root** (enraizamiento por punto medio, no usa outgroup).
3. **Mostrar soportes:** panel izquierdo → *Node Labels* → activar, y en el desplegable elegir el valor que IQ-TREE guardó como etiqueta de nodo (`SH-aLRT/UFBoot`; en la corrida `-b 100` la etiqueta es el bootstrap estándar). Ajustar el tamaño de fuente si se solapan.
4. **Longitudes de rama:** activar *Branch Labels → Branch Lengths* para verlas, o dejar la escala (*Scale Axis*) visible abajo.
5. **Exportar:** **File → Export PDF…** (vectorial, ideal para el informe) o **Export PNG…**; opcionalmente **File → Save Tree…** para guardar el árbol enraizado/anotado en Newick.
6. Para el set3 abrir `resultados/tree_set3.nwk`: al ser una trifurcación no tiene sentido re-enraizar (FigTree lo permite pero el resultado es arbitrario).

## 10. Comandos exactos ejecutados

Todos se corren desde la raíz del repositorio.

```bash
# 1) Conjuntos con etiquetas cortas
python3 src/01_prepare_sets.py

# 2) Alineamientos
mafft --localpair --maxiterate 1000 \
    datos/set9.fasta > resultados/aln_mafft_set9.fasta
mafft --localpair --maxiterate 1000 \
    datos/set3.fasta > resultados/aln_mafft_set3.fasta
bash src/02_ebi_clustalo.sh        # Clustal Omega (REST EBI, outfmt=fa)

# 3) Matrices de identidad / p-distancia y chequeo colapsando gaps
python3 src/03_pairwise.py

# 4) Mosaico de bartmosca + heatmap
MPLCONFIGDIR="$PWD/.mplcache" python3 src/04_mosaico.py
MPLCONFIGDIR="$PWD/.mplcache" python3 src/06_heatmap.py

# 5) IQ-TREE
IQ=tools/iqtree-2.3.6-Linux-intel/bin/iqtree2
$IQ -s resultados/aln_clustalo_set9.fasta -m MFP -B 1000 -alrt 1000 \
    -T 2 -o rhizobium -seed 12345 --prefix resultados/iqtree_set9 -redo
# bootstrap estandar: 100 replicas con modelo fijo (el -b 1000 tardaba >1 h)
$IQ -s resultados/aln_clustalo_set9.fasta -m Dayhoff+G4 -b 100 -T 2 \
    -o rhizobium --prefix resultados/iqtree_set9_boot100 -seed 12345 -redo
# UFBoot 1000 con modelo fijo, guardando las replicas para comparar 100 vs 1000
$IQ -s resultados/aln_clustalo_set9.fasta -m Dayhoff+G4 -B 1000 -T 2 \
    -o rhizobium --boot-trees --prefix resultados/iqtree_set9_B1000 \
    -seed 12345 -redo
# OJO: -B 100 no existe en IQ-TREE 2.3.6 (exige >=1000 replicas)
# set3 con modelo fijo
$IQ -s resultados/aln_clustalo_set3.fasta -m Dayhoff+G4 -T 2 -seed 12345 \
    --prefix resultados/iqtree_set3 -redo

# comparacion UFBoot 100 vs 1000 (submuestreo de las mismas replicas)
python3 src/10_boot_comparacion.py

# 6) Árboles enraizados/sin enraizar y figuras
MPLCONFIGDIR="$PWD/.mplcache" python3 src/05_trees.py

# 7) Métricas y resumen
python3 src/07_metricas.py
python3 src/08_resumen.py
```

El detalle del envío/polling a EBI está en `src/02_ebi_clustalo.sh`; los IDs de trabajo quedaron registrados en `resultados/ebi_jobids.txt`:

| Set | Job ID de EBI |
|---|---|
| set3 | clustalo-R20261028-001614-0565-18522630-p1m |
| set9 | clustalo-R20261028-001625-0783-85075839-p1m |

## 11. Archivos generados

| Archivo | Qué es |
|---|---|
| `datos/set3.fasta` | 3 secuencias con etiquetas cortas |
| `datos/set9.fasta` | 9 secuencias con etiquetas cortas |
| `resultados/aln_clustalo_set{3,9}.fasta` | alineamientos de Clustal Omega (EBI) |
| `resultados/aln_mafft_set{3,9}.fasta` | alineamientos de MAFFT |
| `resultados/pairwise_{clustalo,mafft}_set{3,9}.csv` | tabla de pares: identidad, p-distancia, identidad colapsada, gaps como diferencias |
| `resultados/pairwise_identidad_*.csv`, `pairwise_pdist_*.csv`, `pairwise_identidad_colapsada_*.csv` | matrices cuadradas (CSV) |
| `resultados/bartmosca_ventanas.csv` | 24 ventanas de 15 aa (paso 5) |
| `resultados/bartmosca_regiones.csv` | identidad N-terminal vs resto |
| `resultados/iqtree_set9.*` | corrida ML + ModelFinder + UFBoot1000/SH-aLRT1000 (set9, árbol principal) |
| `resultados/iqtree_set9_boot100.*` | bootstrap estándar 100, modelo fijo Dayhoff+G4 (set9) |
| `resultados/iqtree_set9_B1000.*` | UFBoot 1000 con modelo fijo + `.ufboot` con las 1000 réplicas |
| `resultados/boot_comparacion.csv` | UFBoot 100 vs UFBoot 1000 por nodo |
| `resultados/iqtree_set3.*` | corrida ML del set3 con modelo fijo Dayhoff+G4 |
| `resultados/iqtree_set3_mfp.*` | referencia: corrida previa del set3 con ModelFinder (eligió DCMut) |
| `resultados/tree_set9_rooted.nwk` | árbol set9 enraizado por outgroup (Newick) |
| `resultados/tree_set9_unrooted.nwk` | árbol set9 sin enraizar (Newick) |
| `resultados/tree_set3.nwk` | árbol set3 sin enraizar (Newick) |
| `resultados/soportes_nodos.csv` | tabla nodo → SH-aLRT / UFBoot / bootstrap / rama |
| `resultados/fig_arbol_set9.png`, `fig_arbol_set3.png`, `fig_ventanas_bartmosca.png`, `fig_heatmap_identidad.png` | figuras |
| `resultados/metricas.json` | todos los números clave en JSON |
| `resultados/resumen_analisis.md` | este documento |

## 12. Limitaciones y notas de honestidad

1. **Clustal Omega vs MAFFT:** para set9 coinciden en 131 columnas pero no en el detalle de los gaps; por eso algunas identidades de pares cambian hasta ~2 puntos porcentuales (p. ej. bartmosca–mosca: 75.47 % con ClustalO vs 77.36 % con MAFFT). El chequeo «colapsando gaps» da prácticamente lo mismo con ambos alineadores (columna correspondiente en §4.1), lo que indica que la diferencia es de colocación de gaps y no de homología.
2. **Bootstrap:** con set3 no se puede correr ningún bootstrap (`-b 1000` ni `-b 100`): IQ-TREE exige ≥ 4 secuencias (`ERROR: It makes no sense to perform bootstrap with less than 4 sequences`). En set9, el bootstrap estándar `-b 1000` se abortó por lentitud (≈ 4 s/réplica, > 1 h) y se reemplazó por `-b 100` con modelo fijo. `-B 100` no existe en IQ-TREE 2.3.6 (mínimo 1000), por lo que la comparación 100 vs 1000 se hizo submuestreando las 1000 réplicas de UFBoot guardadas con `--boot-trees` (subconjunto anidado, mismo seed).
3. **Enraizamiento set9:** se usó el outgroup `rhizobium`. La rama del outgroup se partió por la mitad para colgar la raíz; en el `.treefile` original esa rama conserva su longitud completa (0.8734).
4. **Set3:** no se puede enraizar por outgroup con los tres metazoos; el árbol no enraizado es único (trifurcación). No se inventó un enraizamiento.
5. **Soportes:** los nodos con UFBoot ≥ 95 y SH-aLRT ≥ 80 se consideran bien apoyados. El clado de plantas verdes tiene soporte bajo (UFBoot 45, SH-aLRT 28.1), y el de plantas verdes + ectocarpus también (UFBoot 56, SH-aLRT 35.8). `bartmosca+humano` (UFBoot 75, SH-aLRT 77.2) queda por debajo de los umbrales: hay que tomarlo como soporte moderado, no como un clado robusto.
6. La identidad «colapsando gaps» se implementó re-alineando cada par por separado (BLOSUM62) en lugar de comparar las cadenas sin gaps posición a posición, porque los indels desplazan el marco y ese cálculo directo daría valores sin sentido (p. ej. ~9 %).

## 13. Reproducibilidad

Orden de ejecución: `src/01_prepare_sets.py` → `src/02_ebi_clustalo.sh` → `src/03_pairwise.py` → `src/04_mosaico.py` → `src/05_trees.py` → `src/06_heatmap.py` → `src/10_boot_comparacion.py` → `src/07_metricas.py` → `src/08_resumen.py`. Las figuras requieren `MPLCONFIGDIR="$PWD/.mplcache"`. Los números de este resumen se generan desde `resultados/metricas.json` y los CSV, de modo que si se vuelve a correr el pipeline el documento se regenera consistente.
