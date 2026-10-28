# Desafío II — El alineamiento

En este desafío alineé las tres secuencias del punto I y el conjunto completo de 9 secuencias, usado luego para la filogenia del [Desafío III](informe_desafio_III.md). También comparé dos programas para evaluar si las conclusiones dependían del alineador.

## Qué corrí

- **Clustal Omega** a través del servicio REST de EBI, que es el mismo motor que está detrás del servidor web que menciona el TP (<https://www.ebi.ac.uk/Tools/msa/clustalo/>). Parámetros del servidor para proteínas, salida FASTA. Los trabajos quedaron registrados con su ID en [`resultados/ebi_jobids.txt`](resultados/ebi_jobids.txt):

  | Set | Job ID de EBI |
  |---|---|
  | 3 secuencias (bartmosca, humano, mosca) | `clustalo-R20261028-001614-0565-18522630-p1m` |
  | 9 secuencias | `clustalo-R20261028-001625-0783-85075839-p1m` |

- **MAFFT v7.526** (`--localpair --maxiterate 1000`, L-INS-i) como control independiente, para verificar que las conclusiones no dependan del alineador.

Los insumos y los alineamientos están en [`datos/`](datos/) y [`resultados/`](resultados/): `aln_clustalo_set3.fasta`, `aln_mafft_set3.fasta`, `aln_clustalo_set9.fasta`, `aln_mafft_set9.fasta`. Las etiquetas cortas (humano, mosca, bartmosca, …) están mapeadas a los encabezados originales en [`resultados/label_map.tsv`](resultados/label_map.tsv).

## Resultado para las tres secuencias del punto I

El alineamiento de las tres secuencias tiene 109 columnas y 7 gaps, ubicados en los extremos y en la región N-terminal. Las columnas 19 a 23 contienen el motivo de unión al grupo hemo **Cys-X-X-Cys-His** (`CAQCH` en bartmosca y mosca, `CSQCH` en humano). Corresponde a Cys15, Cys18 e His19 de la secuencia completa (Cys14, Cys17 e His18 en la proteína madura). La metionina axial (Met81 completa / Met80 madura) también queda alineada en la misma columna en las tres secuencias. Estos motivos sirven como referencia para revisar la alineación del núcleo.

```
humano    ----MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE
mosca     MGVPAGDVEKGKKLFVQRCAQCHTVEAGGKHKVGPNLHGLIGRKTGQAAGFAYTDANKAKGITWNEDTLFEYLENPKKYIPGTKMIFAGLKKPNERGDLIAYLKSATK-
bartmosca --MGSGDAENGKKIFVQKCAQCHTYEVGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE
```

En el alineamiento se observa que:

- El **núcleo conservado** (desde `GGKHKTG` hasta el final) es prácticamente un calco entre las tres secuencias; el único gap de ese tramo es el terminal de la mosca.
- Las diferencias de bartmosca con el humano se concentran **antes de la columna 30**: el alineador pone 2 residuos extra en bartmosca y resuelve el resto del N-terminal con sustituciones.
- La mosca tiene un gap terminal al final (`…ATK–` frente a `…ATNE`), sin consecuencias para el resto del alineamiento.

## Comparación con MAFFT

| Set | Alineador | Columnas | Gaps | Ident. bartmosca–humano | Ident. bartmosca–mosca |
|---|---|---|---|---|---|
| 3 | Clustal Omega | 109 | 7 | 92,4% | 75,5% |
| 3 | MAFFT (L-INS-i) | 109 | 7 | 92,4% | 77,4% |
| 9 | Clustal Omega | 131 | 178 | 92,4% | 75,5% |
| 9 | MAFFT (L-INS-i) | 131 | 178 | 93,3% | 77,4% |

Los dos alineadores coinciden en la longitud del alineamiento y en el núcleo; difieren en **dónde ponen los gaps** del N-terminal. Eso mueve algunos porcentajes de identidad entre 0 y 2 puntos (bartmosca–mosca: 75,5% con Clustal Omega vs 77,4% con MAFFT) porque cambian qué posiciones se comparan. Cuando se neutraliza el efecto del alineador (se quitan los gaps y se realinea cada par por separado), los valores convergen; ese chequeo está en [`resultados/pairwise_*_set3.csv`](resultados/) y en el resumen técnico.

Para lo que importa en este TP, la conclusión es estable: bartmosca es ~92% idéntica al humano y ~76% a la mosca con cualquiera de los dos alineadores.

## Región N-terminal

El tramo inicial es el único lugar donde el alineamiento es ambiguo, y no por culpa del software: las secuencias son `MGSGDAENG…` (bartmosca) y `MGDVEKG…` (humano), sin una correspondencia única. Clustal Omega resuelve los dos residuos extra al principio; MAFFT los corre dos columnas:

```
Clustal Omega   bartmosca  --MGSGDAENGKKIFVQKCAQCHT
                humano     ----MGDVEKGKKIFIMKCSQCHT
MAFFT           bartmosca  MG--SGDAENGKKIFVQKCAQCHT
                humano     ----MGDVEKGKKIFIMKCSQCHT
```

Cambia dónde se ubican los gaps, no la composición de residuos ni la conclusión. Pero es un recordatorio concreto de lo que dice el TP sobre que solo el alineamiento correcto produce inferencias correctas: en esta región no hay una respuesta única, hay que decidirla a mano (mirando la estructura, el motivo CXXCH y, si se puede, secuencias más cercanas). Lo que no conviene es aceptar el primer alineamiento automático y sacar conclusiones al residuo exacto en esa zona.

## Revisión del alineamiento

Un alineamiento múltiple establece correspondencias entre posiciones de las secuencias. Si esas correspondencias son incorrectas, pueden cambiar las identidades, las distancias y la topología o el soporte del árbol. Por eso revisé los motivos conservados y la ubicación de los gaps antes de usar el alineamiento en los análisis posteriores.

En este caso las comprobaciones que hice fueron:

1. **Motivos funcionales alineados**: CXXCH (Cys14–Cys17–His18) y Met80 en la misma columna en las tres secuencias.
2. **Gaps fuera del núcleo**: los 7 gaps del set de tres están en los extremos o en el N-terminal; ninguno cae dentro del bloque `GGKHKTG…ATNE`.
3. **Dos alineadores independientes** dan la misma longitud y el mismo núcleo; las diferencias quedan confinadas a la región ambigua.
4. **Verificación contra la estructura**: el citocromo c es una proteína globular chica y su alineamiento no debería romper las hélices; el núcleo libre de gaps es consistente con eso.

Si el alineamiento fuera de un set difícil (proteínas con dominios repetidos, indels largos), el orden sería otro: alineamiento estructural o con información de estructura, ajuste de penalidades de gaps, y curaduría manual de las columnas dudosas. Acá no hizo falta llegar a eso.

## Con el set completo (9 secuencias)

El alineamiento de las 9 secuencias tiene 131 columnas y 178 gaps, contra 109 y 7 del set de tres. La diferencia no es casual: al sumar plantas, un alga y una bacteria aparecen indels que en el set de tres no existían (el rhizobium trae una extensión N-terminal de ~23 residuos). El alineamiento sigue siendo razonable en el núcleo, pero es el que hay que mirar con más cuidado si se van a interpretar nodos profundos. IQ-TREE resume el contenido informativo de este alineamiento: **44 sitios parsimonia-informativos** de 131 columnas (55 constantes). Es poco; de ahí que varios nodos del árbol queden con soporte bajo.

Para la filogenia no recorté columnas (no usé Gblocks/trimAl). IQ-TREE maneja los gaps como datos faltantes y el alineamiento completo funcionó; recortar con un método automático es una alternativa razonable, pero cambia el conjunto de sitios y no era necesario para este caso.

## Referencias

- Ashkenazy, H., Sela, I., Levy Karin, E., Landan, G. y Pupko, T. 2019. Multiple sequence alignment averaging improves phylogeny reconstruction. *Systematic Biology* 68(1): 117–130.
- Iantorno, S., Gori, K., Goldman, N., Gil, M. y Dessimoz, C. 2014. Who watches the watchmen? An appraisal of benchmarks for multiple sequence alignment. *Methods in Molecular Biology* 1079: 59–73.
- Katoh, K. y Standley, D.M. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30(4): 772–780.
- Larkin, M.A. *et al.* 2007. Clustal W and Clustal X version 2.0. *Bioinformatics* 23(21): 2947–2948.
- Sievers, F. *et al.* 2011. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. *Molecular Systems Biology* 7: 539.
- Thompson, J.D., Higgins, D.G. y Gibson, T.J. 1994. CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice. *Nucleic Acids Research* 22(22): 4673–4680.
