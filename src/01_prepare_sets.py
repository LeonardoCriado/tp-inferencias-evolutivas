#!/usr/bin/env python3
"""Prepara datos/set3.fasta y datos/set9.fasta con etiquetas cortas, sin espacios.

Entrada: datos/SecuenciasCytocromoC.fasta (9 secuencias de citocromo c).
Salida:
  datos/set3.fasta   -> bartmosca, humano, mosca
  datos/set9.fasta   -> las 9 secuencias (bartmosca + 2 metazoos + 6 referencia)
  resultados/label_map.tsv -> mapeo etiqueta_corta <-> header original <-> accession
"""
import os
import re
from Bio import SeqIO

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "datos")
RESULTS = os.path.join(BASE, "resultados")

# etiqueta corta -> (accession en el header original, subcadena identificadora)
LABELS = [
    ("humano",      "NP_061820.1",  "Homo_sapiens"),
    ("mosca",       "NP_001285984.1", "Drosophila_melanogaster"),
    ("bartmosca",   None,            "bartmosca"),
    ("cebada",      "KAE8767936.1",  "Hordeum_vulgare"),
    ("altramuz",    "KAE9612850.1",  "Lupinus_albus"),
    ("abedul",      "QEE59979.1",    "Betula_platyphylla"),
    ("vigna",       "NP_001363126.1", "Vigna_radiata"),
    ("ectocarpus",  "CBJ31344.1",    "Ectocarpus_siliculosus"),
    ("rhizobium",   "WP_037127040.1", "Rhizobium_sp._CF097"),
]
SET3 = {"bartmosca", "humano", "mosca"}


def main():
    records = list(SeqIO.parse(os.path.join(DATA, "SecuenciasCytocromoC.fasta"), "fasta"))
    print(f"Secuencias leidas: {len(records)}")

    assigned = {}
    used = set()
    for label, acc, tag in LABELS:
        hits = [r for r in records if (acc and acc in r.id) or (tag in r.description) or (tag in r.id)]
        hits = [h for h in hits if h.id not in used]
        if len(hits) != 1:
            raise SystemExit(f"ERROR: no pude asignar '{label}' de forma unica (candidatos={[h.id for h in hits]})")
        used.add(hits[0].id)
        assigned[label] = hits[0]
    if len(used) != len(records):
        raise SystemExit(f"ERROR: {len(records) - len(used)} secuencias quedaron sin etiqueta")

    # orden de set9: metazoos primero, luego bartmosca y referencias
    order9 = ["humano", "mosca", "bartmosca", "cebada", "altramuz", "abedul", "vigna", "ectocarpus", "rhizobium"]
    order3 = ["humano", "mosca", "bartmosca"]

    def dump(path, order):
        with open(path, "w") as fh:
            for lab in order:
                rec = assigned[lab]
                fh.write(f">{lab}\n{str(rec.seq)}\n")
        print(f"escrito {path} ({len(order)} seqs)")

    dump(os.path.join(DATA, "set9.fasta"), order9)
    dump(os.path.join(DATA, "set3.fasta"), order3)

    with open(os.path.join(RESULTS, "label_map.tsv"), "w") as fh:
        fh.write("etiqueta\theader_original\taccession\tlargo\n")
        for lab in order9:
            r = assigned[lab]
            fh.write(f"{lab}\t{r.description}\t{r.id}\t{len(r.seq)}\n")

    # chequeo: el nombre "bartmosca" no debe tener espacios ni caracteres raros
    for lab in order9:
        assert re.fullmatch(r"[a-z]+", lab), lab
    lens = {lab: len(assigned[lab].seq) for lab in order9}
    print("Largos crudos:", lens)


if __name__ == "__main__":
    main()
