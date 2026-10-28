#!/usr/bin/env python3
"""Matrices de identidad y p-distancia para los alineamientos (Clustal Omega y MAFFT).

Definiciones (documentadas en el resumen):
  ident_pairdel  : matches / columnas donde AMBAS secuencias tienen residuo
                   (delecion de pares; identidad estandar del TP).
  ident_colapsada: chequeo "colapsando gaps": se eliminan los gaps de cada
                   secuencia y el par se re-alinea de forma independiente
                   (PairwiseAligner global, BLOSUM62, gap open -10, extend -0.5);
                   la identidad se recalcula con delecion de pares sobre ese
                   alineamiento par a par. Si difiere de ident_pairdel, el valor
                   depende del alineador multiple elegido.
  ident_full     : matches / largo total del alineamiento (los gaps cuentan como
                   diferencia). Cota inferior, util para comparar con la literatura.
  p_dist         : 1 - ident_pairdel (p-distancia no corregida, delecion de pares).
"""
import os
import sys
import itertools
from Bio import AlignIO
from Bio.Align import PairwiseAligner, substitution_matrices
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "resultados")

_ALIGNER = None


def aligner():
    global _ALIGNER
    if _ALIGNER is None:
        a = PairwiseAligner()
        a.mode = "global"
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score = -10.0
        a.extend_gap_score = -0.5
        a.end_gap_score = 0.0
        _ALIGNER = a
    return _ALIGNER


def collapsed_identity(sa, sb):
    """Identidad tras quitar gaps y re-alinear el par de forma independiente."""
    ua = sa.replace("-", "")
    ub = sb.replace("-", "")
    aln = aligner().align(ua, ub)[0]
    a, b = str(aln[0]), str(aln[1])
    both = match = 0
    for x, y in zip(a, b):
        if x != "-" and y != "-":
            both += 1
            if x == y:
                match += 1
    return match / both if both else float("nan")


def pair_metrics(sa, sb):
    n = len(sa)
    assert len(sb) == n
    both = match = 0
    for x, y in zip(sa, sb):
        if x != "-" and y != "-":
            both += 1
            if x == y:
                match += 1
    ident_pairdel = match / both if both else float("nan")
    ident_full = match / n if n else float("nan")
    return dict(n_columns=n, n_posiciones_ambas=both, matches_pairdel=match,
                ident_pairdel=ident_pairdel,
                ident_colapsada=collapsed_identity(sa, sb),
                ident_full=ident_full, p_dist=1 - ident_pairdel,
                gaps_a=sa.count("-"), gaps_b=sb.count("-"))


def write_set(path, tag):
    aln = AlignIO.read(path, "fasta")
    seqs = {r.id: str(r.seq).upper() for r in aln}
    ids = [r.id for r in aln]
    L = aln.get_alignment_length()
    rows = []
    mat_id, mat_pd, mat_col = {}, {}, {}
    for a, b in itertools.combinations(ids, 2):
        d = pair_metrics(seqs[a], seqs[b])
        rows.append(dict(par=f"{a}-{b}", **d))
        for m, key in ((mat_id, "ident_pairdel"), (mat_pd, "p_dist"), (mat_col, "ident_colapsada")):
            m.setdefault(a, {})[b] = round(d[key], 4)
            m.setdefault(b, {})[a] = round(d[key], 4)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, f"pairwise_{tag}.csv"), index=False, float_format="%.4f")

    for name, m in (("identidad", mat_id), ("pdist", mat_pd), ("identidad_colapsada", mat_col)):
        md = pd.DataFrame(m).reindex(index=ids, columns=ids)
        md.index.name = "seq"
        md.to_csv(os.path.join(RES, f"pairwise_{name}_{tag}.csv"), float_format="%.4f")

    gapfree = sum(1 for j in range(L) if all(seqs[i][j] != "-" for i in ids))
    print(f"[{tag}] n={len(ids)} aln_len={L} columnas_sin_gap={gapfree} "
          f"columnas_todo_gap={sum(1 for j in range(L) if all(seqs[i][j]=='-' for i in ids))} pares={len(rows)}")
    return df, L, ids


if __name__ == "__main__":
    tags = []
    for tag in ["clustalo_set3", "mafft_set3", "clustalo_set9", "mafft_set9"]:
        p = os.path.join(RES, f"aln_{tag}.fasta")
        if os.path.exists(p):
            write_set(p, tag)
            tags.append(tag)
        else:
            print(f"FALTA {p}", file=sys.stderr)
    print("listo:", tags)
