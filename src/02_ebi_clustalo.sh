#!/usr/bin/env bash
# Envia alineamientos a Clustal Omega (REST de EBI) y espera el resultado.
# Uso: ./02_ebi_clustalo.sh
set -u
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE"
# El servicio pide un mail de contacto; se puede pisar con la variable EBI_EMAIL.
EMAIL="${EBI_EMAIL:-usuario@example.com}"
REST="https://www.ebi.ac.uk/Tools/services/rest/clustalo"
mkdir -p resultados datos
: > resultados/ebi_jobids.txt

submit () {
  local set="$1" infile="$2" outfile="$3"
  echo "[$set] enviando $infile ..." >&2
  local jobid
  jobid=$(curl -sS --max-time 120 -X POST "$REST/run" \
      --data-urlencode "email=$EMAIL" \
      --data-urlencode "stype=protein" \
      --data-urlencode "outfmt=fa" \
      --data-urlencode "sequence@$infile" 2>>resultados/clustalo_curl.err)
  echo "[$set] jobid=${jobid}" >&2
  if ! [[ "$jobid" =~ ^[a-zA-Z0-9_-]+$ ]]; then
     echo "ERROR: jobid invalido para $set: $jobid" >&2
     echo -e "$set\tFAILED\t-\t$jobid" >> resultados/ebi_jobids.txt
     return 1
  fi
  echo -e "$set\t$jobid\t$infile\t$outfile" >> resultados/ebi_jobids.txt

  local t=0 status=""
  while [ $t -lt 360 ]; do
    status=$(curl -sS --max-time 60 "$REST/status/$jobid" 2>>resultados/clustalo_curl.err)
    echo "[$set] t=${t}s status=$status" >&2
    case "$status" in
      FINISHED) break ;;
      ERROR|FAILURE|NOT_FOUND|CANCELLED)
          echo "ERROR: job $jobid de $set termino en $status" >&2
          curl -sS "$REST/status/$jobid" >&2
          return 2 ;;
    esac
    sleep 5
    t=$((t+5))
  done
  if [ "$status" != "FINISHED" ]; then
    echo "ERROR: timeout esperando $jobid ($set)" >&2
    return 3
  fi
  curl -sS --max-time 120 "$REST/result/$jobid/aln-fasta" -o "$outfile"
  echo "[$set] escrito $outfile" >&2
  head -2 "$outfile" >&2
}

submit set3 datos/set3.fasta resultados/aln_clustalo_set3.fasta
submit set9 datos/set9.fasta resultados/aln_clustalo_set9.fasta
echo "DONE"
