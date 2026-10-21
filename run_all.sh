#!/bin/bash
# End-to-end pipeline. Resumable: every step skips work whose outputs already exist.
#
#   GPUS="0 1 2 3" bash run_all.sh            # all steps
#   GPUS="0 1" STEPS="extract eval" bash run_all.sh
#
# Env: GPUS (CUDA ids, default "0"), ENCODERS, STEPS (download index gate extract eval report),
#      PY (python, default .venv/bin/python), BS (extraction batch, default 256), DECODE (decode threads per
#      extraction process, default 12), ABMIL_JOBS (parallel ABMIL fits, default 4 per GPU), LR_JOBS (default 16).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-.venv/bin/python}
GPUS=${GPUS:-0}
read -r -a GPU_ARR <<< "$GPUS"
NG=${#GPU_ARR[@]}
# slowest first so the GPU queues balance
ENCODERS=${ENCODERS:-"gigapath uni2_h virchow2 phikon_v2 uni phikon"}
STEPS=${STEPS:-"download index gate extract features eval report"}
TIA=${TIA:-.venv-tia/bin/python}   # separate env with tiatoolbox for cells.py
FEATS="feat_tissue feat_cells feat_interp"   # interpretable features: LR + confound probes only
BS=${BS:-256}
DECODE=${DECODE:-12}
ABMIL_JOBS=${ABMIL_JOBS:-$((4 * NG))}
LR_JOBS=${LR_JOBS:-16}
mkdir -p logs results/metrics results/preds
has() { [[ " $STEPS " == *" $1 "* ]]; }

if has download; then
  $PY src/download.py all 2>&1 | tee logs/download.log
fi

if has index; then
  [ -f data/index.parquet ] || $PY src/build_index.py 2>&1 | tee logs/index.log
fi

if has gate; then
  # fp16 autocast must match fp32 (min cosine > 0.999 on 256 patches) on this hardware before we use it
  for e in $ENCODERS; do
    [ -f emb/$e/fp16_check.json ] || CUDA_VISIBLE_DEVICES=${GPU_ARR[0]} $PY src/extract.py --encoder $e --check-fp16
    grep -q '"pass": true' emb/$e/fp16_check.json || { echo "fp16 check failed for $e"; exit 1; }
  done
fi

if has extract; then
  # static round-robin of encoders over GPUs; each GPU runs its queue sequentially
  i=0
  declare -A QUEUE
  for e in $ENCODERS; do g=${GPU_ARR[$((i % NG))]}; QUEUE[$g]="${QUEUE[$g]:-} $e"; i=$((i + 1)); done
  for g in "${GPU_ARR[@]}"; do
    (
      for e in ${QUEUE[$g]:-}; do
        [ -f emb/$e/meta.json ] && continue
        CUDA_VISIBLE_DEVICES=$g $PY src/extract.py --encoder $e --fp16 --bs $BS --workers $DECODE \
          > logs/extract_$e.log 2>&1 || echo "EXTRACT FAILED $e (see logs/extract_$e.log)"
      done
    ) &
  done
  wait
fi

if has features; then
  # count-invariant interpretable features: tissue-type fractions (UNI probe on NCT-CRC-HE-100K) + nucleus densities
  $PY src/download.py nct 2>&1 | tee logs/download_nct.log   # idempotent, 12.5 GB
  [ -f data/nct/cells_val.npz ] || CUDA_VISIBLE_DEVICES=${GPU_ARR[0]} $TIA src/cells.py --set nct_val --fp16
  [ -f emb/feat_cells/meta.json ] || CUDA_VISIBLE_DEVICES=${GPU_ARR[0]} $TIA src/cells.py --set kather --fp16 > logs/cells.log 2>&1
  [ -f data/nct/uni_NCT-CRC-HE-100K.npz ] || CUDA_VISIBLE_DEVICES=${GPU_ARR[0]} $PY src/tissue.py embed-nct
  [ -f emb/uni/quads2x.f16.npy ] || CUDA_VISIBLE_DEVICES=${GPU_ARR[0]} $PY src/tissue.py embed-quads
  [ -f emb/feat_tissue/meta.json ] || $PY src/tissue.py fit 2>&1 | tee logs/tissue_fit.log
  [ -f emb/feat_interp/meta.json ] || $PY src/features.py
fi

if has eval; then
  [ -f emb/phikon_ref/meta.json ] || $PY src/build_ref_emb.py 2>&1 | tee logs/build_ref_emb.log
  [ -f results/metrics/count_stats.csv ] || $PY src/confound.py --counts > /dev/null
  for e in phikon_ref $ENCODERS; do
    [ -f emb/$e/meta.json ] || { echo "skip $e (no embeddings)"; continue; }
    $PY src/pool.py --encoder $e
    [ -f results/metrics/${e}_lr.csv ] || $PY src/eval_lr.py --encoder $e --jobs $LR_JOBS 2>&1 | tee logs/eval_lr_$e.log
    [ -f results/metrics/${e}_abmil.csv ] || CUDA_VISIBLE_DEVICES=$(echo $GPUS | tr ' ' ',') $PY src/abmil.py --encoder $e --jobs $ABMIL_JOBS 2>&1 | tee logs/abmil_$e.log
    $PY src/confound.py --encoder $e > logs/confound_$e.log
  done
  for e in $FEATS; do
    [ -f emb/$e/meta.json ] || { echo "skip $e (no features)"; continue; }
    $PY src/pool.py --encoder $e
    [ -f results/metrics/${e}_lr.csv ] || $PY src/eval_lr.py --encoder $e --jobs $LR_JOBS 2>&1 | tee logs/eval_lr_$e.log
    $PY src/confound.py --encoder $e > logs/confound_$e.log
  done
fi

if has report; then
  $PY src/report.py 2>&1 | tee logs/report.log
fi
