#!/usr/bin/env bash
set -euo pipefail
# --- User knobs (set via env) ---
PARAMS_FILE="${PARAMS_FILE:-input/params.txt}"
QUENCH_BIN="${QUENCH_BIN:-./quenching}"
QHAT_VALUES=(${QHAT_VALUES:-0.051 0.09}) # 0.075
ALPHA_MODES=(${ALPHA_MODES:-0.5 0})        # nonzero=constant, 0=running
FORCE="${FORCE:-0}"                        # 1 to overwrite
DRYRUN="${DRYRUN:-0}"
# e.g., ROOTSNN_LIST="5023 8160" or "8160" or "200 2760 5023 8160"
ROOTSNN_LIST=(${ROOTSNN_LIST:-5023 8160})  # <-- unquoted default to split into two items

read_kv() { awk -v k="$1" '
  /^[[:space:]]*\/\//{next} /^[[:space:]]*$/{next} $1==k{print $2; exit}' "$PARAMS_FILE"; }

CTYPE=$(read_kv collisionType); case "${CTYPE:-1}" in 1) COLL=pPb;; 2) COLL=AB;; 0) COLL=both;; *) COLL=unknown;; esac
PTYPE=$(read_kv particleType); case "${PTYPE:-1}" in 1) PART=JPsi;; 0) PART=Upsilon;; *) PART=Unknown;; esac

# Treat any non-zero as "constant", 0 as "running" (works for 0.5, 1, etc.)
alpha_tag(){
  awk -v a="$1" 'BEGIN{exit (a!=0)?0:1}'; [[ $? -eq 0 ]] && echo constant || echo running
}

# Safely swap rootsnn inside PARAMS_FILE
_patch_rootsnn() {
  local new="$1"
  awk -v rs="$new" '
    BEGIN{done=0}
    /^[[:space:]]*\/\//{print; next}
    /^[[:space:]]*$/{print; next}
    $1=="rootsnn" && !done {print "rootsnn       " rs; done=1; next}
    {print}
    END{if(!done) print "rootsnn       " rs}
  ' "$PARAMS_FILE" > "${PARAMS_FILE}.tmp"
  mv "${PARAMS_FILE}.tmp" "$PARAMS_FILE"
}

run_one() {
  local rootsnn="$1" a="$2" q="$3" tag; tag=$(alpha_tag "$a")
  local base="output_${rootsnn}GeV_${COLL}_alpha_${tag}"
  local dest="${base}/output_qhat0_${q}"
  echo "[INFO] √sNN=$rootsnn, a=$a ($tag), qhat0=$q → $dest"
  [[ $DRYRUN -eq 1 ]] && return 0

  mkdir -p "$dest"
  "$QUENCH_BIN" -alphas "$a" -qhat0 "$q" | tee "${dest}/quenching.log"

  [[ -d output ]] || { echo "ERROR: 'output' dir not produced"; exit 2; }
  shopt -s dotglob && mv output/* "$dest"/ && rmdir output || true

  cp -f "$PARAMS_FILE" "${dest}/params.txt"
  printf "rootsnn:%s\ncollision:%s\nparticle:%s\nalphas:%s\nqhat0:%s\n" \
    "$rootsnn" "$COLL ($CTYPE)" "$PART ($PTYPE)" "$a" "$q" > "${dest}/overrides.txt"
}

# --- MAIN ---
PARAMS_BAK="${PARAMS_FILE}.bak"
cp -f "$PARAMS_FILE" "$PARAMS_BAK"

for R in "${ROOTSNN_LIST[@]}"; do
  echo "[INFO] === Running √sNN=${R} GeV ==="
  _patch_rootsnn "$R"

  for a in "${ALPHA_MODES[@]}"; do
    for q in "${QHAT_VALUES[@]}"; do
      if [[ $FORCE -ne 1 ]]; then
        tag=$(alpha_tag "$a"); base="output_${R}GeV_${COLL}_alpha_${tag}"
        dest="${base}/output_qhat0_${q}"
        if compgen -G "${dest}/${PART}/cent_*_*" >/dev/null 2>&1; then
          echo "[INFO] skip existing ${dest} (FORCE=0)"; continue
        fi
      fi
      run_one "$R" "$a" "$q"
    done
  done

  # restore pristine params before next energy
  cp -f "$PARAMS_BAK" "$PARAMS_FILE"
done

# final restore (paranoid)
cp -f "$PARAMS_BAK" "$PARAMS_FILE"
rm -f "$PARAMS_BAK"

echo "[INFO] all runs done."

