#!/usr/bin/env bash
set -euo pipefail

# --- User knobs (set via env) ---
PARAMS_FILE="${PARAMS_FILE:-input/params.txt}"
QUENCH_BIN="${QUENCH_BIN:-./quenching}"
QHAT_VALUES=(${QHAT_VALUES:-0.05 0.075 0.09})
ALPHA_MODES=(${ALPHA_MODES:-1 0})        # 1=constant, 0=running
FORCE="${FORCE:-0}"                      # 1 to overwrite
DRYRUN="${DRYRUN:-0}"

read_kv() { awk -v k="$1" '
  /^[[:space:]]*\/\//{next} /^[[:space:]]*$/{next} $1==k{print $2; exit}' "$PARAMS_FILE"; }

ROOTSNN=$(read_kv rootsnn)
CTYPE=$(read_kv collisionType); case "${CTYPE:-1}" in 1) COLL=pPb;; 2) COLL=AB;; 0) COLL=both;; *) COLL=unknown;; esac
PTYPE=$(read_kv particleType); case "${PTYPE:-1}" in 1) PART=JPsi;; 0) PART=Upsilon;; *) PART=Unknown;; esac

alpha_tag(){ [[ $1 -eq 1 ]] && echo constant || echo running; }

run_one() {
  local a="$1" q="$2" tag=$(alpha_tag "$a")
  local base="output_${ROOTSNN}GeV_${COLL}_alpha_${tag}"
  local dest="${base}/output_qhat0_${q}"
  echo "[INFO] a=$a ($tag), qhat0=$q → $dest"
  [[ $DRYRUN -eq 1 ]] && return 0

  mkdir -p "$dest"
  "$QUENCH_BIN" -alphas "$a" -qhat0 "$q" | tee "${dest}/quenching.log"

  [[ -d output ]] || { echo "ERROR: 'output' dir not produced"; exit 2; }
  # move contents; keep layout as quenching produces it
  shopt -s dotglob && mv output/* "$dest"/ && rmdir output || true

  cp -f "$PARAMS_FILE" "${dest}/params.txt"
  printf "rootsnn:%s\ncollision:%s\nparticle:%s\nalphas:%s\nqhat0:%s\n" \
    "$ROOTSNN" "$COLL ($CTYPE)" "$PART ($PTYPE)" "$a" "$q" > "${dest}/overrides.txt"
}

for a in "${ALPHA_MODES[@]}"; do
  for q in "${QHAT_VALUES[@]}"; do
    if [[ $FORCE -ne 1 ]]; then
      tag=$(alpha_tag "$a"); base="output_${ROOTSNN}GeV_${COLL}_alpha_${tag}"
      dest="${base}/output_qhat0_${q}"
      if compgen -G "${dest}/${PART}/cent_*_*" >/dev/null 2>&1; then
        echo "[INFO] skip existing ${dest} (FORCE=0)"; continue
      fi
    fi
    run_one "$a" "$q"
  done
done
echo "[INFO] all runs done."

