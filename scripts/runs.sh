#!/usr/bin/env bash
set -euo pipefail

# --------------------------
# User-configurable section
# --------------------------

# L_eff values to sweep (fm)
LEFF_VALS=(13.46 9.55 6.29 3.39)

# Particles: 0 = Upsilon, 1 = J/Psi
PARTICLES=(0 1)

# Alphas cases: 0 = running (4-loop, mu=d*pT), 1 = constant
ALPHAS_CASES=(0 1)

# Path to params file and executable
PARAMS_FILE="./input/params.txt"
EXEC="./quenching"

# Where to collect results (created if missing)
OUTROOT="./runs"

# If your code generates files in place, list patterns to scoop up
# into each run's folder (adjust as needed).
COLLECT_PATTERNS=("*.dat" "*.txt" "*.out" "output*" "*.log")

# --------------------------
# Helpers
# --------------------------

# In-place update of "key  value" pairs in PARAMS_FILE.
# Matches non-comment lines that start with the key.
update_param() {
  local key="$1"
  local val="$2"
  local file="$3"
  # Use sed to replace: (start)(optional spaces)(key)(spaces)(anything) -> key + single space + val
  # Skip lines that start with '//' (comments).
  # Use a temporary file for portability (macOS sed compatibility).
  local tmp
  tmp="$(mktemp)"
  awk -v k="$key" -v v="$val" '
    BEGIN { OFS=" " }
    /^\s*\/\// { print; next }                  # keep comment lines
    {
      # If the first token equals the key, replace the second token with v,
      # and rebuild the line as "key<space>value", dropping trailing comments.
      if ($1 == k) {
        $2 = v
        # Reprint as "key value" exactly; ignore anything beyond $2 on that line.
        print $1, $2
        next
      }
      print
    }
  ' "$file" > "$tmp"
  mv "$tmp" "$file"
}

# Create a descriptive directory for each run.
# Convention: runs/<particle_name>/alpha_<run|const>/L_<value>
run_dir_path() {
  local particle="$1"   # 0 or 1
  local alpha="$2"      # 0 or 1
  local leff="$3"       # float
  local pname="upsilon"
  [[ "$particle" -eq 1 ]] && pname="jpsi"
  local aname="run"
  [[ "$alpha" -eq 1 ]] && aname="const"
  echo "${OUTROOT}/${pname}/alpha_${aname}/L_${leff}"
}

# Copy collected files/patterns into run dir (ignore missing).
collect_outputs() {
  local dest="$1"
  mkdir -p "$dest"
  for pat in "${COLLECT_PATTERNS[@]}"; do
    shopt -s nullglob
    for f in $pat; do
      # Preserve relative folder structure if it's a directory match like 'output*'
      if [[ -d "$f" ]]; then
        cp -r "$f" "$dest/" || true
      else
        cp -p "$f" "$dest/" || true
      fi
    done
  done
}

# --------------------------
# Pre-flight checks
# --------------------------

if [[ ! -x "$EXEC" ]]; then
  echo "ERROR: Executable not found or not executable: $EXEC" >&2
  exit 1
fi

if [[ ! -f "$PARAMS_FILE" ]]; then
  echo "ERROR: Params file not found: $PARAMS_FILE" >&2
  exit 1
fi

mkdir -p "$OUTROOT"

# Backup original params to restore afterwards
PARAMS_BAK="${PARAMS_FILE}.bak.$(date +%Y%m%d_%H%M%S)"
cp -p "$PARAMS_FILE" "$PARAMS_BAK"
echo "Backed up params to: $PARAMS_BAK"

# --------------------------
# Main sweep
# --------------------------

for particle in "${PARTICLES[@]}"; do
  # Set particleType (0: Upsilon, 1: J/Psi)
  update_param "particleType" "$particle" "$PARAMS_FILE"

  for alpha in "${ALPHAS_CASES[@]}"; do
    # Set alphas (0: running, 1: constant)
    update_param "alphas" "$alpha" "$PARAMS_FILE"

    for leff in "${LEFF_VALS[@]}"; do
      # For AB collisions, tie LA and LB to the same effective length
      update_param "lA" "$leff" "$PARAMS_FILE"
      update_param "lB" "$leff" "$PARAMS_FILE"

      # Optionally force collisionType=0 (both) if you want to be explicit:
      # update_param "collisionType" "0" "$PARAMS_FILE"

      # Make a per-run directory and persist a copy of the exact params used
      RUN_DIR="$(run_dir_path "$particle" "$alpha" "$leff")"
      mkdir -p "$RUN_DIR"
      cp -p "$PARAMS_FILE" "${RUN_DIR}/params_used.txt"

      # Run and capture logs
      echo "=== Running: particleType=${particle}, alphas=${alpha}, L_eff=${leff} fm ==="
      LOGFILE="${RUN_DIR}/run.log"
      ERRFILE="${RUN_DIR}/run.err"

      # If your binary supports setting an outTag via env or flags, add that here.
      # Otherwise, we just run in place and collect artifacts after.
      set +e
      "$EXEC" >"$LOGFILE" 2>"$ERRFILE"
      rc=$?
      set -e

      # Collect output artifacts (edit COLLECT_PATTERNS as needed)
      collect_outputs "$RUN_DIR"

      # Quick run status note
      if [[ $rc -eq 0 ]]; then
        echo "✔ Completed: ${RUN_DIR}"
      else
        echo "✖ Failed (rc=$rc): ${RUN_DIR} — check run.err"
      fi

    done
  done
done

# Restore original params
cp -p "$PARAMS_BAK" "$PARAMS_FILE"
echo "Params restored from backup."

echo "All sweeps finished. Results under: $OUTROOT"

