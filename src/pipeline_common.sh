#!/usr/bin/env bash

[[ "${BASH_SOURCE[0]}" != "$0" ]] || { echo "Source this file, do not execute it." >&2; exit 1; }

if [[ -t 1 ]]; then
    _RED=$'\033[0;31m' _GREEN=$'\033[0;32m' _YELLOW=$'\033[1;33m' _NC=$'\033[0m'
else
    _RED='' _GREEN='' _YELLOW='' _NC=''
fi

log()     { echo "[$(date '+%H:%M:%S')] $*"; [[ -n "${LOG:-}" ]] && echo "[$(date '+%H:%M:%S')] $*" >> "$LOG" || true; }
die()     { echo "${_RED}FATAL ERROR:${_NC} $*" >&2; [[ -n "${LOG:-}" ]] && echo "FATAL ERROR: $*" >> "$LOG" || true; exit 1; }
warn()    { echo "${_YELLOW}WARNING:${_NC} $*" >&2;  [[ -n "${LOG:-}" ]] && echo "WARNING: $*" >> "$LOG" || true; }
success() { echo "${_GREEN}OK${_NC} $*";              [[ -n "${LOG:-}" ]] && echo "OK $*" >> "$LOG" || true; }

write_log_header() {
    local label="$1" id="$2" desc="${3:-}"
    [[ -n "${LOG:-}" ]] || return 0
    { flock 200
      printf '\n%s\n' '====================================================================='
      printf '[%s] START %s: %s' "$(date '+%H:%M:%S')" "$label" "$id"
      [[ -n "$desc" ]] && printf ' | %s' "$desc"
      printf '\n%s\n' '---------------------------------------------------------------------'
    } >> "$LOG" 200>> "$LOG"
}

tsv_update() {
    local tsv="$1" acc="$2"; shift 2
    flock -x "${tsv}.lock" python3 src/tsv_updater.py "$tsv" "$acc" "$@"
}

collect_tsv_files() {
    if [[ $# -gt 0 ]]; then
        TSV_FILES=("$@")
    else
        mapfile -t TSV_FILES < <(find "$TSV_DIR" -maxdepth 1 -name "*.tsv" | sort)
    fi
    (( ${#TSV_FILES[@]} > 0 )) || die "No .tsv files found in $TSV_DIR"
}