#!/usr/bin/env bash
# Extract per-NVTX-region GPU/wall time from a Nsight Systems .sqlite export.
# Usage: nvtx_regions.sh file.sqlite [region1 region2 ...]
# With no region list, prints the top 30 regions by total time.
# Requires only sqlite3 (native on macOS); no nsys binary needed.
set -euo pipefail
db="$1"; shift || true
if [ "$#" -gt 0 ]; then
  inlist=$(printf "'%s'," "$@"); inlist="${inlist%,}"
  filter="AND text IN ($inlist)"
else
  filter=""
fi
sqlite3 -header -column "$db" "
  SELECT substr(text,1,30) AS region,
         count(*)               AS n,
         ROUND(SUM(end-start)/1e6, 2) AS total_ms,
         ROUND(AVG(end-start)/1e6, 3) AS avg_ms
  FROM NVTX_EVENTS
  WHERE end IS NOT NULL AND text IS NOT NULL $filter
  GROUP BY text
  ORDER BY total_ms DESC
  LIMIT 30;"
