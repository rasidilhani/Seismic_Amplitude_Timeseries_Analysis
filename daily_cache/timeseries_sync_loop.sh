#!/bin/bash
sync_one() {
  relpath="$1"
  destdir="$DEST_BASE/$(dirname "$relpath")"
  mkdir -p "$destdir"
  if rsync -a "$SRC_BASE/$relpath" "$destdir/" >/dev/null 2>&1; then
    echo "OK   $relpath"
  else
    echo "MISS $relpath"
  fi
}
export -f sync_one
export SRC_BASE='rsync://illslef1@core.geo.vuw.ac.nz/mlvolc/seismic_amplitude_timeseries_out'
export DEST_BASE='/mnt/c/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/daily_cache'
export RSYNC_PASSWORD='RPrs7RYfAUN0cHat'
cat '/mnt/c/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/daily_cache/timeseries_files_from.txt' | xargs -P 8 -I {} bash -c 'sync_one "$@"' _ {}
