#!/bin/bash
# ============================================================================
# Assemble frames/frame_NNN.png into the movie.
#
#   ./make_movie.sh                 # 10 frames per second, ~20 s
#   FPS=8 ./make_movie.sh
#   FRAMES=series/2010/frames OUT=series/outcome_movie_2010.mp4 ./make_movie.sh
#
# FRAMES and OUT may be absolute or relative to this folder. FFMPEG names the
# ffmpeg binary when it is not on the PATH (job_series.sh sets it on the HPC).
#
# One frame per tracking day, played straight through: no frame is held at
# the start or the end, so the larvae are moving from the first second to
# the last. H.264 in yuv420p with faststart is what LinkedIn and most players
# expect.
# ============================================================================

set -euo pipefail

FPS="${FPS:-10}"
OUT="${OUT:-outcome_movie.mp4}"
FRAMES="${FRAMES:-frames}"
FFMPEG="${FFMPEG:-ffmpeg}"

cd "$(dirname "${BASH_SOURCE[0]}")"

"${FFMPEG}" -y -framerate "${FPS}" -i "${FRAMES}/frame_%03d.png" \
    -c:v libx264 -pix_fmt yuv420p -crf 18 -movflags +faststart \
    "${OUT}"

echo "Wrote ${OUT}"
