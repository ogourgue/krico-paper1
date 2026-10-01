#!/bin/bash
# ============================================================================
# Assemble frames/frame_NNN.png into the movie.
#
#   ./make_movie.sh                 # 10 frames per second, ~20 s
#   FPS=8 ./make_movie.sh
#
# One frame per tracking day, played straight through: no frame is held at
# the start or the end, so the larvae are moving from the first second to
# the last. H.264 in yuv420p with faststart is what LinkedIn and most players
# expect.
# ============================================================================

set -euo pipefail

FPS="${FPS:-10}"
OUT="${OUT:-outcome_movie.mp4}"

cd "$(dirname "${BASH_SOURCE[0]}")"

ffmpeg -y -framerate "${FPS}" -i frames/frame_%03d.png \
    -c:v libx264 -pix_fmt yuv420p -crf 18 -movflags +faststart \
    "${OUT}"

echo "Wrote ${OUT}"
