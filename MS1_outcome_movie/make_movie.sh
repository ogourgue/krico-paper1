#!/bin/bash
# ============================================================================
# Assemble frames/frame_NNN.png into the movie.
#
#   ./make_movie.sh                 # 12 frames per second
#   FPS=10 ./make_movie.sh
#
# The first frame is held HOLD_START seconds (release positions) and the last
# HOLD_END seconds (the final outcome map), so the opening and the result can
# be read. H.264 in yuv420p with faststart is what
# LinkedIn and most players expect.
# ============================================================================

set -euo pipefail

FPS="${FPS:-12}"
HOLD_START="${HOLD_START:-1.5}"
HOLD_END="${HOLD_END:-4}"
OUT="${OUT:-outcome_movie.mp4}"

cd "$(dirname "${BASH_SOURCE[0]}")"

ffmpeg -y -framerate "${FPS}" -i frames/frame_%03d.png \
    -vf "tpad=start_mode=clone:start_duration=${HOLD_START}:stop_mode=clone:stop_duration=${HOLD_END}" \
    -c:v libx264 -pix_fmt yuv420p -crf 18 -movflags +faststart \
    "${OUT}"

echo "Wrote ${OUT}"
