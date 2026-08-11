#!/usr/bin/env bash
#
# Download the February and September monthly files of the NOAA/NSIDC Climate
# Data Record of Passive Microwave Sea Ice Concentration, Version 6 (G02202 v6),
# southern hemisphere, into data/sic/.
#
# February is the southern-hemisphere ice minimum, September the maximum.
#
# The downloaded files are gitignored (~22 MB). They are not needed to reproduce
# the figure -- data/aggregated.nc carries the two climatologies -- only to
# re-run aggregate.py from scratch.
#
# sic_provenance.txt IS tracked: it records the DOI, source URL, access date and
# the sha256 of every file, so the exact inputs stay verifiable after NSIDC
# retires the version 6 distribution path.
#
# Access is anonymous HTTPS; no Earthdata account is required.
#
# Usage:
#   ./download_sic.sh

set -euo pipefail

BASE_URL="https://noaadata.apps.nsidc.org/NOAA/G02202_V6/south/monthly"
DOI="10.7265/b18j-z797"
MONTHS="02 09"
YEAR_MIN=1979
YEAR_MAX=2025

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIC_DIR="${HERE}/data/sic"
PROVENANCE="${HERE}/sic_provenance.txt"

INDEX="$(mktemp)"
LIST="$(mktemp)"
trap 'rm -f "$INDEX" "$LIST"' EXIT

mkdir -p "$SIC_DIR"

# --insecure mirrors the --no-check-certificate option documented by NSIDC.
# Drop it if your machine resolves the certificate chain cleanly.
CURL_OPTS=(-fsSL --insecure)

echo "[1/3] Listing ${BASE_URL}"
curl "${CURL_OPTS[@]}" "${BASE_URL}/" > "$INDEX"

# Keep only NetCDF files whose YYYYMM stamp falls in a wanted month of the
# wanted year range. The filename convention has changed across G02202 versions, so
# the year-month stamp is matched rather than a fixed prefix.
grep -oE 'href="[^"]+\.nc"' "$INDEX" \
    | sed -E 's/^href="//; s/"$//' \
    | while read -r name; do
        stamp="$(echo "$name" | grep -oE '(19|20)[0-9]{4}' | head -n 1 || true)"
        [ -z "$stamp" ] && continue
        year="${stamp:0:4}"
        month="${stamp:4:2}"
        case " $MONTHS " in
            *" $month "*) ;;
            *) continue ;;
        esac
        if [ "$year" -ge "$YEAR_MIN" ] && [ "$year" -le "$YEAR_MAX" ]; then
            echo "$name"
        fi
    done | sort -u > "$LIST"

n_found="$(wc -l < "$LIST" | tr -d ' ')"
n_months=$(echo "$MONTHS" | wc -w | tr -d " ")
n_expected=$(((YEAR_MAX - YEAR_MIN + 1) * n_months))
echo "  ${n_found} files matched for months [${MONTHS}] (expected ${n_expected})"
if [ "$n_found" -ne "$n_expected" ]; then
    echo "  WARNING: unexpected count. Inspect the listing before aggregating." >&2
fi

echo "[2/3] Downloading into ${SIC_DIR}"
while read -r name; do
    if [ -f "${SIC_DIR}/${name}" ]; then
        echo "  skip ${name} (already present)"
    else
        echo "  get  ${name}"
        curl "${CURL_OPTS[@]}" -o "${SIC_DIR}/${name}" "${BASE_URL}/${name}"
    fi
done < "$LIST"

echo "[3/3] Writing ${PROVENANCE}"
{
    echo "# Sea-ice concentration inputs for F1 (February sea-ice edge)"
    echo "#"
    echo "# dataset    : NOAA/NSIDC Climate Data Record of Passive Microwave"
    echo "#              Sea Ice Concentration, G02202 Version 6"
    echo "# doi        : ${DOI}"
    echo "# source_url : ${BASE_URL}"
    echo "# subset     : southern hemisphere, monthly means for months [${MONTHS}], ${YEAR_MIN}-${YEAR_MAX}"
    echo "# accessed   : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "# n_files    : ${n_found}"
    echo "#"
    echo "# The files themselves are gitignored (data/sic/). Re-download with"
    echo "# ./download_sic.sh and verify against the checksums below."
    echo "#"
    echo "# sha256  filename"
    (cd "$SIC_DIR" && find . -maxdepth 1 -name '*.nc' -print0 | sort -z \
        | xargs -0 shasum -a 256 | sed 's|\./||')
} > "$PROVENANCE"

echo "Done. ${n_found} files in data/sic/; provenance in sic_provenance.txt"
echo "Next: python aggregate.py"