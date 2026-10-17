#!/usr/bin/env bash
# Downloads the reference papers from E. Engel's page (Uni Frankfurt).
# They are (c) American Physical Society and are not stored in this repo.
set -euo pipefail
cd "$(dirname "$0")"
base=http://www.itp.uni-frankfurt.de/~engel/papers
curl -fsSL -o Kreibich_Gross_Engel_PRA57_138_1998.pdf "$base/kge98.pdf"
curl -fsSL -o Engel_et_al_PRA52_2750_1995.pdf        "$base/ekfmd95.pdf"
curl -fsSL -o Engel_et_al_PRA58_964_1998.pdf         "$base/efkad98.pdf"
echo "papers downloaded to $(pwd)"
