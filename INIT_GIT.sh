#!/usr/bin/env bash
# Run once, locally, to turn this folder into a git repository and make the
# first commit. Then add your remote and push.
set -euo pipefail
cd "$(dirname "$0")"
git init
git add -A
git commit -m "Initial commit: SMR toolkit vs Norwegian general population

- Analysis scripts for DOB, age-at-entry, and three-date cohort inputs
- Exact Lexis person-time expansion; Poisson (Garwood) exact CIs
- SSB table 07902 reference rates (committed) + reproducible fetch script
- Base-R and ggplot2 forest plots
- Same-day follow-up and DOB-from-age handling as documented toggles
- Death-reconciliation utilities, tests, methods text, data dictionary"
echo
echo "Local repository created. Now connect a remote and push, e.g.:"
echo "  git branch -M main"
echo "  git remote add origin git@github.com:<you>/smr-norway.git"
echo "  git push -u origin main"
