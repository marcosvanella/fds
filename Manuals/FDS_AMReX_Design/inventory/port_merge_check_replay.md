# Replay of port_merge_check over upstream history (old text heuristic against the manifest classifier)

Window: non-merge, non-root commits since the survey window start that touch a Source/*.f90 file holding a candidate routine (`docs/inventory/port_kernel_map.csv`), each classified as `COMMIT^..COMMIT`; markers pinned at 36975d7 and translated to each commit through `git diff`. Per-hunk data: `merge_check/results/replay_hunks.csv`. Reproduce: `python3 merge_check/replay.py --out DIR; python3 merge_check/summarize.py DIR/replay_hunks.csv`.

## All commits

- commits 62, hunks in candidate routines 641
- hunks in a loop nest of a routine with no marker (text heuristic kept): 371, of which old=new verdict body-only 213, surface-changing 158
- hunks in a routine that holds marked kernels: 30 (in 9 commits, 13 commit-kernel pairs)

| manifest verdict | old heuristic body-only | old heuristic surface-changing |
|---|---|---|
| body-only | 0 | 0 |
| surface-changing | 0 | 8 |
| outside-kernel | 5 | 17 |

## Excluding bulk commits (more than 500 changed lines in one file)

- commits 58, hunks in candidate routines 603
- hunks in a loop nest of a routine with no marker (text heuristic kept): 355, of which old=new verdict body-only 210, surface-changing 145
- hunks in a routine that holds marked kernels: 30 (in 9 commits, 13 commit-kernel pairs)

| manifest verdict | old heuristic body-only | old heuristic surface-changing |
|---|---|---|
| body-only | 0 | 0 |
| surface-changing | 0 | 8 |
| outside-kernel | 5 | 17 |

## Reading the numbers

- The 25 marked kernels cover 3 files (divg.f90, mass.f90, velo.f90) and a small part of the candidate code, so most hunks fall in unmarked code and stay on the old heuristic (the `unmarked-text-heuristic` rows). The old and new classifier agree there by construction.
- On routines that hold markers the old heuristic treated every hunk in a nest as a question about the nest; the manifest separates hunks inside a marked range, hunks elsewhere in the routine that leave every kernel manifest unchanged (`outside-kernel`: 22 hunks, of which the old heuristic flagged 17 as surface-changing, mostly lines that name a routine, a `DO` or a pointer elsewhere in the routine), and routine edits that do change a kernel manifest.
- Real kernel-surface events found by the replay: commit 0be09053a1 (the U, V, W work arrays gain an extra row, `-1:IBAR+1` style bounds: four kernels change argument bounds although the marked lines were not touched; the old heuristic cannot see this because the edit is in the allocation, outside the marked lines), and commit 87157833f5 (a revert of a flux limiter change: the pinned range of `flux_mw_fix_zz` no longer resolves to a translatable nest in that tree).
- No hunk in the window changed only the translated text of a marked loop, so the replay contains no `body-only` verdict on a marked kernel; that class is exercised by the fixtures in `test/test_merge_check.py` (arithmetic edit, comment edit, callee body edit).
- Limitation of the replay itself: the markers are pinned at 36975d7, the newest commit, and every replayed commit is older. Ranges are translated back through `git diff`, so a kernel that did not exist in that form at the older commit (for example `flux_mw_fix_zz` in 87157833f5) is reported as surface-changing because its range no longer parses as a nest. These rows say "this commit changed the code the marker now sits on", which is the question the daily merge check asks for a forward merge; they are not a count of kernel-surface changes that a forward merge would have met.
