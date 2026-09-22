# Release preparation audit — 2026-09-22

## Scope

Audited the submission-package MATLAB analysis and its consolidated workbook. Historical writing, translation, Word/PDF conversion and figure-preview scripts were excluded from the release. Numerical figure/table source data remain in the workbook; this release does not recreate the final manuscript layout or all figure artwork.

## Changes

1. Default model mode changed from the obsolete `full` path to the manuscript `all4test` path; accepted modes are validated before work begins.
2. Removed obsolete `full/refine/boundary` branches requiring a historical MAT file that the package cannot generate. Kept manuscript modes `all4test`, `fourrefine`, `s3boundary`, `s3test` and a reduced `smoke` mode.
3. Removed unreachable legacy workbook-reading helpers from the segmentation program.
4. Corrected comments that overstated protection against data leakage. Documented record-level splits and workflow selection limits.
5. Added a stage runner and small input/dependency/fit checks. Excluded generated logs, caches and model outputs.
6. Cleared workbook creator/last-editor properties only. All worksheet XML remained byte-identical to the input workbook.
7. Retained third-party headers and added verified libPLS/airPLS licenses. Two functions with unverified redistribution terms are external dependencies.

## Validation scope

MATLAB R2023b static analysis and small runtime checks were performed locally. No parser errors were reported. Remaining analyzer notes include legacy date functions, array growth and unused variables in the original code; these are not being presented as runtime failures. Inputs passed expected dimensions and finite-value checks. Small SPXY, SPA, airPLS, CARS, Random Frog and SVR calls were exercised with the original external dependencies.

The full 240-model comparison and repeated nested-validation grids were not rerun in this release-preparation check. The existing manuscript estimates have not been replaced. Do not interpret a smoke-check pass as certification of scientific validity or exact reproduction of every reported estimate.

No account credentials, author-private email addresses, absolute machine paths or chat instructions were found in the selected public code/data scan. Original third-party public attribution addresses remain intentionally.

The complete 500-resample segmentation stage was rerun. All eight output worksheets matched the original results within relative tolerance 1e-9 and absolute tolerance 1e-10. Retained numerical local functions across eight core files were also compared verbatim and matched the originals.
