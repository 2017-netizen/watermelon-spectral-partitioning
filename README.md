# Watermelon spectral partitioning and SSC prediction

MATLAB code and data for **Spectral partitioning of peel thickness effects and evaluation of wavelength ranges for soluble solids prediction in watermelon**.

## Data and scope

- 391 peel-thickness spectra: original peel (99), 1 mm (150), 5 mm (142).
- 372 SSC spectra from 99 fruits; 151 wavelengths at 2 nm spacing from 650 to 950 nm.
- Cultivar: L600 (*Citrullus lanatus*). Input and figure/table source data are in `02_Data/L600_Data_Consolidated_EN_V02.xlsx`.
- Results use record-level splits. Records from the same fruit may cross folds. Preprocessing/regressor choice was made outside the outer validation loop; CARS is fitted on outer training data and is not reselected inside every inner SVR fold. Results are conditional internal estimates, not independent-fruit validation.
- The workbook does not provide a verified fruit-ID mapping. Do not invent fruit grouping from row numbers. SSC record indices <=150 denote the 2026 competition subset in the original scripts; the remaining indices denote the 2025 market subset.

## Requirements

MATLAB R2023b, Statistics and Machine Learning Toolbox, and Signal Processing Toolbox. Some libPLS branches also use `sumsqr` (Deep Learning Toolbox). The original environment is Windows; other systems have not been runtime-tested. Source files use UTF-8 and output filenames retain Chinese labels.

Read [dependency instructions](docs/DEPENDENCIES.md) before running. Two external functions, `spxy.m` and `projections_qr.m`, are intentionally not redistributed because their redistribution terms were not available with the supplied copies. Add your lawfully obtained copies to the MATLAB path. Their expected SHA-256 hashes are recorded in `docs/dependency_versions.json`.

## Download and install the two external functions

The following original-author links were checked on 2026-09-22:

| Function | Original code / download | Original method reference |
| --- | --- | --- |
| `spxy.m` | [Direct MATLAB source download](https://www.ele.ita.br/~kawakami/spa/spxy.m) from the [authors’ SPA website](https://www.ele.ita.br/~kawakami/spa/) | Galvão et al. (2005), *A method for calibration and validation subset partitioning*, Talanta 67, 736–740. [DOI: 10.1016/j.talanta.2005.03.025](https://doi.org/10.1016/j.talanta.2005.03.025) |
| `projections_qr.m` | [Original-author ZIP download: gui_spa.zip](https://www.ele.ita.br/~kawakami/spa/gui_spa.zip). Extract the file named `projections_qr.m` from this archive. | Araújo et al. (2001), *The successive projections algorithm for variable selection in spectroscopic multicomponent analysis*, Chemometrics and Intelligent Laboratory Systems 57, 65–73. [DOI: 10.1016/S0169-7439(01)00119-8](https://doi.org/10.1016/S0169-7439(01)00119-8) |

Installation:

1. Save the SPXY link as `spxy.m` (not `spxy.m.txt`). If the browser displays source text, use Save As.
2. Download and extract `gui_spa.zip`; take `projections_qr.m`. The ZIP does **not** contain `spxy.m`. The homepage’s latest `SPA_GUI.p` is a GUI, not a replacement for these two functions.
3. Put both files in `03_Code/Dependencies/` for local use, retaining their original headers. These filenames are ignored by this repository’s `.gitignore`; do not add them to a public upload without checking their redistribution terms.
4. From the repository root, run `run_reproduction('check')`. Use `which spxy -all` and `which projections_qr -all` to identify any conflicting MATLAB path entries.

Downloaded copies matched the research copies after removing whitespace. Byte-level SHA-256 values differ because line endings/whitespace differ; research-file hashes are recorded in `docs/dependency_versions.json`. Verified download hashes:

- `spxy.m`: `1d92cb30097eb3b2d77ea0a6e2ca64bd2c37119a4ec1cba4613c687737e8d911`
- `projections_qr.m` extracted from the ZIP: `0c0ae9345d51c9bd48e1c0682bcf98ca8216616a28cfcbab9a8bdd22a030e99e`

The [authors’ GUI manual](https://www.ele.ita.br/~kawakami/spa/SPA_GUI_Manual_v3p3.pdf) describes academic/non-commercial use and redistribution of SPA_GUI. This has not been treated as an MIT license or blanket permission for the two standalone research copies. They remain external dependencies.

If the download links become unavailable, use the [original project homepage](https://www.ele.ita.br/~kawakami/spa/) or contact Roberto Kawakami Harrop Galvão through the [ITA faculty directory](https://www.ele.ita.br/) (`kawakami@ita.br`). Request the specific filenames and applicable terms rather than substituting another implementation without validation.

## Run

Open this repository folder as the MATLAB current folder:

```matlab
run_reproduction('check')         % Input dimensions and small dependency/fit checks
run_reproduction('segmentation')  % 500 bootstrap resamples; candidate partitions
run_reproduction('fixed')         % 240 fixed-split candidates
run_reproduction('refine')        % Full/S3/S4 fine-grid comparisons
run_reproduction('nested')        % Repeated outer/inner CV, full spectrum and S1-S4
run_reproduction('boundary')      % Repeated CV for S3 boundary perturbations
run_reproduction('summaries')     % Requires fixed, refine, nested, boundary outputs
```

`run_reproduction('all')` runs the analysis stages in that order. Full model fitting can take substantial time. Existing files in `04_Model_results` are overwritten by the corresponding stage. Inputs are not modified. Model outputs and local logs are excluded by `.gitignore`.

The 4-range definitions used for SSC comparisons are fixed in the modeling scripts to the manuscript ranges (650–688, 690–746, 748–826, 828–950 nm). The segmentation stage provides the supporting thickness-only partition evidence; it does not automatically optimize the SSC ranges.

The legacy function name `stage2_l600_s3_models` is retained for traceability. Its default now runs `all4test`; obsolete modes depending on an absent historical segmentation file were removed. `smoke` uses the full spectrum and S3 with reduced settings and must not replace manuscript results.

## Validation and provenance

See [audit notes](docs/AUDIT.md). Model arithmetic, manuscript-mode grids, seeds, and supplied data values were preserved. The release check is a small operational test, not a rerun of all reported models. Manuscript drafting tools, personal paths, author correspondence, temporary files and historical manuscript revisions are not part of this repository.

## Citation and reuse

See `CITATION.cff`; add the final release DOI after archiving. The manuscript is not represented as published or accepted. Third-party copyrights and licenses remain in force. The authors' own code/documentation use the MIT license (`LICENSE`); data use CC BY 4.0 (`02_Data/LICENSE.txt`).
