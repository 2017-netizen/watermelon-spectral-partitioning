# Dependencies and attribution

## Bundled without numerical changes

- `airPLS.m`: Zhimin Zhang; BSD 3-Clause. License: `LICENSE_airPLS.txt`. Author source: https://github.com/zmzhang/airPLS .
- libPLS functions `carspls`, `pls_`, `plscv`, `plsnipals`, `pretreat`, `randomfrog_pls`, `tp`, `vip`: Hongdong Li and contributors; license copied from the author's MATLAB Central release: `LICENSE_libPLS.txt`. https://www.mathworks.com/matlabcentral/fileexchange/47767-libpls_1-95-zip .

Original headers, references and attribution addresses are retained. They are software provenance, not the manuscript authors' private contact data. The third-party code is not relicensed as the manuscript authors' code. `dependency_versions.json` records the actual supplied versions by hash; replacing them with another upstream version can change results.

## External, not redistributed

- `spxy.m`: Sample-set Partitioning based on joint X–Y distances. Reference: Galvao et al., Talanta 67 (2005), 736–740, DOI 10.1016/j.talanta.2005.03.025.
- `projections_qr.m`: QR-based projection-chain helper for the successive projections algorithm (SPA).

Verified original-author downloads (2026-09-22):

- SPXY: https://www.ele.ita.br/~kawakami/spa/spxy.m
- QR projection helper: https://www.ele.ita.br/~kawakami/spa/gui_spa.zip (extract `projections_qr.m`; SPXY is not in this ZIP).
- Author homepage: https://www.ele.ita.br/~kawakami/spa/
- GUI manual: https://www.ele.ita.br/~kawakami/spa/SPA_GUI_Manual_v3p3.pdf
- SPA original paper: https://doi.org/10.1016/S0169-7439(01)00119-8

Both downloaded source files match the research copies after removing whitespace. See README for installation and download-file hashes. The GUI manual states academic/non-commercial terms for SPA_GUI; no MIT relicensing or blanket redistribution permission for the two standalone functions is inferred. Keep the external functions out of public uploads until their applicable terms are established.

The test performed for this release used the original local copies. The public package therefore requires these two functions before the complete modeling workflow can run. No unverified substitute implementation has been inserted.
