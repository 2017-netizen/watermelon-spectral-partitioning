# Data dictionary

The workbook is the numerical submission dataset. Only document creator/editor metadata were cleared; worksheet XML is unchanged.

- Thickness_spectra: header in row 1; treatment and record index in columns 1–2; absorbance at 650, 652, ..., 950 nm in columns 3–153. Total 391 records: 99 original peel, 150 at 1 mm and 142 at 5 mm. Absorbance is dimensionless.
- SSC_model_spectra: header in row 1; original record index in column 1; reference SSC in degrees Brix in column 2; the same 151 absorbance wavelengths in columns 3–153. Total 372 records from 99 fruits. Record indices are not fruit identifiers.
- Table1–Table3: manuscript dataset composition, interval properties and validation summary.
- Supp_S1_1–Supp_S1_4: analysis, preprocessing and validation settings.
- Fig3_*–Fig8_*: original numerical figure source tables, not newly refitted results.
- Index and Readme: original workbook descriptions.

The input records are already quality-screened. No additional outlier removal was performed during release preparation. Blank result cells remain blank. Market-fruit SSC references are local, whereas competition-fruit spectra share pooled fruit-level SSC references; see the manuscript limitations. No verified fruit-ID mapping is supplied, so do not infer fruit-grouped validation from row order.
