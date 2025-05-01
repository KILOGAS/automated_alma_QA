# ALMA FITS Data Product QA Pipeline

## Introduction
This repository provides automated quality assurance (QA) tools for FITS data products. It verifies FITS file integrity, WCS validity, header units, S/N map statistics, and physical consistency of derived maps, ensuring robust, reproducible, and standards-compliant data products for scientific analysis.

## User Instructions
1. Edit `config.md` to set your data paths, file patterns, and logging/reporting options.
2. Place your FITS data products in the appropriate directory structure.
3. Run the QA scripts on your files.
4. Review the output report for any flagged issues or errors.

## Configuration
- All paths and file patterns are set in `config.md` (YAML frontmatter).
- `file_patterns` must include:
  - `unmaskedcube`: Full cube, used for detection and total flux (e.g., `{object_id}_co2-1_10.0kmps_12m.image.pbcor.ifumatched.fits`)
  - `maskedcube`: Pruned subcube, used for QA on masked data (e.g., `{object_id}_expanded_pruned_subcube.fits`)
  - Other products: mask, ico, lco, sigma_mol, mmol, etc.
- See the example `config.md` for details.

## QA Checks
- **FITS File Verification**: Ensures FITS files are readable and conform to standards using Astropy's verification machinery.
- **Header Verification**: Checks for required header keywords and FITS-legal formatting.
- **WCS Validation**: Validates WCS presence, required keywords, and instantiation; reports missing or malformed WCS.
- **Unit Parsing & Compliance**: Parses and validates all header unit strings using Astropy units; flags typos, non-standard, or non-FITS-compliant units for all science and error maps.
- **Unmasked Cube Detection**: Checks for detected signal in the unmasked cube. Detection now requires both >5 voxels above SNR threshold and at least 3 consecutive channels above threshold in any pixel.
- **Masked Cube QA**: All QA checks (edge emission, beam, S/N, etc.) are performed on the masked/pruned cube.
- **Masked vs Unmasked Flux**: Compares total integrated flux in masked and unmasked cubes; flags if the ratio differs by >20%.
- **Edge Emission**: Checks for missed emission at cube velocity edges.
- **Beam and Pixel Size**: Confirms round beams and correct pixel/beam units.
- **S/N Map QA**: Computes fractions of pixels with S/N > 3, 5, 10; flags low S/N and outliers. If no S/N map is provided, computes S/N as Ico / Ico_err.
- **Moment 0 & S/N Correspondence**: Verifies spatial correspondence between moment 0 peaks and high S/N regions.
- **Physical Consistency**: Checks scaling between Sigma_mol, L_CO, and moment 0 maps.
- **All Positive Check**: Ensures all values in science maps are positive.
- **Mask Non-Blank Check**: Ensures mask has sufficient non-blank pixels.
- **Percentile Reporting**: Reports 0, 16, 50, 84, 100 percentiles for ICO, LCO, Sigma_mol, and Mmol moment maps.

## Reporting Output
- Each QA function returns a dictionary with pass/fail flags, measured values, and error messages.
- The final report summarizes all checks, including detection, flux comparison, unit compliance, and percentiles, highlighting any failures or warnings for user review.
- A flagged report lists all flagged objects and the tests failed for each, including robust cube detection failures.
- If enabled in config, the report is written to a timestamped file in the specified log path.

---

For more details, see the code and docstrings in each function.

### Install dependencies from requirements.txt:

`>> pip install -r requirements.txt`
