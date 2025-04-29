# Automated QA for KILOGAS FITS Data Products

## Introduction
This repository provides automated quality assurance (QA) tools for KILOGAS ALMA data products. It verifies FITS file integrity, WCS validity, header units, S/N map statistics, and physical consistency of derived maps, ensuring robust, reproducible, and standards-compliant data products for scientific analysis.

## User Instructions
1. Place your FITS data products in the appropriate directory.
2. Run the QA scripts on your files (see example usage below).
3. Review the output report for any flagged issues or errors.

## Example Usage
In main.py
"""
# CONFIGURATION (edit these paths as needed)
SUMMARY_TABLE_PATH = './../DR1_co2-1_10.0kmps_DP_QA0_simple.csv'  # <-- update this to table with KGAS IDs column (1-N)
BASE_DIR = './../products/matched'  # absolute product directory path or relative to this script
"""

## QA Checks
- **FITS File Verification**: Ensures FITS files are readable and conform to standards using Astropy's verification machinery.
- **Header Verification**: Checks for required header keywords and FITS-legal formatting.
- **WCS Validation**: Validates WCS presence, required keywords, and instantiation; reports missing or malformed WCS.
- **Unit Parsing**: Parses and validates all header unit strings using Astropy units; flags typos or non-standard units.
- **Edge Emission**: Checks for missed emission at cube velocity edges.
- **Beam and Pixel Size**: Confirms round beams and correct pixel/beam units.
- **S/N Map QA**: Computes fractions of pixels with S/N > 3, 5, 10; flags low S/N and outliers.
- **Moment 0 & S/N Correspondence**: Verifies spatial correspondence between moment 0 peaks and high S/N regions.
- **Error Map Variation**: Ensures error map variation is <10% across the moment 0 map.
- **Physical Consistency**: Checks scaling between Sigma_mol, L_CO, and moment 0 maps.

## Reporting Output
Each QA function returns a dictionary with pass/fail flags, measured values, and error messages. The final report summarizes all checks, highlighting any failures or warnings for user review.

---

For more details, see the code and docstrings in each function. 