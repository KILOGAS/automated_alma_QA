# Upper Limit Measurements - User Guide

## Overview

The upper limit QA system measures specific parameters from upper limit FITS files and cubes, outputting them to a clean CSV file for analysis.

## Quick Start

```bash
cd /Users/thbrown/kilogas/qa/alma_product_qa/automated_alma_QA
python main_upper_limit.py
```

## Measurements

### From `Ico_ul.fits` Files

For each galaxy and velocity width (10 km/s, 30 km/s):

1. **Minimum value** (`ul_min`)
   - The minimum value across all pixels in the upper limit map
   - Units: K km/s
   - Calculation: `np.nanmin(data)`

2. **Beam major axis** (`ul_bmaj`)
   - BMAJ parameter from FITS header
   - Units: degrees
   - Header keyword: `BMAJ`

3. **Beam minor axis** (`ul_bmin`)
   - BMIN parameter from FITS header
   - Units: degrees
   - Header keyword: `BMIN`

4. **Non-blank pixels** (`ul_nonblank_pixels`)
   - Count of finite (non-NaN) pixels in the map
   - Calculation: `np.sum(np.isfinite(data))`

### From Cube Files (`*image.fits`, NOT pbcor)

For each galaxy and velocity width:

1. **RMS in lowest 300 km/s** (`cube_rms`)
   - Standard deviation in the lowest 300 km/s velocity range
   - Units: K
   - Calculation: `np.nanstd(data[:n_channels])`
   - Where `n_channels = 300 / channel_width`

2. **Beam major axis** (`cube_bmaj`)
   - BMAJ parameter from cube FITS header
   - Units: degrees
   - Header keyword: `BMAJ`

3. **Beam minor axis** (`cube_bmin`)
   - BMIN parameter from cube FITS header
   - Units: degrees
   - Header keyword: `BMIN`

## CSV Output Format

The output CSV contains these columns:

```
object_id,
10_ul_min, 10_ul_bmaj, 10_ul_bmin, 10_ul_nonblank_pixels, 10_cube_rms, 10_cube_bmaj, 10_cube_bmin,
30_ul_min, 30_ul_bmaj, 30_ul_bmin, 30_ul_nonblank_pixels, 30_cube_rms, 30_cube_bmaj, 30_cube_bmin
```

### Column Naming Convention

- Prefix `10_` = 10 km/s velocity width
- Prefix `30_` = 30 km/s velocity width
- `ul_` = from upper limit FITS file
- `cube_` = from cube FITS file

### Example Row

```csv
KGAS1,1.918858,0.000254,0.000211,8358,0.034327,0.000254,0.000211,1.380015,0.000285,0.000250,8358,0.014138,0.000285,0.000250
```

Interpreted as:
- **Galaxy:** KGAS1
- **10 km/s:**
  - UL min: 1.918858 K km/s
  - UL beam: 0.000254° × 0.000211°
  - UL non-blank pixels: 8358
  - Cube RMS: 0.034327 K
  - Cube beam: 0.000254° × 0.000211°
- **30 km/s:**
  - UL min: 1.380015 K km/s
  - UL beam: 0.000285° × 0.000250°
  - UL non-blank pixels: 8358
  - Cube RMS: 0.014138 K
  - Cube beam: 0.000285° × 0.000250°

## File Locations

### Input Files

The system looks for:

**Upper limit files:**
```
{data_root}/{object_id}/{velocity}kms/{object_id}_Ico_K_kms-1_ul.fits
```

**Cube files (tries in order):**
```
{cube_root}/{object_id}/{object_id}_co2-1_{velocity}.0kmps_7m+12m.image.fits
{cube_root}/{object_id}/{object_id}_co2-1_{velocity}.0kmps_12m.image.fits
{cube_root}/{object_id}/{velocity}kms/{object_id}_co2-1_{velocity}.0kmps_7m+12m.image.fits
{cube_root}/{object_id}/{velocity}kms/{object_id}_co2-1_{velocity}.0kmps_12m.image.fits
```

**Note:** The system looks for `*image.fits` files (NOT `*image.pbcor.fits`)

### Output Files

Results are saved to `logs/` directory:

- **CSV:** `upper_limit_qa_YYYYMMDD_HHMMSS.csv` (for analysis)
- **Text:** `upper_limit_qa_YYYYMMDD_HHMMSS.txt` (human-readable)

## Configuration

Edit `config_upper_limit.md`:

```yaml
data_root: ../../../products/v1.1/original/by_galaxy
cube_root: ../../../cubes/v1.0/original
summary_table: ../sample/DR1_co2-1_10.0kmps_DP_QA0_simple.csv
```

- `data_root`: Where upper limit FITS files are located
- `cube_root`: Where cube files are located
- `summary_table`: CSV file with list of galaxy IDs

## Using the CSV Output

### In Excel

1. Open the CSV file in Excel
2. Sort/filter by any column
3. Create scatter plots (e.g., 10 km/s vs 30 km/s values)
4. Calculate ratios or differences

### In Python

```python
import pandas as pd

# Load data
df = pd.read_csv('logs/upper_limit_qa_20251218_152823.csv')

# Compare 10 vs 30 km/s UL minimums
df['ul_ratio'] = df['10_ul_min'] / df['30_ul_min']

# Check beam consistency
df['beam_match_10'] = (df['10_ul_bmaj'] == df['10_cube_bmaj'])

# Find outliers
outliers = df[df['10_cube_rms'] > 0.05]

# Plot
import matplotlib.pyplot as plt
plt.scatter(df['10_ul_min'], df['30_ul_min'])
plt.xlabel('10 km/s UL min')
plt.ylabel('30 km/s UL min')
plt.show()
```

## Detailed Report

In addition to the CSV, a detailed text report is generated showing all measurements for each galaxy:

```
============================================================
Object: KGAS1
============================================================

  10kms:
  --------------------------------------------------------
    FROM ICO_UL.FITS:
      Minimum value:       1.918858
      BMAJ:                0.00025399
      BMIN:                0.00021130
      Non-blank pixels:    8358

    FROM CUBE (*IMAGE.FITS):
      RMS (lowest 300 km/s): 0.03432700 K
      Channels used:       27
      BMAJ:                0.00025399
      BMIN:                0.00021130
```

## Troubleshooting

### "No upper limit files found"

Check that files exist at:
```
{data_root}/{object_id}/10kms/{object_id}_Ico_K_kms-1_ul.fits
{data_root}/{object_id}/30kms/{object_id}_Ico_K_kms-1_ul.fits
```

### "Cube file not found"

The system looks for `*image.fits` files (NOT `*image.pbcor.fits`).

Check that files exist at:
```
{cube_root}/{object_id}/{object_id}_co2-1_10.0kmps_7m+12m.image.fits
{cube_root}/{object_id}/{object_id}_co2-1_30.0kmps_7m+12m.image.fits
```

Or in subdirectories:
```
{cube_root}/{object_id}/10kms/{object_id}_co2-1_10.0kmps_7m+12m.image.fits
```

### Cube RMS values seem wrong

The RMS is calculated from the **lowest 300 km/s** of the cube. This assumes:
- The spectral axis is in velocity units (m/s or km/s)
- The lowest velocity channels are line-free (noise only)
- CDELT3 header keyword exists

If your cube has different characteristics, you may need to adjust the `velocity_range_kms` parameter in the code.

### Beam values differ between UL and cube

This can happen if:
- The upper limit map was smoothed/convolved
- Different arrays were used (7m+12m vs 12m only)
- The cube was processed differently

This is not necessarily an error - just indicates different beam sizes.

## Summary

The system provides a clean CSV with exactly the measurements you need:

✅ **From UL files:** min value, bmaj, bmin, non-blank pixels  
✅ **From cubes:** RMS (lowest 300 km/s), bmaj, bmin  
✅ **For both:** 10 km/s and 30 km/s velocity widths  
✅ **Output:** CSV for analysis + detailed text report  

Run `python main_upper_limit.py` to generate the reports!

