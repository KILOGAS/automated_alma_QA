---
data_root: ../../../products/v0.1/matched/by_galaxy
cube_root: ../../../cubes/v1.0/matched
summary_table: ../sample/DR1_co2-1_10.0kmps_DP_QA0_simple.csv
file_patterns:
  unmaskedcube:
    - "{object_id}_co2-1_10.0kmps_7m+12m.image.pbcor.ifumatched.fits"
    - "{object_id}_co2-1_10.0kmps_12m.image.pbcor.ifumatched.fits"
  maskedcube: "{object_id}_clipped_cube.fits"
  mask: "{object_id}_mask_cube.fits"
  ico: "{object_id}_Ico_K_kms-1.fits"
  lco: "{object_id}_Lco_K_kms-1_pc2.fits"
  sigma_mol: "{object_id}_mmol_pc-2.fits"
  mmol: "{object_id}_mmol_pix-1.fits"
logging:
  log_path: ./../logs
  report_to_file: true
---
# Configuration for Generic FITS QA Pipeline
#
# data_root: Root directory for data products
# summary_table: Path to summary table (CSV or Excel)
# file_patterns: Patterns for each data product, use {object_id} as placeholder
#   unmaskedcube: Full cube, used for detection and total flux
#   maskedcube: Pruned subcube, used for QA on masked data
# logging.log_path: Directory for log/report output (required if logging.enabled)
# logging.report_to_file: If true, write the report to a timestamped file in log_path 