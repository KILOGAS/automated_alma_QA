import logging
from dataclasses import dataclass, field
from typing import List, Dict, Any

@dataclass
class QAGalaxyResult:
    KGAS_ID: str
    checks: Dict[str, Any] = field(default_factory=dict)

@dataclass
class QASummary:
    n_total: int
    n_checked: int
    n_passed: int
    n_flagged: int
    n_skipped: int
    skipped: List[str]
    flagged: List[str]

# --- Data computation functions ---
def compute_qa_summary(kgas_ids, results, skipped):
    n_total = len(kgas_ids)
    n_skipped = len(skipped)
    n_checked = len(results)
    flagged = [r['object_id'] for r in results if (
        r['edge_flag'] or r['flag_round_beam'] or r['flag_kelvin_units'] or r.get('flag_sigma_mol_ico', False) or r.get('flag_lco_ico', False) or r.get('flag_mmol_lco', False) or r.get('flag_cube_detected', False) or r.get('flag_ico_detected', False) or r.get('flag_lco_detected', False) or r.get('flag_lco_gt_ico', False) or r.get('flag_scaling_consistency', False) or
        any(r.get(key, {}).get('fail_unit', False) or r.get(key, {}).get('fits_compliant') is False for key in [
            'cube_unit_check', 'mask_unit_check', 'ico_unit_check', 'lco_unit_check', 'sigma_mol_unit_check', 'mmol_unit_check',
            'ico_err_unit_check', 'lco_err_unit_check', 'sigma_mol_err_unit_check', 'mmol_err_unit_check']) or
        r.get('moment_map_error') or r.get('flag_snr_below5', False)
    )]
    n_flagged = len(flagged)
    n_passed = n_checked - n_flagged
    return QASummary(n_total, n_checked, n_passed, n_flagged, n_skipped, skipped, flagged)

# --- Presentation functions ---
def log_qa_summary(summary: QASummary):
    logging.info(f"\nQA Summary:")
    logging.info(f"  Total KGAS IDs: {summary.n_total}")
    logging.info(f"  Checked: {summary.n_checked}")
    logging.info(f"  Passed: {summary.n_passed}")
    logging.info(f"  Flagged: {summary.n_flagged}")
    logging.info(f"  Skipped (not found): {summary.n_skipped}")
    if summary.skipped:
        logging.info("\nKGAS IDs skipped (cube not found): %s", ' '.join(summary.skipped))
    else:
        logging.info("\nNo KGAS IDs were skipped.")
    if summary.flagged:
        logging.info("\nFlagged Galaxies: %s", ', '.join(summary.flagged))


def log_detailed_report(results: List[Dict[str, Any]]):
    for r in results:
        logging.info(f"\nObject {r['object_id']}:")
        # Detection checks (unmasked and masked)
        if r.get('flag_cube_detected', False):
            logging.warning(f"  [FAIL] Unmasked cube detection: max={r.get('cube_max')}, rms={r.get('cube_rms')}, n_voxels_above={r.get('cube_n_voxels_above')} (min required: {r.get('cube_detection_min_voxels', 5)})")
        else:
            logging.info(f"  [PASS] Unmasked cube detection: max={r.get('cube_max')}, rms={r.get('cube_rms')}, n_voxels_above={r.get('cube_n_voxels_above')} (min required: {r.get('cube_detection_min_voxels', 5)})")
        if r.get('masked_cube_detected', False) is not None:
            if not r.get('masked_cube_detected', False):
                logging.warning(f"  [FAIL] Masked cube detection: max={r.get('masked_cube_max')}, rms={r.get('masked_cube_rms')}, n_voxels_above={r.get('masked_cube_n_voxels_above')} (min required: {r.get('cube_detection_min_voxels', 5)})")
            else:
                logging.info(f"  [PASS] Masked cube detection: max={r.get('masked_cube_max')}, rms={r.get('masked_cube_rms')}, n_voxels_above={r.get('masked_cube_n_voxels_above')} (min required: {r.get('cube_detection_min_voxels', 5)})")
        # Integrated flux comparison
        logging.info(f"  Total integrated flux (masked): {r.get('masked_flux')}")
        logging.info(f"  Total integrated flux (unmasked): {r.get('unmasked_flux')}")
        logging.info(f"  Flux ratio (masked/unmasked): {r.get('flux_ratio')}")
        if r.get('flag_flux_diff', False):
            logging.warning(f"  [FAIL] Major difference in total integrated flux between masked and unmasked cubes!")
        # LCO > ICO
        if r.get('flag_lco_gt_ico', False):
            logging.warning(f"  [FAIL] LCO > ICO: lco_sum={r.get('lco_sum')}, ico_sum={r.get('ico_sum')}, lco_mean={r.get('lco_mean')}, ico_mean={r.get('ico_mean')}")
        else:
            logging.info(f"  [PASS] LCO > ICO: lco_sum={r.get('lco_sum')}, ico_sum={r.get('ico_sum')}, lco_mean={r.get('lco_mean')}, ico_mean={r.get('ico_mean')}")
        # Scaling factor consistency
        if r.get('flag_scaling_consistency', False):
            logging.warning(f"  [FAIL] Scaling factor consistency: scaling_factor={r.get('scaling_factor')}, error_scaling_factor={r.get('error_scaling_factor')}")
        else:
            logging.info(f"  [PASS] Scaling factor consistency: scaling_factor={r.get('scaling_factor')}, error_scaling_factor={r.get('error_scaling_factor')}")
        # Min/max/units
        logging.info(f"  ICO min/max/units: {r.get('ico_min')}/{r.get('ico_max')} [{r.get('ico_units')}]")
        logging.info(f"  LCO min/max/units: {r.get('lco_min')}/{r.get('lco_max')} [{r.get('lco_units')}]")
        logging.info(f"  Sigma_mol min/max/units: {r.get('sigma_mol_min')}/{r.get('sigma_mol_max')} [{r.get('sigma_mol_units')}]")
        logging.info(f"  Mmol min/max/units: {r.get('mmol_min')}/{r.get('mmol_max')} [{r.get('mmol_units')}]")
        # Edge emission
        if r.get('edge_flag', False):
            logging.warning(f"  [FAIL] Edge emission: max_edge={r.get('max_edge')}, rms={r.get('rms')}, threshold={r.get('threshold')}, edge_channels={r.get('edge_channels')}")
        else:
            logging.info(f"  [PASS] Edge emission: max_edge={r.get('max_edge')}, rms={r.get('rms')}, threshold={r.get('threshold')}, edge_channels={r.get('edge_channels')}")
        # Beam roundness
        if r.get('flag_round_beam', False):
            logging.warning(f"  [FAIL] Beam roundness: bmaj={r.get('bmaj')}, bmin={r.get('bmin')}, bpa={r.get('bpa')}")
        else:
            logging.info(f"  [PASS] Beam roundness: bmaj={r.get('bmaj')}, bmin={r.get('bmin')}, bpa={r.get('bpa')}")
        # Kelvin units
        if r.get('flag_kelvin_units', False):
            logging.warning(f"  [FAIL] Units: units={r.get('units')}")
        else:
            logging.info(f"  [PASS] Units: units={r.get('units')}")
        # Pixel and beam size
        logging.info(f"  Pixel size: cdelt1={r.get('cdelt1')}, cdelt2={r.get('cdelt2')}, cunit1={r.get('cunit1')}, cunit2={r.get('cunit2')}")
        logging.info(f"  Beam size: bmaj={r.get('bmaj')}, bmin={r.get('bmin')}, bpa={r.get('bpa')}")
        # RMS
        logging.info(f"  RMS: rms_start={r.get('rms_start')}, rms_end={r.get('rms_end')}, rms_avg={r.get('rms_avg')}")
        # Max intensity
        logging.info(f"  Cube max intensity: {r.get('cube_max')}")
        # Velocity range
        logging.info(f"  Velocity range of non-blank pixels: vmin_chan={r.get('vmin_chan')}, vmax_chan={r.get('vmax_chan')}, v_channels={r.get('v_channels')}")
        # Mask stats
        logging.info(f"  Mask non-blank: {r.get('mask_nonblank')} / {r.get('mask_total')} (frac={r.get('mask_frac')})")
        # Moment map QA
        if 'moment_map_error' in r:
            logging.error(f"  [FAIL] Moment map inspection: {r['moment_map_error']}")
        else:
            logging.info(f"  [PASS] Moment 0 map: mean={r.get('ico_mean')}, min={r.get('ico_min')}, max={r.get('ico_max')}")
        # S/N QA
        if r.get('flag_snr_below5', False):
            logging.warning(f"  [FAIL] S/N QA: min S/N = {r.get('snr_min')}, some pixels S/N < 5")
        else:
            logging.info(f"  [PASS] S/N QA: min S/N = {r.get('snr_min')}")
        logging.info(f"    Fraction S/N > 3: {r.get('frac_snr3')}")
        logging.info(f"    Fraction S/N > 5: {r.get('frac_snr5')}")
        logging.info(f"    Fraction S/N > 10: {r.get('frac_snr10')}")
        # SNR map consistency check
        if r.get('flag_ico_snr_mismatch', False):
            logging.warning(f"  [FAIL] ICO SNR map does not match ICO/ICO_err: {r.get('ico_snr_mismatch_summary', '')}")
        if r.get('flag_lco_snr_mismatch', False):
            logging.warning(f"  [FAIL] LCO SNR map does not match LCO/LCO_err: {r.get('lco_snr_mismatch_summary', '')}")
        if r.get('flag_sigma_mol_snr_mismatch', False):
            logging.warning(f"  [FAIL] Sigma_mol SNR map does not match Sigma_mol/Sigma_mol_err: {r.get('sigma_mol_snr_mismatch_summary', '')}")
        if r.get('flag_mmol_snr_mismatch', False):
            logging.warning(f"  [FAIL] Mmol SNR map does not match Mmol/Mmol_err: {r.get('mmol_snr_mismatch_summary', '')}")
        # Sigma_mol/ICO scaling
        if r.get('flag_sigma_mol_ico', False):
            logging.warning(f"  [FAIL] Sigma_mol/ICO scaling: ratio={r.get('sigma_mol_ico_ratio')}")
        else:
            logging.info(f"  [PASS] Sigma_mol/ICO scaling: ratio={r.get('sigma_mol_ico_ratio')}")
        # LCO/ICO scaling
        if r.get('flag_lco_ico', False):
            logging.warning(f"  [FAIL] LCO/ICO scaling: ratio={r.get('lco_ico_ratio')}")
        else:
            logging.info(f"  [PASS] LCO/ICO scaling: ratio={r.get('lco_ico_ratio')}")
        # Mmol/LCO scaling
        if r.get('flag_mmol_lco', False):
            logging.warning(f"  [FAIL] Mmol/LCO scaling: ratio={r.get('mmol_lco_ratio')}")
        else:
            logging.info(f"  [PASS] Mmol/LCO scaling: ratio={r.get('mmol_lco_ratio')}")
        # WCS validation
        for key, label in [
            ('wcs_masked_cube', 'Masked Cube'),
            ('wcs_unmasked_cube', 'Unmasked Cube'),
            ('wcs_mask', 'Mask'),
            ('wcs_ico', 'ICO'),
            ('wcs_lco', 'LCO'),
            ('wcs_sigma_mol', 'Sigma_mol'),
            ('wcs_mmol', 'Mmol')]:
            wcs = r.get(key, {})
            if wcs.get('wcs_validation_fail', False):
                logging.warning(f"  [FAIL] WCS validation {label}: {wcs.get('wcs_validation')}")
            else:
                logging.info(f"  [PASS] WCS validation {label}")
        # Unit checks
        for key, label in [
            ('masked_cube_unit_check', 'Masked Cube'),
            ('unmasked_cube_unit_check', 'Unmasked Cube'),
            ('mask_unit_check', 'Mask'),
            ('ico_unit_check', 'ICO'),
            ('lco_unit_check', 'LCO'),
            ('sigma_mol_unit_check', 'Sigma_mol'),
            ('mmol_unit_check', 'Mmol'),
            ('ico_err_unit_check', 'ICO error'),
            ('lco_err_unit_check', 'LCO error'),
            ('sigma_mol_err_unit_check', 'Sigma_mol error'),
            ('mmol_err_unit_check', 'Mmol error')]:
            unit_check = r.get(key, {})
            if unit_check.get('fits_compliant') is False:
                logging.warning(f"  [FAIL] {label} unit: Not FITS-compliant: {unit_check.get('fits_compliance_error')}")
            elif unit_check.get('fail_unit', False):
                logging.warning(f"  [FAIL] {label} unit: {unit_check.get('unit')} (expected: {unit_check.get('expected_unit')}) Reason: {unit_check.get('fail_reason')}")
            elif unit_check.get('fits_compliant') is None:
                logging.info(f"  [SKIP] {label} unit: {unit_check.get('fail_reason')}")
            else:
                logging.info(f"  [PASS] {label} unit: {unit_check.get('unit')}")
        # Positivity checks
        for key, label in [
            ('ico_positive_check', 'ICO'),
            ('lco_positive_check', 'LCO'),
            ('sigma_mol_positive_check', 'Sigma_mol'),
            ('mmol_positive_check', 'Mmol')]:
            pos_check = r.get(key, {})
            if pos_check.get('fail_positive', False):
                logging.warning(f"  [FAIL] {label} positivity: min={pos_check.get('min')}, max={pos_check.get('max')}")
            else:
                logging.info(f"  [PASS] {label} positivity: min={pos_check.get('min')}, max={pos_check.get('max')}")
        # Mask non-blank check
        mask_check = r.get('mask_nonblank_check', {})
        if mask_check.get('fail_mask_frac', False):
            logging.warning(f"  [FAIL] Mask non-blank: {mask_check.get('mask_nonblank')} / {mask_check.get('mask_total')} (frac={mask_check.get('mask_frac')})")
        else:
            logging.info(f"  [PASS] Mask non-blank: {mask_check.get('mask_nonblank')} / {mask_check.get('mask_total')} (frac={mask_check.get('mask_frac')})")

# --- Usage in main.py ---
# In your main script, after collecting results:
# summary = compute_qa_summary(kgas_ids, results, skipped)
# log_qa_summary(summary)
# log_detailed_report(results)

# --- Logging setup ---
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s') 