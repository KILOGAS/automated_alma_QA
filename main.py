import os
import numpy as np
from io_utils import load_config, load_summary_table, find_data_files, find_error_map_paths
from qa_checks import (
    check_edge_emission,
    check_beam_and_units,
    measure_pixel_and_beam_size,
    measure_rms_cube_ends,
    measure_cube_max,
    velocity_range_nonblank,
    mask_nonblank_stats,
    inspect_moment_maps,
    compare_sigma_mol_to_ico,
    compare_lco_to_ico,
    compare_mmol_to_lco,
    check_cube_detection,
    check_map_detection,
    check_lco_larger_than_ico,
    check_scaling_factor_consistency,
    get_map_min_max_units,
    run_wcs_validation,
    assert_header_units,
    check_all_positive,
    check_mask_nonblank,
    extract_velocity_axis
)
from reporting import compute_qa_summary, log_qa_summary, log_detailed_report
import logging
from datetime import datetime
import yaml

CONFIG_PATH = './config.md'

def compute_percentiles(fits_path):
    try:
        from astropy.io import fits
        with fits.open(fits_path) as hdul:
            data = hdul[0].data
            percentiles = np.nanpercentile(data, [0, 16, 50, 84, 100])
            return {
                'p0': float(percentiles[0]),
                'p16': float(percentiles[1]),
                'p50': float(percentiles[2]),
                'p84': float(percentiles[3]),
                'p100': float(percentiles[4])
            }
    except Exception as e:
        return {'p0': None, 'p16': None, 'p50': None, 'p84': None, 'p100': None, 'percentile_error': str(e)}

def get_failed_tests(result):
    failed = []
    # Add all flags and unit/positivity/mask checks
    for key, label in [
        ('edge_flag', 'Edge emission'),
        ('flag_round_beam', 'Beam roundness'),
        ('flag_kelvin_units', 'Kelvin units'),
        ('flag_sigma_mol_ico', 'Sigma_mol/ICO scaling'),
        ('flag_lco_ico', 'LCO/ICO scaling'),
        ('flag_mmol_lco', 'Mmol/LCO scaling'),
        ('flag_cube_detected', 'Cube detection'),
        ('flag_ico_detected', 'ICO detection'),
        ('flag_lco_detected', 'LCO detection'),
        ('flag_lco_gt_ico', 'LCO > ICO'),
        ('flag_scaling_consistency', 'Scaling factor consistency'),
        ('flag_snr_below5', 'S/N < 5'),
    ]:
        if result.get(key, False):
            failed.append(label)
    # Unit checks
    for key, label in [
        ('cube_unit_check', 'Cube unit'),
        ('mask_unit_check', 'Mask unit'),
        ('ico_unit_check', 'ICO unit'),
        ('lco_unit_check', 'LCO unit'),
        ('sigma_mol_unit_check', 'Sigma_mol unit'),
        ('mmol_unit_check', 'Mmol unit'),
        ('ico_err_unit_check', 'ICO error unit'),
        ('lco_err_unit_check', 'LCO error unit'),
        ('sigma_mol_err_unit_check', 'Sigma_mol error unit'),
        ('mmol_err_unit_check', 'Mmol error unit')]:
        unit_check = result.get(key, {})
        if unit_check.get('fail_unit', False) or unit_check.get('fits_compliant') is False:
            failed.append(label)
    # Positivity
    for key, label in [
        ('ico_positive_check', 'ICO positivity'),
        ('lco_positive_check', 'LCO positivity'),
        ('sigma_mol_positive_check', 'Sigma_mol positivity'),
        ('mmol_positive_check', 'Mmol positivity')]:
        pos_check = result.get(key, {})
        if pos_check.get('fail_positive', False):
            failed.append(label)
    # Mask
    mask_check = result.get('mask_nonblank_check', {})
    if mask_check.get('fail_mask_frac', False):
        failed.append('Mask non-blank')
    # Moment map error
    if result.get('moment_map_error'):
        failed.append('Moment map error')
    return failed

def write_report_to_file(report_lines, log_path):
    os.makedirs(log_path, exist_ok=True)
    dt = datetime.now().strftime('%Y%m%d_%H%M%S')
    out_path = os.path.join(log_path, f'qa_report_{dt}.txt')
    with open(out_path, 'w') as f:
        for line in report_lines:
            f.write(line + '\n')
    return out_path

def main():
    config = load_config(CONFIG_PATH)
    summary_table_path = config['summary_table']
    data_root = config['data_root']
    object_ids, summary_df = load_summary_table(summary_table_path)
    results = []
    skipped = []
    for object_id in object_ids:
        files = find_data_files(config, object_id)
        unmasked_cube_path = files['unmaskedcube']
        masked_cube_path = files['maskedcube']
        mask_path = files['mask']
        ico_path = files['ico']
        lco_path = files['lco']
        sigma_mol_path = files['sigma_mol']
        mmol_path = files['mmol']
        # Find error map paths
        err_files = find_error_map_paths(files)
        ico_err_path = err_files['ico_err']
        lco_err_path = err_files['lco_err']
        mmol_err_path = err_files['mmol_err']
        sigma_mol_err_path = err_files['sigma_mol_err']
        # Check for all required files (not error maps)
        required_files = [unmasked_cube_path, masked_cube_path, mask_path, ico_path, lco_path, sigma_mol_path]
        missing = [f for f in required_files if f is None or not os.path.exists(f)]
        if missing:
            skipped.append(object_id)
            continue
        # Detection checks (use unmasked cube)
        unmasked_cube_detect_result = check_cube_detection(unmasked_cube_path, threshold_sigma=5, min_voxels=5, min_consecutive_channels=2)
        # Masked cube detection (optional, for completeness)
        masked_cube_detect_result = check_cube_detection(masked_cube_path, threshold_sigma=5, min_voxels=5, min_consecutive_channels=2)
        # LCO > ICO check
        lco_gt_ico_result = check_lco_larger_than_ico(lco_path, ico_path)
        # Detection checks for ICO and LCO
        ico_detect_result = check_map_detection(ico_path)
        lco_detect_result = check_map_detection(lco_path)
        # Scaling factor consistency (ICO/LCO, use error maps if present)
        scale_consistency_result = check_scaling_factor_consistency(
            ico_path, lco_path,
            err1_path=ico_err_path if ico_err_path and os.path.exists(ico_err_path) else None,
            err2_path=lco_err_path if lco_err_path and os.path.exists(lco_err_path) else None
        )
        # Min/max/units for each map
        ico_minmax = get_map_min_max_units(ico_path)
        lco_minmax = get_map_min_max_units(lco_path)
        sigma_mol_minmax = get_map_min_max_units(sigma_mol_path)
        mmol_minmax = get_map_min_max_units(mmol_path)
        # Percentiles for each map
        ico_percentiles = compute_percentiles(ico_path)
        lco_percentiles = compute_percentiles(lco_path)
        sigma_mol_percentiles = compute_percentiles(sigma_mol_path)
        mmol_percentiles = compute_percentiles(mmol_path)
        # Integrated flux for masked and unmasked cubes
        def total_flux(cube_path):
            try:
                from astropy.io import fits
                with fits.open(cube_path) as hdul:
                    data = hdul[0].data
                    return float(np.nansum(data))
            except Exception:
                return None
        masked_flux = total_flux(masked_cube_path)
        unmasked_flux = total_flux(unmasked_cube_path)
        flux_ratio = masked_flux / unmasked_flux if (masked_flux is not None and unmasked_flux not in (None, 0)) else None
        flag_flux_diff = False
        if flux_ratio is not None:
            flag_flux_diff = (flux_ratio < 0.8 or flux_ratio > 1.2)
        # Run all QA checks (use masked cube for QA)
        edge_result = check_edge_emission(masked_cube_path)
        beam_units_result = check_beam_and_units(masked_cube_path)
        pix_beam_result = measure_pixel_and_beam_size(masked_cube_path)
        rms_result = measure_rms_cube_ends(masked_cube_path)
        max_result = measure_cube_max(masked_cube_path)
        vel_range_result = velocity_range_nonblank(masked_cube_path)
        mask_stats = mask_nonblank_stats(mask_path) if os.path.exists(mask_path) else {'mask_nonblank': None, 'mask_total': None, 'mask_frac': None}
        # Moment map QA
        snr_path = None
        if ico_err_path and os.path.exists(ico_err_path):
            import tempfile
            from astropy.io import fits as afits
            with afits.open(ico_path) as hdul_ico, afits.open(ico_err_path) as hdul_err:
                ico_data = hdul_ico[0].data
                err_data = hdul_err[0].data
                snr_data = np.where(err_data != 0, ico_data / err_data, 0)
                with tempfile.NamedTemporaryFile(suffix='.fits', delete=False) as tmp:
                    hdu = afits.PrimaryHDU(snr_data, header=hdul_ico[0].header)
                    hdu.writeto(tmp.name, overwrite=True)
                    snr_path = tmp.name
        moment_map_result = inspect_moment_maps(ico_path, err_path=ico_err_path, snr_path=snr_path)
        sigma_mol_ico_result = compare_sigma_mol_to_ico(sigma_mol_path, ico_path)
        # lco_ico_result = compare_lco_to_ico(lco_path, ico_path, pixel_area_pc2=1.0) # requires pixel area in pc2
        mmol_lco_result = compare_mmol_to_lco(mmol_path, lco_path)
        # WCS validation for all products
        wcs_masked_cube = run_wcs_validation(masked_cube_path)
        wcs_unmasked_cube = run_wcs_validation(unmasked_cube_path)
        wcs_mask = run_wcs_validation(mask_path)
        wcs_ico = run_wcs_validation(ico_path)
        wcs_lco = run_wcs_validation(lco_path)
        wcs_sigma_mol = run_wcs_validation(sigma_mol_path)
        wcs_mmol = run_wcs_validation(mmol_path)
        # Velocity axis info (from masked cube)
        header = None
        try:
            from astropy.io import fits as afits
            with afits.open(masked_cube_path) as hdul:
                header = hdul[0].header
        except Exception:
            header = None
        if header is not None:
            cdelt3, crval3, crpix3, cunit3 = extract_velocity_axis(header)
        else:
            cdelt3 = crval3 = crpix3 = cunit3 = None
        # Header unit checks
        masked_cube_unit_check = assert_header_units(masked_cube_path, expected_unit='K')
        unmasked_cube_unit_check = assert_header_units(unmasked_cube_path, expected_unit='K')
        mask_unit_check = assert_header_units(mask_path, expected_unit='')
        ico_unit_check = assert_header_units(ico_path, expected_unit='K km/s')
        lco_unit_check = assert_header_units(lco_path, expected_unit='K km/s pc^2')
        sigma_mol_unit_check = assert_header_units(sigma_mol_path, expected_unit='Msun / pc2', allow_log=True)
        mmol_unit_check = assert_header_units(mmol_path, expected_unit='Msun', allow_log=True)
        # All positive checks
        ico_positive_check = check_all_positive(ico_path)
        lco_positive_check = check_all_positive(lco_path)
        sigma_mol_positive_check = check_all_positive(sigma_mol_path)
        mmol_positive_check = check_all_positive(mmol_path)
        # Mask non-blank check
        mask_nonblank_check = check_mask_nonblank(mask_path)
        # Error map unit checks
        ico_err_unit_check = assert_header_units(ico_err_path, expected_unit='K km/s') if ico_err_path and os.path.exists(ico_err_path) else {'unit': None, 'expected_unit': 'K km/s', 'fail_unit': None, 'fail_reason': 'File missing', 'fits_compliant': None, 'fits_compliance_error': 'File missing'}
        lco_err_unit_check = assert_header_units(lco_err_path, expected_unit='K km/s pc^2') if lco_err_path and os.path.exists(lco_err_path) else {'unit': None, 'expected_unit': 'K km/s pc^2', 'fail_unit': None, 'fail_reason': 'File missing', 'fits_compliant': None, 'fits_compliance_error': 'File missing'}
        sigma_mol_err_unit_check = assert_header_units(sigma_mol_err_path, expected_unit='Msun / pc2', allow_log=True) if sigma_mol_err_path and os.path.exists(sigma_mol_err_path) else {'unit': None, 'expected_unit': 'Msun / pc2', 'fail_unit': None, 'fail_reason': 'File missing', 'fits_compliant': None, 'fits_compliance_error': 'File missing'}
        mmol_err_unit_check = assert_header_units(mmol_err_path, expected_unit='Msun', allow_log=True) if mmol_err_path and os.path.exists(mmol_err_path) else {'unit': None, 'expected_unit': 'Msun', 'fail_unit': None, 'fail_reason': 'File missing', 'fits_compliant': None, 'fits_compliance_error': 'File missing'}
        # Aggregate all results
        result = {
            'object_id': object_id,
            **unmasked_cube_detect_result,
            'masked_cube_detected': masked_cube_detect_result.get('cube_detected', None),
            'masked_cube_max': masked_cube_detect_result.get('cube_max', None),
            'masked_cube_rms': masked_cube_detect_result.get('cube_rms', None),
            'masked_cube_n_voxels_above': masked_cube_detect_result.get('cube_n_voxels_above', None),
            'cube_detection_min_voxels': masked_cube_detect_result.get('cube_detection_min_voxels', 5),
            'masked_flux': masked_flux,
            'unmasked_flux': unmasked_flux,
            'flux_ratio': flux_ratio,
            'flag_flux_diff': flag_flux_diff,
            **ico_detect_result,
            **lco_detect_result,
            **lco_gt_ico_result,
            **scale_consistency_result,
            'ico_min': ico_minmax['min'], 'ico_max': ico_minmax['max'], 'ico_units': ico_minmax['units'],
            'lco_min': lco_minmax['min'], 'lco_max': lco_minmax['max'], 'lco_units': lco_minmax['units'],
            'sigma_mol_min': sigma_mol_minmax['min'], 'sigma_mol_max': sigma_mol_minmax['max'], 'sigma_mol_units': sigma_mol_minmax['units'],
            'mmol_min': mmol_minmax['min'], 'mmol_max': mmol_minmax['max'], 'mmol_units': mmol_minmax['units'],
            'ico_percentiles': ico_percentiles,
            'lco_percentiles': lco_percentiles,
            'sigma_mol_percentiles': sigma_mol_percentiles,
            'mmol_percentiles': mmol_percentiles,
            'edge_flag': edge_result['flag'],
            **edge_result['details'],
            **beam_units_result,
            **pix_beam_result,
            **rms_result,
            **max_result,
            **vel_range_result,
            **mask_stats,
            **moment_map_result,
            **sigma_mol_ico_result,
            # **lco_ico_result, # SEE ABOVE, requires pixel area in pc2
            **mmol_lco_result,
            'wcs_masked_cube': wcs_masked_cube,
            'wcs_unmasked_cube': wcs_unmasked_cube,
            'wcs_mask': wcs_mask,
            'wcs_ico': wcs_ico,
            'wcs_lco': wcs_lco,
            'wcs_sigma_mol': wcs_sigma_mol,
            'wcs_mmol': wcs_mmol,
            'velaxis_cdelt3': cdelt3,
            'velaxis_crval3': crval3,
            'velaxis_crpix3': crpix3,
            'velaxis_cunit3': cunit3,
            'masked_cube_unit_check': masked_cube_unit_check,
            'unmasked_cube_unit_check': unmasked_cube_unit_check,
            'mask_unit_check': mask_unit_check,
            'ico_unit_check': ico_unit_check,
            'lco_unit_check': lco_unit_check,
            'sigma_mol_unit_check': sigma_mol_unit_check,
            'mmol_unit_check': mmol_unit_check,
            'ico_positive_check': ico_positive_check,
            'lco_positive_check': lco_positive_check,
            'sigma_mol_positive_check': sigma_mol_positive_check,
            'mmol_positive_check': mmol_positive_check,
            'mask_nonblank_check': mask_nonblank_check,
            'ico_err_unit_check': ico_err_unit_check,
            'lco_err_unit_check': lco_err_unit_check,
            'sigma_mol_err_unit_check': sigma_mol_err_unit_check,
            'mmol_err_unit_check': mmol_err_unit_check
        }
        # Ensure all expected flag keys are present
        result['flag_round_beam'] = not beam_units_result.get('round_beam', False)
        result['flag_kelvin_units'] = not beam_units_result.get('kelvin_units', False)
        result['flag_cube_detected'] = not unmasked_cube_detect_result.get('cube_detected', False)
        result['flag_ico_detected'] = not ico_detect_result.get('map_detected', False)
        result['flag_lco_detected'] = not lco_detect_result.get('map_detected', False)
        result['flag_lco_gt_ico'] = lco_gt_ico_result.get('flag_lco_gt_ico', False)
        result['flag_scaling_consistency'] = scale_consistency_result.get('flag_scaling_consistency', False)
        results.append(result)
    # After results are collected
    summary = compute_qa_summary(object_ids, results, skipped)
    log_qa_summary(summary)
    log_detailed_report(results)
    # Flagged report
    flagged_lines = []
    flagged = summary.flagged
    for r in results:
        # Always flag if check_cube_detection failed
        if r.get('flag_cube_detected', False):
            failed_tests = get_failed_tests(r)
            if 'Cube detection' not in failed_tests:
                failed_tests.insert(0, 'Cube detection')
            flagged_lines.append(f"{r['object_id']}: {', '.join(failed_tests)}")
        elif r['object_id'] in flagged:
            failed_tests = get_failed_tests(r)
            flagged_lines.append(f"{r['object_id']}: {', '.join(failed_tests)}")
    if flagged_lines:
        logging.warning("\nFlagged Objects and Failed Tests:")
        for line in flagged_lines:
            logging.warning(line)
    # Write report to file if enabled
    if config['logging']['report_to_file']:
        log_path = config['logging']['log_path']
        dt = datetime.now().strftime('%Y%m%d_%H%M%S')
        out_path = os.path.join(log_path, f'qa_report_{dt}.txt')
        os.makedirs(log_path, exist_ok=True)
        with open(out_path, 'w') as f:
            f.write('QA Summary\n')
            f.write(f"Total objects: {summary.n_total}\nChecked: {summary.n_checked}\nPassed: {summary.n_passed}\nFlagged: {summary.n_flagged}\nSkipped: {summary.n_skipped}\n")
            if summary.skipped:
                f.write(f"Skipped objects: {' '.join(summary.skipped)}\n")
            if summary.flagged:
                f.write(f"Flagged objects: {', '.join(summary.flagged)}\n")
            f.write('\nFlagged Objects and Failed Tests:\n')
            for line in flagged_lines:
                f.write(line + '\n')
            f.write('\nDetailed Report:\n')
            for r in results:
                f.write(f"\nObject {r['object_id']}\n")
                for k, v in r.items():
                    if isinstance(v, dict):
                        f.write(f"  {k}: {v}\n")
                    else:
                        f.write(f"  {k}: {v}\n")
        logging.info(f"Report written to {out_path}")

if __name__ == "__main__":
    main() 