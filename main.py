import os
from io_utils import load_summary_table, find_kgas_files, find_error_map_paths
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
    extract_velocity_axis,
    check_error_map_variation
)
from reporting import print_qa_summary, print_detailed_report, print_flagged_report

# CONFIGURATION (edit these paths as needed)
SUMMARY_TABLE_PATH = './../DR1_co2-1_10.0kmps_DP_QA0_simple.csv'  # <-- update this
BASE_DIR = './../products/matched'  # or relative to this script


def main():
    # kgas_ids are now in 'KGAS{ID}' format (e.g., 'KGAS1')
    kgas_ids, summary_df = load_summary_table(SUMMARY_TABLE_PATH)
    results = []
    skipped = []
    for kgas_id in kgas_ids:
        files = find_kgas_files(BASE_DIR, kgas_id)
        cube_path = files['cube']
        mask_path = files['mask']
        ico_path = files['ico']
        lco_path = files['lco']
        sigma_mol_path = files['sigma_mol']
        mmol_path = files['mmol']  # For now, use sigma_mol as mmol if no separate path
        # Find error map paths
        err_files = find_error_map_paths(files)
        ico_err_path = err_files['ico_err']
        lco_err_path = err_files['lco_err']
        mmol_err_path = err_files['mmol_err']
        sigma_mol_err_path = err_files['sigma_mol_err']
        # Check for all required files (not error maps)
        required_files = [cube_path, mask_path, ico_path, lco_path, sigma_mol_path]
        missing = [f for f in required_files if not os.path.exists(f)]
        if missing:
            # print(f"Skipping {kgas_id}", end="   ")
            skipped.append(kgas_id)
            continue
        # Detection checks
        cube_detect_result = check_cube_detection(cube_path)
        ico_detect_result = check_map_detection(ico_path)
        lco_detect_result = check_map_detection(lco_path)
        # LCO > ICO check
        lco_gt_ico_result = check_lco_larger_than_ico(lco_path, ico_path)
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
        if not os.path.exists(cube_path):
            # print(f"Cube not found for {kgas_id}, skipping.")
            skipped.append(kgas_id)
            continue
        # Run all QA checks
        edge_result = check_edge_emission(cube_path)
        beam_units_result = check_beam_and_units(cube_path)
        pix_beam_result = measure_pixel_and_beam_size(cube_path)
        rms_result = measure_rms_cube_ends(cube_path)
        max_result = measure_cube_max(cube_path)
        vel_range_result = velocity_range_nonblank(cube_path)
        mask_stats = mask_nonblank_stats(mask_path) if os.path.exists(mask_path) else {'mask_nonblank': None, 'mask_total': None, 'mask_frac': None}
        # Moment map QA
        moment_map_result = inspect_moment_maps(ico_path)
        sigma_mol_ico_result = compare_sigma_mol_to_ico(sigma_mol_path, ico_path)
        lco_ico_result = compare_lco_to_ico(lco_path, ico_path, pixel_area_pc2=1.0)
        mmol_lco_result = compare_mmol_to_lco(mmol_path, lco_path)
        # WCS validation for all products
        wcs_cube = run_wcs_validation(cube_path)
        wcs_mask = run_wcs_validation(mask_path)
        wcs_ico = run_wcs_validation(ico_path)
        wcs_lco = run_wcs_validation(lco_path)
        wcs_sigma_mol = run_wcs_validation(sigma_mol_path)
        wcs_mmol = run_wcs_validation(mmol_path)
        # Velocity axis info
        header = None
        try:
            with fits.open(cube_path) as hdul:
                header = hdul[0].header
        except Exception:
            header = None
        if header is not None:
            cdelt3, crval3, crpix3, cunit3 = extract_velocity_axis(header)
        else:
            cdelt3 = crval3 = crpix3 = cunit3 = None

        # Header unit checks
        cube_unit_check = assert_header_units(cube_path, expected_unit='K')
        mask_unit_check = assert_header_units(mask_path, expected_unit='')
        ico_unit_check = assert_header_units(ico_path, expected_unit='K km/s')
        lco_unit_check = assert_header_units(lco_path, expected_unit='K km/s pc2')
        sigma_mol_unit_check = assert_header_units(sigma_mol_path, expected_unit='Msun / pc2', allow_log=False)
        mmol_unit_check = assert_header_units(mmol_path, expected_unit='Msun', allow_log=False)

        # All positive checks
        ico_positive_check = check_all_positive(ico_path)
        lco_positive_check = check_all_positive(lco_path)
        sigma_mol_positive_check = check_all_positive(sigma_mol_path)
        mmol_positive_check = check_all_positive(mmol_path)

        # Mask non-blank check
        mask_nonblank_check = check_mask_nonblank(mask_path)

        # Error map unit checks
        ico_err_unit_check = assert_header_units(ico_err_path, expected_unit='K km / s') if ico_err_path and os.path.exists(ico_err_path) else {'unit': None, 'expected_unit': 'K km / s', 'fail_unit': None, 'fail_reason': 'File missing', 'fits_compliant': None, 'fits_compliance_error': 'File missing'}
        lco_err_unit_check = assert_header_units(lco_err_path, expected_unit='K km / s pc2') if lco_err_path and os.path.exists(lco_err_path) else {'unit': None, 'expected_unit': 'K km / s pc2', 'fail_unit': None, 'fail_reason': 'File missing', 'fits_compliant': None, 'fits_compliance_error': 'File missing'}
        sigma_mol_err_unit_check = assert_header_units(sigma_mol_err_path, expected_unit='Msun / pc2', allow_log=True) if sigma_mol_err_path and os.path.exists(sigma_mol_err_path) else {'unit': None, 'expected_unit': 'Msun / pc2', 'fail_unit': None, 'fail_reason': 'File missing', 'fits_compliant': None, 'fits_compliance_error': 'File missing'}
        mmol_err_unit_check = assert_header_units(mmol_err_path, expected_unit='Msun', allow_log=True) if mmol_err_path and os.path.exists(mmol_err_path) else {'unit': None, 'expected_unit': 'Msun', 'fail_unit': None, 'fail_reason': 'File missing', 'fits_compliant': None, 'fits_compliance_error': 'File missing'}

        # Error map error variation checks
        ico_err_var = check_error_map_variation(ico_err_path) if ico_err_path and os.path.exists(ico_err_path) else {'err_var_frac': None, 'err_var_fail': None, 'err_mean': None, 'err_std': None, 'err_var_error': 'File missing'}
        lco_err_var = check_error_map_variation(lco_err_path) if lco_err_path and os.path.exists(lco_err_path) else {'err_var_frac': None, 'err_var_fail': None, 'err_mean': None, 'err_std': None, 'err_var_error': 'File missing'}
        sigma_mol_err_var = check_error_map_variation(sigma_mol_err_path) if sigma_mol_err_path and os.path.exists(sigma_mol_err_path) else {'err_var_frac': None, 'err_var_fail': None, 'err_mean': None, 'err_std': None, 'err_var_error': 'File missing'}
        mmol_err_var = check_error_map_variation(mmol_err_path) if mmol_err_path and os.path.exists(mmol_err_path) else {'err_var_frac': None, 'err_var_fail': None, 'err_mean': None, 'err_std': None, 'err_var_error': 'File missing'}

        # Aggregate all results
        result = {
            'KGAS_ID': kgas_id,
            **cube_detect_result,
            **ico_detect_result,
            **lco_detect_result,
            **lco_gt_ico_result,
            **scale_consistency_result,
            'ico_min': ico_minmax['min'], 'ico_max': ico_minmax['max'], 'ico_units': ico_minmax['units'],
            'lco_min': lco_minmax['min'], 'lco_max': lco_minmax['max'], 'lco_units': lco_minmax['units'],
            'sigma_mol_min': sigma_mol_minmax['min'], 'sigma_mol_max': sigma_mol_minmax['max'], 'sigma_mol_units': sigma_mol_minmax['units'],
            'mmol_min': mmol_minmax['min'], 'mmol_max': mmol_minmax['max'], 'mmol_units': mmol_minmax['units'],
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
            **lco_ico_result,
            **mmol_lco_result,
            'wcs_cube': wcs_cube,
            'wcs_mask': wcs_mask,
            'wcs_ico': wcs_ico,
            'wcs_lco': wcs_lco,
            'wcs_sigma_mol': wcs_sigma_mol,
            'wcs_mmol': wcs_mmol,
            # New QA checks:
            'velaxis_cdelt3': cdelt3,
            'velaxis_crval3': crval3,
            'velaxis_crpix3': crpix3,
            'velaxis_cunit3': cunit3,
            'cube_unit_check': cube_unit_check,
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
            'mmol_err_unit_check': mmol_err_unit_check,
            'ico_err_var': ico_err_var,
            'lco_err_var': lco_err_var,
            'sigma_mol_err_var': sigma_mol_err_var,
            'mmol_err_var': mmol_err_var
        }
        # Add flags for failed checks
        result['flag_round_beam'] = not beam_units_result['round_beam']
        result['flag_kelvin_units'] = not beam_units_result['kelvin_units']
        # Add flags for new checks
        result['flag_cube_detected'] = not cube_detect_result['cube_detected']
        result['flag_ico_detected'] = not ico_detect_result['map_detected']
        result['flag_lco_detected'] = not lco_detect_result['map_detected']
        result['flag_lco_gt_ico'] = lco_gt_ico_result['flag_lco_gt_ico']
        result['flag_scaling_consistency'] = scale_consistency_result['flag_scaling_consistency']
        # Already included in *_result dicts
        results.append(result)
    # Optionally, save results to CSV
    # import pandas as pd
    # pd.DataFrame(results).to_csv('qa_edge_emission_results.csv', index=False)

    # After results are collected
    print_qa_summary(kgas_ids, results, skipped)
    print_detailed_report(results)
    print_flagged_report(results)

if __name__ == "__main__":
    main() 