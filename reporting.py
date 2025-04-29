def print_qa_summary(kgas_ids, results, skipped):
    n_total = len(kgas_ids)
    n_skipped = len(skipped)
    n_checked = len(results)
    # Count flagged if any fail_unit or not fits_compliant in any unit check
    def any_unit_flag(r):
        for key in [
            'cube_unit_check', 'mask_unit_check', 'ico_unit_check', 'lco_unit_check', 'sigma_mol_unit_check', 'mmol_unit_check',
            'ico_err_unit_check', 'lco_err_unit_check', 'sigma_mol_err_unit_check', 'mmol_err_unit_check']:
            unit_check = r.get(key, {})
            if unit_check.get('fail_unit', False) or unit_check.get('fits_compliant') is False:
                return True
        return False
    def any_err_var_flag(r):
        for key in ['ico_err_var', 'lco_err_var', 'sigma_mol_err_var', 'mmol_err_var']:
            err_var = r.get(key, {})
            if err_var.get('err_var_fail', False):
                return True
        return False
    n_flagged = sum(1 for r in results if r['edge_flag'] or r['flag_round_beam'] or r['flag_kelvin_units'] or r.get('flag_sigma_mol_ico', False) or r.get('flag_lco_ico', False) or r.get('flag_mmol_lco', False) or r.get('flag_cube_detected', False) or r.get('flag_ico_detected', False) or r.get('flag_lco_detected', False) or r.get('flag_lco_gt_ico', False) or r.get('flag_scaling_consistency', False) or any_unit_flag(r) or any_err_var_flag(r))
    n_passed = n_checked - n_flagged
    print(f"\nQA Summary:")
    print(f"  Total KGAS IDs: {n_total}")
    print(f"  Checked: {n_checked}")
    print(f"  Passed: {n_passed}")
    print(f"  Flagged: {n_flagged}")
    print(f"  Skipped (not found): {n_skipped}")
    if skipped:
        print("\nKGAS IDs skipped (cube not found):", ' '.join(skipped))
    else:
        print("\nNo KGAS IDs were skipped.")

def print_detailed_report(results):
    print("\nDetailed QA Report for All Galaxies:")
    for r in results:
        print(f"\nGalaxy {r['KGAS_ID']}:")
        # Detection checks
        if r['flag_cube_detected']:
            print(f"  [FAIL] Cube detection: max={r['cube_max']:.4g}, rms={r['cube_rms']:.4g}")
        else:
            print(f"  [PASS] Cube detection: max={r['cube_max']:.4g}, rms={r['cube_rms']:.4g}")
        if r['flag_ico_detected']:
            print(f"  [FAIL] ICO detection: max={r['map_max']:.4g}, rms={r['map_rms']:.4g}")
        else:
            print(f"  [PASS] ICO detection: max={r['map_max']:.4g}, rms={r['map_rms']:.4g}")
        if r['flag_lco_detected']:
            print(f"  [FAIL] LCO detection: max={r['map_max']:.4g}, rms={r['map_rms']:.4g}")
        else:
            print(f"  [PASS] LCO detection: max={r['map_max']:.4g}, rms={r['map_rms']:.4g}")
        # LCO > ICO
        if r['flag_lco_gt_ico']:
            print(f"  [FAIL] LCO > ICO: lco_sum={r['lco_sum']:.4g}, ico_sum={r['ico_sum']:.4g}, lco_mean={r['lco_mean']:.4g}, ico_mean={r['ico_mean']:.4g}")
        else:
            print(f"  [PASS] LCO > ICO: lco_sum={r['lco_sum']:.4g}, ico_sum={r['ico_sum']:.4g}, lco_mean={r['lco_mean']:.4g}, ico_mean={r['ico_mean']:.4g}")
        # Scaling factor consistency
        if r['flag_scaling_consistency']:
            print(f"  [FAIL] Scaling factor consistency: scaling_factor={r['scaling_factor']:.4g}, error_scaling_factor={r['error_scaling_factor']}")
        else:
            print(f"  [PASS] Scaling factor consistency: scaling_factor={r['scaling_factor']:.4g}, error_scaling_factor={r['error_scaling_factor']}")
        if r['error_scaling_factor'] is None:
            print(f"    [NOTE] Error scaling factor is None: error maps may be missing.")
        # Min/max/units
        print(f"  ICO min/max/units: {r['ico_min']:.4g}/{r['ico_max']:.4g} [{r['ico_units']}]")
        print(f"  LCO min/max/units: {r['lco_min']:.4g}/{r['lco_max']:.4g} [{r['lco_units']}]")
        print(f"  Sigma_mol min/max/units: {r['sigma_mol_min']:.4g}/{r['sigma_mol_max']:.4g} [{r['sigma_mol_units']}]")
        print(f"  Mmol min/max/units: {r['mmol_min']:.4g}/{r['mmol_max']:.4g} [{r['mmol_units']}]")
        # Edge emission
        if r['edge_flag']:
            print(f"  [FAIL] Edge emission: max_edge={r['max_edge']:.4g}, rms={r['rms']:.4g}, threshold={r['threshold']:.4g}, edge_channels={r['edge_channels']}")
        else:
            print(f"  [PASS] Edge emission: max_edge={r['max_edge']:.4g}, rms={r['rms']:.4g}, threshold={r['threshold']:.4g}, edge_channels={r['edge_channels']}")
        # Beam roundness
        if r['flag_round_beam']:
            print(f"  [FAIL] Beam roundness: bmaj={r['bmaj']}, bmin={r['bmin']}, bpa={r['bpa']}")
        else:
            print(f"  [PASS] Beam roundness: bmaj={r['bmaj']}, bmin={r['bmin']}, bpa={r['bpa']}")
        # Kelvin units
        if r['flag_kelvin_units']:
            print(f"  [FAIL] Units: units={r['units']}")
        else:
            print(f"  [PASS] Units: units={r['units']}")
        # Pixel and beam size
        print(f"  Pixel size: cdelt1={r['cdelt1']}, cdelt2={r['cdelt2']}, cunit1={r['cunit1']}, cunit2={r['cunit2']}")
        print(f"  Beam size: bmaj={r['bmaj']}, bmin={r['bmin']}, bpa={r['bpa']}")
        # RMS
        print(f"  RMS: rms_start={r['rms_start']:.4g}, rms_end={r['rms_end']:.4g}, rms_avg={r['rms_avg']:.4g}")
        # Max intensity
        print(f"  Cube max intensity: {r['cube_max']:.4g}")
        # Velocity range
        print(f"  Velocity range of non-blank pixels: vmin_chan={r['vmin_chan']}, vmax_chan={r['vmax_chan']}, v_channels={r['v_channels']}")
        # Mask stats
        print(f"  Mask non-blank: {r['mask_nonblank']} / {r['mask_total']} (frac={r['mask_frac']})")
        # Moment map QA
        if 'moment_map_error' in r:
            print(f"  [FAIL] Moment map inspection: {r['moment_map_error']}")
        else:
            print(f"  [PASS] Moment 0 map: mean={r.get('ico_mean')}, min={r.get('ico_min')}, max={r.get('ico_max')}")
        # Sigma_mol/ICO scaling
        if r.get('flag_sigma_mol_ico', False):
            print(f"  [FAIL] Sigma_mol/ICO scaling: ratio={r.get('sigma_mol_ico_ratio')}")
        else:
            print(f"  [PASS] Sigma_mol/ICO scaling: ratio={r.get('sigma_mol_ico_ratio')}")
        # LCO/ICO scaling
        if r.get('flag_lco_ico', False):
            print(f"  [FAIL] LCO/ICO scaling: ratio={r.get('lco_ico_ratio')}")
        else:
            print(f"  [PASS] LCO/ICO scaling: ratio={r.get('lco_ico_ratio')}")
        # Mmol/LCO scaling
        if r.get('flag_mmol_lco', False):
            print(f"  [FAIL] Mmol/LCO scaling: ratio={r.get('mmol_lco_ratio')}")
        else:
            print(f"  [PASS] Mmol/LCO scaling: ratio={r.get('mmol_lco_ratio')}")
        # WCS validation
        for key, label in [
            ('wcs_cube', 'Cube'),
            ('wcs_mask', 'Mask'),
            ('wcs_ico', 'ICO'),
            ('wcs_lco', 'LCO'),
            ('wcs_sigma_mol', 'Sigma_mol'),
            ('wcs_mmol', 'Mmol')]:
            wcs = r.get(key, {})
            if wcs.get('wcs_validation_fail', False):
                print(f"  [FAIL] WCS validation {label}: {wcs.get('wcs_validation')}")
            else:
                print(f"  [PASS] WCS validation {label}")

        # Velocity axis info
        print(f"  Velocity axis: cdelt3={r.get('velaxis_cdelt3')}, crval3={r.get('velaxis_crval3')}, crpix3={r.get('velaxis_crpix3')}, cunit3={r.get('velaxis_cunit3')}")

        # Unit checks
        for key, label in [
            ('cube_unit_check', 'Cube'),
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
                print(f"  [FAIL] {label} unit: Not FITS-compliant: {unit_check.get('fits_compliance_error')}")
            elif unit_check.get('fail_unit', False):
                print(f"  [FAIL] {label} unit: {unit_check.get('unit')} (expected: {unit_check.get('expected_unit')}) Reason: {unit_check.get('fail_reason')}")
            elif unit_check.get('fits_compliant') is None:
                print(f"  [SKIP] {label} unit: {unit_check.get('fail_reason')}")
            else:
                print(f"  [PASS] {label} unit: {unit_check.get('unit')}")

        # Positivity checks
        for key, label in [
            ('ico_positive_check', 'ICO'),
            ('lco_positive_check', 'LCO'),
            ('sigma_mol_positive_check', 'Sigma_mol'),
            ('mmol_positive_check', 'Mmol')]:
            pos_check = r.get(key, {})
            if pos_check.get('fail_positive', False):
                print(f"  [FAIL] {label} positivity: min={pos_check.get('min')}, max={pos_check.get('max')}")
            else:
                print(f"  [PASS] {label} positivity: min={pos_check.get('min')}, max={pos_check.get('max')}")

        # Mask non-blank check
        mask_check = r.get('mask_nonblank_check', {})
        if mask_check.get('fail_mask_frac', False):
            print(f"  [FAIL] Mask non-blank: {mask_check.get('mask_nonblank')} / {mask_check.get('mask_total')} (frac={mask_check.get('mask_frac')})")
        else:
            print(f"  [PASS] Mask non-blank: {mask_check.get('mask_nonblank')} / {mask_check.get('mask_total')} (frac={mask_check.get('mask_frac')})")

        # Error map error variation checks
        for key, label in [
            ('ico_err_var', 'ICO error'),
            ('lco_err_var', 'LCO error'),
            ('sigma_mol_err_var', 'Sigma_mol error'),
            ('mmol_err_var', 'Mmol error')]:
            err_var = r.get(key, {})
            if err_var.get('err_var_fail') is True:
                print(f"  [FAIL] {label} error variation: std/mean={err_var.get('err_var_frac'):.3g} (std={err_var.get('err_std'):.3g}, mean={err_var.get('err_mean'):.3g})")
            elif err_var.get('err_var_fail') is False:
                print(f"  [PASS] {label} error variation: std/mean={err_var.get('err_var_frac'):.3g} (std={err_var.get('err_std'):.3g}, mean={err_var.get('err_mean'):.3g})")
            else:
                print(f"  [SKIP] {label} error variation: {err_var.get('err_var_error', 'Not checked')}")

def print_flagged_report(results):
    n_flagged = sum(1 for r in results if r['edge_flag'] or r['flag_round_beam'] or r['flag_kelvin_units'] or r.get('flag_sigma_mol_ico', False) or r.get('flag_lco_ico', False) or r.get('flag_mmol_lco', False) or r.get('flag_cube_detected', False) or r.get('flag_ico_detected', False) or r.get('flag_lco_detected', False) or r.get('flag_lco_gt_ico', False) or r.get('flag_scaling_consistency', False))
    if n_flagged > 0:
        print("\nFlagged Galaxies (Any QA Test):")
        for r in results:
            if r['edge_flag'] or r['flag_round_beam'] or r['flag_kelvin_units'] or r.get('flag_sigma_mol_ico', False) or r.get('flag_lco_ico', False) or r.get('flag_mmol_lco', False) or r.get('flag_cube_detected', False) or r.get('flag_ico_detected', False) or r.get('flag_lco_detected', False) or r.get('flag_lco_gt_ico', False) or r.get('flag_scaling_consistency', False):
                print(f"  {r['KGAS_ID']}:")
                if r.get('flag_cube_detected', False):
                    print(f"    Cube detection: max={r['cube_max']:.4g}, rms={r['cube_rms']:.4g}")
                if r.get('flag_ico_detected', False):
                    print(f"    ICO detection: max={r['map_max']:.4g}, rms={r['map_rms']:.4g}")
                if r.get('flag_lco_detected', False):
                    print(f"    LCO detection: max={r['map_max']:.4g}, rms={r['map_rms']:.4g}")
                if r.get('flag_lco_gt_ico', False):
                    print(f"    LCO > ICO: lco_sum={r['lco_sum']:.4g}, ico_sum={r['ico_sum']:.4g}, lco_mean={r['lco_mean']:.4g}, ico_mean={r['ico_mean']:.4g}")
                if r.get('flag_scaling_consistency', False):
                    print(f"    Scaling factor consistency: scaling_factor={r['scaling_factor']:.4g}, error_scaling_factor={r['error_scaling_factor']}")
                if r['edge_flag']:
                    print(f"    Edge emission: max_edge={r['max_edge']:.4g}, rms={r['rms']:.4g}, threshold={r['threshold']:.4g}, edge_channels={r['edge_channels']}")
                if r['flag_round_beam']:
                    print(f"    Beam not round: bmaj={r['bmaj']}, bmin={r['bmin']}, bpa={r['bpa']}")
                if r['flag_kelvin_units']:
                    print(f"    Units not Kelvin: units={r['units']}")
                if r.get('flag_sigma_mol_ico', False):
                    print(f"    Sigma_mol/ICO scaling: ratio={r.get('sigma_mol_ico_ratio')}")
                if r.get('flag_lco_ico', False):
                    print(f"    LCO/ICO scaling: ratio={r.get('lco_ico_ratio')}")
                if r.get('flag_mmol_lco', False):
                    print(f"    Mmol/LCO scaling: ratio={r.get('mmol_lco_ratio')}")
                # New: Unit checks
                for key, label in [
                    ('cube_unit_check', 'Cube'),
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
                        print(f"    [FAIL] {label} unit: Not FITS-compliant: {unit_check.get('fits_compliance_error')}")
                    elif unit_check.get('fail_unit', False):
                        print(f"    [FAIL] {label} unit: {unit_check.get('unit')} (expected: {unit_check.get('expected_unit')}) Reason: {unit_check.get('fail_reason')}")
                    elif unit_check.get('fits_compliant') is None:
                        print(f"    [SKIP] {label} unit: {unit_check.get('fail_reason')}")
                # New: Positivity checks
                for key, label in [
                    ('ico_positive_check', 'ICO'),
                    ('lco_positive_check', 'LCO'),
                    ('sigma_mol_positive_check', 'Sigma_mol'),
                    ('mmol_positive_check', 'Mmol')]:
                    pos_check = r.get(key, {})
                    if pos_check.get('fail_positive', False):
                        print(f"    [FAIL] {label} positivity: min={pos_check.get('min')}, max={pos_check.get('max')}")
                # New: Mask non-blank check
                mask_check = r.get('mask_nonblank_check', {})
                if mask_check.get('fail_mask_frac', False):
                    print(f"    [FAIL] Mask non-blank: {mask_check.get('mask_nonblank')} / {mask_check.get('mask_total')} (frac={mask_check.get('mask_frac')})")
                # New: Error map error variation checks
                for key, label in [
                    ('ico_err_var', 'ICO error'),
                    ('lco_err_var', 'LCO error'),
                    ('sigma_mol_err_var', 'Sigma_mol error'),
                    ('mmol_err_var', 'Mmol error')]:
                    err_var = r.get(key, {})
                    if err_var.get('err_var_fail') is True:
                        print(f"    [FAIL] {label} error variation: std/mean={err_var.get('err_var_frac'):.3g} (std={err_var.get('err_std'):.3g}, mean={err_var.get('err_mean'):.3g})")
                    elif err_var.get('err_var_fail') is False:
                        print(f"    [PASS] {label} error variation: std/mean={err_var.get('err_var_frac'):.3g} (std={err_var.get('err_std'):.3g}, mean={err_var.get('err_mean'):.3g})")
                    else:
                        print(f"    [SKIP] {label} error variation: {err_var.get('err_var_error', 'Not checked')}")
    else:
        print("\nNo galaxies were flagged by the QA tests.") 