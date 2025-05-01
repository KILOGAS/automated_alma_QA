import numpy as np
from astropy.io import fits
from io_utils import read_fits_header, extract_beam_params, extract_pixel_size, extract_units, extract_velocity_axis
import os
from astropy.wcs import validate as wcs_validate
from astropy.wcs import WCS
from astropy.units import Unit, UnitsError
from astropy import units as u

def check_edge_emission(cube_path, edge_channels=3, threshold_sigma=3):
    """
    Check for missed emission at velocity range edges.
    Returns: dict with 'flag', 'max_edge', 'rms', 'threshold', 'details', 'error'.
    """
    result = {}
    try:
        with fits.open(cube_path, do_not_scale_image_data=True, memmap=False) as hdul:
            hdul.verify('exception')
            data = hdul[0].data
            # Assume data shape is (velocity, y, x) or (z, y, x)
            edge_start = data[:edge_channels].flatten()
            edge_end = data[-edge_channels:].flatten()
            edge_data = np.concatenate([edge_start, edge_end])
            rms = np.std(edge_data)
            max_edge = np.max(np.abs(edge_data))
            threshold = threshold_sigma * rms
            flag = max_edge > threshold
            details = {
                'max_edge': float(max_edge),
                'rms': float(rms),
                'threshold': float(threshold),
                'edge_channels': edge_channels
            }
            result.update({'flag': flag, 'details': details})
    except Exception as e:
        result['error'] = f'FITS open/verify failed: {e}'
    return result

def check_beam_and_units(cube_path):
    """Confirm beams are round and units are in Kelvin. Also verify header and units."""
    result = {}
    try:
        header = read_fits_header(cube_path)
        # Re-verify header
        try:
            with fits.open(cube_path) as hdul:
                hdul[0].verify('exception')
        except Exception as e:
            result['header_verify_error'] = str(e)
        bmaj, bmin, bpa = extract_beam_params(header)
        units = extract_units(header)
        round_beam = (bmaj is not None and bmin is not None and abs(bmaj - bmin) < 1e-6)
        kelvin_units = False
        try:
            u = Unit(units, format='fits')
            kelvin_units = u.is_equivalent(Unit('K'))
        except Exception as e:
            result['unit_parse_error'] = str(e)
        result.update({
            'round_beam': round_beam,
            'bmaj': bmaj,
            'bmin': bmin,
            'bpa': bpa,
            'units': units,
            'kelvin_units': kelvin_units
        })
    except Exception as e:
        result['error'] = f'Header/beam/unit check failed: {e}'
    return result

def measure_pixel_and_beam_size(cube_path):
    """Measure pixel size and synthesized beam size from cube header. Verify header and units."""
    result = {}
    try:
        header = read_fits_header(cube_path)
        try:
            with fits.open(cube_path) as hdul:
                hdul[0].verify('exception')
        except Exception as e:
            result['header_verify_error'] = str(e)
        cdelt1, cdelt2, cunit1, cunit2 = extract_pixel_size(header)
        bmaj, bmin, bpa = extract_beam_params(header)
        # Validate units
        try:
            Unit(cunit1, format='fits')
            Unit(cunit2, format='fits')
        except Exception as e:
            result['pixel_unit_parse_error'] = str(e)
        result.update({
            'cdelt1': cdelt1,
            'cdelt2': cdelt2,
            'cunit1': cunit1,
            'cunit2': cunit2,
            'bmaj': bmaj,
            'bmin': bmin,
            'bpa': bpa
        })
    except Exception as e:
        result['error'] = f'Pixel/beam size check failed: {e}'
    return result

def measure_rms_cube_ends(cube_path, edge_channels=3):
    """Measure rms noise at cube ends and average the values."""
    with fits.open(cube_path) as hdul:
        data = hdul[0].data
        edge_start = data[:edge_channels].flatten()
        edge_end = data[-edge_channels:].flatten()
        rms_start = np.std(edge_start)
        rms_end = np.std(edge_end)
        rms_avg = 0.5 * (rms_start + rms_end)
        return {'rms_start': float(rms_start), 'rms_end': float(rms_end), 'rms_avg': float(rms_avg)}

def measure_cube_max(cube_path):
    """Measure maximum intensity value in cube."""
    with fits.open(cube_path) as hdul:
        data = hdul[0].data
        max_val = np.nanmax(data)
        return {'cube_max': float(max_val)}

def velocity_range_nonblank(cube_path):
    """Identify velocity range of non-blank pixels, assuming velocity is 3rd axis (z). Convert channel to velocity using header. Output velocities in km/s if possible."""
    with fits.open(cube_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        # Assume data shape is (x, y, z) or (y, x, z), but velocity is always last axis
        if data.ndim == 3:
            z_axis = 2
        else:
            raise ValueError('Cube data is not 3D')
        nonblank = np.where(np.isfinite(data) & (data > 0))
        if len(nonblank[z_axis]) == 0:
            return {'vmin_chan': None, 'vmax_chan': None, 'v_channels': [], 'vmin': None, 'vmax': None, 'vunit': None}
        vchans = np.unique(nonblank[z_axis])
        vmin_chan = int(np.min(vchans))
        vmax_chan = int(np.max(vchans))
        # Convert channel to velocity using header
        crval3 = header.get('CRVAL3')
        cdelt3 = header.get('CDELT3')
        crpix3 = header.get('CRPIX3')
        cunit3 = header.get('CUNIT3', '').strip().lower()
        vmin = vmax = vunit = None
        if crval3 is not None and cdelt3 is not None and crpix3 is not None:
            vmin = crval3 + (vmin_chan + 1 - crpix3) * cdelt3
            vmax = crval3 + (vmax_chan + 1 - crpix3) * cdelt3
            vunit = cunit3
            # Convert to km/s if in m/s
            if cunit3 == 'm/s':
                vmin /= 1000.0
                vmax /= 1000.0
                vunit = 'km/s'
        return {'vmin_chan': vmin_chan, 'vmax_chan': vmax_chan, 'v_channels': vchans.tolist(), 'vmin': vmin, 'vmax': vmax, 'vunit': vunit}

def mask_nonblank_stats(mask_path):
    """Calculate number and fraction of non-blank pixels in mask cube."""
    with fits.open(mask_path) as hdul:
        data = hdul[0].data
        total = np.prod(data.shape)
        nonblank = np.sum(np.isfinite(data) & (data != 0))
        frac = nonblank / total if total > 0 else 0
        return {'mask_nonblank': int(nonblank), 'mask_total': int(total), 'mask_frac': float(frac)}

def inspect_moment_maps(ico_path, err_path=None, snr_path=None):
    """Inspect moment 0, error, and S/N maps for reasonable relative values, spatial correspondence, error variation, and outlier S/N fraction.
    Returns a dict with S/N fractions, min S/N, and flags for low S/N and outliers.
    """
    result = {}
    try:
        # Moment 0
        with fits.open(ico_path) as hdul:
            hdul[0].verify('exception')
            moment0 = hdul[0].data
        result['ico_mean'] = float(np.nanmean(moment0))
        result['ico_max'] = float(np.nanmax(moment0))
        result['ico_min'] = float(np.nanmin(moment0))
        result['ico_nan'] = int(np.isnan(moment0).sum())
        result['ico_inf'] = int(np.isinf(moment0).sum())
        # Error map
        if err_path and os.path.exists(err_path):
            with fits.open(err_path) as hdul:
                hdul[0].verify('exception')
                err = hdul[0].data
            result['err_mean'] = float(np.nanmean(err))
            result['err_max'] = float(np.nanmax(err))
            result['err_min'] = float(np.nanmin(err))
            result['err_nan'] = int(np.isnan(err).sum())
            result['err_inf'] = int(np.isinf(err).sum())
        # S/N map
        if snr_path and os.path.exists(snr_path):
            with fits.open(snr_path) as hdul:
                hdul[0].verify('exception')
                snr = hdul[0].data
            result['snr_mean'] = float(np.nanmean(snr))
            result['snr_max'] = float(np.nanmax(snr))
            result['snr_min'] = float(np.nanmin(snr))
            result['snr_nan'] = int(np.isnan(snr).sum())
            result['snr_inf'] = int(np.isinf(snr).sum())
            # S/N thresholds
            total_pix = np.isfinite(snr).sum()
            snr3 = np.sum(snr > 3)
            snr5 = np.sum(snr > 5)
            snr10 = np.sum(snr > 10)
            result['frac_snr3'] = snr3 / total_pix if total_pix else 0
            result['frac_snr5'] = snr5 / total_pix if total_pix else 0
            result['frac_snr10'] = snr10 / total_pix if total_pix else 0
            result['flag_snr_below3'] = snr3 < total_pix
            result['flag_snr_below5'] = (np.nanmin(snr) < 5)
            result['flag_snr10_outlier'] = result['frac_snr10'] > 0.02
            # Spatial correspondence: peaks in moment0 and high S/N
            if moment0.shape == snr.shape:
                peak_mask = moment0 > 0.9 * result['ico_max']
                highsnr_mask = snr > 10
                overlap = np.sum(peak_mask & highsnr_mask)
                result['peak_highsnr_overlap'] = int(overlap)
                result['flag_peak_highsnr'] = overlap > 0
            else:
                result['spatial_shape_mismatch'] = f"moment0 shape {moment0.shape} != snr shape {snr.shape}"
    except Exception as e:
        result['moment_map_error'] = str(e)
    return result

def compare_sigma_mol_to_ico(sigma_mol_path, ico_path):
    """Compare Sigma_mol map scaling to moment 0 map (factor 6.2–1.6, inclination dependent)."""
    try:
        with fits.open(sigma_mol_path) as hdul:
            sigma_mol = hdul[0].data
        with fits.open(ico_path) as hdul:
            ico = hdul[0].data
        # Convert log(Sigma_mol) to linear if needed
        if np.nanmax(sigma_mol) < 10:  # assume log10(Msun/pc^2)
            sigma_mol = 10 ** sigma_mol
        ratio = np.nanmean(sigma_mol) / np.nanmean(ico) if np.nanmean(ico) != 0 else np.nan
        flag = not (1.6 <= ratio <= 6.2)
        return {'sigma_mol_ico_ratio': float(ratio), 'flag_sigma_mol_ico': flag}
    except Exception as e:
        return {'sigma_mol_ico_ratio': np.nan, 'flag_sigma_mol_ico': True, 'sigma_mol_ico_error': str(e)}

def compare_lco_to_ico(lco_path, ico_path, pixel_area_pc2=1.0):
    """Compare L_CO map scaling to moment 0 map (should match pixel area scaling)."""
    try:
        with fits.open(lco_path) as hdul:
            lco = hdul[0].data
        with fits.open(ico_path) as hdul:
            ico = hdul[0].data
        ratio = np.nanmean(lco) / (np.nanmean(ico) * pixel_area_pc2) if np.nanmean(ico) != 0 else np.nan
        flag = not (0.9 <= ratio <= 1.1)  # allow 10% tolerance
        return {'lco_ico_ratio': float(ratio), 'flag_lco_ico': flag}
    except Exception as e:
        return {'lco_ico_ratio': np.nan, 'flag_lco_ico': True, 'lco_ico_error': str(e)}

def compare_mmol_to_lco(mmol_path, lco_path):
    """Compare Mmol map scaling to L_CO map (Mmol ≈ 6.2 × L_CO)."""
    try:
        with fits.open(mmol_path) as hdul:
            mmol = hdul[0].data
        with fits.open(lco_path) as hdul:
            lco = hdul[0].data
        # If mmol is log10, convert to linear
        
        if np.nanmax(mmol) < 10:
            mmol = 10 ** mmol
        
        ratio = np.nanmean(mmol) / (6.2 * np.nanmean(lco)) if np.nanmean(lco) != 0 else np.nan
        flag = not (0.9 <= ratio <= 1.1)  # allow 10% tolerance
        
        return {'mmol_lco_ratio': float(ratio), 'flag_mmol_lco': flag}
    except Exception as e:
        return {'mmol_lco_ratio': np.nan, 'flag_mmol_lco': True, 'mmol_lco_error': str(e)}

def check_cube_detection(cube_path, threshold_sigma=3):
    """Check for detected signal in the cube (flag if no detection above threshold)."""
    with fits.open(cube_path) as hdul:
        data = hdul[0].data
        rms = np.nanstd(data)
        max_val = np.nanmax(data)
        detected = max_val > threshold_sigma * rms
        return {'cube_detected': detected, 'cube_max': float(max_val), 'cube_rms': float(rms)}

def check_map_detection(map_path, threshold_sigma=3):
    """Check for detected signal in a moment map (flag if no detection above threshold)."""
    with fits.open(map_path) as hdul:
        data = hdul[0].data
        rms = np.nanstd(data)
        max_val = np.nanmax(data)
        detected = max_val > threshold_sigma * rms
        return {'map_detected': detected, 'map_max': float(max_val), 'map_rms': float(rms)}

def check_lco_larger_than_ico(lco_path, ico_path):
    """Check that LCO is larger than ICO (mean or sum)."""
    with fits.open(lco_path) as hdul:
        lco = hdul[0].data
    with fits.open(ico_path) as hdul:
        ico = hdul[0].data
    lco_sum = np.nansum(lco)
    ico_sum = np.nansum(ico)
    lco_mean = np.nanmean(lco)
    ico_mean = np.nanmean(ico)
    flag = not (lco_sum > ico_sum and lco_mean > ico_mean)
    return {'lco_sum': float(lco_sum), 'ico_sum': float(ico_sum), 'lco_mean': float(lco_mean), 'ico_mean': float(ico_mean), 'flag_lco_gt_ico': flag}

def check_scaling_factor_consistency(map1_path, map2_path, err1_path=None, err2_path=None):
    """Check scaling factor for each pair of maps and their errors is the same (if error maps exist)."""
    result = {}
    with fits.open(map1_path) as hdul:
        map1 = hdul[0].data
    with fits.open(map2_path) as hdul:
        map2 = hdul[0].data
    scale = np.nanmean(map2) / np.nanmean(map1) if np.nanmean(map1) != 0 else np.nan
    result['scaling_factor'] = float(scale)
    if err1_path and err2_path and os.path.exists(err1_path) and os.path.exists(err2_path):
        with fits.open(err1_path) as hdul:
            err1 = hdul[0].data
        with fits.open(err2_path) as hdul:
            err2 = hdul[0].data
        err_scale = np.nanmean(err2) / np.nanmean(err1) if np.nanmean(err1) != 0 else np.nan
        result['error_scaling_factor'] = float(err_scale)
        result['flag_scaling_consistency'] = not np.isclose(scale, err_scale, rtol=0.05)  # 5% tolerance
    else:
        result['error_scaling_factor'] = None
        result['flag_scaling_consistency'] = False
    return result

def get_map_min_max_units(map_path):
    """Return min, max, and units for a FITS map."""
    with fits.open(map_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        min_val = np.nanmin(data)
        max_val = np.nanmax(data)
        units = header.get('BUNIT', None)
        return {'min': float(min_val), 'max': float(max_val), 'units': units}

def assert_fits_unit(unit_str):
    """
    Raise a ValueError if unit_str is not FITS-compliant.
    """
    try:
        u.Unit(unit_str, format='fits', parse_strict='raise')
        return True, None
    except Exception as e:
        return False, str(e)

def assert_header_units(fits_path, expected_unit, allow_log=False):
    """Assert the header unit matches expected. Fail if log units or not as expected. Use astropy.units.Unit for parsing and FITS compliance."""
    result = {}
    try:
        with fits.open(fits_path) as hdul:
            hdul[0].verify('exception')
            header = hdul[0].header
            unit = header.get('BUNIT', '').strip()
            fail = False
            reason = ''
            fits_compliant, fits_compliance_error = assert_fits_unit(unit)
            if not fits_compliant:
                fail = True
                reason = f'Not FITS-compliant: {fits_compliance_error}'
            elif not allow_log and 'log' in unit.lower():
                fail = True
                reason = f'Unit is log: {unit}'
            else:
                try:
                    u = Unit(unit, format='fits')
                    if str(u) != expected_unit:
                        fail = True
                        reason = f'Unit mismatch: {unit} != {expected_unit}'
                except UnitsError as ue:
                    fail = True
                    reason = f'Unit parse error: {ue}'
            result.update({'unit': unit, 'expected_unit': expected_unit, 'fail_unit': fail, 'fail_reason': reason, 'fits_compliant': fits_compliant, 'fits_compliance_error': fits_compliance_error})
    except Exception as e:
        result['error'] = f'Header unit check failed: {e}'
    return result

def check_all_positive(fits_path):
    """Check all values in the FITS product are positive. Report min/max and fail if any <= 0."""
    with fits.open(fits_path) as hdul:
        data = hdul[0].data
        min_val = np.nanmin(data)
        max_val = np.nanmax(data)
        fail = min_val <= 0
        return {'min': float(min_val), 'max': float(max_val), 'fail_positive': fail}

def check_mask_nonblank(mask_path, min_frac=0.01):
    """Assert mask non-blank > min_frac (default 1%)."""
    with fits.open(mask_path) as hdul:
        data = hdul[0].data
        total = np.prod(data.shape)
        nonblank = np.sum(np.isfinite(data) & (data != 0))
        frac = nonblank / total if total > 0 else 0
        fail = frac <= min_frac
        return {'mask_nonblank': int(nonblank), 'mask_total': int(total), 'mask_frac': float(frac), 'fail_mask_frac': fail}

def run_wcs_validation(fits_path):
    """Run astropy.wcs.validate and instantiate WCS to check for required keywords and validity."""
    try:
        results = wcs_validate(fits_path)
        # Try to instantiate WCS and check for required keywords
        with fits.open(fits_path) as hdul:
            header = hdul[0].header
            wcs = WCS(header)
            required = ['CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2']
            missing = [k for k in required if k not in header]
            if missing:
                return {'wcs_validation': str(results), 'wcs_validation_fail': True, 'missing_wcs_keywords': missing}
        return {'wcs_validation': str(results), 'wcs_validation_fail': False}
    except Exception as e:
        return {'wcs_validation': str(e), 'wcs_validation_fail': True}

def check_error_map_variation(err_map_path, threshold=0.2):
    """
    Check if the error map varies by more than the given threshold (default 20%).
    Returns: dict with 'err_var_frac', 'err_var_fail', 'err_mean', 'err_std'.
    """
    try:
        with fits.open(err_map_path) as hdul:
            data = hdul[0].data
            mean = float(np.nanmean(data))
            std = float(np.nanstd(data))
            var_frac = std / mean if mean != 0 else np.nan
            fail = var_frac > threshold
            return {'err_var_frac': var_frac, 'err_var_fail': fail, 'err_mean': mean, 'err_std': std}
    except Exception as e:
        return {'err_var_frac': np.nan, 'err_var_fail': True, 'err_mean': np.nan, 'err_std': np.nan, 'err_var_error': str(e)} 