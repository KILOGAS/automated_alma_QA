"""
Upper Limit QA Module

Measures key parameters from upper limit FITS files and cubes:
- From Ico_ul.fits: minimum value, beam parameters (bmin/bmaj), non-blank pixel count
- From cubes: RMS in lowest 300 km/s, beam parameters (bmin/bmaj)
"""

import os
import numpy as np
from astropy.io import fits


def measure_ul_fits(ul_fits_path):
    """
    Measure parameters from upper limit FITS file.
    
    Args:
        ul_fits_path: Path to Ico_ul FITS file
    
    Returns:
        dict with ul_min, ul_bmaj, ul_bmin, ul_nonblank_pixels
    """
    try:
        with fits.open(ul_fits_path) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            
            # Minimum value in image
            ul_min = float(np.nanmin(data))
            
            # Beam parameters from header
            ul_bmaj = header.get('BMAJ', None)
            ul_bmin = header.get('BMIN', None)
            
            # Number of non-blank (finite, non-NaN) pixels
            ul_nonblank = int(np.sum(np.isfinite(data)))
            
            return {
                'ul_min': ul_min,
                'ul_bmaj': ul_bmaj,
                'ul_bmin': ul_bmin,
                'ul_nonblank_pixels': ul_nonblank,
                'error': None
            }
    except Exception as e:
        return {
            'ul_min': None,
            'ul_bmaj': None,
            'ul_bmin': None,
            'ul_nonblank_pixels': None,
            'error': str(e)
        }


def measure_cube_rms_and_beam(cube_path, velocity_range_kms=300):
    """
    Measure RMS in lowest velocity range and beam parameters from cube.
    
    Args:
        cube_path: Path to cube FITS file (*image.fits, not pbcor)
        velocity_range_kms: Velocity range to use for RMS calculation (default 300 km/s)
    
    Returns:
        dict with cube_rms, cube_bmaj, cube_bmin
    """
    try:
        with fits.open(cube_path) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            
            # Get velocity axis information
            cdelt3 = header.get('CDELT3', None)  # Velocity step in m/s or km/s
            cunit3 = header.get('CUNIT3', '').strip().lower()
            
            if cdelt3 is None:
                return {
                    'cube_rms': None,
                    'cube_bmaj': None,
                    'cube_bmin': None,
                    'error': 'No CDELT3 in header'
                }
            
            # Convert velocity step to km/s if needed
            if cunit3 == 'm/s':
                cdelt3_kms = abs(cdelt3) / 1000.0
            else:
                cdelt3_kms = abs(cdelt3)
            
            # Calculate number of channels for velocity_range_kms
            n_channels = int(velocity_range_kms / cdelt3_kms)
            
            # Use the lowest n_channels (first channels in cube)
            if n_channels > data.shape[0]:
                n_channels = data.shape[0]
            
            lowest_channels = data[:n_channels]
            
            # Calculate RMS (standard deviation)
            cube_rms = float(np.nanstd(lowest_channels))
            
            # Beam parameters from header
            cube_bmaj = header.get('BMAJ', None)
            cube_bmin = header.get('BMIN', None)
            
            return {
                'cube_rms': cube_rms,
                'cube_bmaj': cube_bmaj,
                'cube_bmin': cube_bmin,
                'n_channels_used': n_channels,
                'velocity_range_kms': velocity_range_kms,
                'error': None
            }
    except Exception as e:
        return {
            'cube_rms': None,
            'cube_bmaj': None,
            'cube_bmin': None,
            'error': str(e)
        }


def check_all_upper_limits(object_id, data_root, cube_root=None):
    """
    Check upper limits for all velocity widths (10/30 kms).
    
    Args:
        object_id: Galaxy ID (e.g., 'KGAS1')
        data_root: Root directory for data products (for upper limit files)
        cube_root: Root directory for cubes (for *image.fits files)
    
    Returns:
        dict with results for each configuration
    """
    if cube_root is None:
        cube_root = data_root
    
    results = {}
    
    # Check both velocity widths
    for velocity_width in [10, 30]:
        subdir = f"{velocity_width}kms"
        
        # Paths to files
        ul_pattern = f"{object_id}_Ico_K_kms-1_ul.fits"
        ul_path = os.path.join(data_root, object_id, subdir, ul_pattern)
        
        # Try different cube patterns (*image.fits, not pbcor)
        # Include both regular and ifumatched versions
        cube_patterns = [
            f"{object_id}_co2-1_{velocity_width}.0kmps_7m+12m.image.fits",
            f"{object_id}_co2-1_{velocity_width}.0kmps_12m.image.fits",
            f"{object_id}_co2-1_{velocity_width}kmps_7m+12m.image.fits",
            f"{object_id}_co2-1_{velocity_width}kmps_12m.image.fits",
            f"{object_id}_co2-1_{velocity_width}.0kmps_7m+12m.image.ifumatched.fits",
            f"{object_id}_co2-1_{velocity_width}.0kmps_12m.image.ifumatched.fits",
            f"{object_id}_co2-1_{velocity_width}kmps_7m+12m.image.ifumatched.fits",
            f"{object_id}_co2-1_{velocity_width}kmps_12m.image.ifumatched.fits",
        ]
        
        cube_path = None
        for pattern in cube_patterns:
            # Try in cube_root
            candidate = os.path.join(cube_root, object_id, pattern)
            if os.path.exists(candidate):
                cube_path = candidate
                break
            # Also try in subdir
            candidate = os.path.join(cube_root, object_id, subdir, pattern)
            if os.path.exists(candidate):
                cube_path = candidate
                break
        
        # Check if files exist
        if not os.path.exists(ul_path):
            results[subdir] = {
                'exists': False,
                'ul_path': ul_path,
                'cube_path': cube_path,
                'error': 'Upper limit file not found'
            }
            continue
        
        # Measure from UL file
        ul_measurements = measure_ul_fits(ul_path)
        
        # Measure from cube
        if cube_path and os.path.exists(cube_path):
            cube_measurements = measure_cube_rms_and_beam(cube_path)
        else:
            cube_measurements = {
                'cube_rms': None,
                'cube_bmaj': None,
                'cube_bmin': None,
                'error': 'Cube file not found'
            }
        
        # Combine results
        result = {
            'exists': True,
            'ul_path': ul_path,
            'cube_path': cube_path,
            **ul_measurements,
            **cube_measurements
        }
        
        results[subdir] = result
    
    return results


def format_ul_result_for_table(object_id, results):
    """
    Format upper limit results for tabular display.
    
    Returns:
        dict with one entry per measurement, suitable for DataFrame
    """
    table_row = {'object_id': object_id}
    
    for config, result in results.items():
        # Extract velocity width (10 or 30)
        vel = config.replace('kms', '')
        
        if not result.get('exists', False):
            # No data
            table_row[f'{vel}_ul_min'] = None
            table_row[f'{vel}_ul_bmaj'] = None
            table_row[f'{vel}_ul_bmin'] = None
            table_row[f'{vel}_ul_nonblank_pixels'] = None
            table_row[f'{vel}_cube_rms'] = None
            table_row[f'{vel}_cube_bmaj'] = None
            table_row[f'{vel}_cube_bmin'] = None
        else:
            # From UL file
            table_row[f'{vel}_ul_min'] = result.get('ul_min')
            table_row[f'{vel}_ul_bmaj'] = result.get('ul_bmaj')
            table_row[f'{vel}_ul_bmin'] = result.get('ul_bmin')
            table_row[f'{vel}_ul_nonblank_pixels'] = result.get('ul_nonblank_pixels')
            
            # From cube
            table_row[f'{vel}_cube_rms'] = result.get('cube_rms')
            table_row[f'{vel}_cube_bmaj'] = result.get('cube_bmaj')
            table_row[f'{vel}_cube_bmin'] = result.get('cube_bmin')
    
    return table_row
