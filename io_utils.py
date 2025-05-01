import os
import pandas as pd
from astropy.io import fits

def load_summary_table(summary_path):
    """Load summary table and return list of KGAS IDs to process (as 'KGAS{ID}')."""
    df = pd.read_csv(summary_path) if summary_path.endswith('.csv') else pd.read_excel(summary_path)
    # Ensure KGAS_ID is int, then create formatted string
    df['KGAS_ID_INT'] = df['KGAS_ID']
    df['KGAS_ID'] = df['KGAS_ID_INT'].apply(lambda x: f"KGAS{int(x)}")
    return df['KGAS_ID'].tolist(), df

def find_kgas_files(base_dir, kgas_id):
    """Find relevant FITS files for a given KGAS ID."""
    kgas_dir = os.path.join(base_dir, kgas_id)
    files = {}
    # Cube
    files['cube'] = os.path.join(kgas_dir, f"{kgas_id}_expanded_pruned_subcube.fits")
    # Mask
    files['mask'] = os.path.join(kgas_dir, f"{kgas_id}_mask_cube.fits")
    # Moment maps
    moment_dir = os.path.join(kgas_dir, 'moment_maps')
    files['ico'] = os.path.join(moment_dir, f"{kgas_id}_Ico_K_kms-1.fits")
    files['lco'] = os.path.join(moment_dir, f"{kgas_id}_Lco_K_kms-1_pc2.fits")
    files['sigma_mol'] = os.path.join(moment_dir, f"{kgas_id}_mmol_pc-2.fits")
    files['mmol'] = os.path.join(moment_dir, f"{kgas_id}_mmol_pix-1.fits")
    return files

def read_fits_header(fits_path):
    """Read FITS header and return as a dictionary."""
    with fits.open(fits_path) as hdul:
        header = hdul[0].header
        return dict(header)

def extract_beam_params(header):
    """Extract synthesized beam parameters from FITS header."""
    bmaj = header.get('BMAJ', None)
    bmin = header.get('BMIN', None)
    bpa = header.get('BPA', None)
    return bmaj, bmin, bpa

def extract_pixel_size(header):
    """Extract pixel size in degrees (CDELT1, CDELT2) from FITS header."""
    cdelt1 = header.get('CDELT1', None)
    cdelt2 = header.get('CDELT2', None)
    cunit1 = header.get('CUNIT1', 'deg')
    cunit2 = header.get('CUNIT2', 'deg')
    return cdelt1, cdelt2, cunit1, cunit2

def extract_units(header):
    """Extract data units from FITS header."""
    return header.get('BUNIT', None)

def extract_velocity_axis(header):
    """Extract velocity axis info (CDELT3, CRVAL3, CRPIX3, CUNIT3) from FITS header."""
    cdelt3 = header.get('CDELT3', None)
    crval3 = header.get('CRVAL3', None)
    crpix3 = header.get('CRPIX3', None)
    cunit3 = header.get('CUNIT3', None)
    return cdelt3, crval3, crpix3, cunit3

def find_error_map_paths(files):
    """Given a files dict, return error map paths for ICO, LCO, SIGMA_MOL, MMOL."""
    ico_err = files['ico'].replace('.fits', '_err.fits') if files.get('ico') else None
    lco_err = files['lco'].replace('.fits', '_err.fits') if files.get('lco') else None
    sigma_mol_err = files['sigma_mol'].replace('.fits', '_err.fits') if files.get('sigma_mol') else None
    mmol_err = files['mmol'].replace('.fits', '_err.fits') if files.get('mmol') else None
    return {'ico_err': ico_err, 'lco_err': lco_err, 'sigma_mol_err': sigma_mol_err, 'mmol_err': mmol_err} 