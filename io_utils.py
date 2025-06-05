import os
import pandas as pd
from astropy.io import fits
import yaml

def load_config(config_path):
    """Load YAML config from config.md (YAML frontmatter)."""
    with open(config_path, 'r') as f:
        # Skip lines until first '---'
        while True:
            line = f.readline()
            if not line:
                raise ValueError('No YAML frontmatter found in config file')
            if line.strip() == '---':
                break
        yaml_lines = []
        for line in f:
            if line.strip() == '---':
                break
            yaml_lines.append(line)
        config = yaml.safe_load(''.join(yaml_lines))
    return config

def load_summary_table(summary_path):
    """Load summary table and return list of object IDs to process (as strings)."""
    df = pd.read_csv(summary_path) if summary_path.endswith('.csv') else pd.read_excel(summary_path)
    # Use the first column as object_id if not specified
    if 'object_id' not in df.columns:
        df['object_id'] = df.iloc[:, 0].astype(str)
    return df['object_id'].tolist(), df

def find_data_files(config, object_id):
    """Find relevant FITS files for a given object ID using config patterns. Returns both unmaskedcube and maskedcube if present.
    Supports multiple patterns for a key (e.g., unmaskedcube as a list). Returns the first file found, or None if none exist.
    """
    base_dir = config['data_root']
    patterns = config['file_patterns']
    files = {}
    for key, pattern in patterns.items():
        # Support multiple patterns (list) for a key
        if isinstance(pattern, list):
            found = None
            for pat in pattern:
                rel_path = pat.format(object_id=object_id)
                candidate = os.path.join(base_dir, object_id, rel_path)
                if os.path.exists(candidate):
                    found = candidate
                    break
            files[key] = found
        else:
            rel_path = pattern.format(object_id=object_id)
            candidate = os.path.join(base_dir, object_id, rel_path)
            files[key] = candidate if os.path.exists(candidate) else None

    # If config has a cube_root, use that (and a subfolder {object_id}) for unmaskedcube, otherwise fall back to data_root.
    if "cube_root" in config:
        unmasked_cube_path = os.path.join(config["cube_root"], object_id, object_id + "_co2-1_10.0kmps_7m+12m.image.pbcor.ifumatched.fits")
    else:
        unmasked_cube_path = os.path.join(config["data_root"], object_id, object_id + "_co2-1_10.0kmps_7m+12m.image.pbcor.ifumatched.fits")
    files["unmaskedcube"] = unmasked_cube_path

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

def find_snr_map_path(config, object_id):
    """Return the expected SNR map path for a given object_id, or None if not found."""
    base_dir = config['data_root']
    snr_map = os.path.join(base_dir, object_id, 'moment_maps', f'{object_id}_mom0_SN.fits')
    return snr_map if os.path.exists(snr_map) else None 