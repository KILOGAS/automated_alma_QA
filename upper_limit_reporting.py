"""
Upper Limit QA Reporting Module

Provides clean, tabular reporting for upper limit measurements.
"""

import pandas as pd
import logging


def create_summary_table(results_list):
    """
    Create a pandas DataFrame with upper limit results for all galaxies.
    
    Args:
        results_list: List of dicts with upper limit results for each galaxy
    
    Returns:
        pandas DataFrame with tabular results
    """
    df = pd.DataFrame(results_list)
    
    # Reorder columns for readability
    col_order = ['object_id']
    
    # Add columns for each velocity width
    for vel in ['10', '30']:
        col_order.extend([
            f'{vel}_ul_min',
            f'{vel}_ul_bmaj',
            f'{vel}_ul_bmin',
            f'{vel}_ul_nonblank_pixels',
            f'{vel}_cube_rms',
            f'{vel}_cube_bmaj',
            f'{vel}_cube_bmin',
        ])
    
    # Only include columns that exist
    col_order = [c for c in col_order if c in df.columns]
    df = df[col_order]
    
    return df


def format_detailed_report(object_id, ul_results):
    """
    Format detailed report for a single galaxy.
    
    Args:
        object_id: Galaxy ID
        ul_results: Dict with upper limit results for all configurations
    
    Returns:
        List of strings for report
    """
    lines = []
    lines.append(f"\n{'='*60}")
    lines.append(f"Object: {object_id}")
    lines.append('='*60)
    
    for config, result in sorted(ul_results.items()):
        lines.append(f"\n  {config}:")
        lines.append("  " + "-"*56)
        
        if not result.get('exists', False):
            lines.append(f"    Status: MISSING")
            lines.append(f"    Error: {result.get('error', 'Unknown')}")
            continue
        
        if result.get('error'):
            lines.append(f"    Status: ERROR")
            lines.append(f"    Error: {result['error']}")
            continue
        
        # From Ico_ul.fits
        lines.append("    FROM ICO_UL.FITS:")
        ul_min = result.get('ul_min')
        ul_bmaj = result.get('ul_bmaj')
        ul_bmin = result.get('ul_bmin')
        ul_nonblank = result.get('ul_nonblank_pixels')
        
        if ul_min is not None:
            lines.append(f"      Minimum value:       {ul_min:.6f}")
        else:
            lines.append(f"      Minimum value:       N/A")
        
        if ul_bmaj is not None:
            lines.append(f"      BMAJ:                {ul_bmaj:.8f}")
        else:
            lines.append(f"      BMAJ:                N/A")
        
        if ul_bmin is not None:
            lines.append(f"      BMIN:                {ul_bmin:.8f}")
        else:
            lines.append(f"      BMIN:                N/A")
        
        if ul_nonblank is not None:
            lines.append(f"      Non-blank pixels:    {ul_nonblank}")
        else:
            lines.append(f"      Non-blank pixels:    N/A")
        
        # From cube
        lines.append("\n    FROM CUBE (*IMAGE.FITS):")
        cube_rms = result.get('cube_rms')
        cube_bmaj = result.get('cube_bmaj')
        cube_bmin = result.get('cube_bmin')
        vel_range = result.get('velocity_range_kms', 300)
        n_chan = result.get('n_channels_used', 'N/A')
        
        if cube_rms is not None:
            lines.append(f"      RMS (lowest {vel_range} km/s): {cube_rms:.8f} K")
            lines.append(f"      Channels used:       {n_chan}")
        else:
            lines.append(f"      RMS:                 N/A")
        
        if cube_bmaj is not None:
            lines.append(f"      BMAJ:                {cube_bmaj:.8f}")
        else:
            lines.append(f"      BMAJ:                N/A")
        
        if cube_bmin is not None:
            lines.append(f"      BMIN:                {cube_bmin:.8f}")
        else:
            lines.append(f"      BMIN:                N/A")
    
    return lines


def log_summary_table(df):
    """
    Log summary table in a readable format.
    
    Args:
        df: pandas DataFrame with upper limit results
    """
    logging.info("\n" + "="*80)
    logging.info("UPPER LIMIT MEASUREMENTS SUMMARY")
    logging.info("="*80)
    
    # Format the dataframe for better display
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', 20)
    
    # Create a formatted version with better column names
    display_df = df.copy()
    
    # Rename columns for clarity
    rename_dict = {}
    for col in display_df.columns:
        if col == 'object_id':
            rename_dict[col] = 'Galaxy'
        elif '_ul_min' in col:
            vel = col.split('_')[0]
            rename_dict[col] = f'{vel}_UL_min'
        elif '_ul_bmaj' in col:
            vel = col.split('_')[0]
            rename_dict[col] = f'{vel}_UL_bmaj'
        elif '_ul_bmin' in col:
            vel = col.split('_')[0]
            rename_dict[col] = f'{vel}_UL_bmin'
        elif '_ul_nonblank' in col:
            vel = col.split('_')[0]
            rename_dict[col] = f'{vel}_UL_npix'
        elif '_cube_rms' in col:
            vel = col.split('_')[0]
            rename_dict[col] = f'{vel}_Cube_RMS'
        elif '_cube_bmaj' in col:
            vel = col.split('_')[0]
            rename_dict[col] = f'{vel}_Cube_bmaj'
        elif '_cube_bmin' in col:
            vel = col.split('_')[0]
            rename_dict[col] = f'{vel}_Cube_bmin'
    
    display_df = display_df.rename(columns=rename_dict)
    
    # Format numeric columns
    for col in display_df.columns:
        if 'RMS' in col or 'min' in col:
            display_df[col] = display_df[col].apply(lambda x: f'{x:.6f}' if pd.notna(x) else 'N/A')
        elif 'bmaj' in col or 'bmin' in col:
            display_df[col] = display_df[col].apply(lambda x: f'{x:.8f}' if pd.notna(x) else 'N/A')
        elif 'npix' in col:
            display_df[col] = display_df[col].apply(lambda x: f'{int(x)}' if pd.notna(x) else 'N/A')
    
    logging.info("\n" + display_df.to_string(index=False))
    logging.info("\n" + "="*80)


def write_summary_report(df, all_results, output_path):
    """
    Write comprehensive report to file.
    
    Args:
        df: pandas DataFrame with summary results
        all_results: Dict mapping object_id to detailed results
        output_path: Path to output file
    """
    with open(output_path, 'w') as f:
        # Header
        f.write("="*80 + "\n")
        f.write("UPPER LIMIT MEASUREMENTS REPORT\n")
        f.write("="*80 + "\n\n")
        
        # Summary statistics
        n_total = len(df)
        
        f.write(f"Total galaxies checked: {n_total}\n\n")
        
        # Summary table
        f.write("SUMMARY TABLE\n")
        f.write("-"*80 + "\n")
        f.write(df.to_string(index=False))
        f.write("\n\n")
        
        # Detailed results for each galaxy
        f.write("\n" + "="*80 + "\n")
        f.write("DETAILED RESULTS\n")
        f.write("="*80 + "\n")
        
        for object_id, ul_results in sorted(all_results.items()):
            report_lines = format_detailed_report(object_id, ul_results)
            for line in report_lines:
                f.write(line + "\n")


def export_summary_csv(df, output_path):
    """
    Export summary table to CSV for easy analysis in Excel/Python.
    
    Args:
        df: pandas DataFrame with summary results
        output_path: Path to output CSV file
    """
    df.to_csv(output_path, index=False)
    logging.info(f"CSV exported to: {output_path}")
