"""
Main script for Upper Limit QA checks.

This focused QA tool checks that upper limit products are reasonable by:
1. Estimating RMS from line-free channels in the cube
2. Calculating expected Ico upper limit (90*rms for 30 km/s, 52*rms for 10 km/s)
3. Measuring minimum value in Ico_ul fits files
4. Comparing measured vs expected values

The output is a clean tabular report highlighting any outliers.
"""

import logging
import os
from datetime import datetime

from io_utils import load_config, load_summary_table
from upper_limit_qa import check_all_upper_limits, format_ul_result_for_table
from upper_limit_reporting import (
    create_summary_table,
    log_summary_table,
    write_summary_report,
    export_summary_csv
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)

CONFIG_PATH = "./config.md"


def main():
    """
    Main QA workflow for upper limit checks.
    """
    # Load configuration
    config = load_config(CONFIG_PATH)
    data_root = config['data_root']
    cube_root = config.get('cube_root', data_root)
    summary_table_path = config['summary_table']
    log_path = config['logging']['log_path']
    
    logging.info("="*80)
    logging.info("Starting Upper Limit QA Checks")
    logging.info("="*80)
    logging.info(f"Data root: {data_root}")
    logging.info(f"Cube root: {cube_root}")
    logging.info(f"Summary table: {summary_table_path}")
    
    # Load list of galaxies to check
    object_ids, summary_df = load_summary_table(summary_table_path)
    logging.info(f"Total galaxies to check: {len(object_ids)}")
    
    # Run upper limit checks for all galaxies
    all_results = {}
    table_results = []
    skipped = []
    
    for i, object_id in enumerate(object_ids, 1):
        logging.info(f"\n[{i}/{len(object_ids)}] Checking {object_id}...")
        
        try:
            # Check all 4 configurations (10/30 kms)
            ul_results = check_all_upper_limits(object_id, data_root, cube_root)
            
            # Check if any data exists
            has_data = any(r.get('exists', False) for r in ul_results.values())
            
            if not has_data:
                logging.warning(f"  No upper limit files found for {object_id}")
                skipped.append(object_id)
                continue
            
            # Store results
            all_results[object_id] = ul_results
            
            # Format for table
            table_row = format_ul_result_for_table(object_id, ul_results)
            table_results.append(table_row)
            
            # Quick status update
            flags = []
            for config, result in ul_results.items():
                if result.get('flag_ul_mismatch', False):
                    flags.append(config)
            
            if flags:
                logging.warning(f"  ⚠ Flags raised for: {', '.join(flags)}")
            else:
                logging.info(f"  ✓ All checks passed")
        
        except Exception as e:
            logging.error(f"  Error processing {object_id}: {e}")
            skipped.append(object_id)
    
    # Create summary table
    logging.info("\n" + "="*80)
    logging.info(f"Completed checks: {len(all_results)} galaxies")
    logging.info(f"Skipped: {len(skipped)} galaxies")
    if skipped:
        logging.info(f"Skipped galaxies: {', '.join(skipped)}")
    
    if not table_results:
        logging.warning("No results to display!")
        return
    
    # Create and display summary table
    df = create_summary_table(table_results)
    log_summary_table(df)
    
    # Extract product type from data_root path
    # e.g., "../../../products/v1.1/original/by_galaxy" -> "original"
    product_type = "unknown"
    try:
        path_parts = data_root.split('/')
        for i, part in enumerate(path_parts):
            if 'products' in part and i + 2 < len(path_parts):
                product_type = path_parts[i + 2]
                break
    except:
        pass
    
    # Write detailed report to file
    os.makedirs(log_path, exist_ok=True)
    dt = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Text report
    report_path = os.path.join(log_path, f"upper_limit_qa_10kms_30kms_{product_type}_{dt}.txt")
    write_summary_report(df, all_results, report_path)
    logging.info(f"\nDetailed report written to: {report_path}")
    
    # CSV export
    csv_path = os.path.join(log_path, f"upper_limit_qa_10kms_30kms_{product_type}_{dt}.csv")
    export_summary_csv(df, csv_path)
    
    logging.info("\n" + "="*80)
    logging.info("Upper Limit QA Complete!")
    logging.info("="*80)


if __name__ == "__main__":
    main()


