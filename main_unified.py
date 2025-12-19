"""
Unified QA Script - Runs both Comprehensive QA and Upper Limit QA

This script generates two complementary reports:
1. Comprehensive QA - Detection, data quality, physical consistency
2. Upper Limit QA - Focused tabular validation of upper limits

Run this for a complete quality assessment of your data products.
"""

import logging
import os
from datetime import datetime

import numpy as np
from io_utils import (
    extract_yaml_from_md,
    find_data_files,
    find_error_map_paths,
    find_ico_10kms_30kms,
    load_config,
    load_summary_table,
)
from qa_checks import (
    check_all_positive,
    check_beam_and_units,
    check_cube_detection,
    check_edge_emission,
    check_lco_larger_than_ico,
    check_map_detection,
    check_mask_nonblank,
    check_scaling_factor_consistency,
    compare_ico_10kms_30kms,
    compare_mmol_to_lco,
    get_map_min_max_units,
    mask_nonblank_stats,
    measure_cube_max,
    measure_pixel_and_beam_size,
    measure_rms_cube_ends,
    velocity_range_nonblank,
)
from reporting_clean import (
    compute_qa_summary,
    get_failed_tests,
    log_clean_detailed_report,
    log_clean_summary,
    write_clean_report_to_file,
)
from upper_limit_qa import check_all_upper_limits, format_ul_result_for_table
from upper_limit_reporting import (
    create_summary_table,
    export_summary_csv,
    log_summary_table,
    write_summary_report,
)

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(message)s")

CONFIG_COMPREHENSIVE = "./config.md"
CONFIG_UPPER_LIMIT = "./config.md"


def run_comprehensive_qa(config, object_ids):
    """Run comprehensive QA checks."""
    logging.info("\n" + "╔" + "="*78 + "╗")
    logging.info("║" + " "*20 + "COMPREHENSIVE QA - STARTING" + " "*31 + "║")
    logging.info("╚" + "="*78 + "╝")
    
    data_root = config["data_root"]
    results = []
    skipped = []
    
    for i, object_id in enumerate(object_ids, 1):
        if i % 50 == 0 or i == 1:
            logging.info(f"\n[{i}/{len(object_ids)}] Processing comprehensive QA...")
        
        # Find files
        files = find_data_files(config, object_id)
        unmasked_cube_path = files["unmaskedcube"]
        masked_cube_path = files["maskedcube"]
        mask_path = files["mask"]
        ico_path = files["ico"]
        lco_path = files["lco"]
        sigma_mol_path = files["sigma_mol"]
        mmol_path = files["mmol"]
        err_files = find_error_map_paths(files)
        
        # Check for missing files
        required_files = [
            unmasked_cube_path,
            masked_cube_path,
            mask_path,
            ico_path,
            lco_path,
            sigma_mol_path,
        ]
        missing = [f for f in required_files if f is None or not os.path.exists(f)]
        if missing:
            skipped.append(object_id)
            continue
        
        # Detection checks
        unmasked_cube_detect = check_cube_detection(
            unmasked_cube_path, threshold_sigma=5, min_voxels=5, min_consecutive_channels=2
        )
        ico_detect = check_map_detection(ico_path)
        lco_detect = check_map_detection(lco_path)
        lco_gt_ico = check_lco_larger_than_ico(lco_path, ico_path)
        
        # Physical consistency
        scale_consistency = check_scaling_factor_consistency(
            ico_path,
            lco_path,
            err1_path=err_files.get("ico_err"),
            err2_path=err_files.get("lco_err"),
        )
        mmol_lco_result = compare_mmol_to_lco(mmol_path, lco_path)
        
        # Map statistics
        ico_minmax = get_map_min_max_units(ico_path)
        lco_minmax = get_map_min_max_units(lco_path)
        
        # QA checks
        edge_result = check_edge_emission(masked_cube_path)
        beam_units_result = check_beam_and_units(masked_cube_path)
        pix_beam_result = measure_pixel_and_beam_size(masked_cube_path)
        rms_result = measure_rms_cube_ends(masked_cube_path)
        max_result = measure_cube_max(masked_cube_path)
        vel_range_result = velocity_range_nonblank(masked_cube_path)
        mask_stats = mask_nonblank_stats(mask_path) if os.path.exists(mask_path) else {}
        
        # 10km/s vs 30km/s comparison
        ico_10kms_path, ico_30kms_path = find_ico_10kms_30kms(config, object_id)
        velocity_comparison = {}
        if ico_10kms_path and ico_30kms_path:
            velocity_comparison = compare_ico_10kms_30kms(ico_10kms_path, ico_30kms_path)
        
        # Aggregate results
        result = {
            "object_id": object_id,
            **unmasked_cube_detect,
            **ico_detect,
            **lco_detect,
            **lco_gt_ico,
            **scale_consistency,
            "ico_min": ico_minmax["min"],
            "ico_max": ico_minmax["max"],
            "ico_units": ico_minmax["units"],
            "lco_min": lco_minmax["min"],
            "lco_max": lco_minmax["max"],
            "lco_units": lco_minmax["units"],
            "edge_flag": edge_result["flag"],
            **edge_result["details"],
            **beam_units_result,
            **pix_beam_result,
            **rms_result,
            **max_result,
            **vel_range_result,
            **mask_stats,
            **mmol_lco_result,
            **velocity_comparison,
        }
        
        # Add flags
        result["flag_round_beam"] = not beam_units_result.get("round_beam", False)
        result["flag_kelvin_units"] = not beam_units_result.get("kelvin_units", False)
        result["flag_cube_detected"] = not unmasked_cube_detect.get("cube_detected", False)
        result["flag_ico_detected"] = not ico_detect.get("map_detected", False)
        result["flag_lco_detected"] = not lco_detect.get("map_detected", False)
        result["flag_lco_gt_ico"] = lco_gt_ico.get("flag_lco_gt_ico", False)
        result["flag_scaling_consistency"] = scale_consistency.get("flag_scaling_consistency", False)
        
        results.append(result)
    
    return results, skipped


def run_upper_limit_qa(config, object_ids):
    """Run upper limit QA checks."""
    logging.info("\n" + "╔" + "="*78 + "╗")
    logging.info("║" + " "*21 + "UPPER LIMIT QA - STARTING" + " "*32 + "║")
    logging.info("╚" + "="*78 + "╝")
    
    data_root = config["data_root"]
    cube_root = config.get("cube_root", data_root)
    
    all_results = {}
    table_results = []
    skipped = []
    
    for i, object_id in enumerate(object_ids, 1):
        if i % 50 == 0 or i == 1:
            logging.info(f"\n[{i}/{len(object_ids)}] Processing upper limit QA...")
        
        try:
            ul_results = check_all_upper_limits(object_id, data_root, cube_root)
            has_data = any(r.get("exists", False) for r in ul_results.values())
            
            if not has_data:
                skipped.append(object_id)
                continue
            
            all_results[object_id] = ul_results
            table_row = format_ul_result_for_table(object_id, ul_results)
            table_results.append(table_row)
        
        except Exception as e:
            logging.debug(f"  Error processing {object_id}: {e}")
            skipped.append(object_id)
    
    return all_results, table_results, skipped


def main():
    """Run both comprehensive and upper limit QA."""
    start_time = datetime.now()
    
    # Header
    logging.info("\n" + "╔" + "="*78 + "╗")
    logging.info("║" + " "*26 + "UNIFIED QA SYSTEM" + " "*35 + "║")
    logging.info("║" + " "*18 + "Comprehensive + Upper Limit QA" + " "*29 + "║")
    logging.info("╚" + "="*78 + "╝")
    
    # Load configurations
    config_comp = load_config(CONFIG_COMPREHENSIVE)
    config_ul = load_config(CONFIG_UPPER_LIMIT)
    
    # Load galaxy list
    summary_table_path = config_comp["summary_table"]
    object_ids, summary_df = load_summary_table(summary_table_path)
    
    logging.info(f"\n📊 Total galaxies to check: {len(object_ids)}")
    
    # ========================================================================
    # COMPREHENSIVE QA
    # ========================================================================
    comp_results, comp_skipped = run_comprehensive_qa(config_comp, object_ids)
    comp_summary = compute_qa_summary(object_ids, comp_results, comp_skipped)
    
    logging.info("\n" + "╔" + "="*78 + "╗")
    logging.info("║" + " "*20 + "COMPREHENSIVE QA - COMPLETE" + " "*31 + "║")
    logging.info("╚" + "="*78 + "╝")
    log_clean_summary(comp_summary)
    
    # ========================================================================
    # UPPER LIMIT QA
    # ========================================================================
    ul_all_results, ul_table_results, ul_skipped = run_upper_limit_qa(config_ul, object_ids)
    
    logging.info("\n" + "╔" + "="*78 + "╗")
    logging.info("║" + " "*21 + "UPPER LIMIT QA - COMPLETE" + " "*32 + "║")
    logging.info("╚" + "="*78 + "╝")
    
    if ul_table_results:
        ul_df = create_summary_table(ul_table_results)
        log_summary_table(ul_df)
    else:
        logging.warning("\n⚠ No upper limit results to display")
    
    # ========================================================================
    # SAVE REPORTS
    # ========================================================================
    log_path = config_comp["logging"]["log_path"]
    os.makedirs(log_path, exist_ok=True)
    dt = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Extract product type and velocity from data_root path
    # e.g., "../../../products/v1.1/original/by_galaxy" -> "original"
    product_type = "unknown"
    velocity = config_comp.get("data_subdir", "unknown")
    try:
        path_parts = config_comp["data_root"].split('/')
        for i, part in enumerate(path_parts):
            if 'products' in part and i + 2 < len(path_parts):
                product_type = path_parts[i + 2]
                break
    except:
        pass
    
    # Comprehensive QA report
    comp_report_path = os.path.join(log_path, f"comprehensive_qa_{velocity}_{product_type}_{dt}.txt")
    yaml_config = extract_yaml_from_md(CONFIG_COMPREHENSIVE)
    
    flagged_lines = []
    for r in comp_results:
        if r["object_id"] in comp_summary.flagged:
            failed_tests = get_failed_tests(r)
            flagged_lines.append(f"{r['object_id']}: {', '.join(failed_tests)}")
    
    write_clean_report_to_file(
        comp_summary, comp_results, flagged_lines, comp_report_path, yaml_config
    )
    
    # Upper limit QA reports
    if ul_table_results:
        ul_report_path = os.path.join(log_path, f"upper_limit_qa_10kms_30kms_{product_type}_{dt}.txt")
        ul_csv_path = os.path.join(log_path, f"upper_limit_qa_10kms_30kms_{product_type}_{dt}.csv")
        
        write_summary_report(ul_df, ul_all_results, ul_report_path)
        export_summary_csv(ul_df, ul_csv_path)
    
    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    elapsed = datetime.now() - start_time
    
    logging.info("\n" + "╔" + "="*78 + "╗")
    logging.info("║" + " "*30 + "COMPLETE!" + " "*37 + "║")
    logging.info("╚" + "="*78 + "╝")
    
    logging.info(f"\n⏱  Processing time: {elapsed}")
    logging.info(f"\n📁 Reports saved to: {log_path}/")
    logging.info(f"   • Comprehensive QA:  comprehensive_qa_{velocity}_{product_type}_{dt}.txt")
    if ul_table_results:
        logging.info(f"   • Upper Limit QA:    upper_limit_qa_10kms_30kms_{product_type}_{dt}.txt")
        logging.info(f"   • Upper Limit CSV:   upper_limit_qa_10kms_30kms_{product_type}_{dt}.csv")
    
    logging.info("\n" + "="*80)


if __name__ == "__main__":
    main()


