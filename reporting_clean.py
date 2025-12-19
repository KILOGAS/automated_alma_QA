"""
Clean Reporting Module for Comprehensive QA

Provides visually improved, concise reporting that focuses on key information
while eliminating verbose data dumps.
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List


@dataclass
class QASummary:
    n_total: int
    n_checked: int
    n_passed: int
    n_flagged: int
    n_skipped: int
    skipped: List[str]
    flagged: List[str]


def compute_qa_summary(kgas_ids, results, skipped):
    """Compute summary statistics for QA results."""
    n_total = len(kgas_ids)
    n_skipped = len(skipped)
    n_checked = len(results)
    
    # Identify flagged objects
    flagged = []
    for r in results:
        is_flagged = (
            r.get("edge_flag", False)
            or r.get("flag_round_beam", False)
            or r.get("flag_kelvin_units", False)
            or r.get("flag_cube_detected", False)
            or r.get("flag_ico_detected", False)
            or r.get("flag_lco_detected", False)
            or r.get("flag_lco_gt_ico", False)
            or r.get("flag_scaling_consistency", False)
            or r.get("moment_map_error")
        )
        if is_flagged:
            flagged.append(r["object_id"])
    
    n_flagged = len(flagged)
    n_passed = n_checked - n_flagged
    
    return QASummary(
        n_total, n_checked, n_passed, n_flagged, n_skipped, skipped, flagged
    )


def get_failed_tests(result):
    """Extract list of failed tests for a result."""
    failed = []
    
    # Detection and physical consistency flags
    checks = [
        ("edge_flag", "Edge emission"),
        ("flag_round_beam", "Beam roundness"),
        ("flag_kelvin_units", "Kelvin units"),
        ("flag_lco_ico", "LCO/ICO scaling"),
        ("flag_mmol_lco", "Mmol/LCO scaling"),
        ("flag_cube_detected", "Cube detection"),
        ("flag_ico_detected", "ICO detection"),
        ("flag_lco_detected", "LCO detection"),
        ("flag_lco_gt_ico", "LCO > ICO"),
        ("flag_scaling_consistency", "Scaling consistency"),
    ]
    
    for key, label in checks:
        if result.get(key, False):
            failed.append(label)
    
    # Moment map error
    if result.get("moment_map_error"):
        failed.append("Moment map error")
    
    return failed


def log_clean_summary(summary: QASummary):
    """Log clean summary with visual formatting."""
    logging.info("\n" + "╔" + "="*78 + "╗")
    logging.info("║" + " "*30 + "QA SUMMARY" + " "*38 + "║")
    logging.info("╚" + "="*78 + "╝")
    
    logging.info(f"\n  📊 Total objects:     {summary.n_total}")
    logging.info(f"  ✓  Checked:           {summary.n_checked}")
    logging.info(f"  ✓  Passed:            {summary.n_passed}")
    logging.info(f"  ⚠  Flagged:           {summary.n_flagged}")
    logging.info(f"  ⊘  Skipped:           {summary.n_skipped}")
    
    if summary.skipped:
        logging.info(f"\n  Skipped objects: {', '.join(summary.skipped[:10])}")
        if len(summary.skipped) > 10:
            logging.info(f"    ... and {len(summary.skipped) - 10} more")
    
    if summary.flagged:
        logging.info(f"\n  ⚠ Flagged objects: {', '.join(summary.flagged)}")


def log_clean_detailed_report(results: List[Dict[str, Any]]):
    """Log detailed report with clean, visual formatting."""
    
    for r in results:
        # Header
        logging.info("\n" + "┌" + "─"*78 + "┐")
        logging.info("│ " + f"Object: {r['object_id']}".ljust(77) + "│")
        logging.info("└" + "─"*78 + "┘")
        
        # Overall status
        failed_tests = get_failed_tests(r)
        if failed_tests:
            logging.warning(f"\n  ⚠ Status: FAILED ({len(failed_tests)} issues)")
            logging.warning(f"  Issues: {', '.join(failed_tests)}")
        else:
            logging.info(f"\n  ✓ Status: PASSED")
        
        # Detection Summary
        logging.info("\n  ─── DETECTION ───")
        _log_detection_info(r)
        
        # Data Quality
        logging.info("\n  ─── DATA QUALITY ───")
        _log_data_quality(r)
        
        # Physical Consistency (if detected)
        if r.get('map_detected', False):
            logging.info("\n  ─── PHYSICAL CONSISTENCY ───")
            _log_physical_consistency(r)
        
        # 10km/s vs 30km/s Comparison (if available)
        if r.get('ico_10kms_integrated_intensity') is not None:
            logging.info("\n  ─── VELOCITY WIDTH COMPARISON ───")
            _log_velocity_comparison(r)


def _log_detection_info(r):
    """Log detection information."""
    # Cube detection
    if r.get('flag_cube_detected', False):
        logging.warning(f"  ✗ Cube detection:    FAIL (max={r.get('cube_max'):.4f}, rms={r.get('cube_rms'):.4f})")
    else:
        logging.info(f"  ✓ Cube detection:    PASS (max={r.get('cube_max'):.4f}, rms={r.get('cube_rms'):.4f})")
    
    # Map detections
    if r.get('flag_ico_detected', False):
        logging.warning(f"  ✗ ICO detection:     FAIL")
    else:
        logging.info(f"  ✓ ICO detection:     PASS (max={r.get('ico_max'):.4f} {r.get('ico_units', '')})")
    
    if r.get('flag_lco_detected', False):
        logging.warning(f"  ✗ LCO detection:     FAIL")
    else:
        logging.info(f"  ✓ LCO detection:     PASS (max={r.get('lco_max'):.4f} {r.get('lco_units', '')})")


def _log_data_quality(r):
    """Log data quality metrics."""
    # Edge emission
    edge_pixels = r.get('total_edge_pixels', 0)
    if r.get('edge_flag', False):
        logging.warning(f"  ✗ Edge emission:     {edge_pixels} pixels at edges")
    else:
        logging.info(f"  ✓ Edge emission:     Clean (0 pixels at edges)")
    
    # Beam
    if r.get('flag_round_beam', False):
        logging.warning(f"  ✗ Beam:              Not round (bmaj={r.get('bmaj'):.6f}, bmin={r.get('bmin'):.6f})")
    else:
        logging.info(f"  ✓ Beam:              Round (bmaj={r.get('bmaj'):.6f})")
    
    # RMS
    rms = r.get('rms_avg')
    if rms:
        logging.info(f"  • RMS:               {rms:.6f} K")
    
    # Velocity range
    vmin = r.get('vmin')
    vmax = r.get('vmax')
    if vmin is not None and vmax is not None:
        logging.info(f"  • Velocity range:    {vmin:.1f} to {vmax:.1f} {r.get('vunit', '')}")


def _log_physical_consistency(r):
    """Log physical consistency checks."""
    # LCO > ICO
    if r.get('flag_lco_gt_ico', False):
        logging.warning(f"  ✗ LCO > ICO:         FAIL")
    else:
        logging.info(f"  ✓ LCO > ICO:         PASS")
    
    # Scaling consistency
    scaling = r.get('scaling_factor')
    error_scaling = r.get('error_scaling_factor')
    if r.get('flag_scaling_consistency', False):
        logging.warning(f"  ✗ Scaling:           INCONSISTENT (map={scaling:.2f}, err={error_scaling:.2f})")
    elif scaling is not None and error_scaling is not None:
        logging.info(f"  ✓ Scaling:           Consistent ({scaling:.2f})")
    
    # Mmol/LCO ratio
    mmol_ratio = r.get('mmol_lco_ratio')
    if r.get('flag_mmol_lco', False) and mmol_ratio:
        logging.warning(f"  ✗ Mmol/LCO:          OFF (ratio={mmol_ratio:.3f})")
    elif mmol_ratio:
        logging.info(f"  ✓ Mmol/LCO:          OK (ratio={mmol_ratio:.3f})")


def _log_velocity_comparison(r):
    """Log 10km/s vs 30km/s comparison."""
    int_10 = r.get('ico_10kms_integrated_intensity')
    int_30 = r.get('ico_30kms_integrated_intensity')
    ratio = r.get('ico_10kms_30kms_integrated_intensity_ratio')
    
    if int_10 and int_30:
        logging.info(f"  • 10 km/s Ico:       {int_10:.2f} K km/s arcsec²")
        logging.info(f"  • 30 km/s Ico:       {int_30:.2f} K km/s arcsec²")
        if ratio:
            logging.info(f"  • Ratio (10/30):     {ratio:.3f}")


def write_clean_report_to_file(summary, results, flagged_details, output_path, config_text):
    """Write clean comprehensive report to file."""
    with open(output_path, 'w') as f:
        # Header
        f.write("="*80 + "\n")
        f.write(" "*25 + "COMPREHENSIVE QA REPORT\n")
        f.write("="*80 + "\n\n")
        
        # Configuration
        f.write("CONFIGURATION\n")
        f.write("-"*80 + "\n")
        f.write(config_text + "\n")
        f.write("-"*80 + "\n\n")
        
        # Summary
        f.write("SUMMARY\n")
        f.write("-"*80 + "\n")
        f.write(f"Total objects:   {summary.n_total}\n")
        f.write(f"Checked:         {summary.n_checked}\n")
        f.write(f"Passed:          {summary.n_passed}\n")
        f.write(f"Flagged:         {summary.n_flagged}\n")
        f.write(f"Skipped:         {summary.n_skipped}\n")
        
        if summary.skipped:
            f.write(f"\nSkipped objects: {', '.join(summary.skipped[:20])}\n")
            if len(summary.skipped) > 20:
                f.write(f"... and {len(summary.skipped) - 20} more\n")
        
        if summary.flagged:
            f.write(f"\nFlagged objects: {', '.join(summary.flagged)}\n")
        
        f.write("\n")
        
        # Flagged objects detail
        if flagged_details:
            f.write("\n" + "="*80 + "\n")
            f.write("FLAGGED OBJECTS\n")
            f.write("="*80 + "\n")
            for line in flagged_details:
                f.write(line + "\n")
        
        # Detailed results
        f.write("\n" + "="*80 + "\n")
        f.write("DETAILED RESULTS\n")
        f.write("="*80 + "\n")
        
        for r in results:
            f.write(f"\n{'─'*80}\n")
            f.write(f"Object: {r['object_id']}\n")
            f.write(f"{'─'*80}\n")
            
            # Status
            failed_tests = get_failed_tests(r)
            if failed_tests:
                f.write(f"\nStatus: FAILED ({len(failed_tests)} issues)\n")
                f.write(f"Issues: {', '.join(failed_tests)}\n")
            else:
                f.write(f"\nStatus: PASSED\n")
            
            # Key metrics
            f.write(f"\nDETECTION:\n")
            f.write(f"  Cube max:        {r.get('cube_max', 'N/A')}\n")
            f.write(f"  Cube RMS:        {r.get('cube_rms', 'N/A')}\n")
            f.write(f"  ICO max:         {r.get('ico_max', 'N/A')} {r.get('ico_units', '')}\n")
            f.write(f"  LCO max:         {r.get('lco_max', 'N/A')} {r.get('lco_units', '')}\n")
            
            f.write(f"\nDATA QUALITY:\n")
            f.write(f"  Edge pixels:     {r.get('total_edge_pixels', 0)}\n")
            f.write(f"  Beam (bmaj):     {r.get('bmaj', 'N/A')}\n")
            f.write(f"  RMS avg:         {r.get('rms_avg', 'N/A')}\n")
            
            if r.get('ico_10kms_integrated_intensity'):
                f.write(f"\nVELOCITY COMPARISON:\n")
                f.write(f"  10 km/s Ico:     {r.get('ico_10kms_integrated_intensity', 'N/A')}\n")
                f.write(f"  30 km/s Ico:     {r.get('ico_30kms_integrated_intensity', 'N/A')}\n")
                f.write(f"  Ratio (10/30):   {r.get('ico_10kms_30kms_integrated_intensity_ratio', 'N/A')}\n")


