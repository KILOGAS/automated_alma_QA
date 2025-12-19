"""
Example script showing how to analyze upper limit QA results.

This demonstrates how to use the CSV output for further analysis.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def load_qa_results(csv_path):
    """Load QA results from CSV file."""
    df = pd.read_csv(csv_path)
    return df


def print_summary_statistics(df):
    """Print summary statistics for upper limit ratios."""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    for vel in ['10', '30']:
        ratio_col = f'{vel}_ratio'
        flag_col = f'{vel}_flag'
        
        if ratio_col in df.columns:
            ratios = df[ratio_col].dropna()
            n_flagged = df[flag_col].sum()
            
            print(f"\n{vel} km/s cubes:")
            print(f"  Total galaxies: {len(ratios)}")
            print(f"  Flagged: {n_flagged} ({100*n_flagged/len(ratios):.1f}%)")
            print(f"  Ratio statistics:")
            print(f"    Mean:   {ratios.mean():.3f}")
            print(f"    Median: {ratios.median():.3f}")
            print(f"    Std:    {ratios.std():.3f}")
            print(f"    Min:    {ratios.min():.3f}")
            print(f"    Max:    {ratios.max():.3f}")


def identify_outliers(df, threshold=3):
    """Identify outliers using standard deviation threshold."""
    print("\n" + "="*80)
    print(f"OUTLIERS (>{threshold} sigma from mean)")
    print("="*80)
    
    for vel in ['10', '30']:
        ratio_col = f'{vel}_ratio'
        
        if ratio_col in df.columns:
            ratios = df[ratio_col].dropna()
            mean = ratios.mean()
            std = ratios.std()
            
            # Find outliers
            outlier_mask = np.abs(df[ratio_col] - mean) > threshold * std
            outliers = df[outlier_mask]
            
            if len(outliers) > 0:
                print(f"\n{vel} km/s cubes ({len(outliers)} outliers):")
                for _, row in outliers.iterrows():
                    print(f"  {row['object_id']}: ratio={row[ratio_col]:.3f}")
            else:
                print(f"\n{vel} km/s cubes: No outliers")


def compare_velocity_widths(df):
    """Compare 10 km/s vs 30 km/s ratios."""
    print("\n" + "="*80)
    print("10 km/s vs 30 km/s COMPARISON")
    print("="*80)
    
    # Galaxies with both measurements
    both = df.dropna(subset=['10_ratio', '30_ratio'])
    
    if len(both) > 0:
        print(f"\nGalaxies with both measurements: {len(both)}")
        
        # Calculate correlation
        corr = both['10_ratio'].corr(both['30_ratio'])
        print(f"Correlation: {corr:.3f}")
        
        # Find galaxies with large discrepancies
        both['ratio_diff'] = np.abs(both['10_ratio'] - both['30_ratio'])
        discrepant = both.nlargest(5, 'ratio_diff')
        
        print("\nTop 5 galaxies with largest ratio differences:")
        for _, row in discrepant.iterrows():
            print(f"  {row['object_id']}: 10km={row['10_ratio']:.3f}, 30km={row['30_ratio']:.3f}, diff={row['ratio_diff']:.3f}")
    else:
        print("\nNo galaxies with both measurements")


def plot_ratio_distributions(df, output_path='ratio_distributions.png'):
    """Plot ratio distributions for both velocity widths."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for i, vel in enumerate(['10', '30']):
        ratio_col = f'{vel}_ratio'
        flag_col = f'{vel}_flag'
        
        if ratio_col in df.columns:
            # Get data
            ratios = df[ratio_col].dropna()
            flagged = df[df[flag_col] == True][ratio_col].dropna()
            passed = df[df[flag_col] == False][ratio_col].dropna()
            
            # Plot histogram
            ax = axes[i]
            ax.hist(passed, bins=30, alpha=0.7, label='Passed', color='green')
            ax.hist(flagged, bins=30, alpha=0.7, label='Flagged', color='red')
            
            # Add reference lines
            ax.axvline(0.8, color='orange', linestyle='--', label='Lower threshold')
            ax.axvline(1.2, color='orange', linestyle='--', label='Upper threshold')
            ax.axvline(1.0, color='black', linestyle='-', linewidth=2, label='Expected')
            
            ax.set_xlabel('Measured / Expected Ratio')
            ax.set_ylabel('Number of Galaxies')
            ax.set_title(f'{vel} km/s Cubes')
            ax.legend()
            ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"\nPlot saved to: {output_path}")


def plot_ratio_comparison(df, output_path='ratio_comparison.png'):
    """Plot 10 km/s vs 30 km/s ratios."""
    both = df.dropna(subset=['10_ratio', '30_ratio'])
    
    if len(both) == 0:
        print("\nNot enough data for comparison plot")
        return
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Separate by flag status
    both_flagged = both[(both['10_flag'] == True) | (both['30_flag'] == True)]
    both_passed = both[(both['10_flag'] == False) & (both['30_flag'] == False)]
    
    # Plot
    ax.scatter(both_passed['10_ratio'], both_passed['30_ratio'], 
               alpha=0.6, s=50, label='Passed', color='green')
    ax.scatter(both_flagged['10_ratio'], both_flagged['30_ratio'], 
               alpha=0.6, s=50, label='Flagged', color='red')
    
    # Add reference lines
    ax.axhline(1.0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    ax.axvline(1.0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    ax.plot([0, 5], [0, 5], 'k--', alpha=0.3, label='1:1 line')
    
    ax.set_xlabel('10 km/s Ratio (Measured / Expected)')
    ax.set_ylabel('30 km/s Ratio (Measured / Expected)')
    ax.set_title('Upper Limit Ratio Comparison')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Plot saved to: {output_path}")


def main():
    """Example analysis workflow."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python example_analysis.py <path_to_csv>")
        print("\nExample:")
        print("  python example_analysis.py ../logs/upper_limit_qa_20251217_220344.csv")
        return
    
    csv_path = sys.argv[1]
    
    print("="*80)
    print("UPPER LIMIT QA ANALYSIS")
    print("="*80)
    print(f"Loading: {csv_path}")
    
    # Load data
    df = load_qa_results(csv_path)
    print(f"Loaded {len(df)} galaxies")
    
    # Print summary statistics
    print_summary_statistics(df)
    
    # Identify outliers
    identify_outliers(df, threshold=3)
    
    # Compare velocity widths
    compare_velocity_widths(df)
    
    # Generate plots (requires matplotlib)
    try:
        plot_ratio_distributions(df)
        plot_ratio_comparison(df)
    except Exception as e:
        print(f"\nCould not generate plots: {e}")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)


if __name__ == "__main__":
    main()


