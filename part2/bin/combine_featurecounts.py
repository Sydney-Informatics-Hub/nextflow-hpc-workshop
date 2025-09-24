#!/usr/bin/env python3
import sys
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Combine featureCounts output files and summary files')
    parser.add_argument('--sample', required=True, type=str,
                        help='Sample name')
    parser.add_argument('--counts', nargs='+', required=True, 
                       help='List of featureCounts .txt output files')
    parser.add_argument('--summaries', nargs='+', required=True,
                       help='List of featureCounts .summary files')
    parser.add_argument('--output-counts', required=True,
                       help='Output file for combined counts')
    parser.add_argument('--output-summary', required=True,
                       help='Output file for combined summary')
    
    args = parser.parse_args()
    
    # Combine count files
    count_dfs = []
    for file in args.counts:
        df = pd.read_csv(file, sep='\t', comment='#')
        # Rename the last column to the sample name
        df.columns = list(df.columns[:-1]) + [args.sample]
        count_dfs.append(df)
    
    # Concatenate all count dataframes
    combined_counts = pd.concat(count_dfs, ignore_index=True)
    combined_counts.to_csv(args.output_counts, sep='\t', index=False)
    
    # Combine summary files
    summary_dfs = []
    for file in args.summaries:
        df = pd.read_csv(file, sep='\t', header=0, names=['Status', 'Count'])
        summary_dfs.append(df)
    
    # Sum values across all summary files
    combined_summary = summary_dfs[0].copy()
    for df in summary_dfs[1:]:
        combined_summary['Count'] += df['Count']

    combined_summary.columns = ['Status', args.sample]
    
    combined_summary.to_csv(args.output_summary, sep='\t', index=False, header=True)
    
    print(f"Combined {len(args.counts)} count files into {args.output_counts}")
    print(f"Combined {len(args.summaries)} summary files into {args.output_summary}")

if __name__ == "__main__":
    main()