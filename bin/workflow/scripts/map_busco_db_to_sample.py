#!/usr/bin/env python3

"""
BUSCO Database Mapping Script
Automatically maps samples to their best BUSCO database based on taxonomy
Author: @LiaOb21
Date: 15/04/2026
Original script: https://github.com/LiaOb21/FSP_assembly_benchmarking/blob/main/workflow/scripts/map_busco_db_to_sample.py
"""

import pandas as pd
import argparse
import sys
from pathlib import Path

def normalise_taxon_name(taxon_name, extension):
    """Convert a taxonomic name to BUSCO database format"""
    if pd.isna(taxon_name) or taxon_name == "":
        return None
    return taxon_name.lower() + extension

def find_best_busco_db(sample_row, busco_lineages, extension, fallback_db):
    """Find the best BUSCO database for a sample based on taxonomy hierarchy"""
    hierarchy = ['Family', 'Order', 'Class', 'Phylum']
    
    for level in hierarchy:
        if level in sample_row and not pd.isna(sample_row[level]):
            potential_db = normalise_taxon_name(sample_row[level], extension)
            if potential_db in busco_lineages:
                return potential_db
    
    return fallback_db

def get_sample_busco_db(taxonomy_file, busco_lineages_file, output_file, extension, fallback_db, sample_name):
    """Get BUSCO database for a specific sample"""
    
    # Check if input files exist
    if not Path(taxonomy_file).exists():
        print(f"ERROR: Taxonomy file '{taxonomy_file}' not found!")
        sys.exit(1)
    
    if not Path(busco_lineages_file).exists():
        print(f"ERROR: BUSCO lineages file '{busco_lineages_file}' not found!")
        sys.exit(1)
    
    # Load data
    try:
        with open(busco_lineages_file, 'r') as f:
            busco_lineages = [line.strip() for line in f if line.strip()]
            
        taxonomy_df = pd.read_csv(taxonomy_file, sep='\t')
    
    except Exception as e:
        print(f"ERROR loading data: {e}")
        sys.exit(1)
    
    # Find the specific sample
    sample_row = taxonomy_df[taxonomy_df['Sample'] == sample_name]
    if sample_row.empty:
        print(f"ERROR: Sample '{sample_name}' not found in taxonomy file!")
        sys.exit(1)
    
    best_db = find_best_busco_db(sample_row.iloc[0], busco_lineages, extension, fallback_db)
    print(f"{sample_name}: {best_db}")
    
    # Save just the database name to output file
    with open(output_file, 'w') as f:
        f.write(best_db)
    
    return best_db

def main():
    parser = argparse.ArgumentParser(
        description="Get appropriate BUSCO database for a specific sample based on taxonomy"
    )
    
    parser.add_argument('-t', '--taxonomy', 
                       required=True,
                       help='Taxonomy file (CSV/TSV with Sample, Family, Order, Class, Phylum columns, tab separated)')
    
    parser.add_argument('-b', '--busco-lineages',
                       required=True, 
                       help='File containing available BUSCO lineages (one per line)')
    
    parser.add_argument('-o', '--output',
                       required=True,
                       help='Output file containing just the database name')
    
    parser.add_argument('-s', '--sample',
                       required=True,
                       help='Sample name to process')
    
    parser.add_argument('-f', '--fallback',
                       default='fungi_odb12',
                       help='Fallback database name, i.e. the highest taxonomic rank possible (default: fungi_odb12)')
    
    parser.add_argument('-e', '--extension',
                       default='_odb12',
                       help='BUSCO database extension (default: _odb12)')
    
    args = parser.parse_args()
    
    # Get the database for the specific sample
    get_sample_busco_db(
        args.taxonomy,
        args.busco_lineages, 
        args.output,
        args.extension,
        args.fallback,
        args.sample
    )

if __name__ == "__main__":
    main()