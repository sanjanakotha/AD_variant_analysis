"""
Domain Mapper - Convert protein domain coordinates to genomic BED format

This script maps protein domain coordinates to genomic coordinates using CDS (coding sequence)
information in BED format. It handles both positive and negative strand genes.
"""

import argparse
import glob
import math
import os
import sys
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def split_coords(data: List[int], strand: str) -> List[List[int]]:
    """
    Split a list of genomic coordinates into continuous segments.
    
    Based on: https://stackoverflow.com/questions/3149440/splitting-list-based-on-missing-numbers-in-a-sequence
    
    Args:
        data: List of genomic coordinates
        strand: Strand orientation ('+', '-', '1', '-1', 1, or -1)
        
    Returns:
        List of continuous coordinate segments
    """
    output = []
    if strand in ["-", "-1", -1]:
        data.reverse()
    for k, g in groupby(enumerate(data), lambda i_x: i_x[0] - i_x[1]):
        output.append(list(map(itemgetter(1), g)))
    return output


def parse_domain_coords(coord_string: str) -> List[Tuple[int, int]]:
    """
    Parse domain coordinate string into list of (start, end) tuples.
    
    Args:
        coord_string: Comma-separated coordinate ranges (e.g., "10-20,30-40")
        
    Returns:
        List of (start, end) coordinate tuples
    """
    domains = []
    if coord_string in ["", "NA-NA"]:
        return domains
    
    for coord in coord_string.split(","):
        if coord and "-" in coord:
            start, end = coord.split("-")
            domains.append((int(start), int(end)))
    
    return domains


def read_tf_data(input_file: Path) -> Dict:
    """
    Read transcription factor data from input file.
    
    Args:
        input_file: Path to TF data file
        
    Returns:
        Dictionary with TF information keyed by ENST ID
    """
    tf_data = {}
    
    with open(input_file) as f:
        f.readline()  # Skip header
        for line in f:
            fields = line.rstrip().split("\t")
            
            tf_id = fields[3]
            if tf_id == "Q15583-2":
                continue
            
            ensg = fields[4]
            enst = fields[5]
            
            tf_data[enst] = {
                'tf_id': tf_id,
                'ensg': ensg,
                'enst': enst,
                'length': len(fields[-2]) if len(fields) > 12 else 0,
                'domains': {
                    'DBD': parse_domain_coords(fields[6]),
                    'AD': parse_domain_coords(fields[7]),
                    'RD': parse_domain_coords(fields[8]),
                    'Bif': parse_domain_coords(fields[9]),
                    'IDR': parse_domain_coords(fields[11]) if len(fields) > 11 else []
                }
            }
    
    return tf_data


def read_cds_coords(cds_file: Path) -> Tuple[List[Tuple[int, int]], str, str]:
    """
    Read CDS coordinates from BED format file.
    
    Args:
        cds_file: Path to CDS BED file
        
    Returns:
        Tuple of (coordinates list, chromosome, strand)
    """
    coords = []
    chrom = None
    strand = None
    
    with open(cds_file) as f:
        for line in f:
            fields = line.rstrip().split("\t")
            chrom = fields[0]
            strand = fields[4]
            
            # Handle strand orientation
            if strand in ["-", "-1", -1]:
                start, end = int(fields[2]), int(fields[1]) + 1
            else:
                start, end = int(fields[1]) + 1, int(fields[2])
            
            coords.append((start, end))
    
    # Sort coordinates based on strand
    if strand in ["-", "-1", -1]:
        coords.sort(key=lambda x: -1 * x[0])
    else:
        coords.sort(key=lambda x: 1 * x[0])
    
    return coords, chrom, strand


def build_aa_to_nt_mapping(coords: List[Tuple[int, int]], strand: str) -> Tuple[Dict[int, List[int]], int]:
    """
    Build mapping from amino acid positions to nucleotide positions.
    
    Args:
        coords: List of (start, end) coordinate tuples
        strand: Strand orientation
        
    Returns:
        Tuple of (aa_to_nt dictionary, aa_length)
    """
    # Collect all nucleotide positions
    all_pos = []
    for start, end in coords:
        if strand in ["+", "1", 1]:
            all_pos += list(range(start, end + 1))
        elif strand in ["-", "-1", -1]:
            all_pos += list(range(start, end - 1, -1))
    
    # Calculate amino acid length
    nt_length = len(all_pos)
    aa_length = int(nt_length / 3)
    
    # Map each amino acid to its three nucleotides
    dict_aa_coords = {}
    start_codon = 0
    for i in range(aa_length):
        end_codon = start_codon + 3
        dict_aa_coords[i + 1] = all_pos[start_codon:end_codon]
        start_codon += 3
    
    return dict_aa_coords, aa_length


def map_domain_to_genomic(domain_coords: List[Tuple[int, int]], 
                          aa_to_nt: Dict[int, List[int]], 
                          strand: str,
                          chrom: str,
                          domain_name: str,
                          ensg: str,
                          enst: str,
                          tf_id: str) -> List[str]:
    """
    Map protein domain coordinates to genomic BED format entries.
    
    Args:
        domain_coords: List of (start, end) domain coordinates in AA
        aa_to_nt: Mapping from AA position to nucleotide positions
        strand: Strand orientation
        chrom: Chromosome name
        domain_name: Name of the domain
        ensg: Ensembl gene ID
        enst: Ensembl transcript ID
        tf_id: Transcription factor ID
        
    Returns:
        List of BED format strings
    """
    bed_lines = []
    
    for start_dom, end_dom in domain_coords:
        dom_pos = []
        
        # Collect all nucleotide positions for this domain
        for aa in range(start_dom, end_dom + 1):
            if aa <= len(aa_to_nt):
                dom_pos += aa_to_nt[aa]
            else:
                print(f"WARNING: {tf_id}", file=sys.stderr)
                print(f"Trying to add AA {aa} but protein has length {len(aa_to_nt)}", file=sys.stderr)
        
        # Split into continuous segments
        for cds_dom in split_coords(dom_pos, strand):
            start_bed = cds_dom[0] - 1
            end_bed = cds_dom[-1]
            
            bed_line = f"{chrom}\t{start_bed}\t{end_bed}\t{domain_name}\t{ensg}\t.\t{strand}\t{enst}"
            bed_lines.append(bed_line)
    
    return bed_lines


def process_transcripts(tf_data_file: Path, 
                       cds_directory: Path, 
                       output_directory: Path,
                       verbose: bool = False) -> None:
    """
    Process all transcripts and generate domain BED files.
    
    Args:
        tf_data_file: Path to TF data input file
        cds_directory: Directory containing CDS BED files
        output_directory: Directory for output BED files
        verbose: Print verbose output
    """
    # Create output directory if it doesn't exist
    output_directory.mkdir(parents=True, exist_ok=True)
    
    # Read TF data
    tf_data = read_tf_data(tf_data_file)
    
    # Get available CDS files
    cds_files = {Path(f).name: f for f in glob.glob(str(cds_directory / "*"))}
    
    # Process each transcript
    for enst, data in tf_data.items():
        if enst not in cds_files:
            print(f"No CDS for {enst}", file=sys.stderr)
            continue
        
        if verbose:
            print(f"Processing {data['tf_id']}")
        
        # Read CDS coordinates
        cds_file = cds_directory / enst
        coords, chrom, strand = read_cds_coords(cds_file)
        
        # Build AA to nucleotide mapping
        aa_to_nt, aa_length = build_aa_to_nt_mapping(coords, strand)
        
        # Process each domain type
        output_file = output_directory / data['tf_id']
        with open(output_file, 'w') as out:
            for domain_name, domain_coords in data['domains'].items():
                if not domain_coords:
                    continue
                
                bed_lines = map_domain_to_genomic(
                    domain_coords, aa_to_nt, strand, chrom, 
                    domain_name, data['ensg'], enst, data['tf_id']
                )
                
                for line in bed_lines:
                    out.write(line + "\n")


def main() -> None:
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Map protein domain coordinates to genomic BED format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python domain_mapper.py -i tf_data.txt -c cds_beds/ -o output_domains/
  python domain_mapper.py -i tf_data.txt -c cds_beds/ -o output_domains/ -v
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        type=Path,
        required=True,
        help='Input file containing TF data with domain coordinates'
    )
    
    parser.add_argument(
        '-c', '--cds-directory',
        type=Path,
        required=True,
        help='Directory containing CDS BED format files (named by ENST ID)'
    )
    
    parser.add_argument(
        '-o', '--output-directory',
        type=Path,
        required=True,
        help='Output directory for domain BED files'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.input.exists():
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not args.cds_directory.exists():
        print(f"Error: CDS directory {args.cds_directory} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Process transcripts
    print(f"Processing transcripts...")
    print(f"  Input file: {args.input}")
    print(f"  CDS directory: {args.cds_directory}")
    print(f"  Output directory: {args.output_directory}")
    print()
    
    process_transcripts(
        args.input,
        args.cds_directory,
        args.output_directory,
        args.verbose
    )
    
    # Summary
    print("\n" + "="*70)
    print("COMPLETE!")
    print("="*70)
    
    # Count output files
    output_files = list(args.output_directory.glob('*'))
    output_files = [f for f in output_files if f.is_file()]
    non_empty = [f for f in output_files if f.stat().st_size > 0]
    
    print(f"\nResults:")
    print(f"  Total domain files created: {len(output_files)}")
    print(f"  Files with domains: {len(non_empty)}")
    print(f"  Files without domains: {len(output_files) - len(non_empty)}")
    
    if non_empty:
        print(f"\nExample domain files:")
        for f in non_empty[:5]:
            line_count = sum(1 for _ in open(f))
            print(f"  {f.name}: {line_count} domain regions")
        if len(non_empty) > 5:
            print(f"  ... and {len(non_empty) - 5} more")
    
    print(f"\nOutput location: {args.output_directory}")


if __name__ == "__main__":
    main()