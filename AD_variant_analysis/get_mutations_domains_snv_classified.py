"""
Classify SNVs within protein domains.

This script intersects classified SNVs in CDS regions with protein domain 
regions to identify which SNVs fall within specific domains (DBD, AD, RD, etc.).
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
from tqdm import tqdm


def check_bedtools() -> bool:
    """Check if bedtools is available.
    
    Returns:
        True if bedtools is available, False otherwise
    """
    try:
        subprocess.run(['bedtools', '--version'], 
                      capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def load_protein_mapping(mapping_file: Path) -> Dict[str, str]:
    """
    Load UniProt ID to ENST mapping.
    
    Args:
        mapping_file: Path to CSV/TSV file with uniprotID and ENST columns
        
    Returns:
        Dictionary mapping UniProt IDs to ENST IDs
    """
    try:
        # Try reading as CSV first, then TSV
        try:
            df = pd.read_csv(mapping_file)
        except:
            df = pd.read_csv(mapping_file, sep='\t')
        
        # Validate required columns
        if 'uniprotID' not in df.columns or 'ENST' not in df.columns:
            raise ValueError("Mapping file must contain 'uniprotID' and 'ENST' columns")
        
        # Extract base ENST ID (remove version number)
        df = df[['uniprotID', 'ENST']].copy()
        df['ENST'] = df['ENST'].str.split('.').str[0]
        
        # Create mapping dictionary
        mapping = dict(zip(df['uniprotID'], df['ENST']))
        
        return mapping
        
    except Exception as e:
        print(f"Error loading mapping file: {e}", file=sys.stderr)
        sys.exit(1)


def run_bedtools_intersect(domain_file: Path, classified_snv_file: Path, 
                          output_file: Path) -> Tuple[bool, Optional[str]]:
    """
    Intersect domain BED file with classified SNVs using bedtools.
    
    Args:
        domain_file: Path to domain BED file
        classified_snv_file: Path to classified SNV BED file
        output_file: Path to output intersection file
        
    Returns:
        Tuple of (success: bool, error_message: str or None)
    """
    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        cmd = (
            f"bedtools intersect "
            f"-a {domain_file} "
            f"-b {classified_snv_file} "
            f"-wb > {output_file}"
        )
        
        result = subprocess.run(
            cmd,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            check=False
        )
        
        if result.returncode != 0 and result.stderr:
            stderr_text = result.stderr.decode().strip()
            if stderr_text:
                return False, stderr_text
        
        return True, None
        
    except Exception as e:
        return False, str(e)


def intersect_domains_with_classified_snvs(
    domain_dir: Path,
    classified_snv_dir: Path,
    output_dir: Path,
    mapping: Dict[str, str],
    domain_type: Optional[str] = None,
    max_workers: int = 8
) -> Tuple[int, int, list]:
    """
    Intersect all domain files with their corresponding classified SNV files.
    
    Args:
        domain_dir: Directory containing domain BED files (named by UniProt ID)
        classified_snv_dir: Directory containing classified SNV files (named by ENST)
        output_dir: Output directory for domain-SNV intersections
        mapping: Dictionary mapping UniProt IDs to ENST IDs
        domain_type: Specific domain subdirectory (e.g., 'DBD', 'AD') or None for all
        max_workers: Number of parallel workers
        
    Returns:
        Tuple of (processed count, total variants, errors list)
    """
    # Determine domain directory
    if domain_type and domain_type.lower() != 'none':
        domain_search_dir = domain_dir / domain_type
        output_subdir = output_dir / domain_type
    else:
        domain_search_dir = domain_dir
        output_subdir = output_dir
    
    if not domain_search_dir.exists():
        print(f"Error: Domain directory {domain_search_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    output_subdir.mkdir(parents=True, exist_ok=True)
    
    # Get all domain files
    domain_files = [f for f in domain_search_dir.iterdir() if f.is_file()]
    
    print(f"Found {len(domain_files)} domain files in {domain_search_dir}")
    
    # Process each domain file
    def process_domain_file(domain_file: Path):
        # Get UniProt ID from filename
        uniprot_id = domain_file.stem  # filename without extension
        
        # Look up ENST ID
        if uniprot_id not in mapping:
            return f"No ENST mapping for {uniprot_id}", 0
        
        enst_id = mapping[uniprot_id]
        
        # Find corresponding classified SNV file
        classified_snv_file = classified_snv_dir / f"{enst_id}.bed"
        
        if not classified_snv_file.exists():
            return f"No classified SNV file for {enst_id} (from {uniprot_id})", 0
        
        # Output file named by ENST
        output_file = output_subdir / f"{enst_id}.bed"
        
        # Run intersection
        success, error = run_bedtools_intersect(
            domain_file, classified_snv_file, output_file
        )
        
        if not success:
            return f"Error intersecting {uniprot_id}: {error}", 0
        
        # Count variants in output
        variant_count = 0
        if output_file.exists():
            with open(output_file) as f:
                variant_count = sum(1 for _ in f)
        
        return None, variant_count
    
    # Process in parallel
    errors = []
    processed = 0
    total_variants = 0
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_domain_file, f) for f in domain_files]
        
        for future in tqdm(as_completed(futures), total=len(futures), 
                          desc="Classifying domain SNVs"):
            error, variant_count = future.result()
            if error is None:
                processed += 1
                total_variants += variant_count
            else:
                errors.append(error)
    
    return processed, total_variants, errors


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Classify SNVs within protein domains",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This tool intersects classified SNVs (from classify-snvs) with protein domain 
regions to identify which variants fall within specific domains and their 
functional impact.

Domain files are named by UniProt ID (e.g., Q15583.bed)
Classified SNV files are named by ENST ID (e.g., ENST00000269305.bed)
A mapping file links UniProt IDs to ENST IDs

Examples:
  # Classify SNVs in all domains
  classify-domain-snvs \\
    -d domain_beds/ \\
    -s classified_snvs/ \\
    -o domain_snvs/ \\
    -m protein_mapping.csv

  # Classify SNVs only in DBD domains
  classify-domain-snvs \\
    -d domain_beds/ \\
    -s classified_snvs/ \\
    -o domain_snvs/ \\
    -m protein_mapping.csv \\
    -t DBD

  # Use more workers
  classify-domain-snvs \\
    -d domain_beds/ \\
    -s classified_snvs/ \\
    -o domain_snvs/ \\
    -m protein_mapping.csv \\
    -w 16

Input structure:
  domain_beds/              OR    domain_beds/
  ├── Q15583.bed                  ├── DBD/
  ├── P53_HUMAN.bed               │   ├── Q15583.bed
  └── ...                         │   └── ...
                                  ├── AD/
                                  │   ├── Q15583.bed
                                  │   └── ...
                                  └── ...

Output structure:
  domain_snvs/
  ├── ENST00000269305.bed
  ├── ENST00000357654.bed
  └── ...

  OR (with -t DBD):
  
  domain_snvs/
  └── DBD/
      ├── ENST00000269305.bed
      └── ...

Output format:
  Domain columns + Classified SNV columns (including Syn/Non-Syn/Nonsense)
        """
    )
    
    parser.add_argument(
        '-d', '--domain-dir',
        type=Path,
        required=True,
        help='Directory containing domain BED files (named by UniProt ID)'
    )
    
    parser.add_argument(
        '-s', '--classified-snv-dir',
        type=Path,
        required=True,
        help='Directory containing classified SNV BED files (named by ENST ID). '
             'Typically the output from classify-snvs.'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for domain-SNV intersections'
    )
    
    parser.add_argument(
        '-m', '--mapping',
        type=Path,
        required=True,
        help='CSV/TSV file mapping UniProt IDs to ENST IDs (columns: uniprotID, ENST)'
    )
    
    parser.add_argument(
        '-t', '--domain-type',
        type=str,
        default=None,
        help='Specific domain type subdirectory (e.g., DBD, AD, RD, Bif, IDR). '
             'If not specified, uses all files in domain-dir.'
    )
    
    parser.add_argument(
        '-w', '--workers',
        type=int,
        default=8,
        help='Number of parallel workers (default: 8)'
    )
    
    args = parser.parse_args()
    
    # Validate bedtools
    if not check_bedtools():
        print("Error: bedtools not found. Please install bedtools.", file=sys.stderr)
        print("\nInstallation:", file=sys.stderr)
        print("  Ubuntu/Debian: sudo apt-get install bedtools", file=sys.stderr)
        print("  macOS: brew install bedtools", file=sys.stderr)
        print("  Conda: conda install -c bioconda bedtools", file=sys.stderr)
        sys.exit(1)
    
    # Validate inputs
    if not args.domain_dir.exists():
        print(f"Error: Domain directory {args.domain_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not args.classified_snv_dir.exists():
        print(f"Error: Classified SNV directory {args.classified_snv_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not args.mapping.exists():
        print(f"Error: Mapping file {args.mapping} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Load mapping
    print(f"Loading UniProt to ENST mapping from {args.mapping}")
    mapping = load_protein_mapping(args.mapping)
    print(f"  Loaded {len(mapping)} protein mappings")
    
    # Run intersections
    print(f"\nIntersecting domain files with classified SNVs...")
    if args.domain_type:
        print(f"  Domain type: {args.domain_type}")
    print(f"  Domain directory: {args.domain_dir}")
    print(f"  Classified SNV directory: {args.classified_snv_dir}")
    print(f"  Output directory: {args.output_dir}")
    print()
    
    processed, total_variants, errors = intersect_domains_with_classified_snvs(
        args.domain_dir,
        args.classified_snv_dir,
        args.output_dir,
        mapping,
        args.domain_type,
        args.workers
    )
    
    # Summary
    print("\n" + "="*70)
    print("COMPLETE!")
    print("="*70)
    
    print(f"\nResults:")
    print(f"  Successfully processed: {processed}")
    print(f"  Total SNVs in domains: {total_variants:,}")
    print(f"  Errors/warnings: {len(errors)}")
    
    if errors:
        print(f"\nWarnings ({len(errors)}):")
        for error in errors[:10]:
            print(f"  {error}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more")
    
    # Count output files and classify by variant type
    output_files = list(args.output_dir.rglob('*.bed'))
    non_empty = [f for f in output_files if f.stat().st_size > 0]
    
    print(f"\nOutput:")
    print(f"  Total files: {len(output_files)}")
    print(f"  Files with domain SNVs: {len(non_empty)}")
    print(f"  Files without SNVs: {len(output_files) - len(non_empty)}")
    
    # Classify SNVs by type
    syn_count = 0
    nonsyn_count = 0
    nonsense_count = 0
    
    for output_file in non_empty:
        with open(output_file) as f:
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) > 0:
                    classification = fields[-1]
                    if classification == "Syn":
                        syn_count += 1
                    elif classification == "Non-Syn":
                        nonsyn_count += 1
                    elif classification == "Nonsense":
                        nonsense_count += 1
    
    total_classified = syn_count + nonsyn_count + nonsense_count
    if total_classified > 0:
        print(f"\nSNV Classification in domains:")
        print(f"  Synonymous: {syn_count:,} ({100*syn_count/total_classified:.1f}%)")
        print(f"  Non-synonymous: {nonsyn_count:,} ({100*nonsyn_count/total_classified:.1f}%)")
        print(f"  Nonsense: {nonsense_count:,} ({100*nonsense_count/total_classified:.1f}%)")
    
    if non_empty:
        print(f"\nExample files with domain SNVs:")
        for f in non_empty[:5]:
            variant_count = sum(1 for _ in open(f))
            rel_path = f.relative_to(args.output_dir)
            print(f"  {rel_path}: {variant_count} SNVs")
        if len(non_empty) > 5:
            print(f"  ... and {len(non_empty) - 5} more")
    
    print(f"\nOutput location: {args.output_dir}")


if __name__ == "__main__":
    main()