"""
Intersect protein domain regions with variants.

This script takes domain BED files (mapped to genomic coordinates) and 
intersects them with variant files to identify variants within specific 
protein domains (DBD, AD, RD, Bif, IDR).
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Optional, Tuple

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
        mapping_file: Path to TSV file with uniprotID and ENST columns
        
    Returns:
        Dictionary mapping UniProt IDs to ENST IDs
    """
    try:
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


def run_bedtools_intersect(domain_file: Path, cds_variants_file: Path, 
                          output_file: Path) -> Tuple[bool, Optional[str]]:
    """
    Intersect domain BED file with variants using bedtools.
    
    Args:
        domain_file: Path to domain BED file
        cds_variants_file: Path to CDS variants BED file
        output_file: Path to output intersection file
        
    Returns:
        Tuple of (success: bool, error_message: str or None)
    """
    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        cmd = (
            f"bedtools intersect "
            f"-a {domain_file} "
            f"-b {cds_variants_file} "
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


def intersect_domains_with_variants(
    domain_dir: Path,
    cds_variants_dir: Path,
    output_dir: Path,
    mapping: Dict[str, str],
    domain_type: Optional[str] = None,
    max_workers: int = 8
) -> Tuple[int, list]:
    """
    Intersect all domain files with their corresponding CDS variant files.
    
    Args:
        domain_dir: Directory containing domain BED files (named by UniProt ID)
        cds_variants_dir: Directory containing CDS variant intersection files (named by ENST)
        output_dir: Output directory for domain-variant intersections
        mapping: Dictionary mapping UniProt IDs to ENST IDs
        domain_type: Specific domain subdirectory (e.g., 'DBD', 'AD') or None for all
        max_workers: Number of parallel workers
        
    Returns:
        Tuple of (processed count, errors list)
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
            return f"No ENST mapping for {uniprot_id}"
        
        enst_id = mapping[uniprot_id]
        
        # Find corresponding CDS variants file
        cds_variants_file = cds_variants_dir / f"{enst_id}.bed"
        
        if not cds_variants_file.exists():
            return f"No CDS variants file for {enst_id} (from {uniprot_id})"
        
        # Output file named by ENST
        output_file = output_subdir / f"{enst_id}.bed"
        
        # Run intersection
        success, error = run_bedtools_intersect(domain_file, cds_variants_file, output_file)
        
        if not success:
            return f"Error intersecting {uniprot_id}: {error}"
        
        return None
    
    # Process in parallel
    errors = []
    processed = 0
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_domain_file, f) for f in domain_files]
        
        for future in tqdm(as_completed(futures), total=len(futures), 
                          desc="Intersecting domains"):
            result = future.result()
            if result is None:
                processed += 1
            else:
                errors.append(result)
    
    return processed, errors


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Intersect protein domain regions with genomic variants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This tool intersects domain BED files (genomic coordinates of protein domains)
with CDS variant BED files to identify which variants fall within specific domains.

Domain files are named by UniProt ID (e.g., Q15583.bed)
CDS variant files are named by ENST ID (e.g., ENST00000269305.bed)
A mapping file links UniProt IDs to ENST IDs

Examples:
  # Intersect all domains with CDS variants
  intersect-domain-variants \\
    -d domain_beds/ \\
    -c cds_variants/intersections/ \\
    -o domain_variants/ \\
    -m protein_mapping.txt

  # Intersect only DBD domains
  intersect-domain-variants \\
    -d domain_beds/ \\
    -c cds_variants/intersections/ \\
    -o domain_variants/ \\
    -m protein_mapping.txt \\
    -t DBD

  # Use more workers
  intersect-domain-variants \\
    -d domain_beds/ \\
    -c cds_variants/intersections/ \\
    -o domain_variants/ \\
    -m protein_mapping.txt \\
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
  domain_variants/
  ├── ENST00000269305.bed
  ├── ENST00000357654.bed
  └── ...

  OR (with -t DBD):
  
  domain_variants/
  └── DBD/
      ├── ENST00000269305.bed
      └── ...
        """
    )
    
    parser.add_argument(
        '-d', '--domain-dir',
        type=Path,
        required=True,
        help='Directory containing domain BED files (named by UniProt ID)'
    )
    
    parser.add_argument(
        '-c', '--cds-variants-dir',
        type=Path,
        required=True,
        help='Directory containing CDS variant BED files (named by ENST ID). '
             'Typically the output from intersect-variants.'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for domain-variant intersections'
    )
    
    parser.add_argument(
        '-m', '--mapping',
        type=Path,
        required=True,
        help='TSV file mapping UniProt IDs to ENST IDs (columns: uniprotID, ENST)'
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
    
    if not args.cds_variants_dir.exists():
        print(f"Error: CDS variants directory {args.cds_variants_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not args.mapping.exists():
        print(f"Error: Mapping file {args.mapping} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Load mapping
    print(f"Loading UniProt to ENST mapping from {args.mapping}")
    mapping = load_protein_mapping(args.mapping)
    print(f"  Loaded {len(mapping)} protein mappings")
    
    # Run intersections
    print(f"\nIntersecting domain files with CDS variants...")
    if args.domain_type:
        print(f"  Domain type: {args.domain_type}")
    print(f"  Domain directory: {args.domain_dir}")
    print(f"  CDS variants directory: {args.cds_variants_dir}")
    print(f"  Output directory: {args.output_dir}")
    print()
    
    processed, errors = intersect_domains_with_variants(
        args.domain_dir,
        args.cds_variants_dir,
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
    print(f"  Errors/warnings: {len(errors)}")
    
    if errors:
        print(f"\nWarnings ({len(errors)}):")
        for error in errors[:10]:
            print(f"  {error}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more")
    
    # Count output files with variants
    output_files = list(args.output_dir.rglob('*.bed'))
    non_empty = [f for f in output_files if f.stat().st_size > 0]
    
    print(f"\nOutput:")
    print(f"  Total files: {len(output_files)}")
    print(f"  Files with domain variants: {len(non_empty)}")
    print(f"  Files without variants: {len(output_files) - len(non_empty)}")
    
    if non_empty:
        print(f"\nExample files with domain variants:")
        for f in non_empty[:5]:
            variant_count = sum(1 for _ in open(f))
            rel_path = f.relative_to(args.output_dir)
            print(f"  {rel_path}: {variant_count} variants")
        if len(non_empty) > 5:
            print(f"  ... and {len(non_empty) - 5} more")
    
    print(f"\nOutput location: {args.output_dir}")


if __name__ == "__main__":
    main()