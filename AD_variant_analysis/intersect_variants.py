"""
Intersect genomic variants with CDS BED files.

This script uses bedtools to find overlaps between genomic variants and 
coding sequences (CDS). It automatically sorts all inputs and stores 
sorted versions for reuse.
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional

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


def run_bedtools_sort(input_file: Path, output_file: Path) -> bool:
    """
    Sort a BED file using bedtools.
    
    Args:
        input_file: Path to input BED file
        output_file: Path to output sorted BED file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        cmd = f"bedtools sort -i {input_file} > {output_file}"
        subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error sorting {input_file}: {e.stderr.decode()}", file=sys.stderr)
        return False


def run_bedtools_intersect(cds_file: Path, variants_file: Path, output_file: Path) -> Optional[str]:
    """
    Intersect CDS BED file with variants using bedtools.
    
    Args:
        cds_file: Path to CDS BED file (must be sorted)
        variants_file: Path to variants BED file (must be sorted)
        output_file: Path to output intersection file
        
    Returns:
        None if successful, error message if failed
    """
    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        cmd = (
            f"bedtools intersect "
            f"-a {cds_file} "
            f"-b {variants_file} "
            f"-wb -sorted -nonamecheck "
            f"> {output_file}"
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
            if stderr_text and "Error" in stderr_text:
                return f"Error on {cds_file.name}: {stderr_text}"
        
        return None
        
    except Exception as e:
        return f"Exception on {cds_file.name}: {e}"


def sort_bed_files(bed_files: List[str], input_dir: Path, sorted_dir: Path, 
                   desc: str = "Sorting", max_workers: int = 8) -> int:
    """
    Sort multiple BED files in parallel.
    
    Args:
        bed_files: List of BED file names
        input_dir: Directory containing input BED files
        sorted_dir: Directory for sorted BED files
        desc: Description for progress bar
        max_workers: Number of parallel workers
        
    Returns:
        Number of successfully sorted files
    """
    sorted_dir.mkdir(parents=True, exist_ok=True)
    
    def sort_file(filename):
        input_file = input_dir / filename
        # Ensure .bed extension
        if not filename.endswith('.bed'):
            output_name = f"{filename}.bed"
        else:
            output_name = filename
        output_file = sorted_dir / output_name
        
        # Skip if already sorted and newer than input
        if output_file.exists():
            if output_file.stat().st_mtime > input_file.stat().st_mtime:
                return True  # Already sorted and up-to-date
        
        return run_bedtools_sort(input_file, output_file)
    
    success_count = 0
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(sort_file, f) for f in bed_files]
        for future in tqdm(as_completed(futures), total=len(futures), desc=desc):
            if future.result():
                success_count += 1
    
    return success_count


def intersect_variants(cds_sorted_dir: Path, variants_sorted_file: Path, 
                      output_dir: Path, cds_files: Optional[List[str]] = None, 
                      max_workers: int = 8) -> None:
    """
    Intersect sorted CDS BED files with sorted variants file.
    
    Args:
        cds_sorted_dir: Directory containing sorted CDS BED files
        variants_sorted_file: Path to sorted variants BED file
        output_dir: Directory for output intersection files
        cds_files: List of specific CDS files to process (None = all)
        max_workers: Number of parallel workers
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get list of sorted CDS files
    if cds_files is None:
        cds_files = [f.name for f in cds_sorted_dir.iterdir() if f.is_file()]
    else:
        # Ensure .bed extension
        cds_files = [f"{f}.bed" if not f.endswith('.bed') else f for f in cds_files]
    
    def process_file(filename):
        cds_file = cds_sorted_dir / filename
        
        # Skip if not a file
        if not cds_file.is_file():
            return f"Skipping non-file: {filename}"
        
        # Determine output name
        name = filename.replace('.bed', '')
        output_file = output_dir / f"{name}.bed"
        
        return run_bedtools_intersect(cds_file, variants_sorted_file, output_file)
    
    # Run intersections in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_file, f) for f in cds_files]
        for future in tqdm(as_completed(futures), total=len(futures), 
                          desc="Intersecting variants"):
            result = future.result()
            if result:
                print(result, file=sys.stderr)


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Intersect genomic variants with CDS BED files using bedtools. "
                    "Automatically sorts all inputs and stores sorted files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - sorts everything automatically
  intersect-variants -c cds_beds/ -v variants.bed -o output/

  # Process specific CDS files only
  intersect-variants -c cds_beds/ -v variants.bed -o output/ -f genes.txt

  # Use more parallel workers
  intersect-variants -c cds_beds/ -v variants.bed -o output/ -w 16

Output structure:
  output/
  ├── intersections/          # Final intersection results
  │   ├── ENST00000001.bed
  │   └── ENST00000002.bed
  └── sorted/                 # Sorted input files (reusable)
      ├── cds/                # Sorted CDS files
      │   ├── ENST00000001.bed
      │   └── ENST00000002.bed
      └── variants/           # Sorted variants file
          └── variants.bed

Requirements:
  - bedtools must be installed and in PATH
  - Input files must be in BED format
        """
    )
    
    parser.add_argument(
        '-c', '--cds-dir',
        type=Path,
        required=True,
        help='Directory containing CDS BED files'
    )
    
    parser.add_argument(
        '-v', '--variants',
        type=Path,
        required=True,
        help='Variants BED file'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for all results (intersections and sorted files)'
    )
    
    parser.add_argument(
        '-f', '--files-list',
        type=Path,
        help='Text file with list of CDS files to process (one per line). '
             'If not provided, processes all files in cds-dir.'
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
    if not args.cds_dir.exists():
        print(f"Error: CDS directory {args.cds_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not args.cds_dir.is_dir():
        print(f"Error: {args.cds_dir} is not a directory", file=sys.stderr)
        sys.exit(1)
    
    if not args.variants.exists():
        print(f"Error: Variants file {args.variants} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Read file list if provided
    cds_files_to_process = None
    if args.files_list:
        if not args.files_list.exists():
            print(f"Error: Files list {args.files_list} does not exist", file=sys.stderr)
            sys.exit(1)
        with open(args.files_list) as f:
            cds_files_to_process = [line.strip() for line in f if line.strip()]
        print(f"Will process {len(cds_files_to_process)} CDS files from {args.files_list}")
    else:
        all_cds = [f.name for f in args.cds_dir.iterdir() if f.is_file()]
        print(f"Will process all {len(all_cds)} CDS files from {args.cds_dir}")
    
    # Setup output directories
    sorted_dir = args.output_dir / 'sorted'
    cds_sorted_dir = sorted_dir / 'cds'
    variants_sorted_dir = sorted_dir / 'variants'
    intersections_dir = args.output_dir / 'intersections'
    
    print(f"\nOutput structure:")
    print(f"  Sorted CDS files: {cds_sorted_dir}")
    print(f"  Sorted variants:  {variants_sorted_dir}")
    print(f"  Intersections:    {intersections_dir}")
    print()
    
    # Step 1: Sort CDS files
    print("Step 1: Sorting CDS files...")
    if cds_files_to_process is None:
        cds_files_to_sort = [f.name for f in args.cds_dir.iterdir() if f.is_file()]
    else:
        cds_files_to_sort = cds_files_to_process
    
    sorted_count = sort_bed_files(
        cds_files_to_sort,
        args.cds_dir,
        cds_sorted_dir,
        desc="Sorting CDS files",
        max_workers=args.workers
    )
    print(f"  Sorted {sorted_count}/{len(cds_files_to_sort)} CDS files")
    
    # Step 2: Sort variants
    print("\nStep 2: Sorting variants file...")
    variants_sorted_file = variants_sorted_dir / args.variants.name
    
    # Check if already sorted and up-to-date
    if variants_sorted_file.exists():
        if variants_sorted_file.stat().st_mtime > args.variants.stat().st_mtime:
            print(f"  Variants already sorted (using cached version)")
        else:
            print(f"  Re-sorting variants (source file updated)")
            run_bedtools_sort(args.variants, variants_sorted_file)
    else:
        run_bedtools_sort(args.variants, variants_sorted_file)
        print(f"  Sorted variants file")
    
    # Step 3: Intersect
    print("\nStep 3: Intersecting CDS files with variants...")
    print(f"  Using sorted CDS files from: {cds_sorted_dir}")
    print(f"  Using sorted variants from: {variants_sorted_file}")
    print(f"  Writing intersections to: {intersections_dir}")
    print()
    
    intersect_variants(
        cds_sorted_dir,
        variants_sorted_file,
        intersections_dir,
        cds_files_to_process,
        max_workers=args.workers
    )
    
    # Summary
    print("\n" + "="*70)
    print("COMPLETE!")
    print("="*70)
    
    output_files = list(intersections_dir.glob('*.bed'))
    non_empty = [f for f in output_files if f.stat().st_size > 0]
    
    print(f"\nResults:")
    print(f"  Total intersection files: {len(output_files)}")
    print(f"  Files with variants:      {len(non_empty)}")
    print(f"  Files without variants:   {len(output_files) - len(non_empty)}")
    print(f"\nOutput locations:")
    print(f"  Intersections: {intersections_dir}")
    print(f"  Sorted files:  {sorted_dir} (reusable for future runs)")
    
    if non_empty:
        print(f"\nExample files with variants:")
        for f in non_empty[:5]:
            variant_count = sum(1 for _ in open(f))
            print(f"  {f.name}: {variant_count} variants")
        if len(non_empty) > 5:
            print(f"  ... and {len(non_empty) - 5} more")


if __name__ == "__main__":
    main()