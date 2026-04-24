"""
Classify Single Nucleotide Variants (SNVs) as synonymous, non-synonymous, or nonsense.

This script analyzes SNVs in CDS regions and determines their effect on protein
translation by comparing wild-type and mutant codons.
"""

import argparse
import os
import pickle
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.Seq import Seq
from tqdm import tqdm


# Nucleotide complement mapping for reverse strand
COMPLEMENT = {"C": "G", "A": "T", "T": "A", "G": "C"}


class SNVClassifier:
    """Classify SNVs based on their effect on protein translation."""
    
    def __init__(self, proteins_file: Path, dna_transcripts_file: Path, 
                 sorted_cds_dir: Path):
        """
        Initialize classifier with reference sequences.
        
        Args:
            proteins_file: Pickle file with protein sequences (ENST -> amino acid sequence)
            dna_transcripts_file: Pickle file with DNA transcripts (ENST -> nucleotide sequence)
            sorted_cds_dir: Directory containing sorted CDS BED files
        """
        self.proteins = pickle.load(open(proteins_file, 'rb'))
        self.dna_transcripts = pickle.load(open(dna_transcripts_file, 'rb'))
        self.sorted_cds_dir = sorted_cds_dir
    
    def load_cds_coords(self, enst: str, strand: str) -> List[Tuple[int, int]]:
        """
        Load CDS coordinates from sorted BED file.
        
        Args:
            enst: ENST transcript ID
            strand: Strand orientation (+, -, 1, -1)
            
        Returns:
            List of (start, end) coordinate tuples
        """
        cds_file = self.sorted_cds_dir / f"{enst}.bed"
        
        if not cds_file.exists():
            raise FileNotFoundError(f"CDS file not found: {cds_file}")
        
        coords = []
        with open(cds_file) as f:
            for line in f:
                fields = line.rstrip().split('\t')
                start = int(fields[1]) + 1  # Convert BED (0-based) to 1-based
                end = int(fields[2])
                coords.append((start, end))
        
        # Reverse if on negative strand
        if strand in ["-", "-1", -1]:
            coords.reverse()
        
        return coords
    
    def build_genomic_to_cds_mapping(self, coords: List[Tuple[int, int]], strand: str, 
                                     nt_seq: str) -> Tuple[Dict[int, str], Dict[int, List[int]]]:
        """
        Build mapping from genomic coordinates to CDS positions and codons.
        
        Args:
            coords: List of (start, end) exon coordinates
            strand: Strand orientation
            nt_seq: Nucleotide sequence of CDS
            
        Returns:
            Tuple of (pos_to_nt dict, aa_to_codon dict)
        """
        # Build list of all genomic positions in CDS order
        total_coords = []
        for start, end in coords:
            if strand in ["-", "-1", -1]:
                total_coords += list(range(end, start - 1, -1))
            else:
                total_coords += list(range(start, end + 1))
        
        # Map genomic position to nucleotide
        pos_to_nt = {}
        for i, coord in enumerate(total_coords):
            if i >= len(nt_seq):
                raise ValueError(
                    f"Coordinate mismatch: position {i} exceeds sequence length {len(nt_seq)}"
                )
            pos_to_nt[coord] = nt_seq[i]
        
        # Map amino acid position to codon coordinates
        aa_to_codon = {}
        for i in range(len(total_coords) // 3):
            start_idx = i * 3
            aa_to_codon[i + 1] = total_coords[start_idx:start_idx + 3]
        
        return pos_to_nt, aa_to_codon
    
    def classify_variant(self, mut_pos: int, mut_nt: str, strand: str,
                        enst: str) -> Optional[Tuple[str, str, str]]:
        """
        Classify a single variant.
        
        Args:
            mut_pos: Genomic position of mutation
            mut_nt: Mutant nucleotide
            strand: Strand orientation
            enst: ENST transcript ID
            
        Returns:
            Tuple of (wt_aa, mt_aa, classification) or None if error
        """
        # Validate input
        if len(mut_nt) != 1 or mut_nt == ".":
            return None
        
        # Reverse complement if on negative strand
        if strand in ["-", "-1", -1]:
            mut_nt = COMPLEMENT[mut_nt]
        
        # Get reference sequences
        if enst not in self.dna_transcripts:
            raise ValueError(f"No DNA transcript for {enst}")
        if enst not in self.proteins:
            raise ValueError(f"No protein sequence for {enst}")
        
        nt_seq = self.dna_transcripts[enst]
        protein_seq = self.proteins[enst]
        
        # Load CDS coordinates
        coords = self.load_cds_coords(enst, strand)
        
        # Build mappings
        pos_to_nt, aa_to_codon = self.build_genomic_to_cds_mapping(
            coords, strand, nt_seq
        )
        
        # Find which codon contains the mutation
        codon_coords = None
        aa_pos = None
        for aa_position, codon in aa_to_codon.items():
            if mut_pos in codon:
                aa_pos = aa_position
                codon_coords = codon
                break
        
        if codon_coords is None:
            raise ValueError(f"Mutation position {mut_pos} not found in any codon")
        
        # Build wild-type codon
        wt_codon = {pos: pos_to_nt[pos] for pos in codon_coords}
        
        # Build mutant codon
        mt_codon = dict(wt_codon)
        mt_codon[mut_pos] = mut_nt
        
        # Convert to strings
        wt_codon_str = ''.join(wt_codon[pos] for pos in codon_coords)
        mt_codon_str = ''.join(mt_codon[pos] for pos in codon_coords)
        
        # Translate
        wt_aa = str(Seq(wt_codon_str).translate())
        mt_aa = str(Seq(mt_codon_str).translate())
        
        # Classify
        if wt_aa == mt_aa:
            classification = "Syn"
        elif mt_aa == "*":
            classification = "Nonsense"
        else:
            classification = "Non-Syn"
        
        return wt_aa, mt_aa, classification
    
    def process_variant_file(self, input_file: Path, output_file: Path) -> Tuple[int, int, str]:
        """
        Process a single variant file and classify all SNVs.
        
        Args:
            input_file: Input BED file with variants
            output_file: Output file for classified variants
            
        Returns:
            Tuple of (total_variants, classified_variants, error_message)
        """
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        total = 0
        classified = 0
        
        try:
            with open(input_file) as f_in, open(output_file, 'w') as f_out:
                for line in f_in:
                    total += 1
                    fields = line.rstrip().split('\t')
                    
                    # Extract variant information
                    try:
                        enst = fields[3]
                        strand = fields[4]
                        mut_pos = int(fields[-4])
                        ref_nt = fields[-3]
                        mut_nt = fields[-2]
                        
                        # Validate SNV
                        if len(mut_nt) != 1 or len(ref_nt) != 1 or mut_nt == ".":
                            continue
                        
                        # Classify variant
                        result = self.classify_variant(mut_pos, mut_nt, strand, enst)
                        
                        if result is None:
                            continue
                        
                        wt_aa, mt_aa, classification = result
                        
                        # Write output
                        fields.extend([wt_aa, mt_aa, classification])
                        f_out.write('\t'.join(fields) + '\n')
                        classified += 1
                        
                    except (IndexError, ValueError, KeyError) as e:
                        # Skip problematic variants but continue processing
                        continue
            
            return total, classified, None
            
        except Exception as e:
            return total, classified, str(e)


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Classify SNVs as synonymous, non-synonymous, or nonsense. "
                    "The provided reference files (raw_files/proteins.dat and "
                    "raw_files/dna_transcripts.dat) cover human TFs only.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This tool analyzes single nucleotide variants (SNVs) in CDS regions and 
classifies them based on their effect on protein translation.

NOTE: The reference sequence files supplied with this repository
(raw_files/proteins.dat and raw_files/dna_transcripts.dat) are specific to
human transcription factors (TFs).

Classifications:
  - Syn: Synonymous (no amino acid change)
  - Non-Syn: Non-synonymous (amino acid change)
  - Nonsense: Creates stop codon

Input files are typically from intersect-variants output.

Examples:
  # Basic usage with provided human TF reference files
  classify-snvs \\
    -i cds_variants/intersections/ \\
    -o classified_variants/ \\
    -p raw_files/proteins.dat \\
    -d raw_files/dna_transcripts.dat \\
    -c cds_beds/sorted/

  # Use more workers
  classify-snvs \\
    -i cds_variants/intersections/ \\
    -o classified_variants/ \\
    -p raw_files/proteins.dat \\
    -d raw_files/dna_transcripts.dat \\
    -c cds_beds/sorted/ \\
    -w 16

Output format:
  Input columns + 3 additional columns:
  - WT amino acid
  - Mutant amino acid
  - Classification (Syn/Non-Syn/Nonsense)
        """
    )
    
    parser.add_argument(
        '-i', '--input-dir',
        type=Path,
        required=True,
        help='Directory containing variant BED files (from intersect-variants)'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for classified variants'
    )
    
    parser.add_argument(
        '-p', '--proteins',
        type=Path,
        required=True,
        help='Pickle file with protein sequences (ENST -> amino acid sequence). '
             'Provided for human TFs at raw_files/proteins.dat.'
    )
    
    parser.add_argument(
        '-d', '--dna-transcripts',
        type=Path,
        required=True,
        help='Pickle file with DNA transcripts (ENST -> nucleotide sequence). '
             'Provided for human TFs at raw_files/dna_transcripts.dat.'
    )
    
    parser.add_argument(
        '-c', '--sorted-cds-dir',
        type=Path,
        required=True,
        help='Directory with sorted CDS BED files'
    )
    
    parser.add_argument(
        '-w', '--workers',
        type=int,
        default=8,
        help='Number of parallel workers (default: 8)'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.input_dir.exists():
        print(f"Error: Input directory {args.input_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not args.proteins.exists():
        print(f"Error: Proteins file {args.proteins} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not args.dna_transcripts.exists():
        print(f"Error: DNA transcripts file {args.dna_transcripts} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not args.sorted_cds_dir.exists():
        print(f"Error: Sorted CDS directory {args.sorted_cds_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Initialize classifier
    print("Loading reference sequences...")
    classifier = SNVClassifier(
        args.proteins,
        args.dna_transcripts,
        args.sorted_cds_dir
    )
    print(f"  Loaded {len(classifier.proteins)} protein sequences")
    print(f"  Loaded {len(classifier.dna_transcripts)} DNA transcripts")
    
    # Get input files
    input_files = [f for f in args.input_dir.iterdir() if f.is_file() and f.suffix == '.bed']
    
    if not input_files:
        print(f"Warning: No .bed files found in {args.input_dir}", file=sys.stderr)
        sys.exit(0)
    
    print(f"\nProcessing {len(input_files)} variant files...")
    print(f"  Input directory: {args.input_dir}")
    print(f"  Output directory: {args.output_dir}")
    print()
    
    # Process files in parallel
    def process_file(input_file: Path) -> Tuple[str, int, int, Optional[str]]:
        """Process a single variant file.
        
        Args:
            input_file: Path to input variant file
            
        Returns:
            Tuple of (filename, total variants, classified variants, error or None)
        """
        output_file = args.output_dir / input_file.name
        total, classified, error = classifier.process_variant_file(input_file, output_file)
        
        if error:
            return input_file.name, total, classified, error
        return input_file.name, total, classified, None
    
    total_variants = 0
    total_classified = 0
    errors = []
    
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(process_file, f) for f in input_files]
        
        for future in tqdm(as_completed(futures), total=len(futures), 
                          desc="Processing SNVs"):
            filename, total, classified, error = future.result()
            
            total_variants += total
            total_classified += classified
            
            if error:
                errors.append((filename, error))
    
    # Summary
    print("\n" + "="*70)
    print("COMPLETE!")
    print("="*70)
    
    print(f"\nResults:")
    print(f"  Total variants processed: {total_variants:,}")
    print(f"  Successfully classified: {total_classified:,}")
    print(f"  Classification rate: {100*total_classified/total_variants:.1f}%" if total_variants > 0 else "  Classification rate: N/A")
    
    if errors:
        print(f"\nErrors encountered ({len(errors)} files):")
        for filename, error in errors[:10]:
            print(f"  {filename}: {error}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more")
    
    # Count classifications
    print("\nClassification summary:")
    syn_count = 0
    nonsyn_count = 0
    nonsense_count = 0
    
    for output_file in args.output_dir.glob('*.bed'):
        with open(output_file) as f:
            for line in f:
                classification = line.rstrip().split('\t')[-1]
                if classification == "Syn":
                    syn_count += 1
                elif classification == "Non-Syn":
                    nonsyn_count += 1
                elif classification == "Nonsense":
                    nonsense_count += 1
    
    total_class = syn_count + nonsyn_count + nonsense_count
    if total_class > 0:
        print(f"  Synonymous: {syn_count:,} ({100*syn_count/total_class:.1f}%)")
        print(f"  Non-synonymous: {nonsyn_count:,} ({100*nonsyn_count/total_class:.1f}%)")
        print(f"  Nonsense: {nonsense_count:,} ({100*nonsense_count/total_class:.1f}%)")
    
    print(f"\nOutput location: {args.output_dir}")


if __name__ == "__main__":
    main()