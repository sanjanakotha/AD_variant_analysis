"""
AD Variant Analysis
===================

A toolkit for analyzing single nucleotide variants (SNVs) in transcription
factor protein domains, including activation domains (ADs), DNA-binding domains
(DBDs), repression domains (RDs), and intrinsically disordered regions (IDRs).

Modules
-------
intersect_variants
    Intersect genomic variants with CDS BED files using bedtools.
classify_snvs
    Classify SNVs as synonymous, non-synonymous, or nonsense.
classify_domain_snvs
    Classify SNVs within specific protein domains.
intersect_domains_variants
    Intersect protein domain regions with genomic variants.
get_mutations_domains_snv_classified
    Get mutations in domains with SNV classification.

CLI Commands
------------
After installation (``pip install .``), the following commands are available:

- ``map-domains``             – map protein domain AA coordinates to genomic BED format
- ``intersect-variants``      – intersect variant BED files with CDS regions
- ``classify-snvs``           – classify CDS variants as Syn/Non-Syn/Nonsense
- ``classify-domain-snvs``    – classify SNVs within protein domains
- ``intersect-domain-variants`` – intersect domain BED files with variant files
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("ad-variant-analysis")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = [
    "intersect_variants",
    "classify_snvs",
    "classify_domain_snvs",
    "intersect_domains_variants",
    "get_mutations_domains_snv_classified",
]
