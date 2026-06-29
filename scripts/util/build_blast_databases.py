#!/usr/bin/env python3
"""Build BLAST+ databases from BioCompiler hazard signature database.

This script extracts hazardous sequence motifs from the biosecurity module
and builds BLAST+ nucleotide and protein databases for homology-based
biosecurity screening.

Usage:
    python scripts/util/build_blast_databases.py [--output-dir /path/to/dbs]

Requirements:
    - NCBI BLAST+ must be installed (makeblastdb)
    - biocompiler must be importable

Environment Variables:
    BIOCOMPILER_BLAST_DB_PATH: Default output directory for databases.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def _parse_fasta(filepath: str) -> dict[str, str]:
    """Parse a FASTA file and return {header: sequence} dict.

    Parameters
    ----------
    filepath : str
        Path to a FASTA file.

    Returns
    -------
    dict[str, str]
        Mapping of sequence headers to sequences.
    """
    sequences: dict[str, str] = {}
    current_header: str | None = None
    current_seq: list[str] = []

    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]  # first word after >
                current_seq = []
            else:
                current_seq.append(line)

    if current_header is not None:
        sequences[current_header] = "".join(current_seq)

    return sequences


def main() -> None:
    """Build BLAST+ hazard databases from BioCompiler biosecurity signatures."""
    parser = argparse.ArgumentParser(
        description="Build BLAST+ hazard databases from BioCompiler biosecurity signatures",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory for BLAST databases (default: $BIOCOMPILER_BLAST_DB_PATH or /tmp/biocompiler_blast_db)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force rebuild even if databases exist",
    )
    args = parser.parse_args()

    output_dir = args.output_dir or os.environ.get(
        "BIOCOMPILER_BLAST_DB_PATH", "/tmp/biocompiler_blast_db"
    )

    # Check BLAST+ availability
    from biocompiler.biosecurity.blast_integration import is_blast_available, build_hazard_db

    if not is_blast_available():
        print("ERROR: BLAST+ is not installed. Install NCBI BLAST+ first.")
        print("See: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/")
        sys.exit(1)

    # Extract hazard sequences from the biosecurity signature database
    from biocompiler.biosecurity.blast_integration import _build_hazard_sequences_from_signatures

    protein_seqs = _build_hazard_sequences_from_signatures("protein")
    dna_seqs = _build_hazard_sequences_from_signatures("dna")

    print(f"Found {len(protein_seqs)} protein hazard signatures")
    print(f"Found {len(dna_seqs)} DNA hazard signatures")

    # Check if databases already exist
    prot_db_path = os.path.join(output_dir, "hazard_protein")
    nucl_db_path = os.path.join(output_dir, "hazard_nucl")

    if not args.force:
        # Check if .phr (protein) or .nhr (nucleotide) files exist
        prot_exists = (
            Path(prot_db_path + ".phr").exists()
            or Path(prot_db_path + ".pdb").exists()
        )
        nucl_exists = (
            Path(nucl_db_path + ".nhr").exists()
            or Path(nucl_db_path + ".ndb").exists()
        )
        if prot_exists and nucl_exists:
            print("BLAST databases already exist. Use --force to rebuild.")
            print(f"Databases stored in: {output_dir}")
            return

    # Build protein database
    if protein_seqs:
        print(f"Building protein BLAST database: {prot_db_path}")
        try:
            build_hazard_db(
                protein_seqs, "BioCompiler_Hazard_Protein", prot_db_path, db_type="prot"
            )
            print(f"  -> Protein database built successfully ({len(protein_seqs)} sequences)")
        except RuntimeError as e:
            print(f"  -> ERROR building protein database: {e}")
            sys.exit(1)
    else:
        print("No protein signatures found, skipping protein database.")

    # Build nucleotide database
    if dna_seqs:
        print(f"Building nucleotide BLAST database: {nucl_db_path}")
        try:
            build_hazard_db(
                dna_seqs, "BioCompiler_Hazard_Nucl", nucl_db_path, db_type="nucl"
            )
            print(f"  -> Nucleotide database built successfully ({len(dna_seqs)} sequences)")
        except RuntimeError as e:
            print(f"  -> ERROR building nucleotide database: {e}")
            sys.exit(1)
    else:
        print("No DNA signatures found, skipping nucleotide database.")

    # Also build from full-length reference sequences if available
    ref_seqs_file = os.path.join(output_dir, "reference_sequences.fasta")
    if os.path.exists(ref_seqs_file):
        print(f"Found reference sequences file: {ref_seqs_file}")
        ref_seqs = _parse_fasta(ref_seqs_file)
        if ref_seqs:
            # Determine type based on content (heuristic: if all ACGT, it is DNA)
            sample_seq = next(iter(ref_seqs.values()))
            is_dna = all(c in "ACGTNacgtn" for c in sample_seq)

            if is_dna:
                ref_db_path = os.path.join(output_dir, "hazard_reference_nucl")
                db_type = "nucl"
                db_name = "BioCompiler_Hazard_Reference_Nucl"
            else:
                ref_db_path = os.path.join(output_dir, "hazard_reference_prot")
                db_type = "prot"
                db_name = "BioCompiler_Hazard_Reference_Prot"

            print(f"Building reference BLAST database: {ref_db_path}")
            try:
                build_hazard_db(ref_seqs, db_name, ref_db_path, db_type=db_type)
                print(f"  -> Reference database built successfully ({len(ref_seqs)} sequences)")
            except RuntimeError as e:
                print(f"  -> WARNING: Could not build reference database: {e}")
    else:
        # Check for bundled reference sequences in the package data directory
        package_data_dir = Path(__file__).resolve().parent.parent / "data"
        bundled_ref = package_data_dir / "hazard_reference_sequences.fasta"
        if bundled_ref.exists():
            print(f"Found bundled reference sequences: {bundled_ref}")
            ref_seqs = _parse_fasta(str(bundled_ref))
            if ref_seqs:
                # Separate protein and DNA reference sequences
                prot_refs = {}
                nucl_refs = {}
                for header, seq in ref_seqs.items():
                    is_dna = all(c in "ACGTNacgtn" for c in seq)
                    if is_dna:
                        nucl_refs[header] = seq
                    else:
                        prot_refs[header] = seq

                if prot_refs:
                    ref_prot_db = os.path.join(output_dir, "hazard_reference_prot")
                    print(f"Building reference protein BLAST database: {ref_prot_db}")
                    try:
                        build_hazard_db(
                            prot_refs,
                            "BioCompiler_Hazard_Reference_Prot",
                            ref_prot_db,
                            db_type="prot",
                        )
                        print(f"  -> Reference protein database built ({len(prot_refs)} sequences)")
                    except RuntimeError as e:
                        print(f"  -> WARNING: Could not build reference protein database: {e}")

                if nucl_refs:
                    ref_nucl_db = os.path.join(output_dir, "hazard_reference_nucl")
                    print(f"Building reference nucleotide BLAST database: {ref_nucl_db}")
                    try:
                        build_hazard_db(
                            nucl_refs,
                            "BioCompiler_Hazard_Reference_Nucl",
                            ref_nucl_db,
                            db_type="nucl",
                        )
                        print(f"  -> Reference nucleotide database built ({len(nucl_refs)} sequences)")
                    except RuntimeError as e:
                        print(f"  -> WARNING: Could not build reference nucleotide database: {e}")

    print()
    print("BLAST database build complete!")
    print(f"Databases stored in: {output_dir}")
    print(f"Set BIOCOMPILER_BLAST_DB_PATH={output_dir} to use these databases")


if __name__ == "__main__":
    main()
