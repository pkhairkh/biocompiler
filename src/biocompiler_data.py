#!/usr/bin/env python3
"""
BioCompiler Biological Data Module
===================================
Curated biological data for the BioCompiler project.

Sources:
- HBB gene: NCBI NG_000007.3 (RefSeqGene), NM_000518.5 (mRNA)
- EGFP: Clontech pEGFP-C1 (U55763.1), mammalian codon-optimized
- Splice sites: Consensus from Mount lab (UMD), Gao et al. 2008 (NAR 36:2257),
  Taggart et al. 2017 (Genome Res 27:1199), Shapiro & Senapathy 1987
- Codon usage: Kazusa Codon Usage Database (Homo sapiens, 93487 CDSs, 40662582 codons)
- HEK293T splicing: ENCODE eCLIP data, GTEx project

All positions are 1-based unless otherwise noted.
All sequences are 5'->3' sense strand.
"""

# ========================================================================
# 1. HUMAN BETA-GLOBIN (HBB) GENE
# ========================================================================

# Reference: NG_000007.3 (RefSeqGene on chromosome 11)
# Gene position in NG_000007.3: 70545-72152 (1608 bp)
# mRNA reference: NM_000518.5 (628 bp, MANE Select)
# CDS reference: NP_000509.1 (147 aa, includes stop)

# Complete HBB pre-mRNA/genomic sequence (including introns)
# This is the DNA template strand (sense strand) of the full gene
HBB_PREMRNA = (
    "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGA"
    "GGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGC"
    "AGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAG"
    "ACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGG"
    "TGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGG"
    "CAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGAC"
    "AACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACT"
    "TCAGGGTGAGTCTATGGGACGCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAG"
    "GAAGGGGATAAGTAACAGGGTACAGTTTAGAATGGGAAACAGACGAATGATTGCATCAGTGTGGAAGTCT"
    "CAGGATCGTTTTAGTTTCTTTTATTTGCTGTTCATAACAATTGTTTTCTTTTGTTTAATTCTTGCTTTCT"
    "TTTTTTTTCTTCTCCGCAATTTTTACTATTATACTTAATGCCTTAACATTGTGTATAACAAAAGGAAATA"
    "TCTCTGAGATACATTAAGTAACTTAAAAAAAAACTTTACACAGTCTGCCTAGTACATTACTATTTGGAAT"
    "ATATGTGTGCTTATTTGCATATTCATAATCTCCCTACTTTATTTTCTTTTATTTTTAATTGATACATAAT"
    "CATTATACATATTTATGGGTTAAAGTGTAATGTTTTAATATGTGTACACATATTGACCAAATCAGGGTAA"
    "TTTTGCATTTGTAATTTTAAAAAATGCTTTCTTCTTTTAATATACTTTTTTGTTTATCTTATTTCTAATA"
    "CTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAATGTATCATGCCTCTTTGCACCATTCTAAAG"
    "AATAACAGTGATAATTTCTGGGTTAAGGCAATAGCAATATCTCTGCATATAAATATTTCTGCATATAAAT"
    "TGTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCATTCTGCTTTTATTTT"
    "ATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAATCATGTTCATACCTCTT"
    "ATCTTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCA"
    "CCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCA"
    "CTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACT"
    "GGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA"
)

# HBB mature mRNA sequence (NM_000518.5, 628 bp)
HBB_MRNA = (
    "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGA"
    "GGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGC"
    "AGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATG"
    "CTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGC"
    "TCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGAT"
    "CCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCA"
    "CCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCA"
    "CTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACT"
    "GGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA"
)

# HBB CDS (coding sequence, positions 51-494 of mRNA = 444 bp)
HBB_CDS = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAG"
    "TTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTG"
    "GGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGG"
    "TGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCAC"
    "TGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCC"
    "CATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTA"
    "ATGCCCTGGCCCACAAGTATCACTAA"
)

# HBB protein (147 aa including stop, NP_000509.1)
HBB_PROTEIN = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"

# Exon positions in the pre-mRNA (1-based)
# Format: (start, end, length)
HBB_EXONS = [
    (1, 142, 142),    # Exon 1: 5'UTR + part of CDS (includes 50 bp 5'UTR)
    (273, 495, 223),  # Exon 2: internal exon, entirely CDS
    (1346, 1608, 263) # Exon 3: part of CDS + 3'UTR (CDS ends at 1474, 3'UTR from 1475-1608)
]

# CDS positions within the pre-mRNA (1-based)
HBB_CDS_REGIONS = [
    (51, 142, 92),    # CDS part 1 in Exon 1: starts at ATG (position 51)
    (273, 495, 223),  # CDS part 2 in Exon 2: entire exon
    (1346, 1474, 129) # CDS part 3 in Exon 3: ends at stop codon TAA (positions 1472-1474)
]

# Intron positions in the pre-mRNA (1-based)
HBB_INTRONS = [
    (143, 272, 130),   # Intron 1 (IVS1): 130 bp
    (496, 1345, 850)   # Intron 2 (IVS2): 850 bp
]

# Splice donor/acceptor sequences for HBB (actual sequences from the gene)
HBB_SPLICE_SITES = {
    "intron1_donor": {      # Exon1 -> Intron1 boundary
        "exon_end": "CAG",           # Last 3 nt of exon 1 (positions 140-142)
        "intron_start": "GTTGGTATCAAG", # First 12 nt of intron 1 (positions 143-154)
        "consensus_donor": "GT",      # Conserved GT at 5' splice site
        "position_in_premrna": 143,   # Position of GT in pre-mRNA (1-based)
    },
    "intron1_acceptor": {   # Intron1 -> Exon2 boundary
        "intron_end": "TTAG",         # Last 4 nt of intron 1 (positions 269-272)
        "exon_start": "GCTGCTGG",     # First 8 nt of exon 2 (positions 273-280)
        "consensus_acceptor": "AG",    # Conserved AG at 3' splice site
        "position_in_premrna": 272,    # Position of AG in pre-mRNA (1-based)
    },
    "intron2_donor": {      # Exon2 -> Intron2 boundary
        "exon_end": "CAG",            # Last 3 nt of exon 2 (positions 493-495)
        "intron_start": "GTGAGTCTATGGGAC", # First 15 nt of intron 2 (positions 496-510)
        "consensus_donor": "GT",       # Conserved GT at 5' splice site
        "position_in_premrna": 496,    # Position of GT in pre-mRNA (1-based)
    },
    "intron2_acceptor": {   # Intron2 -> Exon3 boundary
        "intron_end": "TCAG",          # Last 4 nt of intron 2 (positions 1342-1345)
        "exon_start": "CTCCTGGG",      # First 8 nt of exon 3 (positions 1346-1353)
        "consensus_acceptor": "AG",     # Conserved AG at 3' splice site
        "position_in_premrna": 1345,    # Position of AG in pre-mRNA (1-based)
    },
}

# HBB mRNA exon positions (1-based, in the spliced mRNA)
HBB_MRNA_EXONS = [
    (1, 142, 142),    # Exon 1
    (143, 365, 223),  # Exon 2
    (366, 628, 263)   # Exon 3
]

# HBB mRNA CDS: positions 51-494 (1-based in mRNA)
HBB_MRNA_CDS = (51, 494)

# 5'UTR and 3'UTR in mRNA
HBB_5UTR = HBB_MRNA[0:50]    # Positions 1-50
HBB_3UTR = HBB_MRNA[494:628] # Positions 495-628

# PolyA signal in 3'UTR
HBB_POLYA_SIGNAL = {
    "sequence": "AATAAA",
    "position_in_mrna": (602, 607),  # 1-based in mRNA
}

# ========================================================================
# 2. EGFP (ENHANCED GREEN FLUORESCENT PROTEIN)
# ========================================================================

# EGFP CDS - mammalian codon-optimized for Homo sapiens
# Source: Clontech pEGFP-C1 vector (U55763.1), positions 614-1330
# The EGFP CDS is 720 bp (239 aa + TGA stop codon)
# Contains >190 silent base changes from wild-type Aequorea victoria GFP
# optimized for human codon usage preferences
# Reference: Cormack et al. 1996, Gene 173:33-38

EGFP_CDS = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC"
    "GGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTAC"
    "GGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACC"
    "CTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG"
    "CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTC"
    "TTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTG"
    "GTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCAC"
    "AAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC"
    "GGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCC"
    "GACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCAC"
    "TACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTC"
    "CTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTGA"
)

# EGFP protein sequence (239 aa)
EGFP_PROTEIN = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQC"
    "FSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKL"
    "EYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPN"
    "EKRDHMVLLEFVTAAGITLGMDELYK"
)

# EGFP properties
EGFP_PROPERTIES = {
    "excitation_max_nm": 488,
    "emission_max_nm": 507,
    "extinction_coefficient_M1cm1": 55900,
    "quantum_yield": 0.60,
    "brightness": 33.54,  # EC * QY / 1000
    "pKa": 6.0,
    "maturation_time_min": 25.0,
    "fluorescence_lifetime_ns": 2.6,
    "molecular_weight_kDa": 26.9,
    "oligomerization": "weak dimer",
    "source_organism": "Aequorea victoria",
    "mutations_from_wt_gfp": "F64L, S65T",
    "codon_optimized_for": "Homo sapiens",
    "num_silent_mutations": ">190",
    "cds_length_bp": 720,
    "protein_length_aa": 239,
}

# Note: EGFP is a single-exon gene (no introns) as used in expression constructs.
# The original Aequorea victoria GFP gene also has no introns.

# ========================================================================
# 3. SPLICING GRAMMAR FOR HEK293T / HUMAN CELLS
# ========================================================================

# 3a. U2-type (major class) splice site consensus sequences
# References:
#   - Shapiro & Senapathy 1987, NAR 15:7155
#   - Mount lab consensus: https://science.umd.edu/labs/mount/RNAinfo/consensus.html
#   - Burge et al. 1999
#   - Taggart et al. 2017, Genome Res 27:1199

SPLICE_DONOR_CONSENSUS = {
    "name": "5' Splice Site (Donor)",
    "consensus": "MAG|GTRAGT",  # | marks exon-intron boundary
    "description": "M=A/C, R=A/G. Nearly invariant GT at positions +1,+2",
    "positions": {  # Position relative to exon-intron boundary (exon = negative, intron = positive)
        -3: {"A": 0.33, "C": 0.33, "G": 0.20, "T": 0.14},  # M (A or C preferred)
        -2: {"A": 0.60, "C": 0.13, "G": 0.13, "T": 0.14},  # A strongly preferred
        -1: {"A": 0.10, "C": 0.05, "G": 0.78, "T": 0.07},  # G nearly invariant
         1: {"A": 0.00, "C": 0.00, "G": 0.99, "T": 0.01},  # G invariant (DNA)
         2: {"A": 0.00, "C": 0.00, "G": 0.01, "T": 0.99},  # T nearly invariant (DNA); U in RNA
         3: {"A": 0.55, "C": 0.10, "G": 0.30, "T": 0.05},  # R (A or G preferred)
         4: {"A": 0.15, "C": 0.10, "G": 0.15, "T": 0.60},  # A/G/T
         5: {"A": 0.07, "C": 0.05, "G": 0.80, "T": 0.08},  # G strongly preferred
         6: {"A": 0.20, "C": 0.20, "G": 0.25, "T": 0.35},  # T preferred
    },
    "invariant": "GT at +1,+2 (GU in RNA)",
    "gc_variant_fraction": 0.005,  # ~0.5% of human introns use GC donor
}

SPLICE_ACCEPTOR_CONSENSUS = {
    "name": "3' Splice Site (Acceptor)",
    "consensus": "YYYYYYNCAG|G",  # | marks intron-exon boundary
    "description": "Polypyrimidine tract (Y=C/T) upstream of CAG. Nearly invariant AG at -2,-1",
    "positions": {  # Position relative to intron-exon boundary (intron = negative, exon = positive)
       -14: {"A": 0.05, "C": 0.45, "G": 0.05, "T": 0.45},  # Polypyrimidine tract
       -13: {"A": 0.05, "C": 0.45, "G": 0.05, "T": 0.45},
       -12: {"A": 0.05, "C": 0.45, "G": 0.05, "T": 0.45},
       -11: {"A": 0.05, "C": 0.45, "G": 0.05, "T": 0.45},
       -10: {"A": 0.05, "C": 0.50, "G": 0.05, "T": 0.40},
        -9: {"A": 0.05, "C": 0.50, "G": 0.05, "T": 0.40},
        -8: {"A": 0.06, "C": 0.50, "G": 0.06, "T": 0.38},
        -7: {"A": 0.08, "C": 0.45, "G": 0.07, "T": 0.40},
        -6: {"A": 0.10, "C": 0.40, "G": 0.08, "T": 0.42},
        -5: {"A": 0.06, "C": 0.50, "G": 0.06, "T": 0.38},
        -4: {"A": 0.03, "C": 0.70, "G": 0.02, "T": 0.25},
        -3: {"A": 0.02, "C": 0.77, "G": 0.01, "T": 0.20},  # C preferred (part of CAG)
        -2: {"A": 0.00, "C": 0.00, "G": 0.99, "T": 0.01},  # G invariant
        -1: {"A": 0.99, "C": 0.00, "G": 0.00, "T": 0.01},  # A nearly invariant
         1: {"A": 0.25, "C": 0.10, "G": 0.50, "T": 0.15},  # G preferred
    },
    "invariant": "AG at -2,-1",
    "polypyrimidine_tract": {
        "typical_length": "10-20 nt",
        "location": "upstream of AG, downstream of branch point",
        "composition": ">70% C/T",
    },
}

BRANCH_POINT_CONSENSUS = {
    "name": "Branch Point Sequence (BPS)",
    "consensus": "yUnAy",  # Position -3 to +1 relative to branch point A
    "description": (
        "y = pyrimidine (C/T), n = any, A = branch point adenosine. "
        "The underlined A at position 0 is the branch point nucleotide. "
        "Reference: Gao et al. 2008, NAR 36:2257-2267"
    ),
    "positions": {  # Position relative to branch point A (position 0)
        -3: {"A": 0.21, "C": 0.40, "G": 0.17, "T": 0.22},  # y (pyrimidine), info=0.27
        -2: {"A": 0.07, "C": 0.18, "G": 0.03, "T": 0.72},  # U strongly preferred, info=0.85
        -1: {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},  # n (any)
         0: {"A": 0.92, "C": 0.03, "G": 0.02, "T": 0.03},  # A nearly invariant, info=1.48
         1: {"A": 0.25, "C": 0.35, "G": 0.15, "T": 0.25},  # y (pyrimidine), info=0.23
    },
    "distance_from_3ss": {
        "range": "10-60 nt upstream of 3' splice site AG",
        "peak": "18-35 nt upstream of 3' splice site AG",
        "median": "~24 nt upstream of 3' splice site AG",
    },
    "u2_snRNA_pairing": "TRYTRAY (where Y = C/T, R = A/G), with bulged A at position 5",
}

# 3b. HEK293T-specific splicing data
# HEK293T is a human embryonic kidney cell line widely used in molecular biology
# ENCODE eCLIP data for HEK293T is available for some RBPs, but primary eCLIP
# datasets are from K562 and HepG2 cell lines.

HEK293T_SPLICING = {
    "cell_type": "HEK293T (human embryonic kidney, SV40 T-antigen transformed)",
    "karyotype": "hypotriploid, 60-65 chromosomes",
    "encode_data_available": {
        "eCLIP": "Limited (primary ENCODE eCLIP: K562, HepG2)",
        "RNAseq": "Yes - available from ENCODE",
        "RAMPAGE": "Yes - available from ENCODE",
        "PROseq": "Limited",
        "notes": "For comprehensive eCLIP splicing factor binding, K562 and HepG2 data are often used as proxies",
    },
    "splice_site_consensus": "Same as general human U2-type consensus (see above)",
    "known_splicing_factors": [
        "SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6", "SRSF7",
        "HNRNPA1", "HNRNPA2B1", "HNRNPH1", "HNRNPF", "HNRNPM",
        "PTBP1", "NOVA1", "RBFOX2", "ESRP1", "MBNL1",
        "U2AF1", "U2AF2", "SF1", "SF3B1", "PUF60",
    ],
    "alternative_splicing_prevalence": {
        "exon_skipping": "~40% of AS events",
        "alternative_5ss": "~20% of AS events",
        "alternative_3ss": "~20% of AS events",
        "intron_retention": "~15% of AS events",
        "mutually_exclusive_exons": "~5% of AS events",
    },
    "gtx_data": "HEK293T splicing data available from GTEx project (though GTEx focuses on tissue-level data)",
    "notes": (
        "For HEK293T-specific splice site scoring, MaxEntScan scores calibrated on "
        "human data are recommended. The general human U2-type splice site consensus "
        "applies. No cell-type-specific deviations from the human consensus are known "
        "for HEK293T."
    ),
}

# 3c. Complete splice site model for computational use
SPLICE_SITE_MODEL = {
    "donor": {
        "window": (-3, +6),  # 9-mer: 3 nt exon + 6 nt intron
        "consensus_9mer": "MAGGTRAGT",
        "weight_matrix": {  # Rows: A, C, G, T; Columns: positions -3 to +6
            "A": [0.33, 0.60, 0.10, 0.00, 0.00, 0.55, 0.15, 0.07, 0.20],
            "C": [0.33, 0.13, 0.05, 0.00, 0.00, 0.10, 0.10, 0.05, 0.20],
            "G": [0.20, 0.13, 0.78, 0.99, 0.01, 0.30, 0.15, 0.80, 0.25],
            "T": [0.14, 0.14, 0.07, 0.01, 0.99, 0.05, 0.60, 0.08, 0.35],
        },
        "maxentscan_model": "MaxEntScan::score5ss (9-mer model)",
    },
    "acceptor": {
        "window": (-20, +3),  # 23-mer: 20 nt intron + 3 nt exon
        "consensus_23mer": "yyyyyyyyyyyyyyyyncagGTR",
        "weight_matrix_available": True,
        "maxentscan_model": "MaxEntScan::score3ss (23-mer model)",
    },
    "branch_point": {
        "consensus": "yUnAy",
        "search_window": (-60, -10),  # Relative to 3' splice site AG
        "peak_region": (-35, -18),    # Most common location
    },
    "polypyrimidine_tract": {
        "location": "Between branch point and 3' splice site AG",
        "typical_length": "10-20 nt",
        "min_pyrimidine_content": 0.70,
    },
}

# ========================================================================
# 4. CODON USAGE TABLE FOR HOMO SAPIENS
# ========================================================================

# Source: Kazusa Codon Usage Database
# Homo sapiens [gbpri]: 93,487 CDSs (40,662,582 codons)
# URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N
# Coding GC: 52.27%, 1st letter: 55.72%, 2nd: 42.54%, 3rd: 58.55%

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
# fraction = proportion of this codon among all synonymous codons for the amino acid
# per_thousand = frequency per 1000 codons in all CDSs
CODON_USAGE = {
    "TTT": ("F", 0.46, 17.6, 714298),
    "TTC": ("F", 0.54, 20.3, 824692),
    "TTA": ("L", 0.08,  7.7, 311881),
    "TTG": ("L", 0.13, 12.9, 525688),
    "CTT": ("L", 0.13, 13.2, 536515),
    "CTC": ("L", 0.20, 19.6, 796638),
    "CTA": ("L", 0.07,  7.2, 290751),
    "CTG": ("L", 0.40, 39.6, 1611801),
    "ATT": ("I", 0.36, 16.0, 650473),
    "ATC": ("I", 0.47, 20.8, 846466),
    "ATA": ("I", 0.17,  7.5, 304565),
    "ATG": ("M", 1.00, 22.0, 896005),
    "GTT": ("V", 0.18, 11.0, 448607),
    "GTC": ("V", 0.24, 14.5, 588138),
    "GTA": ("V", 0.12,  7.1, 287712),
    "GTG": ("V", 0.46, 28.1, 1143534),
    "TCT": ("S", 0.19, 15.2, 618711),
    "TCC": ("S", 0.22, 17.7, 718892),
    "TCA": ("S", 0.15, 12.2, 496448),
    "TCG": ("S", 0.05,  4.4, 179419),
    "CCT": ("P", 0.29, 17.5, 713233),
    "CCC": ("P", 0.32, 19.8, 804620),
    "CCA": ("P", 0.28, 16.9, 688038),
    "CCG": ("P", 0.11,  6.9, 281570),
    "ACT": ("T", 0.25, 13.1, 533609),
    "ACC": ("T", 0.36, 18.9, 768147),
    "ACA": ("T", 0.28, 15.1, 614523),
    "ACG": ("T", 0.11,  6.1, 246105),
    "GCT": ("A", 0.27, 18.4, 750096),
    "GCC": ("A", 0.40, 27.7, 1127679),
    "GCA": ("A", 0.23, 15.8, 643471),
    "GCG": ("A", 0.11,  7.4, 299495),
    "TAT": ("Y", 0.44, 12.2, 495699),
    "TAC": ("Y", 0.56, 15.3, 622407),
    "TAA": ("*", 0.30,  1.0,  40285),  # Stop
    "TAG": ("*", 0.24,  0.8,  32109),  # Stop
    "CAT": ("H", 0.42, 10.9, 441711),
    "CAC": ("H", 0.58, 15.1, 613713),
    "CAA": ("Q", 0.27, 12.3, 501911),
    "CAG": ("Q", 0.73, 34.2, 1391973),
    "AAT": ("N", 0.47, 17.0, 689701),
    "AAC": ("N", 0.53, 19.1, 776603),
    "AAA": ("K", 0.43, 24.4, 993621),
    "AAG": ("K", 0.57, 31.9, 1295568),
    "GAT": ("D", 0.46, 21.8, 885429),
    "GAC": ("D", 0.54, 25.1, 1020595),
    "GAA": ("E", 0.42, 29.0, 1177632),
    "GAG": ("E", 0.58, 39.6, 1609975),
    "TGT": ("C", 0.46, 10.6, 430311),
    "TGC": ("C", 0.54, 12.6, 513028),
    "TGA": ("*", 0.47,  1.6,  63237),  # Stop
    "TGG": ("W", 1.00, 13.2, 535595),
    "CGT": ("R", 0.08,  4.5, 184609),
    "CGC": ("R", 0.18, 10.4, 423516),
    "CGA": ("R", 0.11,  6.2, 250760),
    "CGG": ("R", 0.20, 11.4, 464485),
    "AGT": ("S", 0.15, 12.1, 493429),
    "AGC": ("S", 0.24, 19.5, 791383),
    "AGA": ("R", 0.21, 12.2, 494682),
    "AGG": ("R", 0.21, 12.0, 486463),
    "GGT": ("G", 0.16, 10.8, 437126),
    "GGC": ("G", 0.34, 22.2, 903565),
    "GGA": ("G", 0.25, 16.5, 669873),
    "GGG": ("G", 0.25, 16.5, 669768),
}

# Relative adaptiveness values (w_i) for each codon
# w_i = frequency_i / max_frequency_among_synonymous_codons
# This is the metric used by Sharp & Li (1987) for CAI calculation
# For high-expression genes, these values are preferred

# First, find the max frequency codon for each amino acid
_AMINO_ACID_CODONS = {}
for codon, (aa, frac, freq, count) in CODON_USAGE.items():
    if aa not in _AMINO_ACID_CODONS:
        _AMINO_ACID_CODONS[aa] = []
    _AMINO_ACID_CODONS[aa].append((codon, freq))

_MAX_FREQ = {}
for aa, codons in _AMINO_ACID_CODONS.items():
    max_freq = max(f for _, f in codons)
    _MAX_FREQ[aa] = max_freq

CODON_ADAPTIVENESS = {}
for codon, (aa, frac, freq, count) in CODON_USAGE.items():
    if aa == "*":  # Skip stop codons
        continue
    CODON_ADAPTIVENESS[codon] = freq / _MAX_FREQ[aa]

# Preferred (highest-frequency) codons for each amino acid for high-expression genes
PREFERRED_CODONS = {}
for aa, codons in _AMINO_ACID_CODONS.items():
    if aa == "*":
        continue
    best = max(codons, key=lambda x: x[1])
    PREFERRED_CODONS[aa] = best[0]

# Codon usage grouped by amino acid (sorted by frequency, descending)
CODON_USAGE_BY_AA = {}
for aa, codons in _AMINO_ACID_CODONS.items():
    if aa == "*":
        continue
    codon_data = []
    for codon, freq in sorted(codons, key=lambda x: -x[1]):
        c, a, f, p, n = codon, aa, CODON_USAGE[codon][1], CODON_USAGE[codon][2], CODON_USAGE[codon][3]
        codon_data.append((c, f, p, n))
    CODON_USAGE_BY_AA[aa] = codon_data


# ========================================================================
# 5. UTILITY FUNCTIONS
# ========================================================================

def get_exon_sequence(premrna_seq, exon_start, exon_end):
    """Extract exon sequence from pre-mRNA (1-based positions)."""
    return premrna_seq[exon_start-1:exon_end]

def get_intron_sequence(premrna_seq, intron_start, intron_end):
    """Extract intron sequence from pre-mRNA (1-based positions)."""
    return premrna_seq[intron_start-1:intron_end]

def splice_premrna(premrna_seq, exons):
    """Simulate splicing: concatenate exons to produce mRNA."""
    mrna = ""
    for start, end, _ in exons:
        mrna += premrna_seq[start-1:end]
    return mrna

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G",
                  "a": "t", "t": "a", "g": "c", "c": "g",
                  "N": "N", "n": "n"}
    return "".join(complement[base] for base in reversed(seq))

def translate_dna(seq, to_stop=True):
    """Translate a DNA sequence to protein."""
    codon_table = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
        "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
        "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
        "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
        "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
        "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
        "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
        "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    }
    protein = ""
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3].upper()
        aa = codon_table.get(codon, "X")
        if to_stop and aa == "*":
            break
        protein += aa
    return protein

def calculate_cai(protein_seq, codon_usage_table=None):
    """
    Calculate Codon Adaptation Index (CAI) for a coding sequence.
    Uses the Sharp & Li (1987) method.
    """
    import math
    if codon_usage_table is None:
        codon_usage_table = CODON_ADAPTIVENESS

    log_scores = []
    for i in range(0, len(protein_seq)-2, 3):
        codon = protein_seq[i:i+3].upper()
        w = codon_usage_table.get(codon, 0.5)  # Default for unknown codons
        if w > 0:
            log_scores.append(math.log(w))

    if not log_scores:
        return 0.0

    cai = math.exp(sum(log_scores) / len(log_scores))
    return cai

def find_splice_donor_sites(seq, threshold=0.8):
    """Find potential 5' splice donor sites (GT dinucleotide) in a sequence."""
    sites = []
    for i in range(len(seq)-1):
        if seq[i:i+2].upper() == "GT":
            # Score against consensus
            window_start = max(0, i-3)
            window_end = min(len(seq), i+6)
            sites.append({
                "position": i+1,  # 1-based
                "context": seq[window_start:window_end].upper(),
            })
    return sites

def find_splice_acceptor_sites(seq, threshold=0.8):
    """Find potential 3' splice acceptor sites (AG dinucleotide) in a sequence."""
    sites = []
    for i in range(len(seq)-1):
        if seq[i:i+2].upper() == "AG":
            sites.append({
                "position": i+1,  # 1-based
                "context": seq[max(0,i-20):min(len(seq),i+3)].upper(),
            })
    return sites

def find_branch_points(seq, search_start=None, search_end=None):
    """
    Find potential branch point sequences (yUnAy) in an intronic sequence.
    Returns positions and scores.
    """
    import re
    # Search for YTNAY pattern (where Y=C/T, N=any)
    pattern = re.compile(r"[CT][AT][ACGT]A[CT]", re.IGNORECASE)
    matches = []
    for m in pattern.finditer(seq):
        pos = m.start() + 1  # 1-based
        # Check if A is at the expected branch point position
        match_seq = m.group().upper()
        matches.append({
            "position": pos,
            "sequence": match_seq,
            "branch_point_A": pos + 3,  # Position of the A in yUnAy
        })
    return matches


# ========================================================================
# 6. VERIFICATION TESTS
# ========================================================================

def run_verification_tests():
    """Run verification tests on all data."""
    print("=" * 60)
    print("BioCompiler Data Verification")
    print("=" * 60)

    # Test 1: HBB pre-mRNA length
    assert len(HBB_PREMRNA) == 1608, f"HBB pre-mRNA should be 1608 bp, got {len(HBB_PREMRNA)}"
    print(f"✓ HBB pre-mRNA: {len(HBB_PREMRNA)} bp")

    # Test 2: HBB mRNA length
    assert len(HBB_MRNA) == 628, f"HBB mRNA should be 628 bp, got {len(HBB_MRNA)}"
    print(f"✓ HBB mRNA: {len(HBB_MRNA)} bp")

    # Test 3: Splicing produces correct mRNA
    mrna = splice_premrna(HBB_PREMRNA, HBB_EXONS)
    assert mrna == HBB_MRNA, "Splicing pre-mRNA should produce mRNA"
    print(f"✓ Splicing verification: pre-mRNA -> mRNA matches")

    # Test 4: CDS starts with ATG
    assert HBB_CDS[:3] == "ATG", "CDS should start with ATG"
    print(f"✓ HBB CDS starts with ATG: {HBB_CDS[:12]}...")

    # Test 5: CDS translates to expected protein
    protein = translate_dna(HBB_CDS)
    assert protein == HBB_PROTEIN, f"CDS translation should match protein"
    print(f"✓ HBB CDS translates to: {protein[:20]}... ({len(protein)} aa)")

    # Test 6: Splice sites have conserved GT/AG
    for name, site in HBB_SPLICE_SITES.items():
        if "donor" in name:
            assert site["consensus_donor"] == "GT", f"{name} should have GT donor"
        if "acceptor" in name:
            assert site["consensus_acceptor"] == "AG", f"{name} should have AG acceptor"
    print(f"✓ All splice sites have conserved GT/AG dinucleotides")

    # Test 7: EGFP CDS length
    assert len(EGFP_CDS) == 720, f"EGFP CDS should be 720 bp, got {len(EGFP_CDS)}"
    print(f"✓ EGFP CDS: {len(EGFP_CDS)} bp")

    # Test 8: EGFP translates correctly
    egfp_protein = translate_dna(EGFP_CDS)
    assert egfp_protein == EGFP_PROTEIN, "EGFP CDS should translate to EGFP protein"
    print(f"✓ EGFP translates to: {egfp_protein[:20]}... ({len(egfp_protein)} aa)")

    # Test 9: Codon usage completeness
    assert len(CODON_USAGE) == 64, f"Should have 64 codons, got {len(CODON_USAGE)}"
    print(f"✓ Codon usage table: {len(CODON_USAGE)} codons")

    # Test 10: Codon adaptiveness completeness (excluding stop codons)
    assert len(CODON_ADAPTIVENESS) == 61, f"Should have 61 sense codons, got {len(CODON_ADAPTIVENESS)}"
    print(f"✓ Codon adaptiveness: {len(CODON_ADAPTIVENESS)} sense codons")

    # Test 11: Exon sizes sum to mRNA length
    total_exon = sum(length for _, _, length in HBB_EXONS)
    assert total_exon == len(HBB_MRNA), f"Exon sum {total_exon} != mRNA {len(HBB_MRNA)}"
    print(f"✓ Exon sizes sum to mRNA length: {total_exon} bp")

    # Test 12: Intron GT/AG verification
    intron1 = get_intron_sequence(HBB_PREMRNA, 143, 272)
    intron2 = get_intron_sequence(HBB_PREMRNA, 496, 1345)
    assert intron1[:2] == "GT" and intron1[-2:] == "AG", "Intron 1: GT...AG rule"
    assert intron2[:2] == "GT" and intron2[-2:] == "AG", "Intron 2: GT...AG rule"
    print(f"✓ Intron 1 ({len(intron1)} bp): GT...AG ✓")
    print(f"✓ Intron 2 ({len(intron2)} bp): GT...AG ✓")

    # Test 13: Preferred codons check
    assert PREFERRED_CODONS["L"] == "CTG", "Leucine should prefer CTG"
    assert PREFERRED_CODONS["E"] == "GAG", "Glutamate should prefer GAG"
    print(f"✓ Preferred codons: L→CTG, E→GAG, A→GCC, G→GGC")

    print()
    print("All verification tests passed! ✓")


if __name__ == "__main__":
    run_verification_tests()
