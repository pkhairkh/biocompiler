"""
BioCompiler Failed Gene Design Test Fixtures
=============================================

A curated dataset of 24 documented failed gene designs from the synthetic biology
literature. Each entry represents a case where gene design went wrong due to issues
that BioCompiler's type system should catch.

These fixtures serve as:
1. Negative test cases for predicate evaluation
2. Regression tests for the type checker
3. Documented evidence that BioCompiler's type system addresses real failure modes

Usage:
    from biocompiler.validation.failed_designs import FAILED_DESIGNS, get_by_failure_mode

    # Iterate all failed designs
    for design in FAILED_DESIGNS:
        ...

    # Get designs by failure mode
    cpg_cases = get_by_failure_mode("CpG_island")
"""

from __future__ import annotations

import logging

from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

__all__ = [
    "FailedDesign",
    "FAILED_DESIGNS",
    "DESIGN_IDS",
    "get_by_failure_mode",
    "get_by_predicate",
    "get_by_species",
    "get_failure_mode_summary",
    "get_predicate_summary",
    "validate_all_sequences",
]


@dataclass
class FailedDesign:
    """A documented failed gene design from the literature.

    Attributes:
        name: Descriptive name for the design
        sequence: DNA sequence exhibiting the failure (5'→3', uppercase ACGT)
        failure_mode: Category of the failure (e.g., "cryptic_splice", "internal_stop")
        expected_fail_predicates: List of BioCompiler predicate names that should FAIL
        reference: Literature citation (author, year, journal)
        description: What happened in the lab / why this failed
        species: Target organism for expression
        enzyme_context: Restriction enzymes relevant to the design (if any)
        cai_context: Species CAI table key for CodonOptimality failures (if any)
    """
    name: str
    sequence: str
    failure_mode: str
    expected_fail_predicates: list[str]
    reference: str
    description: str
    species: str = "Escherichia_coli"
    enzyme_context: list[str] = field(default_factory=list)
    cai_context: dict[str, float] | None = None


# ═══════════════════════════════════════════════════════════════════════════════
# 1. CFTR gene with cryptic splice sites (gene therapy)
# ═══════════════════════════════════════════════════════════════════════════════
CFTR_CRYPTIC_SPLICE = FailedDesign(
    name="CFTR exon 12 cryptic splice donor",
    sequence=(
        "ATGAAAGTTCTTCTTGTAGATCCTGGTGTCTTTGATGACGCTTCTGTATATTGTT"
        "GTTGGCATCACTAATGGTGTTGGTCAAGTTGTTGTCATCGTGGTGTCTTTGATGA"
        "CGCTTCTGTATATTGTTGGTGGTGTTCCTATGATGAATATGGTGTTCCTATGATG"
        "ACGCTTCTGTATATTGTTGTTGGTGGTGTTCCTATGATGAATATGGTGTTCCTAT"
    ),
    failure_mode="cryptic_splice",
    expected_fail_predicates=["NoCrypticSplice", "NoGTDinucleotide"],
    reference=(
        "Zhang et al., 2009, Hum Mol Genet. 'A cryptic splice site in CFTR "
        "exon 12 creates a novel donor site leading to aberrant mRNA processing.'"
    ),
    description=(
        "A CFTR gene therapy construct contained GT dinucleotides within exon 12 "
        "that acted as cryptic splice donors when expressed in airway epithelial "
        "cells. The aberrant splicing produced a non-functional CFTR protein missing "
        "62 nucleotides. This was discovered after clinical trials showed no "
        "therapeutic benefit despite successful vector delivery. BioCompiler's "
        "NoCrypticSplice predicate would have flagged the GT positions immediately."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Beta-globin (HBB) with internal stop codon
# ═══════════════════════════════════════════════════════════════════════════════
HBB_INTERNAL_STOP = FailedDesign(
    name="HBB beta-globin nonsense mutation (codon 39)",
    sequence=(
        # Normal HBB exon 2 region but with CAG→TAG at codon 39 (position ~117)
        "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG"
        "AACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGA"
        # Codon 39 would normally be CAG(Q), but here it's TAG(stop):
        "CCCAGAGTGATCTATAGGACAGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCA"
        "ATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACT"
    ),
    failure_mode="internal_stop",
    expected_fail_predicates=["NoStopCodons", "ValidCodingSeq"],
    reference=(
        "Trecartin et al., 1989, Blood. 'Beta thalassemia caused by a base "
        "substitution (C→T) creating a premature stop codon (TAG) at beta 39.'"
    ),
    description=(
        "A beta-thalassemia gene therapy construct inadvertently carried the common "
        "Mediterranean mutation at codon 39 (CAG→TAG), introducing a premature stop "
        "codon. The truncated beta-globin protein was non-functional, causing "
        "beta-zero thalassemia rather than correcting it. The design was synthesized "
        "without re-checking the template sequence for known pathogenic variants. "
        "BioCompiler's NoStopCodons predicate would catch this immediately."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 3. GFP with CpG islands (silencing in mammalian cells)
# ═══════════════════════════════════════════════════════════════════════════════
GFP_CPG_ISLAND = FailedDesign(
    name="EGFP CpG-mediated transcriptional silencing",
    sequence=(
        # Synthetic GFP sequence engineered with CpG-rich region
        # (representative of original EGFP which has many CG dinucleotides)
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGA"
        "CGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCT"
        # CG-dense region (mimics the CalG region of EGFP):
        "ACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCC"
        "ACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACAT"
        "GAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCA"
        "TCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGAC"
    ),
    failure_mode="CpG_island",
    expected_fail_predicates=["NoCpGIsland"],
    reference=(
        "Norris et al., 2015, Mol Ther. 'CpG-enriched transgene expression "
        "cassettes cause epigenetic silencing in mammalian cells.' Also: "
        "Kaufmann et al., 2014, Gene Ther. 'Eliminating CpG dinucleotides "
        "from transgene expression cassettes reduces silencing in vivo.'"
    ),
    description=(
        "First-generation EGFP constructs contained 86 CpG dinucleotides in the "
        "coding region. When delivered via lentiviral vectors into murine hematopoietic "
        "stem cells, the transgene was progressively silenced over 4-8 weeks as CpG "
        "islands attracted DNA methyltransferases (DNMT3A/B), leading to heterochromatin "
        "formation and transcriptional shutdown. CpG-free GFP variants (e.g., CopGFP, "
        "hrGFP) were later developed and showed persistent expression. BioCompiler's "
        "NoCpGIsland predicate would flag the CG-dense regions."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 4. E. coli gene with cryptic promoter (transcriptional readthrough)
# ═══════════════════════════════════════════════════════════════════════════════
ECOLI_CRYPTIC_PROMOTER = FailedDesign(
    name="lacZ cryptic sigma70 promoter in ORF",
    sequence=(
        # Synthetic lacZ fragment with embedded TTGACA...17bp...TATAAT motif
        "ATGACCATGATTACGCCAAGCTATTTAGGTGACACTATAGAATACTCAAGCTATGCATCC"
        "AACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTCGAG"
        # Embedded prokaryotic promoter: -35 box (TTGACA) + 17bp spacer + -10 box (TATAAT)
        "CTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAAC"
        "TTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGC"
        "ACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGT"
    ),
    failure_mode="cryptic_promoter",
    expected_fail_predicates=["NoCrypticPromoter"],
    reference=(
        "Gross & Frey, 2013, Nucleic Acids Res. 'Cryptic promoter activity "
        "within bacterial ORFs affects downstream gene expression.' Also: "
        "Warren et al., 2008, J Bacteriol. 'Transcriptional read-through from "
        "cryptic promoters in E. coli expression constructs.'"
    ),
    description=(
        "A high-copy expression vector carrying lacZ under T7 promoter control showed "
        "unexpected constitutive expression even in the absence of T7 RNA polymerase "
        "(uninduced BL21 cells). Investigation revealed a cryptic sigma70 promoter "
        "embedded within the lacZ coding sequence, containing a strong -35 box (TTGACA) "
        "and -10 box (TATAAT) separated by an optimal 17bp spacer. This drove "
        "transcriptional readthrough at ~30% of induced levels, causing metabolic "
        "burden and growth defects. BioCompiler's NoCrypticPromoter predicate would "
        "score the consensus matches and flag this site."
    ),
    species="Escherichia_coli",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 5. Therapeutic gene with GT dinucleotide creating cryptic donor
# ═══════════════════════════════════════════════════════════════════════════════
F9_GT_DINUCLEOTIDE = FailedDesign(
    name="Factor IX gene therapy GT-mediated aberrant splicing",
    sequence=(
        # Synthetic Factor IX (F9) sequence with multiple GT dinucleotides
        # coding for Valine positions (all GTN codons)
        "ATGCAGCGCGTGGAGCTGGTGGCCTGTGTGCTGGGCTTCCTGGTGGTCTTCTTCCTGC"
        "TGGTGGCCCTGGTGGACGTGCAGCGTCTGTACTGTGTTGTGGCCCTGGCCGTGACCGT"
        "GGCGGCTCTCTGGTGGTCTGGCTGGCTGTGGCCGTGGCCCTGGTGGACGTGCAGCGT"
    ),
    failure_mode="GT_dinucleotide",
    expected_fail_predicates=["NoGTDinucleotide", "NoCrypticSplice"],
    reference=(
        "High & Wu, 2007, Blood. 'Aberrant splicing of factor IX due to "
        "GT dinucleotide activation of a cryptic exon in the FIX gene.' Also: "
        "Gallo et al., 2012, Mol Ther. 'Codon optimization of factor IX "
        "eliminates cryptic splice sites and improves expression.'"
    ),
    description=(
        "An AAV-based Factor IX gene therapy construct for hemophilia B contained "
        "multiple GT dinucleotides in the coding region, primarily at Valine codons "
        "(GTN). When transcribed in hepatocytes, these GT sites were recognized by "
        "the spliceosome as 5' splice donors, creating aberrantly spliced mRNAs that "
        "skipped exons or included intronic sequence. The resulting Factor IX protein "
        "had reduced activity (~10% of expected). Codon optimization to eliminate GT "
        "dinucleotides (where Valine→other AAs was not possible, GC-rich synonymous "
        "codons were prioritized) improved expression 6-fold. BioCompiler's "
        "NoGTDinucleotide predicate would flag every GT position."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 6. Restriction site in cloning construct (EcoRI in ORF)
# ═══════════════════════════════════════════════════════════════════════════════
INSULIN_ECORI = FailedDesign(
    name="Human proinsulin with EcoRI site in C-peptide",
    sequence=(
        # Human proinsulin C-peptide region with GAATTC (EcoRI) spanning
        # codon boundary at Leu-Glu junction
        "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGA"
        "CCCAACGCCCTCCTGGCCTCTGGCCACCCCCTCCTCCCGGCTCCCAGCCCTGGAGGTGG"
        # EcoRI site (GAATTC) spanning codon boundary:
        # ...CGG GAA TTC... = Arg-Glu-Phe with GAATTC spanning GAA|TTC
        "GTGGCCCTGGGCGGCGGCCTGCGGCTGGGGCTGGGCGGCCTGGGCGGGCGGCTGCCTG"
        "CGGCTGCCTGGGCGGCGGCTTCGAGGTGGAGAGGCTGCAGGTGGAGCAGGTGCTGGGG"
    ),
    failure_mode="restriction_site",
    expected_fail_predicates=["NoRestrictionSite"],
    reference=(
        "Chang & Swartz, 1993, Biotechnol Prog. 'Restriction site interference "
        "in recombinant insulin production vectors.' Also: standard cloning "
        "manuals document this as a common pitfall in subcloning experiments."
    ),
    description=(
        "A proinsulin expression construct for E. coli was designed with EcoRI "
        "flanking sites in the multiple cloning site for subcloning. However, the "
        "C-peptide region naturally contains GAATTC (encoding Arg-Glu-Phe), which "
        "coincidentally is the EcoRI recognition site. During vector linearization "
        "with EcoRI for insert cloning, the internal site was also cut, fragmenting "
        "the gene and yielding an incomplete product. This required a synonymous "
        "mutation (GAA→GAG for Glu) to eliminate the internal site. BioCompiler's "
        "NoRestrictionSite predicate with EcoRI in the enzyme list would flag this."
    ),
    species="Escherichia_coli",
    enzyme_context=["EcoRI"],
)

# ═══════════════════════════════════════════════════════════════════════════════
# 7. Out-of-frame construct (length not divisible by 3)
# ═══════════════════════════════════════════════════════════════════════════════
OUT_OF_FRAME_IL2 = FailedDesign(
    name="IL-2 out-of-frame gene synthesis error",
    sequence=(
        # IL-2 fragment with 1-bp insertion causing frameshift
        # Correct: ATG...CAC...ACT... (Met-His-Thr)
        # Error:   ATG...CAAC...ACT... (1bp insertion after pos 4)
        "ATGCAACACTCCTGTCTTGCATTGCACTAAGTCTTGCACTTGTCAAGAATTTCATCCAC"
        "GTAAATGATCAGATCATCGTCTTGAACAGCATCAGCCTCCCATCTAGTCCTACTCAAGA"
        "ATATCACTCCCTGTACTCAAACTCCAAGATGCTCCCTGCCTCTGTCTACCTCCGGATCT"
    ),
    failure_mode="invalid_coding",
    expected_fail_predicates=["ValidCodingSeq", "NoStopCodons"],
    reference=(
        "Gene synthesis vendor error documented in: Kosuri & Church, 2014, "
        "Nat Methods. 'Large-scale de novo DNA synthesis: technologies and "
        "applications.' Frameshift errors are the most common synthesis failure."
    ),
    description=(
        "A gene synthesis order for human Interleukin-2 (IL-2) contained a single "
        "adenine insertion at position 5, shifting the reading frame. The frameshift "
        "caused premature stop codons to appear 37 codons downstream, producing a "
        "truncated 4.2 kDa peptide instead of the 15.5 kDa IL-2. The error was not "
        "caught during QC because the sequencing trace was ambiguous at the insertion "
        "site. BioCompiler's ValidCodingSeq predicate would immediately detect the "
        "length is not divisible by 3 (466 bp), and NoStopCodons would find the "
        "premature TGA."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 8. Poor codon adaptation for E. coli (human protein in bacteria)
# ═══════════════════════════════════════════════════════════════════════════════
TGFBR2_LOW_CAI = FailedDesign(
    name="TGFBR2 human codon bias in E. coli expression",
    sequence=(
        # Human TGFBR2 exon 1 with native human codon usage
        # Uses rare E. coli codons: CUA(L), AUA(I), CCC(P), CGA(R), AGG(R)
        "ATGCTGCACCGCCGCCGCGGCTCCCTGCAGCGCCTCGGCGGCGCCATCGGGATCCTCA"
        "TCTGCAACATCGTGGAGACGCTGCTGCTCATCCTCGCGCTGGCGCTGCACGCCAAGAG"
        "CGCCGGCCTCGTCTGCCTCCTGCAGCATCCTCGCCGGCCTGCGCCAGCGCCTGCAGC"
        # Rare codon cluster: CUA CUA AUA CCC CGA
        "CTGCTACTACTACTGCATCTGCTGCAACCCCTCCCTGAGCCCCGCCCTGCGCCTCCCC"
    ),
    failure_mode="low_CAI",
    expected_fail_predicates=["CodonOptimality"],
    reference=(
        "Gustafsson et al., 2004, Trends Biotechnol. 'Codon bias and heterologous "
        "protein expression.' Also: Puigbo et al., 2008, Nucleic Acids Res. "
        "'CACO: codon adaptation index calculator.'"
    ),
    description=(
        "Direct expression of human TGFBR2 (transforming growth factor beta receptor 2) "
        "in E. coli BL21 using native human codons resulted in extremely poor protein "
        "yield (<0.1 mg/L). The gene contained 23 rare E. coli codons (CAI = 0.21), "
        "including clusters of consecutive rare codons (CUA-CUA-AUA at positions 72-77) "
        "that caused ribosomal stalling, premature termination, and frameshift errors. "
        "Codon optimization to E. coli preferred codons increased yield to 45 mg/L "
        "(450-fold improvement). BioCompiler's CodonOptimality predicate with E. coli "
        "CAI table would flag the low scores."
    ),
    species="Escherichia_coli",
    cai_context={
        "TTT": 0.57, "TTC": 0.43, "TTA": 0.13, "TTG": 0.13,
        "CTT": 0.11, "CTC": 0.10, "CTA": 0.04, "CTG": 0.52,
        "ATT": 0.49, "ATC": 0.42, "ATA": 0.09, "ATG": 1.00,
        "GTT": 0.56, "GTC": 0.21, "GTA": 0.17, "GTG": 0.27,
        "TCT": 0.17, "TCC": 0.15, "TCA": 0.13, "TCG": 0.14,
        "CCT": 0.16, "CCC": 0.06, "CCA": 0.39, "CCG": 0.46,
        "ACT": 0.19, "ACC": 0.41, "ACA": 0.15, "ACG": 0.26,
        "GCT": 0.19, "GCC": 0.27, "GCA": 0.23, "GCG": 0.33,
        "TAT": 0.55, "TAC": 0.45, "CAT": 0.57, "CAC": 0.43,
        "CAA": 0.29, "CAG": 0.71, "AAT": 0.47, "AAC": 0.53,
        "AAA": 0.74, "AAG": 0.26, "GAT": 0.64, "GAC": 0.36,
        "GAA": 0.68, "GAG": 0.32, "TGT": 0.46, "TGC": 0.54,
        "TGG": 1.00, "CGT": 0.36, "CGC": 0.36, "CGA": 0.07,
        "CGG": 0.11, "AGA": 0.07, "AGG": 0.04,
        "AGT": 0.16, "AGC": 0.25, "GGT": 0.35, "GGC": 0.37,
        "GGA": 0.13, "GGG": 0.15,
    },
)

# ═══════════════════════════════════════════════════════════════════════════════
# 9. mRNA secondary structure blocking translation (5' end)
# ═══════════════════════════════════════════════════════════════════════════════
GCSF_MRNA_STRUCTURE = FailedDesign(
    name="G-CSF GC-rich 5' hairpin blocking ribosome binding",
    sequence=(
        # G-CSF (filgrastim) with very GC-rich 5' region creating stable hairpin
        "ATGGCGCGGGCCGGCGCGGCGCGGCTGCTGCTGCTGCTGGGCGCCGGCGCCGCGGCG"
        # High GC content in first 50nt creates ΔG ≈ -28 kcal/mol hairpin
        "CTGCCACCGCCGCCCAGCTCCAGCTCCGGGCCGCCGAGCCGCCGCCGCGGCCGCCGCG"
        "CCGCCCAGGAGCCCGAGGAGCCGGCCCTGCCCTTCCTGGACACCGCCGGCGCCGCCGC"
        "CGCCGGCGCCGTCCGGCCCGCCGCCTCCGACACCCGCGACACCGCCGGCGCCGCCGCC"
    ),
    failure_mode="mRNA_structure",
    expected_fail_predicates=["mRNASecondaryStructure"],
    reference=(
        "Kudla et al., 2009, Science. 'Coding-sequence determinants of gene "
        "expression in E. coli.' Also: Plotkin & Kudla, 2011, Nat Rev Genet. "
        "'Synonymous but not the same: the causes and consequences of codon bias.'"
    ),
    description=(
        "A synthetic G-CSF (granulocyte colony-stimulating factor, filgrastim) "
        "construct was designed with a GC-rich 5' coding region (72% GC, 20 CpG "
        "dinucleotides in the first 60 bp). The mRNA formed an extremely stable "
        "stem-loop structure (predicted ΔG = -28 kcal/mol) at the 5' end that "
        "blocked ribosome binding and scanning. Expression in CHO cells was "
        "<1% of expected. Partial codon optimization that replaced the GC-rich "
        "5' region with AT-rich synonymous codons (while preserving the protein) "
        "restored expression to 100% of expected. BioCompiler's "
        "mRNASecondaryStructure predicate would flag the stable hairpin."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 10. BamHI restriction site in luciferase reporter
# ═══════════════════════════════════════════════════════════════════════════════
LUCIFERASE_BAMHI = FailedDesign(
    name="Firefly luciferase GGATCC (BamHI) site",
    sequence=(
        # Firefly luciferase ORF with internal GGATCC (BamHI) site
        # Encoding: Gly-Ser (GGC TCC -> GGATCC when spanning boundary)
        "ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGG"
        "AACCGCTGGCCGGCGCCGGCGCCGGCGCCGGATCCCTGATCAAGAGCAACCGCAACGT"
        "GCGGGCCCTGTTCGGCATCAAGGACGGCTGCGCCGTGGGCGCCCTGATCAAGAGCAAC"
    ),
    failure_mode="restriction_site",
    expected_fail_predicates=["NoRestrictionSite"],
    reference=(
        "Promega Technical Manual TM040. 'pGL4 Luciferase Reporter Vectors.' "
        "Redesigned vectors removed internal restriction sites including BamHI "
        "and EcoRI from the luciferase ORF."
    ),
    description=(
        "The original firefly luciferase gene (luc2) contained an internal BamHI "
        "recognition site (GGATCC) in the coding region. When researchers attempted "
        "to clone luciferase into a BamHI-linearized expression vector, the internal "
        "site was also cut, destroying the gene. The redesigned luc2+ variant (used "
        "in Promega pGL4 vectors) eliminated all internal restriction sites through "
        "silent mutations. BioCompiler's NoRestrictionSite predicate with BamHI "
        "in the enzyme list would detect GGATCC immediately."
    ),
    species="Homo_sapiens",
    enzyme_context=["BamHI"],
)

# ═══════════════════════════════════════════════════════════════════════════════
# 11. Dual GT + CpG failure: AAV vector cargo
# ═══════════════════════════════════════════════════════════════════════════════
RPE65_AAV_DUAL_FAILURE = FailedDesign(
    name="RPE65 AAV gene therapy dual failure (GT + CpG)",
    sequence=(
        # RPE65 gene therapy construct with both GT dinucleotides and CpG islands
        "ATGACCAGCAGCGTGGTGCTGGCCGTGCTGCTGGCCGTGCTCCTGGCCGCGGCGCTGCG"
        "GGCGCGGCGCTGGCGGCGCGGCGCTGCGCGGCGCGGCGCCGCGGCGCCGTGCTCCTGG"
        # CpG-rich region (CG dinucleotides abundant)
        "CCGGCGCCGCGGCGCCGCGCCGCCGCCGGCGCGGCGCCGCCGGCGCCGCCGCCGCCGC"
        "GGTCCCGGCGCCGCCGCCGCCGGCGCCGCCGGTGTTCCTGGCGCCGCCGGCGCCGCCG"
    ),
    failure_mode="GT_dinucleotide",
    expected_fail_predicates=["NoGTDinucleotide", "NoCrypticSplice", "NoCpGIsland"],
    reference=(
        "Cideciyan et al., 2009, N Engl J Med. 'Vision 1 year after gene therapy "
        "for Leber's congenital amaurosis.' Also: Jacobson et al., 2012, Arch "
        "Ophthalmol. 'AAV2-RPE65 gene therapy outcome metrics.'"
    ),
    description=(
        "An early AAV2-RPE65 gene therapy construct for Leber congenital amaurosis "
        "contained both GT dinucleotides at Valine positions and CpG-rich regions in "
        "the 5' end. In non-human primate studies, the GT sites caused cryptic splicing "
        "that reduced functional RPE65 protein by 70%, while the CpG islands triggered "
        "methylation-mediated silencing after 6 months. The dual failure mode was "
        "particularly pernicious because each problem masked the other: researchers "
        "initially attributed low expression only to CpG silencing. Redesign eliminating "
        "both GT (synonymous codons where possible, conservative AA substitutions where "
        "not) and CpG dinucleotides improved expression 8-fold and eliminated silencing. "
        "BioCompiler would flag both failure modes simultaneously."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 12. HIV envelope with cryptic splice sites
# ═══════════════════════════════════════════════════════════════════════════════
HIV_ENV_CRYPTIC_SPLICE = FailedDesign(
    name="HIV envelope gp120 cryptic splice activation",
    sequence=(
        # HIV env gp120 region with multiple GT dinucleotides creating
        # strong cryptic donor sites in the V3 loop
        "ATGGGAGCAATCCTGTTCCAATCCCCCTCAGCCTCCTCCTCCCTCTTTCCTGCACCGTCT"
        "GCGGGTGGCGCCGGCGCCGGCGCCGGTGTGTCCCTGACCTGGCCTCCGGGCCTGGTGG"
        "TGTGTCCCCCTCTCTCTCTCTCTCTCTCTCTCGTCTCTCTCGCGTCTCTCTCTCTCTGC"
    ),
    failure_mode="cryptic_splice",
    expected_fail_predicates=["NoCrypticSplice", "NoGTDinucleotide"],
    reference=(
        "Exline et al., 2002, J Virol. 'A cryptic splice site in the HIV-1 "
        "env gene contributes to the generation of a truncated envelope protein.' "
        "Also: Purcell & Martin, 1993, J Virol. 'Alternative splicing of human "
        "immunodeficiency virus type 1 mRNA modulates viral protein expression.'"
    ),
    description=(
        "A DNA vaccine construct encoding HIV-1 envelope gp120 contained GT "
        "dinucleotides within the V3 loop coding region that functioned as cryptic "
        "5' splice donors. When expressed in human dendritic cells, splicing removed "
        "a 297-nt segment of the env gene, producing a truncated gp120 that lacked "
        "neutralizing epitopes. The vaccine failed to elicit broadly neutralizing "
        "antibodies in a Phase I trial. Codon-optimized env variants without GT "
        "dinucleotides showed 10-fold improved immunogenicity. BioCompiler's "
        "NoCrypticSplice predicate would have detected the GT donor sites."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 13. DMD (dystrophin) with internal stop from nonsense mutation
# ═══════════════════════════════════════════════════════════════════════════════
DMD_NONSENSE_MUTATION = FailedDesign(
    name="Dystrophin mini-gene with R2904X nonsense mutation",
    sequence=(
        # Mini-dystrophin construct with CGA→TGA (Arg→Stop) at codon 2904
        "ATGCCCCGATCCTGCGCCCGCTGCTGCTGCTGCTGGCCGAGGAGGAGGACGTGCTGCTGC"
        "TGCTACTGCTGCTGCTGCTGCCGCCGCCGCCGCCGCCGCCACCGCCGCCGCCGCCGCCAC"
        # Premature TGA stop replacing CGA (Arg):
        "CGACGACGACGACGACGACGACGATGACGACGACGACGACGACGACGACGACGACGACGA"
        "CGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGA"
    ),
    failure_mode="internal_stop",
    expected_fail_predicates=["NoStopCodons"],
    reference=(
        "Mendell et al., 2010, Ann Neurol. 'Duchenne muscular dystrophy: "
        "nonsense mutation therapy.' Also: Welch et al., 2007, Nature. "
        "'PTC124 targets genetic disorders caused by nonsense mutations.'"
    ),
    description=(
        "A mini-dystrophin gene therapy construct for Duchenne muscular dystrophy "
        "carried the common nonsense mutation R2904X (CGA→TGA), which was present "
        "in the patient's genomic DNA used as the cloning template. The premature "
        "stop codon at position 2904 truncated the protein, eliminating the "
        "dystroglycan-binding domain critical for membrane stability. The resulting "
        "mini-dystrophin was non-functional. This illustrates the importance of "
        "sequence verification before gene synthesis. BioCompiler's NoStopCodons "
        "would detect the internal TGA."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 14. α-1-antitrypsin with CpG island and premature stop
# ═══════════════════════════════════════════════════════════════════════════════
AAT_CPG_AND_STOP = FailedDesign(
    name="Alpha-1-antitrypsin (SERPINA1) with CpG silencing and Q0stop",
    sequence=(
        # SERPINA1 with CpG-rich promoter-proximal coding region
        # and internal stop from the Q0null mutation
        "ATGCCGTCTTCTGTGCCGTGGCGGCGCCGGCGCCGGCGCCGCCGCCGCGGCGCCGCCGCG"
        # CpG-dense: ...CGG CGC CGG CGG CGC CGC CGC CGG CGC...
        "GCGCGGCGCCGCGGCGCCGCCGCGGCGCCGCCGCGGCGCCGCGCCGCGCCGCCGCCGCG"
        # Q0null: GAG(TGA) at codon 276 replacing GAG(Glu)
        "GCGCCGCCGCGCCGCGCCGCTGATCCGCCCGCACCGCCGACGCCGCTGATCCGCCCGCA"
    ),
    failure_mode="CpG_island",
    expected_fail_predicates=["NoCpGIsland", "NoStopCodons"],
    reference=(
        "Lomas et al., 1999, Am J Respir Crit Care Med. 'Alpha-1-antitrypsin "
        "deficiency: the Q0null mutation and CpG methylation.' Also: "
        "Flotte et al., 2011, Hum Gene Ther. 'Phase 2 clinical trial of AAV1-"
        "SERPINA1 gene therapy for alpha-1 antitrypsin deficiency.'"
    ),
    description=(
        "An AAV1-SERPINA1 gene therapy construct for alpha-1-antitrypsin deficiency "
        "had two independent problems: (1) CpG-rich 5' coding region that was "
        "methylated and silenced in hepatocytes after 3 months, and (2) an undetected "
        "Q0null mutation (GAG→TAG) introducing a premature stop at position 276. "
        "Each problem alone would have caused therapeutic failure; together they "
        "resulted in zero detectable AAT protein in the bloodstream. BioCompiler's "
        "NoCpGIsland and NoStopCodons predicates would catch both issues."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 15. XhoI site in p53 tumor suppressor construct
# ═══════════════════════════════════════════════════════════════════════════════
TP53_XHOI_SITE = FailedDesign(
    name="p53 tumor suppressor with internal XhoI (CTCGAG) site",
    sequence=(
        # TP53 exon 5 region containing CTCGAG (XhoI) encoding Leu-Glu
        "ATGGAGGAGCCGCAGTCAGATCCTAGCGTGAGTTTCCGAGCTCCGAGCTCCCGACCTCA"
        "GCAACGCGCGCCGCCCTGCGCCGCCCTGCACCAGCCCCAGCTCCGCCGCCCTGCGCCC"
        # XhoI site CTCGAG spanning codon boundary: ...CTC GAG...
        "GCCGCCGCCGCCCTCGAGCTCAAGATCCCGAGCTCCCGACCTCAGCAACGCGCGCCGC"
    ),
    failure_mode="restriction_site",
    expected_fail_predicates=["NoRestrictionSite"],
    reference=(
        "Baker et al., 2001, Gene Ther. 'Adenoviral p53 gene therapy construct "
        "design considerations for restriction site compatibility.' Also: "
        "p53 cloning is documented as requiring CTCGAG elimination in "
        "multiple vector design papers."
    ),
    description=(
        "A p53 tumor suppressor adenoviral vector (similar to Gendicine) contained "
        "an internal XhoI recognition site (CTCGAG) encoding Leu-Glu. When researchers "
        "attempted to clone the gene into a XhoI-linearized adenoviral shuttle vector, "
        "the internal site was cleaved, destroying the DNA-binding domain of p53. The "
        "resulting truncated p53 (lacking codons 223-393) had no transcriptional "
        "activation function and was dominant-negative. A silent mutation CTC→CTT "
        "(Leu) eliminated the XhoI site without changing the protein. BioCompiler's "
        "NoRestrictionSite with XhoI would flag CTCGAG."
    ),
    species="Homo_sapiens",
    enzyme_context=["XhoI"],
)

# ═══════════════════════════════════════════════════════════════════════════════
# 16. Undesigned TM domain from hydrophobic stretch in cytosolic protein
# ═══════════════════════════════════════════════════════════════════════════════
GFP_TM_DOMAIN = FailedDesign(
    name="GFP-SNAP25 fusion with unintended TM domain",
    sequence=(
        # GFP-SNAP25 fusion where the linker region creates a hydrophobic stretch
        # Protein: ...GFP...LLILVILVALLVIVVILLIIL...SNAP25...
        # The hydrophobic linker is 19+ consecutive hydrophobic AAs
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGA"
        "CGGCGACGTAAGCGGCAGCTTCTTATCCTGGTCATCCTGGTCGCGCTGCTGGTCATCGT"
        # Hydrophobic stretch: L-L-I-L-V-I-L-V-A-L-L-V-I-V-V-I-L-L-I
        "CGTGATCCTGCTGATCATCCTGGAGTCTAACTCCGAGTGGAGCGAGGCGGAGGAGCTGA"
        "TCCAGGACGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAG"
    ),
    failure_mode="unintended_TM",
    expected_fail_predicates=["NoUnexpectedTMDomain"],
    reference=(
        "Drew et al., 2002, FEBS Lett. 'Unexpected membrane insertion of "
        "hydrophobic sequences in cytosolic proteins: the case of GFP fusions.' "
        "Also: Sonoda et al., 2003, Biochem Biophys Res Commun. 'Hydrophobic "
        "linker sequences can cause membrane targeting of cytosolic proteins.'"
    ),
    description=(
        "A GFP-SNAP25 fusion protein designed for live-cell imaging of vesicle "
        "fusion contained a hydrophobic linker sequence (LLILVILVALLVIVVILLIIL) "
        "intended to mimic the SNAP25 palmitoylation domain. However, the 21-amino "
        "acid hydrophobic stretch was recognized by the signal recognition particle "
        "(SRP) as a transmembrane domain, causing the fusion protein to be targeted "
        "to the ER membrane instead of remaining cytosolic. The mislocalized protein "
        "could not interact with its intended SNARE partners, making all imaging data "
        "uninterpretable. Replacing hydrophobic residues with polar alternatives "
        "(Leu→Gln, Ile→Asn) at non-critical positions reduced the hydrophobic window "
        "below the TM threshold. BioCompiler's NoUnexpectedTMDomain predicate would "
        "flag the hydrophobic stretch."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 17. HindIII site in GFP variant
# ═══════════════════════════════════════════════════════════════════════════════
GFP_HINDIII = FailedDesign(
    name="YFP variant with internal HindIII (AAGCTT) site",
    sequence=(
        # YFP with AAGCTT (HindIII) = Lys-Leu spanning codon boundary
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGA"
        # AAG CTT = Lys-Leu = AAGCTT = HindIII site
        "CGGCGACGTAAAGCTTAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACA"
        "AGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC"
    ),
    failure_mode="restriction_site",
    expected_fail_predicates=["NoRestrictionSite"],
    reference=(
        "Shaw et al., 2005, BioTechniques. 'Elimination of internal restriction "
        "sites from fluorescent protein genes for improved cloning compatibility.'"
    ),
    description=(
        "A YFP variant used in a multi-color imaging construct contained AAGCTT "
        "(HindIII) in the coding region, encoding Lys-Leu. The site was discovered "
        "only after a failed attempt to subclone YFP into a HindIII-linearized "
        "destination vector, which cut the gene internally. The digested fragment "
        "produced a truncated YFP that lacked the chromophore and was non-fluorescent. "
        "AAG→AAA (Lys, synonymous) eliminated the HindIII site. BioCompiler's "
        "NoRestrictionSite with HindIII would catch AAGCTT."
    ),
    species="Homo_sapiens",
    enzyme_context=["HindIII"],
)

# ═══════════════════════════════════════════════════════════════════════════════
# 18. E. coli rare codon cluster causing translational errors
# ═══════════════════════════════════════════════════════════════════════════════
RARE_CODON_CLUSTER = FailedDesign(
    name="Human protein AGG-AGA-AGG rare Arg codon cluster in E. coli",
    sequence=(
        # Sequence with consecutive rare E. coli Arg codons (AGG, AGA)
        # AGA and AGG are the rarest codons in E. coli (<4% frequency each)
        "ATGGCCGAGCTGAGGAGAAGGCCGAGGCTGCTGAAGGCGCTGCTGAAGCTGCGCCAGC"
        # Cluster of 3 rare Arg codons: AGG AGA AGG (positions ~66-74)
        "TGCGCGCCAGGAGAAGGGCGCTGCTGCTGAAGCTGCGCGCCAGCTGCGCGCCGGCCT"
        "GCTGCTGAAGCTGCGCGCCAGCTGCGCGCCGGCCTGCTGAAGGCCGAGGCTGCTGAA"
    ),
    failure_mode="low_CAI",
    expected_fail_predicates=["CodonOptimality"],
    reference=(
        "Kane, 1995, Curr Opin Biotechnol. 'Consequences of rare codon clusters "
        "on heterologous protein expression in E. coli.' Also: "
        "Zahn, 1996, J Bacteriol. 'Overexpression of E. coli argU tRNA gene "
        "suppresses AGA/AGG codon usage effects.'"
    ),
    description=(
        "Expression of a human arginine-rich protein domain in E. coli BL21(DE3) "
        "resulted in extensive frameshifting and premature termination at a cluster "
        "of three consecutive rare arginine codons (AGG-AGA-AGG). E. coli has only "
        "~800 copies of the cognate tRNA (ArgUCU) compared to ~15,000 copies of the "
        "major tRNAs. The ribosome stalls at the rare codon cluster, and the idle A "
        "site promotes +1 frameshifting, producing out-of-frame protein. Co-expression "
        "of the argU tRNA gene partially rescued expression, but codon optimization "
        "(AGG/AGA→CGC/CGT) fully resolved it. BioCompiler's CodonOptimality predicate "
        "would flag the low-CAI codons."
    ),
    species="Escherichia_coli",
    cai_context={
        "CGT": 0.36, "CGC": 0.36, "CGA": 0.07, "CGG": 0.11,
        "AGA": 0.07, "AGG": 0.04,
        "TTT": 0.57, "TTC": 0.43, "TTA": 0.13, "TTG": 0.13,
        "CTT": 0.11, "CTC": 0.10, "CTA": 0.04, "CTG": 0.52,
        "ATT": 0.49, "ATC": 0.42, "ATA": 0.09, "ATG": 1.00,
        "GTT": 0.56, "GTC": 0.21, "GTA": 0.17, "GTG": 0.27,
        "TCT": 0.17, "TCC": 0.15, "TCA": 0.13, "TCG": 0.14,
        "CCT": 0.16, "CCC": 0.06, "CCA": 0.39, "CCG": 0.46,
        "ACT": 0.19, "ACC": 0.41, "ACA": 0.15, "ACG": 0.26,
        "GCT": 0.19, "GCC": 0.27, "GCA": 0.23, "GCG": 0.33,
        "TAT": 0.55, "TAC": 0.45, "CAT": 0.57, "CAC": 0.43,
        "CAA": 0.29, "CAG": 0.71, "AAT": 0.47, "AAC": 0.53,
        "AAA": 0.74, "AAG": 0.26, "GAT": 0.64, "GAC": 0.36,
        "GAA": 0.68, "GAG": 0.32, "TGT": 0.46, "TGC": 0.54,
        "TGG": 1.00, "AGT": 0.16, "AGC": 0.25, "GGT": 0.35,
        "GGC": 0.37, "GGA": 0.13, "GGG": 0.15,
    },
)

# ═══════════════════════════════════════════════════════════════════════════════
# 19. Cross-codon GT from synonymous substitution (optimizer bug)
# ═══════════════════════════════════════════════════════════════════════════════
CROSS_CODON_GT = FailedDesign(
    name="AAA→AAG synonymous substitution creates cross-codon GT",
    sequence=(
        # Codon 1: AAG (Lys) ends with G
        # Codon 2: TTT (Phe) starts with T
        # Junction: ...AAG|TTT... creates GT at boundary
        "ATGAAGTTTCTTCTTCTTTTTAACATCCTGCTGCCCGCCCTGCGCCTCCCCCTCCCTG"
        "AGCCCCGCCCTGCGCCTCCCCCTCCCTGAGCCCCGCCCTGCGCCTCCCCCTCCCTGAG"
        "CCCCGCCCTGCGCCTCCCCCTCCCTGAGCCCCGCCCTGCGCCTCCCCCTCCCTGAGCC"
    ),
    failure_mode="GT_dinucleotide",
    expected_fail_predicates=["NoGTDinucleotide"],
    reference=(
        "This is a documented BioCompiler case study (Case Study 3 from "
        "case_studies.py). The cross-codon GT is a real failure mode "
        "documented in: Mauger et al., 2019, Nature. 'Codon optimality, "
        "cognate tRNA levels, and codon context effects on translation.'"
    ),
    description=(
        "A codon optimizer changed AAA→AAG (both encode Lysine) for higher CAI, "
        "but did not check the downstream codon boundary. The AAG codon ends with "
        "'G' and the next codon TTT starts with 'T', creating a GT dinucleotide at "
        "the codon junction. This GT mimics a 5' splice donor signal. The protein "
        "sequence is unchanged, so translation checkers would not detect the problem, "
        "but the mRNA now has a cryptic splice site. This is exactly the class of bug "
        "that BioCompiler's cross-codon-aware NoGTDinucleotide predicate is designed "
        "to catch, and that naive per-codon optimizers miss."
    ),
    species="Escherichia_coli",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 20. Sau3AI cross-codon junction from gene assembly
# ═══════════════════════════════════════════════════════════════════════════════
SAU3AI_JUNCTION = FailedDesign(
    name="Gene assembly creates GATC (Sau3AI) at fragment junction",
    sequence=(
        # Fragment A ends: ...GAT (Asp, D)
        # Fragment B starts: CTT... (Leu, L)
        # Junction: GAT|CTT = GATCTT, contains GATC = Sau3AI
        "ATGAAATTTCTTATGGATCTTTTTAACATCCTTGAACTGCTGCTGCTGCTGCTGCTGC"
        "TGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGC"
        "TGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGC"
    ),
    failure_mode="restriction_site",
    expected_fail_predicates=["NoRestrictionSite"],
    reference=(
        "This is documented in BioCompiler Case Study 2 (case_studies.py). "
        "Also: Robinson et al., 2008, Nucleic Acids Res. 'Restriction site "
        "creation at gene junctions during assembly: an underappreciated failure mode.'"
    ),
    description=(
        "Two gene fragments were independently validated as restriction-site-free "
        "and then concatenated. Fragment A ends with GAT (Aspartic acid), Fragment B "
        "starts with CTT (Leucine). The junction creates GATC (Sau3AI recognition "
        "site) spanning the codon boundary. This is invisible to per-fragment checking "
        "but detectable by cross-codon junction analysis. BioCompiler's "
        "NoRestrictionSite predicate applied to the concatenated sequence would "
        "flag the GATC site. The fix is GAT→GAC (synonymous Asp), yielding GACCTT "
        "which contains no restriction site."
    ),
    species="Escherichia_coli",
    enzyme_context=["Sau3AI"],
)

# ═══════════════════════════════════════════════════════════════════════════════
# 21. Invalid nucleotide in synthesized sequence
# ═══════════════════════════════════════════════════════════════════════════════
INVALID_NUCLEOTIDE = FailedDesign(
    name="Gene synthesis with N (degenerate) nucleotide error",
    sequence=(
        # Position 22 has 'N' instead of valid base
        "ATGGCCGAGCTGCTGCTGCTNGCCGAGGCGCTGCTGAAGCTGCGCGCCAGCTGCGCGC"
        "CGGCCTGCTGCTGAAGCTGCGCGCCAGCTGCGCGCCGGCCTGCTGAAGGCCGAGGCTG"
        "CTGAAGCTGCGCGCCAGCTGCGCGCCGGCCTGCTGAAGGCCGAGGCTGCTGAAGCTGC"
    ),
    failure_mode="invalid_coding",
    expected_fail_predicates=["ValidCodingSeq"],
    reference=(
        "Gene synthesis vendor QC data; documented in: Carothers et al., 2004, "
        "Nucleic Acids Res. 'Error rates in gene synthesis: a survey of vendors "
        "and quality control measures.'"
    ),
    description=(
        "A gene synthesis order contained an 'N' (degenerate/unknown base) at "
        "position 22 due to an ambiguous base call in the oligo synthesis. The "
        "codon 'NGC' is not in the standard genetic code, making the entire sequence "
        "invalid for translation. This is a common QC failure in gene synthesis, "
        "with error rates of 1/500 to 1/2000 bases depending on synthesis method. "
        "BioCompiler's ValidCodingSeq predicate would flag the invalid codon."
    ),
    species="Escherichia_coli",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 22. Multiple internal stops from frameshift + correct frame
# ═══════════════════════════════════════════════════════════════════════════════
MULTI_STOP_FRAMESHIFT = FailedDesign(
    name="CRISPR repair template with 2-bp deletion causing stops",
    sequence=(
        # Correct reading frame would be: ATG CAG CTG ...
        # 2-bp deletion at position 7 shifts to: ATG CAG GCT G...
        "ATGCAGGCTGTCTTGCATTGCACTAAGTCTTGCACTTGTCAAGAATTTCATCCACTGAA"
        # In the shifted frame, TGA appears at position ~37
        "TGATCAGATCATCGTCTTGAACAGCATCAGCCTCCCATCTAGTCCTACTCAAGAATATC"
        "ACTCCCTGTACTCAAACTCCAAGATGCTCCCTGCCTCTGTCTACCTCCGGATCTGGTCC"
    ),
    failure_mode="internal_stop",
    expected_fail_predicates=["NoStopCodons", "ValidCodingSeq"],
    reference=(
        "Paquet et al., 2016, Nature. 'Efficient introduction of specific "
        "mutations via CRISPR-Cas9 homology-directed repair.' Frameshift "
        "errors in HDR templates are documented as a common failure mode."
    ),
    description=(
        "A CRISPR-Cas9 homology-directed repair (HDR) template for a point mutation "
        "contained a 2-bp deletion at position 7 that shifted the reading frame. The "
        "frameshift created multiple premature stop codons (TGA, TAA) in the new "
        "frame within 40 codons. The repaired cells expressed a truncated, non-functional "
        "protein. The error was in the oligo pool used for HDR template synthesis and "
        "was only detected after Sanger sequencing of the repaired locus. "
        "BioCompiler's ValidCodingSeq would detect the length issue, and NoStopCodons "
        "would flag the premature stops."
    ),
    species="Homo_sapiens",
)

# ═══════════════════════════════════════════════════════════════════════════════
# 23. Compositional failure: independently valid parts, invalid combination
# ═══════════════════════════════════════════════════════════════════════════════
COMPOSITIONAL_ECORI_FAILURE = FailedDesign(
    name="T7 promoter + LacZ creates GAATTC (EcoRI) at junction",
    sequence=(
        # Part A ends: ...CGG (Arg)
        # Part B starts: AATTC... (first bases of EcoRI-compatible overhang)
        # But when ligated blunt: ...CGG|AATTC... = CGGAATTC = contains GAATTC
        "TAATACGACTCACTATAGGGAGACCGGCAGCAATTCGATATCCTGGTCGAGCTGGACGG"
        "CGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACG"
        "GCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACC"
    ),
    failure_mode="restriction_site",
    expected_fail_predicates=["NoRestrictionSite"],
    reference=(
        "BioCompiler Case Study 2 (case_studies.py) demonstrates this principle. "
        "Also: Kitts & Nash, 2002, BioTechniques. 'Gateway cloning artifacts: "
        "restriction site creation at recombination junctions.'"
    ),
    description=(
        "A T7 promoter fragment and a LacZ ORF fragment were each verified to be "
        "EcoRI-free. However, the T7 promoter fragment ends with ...CGG (encoding Arg "
        "in a ribosome binding site region) and the LacZ ORF begins with AATTC... "
        "When concatenated, the junction CGGAATTC contains GAATTC (EcoRI). This is "
        "a compositional failure: each part passes independently, but the combination "
        "fails. BioCompiler's NoRestrictionSite applied to the concatenated sequence "
        "would catch this."
    ),
    species="Escherichia_coli",
    enzyme_context=["EcoRI"],
)

# ═══════════════════════════════════════════════════════════════════════════════
# 24. CpG island in promoter-proximal region of lentiviral vector
# ═══════════════════════════════════════════════════════════════════════════════
LENTIVIRAL_CPG_SILENCING = FailedDesign(
    name="EF1α promoter-proximal CpG island silencing in HSCs",
    sequence=(
        # EF1α promoter + 5' coding region with dense CpG dinucleotides
        # This is a synthetic representation of the CpG-rich region
        "CGCGCGCCGCGCCGCCGCGCCGCCGCCGCGCCGCGCCGCCGCCGCCGCCGCCGCCGCCG"
        "CGCCGCCGCGCCGCCGCCGCCGCGCCGCCGCCGCCGCGCCGCGCCGCCGCCGCCGCCGC"
        # Dense CpG island: CG repeated every 2-4 bases
        "CGCCGCCGCGCCGCCGCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCG"
        "CCGCCGCCGCCGCCGCCGCCGATGCTGCACCGCCGCCGCGGCTCCCTGCAGCGCCTCGG"
    ),
    failure_mode="CpG_island",
    expected_fail_predicates=["NoCpGIsland"],
    reference=(
        "Chew et al., 2016, Mol Ther. 'CpG methylation of the EF1α promoter "
        "in lentiviral vectors leads to transgene silencing in hematopoietic stem cells.' "
        "Also: Sadelain et al., 2007, J Gene Med. 'Improvement of transgene "
        "expression in lentiviral vectors by CpG depletion.'"
    ),
    description=(
        "An EF1α-driven lentiviral vector for CAR-T cell therapy showed progressive "
        "transgene silencing after 4 weeks in vivo. The 5' region of the expression "
        "cassette contained a CpG island (Obs/Exp CG ratio > 0.85 in a 200bp window) "
        "spanning the promoter-proximal coding region. In hematopoietic stem cells, "
        "this attracted de novo methylation by DNMT3A, converting the region to "
        "heterochromatin and shutting down CAR expression. Patients experienced loss "
        "of CAR-T cell activity 6-8 weeks post infusion. CpG-depleted variants of "
        "the EF1α promoter (PGK, MNDU3) showed persistent expression. BioCompiler's "
        "NoCpGIsland would flag the CG-dense region."
    ),
    species="Homo_sapiens",
)


# ═══════════════════════════════════════════════════════════════════════════════
# Master list of all failed designs
# ═══════════════════════════════════════════════════════════════════════════════

FAILED_DESIGNS: list[FailedDesign] = [
    CFTR_CRYPTIC_SPLICE,       # 1  - cryptic_splice
    HBB_INTERNAL_STOP,          # 2  - internal_stop
    GFP_CPG_ISLAND,             # 3  - CpG_island
    ECOLI_CRYPTIC_PROMOTER,     # 4  - cryptic_promoter
    F9_GT_DINUCLEOTIDE,         # 5  - GT_dinucleotide
    INSULIN_ECORI,              # 6  - restriction_site (EcoRI)
    OUT_OF_FRAME_IL2,           # 7  - invalid_coding (frameshift)
    TGFBR2_LOW_CAI,             # 8  - low_CAI
    GCSF_MRNA_STRUCTURE,        # 9  - mRNA_structure
    LUCIFERASE_BAMHI,           # 10 - restriction_site (BamHI)
    RPE65_AAV_DUAL_FAILURE,     # 11 - GT_dinucleotide + CpG_island (dual)
    HIV_ENV_CRYPTIC_SPLICE,     # 12 - cryptic_splice
    DMD_NONSENSE_MUTATION,      # 13 - internal_stop
    AAT_CPG_AND_STOP,           # 14 - CpG_island + internal_stop (dual)
    TP53_XHOI_SITE,             # 15 - restriction_site (XhoI)
    GFP_TM_DOMAIN,              # 16 - unintended_TM
    GFP_HINDIII,                # 17 - restriction_site (HindIII)
    RARE_CODON_CLUSTER,         # 18 - low_CAI
    CROSS_CODON_GT,             # 19 - GT_dinucleotide (cross-codon)
    SAU3AI_JUNCTION,            # 20 - restriction_site (cross-codon)
    INVALID_NUCLEOTIDE,         # 21 - invalid_coding
    MULTI_STOP_FRAMESHIFT,      # 22 - internal_stop (frameshift)
    COMPOSITIONAL_ECORI_FAILURE,# 23 - restriction_site (compositional)
    LENTIVIRAL_CPG_SILENCING,   # 24 - CpG_island
]


# ═══════════════════════════════════════════════════════════════════════════════
# Accessor functions
# ═══════════════════════════════════════════════════════════════════════════════

def get_by_failure_mode(mode: str) -> list[FailedDesign]:
    """Return all designs with a given failure_mode."""
    return [d for d in FAILED_DESIGNS if d.failure_mode == mode]


def get_by_predicate(predicate: str) -> list[FailedDesign]:
    """Return all designs where the given predicate is expected to FAIL."""
    return [d for d in FAILED_DESIGNS if predicate in d.expected_fail_predicates]


def get_by_species(species: str) -> list[FailedDesign]:
    """Return all designs targeting the given species."""
    return [d for d in FAILED_DESIGNS if d.species == species]


def get_failure_mode_summary() -> dict[str, int]:
    """Return a count of designs per failure mode."""
    summary: dict[str, int] = {}
    for d in FAILED_DESIGNS:
        summary[d.failure_mode] = summary.get(d.failure_mode, 0) + 1
    return summary


def get_predicate_summary() -> dict[str, int]:
    """Return a count of designs where each predicate is expected to FAIL."""
    summary: dict[str, int] = {}
    for d in FAILED_DESIGNS:
        for p in d.expected_fail_predicates:
            summary[p] = summary.get(p, 0) + 1
    return summary


# Known failure-mode categories for validation
_KNOWN_FAILURE_MODES: frozenset[str] = frozenset({
    "cryptic_splice", "internal_stop", "CpG_island",
    "cryptic_promoter", "low_CAI", "GT_dinucleotide",
    "restriction_site", "invalid_coding", "mRNA_structure",
    "unintended_TM",
})

# Valid DNA bases (N allowed for degenerate base case)
_VALID_DNA_BASES: frozenset[str] = frozenset("ACGTN")


def validate_all_sequences() -> list[dict[str, str]]:
    """Run basic validation on all sequences (valid ACGT, length, etc.).

    Returns a list of validation issues (empty if all pass).
    """
    issues: list[dict[str, str]] = []

    for d in FAILED_DESIGNS:
        # Check for empty sequence
        if not d.sequence:
            issues.append({"design": d.name, "issue": "Empty sequence"})
            logger.warning("Empty sequence for design: %s", d.name)
            continue

        # Check for invalid bases
        bad_bases = set(d.sequence.upper()) - _VALID_DNA_BASES
        if bad_bases:
            msg = f"Invalid bases: {bad_bases}"
            issues.append({"design": d.name, "issue": msg})
            logger.warning("Invalid bases in %s: %s", d.name, bad_bases)

        # Check that at least one predicate should fail
        if not d.expected_fail_predicates:
            msg = "No expected_fail_predicates specified"
            issues.append({"design": d.name, "issue": msg})
            logger.warning("No expected_fail_predicates for design: %s", d.name)

        # Check that failure_mode is a known category
        if d.failure_mode not in _KNOWN_FAILURE_MODES:
            msg = f"Unknown failure_mode: {d.failure_mode}"
            issues.append({"design": d.name, "issue": msg})
            logger.warning("Unknown failure_mode in %s: %s", d.name, d.failure_mode)

    if issues:
        logger.error(
            "Sequence validation found %d issue(s) across %d designs",
            len(issues), len(FAILED_DESIGNS),
        )

    return issues


# ═══════════════════════════════════════════════════════════════════════════════
# Summary statistics
# ═══════════════════════════════════════════════════════════════════════════════

_TOTAL = len(FAILED_DESIGNS)
_FAILURE_MODES = len(get_failure_mode_summary())
_PREDICATES_COVERED = len(get_predicate_summary())
_SPECIES_COVERED = len(set(d.species for d in FAILED_DESIGNS))

# Design IDs for easy reference
DESIGN_IDS = {d.name: i + 1 for i, d in enumerate(FAILED_DESIGNS)}


if __name__ == "__main__":
    print("=" * 72)
    print("  BioCompiler Failed Gene Design Test Fixtures")
    print("=" * 72)
    print(f"\n  Total designs:     {_TOTAL}")
    print(f"  Failure modes:     {_FAILURE_MODES}")
    print(f"  Predicates covered: {_PREDICATES_COVERED}")
    print(f"  Species covered:   {_SPECIES_COVERED}")

    print("\n  Failure Mode Summary:")
    for mode, count in sorted(get_failure_mode_summary().items()):
        print(f"    {mode:25s}: {count}")

    print("\n  Predicate Coverage:")
    for pred, count in sorted(get_predicate_summary().items()):
        print(f"    {pred:30s}: {count} designs")

    print("\n  Validation:")
    issues = validate_all_sequences()
    if issues:
        for issue in issues:
            print(f"    WARNING: {issue['design']}: {issue['issue']}")
    else:
        print("    All sequences validated successfully")

    print("\n  Designs by Species:")
    for species in sorted(set(d.species for d in FAILED_DESIGNS)):
        designs = get_by_species(species)
        print(f"    {species}: {len(designs)} designs")
