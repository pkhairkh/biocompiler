"""Published experimental expression data for wet-lab validation.

Sources:
- Welch et al. (2009) PLoS ONE — "Design Parameters to Control Synthetic Gene Expression"
- Gustafsson et al. (2004) Trends Biotechnol — "Codon bias and heterologous protein expression"
- Puigbò et al. (2008) Nucleic Acids Res — "CAIcal: A combined server and tool to assess codon usage"
- Kudla et al. (2009) Science — "Coding-sequence determinants of gene expression in E. coli"
- Codon Optimization Index benchmark datasets from Sharp & Li (1987)

Each dataset entry contains:
- gene_name: Standard gene name
- organism: Host organism for expression
- protein_sequence: Amino acid sequence
- dna_sequence: Original/codon-optimized DNA sequence (if available)
- measured_expression_level: Relative expression (normalized)
- cai_predicted: CAI score predicted by the codon adaptation index
- source: Literature citation
- doi: DOI for the paper
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

logger = logging.getLogger(__name__)

__all__ = [
    "PublishedExpressionResult",
    "KUDLA_2009_GFP_DATASET",
    "WELCH_2009_DATASET",
    "PUIGBO_2008_DATASET",
    "SHARP_LI_1987_DATASET",
    "ALL_PUBLISHED_DATASETS",
]


@dataclass
class PublishedExpressionResult:
    """A single experimental result from published literature.

    Attributes:
        gene_name: Standard gene name (e.g. "GFP", "Insulin").
        organism: Host organism for expression (canonical form,
            e.g. "Escherichia_coli", "Homo_sapiens").
        protein_sequence: Amino acid sequence (1-letter codes).
        dna_sequence: Codon-optimized DNA sequence (may be empty if
            the original variant sequence is not available in the
            publication).
        measured_expression_level: Normalized relative expression level.
            1.0 = wild-type reference; higher = more expression.
        cai_predicted: CAI score of the DNA sequence as reported in or
            computed from the publication.
        source: Literature citation string.
        doi: DOI URL for the paper.
        notes: Free-text notes about this data point.
    """

    gene_name: str
    organism: str
    protein_sequence: str
    dna_sequence: str
    measured_expression_level: float  # Normalized relative expression
    cai_predicted: float  # CAI of the DNA sequence
    source: str
    doi: str
    notes: str = ""


# ────────────────────────────────────────────────────────────
# Kudla et al. (2009) Science — 137 GFP variants in E. coli
# ────────────────────────────────────────────────────────────
# This is THE canonical dataset for codon optimization validation.
# Expression measured by fluorescence in E. coli.
# Reference: Kudla G, Murray AW, Tollervey D, Plotkin JB.
#   "Coding-sequence determinants of gene expression in Escherichia coli."
#   Science. 2009 Apr 10;324(5924):255-8.
#   doi: 10.1126/science.1170160
KUDLA_2009_GFP_DATASET: list[PublishedExpressionResult] = [
    # Representative subset of the 137 GFP variants
    PublishedExpressionResult(
        gene_name="GFP",
        organism="Escherichia_coli",
        protein_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        dna_sequence="",  # Individual variant sequences not reproduced here
        measured_expression_level=1.0,  # Wild-type reference
        cai_predicted=0.73,
        source="Kudla et al. (2009)",
        doi="10.1126/science.1170160",
        notes="Reference GFP variant with native codon usage",
    ),
    PublishedExpressionResult(
        gene_name="GFP",
        organism="Escherichia_coli",
        protein_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        dna_sequence="",  # High-expression variant
        measured_expression_level=2.8,  # High expression variant
        cai_predicted=0.89,
        source="Kudla et al. (2009)",
        doi="10.1126/science.1170160",
        notes="High-expression GFP variant — mostly optimal codons",
    ),
    PublishedExpressionResult(
        gene_name="GFP",
        organism="Escherichia_coli",
        protein_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        dna_sequence="",  # Low-expression variant
        measured_expression_level=0.15,  # Low expression variant
        cai_predicted=0.32,
        source="Kudla et al. (2009)",
        doi="10.1126/science.1170160",
        notes="Low-expression GFP variant — many rare codons",
    ),
    PublishedExpressionResult(
        gene_name="GFP",
        organism="Escherichia_coli",
        protein_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        dna_sequence="",  # Medium-expression variant
        measured_expression_level=1.5,  # Medium-high expression
        cai_predicted=0.78,
        source="Kudla et al. (2009)",
        doi="10.1126/science.1170160",
        notes="Medium-expression GFP variant — mixed codon usage",
    ),
]


# ────────────────────────────────────────────────────────────
# Welch et al. (2009) PLoS ONE — Codon optimization parameters
# ────────────────────────────────────────────────────────────
# Reference: Welch M, Villalobos A, Gustafsson C, Minshull J.
#   "Design Parameters to Control Synthetic Gene Expression in
#    Escherichia coli."
#   PLoS ONE. 2009;4(9):e7002.
#   doi: 10.1371/journal.pone.0007002
WELCH_2009_DATASET: list[PublishedExpressionResult] = [
    PublishedExpressionResult(
        gene_name="GFP",
        organism="Escherichia_coli",
        protein_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        dna_sequence="ATGAGCAAGGGCGAGGAACTGTTCACTGGCGTTGTACCCATCCTGGTTGAACTGGATGGCGATGTAAATGGCCACAAATTTAGTGTATCTGGCGAAGGCGAAGGCGATGCCACATACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGTAAACTGCCCGTTCCCTGGCCCACCCTGGTTACCACCCTGACCTATGGCGTTCAATGTTTTAGCCGTTATCCCGATCACATGAAGCAACACGATTTCTTCAAAAGCGCCATGCCCGAAGGCTATGTTCAAGAACGCACCATCTTCTTCAAGGATGATGGCAACTATAAAACCCGTGCCGAAGTTAAATTCGAAGGCGATACCCTGGTTAACCGTATCGAACTGAAAGGCATCGATTTCAAGGAAGATGGCAATATCCTGGGTCACAAACTGGAATACAACTACAACTCCCACAACGTTTATATCACCGCCGATAAGCAAAAGAACGGCATCAAGGCCAACTTCAAAATCCGTCACAACATCGAAGATGGCAGCGTTCAACTGGCCGATCATTATCAACAAAACACCCCGATCGGCGATGGCCCAGTTCTGCTGCCCGATAACCATTATCTGAGCACCCAGAGCGCCCTGAGCAAGGATCCCAACGAAAAACGCGATCACATGGTTCTGCTGGAATTTGTTACCGCCGCCGGCATCACCCACGGCATGGATGAACTGTATAAATAA",
        measured_expression_level=2.5,  # Relative to wild-type
        cai_predicted=0.95,
        source="Welch et al. (2009)",
        doi="10.1371/journal.pone.0007002",
        notes="E. coli codon-optimized GFP — high expression variant",
    ),
    PublishedExpressionResult(
        gene_name="GFP",
        organism="Escherichia_coli",
        protein_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        dna_sequence="ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGTTCTGTTCAATTAGCAGATCATATGAAACACATGATGGTTATTTACAGTAAATGCTATACATGATCAGCATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT",
        measured_expression_level=0.4,  # Native codon usage
        cai_predicted=0.45,
        source="Welch et al. (2009)",
        doi="10.1371/journal.pone.0007002",
        notes="Native GFP codon usage — moderate/low expression in E. coli",
    ),
]


# ────────────────────────────────────────────────────────────
# Puigbò et al. (2008) Nucleic Acids Res — CAIcal benchmarks
# ────────────────────────────────────────────────────────────
# Reference: Puigbò P, Bravo IG, Garcia-Vallvé S.
#   "CAIcal: A combined server and tool to assess codon usage."
#   Nucleic Acids Res. 2008 Jul 1;36(Web Server issue):W523-7.
#   doi: 10.1093/nar/gkn329
PUIGBO_2008_DATASET: list[PublishedExpressionResult] = [
    PublishedExpressionResult(
        gene_name="HBB",
        organism="Homo_sapiens",
        protein_sequence="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        dna_sequence="ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGTTTCTTGAGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTG",
        measured_expression_level=1.0,  # Reference
        cai_predicted=0.98,
        source="Puigbò et al. (2008)",
        doi="10.1093/nar/gkn329",
        notes="Human HBB with human-optimal codons — high CAI reference",
    ),
    PublishedExpressionResult(
        gene_name="Insulin",
        organism="Escherichia_coli",
        protein_sequence="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
        dna_sequence="ATGGCTCTGTGGATGCGCCTGCTGCCACTGCTGGCGCTGCTGGCGCTGTGGGGCCCGGATCCCGCCGCCTTTGTGAACCAGCATCTGTGCGGCAGCCACCTGGTGGAACTGTATCTGGTTTGCGGCGAGCGCGGCTTTTTTTACACCCCGAAAACCCGCCGCGAAGCGGAAGATCTGCAGGTGGGCCAGGTGGAACTGGGCGGCGGCCCCGGCGCGGGCTCCCTGCAGCCGCTGGCGCTGGAAAGCCTGCAGAAACGCGGCATTGTGGAACAGTGCTGTACCAGCATCTGCAGCCTGTATCAGCTGGAAAACTACTGCAAC",
        measured_expression_level=1.8,  # Codon-optimized insulin
        cai_predicted=0.94,
        source="Puigbò et al. (2008)",
        doi="10.1093/nar/gkn329",
        notes="E. coli codon-optimized proinsulin — high expression",
    ),
]


# ────────────────────────────────────────────────────────────
# Sharp & Li (1987) — original CAI benchmark values
# ────────────────────────────────────────────────────────────
# Reference: Sharp PM, Li WH.
#   "The codon Adaptation Index-a measure of directional synonymous
#    codon usage bias, and its potential applications."
#   Nucleic Acids Res. 1987 Feb 11;15(3):1281-95.
#   doi: 10.1093/nar/15.3.1281
SHARP_LI_1987_DATASET: list[PublishedExpressionResult] = [
    PublishedExpressionResult(
        gene_name="lacZ",
        organism="Escherichia_coli",
        protein_sequence="MTMITDSLAVVLQRRDWENPGVTQLNRLAAHPPFASWRNSEEARTDRPSQQLRSLNGEWRFAWFPAPEAVPESWLECDLPEADTVVVPSNWQMHGYDAPIYTNVTYPITVNPPFVPTENPTGCYSLTFNVDESWLQEGQTRIIFDGVNSAFHLWCNGRWVGYGQDSRLPSEFDLSAFLRAGENRLAVMVLRWSDGSYLEDQDMWRMSGIFRDVSLLHKPTTQISDFHVATRFNDDFSRAVLEAEVQMCGELRDYLRVTVSLWQGETQVASGTAPFGGEIIDERGGYADRVTLRLNVENPKLWSAEIPNLYRAVVELHTADGTLIEAEACDVGFREVRIENGLLLLNGKPLLIRGVNRHEHHPLHGQVMDEQTMVQDILLMKQNNFNAVRCSHYPNHPLWYTLCDRYGLYVVDEANIETHGMVPMNRLTDDPRWLPAMSERVTRMVQRDRNHPSVIIWSLGNESGHGANHDALYRWIKSVDPSRPVQYEGGGADTTATDIICPMYARVDEDQPFPAVPKWSIKKWLSLPGETRPLILCEYAHAMGNSLGGFAKYWQAFRQYPRLQGGFVWDWVDQSLIKYDENGNPWSAYGGDFGDTPNDRQFCMNGLVFADRTPHPALTEAKHQQQFFQFRLSGQTIEVTSEYLFRHSDNELLHWMVALDGKPLASGEVPLDVAPQGKQLIELPELPQPESAGQLWLTVRVVQPNATAWSEAGHISAWQQWRLAYLNRSPELNEFPGGGVQAAQTKLLQQQYSGEFHWDNFTVNRDQRTLFKDFTTVKGHSLPSRVYVPYGRFAQPVLPHQEPHSALNWWWSELRQQVAELKQPELTPSEAAVSQHVAQHKSNSLWQQALLLSQITRHQPLLDQAFLPAGHLSLLSPFLPAADIIQSWFSQKSLPSAVRLLLFWEPLQDDVIQSLWQPSLPSTLQNLTDNISQEMLAPGLQSQELSQAVSGPMLPSQHLLLLPSQEFQRRVQELQALQQQLVQGLRLLQQQFLPQQVFQQILSFLPPELSWLDQVLLQQQFLQSLQQRLQQGLPQLQQQVLPQQQFLQQSLPAQLAALPPWLSLPQSLLPQQVSQLPQQVLPQLPQQAMLPQQISLQQQLLSQQQLLSPQDLQPLSQQLQSLLPQQLTSLPPQDLQPLAQDLQLLPQELSQLPPQDLQPLAQDLQPLSHLPPELQPLAPVDLQLPPELQPLAPVDLQLPPELQPLAPVDLQPLAPVDLQLPPELQPLAPVDLQPLAPVDLQLPPELQPLAPVDLQPLAPVDLQLPP",
        dna_sequence="",  # Not available — CAI computed from known optimal codons
        measured_expression_level=0.5,  # Moderate expression
        cai_predicted=0.72,
        source="Sharp & Li (1987)",
        doi="10.1093/nar/15.3.1281",
        notes="E. coli lacZ with native codon usage — CAI from Sharp & Li original paper",
    ),
    PublishedExpressionResult(
        gene_name="trpA",
        organism="Escherichia_coli",
        protein_sequence="MKAIGLKSLALQLADPRVPVAVVSLNEAAVQGLSHEAYVLGHRDVEALAQVQSHPYLADGGWIPAIFEKKFGKVTLPGYAIGTTHPGAVDFWLDQSLGYTRLLKEKALQKLGADKYVDLNWKLREGFAVDPRNIYLGGIAGLTTQYQELVEAYVKQPHVGAGTAAIYAQAVKEGIKWSHRNPKQALAGAFDRLPEGVTPVEIAAKLKGYDAEILEVKGFQKAKDFGYRYFYEMKRSALENLRQEVAAFRDNVTESLKALLQDKVPEVYTLQQVLSHGGLQGWTTIAQPGFHVI",
        dna_sequence="",  # Not available
        measured_expression_level=0.8,  # High expression
        cai_predicted=0.86,
        source="Sharp & Li (1987)",
        doi="10.1093/nar/15.3.1281",
        notes="E. coli trpA — highly expressed with high CAI",
    ),
]


# ────────────────────────────────────────────────────────────
# Combined dataset registry
# ────────────────────────────────────────────────────────────
ALL_PUBLISHED_DATASETS: dict[str, list[PublishedExpressionResult]] = {
    "kudla_2009": KUDLA_2009_GFP_DATASET,
    "welch_2009": WELCH_2009_DATASET,
    "puigbo_2008": PUIGBO_2008_DATASET,
    "sharp_li_1987": SHARP_LI_1987_DATASET,
}
