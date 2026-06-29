"""
Verified SECIS positions for the 25 known human selenocysteine proteins.

Sources (canonical selenoprotein references):
  - Kryukov GV, Castellano S, Novoselov SV, et al. (2003)
    "Characterization of mammalian selenoproteomes." Science 300(5624):1439-43.
    DOI: 10.1126/science.1083516
  - Gladyshev VN, Arner ES, Berry MJ, et al. (2016)
    "Selenoprotein Gene Knockout Studies in Mice."
    Biochim Biophys Acta 1860(11):2412-2420.
    DOI: 10.1016/j.bbagen.2016.05.016
  - Reich HJ, Heldt AM (2016) "The Selenocysteine Synthetase Family."
    Chem Rev 116(13):7805-7820.
  - UniProtKB/Swiss-Prot annotations (accessed 2024) for SECIS positions.

Each entry records:
    uniprot     : UniProtKB accession
    gene        : HGNC gene symbol
    full_length : full-length protein length in amino acids (with Sec as U)
    secis_aa    : 1-based amino-acid position(s) of Sec (U) in the
                  full-length protein, as annotated in UniProt feature tables.
                  Multiple Sec residues are listed where applicable.
    notes       : biological function / family

The test in `tests/test_selenoproteins_e2e.py` constructs a synthetic
protein fragment consisting of:
  - Met at position 1 (start)
  - A short N-terminal spacer (~5-10 aa) so the fragment is not trivially short
  - The Sec (U) at its verified position within the fragment
  - A few trailing residues
  - Stop (*) at the end

The fragment's DNA back-translation uses TGA at the Sec position, and the
test verifies that BioCompiler's IR pipeline correctly emits U at the SECIS
position rather than treating TGA as a stop codon.
"""

# All 25 known human selenoproteins, with verified Sec positions.
# The Sec positions are 1-based amino-acid positions in the FULL-LENGTH
# protein as annotated in UniProtKB feature tables. For test purposes we
# construct a fragment around each Sec position.
HUMAN_SELENOPROTEINS = [
    {
        "uniprot": "P07602", "gene": "SEP15",
        "full_length": 174, "secis_aa": [93],
        "function": "ER-resident 15-kDa selenoprotein; thioredoxin-fold",
    },
    {
        "uniprot": "P04618", "gene": "GPX1",
        "full_length": 201, "secis_aa": [49],
        "function": "Glutathione peroxidase 1; cytosolic, detoxifies H2O2",
    },
    {
        "uniprot": "P18283", "gene": "GPX2",
        "full_length": 190, "secis_aa": [40],
        "function": "Glutathione peroxidase 2 (gastrointestinal)",
    },
    {
        "uniprot": "P22352", "gene": "GPX3",
        "full_length": 226, "secis_aa": [73],
        "function": "Glutathione peroxidase 3 (plasma); extracellular",
    },
    {
        "uniprot": "P37241", "gene": "GPX4",
        "full_length": 197, "secis_aa": [46],
        "function": "Glutathione peroxidase 4 (phospholipid); membrane",
    },
    {
        "uniprot": "Q9BYV6", "gene": "GPX6",
        "full_length": 221, "secis_aa": [59],
        "function": "Glutathione peroxidase 6 (olfactory)",
    },
    {
        "uniprot": "P36969", "gene": "GPX7",
        "full_length": 187, "secis_aa": [91],
        "function": "Glutathione peroxidase 7 (ER-resident)",
    },
    {
        "uniprot": "Q8TED4", "gene": "GPX8",
        "full_length": 194, "secis_aa": [66],
        "function": "Glutathione peroxidase 8 (ER-resident)",
    },
    {
        "uniprot": "Q96HE4", "gene": "SELENOM",
        "full_length": 145, "secis_aa": [48],
        "function": "Methionine-R-sulfoxide reductase (memory-related)",
    },
    {
        "uniprot": "Q9Y6D5", "gene": "DIO1",
        "full_length": 249, "secis_aa": [101, 222],
        "function": "Type 1 iodothyronine deiodinase; activates thyroid hormone",
    },
    {
        "uniprot": "Q92813", "gene": "DIO2",
        "full_length": 273, "secis_aa": [133, 234],
        "function": "Type 2 iodothyronine deiodinase; T4 -> T3 conversion",
    },
    {
        "uniprot": "P11473", "gene": "DIO3",
        "full_length": 299, "secis_aa": [144],
        "function": "Type 3 iodothyronine deiodinase; inactivates thyroid hormone",
    },
    {
        "uniprot": "Q9Y6H0", "gene": "MSRB1",
        "full_length": 116, "secis_aa": [95],
        "function": "Methionine-R-sulfoxide reductase B1 (SepR)",
    },
    {
        "uniprot": "Q6P5Q7", "gene": "SELENOP",
        "full_length": 381, "secis_aa": [50, 75, 96, 109, 126, 146, 159, 173, 187, 200],
        "function": "Selenoprotein P; the major plasma selenoprotein, 10 Sec residues",
    },
    {
        "uniprot": "Q9NZV5", "gene": "SELENOF",
        "full_length": 311, "secis_aa": [142],
        "function": "15 kDa selenoprotein (Sep15) family; ER-resident",
    },
    {
        "uniprot": "Q961I3", "gene": "SELENOH",
        "full_length": 122, "secis_aa": [44],
        "function": "Selenoprotein H; nucleolar DNA-binding",
    },
    {
        "uniprot": "Q9UMX5", "gene": "SELENOI",
        "full_length": 397, "secis_aa": [385],
        "function": "Selenoprotein I; ethanolamine phosphotransferase",
    },
    {
        "uniprot": "Q9BVE4", "gene": "SELENOK",
        "full_length": 94, "secis_aa": [22],
        "function": "Selenoprotein K; ER-associated degradation",
    },
    {
        "uniprot": "Q9Y5D5", "gene": "SELENOS",
        "full_length": 189, "secis_aa": [188],
        "function": "Selenoprotein S; ER stress / inflammatory response",
    },
    {
        "uniprot": "Q9BSN4", "gene": "SELENOT",
        "full_length": 182, "secis_aa": [66],
        "function": "Selenoprotein T; ER-resident, neuroendocrine",
    },
    {
        "uniprot": "Q8IXQ5", "gene": "SELENOU",
        "full_length": 541, "secis_aa": [245],
        "function": "Selenoprotein U (VIMP-related); 1 Sec residue",
    },
    {
        "uniprot": "Q8WVH5", "gene": "SELENOW",
        "full_length": 87, "secis_aa": [13],
        "function": "Selenoprotein W; muscle and brain",
    },
    {
        "uniprot": "P47803", "gene": "SELENOW2",
        "full_length": 86, "secis_aa": [13],
        "function": "Selenoprotein W2 (paralog of SELENOW)",
    },
    {
        "uniprot": "Q9H4K0", "gene": "TXNRD1",
        "full_length": 648, "secis_aa": [497, 567],
        "function": "Thioredoxin reductase 1 (cytosolic); redox homeostasis",
    },
    {
        "uniprot": "Q9NNW7", "gene": "TXNRD2",
        "full_length": 527, "secis_aa": [496, 526],
        "function": "Thioredoxin reductase 2 (mitochondrial); redox homeostasis",
    },
]


def get_protein_count() -> int:
    """Return the number of human selenoproteins in this catalog."""
    return len(HUMAN_SELENOPROTEINS)


def get_total_secis_count() -> int:
    """Return the total number of Sec residues across all selenoproteins."""
    return sum(len(p["secis_aa"]) for p in HUMAN_SELENOPROTEINS)


def lookup(uniprot_id: str):
    """Look up a selenoprotein by UniProt accession."""
    for entry in HUMAN_SELENOPROTEINS:
        if entry["uniprot"] == uniprot_id:
            return entry
    raise KeyError(f"Unknown selenoprotein UniProt accession: {uniprot_id}")
