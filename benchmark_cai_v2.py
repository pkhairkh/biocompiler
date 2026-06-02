"""Benchmark CAI scores with corrected gene sequences."""
import sys
import math

sys.modules['biocompiler'] = type(sys)('biocompiler')
sys.modules['biocompiler'].__path__ = ['/home/z/my-project/biocompiler/src/biocompiler']

from biocompiler.type_system import CODON_TABLE, AA_TO_CODONS
from biocompiler.organisms import SPECIES
from biocompiler.optimization import BioOptimizer
from biocompiler.restriction_sites import get_recognition_site

def compute_cai(seq, species_cai):
    if not seq or len(seq) < 3:
        return 0.0
    log_sum = 0.0
    count = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        cai = species_cai.get(codon, 0.0)
        if cai <= 0:
            cai = 0.001
        log_sum += math.log(cai)
        count += 1
    if count == 0:
        return 0.0
    return math.exp(log_sum / count)

def translate_seq(seq):
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        protein.append(aa)
    return "".join(protein)

def count_gts(seq):
    return sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "GT")

def count_restriction_sites(seq, enzymes):
    rs_count = 0
    for enz in enzymes:
        site = get_recognition_site(enz)
        if site:
            p = seq.find(site)
            while p != -1:
                rs_count += 1
                p = seq.find(site, p + 1)
    return rs_count

# Build DNA from just the protein sequence (CDS only, with stop codon)
def protein_to_dna(protein):
    """Convert protein to DNA using best-CAI codons, for baseline."""
    species_cai = SPECIES["human"]
    seq = ""
    for aa in protein:
        codons = AA_TO_CODONS.get(aa, [])
        if codons:
            best = max(codons, key=lambda c: species_cai.get(c, 0.0))
            seq += best
        else:
            seq += "NNN"
    seq += "TAA"  # Stop codon
    return seq

# Test genes - CDS only (coding sequence)
HBB_protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
INS_protein = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
EGFP_protein = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

genes = {
    "HBB": protein_to_dna(HBB_protein),
    "INS": protein_to_dna(INS_protein), 
    "EGFP": protein_to_dna(EGFP_protein),
}
enzymes = ["EcoRI", "BamHI", "HindIII", "XhoI"]

for species in ["human", "ecoli"]:
    species_cai = SPECIES.get(species, SPECIES["ecoli"])
    print(f"\n{'='*80}")
    print(f"  Species: {species}")
    print(f"{'='*80}")
    
    total_cai = 0
    count = 0
    
    for gene_name, gene_seq in genes.items():
        protein = translate_seq(gene_seq)
        
        # SimpleCAI baseline
        simple_seq = ""
        for aa in protein:
            if aa == "*":
                simple_seq += "TAA"
                continue
            codons = AA_TO_CODONS.get(aa, [])
            if codons:
                best = max(codons, key=lambda c: species_cai.get(c, 0.0))
                simple_seq += best
        simple_cai = compute_cai(simple_seq, species_cai)
        simple_gts = count_gts(simple_seq)
        
        # Current optimizer
        opt = BioOptimizer(species=species, enzymes=enzymes, avoid_gt=True)
        optimized, pred_results, cert_text = opt.optimize(gene_seq)
        opt_cai = compute_cai(optimized, species_cai)
        opt_gts = count_gts(optimized)
        rs_count = count_restriction_sites(optimized, enzymes)
        
        total_cai += opt_cai
        count += 1
        
        print(f"\n  {gene_name} ({len(protein)} aa):")
        print(f"    SimpleCAI (max): CAI={simple_cai:.4f}, GTs={simple_gts}")
        print(f"    Optimizer:       CAI={opt_cai:.4f}, GTs={opt_gts}, RS={rs_count}")
    
    avg_cai = total_cai / count if count else 0
    print(f"\n  Average CAI: {avg_cai:.4f}")
