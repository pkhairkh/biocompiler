"""Benchmark CAI scores - directly imports to avoid __init__.py issues."""
import sys
import math

# Direct imports bypassing broken __init__.py
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

# Test genes
HBB_DNA = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCTCACTGCAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGC"

INS_DNA = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"

EGFP_DNA = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"

genes = {"HBB": HBB_DNA, "INS": INS_DNA, "EGFP": EGFP_DNA}
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
        print(f"    Current opt:     CAI={opt_cai:.4f}, GTs={opt_gts}, RS={rs_count}")
        
        # Analyze GTs
        gt_positions = [i for i in range(len(optimized)-1) if optimized[i:i+2] == "GT"]
        if gt_positions:
            print(f"    GT count: {len(gt_positions)}")
            for pos in gt_positions[:5]:
                codon_start = (pos // 3) * 3
                codon = optimized[codon_start:codon_start+3]
                aa = CODON_TABLE.get(codon, "?")
                is_cross = (pos + 1) % 3 == 0
                if is_cross:
                    next_cs = codon_start + 3
                    next_codon = optimized[next_cs:next_cs+3] if next_cs + 3 <= len(optimized) else "?"
                    next_aa = CODON_TABLE.get(next_codon, "?")
                    print(f"      pos {pos}: cross-codon GT, {codon}({aa})-{next_codon}({next_aa})")
                else:
                    print(f"      pos {pos}: within-codon GT in {codon}({aa})")
    
    avg_cai = total_cai / count if count else 0
    print(f"\n  Average CAI: {avg_cai:.4f}")
