# BioCompiler Full-Spectrum E2E Validation Report

## Overview

This report validates ALL 28 BioCompiler predicates across 5 layers
as a computational substitute for wet-lab testing:

| Layer | Predicates | Count |
|-------|-----------|-------|
| DNA | Codon fidelity, CAI, GC%, restriction sites, splice sites, CpG, GT avoidance, coding validity | 8 |
| mRNA | Cryptic promoters, TM domains, secondary structure, co-translational folding, stability motifs, GC range | 6 |
| Protein Structure | ESMFold confidence, misfolding risk, fold topology, unexpected interactions | 4 |
| Protein Biophysics | Stability, mutations, disulfides, hydrophobic core, solubility, aggregation, charge, hydrophobic stretches | 8 |
| Immunogenicity | Overall immunogenicity, T-cell epitopes, B-cell epitopes, population coverage | 4 |

## Summary by Gene

| Gene | Organism | Length | CAI | GC% | Opt Time | Total Preds | Pass | Uncertain | Fail | Pass Rate |
|------|----------|--------|-----|-----|----------|-------------|------|-----------|------|-----------|
| GFP | E.coli | 238 | 0.9989 | 0.5420 | 4614ms | 30 | 20 | 4 | 2 | 76.9% |
| mCherry | E.coli | 225 | 0.9966 | 0.5452 | 2ms | 30 | 17 | 7 | 2 | 65.4% |
| T4_lysozyme | E.coli | 216 | 1.0000 | 0.5463 | 1ms | 30 | 19 | 7 | 0 | 73.1% |
| groEL | E.coli | 338 | 0.9993 | 0.5306 | 2ms | 30 | 21 | 4 | 1 | 80.8% |
| Insulin | Human | 110 | 1.0000 | 0.6879 | 8ms | 30 | 17 | 4 | 5 | 65.4% |
| HBB | Human | 147 | 1.0000 | 0.6712 | 3ms | 30 | 16 | 5 | 5 | 61.5% |
| EPO | Human | 193 | 1.0000 | 0.6684 | 4ms | 30 | 13 | 5 | 8 | 50.0% |
| GH1 | Human | 217 | 1.0000 | 0.6144 | 8ms | 30 | 15 | 5 | 6 | 57.7% |
| IFNA2 | Human | 188 | 0.9999 | 0.6099 | 5ms | 30 | 16 | 4 | 6 | 61.5% |
| Albumin | Human | 607 | 1.0000 | 0.6216 | 12ms | 30 | 16 | 4 | 6 | 61.5% |
| TDH3 | Yeast | 190 | 0.9807 | 0.2965 | 3ms | 30 | 20 | 3 | 3 | 76.9% |
| PGK1 | Yeast | 240 | 0.9949 | 0.3306 | 3ms | 30 | 18 | 5 | 3 | 69.2% |
| ADH1 | Yeast | 306 | 0.9921 | 0.3007 | 8ms | 30 | 13 | 8 | 5 | 50.0% |
| mIL2 | Mouse | 163 | 1.0000 | 0.5685 | 3ms | 30 | 16 | 6 | 4 | 61.5% |
| mIFNG | Mouse | 292 | 1.0000 | 0.6210 | 5ms | 30 | 18 | 3 | 5 | 69.2% |
| CHO_EPO | CHO | 193 | 0.9900 | 0.6995 | 5ms | 30 | 13 | 5 | 8 | 50.0% |
| CHO_GFP | CHO | 238 | 1.0000 | 0.6275 | 4ms | 30 | 17 | 4 | 5 | 65.4% |

## Aggregate Statistics

- **Total genes tested**: 17
- **Total predicates evaluated**: 510
- **Overall pass rate**: 55.9% (285/510)
- **Uncertain verdicts**: 83 (16.3%)
- **Failed verdicts**: 74 (14.5%)
- **Mean CAI**: 0.9972

## Per-Organism Breakdown

### E. coli
- Genes: 4
- Pass rate: 64.2% (77/120)
- Mean CAI: 0.9987

### Human
- Genes: 6
- Pass rate: 51.7% (93/180)
- Mean CAI: 1.0000

### Yeast
- Genes: 3
- Pass rate: 56.7% (51/90)
- Mean CAI: 0.9892

### Mouse
- Genes: 2
- Pass rate: 56.7% (34/60)
- Mean CAI: 1.0000

### CHO-K1
- Genes: 2
- Pass rate: 50.0% (30/60)
- Mean CAI: 0.9950

## Detailed Predicate Results

### NoStopCodons
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | No internal stop codons |
| mCherry | E.coli | PASS | No internal stop codons |
| T4_lysozyme | E.coli | PASS | No internal stop codons |
| groEL | E.coli | PASS | No internal stop codons |
| Insulin | Human | PASS | No internal stop codons |
| HBB | Human | PASS | No internal stop codons |
| EPO | Human | PASS | No internal stop codons |
| GH1 | Human | PASS | No internal stop codons |
| IFNA2 | Human | PASS | No internal stop codons |
| Albumin | Human | PASS | No internal stop codons |
| TDH3 | Yeast | PASS | No internal stop codons |
| PGK1 | Yeast | PASS | No internal stop codons |
| ADH1 | Yeast | PASS | No internal stop codons |
| mIL2 | Mouse | PASS | No internal stop codons |
| mIFNG | Mouse | PASS | No internal stop codons |
| CHO_EPO | CHO | PASS | No internal stop codons |
| CHO_GFP | CHO | PASS | No internal stop codons |

### NoCrypticSplice
Pass rate: 23.5% (4/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | Cryptic splice check skipped for prokaryotic organism 'Escherichia_coli' |
| mCherry | E.coli | PASS | Cryptic splice check skipped for prokaryotic organism 'Escherichia_coli' |
| T4_lysozyme | E.coli | PASS | Cryptic splice check skipped for prokaryotic organism 'Escherichia_coli' |
| groEL | E.coli | PASS | Cryptic splice check skipped for prokaryotic organism 'Escherichia_coli' |
| Insulin | Human | FAIL | Worst splice score 7.11 at pos 47 |
| HBB | Human | FAIL | Worst splice score 7.11 at pos 44 |
| EPO | Human | FAIL | Worst splice score 6.02 at pos 72 |
| GH1 | Human | UNCERTAIN | Worst splice score 5.92 at pos 233 |
| IFNA2 | Human | FAIL | Worst splice score 6.02 at pos 231 |
| Albumin | Human | FAIL | Worst splice score 6.31 at pos 1448 |
| TDH3 | Yeast | UNCERTAIN | Worst splice score 5.69 at pos 102 |
| PGK1 | Yeast | UNCERTAIN | Worst splice score 5.72 at pos 614 |
| ADH1 | Yeast | UNCERTAIN | Worst splice score 5.50 at pos 285 |
| mIL2 | Mouse | FAIL | Worst splice score 6.15 at pos 467 |
| mIFNG | Mouse | UNCERTAIN | Worst splice score 5.22 at pos 156 |
| CHO_EPO | CHO | FAIL | Worst splice score 6.02 at pos 72 |
| CHO_GFP | CHO | FAIL | Worst splice score 6.02 at pos 576 |

### NoCpGIsland
Pass rate: 58.8% (10/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | CpG island check skipped for prokaryotic organism 'Escherichia_coli' |
| mCherry | E.coli | PASS | CpG island check skipped for prokaryotic organism 'Escherichia_coli' |
| T4_lysozyme | E.coli | PASS | CpG island check skipped for prokaryotic organism 'Escherichia_coli' |
| groEL | E.coli | PASS | CpG island check skipped for prokaryotic organism 'Escherichia_coli' |
| Insulin | Human | PASS | Worst CpG Obs/Exp ratio 0.533 <= 0.6 |
| HBB | Human | FAIL | CpG island at pos 23, Obs/Exp=0.758 > 0.6 |
| EPO | Human | FAIL | CpG island at pos 374, Obs/Exp=0.652 > 0.6 |
| GH1 | Human | FAIL | CpG island at pos 359, Obs/Exp=0.766 > 0.6 |
| IFNA2 | Human | FAIL | CpG island at pos 193, Obs/Exp=0.689 > 0.6 |
| Albumin | Human | FAIL | CpG island at pos 177, Obs/Exp=1.012 > 0.6 |
| TDH3 | Yeast | PASS | Worst CpG Obs/Exp ratio 0.000 <= 0.6 |
| PGK1 | Yeast | PASS | Worst CpG Obs/Exp ratio 0.000 <= 0.6 |
| ADH1 | Yeast | PASS | Worst CpG Obs/Exp ratio 0.000 <= 0.6 |
| mIL2 | Mouse | PASS | Worst CpG Obs/Exp ratio 0.417 <= 0.6 |
| mIFNG | Mouse | PASS | Worst CpG Obs/Exp ratio 0.593 <= 0.6 |
| CHO_EPO | CHO | FAIL | CpG island at pos 378, Obs/Exp=0.703 > 0.6 |
| CHO_GFP | CHO | FAIL | CpG island at pos 195, Obs/Exp=0.984 > 0.6 |

### NoRestrictionSite
Pass rate: 88.2% (15/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | No restriction sites found |
| mCherry | E.coli | FAIL | Restriction sites found at [418] |
| T4_lysozyme | E.coli | PASS | No restriction sites found |
| groEL | E.coli | PASS | No restriction sites found |
| Insulin | Human | PASS | No restriction sites found |
| HBB | Human | PASS | No restriction sites found |
| EPO | Human | PASS | No restriction sites found |
| GH1 | Human | PASS | No restriction sites found |
| IFNA2 | Human | PASS | No restriction sites found |
| Albumin | Human | PASS | No restriction sites found |
| TDH3 | Yeast | PASS | No restriction sites found |
| PGK1 | Yeast | PASS | No restriction sites found |
| ADH1 | Yeast | FAIL | Restriction sites found at [625] |
| mIL2 | Mouse | PASS | No restriction sites found |
| mIFNG | Mouse | PASS | No restriction sites found |
| CHO_EPO | CHO | PASS | No restriction sites found |
| CHO_GFP | CHO | PASS | No restriction sites found |

### NoAvoidableGT
Pass rate: 23.5% (4/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | GT dinucleotide check skipped for prokaryotic organism 'Escherichia_coli' |
| mCherry | E.coli | PASS | GT dinucleotide check skipped for prokaryotic organism 'Escherichia_coli' |
| T4_lysozyme | E.coli | PASS | GT dinucleotide check skipped for prokaryotic organism 'Escherichia_coli' |
| groEL | E.coli | PASS | GT dinucleotide check skipped for prokaryotic organism 'Escherichia_coli' |
| Insulin | Human | FAIL | Avoidable GT dinucleotides at [8, 47, 89, 116, 281, 305]; unavoidable at [75, 10 |
| HBB | Human | FAIL | Avoidable GT dinucleotides at [44, 365, 434]; unavoidable at [3, 33, 54, 60, 69, |
| EPO | Human | FAIL | Avoidable GT dinucleotides at [14, 29, 503, 545]; unavoidable at [6, 72, 111, 20 |
| GH1 | Human | FAIL | Avoidable GT dinucleotides at [47, 104, 167, 200, 233, 350, 512, 566, 620]; unav |
| IFNA2 | Human | FAIL | Avoidable GT dinucleotides at [146, 194, 206, 266, 317, 332, 431, 452, 470]; una |
| Albumin | Human | FAIL | Avoidable GT dinucleotides at [5, 29, 35, 158, 170, 215, 278, 293, 371, 482, 488 |
| TDH3 | Yeast | FAIL | Avoidable GT dinucleotides at [22, 82, 109, 151, 181, 229, 241, 331, 334, 376, 4 |
| PGK1 | Yeast | FAIL | Avoidable GT dinucleotides at [2, 58, 85, 130, 139, 146, 157, 181, 199, 220, 223 |
| ADH1 | Yeast | FAIL | Avoidable GT dinucleotides at [2, 23, 64, 97, 103, 112, 133, 181, 184, 187, 196, |
| mIL2 | Mouse | FAIL | Avoidable GT dinucleotides at [143, 296, 404, 461, 467]; unavoidable at [2, 27,  |
| mIFNG | Mouse | FAIL | Avoidable GT dinucleotides at [632, 692, 824, 872]; unavoidable at [2, 9, 15, 21 |
| CHO_EPO | CHO | FAIL | Avoidable GT dinucleotides at [14, 29, 122, 166, 503, 545]; unavoidable at [6, 7 |
| CHO_GFP | CHO | FAIL | Avoidable GT dinucleotides at [20, 77, 134, 206, 218, 338, 425, 665, 707]; unavo |

### ValidCodingSeq
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | All codons valid |
| mCherry | E.coli | PASS | All codons valid |
| T4_lysozyme | E.coli | PASS | All codons valid |
| groEL | E.coli | PASS | All codons valid |
| Insulin | Human | PASS | All codons valid |
| HBB | Human | PASS | All codons valid |
| EPO | Human | PASS | All codons valid |
| GH1 | Human | PASS | All codons valid |
| IFNA2 | Human | PASS | All codons valid |
| Albumin | Human | PASS | All codons valid |
| TDH3 | Yeast | PASS | All codons valid |
| PGK1 | Yeast | PASS | All codons valid |
| ADH1 | Yeast | PASS | All codons valid |
| mIL2 | Mouse | PASS | All codons valid |
| mIFNG | Mouse | PASS | All codons valid |
| CHO_EPO | CHO | PASS | All codons valid |
| CHO_GFP | CHO | PASS | All codons valid |

### ConservationScore
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | Translation matches original protein |
| mCherry | E.coli | PASS | Translation matches original protein |
| T4_lysozyme | E.coli | PASS | Translation matches original protein |
| groEL | E.coli | PASS | Translation matches original protein |
| Insulin | Human | PASS | Translation matches original protein |
| HBB | Human | PASS | Translation matches original protein |
| EPO | Human | PASS | Translation matches original protein |
| GH1 | Human | PASS | Translation matches original protein |
| IFNA2 | Human | PASS | Translation matches original protein |
| Albumin | Human | PASS | Translation matches original protein |
| TDH3 | Yeast | PASS | Translation matches original protein |
| PGK1 | Yeast | PASS | Translation matches original protein |
| ADH1 | Yeast | PASS | Translation matches original protein |
| mIL2 | Mouse | PASS | Translation matches original protein |
| mIFNG | Mouse | PASS | Translation matches original protein |
| CHO_EPO | CHO | PASS | Translation matches original protein |
| CHO_GFP | CHO | PASS | Translation matches original protein |

### CodonOptimality
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | CAI=0.9989 |
| mCherry | E.coli | PASS | CAI=0.9966 |
| T4_lysozyme | E.coli | PASS | CAI=1.0000 |
| groEL | E.coli | PASS | CAI=0.9993 |
| Insulin | Human | PASS | CAI=1.0000 |
| HBB | Human | PASS | CAI=1.0000 |
| EPO | Human | PASS | CAI=1.0000 |
| GH1 | Human | PASS | CAI=1.0000 |
| IFNA2 | Human | PASS | CAI=0.9999 |
| Albumin | Human | PASS | CAI=1.0000 |
| TDH3 | Yeast | PASS | CAI=0.9807 |
| PGK1 | Yeast | PASS | CAI=0.9949 |
| ADH1 | Yeast | PASS | CAI=0.9921 |
| mIL2 | Mouse | PASS | CAI=1.0000 |
| mIFNG | Mouse | PASS | CAI=1.0000 |
| CHO_EPO | CHO | PASS | CAI=0.9900 |
| CHO_GFP | CHO | PASS | CAI=1.0000 |

### NoCrypticPromoter
Pass rate: 0.0% (0/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | UNCERTAIN | Worst promoter score 0.583 at pos 129 |
| mCherry | E.coli | UNCERTAIN | Worst promoter score 0.583 at pos 42 |
| T4_lysozyme | E.coli | UNCERTAIN | Worst promoter score 0.583 at pos 28 |
| groEL | E.coli | UNCERTAIN | Worst promoter score 0.667 at pos 580 |
| Insulin | Human | FAIL | Worst promoter score 0.750 at pos 115 |
| HBB | Human | FAIL | Worst promoter score 0.750 at pos 282 |
| EPO | Human | FAIL | Worst promoter score 0.762 at pos 494 |
| GH1 | Human | FAIL | Worst promoter score 0.750 at pos 151 |
| IFNA2 | Human | FAIL | Worst promoter score 0.762 at pos 398 |
| Albumin | Human | FAIL | Worst promoter score 0.762 at pos 362 |
| TDH3 | Yeast | FAIL | Worst promoter score 0.917 at pos 93 |
| PGK1 | Yeast | FAIL | Worst promoter score 0.917 at pos 75 |
| ADH1 | Yeast | FAIL | Worst promoter score 0.917 at pos 413 |
| mIL2 | Mouse | FAIL | Worst promoter score 0.762 at pos 1 |
| mIFNG | Mouse | FAIL | Worst promoter score 0.762 at pos 143 |
| CHO_EPO | CHO | FAIL | Worst promoter score 0.762 at pos 511 |
| CHO_GFP | CHO | FAIL | Worst promoter score 0.762 at pos 230 |

### NoUnexpectedTMDomain
Pass rate: 5.9% (1/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.579 at AA pos 217 exceeds 0.5 |
| mCherry | E.coli | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.579 at AA pos 176 exceeds 0.5 |
| T4_lysozyme | E.coli | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.632 at AA pos 86 exceeds 0.57 |
| groEL | E.coli | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.579 at AA pos 4 exceeds 0.578 |
| Insulin | Human | FAIL | TM domain detected: worst hydrophobic fraction 0.789 at AA pos 0 exceeds thresho |
| HBB | Human | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.632 at AA pos 124 exceeds 0.5 |
| EPO | Human | FAIL | TM domain detected: worst hydrophobic fraction 0.737 at AA pos 7 exceeds thresho |
| GH1 | Human | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.632 at AA pos 8 exceeds 0.578 |
| IFNA2 | Human | FAIL | TM domain detected: worst hydrophobic fraction 0.684 at AA pos 0 exceeds thresho |
| Albumin | Human | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.632 at AA pos 0 exceeds 0.578 |
| TDH3 | Yeast | PASS | No TM domain detected (worst hydrophobic fraction 0.474) |
| PGK1 | Yeast | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.579 at AA pos 200 exceeds 0.5 |
| ADH1 | Yeast | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.632 at AA pos 255 exceeds 0.5 |
| mIL2 | Mouse | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.632 at AA pos 0 exceeds 0.578 |
| mIFNG | Mouse | FAIL | TM domain detected: worst hydrophobic fraction 0.737 at AA pos 3 exceeds thresho |
| CHO_EPO | CHO | FAIL | TM domain detected: worst hydrophobic fraction 0.737 at AA pos 7 exceeds thresho |
| CHO_GFP | CHO | UNCERTAIN | Borderline TM domain: worst hydrophobic fraction 0.579 at AA pos 217 exceeds 0.5 |

### mRNASecondaryStructure
Pass rate: 0.0% (0/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | FAIL | Strong mRNA secondary structure (Nussinov fallback): ΔG≈-16.1 kcal/mol <= -15.0 |
| mCherry | E.coli | UNCERTAIN | Moderate mRNA secondary structure (Nussinov fallback): ΔG≈-10.5 kcal/mol <= -10. |
| T4_lysozyme | E.coli | UNCERTAIN | Moderate mRNA secondary structure (Nussinov fallback): ΔG≈-11.9 kcal/mol <= -10. |
| groEL | E.coli | FAIL | Strong mRNA secondary structure (Nussinov fallback): ΔG≈-17.4 kcal/mol <= -15.0 |
| Insulin | Human | UNCERTAIN | Moderate mRNA secondary structure (Nussinov fallback): ΔG≈-14.0 kcal/mol <= -10. |
| HBB | Human | UNCERTAIN | Moderate mRNA secondary structure (Nussinov fallback): ΔG≈-13.9 kcal/mol <= -10. |
| EPO | Human | FAIL | Strong mRNA secondary structure (Nussinov fallback): ΔG≈-16.1 kcal/mol <= -15.0 |
| GH1 | Human | FAIL | Strong mRNA secondary structure (Nussinov fallback): ΔG≈-19.6 kcal/mol <= -15.0 |
| IFNA2 | Human | UNCERTAIN | Moderate mRNA secondary structure (Nussinov fallback): ΔG≈-12.3 kcal/mol <= -10. |
| Albumin | Human | UNCERTAIN | Moderate mRNA secondary structure (Nussinov fallback): ΔG≈-13.4 kcal/mol <= -10. |
| TDH3 | Yeast | FAIL | Strong mRNA secondary structure (Nussinov fallback): ΔG≈-16.9 kcal/mol <= -15.0 |
| PGK1 | Yeast | UNCERTAIN | Moderate mRNA secondary structure (Nussinov fallback): ΔG≈-11.7 kcal/mol <= -10. |
| ADH1 | Yeast | FAIL | Strong mRNA secondary structure (Nussinov fallback): ΔG≈-17.8 kcal/mol <= -15.0 |
| mIL2 | Mouse | FAIL | Strong mRNA secondary structure (Nussinov fallback): ΔG≈-15.9 kcal/mol <= -15.0 |
| mIFNG | Mouse | FAIL | Strong mRNA secondary structure (Nussinov fallback): ΔG≈-24.9 kcal/mol <= -15.0 |
| CHO_EPO | CHO | FAIL | Strong mRNA secondary structure (Nussinov fallback): ΔG≈-41.6 kcal/mol <= -15.0 |
| CHO_GFP | CHO | UNCERTAIN | Moderate mRNA secondary structure (Nussinov fallback): ΔG≈-11.5 kcal/mol <= -10. |

### CoTranslationalFolding
Pass rate: 0.0% (0/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| mCherry | E.coli | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| T4_lysozyme | E.coli | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| groEL | E.coli | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| Insulin | Human | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| HBB | Human | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| EPO | Human | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| GH1 | Human | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| IFNA2 | Human | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| Albumin | Human | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| TDH3 | Yeast | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| PGK1 | Yeast | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| ADH1 | Yeast | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| mIL2 | Mouse | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| mIFNG | Mouse | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| CHO_EPO | CHO | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |
| CHO_GFP | CHO | UNCERTAIN | Error: '<' not supported between instances of 'str' and 'int' |

### MRNAStability
Pass rate: 0.0% (0/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| mCherry | E.coli | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| T4_lysozyme | E.coli | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| groEL | E.coli | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| Insulin | Human | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| HBB | Human | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| EPO | Human | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| GH1 | Human | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| IFNA2 | Human | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| Albumin | Human | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| TDH3 | Yeast | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| PGK1 | Yeast | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| ADH1 | Yeast | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| mIL2 | Mouse | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| mIFNG | Mouse | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| CHO_EPO | CHO | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |
| CHO_GFP | CHO | UNCERTAIN | Error: cannot import name 'scan_instability_motifs' from 'biocompiler.mrna_stabi |

### GCContent
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | GC=0.5420 |
| mCherry | E.coli | PASS | GC=0.5452 |
| T4_lysozyme | E.coli | PASS | GC=0.5463 |
| groEL | E.coli | PASS | GC=0.5306 |
| Insulin | Human | PASS | GC=0.6879 |
| HBB | Human | PASS | GC=0.6712 |
| EPO | Human | PASS | GC=0.6684 |
| GH1 | Human | PASS | GC=0.6144 |
| IFNA2 | Human | PASS | GC=0.6099 |
| Albumin | Human | PASS | GC=0.6216 |
| TDH3 | Yeast | LIKELY_PASS | GC=0.2965 |
| PGK1 | Yeast | PASS | GC=0.3306 |
| ADH1 | Yeast | PASS | GC=0.3007 |
| mIL2 | Mouse | PASS | GC=0.5685 |
| mIFNG | Mouse | PASS | GC=0.6210 |
| CHO_EPO | CHO | PASS | GC=0.6995 |
| CHO_GFP | CHO | PASS | GC=0.6275 |

### StructureConfidence
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | LIKELY_PASS | pLDDT=53.2 (method=heuristic_fallback) |
| mCherry | E.coli | LIKELY_PASS | pLDDT=52.4 (method=heuristic_fallback) |
| T4_lysozyme | E.coli | LIKELY_PASS | pLDDT=54.0 (method=heuristic_fallback) |
| groEL | E.coli | LIKELY_PASS | pLDDT=53.8 (method=heuristic_fallback) |
| Insulin | Human | LIKELY_PASS | pLDDT=53.4 (method=heuristic_fallback) |
| HBB | Human | LIKELY_PASS | pLDDT=53.5 (method=heuristic_fallback) |
| EPO | Human | LIKELY_PASS | pLDDT=53.2 (method=heuristic_fallback) |
| GH1 | Human | LIKELY_PASS | pLDDT=52.8 (method=heuristic_fallback) |
| IFNA2 | Human | LIKELY_PASS | pLDDT=53.7 (method=heuristic_fallback) |
| Albumin | Human | LIKELY_PASS | pLDDT=52.5 (method=heuristic_fallback) |
| TDH3 | Yeast | LIKELY_PASS | pLDDT=53.7 (method=heuristic_fallback) |
| PGK1 | Yeast | LIKELY_PASS | pLDDT=52.8 (method=heuristic_fallback) |
| ADH1 | Yeast | LIKELY_PASS | pLDDT=52.5 (method=heuristic_fallback) |
| mIL2 | Mouse | LIKELY_PASS | pLDDT=54.0 (method=heuristic_fallback) |
| mIFNG | Mouse | LIKELY_PASS | pLDDT=53.6 (method=heuristic_fallback) |
| CHO_EPO | CHO | LIKELY_PASS | pLDDT=53.2 (method=heuristic_fallback) |
| CHO_GFP | CHO | LIKELY_PASS | pLDDT=53.2 (method=heuristic_fallback) |

### NoMisfoldingRisk
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | 0 long hydrophobic stretches (>15aa) |
| mCherry | E.coli | PASS | 0 long hydrophobic stretches (>15aa) |
| T4_lysozyme | E.coli | PASS | 0 long hydrophobic stretches (>15aa) |
| groEL | E.coli | PASS | 0 long hydrophobic stretches (>15aa) |
| Insulin | Human | PASS | 0 long hydrophobic stretches (>15aa) |
| HBB | Human | PASS | 0 long hydrophobic stretches (>15aa) |
| EPO | Human | PASS | 0 long hydrophobic stretches (>15aa) |
| GH1 | Human | PASS | 0 long hydrophobic stretches (>15aa) |
| IFNA2 | Human | PASS | 0 long hydrophobic stretches (>15aa) |
| Albumin | Human | PASS | 0 long hydrophobic stretches (>15aa) |
| TDH3 | Yeast | PASS | 0 long hydrophobic stretches (>15aa) |
| PGK1 | Yeast | PASS | 0 long hydrophobic stretches (>15aa) |
| ADH1 | Yeast | PASS | 0 long hydrophobic stretches (>15aa) |
| mIL2 | Mouse | PASS | 0 long hydrophobic stretches (>15aa) |
| mIFNG | Mouse | PASS | 0 long hydrophobic stretches (>15aa) |
| CHO_EPO | CHO | PASS | 0 long hydrophobic stretches (>15aa) |
| CHO_GFP | CHO | PASS | 0 long hydrophobic stretches (>15aa) |

### CorrectFoldTopology
Pass rate: 94.1% (16/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS |  |
| mCherry | E.coli | UNCERTAIN | Hydrophobic fraction 0.293 is below normal range [0.3, 0.45] -- insufficient hyd |
| T4_lysozyme | E.coli | PASS |  |
| groEL | E.coli | PASS |  |
| Insulin | Human | PASS |  |
| HBB | Human | PASS |  |
| EPO | Human | PASS |  |
| GH1 | Human | PASS |  |
| IFNA2 | Human | PASS |  |
| Albumin | Human | PASS |  |
| TDH3 | Yeast | PASS |  |
| PGK1 | Yeast | PASS |  |
| ADH1 | Yeast | PASS |  |
| mIL2 | Mouse | PASS |  |
| mIFNG | Mouse | PASS |  |
| CHO_EPO | CHO | PASS |  |
| CHO_GFP | CHO | PASS |  |

### NoUnexpectedInteraction
Pass rate: 70.6% (12/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS |  |
| mCherry | E.coli | PASS |  |
| T4_lysozyme | E.coli | UNCERTAIN | High isoelectric point (pI=10.25 > 9.0): protein may precipitate near its pI in  |
| groEL | E.coli | PASS |  |
| Insulin | Human | PASS |  |
| HBB | Human | PASS |  |
| EPO | Human | UNCERTAIN | High isoelectric point (pI=9.86 > 9.0): protein may precipitate near its pI in t |
| GH1 | Human | PASS |  |
| IFNA2 | Human | PASS |  |
| Albumin | Human | PASS |  |
| TDH3 | Yeast | PASS |  |
| PGK1 | Yeast | PASS |  |
| ADH1 | Yeast | UNCERTAIN | High isoelectric point (pI=10.46 > 9.0): protein may precipitate near its pI in  |
| mIL2 | Mouse | UNCERTAIN | High isoelectric point (pI=10.32 > 9.0): protein may precipitate near its pI in  |
| mIFNG | Mouse | PASS |  |
| CHO_EPO | CHO | UNCERTAIN | High isoelectric point (pI=9.86 > 9.0): protein may precipitate near its pI in t |
| CHO_GFP | CHO | PASS |  |

### StableFolding
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS |  |
| mCherry | E.coli | PASS |  |
| T4_lysozyme | E.coli | PASS |  |
| groEL | E.coli | PASS |  |
| Insulin | Human | PASS |  |
| HBB | Human | PASS |  |
| EPO | Human | PASS |  |
| GH1 | Human | PASS |  |
| IFNA2 | Human | PASS |  |
| Albumin | Human | PASS |  |
| TDH3 | Yeast | PASS |  |
| PGK1 | Yeast | PASS |  |
| ADH1 | Yeast | PASS |  |
| mIL2 | Mouse | PASS |  |
| mIFNG | Mouse | PASS |  |
| CHO_EPO | CHO | PASS |  |
| CHO_GFP | CHO | PASS |  |

### NoDestabilizingMutation
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS | No mutations (synonymous codon optimization preserves protein) |
| mCherry | E.coli | PASS | No mutations (synonymous codon optimization preserves protein) |
| T4_lysozyme | E.coli | PASS | No mutations (synonymous codon optimization preserves protein) |
| groEL | E.coli | PASS | No mutations (synonymous codon optimization preserves protein) |
| Insulin | Human | PASS | No mutations (synonymous codon optimization preserves protein) |
| HBB | Human | PASS | No mutations (synonymous codon optimization preserves protein) |
| EPO | Human | PASS | No mutations (synonymous codon optimization preserves protein) |
| GH1 | Human | PASS | No mutations (synonymous codon optimization preserves protein) |
| IFNA2 | Human | PASS | No mutations (synonymous codon optimization preserves protein) |
| Albumin | Human | PASS | No mutations (synonymous codon optimization preserves protein) |
| TDH3 | Yeast | PASS | No mutations (synonymous codon optimization preserves protein) |
| PGK1 | Yeast | PASS | No mutations (synonymous codon optimization preserves protein) |
| ADH1 | Yeast | PASS | No mutations (synonymous codon optimization preserves protein) |
| mIL2 | Mouse | PASS | No mutations (synonymous codon optimization preserves protein) |
| mIFNG | Mouse | PASS | No mutations (synonymous codon optimization preserves protein) |
| CHO_EPO | CHO | PASS | No mutations (synonymous codon optimization preserves protein) |
| CHO_GFP | CHO | PASS | No mutations (synonymous codon optimization preserves protein) |

### DisulfideBondIntegrity
Pass rate: 58.8% (10/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS |  |
| mCherry | E.coli | LIKELY_FAIL | Odd number of cysteines (1): at least one unpaired Cys that may form incorrect d |
| T4_lysozyme | E.coli | PASS |  |
| groEL | E.coli | PASS |  |
| Insulin | Human | PASS |  |
| HBB | Human | PASS |  |
| EPO | Human | LIKELY_FAIL | Odd number of cysteines (5): at least one unpaired Cys that may form incorrect d |
| GH1 | Human | LIKELY_FAIL | Odd number of cysteines (5): at least one unpaired Cys that may form incorrect d |
| IFNA2 | Human | PASS |  |
| Albumin | Human | LIKELY_FAIL | Odd number of cysteines (35): at least one unpaired Cys that may form incorrect  |
| TDH3 | Yeast | PASS |  |
| PGK1 | Yeast | LIKELY_FAIL | Odd number of cysteines (1): at least one unpaired Cys that may form incorrect d |
| ADH1 | Yeast | LIKELY_FAIL | Odd number of cysteines (1): at least one unpaired Cys that may form incorrect d |
| mIL2 | Mouse | PASS |  |
| mIFNG | Mouse | PASS |  |
| CHO_EPO | CHO | LIKELY_FAIL | Odd number of cysteines (5): at least one unpaired Cys that may form incorrect d |
| CHO_GFP | CHO | PASS |  |

### HydrophobicCoreQuality
Pass rate: 94.1% (16/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS |  |
| mCherry | E.coli | UNCERTAIN | Hydrophobic fraction 0.293 is below normal range [0.3, 0.45] -- insufficient hyd |
| T4_lysozyme | E.coli | PASS |  |
| groEL | E.coli | PASS |  |
| Insulin | Human | PASS |  |
| HBB | Human | PASS |  |
| EPO | Human | PASS |  |
| GH1 | Human | PASS |  |
| IFNA2 | Human | PASS |  |
| Albumin | Human | PASS |  |
| TDH3 | Yeast | PASS |  |
| PGK1 | Yeast | PASS |  |
| ADH1 | Yeast | PASS |  |
| mIL2 | Mouse | PASS |  |
| mIFNG | Mouse | PASS |  |
| CHO_EPO | CHO | PASS |  |
| CHO_GFP | CHO | PASS |  |

### SolubleExpression
Pass rate: 58.8% (10/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | LIKELY_PASS |  |
| mCherry | E.coli | LIKELY_PASS |  |
| T4_lysozyme | E.coli | LIKELY_PASS |  |
| groEL | E.coli | LIKELY_PASS |  |
| Insulin | Human | UNCERTAIN | Marginal solubility: CamSol score -0.512 is in the uncertain range [-1.0, 0.0) |
| HBB | Human | UNCERTAIN | Marginal solubility: CamSol score -0.253 is in the uncertain range [-1.0, 0.0) |
| EPO | Human | UNCERTAIN | Marginal solubility: CamSol score -0.233 is in the uncertain range [-1.0, 0.0) |
| GH1 | Human | UNCERTAIN | Marginal solubility: CamSol score -0.140 is in the uncertain range [-1.0, 0.0) |
| IFNA2 | Human | UNCERTAIN | Marginal solubility: CamSol score -0.199 is in the uncertain range [-1.0, 0.0) |
| Albumin | Human | LIKELY_PASS |  |
| TDH3 | Yeast | LIKELY_PASS |  |
| PGK1 | Yeast | LIKELY_PASS |  |
| ADH1 | Yeast | UNCERTAIN | Marginal solubility: CamSol score -0.065 is in the uncertain range [-1.0, 0.0) |
| mIL2 | Mouse | LIKELY_PASS |  |
| mIFNG | Mouse | LIKELY_PASS |  |
| CHO_EPO | CHO | UNCERTAIN | Marginal solubility: CamSol score -0.233 is in the uncertain range [-1.0, 0.0) |
| CHO_GFP | CHO | LIKELY_PASS |  |

### NoAggregationProneRegion
Pass rate: 29.4% (5/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | LIKELY_FAIL | Long aggregation-prone region of 14 residues (max allowed: 5) |
| mCherry | E.coli | LIKELY_PASS | Borderline aggregation-prone region of 7 residues (max allowed: 5) |
| T4_lysozyme | E.coli | LIKELY_PASS | Borderline aggregation-prone region of 6 residues (max allowed: 5) |
| groEL | E.coli | LIKELY_PASS | Borderline aggregation-prone region of 7 residues (max allowed: 5) |
| Insulin | Human | FAIL | Very long aggregation-prone region of 17 residues (max allowed: 5) |
| HBB | Human | LIKELY_FAIL | Long aggregation-prone region of 11 residues (max allowed: 5) |
| EPO | Human | FAIL | Very long aggregation-prone region of 19 residues (max allowed: 5) |
| GH1 | Human | LIKELY_FAIL | Long aggregation-prone region of 11 residues (max allowed: 5) |
| IFNA2 | Human | LIKELY_FAIL | Long aggregation-prone region of 14 residues (max allowed: 5) |
| Albumin | Human | LIKELY_FAIL | Long aggregation-prone region of 14 residues (max allowed: 5) |
| TDH3 | Yeast | PASS |  |
| PGK1 | Yeast | PASS |  |
| ADH1 | Yeast | UNCERTAIN | Aggregation-prone region of 10 residues detected (max allowed: 5) |
| mIL2 | Mouse | UNCERTAIN | Aggregation-prone region of 9 residues detected (max allowed: 5) |
| mIFNG | Mouse | FAIL | Very long aggregation-prone region of 17 residues (max allowed: 5) |
| CHO_EPO | CHO | FAIL | Very long aggregation-prone region of 19 residues (max allowed: 5) |
| CHO_GFP | CHO | LIKELY_FAIL | Long aggregation-prone region of 14 residues (max allowed: 5) |

### ChargeComposition
Pass rate: 70.6% (12/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS |  |
| mCherry | E.coli | PASS |  |
| T4_lysozyme | E.coli | UNCERTAIN | High isoelectric point (pI=10.25 > 9.0): protein may precipitate near its pI in  |
| groEL | E.coli | PASS |  |
| Insulin | Human | PASS |  |
| HBB | Human | PASS |  |
| EPO | Human | UNCERTAIN | High isoelectric point (pI=9.86 > 9.0): protein may precipitate near its pI in t |
| GH1 | Human | PASS |  |
| IFNA2 | Human | PASS |  |
| Albumin | Human | PASS |  |
| TDH3 | Yeast | PASS |  |
| PGK1 | Yeast | PASS |  |
| ADH1 | Yeast | UNCERTAIN | High isoelectric point (pI=10.46 > 9.0): protein may precipitate near its pI in  |
| mIL2 | Mouse | UNCERTAIN | High isoelectric point (pI=10.32 > 9.0): protein may precipitate near its pI in  |
| mIFNG | Mouse | PASS |  |
| CHO_EPO | CHO | UNCERTAIN | High isoelectric point (pI=9.86 > 9.0): protein may precipitate near its pI in t |
| CHO_GFP | CHO | PASS |  |

### NoLongHydrophobicStretch
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS |  |
| mCherry | E.coli | PASS |  |
| T4_lysozyme | E.coli | PASS |  |
| groEL | E.coli | PASS |  |
| Insulin | Human | LIKELY_PASS | Hydrophobic stretch of 8 residues slightly exceeds limit of 7 (borderline) |
| HBB | Human | PASS |  |
| EPO | Human | PASS |  |
| GH1 | Human | PASS |  |
| IFNA2 | Human | LIKELY_PASS | Hydrophobic stretch of 10 residues slightly exceeds limit of 7 (borderline) |
| Albumin | Human | LIKELY_PASS | Hydrophobic stretch of 8 residues slightly exceeds limit of 7 (borderline) |
| TDH3 | Yeast | PASS |  |
| PGK1 | Yeast | PASS |  |
| ADH1 | Yeast | PASS |  |
| mIL2 | Mouse | PASS |  |
| mIFNG | Mouse | PASS |  |
| CHO_EPO | CHO | PASS |  |
| CHO_GFP | CHO | PASS |  |

### LowImmunogenicity
Pass rate: 0.0% (0/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | FAIL | Immunogenicity score 0.872 is high (>0.7) |
| mCherry | E.coli | FAIL | Immunogenicity score 0.900 is high (>0.7) |
| T4_lysozyme | E.coli | FAIL | Immunogenicity score 0.915 is high (>0.7) |
| groEL | E.coli | FAIL | Immunogenicity score 0.905 is high (>0.7) |
| Insulin | Human | FAIL | Immunogenicity score 0.975 is high (>0.7) |
| HBB | Human | FAIL | Immunogenicity score 0.951 is high (>0.7) |
| EPO | Human | FAIL | Immunogenicity score 0.944 is high (>0.7) |
| GH1 | Human | FAIL | Immunogenicity score 0.906 is high (>0.7) |
| IFNA2 | Human | FAIL | Immunogenicity score 0.968 is high (>0.7) |
| Albumin | Human | FAIL | Immunogenicity score 0.955 is high (>0.7) |
| TDH3 | Yeast | FAIL | Immunogenicity score 0.810 is high (>0.7) |
| PGK1 | Yeast | FAIL | Immunogenicity score 0.920 is high (>0.7) |
| ADH1 | Yeast | FAIL | Immunogenicity score 0.881 is high (>0.7) |
| mIL2 | Mouse | FAIL | Immunogenicity score 0.961 is high (>0.7) |
| mIFNG | Mouse | FAIL | Immunogenicity score 0.948 is high (>0.7) |
| CHO_EPO | CHO | FAIL | Immunogenicity score 0.944 is high (>0.7) |
| CHO_GFP | CHO | FAIL | Immunogenicity score 0.872 is high (>0.7) |

### NoStrongTCellEpitope
Pass rate: 0.0% (0/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | FAIL | 356 strong T-cell epitopes found (high immunogenicity risk) |
| mCherry | E.coli | FAIL | 331 strong T-cell epitopes found (high immunogenicity risk) |
| T4_lysozyme | E.coli | FAIL | 317 strong T-cell epitopes found (high immunogenicity risk) |
| groEL | E.coli | FAIL | 504 strong T-cell epitopes found (high immunogenicity risk) |
| Insulin | Human | FAIL | 124 strong T-cell epitopes found (high immunogenicity risk) |
| HBB | Human | FAIL | 189 strong T-cell epitopes found (high immunogenicity risk) |
| EPO | Human | FAIL | 288 strong T-cell epitopes found (high immunogenicity risk) |
| GH1 | Human | FAIL | 313 strong T-cell epitopes found (high immunogenicity risk) |
| IFNA2 | Human | FAIL | 296 strong T-cell epitopes found (high immunogenicity risk) |
| Albumin | Human | FAIL | 914 strong T-cell epitopes found (high immunogenicity risk) |
| TDH3 | Yeast | FAIL | 272 strong T-cell epitopes found (high immunogenicity risk) |
| PGK1 | Yeast | FAIL | 265 strong T-cell epitopes found (high immunogenicity risk) |
| ADH1 | Yeast | FAIL | 375 strong T-cell epitopes found (high immunogenicity risk) |
| mIL2 | Mouse | FAIL | 238 strong T-cell epitopes found (high immunogenicity risk) |
| mIFNG | Mouse | FAIL | 363 strong T-cell epitopes found (high immunogenicity risk) |
| CHO_EPO | CHO | FAIL | 288 strong T-cell epitopes found (high immunogenicity risk) |
| CHO_GFP | CHO | FAIL | 356 strong T-cell epitopes found (high immunogenicity risk) |

### NoDominantBCellEpitope
Pass rate: 0.0% (0/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | FAIL | 11 dominant B-cell epitopes found (high immunogenicity risk) |
| mCherry | E.coli | FAIL | 4 dominant B-cell epitopes found (high immunogenicity risk) |
| T4_lysozyme | E.coli | FAIL | 9 dominant B-cell epitopes found (high immunogenicity risk) |
| groEL | E.coli | FAIL | 5 dominant B-cell epitopes found (high immunogenicity risk) |
| Insulin | Human | FAIL | 6 dominant B-cell epitopes found (high immunogenicity risk) |
| HBB | Human | FAIL | 13 dominant B-cell epitopes found (high immunogenicity risk) |
| EPO | Human | FAIL | 7 dominant B-cell epitopes found (high immunogenicity risk) |
| GH1 | Human | FAIL | 7 dominant B-cell epitopes found (high immunogenicity risk) |
| IFNA2 | Human | FAIL | 9 dominant B-cell epitopes found (high immunogenicity risk) |
| Albumin | Human | FAIL | 15 dominant B-cell epitopes found (high immunogenicity risk) |
| TDH3 | Yeast | FAIL | 13 dominant B-cell epitopes found (high immunogenicity risk) |
| PGK1 | Yeast | FAIL | 8 dominant B-cell epitopes found (high immunogenicity risk) |
| ADH1 | Yeast | FAIL | 10 dominant B-cell epitopes found (high immunogenicity risk) |
| mIL2 | Mouse | FAIL | 5 dominant B-cell epitopes found (high immunogenicity risk) |
| mIFNG | Mouse | FAIL | 15 dominant B-cell epitopes found (high immunogenicity risk) |
| CHO_EPO | CHO | FAIL | 7 dominant B-cell epitopes found (high immunogenicity risk) |
| CHO_GFP | CHO | FAIL | 11 dominant B-cell epitopes found (high immunogenicity risk) |

### PopulationCoverageSafe
Pass rate: 100.0% (17/17)

| Gene | Organism | Verdict | Details |
|------|----------|---------|---------|
| GFP | E.coli | PASS |  |
| mCherry | E.coli | PASS |  |
| T4_lysozyme | E.coli | PASS |  |
| groEL | E.coli | PASS |  |
| Insulin | Human | PASS |  |
| HBB | Human | PASS |  |
| EPO | Human | PASS |  |
| GH1 | Human | PASS |  |
| IFNA2 | Human | PASS |  |
| Albumin | Human | PASS |  |
| TDH3 | Yeast | PASS |  |
| PGK1 | Yeast | PASS |  |
| ADH1 | Yeast | PASS |  |
| mIL2 | Mouse | PASS |  |
| mIFNG | Mouse | PASS |  |
| CHO_EPO | CHO | PASS |  |
| CHO_GFP | CHO | PASS |  |

## Knowledge Gaps (Substitutes for Wet-Lab)

- Aggregation-prone regions identified from intrinsic sequence only.  Structural accessibility correction would refine the prediction.
- Cannot determine if unpaired Cys is buried without structural data.
- Charge composition assessed from sequence alone.  Surface accessibility of charged residues (from structure) would improve solubility prediction.
- Cysteine count is even but spatial pairing cannot be verified without structural data.
- Heuristic fallback used
- Hydrophobic fraction is within normal range, but core burial cannot be verified without structural data.
- Hydrophobic stretches assessed from sequence alone.  Structural context (buried vs. exposed) would clarify whether stretches contribute to aggregation.
- No PDB structure provided; solubility estimated from intrinsic sequence properties only.  Structural correction would improve accuracy.
- No PDB structure provided; stability estimated from sequence composition only.  Structural analysis would improve confidence.
