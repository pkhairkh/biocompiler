"""
Mouse (Mus musculus) Codon Usage Data

Source: Kazusa Codon Usage Database
59,432 CDSs, 24,040,270 codons
Coding GC: 52.25%
"""

# Format: {codon: (amino_acid, fraction, per_thousand, count)}
MOUSE_CODON_USAGE: dict[str, tuple[str, float, float, int]] = {
    "TTT": ("F", 0.46, 17.5, 420837),
    "TTC": ("F", 0.54, 20.6, 494997),
    "TTA": ("L", 0.07, 7.6, 183106),
    "TTG": ("L", 0.13, 13.1, 315351),
    "CTT": ("L", 0.13, 12.6, 302711),
    "CTC": ("L", 0.22, 21.0, 504980),
    "CTA": ("L", 0.07, 7.2, 173960),
    "CTG": ("L", 0.38, 39.6, 952185),
    "ATT": ("I", 0.36, 16.0, 384673),
    "ATC": ("I", 0.48, 21.4, 514619),
    "ATA": ("I", 0.16, 7.3, 174852),
    "ATG": ("M", 1.00, 22.3, 535014),
    "GTT": ("V", 0.19, 11.7, 280603),
    "GTC": ("V", 0.24, 14.9, 358478),
    "GTA": ("V", 0.12, 7.2, 173960),
    "GTG": ("V", 0.46, 28.6, 687491),
    "TCT": ("S", 0.18, 14.6, 351289),
    "TCC": ("S", 0.23, 18.3, 440468),
    "TCA": ("S", 0.15, 12.1, 290522),
    "TCG": ("S", 0.05, 4.3, 103373),
    "CCT": ("P", 0.28, 17.4, 418454),
    "CCC": ("P", 0.33, 20.3, 488352),
    "CCA": ("P", 0.28, 16.9, 406531),
    "CCG": ("P", 0.11, 6.9, 165877),
    "ACT": ("T", 0.25, 13.2, 317071),
    "ACC": ("T", 0.37, 19.4, 467021),
    "ACA": ("T", 0.27, 14.4, 346178),
    "ACG": ("T", 0.10, 5.5, 131823),
    "GCT": ("A", 0.26, 18.4, 442344),
    "GCC": ("A", 0.41, 29.1, 699443),
    "GCA": ("A", 0.22, 15.8, 379664),
    "GCG": ("A", 0.11, 7.6, 182868),
    "TAT": ("Y", 0.45, 12.3, 295743),
    "TAC": ("Y", 0.55, 15.1, 362810),
    "TAA": ("*", 0.30, 1.0, 24040),
    "TAG": ("*", 0.23, 0.8, 19232),
    "CAT": ("H", 0.42, 10.8, 259636),
    "CAC": ("H", 0.58, 14.8, 356197),
    "CAA": ("Q", 0.26, 12.3, 295743),
    "CAG": ("Q", 0.74, 34.9, 839005),
    "AAT": ("N", 0.47, 17.5, 420717),
    "AAC": ("N", 0.53, 19.5, 469287),
    "AAA": ("K", 0.43, 24.2, 581773),
    "AAG": ("K", 0.57, 32.0, 769289),
    "GAT": ("D", 0.46, 22.0, 528886),
    "GAC": ("D", 0.54, 26.2, 629856),
    "GAA": ("E", 0.42, 29.4, 706784),
    "GAG": ("E", 0.58, 40.8, 980430),
    "TGT": ("C", 0.45, 10.5, 252422),
    "TGC": ("C", 0.55, 12.8, 307715),
    "TGA": ("*", 0.47, 1.7, 40868),
    "TGG": ("W", 1.00, 13.3, 319735),
    "CGT": ("R", 0.09, 4.9, 117798),
    "CGC": ("R", 0.19, 11.3, 271656),
    "CGA": ("R", 0.10, 6.3, 151454),
    "CGG": ("R", 0.20, 12.0, 288483),
    "AGT": ("S", 0.16, 12.3, 295743),
    "AGC": ("S", 0.24, 19.1, 459168),
    "AGA": ("R", 0.22, 12.1, 290522),
    "AGG": ("R", 0.21, 12.0, 288483),
    "GGT": ("G", 0.16, 11.0, 264443),
    "GGC": ("G", 0.35, 23.5, 564944),
    "GGA": ("G", 0.25, 16.8, 403876),
    "GGG": ("G", 0.24, 16.3, 391856),
}

# Compute relative adaptiveness
_AA_MAX_FREQ: dict[str, float] = {}
for _codon, (_aa, _frac, _freq, _count) in MOUSE_CODON_USAGE.items():
    if _aa != "*":
        _current = _AA_MAX_FREQ.get(_aa, 0.0)
        if _freq > _current:
            _AA_MAX_FREQ[_aa] = _freq

MOUSE_CODON_ADAPTIVENESS: dict[str, float] = {}
for _codon, (_aa, _frac, _freq, _count) in MOUSE_CODON_USAGE.items():
    if _aa != "*":
        MOUSE_CODON_ADAPTIVENESS[_codon] = _freq / _AA_MAX_FREQ[_aa] if _AA_MAX_FREQ[_aa] > 0 else 0.0

MOUSE_PREFERRED_CODONS: dict[str, str] = {}
_AA_CODONS: dict[str, list[tuple[str, float]]] = {}
for _codon, (_aa, _frac, _freq, _count) in MOUSE_CODON_USAGE.items():
    if _aa != "*":
        _AA_CODONS.setdefault(_aa, []).append((_codon, _freq))

for _aa, _codons in _AA_CODONS.items():
    MOUSE_PREFERRED_CODONS[_aa] = max(_codons, key=lambda x: x[1])[0]
