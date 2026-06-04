"""
Incremental Sequence State — O(1) constraint tracking for gene optimization.

Eliminates redundant full-sequence scans (_count_gts, CG counting) by maintaining
incremental counters that update in O(1) when codons are swapped.

Performance impact: Replaces O(N) scans at 49+ call sites with O(1) updates,
reducing total optimization time by an estimated 30-50%.
"""
from typing import Dict, List, Optional, Tuple, Set
from .type_system import CODON_TABLE, AA_TO_CODONS


class IncrementalSequenceState:
    """Maintains mutable sequence state with O(1) incremental constraint tracking.
    
    Instead of scanning the entire sequence to count GT/CG dinucleotides after
    every codon swap, this class tracks them incrementally. When a codon at
    position `codon_idx` (0-based codon index) is changed, only the boundary
    positions affected by that codon need to be rechecked — at most 4 dinucleotide
    positions: (codon_start-1, codon_start), (codon_start+1, codon_start+2),
    (codon_start+2, codon_start+3), and within-codon positions.
    
    Usage:
        state = IncrementalSequenceState("ATGGTCAAG...")
        # state.gt_count is immediately available
        # state.cg_count is immediately available
        
        # Swap codon at codon index 5 to "GCA"
        state.swap_codon(5, "GCA")
        # state.gt_count and state.cg_count are updated in O(1)
    """
    
    def __init__(self, sequence: str):
        """Initialize from a DNA sequence string.
        
        Pre-computes:
        - gt_count: number of GT dinucleotides
        - cg_count: number of CG dinucleotides  
        - _gt_positions: set of positions where GT occurs
        - _cg_positions: set of positions where CG occurs
        - _seq_list: mutable list of characters
        """
        # Store as mutable list for in-place swaps
        self._seq_list: List[str] = list(sequence)
        self._n: int = len(sequence)
        
        # Pre-compute dinucleotide positions
        self._gt_positions: Set[int] = set()
        self._cg_positions: Set[int] = set()
        
        for i in range(self._n - 1):
            if sequence[i] == 'G' and sequence[i+1] == 'T':
                self._gt_positions.add(i)
            elif sequence[i] == 'C' and sequence[i+1] == 'G':
                self._cg_positions.add(i)
        
        # Cached counts (derived from position sets)
        self.gt_count: int = len(self._gt_positions)
        self.cg_count: int = len(self._cg_positions)
        
        # Pre-compute amino acid for each codon position
        self._codon_aas: List[Optional[str]] = []
        for i in range(0, self._n - 2, 3):
            codon = sequence[i:i+3]
            self._codon_aas.append(CODON_TABLE.get(codon))
        
        # Cached sequence string — invalidated on swap, recomputed on access
        self._sequence_cache: Optional[str] = sequence
    
    @property
    def sequence(self) -> str:
        """Return current sequence as string (cached, O(1) on repeated access)."""
        if self._sequence_cache is None:
            self._sequence_cache = "".join(self._seq_list)
        return self._sequence_cache
    
    @property 
    def num_codons(self) -> int:
        """Number of codons in the sequence."""
        return self._n // 3
    
    def get_codon(self, codon_idx: int) -> str:
        """Get the codon at the given 0-based codon index."""
        start = codon_idx * 3
        return "".join(self._seq_list[start:start+3])
    
    def get_aa(self, codon_idx: int) -> Optional[str]:
        """Get the amino acid at the given codon index."""
        if 0 <= codon_idx < len(self._codon_aas):
            return self._codon_aas[codon_idx]
        return None
    
    def swap_codon(self, codon_idx: int, new_codon: str) -> str:
        """Swap the codon at codon_idx to new_codon, updating counters in O(1).
        
        Returns the OLD codon that was replaced.
        
        Only rechecks the dinucleotide positions affected by the swap:
        - Within the codon: positions codon_start and codon_start+1
        - Left boundary: position codon_start-1 (between prev codon and this one)
        - Right boundary: position codon_start+2 (between this codon and next one)
        
        Total: at most 4 dinucleotide positions to recheck, regardless of sequence length.
        """
        start = codon_idx * 3
        old_codon = "".join(self._seq_list[start:start+3])
        
        # Identify all dinucleotide positions affected by this codon
        # A dinucleotide at position i involves bases[i] and bases[i+1]
        # The codon at [start, start+1, start+2] affects dinucleotides at:
        #   start-1 (prev base + new base[0])
        #   start   (new base[0] + new base[1])  
        #   start+1 (new base[1] + new base[2])
        #   start+2 (new base[2] + next base)
        affected_positions = []
        for pos in [start - 1, start, start + 1, start + 2]:
            if 0 <= pos < self._n - 1:
                affected_positions.append(pos)
        
        # Remove old dinucleotide state at affected positions
        for pos in affected_positions:
            self._gt_positions.discard(pos)
            self._cg_positions.discard(pos)
        
        # Apply the swap
        for k, base in enumerate(new_codon):
            self._seq_list[start + k] = base
        
        # Recompute dinucleotide state at affected positions
        for pos in affected_positions:
            if self._seq_list[pos] == 'G' and self._seq_list[pos + 1] == 'T':
                self._gt_positions.add(pos)
            elif self._seq_list[pos] == 'C' and self._seq_list[pos + 1] == 'G':
                self._cg_positions.add(pos)
        
        # Update cached counts
        self.gt_count = len(self._gt_positions)
        self.cg_count = len(self._cg_positions)
        
        # Update amino acid cache
        if 0 <= codon_idx < len(self._codon_aas):
            self._codon_aas[codon_idx] = CODON_TABLE.get(new_codon)
        
        # Invalidate sequence cache
        self._sequence_cache = None
        
        return old_codon
    
    def try_swap_codon(self, codon_idx: int, new_codon: str, 
                        check_fn=None) -> Tuple[bool, str]:
        """Try swapping a codon, run a check, and rollback if check fails.
        
        This avoids creating seq_list[:] copies — we swap in-place, check,
        and swap back if the check fails.
        
        Args:
            codon_idx: 0-based codon index to swap
            new_codon: The new codon to try
            check_fn: Optional callable(self) -> bool. If it returns False, rollback.
            
        Returns:
            (success, old_codon): Whether the swap was accepted, and what the old codon was.
        """
        old_codon = self.swap_codon(codon_idx, new_codon)
        
        if check_fn is not None and not check_fn(self):
            # Rollback
            self.swap_codon(codon_idx, old_codon)
            return (False, old_codon)
        
        return (True, old_codon)
    
    def gt_positions_list(self) -> List[int]:
        """Return sorted list of GT dinucleotide positions."""
        return sorted(self._gt_positions)
    
    def cg_positions_list(self) -> List[int]:
        """Return sorted list of CG dinucleotide positions."""
        return sorted(self._cg_positions)
    
    def has_gt_at(self, pos: int) -> bool:
        """Check if there's a GT dinucleotide at position pos."""
        return pos in self._gt_positions
    
    def has_cg_at(self, pos: int) -> bool:
        """Check if there's a CG dinucleotide at position pos."""
        return pos in self._cg_positions
    
    def count_gts_in_region(self, start: int, end: int) -> int:
        """Count GT dinucleotides within [start, end) positions."""
        return sum(1 for p in self._gt_positions if start <= p < end)
    
    def count_cgs_in_region(self, start: int, end: int) -> int:
        """Count CG dinucleotides within [start, end) positions."""
        return sum(1 for p in self._cg_positions if start <= p < end)
    
    def boundary_creates_gt(self, codon_idx: int, new_codon: str) -> Tuple[bool, bool]:
        """Check if swapping codon at codon_idx to new_codon would create GT at boundaries.
        
        Returns (left_boundary_gt, right_boundary_gt).
        Does NOT actually perform the swap.
        """
        start = codon_idx * 3
        
        left_gt = False
        right_gt = False
        
        # Left boundary: seq[start-1] + new_codon[0]
        if start > 0 and self._seq_list[start - 1] == 'G' and new_codon[0] == 'T':
            left_gt = True
        
        # Right boundary: new_codon[-1] + seq[start+3]
        if start + 3 < self._n and new_codon[-1] == 'G' and self._seq_list[start + 3] == 'T':
            right_gt = True
        
        return (left_gt, right_gt)
    
    def boundary_creates_cg(self, codon_idx: int, new_codon: str) -> Tuple[bool, bool]:
        """Check if swapping codon at codon_idx to new_codon would create CG at boundaries.
        
        Returns (left_boundary_cg, right_boundary_cg).
        Does NOT actually perform the swap.
        """
        start = codon_idx * 3
        
        left_cg = False
        right_cg = False
        
        if start > 0 and self._seq_list[start - 1] == 'C' and new_codon[0] == 'G':
            left_cg = True
        
        if start + 3 < self._n and new_codon[-1] == 'C' and self._seq_list[start + 3] == 'G':
            right_cg = True
        
        return (left_cg, right_cg)
    
    def would_gt_increase(self, codon_idx: int, new_codon: str) -> bool:
        """Check if swapping to new_codon would increase GT count (without performing swap).
        
        Computes the net GT change by checking only the 4 affected dinucleotide positions.
        """
        start = codon_idx * 3
        old_codon = "".join(self._seq_list[start:start+3])
        
        # Count current GTs at affected positions
        old_gt_at_affected = 0
        for pos in [start - 1, start, start + 1, start + 2]:
            if 0 <= pos < self._n - 1 and pos in self._gt_positions:
                old_gt_at_affected += 1
        
        # Compute new GTs at affected positions
        # Build the 4 bases involved: prev_base, new[0], new[1], new[2], next_base
        prev_base = self._seq_list[start - 1] if start > 0 else ''
        next_base = self._seq_list[start + 3] if start + 3 < self._n else ''
        bases = [prev_base, new_codon[0], new_codon[1], new_codon[2], next_base]
        
        new_gt_at_affected = 0
        for i in range(len(bases) - 1):
            if bases[i] == 'G' and bases[i+1] == 'T':
                pos = start - 1 + i
                if 0 <= pos < self._n - 1:
                    new_gt_at_affected += 1
        
        return new_gt_at_affected > old_gt_at_affected
    
    def would_cg_increase(self, codon_idx: int, new_codon: str) -> bool:
        """Check if swapping to new_codon would increase CG count (without performing swap)."""
        start = codon_idx * 3
        old_codon = "".join(self._seq_list[start:start+3])
        
        old_cg_at_affected = 0
        for pos in [start - 1, start, start + 1, start + 2]:
            if 0 <= pos < self._n - 1 and pos in self._cg_positions:
                old_cg_at_affected += 1
        
        prev_base = self._seq_list[start - 1] if start > 0 else ''
        next_base = self._seq_list[start + 3] if start + 3 < self._n else ''
        bases = [prev_base, new_codon[0], new_codon[1], new_codon[2], next_base]
        
        new_cg_at_affected = 0
        for i in range(len(bases) - 1):
            if bases[i] == 'C' and bases[i+1] == 'G':
                pos = start - 1 + i
                if 0 <= pos < self._n - 1:
                    new_cg_at_affected += 1
        
        return new_cg_at_affected > old_cg_at_affected
    
    def new_gt_positions_after_swap(self, codon_idx: int, new_codon: str) -> Set[int]:
        """Return the set of NEW GT positions that would be created by the swap.
        
        Useful for checking if new GTs are unavoidable.
        """
        start = codon_idx * 3
        old_gt = set(self._gt_positions)
        
        # Simulate the swap
        old_codon = "".join(self._seq_list[start:start+3])
        
        # Remove affected positions from old set
        new_gt = set(old_gt)
        for pos in [start - 1, start, start + 1, start + 2]:
            if 0 <= pos < self._n - 1:
                new_gt.discard(pos)
        
        # Add new GTs at affected positions
        prev_base = self._seq_list[start - 1] if start > 0 else ''
        next_base = self._seq_list[start + 3] if start + 3 < self._n else ''
        bases = [prev_base, new_codon[0], new_codon[1], new_codon[2], next_base]
        
        for i in range(len(bases) - 1):
            if bases[i] == 'G' and bases[i+1] == 'T':
                pos = start - 1 + i
                if 0 <= pos < self._n - 1:
                    new_gt.add(pos)
        
        return new_gt - old_gt


# ────────────────────────────────────────────────────────────
# Pre-computed Codon Cache — avoids repeated sorting in hot paths
# ────────────────────────────────────────────────────────────

class CodonCache:
    """Pre-computed codon information sorted by CAI weight.
    
    Instead of calling sorted(AA_TO_CODONS[aa], key=..., reverse=True) inside
    every iteration of the hot loop, this cache computes it once and reuses.
    
    Usage:
        cache = CodonCache(species_cai_weights)
        sorted_codons = cache.get_sorted_codons('L')  # ['CTG', 'CTC', 'CTT', ...]
    """
    
    def __init__(self, species_cai: Dict[str, float]):
        self._species_cai = species_cai
        self._sorted: Dict[str, List[str]] = {}
        self._best: Dict[str, str] = {}
        self._gt_free: Dict[str, List[str]] = {}
        self._cg_free: Dict[str, List[str]] = {}
        
        for aa, codons in AA_TO_CODONS.items():
            # Sort by CAI, highest first
            sorted_codons = sorted(codons, key=lambda c: species_cai.get(c, 0.0), reverse=True)
            self._sorted[aa] = sorted_codons
            self._best[aa] = sorted_codons[0] if sorted_codons else ""
            
            # Pre-compute GT-free and CG-free subsets
            self._gt_free[aa] = [c for c in sorted_codons if "GT" not in c]
            self._cg_free[aa] = [c for c in sorted_codons if "CG" not in c]
    
    def get_sorted_codons(self, aa: str) -> List[str]:
        """Get codons for an amino acid, sorted by CAI (highest first)."""
        return self._sorted.get(aa, [])
    
    def get_best_codon(self, aa: str) -> str:
        """Get the highest-CAI codon for an amino acid."""
        return self._best.get(aa, "")
    
    def get_cai(self, codon: str) -> float:
        """Get the CAI weight for a codon."""
        return self._species_cai.get(codon, 0.0)
    
    def get_gt_free_codons(self, aa: str) -> List[str]:
        """Get GT-free codons for an amino acid, sorted by CAI (highest first)."""
        return self._gt_free.get(aa, [])
    
    def get_cg_free_codons(self, aa: str) -> List[str]:
        """Get CG-free codons for an amino acid, sorted by CAI (highest first)."""
        return self._cg_free.get(aa, [])


class EnzymeSiteCache:
    """Pre-computed enzyme→recognition_site mapping.
    
    Avoids calling get_recognition_site() inside inner loops.
    """
    
    def __init__(self, enzymes: List[str]):
        from .restriction_sites import get_recognition_site
        self._sites: Dict[str, Optional[str]] = {}
        self._site_lengths: Dict[str, int] = {}
        
        for enzyme in enzymes:
            site = get_recognition_site(enzyme)
            self._sites[enzyme] = site
            if site is not None:
                self._site_lengths[enzyme] = len(site)
    
    def get_site(self, enzyme: str) -> Optional[str]:
        """Get the recognition site for an enzyme."""
        return self._sites.get(enzyme)
    
    def get_site_length(self, enzyme: str) -> int:
        """Get the length of the recognition site for an enzyme."""
        return self._site_lengths.get(enzyme, 0)
    
    def check_any_site_present(self, sequence: str) -> bool:
        """Check if any enzyme's recognition site is present in the sequence."""
        for enzyme, site in self._sites.items():
            if site is not None and site in sequence:
                return True
        return False
    
    def find_sites(self, sequence: str) -> List[Tuple[str, int]]:
        """Find all enzyme sites in the sequence. Returns [(enzyme, position), ...]."""
        results = []
        for enzyme, site in self._sites.items():
            if site is None:
                continue
            pos = 0
            while True:
                pos = sequence.find(site, pos)
                if pos == -1:
                    break
                results.append((enzyme, pos))
                pos += 1
        return results
