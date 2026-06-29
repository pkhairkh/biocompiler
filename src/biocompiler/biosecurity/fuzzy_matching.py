"""
Fuzzy matching utilities for biosecurity screening.

Includes reverse complement, Hamming distance, Levenshtein distance,
and fuzzy match search functions.
"""

from __future__ import annotations

# ─────────────────────────────────────────────────────────────────────────────
# Sequence utilities
# ─────────────────────────────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

_KMER_CACHE: dict = {}

# prior fix / H17: k-mer pre-filter soundness threshold.
# The pigeonhole bound guarantees that two strings of length n with
# Hamming/Levenshtein distance <= d share at least one k-mer of size
# floor(n / (d+1)) (Hamming) or floor((n-d) / (d+1)) (Levenshtein).
# For the default k=4 with max_distance=2, the bound fails on short
# motifs: a length-10 motif under max_distance=2 has floor(10/3) = 3
# (Hamming) guaranteed matching k-mers, so k=4 is UNSOUND and produces
# ~6.7% false negatives (e.g. 2-substitution patterns at positions
# 3,6 or 3,7 break all 7 of the motif's 4-mers).  To stay sound for
# both Hamming and Levenshtein, skip the pre-filter entirely whenever
# the conservative bound k * (max_distance + 1) > len(motif) holds.
# Concretely: with k=4, d=2 this means len(motif) < 12 -> skip.
def _sound_kmer_filter(motif_len: int, k: int, max_distance: int) -> bool:
    """Return True iff the k-mer pre-filter is provably sound for this motif.

    Soundness condition (conservative, valid for both Hamming and
    Levenshtein distance): k * (max_distance + 1) <= motif_len.
    """
    if max_distance < 1:
        # max_distance == 0 means exact match only; any k works because
        # the pre-filter is a strict subset check (no false negatives).
        return True
    return k * (max_distance + 1) <= motif_len
# ─── NUMBA acceleration ─────────────────────────────────────────────
try:
    import numpy as np
    from biocompiler.optimizer.numba_kernels import (
        HAS_NUMBA as _HAS_NUMBA,
        _numba_fuzzy_match_hamming as _numba_fuzzy,
        _numba_hamming_distance as _numba_hamming,
        _numba_levenshtein_within as _numba_lev_within,
        _numba_levenshtein_distance as _numba_lev_dist,
        _numba_has_shared_kmer_fast as _numba_kmer_fast,
    )
    import numpy as _np
except (ImportError, ModuleNotFoundError):
    _HAS_NUMBA = False
    _numba_fuzzy = None
    _numba_hamming = None
    _numba_lev_within = None
    _numba_lev_dist = None
    _numba_kmer_fast = None


def _get_motif_kmers(
    motif: str, k: int = 4, max_distance: int = 0
) -> set | None:
    """Return the set of k-mers for *motif*, or ``None`` to skip the filter.

    Parameters
    ----------
    motif : str
        The motif whose k-mers are wanted.
    k : int
        K-mer length (default 4, matching the historical pre-filter).
    max_distance : int
        Maximum fuzzy distance the caller will accept.  When the
        pigeonhole soundness bound ``k * (max_distance + 1) <= len(motif)``
        is violated, this function returns ``None`` so the caller falls
        back to a direct (unfiltered) Hamming/Levenshtein comparison —
        see :func:`_sound_kmer_filter` and the H17 fix note above.

    Returns
    -------
    set[str] | None
        Set of k-mers, or ``None`` if the pre-filter should be skipped
        (motif shorter than k, or soundness bound violated).
    """
    # Skip pre-filter when soundness bound fails (H17) or motif is too
    # short to extract any k-mers of size k.
    if not _sound_kmer_filter(len(motif), k, max_distance):
        return None
    key = (motif, k)
    if key not in _KMER_CACHE:
        if len(motif) >= k:
            _KMER_CACHE[key] = {motif[j:j+k] for j in range(len(motif) - k + 1)}
        else:
            _KMER_CACHE[key] = None
    return _KMER_CACHE[key]


# ─── Numba k-mer hash cache ─────────────────────────────────────────
# When numba is available, we pre-compute a sorted numpy array of k-mer
# hashes for each motif. This allows O(log n) binary search instead of
# O(n) set lookup in the hot loop.

class _CachedKmerSet:
    """Wraps a set of k-mers with a pre-computed sorted hash array for numba."""
    def __init__(self, kmer_set, k=4):
        self._set = kmer_set
        self._k = k
        if _HAS_NUMBA and kmer_set is not None:
            import numpy as _np_mod
            hashes = []
            for km in kmer_set:
                h = 0
                for ch in km:
                    h = (h * 2654435761 + ord(ch)) % 9223372036854775807  # keep within int64
                hashes.append(h)
            self._hash_array = _np_mod.array(sorted(hashes), dtype=_np_mod.int64) if hashes else _np_mod.zeros(0, dtype=_np_mod.int64)
        else:
            self._hash_array = None
    
    def __contains__(self, km):
        return km in self._set
    
    def __bool__(self):
        return bool(self._set)
    
    def __len__(self):
        return len(self._set)


_kmer_cache: dict = {}

def _get_cached_kmers(motif: str, k: int = 4, max_distance: int = 0):
    """Get cached k-mer set with pre-computed hash array."""
    cache_key = (motif, k, max_distance)
    if cache_key not in _kmer_cache:
        kmers = _get_motif_kmers(motif, k, max_distance)
        if kmers is None:
            _kmer_cache[cache_key] = None
        else:
            _kmer_cache[cache_key] = _CachedKmerSet(kmers, k)
    return _kmer_cache[cache_key]

def _has_shared_kmer(window: str, motif_kmers) -> bool:
    if motif_kmers is None or len(window) < 4:
        return True
    # NUMBA fast path: use pre-computed sorted hash array
    if _HAS_NUMBA and _numba_kmer_fast is not None and hasattr(motif_kmers, '_hash_array') and motif_kmers._hash_array is not None:
        try:
            return _numba_kmer_fast(window.encode('ascii'), motif_kmers._hash_array, 4)
        except Exception:
            pass
    # Pure Python fallback
    if hasattr(motif_kmers, '_set'):
        kmers = motif_kmers._set
    else:
        kmers = motif_kmers
    for j in range(len(window) - 3):
        if window[j:j+4] in kmers:
            return True
    return False

def _build_seq_kmer_index(sequence: str, k: int = 4) -> set:
    """Build a set of all k-mers in the sequence. O(n) once."""
    if len(sequence) < k:
        return set()
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def _motif_could_match(motif_kmers: set | None, seq_kmers: set) -> bool:
    """Quick check: does ANY motif k-mer appear ANYWHERE in the sequence?
    If not, no window can match — skip the motif entirely.
    """
    if motif_kmers is None:
        return True  # can't filter
    # Check if any motif k-mer is in the sequence's k-mer set
    return not motif_kmers.isdisjoint(seq_kmers)


def reverse_complement(dna: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Parameters
    ----------
    dna : str
        DNA sequence consisting of A, C, G, T (uppercase expected).

    Returns
    -------
    str
        The reverse complement.
    """
    return dna.translate(_COMPLEMENT)[::-1]


def _hamming_distance(s1: str, s2: str) -> int:
    """Compute the Hamming distance between two equal-length strings.

    Raises ``ValueError`` if the strings differ in length.
    """
    if len(s1) != len(s2):
        raise ValueError(
            f"Hamming distance requires equal-length strings, "
            f"got {len(s1)} and {len(s2)}"
        )
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def _levenshtein_distance(s1: str, s2: str) -> int:
    """Compute Levenshtein distance. Uses numba when available."""
    if _HAS_NUMBA and _numba_lev_dist is not None:
        try:
            return _numba_lev_dist(s1.encode('ascii'), s2.encode('ascii'))
        except Exception:
            pass
    if len(s1) < len(s2):
        s1, s2 = s2, s1
    previous = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous[j + 1] + 1
            deletions = current[j] + 1
            substitutions = previous[j] + (c1 != c2)
            current.append(min(insertions, deletions, substitutions))
        previous = current

    return previous[-1]


def _levenshtein_within(s1: str, s2: str, max_dist: int) -> bool:
    """Check if Levenshtein distance(s1, s2) <= max_dist. Uses numba when available."""
    if _HAS_NUMBA and _numba_lev_within is not None:
        try:
            return _numba_lev_within(s1.encode('ascii'), s2.encode('ascii'), max_dist)
        except Exception:
            pass
    if abs(len(s1) - len(s2)) > max_dist:
        return False
    if len(s1) < len(s2):
        s1, s2 = s2, s1
    previous = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current = [i + 1]
        row_min = i + 1
        for j, c2 in enumerate(s2):
            insertions = previous[j + 1] + 1
            deletions = current[j] + 1
            substitutions = previous[j] + (c1 != c2)
            val = min(insertions, deletions, substitutions)
            current.append(val)
            if val < row_min:
                row_min = val
        # Early exit: if the minimum value in this row exceeds max_dist,
        # the final distance will also exceed max_dist
        if row_min > max_dist:
            return False
        previous = current

    return previous[-1] <= max_dist


def _fuzzy_match_hamming(
    sequence: str,
    motif: str,
    max_distance: int = 2,
) -> list[tuple[int, int, list[tuple[int, str, str]]]]:
    """Find fuzzy matches using Hamming distance. Uses numba JIT when available."""
    """Find fuzzy matches of *motif* in *sequence* using Hamming distance.

    Only returns matches with distance >= 1 (i.e. excludes exact matches)
    and distance <= *max_distance*.

    Parameters
    ----------
    sequence : str
        The sequence to search in.
    motif : str
        The motif to search for.
    max_distance : int
        Maximum Hamming distance to report (default 2).

    Returns
    -------
    list of (position, distance, substitutions)
        Each element is a tuple ``(pos, dist, subs)`` where *subs* is a
        list of ``(position_in_window, original_char, replacement_char)``
        tuples describing the substitutions.
    """
    mlen = len(motif)
    slen = len(sequence)
    if mlen == 0 or slen < mlen:
        return []

    # ─── NUMBA fast path ───
    if _HAS_NUMBA and _numba_fuzzy is not None:
        try:
            seq_bytes = sequence.encode('ascii')
            motif_bytes = motif.encode('ascii')
            matches = _numba_fuzzy(seq_bytes, motif_bytes, max_distance, 4)
            results = []
            for row in matches:
                pos = int(row[0])
                dist = int(row[1])
                window = sequence[pos:pos + mlen]
                subs = [
                    (j, motif[j], window[j])
                    for j in range(mlen)
                    if motif[j] != window[j]
                ]
                results.append((pos, dist, subs))
            return results
        except Exception:
            pass  # fall through to pure Python

    # ─── Pure Python fallback ───
    motif_kmers = _get_cached_kmers(motif, max_distance=max_distance)
    results: list[tuple[int, int, list[tuple[int, str, str]]]] = []
    for i in range(slen - mlen + 1):
        window = sequence[i : i + mlen]
        if not _has_shared_kmer(window, motif_kmers):
            continue
        dist = _hamming_distance(window, motif)
        if 1 <= dist <= max_distance:
            subs = [
                (j, motif[j], window[j])
                for j in range(mlen)
                if motif[j] != window[j]
            ]
            results.append((i, dist, subs))

    return results


def _fuzzy_match_edit_distance(
    sequence: str,
    motif: str,
    max_distance: int = 1,
) -> list[tuple[int, int]]:
    """Find fuzzy matches of *motif* in *sequence* using Levenshtein distance.

    Uses a sliding-window approach with windows of varying length around
    the motif length to catch insertions and deletions.

    Only returns matches with distance >= 1 and <= *max_distance*.

    Parameters
    ----------
    sequence : str
        The sequence to search in.
    motif : str
        The motif to search for.
    max_distance : int
        Maximum edit distance to report (default 1).

    Returns
    -------
    list of (position, distance)
        Each element is a tuple ``(pos, dist)``.
    """
    mlen = len(motif)
    slen = len(sequence)
    if mlen == 0 or slen == 0:
        return []

    results: list[tuple[int, int]] = []
    # k-mer pre-filter: pre-compute the set of 4-mers in the motif.
    # If a window shares NO 4-mer with the motif, the edit distance is
    # guaranteed > max_distance (for max_distance <= 1, k = 4), so we skip
    # the expensive Levenshtein computation entirely.
    # This reduces Levenshtein calls by ~90% in practice.
    #
    # prior fix / H17: pass max_distance through to _get_motif_kmers so
    # the soundness check can disable the pre-filter for short motifs
    # (pigeonhole bound k*(d+1) > len(motif) violated -> pre-filter
    # would produce false negatives).  In that case the function
    # returns None and _has_shared_kmer falls back to True (no filter),
    # giving a correct but slower direct Levenshtein comparison.
    motif_kmers = _get_cached_kmers(motif, max_distance=max_distance)
    for window_len in range(max(1, mlen - max_distance), mlen + max_distance + 1):
        if window_len > slen:
            continue
        for i in range(slen - window_len + 1):
            window = sequence[i : i + window_len]
            if not _has_shared_kmer(window, motif_kmers):
                continue
            if _levenshtein_within(window, motif, max_distance):
                dist = _levenshtein_distance(window, motif)
                if 1 <= dist <= max_distance:
                    results.append((i, dist))

    # Deduplicate: keep the best distance per position
    best: dict[int, int] = {}
    for pos, dist in results:
        if pos not in best or dist < best[pos]:
            best[pos] = dist

    return sorted(best.items())
