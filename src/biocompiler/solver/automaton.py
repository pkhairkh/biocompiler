"""
Deterministic Finite Automata (DFA) for forbidden DNA patterns.

This module constructs DFAs used by the OR-Tools CP-SAT solver's
``AddAutomaton`` constraint to enforce "no forbidden substring" constraints
across codon boundaries — the KEY cross-codon constraint in the BioCompiler
CSP solver.

The core idea: a restriction site like EcoRI (GAATTC) must not appear anywhere
in the designed DNA sequence, even spanning a codon boundary.  OR-Tools'
``AddAutomaton`` natively encodes "the entire sequence is accepted by this DFA",
so we build a DFA that *accepts* exactly those strings that do **not** contain
any forbidden pattern, and feed it to the solver.

Alphabet encoding (matches ``constants.py``)::

    {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

Key design decisions
--------------------
* This module stores DFAs as ``transition_table[state][letter] = next_state``
  for easy simulation and inspection.  Use :func:`dfa_to_ortools_format` to
  convert to the ``(initial_state, triple_list, final_states)`` format that
  OR-Tools ``CpModel.add_automaton`` requires.
* A state is *accepting* if the DFA can legally end there — i.e. the
  forbidden pattern has **not** been seen.
* The *forbidden* (pattern-matched) state is a **trap**: once entered, the
  automaton never leaves, ensuring any string containing the pattern is
  rejected regardless of what follows.
* For single patterns we use a KMP-style failure function; for multiple
  patterns we use full Aho-Corasick construction, which is more efficient
  than intersecting individual DFAs.

Public API
----------
build_forbidden_pattern_dfa   — DFA for a single forbidden pattern
build_composite_dfa           — DFA for multiple forbidden patterns
build_reverse_complement_dfa  — DFA that also forbids the reverse complement
build_trun_dfa                — DFA that forbids long poly-T runs
dfa_accepts                   — test whether a DFA accepts a sequence
dfa_to_dot                    — export DFA as Graphviz DOT string
dfa_to_ortools_format         — convert to OR-Tools triple-list format
negate_dfa                    — swap accepting / non-accepting states
"""

from __future__ import annotations

import logging
from collections import defaultdict, deque

from biocompiler.shared.constants import BASE_MAP, BASE_REV, reverse_complement

__all__ = [
    "build_forbidden_pattern_dfa",
    "build_composite_dfa",
    "build_reverse_complement_dfa",
    "build_trun_dfa",
    "dfa_accepts",
    "dfa_to_dot",
    "dfa_to_ortools_format",
    "negate_dfa",
]

_logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Internal constants
# ---------------------------------------------------------------------------

_ALPHABET_SIZE = 4  # A=0, C=1, G=2, T=3
_INITIAL_STATE = 0  # Start state for all DFAs in this module
_UNDEFINED_GOTO = -1  # Placeholder for undefined trie edges in Aho-Corasick
_DEFAULT_MAX_T_RUN = 5  # Default maximum allowed consecutive T bases
_T_BASE_IDX: int = BASE_REV["T"]  # Numeric index of T in the alphabet


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _encode_pattern(pattern: str) -> list[int]:
    """Encode a DNA pattern string as a list of integers using ``BASE_REV``.

    Parameters
    ----------
    pattern:
        Uppercase DNA string over {A, C, G, T}.

    Returns
    -------
    list[int]
        Encoded pattern.

    Raises
    ------
    ValueError
        If *pattern* contains characters not in ``BASE_REV``.
    """
    try:
        return [BASE_REV[c] for c in pattern.upper()]
    except KeyError as exc:
        _logger.debug("Invalid base in pattern '%s': %s", pattern, exc)
        raise ValueError(
            f"Invalid base '{exc.args[0]}' in pattern '{pattern}'. "
            f"Allowed: A, C, G, T."
        ) from None


def _compute_kmp_failure(encoded: list[int]) -> list[int]:
    """Compute the KMP failure (partial-match) function.

    ``fail[i]`` is the length of the longest *proper* prefix of
    ``encoded[0:i]`` that is also a suffix of ``encoded[0:i]``.

    Returns a list of length ``len(encoded) + 1`` with ``fail[0] = 0``.
    """
    m: int = len(encoded)
    fail: list[int] = [0] * (m + 1)
    for i in range(2, m + 1):
        k = fail[i - 1]
        while k > 0 and encoded[k] != encoded[i - 1]:
            k = fail[k]
        if encoded[k] == encoded[i - 1]:
            fail[i] = k + 1
    return fail


# ---------------------------------------------------------------------------
# 1. Single-pattern DFA (KMP-style)
# ---------------------------------------------------------------------------

def build_forbidden_pattern_dfa(
    pattern: str,
) -> tuple[list[list[int]], list[int]]:
    """Build a DFA that **accepts** strings *not* containing *pattern*.

    Uses a KMP-style failure function to construct the transition table.
    The last state (``len(pattern)``) is the *forbidden trap* state:
    once entered, the automaton stays there forever.

    Parameters
    ----------
    pattern:
        DNA pattern to forbid (e.g. ``"GAATTC"`` for EcoRI).

    Returns
    -------
    transition_table : list[list[int]]
        ``transition_table[state][letter] = next_state``.
    accepting_states : list[int]
        States where the DFA may legally end (all states except the trap).

    Raises
    ------
    ValueError
        If *pattern* contains characters not in {A, C, G, T}.

    Examples
    --------
    >>> trans, accept = build_forbidden_pattern_dfa("GAATTC")
    >>> dfa_accepts((trans, accept), "GAATTC")
    False
    >>> dfa_accepts((trans, accept), "GAATTG")
    True
    """
    encoded = _encode_pattern(pattern)
    m = len(encoded)

    if m == 0:
        # Nothing is forbidden — accept everything with a single state.
        return [[_INITIAL_STATE] * _ALPHABET_SIZE], [_INITIAL_STATE]

    # KMP failure function
    fail = _compute_kmp_failure(encoded)

    # Build transition table for states 0 .. m
    num_states: int = m + 1
    delta: list[list[int]] = [[_INITIAL_STATE] * _ALPHABET_SIZE for _ in range(num_states)]

    for state in range(m):
        for c in range(_ALPHABET_SIZE):
            if c == encoded[state]:
                # Extend the current match by one character.
                delta[state][c] = state + 1
            else:
                # Follow failure links to find the longest prefix of
                # the pattern that could serve as the current match
                # after seeing character c.
                k: int = state
                while k > 0 and c != encoded[k]:
                    k = fail[k]
                if c == encoded[k]:
                    delta[state][c] = k + 1
                else:
                    delta[state][c] = _INITIAL_STATE

    # Trap state: once the full pattern is matched, stay trapped.
    for c in range(_ALPHABET_SIZE):
        delta[m][c] = m

    # Accepting states: everything except the trap.
    accepting: list[int] = list(range(m))
    return delta, accepting


# ---------------------------------------------------------------------------
# 2. Multi-pattern DFA (Aho-Corasick)
# ---------------------------------------------------------------------------

def build_composite_dfa(
    patterns: list[str],
) -> tuple[list[list[int]], list[int]]:
    """Build a DFA that rejects strings containing **any** of *patterns*.

    Uses Aho-Corasick construction for multiple patterns simultaneously,
    which is more efficient than intersecting individual DFAs.

    Parameters
    ----------
    patterns:
        List of DNA pattern strings to forbid.

    Returns
    -------
    transition_table : list[list[int]]
        ``transition_table[state][letter] = next_state``.
    accepting_states : list[int]
        States where the DFA may legally end.

    Raises
    ------
    ValueError
        If any pattern contains characters not in {A, C, G, T}.

    Notes
    -----
    The Aho-Corasick automaton is built in three phases:

    1. **Trie construction** — insert all patterns into a trie.
    2. **Failure-link computation** — BFS from the root to compute
       failure links and the full transition (delta) table.
    3. **Output propagation** — a state is an *output* state if it
       or any state in its failure chain is the end of a pattern.
       Output states become trap (non-accepting) states.

    Examples
    --------
    >>> trans, accept = build_composite_dfa(["GAATTC", "GGATCC"])
    >>> dfa_accepts((trans, accept), "GAATTC")
    False
    >>> dfa_accepts((trans, accept), "GGATCC")
    False
    >>> dfa_accepts((trans, accept), "GAATTG")
    True
    """
    if not patterns:
        # No patterns to forbid — accept everything.
        return [[_INITIAL_STATE] * _ALPHABET_SIZE], [_INITIAL_STATE]

    # Remove empty patterns (they automatically match everything if kept).
    patterns = [p for p in patterns if p]
    if not patterns:
        return [[_INITIAL_STATE] * _ALPHABET_SIZE], [_INITIAL_STATE]

    # Encode all patterns.
    encoded_patterns = [_encode_pattern(p) for p in patterns]

    # ---- Phase 1: Build trie -------------------------------------------
    # goto[s][c] = next state, or _UNDEFINED_GOTO if undefined.
    goto: list[list[int]] = [[_UNDEFINED_GOTO] * _ALPHABET_SIZE]
    direct_output: set[int] = set()  # states that are the end of a pattern
    num_states: int = 1  # root = _INITIAL_STATE

    for encoded in encoded_patterns:
        current: int = _INITIAL_STATE
        for c in encoded:
            if goto[current][c] == _UNDEFINED_GOTO:
                goto[current][c] = num_states
                goto.append([_UNDEFINED_GOTO] * _ALPHABET_SIZE)
                num_states += 1
            current = goto[current][c]
        direct_output.add(current)

    # ---- Phase 2: BFS to compute failure links & delta ------------------
    fail: list[int] = [_INITIAL_STATE] * num_states
    delta: list[list[int]] = [[_INITIAL_STATE] * _ALPHABET_SIZE for _ in range(num_states)]
    queue: deque[int] = deque()

    # Root transitions: undefined edges loop back to root.
    for c in range(_ALPHABET_SIZE):
        if goto[_INITIAL_STATE][c] != _UNDEFINED_GOTO:
            child: int = goto[_INITIAL_STATE][c]
            fail[child] = _INITIAL_STATE
            delta[_INITIAL_STATE][c] = child
            queue.append(child)
        else:
            delta[_INITIAL_STATE][c] = _INITIAL_STATE

    # BFS: for each trie node, compute failure links and fill delta.
    while queue:
        state: int = queue.popleft()
        for c in range(_ALPHABET_SIZE):
            if goto[state][c] != _UNDEFINED_GOTO:
                child = goto[state][c]
                # Failure link = parent's failure transition on c.
                fail[child] = delta[fail[state]][c]
                delta[state][c] = child
                queue.append(child)
            else:
                # No trie edge: delegate to failure-link's delta.
                delta[state][c] = delta[fail[state]][c]

    # ---- Phase 3: Propagate outputs through failure chain ---------------
    # A state is an *output* if it or any state reachable via failure links
    # from it is a direct-output state.  This ensures that e.g. state
    # representing "ATG" is marked as output when "G" is a forbidden pattern
    # (since fail["ATG"] → state for "G").
    output: set[int] = set(direct_output)

    # Iterate until fixed-point.  Each pass can only add states, so this
    # converges in at most num_states passes.
    changed: bool = True
    while changed:
        changed = False
        for state in range(1, num_states):
            if state not in output and fail[state] in output:
                output.add(state)
                changed = True

    # ---- Phase 4: Make output states trap states ------------------------
    for state in output:
        for c in range(_ALPHABET_SIZE):
            delta[state][c] = state

    # ---- Phase 5: Accepting = non-output states -------------------------
    accepting: list[int] = [s for s in range(num_states) if s not in output]

    return delta, accepting


# ---------------------------------------------------------------------------
# 3. Reverse-complement DFA
# ---------------------------------------------------------------------------

def build_reverse_complement_dfa(
    pattern: str,
) -> tuple[list[list[int]], list[int]]:
    """Build a DFA that rejects both *pattern* and its reverse complement.

    For **palindromic** patterns (where the pattern equals its own reverse
    complement, like the EcoRI site ``GAATTC``), only one direction is
    needed and ``build_forbidden_pattern_dfa`` is used directly.

    Parameters
    ----------
    pattern:
        DNA pattern string (e.g. ``"GAATTC"``).

    Returns
    -------
    transition_table : list[list[int]]
        ``transition_table[state][letter] = next_state``.
    accepting_states : list[int]
        States where the DFA may legally end.

    Examples
    --------
    >>> # EcoRI is palindromic: GAATTC == revcomp(GAATTC)
    >>> trans, acc = build_reverse_complement_dfa("GAATTC")
    >>> dfa_accepts((trans, acc), "GAATTC")
    False

    >>> # BamHI: GGATCC, revcomp = GGATCC (also palindromic)
    >>> trans, acc = build_reverse_complement_dfa("GGATCC")
    >>> dfa_accepts((trans, acc), "GGATCC")
    False

    >>> # Non-palindromic example
    >>> trans, acc = build_reverse_complement_dfa("AAGCTT")
    >>> dfa_accepts((trans, acc), "AAGCTT")   # forward
    False
    >>> dfa_accepts((trans, acc), "AAGCTT")   # check revcomp
    False
    """
    rc: str = reverse_complement(pattern)

    if pattern.upper() == rc.upper():
        # Palindromic — pattern and reverse complement are identical.
        return build_forbidden_pattern_dfa(pattern)

    # Non-palindromic — forbid both directions.
    return build_composite_dfa([pattern, rc])


# ---------------------------------------------------------------------------
# 4. DFA acceptance test
# ---------------------------------------------------------------------------

def dfa_accepts(
    dfa_tuple: tuple[list[list[int]], list[int]],
    sequence: str,
) -> bool:
    """Test whether a DFA accepts a given DNA sequence.

    Simulates the DFA on *sequence* character by character and checks
    whether the final state is accepting.

    Parameters
    ----------
    dfa_tuple:
        A ``(transition_table, accepting_states)`` tuple as returned by
        any ``build_*_dfa`` function in this module.
    sequence:
        DNA sequence string to test.

    Returns
    -------
    bool
        ``True`` if the DFA accepts *sequence*, ``False`` otherwise.

    Examples
    --------
    >>> trans, acc = build_forbidden_pattern_dfa("AAA")
    >>> dfa_accepts((trans, acc), "AAC")
    True
    >>> dfa_accepts((trans, acc), "AAA")
    False
    """
    delta, accepting = dfa_tuple
    accepting_set: set[int] = set(accepting)
    state: int = _INITIAL_STATE  # Start state is always _INITIAL_STATE.

    for char in sequence.upper():
        c: int = BASE_REV[char]
        state = delta[state][c]

    return state in accepting_set


# ---------------------------------------------------------------------------
# 5. DOT export
# ---------------------------------------------------------------------------

def dfa_to_dot(
    dfa_tuple: tuple[list[list[int]], list[int]],
    name: str = "dfa",
) -> str:
    """Export a DFA as a Graphviz DOT string.

    Useful for debugging and visualisation.  Accepting states are drawn as
    double circles; non-accepting states as single circles.  Transitions
    sharing the same ``(source, target)`` pair are merged with comma-
    separated labels.

    Parameters
    ----------
    dfa_tuple:
        A ``(transition_table, accepting_states)`` tuple.
    name:
        Name for the ``digraph`` (used in the DOT header).

    Returns
    -------
    str
        A string in Graphviz DOT format.

    Examples
    --------
    >>> trans, acc = build_forbidden_pattern_dfa("GC")
    >>> dot = dfa_to_dot((trans, acc), name="gc_dfa")
    >>> "digraph gc_dfa" in dot
    True
    """
    delta, accepting = dfa_tuple
    accepting_set: set[int] = set(accepting)
    num_states: int = len(delta)

    lines: list[str] = [
        f"digraph {name} {{",
        "    rankdir=LR;",
        "",
    ]

    # Accepting states → double circle.
    if accepting:
        acc_str = " ".join(str(s) for s in sorted(accepting))
        lines.append(f"    node [shape=doublecircle]; {acc_str};")

    # Non-accepting states → single circle.
    non_accepting: list[int] = sorted(s for s in range(num_states) if s not in accepting_set)
    if non_accepting:
        nacc_str = " ".join(str(s) for s in non_accepting)
        lines.append(f"    node [shape=circle]; {nacc_str};")

    lines.append("")

    # Invisible start node + arrow to state 0.
    lines.append('    "" [shape=point];')
    lines.append(f'    "" -> {_INITIAL_STATE} [label="start"];')
    lines.append("")

    # Merge parallel edges for readability.
    edges: dict[tuple[int, int], list[str]] = defaultdict(list)
    for state in range(num_states):
        for c in range(_ALPHABET_SIZE):
            target: int = delta[state][c]
            edges[(state, target)].append(BASE_MAP[c])

    for (src, dst), labels in sorted(edges.items()):
        merged: str = ",".join(labels)
        lines.append(f"    {src} -> {dst} [label=\"{merged}\"];")

    lines.append("}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 6. DFA negation
# ---------------------------------------------------------------------------

def negate_dfa(
    transition_table: list[list[int]],
    accepting_states: list[int],
) -> tuple[list[list[int]], list[int]]:
    """Negate a DFA by swapping accepting and non-accepting states.

    If the original DFA accepts language *L*, the negated DFA accepts the
    complement *Σ* \\ *L*.  This is useful when you have a DFA that
    **accepts** forbidden patterns and want to invert its verdict.

    .. note::

       This simple negation is correct **only** for *complete* DFAs
       (every state has a defined transition for every symbol).  All DFAs
       produced by this module are complete by construction.

    Parameters
    ----------
    transition_table:
        The transition table (returned **unchanged**).
    accepting_states:
        The current set of accepting states.

    Returns
    -------
    transition_table : list[list[int]]
        Same as input (identity — the topology does not change).
    new_accepting_states : list[int]
        All states that were previously non-accepting.

    Examples
    --------
    >>> trans, acc = build_forbidden_pattern_dfa("AA")
    >>> # Original: accepts strings WITHOUT "AA"
    >>> dfa_accepts((trans, acc), "AC")
    True
    >>> # Negated: accepts strings WITH "AA"
    >>> ntrans, nacc = negate_dfa(trans, acc)
    >>> dfa_accepts((ntrans, nacc), "AAC")
    True
    """
    num_states: int = len(transition_table)
    old_accepting: set[int] = set(accepting_states)
    new_accepting: list[int] = [s for s in range(num_states) if s not in old_accepting]
    return transition_table, new_accepting


# ---------------------------------------------------------------------------
# 7. T-run DFA (forbids max_t+1 consecutive T's)
# ---------------------------------------------------------------------------

def build_trun_dfa(
    max_t: int = _DEFAULT_MAX_T_RUN,
) -> tuple[list[list[int]], list[int]]:
    """Build a DFA that rejects strings with *max_t*+1 or more consecutive T's.

    States 0..max_t represent "have seen this many consecutive T's so far";
    state max_t+1 is the **trap** (forbidden) state.  Any non-T character
    resets the counter to 0.

    Parameters
    ----------
    max_t:
        Maximum allowed number of consecutive T's.  Default is 5, meaning
        6+ consecutive T's are forbidden.

    Returns
    -------
    transition_table : list[list[int]]
        ``transition_table[state][letter] = next_state``.
    accepting_states : list[int]
        States 0..max_t (all except the trap).

    Examples
    --------
    >>> trans, acc = build_trun_dfa(5)
    >>> dfa_accepts((trans, acc), "TTTTT")   # 5 T's — OK
    True
    >>> dfa_accepts((trans, acc), "TTTTTT")  # 6 T's — forbidden
    False
    >>> dfa_accepts((trans, acc), "TTTATTT")  # Broken run — OK
    True
    """
    trap: int = max_t + 1
    n_states: int = trap + 1  # states 0..trap

    delta: list[list[int]] = [[_INITIAL_STATE] * _ALPHABET_SIZE for _ in range(n_states)]

    for state in range(n_states):
        if state == trap:
            # Trap state: stay trapped forever regardless of input.
            for c in range(_ALPHABET_SIZE):
                delta[state][c] = trap
        else:
            for c in range(_ALPHABET_SIZE):
                if c == _T_BASE_IDX:
                    # T base: increment T count (or enter trap).
                    delta[state][c] = state + 1 if state < max_t else trap
                else:
                    # Non-T base: reset counter to 0.
                    delta[state][c] = _INITIAL_STATE

    # Accepting states: everything except the trap.
    accepting: list[int] = list(range(max_t + 1))
    return delta, accepting


# ---------------------------------------------------------------------------
# 8. OR-Tools format conversion
# ---------------------------------------------------------------------------

def dfa_to_ortools_format(
    transition_table: list[list[int]],
    accepting_states: list[int],
) -> tuple[int, list[list[int]], list[int]]:
    """Convert a DFA to OR-Tools ``AddAutomaton`` format.

    OR-Tools' ``CpModel.AddAutomaton`` expects:

    * ``initial_state``: int
    * ``transition_list``: list of ``[from_state, symbol, to_state]``
    * ``final_states``: list of accepting state IDs

    The transition table format used by this module
    (``delta[state][symbol] = next_state``) is converted to the flat
    triple-list format.

    Parameters
    ----------
    transition_table:
        ``delta[state][letter] = next_state`` as produced by any
        ``build_*_dfa`` function.
    accepting_states:
        List of accepting state IDs.

    Returns
    -------
    initial_state : int
        Always 0 (root state for all DFAs in this module).
    transition_list : list[list[int]]
        ``[[from_state, symbol, to_state], ...]``
    final_states : list[int]
        Sorted list of accepting state IDs.

    Examples
    --------
    >>> trans, acc = build_forbidden_pattern_dfa("GC")
    >>> init, tlist, final = dfa_to_ortools_format(trans, acc)
    >>> init
    0
    >>> len(tlist)  # 3 states × 4 symbols = 12 transitions
    12
    """
    transition_list: list[list[int]] = []
    for state in range(len(transition_table)):
        for symbol in range(_ALPHABET_SIZE):
            target: int = transition_table[state][symbol]
            transition_list.append([state, symbol, target])

    return _INITIAL_STATE, transition_list, sorted(accepting_states)
