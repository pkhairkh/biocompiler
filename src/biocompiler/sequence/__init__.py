"""
BioCompiler Sequence Domain
============================

Sequence-level analysis and manipulation — everything that operates
on DNA strings directly.  May import from ``biocompiler.shared`` but
NOT from optimizer/, solver/, or other higher-level subpackages.

Contents:
  - scanner.py           : Multi-DFA motif detection (scan_sequence, gc_content, …)
  - aho_corasick.py      : O(L+M) multi-pattern matching
  - restriction_sites.py : Curated restriction enzyme database
  - iupac.py             : IUPAC ambiguity code resolution
  - pattern_enforcement  : EnforcePattern / AvoidPattern constraints
  - sliding_gc.py        : Sliding-window GC constraint
  - local_gc.py          : Region-specific GC constraints
  - splicing.py          : NDFST splicing engine (isoform computation)
  - maxentscan.py        : MaxEntScan splice site scoring (Yeo & Burge 2004 trained parameters)
  - import_seq.py        : FASTA / GenBank sequence import
"""

# Re-export key symbols for convenient access:
#   from biocompiler.sequence import scan_sequence, score_donor, …

from biocompiler.sequence.scanner import (  # noqa: F401
    validate_dna_sequence,
    gc_content,
    scan_sequence,
)

from biocompiler.sequence.aho_corasick import (  # noqa: F401
    AhoCorasickScanner,
    build_scanner_from_enzymes,
    build_scanner_from_sites,
)

from biocompiler.sequence.restriction_sites import (  # noqa: F401
    RESTRICTION_SITES,
    MIN_SITE_LENGTH,
    get_recognition_site,
    expand_iupac_site,
    get_eliminable_sites,
)

from biocompiler.sequence.iupac import (  # noqa: F401
    IUPAC_DNA,
    resolve_ambiguous,
    is_ambiguous,
    expand_ambiguous,
    has_ambiguous,
    validate_iupac_sequence,
)

from biocompiler.sequence.pattern_enforcement import (  # noqa: F401
    PatternConstraint,
    PatternResult,
    check_pattern,
    check_patterns,
    enforce_pattern,
    enforce_patterns,
    build_avoidance_scanner,
)

from biocompiler.sequence.sliding_gc import (  # noqa: F401
    WindowViolation,
    SlidingGCResult,
    check_sliding_gc,
    fix_sliding_gc_violations,
    evaluate_sliding_gc,
)

from biocompiler.sequence.local_gc import (  # noqa: F401
    LocalGCConstraint,
    LocalGCResult,
    check_local_gc,
    optimize_local_gc,
)

from biocompiler.sequence.splicing import (  # noqa: F401
    maxent_score,
    maxent_score_v2,
    score_splice_sites,
    compute_splice_isoforms,
)

from biocompiler.sequence.maxentscan import (  # noqa: F401
    score_donor,
    score_acceptor,
    scan_splice_sites,
    max_donor_score,
    max_acceptor_score,
    validate_against_published,
    CRYPTIC_SPLICE_THRESHOLD,
)

from biocompiler.sequence.import_seq import (  # noqa: F401
    import_fasta,
    import_genbank,
    import_sequence,
)
