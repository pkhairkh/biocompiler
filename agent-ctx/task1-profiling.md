# Task 1: Eukaryotic Pipeline Profiling Results

## Baseline Performance
| Protein  | AvgTime | CAI    |
|----------|---------|--------|
| GFP      | 12.2ms  | 1.0000 |
| mCherry  | 11.7ms  | 1.0000 |
| Insulin  | 3.1ms   | 1.0000 |
| HBB      | 18.9ms  | 1.0000 |

## Phase Breakdown (HBB)
- Phase 1 (greedy init): 0.05ms — trivial
- Phase 2 (constraint sat): 16.91ms — **THE BOTTLENECK** (89% of time)
- Phase 3 (CAI hill climb): 0.71ms — minor

## Phase 2 Root Cause Analysis
- 50 iterations, 84 "fixes" attempted
- **Same 3 avoidable GT violations persist across ALL iterations** (positions 44, 365, 434)
- **Same 5 CpG violations persist across ALL iterations**
- These are cross-codon GTs that can't be fixed: every swap attempt fails and rolls back
- Total detect time: 2.62ms, Total fix time: 12.57ms
- Fix attempts dominate: each swap-check-rollback cycle is expensive

## Key Issues
1. **No stale violation detection**: Loop doesn't detect that violation set is unchanged
2. **No early termination**: 50 iterations even when violations are unfixable
3. **_compute_cai recomputes max_adapt every call**: Already precomputed as self._max_adapt
4. **CpG fix too conservative**: Requires no GT increase, blocking most CpG fixes
5. **has_any_restriction_site() scans full sequence** after every attempted swap
6. **MaxEntScan scores recomputed**: No caching of 9-mer/23-mer scores
