# DOC-14: SLOT Predicate Proof-Implementation Gap

**Document ID:** DOC-14
**Date:** 2026-06-04
**Status:** Final
**Classification:** Design Documentation

---

## 1. Overview

This document formally describes the proof-implementation gap between the Lean4 formal model and the Python implementation for SLOT (Subject to Limited Oracles and Tools) predicates in BioCompiler. This gap is a **deliberate engineering choice**, not an oversight, and is documented here for full transparency.

## 2. The SLOT Predicate Class

BioCompiler classifies its 33 type predicates into two exclusive categories:

| Category | Count | Can produce PASS? | Formal proof |
|----------|-------|--------------------|--------------|
| Core (deterministic) | 13 | Yes | Full soundness proof |
| SLOT-dependent | 20 | Mode-dependent | Varies by mode |

The 20 SLOT predicates depend on external tools (MaxEntScan, ESMFold, FoldX, CamSol, NetMHC, etc.) whose behavior cannot be formally verified within Lean4:

1. ConservationScore (BLOSUM62 heuristic)
2. NoUnexpectedTMDomain (hydrophobic fraction heuristic)
3. mRNASecondaryStructure (simplified folding model)
4. CoTranslationalFolding (CAI-based heuristic)
5. StructureConfidence (requires ESMFold)
6. NoMisfoldingRisk (requires FoldX)
7. CorrectFoldTopology (requires ESMFold)
8. NoUnexpectedInteraction (requires ESMFold)
9. StableFolding (requires FoldX)
10. NoDestabilizingMutation (requires FoldX)
11. DisulfideBondIntegrity (requires ESMFold)
12. HydrophobicCoreQuality (requires FoldX)
13. SolubleExpression (requires CamSol)
14. NoAggregationProneRegion (requires CamSol)
15. ChargeComposition (requires ExPASy)
16. NoLongHydrophobicStretch (hydrophobicity scanner)
17. LowImmunogenicity (requires NetMHC)
18. NoStrongTCellEpitope (requires NetMHCpan)
19. NoDominantBCellEpitope (requires BepiPred)
20. PopulationCoverageSafe (requires IEDB)

## 3. The Formal Model (Lean4)

### 3.1 Conservative Mode

In the Lean4 model (`proof/BioCompiler/SLOTVerification.lean`), the default behavior for SLOT predicates is **conservative mode**: always return `UNCERTAIN`.

```lean
def evaluateSLOT (mode : SLOTMode) (vctx : VerificationContext)
    (P : TypePredicate) (seq : Sequence) (ctx : CellularContext) : Verdict :=
  match mode with
  | SLOTMode.conservative => UNCERTAIN
  | SLOTMode.verified =>
    if allVCHold vctx (slotVCs P) then PASS else UNCERTAIN
  | SLOTMode.permissive =>
    if anyVCHold vctx (slotVCs P) then PASS else UNCERTAIN
```

The soundness theorem `slot_soundness_conservative` is **vacuously true**:

> If `evaluateWithMode SLOTMode.conservative vctx P seq ctx = PASS`, then `propertyHolds P seq ctx`.

Since conservative mode never returns PASS for SLOT predicates (proven in `conservative_is_safe`), the implication holds trivially. This is the strongest possible formal guarantee — it cannot be wrong because it never makes a positive claim.

### 3.2 Verified Mode

The Lean4 model also defines a verified mode, where SLOT predicates can return PASS when all verification conditions hold. The theorem `slot_soundness_verified` establishes:

> If `evaluateSLOT SLOTMode.verified vctx P seq ctx = PASS`, then `slotPropertySemantics P seq ctx`.

This is backed by the axiom `verification_conditions_imply_property`:

> If all verification conditions for a SLOT predicate hold, then the predicate's semantic property holds.

This axiom is the explicit "social contract" — the only unproven link is the trust that external tools are correct. This is analogous to how the TCB trusts scanner axioms 1-3.

### 3.3 Permissive Mode

The Lean4 model defines permissive mode structurally (PASS if any VC holds), but provides no soundness theorem for it. Permissive mode goes beyond what is formally proven.

## 4. The Python Implementation

The Python implementation (`src/biocompiler/slot_verification.py`) operates in three modes:

### 4.1 CONSERVATIVE Mode

```python
if slot_mode == SLOTMode.CONSERVATIVE:
    return Verdict.UNCERTAIN, VerificationEvidence(...)
```

**Matches the Lean4 model EXACTLY.** The Lean4 soundness proof covers this mode perfectly (vacuously sound).

### 4.2 VERIFIED Mode

```python
if slot_mode == SLOTMode.VERIFIED:
    if tool_available:
        result = run_check(seq, threshold)
        return result.verdict, VerificationEvidence(
            verified=(result.verdict == Verdict.PASS),
            tool_result=result.details,
            threshold_used=threshold,
        )
    else:
        return Verdict.UNCERTAIN, VerificationEvidence(
            tool_available=False, verified=False,
        )
```

**Sound-by-construction:** If the external tool is wrong, the predicate may be wrong, but it will **never** claim PASS without evidence. The tool's output IS the evidence. This corresponds to the Lean4 verified mode under the axiom `verification_conditions_imply_property`.

### 4.3 PERMISSIVE Mode

```python
# Relaxed thresholds
effective_threshold = threshold * 0.85  # or similar relaxation
# UNCERTAIN promoted to PASS
if result.verdict == Verdict.UNCERTAIN:
    return Verdict.PASS, VerificationEvidence(...)
```

**Goes beyond what is formally proven.** Provides practical utility at the cost of formal guarantees. The Python implementation relaxes thresholds and promotes UNCERTAIN to PASS in borderline cases.

## 5. Formal Coverage Summary

```
┌──────────────┬───────────────────────┬──────────────────────────┐
│ Mode         │ Python behavior       │ Lean4 proof coverage     │
├──────────────┼───────────────────────┼──────────────────────────┤
│ CONSERVATIVE │ Always UNCERTAIN      │ Full (vacuously sound)   │
│ VERIFIED     │ PASS with evidence    │ Partial (under axiom)    │
│ PERMISSIVE   │ PASS with weak eviden.│ None (beyond proof)      │
└──────────────┴───────────────────────┴──────────────────────────┘
```

### 5.1 What "Partial (under axiom)" Means for VERIFIED Mode

VERIFIED mode is not "unproven." The Lean4 proof establishes:

1. **PASS implies VCs hold** (proven: `verified_pass_implies_all_vcs`)
2. **VCs hold implies property holds** (axiom: `verification_conditions_imply_property`)

The only unproven link is the axiom that says "if the verification conditions hold, the property actually holds." This axiom is the explicit trust assumption that:
- The external tool is available and functional (`toolAvailable`)
- The tool's score meets the threshold (`scoreAboveThreshold`)
- The tool's structural claim is correct (`structural`)

This is no weaker than the existing TCB's trust in scanner axioms 1-3.

### 5.2 What "None (beyond proof)" Means for PERMISSIVE Mode

PERMISSIVE mode cannot be sound within the Lean4 model because:
- It returns PASS when only *some* VCs hold (not all)
- It relaxes thresholds (e.g., `threshold * 0.85`)
- It promotes UNCERTAIN to PASS

There is no soundness theorem for permissive mode because none can hold in general — the claims it makes are strictly stronger than what the verification conditions support.

## 6. Why This Gap Is Deliberate

### 6.1 The Practical Argument

A system that returns UNCERTAIN for 20 of 33 predicates is formally impeccable but practically useless. Users need:
- **Actionable results**: knowing a sequence is safe (PASS) or unsafe (FAIL), not just "maybe"
- **Risk assessment**: evidence-based PASS is more useful than no-information UNCERTAIN
- **Progressive trust**: start conservative, upgrade to verified as confidence grows

### 6.2 The Soundness-by-Construction Argument

VERIFIED mode is sound-by-construction:
- It never claims PASS without positive evidence from an external tool
- If the tool is wrong, the PASS may be wrong, but the failure mode is well-understood
- The VerificationEvidence object makes the trust chain auditable

### 6.3 The Explicit-Axiom Argument

The gap is not hidden — it is captured by an explicit axiom (`verification_conditions_imply_property`). This is philosophically cleaner than:
- Implicitly trusting tools without documentation
- Pretending the gap doesn't exist
- Weakening the formal model to accommodate practical needs

## 7. Mitigations

The proof-implementation gap is mitigated by:

1. **VerificationEvidence objects**: Every PASS verdict in VERIFIED/PERMISSIVE mode carries an evidence record documenting the tool, threshold, and result. This makes the trust chain auditable.

2. **Property-based tests**: `test_property_predicates.py` (66 tests) and `test_property_three_valued.py` (63 tests) verify that the Python implementation satisfies the same algebraic properties proven in Lean4.

3. **Refinement mapping**: `docs/11-Refinement-Mapping.md` explicitly documents all 7 known refinement gaps, including this one.

4. **Default to CONSERVATIVE**: The default SLOT mode is CONSERVATIVE, requiring explicit opt-in for VERIFIED or PERMISSIVE. Users must consciously choose to go beyond the formally proven boundary.

5. **Progressive strengthening**: The framework supports replacing the axiom `verification_conditions_imply_property` with proofs for specific predicates as concrete tool implementations are verified (as was done for scanner axioms 4-18).

## 8. Related Artifacts

| Artifact | Location | Description |
|----------|----------|-------------|
| Lean4 model | `proof/BioCompiler/SLOTVerification.lean` | Formal SLOT verification with mode support |
| Lean4 independence | `proof/BioCompiler/SLOTIndependence.lean` | SLOT-independence theorems |
| Python implementation | `src/biocompiler/slot_verification.py` | Three-mode SLOT verification |
| Refinement mapping | `docs/11-Refinement-Mapping.md` | Gap #5: SLOT predicates formal UNCERTAIN vs. Python PASS/FAIL |
| Property tests | `tests/test_property_predicates.py` | Python soundness verification |
| Worklog entry | `worklog.md` | Task 1b, Task 1f entries |

## 9. Decision Record

- **Decision**: Allow SLOT predicates to return PASS in VERIFIED and PERMISSIVE modes, despite the Lean4 model only proving soundness for CONSERVATIVE mode.
- **Rationale**: Practical utility outweighs formal completeness. The gap is made explicit via axioms and documented evidence.
- **Alternatives considered**:
  - Always UNCERTAIN: formally sound but practically useless
  - Drop Lean4 model: lose all formal guarantees
  - Prove tool correctness: infeasible for external tools (ESMFold, FoldX, etc.)
- **Risk**: PASS verdicts in VERIFIED/PERMISSIVE mode may be wrong if external tools are wrong
- **Mitigation**: VerificationEvidence audit trail; default CONSERVATIVE mode; progressive strengthening path

---

*This document is referenced from the module docstrings of `SLOTVerification.lean` and `slot_verification.py`.*
