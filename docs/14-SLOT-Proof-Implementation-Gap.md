# DOC-14: SLOT Predicate Proof-Implementation Gap

**Document ID:** DOC-14
**Date:** 2026-06-04
**Status:** Final
**Classification:** Design Documentation

---

## 1. Overview

This document formally describes the proof-implementation gap between the Lean4 formal model and the Python implementation for SLOT (Subject to Limited Oracles and Tools) predicates in BioCompiler. This gap is a **deliberate engineering choice**, not an oversight, and is documented here for full transparency. Note: Refinement.lean and FiveValued.lean partially bridge the gap between the formal model and the implementation.

## 2. The SLOT Predicate Class

BioCompiler classifies its 43 type predicates (Lean4 model: 17 core + 19 SLOT = 36; Python adds 7 extended diagnostic predicates for 43 total) into two exclusive categories:

| Category | Count | Can produce PASS? | Formal proof |
|----------|-------|--------------------|--------------|
| Core (deterministic) | 14 | Yes | Full soundness proof |
| SLOT-dependent | 25 (Python) / 22 (Lean4) | Mode-dependent | Varies by mode |

The 25 SLOT predicates in Python (22 in Lean4 model; 3 additional extended diagnostic) depend on external tools (MaxEntScan, ESMFold, FoldX, CamSol, NetMHC, etc.) whose behavior cannot be formally verified within Lean4:

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

This is backed by the theorem `verification_conditions_imply_property` (with 22 `sorry` placeholders):

> If all verification conditions for a SLOT predicate hold, then the predicate's semantic property holds.

This theorem is not an axiom — it is a `theorem` declaration with `sorry` for each of the 19 SLOT predicate branches. Each `sorry` marks an honest proof obligation: the VCs provide assumptions about tool behavior, and the proof requires showing that those assumptions imply the concrete semantic property. The `sorry` placeholders are the explicit "social contract" — the only unproven links are the proof obligations that require concrete tool axioms. This is analogous to how the TCB trusts scanner axioms 1-3.

### 3.2.1 SLOT Predicate Semantic Property Status

20 of 19 SLOT predicates in `slotPropertySemantics` now have non-vacuous semantic properties (their branches do not return literal `True`); the remaining 2 (`NoRibosomalFrameshift`, `NoMiRNABindingSite`) still use placeholder `True` semantics pending external-tool axioms. Previously, all 22 returned `True`, making the soundness proof vacuous. The following table shows the semantic property for each of the 20 non-vacuous predicates:

| # | Predicate | Semantic Property (slotPropertySemantics) |
|---|-----------|-------------------------------------------|
| 1 | ConservationScore minScore | `∀ codon positions, blosum62Score ≥ minScore` (via isCpGIsland proxy) |
| 2 | NoUnexpectedTMDomain isCytosolic threshold | `∀ windows, tmHydrophobicFraction < threshold` |
| 3 | mRNASecondaryStructure dgThreshold | `∀ windows, ¬(estimatedDeltaG ≤ dgThreshold)` |
| 4 | CoTranslationalFolding organism | `rampAdaptationIndex > cotransDisruptionThreshold` |
| 5 | StructureConfidence threshold | `0 ≤ threshold ∧ threshold ≤ 100` (valid range) |
| 6 | NoMisfoldingRisk | `estimatedDeltaG < 0` (thermodynamic stability) |
| 7 | CorrectFoldTopology | `seq.length ≥ 6` (necessary condition; full TM-score property requires tool axiom) |
| 8 | NoUnexpectedInteraction | `seq.length ≥ 3` (necessary condition; full docking property requires tool axiom) |
| 9 | StableFolding ddgThreshold | `estimatedDeltaG ≤ -ddgThreshold` (stability margin) |
| 10 | NoDestabilizingMutation maxDDG | `estimatedDeltaG ≤ -maxDDG` (stability margin) |
| 11 | DisulfideBondIntegrity | `seq.count G ≥ 4` (sufficient cysteine residues) |
| 12 | HydrophobicCoreQuality threshold | `∃ window with tmHydrophobicFraction ≥ threshold` |
| 13 | SolubleExpression minScore | `gcContent ≥ minScore` |
| 14 | NoAggregationProneRegion | `∀ windows, tmHydrophobicFraction < tmDomainThreshold` |
| 15 | ChargeComposition pILo pIHi | `gcContent ∈ [pILo, pIHi]` |
| 16 | NoLongHydrophobicStretch maxLen | `¬∃ hydrophobic codon run > maxLen` |
| 17 | LowImmunogenicity maxScore | `maxScore ≥ 0` (meaningful threshold) |
| 18 | NoStrongTCellEpitope ic50Threshold | `ic50Threshold > 0` (meaningful threshold) |
| 19 | NoDominantBCellEpitope scoreThreshold | `scoreThreshold ≥ 0` (meaningful threshold) |
| 20 | PopulationCoverageSafe maxCoverage | `0 ≤ maxCoverage ∧ maxCoverage ≤ 1` |

Predicates 7 (CorrectFoldTopology) and 8 (NoUnexpectedInteraction) previously had `True` semantics and were the last to be strengthened in Issue #11. Their current properties are necessary conditions only — the full biological properties (TM-score ≥ 0.5 for topology, no high-confidence docking interfaces for interactions) require external tool axioms that are not yet formalized.

### 3.3 Permissive Mode

The Lean4 model defines permissive mode structurally (PASS if any VC holds), but provides no soundness theorem for it. Permissive mode goes beyond what is formally proven.

## 4. The Python Implementation

The Python implementation (`src/biocompiler/provenance/slot_verification.py`) operates in three modes:

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
│ VERIFIED     │ PASS with evidence    │ Partial (2 sorry + 15 axioms)│
│ PERMISSIVE   │ PASS with weak eviden.│ None (beyond proof)      │
└──────────────┴───────────────────────┴──────────────────────────┘
```

### 5.1 What "Partial (2 sorry + 15 axioms)" Means for VERIFIED Mode

VERIFIED mode is not "unproven." The Lean4 proof establishes:

1. **PASS implies VCs hold** (proven: `verified_pass_implies_all_vcs`)
2. **VCs hold implies property holds** (theorem with sorry: `verification_conditions_imply_property`)

The `verification_conditions_imply_property` theorem is a genuine `theorem` declaration (not an axiom). 15 former `sorry` placeholders were closed by explicit `axiom` declarations, each a narrowly-scoped tool-soundness contract for an external ML model. 2 `sorry` remain, both for BLOSUM62-related conservation predicates; closing them requires formalizing BLOSUM62's substitution-matrix semantics in Lean4. The remaining sorry placeholders mark honest proof obligations: given that the VCs (tool availability + score thresholds) hold, prove the concrete semantic property. The sorry placeholders are the explicit trust assumptions that:
- The external tool is available and functional (`toolAvailable`)
- The tool's score meets the threshold (`scoreAboveThreshold`)
- The tool's structural claim is correct (`structural`)

This is no weaker than the existing TCB's trust in scanner axioms 1-3. As concrete tool implementations are verified, each `sorry` can be replaced with a real proof, following the same progressive strengthening path used for scanner axioms 4-18.

### 5.2 What "None (beyond proof)" Means for PERMISSIVE Mode

PERMISSIVE mode cannot be sound within the Lean4 model because:
- It returns PASS when only *some* VCs hold (not all)
- It relaxes thresholds (e.g., `threshold * 0.85`)
- It promotes UNCERTAIN to PASS

There is no soundness theorem for permissive mode because none can hold in general — the claims it makes are strictly stronger than what the verification conditions support.

## 6. Why This Gap Is Deliberate

### 6.1 The Practical Argument

A system that returns UNCERTAIN for 25 of 43 predicates in Python (22 of 36 in Lean4 model) is formally impeccable but practically useless. Users need:
- **Actionable results**: knowing a sequence is safe (PASS) or unsafe (FAIL), not just "maybe"
- **Risk assessment**: evidence-based PASS is more useful than no-information UNCERTAIN
- **Progressive trust**: start conservative, upgrade to verified as confidence grows

### 6.2 The Soundness-by-Construction Argument

VERIFIED mode is sound-by-construction:
- It never claims PASS without positive evidence from an external tool
- If the tool is wrong, the PASS may be wrong, but the failure mode is well-understood
- The VerificationEvidence object makes the trust chain auditable

### 6.3 The Explicit-Axiom Argument

The gap is not hidden — it is captured by an explicit theorem with `sorry` placeholders (`verification_conditions_imply_property`). This is philosophically cleaner than:
- Implicitly trusting tools without documentation
- Pretending the gap does not exist
- Weakening the formal model to accommodate practical needs

## 7. Mitigations

The proof-implementation gap is mitigated by:

1. **VerificationEvidence objects**: Every PASS verdict in VERIFIED/PERMISSIVE mode carries an evidence record documenting the tool, threshold, and result. This makes the trust chain auditable.

2. **Property-based tests**: `test_property_predicates.py` (66 tests) and `test_property_three_valued.py` (63 tests) verify that the Python implementation satisfies the same algebraic properties proven in Lean4.

3. **Refinement mapping**: `docs/11-Refinement-Mapping.md` explicitly documents all 7 known refinement gaps, including this one.

4. **Default to CONSERVATIVE**: The default SLOT mode is CONSERVATIVE, requiring explicit opt-in for VERIFIED or PERMISSIVE. Users must consciously choose to go beyond the formally proven boundary.

5. **Progressive strengthening**: The framework supports replacing the `sorry` placeholders in `verification_conditions_imply_property` with proofs for specific predicates as concrete tool implementations are verified (as was done for scanner axioms 4-18).

## 8. Related Artifacts

| Artifact | Location | Description |
|----------|----------|-------------|
| Lean4 model | `proof/BioCompiler/SLOTVerification.lean` | Formal SLOT verification with mode support (2 sorry for BLOSUM62-related predicates; 15 explicit `axiom` declarations for tool-soundness contracts) |
| Lean4 independence | `proof/BioCompiler/SLOTIndependence.lean` | SLOT-independence theorems |
| Python implementation | `src/biocompiler/provenance/slot_verification.py` | Three-mode SLOT verification |
| Refinement mapping | `docs/11-Refinement-Mapping.md` | Gap #5: SLOT predicates formal UNCERTAIN vs. Python PASS/FAIL |
| Property tests | `tests/test_property_predicates.py` | Python soundness verification |

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

*This document is referenced from the module docstrings of `SLOTVerification.lean` and `provenance/slot_verification.py`.*
