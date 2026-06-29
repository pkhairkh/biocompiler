# ADR-0015: Biosecurity Sequence Screening

## Status: Accepted

## Date

2026-03-05

## Context

BioCompiler designs DNA sequences that could be physically synthesized by gene synthesis services. Without screening, it could be used to optimize hazardous sequences — toxins, viral proteins, antibiotic resistance markers, and other sequences of concern identified by international biosecurity frameworks.

The dual-use nature of gene design tools is well-documented. A codon optimizer that produces high-CAI, constraint-satisfying sequences for any input protein could inadvertently (or deliberately) be used to enhance the expressibility of dangerous biological agents. Specifically:

1. **Toxin genes**: Sequences encoding potent toxins (e.g., ricin, botulinum, shiga) could be optimized for high expression in target organisms, increasing yield and potency.

2. **Viral proteins**: Capsid, polymerase, and other structural/functional viral proteins could be optimized for expression in heterologous systems, facilitating reverse genetics or subunit production without appropriate biosafety oversight.

3. **Antibiotic resistance markers**: Optimizing resistance genes (e.g., beta-lactamases, vancomycin resistance) for high expression could contribute to the spread of antimicrobial resistance, a WHO-classified global health threat.

4. **Regulated pathogens**: Sequences from Select Agents (US), Category A-C biological agents (EU), and other nationally regulated pathogens could be designed without any checkpoint or review.

International frameworks mandate screening at the gene synthesis stage (IHRA Screening Framework, WHO Laboratory Biorisk Management), but these apply only at synthesis. A design-stage tool has a responsibility to implement upstream screening, as the design output is often the direct input to a synthesis order.

**Alternatives Considered:**

1. **No screening** — Trust the user and rely on synthesis providers to screen. But: not all providers screen consistently; some users synthesize in-house; and the design tool has independent ethical obligations.

2. **Post-hoc screening only** — Screen the output but not the input. But: by the time the output is produced, computational resources have already been spent optimizing a hazardous sequence, and the result exists in the user's environment.

3. **Mandatory pre-optimization screening** — Screen the input protein sequence against known hazardous sequence signatures before any optimization begins. Block critical/high-risk sequences entirely. Warn on medium-risk sequences but allow the user to proceed with explicit acknowledgment.

## Decision

Implement mandatory pre-optimization screening against known hazardous sequence signatures (Alternative 3).

1. **Screening database**: Maintain a curated database of hazardous sequence signatures derived from:
   - Select Agent list protein sequences (US CDC/APHIS)
   - WHO high-consequence pathogen protein sequences
   - Antibiotic resistance gene databases (CARD, ResFinder)
   - Toxin sequences from UniProt annotations
   - Sequences identified in the IHRA Screening Framework

2. **Risk classification**:
   - **Critical**: Exact or near-exact matches (>95% identity over >80% length) to known hazardous sequences. Optimization is **blocked** — the function raises a `BiosecurityError` with details about the matched signature and its risk category.
   - **High-risk**: Strong matches to regulated pathogens or toxin sequences. Optimization is **blocked** — same behavior as Critical.
   - **Medium-risk**: Partial matches or signatures from regulated but lower-concern organisms. Optimization proceeds but emits a `UserWarning` that is included in the result metadata and all export formats.
   - **Low-risk**: Weak matches or distant homology. Optimization proceeds with an informational note in the result metadata.
   - **No match**: No action required. Normal optimization proceeds.

3. **Screening timing**: Screening runs **before** any optimization step. The input protein sequence is checked against the signature database before codon selection, constraint solving, or any transformation. This ensures no computational resources are spent on blocked sequences and no intermediate results exist for blocked inputs.

4. **Biosafety annotations**: All optimization results include a `biosecurity_screen_result` field with:
   - `risk_level`: "none", "low", "medium", "high", or "critical"
   - `matched_signatures`: list of matched signature IDs and descriptions
   - `screening_timestamp`: when the screen was performed
   - `database_version`: version of the signature database used

5. **Database maintenance**: The signature database is versioned and shipped with each release. Updates follow the project's security advisory process. Users can configure a custom database path for institution-specific signatures.

## Consequences

- **Positive**: Prevents misuse of BioCompiler for enhancing hazardous sequences. Aligns with biosecurity best practices recommended by the International Gene Synthesis Consortium (IGSC), the IHRA Screening Framework, and WHO Laboratory Biorisk Management guidance.
- **Positive**: Design-stage screening provides an earlier checkpoint than synthesis-stage screening alone, reducing the risk that hazardous designs reach the synthesis stage.
- **Positive**: Biosafety annotations in all exports ensure that downstream consumers (synthesis providers, lab LIMS systems) can make informed decisions about the sequence.
- **Negative**: May produce false positives for legitimate research — e.g., a researcher studying toxin mechanisms may have a legitimate need to optimize a toxin sequence for structural biology. The medium-risk pathway (warning + acknowledgment) provides an escape valve, but the high-risk block cannot be overridden.
- **Negative**: The signature database requires ongoing maintenance to remain current with evolving biosecurity threat assessments. Stale databases may miss emerging threats; overly aggressive databases may produce excessive false positives.
- **Trade-off**: Speed impact is minimal — pattern matching against a curated database is fast (O(n×m) where n is input length and m is signature count, with BLAST-style indexing for larger databases). Completeness of the signature database vs. false positive rate is the primary trade-off; the 95%/80% threshold was chosen to minimize false positives while catching functional variants.

## References

- International Health Regulations (IHR) Article 44: Biosafety and Biosecurity
- WHO Laboratory Biorisk Management Framework (2024)
- IGSC Harmonized Screening Protocol v2.0
- US Select Agent Regulations (7 CFR Part 331, 9 CFR Part 121, 42 CFR Part 73)
- ADR-0016: Default Safety Measures (biosafety annotations in exports)
