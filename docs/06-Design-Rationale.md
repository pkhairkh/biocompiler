# DOC-06: Design Rationale

| Field | Value |
|---|---|
| **Document ID** | DOC-06 |
| **Version** | 1.0.0-draft |
| **Status** | ROUGH DRAFT |
| **Date** | 2026-05-30 |
| **Prepared By** | BioCompiler Project Team |
| **Reviewed By** | [TBD at baseline review] |
| **Approved By** | [TBD at baseline review] |
| **Standard** | ISO/IEC/IEEE 42010 §7 — Architecture Rationale |

---

## 1. Introduction

### 1.1 Purpose

This document records the reasoning behind every major design decision in the BioCompiler system — not merely what was decided, but *why* it was decided, what alternatives were considered, and what was explicitly rejected. Design rationale serves three distinct audiences. First, it allows future maintainers to evaluate whether the original reasoning still holds when requirements change or new alternatives become available. Second, it provides stakeholders and reviewers with a traceable argument from stakeholder needs through architectural decisions to concrete design choices, enabling informed review rather than blind acceptance. Third, it prevents "decision drift," the common phenomenon where the original justification for a design choice is forgotten and the choice is either retained for the wrong reasons or modified without understanding why it was made.

A design rationale document is distinct from a design document. DOC-03 (Software Design Document) specifies *what* the design is — the algorithms, data structures, invariants, and interfaces. DOC-02 (Software Architecture Document) specifies *how* the system is decomposed into components and how they interact. This document specifies *why* those choices were made, and — equally important — *why alternatives were not chosen*. Every rejected alternative represents an explicit engineering judgment, and recording that judgment prevents later reintroduction of a known-problematic approach under the guise of "improvement."

This document integrates two foundational analyses that shaped the BioCompiler design: the critical analysis of the original "Compiler for Protein Synthesis" concept (DOC-09), which identified nine fatal flaws in the original proposal, and the deterministic methods framework (DOC-10), which established six formal methods for producing deterministic answers from non-deterministic biological systems. Both analyses are summarized here with cross-references to their full treatments.

### 1.2 Scope

This document covers all architectural decisions (AD-01 through AD-07) documented in DOC-02, the critical analysis of the original framework, the deterministic methods approach that replaced probabilistic reasoning, and the novel capabilities that emerge from the interaction of these decisions. It does not repeat the detailed interface contracts (DOC-04), the detailed algorithms (DOC-03), or the verification plan (DOC-05); instead, it provides the *reasoning* that led to those specifications.

### 1.3 References

| Ref ID | Document / Source | Description |
|---|---|---|
| **REF-01** | DOC-01: Software Requirements Specification | Defines the SHALL and SHALL NOT requirements that design decisions must satisfy. |
| **REF-02** | DOC-02: Software Architecture Document | Defines the component decomposition and architectural decisions (AD-01 through AD-07) that this document justifies. |
| **REF-03** | DOC-03: Software Design Document | Detailed algorithms, data structures, and invariants implementing the design decisions. |
| **REF-04** | DOC-09: Critical Analysis of Original Framework | Full treatment of the nine fatal flaws in the original "Compiler for Protein Synthesis" concept. |
| **REF-05** | DOC-10: Deterministic Methods for Non-Deterministic Biology | Full treatment of the six deterministic methods that replace probabilistic approaches. |
| **REF-06** | Cousot & Cousot, "Abstract Interpretation," POPL 1977 | Foundational reference for abstract interpretation as a sound over-approximation technique. |
| **REF-07** | Milner, "A Theory of Type Polymorphism in Programming," JCSS 1978 | Foundational reference for type-based verification: "Well-typed programs don't go wrong." |
| **REF-08** | Angluin, "Learning Regular Sets from Queries and Counterexamples," 1987 | The L* algorithm; referenced as a rejected alternative for grammar induction. |
| **REF-09** | ISO/IEC/IEEE 42010:2022 | Architecture description standard governing the structure of this rationale. |
| **REF-10** | Lattner & Adve, "LLVM: A Compilation Framework for Lifelong Program Analysis & Transformation," CGO 2004 | The LLVM architecture that inspired the BioCompiler pipeline design. |

---

## 2. Background: Why the Original Idea Doesn't Work

### 2.1 The Original Proposal

The original concept, titled "Compiler Framework for Protein Synthesis using Intermediate Representations," proposed treating the entire gene-to-protein pathway as a compilation process. The idea was seductively elegant: just as a software compiler transforms source code through intermediate representations (IRs) into machine code, a "protein compiler" would transform mRNA sequences through biological IRs into functional proteins. The proposal included four key claims: (1) that the gene-to-protein pathway could be modeled as a multi-stage compilation pipeline with typed IRs at each level; (2) that the mapping from mRNA to protein structure could be *learned* via grammar induction (specifically the ADIOS algorithm) from biological sequence data, producing a formal grammar that captures the "syntax" of protein synthesis; (3) that this grammar would support both forward compilation (mRNA → protein) and inverse compilation (protein → mRNA); and (4) that the entire framework would provide deterministic, verifiable guarantees about the compiled output.

The intellectual appeal is obvious. Compilation is one of the most successful engineering paradigms in computer science: it transforms a high-level specification through well-defined stages into a low-level artifact, with formal guarantees at each stage. If biology implements something like compilation, then the full weight of compiler engineering — type systems, formal verification, intermediate representations, optimization passes — could be brought to bear on biological engineering. The LLVM project demonstrated that a well-designed compiler infrastructure can support dozens of front-end languages, optimization passes, and back-end targets through a common IR. The original proposal aspired to an analogous "LLVM for biology."

The problem is that the analogy breaks down at nearly every point where it matters. The following nine flaws are not minor engineering difficulties that can be patched; they are fundamental mismatches between the assumptions of compilation and the realities of molecular biology. Several are fatal, meaning the original proposal cannot be made to work as stated regardless of engineering effort.

### 2.2 Nine Critical Flaws

#### Flaw 1: Protein Folding Is Not String Rewriting

| Attribute | Value |
|---|---|
| **Name** | Protein Folding Is Not String Rewriting |
| **Category** | Representational |
| **Severity** | Fatal |
| **Can Be Fixed?** | No — not within a string-rewriting framework |

Protein folding is a thermodynamic process governed by free-energy minimization in a high-dimensional conformational space. The native structure of a protein is the global free-energy minimum (or a deep local minimum) under physiological conditions, and it depends on the entire amino acid sequence, the solvent environment, temperature, pH, ionic strength, and the presence of chaperones and cofactors. The original proposal treated folding as a string-rewriting production — a grammar rule that transforms an amino acid sequence into a structural description. This is a category error. String rewriting is a purely symbolic operation: it transforms one string into another based on local pattern-matching rules. Folding is a global optimization process in which every residue potentially interacts with every other residue and with the solvent.

Graph grammars — the most sophisticated string-rewriting formalism the original proposal could invoke — can *describe* a folded structure after the fact, but they cannot *predict* it. Describing a known structure as the result of applying a sequence of graph productions is an exercise in post-hoc rationalization, not prediction. To predict folding, the grammar would need to encode the laws of physics (electrostatics, van der Waals forces, hydrogen bonding, hydrophobic effect, entropy of the unfolded state) as production rules, which is not what grammars do. The Anfinsen dogma (the sequence determines the structure) is true in a thermodynamic sense, but "determines" here means "constrains the free-energy landscape," not "specifies through symbolic rules."

This flaw is fatal because it undermines the entire premise of grammar-based folding prediction. No amount of grammar engineering — more productions, more nonterminals, context-sensitive rules — can transform a symbolic system into a physical one. The fix adopted in BioCompiler is to *not attempt folding prediction internally* and instead delegate it to external physics-based and machine-learning tools through the FFI boundary (AD-06).

#### Flaw 2: Post-Translational Modifications Defy Formal Rules

| Attribute | Value |
|---|---|
| **Name** | PTMs Defy Formal Rules |
| **Category** | Representational |
| **Severity** | Fatal |
| **Can Be Fixed?** | No — not within a sequence-grammar framework |

Post-translational modifications (PTMs) — phosphorylation, glycosylation, acetylation, ubiquitination, methylation, and hundreds of others — are chemical modifications applied to proteins after translation. The original proposal needed PTMs to be derivable from the protein sequence alone (or from the sequence plus a formal grammar), because the IR pipeline could only carry sequence-derived information forward. In reality, PTMs are co-determined by the protein sequence and the cellular environment in ways that have no known formal characterization.

Phosphorylation, for example, occurs at serine, threonine, and tyrosine residues, but not at every such residue — only at those recognized by specific kinases that are present in specific cell types under specific signaling conditions. The same protein will be phosphorylated at different sites in a neuron than in a hepatocyte, and at different sites in the same cell type depending on signaling state. Glycosylation patterns depend on the enzyme repertoire of the endoplasmic reticulum and Golgi apparatus, which varies by tissue, developmental stage, and disease state. The mapping from (sequence, cell type) to PTM pattern is not a function — the same sequence in the same cell type can produce different PTM patterns depending on the cell's signaling history.

This flaw is fatal for the same reason as Flaw 1: the information required to determine PTMs is not contained in the sequence and cannot be captured by a grammar over sequences. The original proposal's IR pipeline could not accommodate PTMs without incorporating a model of cellular state that goes far beyond what any formal grammar can represent. The BioCompiler fix is the same as for folding: delegate PTM prediction to external tools through the FFI boundary, and treat PTM predictions as non-deterministic SLOT fields in the IR (AD-06).

#### Flaw 3: ADIOS Learns Statistics, Not Physics

| Attribute | Value |
|---|---|
| **Name** | ADIOS Learns Statistics, Not Physics |
| **Category** | Methodological |
| **Severity** | Critical |
| **Can Be Fixed?** | Partially — by replacing grammar induction with hand-curated grammars from biological knowledge |

The original proposal relied on the ADIOS (Automatic Distillation of Structure) algorithm to learn a grammar for the mRNA-to-protein mapping from biological sequence data. ADIOS is a grammar induction algorithm that extracts statistical regularities from corpora of symbol sequences. It identifies significant patterns (equivalence classes and sequential motifs) by computing statistical significance of n-gram occurrences against a null model. The fundamental problem is that statistical regularities in biological sequences are not the same as physical laws governing protein behavior.

Consider: the genetic code itself is a statistical regularity in biological sequences (certain codons co-occur with certain amino acids), but it is also a physical mapping implemented by the ribosome and tRNA synthetases. The statistical regularity reflects the physical mechanism, but ADIOS cannot distinguish between a statistical regularity that reflects a physical law and one that reflects historical accident or sampling bias. A grammar learned by ADIOS would capture whatever patterns are statistically significant in the training corpus, including artifacts of corpus composition (over-representation of certain gene families, species bias, sequencing errors) and biologically irrelevant correlations (GC content covarying with taxonomic origin, not with functional properties).

This flaw is critical rather than fatal because the *idea* of grammar induction is not inherently wrong — it is the choice of ADIOS and the assumption that statistical learning produces physically meaningful grammars that is the problem. The BioCompiler fix is to abandon grammar induction entirely (REQ-CON-002) and instead construct grammars from curated biological knowledge: GENCODE annotations for splice sites, experimentally validated consensus sequences for regulatory elements, and organism-specific codon usage tables. These are not learned from raw data; they are compiled from the experimental literature by domain experts. This trades the ambition of automatic discovery for the reliability of curated knowledge, which is the correct trade-off for an engineering tool that must provide formal guarantees.

#### Flaw 4: Grammar Induction for CSGs Is Undecidable

| Attribute | Value |
|---|---|
| **Name** | Grammar Induction for Context-Sensitive Languages Is Undecidable |
| **Category** | Methodological |
| **Severity** | Critical |
| **Can Be Fixed?** | Partially — by restricting to regular or context-free approximations, but these are insufficient for biology |

The original proposal needed the learned grammar to be at least context-sensitive, because biological processes exhibit context-dependent behavior (splicing depends on the presence of enhancer and silencer elements at arbitrary distances from the splice site; folding depends on long-range residue interactions). Grammar induction — the problem of learning a grammar from examples of the language it generates — is well-understood for regular languages (Angluin's L* algorithm, 1987, can learn a DFA in polynomial time given a minimally adequate teacher). For context-free languages, the problem becomes significantly harder, though some positive results exist for subclasses. For context-sensitive languages, grammar induction is undecidable in general: there is no algorithm that, given a finite set of positive and negative examples, can determine the context-sensitive grammar that generated them.

This is not an engineering limitation that faster hardware or more data can overcome; it is a fundamental result in formal language theory. The original proposal's reliance on grammar induction for context-sensitive biological grammars is therefore untenable. Even if ADIOS could learn *some* grammar from the data (which it can, as a statistical procedure), there is no guarantee that the learned grammar is the "correct" one — or even that it is in the right class. A grammar that under-approximates the true language will miss valid biological sequences; one that over-approximates will include invalid ones.

The BioCompiler fix is consistent with Flaw 3: abandon grammar induction and use hand-curated grammars. The splicing grammar (AD-03) is explicitly restricted to the class of regular languages (finite-state transducers), which is sufficient for modeling the local consensus sequences and regulatory element thresholds that govern splicing decisions, even though it cannot capture arbitrary long-range dependencies. This is a deliberate approximation that accepts incompleteness in exchange for decidability and formal guarantees.

#### Flaw 5: Determinism Assumption Is Fatal

| Attribute | Value |
|---|---|
| **Name** | Determinism Assumption Is Fatal |
| **Category** | Methodological |
| **Severity** | Fatal |
| **Can Be Fixed?** | Yes — by modeling non-determinism explicitly (as BioCompiler does with NDFSTs and three-valued logic) |

The original proposal assumed that the mapping from mRNA to protein was deterministic — that a given mRNA sequence produces a unique protein product. This assumption is catastrophically wrong. In eukaryotic biology, 95% or more of human multi-exon genes undergo alternative splicing, meaning a single pre-mRNA can produce multiple distinct mRNA isoforms, each encoding a different protein. The "determinism" of the original proposal was achieved by averaging over isoforms or selecting the most common one, which destroys precisely the information that makes the mapping interesting and biologically relevant.

Consider the *DSCAM* gene in *Drosophila*, which can produce 38,016 distinct protein isoforms through alternative splicing. A "deterministic" mapping that averages over these isoforms produces a meaningless chimera; one that selects the most common isoform discards 38,015 biologically real proteins. The assumption of determinism is not merely inaccurate — it fundamentally misrepresents the nature of the biological process. Alternative splicing is not noise to be averaged away; it is a regulated, biologically essential mechanism for proteomic diversity.

This flaw is fatal for the original proposal but fixable in the BioCompiler framework. The fix is twofold: (1) model splicing as a non-deterministic process that produces a *set* of possible isoforms, not a single deterministic output (AD-03, NDFSTs); and (2) use three-valued logic (PASS/FAIL/UNCERTAIN) instead of binary pass/fail, so that the system can honestly represent situations where it cannot determine which isoform will be produced (AD-05). The key insight is that the *computation* is deterministic (the same input always produces the same isoform set) even though the *biological process* is non-deterministic (the same pre-mRNA can produce different isoforms in different cells or at different times).

#### Flaw 6: Grammar Explosion from Context-Sensitivity

| Attribute | Value |
|---|---|
| **Name** | Grammar Explosion from Context-Sensitivity |
| **Category** | Practical |
| **Severity** | Severe |
| **Can Be Fixed?** | Partially — by restricting to regular approximations (finite-state transducers) |

Biological processes involve long-range interactions: the inclusion of an exon depends on enhancer elements that may be thousands of nucleotides away; the folding of a protein depends on interactions between residues that are distant in the primary sequence but proximal in the tertiary structure; post-translational modifications depend on the cellular signaling environment. In a formal grammar, each such long-range dependency requires context-sensitive productions — rules whose applicability depends on the context surrounding the nonterminal being rewritten.

The number of context-sensitive productions required grows exponentially with the number of interacting elements. For a protein of length *n* with pairwise residue interactions, a context-sensitive grammar would need O(n²) productions just to encode the pairwise dependencies; for higher-order interactions (e.g., hydrophobic core formation involving dozens of residues simultaneously), the number of productions grows combinatorially. This is not merely a practical inconvenience — it makes the grammar intractable for any realistic protein, both in terms of specification (who writes millions of grammar rules?) and in terms of computation (parsing context-sensitive grammars is PSPACE-complete in general).

The BioCompiler fix is to restrict the splicing grammar to the regular language class (finite-state transducers), which cannot express arbitrary long-range dependencies but can express the local dependencies that are well-characterized in the splicing literature: splice site consensus sequences, branch point motifs, polypyrimidine tracts, and the presence/absence of nearby enhancer and silencer elements. Long-range dependencies (e.g., exon definition interactions across introns) are approximated by context-dependent thresholds that are parameterized by cellular context. This is a lossy approximation, but it is an honest one: the system explicitly reports UNCERTAIN when the grammar cannot resolve a dependency, rather than silently producing a wrong answer.

#### Flaw 7: IR Levels Are Hardcoded, Not Emergent

| Attribute | Value |
|---|---|
| **Name** | IR Levels Are Hardcoded, Not Emergent |
| **Category** | Conceptual |
| **Severity** | Moderate |
| **Can Be Fixed?** | Yes — by acknowledging that IR levels come from biology, not grammar induction |

The original proposal implied that the IR levels in the compilation pipeline (mRNA, peptide, structure, function) would be *discovered* by grammar induction — that the ADIOS algorithm would naturally segment the compilation process into stages, just as compiler researchers discovered the utility of SSA form and other IR levels through decades of practice. In reality, the IR levels in the BioCompiler pipeline are dictated by the known stages of molecular biology: transcription (DNA → pre-mRNA), splicing (pre-mRNA → mRNA), translation (mRNA → peptide), folding (peptide → structure), and post-translational modification (structure → functional protein). These stages are not discovered by the system; they are imposed by the designer's knowledge of biology.

This is a moderate flaw because it does not invalidate the pipeline architecture — the IR levels are still useful and well-motivated — but it does undermine the claim that grammar induction provides a principled method for discovering the architecture. The BioCompiler fix is simply to acknowledge the source of the IR levels: they come from biological knowledge, not from algorithmic discovery. This is consistent with the broader theme of replacing automatic discovery with curated knowledge (Flaw 3, Flaw 4). The pipeline architecture (AD-01) is still justified by the engineering benefits of modularity, testability, and extensibility — it just isn't justified by the claim of emergent structure.

#### Flaw 8: Inverse Compiler Is Ill-Posed

| Attribute | Value |
|---|---|
| **Name** | Inverse Compiler Is Ill-Posed |
| **Category** | Practical |
| **Severity** | Severe |
| **Can Be Fixed?** | Partially — by reframing as constraint satisfaction, not grammar inversion |

The original proposal promised an "inverse compiler" that, given a target protein, would produce the mRNA sequence that compiles to it. This is the biological analog of decompilation: given machine code, recover the source. The problem is that the forward mapping (mRNA → protein) is many-to-one due to the redundancy of the genetic code (61 sense codons encoding 20 amino acids) and the multiple possible splicing patterns that can produce the same protein. For a protein of length *n*, there are approximately 3^n possible mRNA sequences (since each amino acid is encoded by 2–6 synonymous codons, with an average of about 3). For a typical protein of 300 amino acids, this is approximately 10^143 possible mRNA sequences — a number that exceeds the number of atoms in the observable universe.

The original proposal treated inverse compilation as grammar inversion: given a grammar G and an output string w, find an input string v such that G(v) = w. For context-free grammars, this is already a hard problem; for context-sensitive grammars (which the proposal needed, per Flaw 4, Flaw 6), it is undecidable in general. Even for regular languages, the inverse problem has multiple solutions, and selecting among them requires an objective function — which brings us back to constraint satisfaction, not grammar inversion.

The BioCompiler fix is to reframe the inverse problem as a *constraint satisfaction problem* (AD-04). Instead of "inverting the grammar," the system finds an mRNA sequence that (a) translates to the target protein (hard constraint from the genetic code), (b) splices correctly under the specified cellular context (hard constraint from the splicing grammar), (c) satisfies codon adaptation thresholds (hard or soft constraint), (d) falls within GC content bounds (hard constraint), and (e) avoids specified restriction sites (hard constraint). This is a well-defined CSP with finite variable domains (the synonymous codons at each position) and well-defined constraints. The CSP solver either finds a feasible assignment or reports INFEASIBLE with a minimal unsatisfiable subset (MUS) explaining why. This is not grammar inversion; it is constraint satisfaction — and it is a much more tractable and useful formulation.

#### Flaw 9: Validation Is Circular

| Attribute | Value |
|---|---|
| **Name** | Validation Is Circular |
| **Category** | Methodological |
| **Severity** | Moderate |
| **Can Be Fixed?** | Yes — by using independent validation criteria |

The original proposal's validation strategy was essentially: "The grammar is correct if it parses the training data correctly." This is circular: the grammar was learned from the training data, so of course it parses it correctly. The real question is whether the grammar generalizes to unseen data — whether it correctly predicts the behavior of biological sequences not in the training set. Without independent validation data and well-defined validation criteria, the grammar's accuracy is untestable.

Parsing success (the grammar accepts a sequence) does not imply biological correctness (the sequence is actually processed as the grammar predicts). A grammar that over-approximates the true language will accept sequences that are never produced in biology; one that under-approximates will reject sequences that are biologically valid. Neither type of error is detectable by parsing alone.

The BioCompiler fix is to establish independent validation criteria for each component: the scanner is validated against curated motif databases (REBASE for restriction sites, GENCODE for splice sites); the splicing engine is validated against known isoform sets from GENCODE annotations; the translation engine is validated against the standard genetic code (a known, deterministic mapping); and the type system is validated by soundness arguments and adversarial testing. The guarantee certificates (COMP-07) provide an additional validation layer: they are independently verifiable by a separate checker program, so the correctness of the certificate does not depend on the correctness of the pipeline that produced it.

### 2.3 Summary Table

| # | Flaw | Category | Severity | Can Be Fixed? | BioCompiler Fix |
|---|---|---|---|---|---|
| 1 | Protein Folding Is Not String Rewriting | Representational | Fatal | No (within string-rewriting framework) | Delegate to FFI (AD-06) |
| 2 | PTMs Defy Formal Rules | Representational | Fatal | No (within sequence-grammar framework) | Delegate to FFI (AD-06) |
| 3 | ADIOS Learns Statistics, Not Physics | Methodological | Critical | Partially — use curated grammars | Hand-curated grammars (AD-07) |
| 4 | Grammar Induction for CSGs Is Undecidable | Methodological | Critical | Partially — restrict to regular approximations | NDFSTs, regular grammars (AD-03) |
| 5 | Determinism Assumption Is Fatal | Methodological | Fatal | Yes — model non-determinism explicitly | NDFSTs + three-valued logic (AD-03, AD-05) |
| 6 | Grammar Explosion from Context-Sensitivity | Practical | Severe | Partially — restrict to finite-state | NDFSTs (AD-03) |
| 7 | IR Levels Are Hardcoded, Not Emergent | Conceptual | Moderate | Yes — acknowledge the source | Pipeline with known stages (AD-01) |
| 8 | Inverse Compiler Is Ill-Posed | Practical | Severe | Partially — reframe as CSP | CSP + MUS (AD-04) |
| 9 | Validation Is Circular | Methodological | Moderate | Yes — use independent validation | Soundness arguments + independent verification (AD-05, COMP-07) |

---

## 3. What Remains Valuable

Despite the nine flaws, the original proposal contained genuine insights that survive the critical analysis and form the foundation of the BioCompiler design. These are not salvage operations — forcing broken ideas to work — but rather genuine contributions that the critical analysis validated rather than undermined.

### 3.1 The Pipeline Architecture Insight

The idea of processing gene sequences through a staged pipeline with typed intermediate representations is sound, independent of whether biology "implements" compilation. LLVM demonstrated that a multi-pass pipeline with typed IR enables independent development, testing, and replacement of passes; compositional analysis; and formal verification at each stage. These are engineering benefits that apply regardless of whether the metaphor of "compilation" is biologically accurate. The BioCompiler pipeline (AD-01) adopts this architecture not because genes compile into proteins, but because the engineering benefits of pipelines — modularity, testability, extensibility — are valuable for *any* complex transformation process.

The critical analysis actually *strengthened* the case for the pipeline architecture by showing that it provides the right abstraction for isolating the components that are formalizable (splicing, translation) from those that are not (folding, PTMs). In a monolithic system, the non-formalizable components would contaminate the formalizable ones; in a pipeline with typed IRs, the FFI boundary (AD-06) cleanly separates them, allowing the deterministic core to provide guarantees while the non-deterministic periphery provides enrichment.

### 3.2 Formal Language Theory for Splicing

The application of formal language theory to alternative splicing is genuinely underexplored in the bioinformatics literature, and the original proposal's instinct to use grammars for splicing was correct — it just applied the wrong kind of grammar to the wrong problem. Splicing *is* a language-recognized process: the spliceosome recognizes specific sequence patterns (donor consensus, branch point, acceptor consensus, regulatory elements) and performs a transformation (removing introns, joining exons) based on those patterns. A finite-state transducer is the right formalism for this: it captures the local pattern-matching nature of splice site recognition and the transformational nature of intron removal, and its non-deterministic variant (NDFST) naturally represents the multiple possible splicing outcomes that constitute alternative splicing.

The BioCompiler NDFST (AD-03) is a direct descendant of this insight, refined by the understanding that the grammar must be hand-curated (Flaw 3, Flaw 4), restricted to the regular language class (Flaw 6), and explicitly non-deterministic (Flaw 5). The formal language theory contribution is not the grand unifying grammar that the original proposal imagined, but a specific, well-motivated application of finite-state transducers to a well-characterized biological process.

### 3.3 A Protein IR as Common Infrastructure

The idea of a typed, schema-enforced intermediate representation as a common interchange format for bioinformatics tools is valuable and has no existing equivalent. Bioinformatics currently suffers from a proliferation of ad-hoc file formats (FASTA, FASTQ, GFF3, GenBank, PDB, mmCIF, BED, VCF, SAM/BAM) that are loosely typed, inconsistently parsed, and difficult to compose into pipelines. A typed IR — with schema enforcement (protocol buffers), versioning, and defined invariants — would provide the same infrastructure benefit that LLVM IR provides for compilers: a common language that tools can produce and consume without knowing about each other.

The BioCompiler IR Bus (AD-02) is a realization of this insight. Each IR level (IR-Seq, IR-Peptide, IR-Structure, IR-Circuit) has a defined schema, defined invariants, and defined producers and consumers. Tools that produce or consume these IR levels can interoperate without custom parsers or format conversion scripts. This is a genuinely novel contribution to bioinformatics infrastructure, independent of the compilation metaphor.

### 3.4 The Educational Bridge

The compilation metaphor, even when its limits are clearly understood, has significant pedagogical value. Students of molecular biology who have programming experience can understand the gene-to-protein pathway more easily when it is presented as a compilation process: DNA is "source code," transcription is "preprocessing," splicing is "optimization" (removing comments/introns), translation is "code generation" (producing the "machine code" of amino acids), and folding is "linking" (producing the executable structure). The metaphor is imperfect — as the nine flaws demonstrate — but as a pedagogical bridge, it provides an intuitive framework for understanding a complex process.

The BioCompiler documentation uses the compilation metaphor as a *design pattern* and *pedagogical tool*, not as a theoretical claim about biology. This distinction is critical: REQ-CON-004 explicitly states that the system shall not claim that biology implements compilation. The metaphor is useful for communicating the system's design to software engineers, but it is not a scientific hypothesis.

---

## 4. Why We Avoided Probability: Six Deterministic Methods

### 4.1 The Core Move

The central methodological innovation of BioCompiler is the shift from probabilistic questions to deterministic ones. Probabilistic questions ask "What is the probability that X?" Deterministic questions ask "Is X guaranteed, impossible, or undetermined?" This shift is not merely semantic — it changes the kind of answer the system can provide, the kind of guarantee it can make, and the kind of reasoning it can perform.

The following table illustrates the transformation for four key questions in gene design:

| Probabilistic Question | Deterministic Equivalent |
|---|---|
| "What is P(mRNA → protein)?" | "Is protein p in the set of all possible outputs from mRNA m?" |
| "What is P(splice site functional)?" | "Is this site guaranteed functional / guaranteed non-functional / uncertain?" |
| "Find mRNA maximizing P(correct) × P(expression)" | "Find mRNA satisfying all hard constraints" |
| "What is the distribution of folding outcomes?" | "What is the range of possible folding energies?" |

Each probabilistic question requires a probability model, training data, calibration, and independence assumptions. Each deterministic equivalent requires only the rules of the system (the splicing grammar, the genetic code, the constraint set) and produces a definitive answer or an honest "cannot determine." The deterministic answers compose without independence assumptions (the conjunction of guarantees is a guarantee), whereas probabilistic answers require independence or conditional independence assumptions that are typically unjustified in biological systems.

This does not mean that probabilistic methods are inherently inferior or that BioCompiler rejects them entirely. It means that BioCompiler provides a *different kind of answer* — one that is appropriate for safety-critical gene design, where a guarantee of correctness is more valuable than a probability estimate. When probability is needed (see §7.2), the system should decline to answer rather than give a deterministic answer that misleads.

### 4.2 Method 1: Abstract Interpretation

Abstract interpretation (Cousot & Cousot, 1977) is a theory of sound approximation of program behavior. Instead of executing a program on concrete inputs (which may be too expensive or impossible), abstract interpretation executes the program on *abstract* inputs that represent sets of concrete inputs, producing an *over-approximation* of the concrete behavior. The key property is soundness: if the abstract interpretation says a property holds, it holds for all concrete executions. If it says a property may not hold, it might be a false alarm (the property might actually hold for all concrete executions), but it will never produce a false negative.

In BioCompiler, abstract interpretation is applied to the splicing grammar. Instead of computing the exact set of isoforms for a given pre-mRNA sequence (which may be exponentially large for genes with many alternative exons), the system can compute an over-approximation: a superset of the true isoform set that is guaranteed to include every real isoform. This over-approximation is used for type checking: if the over-approximation satisfies a property (e.g., "no isoform contains a premature stop codon"), then the property is guaranteed for the exact isoform set. If the over-approximation violates a property, the system checks whether the violation is real (present in the exact set) or spurious (introduced by the approximation).

The three-valued logic (PASS/FAIL/UNCERTAIN) maps directly to abstract interpretation: PASS corresponds to "the abstract interpretation proves the property holds," FAIL corresponds to "the concrete behavior definitely violates the property," and UNCERTAIN corresponds to "the abstract interpretation cannot prove the property, but the concrete behavior might still satisfy it." This is precisely the YES/NO/MAYBE of abstract interpretation, instantiated in a biological context.

### 4.3 Method 2: Type Systems

Type systems in programming languages provide static verification without executing the program. Milner's famous slogan — "Well-typed programs don't go wrong" — captures the essence: a type system defines a set of "bad" behaviors (type errors), and any program that passes the type checker is guaranteed not to exhibit those behaviors. The type checker may reject some programs that would actually run correctly (false positives), but it will never accept a program that could go wrong (false negatives).

In BioCompiler, the type system (COMP-05) checks mRNA sequences against biological correctness properties: SpliceCorrect(CellType), NoCrypticSplice, CodonAdapted(Organism, threshold), GCInRange(lo, hi), NoRestrictionSite(EnzymeSet), InFrame, and NoInstabilityMotif. Each type predicate defines a "bad" behavior (incorrect splicing, cryptic splice sites, poor codon adaptation, GC content out of range, restriction site present, frame shift, instability motif present), and any sequence that passes the type checker is guaranteed not to exhibit that behavior. The type checker may reject some sequences that are actually correct (conservatism), but it will never accept a sequence that violates a stated property (soundness).

The type system approach is particularly appropriate for gene design because the "program" (the mRNA sequence) is executed exactly once per cell, and the "type error" (incorrect splicing, frame shift) has irreversible consequences. Static verification before execution is exactly what is needed.

### 4.4 Method 3: Constraint Satisfaction

Constraint satisfaction reformulates the gene design problem as a system of hard constraints over a set of variables with finite domains. Each codon position is a variable whose domain is the set of synonymous codons for the amino acid at that position. The constraints are: splicing correctness (no cryptic splice sites introduced by the chosen codon), codon adaptation (CAI above a threshold), GC content (within a specified range), restriction site absence (no match for any enzyme in the specified set), reading frame consistency, and instability motif absence.

The CSP approach provides three properties that no probabilistic method can match. First, *soundness*: every solution found by the solver satisfies all constraints, by construction. There is no "probably satisfies" — it either satisfies or it doesn't. Second, *completeness*: if a feasible assignment exists, the solver will find one (given sufficient time; the CSP is NP-hard in general, but the domain sizes are small enough that practical instances are tractable). Third, *diagnosis*: when the CSP is infeasible, the solver can compute a Minimal Unsatisfiable Subset (MUS) — the smallest subset of constraints that conflict. This tells the user *exactly which constraints are incompatible*, enabling targeted relaxation of constraints rather than blind trial-and-error.

The contrast with probabilistic optimization (e.g., genetic algorithms with a fitness function combining multiple weighted objectives) is stark. A genetic algorithm can find "good" solutions quickly but provides no guarantees: it cannot prove that no better solution exists, it cannot prove that the solution satisfies all constraints, and it cannot diagnose infeasibility. For safety-critical gene design, these are unacceptable limitations.

### 4.5 Method 4: Non-Deterministic Automata

Non-deterministic finite-state automata (NFAs) and transducers (NDFSTs) model processes that have multiple possible outcomes without assigning probabilities to those outcomes. An NFA accepts a set of strings; an NDFST produces a *set-valued* function mapping each input string to a set of possible output strings. The computation is deterministic (the same input always produces the same set of outputs), but the modeled process is non-deterministic (the "actual" output is some element of the set, but we don't know which one).

This is exactly the right model for alternative splicing. Given a pre-mRNA sequence, the NDFST produces the set of all possible splice isoforms consistent with the splicing grammar and the specified cellular context. We don't know which isoform will be produced in any given cell at any given time — that depends on stochastic splicing factor concentrations — but we *do* know that the produced isoform will be an element of the computed set. This is a guarantee: "the actual isoform is guaranteed to be in this set." It is not a probabilistic statement ("the actual isoform has probability 0.7 of being isoform A"), because we don't have a calibrated probability model for splicing factor concentrations.

The NDFST approach composes without probability. If the splicing NDFST produces isoform set S₁ for gene 1 and isoform set S₂ for gene 2, the combined isoform set for a circuit containing both genes is S₁ × S₂ (Cartesian product). No independence assumption is needed; the set-valued computation is exact. If we used probabilities instead, combining P(isoform for gene 1) and P(isoform for gene 2) would require assuming independence, which is unjustified when both genes are expressed in the same cell and compete for the same splicing machinery.

### 4.6 Method 5: Game Semantics

Game semantics models computation as a game between two players: the system (which makes choices) and the environment (which also makes choices). The system's goal is to satisfy a specification regardless of the environment's choices. This provides worst-case guarantees: even under adversarial conditions, the system's behavior is bounded.

In BioCompiler, game semantics is applied informally to the cellular context. The "system" is the designed mRNA sequence; the "environment" is the cellular context (splicing factor concentrations, kinase activities, chaperone availability). The system makes choices (codon assignment) to satisfy the specification (correct splicing, no cryptic sites, adequate CAI), while the "adversary" (the cellular context) makes choices that could violate the specification (activating a cryptic splice site by upregulating a splicing factor, phosphorylating an unexpected residue). The gene design "wins" if it satisfies the specification under all cellular contexts in the specified range; it "loses" if there exists a context in the range that causes a specification violation.

This framing explains why BioCompiler uses hard constraints rather than soft optimization: a gene design must be correct under all contexts in the specified range, not merely "optimal on average." A design that is 95% likely to splice correctly in the "average" cell is not acceptable if there is a 5% chance of incorrect splicing that produces a toxic protein. Game semantics formalizes this intuition: the design must win against all adversarial contexts, not just the average one.

### 4.7 Method 6: Three-Valued Logic

Three-valued logic (PASS/FAIL/UNCERTAIN) replaces the binary pass/fail of traditional verification with a ternary logic that explicitly represents uncertainty. PASS means the property is guaranteed to hold (there exists a proof). FAIL means the property is guaranteed to be violated (there exists a counterexample). UNCERTAIN means the system cannot determine whether the property holds or is violated (there is neither a proof nor a counterexample accessible to the system).

Three-valued logic has three key advantages over probabilistic scoring. First, it requires no calibration: PASS means PASS, not "PASS with probability 0.95." There is no need to train a model on labeled data to calibrate probability estimates, and no risk of miscalibration leading to false confidence. Second, it composes without independence assumptions: the conjunction of PASS and PASS is PASS (guaranteed to hold), the conjunction of PASS and UNCERTAIN is UNCERTAIN (the second property might be violated), and the conjunction of anything with FAIL is FAIL (the overall specification is violated). Third, it is honest about what the system knows and doesn't know: UNCERTAIN is not a euphemism for "probably fine" but a declaration of ignorance that the user can act on.

The disadvantage of three-valued logic is that UNCERTAIN is less actionable than a probability. A user who receives UNCERTAIN for a property knows that the system cannot guarantee correctness, but does not know whether the property is likely to hold or likely to be violated. This is an inherent limitation of deterministic methods: they cannot provide information they don't have. The BioCompiler design principle is that it is better to be honestly uncertain than to be spuriously precise.

### 4.8 Comparison Table

The following table compares three approaches across key dimensions:

| Dimension | Original Paper (Fake Determinism) | Revised Framework (Probability) | This Document (True Determinism) |
|---|---|---|---|
| **Folding prediction** | Grammar-based (string rewriting) | Probabilistic ML (AlphaFold) | External FFI (AD-06); non-deterministic SLOT |
| **Splicing model** | Deterministic CFG (single isoform) | PCFG or HMM (probability per isoform) | NDFST (set-valued, no probabilities) |
| **Gene design** | Grammar inversion | Multi-objective weighted optimization | CSP with hard constraints (AD-04) |
| **Validation verdict** | Binary (parse/no-parse) | Probability score (0–1) | Three-valued (PASS/FAIL/UNCERTAIN) |
| **Composition** | Sequential (one isoform) | Requires independence assumptions | Set-valued composition (no independence) |
| **Guarantee type** | None (circular validation) | Probabilistic (calibrated?) | Formal (soundness proof) |
| **Diagnosis on failure** | None | Sensitivity analysis | MUS (minimal unsatisfiable subset) |
| **Handling of uncertainty** | Ignored (assumed deterministic) | Quantified (probability estimates) | Explicit (UNCERTAIN verdict) |
| **Calibration required?** | N/A | Yes | No |
| **Training data required?** | Yes (for ADIOS) | Yes (for ML models) | No (curated knowledge) |
| **False positive rate** | Unknown | Depends on calibration | Zero (soundness) |
| **False negative rate** | Unknown | Depends on calibration | Non-zero (conservatism) |

---

## 5. Architectural Decision Records

> **Note:** Each ADR below is also available as a standalone Nygard-format file in the `adr/` directory:
> - [ADR-0001: Pipeline Architecture](adr/ADR-0001.md)
> - [ADR-0002: Protocol Buffers for IR Schemas](adr/ADR-0002.md)
> - [ADR-0003: Non-Deterministic Finite-State Transducers for Splicing](adr/ADR-0003.md)
> - [ADR-0004: Constraint Satisfaction for Gene Design](adr/ADR-0004.md)
> - [ADR-0005: Three-Valued Logic for Verdicts](adr/ADR-0005.md)
> - [ADR-0006: Foreign Function Interface for Folding and PTMs](adr/ADR-0006.md)
> - [ADR-0007: Declarative Grammar Configuration](adr/ADR-0007.md)
>
> The standalone files follow the Michael Nygard ADR format (Title / Status / Date / Context / Decision / Consequences). The sections below provide the same decisions in the ISO 42010 §7 format used throughout this document, with additional narrative context linking each decision to the critical analysis and deterministic methods framework.

### AD-01: Pipeline Architecture

| Field | Value |
|---|---|
| **Decision** | Use a staged pipeline architecture with typed IRs (not monolith, not microservices) |
| **Status** | Accepted |
| **Date** | 2026-05-30 |

**Context:** The BioCompiler system needs to process gene sequences through multiple transformation stages — scanning, splicing, translation, type checking, optimization, and certification — with well-defined interfaces between stages. The stages have different computational characteristics: some are purely symbolic (scanning, translation), some are non-deterministic (splicing), and some require external tools (folding, PTM prediction). The architecture must support independent development, testing, and replacement of stages while maintaining end-to-end correctness guarantees.

**Alternatives Considered:**

(a) **Monolithic tool**: A single program that performs all stages in sequence, with internal data structures passed between functions. This is simpler to implement initially and avoids the overhead of IR serialization and deserialization. However, it couples all stages together, making it impossible to test stages independently, replace a stage with an improved version, or skip stages that are not needed for a particular use case. It also makes it difficult to isolate the non-deterministic stages (FFI) from the deterministic core.

(b) **Microservices**: Each stage runs as an independent service with a network API, communicating via HTTP or gRPC. This provides maximum deployment flexibility and language independence, but it introduces network latency, serialization overhead, and operational complexity (service discovery, health monitoring, fault tolerance) that is disproportionate for a batch processing system that runs on a single machine.

(c) **Pipeline with typed IR** (chosen): Each stage is a pure function that consumes and produces well-typed intermediate representations through an IR Bus. Stages communicate exclusively through IR records; no stage directly accesses another stage's internal state. This provides the modularity and testability of microservices without the operational overhead, and the typing discipline of a monolith without the coupling.

**Rationale:** The pipeline architecture mirrors LLVM's proven multi-pass design (Lattner & Adve, 2004), which demonstrated that a staged pipeline with typed IR enables: (1) independent development of passes by different teams, (2) independent testing of each pass against its input/output contract, (3) replacement of passes without affecting others, (4) compositional analysis where each pass adds information to the IR, and (5) formal verification at each stage via IR invariants. These benefits are directly applicable to BioCompiler. The pipeline architecture also provides the right abstraction for isolating formalizable stages (splicing, translation) from non-formalizable ones (folding, PTMs) behind the FFI boundary, which is essential given Flaw #1 and Flaw #2.

**Consequences:**

- *Positive*: (1) Each stage is independently testable — the scanner can be tested with fixture inputs without running the splicing engine, and vice versa. (2) Stages can be replaced — if a better splicing engine is developed, it can be swapped in without modifying any other stage. (3) IR invariants provide formal checkpoints — if an invariant is violated, the violating stage is immediately identified. (4) The FFI boundary cleanly separates deterministic and non-deterministic computation. (5) New analysis passes can be added as new consumers of existing IR levels without modifying producers.

- *Negative*: (1) More upfront design effort is required for IR schemas — each IR level must be specified before stages can be developed in parallel. (2) IR serialization and deserialization add overhead compared to in-process function calls, though this is negligible for the data sizes involved (sequences of thousands to millions of nucleotides). (3) Schema evolution requires careful version management (see AD-02). (4) The pipeline architecture implies a fixed ordering of stages, which may not be optimal for all use cases — some use cases might benefit from iterative refinement loops between stages, which the pipeline does not naturally support (the optimizer loop is a special case handled by COMP-06).

---

### AD-02: Protocol Buffers for IR Schemas

| Field | Value |
|---|---|
| **Decision** | Define IR schemas in Protocol Buffers (.proto files) |
| **Status** | Accepted |
| **Date** | 2026-05-30 |

**Context:** The IR Bus requires schemas that define the structure of each IR level (IR-Seq, IR-Peptide, IR-Structure, IR-Circuit). These schemas must be enforced at runtime, support backward-compatible evolution, support efficient serialization for persistence and inter-process communication, and ideally support code generation for multiple programming languages.

**Alternatives Considered:**

(a) **JSON Schema only**: JSON is human-readable and widely supported, but JSON Schema is a validation language, not a code generation language. It does not produce typed data classes in any programming language, requiring manual implementation of parsing, validation, and serialization. It also lacks efficient binary serialization — JSON is verbose and slow to parse for large datasets. Furthermore, JSON Schema's typing discipline is weaker than Protocol Buffers: optional fields, repeated fields, and nested messages are clumsier to specify.

(b) **HDF5**: The Hierarchical Data Format is optimized for large numerical arrays and is widely used in scientific computing. However, HDF5 is poorly suited for the structured, record-oriented data that IR levels represent. HDF5's data model is based on multidimensional arrays and groups of arrays, not on typed records with fields. It also lacks code generation, schema evolution, and multi-language support comparable to Protocol Buffers. HDF5 would be appropriate for large sequence matrices but not for the IR metadata and annotations that BioCompiler requires.

(c) **Custom binary format**: Maximum flexibility and potentially the most efficient serialization, but also maximum maintenance burden. A custom format requires writing serializers, deserializers, and validators for every language the system needs to support, and there is no ecosystem of tools (schema registries, compatibility checkers, migration utilities) that mature formats provide.

**Rationale:** Protocol Buffers provide the best combination of schema enforcement, backward compatibility, efficient serialization, and multi-language code generation. Schema enforcement is automatic: the generated Python classes have typed fields and raise errors on schema violations. Backward compatibility is built into the format: new fields can be added without breaking existing consumers (unknown fields are preserved). Efficient binary serialization reduces IR size by 3–10× compared to JSON. Code generation for Python (and potentially C++, Java, Go) eliminates manual parsing and validation. The proto3 dialect (used by BioCompiler) simplifies the schema language by removing required fields and default values, which reduces the risk of schema evolution errors.

**Consequences:**

- *Positive*: (1) Type safety at the language level — generated Python classes have typed fields that mypy can check. (2) Backward-compatible schema evolution — new IR levels can add fields without breaking existing stages. (3) Efficient binary serialization — smaller IR files and faster I/O. (4) Multi-language code generation — the same schema can produce Python, C++, and Java data classes. (5) Schema registry — proto files serve as a single source of truth for IR structure.

- *Negative*: (1) Build step required — `protoc` must be run to generate Python stubs from `.proto` files before the system can be built. This adds complexity to the build process and requires developers to remember to regenerate stubs after schema changes. (2) Proto3 limitations — the proto3 dialect does not support `optional` field modifiers (until recent versions), which makes it harder to distinguish between "field not set" and "field set to default value." This is mitigated by using wrapper messages or `oneof` for optional fields. (3) Limited support for complex data structures — Protocol Buffers do not natively support algebraic data types or recursive types, which are occasionally useful for representing parse trees and derivation traces. These must be approximated with nested messages and `oneof` fields.

---

### AD-03: Non-Deterministic FSTs for Splicing

| Field | Value |
|---|---|
| **Decision** | Model splicing as a Non-Deterministic Finite-State Transducer (NDFST) producing set-valued output (no probabilities) |
| **Status** | Accepted |
| **Date** | 2026-05-30 |

**Context:** The splicing stage must model alternative splicing — the process by which a single pre-mRNA sequence can produce multiple distinct mRNA isoforms through different combinations of exon inclusion and exclusion. The model must capture the essential non-determinism of the process (multiple valid isoforms for the same input) while maintaining deterministic computation (the same input always produces the same isoform set). The model must also support cellular context parameterization (different cell types produce different isoform sets due to different splicing factor concentrations).

**Alternatives Considered:**

(a) **Probabilistic Context-Free Grammar (PCFG)**: Assigns probabilities to production rules, producing a probability distribution over isoforms. This is the standard approach in computational linguistics and has been applied to RNA secondary structure prediction. However, PCFGs require training data to estimate rule probabilities, which are context-dependent (the probability of an exon being included depends on the cellular context, not just the sequence). The probabilities would need to be re-estimated for each cell type, which is infeasible given the sparsity of cell-type-specific isoform quantification data. Furthermore, PCFGs do not provide the guarantees that BioCompiler requires: a 95% probability that an isoform is the only one produced is not a guarantee.

(b) **Hidden Markov Model (HMM)**: Assigns emission and transition probabilities to states, producing a probability distribution over state sequences (isoforms). HMMs are well-understood for gene finding (GENSCAN, AUGUSTUS) but share the same limitations as PCFGs: they require training data, produce probabilities rather than guarantees, and do not compose without independence assumptions. They also have difficulty modeling long-range dependencies in splicing (e.g., exon definition interactions), which require higher-order HMMs with exponential state spaces.

(c) **Neural network**: A sequence-to-sequence model (e.g., transformer) trained on pre-mRNA → isoform pairs. This is the most flexible approach and can potentially capture arbitrary long-range dependencies, but it is a black box: it provides no formal guarantees, no interpretability, and no composability. It also requires large training datasets that do not exist for most cell types, and its predictions cannot be verified independently of the model.

(d) **NDFST** (chosen): A finite-state transducer that maps each input string to a *set* of possible output strings, with no probabilities assigned to the set elements. The computation is deterministic (same input → same isoform set), but the modeled process is non-deterministic (the actual isoform is some element of the set). The set-valued computation captures the essential non-determinism of alternative splicing without requiring probability estimates, training data, or independence assumptions.

**Rationale:** The NDFST is the right formalism for splicing for three reasons. First, it captures the essential non-determinism: alternative splicing produces a *set* of possible isoforms, not a single isoform or a probability distribution over isoforms. Second, it avoids the need for probability estimates: the system computes the set of possible isoforms and checks whether all of them satisfy the desired properties (or whether any of them violate them). Third, it composes: the combined isoform set for a circuit of multiple genes is the Cartesian product of the individual isoform sets, with no independence assumptions. The cellular context parameterizes the NDFST by enabling or disabling transitions based on splicing factor concentration thresholds, which is a deterministic modulation of a non-deterministic machine.

**Consequences:**

- *Positive*: (1) Formal guarantees — the NDFST is sound (every produced isoform satisfies the grammar) and complete (every isoform satisfying the grammar is produced). (2) No training data required — the grammar rules are hand-curated from biological knowledge, not learned from data. (3) Composable — isoform sets compose via Cartesian product without independence assumptions. (4) Cellular context parameterization — thresholds modulate which transitions are active, producing different isoform sets for different cell types. (5) Efficient computation — finite-state transducers can be implemented in O(n × s) time where n is the input length and s is the number of states, using subset construction for the non-deterministic case.

- *Negative*: (1) Cannot rank isoforms by likelihood — the NDFST produces a set, not a distribution, so it cannot tell the user which isoform is most likely to be produced. This limits its usefulness for applications that require quantitative predictions (see §7.2). (2) Conservative — the isoform set is an over-approximation (it may include rare or non-functional isoforms that the grammar cannot rule out). This means that a SpliceCorrect verdict (the target isoform is the only one produced) may be FAIL even if the target isoform is the dominant one, because a rare alternative isoform exists. (3) Limited expressiveness — finite-state transducers cannot capture arbitrary long-range dependencies, so some biologically relevant splicing interactions (e.g., exon definition interactions across introns) are only approximately modeled.

---

### AD-04: Constraint Satisfaction for Gene Design

| Field | Value |
|---|---|
| **Decision** | Use Constraint Satisfaction Problem (CSP) with hard constraints rather than weighted optimization |
| **Status** | Accepted |
| **Date** | 2026-05-30 |

**Context:** The gene design problem is to find an mRNA sequence that translates to a target protein and satisfies multiple biological correctness constraints: splicing correctness (no cryptic splice sites, correct splicing pattern), codon adaptation (CAI above a threshold), GC content (within a specified range), restriction site absence, reading frame consistency, and instability motif absence. The problem is an inverse problem: given a target protein and a set of constraints, find an mRNA sequence that satisfies all constraints.

**Alternatives Considered:**

(a) **Multi-objective optimization**: Define a weighted objective function combining multiple criteria (CAI, GC content, splice site score, etc.) and use an optimizer (e.g., gradient descent, CMA-ES) to find a sequence that maximizes the weighted objective. This is the approach taken by most existing gene design tools (DNAWorks, GeneDesign, OPTIMIZER). The fundamental problem is that weighting factors are arbitrary: there is no principled way to determine whether "CAI ≥ 0.8" is worth 2× or 3× as much as "GC content ∈ [40%, 60%]." Different weighting choices produce different "optimal" designs, and the user has no way to know whether a different weighting would produce a better design for their application. More critically, weighted optimization provides no guarantees: it cannot prove that a solution satisfying all constraints exists, and it cannot prove that no such solution exists.

(b) **Genetic algorithm**: Use a population of candidate sequences, evaluate fitness (a weighted combination of criteria), and evolve the population over generations. Genetic algorithms are widely used in bioinformatics for sequence optimization because they can explore large search spaces without gradient information. However, they provide no guarantees: they cannot prove that a better solution does not exist, they cannot prove that all constraints are satisfied (only that the fitness function assigns a high score), and they cannot diagnose infeasibility (they simply fail to converge). For safety-critical gene design, these are unacceptable limitations.

(c) **Simulated annealing**: Start with a random sequence and iteratively make random changes, accepting improvements and occasionally accepting deteriorations (with decreasing probability). Simulated annealing has the same limitations as genetic algorithms: no guarantees, no completeness, no diagnosis.

(d) **CSP with MUS** (chosen): Formulate the gene design problem as a Constraint Satisfaction Problem with finite variable domains (the synonymous codons at each position) and hard constraints (all constraints must be satisfied exactly). Use a CSP solver (based on AC-3 constraint propagation and backtracking search) to find a feasible assignment or prove that no feasible assignment exists. When infeasible, compute a Minimal Unsatisfiable Subset (MUS) to identify the smallest set of conflicting constraints.

**Rationale:** Safety-critical gene design requires hard guarantees, not soft scores. A gene designed for therapeutic use must be *guaranteed* to splice correctly, not merely *likely* to splice correctly. The CSP formulation provides three properties that no optimization approach can match. First, *soundness*: every solution found by the CSP solver satisfies all constraints by construction. There is no "approximately satisfies" — the solution either satisfies every constraint exactly or it is not returned. Second, *completeness*: if a feasible assignment exists, the solver will find one (given sufficient time; the CSP is NP-hard in general, but biological instances have small domain sizes — 1 to 6 synonymous codons per position — making them tractable). Third, *diagnosis*: the MUS computation tells the user exactly which constraints conflict, enabling targeted constraint relaxation rather than blind parameter tuning.

The CSP formulation also supports the three-valued logic of the type system: a PASS verdict from the type system translates directly to a set of constraints that the CSP solver must satisfy; a FAIL verdict adds constraints (the violated property must be fixed); and an UNCERTAIN verdict identifies constraints that cannot be formally verified. This integration between the type system and the CSP solver is only possible because both operate on hard constraints, not on weighted objectives.

**Consequences:**

- *Positive*: (1) Formal guarantees — every solution satisfies all constraints; no false positives. (2) Diagnosis — MUS identifies the minimal conflicting constraint set, enabling targeted constraint relaxation. (3) Integration with type system — type-check results translate directly to CSP constraints. (4) Determinism — the solver is deterministic (same input → same solution or same INFEASIBLE report). (5) No arbitrary weighting — constraints are hard, not weighted, so the user does not need to choose weighting factors.

- *Negative*: (1) May report INFEASIBLE when a "good enough" solution exists — the hard-constraint formulation is less flexible than weighted optimization. If the user specifies CAI ≥ 0.9 and the maximum achievable CAI (while satisfying all other constraints) is 0.87, the CSP solver reports INFEASIBLE, whereas a weighted optimizer would return the 0.87 solution. This can be partially mitigated by iteratively relaxing constraints guided by the MUS. (2) Scalability — the CSP is NP-hard in general, though biological instances are typically tractable due to small domain sizes and sparse constraint interactions. For genes with many alternative codons and complex splicing constraints, the solver may require significant computation time. (3) Single-objective — the CSP formulation optimizes a single scalar objective (typically CAI) subject to hard constraints. Multi-objective optimization (Pareto-optimal designs trading off CAI against GC content, for example) requires multiple CSP solves with different constraint bounds, which is less elegant than a single multi-objective optimization.

---

### AD-05: Three-Valued Logic for Verdicts

| Field | Value |
|---|---|
| **Decision** | Use PASS/FAIL/UNCERTAIN instead of probabilities or binary pass/fail |
| **Status** | Accepted |
| **Date** | 2026-05-30 |

**Context:** The type system must produce verdicts for each biological correctness property checked against an mRNA sequence. The verdict must honestly represent the system's state of knowledge: guaranteed correct, guaranteed incorrect, or cannot determine. The verdict must also compose: when multiple properties are checked, the combined verdict must accurately represent the combined state of knowledge.

**Alternatives Considered:**

(a) **Probability scores**: Assign a probability to each property (e.g., "P(splicing correct) = 0.93") and use a threshold to convert to pass/fail. This is the approach taken by most bioinformatics tools (MaxEntScan for splice site scoring, NetGene2 for splice site prediction). The fundamental problems are: (1) probability scores require calibration — a score of 0.93 from one model does not mean the same thing as 0.93 from another model, and uncalibrated scores are misleading; (2) probabilities do not compose — P(A and B) ≠ P(A) × P(B) unless A and B are independent, which biological properties typically are not; (3) probability scores give a false sense of precision — a score of 0.93 suggests that the property is "probably" correct, but without calibration, the true probability could be 0.5 or 0.99.

(b) **Binary pass/fail**: Each property is either satisfied or not. This is simpler than three-valued logic but hides uncertainty. A property that cannot be verified (e.g., the splicing grammar cannot determine whether a weak splice site is functional in the specified cell type) would be forced into PASS or FAIL, both of which are misleading. If forced to PASS, the system claims a guarantee it cannot make; if forced to FAIL, the system rejects a design that might be correct. Neither outcome is honest.

(c) **Confidence intervals**: Report a range (e.g., "splicing correctness: 85%–98%") instead of a point estimate. This is more informative than a point probability but still probabilistic — it requires a statistical model, calibration, and independence assumptions for composition. It also does not map cleanly to the yes/no decision that gene design requires: "Is this design acceptable or not?"

**Rationale:** Three-valued logic captures exactly the information that BioCompiler has: PASS means the system can *prove* the property holds (there exists a derivation trace), FAIL means the system can *prove* the property is violated (there exists a counterexample), and UNCERTAIN means the system has *insufficient information* to determine either way. This is not a simplification or an approximation — it is a precise representation of the system's epistemic state. Three-valued logic composes without independence assumptions: PASS ∧ PASS = PASS (both properties are guaranteed), PASS ∧ UNCERTAIN = UNCERTAIN (the second property might not hold), and anything ∧ FAIL = FAIL (at least one property is violated). No calibration is required because the verdicts are not probabilities. No training data is required because the verdicts are computed from the rules of the system, not from statistical models.

The three-valued logic also maps directly to abstract interpretation (Method 1, §4.2): PASS corresponds to the abstract interpretation proving the property, FAIL corresponds to the concrete behavior violating the property, and UNCERTAIN corresponds to the abstract interpretation being too imprecise to determine the property. This connection provides a theoretical foundation for the three-valued logic that goes beyond mere engineering convenience.

**Consequences:**

- *Positive*: (1) Honest about uncertainty — UNCERTAIN is not a euphemism for "probably fine" but a declaration of ignorance that the user can act on. (2) Composable without independence assumptions — the conjunction and disjunction of three-valued verdicts preserve soundness. (3) No calibration required — PASS means guaranteed correct, not "estimated 95% likely to be correct." (4) Maps to abstract interpretation — the three-valued logic has a well-understood theoretical foundation. (5) Supports the CSP solver — PASS constraints are added to the CSP, FAIL constraints trigger the optimizer, and UNCERTAIN constraints are flagged for the user.

- *Negative*: (1) "UNCERTAIN" is less actionable than a probability — a user who receives UNCERTAIN for a property knows that the system cannot guarantee correctness, but does not know whether the property is likely to hold or likely to be violated. The user must decide whether to accept the uncertainty or invest in resolving it (e.g., by providing additional information or using a more precise analysis). (2) May be too conservative — the system may report UNCERTAIN for properties that are "almost certainly" correct, leading to unnecessary rejection of valid designs. This conservatism is intentional (it preserves soundness) but may frustrate users in exploratory research contexts where some risk is acceptable. (3) Does not support quantitative comparison — two designs that both have all-PASS verdicts are indistinguishable, even if one is "more robust" than the other. The CSP's scalar objective (CAI) partially addresses this for the optimization dimension, but not for other dimensions of design quality.

---

### AD-06: FFI for Folding/PTMs

| Field | Value |
|---|---|
| **Decision** | Invoke folding and PTM prediction as external tools through a Foreign Function Interface (FFI), not internal models |
| **Status** | Accepted |
| **Date** | 2026-05-30 |

**Context:** Protein folding and post-translational modification prediction are essential for full gene-to-protein analysis, but they are not formalizable as string transformations (Flaw #1, Flaw #2). The system must incorporate their predictions without compromising the determinism and formal guarantees of the core pipeline. The FFI boundary must cleanly separate the deterministic core (which provides guarantees) from the non-deterministic periphery (which provides enrichment).

**Alternatives Considered:**

(a) **Internal ML model**: Train or integrate a machine learning model (e.g., a transformer-based folding model, a CNN-based PTM predictor) as an internal component of the pipeline. This would eliminate the need for external tool installation and could be optimized for the pipeline's specific use cases. However, it would add a probabilistic component to the core pipeline, violating the determinism principle (REQ-CON-001). It would also require training data, model maintenance, and regular retraining as new data becomes available, which is a significant ongoing engineering burden. Most critically, it would undermine the system's ability to provide formal guarantees: a guarantee that depends on a probabilistic model is not a formal guarantee.

(b) **Internal physics simulation**: Implement a molecular dynamics simulation or energy minimization algorithm as an internal component. This would be the most principled approach (folding is a physical process, so a physics simulation is the most accurate model), but it is computationally infeasible: a single protein folding simulation takes hours to days on specialized hardware, whereas BioCompiler must process sequences in seconds to minutes. Even with recent advances in GPU-accelerated MD, the computational cost is orders of magnitude too high for a tool that must process hundreds of candidate designs.

(c) **FFI to existing tools** (chosen): Use the best-in-class external tools for folding (AlphaFold2/3, ColabFold, RoseTTAFold) and PTM prediction (NetPhos, PhosphoSitePlus, dbPTM, MusiteDeep) through a defined adapter interface. The FFI boundary isolates the non-determinism of these tools from the deterministic core: FFI output is placed in SLOT fields in the IR, which the type system treats as non-deterministic (UNCERTAIN verdicts for properties that depend on FFI output). No guarantee certificate depends on FFI output for its core validity — certificates are based on deterministic properties (splicing, translation, codon usage, GC content, restriction sites) that the core pipeline can verify.

**Rationale:** Folding and PTMs are not formalizable as string transformations (Flaw #1, Flaw #2), so they cannot be part of the deterministic core. The FFI boundary is the right abstraction for this situation: it allows the system to use the best available tools for these tasks while maintaining clean separation between the deterministic core (which provides guarantees) and the non-deterministic periphery (which provides enrichment). The SLOT mechanism (empty fields in the IR that are filled by FFI adapters) ensures that FFI output is clearly labeled as non-deterministic and cannot be confused with deterministic analysis results. The type system's treatment of SLOT-dependent properties as UNCERTAIN ensures that no false guarantees are issued based on FFI output.

**Consequences:**

- *Positive*: (1) Uses best-in-class tools — the system always has access to the latest and most accurate folding and PTM prediction tools without reimplementing them. (2) Clean separation of paradigms — the deterministic core and the non-deterministic periphery are separated by a well-defined boundary with clear semantics. (3) No internal probabilistic models — the core pipeline remains fully deterministic, preserving formal guarantees. (4) Extensible — new external tools can be supported by writing new FFI adapters without modifying the core pipeline. (5) The FFI boundary can be independently tested — each adapter has a defined contract (format_input, invoke, parse_output, validate_output) that can be tested in isolation.

- *Negative*: (1) Dependency on external tools — the system requires users to install and configure external tools (AlphaFold, NetPhos, etc.) to use FFI features. This adds installation complexity and version management burden. The core pipeline operates without FFI tools (REQ-NFR-042), but the enriched analysis requires them. (2) No guarantees on FFI output — the system cannot verify that AlphaFold's folding prediction is correct; it can only verify that the output conforms to the expected schema and invariants. If the folding prediction is wrong, the system will not detect the error (it will report the prediction with its confidence score and flag the dependent properties as UNCERTAIN). (3) Version management — external tools evolve independently, and different versions may produce different outputs for the same input. The FFI adapter records the tool version in the provenance metadata, but cannot ensure reproducibility across tool versions. (4) FFI invocation overhead — each external tool invocation requires subprocess creation, input/output serialization, and output parsing, which adds latency compared to an internal model. This is mitigated by parallelizing FFI calls across isoforms (see §3.3 of DOC-02).

---

### AD-07: Declarative Grammar Configuration

| Field | Value |
|---|---|
| **Decision** | Splicing grammar rules in YAML configuration files, not hardcoded in source code |
| **Status** | Accepted |
| **Date** | 2026-05-30 |

**Context:** The splicing grammar defines the rules that the NDFST uses to parse pre-mRNA sequences: splice donor and acceptor consensus sequences, branch point motifs, polypyrimidine tract requirements, exonic and intronic splicing enhancer and silencer patterns, and cellular context thresholds. These rules are scientific knowledge that evolves as new splice sites are annotated, new regulatory elements are discovered, and new cell-type-specific splicing data becomes available.

**Alternatives Considered:**

(a) **Hardcoded rules**: Embed the splicing grammar directly in the source code as constant data structures (Python dictionaries, lists, or dataclasses). This is the simplest approach: no config file parsing, no runtime validation, no risk of malformed configuration. However, it requires recompilation and redeployment whenever a grammar rule changes, which is unacceptable for a tool used by researchers who need to experiment with different rule sets. It also prevents community contributions — a researcher who discovers a new splicing regulatory element cannot add it to the grammar without modifying the source code.

(b) **Learned rules (grammar induction)**: Use the ADIOS algorithm or another grammar induction method to learn grammar rules from biological sequence data. This was the original proposal's approach, but it is rejected for the reasons documented in Flaw #3 (ADIOS learns statistics, not physics) and Flaw #4 (grammar induction for context-sensitive languages is undecidable). Even if a restricted form of grammar induction could be made to work, the learned grammar would not be interpretable — the researcher could not inspect the rules to understand why the system makes a particular prediction — which defeats the purpose of a tool that must provide transparent, verifiable guarantees.

(c) **Declarative YAML configuration** (chosen): Specify grammar rules in YAML configuration files that are loaded at pipeline initialization time. Each organism has its own configuration file (e.g., `homo_sapiens.yaml`, `e_coli.yaml`, `saccharomyces.yaml`) containing splice site consensus sequences, branch point patterns, polypyrimidine tract requirements, regulatory element definitions, and cell-type-specific thresholds. The configuration is immutable during a pipeline run but can be updated between runs without recompilation.

**Rationale:** Splicing rules are scientific knowledge, not engineering parameters. They are derived from experimental observations (GENCODE annotations, published splicing assays, high-throughput splicing screens) and evolve as new data becomes available. YAML configuration enables updates without recompilation, community contributions (researchers can submit pull requests with updated grammar files for their organism of interest), and transparency (the rules are human-readable and can be inspected to understand the system's behavior). The declarative format also supports validation: the config loader can check that all required fields are present, that thresholds are within valid ranges, and that consensus sequences are valid IUPAC patterns, before the pipeline starts.

The YAML format was chosen over alternatives (JSON, TOML, XML) for its readability: YAML supports comments (essential for documenting the scientific basis of each rule), multi-line strings (useful for consensus sequences and notes), and a clean syntax that non-programmers (biologists) can understand and modify. The main risk — runtime parsing errors due to malformed YAML — is mitigated by strict schema validation at load time: the pipeline fails fast with a clear error message if the configuration is invalid, rather than silently producing wrong results.

**Consequences:**

- *Positive*: (1) Updatable without recompilation — new splicing rules can be added or existing rules modified by editing YAML files, without changing the source code. (2) Community-friendly — researchers can contribute grammar files for their organisms of interest via pull requests, lowering the barrier to community contributions. (3) Transparent — the rules are human-readable, enabling users to understand and audit the system's behavior. (4) Validated at load time — the config loader checks the configuration against a schema before the pipeline starts, preventing runtime errors from malformed rules. (5) Organism-specific — each organism has its own configuration file, enabling organism-specific grammar rules without code changes.

- *Negative*: (1) Runtime parsing overhead — YAML parsing adds a small overhead at pipeline initialization time. This is negligible compared to the computation time of the pipeline stages (seconds to minutes for typical inputs). (2) Configuration validation complexity — the config loader must validate that all required fields are present, that thresholds are within valid ranges, that consensus sequences are valid IUPAC patterns, and that cell-type-specific parameters are consistent with the organism's biology. This validation code must be maintained as the configuration schema evolves. (3) Risk of user error — a user who manually edits a YAML file may introduce errors (typos, invalid values, missing fields) that cause the pipeline to fail or produce incorrect results. The schema validation mitigates this risk but cannot eliminate it entirely. (4) No version control for runtime configurations — the pipeline uses the YAML file as-is at load time; if the file is modified between pipeline runs, the results may differ. This is managed by recording the configuration hash in the provenance metadata of the guarantee certificate.

---

## 6. Novel Capabilities Enabled by These Decisions

The seven architectural decisions, individually motivated by the need to address the nine flaws and adopt deterministic methods, collectively enable capabilities that no existing bioinformatics tool provides. These capabilities are emergent properties of the design — they arise from the interaction of multiple decisions, not from any single decision alone.

### 6.1 Splicing-Aware Gene Design

No existing gene design tool combines codon optimization with splicing grammar enforcement. Current tools (DNAWorks, GeneDesign, OPTIMIZER, IDT Codon Optimization) optimize codon usage for expression (maximizing CAI, controlling GC content, avoiding restriction sites) but treat splicing as an afterthought — they may check for cryptic splice sites after optimization, but they do not *prevent* cryptic splice sites during optimization. The result is that codon-optimized genes frequently introduce cryptic splice sites that cause aberrant splicing in mammalian cells, a well-documented problem in the gene therapy and mRNA vaccine literature.

BioCompiler's combination of NDFST-based splicing analysis (AD-03) and CSP-based optimization (AD-04) enables splicing-aware gene design: the optimizer's constraint set includes splicing correctness constraints (no cryptic splice sites, correct splicing pattern) alongside codon usage constraints (CAI, GC content, restriction sites). The CSP solver finds codon assignments that satisfy all constraints simultaneously, eliminating the need for a manual "optimize-then-check-then-fix" cycle. This capability is enabled by the integration of AD-03 (NDFST provides the splicing constraint function), AD-04 (CSP solver handles the joint constraint set), and AD-05 (three-valued logic provides the constraint language).

### 6.2 Typed IR for Cross-Tool Interoperability

No LLVM-IR equivalent exists for bioinformatics. Current bioinformatics pipelines are assembled from tools that communicate through ad-hoc file formats (FASTA, GFF3, PDB, etc.) with inconsistent typing, no schema enforcement, and no version management. Pipeline integrators spend significant effort writing format conversion scripts and debugging format incompatibilities.

BioCompiler's typed IR (AD-02) provides a common interchange format with schema enforcement (Protocol Buffers), versioning (semver), and defined invariants. Tools that produce or consume IR levels can interoperate without custom parsers or format conversion scripts. A folding tool that produces IR-Structure can be consumed by a type checker that expects IR-Structure, without either tool knowing about the other. This is the same infrastructure benefit that LLVM IR provides for compilers, and it is enabled by AD-01 (pipeline architecture with IR Bus), AD-02 (Protocol Buffers for IR schemas), and AD-06 (FFI adapters produce IR-conformant output).

### 6.3 Compositional Verification of Genetic Circuits

No existing tool performs splicing-aware compositional verification of multi-gene circuits. Current circuit verification tools (Cello, SBOL) focus on transcriptional regulation (promoter-repressor interactions) but do not check for post-transcriptional interactions: splicing interference (cryptic splice sites in one gene affecting another's splicing), RNA-RNA interactions (complementary transcript regions forming dsRNA), or resource competition (ribosome demand exceeding capacity). These interactions can cause circuit failure even when each individual gene is correctly designed.

BioCompiler's compositional verifier (COMP-08) checks four categories of cross-component constraints: promoter conflict, resource competition, splicing interference, and RNA-RNA interaction. Each check produces a three-valued verdict that composes with the per-gene verdicts. The circuit certificate (COMP-07) includes both per-gene certificates and composition check results, enabling independent verification without re-running the pipeline. This capability is enabled by AD-01 (pipeline architecture supporting per-gene and circuit-level analysis), AD-03 (NDFST provides the splicing model for interference detection), AD-05 (three-valued logic composes across genes), and AD-07 (declarative grammar enables organism-specific composition rules).

### 6.4 Grammar-Guided Mutation Landscape Exploration

No existing tool systematically explores the mutation space of a gene using the structure of the splicing grammar to decompose and prune the search space. Current approaches to mutation analysis either enumerate all possible point mutations (exponential in sequence length) or sample random mutations (no coverage guarantees). Neither approach exploits the structure of the splicing grammar to identify which mutations are relevant (those near splice sites or regulatory elements) and which are independent (those in different exons that do not interact).

BioCompiler's mutation explorer (COMP-09) decomposes the mutation space by splicing grammar nonterminals (intra-exonic, splice site, regulatory element), enumerates legal multi-mutation combinations within each category, exploits independence across exons (mutations in different exons can be combined freely), and detects constraint conflicts (mutations that are individually legal but jointly violate the splicing grammar). This systematic exploration provides guarantees that random sampling cannot: every legal mutation within the specified edit distance is found, and every conflict is detected. This capability is enabled by AD-03 (NDFST provides the grammar structure for decomposition) and AD-04 (CSP provides the constraint checking for conflict detection).

### 6.5 Overlapping Reading Frame Analysis

No existing tool computes shared constraint sets for overlapping reading frames in viral and compact genomes. Current tools (Geneious, SnapGene) can visualize overlapping ORFs but cannot compute the coupling between frames (which nucleotide positions affect multiple proteins) or detect constraint conflicts between frames (where the optimal codon for one frame is suboptimal for another).

BioCompiler's ORF analyzer (COMP-10) computes shared constraint sets, classifies positions as high-coupling (affecting multiple proteins) or low-coupling (affecting one protein), and detects constraint conflicts between frames. This information is essential for rational vaccine design (where mutations must not disrupt overlapping reading frames) and for understanding the evolutionary constraints on compact genomes. This capability is enabled by AD-03 (NDFST for splicing analysis of each frame), AD-04 (CSP for joint constraint optimization across frames), and AD-05 (three-valued logic for cross-frame conflict verdicts).

### 6.6 Proof-Carrying Gene Designs

No existing gene design tool provides machine-checkable certificates for its designs. Current tools produce output sequences with summary statistics (CAI, GC content) but no independently verifiable proof that the design satisfies all stated constraints. A regulatory reviewer who receives a gene design from BioCompiler can independently verify the guarantee certificate using the standalone verifier (COMP-07), without trusting the BioCompiler pipeline or its developers. The certificate contains the verified sequence, all type predicates with PASS verdicts and derivation traces, the CSP constraint set and assignment, and provenance metadata. The verifier re-runs each predicate against the sequence and checks that each derivation step follows from the previous one.

This capability is the culmination of the entire design rationale: it depends on AD-01 (pipeline architecture with IR invariants), AD-02 (typed IR with schema enforcement), AD-03 (NDFST for sound splicing analysis), AD-04 (CSP for provably correct optimization), AD-05 (three-valued logic for composable verdicts), and AD-06 (FFI isolation ensuring that core guarantees do not depend on non-deterministic tools). Proof-carrying gene designs are the primary differentiator of BioCompiler relative to all existing gene design tools, and they are only possible because every component of the pipeline provides a formal guarantee that can be independently verified.

---

## 7. Honest Limitations

### 7.1 What This Framework Cannot Do

The deterministic, guarantee-oriented design of BioCompiler comes with inherent limitations that are not engineering deficiencies but fundamental consequences of the chosen approach. Recording these limitations honestly is essential for preventing misuse of the system and for setting appropriate expectations.

**Cannot predict protein structures.** Protein folding is a thermodynamic process (Flaw #1) that is delegated to external tools through the FFI boundary (AD-06). BioCompiler cannot guarantee that a designed mRNA will produce a correctly folded protein; it can only guarantee that the mRNA will be correctly spliced and translated. Users who need folding guarantees must use the FFI tools (AlphaFold, ColabFold) and accept their probabilistic predictions, which BioCompiler treats as non-deterministic (UNCERTAIN verdicts).

**Cannot predict PTM patterns.** Post-translational modifications depend on cellular context in ways that are not formalizable as sequence transformations (Flaw #2). BioCompiler cannot guarantee that a designed protein will not undergo unwanted PTMs; it can only invoke external PTM prediction tools through the FFI boundary and flag the results as non-deterministic.

**Cannot answer probabilistic questions.** The system is designed for deterministic guarantees, not probability estimates. It cannot answer questions like "What fraction of transcripts will include exon 5?" or "What is the probability that this splice site will be used?" These questions require a probabilistic model of splicing, which BioCompiler deliberately does not include (see §4.1). The NDFST produces a *set* of possible isoforms, not a distribution over isoforms.

**Cannot compare designs quantitatively.** Two designs that both receive all-PASS verdicts are indistinguishable from the guarantee perspective, even if one is "more robust" than the other (e.g., has a higher CAI, is further from the GC content boundary, or has weaker cryptic splice sites). The CSP's scalar objective (CAI) partially addresses this for the optimization dimension, but BioCompiler does not provide a general framework for quantitative design comparison.

**May be too conservative for exploratory research.** The three-valued logic (AD-05) and hard-constraint CSP (AD-04) are designed for safety-critical applications where false positives (guaranteeing a property that does not hold) are unacceptable. In exploratory research, where some risk is acceptable and the goal is to generate candidate designs for experimental testing, the conservatism may be counterproductive: the system may reject designs that would work in practice because it cannot formally guarantee their correctness.

### 7.2 When Probability Is Unavoidable

Some biologically important questions are inherently probabilistic and cannot be answered by deterministic methods without misrepresentation. These include:

- **Isoform quantification**: "What fraction of transcripts from this gene will include exon 5?" This requires a probabilistic model of splicing regulation, which depends on splicing factor concentrations that vary across cells and over time.

- **Expression level prediction**: "How much protein will this gene produce?" This depends on transcription rate, mRNA stability, translation rate, and protein degradation rate, all of which are stochastic processes with unknown distributions.

- **Comparative ranking**: "Is design A better than design B?" This requires a scalar metric that combines multiple dimensions of design quality, which in turn requires weighting factors that are inherently subjective.

For questions requiring quantitative predictions, probability is unavoidable. The BioCompiler framework should *explicitly decline to answer these questions* rather than giving deterministic answers that mislead. The three-valued logic supports this: when the system cannot provide a guarantee, it reports UNCERTAIN, which is an honest statement of the system's epistemic limits. A user who needs a probability should use a probabilistic tool (e.g., a splicing quantification model, an expression prediction tool), not BioCompiler. The two tools answer different questions and provide different kinds of information; they are complementary, not competing.

---

## 8. Design Principles Summary

The following table maps each design principle to the architectural decisions that enforce it and the requirements that mandate it. This table serves as a quick reference for understanding how abstract principles are realized in concrete engineering choices.

| Principle | Decisions That Enforce It | Requirements That Mandate It | Description |
|---|---|---|---|
| **Determinism by Construction** | AD-01 (pipeline with pure-function stages), AD-03 (NDFST — deterministic computation of non-deterministic output), AD-04 (CSP solver is deterministic), AD-05 (three-valued logic is deterministic), AD-06 (FFI isolates non-determinism), AD-07 (declarative config is immutable during a run) | REQ-CON-001 (no probabilistic models for internal stages), REQ-NFR-001 (bit-identical output for identical input) | The pipeline produces the same output for the same input on every run. No randomness, no sampling, no probabilistic inference in any internal component. |
| **Soundness by Proof** | AD-03 (NDFST is sound — every produced isoform satisfies the grammar), AD-04 (CSP solver produces only feasible solutions), AD-05 (PASS verdict requires a derivation trace), AD-06 (FFI output is treated as non-deterministic, not as proof) | REQ-NFR-002 (no false PASS verdicts), INV-TYP-01 (soundness invariant for type system) | No PASS verdict is issued for a property that does not hold. Every PASS is backed by a machine-checkable derivation trace. |
| **Separation of Paradigms** | AD-01 (pipeline architecture isolates stages), AD-06 (FFI boundary isolates non-deterministic tools), AD-03 (NDFST isolates non-deterministic biology within a deterministic computation) | REQ-CON-001 (no probabilistic models), REQ-CON-003 (no internal folding/PTM models) | Symbolic (grammar-based) and continuous (ML/physics) computations never mix in the same component. The FFI boundary is the membrane that separates them. |
| **Honesty About Uncertainty** | AD-05 (three-valued logic — PASS/FAIL/UNCERTAIN), AD-03 (NDFST produces sets, not distributions), AD-06 (FFI output flagged as non-deterministic) | REQ-CON-004 (no claim that biology implements compilation), §7.2 (probability is sometimes unavoidable) | The system reports UNCERTAIN when it cannot determine a property, rather than forcing a PASS or FAIL verdict. The system declines to answer probabilistic questions rather than giving misleading deterministic answers. |
| **Composability** | AD-01 (pipeline with typed IR enables composition of stages), AD-03 (NDFST isoform sets compose via Cartesian product), AD-05 (three-valued logic composes via defined algebra), AD-04 (CSP constraints compose additively) | REQ-NFR-007 (compositional verification of circuits) | Multi-gene circuits are verifiable from per-gene certificates without re-running the pipeline. Verdicts compose via a defined algebra without independence assumptions. |
| **Extensibility** | AD-01 (pipeline stages are independently replaceable), AD-02 (Protocol Buffers support schema evolution), AD-06 (new FFI adapters can be added without modifying the core), AD-07 (new grammar rules can be added via config files) | REQ-NFR-009 (extensibility for new organisms, folding algorithms, constraint classes) | New organisms, folding algorithms, and constraint classes can be added without modifying the core pipeline. FFI adapters and YAML configs are the primary extension points. |
| **Verifiability** | AD-05 (derivation traces for every PASS verdict), AD-04 (CSP assignments can be independently verified), AD-02 (typed IR enables independent checking), AD-06 (FFI provenance recorded) | REQ-NFR-004 (independently verifiable guarantee certificates) | Every PASS verdict is accompanied by a certificate that is independently verifiable by a separate checker program without access to the BioCompiler pipeline. |
| **Diagnosis on Failure** | AD-04 (MUS for infeasible CSP), AD-05 (violation identification for FAIL, knowledge gap specification for UNCERTAIN), AD-03 (parse path identification for splicing failures) | REQ-NFR-006 (diagnostic failure reports) | When the system fails (INFEASIBLE, FAIL, UNCERTAIN), it provides a precise diagnosis identifying the root cause, not just a binary "design rejected" message. |

---

## 9. v7.1 Design Decisions: GT-Free Codon Prioritization and CpG Avoidance

Version 7.1 introduces three interrelated design decisions that significantly improve the optimizer's ability to satisfy the NoCrypticSplice and NoCpGIsland predicates. These decisions were driven by empirical observation of optimization failures on real gene sequences and represent a refinement of the greedy optimizer's phase architecture (ADR-0008).

### 9.1 GT-Free Codon Prioritization (ADR-011)

**Problem**: Phase 7 (cryptic splice elimination) was failing on ~80% of genes because it did not prioritize GT-free codons for amino acids that have them. The GT dinucleotide is the core recognition sequence for splice donors. Any codon containing GT is a potential cryptic splice donor. For Valine, ALL four codons (GTT, GTC, GTA, GTG) contain GT — no codon swap can eliminate it. But for Alanine, Glycine, Arginine, and Serine, GT-free alternatives exist and provide a guaranteed path to eliminating the cryptic donor.

**Why naive codon swapping fails**: The previous Phase 7 implementation attempted codon swaps without considering whether GT-free alternatives existed. This led to futile cycling: for a Cysteine position where both codons (TGC, TGT) contain GT, the optimizer would try both and fail, while for an Alanine position, it might miss the obvious fix (swap GCG→GCC or GCA or GCT). The optimizer treated all cryptic donor positions equally, regardless of whether a guaranteed fix was available.

**The 3-strategy approach**:

1. **GT-free codon swap** (highest priority): For amino acids with GT-free alternatives, swap to the highest-CAI GT-free codon. This is a guaranteed elimination — no search, no backtracking. For example, Alanine's GCC (GT-free, highest human CAI) is always preferred over GCG (contains GT).

2. **Context disruption** (Valine only): Since no V codon eliminates GT, the optimizer instead tries to reduce the MaxEntScan score by selecting the V codon that produces the weakest splice site in its 9-mer context, or disrupts the context by swapping neighboring codons.

3. **Accept unrepairable**: Valine positions that remain above threshold are flagged with `gt_mandatory=True` in the type system's derivation, directing the mutagenesis engine to propose V→I substitutions.

**Rationale for priority ordering**: Strategy 1 is tried first because it is guaranteed to work and has minimal CAI impact. Strategy 2 is tried second because it may work but is not guaranteed. Strategy 3 is not an optimizer action — it delegates to the mutagenesis engine, which is a separate concern.

### 9.2 CpG Avoidance Phase (ADR-012)

**Problem**: The NoCpGIsland predicate was failing on many optimized sequences because the optimizer had no phase to avoid CG dinucleotides. High-CAI codon selection in human often favors GC-rich codons (GCC, GGC, CGC), which naturally create CpG dinucleotides both within codons and at codon boundaries.

**Why a dedicated phase is necessary**: CpG avoidance cannot be incorporated into the CAI maximization phase without fundamentally changing the greedy algorithm's objective. The CAI phase optimizes for codon adaptation; adding a secondary CpG avoidance objective would require multi-objective optimization that sacrifices CAI unnecessarily. Many CpG positions can be fixed post-hoc with minimal or zero CAI impact, making a dedicated phase more efficient.

**Phase ordering rationale**: Phase 7.5 (CpG avoidance) runs after Phase 7 (cryptic splice elimination) because CpG disruption must not reintroduce cryptic splice donors. Phase 8.5 (reconciliation) runs after Phase 8 (final GC adjustment) to ensure CpG fixes don't undo restriction site removal or GC content corrections.

**Best-effort nature**: The CpG avoidance phase is best-effort — not all CG dinucleotides can be eliminated without changing the amino acid sequence. Arginine codons (CGN) all contain CG, and the alternative AGA/AGG codons may violate other constraints. The phase reports unrepairable positions rather than silently accepting them.

### 9.3 GT-Mandatory vs Optimizer Weakness Distinction (ADR-013)

**Problem**: The mutagenesis engine was proposing amino acid substitutions for positions where the optimizer should have fixed the problem by choosing a GT-free codon. This conflated two fundamentally different issues: (1) GT-mandatory positions where no codon swap can help (Valine), and (2) optimizer weaknesses where GT-free codons exist but weren't used.

**Why the distinction matters**: If the mutagenesis engine proposes a substitution for an optimizer weakness, it masks the real problem. The optimizer bug remains unfixed, and the protein is unnecessarily modified. By distinguishing GT-mandatory from optimizer weaknesses, the system ensures that:

- Mutagenesis is only proposed for positions where it is truly necessary
- Optimizer bugs are surfaced for repair rather than hidden
- Protein identity is preserved when possible
- The type system's derivation provides diagnostic information about WHY a position fails

**Implementation**: The `is_gt_mandatory(aa)` function classifies each amino acid. Currently, only Valine is GT-mandatory. The `diagnose_optimizer_weakness()` function identifies positions where the optimizer failed to use available GT-free codons, providing `(position, current_codon, gt_free_alternatives)` tuples for debugging. The type system's derivation now includes `gt_mandatory` and `gt_free_alternatives` fields, enriching the diagnostic information available to both automated tools and human reviewers.

### 9.4 Interaction Between v7.1 Decisions

The three v7.1 decisions form a coherent improvement to the optimizer-mutagenesis pipeline:

1. ADR-011 (GT-free codon prioritization) fixes the optimizer for non-Valine amino acids — these positions are now resolved by Phase 7, not by mutagenesis.
2. ADR-012 (CpG avoidance) adds a new optimization capability that was entirely missing — the optimizer now actively disrupts CpG dinucleotides.
3. ADR-013 (GT-mandatory distinction) ensures that the mutagenesis engine only acts on positions that the optimizer truly cannot fix (Valine), rather than masking optimizer bugs with unnecessary protein modifications.

The net effect is: more constraints satisfied by the optimizer, fewer protein modifications by the mutagenesis engine, and clearer diagnostic information when constraints cannot be satisfied.

---

## Appendix A: Relationship to DOC-09 and DOC-10

This document integrates and synthesizes the analyses from DOC-09 (Critical Analysis of Original Framework) and DOC-10 (Deterministic Methods for Non-Deterministic Biology), but it does not replace them. DOC-09 provides the full technical argument for each of the nine flaws, including mathematical proofs where applicable, experimental evidence from the literature, and detailed comparison with existing bioinformatics tools. DOC-10 provides the full formal treatment of each of the six deterministic methods, including formal definitions, soundness proofs, composition theorems, and worked examples. This document provides the *engineering rationale* — why the flaws led to specific architectural decisions, and why the deterministic methods led to specific design choices — which is a different kind of argument that neither DOC-09 nor DOC-10 alone provides.

Readers who need the technical details behind a flaw or a method should consult DOC-09 or DOC-10 directly. Readers who need to understand why a design decision was made should consult this document first and then follow the cross-references to DOC-09 and DOC-10 for deeper analysis.

---

## Appendix B: Terminology Cross-Reference

| Term in This Document | Equivalent in DOC-01 (SRS) | Equivalent in DOC-02 (SAD) | Notes |
|---|---|---|---|
| NDFST | REQ-FUNC-013 through REQ-FUNC-022 | COMP-02 | Non-Deterministic Finite-State Transducer |
| Three-valued logic | REQ-FUNC-030 through REQ-FUNC-039 | COMP-05 | PASS / FAIL / UNCERTAIN |
| CSP | REQ-FUNC-040 through REQ-FUNC-050 | COMP-06 | Constraint Satisfaction Problem |
| MUS | REQ-FUNC-046, REQ-FUNC-047 | COMP-06 | Minimal Unsatisfiable Subset |
| FFI | REQ-FUNC-051 through REQ-FUNC-058 | COMP-04 | Foreign Function Interface |
| IR Bus | REQ-CON-005 through REQ-CON-008 | §2.3 | Intermediate Representation Bus |
| SLOT | REQ-FUNC-053, REQ-FUNC-054 | COMP-04 | Empty field in IR filled by FFI adapter |
| Guarantee Certificate | REQ-FUNC-059 through REQ-FUNC-065 | COMP-07 | Machine-checkable JSON certificate |
| Splicing Grammar | REQ-FUNC-013 through REQ-FUNC-022 | COMP-02 | Grammar rules loaded from YAML config (AD-07) |
| Compositional Verification | REQ-FUNC-066 through REQ-FUNC-073 | COMP-08 | Cross-gene constraint checking |
| Mutation Explorer | REQ-FUNC-074 through REQ-FUNC-080 | COMP-09 | Grammar-guided mutation enumeration |
| ORF Analyzer | REQ-FUNC-081 through REQ-FUNC-088 | COMP-10 | Overlapping reading frame analysis |

---

*End of DOC-06: Design Rationale — Version 1.0.0-draft*
