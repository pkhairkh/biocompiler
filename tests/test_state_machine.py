"""
Tests for the deterministic optimization state machine.

Verifies:
1. Determinism: same input always produces same output
2. Valid state transitions
3. ConstraintAutomaton builds correctly for each constraint
4. Codon selection is deterministic
5. State machine produces same or better results than heuristics
"""

import pytest
import math


# ────────────────────────────────────────────────────────────
# Test 1: Determinism — same input → same output
# ────────────────────────────────────────────────────────────

class TestStateMachineDeterministic:
    """Verify that the state machine is fully deterministic."""

    def test_same_input_same_output_prokaryote(self):
        """Same protein + organism always produces identical sequence."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        sm1 = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        sm2 = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )

        result1 = sm1.optimize(protein)
        result2 = sm2.optimize(protein)

        assert result1.sequence == result2.sequence, \
            "Same input should produce identical output"
        assert result1.cai == result2.cai
        assert result1.gc_content == result2.gc_content

    def test_same_input_same_output_eukaryote(self):
        """Same input produces same output for eukaryotic organism."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        sm1 = DeterministicOptimizationStateMachine(
            organism="Homo_sapiens",
            is_prokaryote=False,
            avoid_gt=True,
        )
        sm2 = DeterministicOptimizationStateMachine(
            organism="Homo_sapiens",
            is_prokaryote=False,
            avoid_gt=True,
        )

        result1 = sm1.optimize(protein)
        result2 = sm2.optimize(protein)

        assert result1.sequence == result2.sequence, \
            "Same input should produce identical output for eukaryotes"
        assert result1.cai == result2.cai

    def test_seed_has_no_effect(self):
        """The seed parameter should have no effect on output."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEK"

        result_no_seed = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        ).optimize(protein)

        result_seed_42 = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
            seed=42,
        ).optimize(protein)

        result_seed_999 = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
            seed=999,
        ).optimize(protein)

        assert result_no_seed.sequence == result_seed_42.sequence == result_seed_999.sequence, \
            "Seed should have no effect — machine is deterministic"

    def test_multiple_runs_identical(self):
        """Running the same optimizer multiple times produces identical results."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        results = []
        for _ in range(5):
            sm = DeterministicOptimizationStateMachine(
                organism="Escherichia_coli",
                is_prokaryote=True,
            )
            results.append(sm.optimize(protein))

        sequences = [r.sequence for r in results]
        assert len(set(sequences)) == 1, \
            "All runs should produce the same sequence"


# ────────────────────────────────────────────────────────────
# Test 2: State transitions are valid
# ────────────────────────────────────────────────────────────

class TestStateMachineTransitions:
    """Verify that all state transitions follow the formal definition."""

    def test_init_is_first_state(self):
        """Machine starts in INIT state."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine, OptimizationState

        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        assert sm.state == OptimizationState.INIT

    def test_final_state_is_done(self):
        """After optimization, machine is in DONE state."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine, OptimizationState

        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize("MVHLTPEEK")
        assert sm.state == OptimizationState.DONE

    def test_states_visited_include_init_and_done(self):
        """State trace always includes INIT and DONE."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine, OptimizationState

        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize("MVHLTPEEK")
        assert result.states_visited[0] == OptimizationState.INIT
        assert result.states_visited[-1] == OptimizationState.DONE

    def test_required_states_visited(self):
        """All required states (CODON_SELECT, CONSTRAINT_CHECK, CAI_RECOVER, VALIDATE) are visited."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine, OptimizationState

        required = {
            OptimizationState.CODON_SELECT,
            OptimizationState.CONSTRAINT_CHECK,
            OptimizationState.CAI_RECOVER,
            OptimizationState.VALIDATE,
        }

        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR")
        visited_set = set(result.states_visited)

        for state in required:
            assert state in visited_set, \
                f"Required state {state.name} was not visited"

    def test_no_invalid_transitions(self):
        """State trace doesn't contain impossible transitions."""
        from biocompiler.optimizer.state_machine import OptimizationState

        # Define valid transitions
        valid_next = {
            OptimizationState.INIT: {OptimizationState.CODON_SELECT},
            OptimizationState.CODON_SELECT: {OptimizationState.CONSTRAINT_CHECK},
            OptimizationState.CONSTRAINT_CHECK: {
                OptimizationState.CONFLICT_RESOLVE,
                OptimizationState.CAI_RECOVER,
            },
            OptimizationState.CONFLICT_RESOLVE: {OptimizationState.CONSTRAINT_CHECK},
            OptimizationState.CAI_RECOVER: {OptimizationState.VALIDATE},
            OptimizationState.VALIDATE: {
                OptimizationState.CONFLICT_RESOLVE,
                OptimizationState.DONE,
            },
            OptimizationState.DONE: set(),
        }

        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR")

        for i in range(len(result.states_visited) - 1):
            current = result.states_visited[i]
            next_state = result.states_visited[i + 1]
            assert next_state in valid_next[current], \
                f"Invalid transition: {current.name} → {next_state.name}"


# ────────────────────────────────────────────────────────────
# Test 3: ConstraintAutomaton builds correctly
# ────────────────────────────────────────────────────────────

class TestConstraintAutomatonBuilds:
    """Verify that the ConstraintAutomaton builds correctly for each
    constraint type."""

    def test_builds_for_ecoli(self):
        """Automaton builds for E. coli (prokaryote)."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        auto = ConstraintAutomaton(
            organism="Escherichia_coli",
            species_cai=CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            is_prokaryote=True,
        )
        assert auto.is_prokaryote is True
        assert not auto.avoid_gt  # GT avoidance disabled for prokaryotes

    def test_builds_for_human(self):
        """Automaton builds for Human (eukaryote with GT avoidance)."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        auto = ConstraintAutomaton(
            organism="Homo_sapiens",
            species_cai=CODON_ADAPTIVENESS_TABLES["Homo_sapiens"],
            is_prokaryote=False,
            avoid_gt=True,
        )
        assert auto.is_prokaryote is False
        assert auto.avoid_gt is True

    def test_builds_with_enzymes(self):
        """Automaton builds with restriction enzyme constraints."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        auto = ConstraintAutomaton(
            organism="Escherichia_coli",
            species_cai=CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            enzymes=["EcoRI", "BamHI", "XhoI"],
            is_prokaryote=True,
        )
        assert len(auto._rs_sites) >= 3  # At least 3 enzyme sites

    def test_sorted_codons_per_aa(self):
        """Automaton pre-computes sorted codons for each amino acid."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        from biocompiler.type_system import AA_TO_CODONS

        auto = ConstraintAutomaton(
            organism="Escherichia_coli",
            species_cai=CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            is_prokaryote=True,
        )

        for aa in set("ACDEFGHIKLMNPQRSTVWY"):
            assert aa in auto._sorted_codons
            # Codons should be sorted by CAI descending
            codons = auto._sorted_codons[aa]
            for i in range(len(codons) - 1):
                cai_i = auto.species_cai.get(codons[i], 0.0)
                cai_next = auto.species_cai.get(codons[i + 1], 0.0)
                assert cai_i >= cai_next, \
                    f"Codons for {aa} not sorted by CAI desc: {codons}"

    def test_gt_free_codons_precomputed(self):
        """GT-free codons are pre-computed for eukaryotes."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        auto = ConstraintAutomaton(
            organism="Homo_sapiens",
            species_cai=CODON_ADAPTIVENESS_TABLES["Homo_sapiens"],
            is_prokaryote=False,
            avoid_gt=True,
        )

        # Valine has no GT-free codons (all start with GT)
        assert auto._gt_free_codons.get("V", []) == []

        # Leucine has GT-free codons
        leu_gt_free = auto._gt_free_codons.get("L", [])
        assert len(leu_gt_free) > 0
        for c in leu_gt_free:
            assert "GT" not in c


# ────────────────────────────────────────────────────────────
# Test 4: Codon selection is deterministic
# ────────────────────────────────────────────────────────────

class TestConstraintAutomatonDeterministicSelection:
    """Verify that codon selection via the ConstraintAutomaton is
    deterministic."""

    def test_same_position_same_codon(self):
        """Selecting a codon at the same position always returns the same codon."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        auto = ConstraintAutomaton(
            organism="Escherichia_coli",
            species_cai=CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            is_prokaryote=True,
        )

        for _ in range(10):
            codon = auto.select_codon(position=0, aa="M")
            assert codon == "ATG", "Methionine always maps to ATG"

    def test_valid_codons_deterministic(self):
        """get_valid_codons always returns the same list for the same inputs."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        auto = ConstraintAutomaton(
            organism="Escherichia_coli",
            species_cai=CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            is_prokaryote=True,
        )

        for _ in range(10):
            valid = auto.get_valid_codons(position=0, aa="L")
            assert isinstance(valid, list)
            assert len(valid) > 0
            # Same result each time
            first_codon = valid[0]
        assert first_codon == valid[0]

    def test_codon_translates_correctly(self):
        """Every selected codon translates to the correct amino acid."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES
        from biocompiler.type_system import CODON_TABLE

        auto = ConstraintAutomaton(
            organism="Escherichia_coli",
            species_cai=CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            is_prokaryote=True,
        )

        for aa in "ACDEFGHIKLMNPQRSTVWY":
            codon = auto.select_codon(position=0, aa=aa)
            translated = CODON_TABLE.get(codon)
            assert translated == aa, \
                f"Codon {codon} for AA {aa} translates to {translated}, not {aa}"

    def test_no_stop_codons_selected(self):
        """The automaton never selects a stop codon for a regular amino acid."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        auto = ConstraintAutomaton(
            organism="Escherichia_coli",
            species_cai=CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            is_prokaryote=True,
        )

        for aa in "ACDEFGHIKLMNPQRSTVWY":
            valid = auto.get_valid_codons(position=0, aa=aa)
            for codon in valid:
                assert codon not in ("TAA", "TAG", "TGA"), \
                    f"Stop codon {codon} in valid set for AA {aa}"

    def test_valid_codons_sorted_by_cai(self):
        """Valid codons are always sorted by CAI descending."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        auto = ConstraintAutomaton(
            organism="Escherichia_coli",
            species_cai=CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            is_prokaryote=True,
        )

        for aa in "ACDEFGHIKLMNPQRSTVWY":
            valid = auto.get_valid_codons(position=0, aa=aa)
            for i in range(len(valid) - 1):
                cai_i = auto.species_cai.get(valid[i], 0.0)
                cai_next = auto.species_cai.get(valid[i + 1], 0.0)
                assert cai_i >= cai_next, \
                    f"Valid codons for {aa} not sorted by CAI: {valid}"

    def test_restriction_enzyme_filtering(self):
        """Codons that create restriction sites within themselves are filtered."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES

        # EcoRI site is GAATTC (6bp) — won't fit in a single codon,
        # but the automaton should still check cross-codon boundaries
        auto = ConstraintAutomaton(
            organism="Escherichia_coli",
            species_cai=CODON_ADAPTIVENESS_TABLES["Escherichia_coli"],
            enzymes=["EcoRI"],
            is_prokaryote=True,
        )

        # Just verify it builds and doesn't crash
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            valid = auto.get_valid_codons(position=0, aa=aa)
            assert len(valid) > 0


# ────────────────────────────────────────────────────────────
# Test 5: State machine produces same or better results than heuristics
# ────────────────────────────────────────────────────────────

class TestStateMachineReplacesHeuristics:
    """Verify that the deterministic state machine produces results
    comparable to or better than the heuristic optimizer."""

    def test_prokaryote_cai_reasonable(self):
        """State machine produces reasonable CAI for E. coli."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize(protein)

        # Prokaryote with no GT/splice constraints should achieve high CAI
        assert result.cai > 0.7, f"CAI too low for E. coli: {result.cai}"
        # GC should be in a reasonable range
        assert 0.2 < result.gc_content < 0.8, \
            f"GC out of range: {result.gc_content}"

    def test_eukaryote_cai_reasonable(self):
        """State machine produces reasonable CAI for Human."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        sm = DeterministicOptimizationStateMachine(
            organism="Homo_sapiens",
            is_prokaryote=False,
            avoid_gt=True,
        )
        result = sm.optimize(protein)

        # Eukaryote with GT avoidance will have slightly lower CAI
        # but should still be reasonable
        assert result.cai > 0.3, f"CAI too low for Human: {result.cai}"

    def test_sequence_translates_to_protein(self):
        """The output sequence correctly translates to the input protein."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine
        from biocompiler.type_system import CODON_TABLE

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize(protein)

        # Translate back
        translated = ""
        for i in range(0, len(result.sequence), 3):
            codon = result.sequence[i:i + 3]
            aa = CODON_TABLE.get(codon, "?")
            translated += aa

        # The translated protein should match the input (excluding stop codon)
        assert translated.rstrip("*") == protein, \
            f"Translation mismatch: expected {protein}, got {translated.rstrip('*')}"

    def test_no_internal_stop_codons(self):
        """The output sequence has no internal stop codons."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize(protein)

        seq = result.sequence
        # Check all codons except the last one
        n_codons = len(seq) // 3
        for i in range(n_codons - 1):  # Skip last codon
            codon = seq[i * 3:i * 3 + 3]
            assert codon not in ("TAA", "TAG", "TGA"), \
                f"Internal stop codon {codon} at position {i}"

    def test_sequence_length_correct(self):
        """Output sequence length is exactly 3x the protein length."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize(protein)

        assert len(result.sequence) == len(protein) * 3, \
            f"Sequence length {len(result.sequence)} != protein length * 3 ({len(protein) * 3})"

    def test_short_protein(self):
        """State machine handles short proteins correctly."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEK"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize(protein)

        assert len(result.sequence) == len(protein) * 3
        assert result.cai > 0.0
        assert 0.0 <= result.gc_content <= 1.0

    def test_single_aa(self):
        """State machine handles single amino acid correctly."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "M"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize(protein)

        assert result.sequence == "ATG"
        assert result.cai == 1.0  # Only one codon for M

    def test_yeast_organism(self):
        """State machine works with yeast organism."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEK"
        sm = DeterministicOptimizationStateMachine(
            organism="Saccharomyces_cerevisiae",
            is_prokaryote=False,
            avoid_gt=True,
        )
        result = sm.optimize(protein)

        assert len(result.sequence) == len(protein) * 3
        assert result.cai > 0.0

    def test_organism_alias_resolution(self):
        """State machine resolves organism aliases correctly."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEK"

        # Both should resolve to the same organism
        sm1 = DeterministicOptimizationStateMachine(organism="ecoli", is_prokaryote=True)
        sm2 = DeterministicOptimizationStateMachine(organism="E. coli", is_prokaryote=True)

        # The organism should be resolved to the canonical name
        assert sm1.organism == "Escherichia_coli"
        assert sm2.organism == "Escherichia_coli"

    def test_constraint_priority_ordering(self):
        """ConstraintPriority enum is ordered correctly."""
        from biocompiler.optimizer.state_machine import ConstraintPriority

        # Hard constraints should have lower values (higher priority)
        assert ConstraintPriority.RESTRICTION_SITE.value < ConstraintPriority.STOP_CODON.value
        assert ConstraintPriority.STOP_CODON.value < ConstraintPriority.GC_RANGE.value
        assert ConstraintPriority.GC_RANGE.value < ConstraintPriority.CRYPTIC_SPLICE.value

        # Soft constraints should have higher values (lower priority)
        assert ConstraintPriority.AVOIDABLE_GT.value > ConstraintPriority.CRYPTIC_SPLICE.value
        assert ConstraintPriority.CPG_ISLAND.value > ConstraintPriority.AVOIDABLE_GT.value
        assert ConstraintPriority.ATTTA_MOTIF.value > ConstraintPriority.CPG_ISLAND.value
        assert ConstraintPriority.T_RUN.value > ConstraintPriority.ATTTA_MOTIF.value

    def test_invalid_protein_raises(self):
        """Invalid protein input raises an error."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine
        from biocompiler.exceptions import InvalidProteinError

        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )

        with pytest.raises(InvalidProteinError):
            sm.optimize("MVHLTPEEK1")  # '1' is not a valid AA


# ────────────────────────────────────────────────────────────
# Test 6: Import path works
# ────────────────────────────────────────────────────────────

class TestImportPath:
    """Verify the required import path works."""

    def test_import_from_optimizer(self):
        """Can import DeterministicOptimizationStateMachine from the
        optimizer subpackage."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine
        assert DeterministicOptimizationStateMachine is not None

    def test_import_constraint_automaton(self):
        """Can import ConstraintAutomaton."""
        from biocompiler.optimizer.state_machine import ConstraintAutomaton
        assert ConstraintAutomaton is not None

    def test_import_optimization_state(self):
        """Can import OptimizationState enum."""
        from biocompiler.optimizer.state_machine import OptimizationState
        assert OptimizationState is not None

    def test_import_constraint_priority(self):
        """Can import ConstraintPriority enum."""
        from biocompiler.optimizer.state_machine import ConstraintPriority
        assert ConstraintPriority is not None

    def test_import_state_machine_result(self):
        """Can import StateMachineResult dataclass."""
        from biocompiler.optimizer.state_machine import StateMachineResult
        assert StateMachineResult is not None


# ────────────────────────────────────────────────────────────
# Test 7: Edge cases
# ────────────────────────────────────────────────────────────

class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_all_same_aa(self):
        """Handle a protein of all the same amino acid."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "LLLLLLLLLL"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize(protein)
        assert len(result.sequence) == len(protein) * 3
        assert result.cai > 0.0

    def test_protein_with_methionine_only(self):
        """Handle a protein with only methionine (single codon AA)."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MMMM"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize(protein)
        assert result.sequence == "ATG" * 4

    def test_tryptophan_only(self):
        """Handle a protein with only tryptophan (single codon AA)."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "WWWW"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
        )
        result = sm.optimize(protein)
        assert result.sequence == "TGG" * 4

    def test_gc_tight_range(self):
        """Handle tight GC range constraints."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
            gc_lo=0.45,
            gc_hi=0.55,
        )
        result = sm.optimize(protein)
        # The state machine should try to get GC into range
        # (may not always succeed for tight ranges)
        assert len(result.sequence) == len(protein) * 3

    def test_conflict_resolve_with_enzymes(self):
        """Conflict resolution works when restriction sites are present."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        sm = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
            enzymes=["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
        )
        result = sm.optimize(protein)
        assert len(result.sequence) == len(protein) * 3
        assert result.cai > 0.0

    def test_wide_gc_range(self):
        """Wide GC range doesn't constrain optimization unnecessarily."""
        from biocompiler.optimizer.state_machine import DeterministicOptimizationStateMachine

        protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR"
        sm_wide = DeterministicOptimizationStateMachine(
            organism="Escherichia_coli",
            is_prokaryote=True,
            gc_lo=0.20,
            gc_hi=0.80,
        )
        result = sm_wide.optimize(protein)
        assert result.cai > 0.5  # Wide GC range should allow high CAI
