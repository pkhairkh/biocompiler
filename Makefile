# ============================================================================
# BioCompiler — Top-level Makefile for POPL Artifact Evaluation
# ============================================================================
#
# Standard targets for building, testing, and verifying the artifact.
# Run `make help` to see all available targets.
#
# Prerequisites:
#   - Python >= 3.10 with pip
#   - Lean4 / Lake (for proof targets)
#   - Docker + Docker Compose (for docker targets)
#   - BLAST+ tools (for blast target)
# ============================================================================

# Default target: show help
.DEFAULT_GOAL := help

# ---------------------------------------------------------------------------
# All targets are phony (no file outputs to track)
# ---------------------------------------------------------------------------
.PHONY: install test proof proof-check-sorry docker docker-up docker-down \
        blast benchmark clean all artifact help

# ---------------------------------------------------------------------------
# install — Install the package in editable mode with dev dependencies
# ---------------------------------------------------------------------------
install:
	pip install -e ".[dev]"

# ---------------------------------------------------------------------------
# test — Run the pytest suite with default options
#   (respects pyproject.toml [tool.pytest.ini_options])
# ---------------------------------------------------------------------------
test:
	pytest

# ---------------------------------------------------------------------------
# proof — Build Lean4 formal proofs via Lake
# ---------------------------------------------------------------------------
proof:
	cd proof && lake build

# ---------------------------------------------------------------------------
# proof-check-sorry — Verify no 'sorry' placeholders remain in proofs
#   Exits with error if any sorry is found (i.e. proof is incomplete).
# ---------------------------------------------------------------------------
proof-check-sorry:
	@echo "Checking for sorry in proof sources..."
	@! grep -r "sorry" proof/BioCompiler/ ; \
	if [ $$? -eq 0 ]; then \
		echo "ERROR: Found 'sorry' in proof sources — proofs are incomplete!"; \
		exit 1; \
	else \
		echo "OK: No 'sorry' found — all proofs are complete."; \
	fi

# ---------------------------------------------------------------------------
# docker — Build the Docker image
# ---------------------------------------------------------------------------
docker:
	docker compose build

# ---------------------------------------------------------------------------
# docker-up — Start the API server in detached mode
# ---------------------------------------------------------------------------
docker-up:
	docker compose up -d

# ---------------------------------------------------------------------------
# docker-down — Stop and remove containers
# ---------------------------------------------------------------------------
docker-down:
	docker compose down

# ---------------------------------------------------------------------------
# blast — Build BLAST databases from reference sequences
# ---------------------------------------------------------------------------
blast:
	python scripts/util/build_blast_databases.py

# ---------------------------------------------------------------------------
# benchmark — Run the benchmarking suite
# ---------------------------------------------------------------------------
benchmark:
	python -m biocompiler.benchmarking.runner

# ---------------------------------------------------------------------------
# clean — Remove build artifacts, caches, and generated files
# ---------------------------------------------------------------------------
clean:
	rm -rf build/ dist/ *.egg-info
	rm -rf src/biocompiler.egg-info
	rm -rf .pytest_cache/ .mypy_cache/ .ruff_cache/
	rm -rf htmlcov/ .coverage coverage.xml
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	cd proof && lake clean 2>/dev/null || true

# ---------------------------------------------------------------------------
# all — Install, test, and verify proofs
# ---------------------------------------------------------------------------
all: install test proof

# ---------------------------------------------------------------------------
# artifact — Full POPL artifact evaluation: install + test + proof + sorry check
# ---------------------------------------------------------------------------
artifact: install test proof proof-check-sorry

# ---------------------------------------------------------------------------
# help — List all available targets
# ---------------------------------------------------------------------------
help:
	@echo "BioCompiler — POPL Artifact Evaluation Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  install            Install package with dev dependencies (pip install -e .[dev])"
	@echo "  test               Run pytest suite"
	@echo "  proof              Build Lean4 formal proofs (lake build)"
	@echo "  proof-check-sorry  Verify no 'sorry' remains in proof sources"
	@echo "  docker             Build Docker image (docker compose build)"
	@echo "  docker-up          Start API server in detached mode"
	@echo "  docker-down        Stop and remove containers"
	@echo "  blast              Build BLAST databases from reference sequences"
	@echo "  benchmark          Run the benchmarking suite"
	@echo "  clean              Remove build artifacts and caches"
	@echo "  all                install + test + proof"
	@echo "  artifact           Full POPL evaluation: install + test + proof + sorry check"
	@echo "  help               Show this help message"
