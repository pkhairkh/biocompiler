# BioCompiler — Production Docker Image
# Multi-stage build for minimal image size

# ─── Build Stage ──────────────────────────────────────────────────
FROM python:3.12-slim AS builder

WORKDIR /build

# Install build dependencies
COPY pyproject.toml ./
COPY src/ src/

# Install package with API optional dependencies
RUN pip install --no-cache-dir --prefix=/install \
    ".[api]" 2>/dev/null || \
    pip install --no-cache-dir --prefix=/install .

# ─── Runtime Stage ────────────────────────────────────────────────
FROM python:3.12-slim AS runtime

LABEL org.opencontainers.image.title="BioCompiler"
LABEL org.opencontainers.image.description="Machine-verified gene design REST API"
LABEL org.opencontainers.image.version="9.2.0"
LABEL org.opencontainers.image.source="https://github.com/pkhairkh/biocompiler"

# Install NCBI BLAST+ for homology-based biosecurity screening
RUN apt-get update && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for security
RUN groupadd -r biocompiler && \
    useradd -r -g biocompiler -d /home/biocompiler -s /sbin/nologin biocompiler && \
    mkdir -p /home/biocompiler/.biocompiler && \
    chown -R biocompiler:biocompiler /home/biocompiler

# Copy installed packages from builder
COPY --from=builder /install /usr/local

# Install z3-solver (via pip, included in dev deps)
RUN pip install --no-cache-dir z3-solver>=4.12.0 || true

# Copy application code
COPY --chown=biocompiler:biocompiler src/ /app/src/
COPY --chown=biocompiler:biocompiler scripts/ /app/scripts/
COPY --chown=biocompiler:biocompiler data/ /app/data/
COPY --chown=biocompiler:biocompiler pyproject.toml /app/

# Install Lean4/elan for proof verification
# Note: This must happen before USER directive since elan installs to /root
COPY proof/ /app/proof/
RUN curl -sSfL https://github.com/leanprover/elan/releases/latest/download/elan-x86_64-unknown-linux-gnu.tar.gz | tar xz \
    && ./elan-init -y --default-toolchain none \
    && rm elan-init \
    && export PATH="$HOME/.elan/bin:$PATH" \
    && cd /app/proof \
    && elan toolchain install leanprover/lean4:v4.30.0 \
    && lake build \
    && cd /app
ENV PATH="/root/.elan/bin:${PATH}"

WORKDIR /app

# Build BLAST databases during image build
RUN python -c "from biocompiler.blast_integration import is_blast_available; print('BLAST+ available:', is_blast_available())" || true
RUN if command -v makeblastdb &> /dev/null; then \
        python scripts/build_blast_databases.py --output-dir /opt/biocompiler/blast_db || true; \
    fi
ENV BIOCOMPILER_BLAST_DB_PATH=/opt/biocompiler/blast_db

# Environment variables with sensible defaults
ENV BIOCOMPILER_HOST=0.0.0.0
ENV BIOCOMPILER_PORT=8000
ENV BIOCOMPILER_API_KEY=""
ENV BIOCOMPILER_RATE_LIMIT=60
ENV BIOCOMPILER_CORS_ORIGINS=*
ENV BIOCOMPILER_DB_PATH=/home/biocompiler/.biocompiler/organisms.db
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Health check
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/health')" || exit 1

# Expose API port
EXPOSE 8000

# Run as non-root user
USER biocompiler

# Start the API server
CMD ["python", "-m", "biocompiler.cli", "serve", \
     "--host", "0.0.0.0", \
     "--port", "8000"]
