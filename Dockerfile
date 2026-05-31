# BioCompiler — Production Docker Image
# Multi-stage build for minimal image size

# ─── Build Stage ──────────────────────────────────────────────────
FROM python:3.12-slim AS builder

WORKDIR /build

# Install build dependencies
COPY pyproject.toml ./
COPY src/ src/

# Install package with all optional dependencies
RUN pip install --no-cache-dir --prefix=/install \
    -e ".[dev,api]" 2>/dev/null || \
    pip install --no-cache-dir --prefix=/install .

# ─── Runtime Stage ────────────────────────────────────────────────
FROM python:3.12-slim AS runtime

LABEL org.opencontainers.image.title="BioCompiler"
LABEL org.opencontainers.image.description="Machine-verified gene design REST API"
LABEL org.opencontainers.image.version="3.1.0"
LABEL org.opencontainers.image.source="https://github.com/pkhairkh/biocompiler"

# Create non-root user for security
RUN groupadd -r biocompiler && \
    useradd -r -g biocompiler -d /home/biocompiler -s /sbin/nologin biocompiler && \
    mkdir -p /home/biocompiler/.biocompiler && \
    chown -R biocompiler:biocompiler /home/biocompiler

# Copy installed packages from builder
COPY --from=builder /install /usr/local

# Copy application code
COPY --chown=biocompiler:biocompiler src/ /app/src/
COPY --chown=biocompiler:biocompiler pyproject.toml /app/

WORKDIR /app

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
