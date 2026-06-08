# REST API Usage Examples

This guide provides practical examples for using the BioCompiler REST API with
`curl` and Python `requests`.

## Getting Started

### Install and Start the Server

```bash
# Install BioCompiler
pip install -e .

# Start the server (with auto-generated API key)
biocompiler serve --port 8000

# Or start without authentication (local development only)
biocompiler serve --port 8000 --no-auth
```

On first startup, the server prints an auto-generated API key:

```
Generated API key: a1b2c3d4e5f6... (save this!)
Key saved to ~/.biocompiler/api_key for reuse across restarts.
```

Set a specific API key for production:

```bash
export BIOCOMPILER_API_KEY="your-secret-key-here"
biocompiler serve --port 8000
```

## Authentication

### Using the API Key

All requests must include the `X-API-Key` header:

```bash
# With curl
curl -X POST http://localhost:8000/optimize \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{"protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE"}'
```

```python
# With Python requests
import requests

API_KEY = "your-secret-key-here"
BASE_URL = "http://localhost:8000"

headers = {
    "X-API-Key": API_KEY,
    "Content-Type": "application/json",
}
```

### Auth Modes

| Mode | Env Variable | Behaviour |
|------|-------------|-----------|
| `required` (default) | `BIOCOMPILER_AUTH_MODE=required` | Unauthenticated → HTTP 401 |
| `optional` | `BIOCOMPILER_AUTH_MODE=optional` | Unauthenticated allowed with warning header |
| `disabled` | `BIOCOMPILER_API_KEY=disabled` | No auth required (dangerous) |

## Health Check

```bash
curl http://localhost:8000/health
```

```python
response = requests.get(f"{BASE_URL}/health", headers=headers)
print(response.json())
```

Expected response:

```json
{
  "status": "healthy",
  "version": "12.0.0",
  "timestamp": "2026-03-04T12:00:00Z",
  "auth_enabled": true,
  "rate_limit_rpm": 60
}
```

## Optimize a Sequence

### Basic Optimization

```bash
curl -X POST http://localhost:8000/optimize \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
    "organism": "Escherichia_coli"
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/optimize",
    headers=headers,
    json={
        "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
        "organism": "Escherichia_coli",
    },
)
result = response.json()
print(f"CAI: {result['cai']:.4f}")
print(f"GC: {result['gc_content']:.4f}")
print(f"Satisfied: {result['satisfied_predicates']}")
print(f"Failed: {result['failed_predicates']}")
```

Expected response (abbreviated):

```json
{
  "sequence": "ATGAGCAAAGGAGAACTGTTCACTGGAG...",
  "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
  "cai": 0.999,
  "gc_content": 0.5280,
  "satisfied_predicates": [
    "no_stop_codons",
    "no_cryptic_splice",
    "no_cpg_island",
    "no_restriction_site",
    "valid_coding_seq",
    "gc_in_range",
    "cai_above_threshold"
  ],
  "failed_predicates": [],
  "fallback_used": false,
  "provenance_id": "opt_abc123",
  "organism_domain": "prokaryote",
  "source_organism": null,
  "therapeutic": false,
  "self_protein": null
}
```

### Optimization for Human (Eukaryote)

```bash
curl -X POST http://localhost:8000/optimize \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
    "organism": "Homo_sapiens",
    "gc_lo": 0.40,
    "gc_hi": 0.60,
    "enzymes": ["EcoRI", "BamHI", "XhoI"]
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/optimize",
    headers=headers,
    json={
        "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
        "organism": "Homo_sapiens",
        "gc_lo": 0.40,
        "gc_hi": 0.60,
        "enzymes": ["EcoRI", "BamHI", "XhoI"],
    },
)
result = response.json()
```

### Therapeutic Protein with Source Organism

```bash
curl -X POST http://localhost:8000/optimize \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
    "organism": "Homo_sapiens",
    "source_organism": "Escherichia_coli",
    "therapeutic": true,
    "strict_mode": true
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/optimize",
    headers=headers,
    json={
        "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
        "organism": "Homo_sapiens",
        "source_organism": "Escherichia_coli",
        "therapeutic": True,
        "strict_mode": True,
    },
)
```

### Organism Short Keys

All organism parameters accept short keys, abbreviated binomials, and display names:

```python
# All of these resolve to Escherichia_coli:
organisms = ["ecoli", "E_coli", "E. coli", "Escherichia_coli"]

for org in organisms:
    response = requests.post(
        f"{BASE_URL}/optimize",
        headers=headers,
        json={"protein": "MSKGEELFTG", "organism": org},
    )
    print(f"{org}: CAI={response.json()['cai']:.4f}")
```

## Check a Sequence

Type-check a DNA sequence against all predicates:

```bash
curl -X POST http://localhost:8000/check \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
    "organism": "Homo_sapiens"
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/check",
    headers=headers,
    json={
        "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
        "organism": "Homo_sapiens",
    },
)
result = response.json()
print(f"Overall verdict: {result['overall_verdict']}")
print(f"GC content: {result['gc_content']:.4f}")
for r in result['results']:
    print(f"  {r['predicate']}: {r['verdict']}")
```

Expected response (abbreviated):

```json
{
  "sequence_length": 39,
  "gc_content": 0.6923,
  "protein": "MVSKGEELFTG",
  "results": [
    {"predicate": "no_stop_codons", "verdict": "PASS", "violation": "", "knowledge_gap": false},
    {"predicate": "gc_in_range", "verdict": "PASS", "violation": "", "knowledge_gap": false},
    {"predicate": "cai_above_threshold", "verdict": "PASS", "violation": "", "knowledge_gap": false}
  ],
  "overall_verdict": "PASS",
  "certificate": { ... }
}
```

## Scan for Issues

```bash
curl -X POST http://localhost:8000/scan \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
    "enzymes": ["EcoRI", "BamHI"],
    "find_orfs": true
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/scan",
    headers=headers,
    json={
        "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
        "enzymes": ["EcoRI", "BamHI"],
        "find_orfs": True,
    },
)
result = response.json()
print(f"Sequence length: {result['sequence_length']}")
print(f"Tokens: {result['tokens']}")
if result.get('orfs'):
    for orf in result['orfs']:
        print(f"  ORF: {orf}")
```

## Export to FASTA

```bash
curl -X POST http://localhost:8000/export/fasta \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
    "identifier": "GFP_optimized",
    "description": "eGFP optimized for E. coli",
    "organism": "Escherichia_coli"
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/export/fasta",
    headers=headers,
    json={
        "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
        "identifier": "GFP_optimized",
        "description": "eGFP optimized for E. coli",
        "organism": "Escherichia_coli",
    },
)
result = response.json()
print(result['content'])
```

Expected response:

```json
{
  "format": "fasta",
  "content": ">GFP_optimized|Escherichia_coli eGFP optimized for E. coli\nATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG\n"
}
```

## Export to GenBank

```bash
curl -X POST http://localhost:8000/export/genbank \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
    "locus_name": "GFP_OPT",
    "definition": "Optimized eGFP gene",
    "organism": "Escherichia_coli",
    "gene_name": "gfp"
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/export/genbank",
    headers=headers,
    json={
        "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
        "locus_name": "GFP_OPT",
        "definition": "Optimized eGFP gene",
        "organism": "Escherichia_coli",
        "gene_name": "gfp",
    },
)
result = response.json()
print(result['content'])
```

## Export to SBOL3

```bash
# XML format
curl -X POST http://localhost:8000/export/sbol3 \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
    "organism": "Escherichia_coli",
    "gene_name": "gfp",
    "format": "sbol3"
  }'

# JSON-LD format
curl -X POST http://localhost:8000/export/sbol3 \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
    "organism": "Escherichia_coli",
    "gene_name": "gfp",
    "format": "sbol3json"
  }'
```

```python
# SBOL3 XML
response = requests.post(
    f"{BASE_URL}/export/sbol3",
    headers=headers,
    json={
        "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
        "organism": "Escherichia_coli",
        "gene_name": "gfp",
        "format": "sbol3",
    },
)
sbol_xml = response.json()["content"]

# SBOL3 JSON-LD
response = requests.post(
    f"{BASE_URL}/export/sbol3",
    headers=headers,
    json={
        "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
        "organism": "Escherichia_coli",
        "gene_name": "gfp",
        "format": "sbol3json",
    },
)
sbol_json = response.json()["content"]
```

See [SBOL3 Export](sbol3_export.md) for detailed SBOL3 documentation.

## Batch Operations

### Batch Type-Check

Check up to 50 sequences in one request:

```bash
curl -X POST http://localhost:8000/batch/check \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "sequences": [
      {"sequence": "ATGGTGAGCAAGGGCGAGG", "organism": "Homo_sapiens"},
      {"sequence": "ATGAGCAAAGGAGAACTGTTC", "organism": "Escherichia_coli"}
    ]
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/batch/check",
    headers=headers,
    json={
        "sequences": [
            {"sequence": "ATGGTGAGCAAGGGCGAGG", "organism": "Homo_sapiens"},
            {"sequence": "ATGAGCAAAGGAGAACTGTTC", "organism": "Escherichia_coli"},
        ]
    },
)
result = response.json()
print(f"Total: {result['summary']['total']}")
print(f"Pass: {result['summary']['pass']}")
print(f"Fail: {result['summary']['fail']}")
```

### Batch Optimize

Optimize up to 20 proteins in one request:

```bash
curl -X POST http://localhost:8000/batch/optimize \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "proteins": [
      {
        "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
        "organism": "Escherichia_coli"
      },
      {
        "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
        "organism": "Homo_sapiens"
      }
    ]
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/batch/optimize",
    headers=headers,
    json={
        "proteins": [
            {
                "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
                "organism": "Escherichia_coli",
            },
            {
                "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
                "organism": "Homo_sapiens",
            },
        ]
    },
)
result = response.json()
print(f"Total: {result['summary']['total']}")
print(f"All satisfied: {result['summary']['all_satisfied']}")
for r in result['results']:
    print(f"  CAI={r['cai']:.4f} GC={r['gc_content']:.4f}")
```

### Fast Batch Optimize (Shared Parameters)

When all proteins target the same organism, use the simplified fast batch endpoint:

```bash
curl -X POST http://localhost:8000/batch/optimize \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "proteins": [
      "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
      "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR"
    ],
    "organism": "Escherichia_coli",
    "gc_lo": 0.45,
    "gc_hi": 0.55
  }'
```

### Batch Export

Export up to 50 sequences in mixed formats:

```bash
curl -X POST http://localhost:8000/batch/export \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "sequences": [
      {
        "sequence": "ATGGTGAGCAAGGGCGAGG",
        "format": "fasta",
        "identifier": "seq1"
      },
      {
        "sequence": "ATGAGCAAAGGAGAACTGTTC",
        "format": "genbank",
        "locus_name": "SEQ2"
      }
    ]
  }'
```

## Protein Analysis Endpoints

### Structure Prediction

```bash
curl -X POST http://localhost:8000/protein/structure/predict \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
    "organism": "Escherichia_coli"
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/protein/structure/predict",
    headers=headers,
    json={
        "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
        "organism": "Escherichia_coli",
    },
)
result = response.json()
print(f"Mean pLDDT: {result['mean_plddt']:.1f}")
print(f"Quality: {result['quality_class']}")
```

### Stability Analysis

```bash
curl -X POST http://localhost:8000/protein/stability/analyze \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
    "organism": "Escherichia_coli"
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/protein/stability/analyze",
    headers=headers,
    json={
        "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
        "organism": "Escherichia_coli",
    },
)
result = response.json()
```

### Solubility Analysis

```bash
curl -X POST http://localhost:8000/protein/solubility/analyze \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{"protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE"}'
```

```python
response = requests.post(
    f"{BASE_URL}/protein/solubility/analyze",
    headers=headers,
    json={"protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE"},
)
```

### Immunogenicity Analysis

```bash
curl -X POST http://localhost:8000/protein/immunogenicity/analyze \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
    "organism": "Homo_sapiens",
    "source_organism": "Escherichia_coli",
    "therapeutic": true
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/protein/immunogenicity/analyze",
    headers=headers,
    json={
        "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
        "organism": "Homo_sapiens",
        "source_organism": "Escherichia_coli",
        "therapeutic": True,
    },
)
```

### Deimmunization

```bash
curl -X POST http://localhost:8000/protein/immunogenicity/deimmunize \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
    "organism": "Homo_sapiens",
    "target_score": 0.2,
    "max_mutations": 5,
    "blosum62_min": 1
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/protein/immunogenicity/deimmunize",
    headers=headers,
    json={
        "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
        "organism": "Homo_sapiens",
        "target_score": 0.2,
        "max_mutations": 5,
        "blosum62_min": 1,
    },
)
```

### Full Protein Assessment

```bash
curl -X POST http://localhost:8000/protein/assessment/full \
  -H "X-API-Key: your-secret-key-here" \
  -H "Content-Type: application/json" \
  -d '{
    "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
    "organism": "Escherichia_coli",
    "run_structure": true,
    "run_stability": true,
    "run_solubility": true,
    "run_immunogenicity": true
  }'
```

```python
response = requests.post(
    f"{BASE_URL}/protein/assessment/full",
    headers=headers,
    json={
        "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
        "organism": "Escherichia_coli",
        "run_structure": True,
        "run_stability": True,
        "run_solubility": True,
        "run_immunogenicity": True,
    },
)
```

## Information Endpoints

### List Supported Organisms

```bash
curl -H "X-API-Key: your-secret-key-here" http://localhost:8000/organisms
```

```python
response = requests.get(f"{BASE_URL}/organisms", headers=headers)
for org in response.json()["organisms"]:
    print(f"  {org}")
```

### List Predicates

```bash
curl -H "X-API-Key: your-secret-key-here" http://localhost:8000/predicates
```

### List Enzymes

```bash
curl -H "X-API-Key: your-secret-key-here" http://localhost:8000/enzymes
```

### Server Info

```bash
curl -H "X-API-Key: your-secret-key-here" http://localhost:8000/info
```

## Provenance Endpoints

### List Provenance Records

```bash
curl -H "X-API-Key: your-secret-key-here" http://localhost:8000/provenance
```

### Get Specific Provenance Record

```bash
curl -H "X-API-Key: your-secret-key-here" http://localhost:8000/provenance/opt_abc123
```

## Error Handling

The API returns standard HTTP status codes:

| Code | Meaning |
|------|---------|
| 200 | Success |
| 400 | Bad request (validation error in input) |
| 401 | Missing or invalid API key |
| 404 | Resource not found |
| 422 | Unprocessable entity (e.g., strict mode optimization failed) |
| 429 | Rate limit exceeded |

Example error response:

```json
{
  "detail": "Invalid or missing API key. Set X-API-Key header."
}
```

```python
response = requests.post(f"{BASE_URL}/optimize", headers=headers, json={...})
if response.status_code == 200:
    result = response.json()
elif response.status_code == 422:
    error = response.json()
    print(f"Optimization failed: {error['detail']}")
elif response.status_code == 429:
    print("Rate limit exceeded. Retry later.")
elif response.status_code == 401:
    print("Invalid API key.")
```

## Rate Limiting

- **Default**: 60 requests per minute per client.
- **Batch**: Each item counts as one request (e.g., 20 items = 20 units).
- **Configuration**: Set `BIOCOMPILER_RATE_LIMIT` environment variable.

If rate-limited, the API returns HTTP 429:

```json
{
  "detail": "Rate limit exceeded: 60 requests/minute. Retry after 60 seconds."
}
```
