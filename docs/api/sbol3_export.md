# SBOL3 Export Documentation

BioCompiler supports export of optimized gene designs in **SBOL3** (Synthetic
Biology Open Language v3) format for interoperability with SBOL-compliant tools
such as Benchling, SynBioHub, and Pigeon.

## Overview

The SBOL3 export module (`biocompiler.export.sbol_export`) provides a **pure-Python**
SBOL3 XML/JSON-LD generator — no external SBOL library dependency is required.
The output conforms to the SBOL3 specification and can be validated with
`sbol-validate` or loaded into any SBOL3-compliant tool.

### What SBOL3 Export Includes

Each SBOL3 document contains:

- **Component** objects for the gene and CDS, with Sequence Ontology (SO)
  role annotations
- **Sequence** objects with the optimized DNA
- **Measure** objects for CAI and GC content (with Ontology of Units of
  Measure / OM units)
- **Activity** and **Plan** elements for provenance (PROV-O vocabulary)
- **Taxonomy** annotations linking to the target organism

### Supported Formats

| Format | Content Type | Extension | Description |
|--------|-------------|-----------|-------------|
| `sbol3` | RDF/XML | `.xml` | SBOL3 RDF/XML (default) |
| `sbol3json` | JSON-LD | `.json` | SBOL3 JSON-LD |

### Namespaces Used

| Prefix | URI |
|--------|-----|
| `sbol3` | `http://sbols.org/v3#` |
| `rdf` | `http://www.w3.org/1999/02/22-rdf-syntax-ns#` |
| `prov` | `http://www.w3.org/ns/prov#` |
| `so` | `http://sequenceontology.org/resource/SO:` |
| `om` | `http://www.ontology-of-units-of-measure.org/resource/om-2/` |
| `dct` | `http://purl.org/dc/terms/` |

### Sequence Ontology Roles

BioCompiler maps biological concepts to SO identifiers:

| BioCompiler Role | SO Identifier | Description |
|-----------------|---------------|-------------|
| `CDS` / `coding_sequence` | SO:0000316 | Coding Sequence |
| `promoter` | SO:0000167 | Promoter |
| `terminator` | SO:0000141 | Terminator |
| `RBS` / `ribosome_binding_site` | SO:0000139 | Ribosome Binding Site |
| `gene` | SO:0000704 | Gene |

## Python API

### Single Gene Export

```python
from biocompiler import optimize_sequence
from biocompiler.export.sbol_export import export_sbol

# Optimize a protein
result = optimize_sequence(
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
    organism="Escherichia_coli",
)

# Export as SBOL3 XML
path = export_sbol(
    optimization_result=result,
    output_path="gfp_optimized.xml",
    gene_name="gfp",
    organism="Escherichia_coli",
)
print(f"SBOL3 document saved to: {path}")

# Export as SBOL3 JSON-LD
path_json = export_sbol(
    optimization_result=result,
    output_path="gfp_optimized.json",
    format="sbol3json",
    gene_name="gfp",
    organism="Escherichia_coli",
)
```

### Collection Export (Multiple Genes)

```python
from biocompiler import optimize_sequence
from biocompiler.export.sbol_export import export_sbol_collection

# Optimize multiple proteins
results = []
for protein, name in [
    ("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE", "gfp"),
    ("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR", "hbb"),
]:
    result = optimize_sequence(protein, organism="Escherichia_coli")
    result.gene_name = name
    results.append(result)

# Export as SBOL3 Collection
path = export_sbol_collection(
    results=results,
    output_path="gene_library.xml",
    collection_name="e_coli_optimization_library",
    organism="Escherichia_coli",
)
```

### Custom Base URI

```python
# Use a custom base URI for SBOL identities
path = export_sbol(
    optimization_result=result,
    output_path="gfp.xml",
    base_uri="https://my-lab.org/sbol3",
    gene_name="gfp",
)
# Component identity will be: https://my-lab.org/sbol3/gfp
```

## REST API

Export via the REST API using the `/export/sbol3` endpoint:

```bash
# SBOL3 XML
curl -X POST http://localhost:8000/export/sbol3 \
  -H "X-API-Key: your-key" \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
    "organism": "Escherichia_coli",
    "gene_name": "gfp",
    "format": "sbol3"
  }'

# SBOL3 JSON-LD
curl -X POST http://localhost:8000/export/sbol3 \
  -H "X-API-Key: your-key" \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTG",
    "organism": "Escherichia_coli",
    "gene_name": "gfp",
    "format": "sbol3json"
  }'
```

## CLI

Export via the CLI `optimize` command with `--format sbol3`:

```bash
biocompiler optimize --input gfp.fasta \
    --organism Escherichia_coli \
    --format sbol3 \
    --output gfp_optimized.xml
```

## Example Input/Output

### Input

Protein sequence for eGFP, optimized for E. coli:

```python
from biocompiler import optimize_sequence

result = optimize_sequence(
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE",
    organism="Escherichia_coli",
)
```

### Output (SBOL3 XML — abbreviated)

```xml
<?xml version='1.0' encoding='utf-8'?>
<rdf:RDF>
  <sbol3:Component rdf:about="https://biocompiler.org/sbol3/gfp">
    <sbol3:displayId>gfp</sbol3:displayId>
    <sbol3:type rdf:resource="http://sequenceontology.org/resource/SO:0000347"/>
    <sbol3:role rdf:resource="http://sequenceontology.org/resource/SO:0000704"/>
    <dct:description>BioCompiler-optimized gene: gfp</dct:description>
    <sbol3:hasSequence rdf:resource="https://biocompiler.org/sbol3/gfp/sequence"/>
    <sbol3:hasMeasure>
      <sbol3:Measure rdf:about="https://biocompiler.org/sbol3/gfp/cai_measure">
        <sbol3:displayId>cai_measure</sbol3:displayId>
        <sbol3:value rdf:datatype="http://www.w3.org/2001/XMLSchema#decimal">0.999000</sbol3:value>
        <sbol3:unit rdf:resource="http://www.ontology-of-units-of-measure.org/resource/om-2/dimensionless"/>
        <dct:title>CAI</dct:title>
      </sbol3:Measure>
    </sbol3:hasMeasure>
    <sbol3:hasMeasure>
      <sbol3:Measure rdf:about="https://biocompiler.org/sbol3/gfp/gc_measure">
        <sbol3:displayId>gc_measure</sbol3:displayId>
        <sbol3:value rdf:datatype="http://www.w3.org/2001/XMLSchema#decimal">0.528000</sbol3:value>
        <sbol3:unit rdf:resource="http://www.ontology-of-units-of-measure.org/resource/om-2/fraction"/>
        <dct:title>GC_content</dct:title>
      </sbol3:Measure>
    </sbol3:hasMeasure>
    <sbol3:wasDerivedFrom rdf:resource="http://identifiers.org/taxonomy/Escherichia coli"/>
    <prov:wasGeneratedBy rdf:resource="https://biocompiler.org/sbol3/Activity/optimization_abc123"/>
    <dct:creator>BioCompiler v0.9.0</dct:creator>
  </sbol3:Component>

  <sbol3:Component rdf:about="https://biocompiler.org/sbol3/gfp_CDS">
    <sbol3:displayId>gfp_CDS</sbol3:displayId>
    <sbol3:type rdf:resource="http://sequenceontology.org/resource/SO:0000347"/>
    <sbol3:role rdf:resource="http://sequenceontology.org/resource/SO:0000316"/>
    <dct:description>Coding sequence for gfp (CAI=0.9990, GC=0.5280)</dct:description>
    <sbol3:hasSequence rdf:resource="https://biocompiler.org/sbol3/gfp_CDS/sequence"/>
    <dct:creator>BioCompiler v0.9.0</dct:creator>
  </sbol3:Component>

  <sbol3:Sequence rdf:about="https://biocompiler.org/sbol3/gfp/sequence">
    <sbol3:displayId>gfp_seq</sbol3:displayId>
    <sbol3:encoding rdf:resource="http://www.chem.qmul.ac.uk/iupac/DNA/"/>
    <sbol3:elements>ATGAGCAAAGGAGAACTGTTCACTGGAG...</sbol3:elements>
  </sbol3:Sequence>

  <prov:Activity rdf:about="https://biocompiler.org/sbol3/Activity/optimization_abc123">
    <sbol3:displayId>optimization_abc123</sbol3:displayId>
    <prov:startedAtTime>2026-03-04T12:00:00Z</prov:startedAtTime>
    <prov:hadPlan rdf:resource="https://biocompiler.org/sbol3/Plan/biocompiler_0.9.0"/>
    <prov:used rdf:resource="http://identifiers.org/taxonomy/Escherichia coli"/>
  </prov:Activity>

  <prov:Plan rdf:about="https://biocompiler.org/sbol3/Plan/biocompiler_0.9.0">
    <sbol3:displayId>biocompiler_0.9.0</sbol3:displayId>
    <dct:title>BioCompiler v0.9.0 codon optimization</dct:title>
  </prov:Plan>
</rdf:RDF>
```

### Output (SBOL3 JSON-LD — abbreviated)

```json
{
  "@context": {
    "sbol3": "http://sbols.org/v3#",
    "prov": "http://www.w3.org/ns/prov#",
    "so": "http://sequenceontology.org/resource/SO:",
    "om": "http://www.ontology-of-units-of-measure.org/resource/om-2/",
    "dct": "http://purl.org/dc/terms/",
    "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
  },
  "@id": "https://biocompiler.org/sbol3",
  "components": [
    {
      "@type": "sbol3:Component",
      "@id": "https://biocompiler.org/sbol3/gfp",
      "sbol3:displayId": "gfp",
      "sbol3:type": "http://sequenceontology.org/resource/SO:0000347",
      "sbol3:role": ["http://sequenceontology.org/resource/SO:0000704"],
      "dct:description": "BioCompiler-optimized gene: gfp",
      "sbol3:hasSequence": "https://biocompiler.org/sbol3/gfp/sequence",
      "sbol3:hasMeasure": [
        {
          "displayId": "cai_measure",
          "value": 0.999,
          "unit": "http://www.ontology-of-units-of-measure.org/resource/om-2/dimensionless",
          "title": "CAI"
        },
        {
          "displayId": "gc_measure",
          "value": 0.528,
          "unit": "http://www.ontology-of-units-of-measure.org/resource/om-2/fraction",
          "title": "GC_content"
        }
      ],
      "dct:creator": "BioCompiler v0.9.0"
    }
  ],
  "sequences": [
    {
      "@type": "sbol3:Sequence",
      "@id": "https://biocompiler.org/sbol3/gfp/sequence",
      "sbol3:displayId": "gfp_seq",
      "sbol3:elements": "ATGAGCAAAGGAGAACTGTTCACTGGAG...",
      "sbol3:encoding": "http://www.chem.qmul.ac.uk/iupac/DNA/"
    }
  ],
  "activities": [
    {
      "@type": "prov:Activity",
      "@id": "https://biocompiler.org/sbol3/Activity/optimization_abc123",
      "sbol3:displayId": "optimization_abc123",
      "prov:startedAtTime": "2026-03-04T12:00:00Z"
    }
  ]
}
```

## Validation

The SBOL3 output can be validated using:

```bash
# Using sbol-validate (if installed)
sbol-validate gfp_optimized.xml

# Using pySBOL3 (if installed)
python -c "
import sbol3
doc = sbol3.Document()
doc.read('gfp_optimized.xml')
report = doc.validate()
print(f'Valid: {len(report.errors) == 0}')
for err in report.errors:
    print(f'  Error: {err}')
"
```

## SBOLComponent Data Class

The internal `SBOLComponent` data class represents a single biological component:

```python
from biocompiler.export.sbol_export import SBOLComponent

comp = SBOLComponent(
    identity="",  # auto-assigned by export_sbol
    display_id="gfp",
    component_type="DNA",  # or "Protein"
    sequence="ATGAGCAAAGGAGAA...",
    roles=["gene"],  # maps to SO:0000704
    description="Optimized GFP gene",
)
```

## Limitations

- The REST endpoint SBOL3 output is minimal (Component + Sequence only), not the full-featured output from `sbol_export.py`. For complete SBOL3 documents with Measures, Activities, Plans, and taxonomy annotations, use the Python API directly.
- The SBOL3 export is a pure-Python generator. For advanced SBOL3 features
  (collections, combinatorial designs, attachments), consider using `pySBOL3`
  to post-process the output.
- The JSON-LD format example shown above is aspirational. The current implementation produces a simplified conversion. For production use, a
  proper RDF library should be used for full RDF serialization.
- Sub-component relationships (e.g., CDS within a gene) are represented as
  separate Component objects. Sequence-level composition constraints are
  not currently modeled.
