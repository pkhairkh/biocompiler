"""BioCompiler Export Subpackage.

Provides SBOL, GenBank, and annotation export/import functionality.

Import order matters: ``sbol3_export`` (canonical SBOL3 implementation)
MUST come after ``sbol_legacy`` so its ``export_sbol3`` symbol wins the
wildcard re-export.  Previously the legacy module was imported last and
shadowed the canonical implementation with its incompatible signature
(see audit issue H24).
"""
from .core import *  # noqa: F401,F403
from biocompiler.export.sbol_export import *  # noqa: F401,F403
from biocompiler.export.sbol_import import *  # noqa: F401,F403
# DEPRECATED: sbol_legacy will be removed in v2.0; use sbol_export or sbol3_export.
# Imported BEFORE sbol3_export so that sbol3_export.export_sbol3 (canonical,
# compatible signature) wins the wildcard re-export — fixes audit issue H24.
from .sbol_legacy import *  # noqa: F401,F403
from biocompiler.export.sbol3_export import *  # noqa: F401,F403
from biocompiler.export.genbank_annotations import *  # noqa: F401,F403
from biocompiler.export.genbank_roundtrip import *  # noqa: F401,F403
from biocompiler.export.annotation import *  # noqa: F401,F403

__all__ = [  # noqa: F405
    "AnnotationDirective", "BC_NS", "BIOPAX_NS", "DCT_NS", "DNAREGION_TYPE",
    "FastaSequenceEntry", "GenBankAnnotationResult", "IUPAC_DNA_ENCODING",
    "OptimizationResult", "RDF_NS", "RestrictionSiteInfo", "RoundTripResult",
    "SBOL3_NS", "SBOLComponent", "SO_CDS", "SO_GENE", "SO_NS",
    "SO_PROMOTER", "SO_TERMINATOR", "SequenceAnnotation",
    "annotate_sequence", "annotate_to_genbank",
    "annotations_to_optimization_params", "compare_sequences",
    "export_batch_fasta", "export_fasta", "export_full_construct",
    "export_genbank", "export_genbank_with_certificate", "export_json",
    "export_multi_fasta", "export_sbol", "export_sbol3",
    "export_sbol_collection", "export_with_annotations",
    "format_biosecurity_report", "generate_sbol3_json",
    "generate_sbol3_xml", "import_sbol", "optimize_from_genbank",
    "parse_annotation_note", "parse_genbank_annotations",
    "sbol_to_genespecs", "verify_annotation_preservation",
    "verify_genbank_roundtrip",
]
