"""
BioCompiler Benchmarking Registry
===================================
Central registry of standard benchmark proteins, target organisms, and
benchmark metrics for tracking optimization performance across releases.

This module provides a single source of truth for:
  - Standard benchmark proteins (EGFP, insulin, BSA, etc.)
  - Target organisms (E. coli, yeast, human, mouse, CHO)
  - Benchmark metrics (CAI, GC%, optimization time, constraint satisfaction rate)
  - Benchmark suite definitions (combinations of proteins x organisms)

Usage::

    from biocompiler.benchmarking.registry import (
        STANDARD_PROTEINS,
        TARGET_ORGANISMS,
        BENCHMARK_METRICS,
        BENCHMARK_SUITES,
        get_suite,
    )

    # Get all proteins for the "core" benchmark suite
    suite = get_suite("core")
    for protein_name, protein_data in suite.proteins.items():
        for organism in suite.organisms:
            ...
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


__all__ = [
    # Protein data
    "StandardProtein",
    "STANDARD_PROTEINS",
    # Organism data
    "TargetOrganism",
    "TARGET_ORGANISMS",
    # Metrics
    "BenchmarkMetric",
    "BENCHMARK_METRICS",
    # Suites
    "BenchmarkSuite",
    "BENCHMARK_SUITES",
    "get_suite",
    "get_all_suites",
]


# ---------------------------------------------------------------------------
# StandardProtein — definition of a benchmark protein
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class StandardProtein:
    """A standard benchmark protein definition.

    Attributes
    ----------
    name : str
        Gene symbol (e.g. ``"EGFP"``, ``"INS"``).
    protein_sequence : str
        Amino acid sequence (one-letter codes).
    description : str
        Human-readable protein name / description.
    uniprot_id : str
        UniProtKB accession, or empty string if synthetic.
    source_organism : str
        Natural source organism of the protein.
    length_category : str
        Size category: ``"small"`` (<150 aa), ``"medium"`` (150-400 aa),
        ``"large"`` (>400 aa), or ``"synthetic"``.
    expected_cai_range : tuple[float, float] | None
        Rough expected CAI range after optimization, or ``None`` if unknown.
    category : str
        Functional category (e.g. ``"reporter"``, ``"therapeutic"``,
        ``"structural"``, ``"enzyme"``, ``"vaccine_antigen"``).
    """

    name: str
    protein_sequence: str
    description: str = ""
    uniprot_id: str = ""
    source_organism: str = ""
    length_category: str = "medium"
    expected_cai_range: tuple[float, float] | None = None
    category: str = "general"

    @property
    def length(self) -> int:
        """Protein length in amino acids."""
        return len(self.protein_sequence)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "name": self.name,
            "protein_sequence": self.protein_sequence,
            "description": self.description,
            "uniprot_id": self.uniprot_id,
            "source_organism": self.source_organism,
            "length_category": self.length_category,
            "length": self.length,
            "expected_cai_range": (
                list(self.expected_cai_range) if self.expected_cai_range else None
            ),
            "category": self.category,
        }


# ---------------------------------------------------------------------------
# TargetOrganism — definition of a target organism for optimization
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TargetOrganism:
    """A target organism for codon optimization benchmarking.

    Attributes
    ----------
    name : str
        Normalised organism name matching BioCompiler internals
        (e.g. ``"Escherichia_coli"``, ``"Homo_sapiens"``).
    display_name : str
        Human-readable name (e.g. ``"E. coli"``, ``"Human"``).
    domain : str
        Taxonomic domain: ``"prokaryote"`` or ``"eukaryote"``.
    typical_gc_range : tuple[float, float]
        Typical GC content range for this organism's genome.
    codon_bias_strength : str
        Qualitative codon bias strength: ``"strong"``, ``"moderate"``, ``"weak"``.
    industrial_relevance : str
        Relevance for industrial applications:
        ``"high"``, ``"medium"``, or ``"low"``.
    """

    name: str
    display_name: str = ""
    domain: str = "eukaryote"
    typical_gc_range: tuple[float, float] = (0.35, 0.65)
    codon_bias_strength: str = "moderate"
    industrial_relevance: str = "medium"

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "name": self.name,
            "display_name": self.display_name,
            "domain": self.domain,
            "typical_gc_range": list(self.typical_gc_range),
            "codon_bias_strength": self.codon_bias_strength,
            "industrial_relevance": self.industrial_relevance,
        }


# ---------------------------------------------------------------------------
# BenchmarkMetric — definition of a tracked benchmark metric
# ---------------------------------------------------------------------------

class MetricDirection(str, Enum):
    """Direction of desirability for a metric."""
    HIGHER_IS_BETTER = "higher_is_better"
    LOWER_IS_BETTER = "lower_is_better"
    CLOSER_TO_TARGET = "closer_to_target"


@dataclass(frozen=True)
class BenchmarkMetric:
    """Definition of a benchmark metric to track across releases.

    Attributes
    ----------
    name : str
        Machine-readable metric name (e.g. ``"cai"``, ``"gc_mean"``).
    display_name : str
        Human-readable name (e.g. ``"CAI"``, ``"Mean GC%"``).
    description : str
        What this metric measures.
    direction : MetricDirection
        Whether higher, lower, or closer-to-target values are better.
    unit : str
        Unit of measurement (e.g. ``""``, ``"%"``, ``"s"``, ``"count"``).
    regression_threshold : float
        Threshold for flagging a regression (as a fraction).
        E.g. 0.10 means a 10% change triggers a regression alert.
    target_value : float | None
        Target value for ``CLOSER_TO_TARGET`` metrics, or ``None``.
    """

    name: str
    display_name: str = ""
    description: str = ""
    direction: MetricDirection = MetricDirection.HIGHER_IS_BETTER
    unit: str = ""
    regression_threshold: float = 0.10
    target_value: float | None = None

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "name": self.name,
            "display_name": self.display_name,
            "description": self.description,
            "direction": self.direction.value,
            "unit": self.unit,
            "regression_threshold": self.regression_threshold,
            "target_value": self.target_value,
        }


# ---------------------------------------------------------------------------
# BenchmarkSuite — a named collection of proteins x organisms
# ---------------------------------------------------------------------------

@dataclass
class BenchmarkSuite:
    """A named benchmark suite combining a set of proteins and target organisms.

    Attributes
    ----------
    name : str
        Suite name (e.g. ``"core"``, ``"full"``, ``"therapeutic"``).
    description : str
        What this suite is designed to test.
    proteins : dict[str, StandardProtein]
        Proteins to benchmark, keyed by gene name.
    organisms : list[TargetOrganism]
        Target organisms for optimization.
    metrics : list[BenchmarkMetric]
        Metrics to collect for each (protein, organism) pair.
    """

    name: str
    description: str = ""
    proteins: dict[str, StandardProtein] = field(default_factory=dict)
    organisms: list[TargetOrganism] = field(default_factory=list)
    metrics: list[BenchmarkMetric] = field(default_factory=list)

    @property
    def num_benchmarks(self) -> int:
        """Total number of (protein, organism) benchmark combinations."""
        return len(self.proteins) * len(self.organisms)

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a JSON-compatible dict."""
        return {
            "name": self.name,
            "description": self.description,
            "num_proteins": len(self.proteins),
            "num_organisms": len(self.organisms),
            "num_benchmarks": self.num_benchmarks,
            "proteins": {k: v.to_dict() for k, v in self.proteins.items()},
            "organisms": [o.to_dict() for o in self.organisms],
            "metrics": [m.to_dict() for m in self.metrics],
        }


# ===========================================================================
# STANDARD PROTEINS
# ===========================================================================

STANDARD_PROTEINS: dict[str, StandardProtein] = {
    # ── Reporters ──
    "EGFP": StandardProtein(
        name="EGFP",
        protein_sequence=(
            "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
            "VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
            "RIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQ"
            "QNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        ),
        description="Enhanced green fluorescent protein; standard reporter",
        uniprot_id="C5MKY7",
        source_organism="Aequorea victoria",
        length_category="medium",
        expected_cai_range=(0.70, 0.95),
        category="reporter",
    ),
    "mCherry": StandardProtein(
        name="mCherry",
        protein_sequence=(
            "MVSKGEEDNMAIIIKFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFA"
            "WDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFI"
            "YKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTY"
            "KAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK"
        ),
        description="mCherry red fluorescent protein; alternative reporter",
        uniprot_id="X5DSL3",
        source_organism="Discosoma sp.",
        length_category="medium",
        expected_cai_range=(0.65, 0.90),
        category="reporter",
    ),

    # ── Therapeutic proteins ──
    "INS": StandardProtein(
        name="INS",
        protein_sequence=(
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQV"
            "GQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        ),
        description="Human insulin precursor; therapeutic: recombinant insulin",
        uniprot_id="P01308",
        source_organism="Homo sapiens",
        length_category="small",
        expected_cai_range=(0.55, 0.85),
        category="therapeutic",
    ),
    "GH1": StandardProtein(
        name="GH1",
        protein_sequence=(
            "MATGSRTSLLLAFGLLCLPWLQEGSAFPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPK"
            "EQKYSFLQNPQTSLCFSESIPTPSNREETQQKSNLELLRISLLLIQSWLEPVQFLRSVFANSLV"
            "YGASDSNVYDLLKDLEEGIQTLMGRLEDGSPRTGQIFKQTYSKFDTNSHNDDALLKNYGLLYCFR"
            "KDMDKVETFLRIVQCRSVEGSCGF"
        ),
        description="Somatotropin (growth hormone); therapeutic: somatropin",
        uniprot_id="P01241",
        source_organism="Homo sapiens",
        length_category="medium",
        expected_cai_range=(0.50, 0.80),
        category="therapeutic",
    ),
    "EPO": StandardProtein(
        name="EPO",
        protein_sequence=(
            "MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNE"
            "NITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSG"
            "LRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGD"
            "R"
        ),
        description="Erythropoietin precursor; therapeutic: epoetin alfa",
        uniprot_id="P01588",
        source_organism="Homo sapiens",
        length_category="medium",
        expected_cai_range=(0.50, 0.80),
        category="therapeutic",
    ),

    # ── Structural proteins ──
    "BSA": StandardProtein(
        name="BSA",
        protein_sequence=(
            "MKWVTFISLLLLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHV"
            "KLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLS"
            "HKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLYYANKYNGVFQEC"
            "CQEDGACLLKNIENPSSSSELSKAGQTALPNLKQAEFAEVSKLVTDLTKTVLTKVHTECCHGDL"
            "LECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDV"
            "CKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHL"
            "VDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERM"
            "PCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHAD"
            "ICTLSEKERQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVST"
            "QTALA"
        ),
        description="Bovine serum albumin; structural protein standard",
        uniprot_id="P02769",
        source_organism="Bos taurus",
        length_category="large",
        expected_cai_range=(0.45, 0.75),
        category="structural",
    ),
    "HBB": StandardProtein(
        name="HBB",
        protein_sequence=(
            "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSQVKGHGKKVAD"
            "ALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFL"
            "ASVSTVLTSKYR"
        ),
        description="Human hemoglobin beta chain; well-characterized self-protein",
        uniprot_id="P68871",
        source_organism="Homo sapiens",
        length_category="small",
        expected_cai_range=(0.55, 0.85),
        category="structural",
    ),

    # ── Enzymes ──
    "LacZ": StandardProtein(
        name="LacZ",
        protein_sequence=(
            "MTMITDSLAVVLQRRDWENPGVTQLNRLAAHPPFASWRNSEEARTDRPSQQLRSLNGEWRFAWFP"
            "APEAVPESWLECDLPEADTVVVPSNWQMHGYDAPIYTNVTYPITVNPPFVPTENPTGCYSLTFN"
            "VDESWLQEGQTRIIFDGVNSAFHLWCNGRWVGYGQDSRLPSEFDLSAFLRAGENRLAVMVLRWS"
            "DGSYLEDQDMWRMSGIFRDVSLLHKPTTQISDFHVATRFNDDFSRAVLEAEVQMCGELRDYLRV"
            "TVSLWQGETQVASGTAPFGGEIIDERGGYADRVTLRLNVENPKLWSAEIPNLYRAVVELHTADGT"
            "LIEAEACDVGFREVRIENGLLLLNGKPLLIRGVNRHEHHPLHGQVMDEQTMVQDILLMKQNNFNA"
            "VRCSHYPNHPLWYTLCDRYGLYVVDEANIETHGMVPMNRLTDDPRWLPAMSERVTRMVQDRKRYT"
            "LEEFVTAEWDESRPAVMEHVKLAPPPHPAEQLVDATLAEVLDRHTTLRGLQAYNAE"
        ),
        description="E. coli beta-galactosidase; classic reporter enzyme",
        uniprot_id="P00722",
        source_organism="Escherichia coli",
        length_category="large",
        expected_cai_range=(0.25, 0.55),
        category="enzyme",
    ),
    "T4_lysozyme": StandardProtein(
        name="T4_lysozyme",
        protein_sequence=(
            "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKAIGRNTNGVITKDEAE"
            "KLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEA"
            "AVNLQSRNFNRPQVNITKDGTSGDSYGAIQFNLRNTDNRIAISTLNGTKAQALKEKYF"
        ),
        description="T4 lysozyme; well-studied model enzyme for protein engineering",
        uniprot_id="P00720",
        source_organism="Enterobacteria phage T4",
        length_category="medium",
        expected_cai_range=(0.40, 0.70),
        category="enzyme",
    ),

    # ── Vaccine antigens ──
    "SARS2_RBD": StandardProtein(
        name="SARS2_RBD",
        protein_sequence=(
            "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSP"
            "TKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYN"
            "YLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFE"
            "LLHAPATVCGPKKSTNLVKNKCVNF"
        ),
        description="SARS-CoV-2 Spike RBD; vaccine target",
        uniprot_id="P0DTC2",
        source_organism="SARS-CoV-2",
        length_category="medium",
        expected_cai_range=(0.50, 0.80),
        category="vaccine_antigen",
    ),
}


# ===========================================================================
# TARGET ORGANISMS
# ===========================================================================

TARGET_ORGANISMS: dict[str, TargetOrganism] = {
    "Escherichia_coli": TargetOrganism(
        name="Escherichia_coli",
        display_name="E. coli",
        domain="prokaryote",
        typical_gc_range=(0.48, 0.55),
        codon_bias_strength="strong",
        industrial_relevance="high",
    ),
    "Saccharomyces_cerevisiae": TargetOrganism(
        name="Saccharomyces_cerevisiae",
        display_name="Yeast (S. cerevisiae)",
        domain="eukaryote",
        typical_gc_range=(0.36, 0.42),
        codon_bias_strength="strong",
        industrial_relevance="high",
    ),
    "Homo_sapiens": TargetOrganism(
        name="Homo_sapiens",
        display_name="Human",
        domain="eukaryote",
        typical_gc_range=(0.38, 0.44),
        codon_bias_strength="moderate",
        industrial_relevance="high",
    ),
    "Mus_musculus": TargetOrganism(
        name="Mus_musculus",
        display_name="Mouse",
        domain="eukaryote",
        typical_gc_range=(0.38, 0.44),
        codon_bias_strength="moderate",
        industrial_relevance="high",
    ),
    "CHO_K1": TargetOrganism(
        name="CHO_K1",
        display_name="CHO (Cricetulus griseus)",
        domain="eukaryote",
        typical_gc_range=(0.38, 0.44),
        codon_bias_strength="moderate",
        industrial_relevance="high",
    ),
}


# ===========================================================================
# BENCHMARK METRICS
# ===========================================================================

BENCHMARK_METRICS: dict[str, BenchmarkMetric] = {
    "cai": BenchmarkMetric(
        name="cai",
        display_name="CAI",
        description="Codon Adaptation Index (Sharp & Li 1987); measures codon usage bias",
        direction=MetricDirection.HIGHER_IS_BETTER,
        unit="",
        regression_threshold=0.05,  # 5% CAI drop = regression
    ),
    "gc_mean": BenchmarkMetric(
        name="gc_mean",
        display_name="Mean GC%",
        description="Mean GC content across the coding sequence",
        direction=MetricDirection.CLOSER_TO_TARGET,
        unit="%",
        regression_threshold=0.10,
        target_value=0.50,
    ),
    "gc_std": BenchmarkMetric(
        name="gc_std",
        display_name="GC Std Dev",
        description="Standard deviation of GC content across sliding windows",
        direction=MetricDirection.LOWER_IS_BETTER,
        unit="",
        regression_threshold=0.15,
    ),
    "optimization_time": BenchmarkMetric(
        name="optimization_time",
        display_name="Optimization Time",
        description="Wall-clock time for the full optimization pipeline",
        direction=MetricDirection.LOWER_IS_BETTER,
        unit="s",
        regression_threshold=0.10,  # 10% slowdown = regression
    ),
    "constraint_satisfaction_rate": BenchmarkMetric(
        name="constraint_satisfaction_rate",
        display_name="Constraint Satisfaction Rate",
        description="Fraction of all constraints satisfied by the optimized sequence",
        direction=MetricDirection.HIGHER_IS_BETTER,
        unit="",
        regression_threshold=0.05,
    ),
    "restriction_site_total": BenchmarkMetric(
        name="restriction_site_total",
        display_name="Restriction Sites",
        description="Total number of restriction enzyme recognition sites found",
        direction=MetricDirection.LOWER_IS_BETTER,
        unit="count",
        regression_threshold=0.20,
    ),
    "cryptic_splice_sites": BenchmarkMetric(
        name="cryptic_splice_sites",
        display_name="Cryptic Splice Sites",
        description="Number of cryptic GT/AG splice-site dinucleotides",
        direction=MetricDirection.LOWER_IS_BETTER,
        unit="count",
        regression_threshold=0.20,
    ),
    "cpg_islands": BenchmarkMetric(
        name="cpg_islands",
        display_name="CpG Islands",
        description="Number of CpG islands detected in the sequence",
        direction=MetricDirection.LOWER_IS_BETTER,
        unit="count",
        regression_threshold=0.20,
    ),
    "mrna_stability": BenchmarkMetric(
        name="mrna_stability",
        display_name="mRNA Stability",
        description="Composite mRNA stability score (T-run + hairpin potential)",
        direction=MetricDirection.HIGHER_IS_BETTER,
        unit="",
        regression_threshold=0.10,
    ),
}


# ===========================================================================
# BENCHMARK SUITES
# ===========================================================================

def _build_core_suite() -> BenchmarkSuite:
    """Build the core benchmark suite — fast, representative cross-section."""
    proteins = {k: STANDARD_PROTEINS[k] for k in [
        "EGFP", "INS", "BSA", "HBB",
    ] if k in STANDARD_PROTEINS}
    organisms = [TARGET_ORGANISMS[k] for k in [
        "Escherichia_coli", "Homo_sapiens",
    ] if k in TARGET_ORGANISMS]
    metrics = [BENCHMARK_METRICS[k] for k in [
        "cai", "gc_mean", "optimization_time", "constraint_satisfaction_rate",
    ] if k in BENCHMARK_METRICS]
    return BenchmarkSuite(
        name="core",
        description=(
            "Core benchmark suite: 4 representative proteins x 2 organisms. "
            "Designed to run in < 5 minutes and catch most regressions."
        ),
        proteins=proteins,
        organisms=organisms,
        metrics=metrics,
    )


def _build_full_suite() -> BenchmarkSuite:
    """Build the full benchmark suite — all proteins x all organisms."""
    metrics = list(BENCHMARK_METRICS.values())
    return BenchmarkSuite(
        name="full",
        description=(
            "Full benchmark suite: all standard proteins x all target organisms "
            "with complete metric tracking. Designed for release validation."
        ),
        proteins=dict(STANDARD_PROTEINS),
        organisms=list(TARGET_ORGANISMS.values()),
        metrics=metrics,
    )


def _build_therapeutic_suite() -> BenchmarkSuite:
    """Build the therapeutic protein benchmark suite."""
    proteins = {k: v for k, v in STANDARD_PROTEINS.items()
                if v.category == "therapeutic"}
    organisms = [TARGET_ORGANISMS[k] for k in [
        "Homo_sapiens", "CHO_K1",
    ] if k in TARGET_ORGANISMS]
    metrics = [BENCHMARK_METRICS[k] for k in [
        "cai", "gc_mean", "optimization_time",
        "constraint_satisfaction_rate", "cryptic_splice_sites", "cpg_islands",
    ] if k in BENCHMARK_METRICS]
    return BenchmarkSuite(
        name="therapeutic",
        description=(
            "Therapeutic protein benchmark: biopharma-relevant proteins "
            "optimized for mammalian expression (human + CHO)."
        ),
        proteins=proteins,
        organisms=organisms,
        metrics=metrics,
    )


def _build_cross_organism_suite() -> BenchmarkSuite:
    """Build the cross-organism benchmark — single protein, all organisms."""
    proteins = {}
    if "EGFP" in STANDARD_PROTEINS:
        proteins["EGFP"] = STANDARD_PROTEINS["EGFP"]
    metrics = [BENCHMARK_METRICS[k] for k in [
        "cai", "gc_mean", "optimization_time", "constraint_satisfaction_rate",
    ] if k in BENCHMARK_METRICS]
    return BenchmarkSuite(
        name="cross_organism",
        description=(
            "Cross-organism benchmark: EGFP optimized for all 5 target "
            "organisms. Tests organism-specific adaptation quality."
        ),
        proteins=proteins,
        organisms=list(TARGET_ORGANISMS.values()),
        metrics=metrics,
    )


BENCHMARK_SUITES: dict[str, BenchmarkSuite] = {
    "core": _build_core_suite(),
    "full": _build_full_suite(),
    "therapeutic": _build_therapeutic_suite(),
    "cross_organism": _build_cross_organism_suite(),
}
"""Registry of named benchmark suites.

Each suite specifies which proteins, organisms, and metrics to use.
"""


def get_suite(name: str) -> BenchmarkSuite:
    """Look up a benchmark suite by name.

    Args:
        name: Suite name (e.g. ``"core"``, ``"full"``).

    Returns:
        BenchmarkSuite instance.

    Raises:
        ValueError: If the suite name is not recognised.
    """
    if name not in BENCHMARK_SUITES:
        available = ", ".join(sorted(BENCHMARK_SUITES.keys()))
        raise ValueError(
            f"Unknown benchmark suite '{name}'. Available: {available}"
        )
    return BENCHMARK_SUITES[name]


def get_all_suites() -> dict[str, BenchmarkSuite]:
    """Return a copy of all registered benchmark suites."""
    return dict(BENCHMARK_SUITES)
