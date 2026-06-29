"""Immunogenicity scoring, MHC binding prediction, and B-cell epitope prediction.

.. warning::

    HONESTY NOTE — FABRICATED PSSM SCORES
    ------------------------------------
    The default immunogenicity PSSMs in this module are APPROXIMATE, not real
    binding data.  The PSSM score values were assigned by hand based on
    published anchor residue preferences (see ``_build_mhc_i_pssms`` and
    ``_build_mhc_ii_pssms`` below — the comments literally say
    "guessed/approximate scores, NOT scores derived from experimental binding
    data").

    Consequently, :func:`compute_immunogenicity` returns an
    :class:`ImmunogenicityResult` whose ``verdict`` field is
    :attr:`~biocompiler.shared.types.Verdict.UNCERTAIN` by default, with
    ``reason="fabricated_scores"`` and ``data_source="guessed_pssm"``.
    The numerical scores are still computed and exposed on the result for
    reference, but they MUST NOT be used to issue a PASS/FAIL verdict.

    To obtain a real PASS/FAIL verdict you must:

    1. Install NetMHCpan (requires DTU license) — a real binding-affinity
       predictor trained on >800,000 IEDB measurements.
    2. Install MHCflurry with downloaded models — a neural-network predictor
       trained on IEDB data.
    3. Call ``compute_immunogenicity(protein, use_real_data=True)``.

    The wrapper predicates in
    :mod:`biocompiler.type_system.predicates` (``check_low_immunogenicity``
    and friends) honour this contract: they return UNCERTAIN whenever the
    underlying ``ImmunogenicityResult.verdict`` is UNCERTAIN, and never
    return PASS on fabricated data.

    This is an honesty measure — we do not return PASS based on fabricated
    scores.

This module consolidates three formerly separate modules:

* **MHC binding** (formerly ``mhc_binding.py``): Predicts peptide-MHC
  binding affinity using position-specific scoring matrices (PSSMs) derived
  from known binding motifs in the Immune Epitope Database (IEDB).
  When ``use_netmhcpan=True``, the NetMHCpan web API is tried first
  for more accurate predictions; the PSSM heuristic serves as fallback.
* **B-cell epitope prediction** (formerly ``epitope.py``): Linear and
  conformational B-cell epitope prediction using multiple classical scales
  and methods (Kolaskar-Tongaonkar, Parker hydrophilicity, Chou-Fasman
  beta-turn, Emini surface accessibility, BepiPred-like composite, and
  conformational epitope prediction from PDB structure).
* **Immunogenicity scoring** (original ``immunogenicity.py``): Combined
  T-cell / B-cell immunogenicity scoring and deimmunization mutation
  suggestions.

Accuracy and Confidence
----------------------
**PSSM-based MHC binding prediction** (default, offline mode):
  - Expected AUC-ROC: 0.60–0.75 for MHC-I binding classification
  - This is significantly below state-of-the-art methods
  - PSSMs capture anchor position preferences but miss subtle
    peptide-MHC interaction features
  - IC50 estimates are rough approximations (log-linear mapping)
  - Binding classification thresholds (50/500/5000 nM) are standard
    but PSSM-derived IC50 values have high uncertainty

**NetMHCpan-based prediction** (when ``use_netmhcpan=True``):
  - Expected AUC-ROC: 0.85–0.95 for MHC-I binding (NetMHCpan 4.1)
  - This is the gold standard for computational MHC binding prediction
  - Requires API connectivity to the NetMHCpan web service
  - Falls back to PSSM if API is unavailable

**B-cell epitope prediction:**
  - Classical scale-based methods (Kolaskar-Tongaonkar, Parker, etc.)
    have typical AUC-ROC of 0.55–0.65
  - Performance varies significantly by epitope type and protein
  - Conformational epitope prediction (when PDB available) is more
    reliable than linear epitope prediction

**Deimmunization mutation suggestions:**
  - Confidence depends on the underlying binding prediction method
  - PSSM-based suggestions: **LOW** confidence
  - NetMHCpan-based suggestions: **MEDIUM-HIGH** confidence
  - Always verify experimentally before clinical use

**Upgrade path:**
  - Replace PSSMs with a neural network-based method for offline use
  - Add MHCflurry as an alternative offline predictor
  - Integrate BepiPred-2.0 for B-cell epitope prediction

  **Confidence levels:**
    - NetMHCpan mode: **HIGH** for MHC-I, **MEDIUM** for MHC-II
    - PSSM mode, strong anchor matches: **MEDIUM**
    - PSSM mode, weak anchor matches: **LOW**
    - B-cell epitope (linear): **LOW**
    - B-cell epitope (conformational with PDB): **MEDIUM**

All predictions are sequence-based heuristics and do not replace
experimental validation. When NetMHCpan integration is enabled
(``use_netmhcpan=True``), predictions use the NetMHCpan 4.1 API,
which provides significantly more accurate binding affinity estimates
than the PSSM-based approach.

References
----------
- Kolaskar & Tongaonkar, FEBS Lett 1990; 276:172-174
- Parker et al., Biochemistry 1986; 25:5424-5432
- Chou & Fasman, Biochemistry 1974; 13:222-245
- Emini et al., J Virol 1985; 55:836-839
- Reynisson et al., Nucleic Acids Res 2020; 48:W449 (NetMHCpan 4.1)
- O'Donnell et al., Bioinformatics 2018; 34:2696 (MHCflurry)
- Jespersen et al., Nucleic Acids Res 2017; 45:W39 (BepiPred-2.0)
"""

from __future__ import annotations

# ----------------------------------------------------------------------
# W8-a refactor: the implementation that formerly lived in this file has
# been split into cohesive sibling modules (``_constants``, ``_pssm``,
# ``_supertypes``, ``_models``, ``mhc_binding``, ``epitopes``,
# ``conformational``, ``scoring``).  This file is now a thin re-export
# shim so that every historical import path keeps working unchanged:
#   from biocompiler.immunogenicity.core import <name>
#   from biocompiler.immunogenicity import <name>
# ----------------------------------------------------------------------

from ._constants import *  # noqa: F401,F403
from ._pssm import *  # noqa: F401,F403
from ._supertypes import *  # noqa: F401,F403
from ._models import *  # noqa: F401,F403
from .mhc_binding import *  # noqa: F401,F403
from .epitopes import *  # noqa: F401,F403
from .conformational import *  # noqa: F401,F403
from .scoring import *  # noqa: F401,F403

# Re-export private helpers / state that are part of the historical
# public contract (imported by tests and sibling modules).
from ._constants import _STANDARD_AA_SET  # noqa: F401
from ._pssm import (  # noqa: F401
    _make_pssm_row, _build_mhc_i_pssms, _build_mhc_ii_pssms,
    _ensure_pssms_built, _get_mhc_i_pssms, _get_mhc_ii_pssms, _pssm_built,
)
from ._supertypes import (  # noqa: F401
    _MHC_SUPERTYPES, _ALLELE_PREFIX_TO_SUPERTYPE, _build_supertype_pssms,
    _ensure_supertype_pssms_built, _get_supertype_pssm,
    _SUPERTYPE_PSSM, _supertype_built,
)
from .mhc_binding import (  # noqa: F401
    _prediction_cache, _is_mhc_ii_allele, _identify_anchor_positions,
)
from .epitopes import (  # noqa: F401
    _validate_protein, _peptide_hydrophobicity_score, _peptide_charge_score,
    _score_peptide_for_allele, _sliding_window_average, _normalize_01,
    _find_regions, _add_region,
)
from .conformational import (  # noqa: F401
    _check_real_binding_data_available, _real_binding_data_cache,
)

# Re-export external names the original core exposed.
from biocompiler.shared.exceptions import ImmunogenicityError  # noqa: F401

__all__ = [
    "ImmunogenicityError",
    "MHCBindingResult",
    "MHCPredictionResult",
    "ImmunogenicityResult",
    "EpitopeRegion",
    "EpitopePredictionResult",
    "TCellEpitopeDict",
    "BCellEpitopeDict",
    "DEFAULT_MHC_I_ALLELES",
    "DEFAULT_MHC_II_ALLELES",
    "POPULATION_COVERAGE",
    "ANTIGENICITY_SCALE",
    "PARKER_SCALE",
    "CHOU_FASMAN_TURN",
    "EMINI_SCALE",
    "ALL_SCALES",
    "MHC_I_PSSM",
    "MHC_II_PSSM",
    "clear_cache",
    "score_peptide_pssm",
    "binding_score_to_ic50",
    "classify_binding",
    "predict_mhc_i_binding",
    "predict_mhc_ii_binding",
    "predict_all",
    "predict_t_cell_epitopes",
    "predict_kolaskar_tongaonkar",
    "predict_parker_hydrophilicity",
    "predict_chou_fasman_beta_turn",
    "predict_eea",
    "predict_bepipred_like",
    "predict_conformational_epitopes",
    "predict_epitopes",
    "compute_surface_accessibility_approx",
    "predict_b_cell_epitopes",
    "compute_immunogenicity",
    "find_deimmunization_mutations",
    "suggest_mutations",
    "compute_epitope_density",
    "compute_immunogenicity_batch",
    "IMMUNOGENICITY_PSSM_AUC_ROC_LOW",
    "IMMUNOGENICITY_PSSM_AUC_ROC_HIGH",
    "IMMUNOGENICITY_NETMHCPAN_AUC_ROC_LOW",
    "IMMUNOGENICITY_NETMHCPAN_AUC_ROC_HIGH",
    "IMMUNOGENICITY_BCELL_AUC_ROC",
    # Binding classification thresholds (nM)
    "IC50_STRONG_BINDER_THRESHOLD",
    "IC50_MODERATE_BINDER_THRESHOLD",
    "IC50_WEAK_BINDER_THRESHOLD",
    # IC50 mapping constants
    "IC50_LOG_INTERCEPT",
    "IC50_LOG_SLOPE",
    # PSSM scoring constants
    "PSSM_UNKNOWN_AA_SCORE",
    "PSSM_CONTRAST_POWER",
    # Hydrophobicity normalization
    "HYDROPHOBICITY_OFFSET",
    "HYDROPHOBICITY_RANGE",
    # Immunogenicity scoring weights
    "T_CELL_WEIGHT",
    "B_CELL_WEIGHT",
    # Immunogenicity classification thresholds
    "IMMUNOGENICITY_LOW_THRESHOLD",
    "IMMUNOGENICITY_HIGH_THRESHOLD",
    # Deimmunization limits
    "MAX_DEIMMUNIZATION_CANDIDATES",
    # MHC-II core peptide length
    "MHC_II_CORE_LENGTH",
    # Conformational epitope constants
    "CONF_EPITOPE_NEIGHBOR_CUTOPT_ANGSTROM",
    "CONF_EPITOPE_MAX_NEIGHBORS",
    # Epitope density constants
    "EPITOPE_DENSITY_CLUSTER_DISTANCE",
    # Offline prediction API
    "PeptideResult",
    "ImmunogenicityPrediction",
    "predict_immunogenicity",
    "scan_peptides",
    "PRECOMPUTED_BINDERS",
]
