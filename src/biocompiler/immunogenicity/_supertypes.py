"""Generic supertype PSSM fallbacks for MHC alleles.

When an allele has no dedicated PSSM, a supertype-family PSSM is used as
a fallback.  Split out of ``core.py`` (W8-a refactor).
"""
from __future__ import annotations

import concurrent.futures
import hashlib
import logging
import math
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, TypedDict

from biocompiler.shared.constants import (
    BLOSUM62,
    DEFAULT_MHC_PEPTIDE_LENGTH,
    HYDROPATHY,
    STANDARD_AAS,
)
from ..engines.base import (
    BaseEngineResult,
    BatchResult,
    EngineTimer,
    MutationResult,
    validate_protein_sequence,
)
from biocompiler.shared.exceptions import ImmunogenicityError
from biocompiler.shared.types import Verdict

logger = logging.getLogger(__name__)
from ._pssm import _make_pssm_row  # noqa: F401
from ._pssm import *  # noqa: F401,F403


# Removed: MHC_I_PREFERENCES, MHC_II_PREFERENCES (backward-compat
# aliases derived from PSSMs — use MHC_I_PSSM / MHC_II_PSSM instead).
# Removed: _DEFAULT_MHC_I_ALLELES, _DEFAULT_MHC_II_ALLELES,
# _DEFAULT_MHC_ALLELES (private backward-compat aliases).


# ═══════════════════════════════════════════════════════════════════════════
# Generic PSSM fallback for alleles without a dedicated PSSM
# ═══════════════════════════════════════════════════════════════════════════

#: MHC supertype families — maps common HLA alleles to their supertype
#: so that alleles without a dedicated PSSM can fall back to a generic
#: supertype PSSM.
_MHC_SUPERTYPES: dict[str, str] = {
    # A2 supertype
    "HLA-A*02:02": "A2", "HLA-A*02:03": "A2", "HLA-A*02:05": "A2",
    "HLA-A*02:06": "A2", "HLA-A*02:07": "A2", "HLA-A*02:09": "A2",
    "HLA-A*02:10": "A2", "HLA-A*02:11": "A2", "HLA-A*02:14": "A2",
    "HLA-A*02:16": "A2", "HLA-A*02:17": "A2", "HLA-A*02:19": "A2",
    "HLA-A*02:20": "A2", "HLA-A*02:24": "A2",
    "HLA-A*68:02": "A2",
    # A1 supertype
    "HLA-A*01:02": "A1", "HLA-A*01:03": "A1",
    "HLA-A*01:04": "A1", "HLA-A*01:06": "A1",
    "HLA-A*36:01": "A1", "HLA-A*80:01": "A1",
    # A3 supertype
    "HLA-A*03:02": "A3", "HLA-A*03:03": "A3",
    "HLA-A*03:04": "A3", "HLA-A*03:06": "A3",
    "HLA-A*11:01": "A3", "HLA-A*11:02": "A3", "HLA-A*11:03": "A3",
    "HLA-A*31:01": "A3", "HLA-A*31:02": "A3",
    "HLA-A*33:01": "A3", "HLA-A*33:03": "A3",
    "HLA-A*68:01": "A3",
    "HLA-A*30:01": "A3", "HLA-A*30:02": "A3",
    # A24 supertype
    "HLA-A*24:03": "A24", "HLA-A*24:07": "A24",
    "HLA-A*24:08": "A24", "HLA-A*24:10": "A24",
    "HLA-A*23:01": "A24", "HLA-A*23:02": "A24",
    "HLA-A*30:03": "A24",
    # B7 supertype
    "HLA-B*07:03": "B7", "HLA-B*07:04": "B7", "HLA-B*07:05": "B7",
    "HLA-B*35:01": "B7", "HLA-B*35:03": "B7", "HLA-B*35:08": "B7",
    "HLA-B*51:01": "B7", "HLA-B*51:02": "B7",
    "HLA-B*53:01": "B7",
    "HLA-B*54:01": "B7", "HLA-B*54:02": "B7",
    "HLA-B*55:01": "B7", "HLA-B*55:02": "B7",
    "HLA-B*56:01": "B7", "HLA-B*56:02": "B7",
    "HLA-B*67:01": "B7",
    "HLA-B*78:01": "B7",
    # B8 supertype
    "HLA-B*08:02": "B8", "HLA-B*08:03": "B8",
    # B44 supertype
    "HLA-B*44:02": "B44", "HLA-B*44:03": "B44",
    "HLA-B*40:01": "B44", "HLA-B*40:02": "B44", "HLA-B*40:06": "B44",
    "HLA-B*45:01": "B44",
    # B58 supertype
    "HLA-B*58:01": "B58", "HLA-B*58:02": "B58",
    "HLA-B*57:01": "B58", "HLA-B*57:02": "B58",
    # B62 supertype
    "HLA-B*62:01": "B62", "HLA-B*62:02": "B62",
    "HLA-B*15:01": "B62", "HLA-B*15:02": "B62", "HLA-B*15:03": "B62",
    "HLA-B*46:01": "B62",
    "HLA-B*52:01": "B62",
    # B27 supertype
    "HLA-B*27:02": "B27", "HLA-B*27:03": "B27", "HLA-B*27:04": "B27",
    "HLA-B*27:05": "B27", "HLA-B*27:06": "B27", "HLA-B*27:07": "B27",
    "HLA-B*38:01": "B27", "HLA-B*39:01": "B27",
    # DR supertypes
    "HLA-DRB1*01:02": "DR1", "HLA-DRB1*01:03": "DR1",
    "HLA-DRB1*04:02": "DR4", "HLA-DRB1*04:03": "DR4",
    "HLA-DRB1*04:04": "DR4", "HLA-DRB1*04:05": "DR4",
    "HLA-DRB1*04:06": "DR4", "HLA-DRB1*04:07": "DR4", "HLA-DRB1*04:08": "DR4",
    "HLA-DRB1*07:02": "DR7", "HLA-DRB1*07:03": "DR7",
    "HLA-DRB1*03:01": "DR3", "HLA-DRB1*03:02": "DR3", "HLA-DRB1*03:03": "DR3",
    "HLA-DRB1*11:01": "DR11", "HLA-DRB1*11:02": "DR11", "HLA-DRB1*11:03": "DR11",
    "HLA-DRB1*11:04": "DR11",
    "HLA-DRB1*13:01": "DR11", "HLA-DRB1*13:02": "DR11", "HLA-DRB1*13:03": "DR11",
    "HLA-DRB1*15:01": "DR15", "HLA-DRB1*15:02": "DR15", "HLA-DRB1*15:03": "DR15",
    "HLA-DRB1*16:01": "DR15", "HLA-DRB1*16:02": "DR15",
    "HLA-DRB1*08:01": "DR8", "HLA-DRB1*08:02": "DR8",
    "HLA-DRB1*09:01": "DR9",
    "HLA-DRB1*10:01": "DR10",
    "HLA-DRB1*12:01": "DR11", "HLA-DRB1*12:02": "DR11",
    "HLA-DRB1*14:01": "DR11",
    "HLA-DRB5*01:01": "DR11", "HLA-DRB5*02:02": "DR15",
}

#: Heuristic supertype mapping based on allele family prefixes.
#: When an exact allele-to-supertype mapping is not found in
#: :data:`_MHC_SUPERTYPES`, this table is consulted by extracting
#: the gene-group prefix (e.g. ``"HLA-A*02"`` from ``"HLA-A*02:14"``)
#: and looking it up here.
_ALLELE_PREFIX_TO_SUPERTYPE: dict[str, str] = {
    # MHC-I A locus
    "HLA-A*01": "A1",
    "HLA-A*02": "A2",
    "HLA-A*03": "A3",
    "HLA-A*11": "A3",
    "HLA-A*23": "A24",
    "HLA-A*24": "A24",
    "HLA-A*30": "A3",
    "HLA-A*31": "A3",
    "HLA-A*33": "A3",
    "HLA-A*36": "A1",
    "HLA-A*68": "A2",  # A*68:01 is A3, A*68:02 is A2; majority rule: A2
    "HLA-A*69": "A2",
    "HLA-A*80": "A1",
    # MHC-I B locus
    "HLA-B*07": "B7",
    "HLA-B*08": "B8",
    "HLA-B*13": "B44",
    "HLA-B*14": "B44",
    "HLA-B*15": "B62",
    "HLA-B*18": "B44",
    "HLA-B*27": "B27",
    "HLA-B*35": "B7",
    "HLA-B*38": "B27",
    "HLA-B*39": "B27",
    "HLA-B*40": "B44",
    "HLA-B*44": "B44",
    "HLA-B*45": "B44",
    "HLA-B*46": "B62",
    "HLA-B*47": "B44",
    "HLA-B*48": "B7",
    "HLA-B*49": "B7",
    "HLA-B*50": "B44",
    "HLA-B*51": "B7",
    "HLA-B*52": "B62",
    "HLA-B*53": "B7",
    "HLA-B*54": "B7",
    "HLA-B*55": "B7",
    "HLA-B*56": "B7",
    "HLA-B*57": "B58",
    "HLA-B*58": "B58",
    "HLA-B*62": "B62",
    "HLA-B*67": "B7",
    "HLA-B*78": "B7",
    # MHC-II DR locus
    "HLA-DRB1*01": "DR1",
    "HLA-DRB1*03": "DR3",
    "HLA-DRB1*04": "DR4",
    "HLA-DRB1*07": "DR7",
    "HLA-DRB1*08": "DR8",
    "HLA-DRB1*09": "DR9",
    "HLA-DRB1*10": "DR10",
    "HLA-DRB1*11": "DR11",
    "HLA-DRB1*12": "DR11",
    "HLA-DRB1*13": "DR11",
    "HLA-DRB1*14": "DR11",
    "HLA-DRB1*15": "DR15",
    "HLA-DRB1*16": "DR15",
}


def _build_supertype_pssms() -> dict[str, list[dict[str, float]]]:
    """Construct generic PSSMs for MHC supertypes.

    These serve as fallbacks when a specific allele does not have its
    own dedicated PSSM but belongs to a known supertype family.

    .. warning::
        These supertype PSSMs are **even less accurate** than the
        allele-specific PSSMs above.  They average over multiple alleles
        within a supertype family and thus lose allele-specific detail.
        For any production use, prefer ``use_netmhcpan=True`` or use
        allele-specific PSSMs (via NetMHCIIpan for MHC-II).
    """
    pssms: dict[str, list[dict[str, float]]] = {}

    # A2 supertype — hydrophobic anchors at P2 and P9
    pssms["A2"] = [
        _make_pssm_row(
            preferred={"L": 1.1, "M": 1.1, "I": 1.1, "V": 1.1, "A": 1.05},
            disfavored={"D": 0.6, "E": 0.6, "K": 0.6, "R": 0.6},
        ),
        _make_pssm_row(
            preferred={"L": 1.8, "M": 1.8, "I": 1.6, "V": 1.6},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.4, "R": 0.4},
            default=0.8,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"V": 1.4, "L": 1.4, "I": 1.2, "A": 1.1},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5, "R": 0.5},
            default=0.8,
        ),
    ]

    # A1 supertype — small/serine at P2, tyrosine at P9
    pssms["A1"] = [
        _make_pssm_row(
            preferred={"A": 1.05, "S": 1.05},
            disfavored={"W": 0.6, "R": 0.6},
        ),
        _make_pssm_row(
            preferred={"T": 1.6, "S": 1.4, "D": 1.3, "E": 1.3},
            disfavored={"L": 0.5, "I": 0.5, "V": 0.6},
            default=0.85,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"Y": 1.6, "F": 1.4},
            disfavored={"K": 0.5, "R": 0.5, "D": 0.6},
            default=0.85,
        ),
    ]

    # A3 supertype — hydrophobic at P2, basic at P9
    pssms["A3"] = [
        _make_pssm_row(
            preferred={"A": 1.05, "S": 1.05},
        ),
        _make_pssm_row(
            preferred={"V": 1.6, "I": 1.6, "L": 1.4, "M": 1.4},
            disfavored={"D": 0.5, "E": 0.5},
            default=0.85,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.8, "R": 1.6, "H": 1.3},
            disfavored={"D": 0.4, "E": 0.4},
            default=0.75,
        ),
    ]

    # A24 supertype — aromatic at P2, hydrophobic at P9
    pssms["A24"] = [
        _make_pssm_row(
            preferred={"Y": 1.1, "F": 1.05},
        ),
        _make_pssm_row(
            preferred={"Y": 1.8, "F": 1.8, "W": 1.5},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.4},
            default=0.75,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"F": 1.4, "L": 1.4, "I": 1.2},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.5},
            default=0.85,
        ),
    ]

    # B7 supertype — proline at P2, hydrophobic at P9
    pssms["B7"] = [
        _make_pssm_row(
            preferred={"A": 1.05, "P": 1.05},
        ),
        _make_pssm_row(
            preferred={"P": 1.8, "A": 1.6},
            disfavored={"D": 0.4, "E": 0.4, "K": 0.4},
            default=0.75,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.4, "I": 1.4, "V": 1.2},
            disfavored={"D": 0.5, "E": 0.5},
            default=0.85,
        ),
    ]

    # B8 supertype — basic at P2, hydrophobic at P9
    pssms["B8"] = [
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.6, "R": 1.6},
            disfavored={"D": 0.4, "E": 0.4, "P": 0.5},
            default=0.85,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.4, "I": 1.2, "V": 1.2},
            disfavored={"D": 0.5, "E": 0.5},
            default=0.85,
        ),
    ]

    # B44 supertype — glutamate at P2, hydrophobic/aromatic at P9
    pssms["B44"] = [
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"E": 1.8, "D": 1.5},
            disfavored={"K": 0.4, "R": 0.4},
            default=0.8,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"Y": 1.5, "F": 1.4, "L": 1.3, "I": 1.2},
            disfavored={"D": 0.5, "E": 0.5},
            default=0.85,
        ),
    ]

    # B58 supertype — small/serine at P2, aromatic/hydrophobic at P9
    pssms["B58"] = [
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"S": 1.5, "A": 1.4, "T": 1.3},
            disfavored={"W": 0.5, "F": 0.6},
            default=0.85,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"W": 1.6, "F": 1.4, "Y": 1.4, "L": 1.2},
            disfavored={"D": 0.5, "E": 0.5},
            default=0.85,
        ),
    ]

    # B62 supertype — glutamine at P2, hydrophobic at P9
    pssms["B62"] = [
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"Q": 1.6, "N": 1.3},
            disfavored={"D": 0.5, "E": 0.5},
            default=0.85,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.4, "I": 1.3, "V": 1.2, "F": 1.3},
            disfavored={"D": 0.5, "E": 0.5},
            default=0.85,
        ),
    ]

    # DR1 supertype — hydrophobic at P1, small at P4, hydrophobic at P6
    pssms["DR1"] = [
        _make_pssm_row(
            preferred={"F": 1.6, "Y": 1.5, "W": 1.4, "L": 1.3, "I": 1.2, "V": 1.2},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.4, "S": 1.3, "T": 1.3, "N": 1.2},
            disfavored={"W": 0.6, "F": 0.7},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.4, "I": 1.3, "V": 1.3, "M": 1.2},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.2, "R": 1.2, "N": 1.1, "Q": 1.1},
        ),
    ]

    # DR4 supertype — hydrophobic at P1, acidic at P4, small at P6
    pssms["DR4"] = [
        _make_pssm_row(
            preferred={"F": 1.6, "Y": 1.5, "W": 1.4, "L": 1.2},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"D": 1.6, "E": 1.4},
            disfavored={"K": 0.5, "R": 0.5},
            default=0.85,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.4, "S": 1.3, "G": 1.2},
            disfavored={"W": 0.6, "F": 0.7},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.3, "I": 1.2, "V": 1.2},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
    ]

    # DR7 supertype — hydrophobic at P1, small at P4, hydrophobic at P6
    pssms["DR7"] = [
        _make_pssm_row(
            preferred={"F": 1.5, "Y": 1.4, "L": 1.3, "I": 1.2, "V": 1.2},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.3, "S": 1.2, "T": 1.2},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.3, "I": 1.2, "V": 1.2, "F": 1.2},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.2, "R": 1.2, "N": 1.1, "Q": 1.1},
        ),
    ]

    # DR3 supertype — similar to DR1 with broader anchor tolerance
    pssms["DR3"] = pssms["DR1"][:]

    # DR11 supertype — similar to DR1 with intermediate P4 preference
    pssms["DR11"] = [
        _make_pssm_row(
            preferred={"F": 1.5, "Y": 1.4, "W": 1.3, "L": 1.2, "I": 1.1},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.3, "S": 1.2, "N": 1.2, "G": 1.1},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.3, "I": 1.2, "V": 1.2},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
    ]

    # DR15 supertype — similar to DR1 with P6 preferring small residues
    pssms["DR15"] = [
        _make_pssm_row(
            preferred={"F": 1.4, "Y": 1.3, "L": 1.2, "V": 1.1},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.3, "S": 1.2},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.3, "S": 1.2, "G": 1.2},
            disfavored={"W": 0.6, "F": 0.7},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.2, "I": 1.2, "V": 1.1},
            default=0.9,
        ),
    ]

    # B27 supertype — arginine at P2, basic/tyrosine at P9
    pssms["B27"] = [
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"R": 1.8, "K": 1.5, "Q": 1.3},
            disfavored={"D": 0.4, "E": 0.4, "P": 0.5},
            default=0.85,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.5, "R": 1.4, "Y": 1.3, "F": 1.2},
            disfavored={"D": 0.5, "E": 0.5},
            default=0.85,
        ),
    ]

    # DR8 supertype — similar to DR1 with broader P4 tolerance
    pssms["DR8"] = [
        _make_pssm_row(
            preferred={"F": 1.4, "Y": 1.3, "L": 1.2, "V": 1.1},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.3, "S": 1.2, "N": 1.1},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.3, "I": 1.2, "V": 1.2},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.2, "R": 1.1, "N": 1.1},
        ),
    ]

    # DR9 supertype — similar to DR4 with hydrophobic at P1
    pssms["DR9"] = [
        _make_pssm_row(
            preferred={"F": 1.5, "Y": 1.4, "W": 1.3, "L": 1.2},
            disfavored={"D": 0.5, "E": 0.5, "K": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"D": 1.5, "E": 1.3, "A": 1.2},
            disfavored={"K": 0.5, "R": 0.5},
            default=0.85,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.3, "S": 1.2, "G": 1.1},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.2, "I": 1.1, "V": 1.1},
            default=0.9,
        ),
    ]

    # DR10 supertype — similar to DR1 with broader P1 tolerance
    pssms["DR10"] = [
        _make_pssm_row(
            preferred={"F": 1.4, "Y": 1.3, "L": 1.2, "I": 1.1, "V": 1.1},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"A": 1.3, "S": 1.2},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"L": 1.3, "I": 1.2, "V": 1.2},
            disfavored={"D": 0.6, "E": 0.6},
            default=0.9,
        ),
        _make_pssm_row(default=1.0),
        _make_pssm_row(default=1.0),
        _make_pssm_row(
            preferred={"K": 1.2, "R": 1.1, "Q": 1.1},
        ),
    ]

    return pssms


_SUPERTYPE_PSSM: dict[str, list[dict[str, float]]] = {}
_supertype_built: bool = False


def _ensure_supertype_pssms_built() -> None:
    """Build supertype PSSMs on first access."""
    global _SUPERTYPE_PSSM, _supertype_built
    if not _supertype_built:
        _SUPERTYPE_PSSM = _build_supertype_pssms()
        _supertype_built = True


def _get_supertype_pssm(allele: str) -> list[dict[str, float]] | None:
    """Get a supertype fallback PSSM for an allele.

    Resolution strategy (in order):

    1. **Exact lookup**: check :data:`_MHC_SUPERTYPES` for a direct
       allele → supertype mapping.
    2. **Prefix heuristic**: extract the gene-group prefix from the
       allele name (e.g. ``"HLA-A*02"`` from ``"HLA-A*02:99"``)
       and look it up in :data:`_ALLELE_PREFIX_TO_SUPERTYPE`.

    Returns None if no supertype mapping can be found by either method.
    """
    _ensure_supertype_pssms_built()

    # Strategy 1: exact allele match
    supertype = _MHC_SUPERTYPES.get(allele)
    if supertype is not None:
        pssm = _SUPERTYPE_PSSM.get(supertype)
        if pssm is not None:
            return pssm

    # Strategy 2: prefix-based heuristic fallback
    # Extract prefix like "HLA-A*02" from "HLA-A*02:99"
    colon_idx = allele.find(":")
    if colon_idx > 0:
        prefix = allele[:colon_idx]
        supertype = _ALLELE_PREFIX_TO_SUPERTYPE.get(prefix)
        if supertype is not None:
            pssm = _SUPERTYPE_PSSM.get(supertype)
            if pssm is not None:
                logger.debug(
                    "Allele %s resolved to supertype %s via prefix %s",
                    allele, supertype, prefix,
                )
                return pssm

    logger.debug("No supertype PSSM fallback found for allele %s", allele)
    return None

__all__ = []  # all names here are private helpers
