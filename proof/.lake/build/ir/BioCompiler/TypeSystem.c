// Lean compiler output
// Module: BioCompiler.TypeSystem
// Imports: public import Init public meta import Init public import BioCompiler.ThreeValued public import BioCompiler.Sequence public import BioCompiler.NDFST public import BioCompiler.Scanners
#include <lean/lean.h>
#if defined(__clang__)
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-label"
#elif defined(__GNUC__) && !defined(__CLANG__)
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-label"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#ifdef __cplusplus
extern "C" {
#endif
uint8_t l_Rat_instDecidableLe(lean_object*, lean_object*);
uint8_t lean_string_dec_eq(lean_object*, lean_object*);
lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg(lean_object*, lean_object*);
lean_object* lp_BioCompiler_BioCompiler_Sequence_gcContent(lean_object*);
uint8_t lp_BioCompiler_BioCompiler_hasAnyRestrictionSite(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_hasInstabilityMotif(lean_object*);
extern lean_object* lp_BioCompiler_BioCompiler_spliceDonorConsensus;
uint8_t lp_BioCompiler_BioCompiler_Sequence_containsPattern(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_isValidCodingSeq(lean_object*);
lean_object* lean_string_length(lean_object*);
lean_object* l_Nat_reprFast(lean_object*);
lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(lean_object*);
lean_object* lean_nat_to_int(lean_object*);
lean_object* l_Std_Format_fill(lean_object*);
lean_object* l_Repr_addAppParen(lean_object*, lean_object*);
lean_object* l_String_quote(lean_object*);
uint8_t lean_nat_dec_le(lean_object*, lean_object*);
uint8_t lean_nat_dec_eq(lean_object*, lean_object*);
lean_object* l_Int_repr(lean_object*);
lean_object* lean_string_append(lean_object*, lean_object*);
uint8_t lean_int_dec_lt(lean_object*, lean_object*);
lean_object* l_Bool_repr___redArg(uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorIdx(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorIdx___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorElim(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SpliceCorrect_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SpliceCorrect_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticSplice_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticSplice_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonAdapted_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonAdapted_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_GCInRange_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_GCInRange_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoRestrictionSite_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoRestrictionSite_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_InFrame_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_InFrame_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoInstabilityMotif_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoInstabilityMotif_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCpGIsland_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCpGIsland_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoGTDinucleotide_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoGTDinucleotide_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoStopCodons_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoStopCodons_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ValidCodingSeq_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ValidCodingSeq_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonOptimality_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonOptimality_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticPromoter_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticPromoter_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ConservationScore_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ConservationScore_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoUnexpectedTMDomain_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoUnexpectedTMDomain_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_mRNASecondaryStructure_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_mRNASecondaryStructure_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CoTranslationalFolding_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CoTranslationalFolding_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_StructureConfidence_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_StructureConfidence_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoMisfoldingRisk_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoMisfoldingRisk_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CorrectFoldTopology_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CorrectFoldTopology_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoUnexpectedInteraction_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoUnexpectedInteraction_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_StableFolding_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_StableFolding_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoDestabilizingMutation_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoDestabilizingMutation_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_DisulfideBondIntegrity_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_DisulfideBondIntegrity_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_HydrophobicCoreQuality_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_HydrophobicCoreQuality_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SolubleExpression_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SolubleExpression_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoAggregationProneRegion_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoAggregationProneRegion_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ChargeComposition_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ChargeComposition_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoLongHydrophobicStretch_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoLongHydrophobicStretch_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_LowImmunogenicity_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_LowImmunogenicity_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoStrongTCellEpitope_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoStrongTCellEpitope_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoDominantBCellEpitope_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoDominantBCellEpitope_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_PopulationCoverageSafe_elim___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_PopulationCoverageSafe_elim(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1_spec__3(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0(lean_object*, lean_object*);
static const lean_string_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 3, .m_capacity = 3, .m_length = 2, .m_data = "[]"};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__1_value;
static const lean_string_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "["};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__2_value;
static const lean_string_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = ","};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__3_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__3_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__4_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__4_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__5_value;
static const lean_string_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "]"};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__6 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__6_value;
static lean_once_cell_t lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7;
static lean_once_cell_t lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__9_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__2_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__9 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__9_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__10_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__6_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__10 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__10_value;
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2___lam__0(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4_spec__6(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___redArg(lean_object*);
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 51, .m_capacity = 51, .m_length = 50, .m_data = "BioCompiler.TypePredicate.NoAggregationProneRegion"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__1_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 49, .m_capacity = 49, .m_length = 48, .m_data = "BioCompiler.TypePredicate.DisulfideBondIntegrity"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__2_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 50, .m_capacity = 50, .m_length = 49, .m_data = "BioCompiler.TypePredicate.NoUnexpectedInteraction"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__5_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 46, .m_capacity = 46, .m_length = 45, .m_data = "BioCompiler.TypePredicate.CorrectFoldTopology"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__6 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__6_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__7_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__6_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__7 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__7_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 43, .m_capacity = 43, .m_length = 42, .m_data = "BioCompiler.TypePredicate.NoMisfoldingRisk"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__8 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__8_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__8_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 42, .m_capacity = 42, .m_length = 41, .m_data = "BioCompiler.TypePredicate.NoCrypticSplice"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__11_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__11 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__11_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__12_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 45, .m_capacity = 45, .m_length = 44, .m_data = "BioCompiler.TypePredicate.NoInstabilityMotif"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__12 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__12_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__13_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__12_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__13 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__13_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 38, .m_capacity = 38, .m_length = 37, .m_data = "BioCompiler.TypePredicate.NoCpGIsland"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 43, .m_capacity = 43, .m_length = 42, .m_data = "BioCompiler.TypePredicate.NoGTDinucleotide"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__17_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__17 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__17_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__18_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 39, .m_capacity = 39, .m_length = 38, .m_data = "BioCompiler.TypePredicate.NoStopCodons"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__18 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__18_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__19_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__18_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__19 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__19_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__20_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 41, .m_capacity = 41, .m_length = 40, .m_data = "BioCompiler.TypePredicate.ValidCodingSeq"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__20 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__20_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__21_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__20_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__21 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__21_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__22_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 40, .m_capacity = 40, .m_length = 39, .m_data = "BioCompiler.TypePredicate.SpliceCorrect"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__22 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__22_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__23_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__22_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__23 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__23_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__24_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__23_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__24 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__24_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__27_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 39, .m_capacity = 39, .m_length = 38, .m_data = "BioCompiler.TypePredicate.CodonAdapted"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__27 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__27_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__28_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__27_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__28 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__28_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__29_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__28_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__29 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__29_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "("};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = " : Rat)/"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__33_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 36, .m_capacity = 36, .m_length = 35, .m_data = "BioCompiler.TypePredicate.GCInRange"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__33 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__33_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__34_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__33_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__34 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__34_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__35_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__34_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__35 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__35_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__36_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 44, .m_capacity = 44, .m_length = 43, .m_data = "BioCompiler.TypePredicate.NoRestrictionSite"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__36 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__36_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__37_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__36_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__37 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__37_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__38_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__37_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__38 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__38_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__39_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 34, .m_capacity = 34, .m_length = 33, .m_data = "BioCompiler.TypePredicate.InFrame"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__39 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__39_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__40_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__39_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__40 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__40_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__41_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__40_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__41 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__41_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__42_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 42, .m_capacity = 42, .m_length = 41, .m_data = "BioCompiler.TypePredicate.CodonOptimality"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__42 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__42_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__43_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__42_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__43 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__43_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__44_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__43_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__44 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__44_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__45_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 44, .m_capacity = 44, .m_length = 43, .m_data = "BioCompiler.TypePredicate.NoCrypticPromoter"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__45 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__45_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__46_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__45_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__46 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__46_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__47_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__46_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__47 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__47_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__48_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 44, .m_capacity = 44, .m_length = 43, .m_data = "BioCompiler.TypePredicate.ConservationScore"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__48 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__48_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__49_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__48_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__49 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__49_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__50_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__49_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__50 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__50_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__51_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 47, .m_capacity = 47, .m_length = 46, .m_data = "BioCompiler.TypePredicate.NoUnexpectedTMDomain"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__51 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__51_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__52_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__51_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__52 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__52_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__53_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__52_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__53 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__53_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__54_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 49, .m_capacity = 49, .m_length = 48, .m_data = "BioCompiler.TypePredicate.mRNASecondaryStructure"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__54 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__54_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__55_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__54_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__55 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__55_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__56_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__55_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__56 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__56_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__57_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 49, .m_capacity = 49, .m_length = 48, .m_data = "BioCompiler.TypePredicate.CoTranslationalFolding"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__57 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__57_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__58_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__57_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__58 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__58_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__59_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__58_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__59 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__59_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__60_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 46, .m_capacity = 46, .m_length = 45, .m_data = "BioCompiler.TypePredicate.StructureConfidence"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__60 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__60_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__61_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__60_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__61 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__61_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__62_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__61_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__62 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__62_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__63_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 40, .m_capacity = 40, .m_length = 39, .m_data = "BioCompiler.TypePredicate.StableFolding"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__63 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__63_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__64_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__63_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__64 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__64_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__65_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__64_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__65 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__65_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__66_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 50, .m_capacity = 50, .m_length = 49, .m_data = "BioCompiler.TypePredicate.NoDestabilizingMutation"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__66 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__66_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__67_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__66_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__67 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__67_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__68_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__67_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__68 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__68_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__69_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 49, .m_capacity = 49, .m_length = 48, .m_data = "BioCompiler.TypePredicate.HydrophobicCoreQuality"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__69 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__69_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__70_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__69_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__70 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__70_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__71_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__70_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__71 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__71_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__72_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 44, .m_capacity = 44, .m_length = 43, .m_data = "BioCompiler.TypePredicate.SolubleExpression"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__72 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__72_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__73_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__72_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__73 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__73_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__74_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__73_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__74 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__74_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__75_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 44, .m_capacity = 44, .m_length = 43, .m_data = "BioCompiler.TypePredicate.ChargeComposition"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__75 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__75_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__76_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__75_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__76 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__76_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__77_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__76_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__77 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__77_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__78_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 51, .m_capacity = 51, .m_length = 50, .m_data = "BioCompiler.TypePredicate.NoLongHydrophobicStretch"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__78 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__78_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__79_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__78_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__79 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__79_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__80_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__79_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__80 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__80_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__81_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 44, .m_capacity = 44, .m_length = 43, .m_data = "BioCompiler.TypePredicate.LowImmunogenicity"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__81 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__81_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__82_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__81_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__82 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__82_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__83_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__82_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__83 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__83_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__84_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 47, .m_capacity = 47, .m_length = 46, .m_data = "BioCompiler.TypePredicate.NoStrongTCellEpitope"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__84 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__84_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__85_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__84_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__85 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__85_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__86_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__85_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__86 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__86_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__87_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 49, .m_capacity = 49, .m_length = 48, .m_data = "BioCompiler.TypePredicate.NoDominantBCellEpitope"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__87 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__87_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__88_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__87_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__88 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__88_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__89_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__88_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__89 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__89_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__90_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 49, .m_capacity = 49, .m_length = 48, .m_data = "BioCompiler.TypePredicate.PopulationCoverageSafe"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__90 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__90_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__91_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__90_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__91 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__91_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__92_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__91_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__92 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__92_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprTypePredicate___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_isSLOT(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_isSLOT___boxed(lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluate___redArg(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluate___redArg___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluate(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluate___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter___redArg(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter___redArg___boxed(lean_object**);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter___boxed(lean_object**);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__1_splitter___redArg(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__1_splitter(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorIdx(lean_object* v_x_1_){
_start:
{
switch(lean_obj_tag(v_x_1_))
{
case 0:
{
lean_object* v___x_2_; 
v___x_2_ = lean_unsigned_to_nat(0u);
return v___x_2_;
}
case 1:
{
lean_object* v___x_3_; 
v___x_3_ = lean_unsigned_to_nat(1u);
return v___x_3_;
}
case 2:
{
lean_object* v___x_4_; 
v___x_4_ = lean_unsigned_to_nat(2u);
return v___x_4_;
}
case 3:
{
lean_object* v___x_5_; 
v___x_5_ = lean_unsigned_to_nat(3u);
return v___x_5_;
}
case 4:
{
lean_object* v___x_6_; 
v___x_6_ = lean_unsigned_to_nat(4u);
return v___x_6_;
}
case 5:
{
lean_object* v___x_7_; 
v___x_7_ = lean_unsigned_to_nat(5u);
return v___x_7_;
}
case 6:
{
lean_object* v___x_8_; 
v___x_8_ = lean_unsigned_to_nat(6u);
return v___x_8_;
}
case 7:
{
lean_object* v___x_9_; 
v___x_9_ = lean_unsigned_to_nat(7u);
return v___x_9_;
}
case 8:
{
lean_object* v___x_10_; 
v___x_10_ = lean_unsigned_to_nat(8u);
return v___x_10_;
}
case 9:
{
lean_object* v___x_11_; 
v___x_11_ = lean_unsigned_to_nat(9u);
return v___x_11_;
}
case 10:
{
lean_object* v___x_12_; 
v___x_12_ = lean_unsigned_to_nat(10u);
return v___x_12_;
}
case 11:
{
lean_object* v___x_13_; 
v___x_13_ = lean_unsigned_to_nat(11u);
return v___x_13_;
}
case 12:
{
lean_object* v___x_14_; 
v___x_14_ = lean_unsigned_to_nat(12u);
return v___x_14_;
}
case 13:
{
lean_object* v___x_15_; 
v___x_15_ = lean_unsigned_to_nat(13u);
return v___x_15_;
}
case 14:
{
lean_object* v___x_16_; 
v___x_16_ = lean_unsigned_to_nat(14u);
return v___x_16_;
}
case 15:
{
lean_object* v___x_17_; 
v___x_17_ = lean_unsigned_to_nat(15u);
return v___x_17_;
}
case 16:
{
lean_object* v___x_18_; 
v___x_18_ = lean_unsigned_to_nat(16u);
return v___x_18_;
}
case 17:
{
lean_object* v___x_19_; 
v___x_19_ = lean_unsigned_to_nat(17u);
return v___x_19_;
}
case 18:
{
lean_object* v___x_20_; 
v___x_20_ = lean_unsigned_to_nat(18u);
return v___x_20_;
}
case 19:
{
lean_object* v___x_21_; 
v___x_21_ = lean_unsigned_to_nat(19u);
return v___x_21_;
}
case 20:
{
lean_object* v___x_22_; 
v___x_22_ = lean_unsigned_to_nat(20u);
return v___x_22_;
}
case 21:
{
lean_object* v___x_23_; 
v___x_23_ = lean_unsigned_to_nat(21u);
return v___x_23_;
}
case 22:
{
lean_object* v___x_24_; 
v___x_24_ = lean_unsigned_to_nat(22u);
return v___x_24_;
}
case 23:
{
lean_object* v___x_25_; 
v___x_25_ = lean_unsigned_to_nat(23u);
return v___x_25_;
}
case 24:
{
lean_object* v___x_26_; 
v___x_26_ = lean_unsigned_to_nat(24u);
return v___x_26_;
}
case 25:
{
lean_object* v___x_27_; 
v___x_27_ = lean_unsigned_to_nat(25u);
return v___x_27_;
}
case 26:
{
lean_object* v___x_28_; 
v___x_28_ = lean_unsigned_to_nat(26u);
return v___x_28_;
}
case 27:
{
lean_object* v___x_29_; 
v___x_29_ = lean_unsigned_to_nat(27u);
return v___x_29_;
}
case 28:
{
lean_object* v___x_30_; 
v___x_30_ = lean_unsigned_to_nat(28u);
return v___x_30_;
}
case 29:
{
lean_object* v___x_31_; 
v___x_31_ = lean_unsigned_to_nat(29u);
return v___x_31_;
}
case 30:
{
lean_object* v___x_32_; 
v___x_32_ = lean_unsigned_to_nat(30u);
return v___x_32_;
}
case 31:
{
lean_object* v___x_33_; 
v___x_33_ = lean_unsigned_to_nat(31u);
return v___x_33_;
}
default: 
{
lean_object* v___x_34_; 
v___x_34_ = lean_unsigned_to_nat(32u);
return v___x_34_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorIdx___boxed(lean_object* v_x_35_){
_start:
{
lean_object* v_res_36_; 
v_res_36_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorIdx(v_x_35_);
lean_dec(v_x_35_);
return v_res_36_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(lean_object* v_t_37_, lean_object* v_k_38_){
_start:
{
switch(lean_obj_tag(v_t_37_))
{
case 0:
{
lean_object* v_cellType_39_; lean_object* v___x_40_; 
v_cellType_39_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_cellType_39_);
lean_dec_ref(v_t_37_);
v___x_40_ = lean_apply_1(v_k_38_, v_cellType_39_);
return v___x_40_;
}
case 2:
{
lean_object* v_organism_41_; lean_object* v_threshold_42_; lean_object* v___x_43_; 
v_organism_41_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_organism_41_);
v_threshold_42_ = lean_ctor_get(v_t_37_, 1);
lean_inc_ref(v_threshold_42_);
lean_dec_ref(v_t_37_);
v___x_43_ = lean_apply_2(v_k_38_, v_organism_41_, v_threshold_42_);
return v___x_43_;
}
case 3:
{
lean_object* v_lo_44_; lean_object* v_hi_45_; lean_object* v___x_46_; 
v_lo_44_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_lo_44_);
v_hi_45_ = lean_ctor_get(v_t_37_, 1);
lean_inc_ref(v_hi_45_);
lean_dec_ref(v_t_37_);
v___x_46_ = lean_apply_2(v_k_38_, v_lo_44_, v_hi_45_);
return v___x_46_;
}
case 4:
{
lean_object* v_enzymeSites_47_; lean_object* v___x_48_; 
v_enzymeSites_47_ = lean_ctor_get(v_t_37_, 0);
lean_inc(v_enzymeSites_47_);
lean_dec_ref(v_t_37_);
v___x_48_ = lean_apply_1(v_k_38_, v_enzymeSites_47_);
return v___x_48_;
}
case 5:
{
lean_object* v_readingFrame_49_; lean_object* v_exonBoundaries_50_; lean_object* v___x_51_; 
v_readingFrame_49_ = lean_ctor_get(v_t_37_, 0);
lean_inc(v_readingFrame_49_);
v_exonBoundaries_50_ = lean_ctor_get(v_t_37_, 1);
lean_inc(v_exonBoundaries_50_);
lean_dec_ref(v_t_37_);
v___x_51_ = lean_apply_2(v_k_38_, v_readingFrame_49_, v_exonBoundaries_50_);
return v___x_51_;
}
case 11:
{
lean_object* v_organism_52_; lean_object* v_threshold_53_; lean_object* v___x_54_; 
v_organism_52_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_organism_52_);
v_threshold_53_ = lean_ctor_get(v_t_37_, 1);
lean_inc_ref(v_threshold_53_);
lean_dec_ref(v_t_37_);
v___x_54_ = lean_apply_2(v_k_38_, v_organism_52_, v_threshold_53_);
return v___x_54_;
}
case 12:
{
lean_object* v_organism_55_; lean_object* v_threshold_56_; lean_object* v___x_57_; 
v_organism_55_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_organism_55_);
v_threshold_56_ = lean_ctor_get(v_t_37_, 1);
lean_inc_ref(v_threshold_56_);
lean_dec_ref(v_t_37_);
v___x_57_ = lean_apply_2(v_k_38_, v_organism_55_, v_threshold_56_);
return v___x_57_;
}
case 13:
{
lean_object* v_minScore_58_; lean_object* v___x_59_; 
v_minScore_58_ = lean_ctor_get(v_t_37_, 0);
lean_inc(v_minScore_58_);
lean_dec_ref(v_t_37_);
v___x_59_ = lean_apply_1(v_k_38_, v_minScore_58_);
return v___x_59_;
}
case 14:
{
uint8_t v_isCytosolic_60_; lean_object* v_threshold_61_; lean_object* v___x_62_; lean_object* v___x_63_; 
v_isCytosolic_60_ = lean_ctor_get_uint8(v_t_37_, sizeof(void*)*1);
v_threshold_61_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_threshold_61_);
lean_dec_ref(v_t_37_);
v___x_62_ = lean_box(v_isCytosolic_60_);
v___x_63_ = lean_apply_2(v_k_38_, v___x_62_, v_threshold_61_);
return v___x_63_;
}
case 15:
{
lean_object* v_dgThreshold_64_; lean_object* v___x_65_; 
v_dgThreshold_64_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_dgThreshold_64_);
lean_dec_ref(v_t_37_);
v___x_65_ = lean_apply_1(v_k_38_, v_dgThreshold_64_);
return v___x_65_;
}
case 16:
{
lean_object* v_organism_66_; lean_object* v___x_67_; 
v_organism_66_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_organism_66_);
lean_dec_ref(v_t_37_);
v___x_67_ = lean_apply_1(v_k_38_, v_organism_66_);
return v___x_67_;
}
case 17:
{
lean_object* v_threshold_68_; lean_object* v___x_69_; 
v_threshold_68_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_threshold_68_);
lean_dec_ref(v_t_37_);
v___x_69_ = lean_apply_1(v_k_38_, v_threshold_68_);
return v___x_69_;
}
case 21:
{
lean_object* v_ddgThreshold_70_; lean_object* v___x_71_; 
v_ddgThreshold_70_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_ddgThreshold_70_);
lean_dec_ref(v_t_37_);
v___x_71_ = lean_apply_1(v_k_38_, v_ddgThreshold_70_);
return v___x_71_;
}
case 22:
{
lean_object* v_maxDDG_72_; lean_object* v___x_73_; 
v_maxDDG_72_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_maxDDG_72_);
lean_dec_ref(v_t_37_);
v___x_73_ = lean_apply_1(v_k_38_, v_maxDDG_72_);
return v___x_73_;
}
case 24:
{
lean_object* v_threshold_74_; lean_object* v___x_75_; 
v_threshold_74_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_threshold_74_);
lean_dec_ref(v_t_37_);
v___x_75_ = lean_apply_1(v_k_38_, v_threshold_74_);
return v___x_75_;
}
case 25:
{
lean_object* v_minScore_76_; lean_object* v___x_77_; 
v_minScore_76_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_minScore_76_);
lean_dec_ref(v_t_37_);
v___x_77_ = lean_apply_1(v_k_38_, v_minScore_76_);
return v___x_77_;
}
case 27:
{
lean_object* v_pILo_78_; lean_object* v_pIHi_79_; lean_object* v___x_80_; 
v_pILo_78_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_pILo_78_);
v_pIHi_79_ = lean_ctor_get(v_t_37_, 1);
lean_inc_ref(v_pIHi_79_);
lean_dec_ref(v_t_37_);
v___x_80_ = lean_apply_2(v_k_38_, v_pILo_78_, v_pIHi_79_);
return v___x_80_;
}
case 28:
{
lean_object* v_maxLen_81_; lean_object* v___x_82_; 
v_maxLen_81_ = lean_ctor_get(v_t_37_, 0);
lean_inc(v_maxLen_81_);
lean_dec_ref(v_t_37_);
v___x_82_ = lean_apply_1(v_k_38_, v_maxLen_81_);
return v___x_82_;
}
case 29:
{
lean_object* v_maxScore_83_; lean_object* v___x_84_; 
v_maxScore_83_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_maxScore_83_);
lean_dec_ref(v_t_37_);
v___x_84_ = lean_apply_1(v_k_38_, v_maxScore_83_);
return v___x_84_;
}
case 30:
{
lean_object* v_ic50Threshold_85_; lean_object* v___x_86_; 
v_ic50Threshold_85_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_ic50Threshold_85_);
lean_dec_ref(v_t_37_);
v___x_86_ = lean_apply_1(v_k_38_, v_ic50Threshold_85_);
return v___x_86_;
}
case 31:
{
lean_object* v_scoreThreshold_87_; lean_object* v___x_88_; 
v_scoreThreshold_87_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_scoreThreshold_87_);
lean_dec_ref(v_t_37_);
v___x_88_ = lean_apply_1(v_k_38_, v_scoreThreshold_87_);
return v___x_88_;
}
case 32:
{
lean_object* v_maxCoverage_89_; lean_object* v___x_90_; 
v_maxCoverage_89_ = lean_ctor_get(v_t_37_, 0);
lean_inc_ref(v_maxCoverage_89_);
lean_dec_ref(v_t_37_);
v___x_90_ = lean_apply_1(v_k_38_, v_maxCoverage_89_);
return v___x_90_;
}
default: 
{
lean_dec(v_t_37_);
return v_k_38_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorElim(lean_object* v_motive_91_, lean_object* v_ctorIdx_92_, lean_object* v_t_93_, lean_object* v_h_94_, lean_object* v_k_95_){
_start:
{
lean_object* v___x_96_; 
v___x_96_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_93_, v_k_95_);
return v___x_96_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___boxed(lean_object* v_motive_97_, lean_object* v_ctorIdx_98_, lean_object* v_t_99_, lean_object* v_h_100_, lean_object* v_k_101_){
_start:
{
lean_object* v_res_102_; 
v_res_102_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim(v_motive_97_, v_ctorIdx_98_, v_t_99_, v_h_100_, v_k_101_);
lean_dec(v_ctorIdx_98_);
return v_res_102_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SpliceCorrect_elim___redArg(lean_object* v_t_103_, lean_object* v_SpliceCorrect_104_){
_start:
{
lean_object* v___x_105_; 
v___x_105_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_103_, v_SpliceCorrect_104_);
return v___x_105_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SpliceCorrect_elim(lean_object* v_motive_106_, lean_object* v_t_107_, lean_object* v_h_108_, lean_object* v_SpliceCorrect_109_){
_start:
{
lean_object* v___x_110_; 
v___x_110_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_107_, v_SpliceCorrect_109_);
return v___x_110_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticSplice_elim___redArg(lean_object* v_t_111_, lean_object* v_NoCrypticSplice_112_){
_start:
{
lean_object* v___x_113_; 
v___x_113_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_111_, v_NoCrypticSplice_112_);
return v___x_113_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticSplice_elim(lean_object* v_motive_114_, lean_object* v_t_115_, lean_object* v_h_116_, lean_object* v_NoCrypticSplice_117_){
_start:
{
lean_object* v___x_118_; 
v___x_118_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_115_, v_NoCrypticSplice_117_);
return v___x_118_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonAdapted_elim___redArg(lean_object* v_t_119_, lean_object* v_CodonAdapted_120_){
_start:
{
lean_object* v___x_121_; 
v___x_121_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_119_, v_CodonAdapted_120_);
return v___x_121_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonAdapted_elim(lean_object* v_motive_122_, lean_object* v_t_123_, lean_object* v_h_124_, lean_object* v_CodonAdapted_125_){
_start:
{
lean_object* v___x_126_; 
v___x_126_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_123_, v_CodonAdapted_125_);
return v___x_126_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_GCInRange_elim___redArg(lean_object* v_t_127_, lean_object* v_GCInRange_128_){
_start:
{
lean_object* v___x_129_; 
v___x_129_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_127_, v_GCInRange_128_);
return v___x_129_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_GCInRange_elim(lean_object* v_motive_130_, lean_object* v_t_131_, lean_object* v_h_132_, lean_object* v_GCInRange_133_){
_start:
{
lean_object* v___x_134_; 
v___x_134_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_131_, v_GCInRange_133_);
return v___x_134_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoRestrictionSite_elim___redArg(lean_object* v_t_135_, lean_object* v_NoRestrictionSite_136_){
_start:
{
lean_object* v___x_137_; 
v___x_137_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_135_, v_NoRestrictionSite_136_);
return v___x_137_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoRestrictionSite_elim(lean_object* v_motive_138_, lean_object* v_t_139_, lean_object* v_h_140_, lean_object* v_NoRestrictionSite_141_){
_start:
{
lean_object* v___x_142_; 
v___x_142_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_139_, v_NoRestrictionSite_141_);
return v___x_142_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_InFrame_elim___redArg(lean_object* v_t_143_, lean_object* v_InFrame_144_){
_start:
{
lean_object* v___x_145_; 
v___x_145_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_143_, v_InFrame_144_);
return v___x_145_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_InFrame_elim(lean_object* v_motive_146_, lean_object* v_t_147_, lean_object* v_h_148_, lean_object* v_InFrame_149_){
_start:
{
lean_object* v___x_150_; 
v___x_150_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_147_, v_InFrame_149_);
return v___x_150_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoInstabilityMotif_elim___redArg(lean_object* v_t_151_, lean_object* v_NoInstabilityMotif_152_){
_start:
{
lean_object* v___x_153_; 
v___x_153_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_151_, v_NoInstabilityMotif_152_);
return v___x_153_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoInstabilityMotif_elim(lean_object* v_motive_154_, lean_object* v_t_155_, lean_object* v_h_156_, lean_object* v_NoInstabilityMotif_157_){
_start:
{
lean_object* v___x_158_; 
v___x_158_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_155_, v_NoInstabilityMotif_157_);
return v___x_158_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCpGIsland_elim___redArg(lean_object* v_t_159_, lean_object* v_NoCpGIsland_160_){
_start:
{
lean_object* v___x_161_; 
v___x_161_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_159_, v_NoCpGIsland_160_);
return v___x_161_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCpGIsland_elim(lean_object* v_motive_162_, lean_object* v_t_163_, lean_object* v_h_164_, lean_object* v_NoCpGIsland_165_){
_start:
{
lean_object* v___x_166_; 
v___x_166_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_163_, v_NoCpGIsland_165_);
return v___x_166_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoGTDinucleotide_elim___redArg(lean_object* v_t_167_, lean_object* v_NoGTDinucleotide_168_){
_start:
{
lean_object* v___x_169_; 
v___x_169_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_167_, v_NoGTDinucleotide_168_);
return v___x_169_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoGTDinucleotide_elim(lean_object* v_motive_170_, lean_object* v_t_171_, lean_object* v_h_172_, lean_object* v_NoGTDinucleotide_173_){
_start:
{
lean_object* v___x_174_; 
v___x_174_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_171_, v_NoGTDinucleotide_173_);
return v___x_174_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoStopCodons_elim___redArg(lean_object* v_t_175_, lean_object* v_NoStopCodons_176_){
_start:
{
lean_object* v___x_177_; 
v___x_177_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_175_, v_NoStopCodons_176_);
return v___x_177_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoStopCodons_elim(lean_object* v_motive_178_, lean_object* v_t_179_, lean_object* v_h_180_, lean_object* v_NoStopCodons_181_){
_start:
{
lean_object* v___x_182_; 
v___x_182_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_179_, v_NoStopCodons_181_);
return v___x_182_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ValidCodingSeq_elim___redArg(lean_object* v_t_183_, lean_object* v_ValidCodingSeq_184_){
_start:
{
lean_object* v___x_185_; 
v___x_185_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_183_, v_ValidCodingSeq_184_);
return v___x_185_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ValidCodingSeq_elim(lean_object* v_motive_186_, lean_object* v_t_187_, lean_object* v_h_188_, lean_object* v_ValidCodingSeq_189_){
_start:
{
lean_object* v___x_190_; 
v___x_190_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_187_, v_ValidCodingSeq_189_);
return v___x_190_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonOptimality_elim___redArg(lean_object* v_t_191_, lean_object* v_CodonOptimality_192_){
_start:
{
lean_object* v___x_193_; 
v___x_193_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_191_, v_CodonOptimality_192_);
return v___x_193_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonOptimality_elim(lean_object* v_motive_194_, lean_object* v_t_195_, lean_object* v_h_196_, lean_object* v_CodonOptimality_197_){
_start:
{
lean_object* v___x_198_; 
v___x_198_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_195_, v_CodonOptimality_197_);
return v___x_198_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticPromoter_elim___redArg(lean_object* v_t_199_, lean_object* v_NoCrypticPromoter_200_){
_start:
{
lean_object* v___x_201_; 
v___x_201_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_199_, v_NoCrypticPromoter_200_);
return v___x_201_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticPromoter_elim(lean_object* v_motive_202_, lean_object* v_t_203_, lean_object* v_h_204_, lean_object* v_NoCrypticPromoter_205_){
_start:
{
lean_object* v___x_206_; 
v___x_206_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_203_, v_NoCrypticPromoter_205_);
return v___x_206_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ConservationScore_elim___redArg(lean_object* v_t_207_, lean_object* v_ConservationScore_208_){
_start:
{
lean_object* v___x_209_; 
v___x_209_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_207_, v_ConservationScore_208_);
return v___x_209_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ConservationScore_elim(lean_object* v_motive_210_, lean_object* v_t_211_, lean_object* v_h_212_, lean_object* v_ConservationScore_213_){
_start:
{
lean_object* v___x_214_; 
v___x_214_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_211_, v_ConservationScore_213_);
return v___x_214_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoUnexpectedTMDomain_elim___redArg(lean_object* v_t_215_, lean_object* v_NoUnexpectedTMDomain_216_){
_start:
{
lean_object* v___x_217_; 
v___x_217_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_215_, v_NoUnexpectedTMDomain_216_);
return v___x_217_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoUnexpectedTMDomain_elim(lean_object* v_motive_218_, lean_object* v_t_219_, lean_object* v_h_220_, lean_object* v_NoUnexpectedTMDomain_221_){
_start:
{
lean_object* v___x_222_; 
v___x_222_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_219_, v_NoUnexpectedTMDomain_221_);
return v___x_222_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_mRNASecondaryStructure_elim___redArg(lean_object* v_t_223_, lean_object* v_mRNASecondaryStructure_224_){
_start:
{
lean_object* v___x_225_; 
v___x_225_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_223_, v_mRNASecondaryStructure_224_);
return v___x_225_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_mRNASecondaryStructure_elim(lean_object* v_motive_226_, lean_object* v_t_227_, lean_object* v_h_228_, lean_object* v_mRNASecondaryStructure_229_){
_start:
{
lean_object* v___x_230_; 
v___x_230_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_227_, v_mRNASecondaryStructure_229_);
return v___x_230_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CoTranslationalFolding_elim___redArg(lean_object* v_t_231_, lean_object* v_CoTranslationalFolding_232_){
_start:
{
lean_object* v___x_233_; 
v___x_233_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_231_, v_CoTranslationalFolding_232_);
return v___x_233_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CoTranslationalFolding_elim(lean_object* v_motive_234_, lean_object* v_t_235_, lean_object* v_h_236_, lean_object* v_CoTranslationalFolding_237_){
_start:
{
lean_object* v___x_238_; 
v___x_238_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_235_, v_CoTranslationalFolding_237_);
return v___x_238_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_StructureConfidence_elim___redArg(lean_object* v_t_239_, lean_object* v_StructureConfidence_240_){
_start:
{
lean_object* v___x_241_; 
v___x_241_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_239_, v_StructureConfidence_240_);
return v___x_241_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_StructureConfidence_elim(lean_object* v_motive_242_, lean_object* v_t_243_, lean_object* v_h_244_, lean_object* v_StructureConfidence_245_){
_start:
{
lean_object* v___x_246_; 
v___x_246_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_243_, v_StructureConfidence_245_);
return v___x_246_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoMisfoldingRisk_elim___redArg(lean_object* v_t_247_, lean_object* v_NoMisfoldingRisk_248_){
_start:
{
lean_object* v___x_249_; 
v___x_249_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_247_, v_NoMisfoldingRisk_248_);
return v___x_249_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoMisfoldingRisk_elim(lean_object* v_motive_250_, lean_object* v_t_251_, lean_object* v_h_252_, lean_object* v_NoMisfoldingRisk_253_){
_start:
{
lean_object* v___x_254_; 
v___x_254_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_251_, v_NoMisfoldingRisk_253_);
return v___x_254_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CorrectFoldTopology_elim___redArg(lean_object* v_t_255_, lean_object* v_CorrectFoldTopology_256_){
_start:
{
lean_object* v___x_257_; 
v___x_257_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_255_, v_CorrectFoldTopology_256_);
return v___x_257_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CorrectFoldTopology_elim(lean_object* v_motive_258_, lean_object* v_t_259_, lean_object* v_h_260_, lean_object* v_CorrectFoldTopology_261_){
_start:
{
lean_object* v___x_262_; 
v___x_262_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_259_, v_CorrectFoldTopology_261_);
return v___x_262_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoUnexpectedInteraction_elim___redArg(lean_object* v_t_263_, lean_object* v_NoUnexpectedInteraction_264_){
_start:
{
lean_object* v___x_265_; 
v___x_265_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_263_, v_NoUnexpectedInteraction_264_);
return v___x_265_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoUnexpectedInteraction_elim(lean_object* v_motive_266_, lean_object* v_t_267_, lean_object* v_h_268_, lean_object* v_NoUnexpectedInteraction_269_){
_start:
{
lean_object* v___x_270_; 
v___x_270_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_267_, v_NoUnexpectedInteraction_269_);
return v___x_270_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_StableFolding_elim___redArg(lean_object* v_t_271_, lean_object* v_StableFolding_272_){
_start:
{
lean_object* v___x_273_; 
v___x_273_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_271_, v_StableFolding_272_);
return v___x_273_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_StableFolding_elim(lean_object* v_motive_274_, lean_object* v_t_275_, lean_object* v_h_276_, lean_object* v_StableFolding_277_){
_start:
{
lean_object* v___x_278_; 
v___x_278_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_275_, v_StableFolding_277_);
return v___x_278_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoDestabilizingMutation_elim___redArg(lean_object* v_t_279_, lean_object* v_NoDestabilizingMutation_280_){
_start:
{
lean_object* v___x_281_; 
v___x_281_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_279_, v_NoDestabilizingMutation_280_);
return v___x_281_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoDestabilizingMutation_elim(lean_object* v_motive_282_, lean_object* v_t_283_, lean_object* v_h_284_, lean_object* v_NoDestabilizingMutation_285_){
_start:
{
lean_object* v___x_286_; 
v___x_286_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_283_, v_NoDestabilizingMutation_285_);
return v___x_286_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_DisulfideBondIntegrity_elim___redArg(lean_object* v_t_287_, lean_object* v_DisulfideBondIntegrity_288_){
_start:
{
lean_object* v___x_289_; 
v___x_289_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_287_, v_DisulfideBondIntegrity_288_);
return v___x_289_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_DisulfideBondIntegrity_elim(lean_object* v_motive_290_, lean_object* v_t_291_, lean_object* v_h_292_, lean_object* v_DisulfideBondIntegrity_293_){
_start:
{
lean_object* v___x_294_; 
v___x_294_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_291_, v_DisulfideBondIntegrity_293_);
return v___x_294_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_HydrophobicCoreQuality_elim___redArg(lean_object* v_t_295_, lean_object* v_HydrophobicCoreQuality_296_){
_start:
{
lean_object* v___x_297_; 
v___x_297_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_295_, v_HydrophobicCoreQuality_296_);
return v___x_297_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_HydrophobicCoreQuality_elim(lean_object* v_motive_298_, lean_object* v_t_299_, lean_object* v_h_300_, lean_object* v_HydrophobicCoreQuality_301_){
_start:
{
lean_object* v___x_302_; 
v___x_302_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_299_, v_HydrophobicCoreQuality_301_);
return v___x_302_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SolubleExpression_elim___redArg(lean_object* v_t_303_, lean_object* v_SolubleExpression_304_){
_start:
{
lean_object* v___x_305_; 
v___x_305_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_303_, v_SolubleExpression_304_);
return v___x_305_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SolubleExpression_elim(lean_object* v_motive_306_, lean_object* v_t_307_, lean_object* v_h_308_, lean_object* v_SolubleExpression_309_){
_start:
{
lean_object* v___x_310_; 
v___x_310_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_307_, v_SolubleExpression_309_);
return v___x_310_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoAggregationProneRegion_elim___redArg(lean_object* v_t_311_, lean_object* v_NoAggregationProneRegion_312_){
_start:
{
lean_object* v___x_313_; 
v___x_313_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_311_, v_NoAggregationProneRegion_312_);
return v___x_313_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoAggregationProneRegion_elim(lean_object* v_motive_314_, lean_object* v_t_315_, lean_object* v_h_316_, lean_object* v_NoAggregationProneRegion_317_){
_start:
{
lean_object* v___x_318_; 
v___x_318_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_315_, v_NoAggregationProneRegion_317_);
return v___x_318_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ChargeComposition_elim___redArg(lean_object* v_t_319_, lean_object* v_ChargeComposition_320_){
_start:
{
lean_object* v___x_321_; 
v___x_321_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_319_, v_ChargeComposition_320_);
return v___x_321_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ChargeComposition_elim(lean_object* v_motive_322_, lean_object* v_t_323_, lean_object* v_h_324_, lean_object* v_ChargeComposition_325_){
_start:
{
lean_object* v___x_326_; 
v___x_326_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_323_, v_ChargeComposition_325_);
return v___x_326_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoLongHydrophobicStretch_elim___redArg(lean_object* v_t_327_, lean_object* v_NoLongHydrophobicStretch_328_){
_start:
{
lean_object* v___x_329_; 
v___x_329_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_327_, v_NoLongHydrophobicStretch_328_);
return v___x_329_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoLongHydrophobicStretch_elim(lean_object* v_motive_330_, lean_object* v_t_331_, lean_object* v_h_332_, lean_object* v_NoLongHydrophobicStretch_333_){
_start:
{
lean_object* v___x_334_; 
v___x_334_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_331_, v_NoLongHydrophobicStretch_333_);
return v___x_334_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_LowImmunogenicity_elim___redArg(lean_object* v_t_335_, lean_object* v_LowImmunogenicity_336_){
_start:
{
lean_object* v___x_337_; 
v___x_337_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_335_, v_LowImmunogenicity_336_);
return v___x_337_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_LowImmunogenicity_elim(lean_object* v_motive_338_, lean_object* v_t_339_, lean_object* v_h_340_, lean_object* v_LowImmunogenicity_341_){
_start:
{
lean_object* v___x_342_; 
v___x_342_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_339_, v_LowImmunogenicity_341_);
return v___x_342_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoStrongTCellEpitope_elim___redArg(lean_object* v_t_343_, lean_object* v_NoStrongTCellEpitope_344_){
_start:
{
lean_object* v___x_345_; 
v___x_345_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_343_, v_NoStrongTCellEpitope_344_);
return v___x_345_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoStrongTCellEpitope_elim(lean_object* v_motive_346_, lean_object* v_t_347_, lean_object* v_h_348_, lean_object* v_NoStrongTCellEpitope_349_){
_start:
{
lean_object* v___x_350_; 
v___x_350_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_347_, v_NoStrongTCellEpitope_349_);
return v___x_350_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoDominantBCellEpitope_elim___redArg(lean_object* v_t_351_, lean_object* v_NoDominantBCellEpitope_352_){
_start:
{
lean_object* v___x_353_; 
v___x_353_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_351_, v_NoDominantBCellEpitope_352_);
return v___x_353_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoDominantBCellEpitope_elim(lean_object* v_motive_354_, lean_object* v_t_355_, lean_object* v_h_356_, lean_object* v_NoDominantBCellEpitope_357_){
_start:
{
lean_object* v___x_358_; 
v___x_358_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_355_, v_NoDominantBCellEpitope_357_);
return v___x_358_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_PopulationCoverageSafe_elim___redArg(lean_object* v_t_359_, lean_object* v_PopulationCoverageSafe_360_){
_start:
{
lean_object* v___x_361_; 
v___x_361_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_359_, v_PopulationCoverageSafe_360_);
return v___x_361_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_PopulationCoverageSafe_elim(lean_object* v_motive_362_, lean_object* v_t_363_, lean_object* v_h_364_, lean_object* v_PopulationCoverageSafe_365_){
_start:
{
lean_object* v___x_366_; 
v___x_366_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_363_, v_PopulationCoverageSafe_365_);
return v___x_366_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1_spec__3(lean_object* v_x_367_, lean_object* v_x_368_, lean_object* v_x_369_){
_start:
{
if (lean_obj_tag(v_x_369_) == 0)
{
lean_dec(v_x_367_);
return v_x_368_;
}
else
{
lean_object* v_head_370_; lean_object* v_tail_371_; lean_object* v___x_373_; uint8_t v_isShared_374_; uint8_t v_isSharedCheck_381_; 
v_head_370_ = lean_ctor_get(v_x_369_, 0);
v_tail_371_ = lean_ctor_get(v_x_369_, 1);
v_isSharedCheck_381_ = !lean_is_exclusive(v_x_369_);
if (v_isSharedCheck_381_ == 0)
{
v___x_373_ = v_x_369_;
v_isShared_374_ = v_isSharedCheck_381_;
goto v_resetjp_372_;
}
else
{
lean_inc(v_tail_371_);
lean_inc(v_head_370_);
lean_dec(v_x_369_);
v___x_373_ = lean_box(0);
v_isShared_374_ = v_isSharedCheck_381_;
goto v_resetjp_372_;
}
v_resetjp_372_:
{
lean_object* v___x_376_; 
lean_inc(v_x_367_);
if (v_isShared_374_ == 0)
{
lean_ctor_set_tag(v___x_373_, 5);
lean_ctor_set(v___x_373_, 1, v_x_367_);
lean_ctor_set(v___x_373_, 0, v_x_368_);
v___x_376_ = v___x_373_;
goto v_reusejp_375_;
}
else
{
lean_object* v_reuseFailAlloc_380_; 
v_reuseFailAlloc_380_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_380_, 0, v_x_368_);
lean_ctor_set(v_reuseFailAlloc_380_, 1, v_x_367_);
v___x_376_ = v_reuseFailAlloc_380_;
goto v_reusejp_375_;
}
v_reusejp_375_:
{
lean_object* v___x_377_; lean_object* v___x_378_; 
v___x_377_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_head_370_);
v___x_378_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_378_, 0, v___x_376_);
lean_ctor_set(v___x_378_, 1, v___x_377_);
v_x_368_ = v___x_378_;
v_x_369_ = v_tail_371_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1(lean_object* v_x_382_, lean_object* v_x_383_, lean_object* v_x_384_){
_start:
{
if (lean_obj_tag(v_x_384_) == 0)
{
lean_dec(v_x_382_);
return v_x_383_;
}
else
{
lean_object* v_head_385_; lean_object* v_tail_386_; lean_object* v___x_388_; uint8_t v_isShared_389_; uint8_t v_isSharedCheck_396_; 
v_head_385_ = lean_ctor_get(v_x_384_, 0);
v_tail_386_ = lean_ctor_get(v_x_384_, 1);
v_isSharedCheck_396_ = !lean_is_exclusive(v_x_384_);
if (v_isSharedCheck_396_ == 0)
{
v___x_388_ = v_x_384_;
v_isShared_389_ = v_isSharedCheck_396_;
goto v_resetjp_387_;
}
else
{
lean_inc(v_tail_386_);
lean_inc(v_head_385_);
lean_dec(v_x_384_);
v___x_388_ = lean_box(0);
v_isShared_389_ = v_isSharedCheck_396_;
goto v_resetjp_387_;
}
v_resetjp_387_:
{
lean_object* v___x_391_; 
lean_inc(v_x_382_);
if (v_isShared_389_ == 0)
{
lean_ctor_set_tag(v___x_388_, 5);
lean_ctor_set(v___x_388_, 1, v_x_382_);
lean_ctor_set(v___x_388_, 0, v_x_383_);
v___x_391_ = v___x_388_;
goto v_reusejp_390_;
}
else
{
lean_object* v_reuseFailAlloc_395_; 
v_reuseFailAlloc_395_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_395_, 0, v_x_383_);
lean_ctor_set(v_reuseFailAlloc_395_, 1, v_x_382_);
v___x_391_ = v_reuseFailAlloc_395_;
goto v_reusejp_390_;
}
v_reusejp_390_:
{
lean_object* v___x_392_; lean_object* v___x_393_; lean_object* v___x_394_; 
v___x_392_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_head_385_);
v___x_393_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_393_, 0, v___x_391_);
lean_ctor_set(v___x_393_, 1, v___x_392_);
v___x_394_ = lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1_spec__3(v_x_382_, v___x_393_, v_tail_386_);
return v___x_394_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0(lean_object* v_x_397_, lean_object* v_x_398_){
_start:
{
if (lean_obj_tag(v_x_397_) == 0)
{
lean_object* v___x_399_; 
lean_dec(v_x_398_);
v___x_399_ = lean_box(0);
return v___x_399_;
}
else
{
lean_object* v_tail_400_; 
v_tail_400_ = lean_ctor_get(v_x_397_, 1);
if (lean_obj_tag(v_tail_400_) == 0)
{
lean_object* v_head_401_; lean_object* v___x_402_; 
lean_dec(v_x_398_);
v_head_401_ = lean_ctor_get(v_x_397_, 0);
lean_inc(v_head_401_);
lean_dec_ref(v_x_397_);
v___x_402_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_head_401_);
return v___x_402_;
}
else
{
lean_object* v_head_403_; lean_object* v___x_404_; lean_object* v___x_405_; 
lean_inc(v_tail_400_);
v_head_403_ = lean_ctor_get(v_x_397_, 0);
lean_inc(v_head_403_);
lean_dec_ref(v_x_397_);
v___x_404_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_head_403_);
v___x_405_ = lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1(v_x_398_, v___x_404_, v_tail_400_);
return v___x_405_;
}
}
}
}
static lean_object* _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7(void){
_start:
{
lean_object* v___x_417_; lean_object* v___x_418_; 
v___x_417_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__2));
v___x_418_ = lean_string_length(v___x_417_);
return v___x_418_;
}
}
static lean_object* _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8(void){
_start:
{
lean_object* v___x_419_; lean_object* v___x_420_; 
v___x_419_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7, &lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7_once, _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7);
v___x_420_ = lean_nat_to_int(v___x_419_);
return v___x_420_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg(lean_object* v_a_425_){
_start:
{
if (lean_obj_tag(v_a_425_) == 0)
{
lean_object* v___x_426_; 
v___x_426_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__1));
return v___x_426_;
}
else
{
lean_object* v___x_427_; lean_object* v___x_428_; lean_object* v___x_429_; lean_object* v___x_430_; lean_object* v___x_431_; lean_object* v___x_432_; lean_object* v___x_433_; lean_object* v___x_434_; uint8_t v___x_435_; lean_object* v___x_436_; 
v___x_427_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__5));
v___x_428_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0(v_a_425_, v___x_427_);
v___x_429_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8, &lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8_once, _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8);
v___x_430_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__9));
v___x_431_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_431_, 0, v___x_430_);
lean_ctor_set(v___x_431_, 1, v___x_428_);
v___x_432_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__10));
v___x_433_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_433_, 0, v___x_431_);
lean_ctor_set(v___x_433_, 1, v___x_432_);
v___x_434_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_434_, 0, v___x_429_);
lean_ctor_set(v___x_434_, 1, v___x_433_);
v___x_435_ = 0;
v___x_436_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_436_, 0, v___x_434_);
lean_ctor_set_uint8(v___x_436_, sizeof(void*)*1, v___x_435_);
return v___x_436_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2___lam__0(lean_object* v___y_437_){
_start:
{
lean_object* v___x_438_; lean_object* v___x_439_; 
v___x_438_ = l_Nat_reprFast(v___y_437_);
v___x_439_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_439_, 0, v___x_438_);
return v___x_439_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4_spec__6(lean_object* v_x_440_, lean_object* v_x_441_, lean_object* v_x_442_){
_start:
{
if (lean_obj_tag(v_x_442_) == 0)
{
lean_dec(v_x_440_);
return v_x_441_;
}
else
{
lean_object* v_head_443_; lean_object* v_tail_444_; lean_object* v___x_446_; uint8_t v_isShared_447_; uint8_t v_isSharedCheck_455_; 
v_head_443_ = lean_ctor_get(v_x_442_, 0);
v_tail_444_ = lean_ctor_get(v_x_442_, 1);
v_isSharedCheck_455_ = !lean_is_exclusive(v_x_442_);
if (v_isSharedCheck_455_ == 0)
{
v___x_446_ = v_x_442_;
v_isShared_447_ = v_isSharedCheck_455_;
goto v_resetjp_445_;
}
else
{
lean_inc(v_tail_444_);
lean_inc(v_head_443_);
lean_dec(v_x_442_);
v___x_446_ = lean_box(0);
v_isShared_447_ = v_isSharedCheck_455_;
goto v_resetjp_445_;
}
v_resetjp_445_:
{
lean_object* v___x_449_; 
lean_inc(v_x_440_);
if (v_isShared_447_ == 0)
{
lean_ctor_set_tag(v___x_446_, 5);
lean_ctor_set(v___x_446_, 1, v_x_440_);
lean_ctor_set(v___x_446_, 0, v_x_441_);
v___x_449_ = v___x_446_;
goto v_reusejp_448_;
}
else
{
lean_object* v_reuseFailAlloc_454_; 
v_reuseFailAlloc_454_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_454_, 0, v_x_441_);
lean_ctor_set(v_reuseFailAlloc_454_, 1, v_x_440_);
v___x_449_ = v_reuseFailAlloc_454_;
goto v_reusejp_448_;
}
v_reusejp_448_:
{
lean_object* v___x_450_; lean_object* v___x_451_; lean_object* v___x_452_; 
v___x_450_ = l_Nat_reprFast(v_head_443_);
v___x_451_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_451_, 0, v___x_450_);
v___x_452_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_452_, 0, v___x_449_);
lean_ctor_set(v___x_452_, 1, v___x_451_);
v_x_441_ = v___x_452_;
v_x_442_ = v_tail_444_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4(lean_object* v_x_456_, lean_object* v_x_457_, lean_object* v_x_458_){
_start:
{
if (lean_obj_tag(v_x_458_) == 0)
{
lean_dec(v_x_456_);
return v_x_457_;
}
else
{
lean_object* v_head_459_; lean_object* v_tail_460_; lean_object* v___x_462_; uint8_t v_isShared_463_; uint8_t v_isSharedCheck_471_; 
v_head_459_ = lean_ctor_get(v_x_458_, 0);
v_tail_460_ = lean_ctor_get(v_x_458_, 1);
v_isSharedCheck_471_ = !lean_is_exclusive(v_x_458_);
if (v_isSharedCheck_471_ == 0)
{
v___x_462_ = v_x_458_;
v_isShared_463_ = v_isSharedCheck_471_;
goto v_resetjp_461_;
}
else
{
lean_inc(v_tail_460_);
lean_inc(v_head_459_);
lean_dec(v_x_458_);
v___x_462_ = lean_box(0);
v_isShared_463_ = v_isSharedCheck_471_;
goto v_resetjp_461_;
}
v_resetjp_461_:
{
lean_object* v___x_465_; 
lean_inc(v_x_456_);
if (v_isShared_463_ == 0)
{
lean_ctor_set_tag(v___x_462_, 5);
lean_ctor_set(v___x_462_, 1, v_x_456_);
lean_ctor_set(v___x_462_, 0, v_x_457_);
v___x_465_ = v___x_462_;
goto v_reusejp_464_;
}
else
{
lean_object* v_reuseFailAlloc_470_; 
v_reuseFailAlloc_470_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_470_, 0, v_x_457_);
lean_ctor_set(v_reuseFailAlloc_470_, 1, v_x_456_);
v___x_465_ = v_reuseFailAlloc_470_;
goto v_reusejp_464_;
}
v_reusejp_464_:
{
lean_object* v___x_466_; lean_object* v___x_467_; lean_object* v___x_468_; lean_object* v___x_469_; 
v___x_466_ = l_Nat_reprFast(v_head_459_);
v___x_467_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_467_, 0, v___x_466_);
v___x_468_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_468_, 0, v___x_465_);
lean_ctor_set(v___x_468_, 1, v___x_467_);
v___x_469_ = lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4_spec__6(v_x_456_, v___x_468_, v_tail_460_);
return v___x_469_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2(lean_object* v_x_472_, lean_object* v_x_473_){
_start:
{
if (lean_obj_tag(v_x_472_) == 0)
{
lean_object* v___x_474_; 
lean_dec(v_x_473_);
v___x_474_ = lean_box(0);
return v___x_474_;
}
else
{
lean_object* v_tail_475_; 
v_tail_475_ = lean_ctor_get(v_x_472_, 1);
if (lean_obj_tag(v_tail_475_) == 0)
{
lean_object* v_head_476_; lean_object* v___x_477_; 
lean_dec(v_x_473_);
v_head_476_ = lean_ctor_get(v_x_472_, 0);
lean_inc(v_head_476_);
lean_dec_ref(v_x_472_);
v___x_477_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2___lam__0(v_head_476_);
return v___x_477_;
}
else
{
lean_object* v_head_478_; lean_object* v___x_479_; lean_object* v___x_480_; 
lean_inc(v_tail_475_);
v_head_478_ = lean_ctor_get(v_x_472_, 0);
lean_inc(v_head_478_);
lean_dec_ref(v_x_472_);
v___x_479_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2___lam__0(v_head_478_);
v___x_480_ = lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4(v_x_473_, v___x_479_, v_tail_475_);
return v___x_480_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___redArg(lean_object* v_a_481_){
_start:
{
if (lean_obj_tag(v_a_481_) == 0)
{
lean_object* v___x_482_; 
v___x_482_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__1));
return v___x_482_;
}
else
{
lean_object* v___x_483_; lean_object* v___x_484_; lean_object* v___x_485_; lean_object* v___x_486_; lean_object* v___x_487_; lean_object* v___x_488_; lean_object* v___x_489_; lean_object* v___x_490_; lean_object* v___x_491_; 
v___x_483_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__5));
v___x_484_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2(v_a_481_, v___x_483_);
v___x_485_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8, &lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8_once, _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8);
v___x_486_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__9));
v___x_487_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_487_, 0, v___x_486_);
lean_ctor_set(v___x_487_, 1, v___x_484_);
v___x_488_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__10));
v___x_489_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_489_, 0, v___x_487_);
lean_ctor_set(v___x_489_, 1, v___x_488_);
v___x_490_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_490_, 0, v___x_485_);
lean_ctor_set(v___x_490_, 1, v___x_489_);
v___x_491_ = l_Std_Format_fill(v___x_490_);
return v___x_491_;
}
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25(void){
_start:
{
lean_object* v___x_531_; lean_object* v___x_532_; 
v___x_531_ = lean_unsigned_to_nat(2u);
v___x_532_ = lean_nat_to_int(v___x_531_);
return v___x_532_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26(void){
_start:
{
lean_object* v___x_533_; lean_object* v___x_534_; 
v___x_533_ = lean_unsigned_to_nat(1u);
v___x_534_ = lean_nat_to_int(v___x_533_);
return v___x_534_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32(void){
_start:
{
lean_object* v___x_543_; lean_object* v___x_544_; 
v___x_543_ = lean_unsigned_to_nat(0u);
v___x_544_ = lean_nat_to_int(v___x_543_);
return v___x_544_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr(lean_object* v_x_665_, lean_object* v_prec_666_){
_start:
{
lean_object* v___y_668_; lean_object* v___y_669_; lean_object* v___y_670_; lean_object* v___y_677_; lean_object* v___y_678_; lean_object* v___y_679_; lean_object* v___y_686_; lean_object* v___y_687_; lean_object* v___y_688_; lean_object* v___y_695_; lean_object* v___y_702_; lean_object* v___y_703_; lean_object* v___y_704_; lean_object* v___y_711_; lean_object* v___y_718_; lean_object* v___y_719_; lean_object* v___y_720_; lean_object* v___y_727_; lean_object* v___y_734_; lean_object* v___y_741_; lean_object* v___y_748_; lean_object* v___y_749_; lean_object* v___y_750_; lean_object* v___y_757_; lean_object* v___y_758_; lean_object* v___y_759_; lean_object* v___y_766_; lean_object* v___y_767_; lean_object* v___y_768_; lean_object* v___y_775_; lean_object* v___y_782_; lean_object* v___y_783_; lean_object* v___y_784_; lean_object* v___y_791_; lean_object* v___y_798_; lean_object* v___y_805_; lean_object* v___y_812_; lean_object* v___y_819_; lean_object* v___y_826_; lean_object* v___y_827_; lean_object* v___y_828_; lean_object* v___y_835_; lean_object* v___y_836_; lean_object* v___y_837_; lean_object* v___y_844_; lean_object* v___y_845_; lean_object* v___y_846_; lean_object* v___y_853_; lean_object* v___y_854_; lean_object* v___y_855_; lean_object* v___y_862_; lean_object* v___y_863_; lean_object* v___y_864_; lean_object* v___y_871_; lean_object* v___y_872_; lean_object* v___y_873_; lean_object* v___y_880_; lean_object* v___y_881_; lean_object* v___y_882_; lean_object* v___y_889_; lean_object* v___y_890_; lean_object* v___y_891_; 
switch(lean_obj_tag(v_x_665_))
{
case 0:
{
lean_object* v_cellType_897_; lean_object* v___x_899_; uint8_t v_isShared_900_; uint8_t v_isSharedCheck_917_; 
v_cellType_897_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_917_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_917_ == 0)
{
v___x_899_ = v_x_665_;
v_isShared_900_ = v_isSharedCheck_917_;
goto v_resetjp_898_;
}
else
{
lean_inc(v_cellType_897_);
lean_dec(v_x_665_);
v___x_899_ = lean_box(0);
v_isShared_900_ = v_isSharedCheck_917_;
goto v_resetjp_898_;
}
v_resetjp_898_:
{
lean_object* v___y_902_; lean_object* v___x_913_; uint8_t v___x_914_; 
v___x_913_ = lean_unsigned_to_nat(1024u);
v___x_914_ = lean_nat_dec_le(v___x_913_, v_prec_666_);
if (v___x_914_ == 0)
{
lean_object* v___x_915_; 
v___x_915_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_902_ = v___x_915_;
goto v___jp_901_;
}
else
{
lean_object* v___x_916_; 
v___x_916_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_902_ = v___x_916_;
goto v___jp_901_;
}
v___jp_901_:
{
lean_object* v___x_903_; lean_object* v___x_904_; lean_object* v___x_906_; 
v___x_903_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__24));
v___x_904_ = l_String_quote(v_cellType_897_);
if (v_isShared_900_ == 0)
{
lean_ctor_set_tag(v___x_899_, 3);
lean_ctor_set(v___x_899_, 0, v___x_904_);
v___x_906_ = v___x_899_;
goto v_reusejp_905_;
}
else
{
lean_object* v_reuseFailAlloc_912_; 
v_reuseFailAlloc_912_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_912_, 0, v___x_904_);
v___x_906_ = v_reuseFailAlloc_912_;
goto v_reusejp_905_;
}
v_reusejp_905_:
{
lean_object* v___x_907_; lean_object* v___x_908_; uint8_t v___x_909_; lean_object* v___x_910_; lean_object* v___x_911_; 
v___x_907_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_907_, 0, v___x_903_);
lean_ctor_set(v___x_907_, 1, v___x_906_);
lean_inc(v___y_902_);
v___x_908_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_908_, 0, v___y_902_);
lean_ctor_set(v___x_908_, 1, v___x_907_);
v___x_909_ = 0;
v___x_910_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_910_, 0, v___x_908_);
lean_ctor_set_uint8(v___x_910_, sizeof(void*)*1, v___x_909_);
v___x_911_ = l_Repr_addAppParen(v___x_910_, v_prec_666_);
return v___x_911_;
}
}
}
}
case 1:
{
lean_object* v___x_918_; uint8_t v___x_919_; 
v___x_918_ = lean_unsigned_to_nat(1024u);
v___x_919_ = lean_nat_dec_le(v___x_918_, v_prec_666_);
if (v___x_919_ == 0)
{
lean_object* v___x_920_; 
v___x_920_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_775_ = v___x_920_;
goto v___jp_774_;
}
else
{
lean_object* v___x_921_; 
v___x_921_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_775_ = v___x_921_;
goto v___jp_774_;
}
}
case 2:
{
lean_object* v_organism_922_; lean_object* v_threshold_923_; lean_object* v___x_925_; uint8_t v_isShared_926_; uint8_t v_isSharedCheck_967_; 
v_organism_922_ = lean_ctor_get(v_x_665_, 0);
v_threshold_923_ = lean_ctor_get(v_x_665_, 1);
v_isSharedCheck_967_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_967_ == 0)
{
v___x_925_ = v_x_665_;
v_isShared_926_ = v_isSharedCheck_967_;
goto v_resetjp_924_;
}
else
{
lean_inc(v_threshold_923_);
lean_inc(v_organism_922_);
lean_dec(v_x_665_);
v___x_925_ = lean_box(0);
v_isShared_926_ = v_isSharedCheck_967_;
goto v_resetjp_924_;
}
v_resetjp_924_:
{
lean_object* v___y_928_; lean_object* v___x_963_; uint8_t v___x_964_; 
v___x_963_ = lean_unsigned_to_nat(1024u);
v___x_964_ = lean_nat_dec_le(v___x_963_, v_prec_666_);
if (v___x_964_ == 0)
{
lean_object* v___x_965_; 
v___x_965_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_928_ = v___x_965_;
goto v___jp_927_;
}
else
{
lean_object* v___x_966_; 
v___x_966_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_928_ = v___x_966_;
goto v___jp_927_;
}
v___jp_927_:
{
lean_object* v_num_929_; lean_object* v_den_930_; lean_object* v___x_932_; uint8_t v_isShared_933_; uint8_t v_isSharedCheck_962_; 
v_num_929_ = lean_ctor_get(v_threshold_923_, 0);
v_den_930_ = lean_ctor_get(v_threshold_923_, 1);
v_isSharedCheck_962_ = !lean_is_exclusive(v_threshold_923_);
if (v_isSharedCheck_962_ == 0)
{
v___x_932_ = v_threshold_923_;
v_isShared_933_ = v_isSharedCheck_962_;
goto v_resetjp_931_;
}
else
{
lean_inc(v_den_930_);
lean_inc(v_num_929_);
lean_dec(v_threshold_923_);
v___x_932_ = lean_box(0);
v_isShared_933_ = v_isSharedCheck_962_;
goto v_resetjp_931_;
}
v_resetjp_931_:
{
lean_object* v___x_934_; lean_object* v___x_935_; lean_object* v___x_936_; lean_object* v___x_937_; lean_object* v___x_939_; 
v___x_934_ = lean_box(1);
v___x_935_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__29));
v___x_936_ = l_String_quote(v_organism_922_);
v___x_937_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_937_, 0, v___x_936_);
if (v_isShared_933_ == 0)
{
lean_ctor_set_tag(v___x_932_, 5);
lean_ctor_set(v___x_932_, 1, v___x_937_);
lean_ctor_set(v___x_932_, 0, v___x_935_);
v___x_939_ = v___x_932_;
goto v_reusejp_938_;
}
else
{
lean_object* v_reuseFailAlloc_961_; 
v_reuseFailAlloc_961_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_961_, 0, v___x_935_);
lean_ctor_set(v_reuseFailAlloc_961_, 1, v___x_937_);
v___x_939_ = v_reuseFailAlloc_961_;
goto v_reusejp_938_;
}
v_reusejp_938_:
{
lean_object* v___x_941_; 
if (v_isShared_926_ == 0)
{
lean_ctor_set_tag(v___x_925_, 5);
lean_ctor_set(v___x_925_, 1, v___x_934_);
lean_ctor_set(v___x_925_, 0, v___x_939_);
v___x_941_ = v___x_925_;
goto v_reusejp_940_;
}
else
{
lean_object* v_reuseFailAlloc_960_; 
v_reuseFailAlloc_960_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_960_, 0, v___x_939_);
lean_ctor_set(v_reuseFailAlloc_960_, 1, v___x_934_);
v___x_941_ = v_reuseFailAlloc_960_;
goto v_reusejp_940_;
}
v_reusejp_940_:
{
lean_object* v___x_942_; uint8_t v___x_943_; 
v___x_942_ = lean_unsigned_to_nat(1u);
v___x_943_ = lean_nat_dec_eq(v_den_930_, v___x_942_);
if (v___x_943_ == 0)
{
lean_object* v___x_944_; lean_object* v___x_945_; lean_object* v___x_946_; lean_object* v___x_947_; lean_object* v___x_948_; lean_object* v___x_949_; lean_object* v___x_950_; lean_object* v___x_951_; 
v___x_944_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_945_ = l_Int_repr(v_num_929_);
lean_dec(v_num_929_);
v___x_946_ = lean_string_append(v___x_944_, v___x_945_);
lean_dec_ref(v___x_945_);
v___x_947_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_948_ = lean_string_append(v___x_946_, v___x_947_);
v___x_949_ = l_Nat_reprFast(v_den_930_);
v___x_950_ = lean_string_append(v___x_948_, v___x_949_);
lean_dec_ref(v___x_949_);
v___x_951_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_951_, 0, v___x_950_);
v___y_766_ = v___y_928_;
v___y_767_ = v___x_941_;
v___y_768_ = v___x_951_;
goto v___jp_765_;
}
else
{
lean_object* v___x_952_; lean_object* v___x_953_; uint8_t v___x_954_; 
lean_dec(v_den_930_);
v___x_952_ = lean_unsigned_to_nat(0u);
v___x_953_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_954_ = lean_int_dec_lt(v_num_929_, v___x_953_);
if (v___x_954_ == 0)
{
lean_object* v___x_955_; lean_object* v___x_956_; 
v___x_955_ = l_Int_repr(v_num_929_);
lean_dec(v_num_929_);
v___x_956_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_956_, 0, v___x_955_);
v___y_766_ = v___y_928_;
v___y_767_ = v___x_941_;
v___y_768_ = v___x_956_;
goto v___jp_765_;
}
else
{
lean_object* v___x_957_; lean_object* v___x_958_; lean_object* v___x_959_; 
v___x_957_ = l_Int_repr(v_num_929_);
lean_dec(v_num_929_);
v___x_958_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_958_, 0, v___x_957_);
v___x_959_ = l_Repr_addAppParen(v___x_958_, v___x_952_);
v___y_766_ = v___y_928_;
v___y_767_ = v___x_941_;
v___y_768_ = v___x_959_;
goto v___jp_765_;
}
}
}
}
}
}
}
}
case 3:
{
lean_object* v_lo_968_; lean_object* v_hi_969_; lean_object* v___x_971_; uint8_t v_isShared_972_; uint8_t v_isSharedCheck_1036_; 
v_lo_968_ = lean_ctor_get(v_x_665_, 0);
v_hi_969_ = lean_ctor_get(v_x_665_, 1);
v_isSharedCheck_1036_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1036_ == 0)
{
v___x_971_ = v_x_665_;
v_isShared_972_ = v_isSharedCheck_1036_;
goto v_resetjp_970_;
}
else
{
lean_inc(v_hi_969_);
lean_inc(v_lo_968_);
lean_dec(v_x_665_);
v___x_971_ = lean_box(0);
v_isShared_972_ = v_isSharedCheck_1036_;
goto v_resetjp_970_;
}
v_resetjp_970_:
{
lean_object* v___y_974_; lean_object* v___y_975_; lean_object* v___y_976_; lean_object* v___y_977_; lean_object* v___y_1009_; lean_object* v___x_1032_; uint8_t v___x_1033_; 
v___x_1032_ = lean_unsigned_to_nat(1024u);
v___x_1033_ = lean_nat_dec_le(v___x_1032_, v_prec_666_);
if (v___x_1033_ == 0)
{
lean_object* v___x_1034_; 
v___x_1034_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1009_ = v___x_1034_;
goto v___jp_1008_;
}
else
{
lean_object* v___x_1035_; 
v___x_1035_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1009_ = v___x_1035_;
goto v___jp_1008_;
}
v___jp_973_:
{
lean_object* v_num_978_; lean_object* v_den_979_; lean_object* v___x_981_; uint8_t v_isShared_982_; uint8_t v_isSharedCheck_1007_; 
v_num_978_ = lean_ctor_get(v_hi_969_, 0);
v_den_979_ = lean_ctor_get(v_hi_969_, 1);
v_isSharedCheck_1007_ = !lean_is_exclusive(v_hi_969_);
if (v_isSharedCheck_1007_ == 0)
{
v___x_981_ = v_hi_969_;
v_isShared_982_ = v_isSharedCheck_1007_;
goto v_resetjp_980_;
}
else
{
lean_inc(v_den_979_);
lean_inc(v_num_978_);
lean_dec(v_hi_969_);
v___x_981_ = lean_box(0);
v_isShared_982_ = v_isSharedCheck_1007_;
goto v_resetjp_980_;
}
v_resetjp_980_:
{
lean_object* v___x_984_; 
lean_inc(v___y_976_);
if (v_isShared_982_ == 0)
{
lean_ctor_set_tag(v___x_981_, 5);
lean_ctor_set(v___x_981_, 1, v___y_977_);
lean_ctor_set(v___x_981_, 0, v___y_976_);
v___x_984_ = v___x_981_;
goto v_reusejp_983_;
}
else
{
lean_object* v_reuseFailAlloc_1006_; 
v_reuseFailAlloc_1006_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1006_, 0, v___y_976_);
lean_ctor_set(v_reuseFailAlloc_1006_, 1, v___y_977_);
v___x_984_ = v_reuseFailAlloc_1006_;
goto v_reusejp_983_;
}
v_reusejp_983_:
{
lean_object* v___x_986_; 
lean_inc(v___y_974_);
if (v_isShared_972_ == 0)
{
lean_ctor_set_tag(v___x_971_, 5);
lean_ctor_set(v___x_971_, 1, v___y_974_);
lean_ctor_set(v___x_971_, 0, v___x_984_);
v___x_986_ = v___x_971_;
goto v_reusejp_985_;
}
else
{
lean_object* v_reuseFailAlloc_1005_; 
v_reuseFailAlloc_1005_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1005_, 0, v___x_984_);
lean_ctor_set(v_reuseFailAlloc_1005_, 1, v___y_974_);
v___x_986_ = v_reuseFailAlloc_1005_;
goto v_reusejp_985_;
}
v_reusejp_985_:
{
lean_object* v___x_987_; uint8_t v___x_988_; 
v___x_987_ = lean_unsigned_to_nat(1u);
v___x_988_ = lean_nat_dec_eq(v_den_979_, v___x_987_);
if (v___x_988_ == 0)
{
lean_object* v___x_989_; lean_object* v___x_990_; lean_object* v___x_991_; lean_object* v___x_992_; lean_object* v___x_993_; lean_object* v___x_994_; lean_object* v___x_995_; lean_object* v___x_996_; 
v___x_989_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_990_ = l_Int_repr(v_num_978_);
lean_dec(v_num_978_);
v___x_991_ = lean_string_append(v___x_989_, v___x_990_);
lean_dec_ref(v___x_990_);
v___x_992_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_993_ = lean_string_append(v___x_991_, v___x_992_);
v___x_994_ = l_Nat_reprFast(v_den_979_);
v___x_995_ = lean_string_append(v___x_993_, v___x_994_);
lean_dec_ref(v___x_994_);
v___x_996_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_996_, 0, v___x_995_);
v___y_782_ = v___y_975_;
v___y_783_ = v___x_986_;
v___y_784_ = v___x_996_;
goto v___jp_781_;
}
else
{
lean_object* v___x_997_; lean_object* v___x_998_; uint8_t v___x_999_; 
lean_dec(v_den_979_);
v___x_997_ = lean_unsigned_to_nat(0u);
v___x_998_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_999_ = lean_int_dec_lt(v_num_978_, v___x_998_);
if (v___x_999_ == 0)
{
lean_object* v___x_1000_; lean_object* v___x_1001_; 
v___x_1000_ = l_Int_repr(v_num_978_);
lean_dec(v_num_978_);
v___x_1001_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1001_, 0, v___x_1000_);
v___y_782_ = v___y_975_;
v___y_783_ = v___x_986_;
v___y_784_ = v___x_1001_;
goto v___jp_781_;
}
else
{
lean_object* v___x_1002_; lean_object* v___x_1003_; lean_object* v___x_1004_; 
v___x_1002_ = l_Int_repr(v_num_978_);
lean_dec(v_num_978_);
v___x_1003_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1003_, 0, v___x_1002_);
v___x_1004_ = l_Repr_addAppParen(v___x_1003_, v___x_997_);
v___y_782_ = v___y_975_;
v___y_783_ = v___x_986_;
v___y_784_ = v___x_1004_;
goto v___jp_781_;
}
}
}
}
}
}
v___jp_1008_:
{
lean_object* v_num_1010_; lean_object* v_den_1011_; lean_object* v___x_1012_; lean_object* v___x_1013_; lean_object* v___x_1014_; uint8_t v___x_1015_; 
v_num_1010_ = lean_ctor_get(v_lo_968_, 0);
lean_inc(v_num_1010_);
v_den_1011_ = lean_ctor_get(v_lo_968_, 1);
lean_inc(v_den_1011_);
lean_dec_ref(v_lo_968_);
v___x_1012_ = lean_box(1);
v___x_1013_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__35));
v___x_1014_ = lean_unsigned_to_nat(1u);
v___x_1015_ = lean_nat_dec_eq(v_den_1011_, v___x_1014_);
if (v___x_1015_ == 0)
{
lean_object* v___x_1016_; lean_object* v___x_1017_; lean_object* v___x_1018_; lean_object* v___x_1019_; lean_object* v___x_1020_; lean_object* v___x_1021_; lean_object* v___x_1022_; lean_object* v___x_1023_; 
v___x_1016_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1017_ = l_Int_repr(v_num_1010_);
lean_dec(v_num_1010_);
v___x_1018_ = lean_string_append(v___x_1016_, v___x_1017_);
lean_dec_ref(v___x_1017_);
v___x_1019_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1020_ = lean_string_append(v___x_1018_, v___x_1019_);
v___x_1021_ = l_Nat_reprFast(v_den_1011_);
v___x_1022_ = lean_string_append(v___x_1020_, v___x_1021_);
lean_dec_ref(v___x_1021_);
v___x_1023_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1023_, 0, v___x_1022_);
v___y_974_ = v___x_1012_;
v___y_975_ = v___y_1009_;
v___y_976_ = v___x_1013_;
v___y_977_ = v___x_1023_;
goto v___jp_973_;
}
else
{
lean_object* v___x_1024_; lean_object* v___x_1025_; uint8_t v___x_1026_; 
lean_dec(v_den_1011_);
v___x_1024_ = lean_unsigned_to_nat(0u);
v___x_1025_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1026_ = lean_int_dec_lt(v_num_1010_, v___x_1025_);
if (v___x_1026_ == 0)
{
lean_object* v___x_1027_; lean_object* v___x_1028_; 
v___x_1027_ = l_Int_repr(v_num_1010_);
lean_dec(v_num_1010_);
v___x_1028_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1028_, 0, v___x_1027_);
v___y_974_ = v___x_1012_;
v___y_975_ = v___y_1009_;
v___y_976_ = v___x_1013_;
v___y_977_ = v___x_1028_;
goto v___jp_973_;
}
else
{
lean_object* v___x_1029_; lean_object* v___x_1030_; lean_object* v___x_1031_; 
v___x_1029_ = l_Int_repr(v_num_1010_);
lean_dec(v_num_1010_);
v___x_1030_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1030_, 0, v___x_1029_);
v___x_1031_ = l_Repr_addAppParen(v___x_1030_, v___x_1024_);
v___y_974_ = v___x_1012_;
v___y_975_ = v___y_1009_;
v___y_976_ = v___x_1013_;
v___y_977_ = v___x_1031_;
goto v___jp_973_;
}
}
}
}
}
case 4:
{
lean_object* v_enzymeSites_1037_; lean_object* v___y_1039_; lean_object* v___x_1047_; uint8_t v___x_1048_; 
v_enzymeSites_1037_ = lean_ctor_get(v_x_665_, 0);
lean_inc(v_enzymeSites_1037_);
lean_dec_ref(v_x_665_);
v___x_1047_ = lean_unsigned_to_nat(1024u);
v___x_1048_ = lean_nat_dec_le(v___x_1047_, v_prec_666_);
if (v___x_1048_ == 0)
{
lean_object* v___x_1049_; 
v___x_1049_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1039_ = v___x_1049_;
goto v___jp_1038_;
}
else
{
lean_object* v___x_1050_; 
v___x_1050_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1039_ = v___x_1050_;
goto v___jp_1038_;
}
v___jp_1038_:
{
lean_object* v___x_1040_; lean_object* v___x_1041_; lean_object* v___x_1042_; lean_object* v___x_1043_; uint8_t v___x_1044_; lean_object* v___x_1045_; lean_object* v___x_1046_; 
v___x_1040_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__38));
v___x_1041_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg(v_enzymeSites_1037_);
v___x_1042_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_1042_, 0, v___x_1040_);
lean_ctor_set(v___x_1042_, 1, v___x_1041_);
lean_inc(v___y_1039_);
v___x_1043_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_1043_, 0, v___y_1039_);
lean_ctor_set(v___x_1043_, 1, v___x_1042_);
v___x_1044_ = 0;
v___x_1045_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_1045_, 0, v___x_1043_);
lean_ctor_set_uint8(v___x_1045_, sizeof(void*)*1, v___x_1044_);
v___x_1046_ = l_Repr_addAppParen(v___x_1045_, v_prec_666_);
return v___x_1046_;
}
}
case 5:
{
lean_object* v_readingFrame_1051_; lean_object* v_exonBoundaries_1052_; lean_object* v___x_1054_; uint8_t v_isShared_1055_; uint8_t v_isSharedCheck_1076_; 
v_readingFrame_1051_ = lean_ctor_get(v_x_665_, 0);
v_exonBoundaries_1052_ = lean_ctor_get(v_x_665_, 1);
v_isSharedCheck_1076_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1076_ == 0)
{
v___x_1054_ = v_x_665_;
v_isShared_1055_ = v_isSharedCheck_1076_;
goto v_resetjp_1053_;
}
else
{
lean_inc(v_exonBoundaries_1052_);
lean_inc(v_readingFrame_1051_);
lean_dec(v_x_665_);
v___x_1054_ = lean_box(0);
v_isShared_1055_ = v_isSharedCheck_1076_;
goto v_resetjp_1053_;
}
v_resetjp_1053_:
{
lean_object* v___y_1057_; lean_object* v___x_1072_; uint8_t v___x_1073_; 
v___x_1072_ = lean_unsigned_to_nat(1024u);
v___x_1073_ = lean_nat_dec_le(v___x_1072_, v_prec_666_);
if (v___x_1073_ == 0)
{
lean_object* v___x_1074_; 
v___x_1074_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1057_ = v___x_1074_;
goto v___jp_1056_;
}
else
{
lean_object* v___x_1075_; 
v___x_1075_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1057_ = v___x_1075_;
goto v___jp_1056_;
}
v___jp_1056_:
{
lean_object* v___x_1058_; lean_object* v___x_1059_; lean_object* v___x_1060_; lean_object* v___x_1061_; lean_object* v___x_1063_; 
v___x_1058_ = lean_box(1);
v___x_1059_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__41));
v___x_1060_ = l_Nat_reprFast(v_readingFrame_1051_);
v___x_1061_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1061_, 0, v___x_1060_);
if (v_isShared_1055_ == 0)
{
lean_ctor_set(v___x_1054_, 1, v___x_1061_);
lean_ctor_set(v___x_1054_, 0, v___x_1059_);
v___x_1063_ = v___x_1054_;
goto v_reusejp_1062_;
}
else
{
lean_object* v_reuseFailAlloc_1071_; 
v_reuseFailAlloc_1071_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1071_, 0, v___x_1059_);
lean_ctor_set(v_reuseFailAlloc_1071_, 1, v___x_1061_);
v___x_1063_ = v_reuseFailAlloc_1071_;
goto v_reusejp_1062_;
}
v_reusejp_1062_:
{
lean_object* v___x_1064_; lean_object* v___x_1065_; lean_object* v___x_1066_; lean_object* v___x_1067_; uint8_t v___x_1068_; lean_object* v___x_1069_; lean_object* v___x_1070_; 
v___x_1064_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_1064_, 0, v___x_1063_);
lean_ctor_set(v___x_1064_, 1, v___x_1058_);
v___x_1065_ = lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___redArg(v_exonBoundaries_1052_);
v___x_1066_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_1066_, 0, v___x_1064_);
lean_ctor_set(v___x_1066_, 1, v___x_1065_);
lean_inc(v___y_1057_);
v___x_1067_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_1067_, 0, v___y_1057_);
lean_ctor_set(v___x_1067_, 1, v___x_1066_);
v___x_1068_ = 0;
v___x_1069_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_1069_, 0, v___x_1067_);
lean_ctor_set_uint8(v___x_1069_, sizeof(void*)*1, v___x_1068_);
v___x_1070_ = l_Repr_addAppParen(v___x_1069_, v_prec_666_);
return v___x_1070_;
}
}
}
}
case 6:
{
lean_object* v___x_1077_; uint8_t v___x_1078_; 
v___x_1077_ = lean_unsigned_to_nat(1024u);
v___x_1078_ = lean_nat_dec_le(v___x_1077_, v_prec_666_);
if (v___x_1078_ == 0)
{
lean_object* v___x_1079_; 
v___x_1079_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_791_ = v___x_1079_;
goto v___jp_790_;
}
else
{
lean_object* v___x_1080_; 
v___x_1080_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_791_ = v___x_1080_;
goto v___jp_790_;
}
}
case 7:
{
lean_object* v___x_1081_; uint8_t v___x_1082_; 
v___x_1081_ = lean_unsigned_to_nat(1024u);
v___x_1082_ = lean_nat_dec_le(v___x_1081_, v_prec_666_);
if (v___x_1082_ == 0)
{
lean_object* v___x_1083_; 
v___x_1083_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_798_ = v___x_1083_;
goto v___jp_797_;
}
else
{
lean_object* v___x_1084_; 
v___x_1084_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_798_ = v___x_1084_;
goto v___jp_797_;
}
}
case 8:
{
lean_object* v___x_1085_; uint8_t v___x_1086_; 
v___x_1085_ = lean_unsigned_to_nat(1024u);
v___x_1086_ = lean_nat_dec_le(v___x_1085_, v_prec_666_);
if (v___x_1086_ == 0)
{
lean_object* v___x_1087_; 
v___x_1087_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_805_ = v___x_1087_;
goto v___jp_804_;
}
else
{
lean_object* v___x_1088_; 
v___x_1088_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_805_ = v___x_1088_;
goto v___jp_804_;
}
}
case 9:
{
lean_object* v___x_1089_; uint8_t v___x_1090_; 
v___x_1089_ = lean_unsigned_to_nat(1024u);
v___x_1090_ = lean_nat_dec_le(v___x_1089_, v_prec_666_);
if (v___x_1090_ == 0)
{
lean_object* v___x_1091_; 
v___x_1091_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_812_ = v___x_1091_;
goto v___jp_811_;
}
else
{
lean_object* v___x_1092_; 
v___x_1092_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_812_ = v___x_1092_;
goto v___jp_811_;
}
}
case 10:
{
lean_object* v___x_1093_; uint8_t v___x_1094_; 
v___x_1093_ = lean_unsigned_to_nat(1024u);
v___x_1094_ = lean_nat_dec_le(v___x_1093_, v_prec_666_);
if (v___x_1094_ == 0)
{
lean_object* v___x_1095_; 
v___x_1095_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_819_ = v___x_1095_;
goto v___jp_818_;
}
else
{
lean_object* v___x_1096_; 
v___x_1096_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_819_ = v___x_1096_;
goto v___jp_818_;
}
}
case 11:
{
lean_object* v_organism_1097_; lean_object* v_threshold_1098_; lean_object* v___x_1100_; uint8_t v_isShared_1101_; uint8_t v_isSharedCheck_1142_; 
v_organism_1097_ = lean_ctor_get(v_x_665_, 0);
v_threshold_1098_ = lean_ctor_get(v_x_665_, 1);
v_isSharedCheck_1142_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1142_ == 0)
{
v___x_1100_ = v_x_665_;
v_isShared_1101_ = v_isSharedCheck_1142_;
goto v_resetjp_1099_;
}
else
{
lean_inc(v_threshold_1098_);
lean_inc(v_organism_1097_);
lean_dec(v_x_665_);
v___x_1100_ = lean_box(0);
v_isShared_1101_ = v_isSharedCheck_1142_;
goto v_resetjp_1099_;
}
v_resetjp_1099_:
{
lean_object* v___y_1103_; lean_object* v___x_1138_; uint8_t v___x_1139_; 
v___x_1138_ = lean_unsigned_to_nat(1024u);
v___x_1139_ = lean_nat_dec_le(v___x_1138_, v_prec_666_);
if (v___x_1139_ == 0)
{
lean_object* v___x_1140_; 
v___x_1140_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1103_ = v___x_1140_;
goto v___jp_1102_;
}
else
{
lean_object* v___x_1141_; 
v___x_1141_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1103_ = v___x_1141_;
goto v___jp_1102_;
}
v___jp_1102_:
{
lean_object* v_num_1104_; lean_object* v_den_1105_; lean_object* v___x_1107_; uint8_t v_isShared_1108_; uint8_t v_isSharedCheck_1137_; 
v_num_1104_ = lean_ctor_get(v_threshold_1098_, 0);
v_den_1105_ = lean_ctor_get(v_threshold_1098_, 1);
v_isSharedCheck_1137_ = !lean_is_exclusive(v_threshold_1098_);
if (v_isSharedCheck_1137_ == 0)
{
v___x_1107_ = v_threshold_1098_;
v_isShared_1108_ = v_isSharedCheck_1137_;
goto v_resetjp_1106_;
}
else
{
lean_inc(v_den_1105_);
lean_inc(v_num_1104_);
lean_dec(v_threshold_1098_);
v___x_1107_ = lean_box(0);
v_isShared_1108_ = v_isSharedCheck_1137_;
goto v_resetjp_1106_;
}
v_resetjp_1106_:
{
lean_object* v___x_1109_; lean_object* v___x_1110_; lean_object* v___x_1111_; lean_object* v___x_1112_; lean_object* v___x_1114_; 
v___x_1109_ = lean_box(1);
v___x_1110_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__44));
v___x_1111_ = l_String_quote(v_organism_1097_);
v___x_1112_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1112_, 0, v___x_1111_);
if (v_isShared_1108_ == 0)
{
lean_ctor_set_tag(v___x_1107_, 5);
lean_ctor_set(v___x_1107_, 1, v___x_1112_);
lean_ctor_set(v___x_1107_, 0, v___x_1110_);
v___x_1114_ = v___x_1107_;
goto v_reusejp_1113_;
}
else
{
lean_object* v_reuseFailAlloc_1136_; 
v_reuseFailAlloc_1136_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1136_, 0, v___x_1110_);
lean_ctor_set(v_reuseFailAlloc_1136_, 1, v___x_1112_);
v___x_1114_ = v_reuseFailAlloc_1136_;
goto v_reusejp_1113_;
}
v_reusejp_1113_:
{
lean_object* v___x_1116_; 
if (v_isShared_1101_ == 0)
{
lean_ctor_set_tag(v___x_1100_, 5);
lean_ctor_set(v___x_1100_, 1, v___x_1109_);
lean_ctor_set(v___x_1100_, 0, v___x_1114_);
v___x_1116_ = v___x_1100_;
goto v_reusejp_1115_;
}
else
{
lean_object* v_reuseFailAlloc_1135_; 
v_reuseFailAlloc_1135_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1135_, 0, v___x_1114_);
lean_ctor_set(v_reuseFailAlloc_1135_, 1, v___x_1109_);
v___x_1116_ = v_reuseFailAlloc_1135_;
goto v_reusejp_1115_;
}
v_reusejp_1115_:
{
lean_object* v___x_1117_; uint8_t v___x_1118_; 
v___x_1117_ = lean_unsigned_to_nat(1u);
v___x_1118_ = lean_nat_dec_eq(v_den_1105_, v___x_1117_);
if (v___x_1118_ == 0)
{
lean_object* v___x_1119_; lean_object* v___x_1120_; lean_object* v___x_1121_; lean_object* v___x_1122_; lean_object* v___x_1123_; lean_object* v___x_1124_; lean_object* v___x_1125_; lean_object* v___x_1126_; 
v___x_1119_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1120_ = l_Int_repr(v_num_1104_);
lean_dec(v_num_1104_);
v___x_1121_ = lean_string_append(v___x_1119_, v___x_1120_);
lean_dec_ref(v___x_1120_);
v___x_1122_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1123_ = lean_string_append(v___x_1121_, v___x_1122_);
v___x_1124_ = l_Nat_reprFast(v_den_1105_);
v___x_1125_ = lean_string_append(v___x_1123_, v___x_1124_);
lean_dec_ref(v___x_1124_);
v___x_1126_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1126_, 0, v___x_1125_);
v___y_826_ = v___x_1116_;
v___y_827_ = v___y_1103_;
v___y_828_ = v___x_1126_;
goto v___jp_825_;
}
else
{
lean_object* v___x_1127_; lean_object* v___x_1128_; uint8_t v___x_1129_; 
lean_dec(v_den_1105_);
v___x_1127_ = lean_unsigned_to_nat(0u);
v___x_1128_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1129_ = lean_int_dec_lt(v_num_1104_, v___x_1128_);
if (v___x_1129_ == 0)
{
lean_object* v___x_1130_; lean_object* v___x_1131_; 
v___x_1130_ = l_Int_repr(v_num_1104_);
lean_dec(v_num_1104_);
v___x_1131_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1131_, 0, v___x_1130_);
v___y_826_ = v___x_1116_;
v___y_827_ = v___y_1103_;
v___y_828_ = v___x_1131_;
goto v___jp_825_;
}
else
{
lean_object* v___x_1132_; lean_object* v___x_1133_; lean_object* v___x_1134_; 
v___x_1132_ = l_Int_repr(v_num_1104_);
lean_dec(v_num_1104_);
v___x_1133_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1133_, 0, v___x_1132_);
v___x_1134_ = l_Repr_addAppParen(v___x_1133_, v___x_1127_);
v___y_826_ = v___x_1116_;
v___y_827_ = v___y_1103_;
v___y_828_ = v___x_1134_;
goto v___jp_825_;
}
}
}
}
}
}
}
}
case 12:
{
lean_object* v_organism_1143_; lean_object* v_threshold_1144_; lean_object* v___x_1146_; uint8_t v_isShared_1147_; uint8_t v_isSharedCheck_1188_; 
v_organism_1143_ = lean_ctor_get(v_x_665_, 0);
v_threshold_1144_ = lean_ctor_get(v_x_665_, 1);
v_isSharedCheck_1188_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1188_ == 0)
{
v___x_1146_ = v_x_665_;
v_isShared_1147_ = v_isSharedCheck_1188_;
goto v_resetjp_1145_;
}
else
{
lean_inc(v_threshold_1144_);
lean_inc(v_organism_1143_);
lean_dec(v_x_665_);
v___x_1146_ = lean_box(0);
v_isShared_1147_ = v_isSharedCheck_1188_;
goto v_resetjp_1145_;
}
v_resetjp_1145_:
{
lean_object* v___y_1149_; lean_object* v___x_1184_; uint8_t v___x_1185_; 
v___x_1184_ = lean_unsigned_to_nat(1024u);
v___x_1185_ = lean_nat_dec_le(v___x_1184_, v_prec_666_);
if (v___x_1185_ == 0)
{
lean_object* v___x_1186_; 
v___x_1186_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1149_ = v___x_1186_;
goto v___jp_1148_;
}
else
{
lean_object* v___x_1187_; 
v___x_1187_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1149_ = v___x_1187_;
goto v___jp_1148_;
}
v___jp_1148_:
{
lean_object* v_num_1150_; lean_object* v_den_1151_; lean_object* v___x_1153_; uint8_t v_isShared_1154_; uint8_t v_isSharedCheck_1183_; 
v_num_1150_ = lean_ctor_get(v_threshold_1144_, 0);
v_den_1151_ = lean_ctor_get(v_threshold_1144_, 1);
v_isSharedCheck_1183_ = !lean_is_exclusive(v_threshold_1144_);
if (v_isSharedCheck_1183_ == 0)
{
v___x_1153_ = v_threshold_1144_;
v_isShared_1154_ = v_isSharedCheck_1183_;
goto v_resetjp_1152_;
}
else
{
lean_inc(v_den_1151_);
lean_inc(v_num_1150_);
lean_dec(v_threshold_1144_);
v___x_1153_ = lean_box(0);
v_isShared_1154_ = v_isSharedCheck_1183_;
goto v_resetjp_1152_;
}
v_resetjp_1152_:
{
lean_object* v___x_1155_; lean_object* v___x_1156_; lean_object* v___x_1157_; lean_object* v___x_1158_; lean_object* v___x_1160_; 
v___x_1155_ = lean_box(1);
v___x_1156_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__47));
v___x_1157_ = l_String_quote(v_organism_1143_);
v___x_1158_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1158_, 0, v___x_1157_);
if (v_isShared_1154_ == 0)
{
lean_ctor_set_tag(v___x_1153_, 5);
lean_ctor_set(v___x_1153_, 1, v___x_1158_);
lean_ctor_set(v___x_1153_, 0, v___x_1156_);
v___x_1160_ = v___x_1153_;
goto v_reusejp_1159_;
}
else
{
lean_object* v_reuseFailAlloc_1182_; 
v_reuseFailAlloc_1182_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1182_, 0, v___x_1156_);
lean_ctor_set(v_reuseFailAlloc_1182_, 1, v___x_1158_);
v___x_1160_ = v_reuseFailAlloc_1182_;
goto v_reusejp_1159_;
}
v_reusejp_1159_:
{
lean_object* v___x_1162_; 
if (v_isShared_1147_ == 0)
{
lean_ctor_set_tag(v___x_1146_, 5);
lean_ctor_set(v___x_1146_, 1, v___x_1155_);
lean_ctor_set(v___x_1146_, 0, v___x_1160_);
v___x_1162_ = v___x_1146_;
goto v_reusejp_1161_;
}
else
{
lean_object* v_reuseFailAlloc_1181_; 
v_reuseFailAlloc_1181_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1181_, 0, v___x_1160_);
lean_ctor_set(v_reuseFailAlloc_1181_, 1, v___x_1155_);
v___x_1162_ = v_reuseFailAlloc_1181_;
goto v_reusejp_1161_;
}
v_reusejp_1161_:
{
lean_object* v___x_1163_; uint8_t v___x_1164_; 
v___x_1163_ = lean_unsigned_to_nat(1u);
v___x_1164_ = lean_nat_dec_eq(v_den_1151_, v___x_1163_);
if (v___x_1164_ == 0)
{
lean_object* v___x_1165_; lean_object* v___x_1166_; lean_object* v___x_1167_; lean_object* v___x_1168_; lean_object* v___x_1169_; lean_object* v___x_1170_; lean_object* v___x_1171_; lean_object* v___x_1172_; 
v___x_1165_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1166_ = l_Int_repr(v_num_1150_);
lean_dec(v_num_1150_);
v___x_1167_ = lean_string_append(v___x_1165_, v___x_1166_);
lean_dec_ref(v___x_1166_);
v___x_1168_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1169_ = lean_string_append(v___x_1167_, v___x_1168_);
v___x_1170_ = l_Nat_reprFast(v_den_1151_);
v___x_1171_ = lean_string_append(v___x_1169_, v___x_1170_);
lean_dec_ref(v___x_1170_);
v___x_1172_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1172_, 0, v___x_1171_);
v___y_757_ = v___x_1162_;
v___y_758_ = v___y_1149_;
v___y_759_ = v___x_1172_;
goto v___jp_756_;
}
else
{
lean_object* v___x_1173_; lean_object* v___x_1174_; uint8_t v___x_1175_; 
lean_dec(v_den_1151_);
v___x_1173_ = lean_unsigned_to_nat(0u);
v___x_1174_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1175_ = lean_int_dec_lt(v_num_1150_, v___x_1174_);
if (v___x_1175_ == 0)
{
lean_object* v___x_1176_; lean_object* v___x_1177_; 
v___x_1176_ = l_Int_repr(v_num_1150_);
lean_dec(v_num_1150_);
v___x_1177_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1177_, 0, v___x_1176_);
v___y_757_ = v___x_1162_;
v___y_758_ = v___y_1149_;
v___y_759_ = v___x_1177_;
goto v___jp_756_;
}
else
{
lean_object* v___x_1178_; lean_object* v___x_1179_; lean_object* v___x_1180_; 
v___x_1178_ = l_Int_repr(v_num_1150_);
lean_dec(v_num_1150_);
v___x_1179_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1179_, 0, v___x_1178_);
v___x_1180_ = l_Repr_addAppParen(v___x_1179_, v___x_1173_);
v___y_757_ = v___x_1162_;
v___y_758_ = v___y_1149_;
v___y_759_ = v___x_1180_;
goto v___jp_756_;
}
}
}
}
}
}
}
}
case 13:
{
lean_object* v_minScore_1189_; lean_object* v___x_1191_; uint8_t v_isShared_1192_; uint8_t v_isSharedCheck_1212_; 
v_minScore_1189_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1212_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1212_ == 0)
{
v___x_1191_ = v_x_665_;
v_isShared_1192_ = v_isSharedCheck_1212_;
goto v_resetjp_1190_;
}
else
{
lean_inc(v_minScore_1189_);
lean_dec(v_x_665_);
v___x_1191_ = lean_box(0);
v_isShared_1192_ = v_isSharedCheck_1212_;
goto v_resetjp_1190_;
}
v_resetjp_1190_:
{
lean_object* v___y_1194_; lean_object* v___x_1208_; uint8_t v___x_1209_; 
v___x_1208_ = lean_unsigned_to_nat(1024u);
v___x_1209_ = lean_nat_dec_le(v___x_1208_, v_prec_666_);
if (v___x_1209_ == 0)
{
lean_object* v___x_1210_; 
v___x_1210_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1194_ = v___x_1210_;
goto v___jp_1193_;
}
else
{
lean_object* v___x_1211_; 
v___x_1211_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1194_ = v___x_1211_;
goto v___jp_1193_;
}
v___jp_1193_:
{
lean_object* v___x_1195_; lean_object* v___x_1196_; uint8_t v___x_1197_; 
v___x_1195_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__50));
v___x_1196_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1197_ = lean_int_dec_lt(v_minScore_1189_, v___x_1196_);
if (v___x_1197_ == 0)
{
lean_object* v___x_1198_; lean_object* v___x_1200_; 
v___x_1198_ = l_Int_repr(v_minScore_1189_);
lean_dec(v_minScore_1189_);
if (v_isShared_1192_ == 0)
{
lean_ctor_set_tag(v___x_1191_, 3);
lean_ctor_set(v___x_1191_, 0, v___x_1198_);
v___x_1200_ = v___x_1191_;
goto v_reusejp_1199_;
}
else
{
lean_object* v_reuseFailAlloc_1201_; 
v_reuseFailAlloc_1201_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1201_, 0, v___x_1198_);
v___x_1200_ = v_reuseFailAlloc_1201_;
goto v_reusejp_1199_;
}
v_reusejp_1199_:
{
v___y_835_ = v___y_1194_;
v___y_836_ = v___x_1195_;
v___y_837_ = v___x_1200_;
goto v___jp_834_;
}
}
else
{
lean_object* v___x_1202_; lean_object* v___x_1203_; lean_object* v___x_1205_; 
v___x_1202_ = lean_unsigned_to_nat(1024u);
v___x_1203_ = l_Int_repr(v_minScore_1189_);
lean_dec(v_minScore_1189_);
if (v_isShared_1192_ == 0)
{
lean_ctor_set_tag(v___x_1191_, 3);
lean_ctor_set(v___x_1191_, 0, v___x_1203_);
v___x_1205_ = v___x_1191_;
goto v_reusejp_1204_;
}
else
{
lean_object* v_reuseFailAlloc_1207_; 
v_reuseFailAlloc_1207_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1207_, 0, v___x_1203_);
v___x_1205_ = v_reuseFailAlloc_1207_;
goto v_reusejp_1204_;
}
v_reusejp_1204_:
{
lean_object* v___x_1206_; 
v___x_1206_ = l_Repr_addAppParen(v___x_1205_, v___x_1202_);
v___y_835_ = v___y_1194_;
v___y_836_ = v___x_1195_;
v___y_837_ = v___x_1206_;
goto v___jp_834_;
}
}
}
}
}
case 14:
{
uint8_t v_isCytosolic_1213_; lean_object* v_threshold_1214_; lean_object* v___y_1216_; lean_object* v___x_1248_; uint8_t v___x_1249_; 
v_isCytosolic_1213_ = lean_ctor_get_uint8(v_x_665_, sizeof(void*)*1);
v_threshold_1214_ = lean_ctor_get(v_x_665_, 0);
lean_inc_ref(v_threshold_1214_);
lean_dec_ref(v_x_665_);
v___x_1248_ = lean_unsigned_to_nat(1024u);
v___x_1249_ = lean_nat_dec_le(v___x_1248_, v_prec_666_);
if (v___x_1249_ == 0)
{
lean_object* v___x_1250_; 
v___x_1250_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1216_ = v___x_1250_;
goto v___jp_1215_;
}
else
{
lean_object* v___x_1251_; 
v___x_1251_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1216_ = v___x_1251_;
goto v___jp_1215_;
}
v___jp_1215_:
{
lean_object* v_num_1217_; lean_object* v_den_1218_; lean_object* v___x_1220_; uint8_t v_isShared_1221_; uint8_t v_isSharedCheck_1247_; 
v_num_1217_ = lean_ctor_get(v_threshold_1214_, 0);
v_den_1218_ = lean_ctor_get(v_threshold_1214_, 1);
v_isSharedCheck_1247_ = !lean_is_exclusive(v_threshold_1214_);
if (v_isSharedCheck_1247_ == 0)
{
v___x_1220_ = v_threshold_1214_;
v_isShared_1221_ = v_isSharedCheck_1247_;
goto v_resetjp_1219_;
}
else
{
lean_inc(v_den_1218_);
lean_inc(v_num_1217_);
lean_dec(v_threshold_1214_);
v___x_1220_ = lean_box(0);
v_isShared_1221_ = v_isSharedCheck_1247_;
goto v_resetjp_1219_;
}
v_resetjp_1219_:
{
lean_object* v___x_1222_; lean_object* v___x_1223_; lean_object* v___x_1224_; lean_object* v___x_1226_; 
v___x_1222_ = lean_box(1);
v___x_1223_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__53));
v___x_1224_ = l_Bool_repr___redArg(v_isCytosolic_1213_);
if (v_isShared_1221_ == 0)
{
lean_ctor_set_tag(v___x_1220_, 5);
lean_ctor_set(v___x_1220_, 1, v___x_1224_);
lean_ctor_set(v___x_1220_, 0, v___x_1223_);
v___x_1226_ = v___x_1220_;
goto v_reusejp_1225_;
}
else
{
lean_object* v_reuseFailAlloc_1246_; 
v_reuseFailAlloc_1246_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1246_, 0, v___x_1223_);
lean_ctor_set(v_reuseFailAlloc_1246_, 1, v___x_1224_);
v___x_1226_ = v_reuseFailAlloc_1246_;
goto v_reusejp_1225_;
}
v_reusejp_1225_:
{
lean_object* v___x_1227_; lean_object* v___x_1228_; uint8_t v___x_1229_; 
v___x_1227_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_1227_, 0, v___x_1226_);
lean_ctor_set(v___x_1227_, 1, v___x_1222_);
v___x_1228_ = lean_unsigned_to_nat(1u);
v___x_1229_ = lean_nat_dec_eq(v_den_1218_, v___x_1228_);
if (v___x_1229_ == 0)
{
lean_object* v___x_1230_; lean_object* v___x_1231_; lean_object* v___x_1232_; lean_object* v___x_1233_; lean_object* v___x_1234_; lean_object* v___x_1235_; lean_object* v___x_1236_; lean_object* v___x_1237_; 
v___x_1230_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1231_ = l_Int_repr(v_num_1217_);
lean_dec(v_num_1217_);
v___x_1232_ = lean_string_append(v___x_1230_, v___x_1231_);
lean_dec_ref(v___x_1231_);
v___x_1233_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1234_ = lean_string_append(v___x_1232_, v___x_1233_);
v___x_1235_ = l_Nat_reprFast(v_den_1218_);
v___x_1236_ = lean_string_append(v___x_1234_, v___x_1235_);
lean_dec_ref(v___x_1235_);
v___x_1237_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1237_, 0, v___x_1236_);
v___y_748_ = v___x_1227_;
v___y_749_ = v___y_1216_;
v___y_750_ = v___x_1237_;
goto v___jp_747_;
}
else
{
lean_object* v___x_1238_; lean_object* v___x_1239_; uint8_t v___x_1240_; 
lean_dec(v_den_1218_);
v___x_1238_ = lean_unsigned_to_nat(0u);
v___x_1239_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1240_ = lean_int_dec_lt(v_num_1217_, v___x_1239_);
if (v___x_1240_ == 0)
{
lean_object* v___x_1241_; lean_object* v___x_1242_; 
v___x_1241_ = l_Int_repr(v_num_1217_);
lean_dec(v_num_1217_);
v___x_1242_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1242_, 0, v___x_1241_);
v___y_748_ = v___x_1227_;
v___y_749_ = v___y_1216_;
v___y_750_ = v___x_1242_;
goto v___jp_747_;
}
else
{
lean_object* v___x_1243_; lean_object* v___x_1244_; lean_object* v___x_1245_; 
v___x_1243_ = l_Int_repr(v_num_1217_);
lean_dec(v_num_1217_);
v___x_1244_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1244_, 0, v___x_1243_);
v___x_1245_ = l_Repr_addAppParen(v___x_1244_, v___x_1238_);
v___y_748_ = v___x_1227_;
v___y_749_ = v___y_1216_;
v___y_750_ = v___x_1245_;
goto v___jp_747_;
}
}
}
}
}
}
case 15:
{
lean_object* v_dgThreshold_1252_; lean_object* v___x_1254_; uint8_t v_isShared_1255_; uint8_t v_isSharedCheck_1289_; 
v_dgThreshold_1252_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1289_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1289_ == 0)
{
v___x_1254_ = v_x_665_;
v_isShared_1255_ = v_isSharedCheck_1289_;
goto v_resetjp_1253_;
}
else
{
lean_inc(v_dgThreshold_1252_);
lean_dec(v_x_665_);
v___x_1254_ = lean_box(0);
v_isShared_1255_ = v_isSharedCheck_1289_;
goto v_resetjp_1253_;
}
v_resetjp_1253_:
{
lean_object* v___y_1257_; lean_object* v___x_1285_; uint8_t v___x_1286_; 
v___x_1285_ = lean_unsigned_to_nat(1024u);
v___x_1286_ = lean_nat_dec_le(v___x_1285_, v_prec_666_);
if (v___x_1286_ == 0)
{
lean_object* v___x_1287_; 
v___x_1287_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1257_ = v___x_1287_;
goto v___jp_1256_;
}
else
{
lean_object* v___x_1288_; 
v___x_1288_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1257_ = v___x_1288_;
goto v___jp_1256_;
}
v___jp_1256_:
{
lean_object* v_num_1258_; lean_object* v_den_1259_; lean_object* v___x_1260_; lean_object* v___x_1261_; uint8_t v___x_1262_; 
v_num_1258_ = lean_ctor_get(v_dgThreshold_1252_, 0);
lean_inc(v_num_1258_);
v_den_1259_ = lean_ctor_get(v_dgThreshold_1252_, 1);
lean_inc(v_den_1259_);
lean_dec_ref(v_dgThreshold_1252_);
v___x_1260_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__56));
v___x_1261_ = lean_unsigned_to_nat(1u);
v___x_1262_ = lean_nat_dec_eq(v_den_1259_, v___x_1261_);
if (v___x_1262_ == 0)
{
lean_object* v___x_1263_; lean_object* v___x_1264_; lean_object* v___x_1265_; lean_object* v___x_1266_; lean_object* v___x_1267_; lean_object* v___x_1268_; lean_object* v___x_1269_; lean_object* v___x_1271_; 
v___x_1263_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1264_ = l_Int_repr(v_num_1258_);
lean_dec(v_num_1258_);
v___x_1265_ = lean_string_append(v___x_1263_, v___x_1264_);
lean_dec_ref(v___x_1264_);
v___x_1266_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1267_ = lean_string_append(v___x_1265_, v___x_1266_);
v___x_1268_ = l_Nat_reprFast(v_den_1259_);
v___x_1269_ = lean_string_append(v___x_1267_, v___x_1268_);
lean_dec_ref(v___x_1268_);
if (v_isShared_1255_ == 0)
{
lean_ctor_set_tag(v___x_1254_, 3);
lean_ctor_set(v___x_1254_, 0, v___x_1269_);
v___x_1271_ = v___x_1254_;
goto v_reusejp_1270_;
}
else
{
lean_object* v_reuseFailAlloc_1272_; 
v_reuseFailAlloc_1272_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1272_, 0, v___x_1269_);
v___x_1271_ = v_reuseFailAlloc_1272_;
goto v_reusejp_1270_;
}
v_reusejp_1270_:
{
v___y_844_ = v___x_1260_;
v___y_845_ = v___y_1257_;
v___y_846_ = v___x_1271_;
goto v___jp_843_;
}
}
else
{
lean_object* v___x_1273_; lean_object* v___x_1274_; uint8_t v___x_1275_; 
lean_dec(v_den_1259_);
v___x_1273_ = lean_unsigned_to_nat(0u);
v___x_1274_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1275_ = lean_int_dec_lt(v_num_1258_, v___x_1274_);
if (v___x_1275_ == 0)
{
lean_object* v___x_1276_; lean_object* v___x_1278_; 
v___x_1276_ = l_Int_repr(v_num_1258_);
lean_dec(v_num_1258_);
if (v_isShared_1255_ == 0)
{
lean_ctor_set_tag(v___x_1254_, 3);
lean_ctor_set(v___x_1254_, 0, v___x_1276_);
v___x_1278_ = v___x_1254_;
goto v_reusejp_1277_;
}
else
{
lean_object* v_reuseFailAlloc_1279_; 
v_reuseFailAlloc_1279_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1279_, 0, v___x_1276_);
v___x_1278_ = v_reuseFailAlloc_1279_;
goto v_reusejp_1277_;
}
v_reusejp_1277_:
{
v___y_844_ = v___x_1260_;
v___y_845_ = v___y_1257_;
v___y_846_ = v___x_1278_;
goto v___jp_843_;
}
}
else
{
lean_object* v___x_1280_; lean_object* v___x_1282_; 
v___x_1280_ = l_Int_repr(v_num_1258_);
lean_dec(v_num_1258_);
if (v_isShared_1255_ == 0)
{
lean_ctor_set_tag(v___x_1254_, 3);
lean_ctor_set(v___x_1254_, 0, v___x_1280_);
v___x_1282_ = v___x_1254_;
goto v_reusejp_1281_;
}
else
{
lean_object* v_reuseFailAlloc_1284_; 
v_reuseFailAlloc_1284_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1284_, 0, v___x_1280_);
v___x_1282_ = v_reuseFailAlloc_1284_;
goto v_reusejp_1281_;
}
v_reusejp_1281_:
{
lean_object* v___x_1283_; 
v___x_1283_ = l_Repr_addAppParen(v___x_1282_, v___x_1273_);
v___y_844_ = v___x_1260_;
v___y_845_ = v___y_1257_;
v___y_846_ = v___x_1283_;
goto v___jp_843_;
}
}
}
}
}
}
case 16:
{
lean_object* v_organism_1290_; lean_object* v___x_1292_; uint8_t v_isShared_1293_; uint8_t v_isSharedCheck_1310_; 
v_organism_1290_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1310_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1310_ == 0)
{
v___x_1292_ = v_x_665_;
v_isShared_1293_ = v_isSharedCheck_1310_;
goto v_resetjp_1291_;
}
else
{
lean_inc(v_organism_1290_);
lean_dec(v_x_665_);
v___x_1292_ = lean_box(0);
v_isShared_1293_ = v_isSharedCheck_1310_;
goto v_resetjp_1291_;
}
v_resetjp_1291_:
{
lean_object* v___y_1295_; lean_object* v___x_1306_; uint8_t v___x_1307_; 
v___x_1306_ = lean_unsigned_to_nat(1024u);
v___x_1307_ = lean_nat_dec_le(v___x_1306_, v_prec_666_);
if (v___x_1307_ == 0)
{
lean_object* v___x_1308_; 
v___x_1308_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1295_ = v___x_1308_;
goto v___jp_1294_;
}
else
{
lean_object* v___x_1309_; 
v___x_1309_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1295_ = v___x_1309_;
goto v___jp_1294_;
}
v___jp_1294_:
{
lean_object* v___x_1296_; lean_object* v___x_1297_; lean_object* v___x_1299_; 
v___x_1296_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__59));
v___x_1297_ = l_String_quote(v_organism_1290_);
if (v_isShared_1293_ == 0)
{
lean_ctor_set_tag(v___x_1292_, 3);
lean_ctor_set(v___x_1292_, 0, v___x_1297_);
v___x_1299_ = v___x_1292_;
goto v_reusejp_1298_;
}
else
{
lean_object* v_reuseFailAlloc_1305_; 
v_reuseFailAlloc_1305_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1305_, 0, v___x_1297_);
v___x_1299_ = v_reuseFailAlloc_1305_;
goto v_reusejp_1298_;
}
v_reusejp_1298_:
{
lean_object* v___x_1300_; lean_object* v___x_1301_; uint8_t v___x_1302_; lean_object* v___x_1303_; lean_object* v___x_1304_; 
v___x_1300_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_1300_, 0, v___x_1296_);
lean_ctor_set(v___x_1300_, 1, v___x_1299_);
lean_inc(v___y_1295_);
v___x_1301_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_1301_, 0, v___y_1295_);
lean_ctor_set(v___x_1301_, 1, v___x_1300_);
v___x_1302_ = 0;
v___x_1303_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_1303_, 0, v___x_1301_);
lean_ctor_set_uint8(v___x_1303_, sizeof(void*)*1, v___x_1302_);
v___x_1304_ = l_Repr_addAppParen(v___x_1303_, v_prec_666_);
return v___x_1304_;
}
}
}
}
case 17:
{
lean_object* v_threshold_1311_; lean_object* v___x_1313_; uint8_t v_isShared_1314_; uint8_t v_isSharedCheck_1348_; 
v_threshold_1311_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1348_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1348_ == 0)
{
v___x_1313_ = v_x_665_;
v_isShared_1314_ = v_isSharedCheck_1348_;
goto v_resetjp_1312_;
}
else
{
lean_inc(v_threshold_1311_);
lean_dec(v_x_665_);
v___x_1313_ = lean_box(0);
v_isShared_1314_ = v_isSharedCheck_1348_;
goto v_resetjp_1312_;
}
v_resetjp_1312_:
{
lean_object* v___y_1316_; lean_object* v___x_1344_; uint8_t v___x_1345_; 
v___x_1344_ = lean_unsigned_to_nat(1024u);
v___x_1345_ = lean_nat_dec_le(v___x_1344_, v_prec_666_);
if (v___x_1345_ == 0)
{
lean_object* v___x_1346_; 
v___x_1346_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1316_ = v___x_1346_;
goto v___jp_1315_;
}
else
{
lean_object* v___x_1347_; 
v___x_1347_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1316_ = v___x_1347_;
goto v___jp_1315_;
}
v___jp_1315_:
{
lean_object* v_num_1317_; lean_object* v_den_1318_; lean_object* v___x_1319_; lean_object* v___x_1320_; uint8_t v___x_1321_; 
v_num_1317_ = lean_ctor_get(v_threshold_1311_, 0);
lean_inc(v_num_1317_);
v_den_1318_ = lean_ctor_get(v_threshold_1311_, 1);
lean_inc(v_den_1318_);
lean_dec_ref(v_threshold_1311_);
v___x_1319_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__62));
v___x_1320_ = lean_unsigned_to_nat(1u);
v___x_1321_ = lean_nat_dec_eq(v_den_1318_, v___x_1320_);
if (v___x_1321_ == 0)
{
lean_object* v___x_1322_; lean_object* v___x_1323_; lean_object* v___x_1324_; lean_object* v___x_1325_; lean_object* v___x_1326_; lean_object* v___x_1327_; lean_object* v___x_1328_; lean_object* v___x_1330_; 
v___x_1322_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1323_ = l_Int_repr(v_num_1317_);
lean_dec(v_num_1317_);
v___x_1324_ = lean_string_append(v___x_1322_, v___x_1323_);
lean_dec_ref(v___x_1323_);
v___x_1325_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1326_ = lean_string_append(v___x_1324_, v___x_1325_);
v___x_1327_ = l_Nat_reprFast(v_den_1318_);
v___x_1328_ = lean_string_append(v___x_1326_, v___x_1327_);
lean_dec_ref(v___x_1327_);
if (v_isShared_1314_ == 0)
{
lean_ctor_set_tag(v___x_1313_, 3);
lean_ctor_set(v___x_1313_, 0, v___x_1328_);
v___x_1330_ = v___x_1313_;
goto v_reusejp_1329_;
}
else
{
lean_object* v_reuseFailAlloc_1331_; 
v_reuseFailAlloc_1331_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1331_, 0, v___x_1328_);
v___x_1330_ = v_reuseFailAlloc_1331_;
goto v_reusejp_1329_;
}
v_reusejp_1329_:
{
v___y_853_ = v___y_1316_;
v___y_854_ = v___x_1319_;
v___y_855_ = v___x_1330_;
goto v___jp_852_;
}
}
else
{
lean_object* v___x_1332_; lean_object* v___x_1333_; uint8_t v___x_1334_; 
lean_dec(v_den_1318_);
v___x_1332_ = lean_unsigned_to_nat(0u);
v___x_1333_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1334_ = lean_int_dec_lt(v_num_1317_, v___x_1333_);
if (v___x_1334_ == 0)
{
lean_object* v___x_1335_; lean_object* v___x_1337_; 
v___x_1335_ = l_Int_repr(v_num_1317_);
lean_dec(v_num_1317_);
if (v_isShared_1314_ == 0)
{
lean_ctor_set_tag(v___x_1313_, 3);
lean_ctor_set(v___x_1313_, 0, v___x_1335_);
v___x_1337_ = v___x_1313_;
goto v_reusejp_1336_;
}
else
{
lean_object* v_reuseFailAlloc_1338_; 
v_reuseFailAlloc_1338_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1338_, 0, v___x_1335_);
v___x_1337_ = v_reuseFailAlloc_1338_;
goto v_reusejp_1336_;
}
v_reusejp_1336_:
{
v___y_853_ = v___y_1316_;
v___y_854_ = v___x_1319_;
v___y_855_ = v___x_1337_;
goto v___jp_852_;
}
}
else
{
lean_object* v___x_1339_; lean_object* v___x_1341_; 
v___x_1339_ = l_Int_repr(v_num_1317_);
lean_dec(v_num_1317_);
if (v_isShared_1314_ == 0)
{
lean_ctor_set_tag(v___x_1313_, 3);
lean_ctor_set(v___x_1313_, 0, v___x_1339_);
v___x_1341_ = v___x_1313_;
goto v_reusejp_1340_;
}
else
{
lean_object* v_reuseFailAlloc_1343_; 
v_reuseFailAlloc_1343_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1343_, 0, v___x_1339_);
v___x_1341_ = v_reuseFailAlloc_1343_;
goto v_reusejp_1340_;
}
v_reusejp_1340_:
{
lean_object* v___x_1342_; 
v___x_1342_ = l_Repr_addAppParen(v___x_1341_, v___x_1332_);
v___y_853_ = v___y_1316_;
v___y_854_ = v___x_1319_;
v___y_855_ = v___x_1342_;
goto v___jp_852_;
}
}
}
}
}
}
case 18:
{
lean_object* v___x_1349_; uint8_t v___x_1350_; 
v___x_1349_ = lean_unsigned_to_nat(1024u);
v___x_1350_ = lean_nat_dec_le(v___x_1349_, v_prec_666_);
if (v___x_1350_ == 0)
{
lean_object* v___x_1351_; 
v___x_1351_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_741_ = v___x_1351_;
goto v___jp_740_;
}
else
{
lean_object* v___x_1352_; 
v___x_1352_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_741_ = v___x_1352_;
goto v___jp_740_;
}
}
case 19:
{
lean_object* v___x_1353_; uint8_t v___x_1354_; 
v___x_1353_ = lean_unsigned_to_nat(1024u);
v___x_1354_ = lean_nat_dec_le(v___x_1353_, v_prec_666_);
if (v___x_1354_ == 0)
{
lean_object* v___x_1355_; 
v___x_1355_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_734_ = v___x_1355_;
goto v___jp_733_;
}
else
{
lean_object* v___x_1356_; 
v___x_1356_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_734_ = v___x_1356_;
goto v___jp_733_;
}
}
case 20:
{
lean_object* v___x_1357_; uint8_t v___x_1358_; 
v___x_1357_ = lean_unsigned_to_nat(1024u);
v___x_1358_ = lean_nat_dec_le(v___x_1357_, v_prec_666_);
if (v___x_1358_ == 0)
{
lean_object* v___x_1359_; 
v___x_1359_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_727_ = v___x_1359_;
goto v___jp_726_;
}
else
{
lean_object* v___x_1360_; 
v___x_1360_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_727_ = v___x_1360_;
goto v___jp_726_;
}
}
case 21:
{
lean_object* v_ddgThreshold_1361_; lean_object* v___x_1363_; uint8_t v_isShared_1364_; uint8_t v_isSharedCheck_1398_; 
v_ddgThreshold_1361_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1398_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1398_ == 0)
{
v___x_1363_ = v_x_665_;
v_isShared_1364_ = v_isSharedCheck_1398_;
goto v_resetjp_1362_;
}
else
{
lean_inc(v_ddgThreshold_1361_);
lean_dec(v_x_665_);
v___x_1363_ = lean_box(0);
v_isShared_1364_ = v_isSharedCheck_1398_;
goto v_resetjp_1362_;
}
v_resetjp_1362_:
{
lean_object* v___y_1366_; lean_object* v___x_1394_; uint8_t v___x_1395_; 
v___x_1394_ = lean_unsigned_to_nat(1024u);
v___x_1395_ = lean_nat_dec_le(v___x_1394_, v_prec_666_);
if (v___x_1395_ == 0)
{
lean_object* v___x_1396_; 
v___x_1396_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1366_ = v___x_1396_;
goto v___jp_1365_;
}
else
{
lean_object* v___x_1397_; 
v___x_1397_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1366_ = v___x_1397_;
goto v___jp_1365_;
}
v___jp_1365_:
{
lean_object* v_num_1367_; lean_object* v_den_1368_; lean_object* v___x_1369_; lean_object* v___x_1370_; uint8_t v___x_1371_; 
v_num_1367_ = lean_ctor_get(v_ddgThreshold_1361_, 0);
lean_inc(v_num_1367_);
v_den_1368_ = lean_ctor_get(v_ddgThreshold_1361_, 1);
lean_inc(v_den_1368_);
lean_dec_ref(v_ddgThreshold_1361_);
v___x_1369_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__65));
v___x_1370_ = lean_unsigned_to_nat(1u);
v___x_1371_ = lean_nat_dec_eq(v_den_1368_, v___x_1370_);
if (v___x_1371_ == 0)
{
lean_object* v___x_1372_; lean_object* v___x_1373_; lean_object* v___x_1374_; lean_object* v___x_1375_; lean_object* v___x_1376_; lean_object* v___x_1377_; lean_object* v___x_1378_; lean_object* v___x_1380_; 
v___x_1372_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1373_ = l_Int_repr(v_num_1367_);
lean_dec(v_num_1367_);
v___x_1374_ = lean_string_append(v___x_1372_, v___x_1373_);
lean_dec_ref(v___x_1373_);
v___x_1375_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1376_ = lean_string_append(v___x_1374_, v___x_1375_);
v___x_1377_ = l_Nat_reprFast(v_den_1368_);
v___x_1378_ = lean_string_append(v___x_1376_, v___x_1377_);
lean_dec_ref(v___x_1377_);
if (v_isShared_1364_ == 0)
{
lean_ctor_set_tag(v___x_1363_, 3);
lean_ctor_set(v___x_1363_, 0, v___x_1378_);
v___x_1380_ = v___x_1363_;
goto v_reusejp_1379_;
}
else
{
lean_object* v_reuseFailAlloc_1381_; 
v_reuseFailAlloc_1381_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1381_, 0, v___x_1378_);
v___x_1380_ = v_reuseFailAlloc_1381_;
goto v_reusejp_1379_;
}
v_reusejp_1379_:
{
v___y_718_ = v___y_1366_;
v___y_719_ = v___x_1369_;
v___y_720_ = v___x_1380_;
goto v___jp_717_;
}
}
else
{
lean_object* v___x_1382_; lean_object* v___x_1383_; uint8_t v___x_1384_; 
lean_dec(v_den_1368_);
v___x_1382_ = lean_unsigned_to_nat(0u);
v___x_1383_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1384_ = lean_int_dec_lt(v_num_1367_, v___x_1383_);
if (v___x_1384_ == 0)
{
lean_object* v___x_1385_; lean_object* v___x_1387_; 
v___x_1385_ = l_Int_repr(v_num_1367_);
lean_dec(v_num_1367_);
if (v_isShared_1364_ == 0)
{
lean_ctor_set_tag(v___x_1363_, 3);
lean_ctor_set(v___x_1363_, 0, v___x_1385_);
v___x_1387_ = v___x_1363_;
goto v_reusejp_1386_;
}
else
{
lean_object* v_reuseFailAlloc_1388_; 
v_reuseFailAlloc_1388_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1388_, 0, v___x_1385_);
v___x_1387_ = v_reuseFailAlloc_1388_;
goto v_reusejp_1386_;
}
v_reusejp_1386_:
{
v___y_718_ = v___y_1366_;
v___y_719_ = v___x_1369_;
v___y_720_ = v___x_1387_;
goto v___jp_717_;
}
}
else
{
lean_object* v___x_1389_; lean_object* v___x_1391_; 
v___x_1389_ = l_Int_repr(v_num_1367_);
lean_dec(v_num_1367_);
if (v_isShared_1364_ == 0)
{
lean_ctor_set_tag(v___x_1363_, 3);
lean_ctor_set(v___x_1363_, 0, v___x_1389_);
v___x_1391_ = v___x_1363_;
goto v_reusejp_1390_;
}
else
{
lean_object* v_reuseFailAlloc_1393_; 
v_reuseFailAlloc_1393_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1393_, 0, v___x_1389_);
v___x_1391_ = v_reuseFailAlloc_1393_;
goto v_reusejp_1390_;
}
v_reusejp_1390_:
{
lean_object* v___x_1392_; 
v___x_1392_ = l_Repr_addAppParen(v___x_1391_, v___x_1382_);
v___y_718_ = v___y_1366_;
v___y_719_ = v___x_1369_;
v___y_720_ = v___x_1392_;
goto v___jp_717_;
}
}
}
}
}
}
case 22:
{
lean_object* v_maxDDG_1399_; lean_object* v___x_1401_; uint8_t v_isShared_1402_; uint8_t v_isSharedCheck_1436_; 
v_maxDDG_1399_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1436_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1436_ == 0)
{
v___x_1401_ = v_x_665_;
v_isShared_1402_ = v_isSharedCheck_1436_;
goto v_resetjp_1400_;
}
else
{
lean_inc(v_maxDDG_1399_);
lean_dec(v_x_665_);
v___x_1401_ = lean_box(0);
v_isShared_1402_ = v_isSharedCheck_1436_;
goto v_resetjp_1400_;
}
v_resetjp_1400_:
{
lean_object* v___y_1404_; lean_object* v___x_1432_; uint8_t v___x_1433_; 
v___x_1432_ = lean_unsigned_to_nat(1024u);
v___x_1433_ = lean_nat_dec_le(v___x_1432_, v_prec_666_);
if (v___x_1433_ == 0)
{
lean_object* v___x_1434_; 
v___x_1434_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1404_ = v___x_1434_;
goto v___jp_1403_;
}
else
{
lean_object* v___x_1435_; 
v___x_1435_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1404_ = v___x_1435_;
goto v___jp_1403_;
}
v___jp_1403_:
{
lean_object* v_num_1405_; lean_object* v_den_1406_; lean_object* v___x_1407_; lean_object* v___x_1408_; uint8_t v___x_1409_; 
v_num_1405_ = lean_ctor_get(v_maxDDG_1399_, 0);
lean_inc(v_num_1405_);
v_den_1406_ = lean_ctor_get(v_maxDDG_1399_, 1);
lean_inc(v_den_1406_);
lean_dec_ref(v_maxDDG_1399_);
v___x_1407_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__68));
v___x_1408_ = lean_unsigned_to_nat(1u);
v___x_1409_ = lean_nat_dec_eq(v_den_1406_, v___x_1408_);
if (v___x_1409_ == 0)
{
lean_object* v___x_1410_; lean_object* v___x_1411_; lean_object* v___x_1412_; lean_object* v___x_1413_; lean_object* v___x_1414_; lean_object* v___x_1415_; lean_object* v___x_1416_; lean_object* v___x_1418_; 
v___x_1410_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1411_ = l_Int_repr(v_num_1405_);
lean_dec(v_num_1405_);
v___x_1412_ = lean_string_append(v___x_1410_, v___x_1411_);
lean_dec_ref(v___x_1411_);
v___x_1413_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1414_ = lean_string_append(v___x_1412_, v___x_1413_);
v___x_1415_ = l_Nat_reprFast(v_den_1406_);
v___x_1416_ = lean_string_append(v___x_1414_, v___x_1415_);
lean_dec_ref(v___x_1415_);
if (v_isShared_1402_ == 0)
{
lean_ctor_set_tag(v___x_1401_, 3);
lean_ctor_set(v___x_1401_, 0, v___x_1416_);
v___x_1418_ = v___x_1401_;
goto v_reusejp_1417_;
}
else
{
lean_object* v_reuseFailAlloc_1419_; 
v_reuseFailAlloc_1419_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1419_, 0, v___x_1416_);
v___x_1418_ = v_reuseFailAlloc_1419_;
goto v_reusejp_1417_;
}
v_reusejp_1417_:
{
v___y_862_ = v___x_1407_;
v___y_863_ = v___y_1404_;
v___y_864_ = v___x_1418_;
goto v___jp_861_;
}
}
else
{
lean_object* v___x_1420_; lean_object* v___x_1421_; uint8_t v___x_1422_; 
lean_dec(v_den_1406_);
v___x_1420_ = lean_unsigned_to_nat(0u);
v___x_1421_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1422_ = lean_int_dec_lt(v_num_1405_, v___x_1421_);
if (v___x_1422_ == 0)
{
lean_object* v___x_1423_; lean_object* v___x_1425_; 
v___x_1423_ = l_Int_repr(v_num_1405_);
lean_dec(v_num_1405_);
if (v_isShared_1402_ == 0)
{
lean_ctor_set_tag(v___x_1401_, 3);
lean_ctor_set(v___x_1401_, 0, v___x_1423_);
v___x_1425_ = v___x_1401_;
goto v_reusejp_1424_;
}
else
{
lean_object* v_reuseFailAlloc_1426_; 
v_reuseFailAlloc_1426_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1426_, 0, v___x_1423_);
v___x_1425_ = v_reuseFailAlloc_1426_;
goto v_reusejp_1424_;
}
v_reusejp_1424_:
{
v___y_862_ = v___x_1407_;
v___y_863_ = v___y_1404_;
v___y_864_ = v___x_1425_;
goto v___jp_861_;
}
}
else
{
lean_object* v___x_1427_; lean_object* v___x_1429_; 
v___x_1427_ = l_Int_repr(v_num_1405_);
lean_dec(v_num_1405_);
if (v_isShared_1402_ == 0)
{
lean_ctor_set_tag(v___x_1401_, 3);
lean_ctor_set(v___x_1401_, 0, v___x_1427_);
v___x_1429_ = v___x_1401_;
goto v_reusejp_1428_;
}
else
{
lean_object* v_reuseFailAlloc_1431_; 
v_reuseFailAlloc_1431_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1431_, 0, v___x_1427_);
v___x_1429_ = v_reuseFailAlloc_1431_;
goto v_reusejp_1428_;
}
v_reusejp_1428_:
{
lean_object* v___x_1430_; 
v___x_1430_ = l_Repr_addAppParen(v___x_1429_, v___x_1420_);
v___y_862_ = v___x_1407_;
v___y_863_ = v___y_1404_;
v___y_864_ = v___x_1430_;
goto v___jp_861_;
}
}
}
}
}
}
case 23:
{
lean_object* v___x_1437_; uint8_t v___x_1438_; 
v___x_1437_ = lean_unsigned_to_nat(1024u);
v___x_1438_ = lean_nat_dec_le(v___x_1437_, v_prec_666_);
if (v___x_1438_ == 0)
{
lean_object* v___x_1439_; 
v___x_1439_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_711_ = v___x_1439_;
goto v___jp_710_;
}
else
{
lean_object* v___x_1440_; 
v___x_1440_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_711_ = v___x_1440_;
goto v___jp_710_;
}
}
case 24:
{
lean_object* v_threshold_1441_; lean_object* v___x_1443_; uint8_t v_isShared_1444_; uint8_t v_isSharedCheck_1478_; 
v_threshold_1441_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1478_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1478_ == 0)
{
v___x_1443_ = v_x_665_;
v_isShared_1444_ = v_isSharedCheck_1478_;
goto v_resetjp_1442_;
}
else
{
lean_inc(v_threshold_1441_);
lean_dec(v_x_665_);
v___x_1443_ = lean_box(0);
v_isShared_1444_ = v_isSharedCheck_1478_;
goto v_resetjp_1442_;
}
v_resetjp_1442_:
{
lean_object* v___y_1446_; lean_object* v___x_1474_; uint8_t v___x_1475_; 
v___x_1474_ = lean_unsigned_to_nat(1024u);
v___x_1475_ = lean_nat_dec_le(v___x_1474_, v_prec_666_);
if (v___x_1475_ == 0)
{
lean_object* v___x_1476_; 
v___x_1476_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1446_ = v___x_1476_;
goto v___jp_1445_;
}
else
{
lean_object* v___x_1477_; 
v___x_1477_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1446_ = v___x_1477_;
goto v___jp_1445_;
}
v___jp_1445_:
{
lean_object* v_num_1447_; lean_object* v_den_1448_; lean_object* v___x_1449_; lean_object* v___x_1450_; uint8_t v___x_1451_; 
v_num_1447_ = lean_ctor_get(v_threshold_1441_, 0);
lean_inc(v_num_1447_);
v_den_1448_ = lean_ctor_get(v_threshold_1441_, 1);
lean_inc(v_den_1448_);
lean_dec_ref(v_threshold_1441_);
v___x_1449_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__71));
v___x_1450_ = lean_unsigned_to_nat(1u);
v___x_1451_ = lean_nat_dec_eq(v_den_1448_, v___x_1450_);
if (v___x_1451_ == 0)
{
lean_object* v___x_1452_; lean_object* v___x_1453_; lean_object* v___x_1454_; lean_object* v___x_1455_; lean_object* v___x_1456_; lean_object* v___x_1457_; lean_object* v___x_1458_; lean_object* v___x_1460_; 
v___x_1452_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1453_ = l_Int_repr(v_num_1447_);
lean_dec(v_num_1447_);
v___x_1454_ = lean_string_append(v___x_1452_, v___x_1453_);
lean_dec_ref(v___x_1453_);
v___x_1455_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1456_ = lean_string_append(v___x_1454_, v___x_1455_);
v___x_1457_ = l_Nat_reprFast(v_den_1448_);
v___x_1458_ = lean_string_append(v___x_1456_, v___x_1457_);
lean_dec_ref(v___x_1457_);
if (v_isShared_1444_ == 0)
{
lean_ctor_set_tag(v___x_1443_, 3);
lean_ctor_set(v___x_1443_, 0, v___x_1458_);
v___x_1460_ = v___x_1443_;
goto v_reusejp_1459_;
}
else
{
lean_object* v_reuseFailAlloc_1461_; 
v_reuseFailAlloc_1461_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1461_, 0, v___x_1458_);
v___x_1460_ = v_reuseFailAlloc_1461_;
goto v_reusejp_1459_;
}
v_reusejp_1459_:
{
v___y_702_ = v___y_1446_;
v___y_703_ = v___x_1449_;
v___y_704_ = v___x_1460_;
goto v___jp_701_;
}
}
else
{
lean_object* v___x_1462_; lean_object* v___x_1463_; uint8_t v___x_1464_; 
lean_dec(v_den_1448_);
v___x_1462_ = lean_unsigned_to_nat(0u);
v___x_1463_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1464_ = lean_int_dec_lt(v_num_1447_, v___x_1463_);
if (v___x_1464_ == 0)
{
lean_object* v___x_1465_; lean_object* v___x_1467_; 
v___x_1465_ = l_Int_repr(v_num_1447_);
lean_dec(v_num_1447_);
if (v_isShared_1444_ == 0)
{
lean_ctor_set_tag(v___x_1443_, 3);
lean_ctor_set(v___x_1443_, 0, v___x_1465_);
v___x_1467_ = v___x_1443_;
goto v_reusejp_1466_;
}
else
{
lean_object* v_reuseFailAlloc_1468_; 
v_reuseFailAlloc_1468_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1468_, 0, v___x_1465_);
v___x_1467_ = v_reuseFailAlloc_1468_;
goto v_reusejp_1466_;
}
v_reusejp_1466_:
{
v___y_702_ = v___y_1446_;
v___y_703_ = v___x_1449_;
v___y_704_ = v___x_1467_;
goto v___jp_701_;
}
}
else
{
lean_object* v___x_1469_; lean_object* v___x_1471_; 
v___x_1469_ = l_Int_repr(v_num_1447_);
lean_dec(v_num_1447_);
if (v_isShared_1444_ == 0)
{
lean_ctor_set_tag(v___x_1443_, 3);
lean_ctor_set(v___x_1443_, 0, v___x_1469_);
v___x_1471_ = v___x_1443_;
goto v_reusejp_1470_;
}
else
{
lean_object* v_reuseFailAlloc_1473_; 
v_reuseFailAlloc_1473_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1473_, 0, v___x_1469_);
v___x_1471_ = v_reuseFailAlloc_1473_;
goto v_reusejp_1470_;
}
v_reusejp_1470_:
{
lean_object* v___x_1472_; 
v___x_1472_ = l_Repr_addAppParen(v___x_1471_, v___x_1462_);
v___y_702_ = v___y_1446_;
v___y_703_ = v___x_1449_;
v___y_704_ = v___x_1472_;
goto v___jp_701_;
}
}
}
}
}
}
case 25:
{
lean_object* v_minScore_1479_; lean_object* v___x_1481_; uint8_t v_isShared_1482_; uint8_t v_isSharedCheck_1516_; 
v_minScore_1479_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1516_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1516_ == 0)
{
v___x_1481_ = v_x_665_;
v_isShared_1482_ = v_isSharedCheck_1516_;
goto v_resetjp_1480_;
}
else
{
lean_inc(v_minScore_1479_);
lean_dec(v_x_665_);
v___x_1481_ = lean_box(0);
v_isShared_1482_ = v_isSharedCheck_1516_;
goto v_resetjp_1480_;
}
v_resetjp_1480_:
{
lean_object* v___y_1484_; lean_object* v___x_1512_; uint8_t v___x_1513_; 
v___x_1512_ = lean_unsigned_to_nat(1024u);
v___x_1513_ = lean_nat_dec_le(v___x_1512_, v_prec_666_);
if (v___x_1513_ == 0)
{
lean_object* v___x_1514_; 
v___x_1514_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1484_ = v___x_1514_;
goto v___jp_1483_;
}
else
{
lean_object* v___x_1515_; 
v___x_1515_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1484_ = v___x_1515_;
goto v___jp_1483_;
}
v___jp_1483_:
{
lean_object* v_num_1485_; lean_object* v_den_1486_; lean_object* v___x_1487_; lean_object* v___x_1488_; uint8_t v___x_1489_; 
v_num_1485_ = lean_ctor_get(v_minScore_1479_, 0);
lean_inc(v_num_1485_);
v_den_1486_ = lean_ctor_get(v_minScore_1479_, 1);
lean_inc(v_den_1486_);
lean_dec_ref(v_minScore_1479_);
v___x_1487_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__74));
v___x_1488_ = lean_unsigned_to_nat(1u);
v___x_1489_ = lean_nat_dec_eq(v_den_1486_, v___x_1488_);
if (v___x_1489_ == 0)
{
lean_object* v___x_1490_; lean_object* v___x_1491_; lean_object* v___x_1492_; lean_object* v___x_1493_; lean_object* v___x_1494_; lean_object* v___x_1495_; lean_object* v___x_1496_; lean_object* v___x_1498_; 
v___x_1490_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1491_ = l_Int_repr(v_num_1485_);
lean_dec(v_num_1485_);
v___x_1492_ = lean_string_append(v___x_1490_, v___x_1491_);
lean_dec_ref(v___x_1491_);
v___x_1493_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1494_ = lean_string_append(v___x_1492_, v___x_1493_);
v___x_1495_ = l_Nat_reprFast(v_den_1486_);
v___x_1496_ = lean_string_append(v___x_1494_, v___x_1495_);
lean_dec_ref(v___x_1495_);
if (v_isShared_1482_ == 0)
{
lean_ctor_set_tag(v___x_1481_, 3);
lean_ctor_set(v___x_1481_, 0, v___x_1496_);
v___x_1498_ = v___x_1481_;
goto v_reusejp_1497_;
}
else
{
lean_object* v_reuseFailAlloc_1499_; 
v_reuseFailAlloc_1499_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1499_, 0, v___x_1496_);
v___x_1498_ = v_reuseFailAlloc_1499_;
goto v_reusejp_1497_;
}
v_reusejp_1497_:
{
v___y_871_ = v___x_1487_;
v___y_872_ = v___y_1484_;
v___y_873_ = v___x_1498_;
goto v___jp_870_;
}
}
else
{
lean_object* v___x_1500_; lean_object* v___x_1501_; uint8_t v___x_1502_; 
lean_dec(v_den_1486_);
v___x_1500_ = lean_unsigned_to_nat(0u);
v___x_1501_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1502_ = lean_int_dec_lt(v_num_1485_, v___x_1501_);
if (v___x_1502_ == 0)
{
lean_object* v___x_1503_; lean_object* v___x_1505_; 
v___x_1503_ = l_Int_repr(v_num_1485_);
lean_dec(v_num_1485_);
if (v_isShared_1482_ == 0)
{
lean_ctor_set_tag(v___x_1481_, 3);
lean_ctor_set(v___x_1481_, 0, v___x_1503_);
v___x_1505_ = v___x_1481_;
goto v_reusejp_1504_;
}
else
{
lean_object* v_reuseFailAlloc_1506_; 
v_reuseFailAlloc_1506_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1506_, 0, v___x_1503_);
v___x_1505_ = v_reuseFailAlloc_1506_;
goto v_reusejp_1504_;
}
v_reusejp_1504_:
{
v___y_871_ = v___x_1487_;
v___y_872_ = v___y_1484_;
v___y_873_ = v___x_1505_;
goto v___jp_870_;
}
}
else
{
lean_object* v___x_1507_; lean_object* v___x_1509_; 
v___x_1507_ = l_Int_repr(v_num_1485_);
lean_dec(v_num_1485_);
if (v_isShared_1482_ == 0)
{
lean_ctor_set_tag(v___x_1481_, 3);
lean_ctor_set(v___x_1481_, 0, v___x_1507_);
v___x_1509_ = v___x_1481_;
goto v_reusejp_1508_;
}
else
{
lean_object* v_reuseFailAlloc_1511_; 
v_reuseFailAlloc_1511_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1511_, 0, v___x_1507_);
v___x_1509_ = v_reuseFailAlloc_1511_;
goto v_reusejp_1508_;
}
v_reusejp_1508_:
{
lean_object* v___x_1510_; 
v___x_1510_ = l_Repr_addAppParen(v___x_1509_, v___x_1500_);
v___y_871_ = v___x_1487_;
v___y_872_ = v___y_1484_;
v___y_873_ = v___x_1510_;
goto v___jp_870_;
}
}
}
}
}
}
case 26:
{
lean_object* v___x_1517_; uint8_t v___x_1518_; 
v___x_1517_ = lean_unsigned_to_nat(1024u);
v___x_1518_ = lean_nat_dec_le(v___x_1517_, v_prec_666_);
if (v___x_1518_ == 0)
{
lean_object* v___x_1519_; 
v___x_1519_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_695_ = v___x_1519_;
goto v___jp_694_;
}
else
{
lean_object* v___x_1520_; 
v___x_1520_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_695_ = v___x_1520_;
goto v___jp_694_;
}
}
case 27:
{
lean_object* v_pILo_1521_; lean_object* v_pIHi_1522_; lean_object* v___x_1524_; uint8_t v_isShared_1525_; uint8_t v_isSharedCheck_1589_; 
v_pILo_1521_ = lean_ctor_get(v_x_665_, 0);
v_pIHi_1522_ = lean_ctor_get(v_x_665_, 1);
v_isSharedCheck_1589_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1589_ == 0)
{
v___x_1524_ = v_x_665_;
v_isShared_1525_ = v_isSharedCheck_1589_;
goto v_resetjp_1523_;
}
else
{
lean_inc(v_pIHi_1522_);
lean_inc(v_pILo_1521_);
lean_dec(v_x_665_);
v___x_1524_ = lean_box(0);
v_isShared_1525_ = v_isSharedCheck_1589_;
goto v_resetjp_1523_;
}
v_resetjp_1523_:
{
lean_object* v___y_1527_; lean_object* v___y_1528_; lean_object* v___y_1529_; lean_object* v___y_1530_; lean_object* v___y_1562_; lean_object* v___x_1585_; uint8_t v___x_1586_; 
v___x_1585_ = lean_unsigned_to_nat(1024u);
v___x_1586_ = lean_nat_dec_le(v___x_1585_, v_prec_666_);
if (v___x_1586_ == 0)
{
lean_object* v___x_1587_; 
v___x_1587_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1562_ = v___x_1587_;
goto v___jp_1561_;
}
else
{
lean_object* v___x_1588_; 
v___x_1588_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1562_ = v___x_1588_;
goto v___jp_1561_;
}
v___jp_1526_:
{
lean_object* v_num_1531_; lean_object* v_den_1532_; lean_object* v___x_1534_; uint8_t v_isShared_1535_; uint8_t v_isSharedCheck_1560_; 
v_num_1531_ = lean_ctor_get(v_pIHi_1522_, 0);
v_den_1532_ = lean_ctor_get(v_pIHi_1522_, 1);
v_isSharedCheck_1560_ = !lean_is_exclusive(v_pIHi_1522_);
if (v_isSharedCheck_1560_ == 0)
{
v___x_1534_ = v_pIHi_1522_;
v_isShared_1535_ = v_isSharedCheck_1560_;
goto v_resetjp_1533_;
}
else
{
lean_inc(v_den_1532_);
lean_inc(v_num_1531_);
lean_dec(v_pIHi_1522_);
v___x_1534_ = lean_box(0);
v_isShared_1535_ = v_isSharedCheck_1560_;
goto v_resetjp_1533_;
}
v_resetjp_1533_:
{
lean_object* v___x_1537_; 
lean_inc(v___y_1529_);
if (v_isShared_1535_ == 0)
{
lean_ctor_set_tag(v___x_1534_, 5);
lean_ctor_set(v___x_1534_, 1, v___y_1530_);
lean_ctor_set(v___x_1534_, 0, v___y_1529_);
v___x_1537_ = v___x_1534_;
goto v_reusejp_1536_;
}
else
{
lean_object* v_reuseFailAlloc_1559_; 
v_reuseFailAlloc_1559_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1559_, 0, v___y_1529_);
lean_ctor_set(v_reuseFailAlloc_1559_, 1, v___y_1530_);
v___x_1537_ = v_reuseFailAlloc_1559_;
goto v_reusejp_1536_;
}
v_reusejp_1536_:
{
lean_object* v___x_1539_; 
lean_inc(v___y_1528_);
if (v_isShared_1525_ == 0)
{
lean_ctor_set_tag(v___x_1524_, 5);
lean_ctor_set(v___x_1524_, 1, v___y_1528_);
lean_ctor_set(v___x_1524_, 0, v___x_1537_);
v___x_1539_ = v___x_1524_;
goto v_reusejp_1538_;
}
else
{
lean_object* v_reuseFailAlloc_1558_; 
v_reuseFailAlloc_1558_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_1558_, 0, v___x_1537_);
lean_ctor_set(v_reuseFailAlloc_1558_, 1, v___y_1528_);
v___x_1539_ = v_reuseFailAlloc_1558_;
goto v_reusejp_1538_;
}
v_reusejp_1538_:
{
lean_object* v___x_1540_; uint8_t v___x_1541_; 
v___x_1540_ = lean_unsigned_to_nat(1u);
v___x_1541_ = lean_nat_dec_eq(v_den_1532_, v___x_1540_);
if (v___x_1541_ == 0)
{
lean_object* v___x_1542_; lean_object* v___x_1543_; lean_object* v___x_1544_; lean_object* v___x_1545_; lean_object* v___x_1546_; lean_object* v___x_1547_; lean_object* v___x_1548_; lean_object* v___x_1549_; 
v___x_1542_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1543_ = l_Int_repr(v_num_1531_);
lean_dec(v_num_1531_);
v___x_1544_ = lean_string_append(v___x_1542_, v___x_1543_);
lean_dec_ref(v___x_1543_);
v___x_1545_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1546_ = lean_string_append(v___x_1544_, v___x_1545_);
v___x_1547_ = l_Nat_reprFast(v_den_1532_);
v___x_1548_ = lean_string_append(v___x_1546_, v___x_1547_);
lean_dec_ref(v___x_1547_);
v___x_1549_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1549_, 0, v___x_1548_);
v___y_686_ = v___x_1539_;
v___y_687_ = v___y_1527_;
v___y_688_ = v___x_1549_;
goto v___jp_685_;
}
else
{
lean_object* v___x_1550_; lean_object* v___x_1551_; uint8_t v___x_1552_; 
lean_dec(v_den_1532_);
v___x_1550_ = lean_unsigned_to_nat(0u);
v___x_1551_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1552_ = lean_int_dec_lt(v_num_1531_, v___x_1551_);
if (v___x_1552_ == 0)
{
lean_object* v___x_1553_; lean_object* v___x_1554_; 
v___x_1553_ = l_Int_repr(v_num_1531_);
lean_dec(v_num_1531_);
v___x_1554_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1554_, 0, v___x_1553_);
v___y_686_ = v___x_1539_;
v___y_687_ = v___y_1527_;
v___y_688_ = v___x_1554_;
goto v___jp_685_;
}
else
{
lean_object* v___x_1555_; lean_object* v___x_1556_; lean_object* v___x_1557_; 
v___x_1555_ = l_Int_repr(v_num_1531_);
lean_dec(v_num_1531_);
v___x_1556_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1556_, 0, v___x_1555_);
v___x_1557_ = l_Repr_addAppParen(v___x_1556_, v___x_1550_);
v___y_686_ = v___x_1539_;
v___y_687_ = v___y_1527_;
v___y_688_ = v___x_1557_;
goto v___jp_685_;
}
}
}
}
}
}
v___jp_1561_:
{
lean_object* v_num_1563_; lean_object* v_den_1564_; lean_object* v___x_1565_; lean_object* v___x_1566_; lean_object* v___x_1567_; uint8_t v___x_1568_; 
v_num_1563_ = lean_ctor_get(v_pILo_1521_, 0);
lean_inc(v_num_1563_);
v_den_1564_ = lean_ctor_get(v_pILo_1521_, 1);
lean_inc(v_den_1564_);
lean_dec_ref(v_pILo_1521_);
v___x_1565_ = lean_box(1);
v___x_1566_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__77));
v___x_1567_ = lean_unsigned_to_nat(1u);
v___x_1568_ = lean_nat_dec_eq(v_den_1564_, v___x_1567_);
if (v___x_1568_ == 0)
{
lean_object* v___x_1569_; lean_object* v___x_1570_; lean_object* v___x_1571_; lean_object* v___x_1572_; lean_object* v___x_1573_; lean_object* v___x_1574_; lean_object* v___x_1575_; lean_object* v___x_1576_; 
v___x_1569_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1570_ = l_Int_repr(v_num_1563_);
lean_dec(v_num_1563_);
v___x_1571_ = lean_string_append(v___x_1569_, v___x_1570_);
lean_dec_ref(v___x_1570_);
v___x_1572_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1573_ = lean_string_append(v___x_1571_, v___x_1572_);
v___x_1574_ = l_Nat_reprFast(v_den_1564_);
v___x_1575_ = lean_string_append(v___x_1573_, v___x_1574_);
lean_dec_ref(v___x_1574_);
v___x_1576_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1576_, 0, v___x_1575_);
v___y_1527_ = v___y_1562_;
v___y_1528_ = v___x_1565_;
v___y_1529_ = v___x_1566_;
v___y_1530_ = v___x_1576_;
goto v___jp_1526_;
}
else
{
lean_object* v___x_1577_; lean_object* v___x_1578_; uint8_t v___x_1579_; 
lean_dec(v_den_1564_);
v___x_1577_ = lean_unsigned_to_nat(0u);
v___x_1578_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1579_ = lean_int_dec_lt(v_num_1563_, v___x_1578_);
if (v___x_1579_ == 0)
{
lean_object* v___x_1580_; lean_object* v___x_1581_; 
v___x_1580_ = l_Int_repr(v_num_1563_);
lean_dec(v_num_1563_);
v___x_1581_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1581_, 0, v___x_1580_);
v___y_1527_ = v___y_1562_;
v___y_1528_ = v___x_1565_;
v___y_1529_ = v___x_1566_;
v___y_1530_ = v___x_1581_;
goto v___jp_1526_;
}
else
{
lean_object* v___x_1582_; lean_object* v___x_1583_; lean_object* v___x_1584_; 
v___x_1582_ = l_Int_repr(v_num_1563_);
lean_dec(v_num_1563_);
v___x_1583_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_1583_, 0, v___x_1582_);
v___x_1584_ = l_Repr_addAppParen(v___x_1583_, v___x_1577_);
v___y_1527_ = v___y_1562_;
v___y_1528_ = v___x_1565_;
v___y_1529_ = v___x_1566_;
v___y_1530_ = v___x_1584_;
goto v___jp_1526_;
}
}
}
}
}
case 28:
{
lean_object* v_maxLen_1590_; lean_object* v___x_1592_; uint8_t v_isShared_1593_; uint8_t v_isSharedCheck_1610_; 
v_maxLen_1590_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1610_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1610_ == 0)
{
v___x_1592_ = v_x_665_;
v_isShared_1593_ = v_isSharedCheck_1610_;
goto v_resetjp_1591_;
}
else
{
lean_inc(v_maxLen_1590_);
lean_dec(v_x_665_);
v___x_1592_ = lean_box(0);
v_isShared_1593_ = v_isSharedCheck_1610_;
goto v_resetjp_1591_;
}
v_resetjp_1591_:
{
lean_object* v___y_1595_; lean_object* v___x_1606_; uint8_t v___x_1607_; 
v___x_1606_ = lean_unsigned_to_nat(1024u);
v___x_1607_ = lean_nat_dec_le(v___x_1606_, v_prec_666_);
if (v___x_1607_ == 0)
{
lean_object* v___x_1608_; 
v___x_1608_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1595_ = v___x_1608_;
goto v___jp_1594_;
}
else
{
lean_object* v___x_1609_; 
v___x_1609_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1595_ = v___x_1609_;
goto v___jp_1594_;
}
v___jp_1594_:
{
lean_object* v___x_1596_; lean_object* v___x_1597_; lean_object* v___x_1599_; 
v___x_1596_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__80));
v___x_1597_ = l_Nat_reprFast(v_maxLen_1590_);
if (v_isShared_1593_ == 0)
{
lean_ctor_set_tag(v___x_1592_, 3);
lean_ctor_set(v___x_1592_, 0, v___x_1597_);
v___x_1599_ = v___x_1592_;
goto v_reusejp_1598_;
}
else
{
lean_object* v_reuseFailAlloc_1605_; 
v_reuseFailAlloc_1605_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1605_, 0, v___x_1597_);
v___x_1599_ = v_reuseFailAlloc_1605_;
goto v_reusejp_1598_;
}
v_reusejp_1598_:
{
lean_object* v___x_1600_; lean_object* v___x_1601_; uint8_t v___x_1602_; lean_object* v___x_1603_; lean_object* v___x_1604_; 
v___x_1600_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_1600_, 0, v___x_1596_);
lean_ctor_set(v___x_1600_, 1, v___x_1599_);
lean_inc(v___y_1595_);
v___x_1601_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_1601_, 0, v___y_1595_);
lean_ctor_set(v___x_1601_, 1, v___x_1600_);
v___x_1602_ = 0;
v___x_1603_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_1603_, 0, v___x_1601_);
lean_ctor_set_uint8(v___x_1603_, sizeof(void*)*1, v___x_1602_);
v___x_1604_ = l_Repr_addAppParen(v___x_1603_, v_prec_666_);
return v___x_1604_;
}
}
}
}
case 29:
{
lean_object* v_maxScore_1611_; lean_object* v___x_1613_; uint8_t v_isShared_1614_; uint8_t v_isSharedCheck_1648_; 
v_maxScore_1611_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1648_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1648_ == 0)
{
v___x_1613_ = v_x_665_;
v_isShared_1614_ = v_isSharedCheck_1648_;
goto v_resetjp_1612_;
}
else
{
lean_inc(v_maxScore_1611_);
lean_dec(v_x_665_);
v___x_1613_ = lean_box(0);
v_isShared_1614_ = v_isSharedCheck_1648_;
goto v_resetjp_1612_;
}
v_resetjp_1612_:
{
lean_object* v___y_1616_; lean_object* v___x_1644_; uint8_t v___x_1645_; 
v___x_1644_ = lean_unsigned_to_nat(1024u);
v___x_1645_ = lean_nat_dec_le(v___x_1644_, v_prec_666_);
if (v___x_1645_ == 0)
{
lean_object* v___x_1646_; 
v___x_1646_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1616_ = v___x_1646_;
goto v___jp_1615_;
}
else
{
lean_object* v___x_1647_; 
v___x_1647_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1616_ = v___x_1647_;
goto v___jp_1615_;
}
v___jp_1615_:
{
lean_object* v_num_1617_; lean_object* v_den_1618_; lean_object* v___x_1619_; lean_object* v___x_1620_; uint8_t v___x_1621_; 
v_num_1617_ = lean_ctor_get(v_maxScore_1611_, 0);
lean_inc(v_num_1617_);
v_den_1618_ = lean_ctor_get(v_maxScore_1611_, 1);
lean_inc(v_den_1618_);
lean_dec_ref(v_maxScore_1611_);
v___x_1619_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__83));
v___x_1620_ = lean_unsigned_to_nat(1u);
v___x_1621_ = lean_nat_dec_eq(v_den_1618_, v___x_1620_);
if (v___x_1621_ == 0)
{
lean_object* v___x_1622_; lean_object* v___x_1623_; lean_object* v___x_1624_; lean_object* v___x_1625_; lean_object* v___x_1626_; lean_object* v___x_1627_; lean_object* v___x_1628_; lean_object* v___x_1630_; 
v___x_1622_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1623_ = l_Int_repr(v_num_1617_);
lean_dec(v_num_1617_);
v___x_1624_ = lean_string_append(v___x_1622_, v___x_1623_);
lean_dec_ref(v___x_1623_);
v___x_1625_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1626_ = lean_string_append(v___x_1624_, v___x_1625_);
v___x_1627_ = l_Nat_reprFast(v_den_1618_);
v___x_1628_ = lean_string_append(v___x_1626_, v___x_1627_);
lean_dec_ref(v___x_1627_);
if (v_isShared_1614_ == 0)
{
lean_ctor_set_tag(v___x_1613_, 3);
lean_ctor_set(v___x_1613_, 0, v___x_1628_);
v___x_1630_ = v___x_1613_;
goto v_reusejp_1629_;
}
else
{
lean_object* v_reuseFailAlloc_1631_; 
v_reuseFailAlloc_1631_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1631_, 0, v___x_1628_);
v___x_1630_ = v_reuseFailAlloc_1631_;
goto v_reusejp_1629_;
}
v_reusejp_1629_:
{
v___y_880_ = v___y_1616_;
v___y_881_ = v___x_1619_;
v___y_882_ = v___x_1630_;
goto v___jp_879_;
}
}
else
{
lean_object* v___x_1632_; lean_object* v___x_1633_; uint8_t v___x_1634_; 
lean_dec(v_den_1618_);
v___x_1632_ = lean_unsigned_to_nat(0u);
v___x_1633_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1634_ = lean_int_dec_lt(v_num_1617_, v___x_1633_);
if (v___x_1634_ == 0)
{
lean_object* v___x_1635_; lean_object* v___x_1637_; 
v___x_1635_ = l_Int_repr(v_num_1617_);
lean_dec(v_num_1617_);
if (v_isShared_1614_ == 0)
{
lean_ctor_set_tag(v___x_1613_, 3);
lean_ctor_set(v___x_1613_, 0, v___x_1635_);
v___x_1637_ = v___x_1613_;
goto v_reusejp_1636_;
}
else
{
lean_object* v_reuseFailAlloc_1638_; 
v_reuseFailAlloc_1638_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1638_, 0, v___x_1635_);
v___x_1637_ = v_reuseFailAlloc_1638_;
goto v_reusejp_1636_;
}
v_reusejp_1636_:
{
v___y_880_ = v___y_1616_;
v___y_881_ = v___x_1619_;
v___y_882_ = v___x_1637_;
goto v___jp_879_;
}
}
else
{
lean_object* v___x_1639_; lean_object* v___x_1641_; 
v___x_1639_ = l_Int_repr(v_num_1617_);
lean_dec(v_num_1617_);
if (v_isShared_1614_ == 0)
{
lean_ctor_set_tag(v___x_1613_, 3);
lean_ctor_set(v___x_1613_, 0, v___x_1639_);
v___x_1641_ = v___x_1613_;
goto v_reusejp_1640_;
}
else
{
lean_object* v_reuseFailAlloc_1643_; 
v_reuseFailAlloc_1643_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1643_, 0, v___x_1639_);
v___x_1641_ = v_reuseFailAlloc_1643_;
goto v_reusejp_1640_;
}
v_reusejp_1640_:
{
lean_object* v___x_1642_; 
v___x_1642_ = l_Repr_addAppParen(v___x_1641_, v___x_1632_);
v___y_880_ = v___y_1616_;
v___y_881_ = v___x_1619_;
v___y_882_ = v___x_1642_;
goto v___jp_879_;
}
}
}
}
}
}
case 30:
{
lean_object* v_ic50Threshold_1649_; lean_object* v___x_1651_; uint8_t v_isShared_1652_; uint8_t v_isSharedCheck_1686_; 
v_ic50Threshold_1649_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1686_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1686_ == 0)
{
v___x_1651_ = v_x_665_;
v_isShared_1652_ = v_isSharedCheck_1686_;
goto v_resetjp_1650_;
}
else
{
lean_inc(v_ic50Threshold_1649_);
lean_dec(v_x_665_);
v___x_1651_ = lean_box(0);
v_isShared_1652_ = v_isSharedCheck_1686_;
goto v_resetjp_1650_;
}
v_resetjp_1650_:
{
lean_object* v___y_1654_; lean_object* v___x_1682_; uint8_t v___x_1683_; 
v___x_1682_ = lean_unsigned_to_nat(1024u);
v___x_1683_ = lean_nat_dec_le(v___x_1682_, v_prec_666_);
if (v___x_1683_ == 0)
{
lean_object* v___x_1684_; 
v___x_1684_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1654_ = v___x_1684_;
goto v___jp_1653_;
}
else
{
lean_object* v___x_1685_; 
v___x_1685_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1654_ = v___x_1685_;
goto v___jp_1653_;
}
v___jp_1653_:
{
lean_object* v_num_1655_; lean_object* v_den_1656_; lean_object* v___x_1657_; lean_object* v___x_1658_; uint8_t v___x_1659_; 
v_num_1655_ = lean_ctor_get(v_ic50Threshold_1649_, 0);
lean_inc(v_num_1655_);
v_den_1656_ = lean_ctor_get(v_ic50Threshold_1649_, 1);
lean_inc(v_den_1656_);
lean_dec_ref(v_ic50Threshold_1649_);
v___x_1657_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__86));
v___x_1658_ = lean_unsigned_to_nat(1u);
v___x_1659_ = lean_nat_dec_eq(v_den_1656_, v___x_1658_);
if (v___x_1659_ == 0)
{
lean_object* v___x_1660_; lean_object* v___x_1661_; lean_object* v___x_1662_; lean_object* v___x_1663_; lean_object* v___x_1664_; lean_object* v___x_1665_; lean_object* v___x_1666_; lean_object* v___x_1668_; 
v___x_1660_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1661_ = l_Int_repr(v_num_1655_);
lean_dec(v_num_1655_);
v___x_1662_ = lean_string_append(v___x_1660_, v___x_1661_);
lean_dec_ref(v___x_1661_);
v___x_1663_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1664_ = lean_string_append(v___x_1662_, v___x_1663_);
v___x_1665_ = l_Nat_reprFast(v_den_1656_);
v___x_1666_ = lean_string_append(v___x_1664_, v___x_1665_);
lean_dec_ref(v___x_1665_);
if (v_isShared_1652_ == 0)
{
lean_ctor_set_tag(v___x_1651_, 3);
lean_ctor_set(v___x_1651_, 0, v___x_1666_);
v___x_1668_ = v___x_1651_;
goto v_reusejp_1667_;
}
else
{
lean_object* v_reuseFailAlloc_1669_; 
v_reuseFailAlloc_1669_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1669_, 0, v___x_1666_);
v___x_1668_ = v_reuseFailAlloc_1669_;
goto v_reusejp_1667_;
}
v_reusejp_1667_:
{
v___y_677_ = v___y_1654_;
v___y_678_ = v___x_1657_;
v___y_679_ = v___x_1668_;
goto v___jp_676_;
}
}
else
{
lean_object* v___x_1670_; lean_object* v___x_1671_; uint8_t v___x_1672_; 
lean_dec(v_den_1656_);
v___x_1670_ = lean_unsigned_to_nat(0u);
v___x_1671_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1672_ = lean_int_dec_lt(v_num_1655_, v___x_1671_);
if (v___x_1672_ == 0)
{
lean_object* v___x_1673_; lean_object* v___x_1675_; 
v___x_1673_ = l_Int_repr(v_num_1655_);
lean_dec(v_num_1655_);
if (v_isShared_1652_ == 0)
{
lean_ctor_set_tag(v___x_1651_, 3);
lean_ctor_set(v___x_1651_, 0, v___x_1673_);
v___x_1675_ = v___x_1651_;
goto v_reusejp_1674_;
}
else
{
lean_object* v_reuseFailAlloc_1676_; 
v_reuseFailAlloc_1676_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1676_, 0, v___x_1673_);
v___x_1675_ = v_reuseFailAlloc_1676_;
goto v_reusejp_1674_;
}
v_reusejp_1674_:
{
v___y_677_ = v___y_1654_;
v___y_678_ = v___x_1657_;
v___y_679_ = v___x_1675_;
goto v___jp_676_;
}
}
else
{
lean_object* v___x_1677_; lean_object* v___x_1679_; 
v___x_1677_ = l_Int_repr(v_num_1655_);
lean_dec(v_num_1655_);
if (v_isShared_1652_ == 0)
{
lean_ctor_set_tag(v___x_1651_, 3);
lean_ctor_set(v___x_1651_, 0, v___x_1677_);
v___x_1679_ = v___x_1651_;
goto v_reusejp_1678_;
}
else
{
lean_object* v_reuseFailAlloc_1681_; 
v_reuseFailAlloc_1681_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1681_, 0, v___x_1677_);
v___x_1679_ = v_reuseFailAlloc_1681_;
goto v_reusejp_1678_;
}
v_reusejp_1678_:
{
lean_object* v___x_1680_; 
v___x_1680_ = l_Repr_addAppParen(v___x_1679_, v___x_1670_);
v___y_677_ = v___y_1654_;
v___y_678_ = v___x_1657_;
v___y_679_ = v___x_1680_;
goto v___jp_676_;
}
}
}
}
}
}
case 31:
{
lean_object* v_scoreThreshold_1687_; lean_object* v___x_1689_; uint8_t v_isShared_1690_; uint8_t v_isSharedCheck_1724_; 
v_scoreThreshold_1687_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1724_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1724_ == 0)
{
v___x_1689_ = v_x_665_;
v_isShared_1690_ = v_isSharedCheck_1724_;
goto v_resetjp_1688_;
}
else
{
lean_inc(v_scoreThreshold_1687_);
lean_dec(v_x_665_);
v___x_1689_ = lean_box(0);
v_isShared_1690_ = v_isSharedCheck_1724_;
goto v_resetjp_1688_;
}
v_resetjp_1688_:
{
lean_object* v___y_1692_; lean_object* v___x_1720_; uint8_t v___x_1721_; 
v___x_1720_ = lean_unsigned_to_nat(1024u);
v___x_1721_ = lean_nat_dec_le(v___x_1720_, v_prec_666_);
if (v___x_1721_ == 0)
{
lean_object* v___x_1722_; 
v___x_1722_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1692_ = v___x_1722_;
goto v___jp_1691_;
}
else
{
lean_object* v___x_1723_; 
v___x_1723_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1692_ = v___x_1723_;
goto v___jp_1691_;
}
v___jp_1691_:
{
lean_object* v_num_1693_; lean_object* v_den_1694_; lean_object* v___x_1695_; lean_object* v___x_1696_; uint8_t v___x_1697_; 
v_num_1693_ = lean_ctor_get(v_scoreThreshold_1687_, 0);
lean_inc(v_num_1693_);
v_den_1694_ = lean_ctor_get(v_scoreThreshold_1687_, 1);
lean_inc(v_den_1694_);
lean_dec_ref(v_scoreThreshold_1687_);
v___x_1695_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__89));
v___x_1696_ = lean_unsigned_to_nat(1u);
v___x_1697_ = lean_nat_dec_eq(v_den_1694_, v___x_1696_);
if (v___x_1697_ == 0)
{
lean_object* v___x_1698_; lean_object* v___x_1699_; lean_object* v___x_1700_; lean_object* v___x_1701_; lean_object* v___x_1702_; lean_object* v___x_1703_; lean_object* v___x_1704_; lean_object* v___x_1706_; 
v___x_1698_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1699_ = l_Int_repr(v_num_1693_);
lean_dec(v_num_1693_);
v___x_1700_ = lean_string_append(v___x_1698_, v___x_1699_);
lean_dec_ref(v___x_1699_);
v___x_1701_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1702_ = lean_string_append(v___x_1700_, v___x_1701_);
v___x_1703_ = l_Nat_reprFast(v_den_1694_);
v___x_1704_ = lean_string_append(v___x_1702_, v___x_1703_);
lean_dec_ref(v___x_1703_);
if (v_isShared_1690_ == 0)
{
lean_ctor_set_tag(v___x_1689_, 3);
lean_ctor_set(v___x_1689_, 0, v___x_1704_);
v___x_1706_ = v___x_1689_;
goto v_reusejp_1705_;
}
else
{
lean_object* v_reuseFailAlloc_1707_; 
v_reuseFailAlloc_1707_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1707_, 0, v___x_1704_);
v___x_1706_ = v_reuseFailAlloc_1707_;
goto v_reusejp_1705_;
}
v_reusejp_1705_:
{
v___y_889_ = v___y_1692_;
v___y_890_ = v___x_1695_;
v___y_891_ = v___x_1706_;
goto v___jp_888_;
}
}
else
{
lean_object* v___x_1708_; lean_object* v___x_1709_; uint8_t v___x_1710_; 
lean_dec(v_den_1694_);
v___x_1708_ = lean_unsigned_to_nat(0u);
v___x_1709_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1710_ = lean_int_dec_lt(v_num_1693_, v___x_1709_);
if (v___x_1710_ == 0)
{
lean_object* v___x_1711_; lean_object* v___x_1713_; 
v___x_1711_ = l_Int_repr(v_num_1693_);
lean_dec(v_num_1693_);
if (v_isShared_1690_ == 0)
{
lean_ctor_set_tag(v___x_1689_, 3);
lean_ctor_set(v___x_1689_, 0, v___x_1711_);
v___x_1713_ = v___x_1689_;
goto v_reusejp_1712_;
}
else
{
lean_object* v_reuseFailAlloc_1714_; 
v_reuseFailAlloc_1714_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1714_, 0, v___x_1711_);
v___x_1713_ = v_reuseFailAlloc_1714_;
goto v_reusejp_1712_;
}
v_reusejp_1712_:
{
v___y_889_ = v___y_1692_;
v___y_890_ = v___x_1695_;
v___y_891_ = v___x_1713_;
goto v___jp_888_;
}
}
else
{
lean_object* v___x_1715_; lean_object* v___x_1717_; 
v___x_1715_ = l_Int_repr(v_num_1693_);
lean_dec(v_num_1693_);
if (v_isShared_1690_ == 0)
{
lean_ctor_set_tag(v___x_1689_, 3);
lean_ctor_set(v___x_1689_, 0, v___x_1715_);
v___x_1717_ = v___x_1689_;
goto v_reusejp_1716_;
}
else
{
lean_object* v_reuseFailAlloc_1719_; 
v_reuseFailAlloc_1719_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1719_, 0, v___x_1715_);
v___x_1717_ = v_reuseFailAlloc_1719_;
goto v_reusejp_1716_;
}
v_reusejp_1716_:
{
lean_object* v___x_1718_; 
v___x_1718_ = l_Repr_addAppParen(v___x_1717_, v___x_1708_);
v___y_889_ = v___y_1692_;
v___y_890_ = v___x_1695_;
v___y_891_ = v___x_1718_;
goto v___jp_888_;
}
}
}
}
}
}
default: 
{
lean_object* v_maxCoverage_1725_; lean_object* v___x_1727_; uint8_t v_isShared_1728_; uint8_t v_isSharedCheck_1762_; 
v_maxCoverage_1725_ = lean_ctor_get(v_x_665_, 0);
v_isSharedCheck_1762_ = !lean_is_exclusive(v_x_665_);
if (v_isSharedCheck_1762_ == 0)
{
v___x_1727_ = v_x_665_;
v_isShared_1728_ = v_isSharedCheck_1762_;
goto v_resetjp_1726_;
}
else
{
lean_inc(v_maxCoverage_1725_);
lean_dec(v_x_665_);
v___x_1727_ = lean_box(0);
v_isShared_1728_ = v_isSharedCheck_1762_;
goto v_resetjp_1726_;
}
v_resetjp_1726_:
{
lean_object* v___y_1730_; lean_object* v___x_1758_; uint8_t v___x_1759_; 
v___x_1758_ = lean_unsigned_to_nat(1024u);
v___x_1759_ = lean_nat_dec_le(v___x_1758_, v_prec_666_);
if (v___x_1759_ == 0)
{
lean_object* v___x_1760_; 
v___x_1760_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25);
v___y_1730_ = v___x_1760_;
goto v___jp_1729_;
}
else
{
lean_object* v___x_1761_; 
v___x_1761_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__26);
v___y_1730_ = v___x_1761_;
goto v___jp_1729_;
}
v___jp_1729_:
{
lean_object* v_num_1731_; lean_object* v_den_1732_; lean_object* v___x_1733_; lean_object* v___x_1734_; uint8_t v___x_1735_; 
v_num_1731_ = lean_ctor_get(v_maxCoverage_1725_, 0);
lean_inc(v_num_1731_);
v_den_1732_ = lean_ctor_get(v_maxCoverage_1725_, 1);
lean_inc(v_den_1732_);
lean_dec_ref(v_maxCoverage_1725_);
v___x_1733_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__92));
v___x_1734_ = lean_unsigned_to_nat(1u);
v___x_1735_ = lean_nat_dec_eq(v_den_1732_, v___x_1734_);
if (v___x_1735_ == 0)
{
lean_object* v___x_1736_; lean_object* v___x_1737_; lean_object* v___x_1738_; lean_object* v___x_1739_; lean_object* v___x_1740_; lean_object* v___x_1741_; lean_object* v___x_1742_; lean_object* v___x_1744_; 
v___x_1736_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__30));
v___x_1737_ = l_Int_repr(v_num_1731_);
lean_dec(v_num_1731_);
v___x_1738_ = lean_string_append(v___x_1736_, v___x_1737_);
lean_dec_ref(v___x_1737_);
v___x_1739_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__31));
v___x_1740_ = lean_string_append(v___x_1738_, v___x_1739_);
v___x_1741_ = l_Nat_reprFast(v_den_1732_);
v___x_1742_ = lean_string_append(v___x_1740_, v___x_1741_);
lean_dec_ref(v___x_1741_);
if (v_isShared_1728_ == 0)
{
lean_ctor_set_tag(v___x_1727_, 3);
lean_ctor_set(v___x_1727_, 0, v___x_1742_);
v___x_1744_ = v___x_1727_;
goto v_reusejp_1743_;
}
else
{
lean_object* v_reuseFailAlloc_1745_; 
v_reuseFailAlloc_1745_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1745_, 0, v___x_1742_);
v___x_1744_ = v_reuseFailAlloc_1745_;
goto v_reusejp_1743_;
}
v_reusejp_1743_:
{
v___y_668_ = v___x_1733_;
v___y_669_ = v___y_1730_;
v___y_670_ = v___x_1744_;
goto v___jp_667_;
}
}
else
{
lean_object* v___x_1746_; lean_object* v___x_1747_; uint8_t v___x_1748_; 
lean_dec(v_den_1732_);
v___x_1746_ = lean_unsigned_to_nat(0u);
v___x_1747_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__32);
v___x_1748_ = lean_int_dec_lt(v_num_1731_, v___x_1747_);
if (v___x_1748_ == 0)
{
lean_object* v___x_1749_; lean_object* v___x_1751_; 
v___x_1749_ = l_Int_repr(v_num_1731_);
lean_dec(v_num_1731_);
if (v_isShared_1728_ == 0)
{
lean_ctor_set_tag(v___x_1727_, 3);
lean_ctor_set(v___x_1727_, 0, v___x_1749_);
v___x_1751_ = v___x_1727_;
goto v_reusejp_1750_;
}
else
{
lean_object* v_reuseFailAlloc_1752_; 
v_reuseFailAlloc_1752_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1752_, 0, v___x_1749_);
v___x_1751_ = v_reuseFailAlloc_1752_;
goto v_reusejp_1750_;
}
v_reusejp_1750_:
{
v___y_668_ = v___x_1733_;
v___y_669_ = v___y_1730_;
v___y_670_ = v___x_1751_;
goto v___jp_667_;
}
}
else
{
lean_object* v___x_1753_; lean_object* v___x_1755_; 
v___x_1753_ = l_Int_repr(v_num_1731_);
lean_dec(v_num_1731_);
if (v_isShared_1728_ == 0)
{
lean_ctor_set_tag(v___x_1727_, 3);
lean_ctor_set(v___x_1727_, 0, v___x_1753_);
v___x_1755_ = v___x_1727_;
goto v_reusejp_1754_;
}
else
{
lean_object* v_reuseFailAlloc_1757_; 
v_reuseFailAlloc_1757_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_1757_, 0, v___x_1753_);
v___x_1755_ = v_reuseFailAlloc_1757_;
goto v_reusejp_1754_;
}
v_reusejp_1754_:
{
lean_object* v___x_1756_; 
v___x_1756_ = l_Repr_addAppParen(v___x_1755_, v___x_1746_);
v___y_668_ = v___x_1733_;
v___y_669_ = v___y_1730_;
v___y_670_ = v___x_1756_;
goto v___jp_667_;
}
}
}
}
}
}
}
v___jp_667_:
{
lean_object* v___x_671_; lean_object* v___x_672_; uint8_t v___x_673_; lean_object* v___x_674_; lean_object* v___x_675_; 
lean_inc(v___y_668_);
v___x_671_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_671_, 0, v___y_668_);
lean_ctor_set(v___x_671_, 1, v___y_670_);
lean_inc(v___y_669_);
v___x_672_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_672_, 0, v___y_669_);
lean_ctor_set(v___x_672_, 1, v___x_671_);
v___x_673_ = 0;
v___x_674_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_674_, 0, v___x_672_);
lean_ctor_set_uint8(v___x_674_, sizeof(void*)*1, v___x_673_);
v___x_675_ = l_Repr_addAppParen(v___x_674_, v_prec_666_);
return v___x_675_;
}
v___jp_676_:
{
lean_object* v___x_680_; lean_object* v___x_681_; uint8_t v___x_682_; lean_object* v___x_683_; lean_object* v___x_684_; 
lean_inc(v___y_678_);
v___x_680_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_680_, 0, v___y_678_);
lean_ctor_set(v___x_680_, 1, v___y_679_);
lean_inc(v___y_677_);
v___x_681_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_681_, 0, v___y_677_);
lean_ctor_set(v___x_681_, 1, v___x_680_);
v___x_682_ = 0;
v___x_683_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_683_, 0, v___x_681_);
lean_ctor_set_uint8(v___x_683_, sizeof(void*)*1, v___x_682_);
v___x_684_ = l_Repr_addAppParen(v___x_683_, v_prec_666_);
return v___x_684_;
}
v___jp_685_:
{
lean_object* v___x_689_; lean_object* v___x_690_; uint8_t v___x_691_; lean_object* v___x_692_; lean_object* v___x_693_; 
v___x_689_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_689_, 0, v___y_686_);
lean_ctor_set(v___x_689_, 1, v___y_688_);
lean_inc(v___y_687_);
v___x_690_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_690_, 0, v___y_687_);
lean_ctor_set(v___x_690_, 1, v___x_689_);
v___x_691_ = 0;
v___x_692_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_692_, 0, v___x_690_);
lean_ctor_set_uint8(v___x_692_, sizeof(void*)*1, v___x_691_);
v___x_693_ = l_Repr_addAppParen(v___x_692_, v_prec_666_);
return v___x_693_;
}
v___jp_694_:
{
lean_object* v___x_696_; lean_object* v___x_697_; uint8_t v___x_698_; lean_object* v___x_699_; lean_object* v___x_700_; 
v___x_696_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__1));
lean_inc(v___y_695_);
v___x_697_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_697_, 0, v___y_695_);
lean_ctor_set(v___x_697_, 1, v___x_696_);
v___x_698_ = 0;
v___x_699_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_699_, 0, v___x_697_);
lean_ctor_set_uint8(v___x_699_, sizeof(void*)*1, v___x_698_);
v___x_700_ = l_Repr_addAppParen(v___x_699_, v_prec_666_);
return v___x_700_;
}
v___jp_701_:
{
lean_object* v___x_705_; lean_object* v___x_706_; uint8_t v___x_707_; lean_object* v___x_708_; lean_object* v___x_709_; 
lean_inc(v___y_703_);
v___x_705_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_705_, 0, v___y_703_);
lean_ctor_set(v___x_705_, 1, v___y_704_);
lean_inc(v___y_702_);
v___x_706_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_706_, 0, v___y_702_);
lean_ctor_set(v___x_706_, 1, v___x_705_);
v___x_707_ = 0;
v___x_708_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_708_, 0, v___x_706_);
lean_ctor_set_uint8(v___x_708_, sizeof(void*)*1, v___x_707_);
v___x_709_ = l_Repr_addAppParen(v___x_708_, v_prec_666_);
return v___x_709_;
}
v___jp_710_:
{
lean_object* v___x_712_; lean_object* v___x_713_; uint8_t v___x_714_; lean_object* v___x_715_; lean_object* v___x_716_; 
v___x_712_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__3));
lean_inc(v___y_711_);
v___x_713_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_713_, 0, v___y_711_);
lean_ctor_set(v___x_713_, 1, v___x_712_);
v___x_714_ = 0;
v___x_715_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_715_, 0, v___x_713_);
lean_ctor_set_uint8(v___x_715_, sizeof(void*)*1, v___x_714_);
v___x_716_ = l_Repr_addAppParen(v___x_715_, v_prec_666_);
return v___x_716_;
}
v___jp_717_:
{
lean_object* v___x_721_; lean_object* v___x_722_; uint8_t v___x_723_; lean_object* v___x_724_; lean_object* v___x_725_; 
lean_inc(v___y_719_);
v___x_721_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_721_, 0, v___y_719_);
lean_ctor_set(v___x_721_, 1, v___y_720_);
lean_inc(v___y_718_);
v___x_722_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_722_, 0, v___y_718_);
lean_ctor_set(v___x_722_, 1, v___x_721_);
v___x_723_ = 0;
v___x_724_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_724_, 0, v___x_722_);
lean_ctor_set_uint8(v___x_724_, sizeof(void*)*1, v___x_723_);
v___x_725_ = l_Repr_addAppParen(v___x_724_, v_prec_666_);
return v___x_725_;
}
v___jp_726_:
{
lean_object* v___x_728_; lean_object* v___x_729_; uint8_t v___x_730_; lean_object* v___x_731_; lean_object* v___x_732_; 
v___x_728_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__5));
lean_inc(v___y_727_);
v___x_729_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_729_, 0, v___y_727_);
lean_ctor_set(v___x_729_, 1, v___x_728_);
v___x_730_ = 0;
v___x_731_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_731_, 0, v___x_729_);
lean_ctor_set_uint8(v___x_731_, sizeof(void*)*1, v___x_730_);
v___x_732_ = l_Repr_addAppParen(v___x_731_, v_prec_666_);
return v___x_732_;
}
v___jp_733_:
{
lean_object* v___x_735_; lean_object* v___x_736_; uint8_t v___x_737_; lean_object* v___x_738_; lean_object* v___x_739_; 
v___x_735_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__7));
lean_inc(v___y_734_);
v___x_736_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_736_, 0, v___y_734_);
lean_ctor_set(v___x_736_, 1, v___x_735_);
v___x_737_ = 0;
v___x_738_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_738_, 0, v___x_736_);
lean_ctor_set_uint8(v___x_738_, sizeof(void*)*1, v___x_737_);
v___x_739_ = l_Repr_addAppParen(v___x_738_, v_prec_666_);
return v___x_739_;
}
v___jp_740_:
{
lean_object* v___x_742_; lean_object* v___x_743_; uint8_t v___x_744_; lean_object* v___x_745_; lean_object* v___x_746_; 
v___x_742_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9));
lean_inc(v___y_741_);
v___x_743_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_743_, 0, v___y_741_);
lean_ctor_set(v___x_743_, 1, v___x_742_);
v___x_744_ = 0;
v___x_745_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_745_, 0, v___x_743_);
lean_ctor_set_uint8(v___x_745_, sizeof(void*)*1, v___x_744_);
v___x_746_ = l_Repr_addAppParen(v___x_745_, v_prec_666_);
return v___x_746_;
}
v___jp_747_:
{
lean_object* v___x_751_; lean_object* v___x_752_; uint8_t v___x_753_; lean_object* v___x_754_; lean_object* v___x_755_; 
v___x_751_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_751_, 0, v___y_748_);
lean_ctor_set(v___x_751_, 1, v___y_750_);
lean_inc(v___y_749_);
v___x_752_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_752_, 0, v___y_749_);
lean_ctor_set(v___x_752_, 1, v___x_751_);
v___x_753_ = 0;
v___x_754_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_754_, 0, v___x_752_);
lean_ctor_set_uint8(v___x_754_, sizeof(void*)*1, v___x_753_);
v___x_755_ = l_Repr_addAppParen(v___x_754_, v_prec_666_);
return v___x_755_;
}
v___jp_756_:
{
lean_object* v___x_760_; lean_object* v___x_761_; uint8_t v___x_762_; lean_object* v___x_763_; lean_object* v___x_764_; 
v___x_760_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_760_, 0, v___y_757_);
lean_ctor_set(v___x_760_, 1, v___y_759_);
lean_inc(v___y_758_);
v___x_761_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_761_, 0, v___y_758_);
lean_ctor_set(v___x_761_, 1, v___x_760_);
v___x_762_ = 0;
v___x_763_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_763_, 0, v___x_761_);
lean_ctor_set_uint8(v___x_763_, sizeof(void*)*1, v___x_762_);
v___x_764_ = l_Repr_addAppParen(v___x_763_, v_prec_666_);
return v___x_764_;
}
v___jp_765_:
{
lean_object* v___x_769_; lean_object* v___x_770_; uint8_t v___x_771_; lean_object* v___x_772_; lean_object* v___x_773_; 
v___x_769_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_769_, 0, v___y_767_);
lean_ctor_set(v___x_769_, 1, v___y_768_);
lean_inc(v___y_766_);
v___x_770_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_770_, 0, v___y_766_);
lean_ctor_set(v___x_770_, 1, v___x_769_);
v___x_771_ = 0;
v___x_772_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_772_, 0, v___x_770_);
lean_ctor_set_uint8(v___x_772_, sizeof(void*)*1, v___x_771_);
v___x_773_ = l_Repr_addAppParen(v___x_772_, v_prec_666_);
return v___x_773_;
}
v___jp_774_:
{
lean_object* v___x_776_; lean_object* v___x_777_; uint8_t v___x_778_; lean_object* v___x_779_; lean_object* v___x_780_; 
v___x_776_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__11));
lean_inc(v___y_775_);
v___x_777_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_777_, 0, v___y_775_);
lean_ctor_set(v___x_777_, 1, v___x_776_);
v___x_778_ = 0;
v___x_779_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_779_, 0, v___x_777_);
lean_ctor_set_uint8(v___x_779_, sizeof(void*)*1, v___x_778_);
v___x_780_ = l_Repr_addAppParen(v___x_779_, v_prec_666_);
return v___x_780_;
}
v___jp_781_:
{
lean_object* v___x_785_; lean_object* v___x_786_; uint8_t v___x_787_; lean_object* v___x_788_; lean_object* v___x_789_; 
v___x_785_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_785_, 0, v___y_783_);
lean_ctor_set(v___x_785_, 1, v___y_784_);
lean_inc(v___y_782_);
v___x_786_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_786_, 0, v___y_782_);
lean_ctor_set(v___x_786_, 1, v___x_785_);
v___x_787_ = 0;
v___x_788_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_788_, 0, v___x_786_);
lean_ctor_set_uint8(v___x_788_, sizeof(void*)*1, v___x_787_);
v___x_789_ = l_Repr_addAppParen(v___x_788_, v_prec_666_);
return v___x_789_;
}
v___jp_790_:
{
lean_object* v___x_792_; lean_object* v___x_793_; uint8_t v___x_794_; lean_object* v___x_795_; lean_object* v___x_796_; 
v___x_792_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__13));
lean_inc(v___y_791_);
v___x_793_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_793_, 0, v___y_791_);
lean_ctor_set(v___x_793_, 1, v___x_792_);
v___x_794_ = 0;
v___x_795_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_795_, 0, v___x_793_);
lean_ctor_set_uint8(v___x_795_, sizeof(void*)*1, v___x_794_);
v___x_796_ = l_Repr_addAppParen(v___x_795_, v_prec_666_);
return v___x_796_;
}
v___jp_797_:
{
lean_object* v___x_799_; lean_object* v___x_800_; uint8_t v___x_801_; lean_object* v___x_802_; lean_object* v___x_803_; 
v___x_799_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15));
lean_inc(v___y_798_);
v___x_800_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_800_, 0, v___y_798_);
lean_ctor_set(v___x_800_, 1, v___x_799_);
v___x_801_ = 0;
v___x_802_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_802_, 0, v___x_800_);
lean_ctor_set_uint8(v___x_802_, sizeof(void*)*1, v___x_801_);
v___x_803_ = l_Repr_addAppParen(v___x_802_, v_prec_666_);
return v___x_803_;
}
v___jp_804_:
{
lean_object* v___x_806_; lean_object* v___x_807_; uint8_t v___x_808_; lean_object* v___x_809_; lean_object* v___x_810_; 
v___x_806_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__17));
lean_inc(v___y_805_);
v___x_807_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_807_, 0, v___y_805_);
lean_ctor_set(v___x_807_, 1, v___x_806_);
v___x_808_ = 0;
v___x_809_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_809_, 0, v___x_807_);
lean_ctor_set_uint8(v___x_809_, sizeof(void*)*1, v___x_808_);
v___x_810_ = l_Repr_addAppParen(v___x_809_, v_prec_666_);
return v___x_810_;
}
v___jp_811_:
{
lean_object* v___x_813_; lean_object* v___x_814_; uint8_t v___x_815_; lean_object* v___x_816_; lean_object* v___x_817_; 
v___x_813_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__19));
lean_inc(v___y_812_);
v___x_814_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_814_, 0, v___y_812_);
lean_ctor_set(v___x_814_, 1, v___x_813_);
v___x_815_ = 0;
v___x_816_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_816_, 0, v___x_814_);
lean_ctor_set_uint8(v___x_816_, sizeof(void*)*1, v___x_815_);
v___x_817_ = l_Repr_addAppParen(v___x_816_, v_prec_666_);
return v___x_817_;
}
v___jp_818_:
{
lean_object* v___x_820_; lean_object* v___x_821_; uint8_t v___x_822_; lean_object* v___x_823_; lean_object* v___x_824_; 
v___x_820_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__21));
lean_inc(v___y_819_);
v___x_821_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_821_, 0, v___y_819_);
lean_ctor_set(v___x_821_, 1, v___x_820_);
v___x_822_ = 0;
v___x_823_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_823_, 0, v___x_821_);
lean_ctor_set_uint8(v___x_823_, sizeof(void*)*1, v___x_822_);
v___x_824_ = l_Repr_addAppParen(v___x_823_, v_prec_666_);
return v___x_824_;
}
v___jp_825_:
{
lean_object* v___x_829_; lean_object* v___x_830_; uint8_t v___x_831_; lean_object* v___x_832_; lean_object* v___x_833_; 
v___x_829_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_829_, 0, v___y_826_);
lean_ctor_set(v___x_829_, 1, v___y_828_);
lean_inc(v___y_827_);
v___x_830_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_830_, 0, v___y_827_);
lean_ctor_set(v___x_830_, 1, v___x_829_);
v___x_831_ = 0;
v___x_832_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_832_, 0, v___x_830_);
lean_ctor_set_uint8(v___x_832_, sizeof(void*)*1, v___x_831_);
v___x_833_ = l_Repr_addAppParen(v___x_832_, v_prec_666_);
return v___x_833_;
}
v___jp_834_:
{
lean_object* v___x_838_; lean_object* v___x_839_; uint8_t v___x_840_; lean_object* v___x_841_; lean_object* v___x_842_; 
lean_inc(v___y_836_);
v___x_838_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_838_, 0, v___y_836_);
lean_ctor_set(v___x_838_, 1, v___y_837_);
lean_inc(v___y_835_);
v___x_839_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_839_, 0, v___y_835_);
lean_ctor_set(v___x_839_, 1, v___x_838_);
v___x_840_ = 0;
v___x_841_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_841_, 0, v___x_839_);
lean_ctor_set_uint8(v___x_841_, sizeof(void*)*1, v___x_840_);
v___x_842_ = l_Repr_addAppParen(v___x_841_, v_prec_666_);
return v___x_842_;
}
v___jp_843_:
{
lean_object* v___x_847_; lean_object* v___x_848_; uint8_t v___x_849_; lean_object* v___x_850_; lean_object* v___x_851_; 
lean_inc(v___y_844_);
v___x_847_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_847_, 0, v___y_844_);
lean_ctor_set(v___x_847_, 1, v___y_846_);
lean_inc(v___y_845_);
v___x_848_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_848_, 0, v___y_845_);
lean_ctor_set(v___x_848_, 1, v___x_847_);
v___x_849_ = 0;
v___x_850_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_850_, 0, v___x_848_);
lean_ctor_set_uint8(v___x_850_, sizeof(void*)*1, v___x_849_);
v___x_851_ = l_Repr_addAppParen(v___x_850_, v_prec_666_);
return v___x_851_;
}
v___jp_852_:
{
lean_object* v___x_856_; lean_object* v___x_857_; uint8_t v___x_858_; lean_object* v___x_859_; lean_object* v___x_860_; 
lean_inc(v___y_854_);
v___x_856_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_856_, 0, v___y_854_);
lean_ctor_set(v___x_856_, 1, v___y_855_);
lean_inc(v___y_853_);
v___x_857_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_857_, 0, v___y_853_);
lean_ctor_set(v___x_857_, 1, v___x_856_);
v___x_858_ = 0;
v___x_859_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_859_, 0, v___x_857_);
lean_ctor_set_uint8(v___x_859_, sizeof(void*)*1, v___x_858_);
v___x_860_ = l_Repr_addAppParen(v___x_859_, v_prec_666_);
return v___x_860_;
}
v___jp_861_:
{
lean_object* v___x_865_; lean_object* v___x_866_; uint8_t v___x_867_; lean_object* v___x_868_; lean_object* v___x_869_; 
lean_inc(v___y_862_);
v___x_865_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_865_, 0, v___y_862_);
lean_ctor_set(v___x_865_, 1, v___y_864_);
lean_inc(v___y_863_);
v___x_866_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_866_, 0, v___y_863_);
lean_ctor_set(v___x_866_, 1, v___x_865_);
v___x_867_ = 0;
v___x_868_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_868_, 0, v___x_866_);
lean_ctor_set_uint8(v___x_868_, sizeof(void*)*1, v___x_867_);
v___x_869_ = l_Repr_addAppParen(v___x_868_, v_prec_666_);
return v___x_869_;
}
v___jp_870_:
{
lean_object* v___x_874_; lean_object* v___x_875_; uint8_t v___x_876_; lean_object* v___x_877_; lean_object* v___x_878_; 
lean_inc(v___y_871_);
v___x_874_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_874_, 0, v___y_871_);
lean_ctor_set(v___x_874_, 1, v___y_873_);
lean_inc(v___y_872_);
v___x_875_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_875_, 0, v___y_872_);
lean_ctor_set(v___x_875_, 1, v___x_874_);
v___x_876_ = 0;
v___x_877_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_877_, 0, v___x_875_);
lean_ctor_set_uint8(v___x_877_, sizeof(void*)*1, v___x_876_);
v___x_878_ = l_Repr_addAppParen(v___x_877_, v_prec_666_);
return v___x_878_;
}
v___jp_879_:
{
lean_object* v___x_883_; lean_object* v___x_884_; uint8_t v___x_885_; lean_object* v___x_886_; lean_object* v___x_887_; 
lean_inc(v___y_881_);
v___x_883_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_883_, 0, v___y_881_);
lean_ctor_set(v___x_883_, 1, v___y_882_);
lean_inc(v___y_880_);
v___x_884_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_884_, 0, v___y_880_);
lean_ctor_set(v___x_884_, 1, v___x_883_);
v___x_885_ = 0;
v___x_886_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_886_, 0, v___x_884_);
lean_ctor_set_uint8(v___x_886_, sizeof(void*)*1, v___x_885_);
v___x_887_ = l_Repr_addAppParen(v___x_886_, v_prec_666_);
return v___x_887_;
}
v___jp_888_:
{
lean_object* v___x_892_; lean_object* v___x_893_; uint8_t v___x_894_; lean_object* v___x_895_; lean_object* v___x_896_; 
lean_inc(v___y_890_);
v___x_892_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_892_, 0, v___y_890_);
lean_ctor_set(v___x_892_, 1, v___y_891_);
lean_inc(v___y_889_);
v___x_893_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_893_, 0, v___y_889_);
lean_ctor_set(v___x_893_, 1, v___x_892_);
v___x_894_ = 0;
v___x_895_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_895_, 0, v___x_893_);
lean_ctor_set_uint8(v___x_895_, sizeof(void*)*1, v___x_894_);
v___x_896_ = l_Repr_addAppParen(v___x_895_, v_prec_666_);
return v___x_896_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___boxed(lean_object* v_x_1763_, lean_object* v_prec_1764_){
_start:
{
lean_object* v_res_1765_; 
v_res_1765_ = lp_BioCompiler_BioCompiler_instReprTypePredicate_repr(v_x_1763_, v_prec_1764_);
lean_dec(v_prec_1764_);
return v_res_1765_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0(lean_object* v_a_1766_, lean_object* v_n_1767_){
_start:
{
lean_object* v___x_1768_; 
v___x_1768_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg(v_a_1766_);
return v___x_1768_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___boxed(lean_object* v_a_1769_, lean_object* v_n_1770_){
_start:
{
lean_object* v_res_1771_; 
v_res_1771_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0(v_a_1769_, v_n_1770_);
lean_dec(v_n_1770_);
return v_res_1771_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1(lean_object* v_a_1772_, lean_object* v_n_1773_){
_start:
{
lean_object* v___x_1774_; 
v___x_1774_ = lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___redArg(v_a_1772_);
return v___x_1774_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___boxed(lean_object* v_a_1775_, lean_object* v_n_1776_){
_start:
{
lean_object* v_res_1777_; 
v_res_1777_ = lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1(v_a_1775_, v_n_1776_);
lean_dec(v_n_1776_);
return v_res_1777_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_isSLOT(lean_object* v_x_1780_){
_start:
{
switch(lean_obj_tag(v_x_1780_))
{
case 0:
{
uint8_t v___x_1781_; 
v___x_1781_ = 0;
return v___x_1781_;
}
case 1:
{
uint8_t v___x_1782_; 
v___x_1782_ = 0;
return v___x_1782_;
}
case 2:
{
uint8_t v___x_1783_; 
v___x_1783_ = 0;
return v___x_1783_;
}
case 3:
{
uint8_t v___x_1784_; 
v___x_1784_ = 0;
return v___x_1784_;
}
case 4:
{
uint8_t v___x_1785_; 
v___x_1785_ = 0;
return v___x_1785_;
}
case 5:
{
uint8_t v___x_1786_; 
v___x_1786_ = 0;
return v___x_1786_;
}
case 6:
{
uint8_t v___x_1787_; 
v___x_1787_ = 0;
return v___x_1787_;
}
case 7:
{
uint8_t v___x_1788_; 
v___x_1788_ = 0;
return v___x_1788_;
}
case 8:
{
uint8_t v___x_1789_; 
v___x_1789_ = 0;
return v___x_1789_;
}
case 9:
{
uint8_t v___x_1790_; 
v___x_1790_ = 0;
return v___x_1790_;
}
case 10:
{
uint8_t v___x_1791_; 
v___x_1791_ = 0;
return v___x_1791_;
}
case 11:
{
uint8_t v___x_1792_; 
v___x_1792_ = 0;
return v___x_1792_;
}
case 12:
{
uint8_t v___x_1793_; 
v___x_1793_ = 0;
return v___x_1793_;
}
case 18:
{
uint8_t v___x_1794_; 
v___x_1794_ = 1;
return v___x_1794_;
}
case 19:
{
uint8_t v___x_1795_; 
v___x_1795_ = 1;
return v___x_1795_;
}
case 20:
{
uint8_t v___x_1796_; 
v___x_1796_ = 1;
return v___x_1796_;
}
case 23:
{
uint8_t v___x_1797_; 
v___x_1797_ = 1;
return v___x_1797_;
}
case 26:
{
uint8_t v___x_1798_; 
v___x_1798_ = 1;
return v___x_1798_;
}
default: 
{
uint8_t v___x_1799_; 
v___x_1799_ = 1;
return v___x_1799_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_isSLOT___boxed(lean_object* v_x_1800_){
_start:
{
uint8_t v_res_1801_; lean_object* v_r_1802_; 
v_res_1801_ = lp_BioCompiler_BioCompiler_isSLOT(v_x_1800_);
lean_dec(v_x_1800_);
v_r_1802_ = lean_box(v_res_1801_);
return v_r_1802_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluate___redArg(lean_object* v_inst_1803_, lean_object* v_inst_1804_, lean_object* v_inst_1805_, lean_object* v_inst_1806_, lean_object* v_inst_1807_, lean_object* v_x_1808_, lean_object* v_x_1809_, lean_object* v_x_1810_){
_start:
{
lean_object* v_org_1812_; lean_object* v_threshold_1813_; lean_object* v_seq_1814_; 
switch(lean_obj_tag(v_x_1808_))
{
case 0:
{
lean_object* v_cellType_1820_; lean_object* v_cellType_1821_; uint8_t v___x_1822_; 
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v_cellType_1820_ = lean_ctor_get(v_x_1808_, 0);
lean_inc_ref(v_cellType_1820_);
lean_dec_ref(v_x_1808_);
v_cellType_1821_ = lean_ctor_get(v_x_1810_, 0);
v___x_1822_ = lean_string_dec_eq(v_cellType_1821_, v_cellType_1820_);
lean_dec_ref(v_cellType_1820_);
if (v___x_1822_ == 0)
{
uint8_t v___x_1823_; 
lean_dec(v_x_1809_);
lean_dec_ref(v_inst_1807_);
v___x_1823_ = 2;
return v___x_1823_;
}
else
{
lean_object* v___x_1824_; 
v___x_1824_ = lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg(v_inst_1807_, v_x_1809_);
lean_dec(v_x_1809_);
if (lean_obj_tag(v___x_1824_) == 1)
{
lean_object* v_tail_1825_; 
v_tail_1825_ = lean_ctor_get(v___x_1824_, 1);
lean_inc(v_tail_1825_);
lean_dec_ref(v___x_1824_);
if (lean_obj_tag(v_tail_1825_) == 0)
{
uint8_t v___x_1826_; 
v___x_1826_ = 0;
return v___x_1826_;
}
else
{
uint8_t v___x_1827_; 
lean_dec(v_tail_1825_);
v___x_1827_ = 1;
return v___x_1827_;
}
}
else
{
uint8_t v___x_1828_; 
lean_dec(v___x_1824_);
v___x_1828_ = 1;
return v___x_1828_;
}
}
}
case 1:
{
lean_object* v_hasCrypticSpliceSite_1829_; lean_object* v_hasBorderlineSpliceSite_1830_; lean_object* v___x_1831_; uint8_t v___x_1832_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
v_hasCrypticSpliceSite_1829_ = lean_ctor_get(v_inst_1803_, 0);
lean_inc_ref(v_hasCrypticSpliceSite_1829_);
v_hasBorderlineSpliceSite_1830_ = lean_ctor_get(v_inst_1803_, 1);
lean_inc_ref(v_hasBorderlineSpliceSite_1830_);
lean_dec_ref(v_inst_1803_);
lean_inc(v_x_1809_);
v___x_1831_ = lean_apply_1(v_hasCrypticSpliceSite_1829_, v_x_1809_);
v___x_1832_ = lean_unbox(v___x_1831_);
if (v___x_1832_ == 0)
{
lean_object* v___x_1833_; uint8_t v___x_1834_; 
v___x_1833_ = lean_apply_1(v_hasBorderlineSpliceSite_1830_, v_x_1809_);
v___x_1834_ = lean_unbox(v___x_1833_);
if (v___x_1834_ == 0)
{
uint8_t v___x_1835_; 
v___x_1835_ = 0;
return v___x_1835_;
}
else
{
uint8_t v___x_1836_; 
v___x_1836_ = 2;
return v___x_1836_;
}
}
else
{
uint8_t v___x_1837_; 
lean_dec_ref(v_hasBorderlineSpliceSite_1830_);
lean_dec(v_x_1809_);
v___x_1837_ = 1;
return v___x_1837_;
}
}
case 2:
{
lean_object* v_organism_1838_; lean_object* v_threshold_1839_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1803_);
v_organism_1838_ = lean_ctor_get(v_x_1808_, 0);
lean_inc_ref(v_organism_1838_);
v_threshold_1839_ = lean_ctor_get(v_x_1808_, 1);
lean_inc_ref(v_threshold_1839_);
lean_dec_ref(v_x_1808_);
v_org_1812_ = v_organism_1838_;
v_threshold_1813_ = v_threshold_1839_;
v_seq_1814_ = v_x_1809_;
goto v___jp_1811_;
}
case 3:
{
lean_object* v_lo_1840_; lean_object* v_hi_1841_; lean_object* v___x_1842_; uint8_t v___x_1843_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v_lo_1840_ = lean_ctor_get(v_x_1808_, 0);
lean_inc_ref(v_lo_1840_);
v_hi_1841_ = lean_ctor_get(v_x_1808_, 1);
lean_inc_ref(v_hi_1841_);
lean_dec_ref(v_x_1808_);
v___x_1842_ = lp_BioCompiler_BioCompiler_Sequence_gcContent(v_x_1809_);
lean_dec(v_x_1809_);
lean_inc_ref(v___x_1842_);
v___x_1843_ = l_Rat_instDecidableLe(v_lo_1840_, v___x_1842_);
if (v___x_1843_ == 0)
{
uint8_t v___x_1844_; 
lean_dec_ref(v___x_1842_);
lean_dec_ref(v_hi_1841_);
v___x_1844_ = 1;
return v___x_1844_;
}
else
{
uint8_t v___x_1845_; 
v___x_1845_ = l_Rat_instDecidableLe(v___x_1842_, v_hi_1841_);
if (v___x_1845_ == 0)
{
uint8_t v___x_1846_; 
v___x_1846_ = 1;
return v___x_1846_;
}
else
{
uint8_t v___x_1847_; 
v___x_1847_ = 0;
return v___x_1847_;
}
}
}
case 4:
{
lean_object* v_enzymeSites_1848_; uint8_t v___x_1849_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v_enzymeSites_1848_ = lean_ctor_get(v_x_1808_, 0);
lean_inc(v_enzymeSites_1848_);
lean_dec_ref(v_x_1808_);
v___x_1849_ = lp_BioCompiler_BioCompiler_hasAnyRestrictionSite(v_x_1809_, v_enzymeSites_1848_);
if (v___x_1849_ == 0)
{
uint8_t v___x_1850_; 
v___x_1850_ = 0;
return v___x_1850_;
}
else
{
uint8_t v___x_1851_; 
v___x_1851_ = 1;
return v___x_1851_;
}
}
case 5:
{
lean_object* v_readingFrame_1852_; lean_object* v_exonBoundaries_1853_; uint8_t v___x_1854_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v_readingFrame_1852_ = lean_ctor_get(v_x_1808_, 0);
lean_inc_n(v_readingFrame_1852_, 2);
v_exonBoundaries_1853_ = lean_ctor_get(v_x_1808_, 1);
lean_inc(v_exonBoundaries_1853_);
lean_dec_ref(v_x_1808_);
v___x_1854_ = lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent(v_exonBoundaries_1853_, v_readingFrame_1852_);
if (v___x_1854_ == 0)
{
uint8_t v___x_1855_; 
lean_dec(v_readingFrame_1852_);
lean_dec(v_x_1809_);
v___x_1855_ = 1;
return v___x_1855_;
}
else
{
uint8_t v___x_1856_; 
v___x_1856_ = lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(v_x_1809_, v_readingFrame_1852_);
lean_dec(v_readingFrame_1852_);
if (v___x_1856_ == 0)
{
uint8_t v___x_1857_; 
v___x_1857_ = 0;
return v___x_1857_;
}
else
{
uint8_t v___x_1858_; 
v___x_1858_ = 1;
return v___x_1858_;
}
}
}
case 6:
{
uint8_t v___x_1859_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v___x_1859_ = lp_BioCompiler_BioCompiler_hasInstabilityMotif(v_x_1809_);
if (v___x_1859_ == 0)
{
uint8_t v___x_1860_; 
v___x_1860_ = 0;
return v___x_1860_;
}
else
{
uint8_t v___x_1861_; 
v___x_1861_ = 1;
return v___x_1861_;
}
}
case 7:
{
lean_object* v___x_1862_; uint8_t v___x_1863_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v___x_1862_ = lean_apply_1(v_inst_1805_, v_x_1809_);
v___x_1863_ = lean_unbox(v___x_1862_);
if (v___x_1863_ == 0)
{
uint8_t v___x_1864_; 
v___x_1864_ = 0;
return v___x_1864_;
}
else
{
uint8_t v___x_1865_; 
v___x_1865_ = 1;
return v___x_1865_;
}
}
case 8:
{
lean_object* v___x_1866_; uint8_t v___x_1867_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v___x_1866_ = lp_BioCompiler_BioCompiler_spliceDonorConsensus;
v___x_1867_ = lp_BioCompiler_BioCompiler_Sequence_containsPattern(v_x_1809_, v___x_1866_);
if (v___x_1867_ == 0)
{
uint8_t v___x_1868_; 
v___x_1868_ = 0;
return v___x_1868_;
}
else
{
uint8_t v___x_1869_; 
v___x_1869_ = 1;
return v___x_1869_;
}
}
case 9:
{
lean_object* v___x_1870_; uint8_t v___x_1871_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v___x_1870_ = lean_unsigned_to_nat(0u);
v___x_1871_ = lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(v_x_1809_, v___x_1870_);
if (v___x_1871_ == 0)
{
uint8_t v___x_1872_; 
v___x_1872_ = 0;
return v___x_1872_;
}
else
{
uint8_t v___x_1873_; 
v___x_1873_ = 1;
return v___x_1873_;
}
}
case 10:
{
uint8_t v___x_1874_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v___x_1874_ = lp_BioCompiler_BioCompiler_isValidCodingSeq(v_x_1809_);
if (v___x_1874_ == 0)
{
uint8_t v___x_1875_; 
v___x_1875_ = 1;
return v___x_1875_;
}
else
{
uint8_t v___x_1876_; 
v___x_1876_ = 0;
return v___x_1876_;
}
}
case 11:
{
lean_object* v_organism_1877_; lean_object* v_threshold_1878_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1803_);
v_organism_1877_ = lean_ctor_get(v_x_1808_, 0);
lean_inc_ref(v_organism_1877_);
v_threshold_1878_ = lean_ctor_get(v_x_1808_, 1);
lean_inc_ref(v_threshold_1878_);
lean_dec_ref(v_x_1808_);
v_org_1812_ = v_organism_1877_;
v_threshold_1813_ = v_threshold_1878_;
v_seq_1814_ = v_x_1809_;
goto v___jp_1811_;
}
case 12:
{
lean_object* v_organism_1879_; lean_object* v_threshold_1880_; lean_object* v_hasCrypticPromoter_1881_; lean_object* v_hasBorderlinePromoter_1882_; lean_object* v___x_1883_; uint8_t v___x_1884_; 
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v_organism_1879_ = lean_ctor_get(v_x_1808_, 0);
lean_inc_ref_n(v_organism_1879_, 2);
v_threshold_1880_ = lean_ctor_get(v_x_1808_, 1);
lean_inc_ref_n(v_threshold_1880_, 2);
lean_dec_ref(v_x_1808_);
v_hasCrypticPromoter_1881_ = lean_ctor_get(v_inst_1806_, 0);
lean_inc_ref(v_hasCrypticPromoter_1881_);
v_hasBorderlinePromoter_1882_ = lean_ctor_get(v_inst_1806_, 1);
lean_inc_ref(v_hasBorderlinePromoter_1882_);
lean_dec_ref(v_inst_1806_);
lean_inc(v_x_1809_);
v___x_1883_ = lean_apply_3(v_hasCrypticPromoter_1881_, v_x_1809_, v_organism_1879_, v_threshold_1880_);
v___x_1884_ = lean_unbox(v___x_1883_);
if (v___x_1884_ == 0)
{
lean_object* v___x_1885_; uint8_t v___x_1886_; 
v___x_1885_ = lean_apply_3(v_hasBorderlinePromoter_1882_, v_x_1809_, v_organism_1879_, v_threshold_1880_);
v___x_1886_ = lean_unbox(v___x_1885_);
if (v___x_1886_ == 0)
{
uint8_t v___x_1887_; 
v___x_1887_ = 0;
return v___x_1887_;
}
else
{
uint8_t v___x_1888_; 
v___x_1888_ = 2;
return v___x_1888_;
}
}
else
{
uint8_t v___x_1889_; 
lean_dec_ref(v_hasBorderlinePromoter_1882_);
lean_dec_ref(v_threshold_1880_);
lean_dec_ref(v_organism_1879_);
lean_dec(v_x_1809_);
v___x_1889_ = 1;
return v___x_1889_;
}
}
default: 
{
uint8_t v___x_1890_; 
lean_dec(v_x_1809_);
lean_dec(v_x_1808_);
lean_dec_ref(v_inst_1807_);
lean_dec_ref(v_inst_1806_);
lean_dec_ref(v_inst_1805_);
lean_dec_ref(v_inst_1804_);
lean_dec_ref(v_inst_1803_);
v___x_1890_ = 2;
return v___x_1890_;
}
}
v___jp_1811_:
{
lean_object* v_computeCAI_1815_; lean_object* v___x_1816_; uint8_t v___x_1817_; 
v_computeCAI_1815_ = lean_ctor_get(v_inst_1804_, 0);
lean_inc_ref(v_computeCAI_1815_);
lean_dec_ref(v_inst_1804_);
v___x_1816_ = lean_apply_2(v_computeCAI_1815_, v_seq_1814_, v_org_1812_);
v___x_1817_ = l_Rat_instDecidableLe(v_threshold_1813_, v___x_1816_);
if (v___x_1817_ == 0)
{
uint8_t v___x_1818_; 
v___x_1818_ = 1;
return v___x_1818_;
}
else
{
uint8_t v___x_1819_; 
v___x_1819_ = 0;
return v___x_1819_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluate___redArg___boxed(lean_object* v_inst_1891_, lean_object* v_inst_1892_, lean_object* v_inst_1893_, lean_object* v_inst_1894_, lean_object* v_inst_1895_, lean_object* v_x_1896_, lean_object* v_x_1897_, lean_object* v_x_1898_){
_start:
{
uint8_t v_res_1899_; lean_object* v_r_1900_; 
v_res_1899_ = lp_BioCompiler_BioCompiler_evaluate___redArg(v_inst_1891_, v_inst_1892_, v_inst_1893_, v_inst_1894_, v_inst_1895_, v_x_1896_, v_x_1897_, v_x_1898_);
lean_dec_ref(v_x_1898_);
v_r_1900_ = lean_box(v_res_1899_);
return v_r_1900_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluate(lean_object* v_inst_1901_, lean_object* v_inst_1902_, lean_object* v_inst_1903_, lean_object* v_inst_1904_, lean_object* v_State_1905_, lean_object* v_inst_1906_, lean_object* v_inst_1907_, lean_object* v_inst_1908_, lean_object* v_x_1909_, lean_object* v_x_1910_, lean_object* v_x_1911_){
_start:
{
uint8_t v___x_1912_; 
v___x_1912_ = lp_BioCompiler_BioCompiler_evaluate___redArg(v_inst_1901_, v_inst_1902_, v_inst_1903_, v_inst_1904_, v_inst_1908_, v_x_1909_, v_x_1910_, v_x_1911_);
return v___x_1912_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluate___boxed(lean_object* v_inst_1913_, lean_object* v_inst_1914_, lean_object* v_inst_1915_, lean_object* v_inst_1916_, lean_object* v_State_1917_, lean_object* v_inst_1918_, lean_object* v_inst_1919_, lean_object* v_inst_1920_, lean_object* v_x_1921_, lean_object* v_x_1922_, lean_object* v_x_1923_){
_start:
{
uint8_t v_res_1924_; lean_object* v_r_1925_; 
v_res_1924_ = lp_BioCompiler_BioCompiler_evaluate(v_inst_1913_, v_inst_1914_, v_inst_1915_, v_inst_1916_, v_State_1917_, v_inst_1918_, v_inst_1919_, v_inst_1920_, v_x_1921_, v_x_1922_, v_x_1923_);
lean_dec_ref(v_x_1923_);
lean_dec(v_inst_1919_);
lean_dec_ref(v_inst_1918_);
v_r_1925_ = lean_box(v_res_1924_);
return v_r_1925_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter___redArg(lean_object* v_x_1926_, lean_object* v_x_1927_, lean_object* v_x_1928_, lean_object* v_h__1_1929_, lean_object* v_h__2_1930_, lean_object* v_h__3_1931_, lean_object* v_h__4_1932_, lean_object* v_h__5_1933_, lean_object* v_h__6_1934_, lean_object* v_h__7_1935_, lean_object* v_h__8_1936_, lean_object* v_h__9_1937_, lean_object* v_h__10_1938_, lean_object* v_h__11_1939_, lean_object* v_h__12_1940_, lean_object* v_h__13_1941_, lean_object* v_h__14_1942_, lean_object* v_h__15_1943_, lean_object* v_h__16_1944_, lean_object* v_h__17_1945_, lean_object* v_h__18_1946_, lean_object* v_h__19_1947_, lean_object* v_h__20_1948_, lean_object* v_h__21_1949_, lean_object* v_h__22_1950_, lean_object* v_h__23_1951_, lean_object* v_h__24_1952_, lean_object* v_h__25_1953_, lean_object* v_h__26_1954_, lean_object* v_h__27_1955_, lean_object* v_h__28_1956_, lean_object* v_h__29_1957_, lean_object* v_h__30_1958_, lean_object* v_h__31_1959_, lean_object* v_h__32_1960_, lean_object* v_h__33_1961_){
_start:
{
switch(lean_obj_tag(v_x_1926_))
{
case 0:
{
lean_object* v_cellType_1962_; lean_object* v___x_1963_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
v_cellType_1962_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_cellType_1962_);
lean_dec_ref(v_x_1926_);
v___x_1963_ = lean_apply_3(v_h__1_1929_, v_cellType_1962_, v_x_1927_, v_x_1928_);
return v___x_1963_;
}
case 1:
{
lean_object* v___x_1964_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__1_1929_);
v___x_1964_ = lean_apply_2(v_h__2_1930_, v_x_1927_, v_x_1928_);
return v___x_1964_;
}
case 2:
{
lean_object* v_organism_1965_; lean_object* v_threshold_1966_; lean_object* v___x_1967_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_organism_1965_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_organism_1965_);
v_threshold_1966_ = lean_ctor_get(v_x_1926_, 1);
lean_inc_ref(v_threshold_1966_);
lean_dec_ref(v_x_1926_);
v___x_1967_ = lean_apply_4(v_h__3_1931_, v_organism_1965_, v_threshold_1966_, v_x_1927_, v_x_1928_);
return v___x_1967_;
}
case 3:
{
lean_object* v_lo_1968_; lean_object* v_hi_1969_; lean_object* v___x_1970_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_lo_1968_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_lo_1968_);
v_hi_1969_ = lean_ctor_get(v_x_1926_, 1);
lean_inc_ref(v_hi_1969_);
lean_dec_ref(v_x_1926_);
v___x_1970_ = lean_apply_4(v_h__4_1932_, v_lo_1968_, v_hi_1969_, v_x_1927_, v_x_1928_);
return v___x_1970_;
}
case 4:
{
lean_object* v_enzymeSites_1971_; lean_object* v___x_1972_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_enzymeSites_1971_ = lean_ctor_get(v_x_1926_, 0);
lean_inc(v_enzymeSites_1971_);
lean_dec_ref(v_x_1926_);
v___x_1972_ = lean_apply_3(v_h__5_1933_, v_enzymeSites_1971_, v_x_1927_, v_x_1928_);
return v___x_1972_;
}
case 5:
{
lean_object* v_readingFrame_1973_; lean_object* v_exonBoundaries_1974_; lean_object* v___x_1975_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_readingFrame_1973_ = lean_ctor_get(v_x_1926_, 0);
lean_inc(v_readingFrame_1973_);
v_exonBoundaries_1974_ = lean_ctor_get(v_x_1926_, 1);
lean_inc(v_exonBoundaries_1974_);
lean_dec_ref(v_x_1926_);
v___x_1975_ = lean_apply_4(v_h__6_1934_, v_readingFrame_1973_, v_exonBoundaries_1974_, v_x_1927_, v_x_1928_);
return v___x_1975_;
}
case 6:
{
lean_object* v___x_1976_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_1976_ = lean_apply_2(v_h__7_1935_, v_x_1927_, v_x_1928_);
return v___x_1976_;
}
case 7:
{
lean_object* v___x_1977_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_1977_ = lean_apply_2(v_h__8_1936_, v_x_1927_, v_x_1928_);
return v___x_1977_;
}
case 8:
{
lean_object* v___x_1978_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_1978_ = lean_apply_2(v_h__9_1937_, v_x_1927_, v_x_1928_);
return v___x_1978_;
}
case 9:
{
lean_object* v___x_1979_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_1979_ = lean_apply_2(v_h__10_1938_, v_x_1927_, v_x_1928_);
return v___x_1979_;
}
case 10:
{
lean_object* v___x_1980_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_1980_ = lean_apply_2(v_h__11_1939_, v_x_1927_, v_x_1928_);
return v___x_1980_;
}
case 11:
{
lean_object* v_organism_1981_; lean_object* v_threshold_1982_; lean_object* v___x_1983_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_organism_1981_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_organism_1981_);
v_threshold_1982_ = lean_ctor_get(v_x_1926_, 1);
lean_inc_ref(v_threshold_1982_);
lean_dec_ref(v_x_1926_);
v___x_1983_ = lean_apply_4(v_h__12_1940_, v_organism_1981_, v_threshold_1982_, v_x_1927_, v_x_1928_);
return v___x_1983_;
}
case 12:
{
lean_object* v_organism_1984_; lean_object* v_threshold_1985_; lean_object* v___x_1986_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_organism_1984_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_organism_1984_);
v_threshold_1985_ = lean_ctor_get(v_x_1926_, 1);
lean_inc_ref(v_threshold_1985_);
lean_dec_ref(v_x_1926_);
v___x_1986_ = lean_apply_4(v_h__13_1941_, v_organism_1984_, v_threshold_1985_, v_x_1927_, v_x_1928_);
return v___x_1986_;
}
case 13:
{
lean_object* v_minScore_1987_; lean_object* v___x_1988_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_minScore_1987_ = lean_ctor_get(v_x_1926_, 0);
lean_inc(v_minScore_1987_);
lean_dec_ref(v_x_1926_);
v___x_1988_ = lean_apply_3(v_h__14_1942_, v_minScore_1987_, v_x_1927_, v_x_1928_);
return v___x_1988_;
}
case 14:
{
uint8_t v_isCytosolic_1989_; lean_object* v_threshold_1990_; lean_object* v___x_1991_; lean_object* v___x_1992_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_isCytosolic_1989_ = lean_ctor_get_uint8(v_x_1926_, sizeof(void*)*1);
v_threshold_1990_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_threshold_1990_);
lean_dec_ref(v_x_1926_);
v___x_1991_ = lean_box(v_isCytosolic_1989_);
v___x_1992_ = lean_apply_4(v_h__15_1943_, v___x_1991_, v_threshold_1990_, v_x_1927_, v_x_1928_);
return v___x_1992_;
}
case 15:
{
lean_object* v_dgThreshold_1993_; lean_object* v___x_1994_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_dgThreshold_1993_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_dgThreshold_1993_);
lean_dec_ref(v_x_1926_);
v___x_1994_ = lean_apply_3(v_h__16_1944_, v_dgThreshold_1993_, v_x_1927_, v_x_1928_);
return v___x_1994_;
}
case 16:
{
lean_object* v_organism_1995_; lean_object* v___x_1996_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_organism_1995_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_organism_1995_);
lean_dec_ref(v_x_1926_);
v___x_1996_ = lean_apply_3(v_h__17_1945_, v_organism_1995_, v_x_1927_, v_x_1928_);
return v___x_1996_;
}
case 17:
{
lean_object* v_threshold_1997_; lean_object* v___x_1998_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_threshold_1997_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_threshold_1997_);
lean_dec_ref(v_x_1926_);
v___x_1998_ = lean_apply_3(v_h__18_1946_, v_threshold_1997_, v_x_1927_, v_x_1928_);
return v___x_1998_;
}
case 18:
{
lean_object* v___x_1999_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_1999_ = lean_apply_2(v_h__19_1947_, v_x_1927_, v_x_1928_);
return v___x_1999_;
}
case 19:
{
lean_object* v___x_2000_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_2000_ = lean_apply_2(v_h__20_1948_, v_x_1927_, v_x_1928_);
return v___x_2000_;
}
case 20:
{
lean_object* v___x_2001_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_2001_ = lean_apply_2(v_h__21_1949_, v_x_1927_, v_x_1928_);
return v___x_2001_;
}
case 21:
{
lean_object* v_ddgThreshold_2002_; lean_object* v___x_2003_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_ddgThreshold_2002_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_ddgThreshold_2002_);
lean_dec_ref(v_x_1926_);
v___x_2003_ = lean_apply_3(v_h__22_1950_, v_ddgThreshold_2002_, v_x_1927_, v_x_1928_);
return v___x_2003_;
}
case 22:
{
lean_object* v_maxDDG_2004_; lean_object* v___x_2005_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_maxDDG_2004_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_maxDDG_2004_);
lean_dec_ref(v_x_1926_);
v___x_2005_ = lean_apply_3(v_h__23_1951_, v_maxDDG_2004_, v_x_1927_, v_x_1928_);
return v___x_2005_;
}
case 23:
{
lean_object* v___x_2006_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_2006_ = lean_apply_2(v_h__24_1952_, v_x_1927_, v_x_1928_);
return v___x_2006_;
}
case 24:
{
lean_object* v_threshold_2007_; lean_object* v___x_2008_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_threshold_2007_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_threshold_2007_);
lean_dec_ref(v_x_1926_);
v___x_2008_ = lean_apply_3(v_h__25_1953_, v_threshold_2007_, v_x_1927_, v_x_1928_);
return v___x_2008_;
}
case 25:
{
lean_object* v_minScore_2009_; lean_object* v___x_2010_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_minScore_2009_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_minScore_2009_);
lean_dec_ref(v_x_1926_);
v___x_2010_ = lean_apply_3(v_h__26_1954_, v_minScore_2009_, v_x_1927_, v_x_1928_);
return v___x_2010_;
}
case 26:
{
lean_object* v___x_2011_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v___x_2011_ = lean_apply_2(v_h__27_1955_, v_x_1927_, v_x_1928_);
return v___x_2011_;
}
case 27:
{
lean_object* v_pILo_2012_; lean_object* v_pIHi_2013_; lean_object* v___x_2014_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_pILo_2012_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_pILo_2012_);
v_pIHi_2013_ = lean_ctor_get(v_x_1926_, 1);
lean_inc_ref(v_pIHi_2013_);
lean_dec_ref(v_x_1926_);
v___x_2014_ = lean_apply_4(v_h__28_1956_, v_pILo_2012_, v_pIHi_2013_, v_x_1927_, v_x_1928_);
return v___x_2014_;
}
case 28:
{
lean_object* v_maxLen_2015_; lean_object* v___x_2016_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_maxLen_2015_ = lean_ctor_get(v_x_1926_, 0);
lean_inc(v_maxLen_2015_);
lean_dec_ref(v_x_1926_);
v___x_2016_ = lean_apply_3(v_h__29_1957_, v_maxLen_2015_, v_x_1927_, v_x_1928_);
return v___x_2016_;
}
case 29:
{
lean_object* v_maxScore_2017_; lean_object* v___x_2018_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_maxScore_2017_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_maxScore_2017_);
lean_dec_ref(v_x_1926_);
v___x_2018_ = lean_apply_3(v_h__30_1958_, v_maxScore_2017_, v_x_1927_, v_x_1928_);
return v___x_2018_;
}
case 30:
{
lean_object* v_ic50Threshold_2019_; lean_object* v___x_2020_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__32_1960_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_ic50Threshold_2019_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_ic50Threshold_2019_);
lean_dec_ref(v_x_1926_);
v___x_2020_ = lean_apply_3(v_h__31_1959_, v_ic50Threshold_2019_, v_x_1927_, v_x_1928_);
return v___x_2020_;
}
case 31:
{
lean_object* v_scoreThreshold_2021_; lean_object* v___x_2022_; 
lean_dec(v_h__33_1961_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_scoreThreshold_2021_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_scoreThreshold_2021_);
lean_dec_ref(v_x_1926_);
v___x_2022_ = lean_apply_3(v_h__32_1960_, v_scoreThreshold_2021_, v_x_1927_, v_x_1928_);
return v___x_2022_;
}
default: 
{
lean_object* v_maxCoverage_2023_; lean_object* v___x_2024_; 
lean_dec(v_h__32_1960_);
lean_dec(v_h__31_1959_);
lean_dec(v_h__30_1958_);
lean_dec(v_h__29_1957_);
lean_dec(v_h__28_1956_);
lean_dec(v_h__27_1955_);
lean_dec(v_h__26_1954_);
lean_dec(v_h__25_1953_);
lean_dec(v_h__24_1952_);
lean_dec(v_h__23_1951_);
lean_dec(v_h__22_1950_);
lean_dec(v_h__21_1949_);
lean_dec(v_h__20_1948_);
lean_dec(v_h__19_1947_);
lean_dec(v_h__18_1946_);
lean_dec(v_h__17_1945_);
lean_dec(v_h__16_1944_);
lean_dec(v_h__15_1943_);
lean_dec(v_h__14_1942_);
lean_dec(v_h__13_1941_);
lean_dec(v_h__12_1940_);
lean_dec(v_h__11_1939_);
lean_dec(v_h__10_1938_);
lean_dec(v_h__9_1937_);
lean_dec(v_h__8_1936_);
lean_dec(v_h__7_1935_);
lean_dec(v_h__6_1934_);
lean_dec(v_h__5_1933_);
lean_dec(v_h__4_1932_);
lean_dec(v_h__3_1931_);
lean_dec(v_h__2_1930_);
lean_dec(v_h__1_1929_);
v_maxCoverage_2023_ = lean_ctor_get(v_x_1926_, 0);
lean_inc_ref(v_maxCoverage_2023_);
lean_dec_ref(v_x_1926_);
v___x_2024_ = lean_apply_3(v_h__33_1961_, v_maxCoverage_2023_, v_x_1927_, v_x_1928_);
return v___x_2024_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter___redArg___boxed(lean_object** _args){
lean_object* v_x_2025_ = _args[0];
lean_object* v_x_2026_ = _args[1];
lean_object* v_x_2027_ = _args[2];
lean_object* v_h__1_2028_ = _args[3];
lean_object* v_h__2_2029_ = _args[4];
lean_object* v_h__3_2030_ = _args[5];
lean_object* v_h__4_2031_ = _args[6];
lean_object* v_h__5_2032_ = _args[7];
lean_object* v_h__6_2033_ = _args[8];
lean_object* v_h__7_2034_ = _args[9];
lean_object* v_h__8_2035_ = _args[10];
lean_object* v_h__9_2036_ = _args[11];
lean_object* v_h__10_2037_ = _args[12];
lean_object* v_h__11_2038_ = _args[13];
lean_object* v_h__12_2039_ = _args[14];
lean_object* v_h__13_2040_ = _args[15];
lean_object* v_h__14_2041_ = _args[16];
lean_object* v_h__15_2042_ = _args[17];
lean_object* v_h__16_2043_ = _args[18];
lean_object* v_h__17_2044_ = _args[19];
lean_object* v_h__18_2045_ = _args[20];
lean_object* v_h__19_2046_ = _args[21];
lean_object* v_h__20_2047_ = _args[22];
lean_object* v_h__21_2048_ = _args[23];
lean_object* v_h__22_2049_ = _args[24];
lean_object* v_h__23_2050_ = _args[25];
lean_object* v_h__24_2051_ = _args[26];
lean_object* v_h__25_2052_ = _args[27];
lean_object* v_h__26_2053_ = _args[28];
lean_object* v_h__27_2054_ = _args[29];
lean_object* v_h__28_2055_ = _args[30];
lean_object* v_h__29_2056_ = _args[31];
lean_object* v_h__30_2057_ = _args[32];
lean_object* v_h__31_2058_ = _args[33];
lean_object* v_h__32_2059_ = _args[34];
lean_object* v_h__33_2060_ = _args[35];
_start:
{
lean_object* v_res_2061_; 
v_res_2061_ = lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter___redArg(v_x_2025_, v_x_2026_, v_x_2027_, v_h__1_2028_, v_h__2_2029_, v_h__3_2030_, v_h__4_2031_, v_h__5_2032_, v_h__6_2033_, v_h__7_2034_, v_h__8_2035_, v_h__9_2036_, v_h__10_2037_, v_h__11_2038_, v_h__12_2039_, v_h__13_2040_, v_h__14_2041_, v_h__15_2042_, v_h__16_2043_, v_h__17_2044_, v_h__18_2045_, v_h__19_2046_, v_h__20_2047_, v_h__21_2048_, v_h__22_2049_, v_h__23_2050_, v_h__24_2051_, v_h__25_2052_, v_h__26_2053_, v_h__27_2054_, v_h__28_2055_, v_h__29_2056_, v_h__30_2057_, v_h__31_2058_, v_h__32_2059_, v_h__33_2060_);
return v_res_2061_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter(lean_object* v_motive_2062_, lean_object* v_x_2063_, lean_object* v_x_2064_, lean_object* v_x_2065_, lean_object* v_h__1_2066_, lean_object* v_h__2_2067_, lean_object* v_h__3_2068_, lean_object* v_h__4_2069_, lean_object* v_h__5_2070_, lean_object* v_h__6_2071_, lean_object* v_h__7_2072_, lean_object* v_h__8_2073_, lean_object* v_h__9_2074_, lean_object* v_h__10_2075_, lean_object* v_h__11_2076_, lean_object* v_h__12_2077_, lean_object* v_h__13_2078_, lean_object* v_h__14_2079_, lean_object* v_h__15_2080_, lean_object* v_h__16_2081_, lean_object* v_h__17_2082_, lean_object* v_h__18_2083_, lean_object* v_h__19_2084_, lean_object* v_h__20_2085_, lean_object* v_h__21_2086_, lean_object* v_h__22_2087_, lean_object* v_h__23_2088_, lean_object* v_h__24_2089_, lean_object* v_h__25_2090_, lean_object* v_h__26_2091_, lean_object* v_h__27_2092_, lean_object* v_h__28_2093_, lean_object* v_h__29_2094_, lean_object* v_h__30_2095_, lean_object* v_h__31_2096_, lean_object* v_h__32_2097_, lean_object* v_h__33_2098_){
_start:
{
switch(lean_obj_tag(v_x_2063_))
{
case 0:
{
lean_object* v_cellType_2099_; lean_object* v___x_2100_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
v_cellType_2099_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_cellType_2099_);
lean_dec_ref(v_x_2063_);
v___x_2100_ = lean_apply_3(v_h__1_2066_, v_cellType_2099_, v_x_2064_, v_x_2065_);
return v___x_2100_;
}
case 1:
{
lean_object* v___x_2101_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__1_2066_);
v___x_2101_ = lean_apply_2(v_h__2_2067_, v_x_2064_, v_x_2065_);
return v___x_2101_;
}
case 2:
{
lean_object* v_organism_2102_; lean_object* v_threshold_2103_; lean_object* v___x_2104_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_organism_2102_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_organism_2102_);
v_threshold_2103_ = lean_ctor_get(v_x_2063_, 1);
lean_inc_ref(v_threshold_2103_);
lean_dec_ref(v_x_2063_);
v___x_2104_ = lean_apply_4(v_h__3_2068_, v_organism_2102_, v_threshold_2103_, v_x_2064_, v_x_2065_);
return v___x_2104_;
}
case 3:
{
lean_object* v_lo_2105_; lean_object* v_hi_2106_; lean_object* v___x_2107_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_lo_2105_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_lo_2105_);
v_hi_2106_ = lean_ctor_get(v_x_2063_, 1);
lean_inc_ref(v_hi_2106_);
lean_dec_ref(v_x_2063_);
v___x_2107_ = lean_apply_4(v_h__4_2069_, v_lo_2105_, v_hi_2106_, v_x_2064_, v_x_2065_);
return v___x_2107_;
}
case 4:
{
lean_object* v_enzymeSites_2108_; lean_object* v___x_2109_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_enzymeSites_2108_ = lean_ctor_get(v_x_2063_, 0);
lean_inc(v_enzymeSites_2108_);
lean_dec_ref(v_x_2063_);
v___x_2109_ = lean_apply_3(v_h__5_2070_, v_enzymeSites_2108_, v_x_2064_, v_x_2065_);
return v___x_2109_;
}
case 5:
{
lean_object* v_readingFrame_2110_; lean_object* v_exonBoundaries_2111_; lean_object* v___x_2112_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_readingFrame_2110_ = lean_ctor_get(v_x_2063_, 0);
lean_inc(v_readingFrame_2110_);
v_exonBoundaries_2111_ = lean_ctor_get(v_x_2063_, 1);
lean_inc(v_exonBoundaries_2111_);
lean_dec_ref(v_x_2063_);
v___x_2112_ = lean_apply_4(v_h__6_2071_, v_readingFrame_2110_, v_exonBoundaries_2111_, v_x_2064_, v_x_2065_);
return v___x_2112_;
}
case 6:
{
lean_object* v___x_2113_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2113_ = lean_apply_2(v_h__7_2072_, v_x_2064_, v_x_2065_);
return v___x_2113_;
}
case 7:
{
lean_object* v___x_2114_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2114_ = lean_apply_2(v_h__8_2073_, v_x_2064_, v_x_2065_);
return v___x_2114_;
}
case 8:
{
lean_object* v___x_2115_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2115_ = lean_apply_2(v_h__9_2074_, v_x_2064_, v_x_2065_);
return v___x_2115_;
}
case 9:
{
lean_object* v___x_2116_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2116_ = lean_apply_2(v_h__10_2075_, v_x_2064_, v_x_2065_);
return v___x_2116_;
}
case 10:
{
lean_object* v___x_2117_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2117_ = lean_apply_2(v_h__11_2076_, v_x_2064_, v_x_2065_);
return v___x_2117_;
}
case 11:
{
lean_object* v_organism_2118_; lean_object* v_threshold_2119_; lean_object* v___x_2120_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_organism_2118_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_organism_2118_);
v_threshold_2119_ = lean_ctor_get(v_x_2063_, 1);
lean_inc_ref(v_threshold_2119_);
lean_dec_ref(v_x_2063_);
v___x_2120_ = lean_apply_4(v_h__12_2077_, v_organism_2118_, v_threshold_2119_, v_x_2064_, v_x_2065_);
return v___x_2120_;
}
case 12:
{
lean_object* v_organism_2121_; lean_object* v_threshold_2122_; lean_object* v___x_2123_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_organism_2121_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_organism_2121_);
v_threshold_2122_ = lean_ctor_get(v_x_2063_, 1);
lean_inc_ref(v_threshold_2122_);
lean_dec_ref(v_x_2063_);
v___x_2123_ = lean_apply_4(v_h__13_2078_, v_organism_2121_, v_threshold_2122_, v_x_2064_, v_x_2065_);
return v___x_2123_;
}
case 13:
{
lean_object* v_minScore_2124_; lean_object* v___x_2125_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_minScore_2124_ = lean_ctor_get(v_x_2063_, 0);
lean_inc(v_minScore_2124_);
lean_dec_ref(v_x_2063_);
v___x_2125_ = lean_apply_3(v_h__14_2079_, v_minScore_2124_, v_x_2064_, v_x_2065_);
return v___x_2125_;
}
case 14:
{
uint8_t v_isCytosolic_2126_; lean_object* v_threshold_2127_; lean_object* v___x_2128_; lean_object* v___x_2129_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_isCytosolic_2126_ = lean_ctor_get_uint8(v_x_2063_, sizeof(void*)*1);
v_threshold_2127_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_threshold_2127_);
lean_dec_ref(v_x_2063_);
v___x_2128_ = lean_box(v_isCytosolic_2126_);
v___x_2129_ = lean_apply_4(v_h__15_2080_, v___x_2128_, v_threshold_2127_, v_x_2064_, v_x_2065_);
return v___x_2129_;
}
case 15:
{
lean_object* v_dgThreshold_2130_; lean_object* v___x_2131_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_dgThreshold_2130_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_dgThreshold_2130_);
lean_dec_ref(v_x_2063_);
v___x_2131_ = lean_apply_3(v_h__16_2081_, v_dgThreshold_2130_, v_x_2064_, v_x_2065_);
return v___x_2131_;
}
case 16:
{
lean_object* v_organism_2132_; lean_object* v___x_2133_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_organism_2132_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_organism_2132_);
lean_dec_ref(v_x_2063_);
v___x_2133_ = lean_apply_3(v_h__17_2082_, v_organism_2132_, v_x_2064_, v_x_2065_);
return v___x_2133_;
}
case 17:
{
lean_object* v_threshold_2134_; lean_object* v___x_2135_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_threshold_2134_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_threshold_2134_);
lean_dec_ref(v_x_2063_);
v___x_2135_ = lean_apply_3(v_h__18_2083_, v_threshold_2134_, v_x_2064_, v_x_2065_);
return v___x_2135_;
}
case 18:
{
lean_object* v___x_2136_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2136_ = lean_apply_2(v_h__19_2084_, v_x_2064_, v_x_2065_);
return v___x_2136_;
}
case 19:
{
lean_object* v___x_2137_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2137_ = lean_apply_2(v_h__20_2085_, v_x_2064_, v_x_2065_);
return v___x_2137_;
}
case 20:
{
lean_object* v___x_2138_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2138_ = lean_apply_2(v_h__21_2086_, v_x_2064_, v_x_2065_);
return v___x_2138_;
}
case 21:
{
lean_object* v_ddgThreshold_2139_; lean_object* v___x_2140_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_ddgThreshold_2139_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_ddgThreshold_2139_);
lean_dec_ref(v_x_2063_);
v___x_2140_ = lean_apply_3(v_h__22_2087_, v_ddgThreshold_2139_, v_x_2064_, v_x_2065_);
return v___x_2140_;
}
case 22:
{
lean_object* v_maxDDG_2141_; lean_object* v___x_2142_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_maxDDG_2141_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_maxDDG_2141_);
lean_dec_ref(v_x_2063_);
v___x_2142_ = lean_apply_3(v_h__23_2088_, v_maxDDG_2141_, v_x_2064_, v_x_2065_);
return v___x_2142_;
}
case 23:
{
lean_object* v___x_2143_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2143_ = lean_apply_2(v_h__24_2089_, v_x_2064_, v_x_2065_);
return v___x_2143_;
}
case 24:
{
lean_object* v_threshold_2144_; lean_object* v___x_2145_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_threshold_2144_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_threshold_2144_);
lean_dec_ref(v_x_2063_);
v___x_2145_ = lean_apply_3(v_h__25_2090_, v_threshold_2144_, v_x_2064_, v_x_2065_);
return v___x_2145_;
}
case 25:
{
lean_object* v_minScore_2146_; lean_object* v___x_2147_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_minScore_2146_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_minScore_2146_);
lean_dec_ref(v_x_2063_);
v___x_2147_ = lean_apply_3(v_h__26_2091_, v_minScore_2146_, v_x_2064_, v_x_2065_);
return v___x_2147_;
}
case 26:
{
lean_object* v___x_2148_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v___x_2148_ = lean_apply_2(v_h__27_2092_, v_x_2064_, v_x_2065_);
return v___x_2148_;
}
case 27:
{
lean_object* v_pILo_2149_; lean_object* v_pIHi_2150_; lean_object* v___x_2151_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_pILo_2149_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_pILo_2149_);
v_pIHi_2150_ = lean_ctor_get(v_x_2063_, 1);
lean_inc_ref(v_pIHi_2150_);
lean_dec_ref(v_x_2063_);
v___x_2151_ = lean_apply_4(v_h__28_2093_, v_pILo_2149_, v_pIHi_2150_, v_x_2064_, v_x_2065_);
return v___x_2151_;
}
case 28:
{
lean_object* v_maxLen_2152_; lean_object* v___x_2153_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_maxLen_2152_ = lean_ctor_get(v_x_2063_, 0);
lean_inc(v_maxLen_2152_);
lean_dec_ref(v_x_2063_);
v___x_2153_ = lean_apply_3(v_h__29_2094_, v_maxLen_2152_, v_x_2064_, v_x_2065_);
return v___x_2153_;
}
case 29:
{
lean_object* v_maxScore_2154_; lean_object* v___x_2155_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_maxScore_2154_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_maxScore_2154_);
lean_dec_ref(v_x_2063_);
v___x_2155_ = lean_apply_3(v_h__30_2095_, v_maxScore_2154_, v_x_2064_, v_x_2065_);
return v___x_2155_;
}
case 30:
{
lean_object* v_ic50Threshold_2156_; lean_object* v___x_2157_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__32_2097_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_ic50Threshold_2156_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_ic50Threshold_2156_);
lean_dec_ref(v_x_2063_);
v___x_2157_ = lean_apply_3(v_h__31_2096_, v_ic50Threshold_2156_, v_x_2064_, v_x_2065_);
return v___x_2157_;
}
case 31:
{
lean_object* v_scoreThreshold_2158_; lean_object* v___x_2159_; 
lean_dec(v_h__33_2098_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_scoreThreshold_2158_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_scoreThreshold_2158_);
lean_dec_ref(v_x_2063_);
v___x_2159_ = lean_apply_3(v_h__32_2097_, v_scoreThreshold_2158_, v_x_2064_, v_x_2065_);
return v___x_2159_;
}
default: 
{
lean_object* v_maxCoverage_2160_; lean_object* v___x_2161_; 
lean_dec(v_h__32_2097_);
lean_dec(v_h__31_2096_);
lean_dec(v_h__30_2095_);
lean_dec(v_h__29_2094_);
lean_dec(v_h__28_2093_);
lean_dec(v_h__27_2092_);
lean_dec(v_h__26_2091_);
lean_dec(v_h__25_2090_);
lean_dec(v_h__24_2089_);
lean_dec(v_h__23_2088_);
lean_dec(v_h__22_2087_);
lean_dec(v_h__21_2086_);
lean_dec(v_h__20_2085_);
lean_dec(v_h__19_2084_);
lean_dec(v_h__18_2083_);
lean_dec(v_h__17_2082_);
lean_dec(v_h__16_2081_);
lean_dec(v_h__15_2080_);
lean_dec(v_h__14_2079_);
lean_dec(v_h__13_2078_);
lean_dec(v_h__12_2077_);
lean_dec(v_h__11_2076_);
lean_dec(v_h__10_2075_);
lean_dec(v_h__9_2074_);
lean_dec(v_h__8_2073_);
lean_dec(v_h__7_2072_);
lean_dec(v_h__6_2071_);
lean_dec(v_h__5_2070_);
lean_dec(v_h__4_2069_);
lean_dec(v_h__3_2068_);
lean_dec(v_h__2_2067_);
lean_dec(v_h__1_2066_);
v_maxCoverage_2160_ = lean_ctor_get(v_x_2063_, 0);
lean_inc_ref(v_maxCoverage_2160_);
lean_dec_ref(v_x_2063_);
v___x_2161_ = lean_apply_3(v_h__33_2098_, v_maxCoverage_2160_, v_x_2064_, v_x_2065_);
return v___x_2161_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter___boxed(lean_object** _args){
lean_object* v_motive_2162_ = _args[0];
lean_object* v_x_2163_ = _args[1];
lean_object* v_x_2164_ = _args[2];
lean_object* v_x_2165_ = _args[3];
lean_object* v_h__1_2166_ = _args[4];
lean_object* v_h__2_2167_ = _args[5];
lean_object* v_h__3_2168_ = _args[6];
lean_object* v_h__4_2169_ = _args[7];
lean_object* v_h__5_2170_ = _args[8];
lean_object* v_h__6_2171_ = _args[9];
lean_object* v_h__7_2172_ = _args[10];
lean_object* v_h__8_2173_ = _args[11];
lean_object* v_h__9_2174_ = _args[12];
lean_object* v_h__10_2175_ = _args[13];
lean_object* v_h__11_2176_ = _args[14];
lean_object* v_h__12_2177_ = _args[15];
lean_object* v_h__13_2178_ = _args[16];
lean_object* v_h__14_2179_ = _args[17];
lean_object* v_h__15_2180_ = _args[18];
lean_object* v_h__16_2181_ = _args[19];
lean_object* v_h__17_2182_ = _args[20];
lean_object* v_h__18_2183_ = _args[21];
lean_object* v_h__19_2184_ = _args[22];
lean_object* v_h__20_2185_ = _args[23];
lean_object* v_h__21_2186_ = _args[24];
lean_object* v_h__22_2187_ = _args[25];
lean_object* v_h__23_2188_ = _args[26];
lean_object* v_h__24_2189_ = _args[27];
lean_object* v_h__25_2190_ = _args[28];
lean_object* v_h__26_2191_ = _args[29];
lean_object* v_h__27_2192_ = _args[30];
lean_object* v_h__28_2193_ = _args[31];
lean_object* v_h__29_2194_ = _args[32];
lean_object* v_h__30_2195_ = _args[33];
lean_object* v_h__31_2196_ = _args[34];
lean_object* v_h__32_2197_ = _args[35];
lean_object* v_h__33_2198_ = _args[36];
_start:
{
lean_object* v_res_2199_; 
v_res_2199_ = lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter(v_motive_2162_, v_x_2163_, v_x_2164_, v_x_2165_, v_h__1_2166_, v_h__2_2167_, v_h__3_2168_, v_h__4_2169_, v_h__5_2170_, v_h__6_2171_, v_h__7_2172_, v_h__8_2173_, v_h__9_2174_, v_h__10_2175_, v_h__11_2176_, v_h__12_2177_, v_h__13_2178_, v_h__14_2179_, v_h__15_2180_, v_h__16_2181_, v_h__17_2182_, v_h__18_2183_, v_h__19_2184_, v_h__20_2185_, v_h__21_2186_, v_h__22_2187_, v_h__23_2188_, v_h__24_2189_, v_h__25_2190_, v_h__26_2191_, v_h__27_2192_, v_h__28_2193_, v_h__29_2194_, v_h__30_2195_, v_h__31_2196_, v_h__32_2197_, v_h__33_2198_);
return v_res_2199_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__1_splitter___redArg(lean_object* v_x_2200_, lean_object* v_h__1_2201_, lean_object* v_h__2_2202_){
_start:
{
if (lean_obj_tag(v_x_2200_) == 1)
{
lean_object* v_tail_2203_; 
v_tail_2203_ = lean_ctor_get(v_x_2200_, 1);
if (lean_obj_tag(v_tail_2203_) == 0)
{
lean_object* v_head_2204_; lean_object* v___x_2205_; 
lean_dec(v_h__2_2202_);
v_head_2204_ = lean_ctor_get(v_x_2200_, 0);
lean_inc(v_head_2204_);
lean_dec_ref(v_x_2200_);
v___x_2205_ = lean_apply_1(v_h__1_2201_, v_head_2204_);
return v___x_2205_;
}
else
{
lean_object* v___x_2206_; 
lean_dec(v_h__1_2201_);
v___x_2206_ = lean_apply_2(v_h__2_2202_, v_x_2200_, lean_box(0));
return v___x_2206_;
}
}
else
{
lean_object* v___x_2207_; 
lean_dec(v_h__1_2201_);
v___x_2207_ = lean_apply_2(v_h__2_2202_, v_x_2200_, lean_box(0));
return v___x_2207_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__1_splitter(lean_object* v_motive_2208_, lean_object* v_x_2209_, lean_object* v_h__1_2210_, lean_object* v_h__2_2211_){
_start:
{
if (lean_obj_tag(v_x_2209_) == 1)
{
lean_object* v_tail_2212_; 
v_tail_2212_ = lean_ctor_get(v_x_2209_, 1);
if (lean_obj_tag(v_tail_2212_) == 0)
{
lean_object* v_head_2213_; lean_object* v___x_2214_; 
lean_dec(v_h__2_2211_);
v_head_2213_ = lean_ctor_get(v_x_2209_, 0);
lean_inc(v_head_2213_);
lean_dec_ref(v_x_2209_);
v___x_2214_ = lean_apply_1(v_h__1_2210_, v_head_2213_);
return v___x_2214_;
}
else
{
lean_object* v___x_2215_; 
lean_dec(v_h__1_2210_);
v___x_2215_ = lean_apply_2(v_h__2_2211_, v_x_2209_, lean_box(0));
return v___x_2215_;
}
}
else
{
lean_object* v___x_2216_; 
lean_dec(v_h__1_2210_);
v___x_2216_ = lean_apply_2(v_h__2_2211_, v_x_2209_, lean_box(0));
return v___x_2216_;
}
}
}
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_ThreeValued(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Sequence(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_NDFST(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Scanners(uint8_t builtin);
static bool _G_initialized = false;
LEAN_EXPORT lean_object* initialize_BioCompiler_BioCompiler_TypeSystem(uint8_t builtin) {
lean_object * res;
if (_G_initialized) return lean_io_result_mk_ok(lean_box(0));
_G_initialized = true;
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_BioCompiler_BioCompiler_ThreeValued(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_BioCompiler_BioCompiler_Sequence(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_BioCompiler_BioCompiler_NDFST(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_BioCompiler_BioCompiler_Scanners(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
