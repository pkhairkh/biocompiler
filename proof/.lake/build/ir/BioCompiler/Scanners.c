// Lean compiler output
// Module: BioCompiler.Scanners
// Imports: public import Init public meta import Init public import BioCompiler.Sequence
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
lean_object* l_Nat_cast___at___00Dyadic_toRat_spec__0(lean_object*);
lean_object* l_Rat_div(lean_object*, lean_object*);
lean_object* l_List_replicateTR___redArg(lean_object*, lean_object*);
lean_object* lean_nat_to_int(lean_object*);
lean_object* l_Nat_reprFast(lean_object*);
lean_object* lean_string_length(lean_object*);
uint8_t lean_nat_dec_eq(lean_object*, lean_object*);
lean_object* l_Int_repr(lean_object*);
lean_object* lean_string_append(lean_object*, lean_object*);
uint8_t lean_int_dec_lt(lean_object*, lean_object*);
lean_object* l_Repr_addAppParen(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_Sequence_matchesAt(lean_object*, lean_object*, lean_object*);
lean_object* l_Rat_ofScientific(lean_object*, uint8_t, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_Sequence_containsPattern(lean_object*, lean_object*);
lean_object* l_String_quote(lean_object*);
lean_object* lp_BioCompiler_BioCompiler_instBEqNucleotide_beq___boxed(lean_object*, lean_object*);
lean_object* l_List_lengthTR___redArg(lean_object*);
lean_object* lean_nat_mod(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(lean_object*, lean_object*);
lean_object* lean_nat_add(lean_object*, lean_object*);
lean_object* l_List_drop___redArg(lean_object*, lean_object*);
lean_object* lean_mk_empty_array_with_capacity(lean_object*);
lean_object* l___private_Init_Data_List_Impl_0__List_takeTR_go___redArg(lean_object*, lean_object*, lean_object*, lean_object*);
lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__0(lean_object*, lean_object*);
lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__1(lean_object*, lean_object*);
lean_object* l_Rat_add(lean_object*, lean_object*);
uint8_t l_Rat_instDecidableLe(lean_object*, lean_object*);
lean_object* l___private_Init_Data_List_Impl_0__List_zipWithTR_go___redArg(lean_object*, lean_object*, lean_object*, lean_object*);
lean_object* l_Rat_mul(lean_object*, lean_object*);
lean_object* l_Rat_neg(lean_object*);
uint8_t l_List_any___redArg(lean_object*, lean_object*);
lean_object* lean_nat_sub(lean_object*, lean_object*);
lean_object* l_List_range(lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasPattern(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasPattern___boxed(lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasAnyRestrictionSite___lam__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasAnyRestrictionSite___lam__0___boxed(lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasAnyRestrictionSite(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasAnyRestrictionSite___boxed(lean_object*, lean_object*);
static const lean_ctor_object lp_BioCompiler_BioCompiler_atttaMotif___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)(((size_t)(0) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_atttaMotif___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_atttaMotif___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(3) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_atttaMotif___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_atttaMotif___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(3) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_atttaMotif___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_atttaMotif___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(3) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__2_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_atttaMotif___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__3_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_atttaMotif___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__3_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_atttaMotif___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__4_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_atttaMotif = (const lean_object*)&lp_BioCompiler_BioCompiler_atttaMotif___closed__4_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_minURichLength;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_tRun(lean_object*);
static lean_once_cell_t lp_BioCompiler_BioCompiler_uRichMotif___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_uRichMotif___closed__0;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_uRichMotif;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasInstabilityMotif(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasInstabilityMotif___boxed(lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_matchesInstabilityAttta(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_matchesInstabilityAttta___boxed(lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_matchesInstabilityURich(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_matchesInstabilityURich___boxed(lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_matchesInstabilityMotif(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_matchesInstabilityMotif___boxed(lean_object*, lean_object*);
static const lean_ctor_object lp_BioCompiler_BioCompiler_spliceDonorConsensus___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(3) << 1) | 1)),((lean_object*)(((size_t)(0) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_spliceDonorConsensus___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_spliceDonorConsensus___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_spliceDonorConsensus___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(2) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_spliceDonorConsensus___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_spliceDonorConsensus___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_spliceDonorConsensus___closed__1_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_spliceDonorConsensus = (const lean_object*)&lp_BioCompiler_BioCompiler_spliceDonorConsensus___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_spliceAcceptorConsensus___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(2) << 1) | 1)),((lean_object*)(((size_t)(0) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_spliceAcceptorConsensus___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_spliceAcceptorConsensus___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_spliceAcceptorConsensus___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_spliceAcceptorConsensus___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_spliceAcceptorConsensus___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_spliceAcceptorConsensus___closed__1_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_spliceAcceptorConsensus = (const lean_object*)&lp_BioCompiler_BioCompiler_spliceAcceptorConsensus___closed__1_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "siteType"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__2_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 5, .m_capacity = 5, .m_length = 4, .m_data = " := "};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__3_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__3_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__5_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__7_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = ","};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__7 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__7_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__7_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__8 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__8_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__9_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "position"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__9 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__9_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__10_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__9_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__10 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__10_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__11_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 3, .m_capacity = 3, .m_length = 2, .m_data = "{ "};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__11 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__11_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__12_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 6, .m_capacity = 6, .m_length = 5, .m_data = "score"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__12 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__12_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__13_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__12_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__13 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__13_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__14_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__14;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__15_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 3, .m_capacity = 3, .m_length = 2, .m_data = " }"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__15 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__15_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__16_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__16;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__18_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__11_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__18 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__18_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__19_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__15_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__19 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__19_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__20_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "("};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__20 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__20_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__21_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = " : Rat)/"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__21 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__21_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch___closed__0_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_crypticThreshold___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_crypticThreshold___closed__0;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_crypticThreshold;
static lean_once_cell_t lp_BioCompiler_BioCompiler_uncertainLoThreshold___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_uncertainLoThreshold___closed__0;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_uncertainLoThreshold;
static const lean_ctor_object lp_BioCompiler_BioCompiler_cpgDinucleotide___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(1) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_spliceAcceptorConsensus___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_cpgDinucleotide___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_cpgDinucleotide___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_cpgDinucleotide = (const lean_object*)&lp_BioCompiler_BioCompiler_cpgDinucleotide___closed__0_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_cpgIslandWindowSize;
static lean_once_cell_t lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__0;
static lean_once_cell_t lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__1_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__1;
static lean_once_cell_t lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__2_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__2;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_cpgIslandGCThreshold;
static lean_once_cell_t lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__0;
static lean_once_cell_t lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1;
static lean_once_cell_t lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__2_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__2;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold;
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_hasCpGIslandConcrete_spec__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_hasCpGIslandConcrete_spec__0___boxed(lean_object*, lean_object*);
static const lean_array_object lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_array_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 246}, .m_size = 0, .m_capacity = 0, .m_data = {}};
static const lean_object* lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instBEqNucleotide_beq___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasCpGIslandConcrete(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___boxed(lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_isValidCodingSeq(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_isValidCodingSeq___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_blosum62Score(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_blosum62Score___boxed(lean_object*, lean_object*);
static const lean_string_object lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "organism"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__3_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprPromoterMatch___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPromoterMatch___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPromoterMatch___closed__0_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_promoterThreshold___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_promoterThreshold___closed__0;
static lean_once_cell_t lp_BioCompiler_BioCompiler_promoterThreshold___closed__1_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_promoterThreshold___closed__1;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_promoterThreshold;
static lean_once_cell_t lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__0;
static lean_once_cell_t lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__1_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__1;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_promoterUncertainThreshold;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__10_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__0_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__1_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 11, .m_capacity = 11, .m_length = 10, .m_data = "windowSize"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__2_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__3_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__4_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__4;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 14, .m_capacity = 14, .m_length = 13, .m_data = "hydroFraction"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__5_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__5_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__6 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__6_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__7_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__7;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprTMDomainMatch___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTMDomainMatch___closed__0_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__0;
static lean_once_cell_t lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__1_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__1;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_tmDomainThreshold;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 7, .m_capacity = 7, .m_length = 6, .m_data = "deltaG"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__1_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__2_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__2;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch___closed__0_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__0;
static lean_once_cell_t lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__1_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__1;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_mrnaStructureThreshold;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 14, .m_capacity = 14, .m_length = 13, .m_data = "codonPosition"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 15, .m_capacity = 15, .m_length = 14, .m_data = "disruptionType"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__5_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__6_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__6;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprFoldingDisruption___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFoldingDisruption___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasPattern(lean_object* v_seq_1_, lean_object* v_pattern_2_){
_start:
{
uint8_t v___x_3_; 
v___x_3_ = lp_BioCompiler_BioCompiler_Sequence_containsPattern(v_seq_1_, v_pattern_2_);
return v___x_3_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasPattern___boxed(lean_object* v_seq_4_, lean_object* v_pattern_5_){
_start:
{
uint8_t v_res_6_; lean_object* v_r_7_; 
v_res_6_ = lp_BioCompiler_BioCompiler_hasPattern(v_seq_4_, v_pattern_5_);
v_r_7_ = lean_box(v_res_6_);
return v_r_7_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasAnyRestrictionSite___lam__0(lean_object* v_seq_8_, lean_object* v_site_9_){
_start:
{
uint8_t v___x_10_; 
v___x_10_ = lp_BioCompiler_BioCompiler_Sequence_containsPattern(v_seq_8_, v_site_9_);
return v___x_10_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasAnyRestrictionSite___lam__0___boxed(lean_object* v_seq_11_, lean_object* v_site_12_){
_start:
{
uint8_t v_res_13_; lean_object* v_r_14_; 
v_res_13_ = lp_BioCompiler_BioCompiler_hasAnyRestrictionSite___lam__0(v_seq_11_, v_site_12_);
v_r_14_ = lean_box(v_res_13_);
return v_r_14_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasAnyRestrictionSite(lean_object* v_seq_15_, lean_object* v_enzymeSites_16_){
_start:
{
lean_object* v___f_17_; uint8_t v___x_18_; 
v___f_17_ = lean_alloc_closure((void*)(lp_BioCompiler_BioCompiler_hasAnyRestrictionSite___lam__0___boxed), 2, 1);
lean_closure_set(v___f_17_, 0, v_seq_15_);
v___x_18_ = l_List_any___redArg(v_enzymeSites_16_, v___f_17_);
return v___x_18_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasAnyRestrictionSite___boxed(lean_object* v_seq_19_, lean_object* v_enzymeSites_20_){
_start:
{
uint8_t v_res_21_; lean_object* v_r_22_; 
v_res_21_ = lp_BioCompiler_BioCompiler_hasAnyRestrictionSite(v_seq_19_, v_enzymeSites_20_);
v_r_22_ = lean_box(v_res_21_);
return v_r_22_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_minURichLength(void){
_start:
{
lean_object* v___x_44_; 
v___x_44_ = lean_unsigned_to_nat(6u);
return v___x_44_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_tRun(lean_object* v_n_45_){
_start:
{
uint8_t v___x_46_; lean_object* v___x_47_; lean_object* v___x_48_; 
v___x_46_ = 3;
v___x_47_ = lean_box(v___x_46_);
v___x_48_ = l_List_replicateTR___redArg(v_n_45_, v___x_47_);
return v___x_48_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_uRichMotif___closed__0(void){
_start:
{
lean_object* v___x_49_; lean_object* v___x_50_; 
v___x_49_ = lean_unsigned_to_nat(6u);
v___x_50_ = lp_BioCompiler_BioCompiler_tRun(v___x_49_);
return v___x_50_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_uRichMotif(void){
_start:
{
lean_object* v___x_51_; 
v___x_51_ = lean_obj_once(&lp_BioCompiler_BioCompiler_uRichMotif___closed__0, &lp_BioCompiler_BioCompiler_uRichMotif___closed__0_once, _init_lp_BioCompiler_BioCompiler_uRichMotif___closed__0);
return v___x_51_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasInstabilityMotif(lean_object* v_seq_52_){
_start:
{
lean_object* v___x_53_; uint8_t v___x_54_; 
v___x_53_ = ((lean_object*)(lp_BioCompiler_BioCompiler_atttaMotif));
lean_inc(v_seq_52_);
v___x_54_ = lp_BioCompiler_BioCompiler_Sequence_containsPattern(v_seq_52_, v___x_53_);
if (v___x_54_ == 0)
{
lean_object* v___x_55_; uint8_t v___x_56_; 
v___x_55_ = lp_BioCompiler_BioCompiler_uRichMotif;
v___x_56_ = lp_BioCompiler_BioCompiler_Sequence_containsPattern(v_seq_52_, v___x_55_);
return v___x_56_;
}
else
{
lean_dec(v_seq_52_);
return v___x_54_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasInstabilityMotif___boxed(lean_object* v_seq_57_){
_start:
{
uint8_t v_res_58_; lean_object* v_r_59_; 
v_res_58_ = lp_BioCompiler_BioCompiler_hasInstabilityMotif(v_seq_57_);
v_r_59_ = lean_box(v_res_58_);
return v_r_59_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_matchesInstabilityAttta(lean_object* v_seq_60_, lean_object* v_pos_61_){
_start:
{
lean_object* v___x_62_; uint8_t v___x_63_; 
v___x_62_ = ((lean_object*)(lp_BioCompiler_BioCompiler_atttaMotif));
v___x_63_ = lp_BioCompiler_BioCompiler_Sequence_matchesAt(v_seq_60_, v___x_62_, v_pos_61_);
return v___x_63_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_matchesInstabilityAttta___boxed(lean_object* v_seq_64_, lean_object* v_pos_65_){
_start:
{
uint8_t v_res_66_; lean_object* v_r_67_; 
v_res_66_ = lp_BioCompiler_BioCompiler_matchesInstabilityAttta(v_seq_64_, v_pos_65_);
lean_dec(v_seq_64_);
v_r_67_ = lean_box(v_res_66_);
return v_r_67_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_matchesInstabilityURich(lean_object* v_seq_68_, lean_object* v_pos_69_){
_start:
{
lean_object* v___x_70_; uint8_t v___x_71_; 
v___x_70_ = lp_BioCompiler_BioCompiler_uRichMotif;
v___x_71_ = lp_BioCompiler_BioCompiler_Sequence_matchesAt(v_seq_68_, v___x_70_, v_pos_69_);
return v___x_71_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_matchesInstabilityURich___boxed(lean_object* v_seq_72_, lean_object* v_pos_73_){
_start:
{
uint8_t v_res_74_; lean_object* v_r_75_; 
v_res_74_ = lp_BioCompiler_BioCompiler_matchesInstabilityURich(v_seq_72_, v_pos_73_);
lean_dec(v_seq_72_);
v_r_75_ = lean_box(v_res_74_);
return v_r_75_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_matchesInstabilityMotif(lean_object* v_seq_76_, lean_object* v_pos_77_){
_start:
{
uint8_t v___x_78_; 
lean_inc(v_pos_77_);
v___x_78_ = lp_BioCompiler_BioCompiler_matchesInstabilityAttta(v_seq_76_, v_pos_77_);
if (v___x_78_ == 0)
{
uint8_t v___x_79_; 
v___x_79_ = lp_BioCompiler_BioCompiler_matchesInstabilityURich(v_seq_76_, v_pos_77_);
return v___x_79_;
}
else
{
lean_dec(v_pos_77_);
return v___x_78_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_matchesInstabilityMotif___boxed(lean_object* v_seq_80_, lean_object* v_pos_81_){
_start:
{
uint8_t v_res_82_; lean_object* v_r_83_; 
v_res_82_ = lp_BioCompiler_BioCompiler_matchesInstabilityMotif(v_seq_80_, v_pos_81_);
lean_dec(v_seq_80_);
v_r_83_ = lean_box(v_res_82_);
return v_r_83_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6(void){
_start:
{
lean_object* v___x_114_; lean_object* v___x_115_; 
v___x_114_ = lean_unsigned_to_nat(12u);
v___x_115_ = lean_nat_to_int(v___x_114_);
return v___x_115_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__14(void){
_start:
{
lean_object* v___x_126_; lean_object* v___x_127_; 
v___x_126_ = lean_unsigned_to_nat(9u);
v___x_127_ = lean_nat_to_int(v___x_126_);
return v___x_127_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__16(void){
_start:
{
lean_object* v___x_129_; lean_object* v___x_130_; 
v___x_129_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__11));
v___x_130_ = lean_string_length(v___x_129_);
return v___x_130_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17(void){
_start:
{
lean_object* v___x_131_; lean_object* v___x_132_; 
v___x_131_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__16, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__16_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__16);
v___x_132_ = lean_nat_to_int(v___x_131_);
return v___x_132_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22(void){
_start:
{
lean_object* v___x_139_; lean_object* v___x_140_; 
v___x_139_ = lean_unsigned_to_nat(0u);
v___x_140_ = lean_nat_to_int(v___x_139_);
return v___x_140_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg(lean_object* v_x_141_){
_start:
{
lean_object* v_siteType_142_; lean_object* v_position_143_; lean_object* v_score_144_; lean_object* v___x_145_; lean_object* v___x_146_; lean_object* v_num_147_; lean_object* v_den_148_; lean_object* v___x_150_; uint8_t v_isShared_151_; uint8_t v_isSharedCheck_209_; 
v_siteType_142_ = lean_ctor_get(v_x_141_, 0);
lean_inc_ref(v_siteType_142_);
v_position_143_ = lean_ctor_get(v_x_141_, 1);
lean_inc(v_position_143_);
v_score_144_ = lean_ctor_get(v_x_141_, 2);
lean_inc_ref(v_score_144_);
lean_dec_ref(v_x_141_);
v___x_145_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4));
v___x_146_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__5));
v_num_147_ = lean_ctor_get(v_score_144_, 0);
v_den_148_ = lean_ctor_get(v_score_144_, 1);
v_isSharedCheck_209_ = !lean_is_exclusive(v_score_144_);
if (v_isSharedCheck_209_ == 0)
{
v___x_150_ = v_score_144_;
v_isShared_151_ = v_isSharedCheck_209_;
goto v_resetjp_149_;
}
else
{
lean_inc(v_den_148_);
lean_inc(v_num_147_);
lean_dec(v_score_144_);
v___x_150_ = lean_box(0);
v_isShared_151_ = v_isSharedCheck_209_;
goto v_resetjp_149_;
}
v_resetjp_149_:
{
lean_object* v___x_152_; lean_object* v___x_153_; lean_object* v___x_154_; lean_object* v___x_156_; 
v___x_152_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6);
v___x_153_ = l_String_quote(v_siteType_142_);
v___x_154_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_154_, 0, v___x_153_);
if (v_isShared_151_ == 0)
{
lean_ctor_set_tag(v___x_150_, 4);
lean_ctor_set(v___x_150_, 1, v___x_154_);
lean_ctor_set(v___x_150_, 0, v___x_152_);
v___x_156_ = v___x_150_;
goto v_reusejp_155_;
}
else
{
lean_object* v_reuseFailAlloc_208_; 
v_reuseFailAlloc_208_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_208_, 0, v___x_152_);
lean_ctor_set(v_reuseFailAlloc_208_, 1, v___x_154_);
v___x_156_ = v_reuseFailAlloc_208_;
goto v_reusejp_155_;
}
v_reusejp_155_:
{
uint8_t v___x_157_; lean_object* v___x_158_; lean_object* v___x_159_; lean_object* v___x_160_; lean_object* v___x_161_; lean_object* v___x_162_; lean_object* v___x_163_; lean_object* v___x_164_; lean_object* v___x_165_; lean_object* v___x_166_; lean_object* v___x_167_; lean_object* v___x_168_; lean_object* v___x_169_; lean_object* v___x_170_; lean_object* v___x_171_; lean_object* v___x_172_; lean_object* v___x_173_; lean_object* v___x_174_; lean_object* v___x_175_; lean_object* v___x_176_; lean_object* v___x_177_; lean_object* v___y_179_; lean_object* v___x_190_; uint8_t v___x_191_; 
v___x_157_ = 0;
v___x_158_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_158_, 0, v___x_156_);
lean_ctor_set_uint8(v___x_158_, sizeof(void*)*1, v___x_157_);
v___x_159_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_159_, 0, v___x_146_);
lean_ctor_set(v___x_159_, 1, v___x_158_);
v___x_160_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__8));
v___x_161_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_161_, 0, v___x_159_);
lean_ctor_set(v___x_161_, 1, v___x_160_);
v___x_162_ = lean_box(1);
v___x_163_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_163_, 0, v___x_161_);
lean_ctor_set(v___x_163_, 1, v___x_162_);
v___x_164_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__10));
v___x_165_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_165_, 0, v___x_163_);
lean_ctor_set(v___x_165_, 1, v___x_164_);
v___x_166_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_166_, 0, v___x_165_);
lean_ctor_set(v___x_166_, 1, v___x_145_);
v___x_167_ = l_Nat_reprFast(v_position_143_);
v___x_168_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_168_, 0, v___x_167_);
v___x_169_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_169_, 0, v___x_152_);
lean_ctor_set(v___x_169_, 1, v___x_168_);
v___x_170_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_170_, 0, v___x_169_);
lean_ctor_set_uint8(v___x_170_, sizeof(void*)*1, v___x_157_);
v___x_171_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_171_, 0, v___x_166_);
lean_ctor_set(v___x_171_, 1, v___x_170_);
v___x_172_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_172_, 0, v___x_171_);
lean_ctor_set(v___x_172_, 1, v___x_160_);
v___x_173_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_173_, 0, v___x_172_);
lean_ctor_set(v___x_173_, 1, v___x_162_);
v___x_174_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__13));
v___x_175_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_175_, 0, v___x_173_);
lean_ctor_set(v___x_175_, 1, v___x_174_);
v___x_176_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_176_, 0, v___x_175_);
lean_ctor_set(v___x_176_, 1, v___x_145_);
v___x_177_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__14, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__14_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__14);
v___x_190_ = lean_unsigned_to_nat(1u);
v___x_191_ = lean_nat_dec_eq(v_den_148_, v___x_190_);
if (v___x_191_ == 0)
{
lean_object* v___x_192_; lean_object* v___x_193_; lean_object* v___x_194_; lean_object* v___x_195_; lean_object* v___x_196_; lean_object* v___x_197_; lean_object* v___x_198_; lean_object* v___x_199_; 
v___x_192_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__20));
v___x_193_ = l_Int_repr(v_num_147_);
lean_dec(v_num_147_);
v___x_194_ = lean_string_append(v___x_192_, v___x_193_);
lean_dec_ref(v___x_193_);
v___x_195_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__21));
v___x_196_ = lean_string_append(v___x_194_, v___x_195_);
v___x_197_ = l_Nat_reprFast(v_den_148_);
v___x_198_ = lean_string_append(v___x_196_, v___x_197_);
lean_dec_ref(v___x_197_);
v___x_199_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_199_, 0, v___x_198_);
v___y_179_ = v___x_199_;
goto v___jp_178_;
}
else
{
lean_object* v___x_200_; lean_object* v___x_201_; uint8_t v___x_202_; 
lean_dec(v_den_148_);
v___x_200_ = lean_unsigned_to_nat(0u);
v___x_201_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22);
v___x_202_ = lean_int_dec_lt(v_num_147_, v___x_201_);
if (v___x_202_ == 0)
{
lean_object* v___x_203_; lean_object* v___x_204_; 
v___x_203_ = l_Int_repr(v_num_147_);
lean_dec(v_num_147_);
v___x_204_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_204_, 0, v___x_203_);
v___y_179_ = v___x_204_;
goto v___jp_178_;
}
else
{
lean_object* v___x_205_; lean_object* v___x_206_; lean_object* v___x_207_; 
v___x_205_ = l_Int_repr(v_num_147_);
lean_dec(v_num_147_);
v___x_206_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_206_, 0, v___x_205_);
v___x_207_ = l_Repr_addAppParen(v___x_206_, v___x_200_);
v___y_179_ = v___x_207_;
goto v___jp_178_;
}
}
v___jp_178_:
{
lean_object* v___x_180_; lean_object* v___x_181_; lean_object* v___x_182_; lean_object* v___x_183_; lean_object* v___x_184_; lean_object* v___x_185_; lean_object* v___x_186_; lean_object* v___x_187_; lean_object* v___x_188_; lean_object* v___x_189_; 
v___x_180_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_180_, 0, v___x_177_);
lean_ctor_set(v___x_180_, 1, v___y_179_);
v___x_181_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_181_, 0, v___x_180_);
lean_ctor_set_uint8(v___x_181_, sizeof(void*)*1, v___x_157_);
v___x_182_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_182_, 0, v___x_176_);
lean_ctor_set(v___x_182_, 1, v___x_181_);
v___x_183_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17);
v___x_184_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__18));
v___x_185_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_185_, 0, v___x_184_);
lean_ctor_set(v___x_185_, 1, v___x_182_);
v___x_186_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__19));
v___x_187_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_187_, 0, v___x_185_);
lean_ctor_set(v___x_187_, 1, v___x_186_);
v___x_188_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_188_, 0, v___x_183_);
lean_ctor_set(v___x_188_, 1, v___x_187_);
v___x_189_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_189_, 0, v___x_188_);
lean_ctor_set_uint8(v___x_189_, sizeof(void*)*1, v___x_157_);
return v___x_189_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr(lean_object* v_x_210_, lean_object* v_prec_211_){
_start:
{
lean_object* v___x_212_; 
v___x_212_ = lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg(v_x_210_);
return v___x_212_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___boxed(lean_object* v_x_213_, lean_object* v_prec_214_){
_start:
{
lean_object* v_res_215_; 
v_res_215_ = lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr(v_x_213_, v_prec_214_);
lean_dec(v_prec_214_);
return v_res_215_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_crypticThreshold___closed__0(void){
_start:
{
lean_object* v___x_218_; uint8_t v___x_219_; lean_object* v___x_220_; lean_object* v___x_221_; 
v___x_218_ = lean_unsigned_to_nat(1u);
v___x_219_ = 1;
v___x_220_ = lean_unsigned_to_nat(30u);
v___x_221_ = l_Rat_ofScientific(v___x_220_, v___x_219_, v___x_218_);
return v___x_221_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_crypticThreshold(void){
_start:
{
lean_object* v___x_222_; 
v___x_222_ = lean_obj_once(&lp_BioCompiler_BioCompiler_crypticThreshold___closed__0, &lp_BioCompiler_BioCompiler_crypticThreshold___closed__0_once, _init_lp_BioCompiler_BioCompiler_crypticThreshold___closed__0);
return v___x_222_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_uncertainLoThreshold___closed__0(void){
_start:
{
lean_object* v___x_223_; uint8_t v___x_224_; lean_object* v___x_225_; lean_object* v___x_226_; 
v___x_223_ = lean_unsigned_to_nat(1u);
v___x_224_ = 1;
v___x_225_ = lean_unsigned_to_nat(15u);
v___x_226_ = l_Rat_ofScientific(v___x_225_, v___x_224_, v___x_223_);
return v___x_226_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_uncertainLoThreshold(void){
_start:
{
lean_object* v___x_227_; 
v___x_227_ = lean_obj_once(&lp_BioCompiler_BioCompiler_uncertainLoThreshold___closed__0, &lp_BioCompiler_BioCompiler_uncertainLoThreshold___closed__0_once, _init_lp_BioCompiler_BioCompiler_uncertainLoThreshold___closed__0);
return v___x_227_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_cpgIslandWindowSize(void){
_start:
{
lean_object* v___x_233_; 
v___x_233_ = lean_unsigned_to_nat(200u);
return v___x_233_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__0(void){
_start:
{
lean_object* v___x_234_; lean_object* v___x_235_; 
v___x_234_ = lean_unsigned_to_nat(6u);
v___x_235_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_234_);
return v___x_235_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__1(void){
_start:
{
lean_object* v___x_236_; lean_object* v___x_237_; 
v___x_236_ = lean_unsigned_to_nat(10u);
v___x_237_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_236_);
return v___x_237_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__2(void){
_start:
{
lean_object* v___x_238_; lean_object* v___x_239_; lean_object* v___x_240_; 
v___x_238_ = lean_obj_once(&lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__1, &lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__1_once, _init_lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__1);
v___x_239_ = lean_obj_once(&lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__0, &lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__0_once, _init_lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__0);
v___x_240_ = l_Rat_div(v___x_239_, v___x_238_);
return v___x_240_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_cpgIslandGCThreshold(void){
_start:
{
lean_object* v___x_241_; 
v___x_241_ = lean_obj_once(&lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__2, &lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__2_once, _init_lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__2);
return v___x_241_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__0(void){
_start:
{
lean_object* v___x_242_; lean_object* v___x_243_; 
v___x_242_ = lean_unsigned_to_nat(65u);
v___x_243_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_242_);
return v___x_243_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1(void){
_start:
{
lean_object* v___x_244_; lean_object* v___x_245_; 
v___x_244_ = lean_unsigned_to_nat(100u);
v___x_245_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_244_);
return v___x_245_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__2(void){
_start:
{
lean_object* v___x_246_; lean_object* v___x_247_; lean_object* v___x_248_; 
v___x_246_ = lean_obj_once(&lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1, &lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1_once, _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1);
v___x_247_ = lean_obj_once(&lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__0, &lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__0_once, _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__0);
v___x_248_ = l_Rat_div(v___x_247_, v___x_246_);
return v___x_248_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold(void){
_start:
{
lean_object* v___x_249_; 
v___x_249_ = lean_obj_once(&lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__2, &lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__2_once, _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__2);
return v___x_249_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_hasCpGIslandConcrete_spec__0(lean_object* v_a_250_, lean_object* v_a_251_){
_start:
{
if (lean_obj_tag(v_a_250_) == 0)
{
return v_a_251_;
}
else
{
lean_object* v_head_252_; uint8_t v___x_253_; 
v_head_252_ = lean_ctor_get(v_a_250_, 0);
v___x_253_ = lean_unbox(v_head_252_);
if (v___x_253_ == 0)
{
lean_object* v_tail_254_; 
v_tail_254_ = lean_ctor_get(v_a_250_, 1);
v_a_250_ = v_tail_254_;
goto _start;
}
else
{
lean_object* v_tail_256_; lean_object* v___x_257_; lean_object* v___x_258_; 
v_tail_256_ = lean_ctor_get(v_a_250_, 1);
v___x_257_ = lean_unsigned_to_nat(1u);
v___x_258_ = lean_nat_add(v_a_251_, v___x_257_);
lean_dec(v_a_251_);
v_a_250_ = v_tail_256_;
v_a_251_ = v___x_258_;
goto _start;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_hasCpGIslandConcrete_spec__0___boxed(lean_object* v_a_260_, lean_object* v_a_261_){
_start:
{
lean_object* v_res_262_; 
v_res_262_ = lp_BioCompiler_List_countP_go___at___00BioCompiler_hasCpGIslandConcrete_spec__0(v_a_260_, v_a_261_);
lean_dec(v_a_260_);
return v_res_262_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0(lean_object* v_seq_265_, lean_object* v___x_266_, lean_object* v___x_267_, lean_object* v___f_268_, lean_object* v_pos_269_){
_start:
{
lean_object* v___x_270_; lean_object* v___x_271_; lean_object* v___x_272_; lean_object* v_window_273_; lean_object* v___x_274_; lean_object* v___x_275_; lean_object* v___x_276_; lean_object* v___x_277_; lean_object* v___x_278_; lean_object* v___x_279_; lean_object* v___x_280_; lean_object* v_gc_281_; lean_object* v___x_282_; uint8_t v___x_283_; 
v___x_270_ = l_List_drop___redArg(v_pos_269_, v_seq_265_);
v___x_271_ = lean_unsigned_to_nat(0u);
v___x_272_ = ((lean_object*)(lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0___closed__0));
lean_inc(v___x_270_);
v_window_273_ = l___private_Init_Data_List_Impl_0__List_takeTR_go___redArg(v___x_270_, v___x_270_, v___x_266_, v___x_272_);
lean_dec(v___x_270_);
v___x_274_ = lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__0(v_window_273_, v___x_271_);
v___x_275_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_274_);
v___x_276_ = lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__1(v_window_273_, v___x_271_);
v___x_277_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_276_);
lean_inc_ref(v___x_277_);
lean_inc_ref(v___x_275_);
v___x_278_ = l_Rat_add(v___x_275_, v___x_277_);
v___x_279_ = l_List_lengthTR___redArg(v_window_273_);
v___x_280_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_279_);
lean_inc_ref(v___x_280_);
v_gc_281_ = l_Rat_div(v___x_278_, v___x_280_);
lean_dec_ref(v___x_278_);
v___x_282_ = lp_BioCompiler_BioCompiler_cpgIslandGCThreshold;
v___x_283_ = l_Rat_instDecidableLe(v___x_282_, v_gc_281_);
if (v___x_283_ == 0)
{
lean_dec_ref(v___x_280_);
lean_dec_ref(v___x_277_);
lean_dec_ref(v___x_275_);
lean_dec(v_window_273_);
lean_dec_ref(v___f_268_);
lean_dec(v___x_267_);
return v___x_283_;
}
else
{
lean_object* v___x_284_; lean_object* v___x_285_; lean_object* v_cpgCount_286_; lean_object* v___x_287_; lean_object* v___x_288_; lean_object* v___x_289_; lean_object* v___x_290_; lean_object* v___x_291_; uint8_t v___x_292_; 
v___x_284_ = l_List_drop___redArg(v___x_267_, v_window_273_);
v___x_285_ = l___private_Init_Data_List_Impl_0__List_zipWithTR_go___redArg(v___f_268_, v_window_273_, v___x_284_, v___x_272_);
v_cpgCount_286_ = lp_BioCompiler_List_countP_go___at___00BioCompiler_hasCpGIslandConcrete_spec__0(v___x_285_, v___x_271_);
lean_dec(v___x_285_);
v___x_287_ = lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold;
v___x_288_ = l_Rat_mul(v___x_287_, v___x_277_);
v___x_289_ = l_Rat_mul(v___x_288_, v___x_275_);
lean_dec_ref(v___x_288_);
v___x_290_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v_cpgCount_286_);
v___x_291_ = l_Rat_mul(v___x_290_, v___x_280_);
lean_dec_ref(v___x_290_);
v___x_292_ = l_Rat_instDecidableLe(v___x_289_, v___x_291_);
return v___x_292_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0___boxed(lean_object* v_seq_293_, lean_object* v___x_294_, lean_object* v___x_295_, lean_object* v___f_296_, lean_object* v_pos_297_){
_start:
{
uint8_t v_res_298_; lean_object* v_r_299_; 
v_res_298_ = lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0(v_seq_293_, v___x_294_, v___x_295_, v___f_296_, v_pos_297_);
lean_dec(v_seq_293_);
v_r_299_ = lean_box(v_res_298_);
return v_r_299_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_hasCpGIslandConcrete(lean_object* v_seq_301_){
_start:
{
lean_object* v___f_302_; lean_object* v___x_303_; lean_object* v___x_304_; lean_object* v___x_305_; lean_object* v___x_306_; lean_object* v___f_307_; lean_object* v___x_308_; lean_object* v___x_309_; uint8_t v___x_310_; 
v___f_302_ = ((lean_object*)(lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___closed__0));
v___x_303_ = l_List_lengthTR___redArg(v_seq_301_);
v___x_304_ = lean_unsigned_to_nat(1u);
v___x_305_ = lean_nat_add(v___x_303_, v___x_304_);
lean_dec(v___x_303_);
v___x_306_ = lean_unsigned_to_nat(200u);
v___f_307_ = lean_alloc_closure((void*)(lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___lam__0___boxed), 5, 4);
lean_closure_set(v___f_307_, 0, v_seq_301_);
lean_closure_set(v___f_307_, 1, v___x_306_);
lean_closure_set(v___f_307_, 2, v___x_304_);
lean_closure_set(v___f_307_, 3, v___f_302_);
v___x_308_ = lean_nat_sub(v___x_305_, v___x_306_);
lean_dec(v___x_305_);
v___x_309_ = l_List_range(v___x_308_);
v___x_310_ = l_List_any___redArg(v___x_309_, v___f_307_);
return v___x_310_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_hasCpGIslandConcrete___boxed(lean_object* v_seq_311_){
_start:
{
uint8_t v_res_312_; lean_object* v_r_313_; 
v_res_312_ = lp_BioCompiler_BioCompiler_hasCpGIslandConcrete(v_seq_311_);
v_r_313_ = lean_box(v_res_312_);
return v_r_313_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_isValidCodingSeq(lean_object* v_seq_314_){
_start:
{
lean_object* v___x_315_; lean_object* v___x_316_; lean_object* v___x_317_; lean_object* v___x_318_; uint8_t v___x_319_; 
v___x_315_ = l_List_lengthTR___redArg(v_seq_314_);
v___x_316_ = lean_unsigned_to_nat(3u);
v___x_317_ = lean_nat_mod(v___x_315_, v___x_316_);
lean_dec(v___x_315_);
v___x_318_ = lean_unsigned_to_nat(0u);
v___x_319_ = lean_nat_dec_eq(v___x_317_, v___x_318_);
lean_dec(v___x_317_);
if (v___x_319_ == 0)
{
lean_dec(v_seq_314_);
return v___x_319_;
}
else
{
uint8_t v___x_320_; 
v___x_320_ = lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(v_seq_314_, v___x_318_);
if (v___x_320_ == 0)
{
return v___x_319_;
}
else
{
uint8_t v___x_321_; 
v___x_321_ = 0;
return v___x_321_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_isValidCodingSeq___boxed(lean_object* v_seq_322_){
_start:
{
uint8_t v_res_323_; lean_object* v_r_324_; 
v_res_323_ = lp_BioCompiler_BioCompiler_isValidCodingSeq(v_seq_322_);
v_r_324_ = lean_box(v_res_323_);
return v_r_324_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_blosum62Score(lean_object* v_aa1_325_, lean_object* v_aa2_326_){
_start:
{
lean_object* v___x_327_; 
v___x_327_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22);
return v___x_327_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_blosum62Score___boxed(lean_object* v_aa1_328_, lean_object* v_aa2_329_){
_start:
{
lean_object* v_res_330_; 
v_res_330_ = lp_BioCompiler_BioCompiler_blosum62Score(v_aa1_328_, v_aa2_329_);
lean_dec_ref(v_aa2_329_);
lean_dec_ref(v_aa1_328_);
return v_res_330_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg(lean_object* v_x_340_){
_start:
{
lean_object* v_organism_341_; lean_object* v_position_342_; lean_object* v_score_343_; lean_object* v___x_344_; lean_object* v___x_345_; lean_object* v_num_346_; lean_object* v_den_347_; lean_object* v___x_349_; uint8_t v_isShared_350_; uint8_t v_isSharedCheck_408_; 
v_organism_341_ = lean_ctor_get(v_x_340_, 0);
lean_inc_ref(v_organism_341_);
v_position_342_ = lean_ctor_get(v_x_340_, 1);
lean_inc(v_position_342_);
v_score_343_ = lean_ctor_get(v_x_340_, 2);
lean_inc_ref(v_score_343_);
lean_dec_ref(v_x_340_);
v___x_344_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4));
v___x_345_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg___closed__3));
v_num_346_ = lean_ctor_get(v_score_343_, 0);
v_den_347_ = lean_ctor_get(v_score_343_, 1);
v_isSharedCheck_408_ = !lean_is_exclusive(v_score_343_);
if (v_isSharedCheck_408_ == 0)
{
v___x_349_ = v_score_343_;
v_isShared_350_ = v_isSharedCheck_408_;
goto v_resetjp_348_;
}
else
{
lean_inc(v_den_347_);
lean_inc(v_num_346_);
lean_dec(v_score_343_);
v___x_349_ = lean_box(0);
v_isShared_350_ = v_isSharedCheck_408_;
goto v_resetjp_348_;
}
v_resetjp_348_:
{
lean_object* v___x_351_; lean_object* v___x_352_; lean_object* v___x_353_; lean_object* v___x_355_; 
v___x_351_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6);
v___x_352_ = l_String_quote(v_organism_341_);
v___x_353_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_353_, 0, v___x_352_);
if (v_isShared_350_ == 0)
{
lean_ctor_set_tag(v___x_349_, 4);
lean_ctor_set(v___x_349_, 1, v___x_353_);
lean_ctor_set(v___x_349_, 0, v___x_351_);
v___x_355_ = v___x_349_;
goto v_reusejp_354_;
}
else
{
lean_object* v_reuseFailAlloc_407_; 
v_reuseFailAlloc_407_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_407_, 0, v___x_351_);
lean_ctor_set(v_reuseFailAlloc_407_, 1, v___x_353_);
v___x_355_ = v_reuseFailAlloc_407_;
goto v_reusejp_354_;
}
v_reusejp_354_:
{
uint8_t v___x_356_; lean_object* v___x_357_; lean_object* v___x_358_; lean_object* v___x_359_; lean_object* v___x_360_; lean_object* v___x_361_; lean_object* v___x_362_; lean_object* v___x_363_; lean_object* v___x_364_; lean_object* v___x_365_; lean_object* v___x_366_; lean_object* v___x_367_; lean_object* v___x_368_; lean_object* v___x_369_; lean_object* v___x_370_; lean_object* v___x_371_; lean_object* v___x_372_; lean_object* v___x_373_; lean_object* v___x_374_; lean_object* v___x_375_; lean_object* v___x_376_; lean_object* v___y_378_; lean_object* v___x_389_; uint8_t v___x_390_; 
v___x_356_ = 0;
v___x_357_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_357_, 0, v___x_355_);
lean_ctor_set_uint8(v___x_357_, sizeof(void*)*1, v___x_356_);
v___x_358_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_358_, 0, v___x_345_);
lean_ctor_set(v___x_358_, 1, v___x_357_);
v___x_359_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__8));
v___x_360_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_360_, 0, v___x_358_);
lean_ctor_set(v___x_360_, 1, v___x_359_);
v___x_361_ = lean_box(1);
v___x_362_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_362_, 0, v___x_360_);
lean_ctor_set(v___x_362_, 1, v___x_361_);
v___x_363_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__10));
v___x_364_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_364_, 0, v___x_362_);
lean_ctor_set(v___x_364_, 1, v___x_363_);
v___x_365_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_365_, 0, v___x_364_);
lean_ctor_set(v___x_365_, 1, v___x_344_);
v___x_366_ = l_Nat_reprFast(v_position_342_);
v___x_367_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_367_, 0, v___x_366_);
v___x_368_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_368_, 0, v___x_351_);
lean_ctor_set(v___x_368_, 1, v___x_367_);
v___x_369_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_369_, 0, v___x_368_);
lean_ctor_set_uint8(v___x_369_, sizeof(void*)*1, v___x_356_);
v___x_370_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_370_, 0, v___x_365_);
lean_ctor_set(v___x_370_, 1, v___x_369_);
v___x_371_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_371_, 0, v___x_370_);
lean_ctor_set(v___x_371_, 1, v___x_359_);
v___x_372_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_372_, 0, v___x_371_);
lean_ctor_set(v___x_372_, 1, v___x_361_);
v___x_373_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__13));
v___x_374_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_374_, 0, v___x_372_);
lean_ctor_set(v___x_374_, 1, v___x_373_);
v___x_375_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_375_, 0, v___x_374_);
lean_ctor_set(v___x_375_, 1, v___x_344_);
v___x_376_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__14, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__14_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__14);
v___x_389_ = lean_unsigned_to_nat(1u);
v___x_390_ = lean_nat_dec_eq(v_den_347_, v___x_389_);
if (v___x_390_ == 0)
{
lean_object* v___x_391_; lean_object* v___x_392_; lean_object* v___x_393_; lean_object* v___x_394_; lean_object* v___x_395_; lean_object* v___x_396_; lean_object* v___x_397_; lean_object* v___x_398_; 
v___x_391_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__20));
v___x_392_ = l_Int_repr(v_num_346_);
lean_dec(v_num_346_);
v___x_393_ = lean_string_append(v___x_391_, v___x_392_);
lean_dec_ref(v___x_392_);
v___x_394_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__21));
v___x_395_ = lean_string_append(v___x_393_, v___x_394_);
v___x_396_ = l_Nat_reprFast(v_den_347_);
v___x_397_ = lean_string_append(v___x_395_, v___x_396_);
lean_dec_ref(v___x_396_);
v___x_398_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_398_, 0, v___x_397_);
v___y_378_ = v___x_398_;
goto v___jp_377_;
}
else
{
lean_object* v___x_399_; lean_object* v___x_400_; uint8_t v___x_401_; 
lean_dec(v_den_347_);
v___x_399_ = lean_unsigned_to_nat(0u);
v___x_400_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22);
v___x_401_ = lean_int_dec_lt(v_num_346_, v___x_400_);
if (v___x_401_ == 0)
{
lean_object* v___x_402_; lean_object* v___x_403_; 
v___x_402_ = l_Int_repr(v_num_346_);
lean_dec(v_num_346_);
v___x_403_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_403_, 0, v___x_402_);
v___y_378_ = v___x_403_;
goto v___jp_377_;
}
else
{
lean_object* v___x_404_; lean_object* v___x_405_; lean_object* v___x_406_; 
v___x_404_ = l_Int_repr(v_num_346_);
lean_dec(v_num_346_);
v___x_405_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_405_, 0, v___x_404_);
v___x_406_ = l_Repr_addAppParen(v___x_405_, v___x_399_);
v___y_378_ = v___x_406_;
goto v___jp_377_;
}
}
v___jp_377_:
{
lean_object* v___x_379_; lean_object* v___x_380_; lean_object* v___x_381_; lean_object* v___x_382_; lean_object* v___x_383_; lean_object* v___x_384_; lean_object* v___x_385_; lean_object* v___x_386_; lean_object* v___x_387_; lean_object* v___x_388_; 
v___x_379_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_379_, 0, v___x_376_);
lean_ctor_set(v___x_379_, 1, v___y_378_);
v___x_380_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_380_, 0, v___x_379_);
lean_ctor_set_uint8(v___x_380_, sizeof(void*)*1, v___x_356_);
v___x_381_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_381_, 0, v___x_375_);
lean_ctor_set(v___x_381_, 1, v___x_380_);
v___x_382_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17);
v___x_383_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__18));
v___x_384_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_384_, 0, v___x_383_);
lean_ctor_set(v___x_384_, 1, v___x_381_);
v___x_385_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__19));
v___x_386_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_386_, 0, v___x_384_);
lean_ctor_set(v___x_386_, 1, v___x_385_);
v___x_387_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_387_, 0, v___x_382_);
lean_ctor_set(v___x_387_, 1, v___x_386_);
v___x_388_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_388_, 0, v___x_387_);
lean_ctor_set_uint8(v___x_388_, sizeof(void*)*1, v___x_356_);
return v___x_388_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr(lean_object* v_x_409_, lean_object* v_prec_410_){
_start:
{
lean_object* v___x_411_; 
v___x_411_ = lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___redArg(v_x_409_);
return v___x_411_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr___boxed(lean_object* v_x_412_, lean_object* v_prec_413_){
_start:
{
lean_object* v_res_414_; 
v_res_414_ = lp_BioCompiler_BioCompiler_instReprPromoterMatch_repr(v_x_412_, v_prec_413_);
lean_dec(v_prec_413_);
return v_res_414_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_promoterThreshold___closed__0(void){
_start:
{
lean_object* v___x_417_; lean_object* v___x_418_; 
v___x_417_ = lean_unsigned_to_nat(7u);
v___x_418_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_417_);
return v___x_418_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_promoterThreshold___closed__1(void){
_start:
{
lean_object* v___x_419_; lean_object* v___x_420_; lean_object* v___x_421_; 
v___x_419_ = lean_obj_once(&lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__1, &lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__1_once, _init_lp_BioCompiler_BioCompiler_cpgIslandGCThreshold___closed__1);
v___x_420_ = lean_obj_once(&lp_BioCompiler_BioCompiler_promoterThreshold___closed__0, &lp_BioCompiler_BioCompiler_promoterThreshold___closed__0_once, _init_lp_BioCompiler_BioCompiler_promoterThreshold___closed__0);
v___x_421_ = l_Rat_div(v___x_420_, v___x_419_);
return v___x_421_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_promoterThreshold(void){
_start:
{
lean_object* v___x_422_; 
v___x_422_ = lean_obj_once(&lp_BioCompiler_BioCompiler_promoterThreshold___closed__1, &lp_BioCompiler_BioCompiler_promoterThreshold___closed__1_once, _init_lp_BioCompiler_BioCompiler_promoterThreshold___closed__1);
return v___x_422_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__0(void){
_start:
{
lean_object* v___x_423_; lean_object* v___x_424_; 
v___x_423_ = lean_unsigned_to_nat(56u);
v___x_424_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_423_);
return v___x_424_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__1(void){
_start:
{
lean_object* v___x_425_; lean_object* v___x_426_; lean_object* v___x_427_; 
v___x_425_ = lean_obj_once(&lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1, &lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1_once, _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1);
v___x_426_ = lean_obj_once(&lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__0, &lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__0_once, _init_lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__0);
v___x_427_ = l_Rat_div(v___x_426_, v___x_425_);
return v___x_427_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_promoterUncertainThreshold(void){
_start:
{
lean_object* v___x_428_; 
v___x_428_ = lean_obj_once(&lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__1, &lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__1_once, _init_lp_BioCompiler_BioCompiler_promoterUncertainThreshold___closed__1);
return v___x_428_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__4(void){
_start:
{
lean_object* v___x_438_; lean_object* v___x_439_; 
v___x_438_ = lean_unsigned_to_nat(14u);
v___x_439_ = lean_nat_to_int(v___x_438_);
return v___x_439_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__7(void){
_start:
{
lean_object* v___x_443_; lean_object* v___x_444_; 
v___x_443_ = lean_unsigned_to_nat(17u);
v___x_444_ = lean_nat_to_int(v___x_443_);
return v___x_444_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg(lean_object* v_x_445_){
_start:
{
lean_object* v_position_446_; lean_object* v_windowSize_447_; lean_object* v_hydroFraction_448_; lean_object* v___x_449_; lean_object* v___x_450_; lean_object* v_num_451_; lean_object* v_den_452_; lean_object* v___x_454_; uint8_t v_isShared_455_; uint8_t v_isSharedCheck_514_; 
v_position_446_ = lean_ctor_get(v_x_445_, 0);
lean_inc(v_position_446_);
v_windowSize_447_ = lean_ctor_get(v_x_445_, 1);
lean_inc(v_windowSize_447_);
v_hydroFraction_448_ = lean_ctor_get(v_x_445_, 2);
lean_inc_ref(v_hydroFraction_448_);
lean_dec_ref(v_x_445_);
v___x_449_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4));
v___x_450_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__1));
v_num_451_ = lean_ctor_get(v_hydroFraction_448_, 0);
v_den_452_ = lean_ctor_get(v_hydroFraction_448_, 1);
v_isSharedCheck_514_ = !lean_is_exclusive(v_hydroFraction_448_);
if (v_isSharedCheck_514_ == 0)
{
v___x_454_ = v_hydroFraction_448_;
v_isShared_455_ = v_isSharedCheck_514_;
goto v_resetjp_453_;
}
else
{
lean_inc(v_den_452_);
lean_inc(v_num_451_);
lean_dec(v_hydroFraction_448_);
v___x_454_ = lean_box(0);
v_isShared_455_ = v_isSharedCheck_514_;
goto v_resetjp_453_;
}
v_resetjp_453_:
{
lean_object* v___x_456_; lean_object* v___x_457_; lean_object* v___x_458_; lean_object* v___x_460_; 
v___x_456_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6);
v___x_457_ = l_Nat_reprFast(v_position_446_);
v___x_458_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_458_, 0, v___x_457_);
if (v_isShared_455_ == 0)
{
lean_ctor_set_tag(v___x_454_, 4);
lean_ctor_set(v___x_454_, 1, v___x_458_);
lean_ctor_set(v___x_454_, 0, v___x_456_);
v___x_460_ = v___x_454_;
goto v_reusejp_459_;
}
else
{
lean_object* v_reuseFailAlloc_513_; 
v_reuseFailAlloc_513_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_513_, 0, v___x_456_);
lean_ctor_set(v_reuseFailAlloc_513_, 1, v___x_458_);
v___x_460_ = v_reuseFailAlloc_513_;
goto v_reusejp_459_;
}
v_reusejp_459_:
{
uint8_t v___x_461_; lean_object* v___x_462_; lean_object* v___x_463_; lean_object* v___x_464_; lean_object* v___x_465_; lean_object* v___x_466_; lean_object* v___x_467_; lean_object* v___x_468_; lean_object* v___x_469_; lean_object* v___x_470_; lean_object* v___x_471_; lean_object* v___x_472_; lean_object* v___x_473_; lean_object* v___x_474_; lean_object* v___x_475_; lean_object* v___x_476_; lean_object* v___x_477_; lean_object* v___x_478_; lean_object* v___x_479_; lean_object* v___x_480_; lean_object* v___x_481_; lean_object* v___x_482_; lean_object* v___y_484_; lean_object* v___x_495_; uint8_t v___x_496_; 
v___x_461_ = 0;
v___x_462_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_462_, 0, v___x_460_);
lean_ctor_set_uint8(v___x_462_, sizeof(void*)*1, v___x_461_);
v___x_463_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_463_, 0, v___x_450_);
lean_ctor_set(v___x_463_, 1, v___x_462_);
v___x_464_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__8));
v___x_465_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_465_, 0, v___x_463_);
lean_ctor_set(v___x_465_, 1, v___x_464_);
v___x_466_ = lean_box(1);
v___x_467_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_467_, 0, v___x_465_);
lean_ctor_set(v___x_467_, 1, v___x_466_);
v___x_468_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__3));
v___x_469_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_469_, 0, v___x_467_);
lean_ctor_set(v___x_469_, 1, v___x_468_);
v___x_470_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_470_, 0, v___x_469_);
lean_ctor_set(v___x_470_, 1, v___x_449_);
v___x_471_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__4, &lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__4_once, _init_lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__4);
v___x_472_ = l_Nat_reprFast(v_windowSize_447_);
v___x_473_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_473_, 0, v___x_472_);
v___x_474_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_474_, 0, v___x_471_);
lean_ctor_set(v___x_474_, 1, v___x_473_);
v___x_475_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_475_, 0, v___x_474_);
lean_ctor_set_uint8(v___x_475_, sizeof(void*)*1, v___x_461_);
v___x_476_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_476_, 0, v___x_470_);
lean_ctor_set(v___x_476_, 1, v___x_475_);
v___x_477_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_477_, 0, v___x_476_);
lean_ctor_set(v___x_477_, 1, v___x_464_);
v___x_478_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_478_, 0, v___x_477_);
lean_ctor_set(v___x_478_, 1, v___x_466_);
v___x_479_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__6));
v___x_480_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_480_, 0, v___x_478_);
lean_ctor_set(v___x_480_, 1, v___x_479_);
v___x_481_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_481_, 0, v___x_480_);
lean_ctor_set(v___x_481_, 1, v___x_449_);
v___x_482_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__7);
v___x_495_ = lean_unsigned_to_nat(1u);
v___x_496_ = lean_nat_dec_eq(v_den_452_, v___x_495_);
if (v___x_496_ == 0)
{
lean_object* v___x_497_; lean_object* v___x_498_; lean_object* v___x_499_; lean_object* v___x_500_; lean_object* v___x_501_; lean_object* v___x_502_; lean_object* v___x_503_; lean_object* v___x_504_; 
v___x_497_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__20));
v___x_498_ = l_Int_repr(v_num_451_);
lean_dec(v_num_451_);
v___x_499_ = lean_string_append(v___x_497_, v___x_498_);
lean_dec_ref(v___x_498_);
v___x_500_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__21));
v___x_501_ = lean_string_append(v___x_499_, v___x_500_);
v___x_502_ = l_Nat_reprFast(v_den_452_);
v___x_503_ = lean_string_append(v___x_501_, v___x_502_);
lean_dec_ref(v___x_502_);
v___x_504_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_504_, 0, v___x_503_);
v___y_484_ = v___x_504_;
goto v___jp_483_;
}
else
{
lean_object* v___x_505_; lean_object* v___x_506_; uint8_t v___x_507_; 
lean_dec(v_den_452_);
v___x_505_ = lean_unsigned_to_nat(0u);
v___x_506_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22);
v___x_507_ = lean_int_dec_lt(v_num_451_, v___x_506_);
if (v___x_507_ == 0)
{
lean_object* v___x_508_; lean_object* v___x_509_; 
v___x_508_ = l_Int_repr(v_num_451_);
lean_dec(v_num_451_);
v___x_509_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_509_, 0, v___x_508_);
v___y_484_ = v___x_509_;
goto v___jp_483_;
}
else
{
lean_object* v___x_510_; lean_object* v___x_511_; lean_object* v___x_512_; 
v___x_510_ = l_Int_repr(v_num_451_);
lean_dec(v_num_451_);
v___x_511_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_511_, 0, v___x_510_);
v___x_512_ = l_Repr_addAppParen(v___x_511_, v___x_505_);
v___y_484_ = v___x_512_;
goto v___jp_483_;
}
}
v___jp_483_:
{
lean_object* v___x_485_; lean_object* v___x_486_; lean_object* v___x_487_; lean_object* v___x_488_; lean_object* v___x_489_; lean_object* v___x_490_; lean_object* v___x_491_; lean_object* v___x_492_; lean_object* v___x_493_; lean_object* v___x_494_; 
v___x_485_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_485_, 0, v___x_482_);
lean_ctor_set(v___x_485_, 1, v___y_484_);
v___x_486_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_486_, 0, v___x_485_);
lean_ctor_set_uint8(v___x_486_, sizeof(void*)*1, v___x_461_);
v___x_487_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_487_, 0, v___x_481_);
lean_ctor_set(v___x_487_, 1, v___x_486_);
v___x_488_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17);
v___x_489_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__18));
v___x_490_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_490_, 0, v___x_489_);
lean_ctor_set(v___x_490_, 1, v___x_487_);
v___x_491_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__19));
v___x_492_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_492_, 0, v___x_490_);
lean_ctor_set(v___x_492_, 1, v___x_491_);
v___x_493_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_493_, 0, v___x_488_);
lean_ctor_set(v___x_493_, 1, v___x_492_);
v___x_494_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_494_, 0, v___x_493_);
lean_ctor_set_uint8(v___x_494_, sizeof(void*)*1, v___x_461_);
return v___x_494_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr(lean_object* v_x_515_, lean_object* v_prec_516_){
_start:
{
lean_object* v___x_517_; 
v___x_517_ = lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg(v_x_515_);
return v___x_517_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___boxed(lean_object* v_x_518_, lean_object* v_prec_519_){
_start:
{
lean_object* v_res_520_; 
v_res_520_ = lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr(v_x_518_, v_prec_519_);
lean_dec(v_prec_519_);
return v_res_520_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__0(void){
_start:
{
lean_object* v___x_523_; lean_object* v___x_524_; 
v___x_523_ = lean_unsigned_to_nat(68u);
v___x_524_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_523_);
return v___x_524_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__1(void){
_start:
{
lean_object* v___x_525_; lean_object* v___x_526_; lean_object* v___x_527_; 
v___x_525_ = lean_obj_once(&lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1, &lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1_once, _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold___closed__1);
v___x_526_ = lean_obj_once(&lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__0, &lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__0_once, _init_lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__0);
v___x_527_ = l_Rat_div(v___x_526_, v___x_525_);
return v___x_527_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_tmDomainThreshold(void){
_start:
{
lean_object* v___x_528_; 
v___x_528_ = lean_obj_once(&lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__1, &lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__1_once, _init_lp_BioCompiler_BioCompiler_tmDomainThreshold___closed__1);
return v___x_528_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__2(void){
_start:
{
lean_object* v___x_532_; lean_object* v___x_533_; 
v___x_532_ = lean_unsigned_to_nat(10u);
v___x_533_ = lean_nat_to_int(v___x_532_);
return v___x_533_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg(lean_object* v_x_534_){
_start:
{
lean_object* v_deltaG_535_; lean_object* v_position_536_; lean_object* v___x_538_; uint8_t v_isShared_539_; uint8_t v_isSharedCheck_597_; 
v_deltaG_535_ = lean_ctor_get(v_x_534_, 1);
v_position_536_ = lean_ctor_get(v_x_534_, 0);
v_isSharedCheck_597_ = !lean_is_exclusive(v_x_534_);
if (v_isSharedCheck_597_ == 0)
{
v___x_538_ = v_x_534_;
v_isShared_539_ = v_isSharedCheck_597_;
goto v_resetjp_537_;
}
else
{
lean_inc(v_deltaG_535_);
lean_inc(v_position_536_);
lean_dec(v_x_534_);
v___x_538_ = lean_box(0);
v_isShared_539_ = v_isSharedCheck_597_;
goto v_resetjp_537_;
}
v_resetjp_537_:
{
lean_object* v_num_540_; lean_object* v_den_541_; lean_object* v___x_543_; uint8_t v_isShared_544_; uint8_t v_isSharedCheck_596_; 
v_num_540_ = lean_ctor_get(v_deltaG_535_, 0);
v_den_541_ = lean_ctor_get(v_deltaG_535_, 1);
v_isSharedCheck_596_ = !lean_is_exclusive(v_deltaG_535_);
if (v_isSharedCheck_596_ == 0)
{
v___x_543_ = v_deltaG_535_;
v_isShared_544_ = v_isSharedCheck_596_;
goto v_resetjp_542_;
}
else
{
lean_inc(v_den_541_);
lean_inc(v_num_540_);
lean_dec(v_deltaG_535_);
v___x_543_ = lean_box(0);
v_isShared_544_ = v_isSharedCheck_596_;
goto v_resetjp_542_;
}
v_resetjp_542_:
{
lean_object* v___x_545_; lean_object* v___x_546_; lean_object* v___x_547_; lean_object* v___x_548_; lean_object* v___x_549_; lean_object* v___x_551_; 
v___x_545_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4));
v___x_546_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__1));
v___x_547_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__6);
v___x_548_ = l_Nat_reprFast(v_position_536_);
v___x_549_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_549_, 0, v___x_548_);
if (v_isShared_544_ == 0)
{
lean_ctor_set_tag(v___x_543_, 4);
lean_ctor_set(v___x_543_, 1, v___x_549_);
lean_ctor_set(v___x_543_, 0, v___x_547_);
v___x_551_ = v___x_543_;
goto v_reusejp_550_;
}
else
{
lean_object* v_reuseFailAlloc_595_; 
v_reuseFailAlloc_595_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_595_, 0, v___x_547_);
lean_ctor_set(v_reuseFailAlloc_595_, 1, v___x_549_);
v___x_551_ = v_reuseFailAlloc_595_;
goto v_reusejp_550_;
}
v_reusejp_550_:
{
uint8_t v___x_552_; lean_object* v___x_553_; lean_object* v___x_555_; 
v___x_552_ = 0;
v___x_553_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_553_, 0, v___x_551_);
lean_ctor_set_uint8(v___x_553_, sizeof(void*)*1, v___x_552_);
if (v_isShared_539_ == 0)
{
lean_ctor_set_tag(v___x_538_, 5);
lean_ctor_set(v___x_538_, 1, v___x_553_);
lean_ctor_set(v___x_538_, 0, v___x_546_);
v___x_555_ = v___x_538_;
goto v_reusejp_554_;
}
else
{
lean_object* v_reuseFailAlloc_594_; 
v_reuseFailAlloc_594_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_594_, 0, v___x_546_);
lean_ctor_set(v_reuseFailAlloc_594_, 1, v___x_553_);
v___x_555_ = v_reuseFailAlloc_594_;
goto v_reusejp_554_;
}
v_reusejp_554_:
{
lean_object* v___x_556_; lean_object* v___x_557_; lean_object* v___x_558_; lean_object* v___x_559_; lean_object* v___x_560_; lean_object* v___x_561_; lean_object* v___x_562_; lean_object* v___x_563_; lean_object* v___y_565_; lean_object* v___x_576_; uint8_t v___x_577_; 
v___x_556_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__8));
v___x_557_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_557_, 0, v___x_555_);
lean_ctor_set(v___x_557_, 1, v___x_556_);
v___x_558_ = lean_box(1);
v___x_559_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_559_, 0, v___x_557_);
lean_ctor_set(v___x_559_, 1, v___x_558_);
v___x_560_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__1));
v___x_561_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_561_, 0, v___x_559_);
lean_ctor_set(v___x_561_, 1, v___x_560_);
v___x_562_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_562_, 0, v___x_561_);
lean_ctor_set(v___x_562_, 1, v___x_545_);
v___x_563_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__2, &lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__2_once, _init_lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg___closed__2);
v___x_576_ = lean_unsigned_to_nat(1u);
v___x_577_ = lean_nat_dec_eq(v_den_541_, v___x_576_);
if (v___x_577_ == 0)
{
lean_object* v___x_578_; lean_object* v___x_579_; lean_object* v___x_580_; lean_object* v___x_581_; lean_object* v___x_582_; lean_object* v___x_583_; lean_object* v___x_584_; lean_object* v___x_585_; 
v___x_578_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__20));
v___x_579_ = l_Int_repr(v_num_540_);
lean_dec(v_num_540_);
v___x_580_ = lean_string_append(v___x_578_, v___x_579_);
lean_dec_ref(v___x_579_);
v___x_581_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__21));
v___x_582_ = lean_string_append(v___x_580_, v___x_581_);
v___x_583_ = l_Nat_reprFast(v_den_541_);
v___x_584_ = lean_string_append(v___x_582_, v___x_583_);
lean_dec_ref(v___x_583_);
v___x_585_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_585_, 0, v___x_584_);
v___y_565_ = v___x_585_;
goto v___jp_564_;
}
else
{
lean_object* v___x_586_; lean_object* v___x_587_; uint8_t v___x_588_; 
lean_dec(v_den_541_);
v___x_586_ = lean_unsigned_to_nat(0u);
v___x_587_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__22);
v___x_588_ = lean_int_dec_lt(v_num_540_, v___x_587_);
if (v___x_588_ == 0)
{
lean_object* v___x_589_; lean_object* v___x_590_; 
v___x_589_ = l_Int_repr(v_num_540_);
lean_dec(v_num_540_);
v___x_590_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_590_, 0, v___x_589_);
v___y_565_ = v___x_590_;
goto v___jp_564_;
}
else
{
lean_object* v___x_591_; lean_object* v___x_592_; lean_object* v___x_593_; 
v___x_591_ = l_Int_repr(v_num_540_);
lean_dec(v_num_540_);
v___x_592_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_592_, 0, v___x_591_);
v___x_593_ = l_Repr_addAppParen(v___x_592_, v___x_586_);
v___y_565_ = v___x_593_;
goto v___jp_564_;
}
}
v___jp_564_:
{
lean_object* v___x_566_; lean_object* v___x_567_; lean_object* v___x_568_; lean_object* v___x_569_; lean_object* v___x_570_; lean_object* v___x_571_; lean_object* v___x_572_; lean_object* v___x_573_; lean_object* v___x_574_; lean_object* v___x_575_; 
v___x_566_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_566_, 0, v___x_563_);
lean_ctor_set(v___x_566_, 1, v___y_565_);
v___x_567_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_567_, 0, v___x_566_);
lean_ctor_set_uint8(v___x_567_, sizeof(void*)*1, v___x_552_);
v___x_568_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_568_, 0, v___x_562_);
lean_ctor_set(v___x_568_, 1, v___x_567_);
v___x_569_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17);
v___x_570_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__18));
v___x_571_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_571_, 0, v___x_570_);
lean_ctor_set(v___x_571_, 1, v___x_568_);
v___x_572_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__19));
v___x_573_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_573_, 0, v___x_571_);
lean_ctor_set(v___x_573_, 1, v___x_572_);
v___x_574_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_574_, 0, v___x_569_);
lean_ctor_set(v___x_574_, 1, v___x_573_);
v___x_575_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_575_, 0, v___x_574_);
lean_ctor_set_uint8(v___x_575_, sizeof(void*)*1, v___x_552_);
return v___x_575_;
}
}
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr(lean_object* v_x_598_, lean_object* v_prec_599_){
_start:
{
lean_object* v___x_600_; 
v___x_600_ = lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___redArg(v_x_598_);
return v___x_600_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr___boxed(lean_object* v_x_601_, lean_object* v_prec_602_){
_start:
{
lean_object* v_res_603_; 
v_res_603_ = lp_BioCompiler_BioCompiler_instReprStructureStabilityMatch_repr(v_x_601_, v_prec_602_);
lean_dec(v_prec_602_);
return v_res_603_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__0(void){
_start:
{
lean_object* v___x_606_; lean_object* v___x_607_; 
v___x_606_ = lean_unsigned_to_nat(15u);
v___x_607_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_606_);
return v___x_607_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__1(void){
_start:
{
lean_object* v___x_608_; lean_object* v___x_609_; 
v___x_608_ = lean_obj_once(&lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__0, &lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__0_once, _init_lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__0);
v___x_609_ = l_Rat_neg(v___x_608_);
return v___x_609_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_mrnaStructureThreshold(void){
_start:
{
lean_object* v___x_610_; 
v___x_610_ = lean_obj_once(&lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__1, &lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__1_once, _init_lp_BioCompiler_BioCompiler_mrnaStructureThreshold___closed__1);
return v___x_610_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__6(void){
_start:
{
lean_object* v___x_623_; lean_object* v___x_624_; 
v___x_623_ = lean_unsigned_to_nat(18u);
v___x_624_ = lean_nat_to_int(v___x_623_);
return v___x_624_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg(lean_object* v_x_625_){
_start:
{
lean_object* v_codonPosition_626_; lean_object* v_disruptionType_627_; lean_object* v___x_629_; uint8_t v_isShared_630_; uint8_t v_isSharedCheck_662_; 
v_codonPosition_626_ = lean_ctor_get(v_x_625_, 0);
v_disruptionType_627_ = lean_ctor_get(v_x_625_, 1);
v_isSharedCheck_662_ = !lean_is_exclusive(v_x_625_);
if (v_isSharedCheck_662_ == 0)
{
v___x_629_ = v_x_625_;
v_isShared_630_ = v_isSharedCheck_662_;
goto v_resetjp_628_;
}
else
{
lean_inc(v_disruptionType_627_);
lean_inc(v_codonPosition_626_);
lean_dec(v_x_625_);
v___x_629_ = lean_box(0);
v_isShared_630_ = v_isSharedCheck_662_;
goto v_resetjp_628_;
}
v_resetjp_628_:
{
lean_object* v___x_631_; lean_object* v___x_632_; lean_object* v___x_633_; lean_object* v___x_634_; lean_object* v___x_635_; lean_object* v___x_637_; 
v___x_631_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__4));
v___x_632_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__3));
v___x_633_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprTMDomainMatch_repr___redArg___closed__7);
v___x_634_ = l_Nat_reprFast(v_codonPosition_626_);
v___x_635_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_635_, 0, v___x_634_);
if (v_isShared_630_ == 0)
{
lean_ctor_set_tag(v___x_629_, 4);
lean_ctor_set(v___x_629_, 1, v___x_635_);
lean_ctor_set(v___x_629_, 0, v___x_633_);
v___x_637_ = v___x_629_;
goto v_reusejp_636_;
}
else
{
lean_object* v_reuseFailAlloc_661_; 
v_reuseFailAlloc_661_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_661_, 0, v___x_633_);
lean_ctor_set(v_reuseFailAlloc_661_, 1, v___x_635_);
v___x_637_ = v_reuseFailAlloc_661_;
goto v_reusejp_636_;
}
v_reusejp_636_:
{
uint8_t v___x_638_; lean_object* v___x_639_; lean_object* v___x_640_; lean_object* v___x_641_; lean_object* v___x_642_; lean_object* v___x_643_; lean_object* v___x_644_; lean_object* v___x_645_; lean_object* v___x_646_; lean_object* v___x_647_; lean_object* v___x_648_; lean_object* v___x_649_; lean_object* v___x_650_; lean_object* v___x_651_; lean_object* v___x_652_; lean_object* v___x_653_; lean_object* v___x_654_; lean_object* v___x_655_; lean_object* v___x_656_; lean_object* v___x_657_; lean_object* v___x_658_; lean_object* v___x_659_; lean_object* v___x_660_; 
v___x_638_ = 0;
v___x_639_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_639_, 0, v___x_637_);
lean_ctor_set_uint8(v___x_639_, sizeof(void*)*1, v___x_638_);
v___x_640_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_640_, 0, v___x_632_);
lean_ctor_set(v___x_640_, 1, v___x_639_);
v___x_641_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__8));
v___x_642_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_642_, 0, v___x_640_);
lean_ctor_set(v___x_642_, 1, v___x_641_);
v___x_643_ = lean_box(1);
v___x_644_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_644_, 0, v___x_642_);
lean_ctor_set(v___x_644_, 1, v___x_643_);
v___x_645_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__5));
v___x_646_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_646_, 0, v___x_644_);
lean_ctor_set(v___x_646_, 1, v___x_645_);
v___x_647_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_647_, 0, v___x_646_);
lean_ctor_set(v___x_647_, 1, v___x_631_);
v___x_648_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__6, &lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg___closed__6);
v___x_649_ = l_String_quote(v_disruptionType_627_);
v___x_650_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_650_, 0, v___x_649_);
v___x_651_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_651_, 0, v___x_648_);
lean_ctor_set(v___x_651_, 1, v___x_650_);
v___x_652_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_652_, 0, v___x_651_);
lean_ctor_set_uint8(v___x_652_, sizeof(void*)*1, v___x_638_);
v___x_653_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_653_, 0, v___x_647_);
lean_ctor_set(v___x_653_, 1, v___x_652_);
v___x_654_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17, &lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__17);
v___x_655_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__18));
v___x_656_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_656_, 0, v___x_655_);
lean_ctor_set(v___x_656_, 1, v___x_653_);
v___x_657_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceSiteMatch_repr___redArg___closed__19));
v___x_658_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_658_, 0, v___x_656_);
lean_ctor_set(v___x_658_, 1, v___x_657_);
v___x_659_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_659_, 0, v___x_654_);
lean_ctor_set(v___x_659_, 1, v___x_658_);
v___x_660_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_660_, 0, v___x_659_);
lean_ctor_set_uint8(v___x_660_, sizeof(void*)*1, v___x_638_);
return v___x_660_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr(lean_object* v_x_663_, lean_object* v_prec_664_){
_start:
{
lean_object* v___x_665_; 
v___x_665_ = lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___redArg(v_x_663_);
return v___x_665_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr___boxed(lean_object* v_x_666_, lean_object* v_prec_667_){
_start:
{
lean_object* v_res_668_; 
v_res_668_ = lp_BioCompiler_BioCompiler_instReprFoldingDisruption_repr(v_x_666_, v_prec_667_);
lean_dec(v_prec_667_);
return v_res_668_;
}
}
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Sequence(uint8_t builtin);
static bool _G_initialized = false;
LEAN_EXPORT lean_object* initialize_BioCompiler_BioCompiler_Scanners(uint8_t builtin) {
lean_object * res;
if (_G_initialized) return lean_io_result_mk_ok(lean_box(0));
_G_initialized = true;
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_BioCompiler_BioCompiler_Sequence(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
lp_BioCompiler_BioCompiler_minURichLength = _init_lp_BioCompiler_BioCompiler_minURichLength();
lean_mark_persistent(lp_BioCompiler_BioCompiler_minURichLength);
lp_BioCompiler_BioCompiler_uRichMotif = _init_lp_BioCompiler_BioCompiler_uRichMotif();
lean_mark_persistent(lp_BioCompiler_BioCompiler_uRichMotif);
lp_BioCompiler_BioCompiler_crypticThreshold = _init_lp_BioCompiler_BioCompiler_crypticThreshold();
lean_mark_persistent(lp_BioCompiler_BioCompiler_crypticThreshold);
lp_BioCompiler_BioCompiler_uncertainLoThreshold = _init_lp_BioCompiler_BioCompiler_uncertainLoThreshold();
lean_mark_persistent(lp_BioCompiler_BioCompiler_uncertainLoThreshold);
lp_BioCompiler_BioCompiler_cpgIslandWindowSize = _init_lp_BioCompiler_BioCompiler_cpgIslandWindowSize();
lean_mark_persistent(lp_BioCompiler_BioCompiler_cpgIslandWindowSize);
lp_BioCompiler_BioCompiler_cpgIslandGCThreshold = _init_lp_BioCompiler_BioCompiler_cpgIslandGCThreshold();
lean_mark_persistent(lp_BioCompiler_BioCompiler_cpgIslandGCThreshold);
lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold = _init_lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold();
lean_mark_persistent(lp_BioCompiler_BioCompiler_cpgIslandObsExpThreshold);
lp_BioCompiler_BioCompiler_promoterThreshold = _init_lp_BioCompiler_BioCompiler_promoterThreshold();
lean_mark_persistent(lp_BioCompiler_BioCompiler_promoterThreshold);
lp_BioCompiler_BioCompiler_promoterUncertainThreshold = _init_lp_BioCompiler_BioCompiler_promoterUncertainThreshold();
lean_mark_persistent(lp_BioCompiler_BioCompiler_promoterUncertainThreshold);
lp_BioCompiler_BioCompiler_tmDomainThreshold = _init_lp_BioCompiler_BioCompiler_tmDomainThreshold();
lean_mark_persistent(lp_BioCompiler_BioCompiler_tmDomainThreshold);
lp_BioCompiler_BioCompiler_mrnaStructureThreshold = _init_lp_BioCompiler_BioCompiler_mrnaStructureThreshold();
lean_mark_persistent(lp_BioCompiler_BioCompiler_mrnaStructureThreshold);
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
