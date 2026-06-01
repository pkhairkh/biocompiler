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
lean_object* l_List_replicateTR___redArg(lean_object*, lean_object*);
lean_object* l_Nat_cast___at___00Dyadic_toRat_spec__0(lean_object*);
lean_object* l_Rat_div(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_Sequence_matchesAt(lean_object*, lean_object*, lean_object*);
lean_object* l_Rat_ofScientific(lean_object*, uint8_t, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_Sequence_containsPattern(lean_object*, lean_object*);
lean_object* lean_nat_to_int(lean_object*);
lean_object* l_String_quote(lean_object*);
lean_object* l_Nat_reprFast(lean_object*);
lean_object* lean_string_length(lean_object*);
uint8_t lean_nat_dec_eq(lean_object*, lean_object*);
lean_object* l_Int_repr(lean_object*);
lean_object* lean_string_append(lean_object*, lean_object*);
uint8_t lean_int_dec_lt(lean_object*, lean_object*);
lean_object* l_Repr_addAppParen(lean_object*, lean_object*);
lean_object* lp_BioCompiler_BioCompiler_instBEqNucleotide_beq___boxed(lean_object*, lean_object*);
lean_object* lean_nat_add(lean_object*, lean_object*);
lean_object* l_List_drop___redArg(lean_object*, lean_object*);
lean_object* lean_mk_empty_array_with_capacity(lean_object*);
lean_object* l___private_Init_Data_List_Impl_0__List_takeTR_go___redArg(lean_object*, lean_object*, lean_object*, lean_object*);
lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__0(lean_object*, lean_object*);
lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__1(lean_object*, lean_object*);
lean_object* l_Rat_add(lean_object*, lean_object*);
lean_object* l_List_lengthTR___redArg(lean_object*);
uint8_t l_Rat_instDecidableLe(lean_object*, lean_object*);
lean_object* l___private_Init_Data_List_Impl_0__List_zipWithTR_go___redArg(lean_object*, lean_object*, lean_object*, lean_object*);
lean_object* l_Rat_mul(lean_object*, lean_object*);
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
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
