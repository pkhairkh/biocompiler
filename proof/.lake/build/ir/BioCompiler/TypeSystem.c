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
uint8_t lean_string_dec_eq(lean_object*, lean_object*);
lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg(lean_object*, lean_object*);
uint8_t l_Rat_instDecidableLe(lean_object*, lean_object*);
lean_object* lp_BioCompiler_BioCompiler_Sequence_gcContent(lean_object*);
uint8_t lp_BioCompiler_BioCompiler_hasAnyRestrictionSite(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_hasInstabilityMotif(lean_object*);
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
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 42, .m_capacity = 42, .m_length = 41, .m_data = "BioCompiler.TypePredicate.NoCrypticSplice"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__1_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 45, .m_capacity = 45, .m_length = 44, .m_data = "BioCompiler.TypePredicate.NoInstabilityMotif"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__2_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 38, .m_capacity = 38, .m_length = 37, .m_data = "BioCompiler.TypePredicate.NoCpGIsland"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__5_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 40, .m_capacity = 40, .m_length = 39, .m_data = "BioCompiler.TypePredicate.SpliceCorrect"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__6 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__6_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__7_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__6_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__7 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__7_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__7_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__8 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__8_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__11_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 39, .m_capacity = 39, .m_length = 38, .m_data = "BioCompiler.TypePredicate.CodonAdapted"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__11 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__11_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__12_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__11_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__12 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__12_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__13_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__12_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__13 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__13_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "("};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = " : Rat)/"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__17_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 36, .m_capacity = 36, .m_length = 35, .m_data = "BioCompiler.TypePredicate.GCInRange"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__17 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__17_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__18_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__17_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__18 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__18_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__19_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__18_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__19 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__19_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__20_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 44, .m_capacity = 44, .m_length = 43, .m_data = "BioCompiler.TypePredicate.NoRestrictionSite"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__20 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__20_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__21_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__20_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__21 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__21_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__22_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__21_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__22 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__22_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__23_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 34, .m_capacity = 34, .m_length = 33, .m_data = "BioCompiler.TypePredicate.InFrame"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__23 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__23_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__24_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__23_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__24 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__24_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__24_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprTypePredicate___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprTypePredicate___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluate___redArg(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluate___redArg___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluate(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluate___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter___redArg(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
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
default: 
{
lean_object* v___x_9_; 
v___x_9_ = lean_unsigned_to_nat(7u);
return v___x_9_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorIdx___boxed(lean_object* v_x_10_){
_start:
{
lean_object* v_res_11_; 
v_res_11_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorIdx(v_x_10_);
lean_dec(v_x_10_);
return v_res_11_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(lean_object* v_t_12_, lean_object* v_k_13_){
_start:
{
switch(lean_obj_tag(v_t_12_))
{
case 0:
{
lean_object* v_cellType_14_; lean_object* v___x_15_; 
v_cellType_14_ = lean_ctor_get(v_t_12_, 0);
lean_inc_ref(v_cellType_14_);
lean_dec_ref(v_t_12_);
v___x_15_ = lean_apply_1(v_k_13_, v_cellType_14_);
return v___x_15_;
}
case 2:
{
lean_object* v_organism_16_; lean_object* v_threshold_17_; lean_object* v___x_18_; 
v_organism_16_ = lean_ctor_get(v_t_12_, 0);
lean_inc_ref(v_organism_16_);
v_threshold_17_ = lean_ctor_get(v_t_12_, 1);
lean_inc_ref(v_threshold_17_);
lean_dec_ref(v_t_12_);
v___x_18_ = lean_apply_2(v_k_13_, v_organism_16_, v_threshold_17_);
return v___x_18_;
}
case 3:
{
lean_object* v_lo_19_; lean_object* v_hi_20_; lean_object* v___x_21_; 
v_lo_19_ = lean_ctor_get(v_t_12_, 0);
lean_inc_ref(v_lo_19_);
v_hi_20_ = lean_ctor_get(v_t_12_, 1);
lean_inc_ref(v_hi_20_);
lean_dec_ref(v_t_12_);
v___x_21_ = lean_apply_2(v_k_13_, v_lo_19_, v_hi_20_);
return v___x_21_;
}
case 4:
{
lean_object* v_enzymeSites_22_; lean_object* v___x_23_; 
v_enzymeSites_22_ = lean_ctor_get(v_t_12_, 0);
lean_inc(v_enzymeSites_22_);
lean_dec_ref(v_t_12_);
v___x_23_ = lean_apply_1(v_k_13_, v_enzymeSites_22_);
return v___x_23_;
}
case 5:
{
lean_object* v_readingFrame_24_; lean_object* v_exonBoundaries_25_; lean_object* v___x_26_; 
v_readingFrame_24_ = lean_ctor_get(v_t_12_, 0);
lean_inc(v_readingFrame_24_);
v_exonBoundaries_25_ = lean_ctor_get(v_t_12_, 1);
lean_inc(v_exonBoundaries_25_);
lean_dec_ref(v_t_12_);
v___x_26_ = lean_apply_2(v_k_13_, v_readingFrame_24_, v_exonBoundaries_25_);
return v___x_26_;
}
default: 
{
lean_dec(v_t_12_);
return v_k_13_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorElim(lean_object* v_motive_27_, lean_object* v_ctorIdx_28_, lean_object* v_t_29_, lean_object* v_h_30_, lean_object* v_k_31_){
_start:
{
lean_object* v___x_32_; 
v___x_32_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_29_, v_k_31_);
return v___x_32_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___boxed(lean_object* v_motive_33_, lean_object* v_ctorIdx_34_, lean_object* v_t_35_, lean_object* v_h_36_, lean_object* v_k_37_){
_start:
{
lean_object* v_res_38_; 
v_res_38_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim(v_motive_33_, v_ctorIdx_34_, v_t_35_, v_h_36_, v_k_37_);
lean_dec(v_ctorIdx_34_);
return v_res_38_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SpliceCorrect_elim___redArg(lean_object* v_t_39_, lean_object* v_SpliceCorrect_40_){
_start:
{
lean_object* v___x_41_; 
v___x_41_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_39_, v_SpliceCorrect_40_);
return v___x_41_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_SpliceCorrect_elim(lean_object* v_motive_42_, lean_object* v_t_43_, lean_object* v_h_44_, lean_object* v_SpliceCorrect_45_){
_start:
{
lean_object* v___x_46_; 
v___x_46_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_43_, v_SpliceCorrect_45_);
return v___x_46_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticSplice_elim___redArg(lean_object* v_t_47_, lean_object* v_NoCrypticSplice_48_){
_start:
{
lean_object* v___x_49_; 
v___x_49_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_47_, v_NoCrypticSplice_48_);
return v___x_49_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCrypticSplice_elim(lean_object* v_motive_50_, lean_object* v_t_51_, lean_object* v_h_52_, lean_object* v_NoCrypticSplice_53_){
_start:
{
lean_object* v___x_54_; 
v___x_54_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_51_, v_NoCrypticSplice_53_);
return v___x_54_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonAdapted_elim___redArg(lean_object* v_t_55_, lean_object* v_CodonAdapted_56_){
_start:
{
lean_object* v___x_57_; 
v___x_57_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_55_, v_CodonAdapted_56_);
return v___x_57_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_CodonAdapted_elim(lean_object* v_motive_58_, lean_object* v_t_59_, lean_object* v_h_60_, lean_object* v_CodonAdapted_61_){
_start:
{
lean_object* v___x_62_; 
v___x_62_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_59_, v_CodonAdapted_61_);
return v___x_62_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_GCInRange_elim___redArg(lean_object* v_t_63_, lean_object* v_GCInRange_64_){
_start:
{
lean_object* v___x_65_; 
v___x_65_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_63_, v_GCInRange_64_);
return v___x_65_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_GCInRange_elim(lean_object* v_motive_66_, lean_object* v_t_67_, lean_object* v_h_68_, lean_object* v_GCInRange_69_){
_start:
{
lean_object* v___x_70_; 
v___x_70_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_67_, v_GCInRange_69_);
return v___x_70_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoRestrictionSite_elim___redArg(lean_object* v_t_71_, lean_object* v_NoRestrictionSite_72_){
_start:
{
lean_object* v___x_73_; 
v___x_73_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_71_, v_NoRestrictionSite_72_);
return v___x_73_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoRestrictionSite_elim(lean_object* v_motive_74_, lean_object* v_t_75_, lean_object* v_h_76_, lean_object* v_NoRestrictionSite_77_){
_start:
{
lean_object* v___x_78_; 
v___x_78_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_75_, v_NoRestrictionSite_77_);
return v___x_78_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_InFrame_elim___redArg(lean_object* v_t_79_, lean_object* v_InFrame_80_){
_start:
{
lean_object* v___x_81_; 
v___x_81_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_79_, v_InFrame_80_);
return v___x_81_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_InFrame_elim(lean_object* v_motive_82_, lean_object* v_t_83_, lean_object* v_h_84_, lean_object* v_InFrame_85_){
_start:
{
lean_object* v___x_86_; 
v___x_86_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_83_, v_InFrame_85_);
return v___x_86_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoInstabilityMotif_elim___redArg(lean_object* v_t_87_, lean_object* v_NoInstabilityMotif_88_){
_start:
{
lean_object* v___x_89_; 
v___x_89_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_87_, v_NoInstabilityMotif_88_);
return v___x_89_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoInstabilityMotif_elim(lean_object* v_motive_90_, lean_object* v_t_91_, lean_object* v_h_92_, lean_object* v_NoInstabilityMotif_93_){
_start:
{
lean_object* v___x_94_; 
v___x_94_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_91_, v_NoInstabilityMotif_93_);
return v___x_94_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCpGIsland_elim___redArg(lean_object* v_t_95_, lean_object* v_NoCpGIsland_96_){
_start:
{
lean_object* v___x_97_; 
v___x_97_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_95_, v_NoCpGIsland_96_);
return v___x_97_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_TypePredicate_NoCpGIsland_elim(lean_object* v_motive_98_, lean_object* v_t_99_, lean_object* v_h_100_, lean_object* v_NoCpGIsland_101_){
_start:
{
lean_object* v___x_102_; 
v___x_102_ = lp_BioCompiler_BioCompiler_TypePredicate_ctorElim___redArg(v_t_99_, v_NoCpGIsland_101_);
return v___x_102_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1_spec__3(lean_object* v_x_103_, lean_object* v_x_104_, lean_object* v_x_105_){
_start:
{
if (lean_obj_tag(v_x_105_) == 0)
{
lean_dec(v_x_103_);
return v_x_104_;
}
else
{
lean_object* v_head_106_; lean_object* v_tail_107_; lean_object* v___x_109_; uint8_t v_isShared_110_; uint8_t v_isSharedCheck_117_; 
v_head_106_ = lean_ctor_get(v_x_105_, 0);
v_tail_107_ = lean_ctor_get(v_x_105_, 1);
v_isSharedCheck_117_ = !lean_is_exclusive(v_x_105_);
if (v_isSharedCheck_117_ == 0)
{
v___x_109_ = v_x_105_;
v_isShared_110_ = v_isSharedCheck_117_;
goto v_resetjp_108_;
}
else
{
lean_inc(v_tail_107_);
lean_inc(v_head_106_);
lean_dec(v_x_105_);
v___x_109_ = lean_box(0);
v_isShared_110_ = v_isSharedCheck_117_;
goto v_resetjp_108_;
}
v_resetjp_108_:
{
lean_object* v___x_112_; 
lean_inc(v_x_103_);
if (v_isShared_110_ == 0)
{
lean_ctor_set_tag(v___x_109_, 5);
lean_ctor_set(v___x_109_, 1, v_x_103_);
lean_ctor_set(v___x_109_, 0, v_x_104_);
v___x_112_ = v___x_109_;
goto v_reusejp_111_;
}
else
{
lean_object* v_reuseFailAlloc_116_; 
v_reuseFailAlloc_116_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_116_, 0, v_x_104_);
lean_ctor_set(v_reuseFailAlloc_116_, 1, v_x_103_);
v___x_112_ = v_reuseFailAlloc_116_;
goto v_reusejp_111_;
}
v_reusejp_111_:
{
lean_object* v___x_113_; lean_object* v___x_114_; 
v___x_113_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_head_106_);
v___x_114_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_114_, 0, v___x_112_);
lean_ctor_set(v___x_114_, 1, v___x_113_);
v_x_104_ = v___x_114_;
v_x_105_ = v_tail_107_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1(lean_object* v_x_118_, lean_object* v_x_119_, lean_object* v_x_120_){
_start:
{
if (lean_obj_tag(v_x_120_) == 0)
{
lean_dec(v_x_118_);
return v_x_119_;
}
else
{
lean_object* v_head_121_; lean_object* v_tail_122_; lean_object* v___x_124_; uint8_t v_isShared_125_; uint8_t v_isSharedCheck_132_; 
v_head_121_ = lean_ctor_get(v_x_120_, 0);
v_tail_122_ = lean_ctor_get(v_x_120_, 1);
v_isSharedCheck_132_ = !lean_is_exclusive(v_x_120_);
if (v_isSharedCheck_132_ == 0)
{
v___x_124_ = v_x_120_;
v_isShared_125_ = v_isSharedCheck_132_;
goto v_resetjp_123_;
}
else
{
lean_inc(v_tail_122_);
lean_inc(v_head_121_);
lean_dec(v_x_120_);
v___x_124_ = lean_box(0);
v_isShared_125_ = v_isSharedCheck_132_;
goto v_resetjp_123_;
}
v_resetjp_123_:
{
lean_object* v___x_127_; 
lean_inc(v_x_118_);
if (v_isShared_125_ == 0)
{
lean_ctor_set_tag(v___x_124_, 5);
lean_ctor_set(v___x_124_, 1, v_x_118_);
lean_ctor_set(v___x_124_, 0, v_x_119_);
v___x_127_ = v___x_124_;
goto v_reusejp_126_;
}
else
{
lean_object* v_reuseFailAlloc_131_; 
v_reuseFailAlloc_131_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_131_, 0, v_x_119_);
lean_ctor_set(v_reuseFailAlloc_131_, 1, v_x_118_);
v___x_127_ = v_reuseFailAlloc_131_;
goto v_reusejp_126_;
}
v_reusejp_126_:
{
lean_object* v___x_128_; lean_object* v___x_129_; lean_object* v___x_130_; 
v___x_128_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_head_121_);
v___x_129_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_129_, 0, v___x_127_);
lean_ctor_set(v___x_129_, 1, v___x_128_);
v___x_130_ = lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1_spec__3(v_x_118_, v___x_129_, v_tail_122_);
return v___x_130_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0(lean_object* v_x_133_, lean_object* v_x_134_){
_start:
{
if (lean_obj_tag(v_x_133_) == 0)
{
lean_object* v___x_135_; 
lean_dec(v_x_134_);
v___x_135_ = lean_box(0);
return v___x_135_;
}
else
{
lean_object* v_tail_136_; 
v_tail_136_ = lean_ctor_get(v_x_133_, 1);
if (lean_obj_tag(v_tail_136_) == 0)
{
lean_object* v_head_137_; lean_object* v___x_138_; 
lean_dec(v_x_134_);
v_head_137_ = lean_ctor_get(v_x_133_, 0);
lean_inc(v_head_137_);
lean_dec_ref(v_x_133_);
v___x_138_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_head_137_);
return v___x_138_;
}
else
{
lean_object* v_head_139_; lean_object* v___x_140_; lean_object* v___x_141_; 
lean_inc(v_tail_136_);
v_head_139_ = lean_ctor_get(v_x_133_, 0);
lean_inc(v_head_139_);
lean_dec_ref(v_x_133_);
v___x_140_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_head_139_);
v___x_141_ = lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0_spec__1(v_x_134_, v___x_140_, v_tail_136_);
return v___x_141_;
}
}
}
}
static lean_object* _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7(void){
_start:
{
lean_object* v___x_153_; lean_object* v___x_154_; 
v___x_153_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__2));
v___x_154_ = lean_string_length(v___x_153_);
return v___x_154_;
}
}
static lean_object* _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8(void){
_start:
{
lean_object* v___x_155_; lean_object* v___x_156_; 
v___x_155_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7, &lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7_once, _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__7);
v___x_156_ = lean_nat_to_int(v___x_155_);
return v___x_156_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg(lean_object* v_a_161_){
_start:
{
if (lean_obj_tag(v_a_161_) == 0)
{
lean_object* v___x_162_; 
v___x_162_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__1));
return v___x_162_;
}
else
{
lean_object* v___x_163_; lean_object* v___x_164_; lean_object* v___x_165_; lean_object* v___x_166_; lean_object* v___x_167_; lean_object* v___x_168_; lean_object* v___x_169_; lean_object* v___x_170_; uint8_t v___x_171_; lean_object* v___x_172_; 
v___x_163_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__5));
v___x_164_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0_spec__0(v_a_161_, v___x_163_);
v___x_165_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8, &lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8_once, _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8);
v___x_166_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__9));
v___x_167_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_167_, 0, v___x_166_);
lean_ctor_set(v___x_167_, 1, v___x_164_);
v___x_168_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__10));
v___x_169_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_169_, 0, v___x_167_);
lean_ctor_set(v___x_169_, 1, v___x_168_);
v___x_170_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_170_, 0, v___x_165_);
lean_ctor_set(v___x_170_, 1, v___x_169_);
v___x_171_ = 0;
v___x_172_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_172_, 0, v___x_170_);
lean_ctor_set_uint8(v___x_172_, sizeof(void*)*1, v___x_171_);
return v___x_172_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2___lam__0(lean_object* v___y_173_){
_start:
{
lean_object* v___x_174_; lean_object* v___x_175_; 
v___x_174_ = l_Nat_reprFast(v___y_173_);
v___x_175_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_175_, 0, v___x_174_);
return v___x_175_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4_spec__6(lean_object* v_x_176_, lean_object* v_x_177_, lean_object* v_x_178_){
_start:
{
if (lean_obj_tag(v_x_178_) == 0)
{
lean_dec(v_x_176_);
return v_x_177_;
}
else
{
lean_object* v_head_179_; lean_object* v_tail_180_; lean_object* v___x_182_; uint8_t v_isShared_183_; uint8_t v_isSharedCheck_191_; 
v_head_179_ = lean_ctor_get(v_x_178_, 0);
v_tail_180_ = lean_ctor_get(v_x_178_, 1);
v_isSharedCheck_191_ = !lean_is_exclusive(v_x_178_);
if (v_isSharedCheck_191_ == 0)
{
v___x_182_ = v_x_178_;
v_isShared_183_ = v_isSharedCheck_191_;
goto v_resetjp_181_;
}
else
{
lean_inc(v_tail_180_);
lean_inc(v_head_179_);
lean_dec(v_x_178_);
v___x_182_ = lean_box(0);
v_isShared_183_ = v_isSharedCheck_191_;
goto v_resetjp_181_;
}
v_resetjp_181_:
{
lean_object* v___x_185_; 
lean_inc(v_x_176_);
if (v_isShared_183_ == 0)
{
lean_ctor_set_tag(v___x_182_, 5);
lean_ctor_set(v___x_182_, 1, v_x_176_);
lean_ctor_set(v___x_182_, 0, v_x_177_);
v___x_185_ = v___x_182_;
goto v_reusejp_184_;
}
else
{
lean_object* v_reuseFailAlloc_190_; 
v_reuseFailAlloc_190_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_190_, 0, v_x_177_);
lean_ctor_set(v_reuseFailAlloc_190_, 1, v_x_176_);
v___x_185_ = v_reuseFailAlloc_190_;
goto v_reusejp_184_;
}
v_reusejp_184_:
{
lean_object* v___x_186_; lean_object* v___x_187_; lean_object* v___x_188_; 
v___x_186_ = l_Nat_reprFast(v_head_179_);
v___x_187_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_187_, 0, v___x_186_);
v___x_188_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_188_, 0, v___x_185_);
lean_ctor_set(v___x_188_, 1, v___x_187_);
v_x_177_ = v___x_188_;
v_x_178_ = v_tail_180_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4(lean_object* v_x_192_, lean_object* v_x_193_, lean_object* v_x_194_){
_start:
{
if (lean_obj_tag(v_x_194_) == 0)
{
lean_dec(v_x_192_);
return v_x_193_;
}
else
{
lean_object* v_head_195_; lean_object* v_tail_196_; lean_object* v___x_198_; uint8_t v_isShared_199_; uint8_t v_isSharedCheck_207_; 
v_head_195_ = lean_ctor_get(v_x_194_, 0);
v_tail_196_ = lean_ctor_get(v_x_194_, 1);
v_isSharedCheck_207_ = !lean_is_exclusive(v_x_194_);
if (v_isSharedCheck_207_ == 0)
{
v___x_198_ = v_x_194_;
v_isShared_199_ = v_isSharedCheck_207_;
goto v_resetjp_197_;
}
else
{
lean_inc(v_tail_196_);
lean_inc(v_head_195_);
lean_dec(v_x_194_);
v___x_198_ = lean_box(0);
v_isShared_199_ = v_isSharedCheck_207_;
goto v_resetjp_197_;
}
v_resetjp_197_:
{
lean_object* v___x_201_; 
lean_inc(v_x_192_);
if (v_isShared_199_ == 0)
{
lean_ctor_set_tag(v___x_198_, 5);
lean_ctor_set(v___x_198_, 1, v_x_192_);
lean_ctor_set(v___x_198_, 0, v_x_193_);
v___x_201_ = v___x_198_;
goto v_reusejp_200_;
}
else
{
lean_object* v_reuseFailAlloc_206_; 
v_reuseFailAlloc_206_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_206_, 0, v_x_193_);
lean_ctor_set(v_reuseFailAlloc_206_, 1, v_x_192_);
v___x_201_ = v_reuseFailAlloc_206_;
goto v_reusejp_200_;
}
v_reusejp_200_:
{
lean_object* v___x_202_; lean_object* v___x_203_; lean_object* v___x_204_; lean_object* v___x_205_; 
v___x_202_ = l_Nat_reprFast(v_head_195_);
v___x_203_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_203_, 0, v___x_202_);
v___x_204_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_204_, 0, v___x_201_);
lean_ctor_set(v___x_204_, 1, v___x_203_);
v___x_205_ = lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4_spec__6(v_x_192_, v___x_204_, v_tail_196_);
return v___x_205_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2(lean_object* v_x_208_, lean_object* v_x_209_){
_start:
{
if (lean_obj_tag(v_x_208_) == 0)
{
lean_object* v___x_210_; 
lean_dec(v_x_209_);
v___x_210_ = lean_box(0);
return v___x_210_;
}
else
{
lean_object* v_tail_211_; 
v_tail_211_ = lean_ctor_get(v_x_208_, 1);
if (lean_obj_tag(v_tail_211_) == 0)
{
lean_object* v_head_212_; lean_object* v___x_213_; 
lean_dec(v_x_209_);
v_head_212_ = lean_ctor_get(v_x_208_, 0);
lean_inc(v_head_212_);
lean_dec_ref(v_x_208_);
v___x_213_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2___lam__0(v_head_212_);
return v___x_213_;
}
else
{
lean_object* v_head_214_; lean_object* v___x_215_; lean_object* v___x_216_; 
lean_inc(v_tail_211_);
v_head_214_ = lean_ctor_get(v_x_208_, 0);
lean_inc(v_head_214_);
lean_dec_ref(v_x_208_);
v___x_215_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2___lam__0(v_head_214_);
v___x_216_ = lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2_spec__4(v_x_209_, v___x_215_, v_tail_211_);
return v___x_216_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___redArg(lean_object* v_a_217_){
_start:
{
if (lean_obj_tag(v_a_217_) == 0)
{
lean_object* v___x_218_; 
v___x_218_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__1));
return v___x_218_;
}
else
{
lean_object* v___x_219_; lean_object* v___x_220_; lean_object* v___x_221_; lean_object* v___x_222_; lean_object* v___x_223_; lean_object* v___x_224_; lean_object* v___x_225_; lean_object* v___x_226_; lean_object* v___x_227_; 
v___x_219_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__5));
v___x_220_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1_spec__2(v_a_217_, v___x_219_);
v___x_221_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8, &lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8_once, _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__8);
v___x_222_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__9));
v___x_223_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_223_, 0, v___x_222_);
lean_ctor_set(v___x_223_, 1, v___x_220_);
v___x_224_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg___closed__10));
v___x_225_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_225_, 0, v___x_223_);
lean_ctor_set(v___x_225_, 1, v___x_224_);
v___x_226_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_226_, 0, v___x_221_);
lean_ctor_set(v___x_226_, 1, v___x_225_);
v___x_227_ = l_Std_Format_fill(v___x_226_);
return v___x_227_;
}
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9(void){
_start:
{
lean_object* v___x_243_; lean_object* v___x_244_; 
v___x_243_ = lean_unsigned_to_nat(2u);
v___x_244_ = lean_nat_to_int(v___x_243_);
return v___x_244_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10(void){
_start:
{
lean_object* v___x_245_; lean_object* v___x_246_; 
v___x_245_ = lean_unsigned_to_nat(1u);
v___x_246_ = lean_nat_to_int(v___x_245_);
return v___x_246_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16(void){
_start:
{
lean_object* v___x_255_; lean_object* v___x_256_; 
v___x_255_ = lean_unsigned_to_nat(0u);
v___x_256_ = lean_nat_to_int(v___x_255_);
return v___x_256_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr(lean_object* v_x_275_, lean_object* v_prec_276_){
_start:
{
lean_object* v___y_278_; lean_object* v___y_279_; lean_object* v___y_280_; lean_object* v___y_287_; lean_object* v___y_294_; lean_object* v___y_295_; lean_object* v___y_296_; lean_object* v___y_303_; lean_object* v___y_310_; 
switch(lean_obj_tag(v_x_275_))
{
case 0:
{
lean_object* v_cellType_316_; lean_object* v___x_318_; uint8_t v_isShared_319_; uint8_t v_isSharedCheck_336_; 
v_cellType_316_ = lean_ctor_get(v_x_275_, 0);
v_isSharedCheck_336_ = !lean_is_exclusive(v_x_275_);
if (v_isSharedCheck_336_ == 0)
{
v___x_318_ = v_x_275_;
v_isShared_319_ = v_isSharedCheck_336_;
goto v_resetjp_317_;
}
else
{
lean_inc(v_cellType_316_);
lean_dec(v_x_275_);
v___x_318_ = lean_box(0);
v_isShared_319_ = v_isSharedCheck_336_;
goto v_resetjp_317_;
}
v_resetjp_317_:
{
lean_object* v___y_321_; lean_object* v___x_332_; uint8_t v___x_333_; 
v___x_332_ = lean_unsigned_to_nat(1024u);
v___x_333_ = lean_nat_dec_le(v___x_332_, v_prec_276_);
if (v___x_333_ == 0)
{
lean_object* v___x_334_; 
v___x_334_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9);
v___y_321_ = v___x_334_;
goto v___jp_320_;
}
else
{
lean_object* v___x_335_; 
v___x_335_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10);
v___y_321_ = v___x_335_;
goto v___jp_320_;
}
v___jp_320_:
{
lean_object* v___x_322_; lean_object* v___x_323_; lean_object* v___x_325_; 
v___x_322_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__8));
v___x_323_ = l_String_quote(v_cellType_316_);
if (v_isShared_319_ == 0)
{
lean_ctor_set_tag(v___x_318_, 3);
lean_ctor_set(v___x_318_, 0, v___x_323_);
v___x_325_ = v___x_318_;
goto v_reusejp_324_;
}
else
{
lean_object* v_reuseFailAlloc_331_; 
v_reuseFailAlloc_331_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_331_, 0, v___x_323_);
v___x_325_ = v_reuseFailAlloc_331_;
goto v_reusejp_324_;
}
v_reusejp_324_:
{
lean_object* v___x_326_; lean_object* v___x_327_; uint8_t v___x_328_; lean_object* v___x_329_; lean_object* v___x_330_; 
v___x_326_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_326_, 0, v___x_322_);
lean_ctor_set(v___x_326_, 1, v___x_325_);
lean_inc(v___y_321_);
v___x_327_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_327_, 0, v___y_321_);
lean_ctor_set(v___x_327_, 1, v___x_326_);
v___x_328_ = 0;
v___x_329_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_329_, 0, v___x_327_);
lean_ctor_set_uint8(v___x_329_, sizeof(void*)*1, v___x_328_);
v___x_330_ = l_Repr_addAppParen(v___x_329_, v_prec_276_);
return v___x_330_;
}
}
}
}
case 1:
{
lean_object* v___x_337_; uint8_t v___x_338_; 
v___x_337_ = lean_unsigned_to_nat(1024u);
v___x_338_ = lean_nat_dec_le(v___x_337_, v_prec_276_);
if (v___x_338_ == 0)
{
lean_object* v___x_339_; 
v___x_339_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9);
v___y_287_ = v___x_339_;
goto v___jp_286_;
}
else
{
lean_object* v___x_340_; 
v___x_340_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10);
v___y_287_ = v___x_340_;
goto v___jp_286_;
}
}
case 2:
{
lean_object* v_organism_341_; lean_object* v_threshold_342_; lean_object* v___x_344_; uint8_t v_isShared_345_; uint8_t v_isSharedCheck_386_; 
v_organism_341_ = lean_ctor_get(v_x_275_, 0);
v_threshold_342_ = lean_ctor_get(v_x_275_, 1);
v_isSharedCheck_386_ = !lean_is_exclusive(v_x_275_);
if (v_isSharedCheck_386_ == 0)
{
v___x_344_ = v_x_275_;
v_isShared_345_ = v_isSharedCheck_386_;
goto v_resetjp_343_;
}
else
{
lean_inc(v_threshold_342_);
lean_inc(v_organism_341_);
lean_dec(v_x_275_);
v___x_344_ = lean_box(0);
v_isShared_345_ = v_isSharedCheck_386_;
goto v_resetjp_343_;
}
v_resetjp_343_:
{
lean_object* v___y_347_; lean_object* v___x_382_; uint8_t v___x_383_; 
v___x_382_ = lean_unsigned_to_nat(1024u);
v___x_383_ = lean_nat_dec_le(v___x_382_, v_prec_276_);
if (v___x_383_ == 0)
{
lean_object* v___x_384_; 
v___x_384_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9);
v___y_347_ = v___x_384_;
goto v___jp_346_;
}
else
{
lean_object* v___x_385_; 
v___x_385_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10);
v___y_347_ = v___x_385_;
goto v___jp_346_;
}
v___jp_346_:
{
lean_object* v_num_348_; lean_object* v_den_349_; lean_object* v___x_351_; uint8_t v_isShared_352_; uint8_t v_isSharedCheck_381_; 
v_num_348_ = lean_ctor_get(v_threshold_342_, 0);
v_den_349_ = lean_ctor_get(v_threshold_342_, 1);
v_isSharedCheck_381_ = !lean_is_exclusive(v_threshold_342_);
if (v_isSharedCheck_381_ == 0)
{
v___x_351_ = v_threshold_342_;
v_isShared_352_ = v_isSharedCheck_381_;
goto v_resetjp_350_;
}
else
{
lean_inc(v_den_349_);
lean_inc(v_num_348_);
lean_dec(v_threshold_342_);
v___x_351_ = lean_box(0);
v_isShared_352_ = v_isSharedCheck_381_;
goto v_resetjp_350_;
}
v_resetjp_350_:
{
lean_object* v___x_353_; lean_object* v___x_354_; lean_object* v___x_355_; lean_object* v___x_356_; lean_object* v___x_358_; 
v___x_353_ = lean_box(1);
v___x_354_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__13));
v___x_355_ = l_String_quote(v_organism_341_);
v___x_356_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_356_, 0, v___x_355_);
if (v_isShared_352_ == 0)
{
lean_ctor_set_tag(v___x_351_, 5);
lean_ctor_set(v___x_351_, 1, v___x_356_);
lean_ctor_set(v___x_351_, 0, v___x_354_);
v___x_358_ = v___x_351_;
goto v_reusejp_357_;
}
else
{
lean_object* v_reuseFailAlloc_380_; 
v_reuseFailAlloc_380_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_380_, 0, v___x_354_);
lean_ctor_set(v_reuseFailAlloc_380_, 1, v___x_356_);
v___x_358_ = v_reuseFailAlloc_380_;
goto v_reusejp_357_;
}
v_reusejp_357_:
{
lean_object* v___x_360_; 
if (v_isShared_345_ == 0)
{
lean_ctor_set_tag(v___x_344_, 5);
lean_ctor_set(v___x_344_, 1, v___x_353_);
lean_ctor_set(v___x_344_, 0, v___x_358_);
v___x_360_ = v___x_344_;
goto v_reusejp_359_;
}
else
{
lean_object* v_reuseFailAlloc_379_; 
v_reuseFailAlloc_379_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_379_, 0, v___x_358_);
lean_ctor_set(v_reuseFailAlloc_379_, 1, v___x_353_);
v___x_360_ = v_reuseFailAlloc_379_;
goto v_reusejp_359_;
}
v_reusejp_359_:
{
lean_object* v___x_361_; uint8_t v___x_362_; 
v___x_361_ = lean_unsigned_to_nat(1u);
v___x_362_ = lean_nat_dec_eq(v_den_349_, v___x_361_);
if (v___x_362_ == 0)
{
lean_object* v___x_363_; lean_object* v___x_364_; lean_object* v___x_365_; lean_object* v___x_366_; lean_object* v___x_367_; lean_object* v___x_368_; lean_object* v___x_369_; lean_object* v___x_370_; 
v___x_363_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14));
v___x_364_ = l_Int_repr(v_num_348_);
lean_dec(v_num_348_);
v___x_365_ = lean_string_append(v___x_363_, v___x_364_);
lean_dec_ref(v___x_364_);
v___x_366_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15));
v___x_367_ = lean_string_append(v___x_365_, v___x_366_);
v___x_368_ = l_Nat_reprFast(v_den_349_);
v___x_369_ = lean_string_append(v___x_367_, v___x_368_);
lean_dec_ref(v___x_368_);
v___x_370_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_370_, 0, v___x_369_);
v___y_278_ = v___x_360_;
v___y_279_ = v___y_347_;
v___y_280_ = v___x_370_;
goto v___jp_277_;
}
else
{
lean_object* v___x_371_; lean_object* v___x_372_; uint8_t v___x_373_; 
lean_dec(v_den_349_);
v___x_371_ = lean_unsigned_to_nat(0u);
v___x_372_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16);
v___x_373_ = lean_int_dec_lt(v_num_348_, v___x_372_);
if (v___x_373_ == 0)
{
lean_object* v___x_374_; lean_object* v___x_375_; 
v___x_374_ = l_Int_repr(v_num_348_);
lean_dec(v_num_348_);
v___x_375_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_375_, 0, v___x_374_);
v___y_278_ = v___x_360_;
v___y_279_ = v___y_347_;
v___y_280_ = v___x_375_;
goto v___jp_277_;
}
else
{
lean_object* v___x_376_; lean_object* v___x_377_; lean_object* v___x_378_; 
v___x_376_ = l_Int_repr(v_num_348_);
lean_dec(v_num_348_);
v___x_377_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_377_, 0, v___x_376_);
v___x_378_ = l_Repr_addAppParen(v___x_377_, v___x_371_);
v___y_278_ = v___x_360_;
v___y_279_ = v___y_347_;
v___y_280_ = v___x_378_;
goto v___jp_277_;
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
lean_object* v_lo_387_; lean_object* v_hi_388_; lean_object* v___x_390_; uint8_t v_isShared_391_; uint8_t v_isSharedCheck_455_; 
v_lo_387_ = lean_ctor_get(v_x_275_, 0);
v_hi_388_ = lean_ctor_get(v_x_275_, 1);
v_isSharedCheck_455_ = !lean_is_exclusive(v_x_275_);
if (v_isSharedCheck_455_ == 0)
{
v___x_390_ = v_x_275_;
v_isShared_391_ = v_isSharedCheck_455_;
goto v_resetjp_389_;
}
else
{
lean_inc(v_hi_388_);
lean_inc(v_lo_387_);
lean_dec(v_x_275_);
v___x_390_ = lean_box(0);
v_isShared_391_ = v_isSharedCheck_455_;
goto v_resetjp_389_;
}
v_resetjp_389_:
{
lean_object* v___y_393_; lean_object* v___y_394_; lean_object* v___y_395_; lean_object* v___y_396_; lean_object* v___y_428_; lean_object* v___x_451_; uint8_t v___x_452_; 
v___x_451_ = lean_unsigned_to_nat(1024u);
v___x_452_ = lean_nat_dec_le(v___x_451_, v_prec_276_);
if (v___x_452_ == 0)
{
lean_object* v___x_453_; 
v___x_453_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9);
v___y_428_ = v___x_453_;
goto v___jp_427_;
}
else
{
lean_object* v___x_454_; 
v___x_454_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10);
v___y_428_ = v___x_454_;
goto v___jp_427_;
}
v___jp_392_:
{
lean_object* v_num_397_; lean_object* v_den_398_; lean_object* v___x_400_; uint8_t v_isShared_401_; uint8_t v_isSharedCheck_426_; 
v_num_397_ = lean_ctor_get(v_hi_388_, 0);
v_den_398_ = lean_ctor_get(v_hi_388_, 1);
v_isSharedCheck_426_ = !lean_is_exclusive(v_hi_388_);
if (v_isSharedCheck_426_ == 0)
{
v___x_400_ = v_hi_388_;
v_isShared_401_ = v_isSharedCheck_426_;
goto v_resetjp_399_;
}
else
{
lean_inc(v_den_398_);
lean_inc(v_num_397_);
lean_dec(v_hi_388_);
v___x_400_ = lean_box(0);
v_isShared_401_ = v_isSharedCheck_426_;
goto v_resetjp_399_;
}
v_resetjp_399_:
{
lean_object* v___x_403_; 
lean_inc(v___y_394_);
if (v_isShared_401_ == 0)
{
lean_ctor_set_tag(v___x_400_, 5);
lean_ctor_set(v___x_400_, 1, v___y_396_);
lean_ctor_set(v___x_400_, 0, v___y_394_);
v___x_403_ = v___x_400_;
goto v_reusejp_402_;
}
else
{
lean_object* v_reuseFailAlloc_425_; 
v_reuseFailAlloc_425_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_425_, 0, v___y_394_);
lean_ctor_set(v_reuseFailAlloc_425_, 1, v___y_396_);
v___x_403_ = v_reuseFailAlloc_425_;
goto v_reusejp_402_;
}
v_reusejp_402_:
{
lean_object* v___x_405_; 
lean_inc(v___y_395_);
if (v_isShared_391_ == 0)
{
lean_ctor_set_tag(v___x_390_, 5);
lean_ctor_set(v___x_390_, 1, v___y_395_);
lean_ctor_set(v___x_390_, 0, v___x_403_);
v___x_405_ = v___x_390_;
goto v_reusejp_404_;
}
else
{
lean_object* v_reuseFailAlloc_424_; 
v_reuseFailAlloc_424_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_424_, 0, v___x_403_);
lean_ctor_set(v_reuseFailAlloc_424_, 1, v___y_395_);
v___x_405_ = v_reuseFailAlloc_424_;
goto v_reusejp_404_;
}
v_reusejp_404_:
{
lean_object* v___x_406_; uint8_t v___x_407_; 
v___x_406_ = lean_unsigned_to_nat(1u);
v___x_407_ = lean_nat_dec_eq(v_den_398_, v___x_406_);
if (v___x_407_ == 0)
{
lean_object* v___x_408_; lean_object* v___x_409_; lean_object* v___x_410_; lean_object* v___x_411_; lean_object* v___x_412_; lean_object* v___x_413_; lean_object* v___x_414_; lean_object* v___x_415_; 
v___x_408_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14));
v___x_409_ = l_Int_repr(v_num_397_);
lean_dec(v_num_397_);
v___x_410_ = lean_string_append(v___x_408_, v___x_409_);
lean_dec_ref(v___x_409_);
v___x_411_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15));
v___x_412_ = lean_string_append(v___x_410_, v___x_411_);
v___x_413_ = l_Nat_reprFast(v_den_398_);
v___x_414_ = lean_string_append(v___x_412_, v___x_413_);
lean_dec_ref(v___x_413_);
v___x_415_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_415_, 0, v___x_414_);
v___y_294_ = v___y_393_;
v___y_295_ = v___x_405_;
v___y_296_ = v___x_415_;
goto v___jp_293_;
}
else
{
lean_object* v___x_416_; lean_object* v___x_417_; uint8_t v___x_418_; 
lean_dec(v_den_398_);
v___x_416_ = lean_unsigned_to_nat(0u);
v___x_417_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16);
v___x_418_ = lean_int_dec_lt(v_num_397_, v___x_417_);
if (v___x_418_ == 0)
{
lean_object* v___x_419_; lean_object* v___x_420_; 
v___x_419_ = l_Int_repr(v_num_397_);
lean_dec(v_num_397_);
v___x_420_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_420_, 0, v___x_419_);
v___y_294_ = v___y_393_;
v___y_295_ = v___x_405_;
v___y_296_ = v___x_420_;
goto v___jp_293_;
}
else
{
lean_object* v___x_421_; lean_object* v___x_422_; lean_object* v___x_423_; 
v___x_421_ = l_Int_repr(v_num_397_);
lean_dec(v_num_397_);
v___x_422_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_422_, 0, v___x_421_);
v___x_423_ = l_Repr_addAppParen(v___x_422_, v___x_416_);
v___y_294_ = v___y_393_;
v___y_295_ = v___x_405_;
v___y_296_ = v___x_423_;
goto v___jp_293_;
}
}
}
}
}
}
v___jp_427_:
{
lean_object* v_num_429_; lean_object* v_den_430_; lean_object* v___x_431_; lean_object* v___x_432_; lean_object* v___x_433_; uint8_t v___x_434_; 
v_num_429_ = lean_ctor_get(v_lo_387_, 0);
lean_inc(v_num_429_);
v_den_430_ = lean_ctor_get(v_lo_387_, 1);
lean_inc(v_den_430_);
lean_dec_ref(v_lo_387_);
v___x_431_ = lean_box(1);
v___x_432_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__19));
v___x_433_ = lean_unsigned_to_nat(1u);
v___x_434_ = lean_nat_dec_eq(v_den_430_, v___x_433_);
if (v___x_434_ == 0)
{
lean_object* v___x_435_; lean_object* v___x_436_; lean_object* v___x_437_; lean_object* v___x_438_; lean_object* v___x_439_; lean_object* v___x_440_; lean_object* v___x_441_; lean_object* v___x_442_; 
v___x_435_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__14));
v___x_436_ = l_Int_repr(v_num_429_);
lean_dec(v_num_429_);
v___x_437_ = lean_string_append(v___x_435_, v___x_436_);
lean_dec_ref(v___x_436_);
v___x_438_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__15));
v___x_439_ = lean_string_append(v___x_437_, v___x_438_);
v___x_440_ = l_Nat_reprFast(v_den_430_);
v___x_441_ = lean_string_append(v___x_439_, v___x_440_);
lean_dec_ref(v___x_440_);
v___x_442_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_442_, 0, v___x_441_);
v___y_393_ = v___y_428_;
v___y_394_ = v___x_432_;
v___y_395_ = v___x_431_;
v___y_396_ = v___x_442_;
goto v___jp_392_;
}
else
{
lean_object* v___x_443_; lean_object* v___x_444_; uint8_t v___x_445_; 
lean_dec(v_den_430_);
v___x_443_ = lean_unsigned_to_nat(0u);
v___x_444_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__16);
v___x_445_ = lean_int_dec_lt(v_num_429_, v___x_444_);
if (v___x_445_ == 0)
{
lean_object* v___x_446_; lean_object* v___x_447_; 
v___x_446_ = l_Int_repr(v_num_429_);
lean_dec(v_num_429_);
v___x_447_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_447_, 0, v___x_446_);
v___y_393_ = v___y_428_;
v___y_394_ = v___x_432_;
v___y_395_ = v___x_431_;
v___y_396_ = v___x_447_;
goto v___jp_392_;
}
else
{
lean_object* v___x_448_; lean_object* v___x_449_; lean_object* v___x_450_; 
v___x_448_ = l_Int_repr(v_num_429_);
lean_dec(v_num_429_);
v___x_449_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_449_, 0, v___x_448_);
v___x_450_ = l_Repr_addAppParen(v___x_449_, v___x_443_);
v___y_393_ = v___y_428_;
v___y_394_ = v___x_432_;
v___y_395_ = v___x_431_;
v___y_396_ = v___x_450_;
goto v___jp_392_;
}
}
}
}
}
case 4:
{
lean_object* v_enzymeSites_456_; lean_object* v___y_458_; lean_object* v___x_466_; uint8_t v___x_467_; 
v_enzymeSites_456_ = lean_ctor_get(v_x_275_, 0);
lean_inc(v_enzymeSites_456_);
lean_dec_ref(v_x_275_);
v___x_466_ = lean_unsigned_to_nat(1024u);
v___x_467_ = lean_nat_dec_le(v___x_466_, v_prec_276_);
if (v___x_467_ == 0)
{
lean_object* v___x_468_; 
v___x_468_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9);
v___y_458_ = v___x_468_;
goto v___jp_457_;
}
else
{
lean_object* v___x_469_; 
v___x_469_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10);
v___y_458_ = v___x_469_;
goto v___jp_457_;
}
v___jp_457_:
{
lean_object* v___x_459_; lean_object* v___x_460_; lean_object* v___x_461_; lean_object* v___x_462_; uint8_t v___x_463_; lean_object* v___x_464_; lean_object* v___x_465_; 
v___x_459_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__22));
v___x_460_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg(v_enzymeSites_456_);
v___x_461_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_461_, 0, v___x_459_);
lean_ctor_set(v___x_461_, 1, v___x_460_);
lean_inc(v___y_458_);
v___x_462_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_462_, 0, v___y_458_);
lean_ctor_set(v___x_462_, 1, v___x_461_);
v___x_463_ = 0;
v___x_464_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_464_, 0, v___x_462_);
lean_ctor_set_uint8(v___x_464_, sizeof(void*)*1, v___x_463_);
v___x_465_ = l_Repr_addAppParen(v___x_464_, v_prec_276_);
return v___x_465_;
}
}
case 5:
{
lean_object* v_readingFrame_470_; lean_object* v_exonBoundaries_471_; lean_object* v___x_473_; uint8_t v_isShared_474_; uint8_t v_isSharedCheck_495_; 
v_readingFrame_470_ = lean_ctor_get(v_x_275_, 0);
v_exonBoundaries_471_ = lean_ctor_get(v_x_275_, 1);
v_isSharedCheck_495_ = !lean_is_exclusive(v_x_275_);
if (v_isSharedCheck_495_ == 0)
{
v___x_473_ = v_x_275_;
v_isShared_474_ = v_isSharedCheck_495_;
goto v_resetjp_472_;
}
else
{
lean_inc(v_exonBoundaries_471_);
lean_inc(v_readingFrame_470_);
lean_dec(v_x_275_);
v___x_473_ = lean_box(0);
v_isShared_474_ = v_isSharedCheck_495_;
goto v_resetjp_472_;
}
v_resetjp_472_:
{
lean_object* v___y_476_; lean_object* v___x_491_; uint8_t v___x_492_; 
v___x_491_ = lean_unsigned_to_nat(1024u);
v___x_492_ = lean_nat_dec_le(v___x_491_, v_prec_276_);
if (v___x_492_ == 0)
{
lean_object* v___x_493_; 
v___x_493_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9);
v___y_476_ = v___x_493_;
goto v___jp_475_;
}
else
{
lean_object* v___x_494_; 
v___x_494_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10);
v___y_476_ = v___x_494_;
goto v___jp_475_;
}
v___jp_475_:
{
lean_object* v___x_477_; lean_object* v___x_478_; lean_object* v___x_479_; lean_object* v___x_480_; lean_object* v___x_482_; 
v___x_477_ = lean_box(1);
v___x_478_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__25));
v___x_479_ = l_Nat_reprFast(v_readingFrame_470_);
v___x_480_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_480_, 0, v___x_479_);
if (v_isShared_474_ == 0)
{
lean_ctor_set(v___x_473_, 1, v___x_480_);
lean_ctor_set(v___x_473_, 0, v___x_478_);
v___x_482_ = v___x_473_;
goto v_reusejp_481_;
}
else
{
lean_object* v_reuseFailAlloc_490_; 
v_reuseFailAlloc_490_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_490_, 0, v___x_478_);
lean_ctor_set(v_reuseFailAlloc_490_, 1, v___x_480_);
v___x_482_ = v_reuseFailAlloc_490_;
goto v_reusejp_481_;
}
v_reusejp_481_:
{
lean_object* v___x_483_; lean_object* v___x_484_; lean_object* v___x_485_; lean_object* v___x_486_; uint8_t v___x_487_; lean_object* v___x_488_; lean_object* v___x_489_; 
v___x_483_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_483_, 0, v___x_482_);
lean_ctor_set(v___x_483_, 1, v___x_477_);
v___x_484_ = lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___redArg(v_exonBoundaries_471_);
v___x_485_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_485_, 0, v___x_483_);
lean_ctor_set(v___x_485_, 1, v___x_484_);
lean_inc(v___y_476_);
v___x_486_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_486_, 0, v___y_476_);
lean_ctor_set(v___x_486_, 1, v___x_485_);
v___x_487_ = 0;
v___x_488_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_488_, 0, v___x_486_);
lean_ctor_set_uint8(v___x_488_, sizeof(void*)*1, v___x_487_);
v___x_489_ = l_Repr_addAppParen(v___x_488_, v_prec_276_);
return v___x_489_;
}
}
}
}
case 6:
{
lean_object* v___x_496_; uint8_t v___x_497_; 
v___x_496_ = lean_unsigned_to_nat(1024u);
v___x_497_ = lean_nat_dec_le(v___x_496_, v_prec_276_);
if (v___x_497_ == 0)
{
lean_object* v___x_498_; 
v___x_498_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9);
v___y_303_ = v___x_498_;
goto v___jp_302_;
}
else
{
lean_object* v___x_499_; 
v___x_499_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10);
v___y_303_ = v___x_499_;
goto v___jp_302_;
}
}
default: 
{
lean_object* v___x_500_; uint8_t v___x_501_; 
v___x_500_ = lean_unsigned_to_nat(1024u);
v___x_501_ = lean_nat_dec_le(v___x_500_, v_prec_276_);
if (v___x_501_ == 0)
{
lean_object* v___x_502_; 
v___x_502_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__9);
v___y_310_ = v___x_502_;
goto v___jp_309_;
}
else
{
lean_object* v___x_503_; 
v___x_503_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10, &lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__10);
v___y_310_ = v___x_503_;
goto v___jp_309_;
}
}
}
v___jp_277_:
{
lean_object* v___x_281_; lean_object* v___x_282_; uint8_t v___x_283_; lean_object* v___x_284_; lean_object* v___x_285_; 
v___x_281_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_281_, 0, v___y_278_);
lean_ctor_set(v___x_281_, 1, v___y_280_);
lean_inc(v___y_279_);
v___x_282_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_282_, 0, v___y_279_);
lean_ctor_set(v___x_282_, 1, v___x_281_);
v___x_283_ = 0;
v___x_284_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_284_, 0, v___x_282_);
lean_ctor_set_uint8(v___x_284_, sizeof(void*)*1, v___x_283_);
v___x_285_ = l_Repr_addAppParen(v___x_284_, v_prec_276_);
return v___x_285_;
}
v___jp_286_:
{
lean_object* v___x_288_; lean_object* v___x_289_; uint8_t v___x_290_; lean_object* v___x_291_; lean_object* v___x_292_; 
v___x_288_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__1));
lean_inc(v___y_287_);
v___x_289_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_289_, 0, v___y_287_);
lean_ctor_set(v___x_289_, 1, v___x_288_);
v___x_290_ = 0;
v___x_291_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_291_, 0, v___x_289_);
lean_ctor_set_uint8(v___x_291_, sizeof(void*)*1, v___x_290_);
v___x_292_ = l_Repr_addAppParen(v___x_291_, v_prec_276_);
return v___x_292_;
}
v___jp_293_:
{
lean_object* v___x_297_; lean_object* v___x_298_; uint8_t v___x_299_; lean_object* v___x_300_; lean_object* v___x_301_; 
v___x_297_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_297_, 0, v___y_295_);
lean_ctor_set(v___x_297_, 1, v___y_296_);
lean_inc(v___y_294_);
v___x_298_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_298_, 0, v___y_294_);
lean_ctor_set(v___x_298_, 1, v___x_297_);
v___x_299_ = 0;
v___x_300_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_300_, 0, v___x_298_);
lean_ctor_set_uint8(v___x_300_, sizeof(void*)*1, v___x_299_);
v___x_301_ = l_Repr_addAppParen(v___x_300_, v_prec_276_);
return v___x_301_;
}
v___jp_302_:
{
lean_object* v___x_304_; lean_object* v___x_305_; uint8_t v___x_306_; lean_object* v___x_307_; lean_object* v___x_308_; 
v___x_304_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__3));
lean_inc(v___y_303_);
v___x_305_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_305_, 0, v___y_303_);
lean_ctor_set(v___x_305_, 1, v___x_304_);
v___x_306_ = 0;
v___x_307_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_307_, 0, v___x_305_);
lean_ctor_set_uint8(v___x_307_, sizeof(void*)*1, v___x_306_);
v___x_308_ = l_Repr_addAppParen(v___x_307_, v_prec_276_);
return v___x_308_;
}
v___jp_309_:
{
lean_object* v___x_311_; lean_object* v___x_312_; uint8_t v___x_313_; lean_object* v___x_314_; lean_object* v___x_315_; 
v___x_311_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___closed__5));
lean_inc(v___y_310_);
v___x_312_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_312_, 0, v___y_310_);
lean_ctor_set(v___x_312_, 1, v___x_311_);
v___x_313_ = 0;
v___x_314_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_314_, 0, v___x_312_);
lean_ctor_set_uint8(v___x_314_, sizeof(void*)*1, v___x_313_);
v___x_315_ = l_Repr_addAppParen(v___x_314_, v_prec_276_);
return v___x_315_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprTypePredicate_repr___boxed(lean_object* v_x_504_, lean_object* v_prec_505_){
_start:
{
lean_object* v_res_506_; 
v_res_506_ = lp_BioCompiler_BioCompiler_instReprTypePredicate_repr(v_x_504_, v_prec_505_);
lean_dec(v_prec_505_);
return v_res_506_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0(lean_object* v_a_507_, lean_object* v_n_508_){
_start:
{
lean_object* v___x_509_; 
v___x_509_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___redArg(v_a_507_);
return v___x_509_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0___boxed(lean_object* v_a_510_, lean_object* v_n_511_){
_start:
{
lean_object* v_res_512_; 
v_res_512_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprTypePredicate_repr_spec__0(v_a_510_, v_n_511_);
lean_dec(v_n_511_);
return v_res_512_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1(lean_object* v_a_513_, lean_object* v_n_514_){
_start:
{
lean_object* v___x_515_; 
v___x_515_ = lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___redArg(v_a_513_);
return v___x_515_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1___boxed(lean_object* v_a_516_, lean_object* v_n_517_){
_start:
{
lean_object* v_res_518_; 
v_res_518_ = lp_BioCompiler_List_repr_x27___at___00BioCompiler_instReprTypePredicate_repr_spec__1(v_a_516_, v_n_517_);
lean_dec(v_n_517_);
return v_res_518_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluate___redArg(lean_object* v_inst_521_, lean_object* v_inst_522_, lean_object* v_inst_523_, lean_object* v_inst_524_, lean_object* v_x_525_, lean_object* v_x_526_, lean_object* v_x_527_){
_start:
{
switch(lean_obj_tag(v_x_525_))
{
case 0:
{
lean_object* v_cellType_528_; lean_object* v_cellType_529_; uint8_t v___x_530_; 
lean_dec_ref(v_inst_523_);
lean_dec_ref(v_inst_522_);
lean_dec_ref(v_inst_521_);
v_cellType_528_ = lean_ctor_get(v_x_525_, 0);
lean_inc_ref(v_cellType_528_);
lean_dec_ref(v_x_525_);
v_cellType_529_ = lean_ctor_get(v_x_527_, 0);
v___x_530_ = lean_string_dec_eq(v_cellType_529_, v_cellType_528_);
lean_dec_ref(v_cellType_528_);
if (v___x_530_ == 0)
{
uint8_t v___x_531_; 
lean_dec(v_x_526_);
lean_dec_ref(v_inst_524_);
v___x_531_ = 2;
return v___x_531_;
}
else
{
lean_object* v___x_532_; 
v___x_532_ = lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg(v_inst_524_, v_x_526_);
lean_dec(v_x_526_);
if (lean_obj_tag(v___x_532_) == 1)
{
lean_object* v_tail_533_; 
v_tail_533_ = lean_ctor_get(v___x_532_, 1);
lean_inc(v_tail_533_);
lean_dec_ref(v___x_532_);
if (lean_obj_tag(v_tail_533_) == 0)
{
uint8_t v___x_534_; 
v___x_534_ = 0;
return v___x_534_;
}
else
{
uint8_t v___x_535_; 
lean_dec(v_tail_533_);
v___x_535_ = 1;
return v___x_535_;
}
}
else
{
uint8_t v___x_536_; 
lean_dec(v___x_532_);
v___x_536_ = 1;
return v___x_536_;
}
}
}
case 1:
{
lean_object* v_hasCrypticSpliceSite_537_; lean_object* v_hasBorderlineSpliceSite_538_; lean_object* v___x_539_; uint8_t v___x_540_; 
lean_dec_ref(v_inst_524_);
lean_dec_ref(v_inst_523_);
lean_dec_ref(v_inst_522_);
v_hasCrypticSpliceSite_537_ = lean_ctor_get(v_inst_521_, 0);
lean_inc_ref(v_hasCrypticSpliceSite_537_);
v_hasBorderlineSpliceSite_538_ = lean_ctor_get(v_inst_521_, 1);
lean_inc_ref(v_hasBorderlineSpliceSite_538_);
lean_dec_ref(v_inst_521_);
lean_inc(v_x_526_);
v___x_539_ = lean_apply_1(v_hasCrypticSpliceSite_537_, v_x_526_);
v___x_540_ = lean_unbox(v___x_539_);
if (v___x_540_ == 0)
{
lean_object* v___x_541_; uint8_t v___x_542_; 
v___x_541_ = lean_apply_1(v_hasBorderlineSpliceSite_538_, v_x_526_);
v___x_542_ = lean_unbox(v___x_541_);
if (v___x_542_ == 0)
{
uint8_t v___x_543_; 
v___x_543_ = 0;
return v___x_543_;
}
else
{
uint8_t v___x_544_; 
v___x_544_ = 2;
return v___x_544_;
}
}
else
{
uint8_t v___x_545_; 
lean_dec_ref(v_hasBorderlineSpliceSite_538_);
lean_dec(v_x_526_);
v___x_545_ = 1;
return v___x_545_;
}
}
case 2:
{
lean_object* v_organism_546_; lean_object* v_threshold_547_; lean_object* v___x_548_; uint8_t v___x_549_; 
lean_dec_ref(v_inst_524_);
lean_dec_ref(v_inst_523_);
lean_dec_ref(v_inst_521_);
v_organism_546_ = lean_ctor_get(v_x_525_, 0);
lean_inc_ref(v_organism_546_);
v_threshold_547_ = lean_ctor_get(v_x_525_, 1);
lean_inc_ref(v_threshold_547_);
lean_dec_ref(v_x_525_);
v___x_548_ = lean_apply_2(v_inst_522_, v_x_526_, v_organism_546_);
v___x_549_ = l_Rat_instDecidableLe(v_threshold_547_, v___x_548_);
if (v___x_549_ == 0)
{
uint8_t v___x_550_; 
v___x_550_ = 1;
return v___x_550_;
}
else
{
uint8_t v___x_551_; 
v___x_551_ = 0;
return v___x_551_;
}
}
case 3:
{
lean_object* v_lo_552_; lean_object* v_hi_553_; lean_object* v___x_554_; uint8_t v___x_555_; 
lean_dec_ref(v_inst_524_);
lean_dec_ref(v_inst_523_);
lean_dec_ref(v_inst_522_);
lean_dec_ref(v_inst_521_);
v_lo_552_ = lean_ctor_get(v_x_525_, 0);
lean_inc_ref(v_lo_552_);
v_hi_553_ = lean_ctor_get(v_x_525_, 1);
lean_inc_ref(v_hi_553_);
lean_dec_ref(v_x_525_);
v___x_554_ = lp_BioCompiler_BioCompiler_Sequence_gcContent(v_x_526_);
lean_dec(v_x_526_);
lean_inc_ref(v___x_554_);
v___x_555_ = l_Rat_instDecidableLe(v_lo_552_, v___x_554_);
if (v___x_555_ == 0)
{
uint8_t v___x_556_; 
lean_dec_ref(v___x_554_);
lean_dec_ref(v_hi_553_);
v___x_556_ = 1;
return v___x_556_;
}
else
{
uint8_t v___x_557_; 
v___x_557_ = l_Rat_instDecidableLe(v___x_554_, v_hi_553_);
if (v___x_557_ == 0)
{
uint8_t v___x_558_; 
v___x_558_ = 1;
return v___x_558_;
}
else
{
uint8_t v___x_559_; 
v___x_559_ = 0;
return v___x_559_;
}
}
}
case 4:
{
lean_object* v_enzymeSites_560_; uint8_t v___x_561_; 
lean_dec_ref(v_inst_524_);
lean_dec_ref(v_inst_523_);
lean_dec_ref(v_inst_522_);
lean_dec_ref(v_inst_521_);
v_enzymeSites_560_ = lean_ctor_get(v_x_525_, 0);
lean_inc(v_enzymeSites_560_);
lean_dec_ref(v_x_525_);
v___x_561_ = lp_BioCompiler_BioCompiler_hasAnyRestrictionSite(v_x_526_, v_enzymeSites_560_);
if (v___x_561_ == 0)
{
uint8_t v___x_562_; 
v___x_562_ = 0;
return v___x_562_;
}
else
{
uint8_t v___x_563_; 
v___x_563_ = 1;
return v___x_563_;
}
}
case 5:
{
lean_object* v_readingFrame_564_; lean_object* v_exonBoundaries_565_; uint8_t v___x_566_; 
lean_dec_ref(v_inst_524_);
lean_dec_ref(v_inst_523_);
lean_dec_ref(v_inst_522_);
lean_dec_ref(v_inst_521_);
v_readingFrame_564_ = lean_ctor_get(v_x_525_, 0);
lean_inc_n(v_readingFrame_564_, 2);
v_exonBoundaries_565_ = lean_ctor_get(v_x_525_, 1);
lean_inc(v_exonBoundaries_565_);
lean_dec_ref(v_x_525_);
v___x_566_ = lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent(v_exonBoundaries_565_, v_readingFrame_564_);
if (v___x_566_ == 0)
{
uint8_t v___x_567_; 
lean_dec(v_readingFrame_564_);
lean_dec(v_x_526_);
v___x_567_ = 1;
return v___x_567_;
}
else
{
uint8_t v___x_568_; 
v___x_568_ = lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(v_x_526_, v_readingFrame_564_);
lean_dec(v_readingFrame_564_);
if (v___x_568_ == 0)
{
uint8_t v___x_569_; 
v___x_569_ = 0;
return v___x_569_;
}
else
{
uint8_t v___x_570_; 
v___x_570_ = 1;
return v___x_570_;
}
}
}
case 6:
{
uint8_t v___x_571_; 
lean_dec_ref(v_inst_524_);
lean_dec_ref(v_inst_523_);
lean_dec_ref(v_inst_522_);
lean_dec_ref(v_inst_521_);
v___x_571_ = lp_BioCompiler_BioCompiler_hasInstabilityMotif(v_x_526_);
if (v___x_571_ == 0)
{
uint8_t v___x_572_; 
v___x_572_ = 0;
return v___x_572_;
}
else
{
uint8_t v___x_573_; 
v___x_573_ = 1;
return v___x_573_;
}
}
default: 
{
lean_object* v___x_574_; uint8_t v___x_575_; 
lean_dec_ref(v_inst_524_);
lean_dec_ref(v_inst_522_);
lean_dec_ref(v_inst_521_);
v___x_574_ = lean_apply_1(v_inst_523_, v_x_526_);
v___x_575_ = lean_unbox(v___x_574_);
if (v___x_575_ == 0)
{
uint8_t v___x_576_; 
v___x_576_ = 0;
return v___x_576_;
}
else
{
uint8_t v___x_577_; 
v___x_577_ = 1;
return v___x_577_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluate___redArg___boxed(lean_object* v_inst_578_, lean_object* v_inst_579_, lean_object* v_inst_580_, lean_object* v_inst_581_, lean_object* v_x_582_, lean_object* v_x_583_, lean_object* v_x_584_){
_start:
{
uint8_t v_res_585_; lean_object* v_r_586_; 
v_res_585_ = lp_BioCompiler_BioCompiler_evaluate___redArg(v_inst_578_, v_inst_579_, v_inst_580_, v_inst_581_, v_x_582_, v_x_583_, v_x_584_);
lean_dec_ref(v_x_584_);
v_r_586_ = lean_box(v_res_585_);
return v_r_586_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluate(lean_object* v_inst_587_, lean_object* v_inst_588_, lean_object* v_inst_589_, lean_object* v_State_590_, lean_object* v_inst_591_, lean_object* v_inst_592_, lean_object* v_inst_593_, lean_object* v_x_594_, lean_object* v_x_595_, lean_object* v_x_596_){
_start:
{
uint8_t v___x_597_; 
v___x_597_ = lp_BioCompiler_BioCompiler_evaluate___redArg(v_inst_587_, v_inst_588_, v_inst_589_, v_inst_593_, v_x_594_, v_x_595_, v_x_596_);
return v___x_597_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluate___boxed(lean_object* v_inst_598_, lean_object* v_inst_599_, lean_object* v_inst_600_, lean_object* v_State_601_, lean_object* v_inst_602_, lean_object* v_inst_603_, lean_object* v_inst_604_, lean_object* v_x_605_, lean_object* v_x_606_, lean_object* v_x_607_){
_start:
{
uint8_t v_res_608_; lean_object* v_r_609_; 
v_res_608_ = lp_BioCompiler_BioCompiler_evaluate(v_inst_598_, v_inst_599_, v_inst_600_, v_State_601_, v_inst_602_, v_inst_603_, v_inst_604_, v_x_605_, v_x_606_, v_x_607_);
lean_dec_ref(v_x_607_);
lean_dec(v_inst_603_);
lean_dec_ref(v_inst_602_);
v_r_609_ = lean_box(v_res_608_);
return v_r_609_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter___redArg(lean_object* v_x_610_, lean_object* v_x_611_, lean_object* v_x_612_, lean_object* v_h__1_613_, lean_object* v_h__2_614_, lean_object* v_h__3_615_, lean_object* v_h__4_616_, lean_object* v_h__5_617_, lean_object* v_h__6_618_, lean_object* v_h__7_619_, lean_object* v_h__8_620_){
_start:
{
switch(lean_obj_tag(v_x_610_))
{
case 0:
{
lean_object* v_cellType_621_; lean_object* v___x_622_; 
lean_dec(v_h__8_620_);
lean_dec(v_h__7_619_);
lean_dec(v_h__6_618_);
lean_dec(v_h__5_617_);
lean_dec(v_h__4_616_);
lean_dec(v_h__3_615_);
lean_dec(v_h__2_614_);
v_cellType_621_ = lean_ctor_get(v_x_610_, 0);
lean_inc_ref(v_cellType_621_);
lean_dec_ref(v_x_610_);
v___x_622_ = lean_apply_3(v_h__1_613_, v_cellType_621_, v_x_611_, v_x_612_);
return v___x_622_;
}
case 1:
{
lean_object* v___x_623_; 
lean_dec(v_h__8_620_);
lean_dec(v_h__7_619_);
lean_dec(v_h__6_618_);
lean_dec(v_h__5_617_);
lean_dec(v_h__4_616_);
lean_dec(v_h__3_615_);
lean_dec(v_h__1_613_);
v___x_623_ = lean_apply_2(v_h__2_614_, v_x_611_, v_x_612_);
return v___x_623_;
}
case 2:
{
lean_object* v_organism_624_; lean_object* v_threshold_625_; lean_object* v___x_626_; 
lean_dec(v_h__8_620_);
lean_dec(v_h__7_619_);
lean_dec(v_h__6_618_);
lean_dec(v_h__5_617_);
lean_dec(v_h__4_616_);
lean_dec(v_h__2_614_);
lean_dec(v_h__1_613_);
v_organism_624_ = lean_ctor_get(v_x_610_, 0);
lean_inc_ref(v_organism_624_);
v_threshold_625_ = lean_ctor_get(v_x_610_, 1);
lean_inc_ref(v_threshold_625_);
lean_dec_ref(v_x_610_);
v___x_626_ = lean_apply_4(v_h__3_615_, v_organism_624_, v_threshold_625_, v_x_611_, v_x_612_);
return v___x_626_;
}
case 3:
{
lean_object* v_lo_627_; lean_object* v_hi_628_; lean_object* v___x_629_; 
lean_dec(v_h__8_620_);
lean_dec(v_h__7_619_);
lean_dec(v_h__6_618_);
lean_dec(v_h__5_617_);
lean_dec(v_h__3_615_);
lean_dec(v_h__2_614_);
lean_dec(v_h__1_613_);
v_lo_627_ = lean_ctor_get(v_x_610_, 0);
lean_inc_ref(v_lo_627_);
v_hi_628_ = lean_ctor_get(v_x_610_, 1);
lean_inc_ref(v_hi_628_);
lean_dec_ref(v_x_610_);
v___x_629_ = lean_apply_4(v_h__4_616_, v_lo_627_, v_hi_628_, v_x_611_, v_x_612_);
return v___x_629_;
}
case 4:
{
lean_object* v_enzymeSites_630_; lean_object* v___x_631_; 
lean_dec(v_h__8_620_);
lean_dec(v_h__7_619_);
lean_dec(v_h__6_618_);
lean_dec(v_h__4_616_);
lean_dec(v_h__3_615_);
lean_dec(v_h__2_614_);
lean_dec(v_h__1_613_);
v_enzymeSites_630_ = lean_ctor_get(v_x_610_, 0);
lean_inc(v_enzymeSites_630_);
lean_dec_ref(v_x_610_);
v___x_631_ = lean_apply_3(v_h__5_617_, v_enzymeSites_630_, v_x_611_, v_x_612_);
return v___x_631_;
}
case 5:
{
lean_object* v_readingFrame_632_; lean_object* v_exonBoundaries_633_; lean_object* v___x_634_; 
lean_dec(v_h__8_620_);
lean_dec(v_h__7_619_);
lean_dec(v_h__5_617_);
lean_dec(v_h__4_616_);
lean_dec(v_h__3_615_);
lean_dec(v_h__2_614_);
lean_dec(v_h__1_613_);
v_readingFrame_632_ = lean_ctor_get(v_x_610_, 0);
lean_inc(v_readingFrame_632_);
v_exonBoundaries_633_ = lean_ctor_get(v_x_610_, 1);
lean_inc(v_exonBoundaries_633_);
lean_dec_ref(v_x_610_);
v___x_634_ = lean_apply_4(v_h__6_618_, v_readingFrame_632_, v_exonBoundaries_633_, v_x_611_, v_x_612_);
return v___x_634_;
}
case 6:
{
lean_object* v___x_635_; 
lean_dec(v_h__8_620_);
lean_dec(v_h__6_618_);
lean_dec(v_h__5_617_);
lean_dec(v_h__4_616_);
lean_dec(v_h__3_615_);
lean_dec(v_h__2_614_);
lean_dec(v_h__1_613_);
v___x_635_ = lean_apply_2(v_h__7_619_, v_x_611_, v_x_612_);
return v___x_635_;
}
default: 
{
lean_object* v___x_636_; 
lean_dec(v_h__7_619_);
lean_dec(v_h__6_618_);
lean_dec(v_h__5_617_);
lean_dec(v_h__4_616_);
lean_dec(v_h__3_615_);
lean_dec(v_h__2_614_);
lean_dec(v_h__1_613_);
v___x_636_ = lean_apply_2(v_h__8_620_, v_x_611_, v_x_612_);
return v___x_636_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__5_splitter(lean_object* v_motive_637_, lean_object* v_x_638_, lean_object* v_x_639_, lean_object* v_x_640_, lean_object* v_h__1_641_, lean_object* v_h__2_642_, lean_object* v_h__3_643_, lean_object* v_h__4_644_, lean_object* v_h__5_645_, lean_object* v_h__6_646_, lean_object* v_h__7_647_, lean_object* v_h__8_648_){
_start:
{
switch(lean_obj_tag(v_x_638_))
{
case 0:
{
lean_object* v_cellType_649_; lean_object* v___x_650_; 
lean_dec(v_h__8_648_);
lean_dec(v_h__7_647_);
lean_dec(v_h__6_646_);
lean_dec(v_h__5_645_);
lean_dec(v_h__4_644_);
lean_dec(v_h__3_643_);
lean_dec(v_h__2_642_);
v_cellType_649_ = lean_ctor_get(v_x_638_, 0);
lean_inc_ref(v_cellType_649_);
lean_dec_ref(v_x_638_);
v___x_650_ = lean_apply_3(v_h__1_641_, v_cellType_649_, v_x_639_, v_x_640_);
return v___x_650_;
}
case 1:
{
lean_object* v___x_651_; 
lean_dec(v_h__8_648_);
lean_dec(v_h__7_647_);
lean_dec(v_h__6_646_);
lean_dec(v_h__5_645_);
lean_dec(v_h__4_644_);
lean_dec(v_h__3_643_);
lean_dec(v_h__1_641_);
v___x_651_ = lean_apply_2(v_h__2_642_, v_x_639_, v_x_640_);
return v___x_651_;
}
case 2:
{
lean_object* v_organism_652_; lean_object* v_threshold_653_; lean_object* v___x_654_; 
lean_dec(v_h__8_648_);
lean_dec(v_h__7_647_);
lean_dec(v_h__6_646_);
lean_dec(v_h__5_645_);
lean_dec(v_h__4_644_);
lean_dec(v_h__2_642_);
lean_dec(v_h__1_641_);
v_organism_652_ = lean_ctor_get(v_x_638_, 0);
lean_inc_ref(v_organism_652_);
v_threshold_653_ = lean_ctor_get(v_x_638_, 1);
lean_inc_ref(v_threshold_653_);
lean_dec_ref(v_x_638_);
v___x_654_ = lean_apply_4(v_h__3_643_, v_organism_652_, v_threshold_653_, v_x_639_, v_x_640_);
return v___x_654_;
}
case 3:
{
lean_object* v_lo_655_; lean_object* v_hi_656_; lean_object* v___x_657_; 
lean_dec(v_h__8_648_);
lean_dec(v_h__7_647_);
lean_dec(v_h__6_646_);
lean_dec(v_h__5_645_);
lean_dec(v_h__3_643_);
lean_dec(v_h__2_642_);
lean_dec(v_h__1_641_);
v_lo_655_ = lean_ctor_get(v_x_638_, 0);
lean_inc_ref(v_lo_655_);
v_hi_656_ = lean_ctor_get(v_x_638_, 1);
lean_inc_ref(v_hi_656_);
lean_dec_ref(v_x_638_);
v___x_657_ = lean_apply_4(v_h__4_644_, v_lo_655_, v_hi_656_, v_x_639_, v_x_640_);
return v___x_657_;
}
case 4:
{
lean_object* v_enzymeSites_658_; lean_object* v___x_659_; 
lean_dec(v_h__8_648_);
lean_dec(v_h__7_647_);
lean_dec(v_h__6_646_);
lean_dec(v_h__4_644_);
lean_dec(v_h__3_643_);
lean_dec(v_h__2_642_);
lean_dec(v_h__1_641_);
v_enzymeSites_658_ = lean_ctor_get(v_x_638_, 0);
lean_inc(v_enzymeSites_658_);
lean_dec_ref(v_x_638_);
v___x_659_ = lean_apply_3(v_h__5_645_, v_enzymeSites_658_, v_x_639_, v_x_640_);
return v___x_659_;
}
case 5:
{
lean_object* v_readingFrame_660_; lean_object* v_exonBoundaries_661_; lean_object* v___x_662_; 
lean_dec(v_h__8_648_);
lean_dec(v_h__7_647_);
lean_dec(v_h__5_645_);
lean_dec(v_h__4_644_);
lean_dec(v_h__3_643_);
lean_dec(v_h__2_642_);
lean_dec(v_h__1_641_);
v_readingFrame_660_ = lean_ctor_get(v_x_638_, 0);
lean_inc(v_readingFrame_660_);
v_exonBoundaries_661_ = lean_ctor_get(v_x_638_, 1);
lean_inc(v_exonBoundaries_661_);
lean_dec_ref(v_x_638_);
v___x_662_ = lean_apply_4(v_h__6_646_, v_readingFrame_660_, v_exonBoundaries_661_, v_x_639_, v_x_640_);
return v___x_662_;
}
case 6:
{
lean_object* v___x_663_; 
lean_dec(v_h__8_648_);
lean_dec(v_h__6_646_);
lean_dec(v_h__5_645_);
lean_dec(v_h__4_644_);
lean_dec(v_h__3_643_);
lean_dec(v_h__2_642_);
lean_dec(v_h__1_641_);
v___x_663_ = lean_apply_2(v_h__7_647_, v_x_639_, v_x_640_);
return v___x_663_;
}
default: 
{
lean_object* v___x_664_; 
lean_dec(v_h__7_647_);
lean_dec(v_h__6_646_);
lean_dec(v_h__5_645_);
lean_dec(v_h__4_644_);
lean_dec(v_h__3_643_);
lean_dec(v_h__2_642_);
lean_dec(v_h__1_641_);
v___x_664_ = lean_apply_2(v_h__8_648_, v_x_639_, v_x_640_);
return v___x_664_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__1_splitter___redArg(lean_object* v_x_665_, lean_object* v_h__1_666_, lean_object* v_h__2_667_){
_start:
{
if (lean_obj_tag(v_x_665_) == 1)
{
lean_object* v_tail_668_; 
v_tail_668_ = lean_ctor_get(v_x_665_, 1);
if (lean_obj_tag(v_tail_668_) == 0)
{
lean_object* v_head_669_; lean_object* v___x_670_; 
lean_dec(v_h__2_667_);
v_head_669_ = lean_ctor_get(v_x_665_, 0);
lean_inc(v_head_669_);
lean_dec_ref(v_x_665_);
v___x_670_ = lean_apply_1(v_h__1_666_, v_head_669_);
return v___x_670_;
}
else
{
lean_object* v___x_671_; 
lean_dec(v_h__1_666_);
v___x_671_ = lean_apply_2(v_h__2_667_, v_x_665_, lean_box(0));
return v___x_671_;
}
}
else
{
lean_object* v___x_672_; 
lean_dec(v_h__1_666_);
v___x_672_ = lean_apply_2(v_h__2_667_, v_x_665_, lean_box(0));
return v___x_672_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_TypeSystem_0__BioCompiler_evaluate_match__1_splitter(lean_object* v_motive_673_, lean_object* v_x_674_, lean_object* v_h__1_675_, lean_object* v_h__2_676_){
_start:
{
if (lean_obj_tag(v_x_674_) == 1)
{
lean_object* v_tail_677_; 
v_tail_677_ = lean_ctor_get(v_x_674_, 1);
if (lean_obj_tag(v_tail_677_) == 0)
{
lean_object* v_head_678_; lean_object* v___x_679_; 
lean_dec(v_h__2_676_);
v_head_678_ = lean_ctor_get(v_x_674_, 0);
lean_inc(v_head_678_);
lean_dec_ref(v_x_674_);
v___x_679_ = lean_apply_1(v_h__1_675_, v_head_678_);
return v___x_679_;
}
else
{
lean_object* v___x_680_; 
lean_dec(v_h__1_675_);
v___x_680_ = lean_apply_2(v_h__2_676_, v_x_674_, lean_box(0));
return v___x_680_;
}
}
else
{
lean_object* v___x_681_; 
lean_dec(v_h__1_675_);
v___x_681_ = lean_apply_2(v_h__2_676_, v_x_674_, lean_box(0));
return v___x_681_;
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
