// Lean compiler output
// Module: BioCompiler.NDFST
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
lean_object* lean_mk_empty_array_with_capacity(lean_object*);
lean_object* lean_array_to_list(lean_object*);
lean_object* l_List_reverse___redArg(lean_object*);
lean_object* l_List_appendTR___redArg(lean_object*, lean_object*);
lean_object* l_List_foldl___at___00Array_appendList_spec__0___redArg(lean_object*, lean_object*);
lean_object* lean_array_push(lean_object*, lean_object*);
uint8_t lp_BioCompiler_BioCompiler_instBEqNucleotide_beq(uint8_t, uint8_t);
lean_object* l_List_eraseDupsBy___redArg(lean_object*, lean_object*);
lean_object* lean_nat_to_int(lean_object*);
lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr(uint8_t, lean_object*);
lean_object* lean_string_length(lean_object*);
lean_object* l_String_quote(lean_object*);
uint8_t lean_nat_dec_eq(lean_object*, lean_object*);
lean_object* l_Int_repr(lean_object*);
lean_object* lean_string_append(lean_object*, lean_object*);
lean_object* l_Nat_reprFast(lean_object*);
uint8_t lean_int_dec_lt(lean_object*, lean_object*);
lean_object* l_Repr_addAppParen(lean_object*, lean_object*);
lean_object* l_Std_Format_joinSep___at___00Lean_Syntax_formatStxAux_spec__2(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_ndfstStep_spec__0___redArg(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1___redArg(lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1___redArg___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
static const lean_array_object lp_BioCompiler_BioCompiler_ndfstStep___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_array_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 246}, .m_size = 0, .m_capacity = 0, .m_data = {}};
static const lean_object* lp_BioCompiler_BioCompiler_ndfstStep___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_ndfstStep___redArg___closed__0_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstStep___redArg(lean_object*, lean_object*, uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstStep___redArg___boxed(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstStep(lean_object*, lean_object*, lean_object*, uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstStep___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_ndfstStep_spec__0(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1(lean_object*, lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0___redArg(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0___redArg___boxed(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstRun___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstRun___redArg___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstRun(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstRun___boxed(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_filterMapTR_go___at___00BioCompiler_ndfstOutputSet_spec__0___redArg(lean_object*, lean_object*, lean_object*);
static const lean_array_object lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_array_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 246}, .m_size = 0, .m_capacity = 0, .m_data = {}};
static const lean_object* lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg___closed__0_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstOutputSet(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstOutputSet___boxed(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_filterMapTR_go___at___00BioCompiler_ndfstOutputSet_spec__0(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_List_beq___at___00List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0_spec__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_beq___at___00List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0_spec__0___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_List_beq___at___00List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0_spec__0___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0___closed__0 = (const lean_object*)&lp_BioCompiler_List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0___closed__0_value;
LEAN_EXPORT lean_object* lp_BioCompiler_List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___boxed(lean_object*, lean_object*, lean_object*);
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "cellType"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__2_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 5, .m_capacity = 5, .m_length = 4, .m_data = " := "};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__3_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__3_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__5_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 3, .m_capacity = 3, .m_length = 2, .m_data = "{ "};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__6 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__6_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__7_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__7;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = ","};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__8 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__8_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__9_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__8_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__9 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__9_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__10_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 13, .m_capacity = 13, .m_length = 12, .m_data = "eseThreshold"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__10 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__10_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__11_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__10_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__11 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__11_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__12_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__12;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__13_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 3, .m_capacity = 3, .m_length = 2, .m_data = " }"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__13 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__13_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__14_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__14;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__15_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__15;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__16_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__6_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__16 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__16_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__17_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__13_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__17 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__17_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__18_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 13, .m_capacity = 13, .m_length = 12, .m_data = "issThreshold"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__18 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__18_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__19_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__18_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__19 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__19_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__20_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "("};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__20 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__20_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__21_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = " : Rat)/"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__21 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__21_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__23_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 13, .m_capacity = 13, .m_length = 12, .m_data = "iseThreshold"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__23 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__23_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__24_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__23_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__24 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__24_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__25_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 13, .m_capacity = 13, .m_length = 12, .m_data = "essThreshold"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__25 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__25_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__26_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__25_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__26 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__26_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprCellularContext___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprCellularContext_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext___closed__0_value;
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0_spec__1_spec__3(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0_spec__1(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0___lam__0(uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0___lam__0___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0(lean_object*, lean_object*);
static const lean_string_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 3, .m_capacity = 3, .m_length = 2, .m_data = "[]"};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__1_value;
static const lean_string_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "["};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__9_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__3_value;
static const lean_string_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "]"};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__4_value;
static lean_once_cell_t lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__5_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__5;
static lean_once_cell_t lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__6_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__6;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__7_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__2_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__7 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__7_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__8 = (const lean_object*)&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__8_value;
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(lean_object*);
static const lean_string_object lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = ")"};
static const lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__0_value;
static lean_once_cell_t lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__1_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__1;
static lean_once_cell_t lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__2_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__2;
static const lean_ctor_object lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__20_value)}};
static const lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__3_value;
static const lean_ctor_object lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__4_value;
LEAN_EXPORT lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__3_spec__5_spec__7(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__3_spec__5(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__3(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1___redArg(lean_object*);
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "sequence"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 15, .m_capacity = 15, .m_length = 14, .m_data = "exonBoundaries"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__5_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__6_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__6;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprSpliceIsoform___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSpliceIsoform___closed__0_value;
LEAN_EXPORT lean_object* lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_ndfstStep_spec__0___redArg(lean_object* v_snd_1_, lean_object* v_a_2_, lean_object* v_a_3_){
_start:
{
if (lean_obj_tag(v_a_2_) == 0)
{
lean_object* v___x_4_; 
lean_dec(v_snd_1_);
v___x_4_ = l_List_reverse___redArg(v_a_3_);
return v___x_4_;
}
else
{
lean_object* v_head_5_; lean_object* v_tail_6_; lean_object* v___x_8_; uint8_t v_isShared_9_; uint8_t v_isSharedCheck_24_; 
v_head_5_ = lean_ctor_get(v_a_2_, 0);
v_tail_6_ = lean_ctor_get(v_a_2_, 1);
v_isSharedCheck_24_ = !lean_is_exclusive(v_a_2_);
if (v_isSharedCheck_24_ == 0)
{
v___x_8_ = v_a_2_;
v_isShared_9_ = v_isSharedCheck_24_;
goto v_resetjp_7_;
}
else
{
lean_inc(v_tail_6_);
lean_inc(v_head_5_);
lean_dec(v_a_2_);
v___x_8_ = lean_box(0);
v_isShared_9_ = v_isSharedCheck_24_;
goto v_resetjp_7_;
}
v_resetjp_7_:
{
lean_object* v_fst_10_; lean_object* v_snd_11_; lean_object* v___x_13_; uint8_t v_isShared_14_; uint8_t v_isSharedCheck_23_; 
v_fst_10_ = lean_ctor_get(v_head_5_, 0);
v_snd_11_ = lean_ctor_get(v_head_5_, 1);
v_isSharedCheck_23_ = !lean_is_exclusive(v_head_5_);
if (v_isSharedCheck_23_ == 0)
{
v___x_13_ = v_head_5_;
v_isShared_14_ = v_isSharedCheck_23_;
goto v_resetjp_12_;
}
else
{
lean_inc(v_snd_11_);
lean_inc(v_fst_10_);
lean_dec(v_head_5_);
v___x_13_ = lean_box(0);
v_isShared_14_ = v_isSharedCheck_23_;
goto v_resetjp_12_;
}
v_resetjp_12_:
{
lean_object* v___x_15_; lean_object* v___x_17_; 
lean_inc(v_snd_1_);
v___x_15_ = l_List_appendTR___redArg(v_snd_1_, v_snd_11_);
if (v_isShared_14_ == 0)
{
lean_ctor_set(v___x_13_, 1, v___x_15_);
v___x_17_ = v___x_13_;
goto v_reusejp_16_;
}
else
{
lean_object* v_reuseFailAlloc_22_; 
v_reuseFailAlloc_22_ = lean_alloc_ctor(0, 2, 0);
lean_ctor_set(v_reuseFailAlloc_22_, 0, v_fst_10_);
lean_ctor_set(v_reuseFailAlloc_22_, 1, v___x_15_);
v___x_17_ = v_reuseFailAlloc_22_;
goto v_reusejp_16_;
}
v_reusejp_16_:
{
lean_object* v___x_19_; 
if (v_isShared_9_ == 0)
{
lean_ctor_set(v___x_8_, 1, v_a_3_);
lean_ctor_set(v___x_8_, 0, v___x_17_);
v___x_19_ = v___x_8_;
goto v_reusejp_18_;
}
else
{
lean_object* v_reuseFailAlloc_21_; 
v_reuseFailAlloc_21_ = lean_alloc_ctor(1, 2, 0);
lean_ctor_set(v_reuseFailAlloc_21_, 0, v___x_17_);
lean_ctor_set(v_reuseFailAlloc_21_, 1, v_a_3_);
v___x_19_ = v_reuseFailAlloc_21_;
goto v_reusejp_18_;
}
v_reusejp_18_:
{
v_a_2_ = v_tail_6_;
v_a_3_ = v___x_19_;
goto _start;
}
}
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1___redArg(lean_object* v_ndfst_25_, uint8_t v_symbol_26_, lean_object* v_a_27_, lean_object* v_a_28_){
_start:
{
if (lean_obj_tag(v_a_27_) == 0)
{
lean_object* v___x_29_; 
lean_dec_ref(v_ndfst_25_);
v___x_29_ = lean_array_to_list(v_a_28_);
return v___x_29_;
}
else
{
lean_object* v_head_30_; lean_object* v_tail_31_; lean_object* v_fst_32_; lean_object* v_snd_33_; lean_object* v_transition_34_; lean_object* v___x_35_; lean_object* v___x_36_; lean_object* v___x_37_; lean_object* v___x_38_; lean_object* v___x_39_; 
v_head_30_ = lean_ctor_get(v_a_27_, 0);
lean_inc(v_head_30_);
v_tail_31_ = lean_ctor_get(v_a_27_, 1);
lean_inc(v_tail_31_);
lean_dec_ref(v_a_27_);
v_fst_32_ = lean_ctor_get(v_head_30_, 0);
lean_inc(v_fst_32_);
v_snd_33_ = lean_ctor_get(v_head_30_, 1);
lean_inc(v_snd_33_);
lean_dec(v_head_30_);
v_transition_34_ = lean_ctor_get(v_ndfst_25_, 0);
v___x_35_ = lean_box(v_symbol_26_);
lean_inc_ref(v_transition_34_);
v___x_36_ = lean_apply_2(v_transition_34_, v_fst_32_, v___x_35_);
v___x_37_ = lean_box(0);
v___x_38_ = lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_ndfstStep_spec__0___redArg(v_snd_33_, v___x_36_, v___x_37_);
v___x_39_ = l_List_foldl___at___00Array_appendList_spec__0___redArg(v_a_28_, v___x_38_);
v_a_27_ = v_tail_31_;
v_a_28_ = v___x_39_;
goto _start;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1___redArg___boxed(lean_object* v_ndfst_41_, lean_object* v_symbol_42_, lean_object* v_a_43_, lean_object* v_a_44_){
_start:
{
uint8_t v_symbol_boxed_45_; lean_object* v_res_46_; 
v_symbol_boxed_45_ = lean_unbox(v_symbol_42_);
v_res_46_ = lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1___redArg(v_ndfst_41_, v_symbol_boxed_45_, v_a_43_, v_a_44_);
return v_res_46_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstStep___redArg(lean_object* v_ndfst_49_, lean_object* v_current_50_, uint8_t v_symbol_51_){
_start:
{
lean_object* v___x_52_; lean_object* v___x_53_; 
v___x_52_ = ((lean_object*)(lp_BioCompiler_BioCompiler_ndfstStep___redArg___closed__0));
v___x_53_ = lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1___redArg(v_ndfst_49_, v_symbol_51_, v_current_50_, v___x_52_);
return v___x_53_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstStep___redArg___boxed(lean_object* v_ndfst_54_, lean_object* v_current_55_, lean_object* v_symbol_56_){
_start:
{
uint8_t v_symbol_boxed_57_; lean_object* v_res_58_; 
v_symbol_boxed_57_ = lean_unbox(v_symbol_56_);
v_res_58_ = lp_BioCompiler_BioCompiler_ndfstStep___redArg(v_ndfst_54_, v_current_55_, v_symbol_boxed_57_);
return v_res_58_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstStep(lean_object* v_State_59_, lean_object* v_ndfst_60_, lean_object* v_current_61_, uint8_t v_symbol_62_){
_start:
{
lean_object* v___x_63_; 
v___x_63_ = lp_BioCompiler_BioCompiler_ndfstStep___redArg(v_ndfst_60_, v_current_61_, v_symbol_62_);
return v___x_63_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstStep___boxed(lean_object* v_State_64_, lean_object* v_ndfst_65_, lean_object* v_current_66_, lean_object* v_symbol_67_){
_start:
{
uint8_t v_symbol_boxed_68_; lean_object* v_res_69_; 
v_symbol_boxed_68_ = lean_unbox(v_symbol_67_);
v_res_69_ = lp_BioCompiler_BioCompiler_ndfstStep(v_State_64_, v_ndfst_65_, v_current_66_, v_symbol_boxed_68_);
return v_res_69_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_ndfstStep_spec__0(lean_object* v_State_70_, lean_object* v_snd_71_, lean_object* v_a_72_, lean_object* v_a_73_){
_start:
{
lean_object* v___x_74_; 
v___x_74_ = lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_ndfstStep_spec__0___redArg(v_snd_71_, v_a_72_, v_a_73_);
return v___x_74_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1(lean_object* v_State_75_, lean_object* v_ndfst_76_, uint8_t v_symbol_77_, lean_object* v_a_78_, lean_object* v_a_79_){
_start:
{
lean_object* v___x_80_; 
v___x_80_ = lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1___redArg(v_ndfst_76_, v_symbol_77_, v_a_78_, v_a_79_);
return v___x_80_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1___boxed(lean_object* v_State_81_, lean_object* v_ndfst_82_, lean_object* v_symbol_83_, lean_object* v_a_84_, lean_object* v_a_85_){
_start:
{
uint8_t v_symbol_boxed_86_; lean_object* v_res_87_; 
v_symbol_boxed_86_ = lean_unbox(v_symbol_83_);
v_res_87_ = lp_BioCompiler___private_Init_Data_List_Impl_0__List_flatMapTR_go___at___00BioCompiler_ndfstStep_spec__1(v_State_81_, v_ndfst_82_, v_symbol_boxed_86_, v_a_84_, v_a_85_);
return v_res_87_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0___redArg(lean_object* v_ndfst_88_, lean_object* v_x_89_, lean_object* v_x_90_){
_start:
{
if (lean_obj_tag(v_x_90_) == 0)
{
lean_dec_ref(v_ndfst_88_);
return v_x_89_;
}
else
{
lean_object* v_head_91_; lean_object* v_tail_92_; uint8_t v___x_93_; lean_object* v___x_94_; 
v_head_91_ = lean_ctor_get(v_x_90_, 0);
v_tail_92_ = lean_ctor_get(v_x_90_, 1);
v___x_93_ = lean_unbox(v_head_91_);
lean_inc_ref(v_ndfst_88_);
v___x_94_ = lp_BioCompiler_BioCompiler_ndfstStep___redArg(v_ndfst_88_, v_x_89_, v___x_93_);
v_x_89_ = v___x_94_;
v_x_90_ = v_tail_92_;
goto _start;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0___redArg___boxed(lean_object* v_ndfst_96_, lean_object* v_x_97_, lean_object* v_x_98_){
_start:
{
lean_object* v_res_99_; 
v_res_99_ = lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0___redArg(v_ndfst_96_, v_x_97_, v_x_98_);
lean_dec(v_x_98_);
return v_res_99_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstRun___redArg(lean_object* v_ndfst_100_, lean_object* v_input_101_){
_start:
{
lean_object* v_initial_102_; lean_object* v___x_103_; lean_object* v___x_104_; lean_object* v___x_105_; lean_object* v___x_106_; 
v_initial_102_ = lean_ctor_get(v_ndfst_100_, 1);
v___x_103_ = lean_box(0);
lean_inc(v_initial_102_);
v___x_104_ = lean_alloc_ctor(0, 2, 0);
lean_ctor_set(v___x_104_, 0, v_initial_102_);
lean_ctor_set(v___x_104_, 1, v___x_103_);
v___x_105_ = lean_alloc_ctor(1, 2, 0);
lean_ctor_set(v___x_105_, 0, v___x_104_);
lean_ctor_set(v___x_105_, 1, v___x_103_);
v___x_106_ = lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0___redArg(v_ndfst_100_, v___x_105_, v_input_101_);
return v___x_106_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstRun___redArg___boxed(lean_object* v_ndfst_107_, lean_object* v_input_108_){
_start:
{
lean_object* v_res_109_; 
v_res_109_ = lp_BioCompiler_BioCompiler_ndfstRun___redArg(v_ndfst_107_, v_input_108_);
lean_dec(v_input_108_);
return v_res_109_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstRun(lean_object* v_State_110_, lean_object* v_ndfst_111_, lean_object* v_input_112_){
_start:
{
lean_object* v___x_113_; 
v___x_113_ = lp_BioCompiler_BioCompiler_ndfstRun___redArg(v_ndfst_111_, v_input_112_);
return v___x_113_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstRun___boxed(lean_object* v_State_114_, lean_object* v_ndfst_115_, lean_object* v_input_116_){
_start:
{
lean_object* v_res_117_; 
v_res_117_ = lp_BioCompiler_BioCompiler_ndfstRun(v_State_114_, v_ndfst_115_, v_input_116_);
lean_dec(v_input_116_);
return v_res_117_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0(lean_object* v_State_118_, lean_object* v_ndfst_119_, lean_object* v_x_120_, lean_object* v_x_121_){
_start:
{
lean_object* v___x_122_; 
v___x_122_ = lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0___redArg(v_ndfst_119_, v_x_120_, v_x_121_);
return v___x_122_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0___boxed(lean_object* v_State_123_, lean_object* v_ndfst_124_, lean_object* v_x_125_, lean_object* v_x_126_){
_start:
{
lean_object* v_res_127_; 
v_res_127_ = lp_BioCompiler_List_foldl___at___00BioCompiler_ndfstRun_spec__0(v_State_123_, v_ndfst_124_, v_x_125_, v_x_126_);
lean_dec(v_x_126_);
return v_res_127_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_filterMapTR_go___at___00BioCompiler_ndfstOutputSet_spec__0___redArg(lean_object* v_ndfst_128_, lean_object* v_a_129_, lean_object* v_a_130_){
_start:
{
if (lean_obj_tag(v_a_129_) == 0)
{
lean_object* v___x_131_; 
lean_dec_ref(v_ndfst_128_);
v___x_131_ = lean_array_to_list(v_a_130_);
return v___x_131_;
}
else
{
lean_object* v_head_132_; lean_object* v_tail_133_; lean_object* v_fst_134_; lean_object* v_snd_135_; lean_object* v_accepting_136_; lean_object* v___x_137_; uint8_t v___x_138_; 
v_head_132_ = lean_ctor_get(v_a_129_, 0);
lean_inc(v_head_132_);
v_tail_133_ = lean_ctor_get(v_a_129_, 1);
lean_inc(v_tail_133_);
lean_dec_ref(v_a_129_);
v_fst_134_ = lean_ctor_get(v_head_132_, 0);
lean_inc(v_fst_134_);
v_snd_135_ = lean_ctor_get(v_head_132_, 1);
lean_inc(v_snd_135_);
lean_dec(v_head_132_);
v_accepting_136_ = lean_ctor_get(v_ndfst_128_, 2);
lean_inc_ref(v_accepting_136_);
v___x_137_ = lean_apply_1(v_accepting_136_, v_fst_134_);
v___x_138_ = lean_unbox(v___x_137_);
if (v___x_138_ == 0)
{
lean_dec(v_snd_135_);
v_a_129_ = v_tail_133_;
goto _start;
}
else
{
lean_object* v___x_140_; 
v___x_140_ = lean_array_push(v_a_130_, v_snd_135_);
v_a_129_ = v_tail_133_;
v_a_130_ = v___x_140_;
goto _start;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg(lean_object* v_ndfst_144_, lean_object* v_input_145_){
_start:
{
lean_object* v___x_146_; lean_object* v___x_147_; lean_object* v___x_148_; 
lean_inc_ref(v_ndfst_144_);
v___x_146_ = lp_BioCompiler_BioCompiler_ndfstRun___redArg(v_ndfst_144_, v_input_145_);
v___x_147_ = ((lean_object*)(lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg___closed__0));
v___x_148_ = lp_BioCompiler_List_filterMapTR_go___at___00BioCompiler_ndfstOutputSet_spec__0___redArg(v_ndfst_144_, v___x_146_, v___x_147_);
return v___x_148_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg___boxed(lean_object* v_ndfst_149_, lean_object* v_input_150_){
_start:
{
lean_object* v_res_151_; 
v_res_151_ = lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg(v_ndfst_149_, v_input_150_);
lean_dec(v_input_150_);
return v_res_151_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstOutputSet(lean_object* v_State_152_, lean_object* v_ndfst_153_, lean_object* v_input_154_){
_start:
{
lean_object* v___x_155_; 
v___x_155_ = lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg(v_ndfst_153_, v_input_154_);
return v___x_155_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstOutputSet___boxed(lean_object* v_State_156_, lean_object* v_ndfst_157_, lean_object* v_input_158_){
_start:
{
lean_object* v_res_159_; 
v_res_159_ = lp_BioCompiler_BioCompiler_ndfstOutputSet(v_State_156_, v_ndfst_157_, v_input_158_);
lean_dec(v_input_158_);
return v_res_159_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_filterMapTR_go___at___00BioCompiler_ndfstOutputSet_spec__0(lean_object* v_State_160_, lean_object* v_ndfst_161_, lean_object* v_a_162_, lean_object* v_a_163_){
_start:
{
lean_object* v___x_164_; 
v___x_164_ = lp_BioCompiler_List_filterMapTR_go___at___00BioCompiler_ndfstOutputSet_spec__0___redArg(v_ndfst_161_, v_a_162_, v_a_163_);
return v___x_164_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_List_beq___at___00List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0_spec__0(lean_object* v_x_165_, lean_object* v_x_166_){
_start:
{
if (lean_obj_tag(v_x_165_) == 0)
{
if (lean_obj_tag(v_x_166_) == 0)
{
uint8_t v___x_167_; 
v___x_167_ = 1;
return v___x_167_;
}
else
{
uint8_t v___x_168_; 
v___x_168_ = 0;
return v___x_168_;
}
}
else
{
if (lean_obj_tag(v_x_166_) == 0)
{
uint8_t v___x_169_; 
v___x_169_ = 0;
return v___x_169_;
}
else
{
lean_object* v_head_170_; lean_object* v_tail_171_; lean_object* v_head_172_; lean_object* v_tail_173_; uint8_t v___x_174_; uint8_t v___x_175_; uint8_t v___x_176_; 
v_head_170_ = lean_ctor_get(v_x_165_, 0);
v_tail_171_ = lean_ctor_get(v_x_165_, 1);
v_head_172_ = lean_ctor_get(v_x_166_, 0);
v_tail_173_ = lean_ctor_get(v_x_166_, 1);
v___x_174_ = lean_unbox(v_head_170_);
v___x_175_ = lean_unbox(v_head_172_);
v___x_176_ = lp_BioCompiler_BioCompiler_instBEqNucleotide_beq(v___x_174_, v___x_175_);
if (v___x_176_ == 0)
{
return v___x_176_;
}
else
{
v_x_165_ = v_tail_171_;
v_x_166_ = v_tail_173_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_beq___at___00List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0_spec__0___boxed(lean_object* v_x_178_, lean_object* v_x_179_){
_start:
{
uint8_t v_res_180_; lean_object* v_r_181_; 
v_res_180_ = lp_BioCompiler_List_beq___at___00List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0_spec__0(v_x_178_, v_x_179_);
lean_dec(v_x_179_);
lean_dec(v_x_178_);
v_r_181_ = lean_box(v_res_180_);
return v_r_181_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0(lean_object* v_as_183_){
_start:
{
lean_object* v___f_184_; lean_object* v___x_185_; 
v___f_184_ = ((lean_object*)(lp_BioCompiler_List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0___closed__0));
v___x_185_ = l_List_eraseDupsBy___redArg(v___f_184_, v_as_183_);
return v___x_185_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg(lean_object* v_ndfst_186_, lean_object* v_input_187_){
_start:
{
lean_object* v___x_188_; lean_object* v___x_189_; 
v___x_188_ = lp_BioCompiler_BioCompiler_ndfstOutputSet___redArg(v_ndfst_186_, v_input_187_);
v___x_189_ = lp_BioCompiler_List_eraseDups___at___00BioCompiler_ndfstUniqueOutputSet_spec__0(v___x_188_);
return v___x_189_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg___boxed(lean_object* v_ndfst_190_, lean_object* v_input_191_){
_start:
{
lean_object* v_res_192_; 
v_res_192_ = lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg(v_ndfst_190_, v_input_191_);
lean_dec(v_input_191_);
return v_res_192_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet(lean_object* v_State_193_, lean_object* v_ndfst_194_, lean_object* v_input_195_){
_start:
{
lean_object* v___x_196_; 
v___x_196_ = lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___redArg(v_ndfst_194_, v_input_195_);
return v___x_196_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet___boxed(lean_object* v_State_197_, lean_object* v_ndfst_198_, lean_object* v_input_199_){
_start:
{
lean_object* v_res_200_; 
v_res_200_ = lp_BioCompiler_BioCompiler_ndfstUniqueOutputSet(v_State_197_, v_ndfst_198_, v_input_199_);
lean_dec(v_input_199_);
return v_res_200_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__7(void){
_start:
{
lean_object* v___x_214_; lean_object* v___x_215_; 
v___x_214_ = lean_unsigned_to_nat(12u);
v___x_215_ = lean_nat_to_int(v___x_214_);
return v___x_215_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__12(void){
_start:
{
lean_object* v___x_222_; lean_object* v___x_223_; 
v___x_222_ = lean_unsigned_to_nat(16u);
v___x_223_ = lean_nat_to_int(v___x_222_);
return v___x_223_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__14(void){
_start:
{
lean_object* v___x_225_; lean_object* v___x_226_; 
v___x_225_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__6));
v___x_226_ = lean_string_length(v___x_225_);
return v___x_226_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__15(void){
_start:
{
lean_object* v___x_227_; lean_object* v___x_228_; 
v___x_227_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__14, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__14_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__14);
v___x_228_ = lean_nat_to_int(v___x_227_);
return v___x_228_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22(void){
_start:
{
lean_object* v___x_238_; lean_object* v___x_239_; 
v___x_238_ = lean_unsigned_to_nat(0u);
v___x_239_ = lean_nat_to_int(v___x_238_);
return v___x_239_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg(lean_object* v_x_246_){
_start:
{
lean_object* v_eseThreshold_247_; lean_object* v_cellType_248_; lean_object* v_essThreshold_249_; lean_object* v_iseThreshold_250_; lean_object* v_issThreshold_251_; lean_object* v_num_252_; lean_object* v_den_253_; lean_object* v___x_255_; uint8_t v_isShared_256_; uint8_t v_isSharedCheck_417_; 
v_eseThreshold_247_ = lean_ctor_get(v_x_246_, 1);
lean_inc_ref(v_eseThreshold_247_);
v_cellType_248_ = lean_ctor_get(v_x_246_, 0);
lean_inc_ref(v_cellType_248_);
v_essThreshold_249_ = lean_ctor_get(v_x_246_, 2);
lean_inc_ref(v_essThreshold_249_);
v_iseThreshold_250_ = lean_ctor_get(v_x_246_, 3);
lean_inc_ref(v_iseThreshold_250_);
v_issThreshold_251_ = lean_ctor_get(v_x_246_, 4);
lean_inc_ref(v_issThreshold_251_);
lean_dec_ref(v_x_246_);
v_num_252_ = lean_ctor_get(v_eseThreshold_247_, 0);
v_den_253_ = lean_ctor_get(v_eseThreshold_247_, 1);
v_isSharedCheck_417_ = !lean_is_exclusive(v_eseThreshold_247_);
if (v_isSharedCheck_417_ == 0)
{
v___x_255_ = v_eseThreshold_247_;
v_isShared_256_ = v_isSharedCheck_417_;
goto v_resetjp_254_;
}
else
{
lean_inc(v_den_253_);
lean_inc(v_num_252_);
lean_dec(v_eseThreshold_247_);
v___x_255_ = lean_box(0);
v_isShared_256_ = v_isSharedCheck_417_;
goto v_resetjp_254_;
}
v_resetjp_254_:
{
lean_object* v___x_257_; lean_object* v___x_258_; lean_object* v___x_259_; lean_object* v___x_260_; lean_object* v___x_261_; lean_object* v___x_263_; 
v___x_257_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__4));
v___x_258_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__5));
v___x_259_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__7);
v___x_260_ = l_String_quote(v_cellType_248_);
v___x_261_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_261_, 0, v___x_260_);
if (v_isShared_256_ == 0)
{
lean_ctor_set_tag(v___x_255_, 4);
lean_ctor_set(v___x_255_, 1, v___x_261_);
lean_ctor_set(v___x_255_, 0, v___x_259_);
v___x_263_ = v___x_255_;
goto v_reusejp_262_;
}
else
{
lean_object* v_reuseFailAlloc_416_; 
v_reuseFailAlloc_416_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_416_, 0, v___x_259_);
lean_ctor_set(v_reuseFailAlloc_416_, 1, v___x_261_);
v___x_263_ = v_reuseFailAlloc_416_;
goto v_reusejp_262_;
}
v_reusejp_262_:
{
uint8_t v___x_264_; lean_object* v___x_265_; lean_object* v___x_266_; lean_object* v___x_267_; lean_object* v___x_268_; lean_object* v___x_269_; lean_object* v___x_270_; lean_object* v___x_271_; lean_object* v___x_272_; lean_object* v___x_273_; lean_object* v___x_274_; lean_object* v___y_276_; lean_object* v___y_277_; lean_object* v___y_289_; lean_object* v___y_290_; lean_object* v___y_326_; lean_object* v___y_327_; lean_object* v___y_363_; lean_object* v___x_398_; uint8_t v___x_399_; 
v___x_264_ = 0;
v___x_265_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_265_, 0, v___x_263_);
lean_ctor_set_uint8(v___x_265_, sizeof(void*)*1, v___x_264_);
v___x_266_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_266_, 0, v___x_258_);
lean_ctor_set(v___x_266_, 1, v___x_265_);
v___x_267_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__9));
v___x_268_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_268_, 0, v___x_266_);
lean_ctor_set(v___x_268_, 1, v___x_267_);
v___x_269_ = lean_box(1);
v___x_270_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_270_, 0, v___x_268_);
lean_ctor_set(v___x_270_, 1, v___x_269_);
v___x_271_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__11));
v___x_272_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_272_, 0, v___x_270_);
lean_ctor_set(v___x_272_, 1, v___x_271_);
v___x_273_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_273_, 0, v___x_272_);
lean_ctor_set(v___x_273_, 1, v___x_257_);
v___x_274_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__12, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__12_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__12);
v___x_398_ = lean_unsigned_to_nat(1u);
v___x_399_ = lean_nat_dec_eq(v_den_253_, v___x_398_);
if (v___x_399_ == 0)
{
lean_object* v___x_400_; lean_object* v___x_401_; lean_object* v___x_402_; lean_object* v___x_403_; lean_object* v___x_404_; lean_object* v___x_405_; lean_object* v___x_406_; lean_object* v___x_407_; 
v___x_400_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__20));
v___x_401_ = l_Int_repr(v_num_252_);
lean_dec(v_num_252_);
v___x_402_ = lean_string_append(v___x_400_, v___x_401_);
lean_dec_ref(v___x_401_);
v___x_403_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__21));
v___x_404_ = lean_string_append(v___x_402_, v___x_403_);
v___x_405_ = l_Nat_reprFast(v_den_253_);
v___x_406_ = lean_string_append(v___x_404_, v___x_405_);
lean_dec_ref(v___x_405_);
v___x_407_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_407_, 0, v___x_406_);
v___y_363_ = v___x_407_;
goto v___jp_362_;
}
else
{
lean_object* v___x_408_; lean_object* v___x_409_; uint8_t v___x_410_; 
lean_dec(v_den_253_);
v___x_408_ = lean_unsigned_to_nat(0u);
v___x_409_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22);
v___x_410_ = lean_int_dec_lt(v_num_252_, v___x_409_);
if (v___x_410_ == 0)
{
lean_object* v___x_411_; lean_object* v___x_412_; 
v___x_411_ = l_Int_repr(v_num_252_);
lean_dec(v_num_252_);
v___x_412_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_412_, 0, v___x_411_);
v___y_363_ = v___x_412_;
goto v___jp_362_;
}
else
{
lean_object* v___x_413_; lean_object* v___x_414_; lean_object* v___x_415_; 
v___x_413_ = l_Int_repr(v_num_252_);
lean_dec(v_num_252_);
v___x_414_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_414_, 0, v___x_413_);
v___x_415_ = l_Repr_addAppParen(v___x_414_, v___x_408_);
v___y_363_ = v___x_415_;
goto v___jp_362_;
}
}
v___jp_275_:
{
lean_object* v___x_278_; lean_object* v___x_279_; lean_object* v___x_280_; lean_object* v___x_281_; lean_object* v___x_282_; lean_object* v___x_283_; lean_object* v___x_284_; lean_object* v___x_285_; lean_object* v___x_286_; lean_object* v___x_287_; 
v___x_278_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_278_, 0, v___x_274_);
lean_ctor_set(v___x_278_, 1, v___y_277_);
v___x_279_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_279_, 0, v___x_278_);
lean_ctor_set_uint8(v___x_279_, sizeof(void*)*1, v___x_264_);
v___x_280_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_280_, 0, v___y_276_);
lean_ctor_set(v___x_280_, 1, v___x_279_);
v___x_281_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__15, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__15_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__15);
v___x_282_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__16));
v___x_283_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_283_, 0, v___x_282_);
lean_ctor_set(v___x_283_, 1, v___x_280_);
v___x_284_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__17));
v___x_285_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_285_, 0, v___x_283_);
lean_ctor_set(v___x_285_, 1, v___x_284_);
v___x_286_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_286_, 0, v___x_281_);
lean_ctor_set(v___x_286_, 1, v___x_285_);
v___x_287_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_287_, 0, v___x_286_);
lean_ctor_set_uint8(v___x_287_, sizeof(void*)*1, v___x_264_);
return v___x_287_;
}
v___jp_288_:
{
lean_object* v_num_291_; lean_object* v_den_292_; lean_object* v___x_294_; uint8_t v_isShared_295_; uint8_t v_isSharedCheck_324_; 
v_num_291_ = lean_ctor_get(v_issThreshold_251_, 0);
v_den_292_ = lean_ctor_get(v_issThreshold_251_, 1);
v_isSharedCheck_324_ = !lean_is_exclusive(v_issThreshold_251_);
if (v_isSharedCheck_324_ == 0)
{
v___x_294_ = v_issThreshold_251_;
v_isShared_295_ = v_isSharedCheck_324_;
goto v_resetjp_293_;
}
else
{
lean_inc(v_den_292_);
lean_inc(v_num_291_);
lean_dec(v_issThreshold_251_);
v___x_294_ = lean_box(0);
v_isShared_295_ = v_isSharedCheck_324_;
goto v_resetjp_293_;
}
v_resetjp_293_:
{
lean_object* v___x_297_; 
if (v_isShared_295_ == 0)
{
lean_ctor_set_tag(v___x_294_, 4);
lean_ctor_set(v___x_294_, 1, v___y_290_);
lean_ctor_set(v___x_294_, 0, v___x_274_);
v___x_297_ = v___x_294_;
goto v_reusejp_296_;
}
else
{
lean_object* v_reuseFailAlloc_323_; 
v_reuseFailAlloc_323_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_323_, 0, v___x_274_);
lean_ctor_set(v_reuseFailAlloc_323_, 1, v___y_290_);
v___x_297_ = v_reuseFailAlloc_323_;
goto v_reusejp_296_;
}
v_reusejp_296_:
{
lean_object* v___x_298_; lean_object* v___x_299_; lean_object* v___x_300_; lean_object* v___x_301_; lean_object* v___x_302_; lean_object* v___x_303_; lean_object* v___x_304_; lean_object* v___x_305_; uint8_t v___x_306_; 
v___x_298_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_298_, 0, v___x_297_);
lean_ctor_set_uint8(v___x_298_, sizeof(void*)*1, v___x_264_);
v___x_299_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_299_, 0, v___y_289_);
lean_ctor_set(v___x_299_, 1, v___x_298_);
v___x_300_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_300_, 0, v___x_299_);
lean_ctor_set(v___x_300_, 1, v___x_267_);
v___x_301_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_301_, 0, v___x_300_);
lean_ctor_set(v___x_301_, 1, v___x_269_);
v___x_302_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__19));
v___x_303_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_303_, 0, v___x_301_);
lean_ctor_set(v___x_303_, 1, v___x_302_);
v___x_304_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_304_, 0, v___x_303_);
lean_ctor_set(v___x_304_, 1, v___x_257_);
v___x_305_ = lean_unsigned_to_nat(1u);
v___x_306_ = lean_nat_dec_eq(v_den_292_, v___x_305_);
if (v___x_306_ == 0)
{
lean_object* v___x_307_; lean_object* v___x_308_; lean_object* v___x_309_; lean_object* v___x_310_; lean_object* v___x_311_; lean_object* v___x_312_; lean_object* v___x_313_; lean_object* v___x_314_; 
v___x_307_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__20));
v___x_308_ = l_Int_repr(v_num_291_);
lean_dec(v_num_291_);
v___x_309_ = lean_string_append(v___x_307_, v___x_308_);
lean_dec_ref(v___x_308_);
v___x_310_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__21));
v___x_311_ = lean_string_append(v___x_309_, v___x_310_);
v___x_312_ = l_Nat_reprFast(v_den_292_);
v___x_313_ = lean_string_append(v___x_311_, v___x_312_);
lean_dec_ref(v___x_312_);
v___x_314_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_314_, 0, v___x_313_);
v___y_276_ = v___x_304_;
v___y_277_ = v___x_314_;
goto v___jp_275_;
}
else
{
lean_object* v___x_315_; lean_object* v___x_316_; uint8_t v___x_317_; 
lean_dec(v_den_292_);
v___x_315_ = lean_unsigned_to_nat(0u);
v___x_316_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22);
v___x_317_ = lean_int_dec_lt(v_num_291_, v___x_316_);
if (v___x_317_ == 0)
{
lean_object* v___x_318_; lean_object* v___x_319_; 
v___x_318_ = l_Int_repr(v_num_291_);
lean_dec(v_num_291_);
v___x_319_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_319_, 0, v___x_318_);
v___y_276_ = v___x_304_;
v___y_277_ = v___x_319_;
goto v___jp_275_;
}
else
{
lean_object* v___x_320_; lean_object* v___x_321_; lean_object* v___x_322_; 
v___x_320_ = l_Int_repr(v_num_291_);
lean_dec(v_num_291_);
v___x_321_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_321_, 0, v___x_320_);
v___x_322_ = l_Repr_addAppParen(v___x_321_, v___x_315_);
v___y_276_ = v___x_304_;
v___y_277_ = v___x_322_;
goto v___jp_275_;
}
}
}
}
}
v___jp_325_:
{
lean_object* v_num_328_; lean_object* v_den_329_; lean_object* v___x_331_; uint8_t v_isShared_332_; uint8_t v_isSharedCheck_361_; 
v_num_328_ = lean_ctor_get(v_iseThreshold_250_, 0);
v_den_329_ = lean_ctor_get(v_iseThreshold_250_, 1);
v_isSharedCheck_361_ = !lean_is_exclusive(v_iseThreshold_250_);
if (v_isSharedCheck_361_ == 0)
{
v___x_331_ = v_iseThreshold_250_;
v_isShared_332_ = v_isSharedCheck_361_;
goto v_resetjp_330_;
}
else
{
lean_inc(v_den_329_);
lean_inc(v_num_328_);
lean_dec(v_iseThreshold_250_);
v___x_331_ = lean_box(0);
v_isShared_332_ = v_isSharedCheck_361_;
goto v_resetjp_330_;
}
v_resetjp_330_:
{
lean_object* v___x_334_; 
if (v_isShared_332_ == 0)
{
lean_ctor_set_tag(v___x_331_, 4);
lean_ctor_set(v___x_331_, 1, v___y_327_);
lean_ctor_set(v___x_331_, 0, v___x_274_);
v___x_334_ = v___x_331_;
goto v_reusejp_333_;
}
else
{
lean_object* v_reuseFailAlloc_360_; 
v_reuseFailAlloc_360_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_360_, 0, v___x_274_);
lean_ctor_set(v_reuseFailAlloc_360_, 1, v___y_327_);
v___x_334_ = v_reuseFailAlloc_360_;
goto v_reusejp_333_;
}
v_reusejp_333_:
{
lean_object* v___x_335_; lean_object* v___x_336_; lean_object* v___x_337_; lean_object* v___x_338_; lean_object* v___x_339_; lean_object* v___x_340_; lean_object* v___x_341_; lean_object* v___x_342_; uint8_t v___x_343_; 
v___x_335_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_335_, 0, v___x_334_);
lean_ctor_set_uint8(v___x_335_, sizeof(void*)*1, v___x_264_);
v___x_336_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_336_, 0, v___y_326_);
lean_ctor_set(v___x_336_, 1, v___x_335_);
v___x_337_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_337_, 0, v___x_336_);
lean_ctor_set(v___x_337_, 1, v___x_267_);
v___x_338_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_338_, 0, v___x_337_);
lean_ctor_set(v___x_338_, 1, v___x_269_);
v___x_339_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__24));
v___x_340_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_340_, 0, v___x_338_);
lean_ctor_set(v___x_340_, 1, v___x_339_);
v___x_341_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_341_, 0, v___x_340_);
lean_ctor_set(v___x_341_, 1, v___x_257_);
v___x_342_ = lean_unsigned_to_nat(1u);
v___x_343_ = lean_nat_dec_eq(v_den_329_, v___x_342_);
if (v___x_343_ == 0)
{
lean_object* v___x_344_; lean_object* v___x_345_; lean_object* v___x_346_; lean_object* v___x_347_; lean_object* v___x_348_; lean_object* v___x_349_; lean_object* v___x_350_; lean_object* v___x_351_; 
v___x_344_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__20));
v___x_345_ = l_Int_repr(v_num_328_);
lean_dec(v_num_328_);
v___x_346_ = lean_string_append(v___x_344_, v___x_345_);
lean_dec_ref(v___x_345_);
v___x_347_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__21));
v___x_348_ = lean_string_append(v___x_346_, v___x_347_);
v___x_349_ = l_Nat_reprFast(v_den_329_);
v___x_350_ = lean_string_append(v___x_348_, v___x_349_);
lean_dec_ref(v___x_349_);
v___x_351_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_351_, 0, v___x_350_);
v___y_289_ = v___x_341_;
v___y_290_ = v___x_351_;
goto v___jp_288_;
}
else
{
lean_object* v___x_352_; lean_object* v___x_353_; uint8_t v___x_354_; 
lean_dec(v_den_329_);
v___x_352_ = lean_unsigned_to_nat(0u);
v___x_353_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22);
v___x_354_ = lean_int_dec_lt(v_num_328_, v___x_353_);
if (v___x_354_ == 0)
{
lean_object* v___x_355_; lean_object* v___x_356_; 
v___x_355_ = l_Int_repr(v_num_328_);
lean_dec(v_num_328_);
v___x_356_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_356_, 0, v___x_355_);
v___y_289_ = v___x_341_;
v___y_290_ = v___x_356_;
goto v___jp_288_;
}
else
{
lean_object* v___x_357_; lean_object* v___x_358_; lean_object* v___x_359_; 
v___x_357_ = l_Int_repr(v_num_328_);
lean_dec(v_num_328_);
v___x_358_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_358_, 0, v___x_357_);
v___x_359_ = l_Repr_addAppParen(v___x_358_, v___x_352_);
v___y_289_ = v___x_341_;
v___y_290_ = v___x_359_;
goto v___jp_288_;
}
}
}
}
}
v___jp_362_:
{
lean_object* v_num_364_; lean_object* v_den_365_; lean_object* v___x_367_; uint8_t v_isShared_368_; uint8_t v_isSharedCheck_397_; 
v_num_364_ = lean_ctor_get(v_essThreshold_249_, 0);
v_den_365_ = lean_ctor_get(v_essThreshold_249_, 1);
v_isSharedCheck_397_ = !lean_is_exclusive(v_essThreshold_249_);
if (v_isSharedCheck_397_ == 0)
{
v___x_367_ = v_essThreshold_249_;
v_isShared_368_ = v_isSharedCheck_397_;
goto v_resetjp_366_;
}
else
{
lean_inc(v_den_365_);
lean_inc(v_num_364_);
lean_dec(v_essThreshold_249_);
v___x_367_ = lean_box(0);
v_isShared_368_ = v_isSharedCheck_397_;
goto v_resetjp_366_;
}
v_resetjp_366_:
{
lean_object* v___x_370_; 
if (v_isShared_368_ == 0)
{
lean_ctor_set_tag(v___x_367_, 4);
lean_ctor_set(v___x_367_, 1, v___y_363_);
lean_ctor_set(v___x_367_, 0, v___x_274_);
v___x_370_ = v___x_367_;
goto v_reusejp_369_;
}
else
{
lean_object* v_reuseFailAlloc_396_; 
v_reuseFailAlloc_396_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_396_, 0, v___x_274_);
lean_ctor_set(v_reuseFailAlloc_396_, 1, v___y_363_);
v___x_370_ = v_reuseFailAlloc_396_;
goto v_reusejp_369_;
}
v_reusejp_369_:
{
lean_object* v___x_371_; lean_object* v___x_372_; lean_object* v___x_373_; lean_object* v___x_374_; lean_object* v___x_375_; lean_object* v___x_376_; lean_object* v___x_377_; lean_object* v___x_378_; uint8_t v___x_379_; 
v___x_371_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_371_, 0, v___x_370_);
lean_ctor_set_uint8(v___x_371_, sizeof(void*)*1, v___x_264_);
v___x_372_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_372_, 0, v___x_273_);
lean_ctor_set(v___x_372_, 1, v___x_371_);
v___x_373_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_373_, 0, v___x_372_);
lean_ctor_set(v___x_373_, 1, v___x_267_);
v___x_374_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_374_, 0, v___x_373_);
lean_ctor_set(v___x_374_, 1, v___x_269_);
v___x_375_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__26));
v___x_376_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_376_, 0, v___x_374_);
lean_ctor_set(v___x_376_, 1, v___x_375_);
v___x_377_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_377_, 0, v___x_376_);
lean_ctor_set(v___x_377_, 1, v___x_257_);
v___x_378_ = lean_unsigned_to_nat(1u);
v___x_379_ = lean_nat_dec_eq(v_den_365_, v___x_378_);
if (v___x_379_ == 0)
{
lean_object* v___x_380_; lean_object* v___x_381_; lean_object* v___x_382_; lean_object* v___x_383_; lean_object* v___x_384_; lean_object* v___x_385_; lean_object* v___x_386_; lean_object* v___x_387_; 
v___x_380_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__20));
v___x_381_ = l_Int_repr(v_num_364_);
lean_dec(v_num_364_);
v___x_382_ = lean_string_append(v___x_380_, v___x_381_);
lean_dec_ref(v___x_381_);
v___x_383_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__21));
v___x_384_ = lean_string_append(v___x_382_, v___x_383_);
v___x_385_ = l_Nat_reprFast(v_den_365_);
v___x_386_ = lean_string_append(v___x_384_, v___x_385_);
lean_dec_ref(v___x_385_);
v___x_387_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_387_, 0, v___x_386_);
v___y_326_ = v___x_377_;
v___y_327_ = v___x_387_;
goto v___jp_325_;
}
else
{
lean_object* v___x_388_; lean_object* v___x_389_; uint8_t v___x_390_; 
lean_dec(v_den_365_);
v___x_388_ = lean_unsigned_to_nat(0u);
v___x_389_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__22);
v___x_390_ = lean_int_dec_lt(v_num_364_, v___x_389_);
if (v___x_390_ == 0)
{
lean_object* v___x_391_; lean_object* v___x_392_; 
v___x_391_ = l_Int_repr(v_num_364_);
lean_dec(v_num_364_);
v___x_392_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_392_, 0, v___x_391_);
v___y_326_ = v___x_377_;
v___y_327_ = v___x_392_;
goto v___jp_325_;
}
else
{
lean_object* v___x_393_; lean_object* v___x_394_; lean_object* v___x_395_; 
v___x_393_ = l_Int_repr(v_num_364_);
lean_dec(v_num_364_);
v___x_394_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_394_, 0, v___x_393_);
v___x_395_ = l_Repr_addAppParen(v___x_394_, v___x_388_);
v___y_326_ = v___x_377_;
v___y_327_ = v___x_395_;
goto v___jp_325_;
}
}
}
}
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr(lean_object* v_x_418_, lean_object* v_prec_419_){
_start:
{
lean_object* v___x_420_; 
v___x_420_ = lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg(v_x_418_);
return v___x_420_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprCellularContext_repr___boxed(lean_object* v_x_421_, lean_object* v_prec_422_){
_start:
{
lean_object* v_res_423_; 
v_res_423_ = lp_BioCompiler_BioCompiler_instReprCellularContext_repr(v_x_421_, v_prec_422_);
lean_dec(v_prec_422_);
return v_res_423_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0_spec__1_spec__3(lean_object* v_x_426_, lean_object* v_x_427_, lean_object* v_x_428_){
_start:
{
if (lean_obj_tag(v_x_428_) == 0)
{
lean_dec(v_x_426_);
return v_x_427_;
}
else
{
lean_object* v_head_429_; lean_object* v_tail_430_; lean_object* v___x_432_; uint8_t v_isShared_433_; uint8_t v_isSharedCheck_442_; 
v_head_429_ = lean_ctor_get(v_x_428_, 0);
v_tail_430_ = lean_ctor_get(v_x_428_, 1);
v_isSharedCheck_442_ = !lean_is_exclusive(v_x_428_);
if (v_isSharedCheck_442_ == 0)
{
v___x_432_ = v_x_428_;
v_isShared_433_ = v_isSharedCheck_442_;
goto v_resetjp_431_;
}
else
{
lean_inc(v_tail_430_);
lean_inc(v_head_429_);
lean_dec(v_x_428_);
v___x_432_ = lean_box(0);
v_isShared_433_ = v_isSharedCheck_442_;
goto v_resetjp_431_;
}
v_resetjp_431_:
{
lean_object* v___x_435_; 
lean_inc(v_x_426_);
if (v_isShared_433_ == 0)
{
lean_ctor_set_tag(v___x_432_, 5);
lean_ctor_set(v___x_432_, 1, v_x_426_);
lean_ctor_set(v___x_432_, 0, v_x_427_);
v___x_435_ = v___x_432_;
goto v_reusejp_434_;
}
else
{
lean_object* v_reuseFailAlloc_441_; 
v_reuseFailAlloc_441_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_441_, 0, v_x_427_);
lean_ctor_set(v_reuseFailAlloc_441_, 1, v_x_426_);
v___x_435_ = v_reuseFailAlloc_441_;
goto v_reusejp_434_;
}
v_reusejp_434_:
{
lean_object* v___x_436_; uint8_t v___x_437_; lean_object* v___x_438_; lean_object* v___x_439_; 
v___x_436_ = lean_unsigned_to_nat(0u);
v___x_437_ = lean_unbox(v_head_429_);
lean_dec(v_head_429_);
v___x_438_ = lp_BioCompiler_BioCompiler_instReprNucleotide_repr(v___x_437_, v___x_436_);
v___x_439_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_439_, 0, v___x_435_);
lean_ctor_set(v___x_439_, 1, v___x_438_);
v_x_427_ = v___x_439_;
v_x_428_ = v_tail_430_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0_spec__1(lean_object* v_x_443_, lean_object* v_x_444_, lean_object* v_x_445_){
_start:
{
if (lean_obj_tag(v_x_445_) == 0)
{
lean_dec(v_x_443_);
return v_x_444_;
}
else
{
lean_object* v_head_446_; lean_object* v_tail_447_; lean_object* v___x_449_; uint8_t v_isShared_450_; uint8_t v_isSharedCheck_459_; 
v_head_446_ = lean_ctor_get(v_x_445_, 0);
v_tail_447_ = lean_ctor_get(v_x_445_, 1);
v_isSharedCheck_459_ = !lean_is_exclusive(v_x_445_);
if (v_isSharedCheck_459_ == 0)
{
v___x_449_ = v_x_445_;
v_isShared_450_ = v_isSharedCheck_459_;
goto v_resetjp_448_;
}
else
{
lean_inc(v_tail_447_);
lean_inc(v_head_446_);
lean_dec(v_x_445_);
v___x_449_ = lean_box(0);
v_isShared_450_ = v_isSharedCheck_459_;
goto v_resetjp_448_;
}
v_resetjp_448_:
{
lean_object* v___x_452_; 
lean_inc(v_x_443_);
if (v_isShared_450_ == 0)
{
lean_ctor_set_tag(v___x_449_, 5);
lean_ctor_set(v___x_449_, 1, v_x_443_);
lean_ctor_set(v___x_449_, 0, v_x_444_);
v___x_452_ = v___x_449_;
goto v_reusejp_451_;
}
else
{
lean_object* v_reuseFailAlloc_458_; 
v_reuseFailAlloc_458_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_458_, 0, v_x_444_);
lean_ctor_set(v_reuseFailAlloc_458_, 1, v_x_443_);
v___x_452_ = v_reuseFailAlloc_458_;
goto v_reusejp_451_;
}
v_reusejp_451_:
{
lean_object* v___x_453_; uint8_t v___x_454_; lean_object* v___x_455_; lean_object* v___x_456_; lean_object* v___x_457_; 
v___x_453_ = lean_unsigned_to_nat(0u);
v___x_454_ = lean_unbox(v_head_446_);
lean_dec(v_head_446_);
v___x_455_ = lp_BioCompiler_BioCompiler_instReprNucleotide_repr(v___x_454_, v___x_453_);
v___x_456_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_456_, 0, v___x_452_);
lean_ctor_set(v___x_456_, 1, v___x_455_);
v___x_457_ = lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0_spec__1_spec__3(v_x_443_, v___x_456_, v_tail_447_);
return v___x_457_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0___lam__0(uint8_t v___y_460_){
_start:
{
lean_object* v___x_461_; lean_object* v___x_462_; 
v___x_461_ = lean_unsigned_to_nat(0u);
v___x_462_ = lp_BioCompiler_BioCompiler_instReprNucleotide_repr(v___y_460_, v___x_461_);
return v___x_462_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0___lam__0___boxed(lean_object* v___y_463_){
_start:
{
uint8_t v___y_611__boxed_464_; lean_object* v_res_465_; 
v___y_611__boxed_464_ = lean_unbox(v___y_463_);
v_res_465_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0___lam__0(v___y_611__boxed_464_);
return v_res_465_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0(lean_object* v_x_466_, lean_object* v_x_467_){
_start:
{
if (lean_obj_tag(v_x_466_) == 0)
{
lean_object* v___x_468_; 
lean_dec(v_x_467_);
v___x_468_ = lean_box(0);
return v___x_468_;
}
else
{
lean_object* v_tail_469_; 
v_tail_469_ = lean_ctor_get(v_x_466_, 1);
if (lean_obj_tag(v_tail_469_) == 0)
{
lean_object* v_head_470_; uint8_t v___x_471_; lean_object* v___x_472_; 
lean_dec(v_x_467_);
v_head_470_ = lean_ctor_get(v_x_466_, 0);
lean_inc(v_head_470_);
lean_dec_ref(v_x_466_);
v___x_471_ = lean_unbox(v_head_470_);
lean_dec(v_head_470_);
v___x_472_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0___lam__0(v___x_471_);
return v___x_472_;
}
else
{
lean_object* v_head_473_; uint8_t v___x_474_; lean_object* v___x_475_; lean_object* v___x_476_; 
lean_inc(v_tail_469_);
v_head_473_ = lean_ctor_get(v_x_466_, 0);
lean_inc(v_head_473_);
lean_dec_ref(v_x_466_);
v___x_474_ = lean_unbox(v_head_473_);
lean_dec(v_head_473_);
v___x_475_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0___lam__0(v___x_474_);
v___x_476_ = lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0_spec__1(v_x_467_, v___x_475_, v_tail_469_);
return v___x_476_;
}
}
}
}
static lean_object* _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__5(void){
_start:
{
lean_object* v___x_485_; lean_object* v___x_486_; 
v___x_485_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__2));
v___x_486_ = lean_string_length(v___x_485_);
return v___x_486_;
}
}
static lean_object* _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__6(void){
_start:
{
lean_object* v___x_487_; lean_object* v___x_488_; 
v___x_487_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__5, &lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__5_once, _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__5);
v___x_488_ = lean_nat_to_int(v___x_487_);
return v___x_488_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(lean_object* v_a_493_){
_start:
{
if (lean_obj_tag(v_a_493_) == 0)
{
lean_object* v___x_494_; 
v___x_494_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__1));
return v___x_494_;
}
else
{
lean_object* v___x_495_; lean_object* v___x_496_; lean_object* v___x_497_; lean_object* v___x_498_; lean_object* v___x_499_; lean_object* v___x_500_; lean_object* v___x_501_; lean_object* v___x_502_; uint8_t v___x_503_; lean_object* v___x_504_; 
v___x_495_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__3));
v___x_496_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0_spec__0(v_a_493_, v___x_495_);
v___x_497_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__6, &lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__6_once, _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__6);
v___x_498_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__7));
v___x_499_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_499_, 0, v___x_498_);
lean_ctor_set(v___x_499_, 1, v___x_496_);
v___x_500_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__8));
v___x_501_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_501_, 0, v___x_499_);
lean_ctor_set(v___x_501_, 1, v___x_500_);
v___x_502_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_502_, 0, v___x_497_);
lean_ctor_set(v___x_502_, 1, v___x_501_);
v___x_503_ = 0;
v___x_504_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_504_, 0, v___x_502_);
lean_ctor_set_uint8(v___x_504_, sizeof(void*)*1, v___x_503_);
return v___x_504_;
}
}
}
static lean_object* _init_lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__1(void){
_start:
{
lean_object* v___x_506_; lean_object* v___x_507_; 
v___x_506_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__20));
v___x_507_ = lean_string_length(v___x_506_);
return v___x_507_;
}
}
static lean_object* _init_lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__2(void){
_start:
{
lean_object* v___x_508_; lean_object* v___x_509_; 
v___x_508_ = lean_obj_once(&lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__1, &lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__1_once, _init_lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__1);
v___x_509_ = lean_nat_to_int(v___x_508_);
return v___x_509_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg(lean_object* v_x_514_){
_start:
{
lean_object* v_fst_515_; lean_object* v_snd_516_; lean_object* v___x_518_; uint8_t v_isShared_519_; uint8_t v_isSharedCheck_540_; 
v_fst_515_ = lean_ctor_get(v_x_514_, 0);
v_snd_516_ = lean_ctor_get(v_x_514_, 1);
v_isSharedCheck_540_ = !lean_is_exclusive(v_x_514_);
if (v_isSharedCheck_540_ == 0)
{
v___x_518_ = v_x_514_;
v_isShared_519_ = v_isSharedCheck_540_;
goto v_resetjp_517_;
}
else
{
lean_inc(v_snd_516_);
lean_inc(v_fst_515_);
lean_dec(v_x_514_);
v___x_518_ = lean_box(0);
v_isShared_519_ = v_isSharedCheck_540_;
goto v_resetjp_517_;
}
v_resetjp_517_:
{
lean_object* v___x_520_; lean_object* v___x_521_; lean_object* v___x_522_; lean_object* v___x_524_; 
v___x_520_ = l_Nat_reprFast(v_fst_515_);
v___x_521_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_521_, 0, v___x_520_);
v___x_522_ = lean_box(0);
if (v_isShared_519_ == 0)
{
lean_ctor_set_tag(v___x_518_, 1);
lean_ctor_set(v___x_518_, 1, v___x_522_);
lean_ctor_set(v___x_518_, 0, v___x_521_);
v___x_524_ = v___x_518_;
goto v_reusejp_523_;
}
else
{
lean_object* v_reuseFailAlloc_539_; 
v_reuseFailAlloc_539_ = lean_alloc_ctor(1, 2, 0);
lean_ctor_set(v_reuseFailAlloc_539_, 0, v___x_521_);
lean_ctor_set(v_reuseFailAlloc_539_, 1, v___x_522_);
v___x_524_ = v_reuseFailAlloc_539_;
goto v_reusejp_523_;
}
v_reusejp_523_:
{
lean_object* v___x_525_; lean_object* v___x_526_; lean_object* v___x_527_; lean_object* v___x_528_; lean_object* v___x_529_; lean_object* v___x_530_; lean_object* v___x_531_; lean_object* v___x_532_; lean_object* v___x_533_; lean_object* v___x_534_; lean_object* v___x_535_; lean_object* v___x_536_; uint8_t v___x_537_; lean_object* v___x_538_; 
v___x_525_ = l_Nat_reprFast(v_snd_516_);
v___x_526_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_526_, 0, v___x_525_);
v___x_527_ = lean_alloc_ctor(1, 2, 0);
lean_ctor_set(v___x_527_, 0, v___x_526_);
lean_ctor_set(v___x_527_, 1, v___x_524_);
v___x_528_ = l_List_reverse___redArg(v___x_527_);
v___x_529_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__3));
v___x_530_ = l_Std_Format_joinSep___at___00Lean_Syntax_formatStxAux_spec__2(v___x_528_, v___x_529_);
v___x_531_ = lean_obj_once(&lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__2, &lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__2_once, _init_lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__2);
v___x_532_ = ((lean_object*)(lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__3));
v___x_533_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_533_, 0, v___x_532_);
lean_ctor_set(v___x_533_, 1, v___x_530_);
v___x_534_ = ((lean_object*)(lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg___closed__4));
v___x_535_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_535_, 0, v___x_533_);
lean_ctor_set(v___x_535_, 1, v___x_534_);
v___x_536_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_536_, 0, v___x_531_);
lean_ctor_set(v___x_536_, 1, v___x_535_);
v___x_537_ = 0;
v___x_538_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_538_, 0, v___x_536_);
lean_ctor_set_uint8(v___x_538_, sizeof(void*)*1, v___x_537_);
return v___x_538_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__3_spec__5_spec__7(lean_object* v_x_541_, lean_object* v_x_542_, lean_object* v_x_543_){
_start:
{
if (lean_obj_tag(v_x_543_) == 0)
{
lean_dec(v_x_541_);
return v_x_542_;
}
else
{
lean_object* v_head_544_; lean_object* v_tail_545_; lean_object* v___x_547_; uint8_t v_isShared_548_; uint8_t v_isSharedCheck_555_; 
v_head_544_ = lean_ctor_get(v_x_543_, 0);
v_tail_545_ = lean_ctor_get(v_x_543_, 1);
v_isSharedCheck_555_ = !lean_is_exclusive(v_x_543_);
if (v_isSharedCheck_555_ == 0)
{
v___x_547_ = v_x_543_;
v_isShared_548_ = v_isSharedCheck_555_;
goto v_resetjp_546_;
}
else
{
lean_inc(v_tail_545_);
lean_inc(v_head_544_);
lean_dec(v_x_543_);
v___x_547_ = lean_box(0);
v_isShared_548_ = v_isSharedCheck_555_;
goto v_resetjp_546_;
}
v_resetjp_546_:
{
lean_object* v___x_550_; 
lean_inc(v_x_541_);
if (v_isShared_548_ == 0)
{
lean_ctor_set_tag(v___x_547_, 5);
lean_ctor_set(v___x_547_, 1, v_x_541_);
lean_ctor_set(v___x_547_, 0, v_x_542_);
v___x_550_ = v___x_547_;
goto v_reusejp_549_;
}
else
{
lean_object* v_reuseFailAlloc_554_; 
v_reuseFailAlloc_554_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_554_, 0, v_x_542_);
lean_ctor_set(v_reuseFailAlloc_554_, 1, v_x_541_);
v___x_550_ = v_reuseFailAlloc_554_;
goto v_reusejp_549_;
}
v_reusejp_549_:
{
lean_object* v___x_551_; lean_object* v___x_552_; 
v___x_551_ = lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg(v_head_544_);
v___x_552_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_552_, 0, v___x_550_);
lean_ctor_set(v___x_552_, 1, v___x_551_);
v_x_542_ = v___x_552_;
v_x_543_ = v_tail_545_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__3_spec__5(lean_object* v_x_556_, lean_object* v_x_557_, lean_object* v_x_558_){
_start:
{
if (lean_obj_tag(v_x_558_) == 0)
{
lean_dec(v_x_556_);
return v_x_557_;
}
else
{
lean_object* v_head_559_; lean_object* v_tail_560_; lean_object* v___x_562_; uint8_t v_isShared_563_; uint8_t v_isSharedCheck_570_; 
v_head_559_ = lean_ctor_get(v_x_558_, 0);
v_tail_560_ = lean_ctor_get(v_x_558_, 1);
v_isSharedCheck_570_ = !lean_is_exclusive(v_x_558_);
if (v_isSharedCheck_570_ == 0)
{
v___x_562_ = v_x_558_;
v_isShared_563_ = v_isSharedCheck_570_;
goto v_resetjp_561_;
}
else
{
lean_inc(v_tail_560_);
lean_inc(v_head_559_);
lean_dec(v_x_558_);
v___x_562_ = lean_box(0);
v_isShared_563_ = v_isSharedCheck_570_;
goto v_resetjp_561_;
}
v_resetjp_561_:
{
lean_object* v___x_565_; 
lean_inc(v_x_556_);
if (v_isShared_563_ == 0)
{
lean_ctor_set_tag(v___x_562_, 5);
lean_ctor_set(v___x_562_, 1, v_x_556_);
lean_ctor_set(v___x_562_, 0, v_x_557_);
v___x_565_ = v___x_562_;
goto v_reusejp_564_;
}
else
{
lean_object* v_reuseFailAlloc_569_; 
v_reuseFailAlloc_569_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_569_, 0, v_x_557_);
lean_ctor_set(v_reuseFailAlloc_569_, 1, v_x_556_);
v___x_565_ = v_reuseFailAlloc_569_;
goto v_reusejp_564_;
}
v_reusejp_564_:
{
lean_object* v___x_566_; lean_object* v___x_567_; lean_object* v___x_568_; 
v___x_566_ = lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg(v_head_559_);
v___x_567_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_567_, 0, v___x_565_);
lean_ctor_set(v___x_567_, 1, v___x_566_);
v___x_568_ = lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__3_spec__5_spec__7(v_x_556_, v___x_567_, v_tail_560_);
return v___x_568_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__3(lean_object* v_x_571_, lean_object* v_x_572_){
_start:
{
if (lean_obj_tag(v_x_571_) == 0)
{
lean_object* v___x_573_; 
lean_dec(v_x_572_);
v___x_573_ = lean_box(0);
return v___x_573_;
}
else
{
lean_object* v_tail_574_; 
v_tail_574_ = lean_ctor_get(v_x_571_, 1);
if (lean_obj_tag(v_tail_574_) == 0)
{
lean_object* v_head_575_; lean_object* v___x_576_; 
lean_dec(v_x_572_);
v_head_575_ = lean_ctor_get(v_x_571_, 0);
lean_inc(v_head_575_);
lean_dec_ref(v_x_571_);
v___x_576_ = lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg(v_head_575_);
return v___x_576_;
}
else
{
lean_object* v_head_577_; lean_object* v___x_578_; lean_object* v___x_579_; 
lean_inc(v_tail_574_);
v_head_577_ = lean_ctor_get(v_x_571_, 0);
lean_inc(v_head_577_);
lean_dec_ref(v_x_571_);
v___x_578_ = lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg(v_head_577_);
v___x_579_ = lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__3_spec__5(v_x_572_, v___x_578_, v_tail_574_);
return v___x_579_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1___redArg(lean_object* v_a_580_){
_start:
{
if (lean_obj_tag(v_a_580_) == 0)
{
lean_object* v___x_581_; 
v___x_581_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__1));
return v___x_581_;
}
else
{
lean_object* v___x_582_; lean_object* v___x_583_; lean_object* v___x_584_; lean_object* v___x_585_; lean_object* v___x_586_; lean_object* v___x_587_; lean_object* v___x_588_; lean_object* v___x_589_; uint8_t v___x_590_; lean_object* v___x_591_; 
v___x_582_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__3));
v___x_583_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__3(v_a_580_, v___x_582_);
v___x_584_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__6, &lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__6_once, _init_lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__6);
v___x_585_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__7));
v___x_586_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_586_, 0, v___x_585_);
lean_ctor_set(v___x_586_, 1, v___x_583_);
v___x_587_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg___closed__8));
v___x_588_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_588_, 0, v___x_586_);
lean_ctor_set(v___x_588_, 1, v___x_587_);
v___x_589_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_589_, 0, v___x_584_);
lean_ctor_set(v___x_589_, 1, v___x_588_);
v___x_590_ = 0;
v___x_591_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_591_, 0, v___x_589_);
lean_ctor_set_uint8(v___x_591_, sizeof(void*)*1, v___x_590_);
return v___x_591_;
}
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__6(void){
_start:
{
lean_object* v___x_604_; lean_object* v___x_605_; 
v___x_604_ = lean_unsigned_to_nat(18u);
v___x_605_ = lean_nat_to_int(v___x_604_);
return v___x_605_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg(lean_object* v_x_606_){
_start:
{
lean_object* v_sequence_607_; lean_object* v_exonBoundaries_608_; lean_object* v___x_610_; uint8_t v_isShared_611_; uint8_t v_isSharedCheck_641_; 
v_sequence_607_ = lean_ctor_get(v_x_606_, 0);
v_exonBoundaries_608_ = lean_ctor_get(v_x_606_, 1);
v_isSharedCheck_641_ = !lean_is_exclusive(v_x_606_);
if (v_isSharedCheck_641_ == 0)
{
v___x_610_ = v_x_606_;
v_isShared_611_ = v_isSharedCheck_641_;
goto v_resetjp_609_;
}
else
{
lean_inc(v_exonBoundaries_608_);
lean_inc(v_sequence_607_);
lean_dec(v_x_606_);
v___x_610_ = lean_box(0);
v_isShared_611_ = v_isSharedCheck_641_;
goto v_resetjp_609_;
}
v_resetjp_609_:
{
lean_object* v___x_612_; lean_object* v___x_613_; lean_object* v___x_614_; lean_object* v___x_615_; lean_object* v___x_617_; 
v___x_612_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__4));
v___x_613_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__3));
v___x_614_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__7);
v___x_615_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_sequence_607_);
if (v_isShared_611_ == 0)
{
lean_ctor_set_tag(v___x_610_, 4);
lean_ctor_set(v___x_610_, 1, v___x_615_);
lean_ctor_set(v___x_610_, 0, v___x_614_);
v___x_617_ = v___x_610_;
goto v_reusejp_616_;
}
else
{
lean_object* v_reuseFailAlloc_640_; 
v_reuseFailAlloc_640_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_640_, 0, v___x_614_);
lean_ctor_set(v_reuseFailAlloc_640_, 1, v___x_615_);
v___x_617_ = v_reuseFailAlloc_640_;
goto v_reusejp_616_;
}
v_reusejp_616_:
{
uint8_t v___x_618_; lean_object* v___x_619_; lean_object* v___x_620_; lean_object* v___x_621_; lean_object* v___x_622_; lean_object* v___x_623_; lean_object* v___x_624_; lean_object* v___x_625_; lean_object* v___x_626_; lean_object* v___x_627_; lean_object* v___x_628_; lean_object* v___x_629_; lean_object* v___x_630_; lean_object* v___x_631_; lean_object* v___x_632_; lean_object* v___x_633_; lean_object* v___x_634_; lean_object* v___x_635_; lean_object* v___x_636_; lean_object* v___x_637_; lean_object* v___x_638_; lean_object* v___x_639_; 
v___x_618_ = 0;
v___x_619_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_619_, 0, v___x_617_);
lean_ctor_set_uint8(v___x_619_, sizeof(void*)*1, v___x_618_);
v___x_620_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_620_, 0, v___x_613_);
lean_ctor_set(v___x_620_, 1, v___x_619_);
v___x_621_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__9));
v___x_622_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_622_, 0, v___x_620_);
lean_ctor_set(v___x_622_, 1, v___x_621_);
v___x_623_ = lean_box(1);
v___x_624_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_624_, 0, v___x_622_);
lean_ctor_set(v___x_624_, 1, v___x_623_);
v___x_625_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__5));
v___x_626_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_626_, 0, v___x_624_);
lean_ctor_set(v___x_626_, 1, v___x_625_);
v___x_627_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_627_, 0, v___x_626_);
lean_ctor_set(v___x_627_, 1, v___x_612_);
v___x_628_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__6, &lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg___closed__6);
v___x_629_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1___redArg(v_exonBoundaries_608_);
v___x_630_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_630_, 0, v___x_628_);
lean_ctor_set(v___x_630_, 1, v___x_629_);
v___x_631_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_631_, 0, v___x_630_);
lean_ctor_set_uint8(v___x_631_, sizeof(void*)*1, v___x_618_);
v___x_632_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_632_, 0, v___x_627_);
lean_ctor_set(v___x_632_, 1, v___x_631_);
v___x_633_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__15, &lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__15_once, _init_lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__15);
v___x_634_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__16));
v___x_635_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_635_, 0, v___x_634_);
lean_ctor_set(v___x_635_, 1, v___x_632_);
v___x_636_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprCellularContext_repr___redArg___closed__17));
v___x_637_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_637_, 0, v___x_635_);
lean_ctor_set(v___x_637_, 1, v___x_636_);
v___x_638_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_638_, 0, v___x_633_);
lean_ctor_set(v___x_638_, 1, v___x_637_);
v___x_639_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_639_, 0, v___x_638_);
lean_ctor_set_uint8(v___x_639_, sizeof(void*)*1, v___x_618_);
return v___x_639_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr(lean_object* v_x_642_, lean_object* v_prec_643_){
_start:
{
lean_object* v___x_644_; 
v___x_644_ = lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___redArg(v_x_642_);
return v___x_644_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr___boxed(lean_object* v_x_645_, lean_object* v_prec_646_){
_start:
{
lean_object* v_res_647_; 
v_res_647_ = lp_BioCompiler_BioCompiler_instReprSpliceIsoform_repr(v_x_645_, v_prec_646_);
lean_dec(v_prec_646_);
return v_res_647_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0(lean_object* v_a_648_, lean_object* v_n_649_){
_start:
{
lean_object* v___x_650_; 
v___x_650_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_a_648_);
return v___x_650_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___boxed(lean_object* v_a_651_, lean_object* v_n_652_){
_start:
{
lean_object* v_res_653_; 
v_res_653_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0(v_a_651_, v_n_652_);
lean_dec(v_n_652_);
return v_res_653_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1(lean_object* v_a_654_, lean_object* v_n_655_){
_start:
{
lean_object* v___x_656_; 
v___x_656_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1___redArg(v_a_654_);
return v___x_656_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1___boxed(lean_object* v_a_657_, lean_object* v_n_658_){
_start:
{
lean_object* v_res_659_; 
v_res_659_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1(v_a_657_, v_n_658_);
lean_dec(v_n_658_);
return v_res_659_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2(lean_object* v_x_660_, lean_object* v_x_661_){
_start:
{
lean_object* v___x_662_; 
v___x_662_ = lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___redArg(v_x_660_);
return v___x_662_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2___boxed(lean_object* v_x_663_, lean_object* v_x_664_){
_start:
{
lean_object* v_res_665_; 
v_res_665_ = lp_BioCompiler_Prod_repr___at___00List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__1_spec__2(v_x_663_, v_x_664_);
lean_dec(v_x_664_);
return v_res_665_;
}
}
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Sequence(uint8_t builtin);
static bool _G_initialized = false;
LEAN_EXPORT lean_object* initialize_BioCompiler_BioCompiler_NDFST(uint8_t builtin) {
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
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
