// Lean compiler output
// Module: BioCompiler.SLOTIndependence
// Imports: public import Init public meta import Init public import BioCompiler.ThreeValued public import BioCompiler.Sequence public import BioCompiler.NDFST public import BioCompiler.Scanners public import BioCompiler.TypeSystem public import BioCompiler.Compositional public import BioCompiler.Certificates
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
lean_object* lean_nat_to_int(lean_object*);
lean_object* l_Nat_reprFast(lean_object*);
lean_object* l_Float_repr(double, lean_object*);
lean_object* lean_string_length(lean_object*);
lean_object* l_String_quote(lean_object*);
lean_object* lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(lean_object*);
lean_object* l_Repr_addAppParen(lean_object*, lean_object*);
uint8_t lean_nat_dec_eq(lean_object*, lean_object*);
lean_object* l_Int_repr(lean_object*);
lean_object* lean_string_append(lean_object*, lean_object*);
uint8_t lean_int_dec_lt(lean_object*, lean_object*);
uint8_t lean_nat_dec_le(lean_object*, lean_object*);
uint8_t l_Rat_instDecidableLe(lean_object*, lean_object*);
static const lean_string_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 3, .m_capacity = 3, .m_length = 2, .m_data = "{ "};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__0_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "atomName"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__2_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 5, .m_capacity = 5, .m_length = 4, .m_data = " := "};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__3_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__6 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__6_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = ","};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__8 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__8_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__9_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__8_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__9 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__9_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__10_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "x"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__10 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__10_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__11_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__10_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__11 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__11_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__12_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__12;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__13_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "y"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__13 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__13_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__14_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__13_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__14 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__14_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__15_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "z"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__15 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__15_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__16_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__15_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__16 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__16_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__17_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 3, .m_capacity = 3, .m_length = 2, .m_data = " }"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__17 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__17_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__18_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__18;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__20_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__20 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__20_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__21_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__17_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__21 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__21_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprAtomCoordinate___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate___closed__0_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "residueI"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "residueJ"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__5_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "paeValue"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__6 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__6_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__7_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__6_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__7 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__7_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprPAEEntry___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPAEEntry___closed__0_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 8, .m_capacity = 8, .m_length = 7, .m_data = "ptmType"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__3_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__4_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__4;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 16, .m_capacity = 16, .m_length = 15, .m_data = "residuePosition"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__5_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__5_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__6 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__6_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__7_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__7;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 6, .m_capacity = 6, .m_length = 5, .m_data = "score"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__8 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__8_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__9_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__8_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__9 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__9_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__10_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__10;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprPTMSiteEntry___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry___closed__0_value;
static const lean_string_object lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 5, .m_capacity = 5, .m_length = 4, .m_data = "none"};
static const lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__0 = (const lean_object*)&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__0_value;
static const lean_ctor_object lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__0_value)}};
static const lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__1 = (const lean_object*)&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__1_value;
static const lean_string_object lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 6, .m_capacity = 6, .m_length = 5, .m_data = "some "};
static const lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__2 = (const lean_object*)&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__2_value;
static const lean_ctor_object lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__2_value)}};
static const lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__3 = (const lean_object*)&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__3_value;
static const lean_string_object lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "("};
static const lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__4 = (const lean_object*)&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__4_value;
static const lean_string_object lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = " : Rat)/"};
static const lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__5 = (const lean_object*)&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__5_value;
static lean_once_cell_t lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__6_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__6;
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3_spec__5_spec__8_spec__11(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3_spec__5_spec__8(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3_spec__5(lean_object*, lean_object*);
static const lean_string_object lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 3, .m_capacity = 3, .m_length = 2, .m_data = "[]"};
static const lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__1_value;
static const lean_string_object lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "["};
static const lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__9_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__3_value;
static const lean_string_object lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 2, .m_capacity = 2, .m_length = 1, .m_data = "]"};
static const lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__4_value;
static lean_once_cell_t lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__5_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__5;
static lean_once_cell_t lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__7_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__2_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__7 = (const lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__7_value;
static const lean_ctor_object lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__8 = (const lean_object*)&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__8_value;
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0_spec__2_spec__5_spec__8(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0_spec__2_spec__5(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0_spec__2(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5_spec__8_spec__11_spec__14(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5_spec__8_spec__11(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5_spec__8(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3___boxed(lean_object*, lean_object*);
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 16, .m_capacity = 16, .m_length = 15, .m_data = "atomCoordinates"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 10, .m_capacity = 10, .m_length = 9, .m_data = "meanPLDDT"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__5_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__6_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__6;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__7_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 10, .m_capacity = 10, .m_length = 9, .m_data = "paeMatrix"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__7 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__7_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__7_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__8 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__8_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__9_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "ptmSites"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__9 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__9_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__10_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__9_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__10 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__10_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprSLOTValues___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprSLOTValues___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_emptySLOTS___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*4 + 0, .m_other = 4, .m_tag = 0}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)(((size_t)(0) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_emptySLOTS___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_emptySLOTS___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_emptySLOTS = (const lean_object*)&lp_BioCompiler_BioCompiler_emptySLOTS___closed__0_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 9, .m_capacity = 9, .m_length = 8, .m_data = "sequence"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 6, .m_capacity = 6, .m_length = 5, .m_data = "slots"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__5_value;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprIRRecord___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprIRRecord_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprIRRecord___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_isCorePredicate(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_isCorePredicate___boxed(lean_object*);
static const lean_string_object lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 53, .m_capacity = 53, .m_length = 52, .m_data = "BioCompiler.FFIDependentPredicate.StructureConfident"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 5}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__1_value),((lean_object*)(((size_t)(1) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__2_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__3_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__3;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__4_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__4;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluateFFIDependent(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluateFFIDependent___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_SLOTIndependence_0__BioCompiler_evaluateFFIDependent_match__3_splitter___redArg(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_SLOTIndependence_0__BioCompiler_evaluateFFIDependent_match__3_splitter(lean_object*, lean_object*, lean_object*, lean_object*);
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7(void){
_start:
{
lean_object* v___x_14_; lean_object* v___x_15_; 
v___x_14_ = lean_unsigned_to_nat(12u);
v___x_15_ = lean_nat_to_int(v___x_14_);
return v___x_15_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__12(void){
_start:
{
lean_object* v___x_22_; lean_object* v___x_23_; 
v___x_22_ = lean_unsigned_to_nat(5u);
v___x_23_ = lean_nat_to_int(v___x_22_);
return v___x_23_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__18(void){
_start:
{
lean_object* v___x_31_; lean_object* v___x_32_; 
v___x_31_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__0));
v___x_32_ = lean_string_length(v___x_31_);
return v___x_32_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19(void){
_start:
{
lean_object* v___x_33_; lean_object* v___x_34_; 
v___x_33_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__18, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__18_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__18);
v___x_34_ = lean_nat_to_int(v___x_33_);
return v___x_34_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg(lean_object* v_x_39_){
_start:
{
lean_object* v_atomName_40_; double v_x_41_; double v_y_42_; double v_z_43_; lean_object* v___x_44_; lean_object* v___x_45_; lean_object* v___x_46_; lean_object* v___x_47_; lean_object* v___x_48_; lean_object* v___x_49_; uint8_t v___x_50_; lean_object* v___x_51_; lean_object* v___x_52_; lean_object* v___x_53_; lean_object* v___x_54_; lean_object* v___x_55_; lean_object* v___x_56_; lean_object* v___x_57_; lean_object* v___x_58_; lean_object* v___x_59_; lean_object* v___x_60_; lean_object* v___x_61_; lean_object* v___x_62_; lean_object* v___x_63_; lean_object* v___x_64_; lean_object* v___x_65_; lean_object* v___x_66_; lean_object* v___x_67_; lean_object* v___x_68_; lean_object* v___x_69_; lean_object* v___x_70_; lean_object* v___x_71_; lean_object* v___x_72_; lean_object* v___x_73_; lean_object* v___x_74_; lean_object* v___x_75_; lean_object* v___x_76_; lean_object* v___x_77_; lean_object* v___x_78_; lean_object* v___x_79_; lean_object* v___x_80_; lean_object* v___x_81_; lean_object* v___x_82_; lean_object* v___x_83_; lean_object* v___x_84_; lean_object* v___x_85_; lean_object* v___x_86_; lean_object* v___x_87_; lean_object* v___x_88_; lean_object* v___x_89_; lean_object* v___x_90_; 
v_atomName_40_ = lean_ctor_get(v_x_39_, 0);
lean_inc_ref(v_atomName_40_);
v_x_41_ = lean_ctor_get_float(v_x_39_, sizeof(void*)*1);
v_y_42_ = lean_ctor_get_float(v_x_39_, sizeof(void*)*1 + 8);
v_z_43_ = lean_ctor_get_float(v_x_39_, sizeof(void*)*1 + 16);
lean_dec_ref(v_x_39_);
v___x_44_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5));
v___x_45_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__6));
v___x_46_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7);
v___x_47_ = l_String_quote(v_atomName_40_);
v___x_48_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_48_, 0, v___x_47_);
v___x_49_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_49_, 0, v___x_46_);
lean_ctor_set(v___x_49_, 1, v___x_48_);
v___x_50_ = 0;
v___x_51_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_51_, 0, v___x_49_);
lean_ctor_set_uint8(v___x_51_, sizeof(void*)*1, v___x_50_);
v___x_52_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_52_, 0, v___x_45_);
lean_ctor_set(v___x_52_, 1, v___x_51_);
v___x_53_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__9));
v___x_54_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_54_, 0, v___x_52_);
lean_ctor_set(v___x_54_, 1, v___x_53_);
v___x_55_ = lean_box(1);
v___x_56_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_56_, 0, v___x_54_);
lean_ctor_set(v___x_56_, 1, v___x_55_);
v___x_57_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__11));
v___x_58_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_58_, 0, v___x_56_);
lean_ctor_set(v___x_58_, 1, v___x_57_);
v___x_59_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_59_, 0, v___x_58_);
lean_ctor_set(v___x_59_, 1, v___x_44_);
v___x_60_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__12, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__12_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__12);
v___x_61_ = lean_unsigned_to_nat(0u);
v___x_62_ = l_Float_repr(v_x_41_, v___x_61_);
v___x_63_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_63_, 0, v___x_60_);
lean_ctor_set(v___x_63_, 1, v___x_62_);
v___x_64_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_64_, 0, v___x_63_);
lean_ctor_set_uint8(v___x_64_, sizeof(void*)*1, v___x_50_);
v___x_65_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_65_, 0, v___x_59_);
lean_ctor_set(v___x_65_, 1, v___x_64_);
v___x_66_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_66_, 0, v___x_65_);
lean_ctor_set(v___x_66_, 1, v___x_53_);
v___x_67_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_67_, 0, v___x_66_);
lean_ctor_set(v___x_67_, 1, v___x_55_);
v___x_68_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__14));
v___x_69_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_69_, 0, v___x_67_);
lean_ctor_set(v___x_69_, 1, v___x_68_);
v___x_70_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_70_, 0, v___x_69_);
lean_ctor_set(v___x_70_, 1, v___x_44_);
v___x_71_ = l_Float_repr(v_y_42_, v___x_61_);
v___x_72_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_72_, 0, v___x_60_);
lean_ctor_set(v___x_72_, 1, v___x_71_);
v___x_73_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_73_, 0, v___x_72_);
lean_ctor_set_uint8(v___x_73_, sizeof(void*)*1, v___x_50_);
v___x_74_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_74_, 0, v___x_70_);
lean_ctor_set(v___x_74_, 1, v___x_73_);
v___x_75_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_75_, 0, v___x_74_);
lean_ctor_set(v___x_75_, 1, v___x_53_);
v___x_76_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_76_, 0, v___x_75_);
lean_ctor_set(v___x_76_, 1, v___x_55_);
v___x_77_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__16));
v___x_78_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_78_, 0, v___x_76_);
lean_ctor_set(v___x_78_, 1, v___x_77_);
v___x_79_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_79_, 0, v___x_78_);
lean_ctor_set(v___x_79_, 1, v___x_44_);
v___x_80_ = l_Float_repr(v_z_43_, v___x_61_);
v___x_81_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_81_, 0, v___x_60_);
lean_ctor_set(v___x_81_, 1, v___x_80_);
v___x_82_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_82_, 0, v___x_81_);
lean_ctor_set_uint8(v___x_82_, sizeof(void*)*1, v___x_50_);
v___x_83_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_83_, 0, v___x_79_);
lean_ctor_set(v___x_83_, 1, v___x_82_);
v___x_84_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19);
v___x_85_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__20));
v___x_86_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_86_, 0, v___x_85_);
lean_ctor_set(v___x_86_, 1, v___x_83_);
v___x_87_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__21));
v___x_88_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_88_, 0, v___x_86_);
lean_ctor_set(v___x_88_, 1, v___x_87_);
v___x_89_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_89_, 0, v___x_84_);
lean_ctor_set(v___x_89_, 1, v___x_88_);
v___x_90_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_90_, 0, v___x_89_);
lean_ctor_set_uint8(v___x_90_, sizeof(void*)*1, v___x_50_);
return v___x_90_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr(lean_object* v_x_91_, lean_object* v_prec_92_){
_start:
{
lean_object* v___x_93_; 
v___x_93_ = lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg(v_x_91_);
return v___x_93_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___boxed(lean_object* v_x_94_, lean_object* v_prec_95_){
_start:
{
lean_object* v_res_96_; 
v_res_96_ = lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr(v_x_94_, v_prec_95_);
lean_dec(v_prec_95_);
return v_res_96_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg(lean_object* v_x_114_){
_start:
{
lean_object* v_residueI_115_; lean_object* v_residueJ_116_; double v_paeValue_117_; lean_object* v___x_118_; lean_object* v___x_119_; lean_object* v___x_120_; lean_object* v___x_121_; lean_object* v___x_122_; lean_object* v___x_123_; uint8_t v___x_124_; lean_object* v___x_125_; lean_object* v___x_126_; lean_object* v___x_127_; lean_object* v___x_128_; lean_object* v___x_129_; lean_object* v___x_130_; lean_object* v___x_131_; lean_object* v___x_132_; lean_object* v___x_133_; lean_object* v___x_134_; lean_object* v___x_135_; lean_object* v___x_136_; lean_object* v___x_137_; lean_object* v___x_138_; lean_object* v___x_139_; lean_object* v___x_140_; lean_object* v___x_141_; lean_object* v___x_142_; lean_object* v___x_143_; lean_object* v___x_144_; lean_object* v___x_145_; lean_object* v___x_146_; lean_object* v___x_147_; lean_object* v___x_148_; lean_object* v___x_149_; lean_object* v___x_150_; lean_object* v___x_151_; lean_object* v___x_152_; lean_object* v___x_153_; lean_object* v___x_154_; lean_object* v___x_155_; 
v_residueI_115_ = lean_ctor_get(v_x_114_, 0);
lean_inc(v_residueI_115_);
v_residueJ_116_ = lean_ctor_get(v_x_114_, 1);
lean_inc(v_residueJ_116_);
v_paeValue_117_ = lean_ctor_get_float(v_x_114_, sizeof(void*)*2);
lean_dec_ref(v_x_114_);
v___x_118_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5));
v___x_119_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__3));
v___x_120_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7);
v___x_121_ = l_Nat_reprFast(v_residueI_115_);
v___x_122_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_122_, 0, v___x_121_);
v___x_123_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_123_, 0, v___x_120_);
lean_ctor_set(v___x_123_, 1, v___x_122_);
v___x_124_ = 0;
v___x_125_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_125_, 0, v___x_123_);
lean_ctor_set_uint8(v___x_125_, sizeof(void*)*1, v___x_124_);
v___x_126_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_126_, 0, v___x_119_);
lean_ctor_set(v___x_126_, 1, v___x_125_);
v___x_127_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__9));
v___x_128_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_128_, 0, v___x_126_);
lean_ctor_set(v___x_128_, 1, v___x_127_);
v___x_129_ = lean_box(1);
v___x_130_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_130_, 0, v___x_128_);
lean_ctor_set(v___x_130_, 1, v___x_129_);
v___x_131_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__5));
v___x_132_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_132_, 0, v___x_130_);
lean_ctor_set(v___x_132_, 1, v___x_131_);
v___x_133_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_133_, 0, v___x_132_);
lean_ctor_set(v___x_133_, 1, v___x_118_);
v___x_134_ = l_Nat_reprFast(v_residueJ_116_);
v___x_135_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_135_, 0, v___x_134_);
v___x_136_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_136_, 0, v___x_120_);
lean_ctor_set(v___x_136_, 1, v___x_135_);
v___x_137_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_137_, 0, v___x_136_);
lean_ctor_set_uint8(v___x_137_, sizeof(void*)*1, v___x_124_);
v___x_138_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_138_, 0, v___x_133_);
lean_ctor_set(v___x_138_, 1, v___x_137_);
v___x_139_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_139_, 0, v___x_138_);
lean_ctor_set(v___x_139_, 1, v___x_127_);
v___x_140_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_140_, 0, v___x_139_);
lean_ctor_set(v___x_140_, 1, v___x_129_);
v___x_141_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg___closed__7));
v___x_142_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_142_, 0, v___x_140_);
lean_ctor_set(v___x_142_, 1, v___x_141_);
v___x_143_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_143_, 0, v___x_142_);
lean_ctor_set(v___x_143_, 1, v___x_118_);
v___x_144_ = lean_unsigned_to_nat(0u);
v___x_145_ = l_Float_repr(v_paeValue_117_, v___x_144_);
v___x_146_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_146_, 0, v___x_120_);
lean_ctor_set(v___x_146_, 1, v___x_145_);
v___x_147_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_147_, 0, v___x_146_);
lean_ctor_set_uint8(v___x_147_, sizeof(void*)*1, v___x_124_);
v___x_148_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_148_, 0, v___x_143_);
lean_ctor_set(v___x_148_, 1, v___x_147_);
v___x_149_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19);
v___x_150_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__20));
v___x_151_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_151_, 0, v___x_150_);
lean_ctor_set(v___x_151_, 1, v___x_148_);
v___x_152_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__21));
v___x_153_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_153_, 0, v___x_151_);
lean_ctor_set(v___x_153_, 1, v___x_152_);
v___x_154_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_154_, 0, v___x_149_);
lean_ctor_set(v___x_154_, 1, v___x_153_);
v___x_155_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_155_, 0, v___x_154_);
lean_ctor_set_uint8(v___x_155_, sizeof(void*)*1, v___x_124_);
return v___x_155_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr(lean_object* v_x_156_, lean_object* v_prec_157_){
_start:
{
lean_object* v___x_158_; 
v___x_158_ = lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg(v_x_156_);
return v___x_158_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___boxed(lean_object* v_x_159_, lean_object* v_prec_160_){
_start:
{
lean_object* v_res_161_; 
v_res_161_ = lp_BioCompiler_BioCompiler_instReprPAEEntry_repr(v_x_159_, v_prec_160_);
lean_dec(v_prec_160_);
return v_res_161_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__4(void){
_start:
{
lean_object* v___x_173_; lean_object* v___x_174_; 
v___x_173_ = lean_unsigned_to_nat(11u);
v___x_174_ = lean_nat_to_int(v___x_173_);
return v___x_174_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__7(void){
_start:
{
lean_object* v___x_178_; lean_object* v___x_179_; 
v___x_178_ = lean_unsigned_to_nat(19u);
v___x_179_ = lean_nat_to_int(v___x_178_);
return v___x_179_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__10(void){
_start:
{
lean_object* v___x_183_; lean_object* v___x_184_; 
v___x_183_ = lean_unsigned_to_nat(9u);
v___x_184_ = lean_nat_to_int(v___x_183_);
return v___x_184_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg(lean_object* v_x_185_){
_start:
{
lean_object* v_ptmType_186_; lean_object* v_residuePosition_187_; double v_score_188_; lean_object* v___x_189_; lean_object* v___x_190_; lean_object* v___x_191_; lean_object* v___x_192_; lean_object* v___x_193_; lean_object* v___x_194_; uint8_t v___x_195_; lean_object* v___x_196_; lean_object* v___x_197_; lean_object* v___x_198_; lean_object* v___x_199_; lean_object* v___x_200_; lean_object* v___x_201_; lean_object* v___x_202_; lean_object* v___x_203_; lean_object* v___x_204_; lean_object* v___x_205_; lean_object* v___x_206_; lean_object* v___x_207_; lean_object* v___x_208_; lean_object* v___x_209_; lean_object* v___x_210_; lean_object* v___x_211_; lean_object* v___x_212_; lean_object* v___x_213_; lean_object* v___x_214_; lean_object* v___x_215_; lean_object* v___x_216_; lean_object* v___x_217_; lean_object* v___x_218_; lean_object* v___x_219_; lean_object* v___x_220_; lean_object* v___x_221_; lean_object* v___x_222_; lean_object* v___x_223_; lean_object* v___x_224_; lean_object* v___x_225_; lean_object* v___x_226_; lean_object* v___x_227_; lean_object* v___x_228_; 
v_ptmType_186_ = lean_ctor_get(v_x_185_, 0);
lean_inc_ref(v_ptmType_186_);
v_residuePosition_187_ = lean_ctor_get(v_x_185_, 1);
lean_inc(v_residuePosition_187_);
v_score_188_ = lean_ctor_get_float(v_x_185_, sizeof(void*)*2);
lean_dec_ref(v_x_185_);
v___x_189_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5));
v___x_190_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__3));
v___x_191_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__4, &lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__4_once, _init_lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__4);
v___x_192_ = l_String_quote(v_ptmType_186_);
v___x_193_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_193_, 0, v___x_192_);
v___x_194_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_194_, 0, v___x_191_);
lean_ctor_set(v___x_194_, 1, v___x_193_);
v___x_195_ = 0;
v___x_196_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_196_, 0, v___x_194_);
lean_ctor_set_uint8(v___x_196_, sizeof(void*)*1, v___x_195_);
v___x_197_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_197_, 0, v___x_190_);
lean_ctor_set(v___x_197_, 1, v___x_196_);
v___x_198_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__9));
v___x_199_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_199_, 0, v___x_197_);
lean_ctor_set(v___x_199_, 1, v___x_198_);
v___x_200_ = lean_box(1);
v___x_201_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_201_, 0, v___x_199_);
lean_ctor_set(v___x_201_, 1, v___x_200_);
v___x_202_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__6));
v___x_203_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_203_, 0, v___x_201_);
lean_ctor_set(v___x_203_, 1, v___x_202_);
v___x_204_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_204_, 0, v___x_203_);
lean_ctor_set(v___x_204_, 1, v___x_189_);
v___x_205_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__7);
v___x_206_ = l_Nat_reprFast(v_residuePosition_187_);
v___x_207_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_207_, 0, v___x_206_);
v___x_208_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_208_, 0, v___x_205_);
lean_ctor_set(v___x_208_, 1, v___x_207_);
v___x_209_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_209_, 0, v___x_208_);
lean_ctor_set_uint8(v___x_209_, sizeof(void*)*1, v___x_195_);
v___x_210_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_210_, 0, v___x_204_);
lean_ctor_set(v___x_210_, 1, v___x_209_);
v___x_211_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_211_, 0, v___x_210_);
lean_ctor_set(v___x_211_, 1, v___x_198_);
v___x_212_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_212_, 0, v___x_211_);
lean_ctor_set(v___x_212_, 1, v___x_200_);
v___x_213_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__9));
v___x_214_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_214_, 0, v___x_212_);
lean_ctor_set(v___x_214_, 1, v___x_213_);
v___x_215_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_215_, 0, v___x_214_);
lean_ctor_set(v___x_215_, 1, v___x_189_);
v___x_216_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__10, &lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__10);
v___x_217_ = lean_unsigned_to_nat(0u);
v___x_218_ = l_Float_repr(v_score_188_, v___x_217_);
v___x_219_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_219_, 0, v___x_216_);
lean_ctor_set(v___x_219_, 1, v___x_218_);
v___x_220_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_220_, 0, v___x_219_);
lean_ctor_set_uint8(v___x_220_, sizeof(void*)*1, v___x_195_);
v___x_221_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_221_, 0, v___x_215_);
lean_ctor_set(v___x_221_, 1, v___x_220_);
v___x_222_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19);
v___x_223_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__20));
v___x_224_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_224_, 0, v___x_223_);
lean_ctor_set(v___x_224_, 1, v___x_221_);
v___x_225_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__21));
v___x_226_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_226_, 0, v___x_224_);
lean_ctor_set(v___x_226_, 1, v___x_225_);
v___x_227_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_227_, 0, v___x_222_);
lean_ctor_set(v___x_227_, 1, v___x_226_);
v___x_228_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_228_, 0, v___x_227_);
lean_ctor_set_uint8(v___x_228_, sizeof(void*)*1, v___x_195_);
return v___x_228_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr(lean_object* v_x_229_, lean_object* v_prec_230_){
_start:
{
lean_object* v___x_231_; 
v___x_231_ = lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg(v_x_229_);
return v___x_231_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___boxed(lean_object* v_x_232_, lean_object* v_prec_233_){
_start:
{
lean_object* v_res_234_; 
v_res_234_ = lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr(v_x_232_, v_prec_233_);
lean_dec(v_prec_233_);
return v_res_234_;
}
}
static lean_object* _init_lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__6(void){
_start:
{
lean_object* v___x_245_; lean_object* v___x_246_; 
v___x_245_ = lean_unsigned_to_nat(0u);
v___x_246_ = lean_nat_to_int(v___x_245_);
return v___x_246_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1(lean_object* v_x_247_, lean_object* v_x_248_){
_start:
{
if (lean_obj_tag(v_x_247_) == 0)
{
lean_object* v___x_249_; 
v___x_249_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__1));
return v___x_249_;
}
else
{
lean_object* v_val_250_; lean_object* v___x_252_; uint8_t v_isShared_253_; uint8_t v_isSharedCheck_291_; 
v_val_250_ = lean_ctor_get(v_x_247_, 0);
v_isSharedCheck_291_ = !lean_is_exclusive(v_x_247_);
if (v_isSharedCheck_291_ == 0)
{
v___x_252_ = v_x_247_;
v_isShared_253_ = v_isSharedCheck_291_;
goto v_resetjp_251_;
}
else
{
lean_inc(v_val_250_);
lean_dec(v_x_247_);
v___x_252_ = lean_box(0);
v_isShared_253_ = v_isSharedCheck_291_;
goto v_resetjp_251_;
}
v_resetjp_251_:
{
lean_object* v_num_254_; lean_object* v_den_255_; lean_object* v___x_257_; uint8_t v_isShared_258_; uint8_t v_isSharedCheck_290_; 
v_num_254_ = lean_ctor_get(v_val_250_, 0);
v_den_255_ = lean_ctor_get(v_val_250_, 1);
v_isSharedCheck_290_ = !lean_is_exclusive(v_val_250_);
if (v_isSharedCheck_290_ == 0)
{
v___x_257_ = v_val_250_;
v_isShared_258_ = v_isSharedCheck_290_;
goto v_resetjp_256_;
}
else
{
lean_inc(v_den_255_);
lean_inc(v_num_254_);
lean_dec(v_val_250_);
v___x_257_ = lean_box(0);
v_isShared_258_ = v_isSharedCheck_290_;
goto v_resetjp_256_;
}
v_resetjp_256_:
{
lean_object* v___x_259_; lean_object* v___y_261_; lean_object* v___x_266_; uint8_t v___x_267_; 
v___x_259_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__3));
v___x_266_ = lean_unsigned_to_nat(1u);
v___x_267_ = lean_nat_dec_eq(v_den_255_, v___x_266_);
if (v___x_267_ == 0)
{
lean_object* v___x_268_; lean_object* v___x_269_; lean_object* v___x_270_; lean_object* v___x_271_; lean_object* v___x_272_; lean_object* v___x_273_; lean_object* v___x_274_; lean_object* v___x_276_; 
v___x_268_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__4));
v___x_269_ = l_Int_repr(v_num_254_);
lean_dec(v_num_254_);
v___x_270_ = lean_string_append(v___x_268_, v___x_269_);
lean_dec_ref(v___x_269_);
v___x_271_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__5));
v___x_272_ = lean_string_append(v___x_270_, v___x_271_);
v___x_273_ = l_Nat_reprFast(v_den_255_);
v___x_274_ = lean_string_append(v___x_272_, v___x_273_);
lean_dec_ref(v___x_273_);
if (v_isShared_253_ == 0)
{
lean_ctor_set_tag(v___x_252_, 3);
lean_ctor_set(v___x_252_, 0, v___x_274_);
v___x_276_ = v___x_252_;
goto v_reusejp_275_;
}
else
{
lean_object* v_reuseFailAlloc_277_; 
v_reuseFailAlloc_277_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_277_, 0, v___x_274_);
v___x_276_ = v_reuseFailAlloc_277_;
goto v_reusejp_275_;
}
v_reusejp_275_:
{
v___y_261_ = v___x_276_;
goto v___jp_260_;
}
}
else
{
lean_object* v___x_278_; lean_object* v___x_279_; uint8_t v___x_280_; 
lean_dec(v_den_255_);
v___x_278_ = lean_unsigned_to_nat(0u);
v___x_279_ = lean_obj_once(&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__6, &lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__6_once, _init_lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__6);
v___x_280_ = lean_int_dec_lt(v_num_254_, v___x_279_);
if (v___x_280_ == 0)
{
lean_object* v___x_281_; lean_object* v___x_283_; 
v___x_281_ = l_Int_repr(v_num_254_);
lean_dec(v_num_254_);
if (v_isShared_253_ == 0)
{
lean_ctor_set_tag(v___x_252_, 3);
lean_ctor_set(v___x_252_, 0, v___x_281_);
v___x_283_ = v___x_252_;
goto v_reusejp_282_;
}
else
{
lean_object* v_reuseFailAlloc_284_; 
v_reuseFailAlloc_284_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_284_, 0, v___x_281_);
v___x_283_ = v_reuseFailAlloc_284_;
goto v_reusejp_282_;
}
v_reusejp_282_:
{
v___y_261_ = v___x_283_;
goto v___jp_260_;
}
}
else
{
lean_object* v___x_285_; lean_object* v___x_287_; 
v___x_285_ = l_Int_repr(v_num_254_);
lean_dec(v_num_254_);
if (v_isShared_253_ == 0)
{
lean_ctor_set_tag(v___x_252_, 3);
lean_ctor_set(v___x_252_, 0, v___x_285_);
v___x_287_ = v___x_252_;
goto v_reusejp_286_;
}
else
{
lean_object* v_reuseFailAlloc_289_; 
v_reuseFailAlloc_289_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v_reuseFailAlloc_289_, 0, v___x_285_);
v___x_287_ = v_reuseFailAlloc_289_;
goto v_reusejp_286_;
}
v_reusejp_286_:
{
lean_object* v___x_288_; 
v___x_288_ = l_Repr_addAppParen(v___x_287_, v___x_278_);
v___y_261_ = v___x_288_;
goto v___jp_260_;
}
}
}
v___jp_260_:
{
lean_object* v___x_263_; 
if (v_isShared_258_ == 0)
{
lean_ctor_set_tag(v___x_257_, 5);
lean_ctor_set(v___x_257_, 1, v___y_261_);
lean_ctor_set(v___x_257_, 0, v___x_259_);
v___x_263_ = v___x_257_;
goto v_reusejp_262_;
}
else
{
lean_object* v_reuseFailAlloc_265_; 
v_reuseFailAlloc_265_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_265_, 0, v___x_259_);
lean_ctor_set(v_reuseFailAlloc_265_, 1, v___y_261_);
v___x_263_ = v_reuseFailAlloc_265_;
goto v_reusejp_262_;
}
v_reusejp_262_:
{
lean_object* v___x_264_; 
v___x_264_ = l_Repr_addAppParen(v___x_263_, v_x_248_);
return v___x_264_;
}
}
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___boxed(lean_object* v_x_292_, lean_object* v_x_293_){
_start:
{
lean_object* v_res_294_; 
v_res_294_ = lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1(v_x_292_, v_x_293_);
lean_dec(v_x_293_);
return v_res_294_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3_spec__5_spec__8_spec__11(lean_object* v_x_295_, lean_object* v_x_296_, lean_object* v_x_297_){
_start:
{
if (lean_obj_tag(v_x_297_) == 0)
{
lean_dec(v_x_295_);
return v_x_296_;
}
else
{
lean_object* v_head_298_; lean_object* v_tail_299_; lean_object* v___x_301_; uint8_t v_isShared_302_; uint8_t v_isSharedCheck_309_; 
v_head_298_ = lean_ctor_get(v_x_297_, 0);
v_tail_299_ = lean_ctor_get(v_x_297_, 1);
v_isSharedCheck_309_ = !lean_is_exclusive(v_x_297_);
if (v_isSharedCheck_309_ == 0)
{
v___x_301_ = v_x_297_;
v_isShared_302_ = v_isSharedCheck_309_;
goto v_resetjp_300_;
}
else
{
lean_inc(v_tail_299_);
lean_inc(v_head_298_);
lean_dec(v_x_297_);
v___x_301_ = lean_box(0);
v_isShared_302_ = v_isSharedCheck_309_;
goto v_resetjp_300_;
}
v_resetjp_300_:
{
lean_object* v___x_304_; 
lean_inc(v_x_295_);
if (v_isShared_302_ == 0)
{
lean_ctor_set_tag(v___x_301_, 5);
lean_ctor_set(v___x_301_, 1, v_x_295_);
lean_ctor_set(v___x_301_, 0, v_x_296_);
v___x_304_ = v___x_301_;
goto v_reusejp_303_;
}
else
{
lean_object* v_reuseFailAlloc_308_; 
v_reuseFailAlloc_308_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_308_, 0, v_x_296_);
lean_ctor_set(v_reuseFailAlloc_308_, 1, v_x_295_);
v___x_304_ = v_reuseFailAlloc_308_;
goto v_reusejp_303_;
}
v_reusejp_303_:
{
lean_object* v___x_305_; lean_object* v___x_306_; 
v___x_305_ = lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg(v_head_298_);
v___x_306_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_306_, 0, v___x_304_);
lean_ctor_set(v___x_306_, 1, v___x_305_);
v_x_296_ = v___x_306_;
v_x_297_ = v_tail_299_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3_spec__5_spec__8(lean_object* v_x_310_, lean_object* v_x_311_, lean_object* v_x_312_){
_start:
{
if (lean_obj_tag(v_x_312_) == 0)
{
lean_dec(v_x_310_);
return v_x_311_;
}
else
{
lean_object* v_head_313_; lean_object* v_tail_314_; lean_object* v___x_316_; uint8_t v_isShared_317_; uint8_t v_isSharedCheck_324_; 
v_head_313_ = lean_ctor_get(v_x_312_, 0);
v_tail_314_ = lean_ctor_get(v_x_312_, 1);
v_isSharedCheck_324_ = !lean_is_exclusive(v_x_312_);
if (v_isSharedCheck_324_ == 0)
{
v___x_316_ = v_x_312_;
v_isShared_317_ = v_isSharedCheck_324_;
goto v_resetjp_315_;
}
else
{
lean_inc(v_tail_314_);
lean_inc(v_head_313_);
lean_dec(v_x_312_);
v___x_316_ = lean_box(0);
v_isShared_317_ = v_isSharedCheck_324_;
goto v_resetjp_315_;
}
v_resetjp_315_:
{
lean_object* v___x_319_; 
lean_inc(v_x_310_);
if (v_isShared_317_ == 0)
{
lean_ctor_set_tag(v___x_316_, 5);
lean_ctor_set(v___x_316_, 1, v_x_310_);
lean_ctor_set(v___x_316_, 0, v_x_311_);
v___x_319_ = v___x_316_;
goto v_reusejp_318_;
}
else
{
lean_object* v_reuseFailAlloc_323_; 
v_reuseFailAlloc_323_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_323_, 0, v_x_311_);
lean_ctor_set(v_reuseFailAlloc_323_, 1, v_x_310_);
v___x_319_ = v_reuseFailAlloc_323_;
goto v_reusejp_318_;
}
v_reusejp_318_:
{
lean_object* v___x_320_; lean_object* v___x_321_; lean_object* v___x_322_; 
v___x_320_ = lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg(v_head_313_);
v___x_321_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_321_, 0, v___x_319_);
lean_ctor_set(v___x_321_, 1, v___x_320_);
v___x_322_ = lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3_spec__5_spec__8_spec__11(v_x_310_, v___x_321_, v_tail_314_);
return v___x_322_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3_spec__5(lean_object* v_x_325_, lean_object* v_x_326_){
_start:
{
if (lean_obj_tag(v_x_325_) == 0)
{
lean_object* v___x_327_; 
lean_dec(v_x_326_);
v___x_327_ = lean_box(0);
return v___x_327_;
}
else
{
lean_object* v_tail_328_; 
v_tail_328_ = lean_ctor_get(v_x_325_, 1);
if (lean_obj_tag(v_tail_328_) == 0)
{
lean_object* v_head_329_; lean_object* v___x_330_; 
lean_dec(v_x_326_);
v_head_329_ = lean_ctor_get(v_x_325_, 0);
lean_inc(v_head_329_);
lean_dec_ref(v_x_325_);
v___x_330_ = lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg(v_head_329_);
return v___x_330_;
}
else
{
lean_object* v_head_331_; lean_object* v___x_332_; lean_object* v___x_333_; 
lean_inc(v_tail_328_);
v_head_331_ = lean_ctor_get(v_x_325_, 0);
lean_inc(v_head_331_);
lean_dec_ref(v_x_325_);
v___x_332_ = lp_BioCompiler_BioCompiler_instReprPAEEntry_repr___redArg(v_head_331_);
v___x_333_ = lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3_spec__5_spec__8(v_x_326_, v___x_332_, v_tail_328_);
return v___x_333_;
}
}
}
}
static lean_object* _init_lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__5(void){
_start:
{
lean_object* v___x_342_; lean_object* v___x_343_; 
v___x_342_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__2));
v___x_343_ = lean_string_length(v___x_342_);
return v___x_343_;
}
}
static lean_object* _init_lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6(void){
_start:
{
lean_object* v___x_344_; lean_object* v___x_345_; 
v___x_344_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__5, &lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__5_once, _init_lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__5);
v___x_345_ = lean_nat_to_int(v___x_344_);
return v___x_345_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg(lean_object* v_a_350_){
_start:
{
if (lean_obj_tag(v_a_350_) == 0)
{
lean_object* v___x_351_; 
v___x_351_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__1));
return v___x_351_;
}
else
{
lean_object* v___x_352_; lean_object* v___x_353_; lean_object* v___x_354_; lean_object* v___x_355_; lean_object* v___x_356_; lean_object* v___x_357_; lean_object* v___x_358_; lean_object* v___x_359_; uint8_t v___x_360_; lean_object* v___x_361_; 
v___x_352_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__3));
v___x_353_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3_spec__5(v_a_350_, v___x_352_);
v___x_354_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6, &lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6_once, _init_lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6);
v___x_355_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__7));
v___x_356_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_356_, 0, v___x_355_);
lean_ctor_set(v___x_356_, 1, v___x_353_);
v___x_357_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__8));
v___x_358_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_358_, 0, v___x_356_);
lean_ctor_set(v___x_358_, 1, v___x_357_);
v___x_359_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_359_, 0, v___x_354_);
lean_ctor_set(v___x_359_, 1, v___x_358_);
v___x_360_ = 0;
v___x_361_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_361_, 0, v___x_359_);
lean_ctor_set_uint8(v___x_361_, sizeof(void*)*1, v___x_360_);
return v___x_361_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2(lean_object* v_x_362_, lean_object* v_x_363_){
_start:
{
if (lean_obj_tag(v_x_362_) == 0)
{
lean_object* v___x_364_; 
v___x_364_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__1));
return v___x_364_;
}
else
{
lean_object* v_val_365_; lean_object* v___x_366_; lean_object* v___x_367_; lean_object* v___x_368_; lean_object* v___x_369_; 
v_val_365_ = lean_ctor_get(v_x_362_, 0);
lean_inc(v_val_365_);
lean_dec_ref(v_x_362_);
v___x_366_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__3));
v___x_367_ = lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg(v_val_365_);
v___x_368_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_368_, 0, v___x_366_);
lean_ctor_set(v___x_368_, 1, v___x_367_);
v___x_369_ = l_Repr_addAppParen(v___x_368_, v_x_363_);
return v___x_369_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2___boxed(lean_object* v_x_370_, lean_object* v_x_371_){
_start:
{
lean_object* v_res_372_; 
v_res_372_ = lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2(v_x_370_, v_x_371_);
lean_dec(v_x_371_);
return v_res_372_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0_spec__2_spec__5_spec__8(lean_object* v_x_373_, lean_object* v_x_374_, lean_object* v_x_375_){
_start:
{
if (lean_obj_tag(v_x_375_) == 0)
{
lean_dec(v_x_373_);
return v_x_374_;
}
else
{
lean_object* v_head_376_; lean_object* v_tail_377_; lean_object* v___x_379_; uint8_t v_isShared_380_; uint8_t v_isSharedCheck_387_; 
v_head_376_ = lean_ctor_get(v_x_375_, 0);
v_tail_377_ = lean_ctor_get(v_x_375_, 1);
v_isSharedCheck_387_ = !lean_is_exclusive(v_x_375_);
if (v_isSharedCheck_387_ == 0)
{
v___x_379_ = v_x_375_;
v_isShared_380_ = v_isSharedCheck_387_;
goto v_resetjp_378_;
}
else
{
lean_inc(v_tail_377_);
lean_inc(v_head_376_);
lean_dec(v_x_375_);
v___x_379_ = lean_box(0);
v_isShared_380_ = v_isSharedCheck_387_;
goto v_resetjp_378_;
}
v_resetjp_378_:
{
lean_object* v___x_382_; 
lean_inc(v_x_373_);
if (v_isShared_380_ == 0)
{
lean_ctor_set_tag(v___x_379_, 5);
lean_ctor_set(v___x_379_, 1, v_x_373_);
lean_ctor_set(v___x_379_, 0, v_x_374_);
v___x_382_ = v___x_379_;
goto v_reusejp_381_;
}
else
{
lean_object* v_reuseFailAlloc_386_; 
v_reuseFailAlloc_386_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_386_, 0, v_x_374_);
lean_ctor_set(v_reuseFailAlloc_386_, 1, v_x_373_);
v___x_382_ = v_reuseFailAlloc_386_;
goto v_reusejp_381_;
}
v_reusejp_381_:
{
lean_object* v___x_383_; lean_object* v___x_384_; 
v___x_383_ = lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg(v_head_376_);
v___x_384_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_384_, 0, v___x_382_);
lean_ctor_set(v___x_384_, 1, v___x_383_);
v_x_374_ = v___x_384_;
v_x_375_ = v_tail_377_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0_spec__2_spec__5(lean_object* v_x_388_, lean_object* v_x_389_, lean_object* v_x_390_){
_start:
{
if (lean_obj_tag(v_x_390_) == 0)
{
lean_dec(v_x_388_);
return v_x_389_;
}
else
{
lean_object* v_head_391_; lean_object* v_tail_392_; lean_object* v___x_394_; uint8_t v_isShared_395_; uint8_t v_isSharedCheck_402_; 
v_head_391_ = lean_ctor_get(v_x_390_, 0);
v_tail_392_ = lean_ctor_get(v_x_390_, 1);
v_isSharedCheck_402_ = !lean_is_exclusive(v_x_390_);
if (v_isSharedCheck_402_ == 0)
{
v___x_394_ = v_x_390_;
v_isShared_395_ = v_isSharedCheck_402_;
goto v_resetjp_393_;
}
else
{
lean_inc(v_tail_392_);
lean_inc(v_head_391_);
lean_dec(v_x_390_);
v___x_394_ = lean_box(0);
v_isShared_395_ = v_isSharedCheck_402_;
goto v_resetjp_393_;
}
v_resetjp_393_:
{
lean_object* v___x_397_; 
lean_inc(v_x_388_);
if (v_isShared_395_ == 0)
{
lean_ctor_set_tag(v___x_394_, 5);
lean_ctor_set(v___x_394_, 1, v_x_388_);
lean_ctor_set(v___x_394_, 0, v_x_389_);
v___x_397_ = v___x_394_;
goto v_reusejp_396_;
}
else
{
lean_object* v_reuseFailAlloc_401_; 
v_reuseFailAlloc_401_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_401_, 0, v_x_389_);
lean_ctor_set(v_reuseFailAlloc_401_, 1, v_x_388_);
v___x_397_ = v_reuseFailAlloc_401_;
goto v_reusejp_396_;
}
v_reusejp_396_:
{
lean_object* v___x_398_; lean_object* v___x_399_; lean_object* v___x_400_; 
v___x_398_ = lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg(v_head_391_);
v___x_399_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_399_, 0, v___x_397_);
lean_ctor_set(v___x_399_, 1, v___x_398_);
v___x_400_ = lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0_spec__2_spec__5_spec__8(v_x_388_, v___x_399_, v_tail_392_);
return v___x_400_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0_spec__2(lean_object* v_x_403_, lean_object* v_x_404_){
_start:
{
if (lean_obj_tag(v_x_403_) == 0)
{
lean_object* v___x_405_; 
lean_dec(v_x_404_);
v___x_405_ = lean_box(0);
return v___x_405_;
}
else
{
lean_object* v_tail_406_; 
v_tail_406_ = lean_ctor_get(v_x_403_, 1);
if (lean_obj_tag(v_tail_406_) == 0)
{
lean_object* v_head_407_; lean_object* v___x_408_; 
lean_dec(v_x_404_);
v_head_407_ = lean_ctor_get(v_x_403_, 0);
lean_inc(v_head_407_);
lean_dec_ref(v_x_403_);
v___x_408_ = lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg(v_head_407_);
return v___x_408_;
}
else
{
lean_object* v_head_409_; lean_object* v___x_410_; lean_object* v___x_411_; 
lean_inc(v_tail_406_);
v_head_409_ = lean_ctor_get(v_x_403_, 0);
lean_inc(v_head_409_);
lean_dec_ref(v_x_403_);
v___x_410_ = lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg(v_head_409_);
v___x_411_ = lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0_spec__2_spec__5(v_x_404_, v___x_410_, v_tail_406_);
return v___x_411_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0___redArg(lean_object* v_a_412_){
_start:
{
if (lean_obj_tag(v_a_412_) == 0)
{
lean_object* v___x_413_; 
v___x_413_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__1));
return v___x_413_;
}
else
{
lean_object* v___x_414_; lean_object* v___x_415_; lean_object* v___x_416_; lean_object* v___x_417_; lean_object* v___x_418_; lean_object* v___x_419_; lean_object* v___x_420_; lean_object* v___x_421_; uint8_t v___x_422_; lean_object* v___x_423_; 
v___x_414_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__3));
v___x_415_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0_spec__2(v_a_412_, v___x_414_);
v___x_416_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6, &lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6_once, _init_lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6);
v___x_417_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__7));
v___x_418_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_418_, 0, v___x_417_);
lean_ctor_set(v___x_418_, 1, v___x_415_);
v___x_419_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__8));
v___x_420_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_420_, 0, v___x_418_);
lean_ctor_set(v___x_420_, 1, v___x_419_);
v___x_421_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_421_, 0, v___x_416_);
lean_ctor_set(v___x_421_, 1, v___x_420_);
v___x_422_ = 0;
v___x_423_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_423_, 0, v___x_421_);
lean_ctor_set_uint8(v___x_423_, sizeof(void*)*1, v___x_422_);
return v___x_423_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0(lean_object* v_x_424_, lean_object* v_x_425_){
_start:
{
if (lean_obj_tag(v_x_424_) == 0)
{
lean_object* v___x_426_; 
v___x_426_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__1));
return v___x_426_;
}
else
{
lean_object* v_val_427_; lean_object* v___x_428_; lean_object* v___x_429_; lean_object* v___x_430_; lean_object* v___x_431_; 
v_val_427_ = lean_ctor_get(v_x_424_, 0);
lean_inc(v_val_427_);
lean_dec_ref(v_x_424_);
v___x_428_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__3));
v___x_429_ = lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0___redArg(v_val_427_);
v___x_430_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_430_, 0, v___x_428_);
lean_ctor_set(v___x_430_, 1, v___x_429_);
v___x_431_ = l_Repr_addAppParen(v___x_430_, v_x_425_);
return v___x_431_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0___boxed(lean_object* v_x_432_, lean_object* v_x_433_){
_start:
{
lean_object* v_res_434_; 
v_res_434_ = lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0(v_x_432_, v_x_433_);
lean_dec(v_x_433_);
return v_res_434_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5_spec__8_spec__11_spec__14(lean_object* v_x_435_, lean_object* v_x_436_, lean_object* v_x_437_){
_start:
{
if (lean_obj_tag(v_x_437_) == 0)
{
lean_dec(v_x_435_);
return v_x_436_;
}
else
{
lean_object* v_head_438_; lean_object* v_tail_439_; lean_object* v___x_441_; uint8_t v_isShared_442_; uint8_t v_isSharedCheck_449_; 
v_head_438_ = lean_ctor_get(v_x_437_, 0);
v_tail_439_ = lean_ctor_get(v_x_437_, 1);
v_isSharedCheck_449_ = !lean_is_exclusive(v_x_437_);
if (v_isSharedCheck_449_ == 0)
{
v___x_441_ = v_x_437_;
v_isShared_442_ = v_isSharedCheck_449_;
goto v_resetjp_440_;
}
else
{
lean_inc(v_tail_439_);
lean_inc(v_head_438_);
lean_dec(v_x_437_);
v___x_441_ = lean_box(0);
v_isShared_442_ = v_isSharedCheck_449_;
goto v_resetjp_440_;
}
v_resetjp_440_:
{
lean_object* v___x_444_; 
lean_inc(v_x_435_);
if (v_isShared_442_ == 0)
{
lean_ctor_set_tag(v___x_441_, 5);
lean_ctor_set(v___x_441_, 1, v_x_435_);
lean_ctor_set(v___x_441_, 0, v_x_436_);
v___x_444_ = v___x_441_;
goto v_reusejp_443_;
}
else
{
lean_object* v_reuseFailAlloc_448_; 
v_reuseFailAlloc_448_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_448_, 0, v_x_436_);
lean_ctor_set(v_reuseFailAlloc_448_, 1, v_x_435_);
v___x_444_ = v_reuseFailAlloc_448_;
goto v_reusejp_443_;
}
v_reusejp_443_:
{
lean_object* v___x_445_; lean_object* v___x_446_; 
v___x_445_ = lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg(v_head_438_);
v___x_446_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_446_, 0, v___x_444_);
lean_ctor_set(v___x_446_, 1, v___x_445_);
v_x_436_ = v___x_446_;
v_x_437_ = v_tail_439_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5_spec__8_spec__11(lean_object* v_x_450_, lean_object* v_x_451_, lean_object* v_x_452_){
_start:
{
if (lean_obj_tag(v_x_452_) == 0)
{
lean_dec(v_x_450_);
return v_x_451_;
}
else
{
lean_object* v_head_453_; lean_object* v_tail_454_; lean_object* v___x_456_; uint8_t v_isShared_457_; uint8_t v_isSharedCheck_464_; 
v_head_453_ = lean_ctor_get(v_x_452_, 0);
v_tail_454_ = lean_ctor_get(v_x_452_, 1);
v_isSharedCheck_464_ = !lean_is_exclusive(v_x_452_);
if (v_isSharedCheck_464_ == 0)
{
v___x_456_ = v_x_452_;
v_isShared_457_ = v_isSharedCheck_464_;
goto v_resetjp_455_;
}
else
{
lean_inc(v_tail_454_);
lean_inc(v_head_453_);
lean_dec(v_x_452_);
v___x_456_ = lean_box(0);
v_isShared_457_ = v_isSharedCheck_464_;
goto v_resetjp_455_;
}
v_resetjp_455_:
{
lean_object* v___x_459_; 
lean_inc(v_x_450_);
if (v_isShared_457_ == 0)
{
lean_ctor_set_tag(v___x_456_, 5);
lean_ctor_set(v___x_456_, 1, v_x_450_);
lean_ctor_set(v___x_456_, 0, v_x_451_);
v___x_459_ = v___x_456_;
goto v_reusejp_458_;
}
else
{
lean_object* v_reuseFailAlloc_463_; 
v_reuseFailAlloc_463_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v_reuseFailAlloc_463_, 0, v_x_451_);
lean_ctor_set(v_reuseFailAlloc_463_, 1, v_x_450_);
v___x_459_ = v_reuseFailAlloc_463_;
goto v_reusejp_458_;
}
v_reusejp_458_:
{
lean_object* v___x_460_; lean_object* v___x_461_; lean_object* v___x_462_; 
v___x_460_ = lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg(v_head_453_);
v___x_461_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_461_, 0, v___x_459_);
lean_ctor_set(v___x_461_, 1, v___x_460_);
v___x_462_ = lp_BioCompiler_List_foldl___at___00List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5_spec__8_spec__11_spec__14(v_x_450_, v___x_461_, v_tail_454_);
return v___x_462_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5_spec__8(lean_object* v_x_465_, lean_object* v_x_466_){
_start:
{
if (lean_obj_tag(v_x_465_) == 0)
{
lean_object* v___x_467_; 
lean_dec(v_x_466_);
v___x_467_ = lean_box(0);
return v___x_467_;
}
else
{
lean_object* v_tail_468_; 
v_tail_468_ = lean_ctor_get(v_x_465_, 1);
if (lean_obj_tag(v_tail_468_) == 0)
{
lean_object* v_head_469_; lean_object* v___x_470_; 
lean_dec(v_x_466_);
v_head_469_ = lean_ctor_get(v_x_465_, 0);
lean_inc(v_head_469_);
lean_dec_ref(v_x_465_);
v___x_470_ = lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg(v_head_469_);
return v___x_470_;
}
else
{
lean_object* v_head_471_; lean_object* v___x_472_; lean_object* v___x_473_; 
lean_inc(v_tail_468_);
v_head_471_ = lean_ctor_get(v_x_465_, 0);
lean_inc(v_head_471_);
lean_dec_ref(v_x_465_);
v___x_472_ = lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg(v_head_471_);
v___x_473_ = lp_BioCompiler_List_foldl___at___00Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5_spec__8_spec__11(v_x_466_, v___x_472_, v_tail_468_);
return v___x_473_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5___redArg(lean_object* v_a_474_){
_start:
{
if (lean_obj_tag(v_a_474_) == 0)
{
lean_object* v___x_475_; 
v___x_475_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__1));
return v___x_475_;
}
else
{
lean_object* v___x_476_; lean_object* v___x_477_; lean_object* v___x_478_; lean_object* v___x_479_; lean_object* v___x_480_; lean_object* v___x_481_; lean_object* v___x_482_; lean_object* v___x_483_; uint8_t v___x_484_; lean_object* v___x_485_; 
v___x_476_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__3));
v___x_477_ = lp_BioCompiler_Std_Format_joinSep___at___00List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5_spec__8(v_a_474_, v___x_476_);
v___x_478_ = lean_obj_once(&lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6, &lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6_once, _init_lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__6);
v___x_479_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__7));
v___x_480_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_480_, 0, v___x_479_);
lean_ctor_set(v___x_480_, 1, v___x_477_);
v___x_481_ = ((lean_object*)(lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg___closed__8));
v___x_482_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_482_, 0, v___x_480_);
lean_ctor_set(v___x_482_, 1, v___x_481_);
v___x_483_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_483_, 0, v___x_478_);
lean_ctor_set(v___x_483_, 1, v___x_482_);
v___x_484_ = 0;
v___x_485_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_485_, 0, v___x_483_);
lean_ctor_set_uint8(v___x_485_, sizeof(void*)*1, v___x_484_);
return v___x_485_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3(lean_object* v_x_486_, lean_object* v_x_487_){
_start:
{
if (lean_obj_tag(v_x_486_) == 0)
{
lean_object* v___x_488_; 
v___x_488_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__1));
return v___x_488_;
}
else
{
lean_object* v_val_489_; lean_object* v___x_490_; lean_object* v___x_491_; lean_object* v___x_492_; lean_object* v___x_493_; 
v_val_489_ = lean_ctor_get(v_x_486_, 0);
lean_inc(v_val_489_);
lean_dec_ref(v_x_486_);
v___x_490_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__3));
v___x_491_ = lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5___redArg(v_val_489_);
v___x_492_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_492_, 0, v___x_490_);
lean_ctor_set(v___x_492_, 1, v___x_491_);
v___x_493_ = l_Repr_addAppParen(v___x_492_, v_x_487_);
return v___x_493_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3___boxed(lean_object* v_x_494_, lean_object* v_x_495_){
_start:
{
lean_object* v_res_496_; 
v_res_496_ = lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3(v_x_494_, v_x_495_);
lean_dec(v_x_495_);
return v_res_496_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__6(void){
_start:
{
lean_object* v___x_509_; lean_object* v___x_510_; 
v___x_509_ = lean_unsigned_to_nat(13u);
v___x_510_ = lean_nat_to_int(v___x_509_);
return v___x_510_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg(lean_object* v_x_517_){
_start:
{
lean_object* v_atomCoordinates_518_; lean_object* v_meanPLDDT_519_; lean_object* v_paeMatrix_520_; lean_object* v_ptmSites_521_; lean_object* v___x_522_; lean_object* v___x_523_; lean_object* v___x_524_; lean_object* v___x_525_; lean_object* v___x_526_; lean_object* v___x_527_; uint8_t v___x_528_; lean_object* v___x_529_; lean_object* v___x_530_; lean_object* v___x_531_; lean_object* v___x_532_; lean_object* v___x_533_; lean_object* v___x_534_; lean_object* v___x_535_; lean_object* v___x_536_; lean_object* v___x_537_; lean_object* v___x_538_; lean_object* v___x_539_; lean_object* v___x_540_; lean_object* v___x_541_; lean_object* v___x_542_; lean_object* v___x_543_; lean_object* v___x_544_; lean_object* v___x_545_; lean_object* v___x_546_; lean_object* v___x_547_; lean_object* v___x_548_; lean_object* v___x_549_; lean_object* v___x_550_; lean_object* v___x_551_; lean_object* v___x_552_; lean_object* v___x_553_; lean_object* v___x_554_; lean_object* v___x_555_; lean_object* v___x_556_; lean_object* v___x_557_; lean_object* v___x_558_; lean_object* v___x_559_; lean_object* v___x_560_; lean_object* v___x_561_; lean_object* v___x_562_; lean_object* v___x_563_; lean_object* v___x_564_; lean_object* v___x_565_; lean_object* v___x_566_; lean_object* v___x_567_; lean_object* v___x_568_; 
v_atomCoordinates_518_ = lean_ctor_get(v_x_517_, 0);
lean_inc(v_atomCoordinates_518_);
v_meanPLDDT_519_ = lean_ctor_get(v_x_517_, 1);
lean_inc(v_meanPLDDT_519_);
v_paeMatrix_520_ = lean_ctor_get(v_x_517_, 2);
lean_inc(v_paeMatrix_520_);
v_ptmSites_521_ = lean_ctor_get(v_x_517_, 3);
lean_inc(v_ptmSites_521_);
lean_dec_ref(v_x_517_);
v___x_522_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5));
v___x_523_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__3));
v___x_524_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__7);
v___x_525_ = lean_unsigned_to_nat(0u);
v___x_526_ = lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0(v_atomCoordinates_518_, v___x_525_);
v___x_527_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_527_, 0, v___x_524_);
lean_ctor_set(v___x_527_, 1, v___x_526_);
v___x_528_ = 0;
v___x_529_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_529_, 0, v___x_527_);
lean_ctor_set_uint8(v___x_529_, sizeof(void*)*1, v___x_528_);
v___x_530_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_530_, 0, v___x_523_);
lean_ctor_set(v___x_530_, 1, v___x_529_);
v___x_531_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__9));
v___x_532_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_532_, 0, v___x_530_);
lean_ctor_set(v___x_532_, 1, v___x_531_);
v___x_533_ = lean_box(1);
v___x_534_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_534_, 0, v___x_532_);
lean_ctor_set(v___x_534_, 1, v___x_533_);
v___x_535_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__5));
v___x_536_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_536_, 0, v___x_534_);
lean_ctor_set(v___x_536_, 1, v___x_535_);
v___x_537_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_537_, 0, v___x_536_);
lean_ctor_set(v___x_537_, 1, v___x_522_);
v___x_538_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__6, &lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__6);
v___x_539_ = lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1(v_meanPLDDT_519_, v___x_525_);
v___x_540_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_540_, 0, v___x_538_);
lean_ctor_set(v___x_540_, 1, v___x_539_);
v___x_541_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_541_, 0, v___x_540_);
lean_ctor_set_uint8(v___x_541_, sizeof(void*)*1, v___x_528_);
v___x_542_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_542_, 0, v___x_537_);
lean_ctor_set(v___x_542_, 1, v___x_541_);
v___x_543_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_543_, 0, v___x_542_);
lean_ctor_set(v___x_543_, 1, v___x_531_);
v___x_544_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_544_, 0, v___x_543_);
lean_ctor_set(v___x_544_, 1, v___x_533_);
v___x_545_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__8));
v___x_546_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_546_, 0, v___x_544_);
lean_ctor_set(v___x_546_, 1, v___x_545_);
v___x_547_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_547_, 0, v___x_546_);
lean_ctor_set(v___x_547_, 1, v___x_522_);
v___x_548_ = lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2(v_paeMatrix_520_, v___x_525_);
v___x_549_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_549_, 0, v___x_538_);
lean_ctor_set(v___x_549_, 1, v___x_548_);
v___x_550_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_550_, 0, v___x_549_);
lean_ctor_set_uint8(v___x_550_, sizeof(void*)*1, v___x_528_);
v___x_551_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_551_, 0, v___x_547_);
lean_ctor_set(v___x_551_, 1, v___x_550_);
v___x_552_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_552_, 0, v___x_551_);
lean_ctor_set(v___x_552_, 1, v___x_531_);
v___x_553_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_553_, 0, v___x_552_);
lean_ctor_set(v___x_553_, 1, v___x_533_);
v___x_554_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg___closed__10));
v___x_555_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_555_, 0, v___x_553_);
lean_ctor_set(v___x_555_, 1, v___x_554_);
v___x_556_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_556_, 0, v___x_555_);
lean_ctor_set(v___x_556_, 1, v___x_522_);
v___x_557_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7);
v___x_558_ = lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3(v_ptmSites_521_, v___x_525_);
v___x_559_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_559_, 0, v___x_557_);
lean_ctor_set(v___x_559_, 1, v___x_558_);
v___x_560_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_560_, 0, v___x_559_);
lean_ctor_set_uint8(v___x_560_, sizeof(void*)*1, v___x_528_);
v___x_561_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_561_, 0, v___x_556_);
lean_ctor_set(v___x_561_, 1, v___x_560_);
v___x_562_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19);
v___x_563_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__20));
v___x_564_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_564_, 0, v___x_563_);
lean_ctor_set(v___x_564_, 1, v___x_561_);
v___x_565_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__21));
v___x_566_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_566_, 0, v___x_564_);
lean_ctor_set(v___x_566_, 1, v___x_565_);
v___x_567_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_567_, 0, v___x_562_);
lean_ctor_set(v___x_567_, 1, v___x_566_);
v___x_568_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_568_, 0, v___x_567_);
lean_ctor_set_uint8(v___x_568_, sizeof(void*)*1, v___x_528_);
return v___x_568_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr(lean_object* v_x_569_, lean_object* v_prec_570_){
_start:
{
lean_object* v___x_571_; 
v___x_571_ = lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg(v_x_569_);
return v___x_571_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___boxed(lean_object* v_x_572_, lean_object* v_prec_573_){
_start:
{
lean_object* v_res_574_; 
v_res_574_ = lp_BioCompiler_BioCompiler_instReprSLOTValues_repr(v_x_572_, v_prec_573_);
lean_dec(v_prec_573_);
return v_res_574_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0(lean_object* v_a_575_, lean_object* v_n_576_){
_start:
{
lean_object* v___x_577_; 
v___x_577_ = lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0___redArg(v_a_575_);
return v___x_577_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0___boxed(lean_object* v_a_578_, lean_object* v_n_579_){
_start:
{
lean_object* v_res_580_; 
v_res_580_ = lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__0_spec__0(v_a_578_, v_n_579_);
lean_dec(v_n_579_);
return v_res_580_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3(lean_object* v_a_581_, lean_object* v_n_582_){
_start:
{
lean_object* v___x_583_; 
v___x_583_ = lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___redArg(v_a_581_);
return v___x_583_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3___boxed(lean_object* v_a_584_, lean_object* v_n_585_){
_start:
{
lean_object* v_res_586_; 
v_res_586_ = lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__2_spec__3(v_a_584_, v_n_585_);
lean_dec(v_n_585_);
return v_res_586_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5(lean_object* v_a_587_, lean_object* v_n_588_){
_start:
{
lean_object* v___x_589_; 
v___x_589_ = lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5___redArg(v_a_587_);
return v___x_589_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5___boxed(lean_object* v_a_590_, lean_object* v_n_591_){
_start:
{
lean_object* v_res_592_; 
v_res_592_ = lp_BioCompiler_List_repr___at___00Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__3_spec__5(v_a_590_, v_n_591_);
lean_dec(v_n_591_);
return v_res_592_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg(lean_object* v_x_610_){
_start:
{
lean_object* v_sequence_611_; lean_object* v_slots_612_; lean_object* v___x_614_; uint8_t v_isShared_615_; uint8_t v_isSharedCheck_645_; 
v_sequence_611_ = lean_ctor_get(v_x_610_, 0);
v_slots_612_ = lean_ctor_get(v_x_610_, 1);
v_isSharedCheck_645_ = !lean_is_exclusive(v_x_610_);
if (v_isSharedCheck_645_ == 0)
{
v___x_614_ = v_x_610_;
v_isShared_615_ = v_isSharedCheck_645_;
goto v_resetjp_613_;
}
else
{
lean_inc(v_slots_612_);
lean_inc(v_sequence_611_);
lean_dec(v_x_610_);
v___x_614_ = lean_box(0);
v_isShared_615_ = v_isSharedCheck_645_;
goto v_resetjp_613_;
}
v_resetjp_613_:
{
lean_object* v___x_616_; lean_object* v___x_617_; lean_object* v___x_618_; lean_object* v___x_619_; lean_object* v___x_621_; 
v___x_616_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__5));
v___x_617_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__3));
v___x_618_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__7);
v___x_619_ = lp_BioCompiler_List_repr___at___00BioCompiler_instReprSpliceIsoform_repr_spec__0___redArg(v_sequence_611_);
if (v_isShared_615_ == 0)
{
lean_ctor_set_tag(v___x_614_, 4);
lean_ctor_set(v___x_614_, 1, v___x_619_);
lean_ctor_set(v___x_614_, 0, v___x_618_);
v___x_621_ = v___x_614_;
goto v_reusejp_620_;
}
else
{
lean_object* v_reuseFailAlloc_644_; 
v_reuseFailAlloc_644_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v_reuseFailAlloc_644_, 0, v___x_618_);
lean_ctor_set(v_reuseFailAlloc_644_, 1, v___x_619_);
v___x_621_ = v_reuseFailAlloc_644_;
goto v_reusejp_620_;
}
v_reusejp_620_:
{
uint8_t v___x_622_; lean_object* v___x_623_; lean_object* v___x_624_; lean_object* v___x_625_; lean_object* v___x_626_; lean_object* v___x_627_; lean_object* v___x_628_; lean_object* v___x_629_; lean_object* v___x_630_; lean_object* v___x_631_; lean_object* v___x_632_; lean_object* v___x_633_; lean_object* v___x_634_; lean_object* v___x_635_; lean_object* v___x_636_; lean_object* v___x_637_; lean_object* v___x_638_; lean_object* v___x_639_; lean_object* v___x_640_; lean_object* v___x_641_; lean_object* v___x_642_; lean_object* v___x_643_; 
v___x_622_ = 0;
v___x_623_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_623_, 0, v___x_621_);
lean_ctor_set_uint8(v___x_623_, sizeof(void*)*1, v___x_622_);
v___x_624_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_624_, 0, v___x_617_);
lean_ctor_set(v___x_624_, 1, v___x_623_);
v___x_625_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__9));
v___x_626_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_626_, 0, v___x_624_);
lean_ctor_set(v___x_626_, 1, v___x_625_);
v___x_627_ = lean_box(1);
v___x_628_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_628_, 0, v___x_626_);
lean_ctor_set(v___x_628_, 1, v___x_627_);
v___x_629_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg___closed__5));
v___x_630_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_630_, 0, v___x_628_);
lean_ctor_set(v___x_630_, 1, v___x_629_);
v___x_631_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_631_, 0, v___x_630_);
lean_ctor_set(v___x_631_, 1, v___x_616_);
v___x_632_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__10, &lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__10_once, _init_lp_BioCompiler_BioCompiler_instReprPTMSiteEntry_repr___redArg___closed__10);
v___x_633_ = lp_BioCompiler_BioCompiler_instReprSLOTValues_repr___redArg(v_slots_612_);
v___x_634_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_634_, 0, v___x_632_);
lean_ctor_set(v___x_634_, 1, v___x_633_);
v___x_635_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_635_, 0, v___x_634_);
lean_ctor_set_uint8(v___x_635_, sizeof(void*)*1, v___x_622_);
v___x_636_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_636_, 0, v___x_631_);
lean_ctor_set(v___x_636_, 1, v___x_635_);
v___x_637_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19, &lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19_once, _init_lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__19);
v___x_638_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__20));
v___x_639_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_639_, 0, v___x_638_);
lean_ctor_set(v___x_639_, 1, v___x_636_);
v___x_640_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprAtomCoordinate_repr___redArg___closed__21));
v___x_641_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_641_, 0, v___x_639_);
lean_ctor_set(v___x_641_, 1, v___x_640_);
v___x_642_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_642_, 0, v___x_637_);
lean_ctor_set(v___x_642_, 1, v___x_641_);
v___x_643_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_643_, 0, v___x_642_);
lean_ctor_set_uint8(v___x_643_, sizeof(void*)*1, v___x_622_);
return v___x_643_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr(lean_object* v_x_646_, lean_object* v_prec_647_){
_start:
{
lean_object* v___x_648_; 
v___x_648_ = lp_BioCompiler_BioCompiler_instReprIRRecord_repr___redArg(v_x_646_);
return v___x_648_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprIRRecord_repr___boxed(lean_object* v_x_649_, lean_object* v_prec_650_){
_start:
{
lean_object* v_res_651_; 
v_res_651_ = lp_BioCompiler_BioCompiler_instReprIRRecord_repr(v_x_649_, v_prec_650_);
lean_dec(v_prec_650_);
return v_res_651_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_isCorePredicate(lean_object* v_x_654_){
_start:
{
switch(lean_obj_tag(v_x_654_))
{
case 1:
{
uint8_t v___x_655_; 
v___x_655_ = 1;
return v___x_655_;
}
case 6:
{
uint8_t v___x_656_; 
v___x_656_ = 1;
return v___x_656_;
}
case 7:
{
uint8_t v___x_657_; 
v___x_657_ = 1;
return v___x_657_;
}
default: 
{
uint8_t v___x_658_; 
v___x_658_ = 1;
return v___x_658_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_isCorePredicate___boxed(lean_object* v_x_659_){
_start:
{
uint8_t v_res_660_; lean_object* v_r_661_; 
v_res_660_ = lp_BioCompiler_BioCompiler_isCorePredicate(v_x_659_);
lean_dec(v_x_659_);
v_r_661_ = lean_box(v_res_660_);
return v_r_661_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__3(void){
_start:
{
lean_object* v___x_668_; lean_object* v___x_669_; 
v___x_668_ = lean_unsigned_to_nat(2u);
v___x_669_ = lean_nat_to_int(v___x_668_);
return v___x_669_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__4(void){
_start:
{
lean_object* v___x_670_; lean_object* v___x_671_; 
v___x_670_ = lean_unsigned_to_nat(1u);
v___x_671_ = lean_nat_to_int(v___x_670_);
return v___x_671_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr(lean_object* v_x_672_, lean_object* v_prec_673_){
_start:
{
lean_object* v___y_675_; lean_object* v___y_676_; lean_object* v___y_677_; lean_object* v___y_684_; lean_object* v___x_706_; uint8_t v___x_707_; 
v___x_706_ = lean_unsigned_to_nat(1024u);
v___x_707_ = lean_nat_dec_le(v___x_706_, v_prec_673_);
if (v___x_707_ == 0)
{
lean_object* v___x_708_; 
v___x_708_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__3, &lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__3_once, _init_lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__3);
v___y_684_ = v___x_708_;
goto v___jp_683_;
}
else
{
lean_object* v___x_709_; 
v___x_709_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__4, &lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__4_once, _init_lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__4);
v___y_684_ = v___x_709_;
goto v___jp_683_;
}
v___jp_674_:
{
lean_object* v___x_678_; lean_object* v___x_679_; uint8_t v___x_680_; lean_object* v___x_681_; lean_object* v___x_682_; 
lean_inc(v___y_675_);
v___x_678_ = lean_alloc_ctor(5, 2, 0);
lean_ctor_set(v___x_678_, 0, v___y_675_);
lean_ctor_set(v___x_678_, 1, v___y_677_);
lean_inc(v___y_676_);
v___x_679_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_679_, 0, v___y_676_);
lean_ctor_set(v___x_679_, 1, v___x_678_);
v___x_680_ = 0;
v___x_681_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_681_, 0, v___x_679_);
lean_ctor_set_uint8(v___x_681_, sizeof(void*)*1, v___x_680_);
v___x_682_ = l_Repr_addAppParen(v___x_681_, v_prec_673_);
return v___x_682_;
}
v___jp_683_:
{
lean_object* v_num_685_; lean_object* v_den_686_; lean_object* v___x_687_; lean_object* v___x_688_; uint8_t v___x_689_; 
v_num_685_ = lean_ctor_get(v_x_672_, 0);
lean_inc(v_num_685_);
v_den_686_ = lean_ctor_get(v_x_672_, 1);
lean_inc(v_den_686_);
lean_dec_ref(v_x_672_);
v___x_687_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___closed__2));
v___x_688_ = lean_unsigned_to_nat(1u);
v___x_689_ = lean_nat_dec_eq(v_den_686_, v___x_688_);
if (v___x_689_ == 0)
{
lean_object* v___x_690_; lean_object* v___x_691_; lean_object* v___x_692_; lean_object* v___x_693_; lean_object* v___x_694_; lean_object* v___x_695_; lean_object* v___x_696_; lean_object* v___x_697_; 
v___x_690_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__4));
v___x_691_ = l_Int_repr(v_num_685_);
lean_dec(v_num_685_);
v___x_692_ = lean_string_append(v___x_690_, v___x_691_);
lean_dec_ref(v___x_691_);
v___x_693_ = ((lean_object*)(lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__5));
v___x_694_ = lean_string_append(v___x_692_, v___x_693_);
v___x_695_ = l_Nat_reprFast(v_den_686_);
v___x_696_ = lean_string_append(v___x_694_, v___x_695_);
lean_dec_ref(v___x_695_);
v___x_697_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_697_, 0, v___x_696_);
v___y_675_ = v___x_687_;
v___y_676_ = v___y_684_;
v___y_677_ = v___x_697_;
goto v___jp_674_;
}
else
{
lean_object* v___x_698_; lean_object* v___x_699_; uint8_t v___x_700_; 
lean_dec(v_den_686_);
v___x_698_ = lean_unsigned_to_nat(0u);
v___x_699_ = lean_obj_once(&lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__6, &lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__6_once, _init_lp_BioCompiler_Option_repr___at___00BioCompiler_instReprSLOTValues_repr_spec__1___closed__6);
v___x_700_ = lean_int_dec_lt(v_num_685_, v___x_699_);
if (v___x_700_ == 0)
{
lean_object* v___x_701_; lean_object* v___x_702_; 
v___x_701_ = l_Int_repr(v_num_685_);
lean_dec(v_num_685_);
v___x_702_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_702_, 0, v___x_701_);
v___y_675_ = v___x_687_;
v___y_676_ = v___y_684_;
v___y_677_ = v___x_702_;
goto v___jp_674_;
}
else
{
lean_object* v___x_703_; lean_object* v___x_704_; lean_object* v___x_705_; 
v___x_703_ = l_Int_repr(v_num_685_);
lean_dec(v_num_685_);
v___x_704_ = lean_alloc_ctor(3, 1, 0);
lean_ctor_set(v___x_704_, 0, v___x_703_);
v___x_705_ = l_Repr_addAppParen(v___x_704_, v___x_698_);
v___y_675_ = v___x_687_;
v___y_676_ = v___y_684_;
v___y_677_ = v___x_705_;
goto v___jp_674_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr___boxed(lean_object* v_x_710_, lean_object* v_prec_711_){
_start:
{
lean_object* v_res_712_; 
v_res_712_ = lp_BioCompiler_BioCompiler_instReprFFIDependentPredicate_repr(v_x_710_, v_prec_711_);
lean_dec(v_prec_711_);
return v_res_712_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluateFFIDependent(lean_object* v_x_715_, lean_object* v_x_716_){
_start:
{
lean_object* v_meanPLDDT_717_; 
v_meanPLDDT_717_ = lean_ctor_get(v_x_716_, 1);
lean_inc(v_meanPLDDT_717_);
lean_dec_ref(v_x_716_);
if (lean_obj_tag(v_meanPLDDT_717_) == 0)
{
uint8_t v___x_718_; 
lean_dec_ref(v_x_715_);
v___x_718_ = 2;
return v___x_718_;
}
else
{
lean_object* v_val_719_; uint8_t v___x_720_; 
v_val_719_ = lean_ctor_get(v_meanPLDDT_717_, 0);
lean_inc(v_val_719_);
lean_dec_ref(v_meanPLDDT_717_);
v___x_720_ = l_Rat_instDecidableLe(v_x_715_, v_val_719_);
if (v___x_720_ == 0)
{
uint8_t v___x_721_; 
v___x_721_ = 1;
return v___x_721_;
}
else
{
uint8_t v___x_722_; 
v___x_722_ = 2;
return v___x_722_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluateFFIDependent___boxed(lean_object* v_x_723_, lean_object* v_x_724_){
_start:
{
uint8_t v_res_725_; lean_object* v_r_726_; 
v_res_725_ = lp_BioCompiler_BioCompiler_evaluateFFIDependent(v_x_723_, v_x_724_);
v_r_726_ = lean_box(v_res_725_);
return v_r_726_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_SLOTIndependence_0__BioCompiler_evaluateFFIDependent_match__3_splitter___redArg(lean_object* v_x_727_, lean_object* v_x_728_, lean_object* v_h__1_729_){
_start:
{
lean_object* v___x_730_; 
v___x_730_ = lean_apply_2(v_h__1_729_, v_x_727_, v_x_728_);
return v___x_730_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_SLOTIndependence_0__BioCompiler_evaluateFFIDependent_match__3_splitter(lean_object* v_motive_731_, lean_object* v_x_732_, lean_object* v_x_733_, lean_object* v_h__1_734_){
_start:
{
lean_object* v___x_735_; 
v___x_735_ = lean_apply_2(v_h__1_734_, v_x_732_, v_x_733_);
return v___x_735_;
}
}
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_ThreeValued(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Sequence(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_NDFST(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Scanners(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_TypeSystem(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Compositional(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Certificates(uint8_t builtin);
static bool _G_initialized = false;
LEAN_EXPORT lean_object* initialize_BioCompiler_BioCompiler_SLOTIndependence(uint8_t builtin) {
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
res = initialize_BioCompiler_BioCompiler_TypeSystem(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_BioCompiler_BioCompiler_Compositional(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_BioCompiler_BioCompiler_Certificates(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
