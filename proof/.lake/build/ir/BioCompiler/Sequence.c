// Lean compiler output
// Module: BioCompiler.Sequence
// Imports: public import Init public meta import Init
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
lean_object* l_List_lengthTR___redArg(lean_object*);
lean_object* lean_nat_add(lean_object*, lean_object*);
uint8_t lean_nat_dec_le(lean_object*, lean_object*);
uint8_t lean_nat_dec_eq(lean_object*, lean_object*);
lean_object* l_List_drop___redArg(lean_object*, lean_object*);
lean_object* lean_mk_empty_array_with_capacity(lean_object*);
lean_object* l___private_Init_Data_List_Impl_0__List_takeTR_go___redArg(lean_object*, lean_object*, lean_object*, lean_object*);
uint8_t l_instDecidableEqList___redArg(lean_object*, lean_object*, lean_object*);
lean_object* lean_nat_mod(lean_object*, lean_object*);
lean_object* l_Nat_cast___at___00Dyadic_toRat_spec__0(lean_object*);
uint8_t l_List_any___redArg(lean_object*, lean_object*);
uint8_t l_List_all___redArg(lean_object*, lean_object*);
uint8_t lean_nat_dec_lt(lean_object*, lean_object*);
lean_object* lean_nat_sub(lean_object*, lean_object*);
lean_object* l_List_range(lean_object*);
lean_object* l_Repr_addAppParen(lean_object*, lean_object*);
lean_object* lean_nat_to_int(lean_object*);
uint8_t lean_nat_dec_le(lean_object*, lean_object*);
lean_object* lean_nat_div(lean_object*, lean_object*);
lean_object* l_List_reverse___redArg(lean_object*);
lean_object* lean_nat_mul(lean_object*, lean_object*);
lean_object* l_Rat_div(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx(uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_toCtorIdx(uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_toCtorIdx___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorElim___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorElim___redArg___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorElim(lean_object*, lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorElim___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_A_elim___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_A_elim___redArg___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_A_elim(lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_A_elim___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_C_elim___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_C_elim___redArg___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_C_elim(lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_C_elim___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_G_elim___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_G_elim___redArg___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_G_elim(lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_G_elim___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_T_elim___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_T_elim___redArg___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_T_elim(lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_T_elim___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Nucleotide_ofNat(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ofNat___boxed(lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instDecidableEqNucleotide(uint8_t, uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instDecidableEqNucleotide___boxed(lean_object*, lean_object*);
static const lean_string_object lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 25, .m_capacity = 25, .m_length = 24, .m_data = "BioCompiler.Nucleotide.A"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__1_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 25, .m_capacity = 25, .m_length = 24, .m_data = "BioCompiler.Nucleotide.C"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__2_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 25, .m_capacity = 25, .m_length = 24, .m_data = "BioCompiler.Nucleotide.G"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__5_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 25, .m_capacity = 25, .m_length = 24, .m_data = "BioCompiler.Nucleotide.T"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__6 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__6_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__7_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__6_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__7 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__7_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr(uint8_t, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprNucleotide___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprNucleotide_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprNucleotide___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instBEqNucleotide_beq(uint8_t, uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instBEqNucleotide_beq___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instBEqNucleotide___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instBEqNucleotide_beq___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instBEqNucleotide___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instBEqNucleotide___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instBEqNucleotide = (const lean_object*)&lp_BioCompiler_BioCompiler_instBEqNucleotide___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instInhabitedNucleotide_default;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instInhabitedNucleotide;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_empty;
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__0___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__1(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__1___boxed(lean_object*, lean_object*);
static lean_once_cell_t lp_BioCompiler_BioCompiler_Sequence_gcContent___closed__0_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_Sequence_gcContent___closed__0;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_gcContent(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_gcContent___boxed(lean_object*);
static const lean_array_object lp_BioCompiler_BioCompiler_Sequence_matchesAt___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_array_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 246}, .m_size = 0, .m_capacity = 0, .m_data = {}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_matchesAt___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_matchesAt___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_matchesAt(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_matchesAt___boxed(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_containsPattern___lam__0(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_containsPattern___lam__0___boxed(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_containsPattern(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_containsPattern___boxed(lean_object*, lean_object*);
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)(((size_t)(0) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__1_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(3) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__1_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(2) << 1) | 1)),((lean_object*)(((size_t)(0) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__3_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(0) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__3_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(3) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__5_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__6_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(2) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__6 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__6_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__7_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)(((size_t)(3) << 1) | 1)),((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__6_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__7 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__7_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__8_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__7_value),((lean_object*)(((size_t)(0) << 1) | 1))}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__8 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__8_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__9_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__5_value),((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__8_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__9 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__9_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__10_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*2 + 0, .m_other = 2, .m_tag = 1}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__2_value),((lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__9_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__10 = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__10_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_Sequence_stopCodons = (const lean_object*)&lp_BioCompiler_BioCompiler_Sequence_stopCodons___closed__10_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_isStopCodon___lam__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_isStopCodon___lam__0___boxed(lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_isStopCodon(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_isStopCodon___boxed(lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop___lam__0(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop___lam__0___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_Sequence_hasPrematureStop_spec__0(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_Sequence_hasPrematureStop_spec__0___boxed(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop___boxed(lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent___lam__0(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent___lam__0___boxed(lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent___boxed(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx(uint8_t v_x_1_){
_start:
{
switch(v_x_1_)
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
default: 
{
lean_object* v___x_5_; 
v___x_5_ = lean_unsigned_to_nat(3u);
return v___x_5_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx___boxed(lean_object* v_x_6_){
_start:
{
uint8_t v_x_boxed_7_; lean_object* v_res_8_; 
v_x_boxed_7_ = lean_unbox(v_x_6_);
v_res_8_ = lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx(v_x_boxed_7_);
return v_res_8_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_toCtorIdx(uint8_t v_x_9_){
_start:
{
lean_object* v___x_10_; 
v___x_10_ = lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx(v_x_9_);
return v___x_10_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_toCtorIdx___boxed(lean_object* v_x_11_){
_start:
{
uint8_t v_x_4__boxed_12_; lean_object* v_res_13_; 
v_x_4__boxed_12_ = lean_unbox(v_x_11_);
v_res_13_ = lp_BioCompiler_BioCompiler_Nucleotide_toCtorIdx(v_x_4__boxed_12_);
return v_res_13_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorElim___redArg(lean_object* v_k_14_){
_start:
{
lean_inc(v_k_14_);
return v_k_14_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorElim___redArg___boxed(lean_object* v_k_15_){
_start:
{
lean_object* v_res_16_; 
v_res_16_ = lp_BioCompiler_BioCompiler_Nucleotide_ctorElim___redArg(v_k_15_);
lean_dec(v_k_15_);
return v_res_16_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorElim(lean_object* v_motive_17_, lean_object* v_ctorIdx_18_, uint8_t v_t_19_, lean_object* v_h_20_, lean_object* v_k_21_){
_start:
{
lean_inc(v_k_21_);
return v_k_21_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ctorElim___boxed(lean_object* v_motive_22_, lean_object* v_ctorIdx_23_, lean_object* v_t_24_, lean_object* v_h_25_, lean_object* v_k_26_){
_start:
{
uint8_t v_t_boxed_27_; lean_object* v_res_28_; 
v_t_boxed_27_ = lean_unbox(v_t_24_);
v_res_28_ = lp_BioCompiler_BioCompiler_Nucleotide_ctorElim(v_motive_22_, v_ctorIdx_23_, v_t_boxed_27_, v_h_25_, v_k_26_);
lean_dec(v_k_26_);
lean_dec(v_ctorIdx_23_);
return v_res_28_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_A_elim___redArg(lean_object* v_A_29_){
_start:
{
lean_inc(v_A_29_);
return v_A_29_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_A_elim___redArg___boxed(lean_object* v_A_30_){
_start:
{
lean_object* v_res_31_; 
v_res_31_ = lp_BioCompiler_BioCompiler_Nucleotide_A_elim___redArg(v_A_30_);
lean_dec(v_A_30_);
return v_res_31_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_A_elim(lean_object* v_motive_32_, uint8_t v_t_33_, lean_object* v_h_34_, lean_object* v_A_35_){
_start:
{
lean_inc(v_A_35_);
return v_A_35_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_A_elim___boxed(lean_object* v_motive_36_, lean_object* v_t_37_, lean_object* v_h_38_, lean_object* v_A_39_){
_start:
{
uint8_t v_t_boxed_40_; lean_object* v_res_41_; 
v_t_boxed_40_ = lean_unbox(v_t_37_);
v_res_41_ = lp_BioCompiler_BioCompiler_Nucleotide_A_elim(v_motive_36_, v_t_boxed_40_, v_h_38_, v_A_39_);
lean_dec(v_A_39_);
return v_res_41_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_C_elim___redArg(lean_object* v_C_42_){
_start:
{
lean_inc(v_C_42_);
return v_C_42_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_C_elim___redArg___boxed(lean_object* v_C_43_){
_start:
{
lean_object* v_res_44_; 
v_res_44_ = lp_BioCompiler_BioCompiler_Nucleotide_C_elim___redArg(v_C_43_);
lean_dec(v_C_43_);
return v_res_44_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_C_elim(lean_object* v_motive_45_, uint8_t v_t_46_, lean_object* v_h_47_, lean_object* v_C_48_){
_start:
{
lean_inc(v_C_48_);
return v_C_48_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_C_elim___boxed(lean_object* v_motive_49_, lean_object* v_t_50_, lean_object* v_h_51_, lean_object* v_C_52_){
_start:
{
uint8_t v_t_boxed_53_; lean_object* v_res_54_; 
v_t_boxed_53_ = lean_unbox(v_t_50_);
v_res_54_ = lp_BioCompiler_BioCompiler_Nucleotide_C_elim(v_motive_49_, v_t_boxed_53_, v_h_51_, v_C_52_);
lean_dec(v_C_52_);
return v_res_54_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_G_elim___redArg(lean_object* v_G_55_){
_start:
{
lean_inc(v_G_55_);
return v_G_55_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_G_elim___redArg___boxed(lean_object* v_G_56_){
_start:
{
lean_object* v_res_57_; 
v_res_57_ = lp_BioCompiler_BioCompiler_Nucleotide_G_elim___redArg(v_G_56_);
lean_dec(v_G_56_);
return v_res_57_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_G_elim(lean_object* v_motive_58_, uint8_t v_t_59_, lean_object* v_h_60_, lean_object* v_G_61_){
_start:
{
lean_inc(v_G_61_);
return v_G_61_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_G_elim___boxed(lean_object* v_motive_62_, lean_object* v_t_63_, lean_object* v_h_64_, lean_object* v_G_65_){
_start:
{
uint8_t v_t_boxed_66_; lean_object* v_res_67_; 
v_t_boxed_66_ = lean_unbox(v_t_63_);
v_res_67_ = lp_BioCompiler_BioCompiler_Nucleotide_G_elim(v_motive_62_, v_t_boxed_66_, v_h_64_, v_G_65_);
lean_dec(v_G_65_);
return v_res_67_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_T_elim___redArg(lean_object* v_T_68_){
_start:
{
lean_inc(v_T_68_);
return v_T_68_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_T_elim___redArg___boxed(lean_object* v_T_69_){
_start:
{
lean_object* v_res_70_; 
v_res_70_ = lp_BioCompiler_BioCompiler_Nucleotide_T_elim___redArg(v_T_69_);
lean_dec(v_T_69_);
return v_res_70_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_T_elim(lean_object* v_motive_71_, uint8_t v_t_72_, lean_object* v_h_73_, lean_object* v_T_74_){
_start:
{
lean_inc(v_T_74_);
return v_T_74_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_T_elim___boxed(lean_object* v_motive_75_, lean_object* v_t_76_, lean_object* v_h_77_, lean_object* v_T_78_){
_start:
{
uint8_t v_t_boxed_79_; lean_object* v_res_80_; 
v_t_boxed_79_ = lean_unbox(v_t_76_);
v_res_80_ = lp_BioCompiler_BioCompiler_Nucleotide_T_elim(v_motive_75_, v_t_boxed_79_, v_h_77_, v_T_78_);
lean_dec(v_T_78_);
return v_res_80_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Nucleotide_ofNat(lean_object* v_n_81_){
_start:
{
lean_object* v___x_82_; uint8_t v___x_83_; 
v___x_82_ = lean_unsigned_to_nat(1u);
v___x_83_ = lean_nat_dec_le(v_n_81_, v___x_82_);
if (v___x_83_ == 0)
{
lean_object* v___x_84_; uint8_t v___x_85_; 
v___x_84_ = lean_unsigned_to_nat(2u);
v___x_85_ = lean_nat_dec_le(v_n_81_, v___x_84_);
if (v___x_85_ == 0)
{
uint8_t v___x_86_; 
v___x_86_ = 3;
return v___x_86_;
}
else
{
uint8_t v___x_87_; 
v___x_87_ = 2;
return v___x_87_;
}
}
else
{
lean_object* v___x_88_; uint8_t v___x_89_; 
v___x_88_ = lean_unsigned_to_nat(0u);
v___x_89_ = lean_nat_dec_le(v_n_81_, v___x_88_);
if (v___x_89_ == 0)
{
uint8_t v___x_90_; 
v___x_90_ = 1;
return v___x_90_;
}
else
{
uint8_t v___x_91_; 
v___x_91_ = 0;
return v___x_91_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Nucleotide_ofNat___boxed(lean_object* v_n_92_){
_start:
{
uint8_t v_res_93_; lean_object* v_r_94_; 
v_res_93_ = lp_BioCompiler_BioCompiler_Nucleotide_ofNat(v_n_92_);
lean_dec(v_n_92_);
v_r_94_ = lean_box(v_res_93_);
return v_r_94_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instDecidableEqNucleotide(uint8_t v_x_95_, uint8_t v_y_96_){
_start:
{
lean_object* v___x_97_; lean_object* v___x_98_; uint8_t v___x_99_; 
v___x_97_ = lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx(v_x_95_);
v___x_98_ = lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx(v_y_96_);
v___x_99_ = lean_nat_dec_eq(v___x_97_, v___x_98_);
lean_dec(v___x_98_);
lean_dec(v___x_97_);
return v___x_99_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instDecidableEqNucleotide___boxed(lean_object* v_x_100_, lean_object* v_y_101_){
_start:
{
uint8_t v_x_13__boxed_102_; uint8_t v_y_14__boxed_103_; uint8_t v_res_104_; lean_object* v_r_105_; 
v_x_13__boxed_102_ = lean_unbox(v_x_100_);
v_y_14__boxed_103_ = lean_unbox(v_y_101_);
v_res_104_ = lp_BioCompiler_BioCompiler_instDecidableEqNucleotide(v_x_13__boxed_102_, v_y_14__boxed_103_);
v_r_105_ = lean_box(v_res_104_);
return v_r_105_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8(void){
_start:
{
lean_object* v___x_118_; lean_object* v___x_119_; 
v___x_118_ = lean_unsigned_to_nat(2u);
v___x_119_ = lean_nat_to_int(v___x_118_);
return v___x_119_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9(void){
_start:
{
lean_object* v___x_120_; lean_object* v___x_121_; 
v___x_120_ = lean_unsigned_to_nat(1u);
v___x_121_ = lean_nat_to_int(v___x_120_);
return v___x_121_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr(uint8_t v_x_122_, lean_object* v_prec_123_){
_start:
{
lean_object* v___y_125_; lean_object* v___y_132_; lean_object* v___y_139_; lean_object* v___y_146_; 
switch(v_x_122_)
{
case 0:
{
lean_object* v___x_152_; uint8_t v___x_153_; 
v___x_152_ = lean_unsigned_to_nat(1024u);
v___x_153_ = lean_nat_dec_le(v___x_152_, v_prec_123_);
if (v___x_153_ == 0)
{
lean_object* v___x_154_; 
v___x_154_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8, &lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8_once, _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8);
v___y_125_ = v___x_154_;
goto v___jp_124_;
}
else
{
lean_object* v___x_155_; 
v___x_155_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9);
v___y_125_ = v___x_155_;
goto v___jp_124_;
}
}
case 1:
{
lean_object* v___x_156_; uint8_t v___x_157_; 
v___x_156_ = lean_unsigned_to_nat(1024u);
v___x_157_ = lean_nat_dec_le(v___x_156_, v_prec_123_);
if (v___x_157_ == 0)
{
lean_object* v___x_158_; 
v___x_158_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8, &lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8_once, _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8);
v___y_132_ = v___x_158_;
goto v___jp_131_;
}
else
{
lean_object* v___x_159_; 
v___x_159_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9);
v___y_132_ = v___x_159_;
goto v___jp_131_;
}
}
case 2:
{
lean_object* v___x_160_; uint8_t v___x_161_; 
v___x_160_ = lean_unsigned_to_nat(1024u);
v___x_161_ = lean_nat_dec_le(v___x_160_, v_prec_123_);
if (v___x_161_ == 0)
{
lean_object* v___x_162_; 
v___x_162_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8, &lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8_once, _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8);
v___y_139_ = v___x_162_;
goto v___jp_138_;
}
else
{
lean_object* v___x_163_; 
v___x_163_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9);
v___y_139_ = v___x_163_;
goto v___jp_138_;
}
}
default: 
{
lean_object* v___x_164_; uint8_t v___x_165_; 
v___x_164_ = lean_unsigned_to_nat(1024u);
v___x_165_ = lean_nat_dec_le(v___x_164_, v_prec_123_);
if (v___x_165_ == 0)
{
lean_object* v___x_166_; 
v___x_166_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8, &lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8_once, _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__8);
v___y_146_ = v___x_166_;
goto v___jp_145_;
}
else
{
lean_object* v___x_167_; 
v___x_167_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9, &lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9_once, _init_lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__9);
v___y_146_ = v___x_167_;
goto v___jp_145_;
}
}
}
v___jp_124_:
{
lean_object* v___x_126_; lean_object* v___x_127_; uint8_t v___x_128_; lean_object* v___x_129_; lean_object* v___x_130_; 
v___x_126_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__1));
lean_inc(v___y_125_);
v___x_127_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_127_, 0, v___y_125_);
lean_ctor_set(v___x_127_, 1, v___x_126_);
v___x_128_ = 0;
v___x_129_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_129_, 0, v___x_127_);
lean_ctor_set_uint8(v___x_129_, sizeof(void*)*1, v___x_128_);
v___x_130_ = l_Repr_addAppParen(v___x_129_, v_prec_123_);
return v___x_130_;
}
v___jp_131_:
{
lean_object* v___x_133_; lean_object* v___x_134_; uint8_t v___x_135_; lean_object* v___x_136_; lean_object* v___x_137_; 
v___x_133_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__3));
lean_inc(v___y_132_);
v___x_134_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_134_, 0, v___y_132_);
lean_ctor_set(v___x_134_, 1, v___x_133_);
v___x_135_ = 0;
v___x_136_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_136_, 0, v___x_134_);
lean_ctor_set_uint8(v___x_136_, sizeof(void*)*1, v___x_135_);
v___x_137_ = l_Repr_addAppParen(v___x_136_, v_prec_123_);
return v___x_137_;
}
v___jp_138_:
{
lean_object* v___x_140_; lean_object* v___x_141_; uint8_t v___x_142_; lean_object* v___x_143_; lean_object* v___x_144_; 
v___x_140_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__5));
lean_inc(v___y_139_);
v___x_141_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_141_, 0, v___y_139_);
lean_ctor_set(v___x_141_, 1, v___x_140_);
v___x_142_ = 0;
v___x_143_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_143_, 0, v___x_141_);
lean_ctor_set_uint8(v___x_143_, sizeof(void*)*1, v___x_142_);
v___x_144_ = l_Repr_addAppParen(v___x_143_, v_prec_123_);
return v___x_144_;
}
v___jp_145_:
{
lean_object* v___x_147_; lean_object* v___x_148_; uint8_t v___x_149_; lean_object* v___x_150_; lean_object* v___x_151_; 
v___x_147_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprNucleotide_repr___closed__7));
lean_inc(v___y_146_);
v___x_148_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_148_, 0, v___y_146_);
lean_ctor_set(v___x_148_, 1, v___x_147_);
v___x_149_ = 0;
v___x_150_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_150_, 0, v___x_148_);
lean_ctor_set_uint8(v___x_150_, sizeof(void*)*1, v___x_149_);
v___x_151_ = l_Repr_addAppParen(v___x_150_, v_prec_123_);
return v___x_151_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprNucleotide_repr___boxed(lean_object* v_x_168_, lean_object* v_prec_169_){
_start:
{
uint8_t v_x_233__boxed_170_; lean_object* v_res_171_; 
v_x_233__boxed_170_ = lean_unbox(v_x_168_);
v_res_171_ = lp_BioCompiler_BioCompiler_instReprNucleotide_repr(v_x_233__boxed_170_, v_prec_169_);
lean_dec(v_prec_169_);
return v_res_171_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instBEqNucleotide_beq(uint8_t v_x_174_, uint8_t v_y_175_){
_start:
{
lean_object* v___x_176_; lean_object* v___x_177_; uint8_t v___x_178_; 
v___x_176_ = lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx(v_x_174_);
v___x_177_ = lp_BioCompiler_BioCompiler_Nucleotide_ctorIdx(v_y_175_);
v___x_178_ = lean_nat_dec_eq(v___x_176_, v___x_177_);
lean_dec(v___x_177_);
lean_dec(v___x_176_);
return v___x_178_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instBEqNucleotide_beq___boxed(lean_object* v_x_179_, lean_object* v_y_180_){
_start:
{
uint8_t v_x_17__boxed_181_; uint8_t v_y_18__boxed_182_; uint8_t v_res_183_; lean_object* v_r_184_; 
v_x_17__boxed_181_ = lean_unbox(v_x_179_);
v_y_18__boxed_182_ = lean_unbox(v_y_180_);
v_res_183_ = lp_BioCompiler_BioCompiler_instBEqNucleotide_beq(v_x_17__boxed_181_, v_y_18__boxed_182_);
v_r_184_ = lean_box(v_res_183_);
return v_r_184_;
}
}
static uint8_t _init_lp_BioCompiler_BioCompiler_instInhabitedNucleotide_default(void){
_start:
{
uint8_t v___x_187_; 
v___x_187_ = 0;
return v___x_187_;
}
}
static uint8_t _init_lp_BioCompiler_BioCompiler_instInhabitedNucleotide(void){
_start:
{
uint8_t v___x_188_; 
v___x_188_ = 0;
return v___x_188_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_Sequence_empty(void){
_start:
{
lean_object* v___x_189_; 
v___x_189_ = lean_box(0);
return v___x_189_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__0(lean_object* v_a_190_, lean_object* v_a_191_){
_start:
{
if (lean_obj_tag(v_a_190_) == 0)
{
return v_a_191_;
}
else
{
lean_object* v_head_192_; lean_object* v_tail_193_; uint8_t v___x_194_; uint8_t v___x_195_; uint8_t v___x_196_; 
v_head_192_ = lean_ctor_get(v_a_190_, 0);
v_tail_193_ = lean_ctor_get(v_a_190_, 1);
v___x_194_ = 2;
v___x_195_ = lean_unbox(v_head_192_);
v___x_196_ = lp_BioCompiler_BioCompiler_instBEqNucleotide_beq(v___x_195_, v___x_194_);
if (v___x_196_ == 0)
{
v_a_190_ = v_tail_193_;
goto _start;
}
else
{
lean_object* v___x_198_; lean_object* v___x_199_; 
v___x_198_ = lean_unsigned_to_nat(1u);
v___x_199_ = lean_nat_add(v_a_191_, v___x_198_);
lean_dec(v_a_191_);
v_a_190_ = v_tail_193_;
v_a_191_ = v___x_199_;
goto _start;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__0___boxed(lean_object* v_a_201_, lean_object* v_a_202_){
_start:
{
lean_object* v_res_203_; 
v_res_203_ = lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__0(v_a_201_, v_a_202_);
lean_dec(v_a_201_);
return v_res_203_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__1(lean_object* v_a_204_, lean_object* v_a_205_){
_start:
{
if (lean_obj_tag(v_a_204_) == 0)
{
return v_a_205_;
}
else
{
lean_object* v_head_206_; lean_object* v_tail_207_; uint8_t v___x_208_; uint8_t v___x_209_; uint8_t v___x_210_; 
v_head_206_ = lean_ctor_get(v_a_204_, 0);
v_tail_207_ = lean_ctor_get(v_a_204_, 1);
v___x_208_ = 1;
v___x_209_ = lean_unbox(v_head_206_);
v___x_210_ = lp_BioCompiler_BioCompiler_instBEqNucleotide_beq(v___x_209_, v___x_208_);
if (v___x_210_ == 0)
{
v_a_204_ = v_tail_207_;
goto _start;
}
else
{
lean_object* v___x_212_; lean_object* v___x_213_; 
v___x_212_ = lean_unsigned_to_nat(1u);
v___x_213_ = lean_nat_add(v_a_205_, v___x_212_);
lean_dec(v_a_205_);
v_a_204_ = v_tail_207_;
v_a_205_ = v___x_213_;
goto _start;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__1___boxed(lean_object* v_a_215_, lean_object* v_a_216_){
_start:
{
lean_object* v_res_217_; 
v_res_217_ = lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__1(v_a_215_, v_a_216_);
lean_dec(v_a_215_);
return v_res_217_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_Sequence_gcContent___closed__0(void){
_start:
{
lean_object* v___x_218_; lean_object* v___x_219_; 
v___x_218_ = lean_unsigned_to_nat(0u);
v___x_219_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_218_);
return v___x_219_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_gcContent(lean_object* v_seq_220_){
_start:
{
lean_object* v___x_221_; lean_object* v___x_222_; uint8_t v___x_223_; 
v___x_221_ = l_List_lengthTR___redArg(v_seq_220_);
v___x_222_ = lean_unsigned_to_nat(0u);
v___x_223_ = lean_nat_dec_eq(v___x_221_, v___x_222_);
if (v___x_223_ == 0)
{
lean_object* v___x_224_; lean_object* v___x_225_; lean_object* v___x_226_; lean_object* v___x_227_; lean_object* v___x_228_; lean_object* v___x_229_; 
v___x_224_ = lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__0(v_seq_220_, v___x_222_);
v___x_225_ = lp_BioCompiler_List_countP_go___at___00BioCompiler_Sequence_gcContent_spec__1(v_seq_220_, v___x_222_);
v___x_226_ = lean_nat_add(v___x_224_, v___x_225_);
lean_dec(v___x_225_);
lean_dec(v___x_224_);
v___x_227_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_226_);
v___x_228_ = l_Nat_cast___at___00Dyadic_toRat_spec__0(v___x_221_);
v___x_229_ = l_Rat_div(v___x_227_, v___x_228_);
lean_dec_ref(v___x_227_);
return v___x_229_;
}
else
{
lean_object* v___x_230_; 
lean_dec(v___x_221_);
v___x_230_ = lean_obj_once(&lp_BioCompiler_BioCompiler_Sequence_gcContent___closed__0, &lp_BioCompiler_BioCompiler_Sequence_gcContent___closed__0_once, _init_lp_BioCompiler_BioCompiler_Sequence_gcContent___closed__0);
return v___x_230_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_gcContent___boxed(lean_object* v_seq_231_){
_start:
{
lean_object* v_res_232_; 
v_res_232_ = lp_BioCompiler_BioCompiler_Sequence_gcContent(v_seq_231_);
lean_dec(v_seq_231_);
return v_res_232_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_matchesAt(lean_object* v_seq_235_, lean_object* v_pattern_236_, lean_object* v_pos_237_){
_start:
{
lean_object* v___x_238_; lean_object* v___x_239_; lean_object* v___x_240_; uint8_t v___x_241_; 
v___x_238_ = l_List_lengthTR___redArg(v_pattern_236_);
v___x_239_ = lean_nat_add(v_pos_237_, v___x_238_);
v___x_240_ = l_List_lengthTR___redArg(v_seq_235_);
v___x_241_ = lean_nat_dec_le(v___x_239_, v___x_240_);
lean_dec(v___x_240_);
lean_dec(v___x_239_);
if (v___x_241_ == 0)
{
lean_dec(v___x_238_);
lean_dec(v_pos_237_);
lean_dec(v_pattern_236_);
return v___x_241_;
}
else
{
lean_object* v___x_242_; lean_object* v___x_243_; lean_object* v___x_244_; lean_object* v___x_245_; uint8_t v___x_246_; 
v___x_242_ = lean_alloc_closure((void*)(lp_BioCompiler_BioCompiler_instDecidableEqNucleotide___boxed), 2, 0);
v___x_243_ = l_List_drop___redArg(v_pos_237_, v_seq_235_);
v___x_244_ = ((lean_object*)(lp_BioCompiler_BioCompiler_Sequence_matchesAt___closed__0));
lean_inc(v___x_243_);
v___x_245_ = l___private_Init_Data_List_Impl_0__List_takeTR_go___redArg(v___x_243_, v___x_243_, v___x_238_, v___x_244_);
lean_dec(v___x_243_);
v___x_246_ = l_instDecidableEqList___redArg(v___x_242_, v___x_245_, v_pattern_236_);
return v___x_246_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_matchesAt___boxed(lean_object* v_seq_247_, lean_object* v_pattern_248_, lean_object* v_pos_249_){
_start:
{
uint8_t v_res_250_; lean_object* v_r_251_; 
v_res_250_ = lp_BioCompiler_BioCompiler_Sequence_matchesAt(v_seq_247_, v_pattern_248_, v_pos_249_);
lean_dec(v_seq_247_);
v_r_251_ = lean_box(v_res_250_);
return v_r_251_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_containsPattern___lam__0(lean_object* v_seq_252_, lean_object* v_pattern_253_, lean_object* v_pos_254_){
_start:
{
uint8_t v___x_255_; 
v___x_255_ = lp_BioCompiler_BioCompiler_Sequence_matchesAt(v_seq_252_, v_pattern_253_, v_pos_254_);
return v___x_255_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_containsPattern___lam__0___boxed(lean_object* v_seq_256_, lean_object* v_pattern_257_, lean_object* v_pos_258_){
_start:
{
uint8_t v_res_259_; lean_object* v_r_260_; 
v_res_259_ = lp_BioCompiler_BioCompiler_Sequence_containsPattern___lam__0(v_seq_256_, v_pattern_257_, v_pos_258_);
lean_dec(v_seq_256_);
v_r_260_ = lean_box(v_res_259_);
return v_r_260_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_containsPattern(lean_object* v_seq_261_, lean_object* v_pattern_262_){
_start:
{
lean_object* v___x_263_; lean_object* v___x_264_; uint8_t v___x_265_; 
v___x_263_ = l_List_lengthTR___redArg(v_pattern_262_);
v___x_264_ = lean_unsigned_to_nat(0u);
v___x_265_ = lean_nat_dec_eq(v___x_263_, v___x_264_);
if (v___x_265_ == 0)
{
lean_object* v___x_266_; uint8_t v___x_267_; 
v___x_266_ = l_List_lengthTR___redArg(v_seq_261_);
v___x_267_ = lean_nat_dec_lt(v___x_266_, v___x_263_);
if (v___x_267_ == 0)
{
lean_object* v___f_268_; lean_object* v_maxPos_269_; lean_object* v___x_270_; lean_object* v___x_271_; lean_object* v___x_272_; uint8_t v___x_273_; 
v___f_268_ = lean_alloc_closure((void*)(lp_BioCompiler_BioCompiler_Sequence_containsPattern___lam__0___boxed), 3, 2);
lean_closure_set(v___f_268_, 0, v_seq_261_);
lean_closure_set(v___f_268_, 1, v_pattern_262_);
v_maxPos_269_ = lean_nat_sub(v___x_266_, v___x_263_);
lean_dec(v___x_263_);
lean_dec(v___x_266_);
v___x_270_ = lean_unsigned_to_nat(1u);
v___x_271_ = lean_nat_add(v_maxPos_269_, v___x_270_);
lean_dec(v_maxPos_269_);
v___x_272_ = l_List_range(v___x_271_);
v___x_273_ = l_List_any___redArg(v___x_272_, v___f_268_);
return v___x_273_;
}
else
{
lean_dec(v___x_266_);
lean_dec(v___x_263_);
lean_dec(v_pattern_262_);
lean_dec(v_seq_261_);
return v___x_265_;
}
}
else
{
lean_dec(v___x_263_);
lean_dec(v_pattern_262_);
lean_dec(v_seq_261_);
return v___x_265_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_containsPattern___boxed(lean_object* v_seq_274_, lean_object* v_pattern_275_){
_start:
{
uint8_t v_res_276_; lean_object* v_r_277_; 
v_res_276_ = lp_BioCompiler_BioCompiler_Sequence_containsPattern(v_seq_274_, v_pattern_275_);
v_r_277_ = lean_box(v_res_276_);
return v_r_277_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_isStopCodon___lam__0(lean_object* v_codon_320_, lean_object* v_c_321_){
_start:
{
lean_object* v___x_322_; uint8_t v___x_323_; 
v___x_322_ = lean_alloc_closure((void*)(lp_BioCompiler_BioCompiler_instDecidableEqNucleotide___boxed), 2, 0);
v___x_323_ = l_instDecidableEqList___redArg(v___x_322_, v_c_321_, v_codon_320_);
return v___x_323_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_isStopCodon___lam__0___boxed(lean_object* v_codon_324_, lean_object* v_c_325_){
_start:
{
uint8_t v_res_326_; lean_object* v_r_327_; 
v_res_326_ = lp_BioCompiler_BioCompiler_Sequence_isStopCodon___lam__0(v_codon_324_, v_c_325_);
v_r_327_ = lean_box(v_res_326_);
return v_r_327_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_isStopCodon(lean_object* v_codon_328_){
_start:
{
lean_object* v___f_329_; lean_object* v___x_330_; uint8_t v___x_331_; 
v___f_329_ = lean_alloc_closure((void*)(lp_BioCompiler_BioCompiler_Sequence_isStopCodon___lam__0___boxed), 2, 1);
lean_closure_set(v___f_329_, 0, v_codon_328_);
v___x_330_ = ((lean_object*)(lp_BioCompiler_BioCompiler_Sequence_stopCodons));
v___x_331_ = l_List_any___redArg(v___x_330_, v___f_329_);
return v___x_331_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_isStopCodon___boxed(lean_object* v_codon_332_){
_start:
{
uint8_t v_res_333_; lean_object* v_r_334_; 
v_res_333_ = lp_BioCompiler_BioCompiler_Sequence_isStopCodon(v_codon_332_);
v_r_334_ = lean_box(v_res_333_);
return v_r_334_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop___lam__0(lean_object* v___x_335_, lean_object* v___x_336_, lean_object* v_seq_337_, lean_object* v_start_338_){
_start:
{
lean_object* v___x_339_; uint8_t v___x_340_; 
v___x_339_ = lean_nat_add(v_start_338_, v___x_335_);
v___x_340_ = lean_nat_dec_le(v___x_339_, v___x_336_);
lean_dec(v___x_339_);
if (v___x_340_ == 0)
{
lean_dec(v_start_338_);
lean_dec(v___x_335_);
return v___x_340_;
}
else
{
lean_object* v___x_341_; lean_object* v___x_342_; lean_object* v___x_343_; uint8_t v___x_344_; 
v___x_341_ = l_List_drop___redArg(v_start_338_, v_seq_337_);
v___x_342_ = ((lean_object*)(lp_BioCompiler_BioCompiler_Sequence_matchesAt___closed__0));
lean_inc(v___x_341_);
v___x_343_ = l___private_Init_Data_List_Impl_0__List_takeTR_go___redArg(v___x_341_, v___x_341_, v___x_335_, v___x_342_);
lean_dec(v___x_341_);
v___x_344_ = lp_BioCompiler_BioCompiler_Sequence_isStopCodon(v___x_343_);
return v___x_344_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop___lam__0___boxed(lean_object* v___x_345_, lean_object* v___x_346_, lean_object* v_seq_347_, lean_object* v_start_348_){
_start:
{
uint8_t v_res_349_; lean_object* v_r_350_; 
v_res_349_ = lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop___lam__0(v___x_345_, v___x_346_, v_seq_347_, v_start_348_);
lean_dec(v_seq_347_);
lean_dec(v___x_346_);
v_r_350_ = lean_box(v_res_349_);
return v_r_350_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_Sequence_hasPrematureStop_spec__0(lean_object* v_readingFrame_351_, lean_object* v_a_352_, lean_object* v_a_353_){
_start:
{
if (lean_obj_tag(v_a_352_) == 0)
{
lean_object* v___x_354_; 
v___x_354_ = l_List_reverse___redArg(v_a_353_);
return v___x_354_;
}
else
{
lean_object* v_head_355_; lean_object* v_tail_356_; lean_object* v___x_358_; uint8_t v_isShared_359_; uint8_t v_isSharedCheck_367_; 
v_head_355_ = lean_ctor_get(v_a_352_, 0);
v_tail_356_ = lean_ctor_get(v_a_352_, 1);
v_isSharedCheck_367_ = !lean_is_exclusive(v_a_352_);
if (v_isSharedCheck_367_ == 0)
{
v___x_358_ = v_a_352_;
v_isShared_359_ = v_isSharedCheck_367_;
goto v_resetjp_357_;
}
else
{
lean_inc(v_tail_356_);
lean_inc(v_head_355_);
lean_dec(v_a_352_);
v___x_358_ = lean_box(0);
v_isShared_359_ = v_isSharedCheck_367_;
goto v_resetjp_357_;
}
v_resetjp_357_:
{
lean_object* v___x_360_; lean_object* v___x_361_; lean_object* v___x_362_; lean_object* v___x_364_; 
v___x_360_ = lean_unsigned_to_nat(3u);
v___x_361_ = lean_nat_mul(v___x_360_, v_head_355_);
lean_dec(v_head_355_);
v___x_362_ = lean_nat_add(v_readingFrame_351_, v___x_361_);
lean_dec(v___x_361_);
if (v_isShared_359_ == 0)
{
lean_ctor_set(v___x_358_, 1, v_a_353_);
lean_ctor_set(v___x_358_, 0, v___x_362_);
v___x_364_ = v___x_358_;
goto v_reusejp_363_;
}
else
{
lean_object* v_reuseFailAlloc_366_; 
v_reuseFailAlloc_366_ = lean_alloc_ctor(1, 2, 0);
lean_ctor_set(v_reuseFailAlloc_366_, 0, v___x_362_);
lean_ctor_set(v_reuseFailAlloc_366_, 1, v_a_353_);
v___x_364_ = v_reuseFailAlloc_366_;
goto v_reusejp_363_;
}
v_reusejp_363_:
{
v_a_352_ = v_tail_356_;
v_a_353_ = v___x_364_;
goto _start;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_Sequence_hasPrematureStop_spec__0___boxed(lean_object* v_readingFrame_368_, lean_object* v_a_369_, lean_object* v_a_370_){
_start:
{
lean_object* v_res_371_; 
v_res_371_ = lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_Sequence_hasPrematureStop_spec__0(v_readingFrame_368_, v_a_369_, v_a_370_);
lean_dec(v_readingFrame_368_);
return v_res_371_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(lean_object* v_seq_372_, lean_object* v_readingFrame_373_){
_start:
{
lean_object* v___x_374_; lean_object* v___x_375_; lean_object* v___x_376_; lean_object* v___f_377_; lean_object* v___x_378_; lean_object* v___x_379_; lean_object* v___x_380_; lean_object* v_codonStarts_381_; uint8_t v___x_382_; 
v___x_374_ = l_List_lengthTR___redArg(v_seq_372_);
v___x_375_ = lean_nat_sub(v___x_374_, v_readingFrame_373_);
v___x_376_ = lean_unsigned_to_nat(3u);
v___f_377_ = lean_alloc_closure((void*)(lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop___lam__0___boxed), 4, 3);
lean_closure_set(v___f_377_, 0, v___x_376_);
lean_closure_set(v___f_377_, 1, v___x_374_);
lean_closure_set(v___f_377_, 2, v_seq_372_);
v___x_378_ = lean_nat_div(v___x_375_, v___x_376_);
lean_dec(v___x_375_);
v___x_379_ = l_List_range(v___x_378_);
v___x_380_ = lean_box(0);
v_codonStarts_381_ = lp_BioCompiler_List_mapTR_loop___at___00BioCompiler_Sequence_hasPrematureStop_spec__0(v_readingFrame_373_, v___x_379_, v___x_380_);
v___x_382_ = l_List_any___redArg(v_codonStarts_381_, v___f_377_);
return v___x_382_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop___boxed(lean_object* v_seq_383_, lean_object* v_readingFrame_384_){
_start:
{
uint8_t v_res_385_; lean_object* v_r_386_; 
v_res_385_ = lp_BioCompiler_BioCompiler_Sequence_hasPrematureStop(v_seq_383_, v_readingFrame_384_);
lean_dec(v_readingFrame_384_);
v_r_386_ = lean_box(v_res_385_);
return v_r_386_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent___lam__0(lean_object* v_readingFrame_387_, lean_object* v_pos_388_){
_start:
{
lean_object* v___x_389_; lean_object* v___x_390_; lean_object* v___x_391_; uint8_t v___x_392_; 
v___x_389_ = lean_unsigned_to_nat(3u);
v___x_390_ = lean_nat_mod(v_pos_388_, v___x_389_);
v___x_391_ = lean_nat_mod(v_readingFrame_387_, v___x_389_);
v___x_392_ = lean_nat_dec_eq(v___x_390_, v___x_391_);
lean_dec(v___x_391_);
lean_dec(v___x_390_);
return v___x_392_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent___lam__0___boxed(lean_object* v_readingFrame_393_, lean_object* v_pos_394_){
_start:
{
uint8_t v_res_395_; lean_object* v_r_396_; 
v_res_395_ = lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent___lam__0(v_readingFrame_393_, v_pos_394_);
lean_dec(v_pos_394_);
lean_dec(v_readingFrame_393_);
v_r_396_ = lean_box(v_res_395_);
return v_r_396_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent(lean_object* v_boundaries_397_, lean_object* v_readingFrame_398_){
_start:
{
lean_object* v___f_399_; uint8_t v___x_400_; 
v___f_399_ = lean_alloc_closure((void*)(lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent___lam__0___boxed), 2, 1);
lean_closure_set(v___f_399_, 0, v_readingFrame_398_);
v___x_400_ = l_List_all___redArg(v_boundaries_397_, v___f_399_);
return v___x_400_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent___boxed(lean_object* v_boundaries_401_, lean_object* v_readingFrame_402_){
_start:
{
uint8_t v_res_403_; lean_object* v_r_404_; 
v_res_403_ = lp_BioCompiler_BioCompiler_Sequence_readingFrameConsistent(v_boundaries_401_, v_readingFrame_402_);
v_r_404_ = lean_box(v_res_403_);
return v_r_404_;
}
}
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_Init(uint8_t builtin);
static bool _G_initialized = false;
LEAN_EXPORT lean_object* initialize_BioCompiler_BioCompiler_Sequence(uint8_t builtin) {
lean_object * res;
if (_G_initialized) return lean_io_result_mk_ok(lean_box(0));
_G_initialized = true;
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
lp_BioCompiler_BioCompiler_instInhabitedNucleotide_default = _init_lp_BioCompiler_BioCompiler_instInhabitedNucleotide_default();
lp_BioCompiler_BioCompiler_instInhabitedNucleotide = _init_lp_BioCompiler_BioCompiler_instInhabitedNucleotide();
lp_BioCompiler_BioCompiler_Sequence_empty = _init_lp_BioCompiler_BioCompiler_Sequence_empty();
lean_mark_persistent(lp_BioCompiler_BioCompiler_Sequence_empty);
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
