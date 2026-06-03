// Lean compiler output
// Module: BioCompiler.ThreeValued
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
uint8_t lean_nat_dec_le(lean_object*, lean_object*);
uint8_t lean_nat_dec_eq(lean_object*, lean_object*);
lean_object* lean_nat_to_int(lean_object*);
lean_object* l_Repr_addAppParen(lean_object*, lean_object*);
uint8_t lean_nat_dec_le(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorIdx(uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorIdx___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_toCtorIdx(uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_toCtorIdx___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorElim___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorElim___redArg___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorElim(lean_object*, lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorElim___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_PASS_elim___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_PASS_elim___redArg___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_PASS_elim(lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_PASS_elim___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_FAIL_elim___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_FAIL_elim___redArg___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_FAIL_elim(lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_FAIL_elim___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim___redArg(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim___redArg___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim(lean_object*, uint8_t, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Verdict_ofNat(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ofNat___boxed(lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instDecidableEqVerdict(uint8_t, uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instDecidableEqVerdict___boxed(lean_object*, lean_object*);
static const lean_string_object lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 25, .m_capacity = 25, .m_length = 24, .m_data = "BioCompiler.Verdict.PASS"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__0_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__1_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__0_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__1 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__1_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__2_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 25, .m_capacity = 25, .m_length = 24, .m_data = "BioCompiler.Verdict.FAIL"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__2 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__2_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__3_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__2_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__3 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__3_value;
static const lean_string_object lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__4_value = {.m_header = {.m_rc = 0, .m_cs_sz = 0, .m_other = 0, .m_tag = 249}, .m_size = 30, .m_capacity = 30, .m_length = 29, .m_data = "BioCompiler.Verdict.UNCERTAIN"};
static const lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__4 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__4_value;
static const lean_ctor_object lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__5_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_ctor_object) + sizeof(void*)*1 + 0, .m_other = 1, .m_tag = 3}, .m_objs = {((lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__4_value)}};
static const lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__5 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__5_value;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6;
static lean_once_cell_t lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7_once = LEAN_ONCE_CELL_INITIALIZER;
static lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7;
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr(uint8_t, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instReprVerdict___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instReprVerdict_repr___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instReprVerdict___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instReprVerdict = (const lean_object*)&lp_BioCompiler_BioCompiler_instReprVerdict___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instBEqVerdict_beq(uint8_t, uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instBEqVerdict_beq___boxed(lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_instBEqVerdict___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_instBEqVerdict_beq___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_instBEqVerdict___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_instBEqVerdict___closed__0_value;
LEAN_EXPORT const lean_object* lp_BioCompiler_BioCompiler_instBEqVerdict = (const lean_object*)&lp_BioCompiler_BioCompiler_instBEqVerdict___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Verdict_and(uint8_t, uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_and___boxed(lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Verdict_or(uint8_t, uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_or___boxed(lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Verdict_not(uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_not___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter___redArg(uint8_t, uint8_t, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter___redArg___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter(lean_object*, uint8_t, uint8_t, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ordering(uint8_t);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ordering___boxed(lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter___redArg(uint8_t, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter___redArg___boxed(lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter(lean_object*, uint8_t, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorIdx(uint8_t v_x_1_){
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
default: 
{
lean_object* v___x_4_; 
v___x_4_ = lean_unsigned_to_nat(2u);
return v___x_4_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorIdx___boxed(lean_object* v_x_5_){
_start:
{
uint8_t v_x_boxed_6_; lean_object* v_res_7_; 
v_x_boxed_6_ = lean_unbox(v_x_5_);
v_res_7_ = lp_BioCompiler_BioCompiler_Verdict_ctorIdx(v_x_boxed_6_);
return v_res_7_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_toCtorIdx(uint8_t v_x_8_){
_start:
{
lean_object* v___x_9_; 
v___x_9_ = lp_BioCompiler_BioCompiler_Verdict_ctorIdx(v_x_8_);
return v___x_9_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_toCtorIdx___boxed(lean_object* v_x_10_){
_start:
{
uint8_t v_x_4__boxed_11_; lean_object* v_res_12_; 
v_x_4__boxed_11_ = lean_unbox(v_x_10_);
v_res_12_ = lp_BioCompiler_BioCompiler_Verdict_toCtorIdx(v_x_4__boxed_11_);
return v_res_12_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorElim___redArg(lean_object* v_k_13_){
_start:
{
lean_inc(v_k_13_);
return v_k_13_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorElim___redArg___boxed(lean_object* v_k_14_){
_start:
{
lean_object* v_res_15_; 
v_res_15_ = lp_BioCompiler_BioCompiler_Verdict_ctorElim___redArg(v_k_14_);
lean_dec(v_k_14_);
return v_res_15_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorElim(lean_object* v_motive_16_, lean_object* v_ctorIdx_17_, uint8_t v_t_18_, lean_object* v_h_19_, lean_object* v_k_20_){
_start:
{
lean_inc(v_k_20_);
return v_k_20_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ctorElim___boxed(lean_object* v_motive_21_, lean_object* v_ctorIdx_22_, lean_object* v_t_23_, lean_object* v_h_24_, lean_object* v_k_25_){
_start:
{
uint8_t v_t_boxed_26_; lean_object* v_res_27_; 
v_t_boxed_26_ = lean_unbox(v_t_23_);
v_res_27_ = lp_BioCompiler_BioCompiler_Verdict_ctorElim(v_motive_21_, v_ctorIdx_22_, v_t_boxed_26_, v_h_24_, v_k_25_);
lean_dec(v_k_25_);
lean_dec(v_ctorIdx_22_);
return v_res_27_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_PASS_elim___redArg(lean_object* v_PASS_28_){
_start:
{
lean_inc(v_PASS_28_);
return v_PASS_28_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_PASS_elim___redArg___boxed(lean_object* v_PASS_29_){
_start:
{
lean_object* v_res_30_; 
v_res_30_ = lp_BioCompiler_BioCompiler_Verdict_PASS_elim___redArg(v_PASS_29_);
lean_dec(v_PASS_29_);
return v_res_30_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_PASS_elim(lean_object* v_motive_31_, uint8_t v_t_32_, lean_object* v_h_33_, lean_object* v_PASS_34_){
_start:
{
lean_inc(v_PASS_34_);
return v_PASS_34_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_PASS_elim___boxed(lean_object* v_motive_35_, lean_object* v_t_36_, lean_object* v_h_37_, lean_object* v_PASS_38_){
_start:
{
uint8_t v_t_boxed_39_; lean_object* v_res_40_; 
v_t_boxed_39_ = lean_unbox(v_t_36_);
v_res_40_ = lp_BioCompiler_BioCompiler_Verdict_PASS_elim(v_motive_35_, v_t_boxed_39_, v_h_37_, v_PASS_38_);
lean_dec(v_PASS_38_);
return v_res_40_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_FAIL_elim___redArg(lean_object* v_FAIL_41_){
_start:
{
lean_inc(v_FAIL_41_);
return v_FAIL_41_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_FAIL_elim___redArg___boxed(lean_object* v_FAIL_42_){
_start:
{
lean_object* v_res_43_; 
v_res_43_ = lp_BioCompiler_BioCompiler_Verdict_FAIL_elim___redArg(v_FAIL_42_);
lean_dec(v_FAIL_42_);
return v_res_43_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_FAIL_elim(lean_object* v_motive_44_, uint8_t v_t_45_, lean_object* v_h_46_, lean_object* v_FAIL_47_){
_start:
{
lean_inc(v_FAIL_47_);
return v_FAIL_47_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_FAIL_elim___boxed(lean_object* v_motive_48_, lean_object* v_t_49_, lean_object* v_h_50_, lean_object* v_FAIL_51_){
_start:
{
uint8_t v_t_boxed_52_; lean_object* v_res_53_; 
v_t_boxed_52_ = lean_unbox(v_t_49_);
v_res_53_ = lp_BioCompiler_BioCompiler_Verdict_FAIL_elim(v_motive_48_, v_t_boxed_52_, v_h_50_, v_FAIL_51_);
lean_dec(v_FAIL_51_);
return v_res_53_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim___redArg(lean_object* v_UNCERTAIN_54_){
_start:
{
lean_inc(v_UNCERTAIN_54_);
return v_UNCERTAIN_54_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim___redArg___boxed(lean_object* v_UNCERTAIN_55_){
_start:
{
lean_object* v_res_56_; 
v_res_56_ = lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim___redArg(v_UNCERTAIN_55_);
lean_dec(v_UNCERTAIN_55_);
return v_res_56_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim(lean_object* v_motive_57_, uint8_t v_t_58_, lean_object* v_h_59_, lean_object* v_UNCERTAIN_60_){
_start:
{
lean_inc(v_UNCERTAIN_60_);
return v_UNCERTAIN_60_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim___boxed(lean_object* v_motive_61_, lean_object* v_t_62_, lean_object* v_h_63_, lean_object* v_UNCERTAIN_64_){
_start:
{
uint8_t v_t_boxed_65_; lean_object* v_res_66_; 
v_t_boxed_65_ = lean_unbox(v_t_62_);
v_res_66_ = lp_BioCompiler_BioCompiler_Verdict_UNCERTAIN_elim(v_motive_61_, v_t_boxed_65_, v_h_63_, v_UNCERTAIN_64_);
lean_dec(v_UNCERTAIN_64_);
return v_res_66_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Verdict_ofNat(lean_object* v_n_67_){
_start:
{
lean_object* v___x_68_; uint8_t v___x_69_; 
v___x_68_ = lean_unsigned_to_nat(0u);
v___x_69_ = lean_nat_dec_le(v_n_67_, v___x_68_);
if (v___x_69_ == 0)
{
lean_object* v___x_70_; uint8_t v___x_71_; 
v___x_70_ = lean_unsigned_to_nat(1u);
v___x_71_ = lean_nat_dec_le(v_n_67_, v___x_70_);
if (v___x_71_ == 0)
{
uint8_t v___x_72_; 
v___x_72_ = 2;
return v___x_72_;
}
else
{
uint8_t v___x_73_; 
v___x_73_ = 1;
return v___x_73_;
}
}
else
{
uint8_t v___x_74_; 
v___x_74_ = 0;
return v___x_74_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ofNat___boxed(lean_object* v_n_75_){
_start:
{
uint8_t v_res_76_; lean_object* v_r_77_; 
v_res_76_ = lp_BioCompiler_BioCompiler_Verdict_ofNat(v_n_75_);
lean_dec(v_n_75_);
v_r_77_ = lean_box(v_res_76_);
return v_r_77_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instDecidableEqVerdict(uint8_t v_x_78_, uint8_t v_y_79_){
_start:
{
lean_object* v___x_80_; lean_object* v___x_81_; uint8_t v___x_82_; 
v___x_80_ = lp_BioCompiler_BioCompiler_Verdict_ctorIdx(v_x_78_);
v___x_81_ = lp_BioCompiler_BioCompiler_Verdict_ctorIdx(v_y_79_);
v___x_82_ = lean_nat_dec_eq(v___x_80_, v___x_81_);
lean_dec(v___x_81_);
lean_dec(v___x_80_);
return v___x_82_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instDecidableEqVerdict___boxed(lean_object* v_x_83_, lean_object* v_y_84_){
_start:
{
uint8_t v_x_13__boxed_85_; uint8_t v_y_14__boxed_86_; uint8_t v_res_87_; lean_object* v_r_88_; 
v_x_13__boxed_85_ = lean_unbox(v_x_83_);
v_y_14__boxed_86_ = lean_unbox(v_y_84_);
v_res_87_ = lp_BioCompiler_BioCompiler_instDecidableEqVerdict(v_x_13__boxed_85_, v_y_14__boxed_86_);
v_r_88_ = lean_box(v_res_87_);
return v_r_88_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6(void){
_start:
{
lean_object* v___x_98_; lean_object* v___x_99_; 
v___x_98_ = lean_unsigned_to_nat(2u);
v___x_99_ = lean_nat_to_int(v___x_98_);
return v___x_99_;
}
}
static lean_object* _init_lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7(void){
_start:
{
lean_object* v___x_100_; lean_object* v___x_101_; 
v___x_100_ = lean_unsigned_to_nat(1u);
v___x_101_ = lean_nat_to_int(v___x_100_);
return v___x_101_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr(uint8_t v_x_102_, lean_object* v_prec_103_){
_start:
{
lean_object* v___y_105_; lean_object* v___y_112_; lean_object* v___y_119_; 
switch(v_x_102_)
{
case 0:
{
lean_object* v___x_125_; uint8_t v___x_126_; 
v___x_125_ = lean_unsigned_to_nat(1024u);
v___x_126_ = lean_nat_dec_le(v___x_125_, v_prec_103_);
if (v___x_126_ == 0)
{
lean_object* v___x_127_; 
v___x_127_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6, &lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6);
v___y_105_ = v___x_127_;
goto v___jp_104_;
}
else
{
lean_object* v___x_128_; 
v___x_128_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7, &lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7);
v___y_105_ = v___x_128_;
goto v___jp_104_;
}
}
case 1:
{
lean_object* v___x_129_; uint8_t v___x_130_; 
v___x_129_ = lean_unsigned_to_nat(1024u);
v___x_130_ = lean_nat_dec_le(v___x_129_, v_prec_103_);
if (v___x_130_ == 0)
{
lean_object* v___x_131_; 
v___x_131_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6, &lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6);
v___y_112_ = v___x_131_;
goto v___jp_111_;
}
else
{
lean_object* v___x_132_; 
v___x_132_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7, &lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7);
v___y_112_ = v___x_132_;
goto v___jp_111_;
}
}
default: 
{
lean_object* v___x_133_; uint8_t v___x_134_; 
v___x_133_ = lean_unsigned_to_nat(1024u);
v___x_134_ = lean_nat_dec_le(v___x_133_, v_prec_103_);
if (v___x_134_ == 0)
{
lean_object* v___x_135_; 
v___x_135_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6, &lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6_once, _init_lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__6);
v___y_119_ = v___x_135_;
goto v___jp_118_;
}
else
{
lean_object* v___x_136_; 
v___x_136_ = lean_obj_once(&lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7, &lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7_once, _init_lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__7);
v___y_119_ = v___x_136_;
goto v___jp_118_;
}
}
}
v___jp_104_:
{
lean_object* v___x_106_; lean_object* v___x_107_; uint8_t v___x_108_; lean_object* v___x_109_; lean_object* v___x_110_; 
v___x_106_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__1));
lean_inc(v___y_105_);
v___x_107_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_107_, 0, v___y_105_);
lean_ctor_set(v___x_107_, 1, v___x_106_);
v___x_108_ = 0;
v___x_109_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_109_, 0, v___x_107_);
lean_ctor_set_uint8(v___x_109_, sizeof(void*)*1, v___x_108_);
v___x_110_ = l_Repr_addAppParen(v___x_109_, v_prec_103_);
return v___x_110_;
}
v___jp_111_:
{
lean_object* v___x_113_; lean_object* v___x_114_; uint8_t v___x_115_; lean_object* v___x_116_; lean_object* v___x_117_; 
v___x_113_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__3));
lean_inc(v___y_112_);
v___x_114_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_114_, 0, v___y_112_);
lean_ctor_set(v___x_114_, 1, v___x_113_);
v___x_115_ = 0;
v___x_116_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_116_, 0, v___x_114_);
lean_ctor_set_uint8(v___x_116_, sizeof(void*)*1, v___x_115_);
v___x_117_ = l_Repr_addAppParen(v___x_116_, v_prec_103_);
return v___x_117_;
}
v___jp_118_:
{
lean_object* v___x_120_; lean_object* v___x_121_; uint8_t v___x_122_; lean_object* v___x_123_; lean_object* v___x_124_; 
v___x_120_ = ((lean_object*)(lp_BioCompiler_BioCompiler_instReprVerdict_repr___closed__5));
lean_inc(v___y_119_);
v___x_121_ = lean_alloc_ctor(4, 2, 0);
lean_ctor_set(v___x_121_, 0, v___y_119_);
lean_ctor_set(v___x_121_, 1, v___x_120_);
v___x_122_ = 0;
v___x_123_ = lean_alloc_ctor(6, 1, 1);
lean_ctor_set(v___x_123_, 0, v___x_121_);
lean_ctor_set_uint8(v___x_123_, sizeof(void*)*1, v___x_122_);
v___x_124_ = l_Repr_addAppParen(v___x_123_, v_prec_103_);
return v___x_124_;
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instReprVerdict_repr___boxed(lean_object* v_x_137_, lean_object* v_prec_138_){
_start:
{
uint8_t v_x_177__boxed_139_; lean_object* v_res_140_; 
v_x_177__boxed_139_ = lean_unbox(v_x_137_);
v_res_140_ = lp_BioCompiler_BioCompiler_instReprVerdict_repr(v_x_177__boxed_139_, v_prec_138_);
lean_dec(v_prec_138_);
return v_res_140_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_instBEqVerdict_beq(uint8_t v_x_143_, uint8_t v_y_144_){
_start:
{
lean_object* v___x_145_; lean_object* v___x_146_; uint8_t v___x_147_; 
v___x_145_ = lp_BioCompiler_BioCompiler_Verdict_ctorIdx(v_x_143_);
v___x_146_ = lp_BioCompiler_BioCompiler_Verdict_ctorIdx(v_y_144_);
v___x_147_ = lean_nat_dec_eq(v___x_145_, v___x_146_);
lean_dec(v___x_146_);
lean_dec(v___x_145_);
return v___x_147_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_instBEqVerdict_beq___boxed(lean_object* v_x_148_, lean_object* v_y_149_){
_start:
{
uint8_t v_x_17__boxed_150_; uint8_t v_y_18__boxed_151_; uint8_t v_res_152_; lean_object* v_r_153_; 
v_x_17__boxed_150_ = lean_unbox(v_x_148_);
v_y_18__boxed_151_ = lean_unbox(v_y_149_);
v_res_152_ = lp_BioCompiler_BioCompiler_instBEqVerdict_beq(v_x_17__boxed_150_, v_y_18__boxed_151_);
v_r_153_ = lean_box(v_res_152_);
return v_r_153_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Verdict_and(uint8_t v_x_156_, uint8_t v_x_157_){
_start:
{
switch(v_x_156_)
{
case 0:
{
return v_x_157_;
}
case 1:
{
if (v_x_157_ == 1)
{
return v_x_157_;
}
else
{
return v_x_156_;
}
}
default: 
{
if (v_x_157_ == 0)
{
return v_x_156_;
}
else
{
return v_x_157_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_and___boxed(lean_object* v_x_158_, lean_object* v_x_159_){
_start:
{
uint8_t v_x_59__boxed_160_; uint8_t v_x_60__boxed_161_; uint8_t v_res_162_; lean_object* v_r_163_; 
v_x_59__boxed_160_ = lean_unbox(v_x_158_);
v_x_60__boxed_161_ = lean_unbox(v_x_159_);
v_res_162_ = lp_BioCompiler_BioCompiler_Verdict_and(v_x_59__boxed_160_, v_x_60__boxed_161_);
v_r_163_ = lean_box(v_res_162_);
return v_r_163_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Verdict_or(uint8_t v_x_164_, uint8_t v_x_165_){
_start:
{
switch(v_x_164_)
{
case 0:
{
if (v_x_165_ == 0)
{
return v_x_165_;
}
else
{
return v_x_164_;
}
}
case 1:
{
return v_x_165_;
}
default: 
{
if (v_x_165_ == 1)
{
return v_x_164_;
}
else
{
return v_x_165_;
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_or___boxed(lean_object* v_x_166_, lean_object* v_x_167_){
_start:
{
uint8_t v_x_59__boxed_168_; uint8_t v_x_60__boxed_169_; uint8_t v_res_170_; lean_object* v_r_171_; 
v_x_59__boxed_168_ = lean_unbox(v_x_166_);
v_x_60__boxed_169_ = lean_unbox(v_x_167_);
v_res_170_ = lp_BioCompiler_BioCompiler_Verdict_or(v_x_59__boxed_168_, v_x_60__boxed_169_);
v_r_171_ = lean_box(v_res_170_);
return v_r_171_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_Verdict_not(uint8_t v_x_172_){
_start:
{
switch(v_x_172_)
{
case 0:
{
uint8_t v___x_173_; 
v___x_173_ = 1;
return v___x_173_;
}
case 1:
{
uint8_t v___x_174_; 
v___x_174_ = 0;
return v___x_174_;
}
default: 
{
return v_x_172_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_not___boxed(lean_object* v_x_175_){
_start:
{
uint8_t v_x_25__boxed_176_; uint8_t v_res_177_; lean_object* v_r_178_; 
v_x_25__boxed_176_ = lean_unbox(v_x_175_);
v_res_177_ = lp_BioCompiler_BioCompiler_Verdict_not(v_x_25__boxed_176_);
v_r_178_ = lean_box(v_res_177_);
return v_r_178_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter___redArg(uint8_t v_x_179_, uint8_t v_x_180_, lean_object* v_h__1_181_, lean_object* v_h__2_182_, lean_object* v_h__3_183_, lean_object* v_h__4_184_, lean_object* v_h__5_185_, lean_object* v_h__6_186_, lean_object* v_h__7_187_, lean_object* v_h__8_188_, lean_object* v_h__9_189_){
_start:
{
switch(v_x_179_)
{
case 0:
{
lean_dec(v_h__9_189_);
lean_dec(v_h__8_188_);
lean_dec(v_h__7_187_);
lean_dec(v_h__6_186_);
lean_dec(v_h__5_185_);
lean_dec(v_h__4_184_);
switch(v_x_180_)
{
case 0:
{
lean_object* v___x_190_; lean_object* v___x_191_; 
lean_dec(v_h__3_183_);
lean_dec(v_h__2_182_);
v___x_190_ = lean_box(0);
v___x_191_ = lean_apply_1(v_h__1_181_, v___x_190_);
return v___x_191_;
}
case 1:
{
lean_object* v___x_192_; lean_object* v___x_193_; 
lean_dec(v_h__2_182_);
lean_dec(v_h__1_181_);
v___x_192_ = lean_box(0);
v___x_193_ = lean_apply_1(v_h__3_183_, v___x_192_);
return v___x_193_;
}
default: 
{
lean_object* v___x_194_; lean_object* v___x_195_; 
lean_dec(v_h__3_183_);
lean_dec(v_h__1_181_);
v___x_194_ = lean_box(0);
v___x_195_ = lean_apply_1(v_h__2_182_, v___x_194_);
return v___x_195_;
}
}
}
case 1:
{
lean_dec(v_h__6_186_);
lean_dec(v_h__5_185_);
lean_dec(v_h__4_184_);
lean_dec(v_h__3_183_);
lean_dec(v_h__2_182_);
lean_dec(v_h__1_181_);
switch(v_x_180_)
{
case 0:
{
lean_object* v___x_196_; lean_object* v___x_197_; 
lean_dec(v_h__9_189_);
lean_dec(v_h__8_188_);
v___x_196_ = lean_box(0);
v___x_197_ = lean_apply_1(v_h__7_187_, v___x_196_);
return v___x_197_;
}
case 1:
{
lean_object* v___x_198_; lean_object* v___x_199_; 
lean_dec(v_h__8_188_);
lean_dec(v_h__7_187_);
v___x_198_ = lean_box(0);
v___x_199_ = lean_apply_1(v_h__9_189_, v___x_198_);
return v___x_199_;
}
default: 
{
lean_object* v___x_200_; lean_object* v___x_201_; 
lean_dec(v_h__9_189_);
lean_dec(v_h__7_187_);
v___x_200_ = lean_box(0);
v___x_201_ = lean_apply_1(v_h__8_188_, v___x_200_);
return v___x_201_;
}
}
}
default: 
{
lean_dec(v_h__9_189_);
lean_dec(v_h__8_188_);
lean_dec(v_h__7_187_);
lean_dec(v_h__3_183_);
lean_dec(v_h__2_182_);
lean_dec(v_h__1_181_);
switch(v_x_180_)
{
case 0:
{
lean_object* v___x_202_; lean_object* v___x_203_; 
lean_dec(v_h__6_186_);
lean_dec(v_h__5_185_);
v___x_202_ = lean_box(0);
v___x_203_ = lean_apply_1(v_h__4_184_, v___x_202_);
return v___x_203_;
}
case 1:
{
lean_object* v___x_204_; lean_object* v___x_205_; 
lean_dec(v_h__5_185_);
lean_dec(v_h__4_184_);
v___x_204_ = lean_box(0);
v___x_205_ = lean_apply_1(v_h__6_186_, v___x_204_);
return v___x_205_;
}
default: 
{
lean_object* v___x_206_; lean_object* v___x_207_; 
lean_dec(v_h__6_186_);
lean_dec(v_h__4_184_);
v___x_206_ = lean_box(0);
v___x_207_ = lean_apply_1(v_h__5_185_, v___x_206_);
return v___x_207_;
}
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter___redArg___boxed(lean_object* v_x_208_, lean_object* v_x_209_, lean_object* v_h__1_210_, lean_object* v_h__2_211_, lean_object* v_h__3_212_, lean_object* v_h__4_213_, lean_object* v_h__5_214_, lean_object* v_h__6_215_, lean_object* v_h__7_216_, lean_object* v_h__8_217_, lean_object* v_h__9_218_){
_start:
{
uint8_t v_x_101__boxed_219_; uint8_t v_x_102__boxed_220_; lean_object* v_res_221_; 
v_x_101__boxed_219_ = lean_unbox(v_x_208_);
v_x_102__boxed_220_ = lean_unbox(v_x_209_);
v_res_221_ = lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter___redArg(v_x_101__boxed_219_, v_x_102__boxed_220_, v_h__1_210_, v_h__2_211_, v_h__3_212_, v_h__4_213_, v_h__5_214_, v_h__6_215_, v_h__7_216_, v_h__8_217_, v_h__9_218_);
return v_res_221_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter(lean_object* v_motive_222_, uint8_t v_x_223_, uint8_t v_x_224_, lean_object* v_h__1_225_, lean_object* v_h__2_226_, lean_object* v_h__3_227_, lean_object* v_h__4_228_, lean_object* v_h__5_229_, lean_object* v_h__6_230_, lean_object* v_h__7_231_, lean_object* v_h__8_232_, lean_object* v_h__9_233_){
_start:
{
switch(v_x_223_)
{
case 0:
{
lean_dec(v_h__9_233_);
lean_dec(v_h__8_232_);
lean_dec(v_h__7_231_);
lean_dec(v_h__6_230_);
lean_dec(v_h__5_229_);
lean_dec(v_h__4_228_);
switch(v_x_224_)
{
case 0:
{
lean_object* v___x_234_; lean_object* v___x_235_; 
lean_dec(v_h__3_227_);
lean_dec(v_h__2_226_);
v___x_234_ = lean_box(0);
v___x_235_ = lean_apply_1(v_h__1_225_, v___x_234_);
return v___x_235_;
}
case 1:
{
lean_object* v___x_236_; lean_object* v___x_237_; 
lean_dec(v_h__2_226_);
lean_dec(v_h__1_225_);
v___x_236_ = lean_box(0);
v___x_237_ = lean_apply_1(v_h__3_227_, v___x_236_);
return v___x_237_;
}
default: 
{
lean_object* v___x_238_; lean_object* v___x_239_; 
lean_dec(v_h__3_227_);
lean_dec(v_h__1_225_);
v___x_238_ = lean_box(0);
v___x_239_ = lean_apply_1(v_h__2_226_, v___x_238_);
return v___x_239_;
}
}
}
case 1:
{
lean_dec(v_h__6_230_);
lean_dec(v_h__5_229_);
lean_dec(v_h__4_228_);
lean_dec(v_h__3_227_);
lean_dec(v_h__2_226_);
lean_dec(v_h__1_225_);
switch(v_x_224_)
{
case 0:
{
lean_object* v___x_240_; lean_object* v___x_241_; 
lean_dec(v_h__9_233_);
lean_dec(v_h__8_232_);
v___x_240_ = lean_box(0);
v___x_241_ = lean_apply_1(v_h__7_231_, v___x_240_);
return v___x_241_;
}
case 1:
{
lean_object* v___x_242_; lean_object* v___x_243_; 
lean_dec(v_h__8_232_);
lean_dec(v_h__7_231_);
v___x_242_ = lean_box(0);
v___x_243_ = lean_apply_1(v_h__9_233_, v___x_242_);
return v___x_243_;
}
default: 
{
lean_object* v___x_244_; lean_object* v___x_245_; 
lean_dec(v_h__9_233_);
lean_dec(v_h__7_231_);
v___x_244_ = lean_box(0);
v___x_245_ = lean_apply_1(v_h__8_232_, v___x_244_);
return v___x_245_;
}
}
}
default: 
{
lean_dec(v_h__9_233_);
lean_dec(v_h__8_232_);
lean_dec(v_h__7_231_);
lean_dec(v_h__3_227_);
lean_dec(v_h__2_226_);
lean_dec(v_h__1_225_);
switch(v_x_224_)
{
case 0:
{
lean_object* v___x_246_; lean_object* v___x_247_; 
lean_dec(v_h__6_230_);
lean_dec(v_h__5_229_);
v___x_246_ = lean_box(0);
v___x_247_ = lean_apply_1(v_h__4_228_, v___x_246_);
return v___x_247_;
}
case 1:
{
lean_object* v___x_248_; lean_object* v___x_249_; 
lean_dec(v_h__5_229_);
lean_dec(v_h__4_228_);
v___x_248_ = lean_box(0);
v___x_249_ = lean_apply_1(v_h__6_230_, v___x_248_);
return v___x_249_;
}
default: 
{
lean_object* v___x_250_; lean_object* v___x_251_; 
lean_dec(v_h__6_230_);
lean_dec(v_h__4_228_);
v___x_250_ = lean_box(0);
v___x_251_ = lean_apply_1(v_h__5_229_, v___x_250_);
return v___x_251_;
}
}
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter___boxed(lean_object* v_motive_252_, lean_object* v_x_253_, lean_object* v_x_254_, lean_object* v_h__1_255_, lean_object* v_h__2_256_, lean_object* v_h__3_257_, lean_object* v_h__4_258_, lean_object* v_h__5_259_, lean_object* v_h__6_260_, lean_object* v_h__7_261_, lean_object* v_h__8_262_, lean_object* v_h__9_263_){
_start:
{
uint8_t v_x_143__boxed_264_; uint8_t v_x_144__boxed_265_; lean_object* v_res_266_; 
v_x_143__boxed_264_ = lean_unbox(v_x_253_);
v_x_144__boxed_265_ = lean_unbox(v_x_254_);
v_res_266_ = lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_and_match__1_splitter(v_motive_252_, v_x_143__boxed_264_, v_x_144__boxed_265_, v_h__1_255_, v_h__2_256_, v_h__3_257_, v_h__4_258_, v_h__5_259_, v_h__6_260_, v_h__7_261_, v_h__8_262_, v_h__9_263_);
return v_res_266_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ordering(uint8_t v_x_267_){
_start:
{
switch(v_x_267_)
{
case 0:
{
lean_object* v___x_268_; 
v___x_268_ = lean_unsigned_to_nat(2u);
return v___x_268_;
}
case 1:
{
lean_object* v___x_269_; 
v___x_269_ = lean_unsigned_to_nat(0u);
return v___x_269_;
}
default: 
{
lean_object* v___x_270_; 
v___x_270_ = lean_unsigned_to_nat(1u);
return v___x_270_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_Verdict_ordering___boxed(lean_object* v_x_271_){
_start:
{
uint8_t v_x_34__boxed_272_; lean_object* v_res_273_; 
v_x_34__boxed_272_ = lean_unbox(v_x_271_);
v_res_273_ = lp_BioCompiler_BioCompiler_Verdict_ordering(v_x_34__boxed_272_);
return v_res_273_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter___redArg(uint8_t v_x_274_, lean_object* v_h__1_275_, lean_object* v_h__2_276_, lean_object* v_h__3_277_){
_start:
{
switch(v_x_274_)
{
case 0:
{
lean_object* v___x_278_; lean_object* v___x_279_; 
lean_dec(v_h__3_277_);
lean_dec(v_h__2_276_);
v___x_278_ = lean_box(0);
v___x_279_ = lean_apply_1(v_h__1_275_, v___x_278_);
return v___x_279_;
}
case 1:
{
lean_object* v___x_280_; lean_object* v___x_281_; 
lean_dec(v_h__2_276_);
lean_dec(v_h__1_275_);
v___x_280_ = lean_box(0);
v___x_281_ = lean_apply_1(v_h__3_277_, v___x_280_);
return v___x_281_;
}
default: 
{
lean_object* v___x_282_; lean_object* v___x_283_; 
lean_dec(v_h__3_277_);
lean_dec(v_h__1_275_);
v___x_282_ = lean_box(0);
v___x_283_ = lean_apply_1(v_h__2_276_, v___x_282_);
return v___x_283_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter___redArg___boxed(lean_object* v_x_284_, lean_object* v_h__1_285_, lean_object* v_h__2_286_, lean_object* v_h__3_287_){
_start:
{
uint8_t v_x_36__boxed_288_; lean_object* v_res_289_; 
v_x_36__boxed_288_ = lean_unbox(v_x_284_);
v_res_289_ = lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter___redArg(v_x_36__boxed_288_, v_h__1_285_, v_h__2_286_, v_h__3_287_);
return v_res_289_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter(lean_object* v_motive_290_, uint8_t v_x_291_, lean_object* v_h__1_292_, lean_object* v_h__2_293_, lean_object* v_h__3_294_){
_start:
{
switch(v_x_291_)
{
case 0:
{
lean_object* v___x_295_; lean_object* v___x_296_; 
lean_dec(v_h__3_294_);
lean_dec(v_h__2_293_);
v___x_295_ = lean_box(0);
v___x_296_ = lean_apply_1(v_h__1_292_, v___x_295_);
return v___x_296_;
}
case 1:
{
lean_object* v___x_297_; lean_object* v___x_298_; 
lean_dec(v_h__2_293_);
lean_dec(v_h__1_292_);
v___x_297_ = lean_box(0);
v___x_298_ = lean_apply_1(v_h__3_294_, v___x_297_);
return v___x_298_;
}
default: 
{
lean_object* v___x_299_; lean_object* v___x_300_; 
lean_dec(v_h__3_294_);
lean_dec(v_h__1_292_);
v___x_299_ = lean_box(0);
v___x_300_ = lean_apply_1(v_h__2_293_, v___x_299_);
return v___x_300_;
}
}
}
}
LEAN_EXPORT lean_object* lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter___boxed(lean_object* v_motive_301_, lean_object* v_x_302_, lean_object* v_h__1_303_, lean_object* v_h__2_304_, lean_object* v_h__3_305_){
_start:
{
uint8_t v_x_51__boxed_306_; lean_object* v_res_307_; 
v_x_51__boxed_306_ = lean_unbox(v_x_302_);
v_res_307_ = lp_BioCompiler___private_BioCompiler_ThreeValued_0__BioCompiler_Verdict_ordering_match__1_splitter(v_motive_301_, v_x_51__boxed_306_, v_h__1_303_, v_h__2_304_, v_h__3_305_);
return v_res_307_;
}
}
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_Init(uint8_t builtin);
static bool _G_initialized = false;
LEAN_EXPORT lean_object* initialize_BioCompiler_BioCompiler_ThreeValued(uint8_t builtin) {
lean_object * res;
if (_G_initialized) return lean_io_result_mk_ok(lean_box(0));
_G_initialized = true;
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
