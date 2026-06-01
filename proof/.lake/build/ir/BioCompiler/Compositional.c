// Lean compiler output
// Module: BioCompiler.Compositional
// Imports: public import Init public meta import Init public import BioCompiler.ThreeValued public import BioCompiler.Sequence public import BioCompiler.NDFST public import BioCompiler.Scanners public import BioCompiler.TypeSystem
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
uint8_t lp_BioCompiler_BioCompiler_evaluate___redArg(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
lean_object* lp_BioCompiler_BioCompiler_Verdict_and___boxed(lean_object*, lean_object*);
lean_object* l_List_mapTR_loop___redArg(lean_object*, lean_object*, lean_object*);
lean_object* l_List_foldl___redArg(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluateAll___redArg___lam__0(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluateAll___redArg___lam__0___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
static const lean_closure_object lp_BioCompiler_BioCompiler_evaluateAll___redArg___closed__0_value = {.m_header = {.m_rc = 0, .m_cs_sz = sizeof(lean_closure_object) + sizeof(void*)*0, .m_other = 0, .m_tag = 245}, .m_fun = (void*)lp_BioCompiler_BioCompiler_Verdict_and___boxed, .m_arity = 2, .m_num_fixed = 0, .m_objs = {} };
static const lean_object* lp_BioCompiler_BioCompiler_evaluateAll___redArg___closed__0 = (const lean_object*)&lp_BioCompiler_BioCompiler_evaluateAll___redArg___closed__0_value;
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluateAll___redArg(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluateAll___redArg___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluateAll(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluateAll___boxed(lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*, lean_object*);
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluateAll___redArg___lam__0(lean_object* v_inst__splice_1_, lean_object* v_inst__cai_2_, lean_object* v_inst__cpg_3_, lean_object* v_inst__ndfst_4_, lean_object* v_seq_5_, lean_object* v_ctx_6_, lean_object* v_P_7_){
_start:
{
uint8_t v___x_8_; 
v___x_8_ = lp_BioCompiler_BioCompiler_evaluate___redArg(v_inst__splice_1_, v_inst__cai_2_, v_inst__cpg_3_, v_inst__ndfst_4_, v_P_7_, v_seq_5_, v_ctx_6_);
return v___x_8_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluateAll___redArg___lam__0___boxed(lean_object* v_inst__splice_9_, lean_object* v_inst__cai_10_, lean_object* v_inst__cpg_11_, lean_object* v_inst__ndfst_12_, lean_object* v_seq_13_, lean_object* v_ctx_14_, lean_object* v_P_15_){
_start:
{
uint8_t v_res_16_; lean_object* v_r_17_; 
v_res_16_ = lp_BioCompiler_BioCompiler_evaluateAll___redArg___lam__0(v_inst__splice_9_, v_inst__cai_10_, v_inst__cpg_11_, v_inst__ndfst_12_, v_seq_13_, v_ctx_14_, v_P_15_);
lean_dec_ref(v_ctx_14_);
v_r_17_ = lean_box(v_res_16_);
return v_r_17_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluateAll___redArg(lean_object* v_inst__splice_19_, lean_object* v_inst__cai_20_, lean_object* v_inst__cpg_21_, lean_object* v_inst__ndfst_22_, lean_object* v_predicates_23_, lean_object* v_seq_24_, lean_object* v_ctx_25_){
_start:
{
lean_object* v___f_26_; lean_object* v___x_27_; uint8_t v___x_28_; lean_object* v___x_29_; lean_object* v___x_30_; lean_object* v___x_31_; lean_object* v___x_32_; uint8_t v___x_33_; 
v___f_26_ = lean_alloc_closure((void*)(lp_BioCompiler_BioCompiler_evaluateAll___redArg___lam__0___boxed), 7, 6);
lean_closure_set(v___f_26_, 0, v_inst__splice_19_);
lean_closure_set(v___f_26_, 1, v_inst__cai_20_);
lean_closure_set(v___f_26_, 2, v_inst__cpg_21_);
lean_closure_set(v___f_26_, 3, v_inst__ndfst_22_);
lean_closure_set(v___f_26_, 4, v_seq_24_);
lean_closure_set(v___f_26_, 5, v_ctx_25_);
v___x_27_ = ((lean_object*)(lp_BioCompiler_BioCompiler_evaluateAll___redArg___closed__0));
v___x_28_ = 0;
v___x_29_ = lean_box(0);
v___x_30_ = l_List_mapTR_loop___redArg(v___f_26_, v_predicates_23_, v___x_29_);
v___x_31_ = lean_box(v___x_28_);
v___x_32_ = l_List_foldl___redArg(v___x_27_, v___x_31_, v___x_30_);
v___x_33_ = lean_unbox(v___x_32_);
lean_dec(v___x_32_);
return v___x_33_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluateAll___redArg___boxed(lean_object* v_inst__splice_34_, lean_object* v_inst__cai_35_, lean_object* v_inst__cpg_36_, lean_object* v_inst__ndfst_37_, lean_object* v_predicates_38_, lean_object* v_seq_39_, lean_object* v_ctx_40_){
_start:
{
uint8_t v_res_41_; lean_object* v_r_42_; 
v_res_41_ = lp_BioCompiler_BioCompiler_evaluateAll___redArg(v_inst__splice_34_, v_inst__cai_35_, v_inst__cpg_36_, v_inst__ndfst_37_, v_predicates_38_, v_seq_39_, v_ctx_40_);
v_r_42_ = lean_box(v_res_41_);
return v_r_42_;
}
}
LEAN_EXPORT uint8_t lp_BioCompiler_BioCompiler_evaluateAll(lean_object* v_inst__splice_43_, lean_object* v_inst__cai_44_, lean_object* v_inst__cpg_45_, lean_object* v_State_46_, lean_object* v_inst__dec_47_, lean_object* v_inst__inhab_48_, lean_object* v_inst__ndfst_49_, lean_object* v_predicates_50_, lean_object* v_seq_51_, lean_object* v_ctx_52_){
_start:
{
uint8_t v___x_53_; 
v___x_53_ = lp_BioCompiler_BioCompiler_evaluateAll___redArg(v_inst__splice_43_, v_inst__cai_44_, v_inst__cpg_45_, v_inst__ndfst_49_, v_predicates_50_, v_seq_51_, v_ctx_52_);
return v___x_53_;
}
}
LEAN_EXPORT lean_object* lp_BioCompiler_BioCompiler_evaluateAll___boxed(lean_object* v_inst__splice_54_, lean_object* v_inst__cai_55_, lean_object* v_inst__cpg_56_, lean_object* v_State_57_, lean_object* v_inst__dec_58_, lean_object* v_inst__inhab_59_, lean_object* v_inst__ndfst_60_, lean_object* v_predicates_61_, lean_object* v_seq_62_, lean_object* v_ctx_63_){
_start:
{
uint8_t v_res_64_; lean_object* v_r_65_; 
v_res_64_ = lp_BioCompiler_BioCompiler_evaluateAll(v_inst__splice_54_, v_inst__cai_55_, v_inst__cpg_56_, v_State_57_, v_inst__dec_58_, v_inst__inhab_59_, v_inst__ndfst_60_, v_predicates_61_, v_seq_62_, v_ctx_63_);
lean_dec(v_inst__inhab_59_);
lean_dec_ref(v_inst__dec_58_);
v_r_65_ = lean_box(v_res_64_);
return v_r_65_;
}
}
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_ThreeValued(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Sequence(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_NDFST(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Scanners(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_TypeSystem(uint8_t builtin);
static bool _G_initialized = false;
LEAN_EXPORT lean_object* initialize_BioCompiler_BioCompiler_Compositional(uint8_t builtin) {
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
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
