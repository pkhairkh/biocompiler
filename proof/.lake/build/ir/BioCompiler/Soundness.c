// Lean compiler output
// Module: BioCompiler.Soundness
// Imports: public import Init public meta import Init public import BioCompiler.ThreeValued public import BioCompiler.Sequence public import BioCompiler.NDFST public import BioCompiler.Scanners public import BioCompiler.TypeSystem public import BioCompiler.Compositional public import BioCompiler.Certificates public import BioCompiler.SLOTIndependence
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
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_ThreeValued(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Sequence(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_NDFST(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Scanners(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_TypeSystem(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Compositional(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_Certificates(uint8_t builtin);
lean_object* initialize_BioCompiler_BioCompiler_SLOTIndependence(uint8_t builtin);
static bool _G_initialized = false;
LEAN_EXPORT lean_object* initialize_BioCompiler_BioCompiler_Soundness(uint8_t builtin) {
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
res = initialize_BioCompiler_BioCompiler_SLOTIndependence(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
