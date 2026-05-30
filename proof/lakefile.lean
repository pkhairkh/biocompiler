import Lake
open Lake DSL

package BioCompiler where
  leanOptions := #[⟨`autoImplicit, false⟩]

@[default_target]
lean_lib BioCompiler where
  roots := #[`BioCompiler.Soundness]
