import Lake
open Lake DSL

package BioCompiler where
  moreLeanArgs := #["-DautoImplicit=false"]

@[default_target]
lean_lib BioCompiler where
  roots := #[`BioCompiler.Soundness, `BioCompiler.ThreeValued, `BioCompiler.Sequence, `BioCompiler.NDFST, `BioCompiler.Scanners, `BioCompiler.TypeSystem, `BioCompiler.Compositional, `BioCompiler.Certificates, `BioCompiler.SLOTIndependence, `BioCompiler.Mutagenesis, `BioCompiler.SplicingResolution]
