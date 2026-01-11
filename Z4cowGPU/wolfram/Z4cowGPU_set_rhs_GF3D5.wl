(* ::Package:: *)

(* Z4cowGPU_set_rhs.wl *)

(* (c) Liwei Ji, 02/2025 *)

(******************)
(* Configurations *)
(******************)

Needs["xAct`xCoba`", FileNameJoin[{Environment["GENERATO"], "src/Generato.wl"}]]

SetPVerbose[False];

SetPrintDate[False];

SetGridPointIndex["[[ijk]]"];

SetTilePointIndex["[[ijk5]]"];

(*SetUseLetterForTensorComponet[True];*)

SetTempVariableType["auto"];

DefManifold[M3, 3, IndexRange[a, z]];

DefChart[cart, M3, {1, 2, 3}, {X[], Y[], Z[]}, ChartColor -> Blue];


(**********************************)
(* Define Variables and Equations *)
(**********************************)

<<wl/Z4c_vars.wl

<<wl/Z4c_rhs.wl

Module[{Mat, invMat},
  Mat =
    Table[gamt[{ii, -cart}, {jj, -cart}] // ToValues, {ii, 1, 3}, {jj, 1, 3}];
  invMat = Inverse[Mat] /. {1 / Det[Mat] -> 1}; (* det(gamt) = 1 enforced *)
  (* SetEQNDelayed[detinvgamt[], 1 / Det[Mat] // Simplify]; *)
  SetEQNDelayed[invgamt[i_, j_], invMat[[i[[1]], j[[1]]]] // Simplify]
];


(******************)
(* Print to Files *)
(******************)

SetOutputFile[FileNameJoin[{Directory[], "Z4cowGPU_set_rhs_GF3D5.hxx"}]];

SetMainPrint[
  (* Initialize grid function names *)
  PrintInitializations[{Mode -> "MainOut"}, dtEvolVarlist];
  PrintInitializations[{Mode -> "MainIn"}, Drop[TmunuVarlist, 1]];
  PrintInitializations[{Mode -> "MainIn"}, Delete[EvolVarlist, {{1}, {-3}}]];
  pr[];
  PrintInitializations[{Mode -> "Derivs", DerivsOrder -> 1,
                        StorageType -> "Tile", TensorType -> "Vect"},
                       dEvolVarlist];
  PrintInitializations[{Mode -> "Derivs", DerivsOrder -> 2,
                        StorageType -> "Tile", TensorType -> "Smat"},
                       ddEvolVarlist];
  pr[];


  (* Loops *)
  pr["noinline([&]() __attribute__((__flatten__, __hot__)) {"];
  pr["grid.loop_int_device<0, 0, 0>("];
  pr["  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {"];
  pr["const int ijk = layout2.linear(p.i, p.j, p.k);"];
  pr["const int ijk5 = layout5.linear(p.i, p.j, p.k);"];
  pr[];

  pr["const CCTK_REAL ceta = calceta(p.x, p.y, p.z);"];
  pr[];

  PrintEquations[{Mode -> "Temp"}, IntermediateVarlist];
  PrintEquations[{Mode -> "Temp"}, DDVarlist];
  PrintEquations[{Mode -> "Temp"}, RVarlist];
  PrintEquations[{Mode -> "Temp"}, MatterVarlist];
  pr[];

  PrintEquations[{Mode -> "MainOut"}, dtEvolVarlist];
  pr[];

  pr["});"];
  pr["});"];
];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetXGPU.wl"}]];
