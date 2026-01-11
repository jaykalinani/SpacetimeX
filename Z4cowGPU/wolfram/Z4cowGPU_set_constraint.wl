(* ::Package:: *)

(* Z4cowGPU_set_constraint.wl *)

(* (c) Liwei Ji, 02/2025 *)

(******************)
(* Configurations *)
(******************)

Needs["xAct`xCoba`", FileNameJoin[{Environment["GENERATO"], "src/Generato.wl"}]]

SetPVerbose[False];

SetPrintDate[False];

SetGridPointIndex["[[ijk]]"];

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

SetOutputFile[FileNameJoin[{Directory[], "Z4cowGPU_set_constraint.hxx"}]];

SetMainPrint[
  (* Initialize grid function names *)
  PrintInitializations[{Mode -> "MainOut"}, Drop[ConstraintVarlist, {-2}]];
  PrintInitializations[{Mode -> "MainIn"}, Drop[TmunuVarlist, 1]];
  PrintInitializations[{Mode -> "MainIn"}, Delete[EvolVarlist, {{1}, {-3}}]];
  pr[];

  (* Loops *)
  pr["noinline([&]() __attribute__((__flatten__, __hot__)) {"];
  pr["grid.loop_int_device<0, 0, 0>("];
  pr["  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {"];
  pr["const int ijk = layout2.linear(p.i, p.j, p.k);"];
  pr[];

  PrintInitializations[{Mode -> "Derivs", DerivsOrder -> 1},
                       Drop[dEvolVarlist, -2]];
  PrintInitializations[{Mode -> "Derivs", DerivsOrder -> 2},
                       Drop[ddEvolVarlist, -2]];
  pr[];

  PrintEquations[{Mode -> "Temp"},
                 Drop[Drop[IntermediateVarlist, {4}], {-4, -2}]];
  PrintEquations[{Mode -> "Temp"}, Drop[DDVarlist, -1]];
  PrintEquations[{Mode -> "Temp"}, RVarlist];
  PrintEquations[{Mode -> "Temp"}, Drop[MatterVarlist, -2]];
  PrintEquations[{Mode -> "Temp"}, dAtUUVarlist];
  pr[];

  PrintEquations[{Mode -> "MainOut"}, ConstraintVarlist];
  pr[];

  pr["});"];
  pr["});"];
];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetXGPU.wl"}]];
