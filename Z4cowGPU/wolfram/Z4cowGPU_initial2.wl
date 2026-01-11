(* ::Package:: *)

(* Z4cowGPU_initial2.wl *)

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

<<wl/ADM_vars.wl

<<wl/Z4c_vars.wl

<<wl/Z4c_rhs.wl

<<wl/Z4cInADM_vars.wl

<<wl/Z4cInADM_rhs.wl

Module[{Mat, invMat},
  Mat =
    Table[gamt[{ii, -cart}, {jj, -cart}] // ToValues, {ii, 1, 3}, {jj, 1, 3}];
  invMat = Inverse[Mat] /. {1 / Det[Mat] -> 1}; (* det(gamt) = 1 enforced *)
  (*SetEQNDelayed[invdetgamt[], 1 / Det[Mat] // Simplify];*)
  SetEQNDelayed[invgamt[i_, j_], invMat[[i[[1]], j[[1]]]] // Simplify]
];


(******************)
(* Print to Files *)
(******************)

SetOutputFile[FileNameJoin[{Directory[], "Z4cowGPU_initial2.hxx"}]];

SetMainPrint[
  (* Initialize grid function names *)
  PrintInitializations[{Mode -> "MainOut"}, Take[EvolVarlist, {5}]];
  PrintInitializations[{Mode -> "MainIn"}, Take[EvolVarlist, {2}]];
  pr[];

  (* Loops *)
  pr["noinline([&]() __attribute__((__flatten__, __hot__)) {"];
  pr["grid.loop_int_device<0, 0, 0>("];
  pr["  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {"];
  pr["const int ijk = layout2.linear(p.i, p.j, p.k);"];
  pr[];

  (* dgamt *)
  PrintInitializations[{Mode -> "Derivs", DerivsOrder -> 1},
                       Take[dEvolVarlist, {2}]];
  pr[];

  (* invgamt, GtDDD *)
  PrintEquations[{Mode -> "Temp"}, Extract[IntermediateVarlist, {{2}, {5}}]];
  pr[];

  (* trGt *)
  PrintEquations[{Mode -> "MainOut"}, Take[EvolVarlist, {5}]];
  pr[];

  pr["});"];
  pr["});"];
];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetXGPU.wl"}]];
