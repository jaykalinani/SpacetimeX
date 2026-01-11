(* ::Package:: *)

(* Z4cowGPU_enforce.wl *)

(* (c) Liwei Ji, 07/2024 *)

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

<<wl/Enforce_vars.wl

<<wl/Enforce_rhs.wl

Module[{Mat},
  Mat = Table[gamtold[{ii, -cart}, {jj, -cart}] // ToValues,
              {ii, 1, 3}, {jj, 1, 3}];
  SetEQNDelayed[invdetgamtold[], 1 / Det[Mat] // Simplify];
];

Module[{Mat, invMat},
  Mat =
    Table[gamt[{ii, -cart}, {jj, -cart}] // ToValues, {ii, 1, 3}, {jj, 1, 3}];
  invMat = Inverse[Mat] /. {1 / Det[Mat] -> 1}; (* det(gamt) = 1 enforced *)
  (* SetEQNDelayed[invdetgamt[], 1 / Det[Mat] // Simplify]; *)
  SetEQNDelayed[invgamt[i_, j_], invMat[[i[[1]], j[[1]]]] // Simplify]
];


(******************)
(* Print to Files *)
(******************)

SetOutputFile[FileNameJoin[{Directory[], "Z4cowGPU_enforce.hxx"}]];

SetMainPrint[
  (* Initialize grid function names *)
  PrintInitializations[{Mode -> "MainOut"},
                       Extract[EvolVarlist, {{2}, {4}, {-2}}]];
  pr[];

  (* Loops *)
  pr["noinline([&]() __attribute__((__flatten__, __hot__)) {"];
  pr["grid.loop_all_device<0, 0, 0>("];
  pr["  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {"];
  pr["const int ijk = layout2.linear(p.i, p.j, p.k);"];
  pr[];

  (* Wold, gamtold, exAtold, alphaold *)
  PrintEquations[{Mode -> "Temp"}, Drop[EnforceVarlist, -3]];

  pr["W[ijk] = fmax(cWfloor, Wold);"];
  pr["alpha[ijk] = fmax(calphafloor, alphaold);"];
  pr[];

  (* invdetgamtold *)
  PrintEquations[{Mode -> "Temp"}, Take[EnforceVarlist, {-3}]];

  pr["const auto Wfac = cbrt(sqrt(invdetgamtold));"];
  pr[];

  (* gamt *)
  PrintEquations[{Mode -> "MainOut"}, Take[EvolVarlist, {2}]];
  pr[];

  (* invgamt *)
  PrintEquations[{Mode -> "Temp"}, Take[IntermediateVarlist, {2}]];
  (* trexAtold *)
  PrintEquations[{Mode -> "Temp"}, Take[EnforceVarlist, {-1}]];
  pr[];

  (* exAt *)
  PrintEquations[{Mode -> "MainOut"}, Take[EvolVarlist, {4}]];
  pr[];

  pr["});"];
  pr["});"];
];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetXGPU.wl"}]];
