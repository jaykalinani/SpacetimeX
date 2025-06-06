(* ::Package:: *)

(* WeylScalo_set_Psi4.wl *)

(* (c) Liwei Ji, 05/2025 *)

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

DefMetric[1, gam[-i, -j], cd, SymbolOfCovD -> {"|", "D"}, PrintAs -> "\[Gamma]"];

MetricInBasis[gam, -cart,
  {
    {gam11[[ijk]], gam12[[ijk]], gam13[[ijk]]},
    {gam12[[ijk]], gam22[[ijk]], gam23[[ijk]]},
    {gam13[[ijk]], gam23[[ijk]], gam33[[ijk]]}
  }];

ComponentValue[X[], coX];
ComponentValue[Y[], coY];
ComponentValue[Z[], coZ];


(**********************************)
(* Define Variables and Equations *)
(**********************************)

<<wl/WeylScal_vars.wl

<<wl/WeylScal_rhs.wl


(******************)
(* Print to Files *)
(******************)

SetOutputFile[FileNameJoin[{Directory[], "WeylScalo_set_Psi4_GF3D5.hxx"}]];

SetMainPrint[
  (* Initialize grid function names *)
  PrintInitializations[{Mode -> "MainOut"}, Psi4Varlist];
  PrintInitializations[{Mode -> "MainIn"}, ADMVarlist];
  pr["const auto &gam11 = gf_gam[0];"];
  pr["const auto &gam12 = gf_gam[1];"];
  pr["const auto &gam13 = gf_gam[2];"];
  pr["const auto &gam22 = gf_gam[3];"];
  pr["const auto &gam23 = gf_gam[4];"];
  pr["const auto &gam33 = gf_gam[5];"];
  pr[];
  PrintInitializations[{Mode -> "Derivs", DerivsOrder -> 1,
                       StorageType -> "Tile", TensorType -> "Vect"},
                       dADMVarlist];
  PrintInitializations[{Mode -> "Derivs", DerivsOrder -> 2,
                       StorageType -> "Tile", TensorType -> "Smat"},
                       ddADMVarlist];
  pr[];


  (* Loops *)
  pr["noinline([&]() __attribute__((__flatten__, __hot__)) {"];
  pr["grid.loop_int_device<0, 0, 0>("];
  pr["  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {"];
  pr["const int ijk = layout2.linear(p.i, p.j, p.k);"];
  pr["const int ijk5 = layout5.linear(p.i, p.j, p.k);"];
  pr[];
  pr["const auto coX = p.x;"];
  pr["const auto coY = p.y;"];
  pr["const auto coZ = p.z;"];
  pr[];

  PrintEquations[{Mode -> "Temp"}, IntermediateVarlist];
  pr[];
  PrintEquations[{Mode -> "Temp"}, EBVarlist];
  pr[];
  PrintEquations[{Mode -> "Temp"}, TetradVarlist];
  pr[];
  PrintEquations[{Mode -> "Main"}, Psi4Varlist];
  pr[];
  pr["  });"];
  pr["});"];
];

Import[FileNameJoin[{Environment["GENERATO"], "codes/CarpetXGPU.wl"}]];
