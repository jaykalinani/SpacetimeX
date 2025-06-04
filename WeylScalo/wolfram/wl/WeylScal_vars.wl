(* ::Package:: *)

(* WeylScal_vars.wl *)

(* (c) Liwei Ji, 05/2025 *)

(*******************)

(* Input Variables *)

(*******************)

ADMVarlist =
  GridTensors[
    (*{gam[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[Gamma]"},*)
    {exK[-i, -j], Symmetric[{-i, -j}], PrintAs -> "K"}
  ];

dADMVarlist =
  TileTensors[
    {dgam[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]\[Gamma]"},
    {dexK[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[PartialD]K"}
  ];

ddADMVarlist =
  TileTensors[
    {ddgam[-k, -l, -i, -j], GenSet[Cycles[{1, 2}], Cycles[{3, 4}]], PrintAs -> "\[PartialD]\[PartialD]\[Gamma]"}
  ];

(********************)

(* Output Variables *)

(********************)

Psi4Varlist =
  GridTensors[
    {Psi4real[], PrintAs -> "\!\(\*SuperscriptBox[SubscriptBox[\(\[CapitalPsi]\), \(4\)], \((r)\)]\)"},
    {Psi4imag[], PrintAs -> "\!\(\*SuperscriptBox[SubscriptBox[\(\[CapitalPsi]\), \(4\)], \((i)\)]\)"}
  ];

(**************************)

(* Intermediate Variables *)

(**************************)

IntermediateVarlist =
  TempTensors[
    {detgam[], PrintAs -> "\[Gamma]"},
    {detinvgam[], PrintAs -> "1/\[Gamma]"},
    {sqrtdetgam[], PrintAs -> "\!\(\*SqrtBox[\(\[Gamma]\)]\)"},
    {sqrtdetinvgam[], PrintAs -> "1/\!\(\*SqrtBox[\(\[Gamma]\)]\)"},
    {invgam[i, j], Symmetric[{i, j}], PrintAs -> "\[Gamma]"},
    {GamDDD [-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[CapitalGamma]"},
    {Gam     [k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "\[CapitalGamma]"},
    {Ricc[-i, -j], Symmetric[{-i, -j}], PrintAs -> "R"},
    {trexK[], PrintAs -> "K"},
    {DexK[-k, -i, -j], Symmetric[{-i, -j}], PrintAs -> "DK"}
  ];

EBVarlist =
  TempTensors[
    {Epart[-i, -j], Symmetric[{-i, -j}], PrintAs -> "E"},
    {Bpart[-i, -j], Symmetric[{-i, -j}], PrintAs -> "B"}
  ];

GaussCodazziVarlist =
  TempTensors[
    {Rsca[], PrintAs -> "R"},
    {Riem[-i, -j, -k, -l], GenSet[-Cycles[{1, 2}],-Cycles[{3, 4}]], PrintAs -> "R"},
    {perpR[-i, -j, -k, -l], GenSet[-Cycles[{1, 2}],-Cycles[{3, 4}]], PrintAs -> "\[Perpendicular]R"},
    {perpRn[-j, -k, -l], Antisymmetric[{-k, -l}], PrintAs -> "\[Perpendicular]\!\(\*SubscriptBox[\(R\), \(n\)]\)"},
    {perpRnn[-i, -j], Symmetric[{-i, -j}], PrintAs -> "\[Perpendicular]\!\(\*SubscriptBox[\(R\), \(nn\)]\)"}
  ];

TetradVarlist = 
  TempTensors[
    {uvec[i], PrintAs -> "u"},
    {vvec[i], PrintAs -> "v"},
    {wvec[i], PrintAs -> "w"},
    {inneruutmp[], PrintAs -> "u\[CenterDot]u"},
    {inneruu[], PrintAs -> "u\[CenterDot]u"},
    {Uvec[i], PrintAs -> "U"},
    {innerUv[], PrintAs -> "U\[CenterDot]v"},
    {Vtmp[i]},
    {innerVVtmp[], PrintAs -> "V\[CenterDot]V"},
    {innerVV[], PrintAs -> "V\[CenterDot]V"},
    {Vvec[i], PrintAs -> "V"},
    {innerUw[], PrintAs -> "U\[CenterDot]w"},
    {innerVw[], PrintAs -> "V\[CenterDot]w"},
    {Wtmp[i]},
    {innerWWtmp[], PrintAs -> "W\[CenterDot]W"},
    {innerWW[], PrintAs -> "W\[CenterDot]W"},
    {Wvec[i], PrintAs -> "W"},
    (* tetrad vectors *)
    {mbmbreal[i, j], Symmetric[{i, j}], PrintAs -> "r\!\(\*OverscriptBox[\(m\), \(_\)]\)\!\(\*OverscriptBox[\(m\), \(_\)]\)"},
    {mbmbimag[i, j], Symmetric[{i, j}], PrintAs -> "i\!\(\*OverscriptBox[\(m\), \(_\)]\)\!\(\*OverscriptBox[\(m\), \(_\)]\)"}
  ];

(************)

(* Constant *)

(************)

DefConstantSymbol[oosqrt2, PrintAs -> "1/\!\(\*SqrtBox[\(2\)]\)"];

DefConstantSymbol[ncutoff, PrintAs -> "\!\(\*SuperscriptBox[\(n\), \(cut\)]\)"];
