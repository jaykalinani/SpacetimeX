(* ::Package:: *)

(* WeylScal_rhs.wl *)

(* (c) Liwei Ji, 05/2025 *)

(***********)

(* Metrics *)

(***********)

Module[{Mat, invMat},
  Mat = Table[gam[{ii, -GetDefaultChart[]}, {jj, -GetDefaultChart[]}] // ToValues, {ii, 1, 3}, {jj, 1, 3}];
  invMat = Inverse[Mat] /. {1 / Det[Mat] -> detinvgam};
  SetEQNDelayed[detgam[], Det[Mat] // Simplify];
  SetEQNDelayed[detinvgam[], 1 / detgam[]];
  SetEQNDelayed[sqrtdetgam[], (detgam[]) ^ (1/2)];
  SetEQNDelayed[sqrtdetinvgam[], 1 / sqrtdetgam[]];
  SetEQNDelayed[invgam[i_, j_], invMat[[i[[1]], j[[1]]]] // Simplify]
];

(****************)

(* Intermediate *)

(****************)

SetEQN[GamDDD[k_, i_, j_], 1/2 (dgam[i, j, k] + dgam[j, k, i] - dgam[k, i, j])];

SetEQN[Gam[k_, i_, j_], invgam[k, l] GamDDD[-l, i, j]];

SetEQN[Ricc[i_, j_], invgam[k, l] (1/2 (-ddgam[-k, -l, i, j] - ddgam[i, j, -k, -l] + ddgam[i, -k, j, -l] + ddgam[j, -k, i, -l]) + Gam[m, i, -k] GamDDD[-m, j, -l] - Gam[m, i, j] GamDDD[-m, -k, -l])];

SetEQN[trexK[], invgam[k, l] exK[-k, -l]];

SetEQN[DexK[k_, i_, j_], dexK[k, i, j] - Gam[l, k, i] exK[-l, j] - Gam[l, k, j] exK[-l, i]];

(********************************)

(* From Gauss-Codazzi equations *)

(********************************)

SetEQN[Rsca[], invgam[k, l] Ricc[-k, -l]];

SetEQN[Riem[i_, j_, k_, l_], gam[i, k] Ricc[j, l] + gam[j, l] Ricc[i, k] - 1/2 Rsca[] gam[i, k] gam[j, l] - gam[i, l] Ricc[j, k] - gam[j, k] Ricc[i, l] + 1/2 Rsca[] gam[i, l] gam[j, k]];

SetEQN[perpR[i_, j_, k_, l_], Riem[i, j, k, l] + exK[i, k] exK[j, l] - exK[i, l] exK[j, k]];

SetEQN[perpRn[j_, k_, l_], DexK[k, j, l] - DexK[l, j, k]];

SetEQN[perpRnn[i_, j_], Ricc[i, j] + trexK[] exK[i, j] - exK[i, -k] exK[j, -l] invgam[k, l]];

(********************************)

(* E and B parts of Weyl tensor *)

(********************************)

SetEQNDelayed[Epart[i_, j_], Ricc[i, j] + trexK[] exK[i, j] - exK[i, -k] exK[j, -l] invgam[k, l]];

SetEQNDelayed[Bpart[i_, j_], (gam[i, -m] epsilongam[m, k, l] DexK[-k, -l, j] /. epsilonToetaUp[gam, GetDefaultChart[]]) /. {1 / Sqrt[Detgamcart[]] -> sqrtdetinvgam[]}]; (* here we use delayed version to prevent contracting gam and epsilongam *)

(**********)

(* Tetrad *)

(**********)

SetEQNDelayed[
  uvec[i_]
  ,
  Switch[i[[1]],
    1,
      -Y[]
    ,
    2,
      X[]
    ,
    3,
      0
    ,
    _,
      Missing["InvalidIndex", i[[1]]]
  ]
];

SetEQNDelayed[
  vvec[i_]
  ,
  Switch[i[[1]],
    1,
      X[]
    ,
    2,
      Y[]
    ,
    3,
      Z[]
    ,
    _,
      Missing["InvalidIndex", i[[1]]]
  ]
];

SetEQN[wvec[i_], (invgam[i, j] epsilongam[-j, -k, -l] uvec[k] vvec[l]) /. epsilonToetaDown[gam, GetDefaultChart[]] /. {Sqrt[Detgamcart[]] -> sqrtdetgam[]}];

(* Gram-Schmidt orthonormalization *)

SetEQN[inneruutmp[], gam[-k, -l] uvec[k] uvec[l]];

SetEQN[inneruu[], Max[ncutoff, inneruutmp[]]];

SetEQN[Uvec[i_], uvec[i] / (inneruu[]) ^ (1/2)];

SetEQN[innerUv[], gam[-k, -l] Uvec[k] vvec[l]];

SetEQN[Vtmp[i_], vvec[i] - innerUv[] Uvec[i]];

SetEQN[innerVVtmp[], gam[-k, -l] Vtmp[k] Vtmp[l]];

SetEQN[innerVV[], Max[ncutoff, innerVVtmp[]]];

SetEQN[Vvec[i_], Vtmp[i] / (innerVV[]) ^ (1/2)];

SetEQN[innerUw[], gam[-k, -l] Uvec[k] wvec[l]];

SetEQN[innerVw[], gam[-k, -l] Vvec[k] wvec[l]];

SetEQN[Wtmp[i_], wvec[i] - innerUw[] Uvec[i] - innerVw[] Vvec[i]];

SetEQN[innerWWtmp[], gam[-k, -l] Wtmp[k] Wtmp[l]];

SetEQN[innerWW[], Max[ncutoff, innerWWtmp[]]];

SetEQN[Wvec[i_], Wtmp[i] / (innerWW[]) ^ (1/2)];

(* The tetrad vectors *)

SetEQN[mbmbreal[i_, j_], (Wvec[i] Wvec[j] - Uvec[i] Uvec[j]) / 2];

SetEQN[mbmbimag[i_, j_], -(Uvec[i] Wvec[j] + Wvec[i] Uvec[j]) / 2];

(* Psi4 *)

SetEQN[Psi4real[], Epart[-i, -j] mbmbreal[i, j] + Bpart[-i, -j] mbmbimag[i, j]];

SetEQN[Psi4imag[], Epart[-i, -j] mbmbimag[i, j] - Bpart[-i, -j] mbmbreal[i, j]];

(* Psi4 using Gauss-Codazzi *)

SetEQN[{SuffixName -> "GaussCodazzi"}, Psi4real[], (perpR[-i, -j, -k, -l] Vvec[i] Vvec[k] - 2 perpRn[-j, -k, -l] Vvec[k] + perpRnn[-j, -l]) mbmbreal[j, l] / 2];

SetEQN[{SuffixName -> "GaussCodazzi"}, Psi4imag[], (perpR[-i, -j, -k, -l] Vvec[i] Vvec[k] - 2 perpRn[-j, -k, -l] Vvec[k] + perpRnn[-j, -l]) mbmbimag[j, l] / 2];
