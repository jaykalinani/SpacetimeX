/* WeylScalo_set_Psi4_GaussCodazzi_GF3D5.hxx */
/* Produced with Generato */

#ifndef WEYLSCALO_SET_PSI4_GAUSSCODAZZI_GF3D5_HXX
#define WEYLSCALO_SET_PSI4_GAUSSCODAZZI_GF3D5_HXX

const auto Psi4real = gf_Psi4real;
const auto Psi4imag = gf_Psi4imag;
const auto exK11 = gf_exK[0];
const auto exK12 = gf_exK[1];
const auto exK13 = gf_exK[2];
const auto exK22 = gf_exK[3];
const auto exK23 = gf_exK[4];
const auto exK33 = gf_exK[5];
const auto &gam11 = gf_gam[0];
const auto &gam12 = gf_gam[1];
const auto &gam13 = gf_gam[2];
const auto &gam22 = gf_gam[3];
const auto &gam23 = gf_gam[4];
const auto &gam33 = gf_gam[5];

const auto dgam111 = tl_dgam[0][0].ptr;
const auto dgam112 = tl_dgam[1][0].ptr;
const auto dgam113 = tl_dgam[2][0].ptr;
const auto dgam122 = tl_dgam[3][0].ptr;
const auto dgam123 = tl_dgam[4][0].ptr;
const auto dgam133 = tl_dgam[5][0].ptr;
const auto dgam211 = tl_dgam[0][1].ptr;
const auto dgam212 = tl_dgam[1][1].ptr;
const auto dgam213 = tl_dgam[2][1].ptr;
const auto dgam222 = tl_dgam[3][1].ptr;
const auto dgam223 = tl_dgam[4][1].ptr;
const auto dgam233 = tl_dgam[5][1].ptr;
const auto dgam311 = tl_dgam[0][2].ptr;
const auto dgam312 = tl_dgam[1][2].ptr;
const auto dgam313 = tl_dgam[2][2].ptr;
const auto dgam322 = tl_dgam[3][2].ptr;
const auto dgam323 = tl_dgam[4][2].ptr;
const auto dgam333 = tl_dgam[5][2].ptr;
const auto dexK111 = tl_dexK[0][0].ptr;
const auto dexK112 = tl_dexK[1][0].ptr;
const auto dexK113 = tl_dexK[2][0].ptr;
const auto dexK122 = tl_dexK[3][0].ptr;
const auto dexK123 = tl_dexK[4][0].ptr;
const auto dexK133 = tl_dexK[5][0].ptr;
const auto dexK211 = tl_dexK[0][1].ptr;
const auto dexK212 = tl_dexK[1][1].ptr;
const auto dexK213 = tl_dexK[2][1].ptr;
const auto dexK222 = tl_dexK[3][1].ptr;
const auto dexK223 = tl_dexK[4][1].ptr;
const auto dexK233 = tl_dexK[5][1].ptr;
const auto dexK311 = tl_dexK[0][2].ptr;
const auto dexK312 = tl_dexK[1][2].ptr;
const auto dexK313 = tl_dexK[2][2].ptr;
const auto dexK322 = tl_dexK[3][2].ptr;
const auto dexK323 = tl_dexK[4][2].ptr;
const auto dexK333 = tl_dexK[5][2].ptr;
const auto ddgam1111 = tl_ddgam[0][0].ptr;
const auto ddgam1112 = tl_ddgam[1][0].ptr;
const auto ddgam1113 = tl_ddgam[2][0].ptr;
const auto ddgam1122 = tl_ddgam[3][0].ptr;
const auto ddgam1123 = tl_ddgam[4][0].ptr;
const auto ddgam1133 = tl_ddgam[5][0].ptr;
const auto ddgam1211 = tl_ddgam[0][1].ptr;
const auto ddgam1212 = tl_ddgam[1][1].ptr;
const auto ddgam1213 = tl_ddgam[2][1].ptr;
const auto ddgam1222 = tl_ddgam[3][1].ptr;
const auto ddgam1223 = tl_ddgam[4][1].ptr;
const auto ddgam1233 = tl_ddgam[5][1].ptr;
const auto ddgam1311 = tl_ddgam[0][2].ptr;
const auto ddgam1312 = tl_ddgam[1][2].ptr;
const auto ddgam1313 = tl_ddgam[2][2].ptr;
const auto ddgam1322 = tl_ddgam[3][2].ptr;
const auto ddgam1323 = tl_ddgam[4][2].ptr;
const auto ddgam1333 = tl_ddgam[5][2].ptr;
const auto ddgam2211 = tl_ddgam[0][3].ptr;
const auto ddgam2212 = tl_ddgam[1][3].ptr;
const auto ddgam2213 = tl_ddgam[2][3].ptr;
const auto ddgam2222 = tl_ddgam[3][3].ptr;
const auto ddgam2223 = tl_ddgam[4][3].ptr;
const auto ddgam2233 = tl_ddgam[5][3].ptr;
const auto ddgam2311 = tl_ddgam[0][4].ptr;
const auto ddgam2312 = tl_ddgam[1][4].ptr;
const auto ddgam2313 = tl_ddgam[2][4].ptr;
const auto ddgam2322 = tl_ddgam[3][4].ptr;
const auto ddgam2323 = tl_ddgam[4][4].ptr;
const auto ddgam2333 = tl_ddgam[5][4].ptr;
const auto ddgam3311 = tl_ddgam[0][5].ptr;
const auto ddgam3312 = tl_ddgam[1][5].ptr;
const auto ddgam3313 = tl_ddgam[2][5].ptr;
const auto ddgam3322 = tl_ddgam[3][5].ptr;
const auto ddgam3323 = tl_ddgam[4][5].ptr;
const auto ddgam3333 = tl_ddgam[5][5].ptr;

noinline([&]() __attribute__((__flatten__, __hot__)) {
grid.loop_int_device<0, 0, 0>(
  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
const int ijk = layout2.linear(p.i, p.j, p.k);
const int ijk5 = layout5.linear(p.i, p.j, p.k);

const auto coX = p.x;
const auto coY = p.y;
const auto coZ = p.z;

const auto
detgam
=
-(Power(gam13[ijk],2)*gam22[ijk]) + 2*gam12[ijk]*gam13[ijk]*gam23[ijk] -
  Power(gam12[ijk],2)*gam33[ijk] +
  gam11[ijk]*(-Power(gam23[ijk],2) + gam22[ijk]*gam33[ijk])
;

const auto
detinvgam
=
1/detgam
;

const auto
sqrtdetgam
=
Sqrt(detgam)
;

const auto
sqrtdetinvgam
=
1/sqrtdetgam
;

const auto
invgam11
=
detinvgam*(-Power(gam23[ijk],2) + gam22[ijk]*gam33[ijk])
;

const auto
invgam12
=
detinvgam*(gam13[ijk]*gam23[ijk] - gam12[ijk]*gam33[ijk])
;

const auto
invgam13
=
detinvgam*(-(gam13[ijk]*gam22[ijk]) + gam12[ijk]*gam23[ijk])
;

const auto
invgam22
=
detinvgam*(-Power(gam13[ijk],2) + gam11[ijk]*gam33[ijk])
;

const auto
invgam23
=
detinvgam*(gam12[ijk]*gam13[ijk] - gam11[ijk]*gam23[ijk])
;

const auto
invgam33
=
detinvgam*(-Power(gam12[ijk],2) + gam11[ijk]*gam22[ijk])
;

const auto
GamDDD111
=
dgam111[ijk5]/2.
;

const auto
GamDDD112
=
dgam211[ijk5]/2.
;

const auto
GamDDD113
=
dgam311[ijk5]/2.
;

const auto
GamDDD122
=
-0.5*dgam122[ijk5] + dgam212[ijk5]
;

const auto
GamDDD123
=
(-dgam123[ijk5] + dgam213[ijk5] + dgam312[ijk5])/2.
;

const auto
GamDDD133
=
-0.5*dgam133[ijk5] + dgam313[ijk5]
;

const auto
GamDDD211
=
dgam112[ijk5] - dgam211[ijk5]/2.
;

const auto
GamDDD212
=
dgam122[ijk5]/2.
;

const auto
GamDDD213
=
(dgam123[ijk5] - dgam213[ijk5] + dgam312[ijk5])/2.
;

const auto
GamDDD222
=
dgam222[ijk5]/2.
;

const auto
GamDDD223
=
dgam322[ijk5]/2.
;

const auto
GamDDD233
=
-0.5*dgam233[ijk5] + dgam323[ijk5]
;

const auto
GamDDD311
=
dgam113[ijk5] - dgam311[ijk5]/2.
;

const auto
GamDDD312
=
(dgam123[ijk5] + dgam213[ijk5] - dgam312[ijk5])/2.
;

const auto
GamDDD313
=
dgam133[ijk5]/2.
;

const auto
GamDDD322
=
dgam223[ijk5] - dgam322[ijk5]/2.
;

const auto
GamDDD323
=
dgam233[ijk5]/2.
;

const auto
GamDDD333
=
dgam333[ijk5]/2.
;

const auto
Gam111
=
GamDDD111*invgam11 + GamDDD211*invgam12 + GamDDD311*invgam13
;

const auto
Gam112
=
GamDDD112*invgam11 + GamDDD212*invgam12 + GamDDD312*invgam13
;

const auto
Gam113
=
GamDDD113*invgam11 + GamDDD213*invgam12 + GamDDD313*invgam13
;

const auto
Gam122
=
GamDDD122*invgam11 + GamDDD222*invgam12 + GamDDD322*invgam13
;

const auto
Gam123
=
GamDDD123*invgam11 + GamDDD223*invgam12 + GamDDD323*invgam13
;

const auto
Gam133
=
GamDDD133*invgam11 + GamDDD233*invgam12 + GamDDD333*invgam13
;

const auto
Gam211
=
GamDDD111*invgam12 + GamDDD211*invgam22 + GamDDD311*invgam23
;

const auto
Gam212
=
GamDDD112*invgam12 + GamDDD212*invgam22 + GamDDD312*invgam23
;

const auto
Gam213
=
GamDDD113*invgam12 + GamDDD213*invgam22 + GamDDD313*invgam23
;

const auto
Gam222
=
GamDDD122*invgam12 + GamDDD222*invgam22 + GamDDD322*invgam23
;

const auto
Gam223
=
GamDDD123*invgam12 + GamDDD223*invgam22 + GamDDD323*invgam23
;

const auto
Gam233
=
GamDDD133*invgam12 + GamDDD233*invgam22 + GamDDD333*invgam23
;

const auto
Gam311
=
GamDDD111*invgam13 + GamDDD211*invgam23 + GamDDD311*invgam33
;

const auto
Gam312
=
GamDDD112*invgam13 + GamDDD212*invgam23 + GamDDD312*invgam33
;

const auto
Gam313
=
GamDDD113*invgam13 + GamDDD213*invgam23 + GamDDD313*invgam33
;

const auto
Gam322
=
GamDDD122*invgam13 + GamDDD222*invgam23 + GamDDD322*invgam33
;

const auto
Gam323
=
GamDDD123*invgam13 + GamDDD223*invgam23 + GamDDD323*invgam33
;

const auto
Gam333
=
GamDDD133*invgam13 + GamDDD233*invgam23 + GamDDD333*invgam33
;

const auto
Ricc11
=
Gam112*GamDDD111*invgam12 - Gam111*GamDDD112*invgam12 +
  Gam212*GamDDD211*invgam12 - Gam211*GamDDD212*invgam12 +
  Gam312*GamDDD311*invgam12 - Gam311*GamDDD312*invgam12 +
  Gam113*GamDDD111*invgam13 - Gam111*GamDDD113*invgam13 +
  Gam213*GamDDD211*invgam13 - Gam211*GamDDD213*invgam13 +
  Gam313*GamDDD311*invgam13 - Gam311*GamDDD313*invgam13 +
  Gam112*GamDDD112*invgam22 - Gam111*GamDDD122*invgam22 +
  Gam212*GamDDD212*invgam22 - Gam211*GamDDD222*invgam22 +
  Gam312*GamDDD312*invgam22 - Gam311*GamDDD322*invgam22 +
  Gam113*GamDDD112*invgam23 + Gam112*GamDDD113*invgam23 -
  2*Gam111*GamDDD123*invgam23 + Gam213*GamDDD212*invgam23 +
  Gam212*GamDDD213*invgam23 - 2*Gam211*GamDDD223*invgam23 +
  Gam313*GamDDD312*invgam23 + Gam312*GamDDD313*invgam23 -
  2*Gam311*GamDDD323*invgam23 + Gam113*GamDDD113*invgam33 -
  Gam111*GamDDD133*invgam33 + Gam213*GamDDD213*invgam33 -
  Gam211*GamDDD233*invgam33 + Gam313*GamDDD313*invgam33 -
  Gam311*GamDDD333*invgam33 - (invgam22*ddgam1122[ijk5])/2. -
  invgam23*ddgam1123[ijk5] - (invgam33*ddgam1133[ijk5])/2. +
  invgam22*ddgam1212[ijk5] + invgam23*ddgam1213[ijk5] +
  invgam23*ddgam1312[ijk5] + invgam33*ddgam1313[ijk5] -
  (invgam22*ddgam2211[ijk5])/2. - invgam23*ddgam2311[ijk5] -
  (invgam33*ddgam3311[ijk5])/2.
;

const auto
Ricc12
=
(-2*Gam112*GamDDD111*invgam11 + 2*Gam111*GamDDD112*invgam11 -
    2*Gam212*GamDDD211*invgam11 + 2*Gam211*GamDDD212*invgam11 -
    2*Gam312*GamDDD311*invgam11 + 2*Gam311*GamDDD312*invgam11 -
    2*Gam112*GamDDD112*invgam12 + 2*Gam111*GamDDD122*invgam12 -
    2*Gam212*GamDDD212*invgam12 + 2*Gam211*GamDDD222*invgam12 -
    2*Gam312*GamDDD312*invgam12 + 2*Gam311*GamDDD322*invgam12 +
    2*Gam113*GamDDD112*invgam13 - 4*Gam112*GamDDD113*invgam13 +
    2*Gam111*GamDDD123*invgam13 + 2*Gam213*GamDDD212*invgam13 -
    4*Gam212*GamDDD213*invgam13 + 2*Gam211*GamDDD223*invgam13 +
    2*Gam313*GamDDD312*invgam13 - 4*Gam312*GamDDD313*invgam13 +
    2*Gam311*GamDDD323*invgam13 + 2*Gam113*GamDDD122*invgam23 -
    2*Gam112*GamDDD123*invgam23 + 2*Gam213*GamDDD222*invgam23 -
    2*Gam212*GamDDD223*invgam23 + 2*Gam313*GamDDD322*invgam23 -
    2*Gam312*GamDDD323*invgam23 + 2*Gam113*GamDDD123*invgam33 -
    2*Gam112*GamDDD133*invgam33 + 2*Gam213*GamDDD223*invgam33 -
    2*Gam212*GamDDD233*invgam33 + 2*Gam313*GamDDD323*invgam33 -
    2*Gam312*GamDDD333*invgam33 + invgam12*ddgam1122[ijk5] +
    invgam13*ddgam1123[ijk5] - 2*invgam12*ddgam1212[ijk5] -
    invgam13*ddgam1213[ijk5] - invgam23*ddgam1223[ijk5] -
    invgam33*ddgam1233[ijk5] - invgam13*ddgam1312[ijk5] +
    invgam23*ddgam1322[ijk5] + invgam33*ddgam1323[ijk5] +
    invgam12*ddgam2211[ijk5] + invgam23*ddgam2213[ijk5] +
    invgam13*ddgam2311[ijk5] - invgam23*ddgam2312[ijk5] +
    invgam33*ddgam2313[ijk5] - invgam33*ddgam3312[ijk5])/2.
;

const auto
Ricc13
=
(-2*Gam113*GamDDD111*invgam11 + 2*Gam111*GamDDD113*invgam11 -
    2*Gam213*GamDDD211*invgam11 + 2*Gam211*GamDDD213*invgam11 -
    2*Gam313*GamDDD311*invgam11 + 2*Gam311*GamDDD313*invgam11 -
    4*Gam113*GamDDD112*invgam12 + 2*Gam112*GamDDD113*invgam12 +
    2*Gam111*GamDDD123*invgam12 - 4*Gam213*GamDDD212*invgam12 +
    2*Gam212*GamDDD213*invgam12 + 2*Gam211*GamDDD223*invgam12 -
    4*Gam313*GamDDD312*invgam12 + 2*Gam312*GamDDD313*invgam12 +
    2*Gam311*GamDDD323*invgam12 - 2*Gam113*GamDDD113*invgam13 +
    2*Gam111*GamDDD133*invgam13 - 2*Gam213*GamDDD213*invgam13 +
    2*Gam211*GamDDD233*invgam13 - 2*Gam313*GamDDD313*invgam13 +
    2*Gam311*GamDDD333*invgam13 - 2*Gam113*GamDDD122*invgam22 +
    2*Gam112*GamDDD123*invgam22 - 2*Gam213*GamDDD222*invgam22 +
    2*Gam212*GamDDD223*invgam22 - 2*Gam313*GamDDD322*invgam22 +
    2*Gam312*GamDDD323*invgam22 - 2*Gam113*GamDDD123*invgam23 +
    2*Gam112*GamDDD133*invgam23 - 2*Gam213*GamDDD223*invgam23 +
    2*Gam212*GamDDD233*invgam23 - 2*Gam313*GamDDD323*invgam23 +
    2*Gam312*GamDDD333*invgam23 + invgam12*ddgam1123[ijk5] +
    invgam13*ddgam1133[ijk5] - invgam12*ddgam1213[ijk5] +
    invgam22*ddgam1223[ijk5] + invgam23*ddgam1233[ijk5] -
    invgam12*ddgam1312[ijk5] - 2*invgam13*ddgam1313[ijk5] -
    invgam22*ddgam1322[ijk5] - invgam23*ddgam1323[ijk5] -
    invgam22*ddgam2213[ijk5] + invgam12*ddgam2311[ijk5] +
    invgam22*ddgam2312[ijk5] - invgam23*ddgam2313[ijk5] +
    invgam13*ddgam3311[ijk5] + invgam23*ddgam3312[ijk5])/2.
;

const auto
Ricc22
=
-(Gam122*GamDDD111*invgam11) + Gam112*GamDDD112*invgam11 -
  Gam222*GamDDD211*invgam11 + Gam212*GamDDD212*invgam11 -
  Gam322*GamDDD311*invgam11 + Gam312*GamDDD312*invgam11 -
  Gam122*GamDDD112*invgam12 + Gam112*GamDDD122*invgam12 -
  Gam222*GamDDD212*invgam12 + Gam212*GamDDD222*invgam12 -
  Gam322*GamDDD312*invgam12 + Gam312*GamDDD322*invgam12 +
  Gam123*GamDDD112*invgam13 - 2*Gam122*GamDDD113*invgam13 +
  Gam112*GamDDD123*invgam13 + Gam223*GamDDD212*invgam13 -
  2*Gam222*GamDDD213*invgam13 + Gam212*GamDDD223*invgam13 +
  Gam323*GamDDD312*invgam13 - 2*Gam322*GamDDD313*invgam13 +
  Gam312*GamDDD323*invgam13 + Gam123*GamDDD122*invgam23 -
  Gam122*GamDDD123*invgam23 + Gam223*GamDDD222*invgam23 -
  Gam222*GamDDD223*invgam23 + Gam323*GamDDD322*invgam23 -
  Gam322*GamDDD323*invgam23 + Gam123*GamDDD123*invgam33 -
  Gam122*GamDDD133*invgam33 + Gam223*GamDDD223*invgam33 -
  Gam222*GamDDD233*invgam33 + Gam323*GamDDD323*invgam33 -
  Gam322*GamDDD333*invgam33 - (invgam11*ddgam1122[ijk5])/2. +
  invgam11*ddgam1212[ijk5] + invgam13*ddgam1223[ijk5] -
  invgam13*ddgam1322[ijk5] - (invgam11*ddgam2211[ijk5])/2. -
  invgam13*ddgam2213[ijk5] - (invgam33*ddgam2233[ijk5])/2. +
  invgam13*ddgam2312[ijk5] + invgam33*ddgam2323[ijk5] -
  (invgam33*ddgam3322[ijk5])/2.
;

const auto
Ricc23
=
(-2*Gam123*GamDDD111*invgam11 + 2*Gam112*GamDDD113*invgam11 -
    2*Gam223*GamDDD211*invgam11 + 2*Gam212*GamDDD213*invgam11 -
    2*Gam323*GamDDD311*invgam11 + 2*Gam312*GamDDD313*invgam11 -
    4*Gam123*GamDDD112*invgam12 + 2*Gam122*GamDDD113*invgam12 +
    2*Gam112*GamDDD123*invgam12 - 4*Gam223*GamDDD212*invgam12 +
    2*Gam222*GamDDD213*invgam12 + 2*Gam212*GamDDD223*invgam12 -
    4*Gam323*GamDDD312*invgam12 + 2*Gam322*GamDDD313*invgam12 +
    2*Gam312*GamDDD323*invgam12 - 2*Gam123*GamDDD113*invgam13 +
    2*Gam112*GamDDD133*invgam13 - 2*Gam223*GamDDD213*invgam13 +
    2*Gam212*GamDDD233*invgam13 - 2*Gam323*GamDDD313*invgam13 +
    2*Gam312*GamDDD333*invgam13 - 2*Gam123*GamDDD122*invgam22 +
    2*Gam122*GamDDD123*invgam22 - 2*Gam223*GamDDD222*invgam22 +
    2*Gam222*GamDDD223*invgam22 - 2*Gam323*GamDDD322*invgam22 +
    2*Gam322*GamDDD323*invgam22 - 2*Gam123*GamDDD123*invgam23 +
    2*Gam122*GamDDD133*invgam23 - 2*Gam223*GamDDD223*invgam23 +
    2*Gam222*GamDDD233*invgam23 - 2*Gam323*GamDDD323*invgam23 +
    2*Gam322*GamDDD333*invgam23 - invgam11*ddgam1123[ijk5] +
    invgam11*ddgam1213[ijk5] - invgam12*ddgam1223[ijk5] +
    invgam13*ddgam1233[ijk5] + invgam11*ddgam1312[ijk5] +
    invgam12*ddgam1322[ijk5] - invgam13*ddgam1323[ijk5] +
    invgam12*ddgam2213[ijk5] + invgam23*ddgam2233[ijk5] -
    invgam11*ddgam2311[ijk5] - invgam12*ddgam2312[ijk5] -
    invgam13*ddgam2313[ijk5] - 2*invgam23*ddgam2323[ijk5] +
    invgam13*ddgam3312[ijk5] + invgam23*ddgam3322[ijk5])/2.
;

const auto
Ricc33
=
-(Gam133*GamDDD111*invgam11) + Gam113*GamDDD113*invgam11 -
  Gam233*GamDDD211*invgam11 + Gam213*GamDDD213*invgam11 -
  Gam333*GamDDD311*invgam11 + Gam313*GamDDD313*invgam11 -
  2*Gam133*GamDDD112*invgam12 + Gam123*GamDDD113*invgam12 +
  Gam113*GamDDD123*invgam12 - 2*Gam233*GamDDD212*invgam12 +
  Gam223*GamDDD213*invgam12 + Gam213*GamDDD223*invgam12 -
  2*Gam333*GamDDD312*invgam12 + Gam323*GamDDD313*invgam12 +
  Gam313*GamDDD323*invgam12 - Gam133*GamDDD113*invgam13 +
  Gam113*GamDDD133*invgam13 - Gam233*GamDDD213*invgam13 +
  Gam213*GamDDD233*invgam13 - Gam333*GamDDD313*invgam13 +
  Gam313*GamDDD333*invgam13 - Gam133*GamDDD122*invgam22 +
  Gam123*GamDDD123*invgam22 - Gam233*GamDDD222*invgam22 +
  Gam223*GamDDD223*invgam22 - Gam333*GamDDD322*invgam22 +
  Gam323*GamDDD323*invgam22 - Gam133*GamDDD123*invgam23 +
  Gam123*GamDDD133*invgam23 - Gam233*GamDDD223*invgam23 +
  Gam223*GamDDD233*invgam23 - Gam333*GamDDD323*invgam23 +
  Gam323*GamDDD333*invgam23 - (invgam11*ddgam1133[ijk5])/2. -
  invgam12*ddgam1233[ijk5] + invgam11*ddgam1313[ijk5] +
  invgam12*ddgam1323[ijk5] - (invgam22*ddgam2233[ijk5])/2. +
  invgam12*ddgam2313[ijk5] + invgam22*ddgam2323[ijk5] -
  (invgam11*ddgam3311[ijk5])/2. - invgam12*ddgam3312[ijk5] -
  (invgam22*ddgam3322[ijk5])/2.
;

const auto
trexK
=
invgam11*exK11[ijk] + 2*invgam12*exK12[ijk] + 2*invgam13*exK13[ijk] +
  invgam22*exK22[ijk] + 2*invgam23*exK23[ijk] + invgam33*exK33[ijk]
;

const auto
DexK111
=
dexK111[ijk5] - 2*(Gam111*exK11[ijk] + Gam211*exK12[ijk] + Gam311*exK13[ijk])
;

const auto
DexK112
=
dexK112[ijk5] - Gam112*exK11[ijk] - Gam111*exK12[ijk] - Gam212*exK12[ijk] -
  Gam312*exK13[ijk] - Gam211*exK22[ijk] - Gam311*exK23[ijk]
;

const auto
DexK113
=
dexK113[ijk5] - Gam113*exK11[ijk] - Gam213*exK12[ijk] - Gam111*exK13[ijk] -
  Gam313*exK13[ijk] - Gam211*exK23[ijk] - Gam311*exK33[ijk]
;

const auto
DexK122
=
dexK122[ijk5] - 2*(Gam112*exK12[ijk] + Gam212*exK22[ijk] + Gam312*exK23[ijk])
;

const auto
DexK123
=
dexK123[ijk5] - Gam113*exK12[ijk] - Gam112*exK13[ijk] - Gam213*exK22[ijk] -
  Gam212*exK23[ijk] - Gam313*exK23[ijk] - Gam312*exK33[ijk]
;

const auto
DexK133
=
dexK133[ijk5] - 2*(Gam113*exK13[ijk] + Gam213*exK23[ijk] + Gam313*exK33[ijk])
;

const auto
DexK211
=
dexK211[ijk5] - 2*(Gam112*exK11[ijk] + Gam212*exK12[ijk] + Gam312*exK13[ijk])
;

const auto
DexK212
=
dexK212[ijk5] - Gam122*exK11[ijk] - Gam112*exK12[ijk] - Gam222*exK12[ijk] -
  Gam322*exK13[ijk] - Gam212*exK22[ijk] - Gam312*exK23[ijk]
;

const auto
DexK213
=
dexK213[ijk5] - Gam123*exK11[ijk] - Gam223*exK12[ijk] - Gam112*exK13[ijk] -
  Gam323*exK13[ijk] - Gam212*exK23[ijk] - Gam312*exK33[ijk]
;

const auto
DexK222
=
dexK222[ijk5] - 2*(Gam122*exK12[ijk] + Gam222*exK22[ijk] + Gam322*exK23[ijk])
;

const auto
DexK223
=
dexK223[ijk5] - Gam123*exK12[ijk] - Gam122*exK13[ijk] - Gam223*exK22[ijk] -
  Gam222*exK23[ijk] - Gam323*exK23[ijk] - Gam322*exK33[ijk]
;

const auto
DexK233
=
dexK233[ijk5] - 2*(Gam123*exK13[ijk] + Gam223*exK23[ijk] + Gam323*exK33[ijk])
;

const auto
DexK311
=
dexK311[ijk5] - 2*(Gam113*exK11[ijk] + Gam213*exK12[ijk] + Gam313*exK13[ijk])
;

const auto
DexK312
=
dexK312[ijk5] - Gam123*exK11[ijk] - Gam113*exK12[ijk] - Gam223*exK12[ijk] -
  Gam323*exK13[ijk] - Gam213*exK22[ijk] - Gam313*exK23[ijk]
;

const auto
DexK313
=
dexK313[ijk5] - Gam133*exK11[ijk] - Gam233*exK12[ijk] - Gam113*exK13[ijk] -
  Gam333*exK13[ijk] - Gam213*exK23[ijk] - Gam313*exK33[ijk]
;

const auto
DexK322
=
dexK322[ijk5] - 2*(Gam123*exK12[ijk] + Gam223*exK22[ijk] + Gam323*exK23[ijk])
;

const auto
DexK323
=
dexK323[ijk5] - Gam133*exK12[ijk] - Gam123*exK13[ijk] - Gam233*exK22[ijk] -
  Gam223*exK23[ijk] - Gam333*exK23[ijk] - Gam323*exK33[ijk]
;

const auto
DexK333
=
dexK333[ijk5] - 2*(Gam133*exK13[ijk] + Gam233*exK23[ijk] + Gam333*exK33[ijk])
;


const auto
Rsca
=
invgam11*Ricc11 + 2*invgam12*Ricc12 + 2*invgam13*Ricc13 + invgam22*Ricc22 +
  2*invgam23*Ricc23 + invgam33*Ricc33
;

const auto
Riem1212
=
-2*Ricc12*gam12[ijk] + (Rsca*Power(gam12[ijk],2))/2. + Ricc11*gam22[ijk] +
  gam11[ijk]*(Ricc22 - (Rsca*gam22[ijk])/2.)
;

const auto
Riem1213
=
-(Ricc12*gam13[ijk]) + gam12[ijk]*(-Ricc13 + (Rsca*gam13[ijk])/2.) +
  Ricc11*gam23[ijk] + gam11[ijk]*(Ricc23 - (Rsca*gam23[ijk])/2.)
;

const auto
Riem1223
=
-(Ricc13*gam22[ijk]) + gam13[ijk]*(-Ricc22 + (Rsca*gam22[ijk])/2.) +
  Ricc12*gam23[ijk] + gam12[ijk]*(Ricc23 - (Rsca*gam23[ijk])/2.)
;

const auto
Riem1312
=
-(Ricc12*gam13[ijk]) + gam12[ijk]*(-Ricc13 + (Rsca*gam13[ijk])/2.) +
  Ricc11*gam23[ijk] + gam11[ijk]*(Ricc23 - (Rsca*gam23[ijk])/2.)
;

const auto
Riem1313
=
-2*Ricc13*gam13[ijk] + (Rsca*Power(gam13[ijk],2))/2. + Ricc11*gam33[ijk] +
  gam11[ijk]*(Ricc33 - (Rsca*gam33[ijk])/2.)
;

const auto
Riem1323
=
-(Ricc13*gam23[ijk]) + gam13[ijk]*(-Ricc23 + (Rsca*gam23[ijk])/2.) +
  Ricc12*gam33[ijk] + gam12[ijk]*(Ricc33 - (Rsca*gam33[ijk])/2.)
;

const auto
Riem2312
=
-(Ricc13*gam22[ijk]) + gam13[ijk]*(-Ricc22 + (Rsca*gam22[ijk])/2.) +
  Ricc12*gam23[ijk] + gam12[ijk]*(Ricc23 - (Rsca*gam23[ijk])/2.)
;

const auto
Riem2313
=
-(Ricc13*gam23[ijk]) + gam13[ijk]*(-Ricc23 + (Rsca*gam23[ijk])/2.) +
  Ricc12*gam33[ijk] + gam12[ijk]*(Ricc33 - (Rsca*gam33[ijk])/2.)
;

const auto
Riem2323
=
-2*Ricc23*gam23[ijk] + (Rsca*Power(gam23[ijk],2))/2. + Ricc22*gam33[ijk] +
  gam22[ijk]*(Ricc33 - (Rsca*gam33[ijk])/2.)
;

const auto
perpR1212
=
Riem1212 - Power(exK12[ijk],2) + exK11[ijk]*exK22[ijk]
;

const auto
perpR1213
=
Riem1213 - exK12[ijk]*exK13[ijk] + exK11[ijk]*exK23[ijk]
;

const auto
perpR1223
=
Riem1223 - exK13[ijk]*exK22[ijk] + exK12[ijk]*exK23[ijk]
;

const auto
perpR1312
=
Riem1312 - exK12[ijk]*exK13[ijk] + exK11[ijk]*exK23[ijk]
;

const auto
perpR1313
=
Riem1313 - Power(exK13[ijk],2) + exK11[ijk]*exK33[ijk]
;

const auto
perpR1323
=
Riem1323 - exK13[ijk]*exK23[ijk] + exK12[ijk]*exK33[ijk]
;

const auto
perpR2312
=
Riem2312 - exK13[ijk]*exK22[ijk] + exK12[ijk]*exK23[ijk]
;

const auto
perpR2313
=
Riem2313 - exK13[ijk]*exK23[ijk] + exK12[ijk]*exK33[ijk]
;

const auto
perpR2323
=
Riem2323 - Power(exK23[ijk],2) + exK22[ijk]*exK33[ijk]
;

const auto
perpRn112
=
DexK112 - DexK211
;

const auto
perpRn113
=
DexK113 - DexK311
;

const auto
perpRn123
=
DexK213 - DexK312
;

const auto
perpRn212
=
DexK122 - DexK212
;

const auto
perpRn213
=
DexK123 - DexK312
;

const auto
perpRn223
=
DexK223 - DexK322
;

const auto
perpRn312
=
DexK123 - DexK213
;

const auto
perpRn313
=
DexK133 - DexK313
;

const auto
perpRn323
=
DexK233 - DexK323
;

const auto
perpRnn11
=
Ricc11 - invgam11*Power(exK11[ijk],2) - invgam22*Power(exK12[ijk],2) -
  2*invgam23*exK12[ijk]*exK13[ijk] - invgam33*Power(exK13[ijk],2) +
  exK11[ijk]*(trexK - 2*invgam12*exK12[ijk] - 2*invgam13*exK13[ijk])
;

const auto
perpRnn12
=
Ricc12 - invgam12*Power(exK12[ijk],2) - invgam23*exK13[ijk]*exK22[ijk] -
  invgam33*exK13[ijk]*exK23[ijk] -
  exK11[ijk]*(invgam12*exK22[ijk] + invgam13*exK23[ijk]) +
  exK12[ijk]*(trexK - invgam11*exK11[ijk] - invgam13*exK13[ijk] -
     invgam22*exK22[ijk] - invgam23*exK23[ijk])
;

const auto
perpRnn13
=
Ricc13 - invgam13*Power(exK13[ijk],2) - invgam22*exK12[ijk]*exK23[ijk] -
  invgam23*exK12[ijk]*exK33[ijk] -
  exK11[ijk]*(invgam12*exK23[ijk] + invgam13*exK33[ijk]) +
  exK13[ijk]*(trexK - invgam11*exK11[ijk] - invgam12*exK12[ijk] -
     invgam23*exK23[ijk] - invgam33*exK33[ijk])
;

const auto
perpRnn22
=
Ricc22 - invgam11*Power(exK12[ijk],2) - invgam22*Power(exK22[ijk],2) -
  invgam33*Power(exK23[ijk],2) -
  2*exK12[ijk]*(invgam12*exK22[ijk] + invgam13*exK23[ijk]) +
  exK22[ijk]*(trexK - 2*invgam23*exK23[ijk])
;

const auto
perpRnn23
=
Ricc23 + trexK*exK23[ijk] - invgam22*exK22[ijk]*exK23[ijk] -
  invgam23*Power(exK23[ijk],2) -
  exK13[ijk]*(invgam12*exK22[ijk] + invgam13*exK23[ijk]) -
  invgam23*exK22[ijk]*exK33[ijk] - invgam33*exK23[ijk]*exK33[ijk] -
  exK12[ijk]*(invgam11*exK13[ijk] + invgam12*exK23[ijk] +
     invgam13*exK33[ijk])
;

const auto
perpRnn33
=
Ricc33 - invgam11*Power(exK13[ijk],2) - invgam22*Power(exK23[ijk],2) +
  trexK*exK33[ijk] - 2*invgam23*exK23[ijk]*exK33[ijk] -
  invgam33*Power(exK33[ijk],2) -
  2*exK13[ijk]*(invgam12*exK23[ijk] + invgam13*exK33[ijk])
;


const auto
uvec1
=
-coY
;

const auto
uvec2
=
coX
;

const auto
uvec3
=
0
;

const auto
vvec1
=
coX
;

const auto
vvec2
=
coY
;

const auto
vvec3
=
coZ
;

const auto
wvec1
=
sqrtdetgam*(-(invgam13*uvec2*vvec1) + invgam12*uvec3*vvec1 +
    invgam13*uvec1*vvec2 - invgam11*uvec3*vvec2 - invgam12*uvec1*vvec3 +
    invgam11*uvec2*vvec3)
;

const auto
wvec2
=
sqrtdetgam*(-(invgam23*uvec2*vvec1) + invgam22*uvec3*vvec1 +
    invgam23*uvec1*vvec2 - invgam12*uvec3*vvec2 - invgam22*uvec1*vvec3 +
    invgam12*uvec2*vvec3)
;

const auto
wvec3
=
sqrtdetgam*(-(invgam33*uvec2*vvec1) + invgam23*uvec3*vvec1 +
    invgam33*uvec1*vvec2 - invgam13*uvec3*vvec2 - invgam23*uvec1*vvec3 +
    invgam13*uvec2*vvec3)
;

const auto
inneruutmp
=
Power(uvec1,2)*gam11[ijk] + 2*uvec1*uvec2*gam12[ijk] +
  2*uvec1*uvec3*gam13[ijk] + Power(uvec2,2)*gam22[ijk] +
  2*uvec2*uvec3*gam23[ijk] + Power(uvec3,2)*gam33[ijk]
;

const auto
inneruu
=
Max(inneruutmp,ncutoff)
;

const auto
Uvec1
=
uvec1/Sqrt(inneruu)
;

const auto
Uvec2
=
uvec2/Sqrt(inneruu)
;

const auto
Uvec3
=
uvec3/Sqrt(inneruu)
;

const auto
innerUv
=
Uvec1*vvec1*gam11[ijk] + (Uvec2*vvec1 + Uvec1*vvec2)*gam12[ijk] +
  Uvec3*vvec1*gam13[ijk] + Uvec1*vvec3*gam13[ijk] + Uvec2*vvec2*gam22[ijk] +
  Uvec3*vvec2*gam23[ijk] + Uvec2*vvec3*gam23[ijk] + Uvec3*vvec3*gam33[ijk]
;

const auto
Vtmp1
=
-(innerUv*Uvec1) + vvec1
;

const auto
Vtmp2
=
-(innerUv*Uvec2) + vvec2
;

const auto
Vtmp3
=
-(innerUv*Uvec3) + vvec3
;

const auto
innerVVtmp
=
Power(Vtmp1,2)*gam11[ijk] + 2*Vtmp1*Vtmp2*gam12[ijk] +
  2*Vtmp1*Vtmp3*gam13[ijk] + Power(Vtmp2,2)*gam22[ijk] +
  2*Vtmp2*Vtmp3*gam23[ijk] + Power(Vtmp3,2)*gam33[ijk]
;

const auto
innerVV
=
Max(innerVVtmp,ncutoff)
;

const auto
Vvec1
=
Vtmp1/Sqrt(innerVV)
;

const auto
Vvec2
=
Vtmp2/Sqrt(innerVV)
;

const auto
Vvec3
=
Vtmp3/Sqrt(innerVV)
;

const auto
innerUw
=
Uvec1*wvec1*gam11[ijk] + (Uvec2*wvec1 + Uvec1*wvec2)*gam12[ijk] +
  Uvec3*wvec1*gam13[ijk] + Uvec1*wvec3*gam13[ijk] + Uvec2*wvec2*gam22[ijk] +
  Uvec3*wvec2*gam23[ijk] + Uvec2*wvec3*gam23[ijk] + Uvec3*wvec3*gam33[ijk]
;

const auto
innerVw
=
Vvec1*wvec1*gam11[ijk] + (Vvec2*wvec1 + Vvec1*wvec2)*gam12[ijk] +
  Vvec3*wvec1*gam13[ijk] + Vvec1*wvec3*gam13[ijk] + Vvec2*wvec2*gam22[ijk] +
  Vvec3*wvec2*gam23[ijk] + Vvec2*wvec3*gam23[ijk] + Vvec3*wvec3*gam33[ijk]
;

const auto
Wtmp1
=
-(innerUw*Uvec1) - innerVw*Vvec1 + wvec1
;

const auto
Wtmp2
=
-(innerUw*Uvec2) - innerVw*Vvec2 + wvec2
;

const auto
Wtmp3
=
-(innerUw*Uvec3) - innerVw*Vvec3 + wvec3
;

const auto
innerWWtmp
=
Power(Wtmp1,2)*gam11[ijk] + 2*Wtmp1*Wtmp2*gam12[ijk] +
  2*Wtmp1*Wtmp3*gam13[ijk] + Power(Wtmp2,2)*gam22[ijk] +
  2*Wtmp2*Wtmp3*gam23[ijk] + Power(Wtmp3,2)*gam33[ijk]
;

const auto
innerWW
=
Max(innerWWtmp,ncutoff)
;

const auto
Wvec1
=
Wtmp1/Sqrt(innerWW)
;

const auto
Wvec2
=
Wtmp2/Sqrt(innerWW)
;

const auto
Wvec3
=
Wtmp3/Sqrt(innerWW)
;

const auto
mbmbreal11
=
(-Power(Uvec1,2) + Power(Wvec1,2))/2.
;

const auto
mbmbreal12
=
(-(Uvec1*Uvec2) + Wvec1*Wvec2)/2.
;

const auto
mbmbreal13
=
(-(Uvec1*Uvec3) + Wvec1*Wvec3)/2.
;

const auto
mbmbreal22
=
(-Power(Uvec2,2) + Power(Wvec2,2))/2.
;

const auto
mbmbreal23
=
(-(Uvec2*Uvec3) + Wvec2*Wvec3)/2.
;

const auto
mbmbreal33
=
(-Power(Uvec3,2) + Power(Wvec3,2))/2.
;

const auto
mbmbimag11
=
-(Uvec1*Wvec1)
;

const auto
mbmbimag12
=
(-(Uvec2*Wvec1) - Uvec1*Wvec2)/2.
;

const auto
mbmbimag13
=
(-(Uvec3*Wvec1) - Uvec1*Wvec3)/2.
;

const auto
mbmbimag22
=
-(Uvec2*Wvec2)
;

const auto
mbmbimag23
=
(-(Uvec3*Wvec2) - Uvec2*Wvec3)/2.
;

const auto
mbmbimag33
=
-(Uvec3*Wvec3)
;


Psi4real[ijk]
=
(2*mbmbreal13*perpRnn13 + mbmbreal22*perpRnn22 + 2*mbmbreal23*perpRnn23 +
    mbmbreal33*perpRnn33 - 2*mbmbreal13*perpRn113*Vvec1 -
    2*mbmbreal22*perpRn212*Vvec1 - 2*mbmbreal23*perpRn213*Vvec1 -
    2*mbmbreal23*perpRn312*Vvec1 - 2*mbmbreal33*perpRn313*Vvec1 +
    mbmbreal22*perpR1212*Power(Vvec1,2) +
    mbmbreal23*perpR1213*Power(Vvec1,2) +
    mbmbreal23*perpR1312*Power(Vvec1,2) +
    mbmbreal33*perpR1313*Power(Vvec1,2) - 2*mbmbreal13*perpRn123*Vvec2 -
    2*mbmbreal23*perpRn223*Vvec2 + 2*mbmbreal13*perpRn312*Vvec2 -
    2*mbmbreal33*perpRn323*Vvec2 - mbmbreal13*perpR1213*Vvec1*Vvec2 +
    mbmbreal23*perpR1223*Vvec1*Vvec2 - mbmbreal13*perpR1312*Vvec1*Vvec2 +
    mbmbreal33*perpR1323*Vvec1*Vvec2 + mbmbreal23*perpR2312*Vvec1*Vvec2 +
    mbmbreal33*perpR2313*Vvec1*Vvec2 -
    mbmbreal13*perpR1223*Power(Vvec2,2) -
    mbmbreal13*perpR2312*Power(Vvec2,2) +
    mbmbreal33*perpR2323*Power(Vvec2,2) + 2*mbmbreal22*perpRn223*Vvec3 +
    2*mbmbreal13*perpRn313*Vvec3 + 2*mbmbreal23*perpRn323*Vvec3 -
    mbmbreal22*perpR1223*Vvec1*Vvec3 - 2*mbmbreal13*perpR1313*Vvec1*Vvec3 -
    mbmbreal23*perpR1323*Vvec1*Vvec3 - mbmbreal22*perpR2312*Vvec1*Vvec3 -
    mbmbreal23*perpR2313*Vvec1*Vvec3 - mbmbreal13*perpR1323*Vvec2*Vvec3 -
    mbmbreal13*perpR2313*Vvec2*Vvec3 - 2*mbmbreal23*perpR2323*Vvec2*Vvec3 +
    mbmbreal22*perpR2323*Power(Vvec3,2) +
    mbmbreal11*(perpRnn11 + 2*perpRn112*Vvec2 + perpR1212*Power(Vvec2,2) +
       2*perpRn113*Vvec3 + perpR1213*Vvec2*Vvec3 + perpR1312*Vvec2*Vvec3 +
       perpR1313*Power(Vvec3,2)) +
    mbmbreal12*(2*perpRnn12 - 2*perpRn112*Vvec1 + 2*perpRn212*Vvec2 -
       2*perpR1212*Vvec1*Vvec2 + 2*perpRn123*Vvec3 + 2*perpRn213*Vvec3 -
       perpR1213*Vvec1*Vvec3 - perpR1312*Vvec1*Vvec3 +
       perpR1223*Vvec2*Vvec3 + perpR2312*Vvec2*Vvec3 +
       perpR1323*Power(Vvec3,2) + perpR2313*Power(Vvec3,2)))/2.
;

Psi4imag[ijk]
=
(2*mbmbimag13*perpRnn13 + mbmbimag22*perpRnn22 + 2*mbmbimag23*perpRnn23 +
    mbmbimag33*perpRnn33 - 2*mbmbimag13*perpRn113*Vvec1 -
    2*mbmbimag22*perpRn212*Vvec1 - 2*mbmbimag23*perpRn213*Vvec1 -
    2*mbmbimag23*perpRn312*Vvec1 - 2*mbmbimag33*perpRn313*Vvec1 +
    mbmbimag22*perpR1212*Power(Vvec1,2) +
    mbmbimag23*perpR1213*Power(Vvec1,2) +
    mbmbimag23*perpR1312*Power(Vvec1,2) +
    mbmbimag33*perpR1313*Power(Vvec1,2) - 2*mbmbimag13*perpRn123*Vvec2 -
    2*mbmbimag23*perpRn223*Vvec2 + 2*mbmbimag13*perpRn312*Vvec2 -
    2*mbmbimag33*perpRn323*Vvec2 - mbmbimag13*perpR1213*Vvec1*Vvec2 +
    mbmbimag23*perpR1223*Vvec1*Vvec2 - mbmbimag13*perpR1312*Vvec1*Vvec2 +
    mbmbimag33*perpR1323*Vvec1*Vvec2 + mbmbimag23*perpR2312*Vvec1*Vvec2 +
    mbmbimag33*perpR2313*Vvec1*Vvec2 -
    mbmbimag13*perpR1223*Power(Vvec2,2) -
    mbmbimag13*perpR2312*Power(Vvec2,2) +
    mbmbimag33*perpR2323*Power(Vvec2,2) + 2*mbmbimag22*perpRn223*Vvec3 +
    2*mbmbimag13*perpRn313*Vvec3 + 2*mbmbimag23*perpRn323*Vvec3 -
    mbmbimag22*perpR1223*Vvec1*Vvec3 - 2*mbmbimag13*perpR1313*Vvec1*Vvec3 -
    mbmbimag23*perpR1323*Vvec1*Vvec3 - mbmbimag22*perpR2312*Vvec1*Vvec3 -
    mbmbimag23*perpR2313*Vvec1*Vvec3 - mbmbimag13*perpR1323*Vvec2*Vvec3 -
    mbmbimag13*perpR2313*Vvec2*Vvec3 - 2*mbmbimag23*perpR2323*Vvec2*Vvec3 +
    mbmbimag22*perpR2323*Power(Vvec3,2) +
    mbmbimag11*(perpRnn11 + 2*perpRn112*Vvec2 + perpR1212*Power(Vvec2,2) +
       2*perpRn113*Vvec3 + perpR1213*Vvec2*Vvec3 + perpR1312*Vvec2*Vvec3 +
       perpR1313*Power(Vvec3,2)) +
    mbmbimag12*(2*perpRnn12 - 2*perpRn112*Vvec1 + 2*perpRn212*Vvec2 -
       2*perpR1212*Vvec1*Vvec2 + 2*perpRn123*Vvec3 + 2*perpRn213*Vvec3 -
       perpR1213*Vvec1*Vvec3 - perpR1312*Vvec1*Vvec3 +
       perpR1223*Vvec2*Vvec3 + perpR2312*Vvec2*Vvec3 +
       perpR1323*Power(Vvec3,2) + perpR2313*Power(Vvec3,2)))/2.
;


  });
});

#endif // #ifndef WEYLSCALO_SET_PSI4_GAUSSCODAZZI_GF3D5_HXX

/* WeylScalo_set_Psi4_GaussCodazzi_GF3D5.hxx */
