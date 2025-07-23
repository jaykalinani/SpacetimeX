/* Z4cowGPU_set_rhs_GF3D5.hxx */
/* Produced with Generato */

#ifndef Z4COWGPU_SET_RHS_GF3D5_HXX
#define Z4COWGPU_SET_RHS_GF3D5_HXX

const auto dtW = gf_dtW;
const auto dtgamt11 = gf_dtgamt[0];
const auto dtgamt12 = gf_dtgamt[1];
const auto dtgamt13 = gf_dtgamt[2];
const auto dtgamt22 = gf_dtgamt[3];
const auto dtgamt23 = gf_dtgamt[4];
const auto dtgamt33 = gf_dtgamt[5];
const auto dtexKh = gf_dtexKh;
const auto dtexAt11 = gf_dtexAt[0];
const auto dtexAt12 = gf_dtexAt[1];
const auto dtexAt13 = gf_dtexAt[2];
const auto dtexAt22 = gf_dtexAt[3];
const auto dtexAt23 = gf_dtexAt[4];
const auto dtexAt33 = gf_dtexAt[5];
const auto dttrGt1 = gf_dttrGt[0];
const auto dttrGt2 = gf_dttrGt[1];
const auto dttrGt3 = gf_dttrGt[2];
const auto dtTheta = gf_dtTheta;
const auto dtalpha = gf_dtalpha;
const auto dtbeta1 = gf_dtbeta[0];
const auto dtbeta2 = gf_dtbeta[1];
const auto dtbeta3 = gf_dtbeta[2];
const auto eTt1 = gf_eTt[0];
const auto eTt2 = gf_eTt[1];
const auto eTt3 = gf_eTt[2];
const auto eT11 = gf_eT[0];
const auto eT12 = gf_eT[1];
const auto eT13 = gf_eT[2];
const auto eT22 = gf_eT[3];
const auto eT23 = gf_eT[4];
const auto eT33 = gf_eT[5];
const auto gamt11 = gf_gamt[0];
const auto gamt12 = gf_gamt[1];
const auto gamt13 = gf_gamt[2];
const auto gamt22 = gf_gamt[3];
const auto gamt23 = gf_gamt[4];
const auto gamt33 = gf_gamt[5];
const auto exKh = gf_exKh;
const auto exAt11 = gf_exAt[0];
const auto exAt12 = gf_exAt[1];
const auto exAt13 = gf_exAt[2];
const auto exAt22 = gf_exAt[3];
const auto exAt23 = gf_exAt[4];
const auto exAt33 = gf_exAt[5];
const auto trGt1 = gf_trGt[0];
const auto trGt2 = gf_trGt[1];
const auto trGt3 = gf_trGt[2];
const auto alpha = gf_alpha;
const auto beta1 = gf_beta[0];
const auto beta2 = gf_beta[1];
const auto beta3 = gf_beta[2];

const auto dW1 = tl_dW[0].ptr;
const auto dW2 = tl_dW[1].ptr;
const auto dW3 = tl_dW[2].ptr;
const auto dgamt111 = tl_dgamt[0][0].ptr;
const auto dgamt112 = tl_dgamt[1][0].ptr;
const auto dgamt113 = tl_dgamt[2][0].ptr;
const auto dgamt122 = tl_dgamt[3][0].ptr;
const auto dgamt123 = tl_dgamt[4][0].ptr;
const auto dgamt133 = tl_dgamt[5][0].ptr;
const auto dgamt211 = tl_dgamt[0][1].ptr;
const auto dgamt212 = tl_dgamt[1][1].ptr;
const auto dgamt213 = tl_dgamt[2][1].ptr;
const auto dgamt222 = tl_dgamt[3][1].ptr;
const auto dgamt223 = tl_dgamt[4][1].ptr;
const auto dgamt233 = tl_dgamt[5][1].ptr;
const auto dgamt311 = tl_dgamt[0][2].ptr;
const auto dgamt312 = tl_dgamt[1][2].ptr;
const auto dgamt313 = tl_dgamt[2][2].ptr;
const auto dgamt322 = tl_dgamt[3][2].ptr;
const auto dgamt323 = tl_dgamt[4][2].ptr;
const auto dgamt333 = tl_dgamt[5][2].ptr;
const auto dexKh1 = tl_dexKh[0].ptr;
const auto dexKh2 = tl_dexKh[1].ptr;
const auto dexKh3 = tl_dexKh[2].ptr;
const auto dexAt111 = tl_dexAt[0][0].ptr;
const auto dexAt112 = tl_dexAt[1][0].ptr;
const auto dexAt113 = tl_dexAt[2][0].ptr;
const auto dexAt122 = tl_dexAt[3][0].ptr;
const auto dexAt123 = tl_dexAt[4][0].ptr;
const auto dexAt133 = tl_dexAt[5][0].ptr;
const auto dexAt211 = tl_dexAt[0][1].ptr;
const auto dexAt212 = tl_dexAt[1][1].ptr;
const auto dexAt213 = tl_dexAt[2][1].ptr;
const auto dexAt222 = tl_dexAt[3][1].ptr;
const auto dexAt223 = tl_dexAt[4][1].ptr;
const auto dexAt233 = tl_dexAt[5][1].ptr;
const auto dexAt311 = tl_dexAt[0][2].ptr;
const auto dexAt312 = tl_dexAt[1][2].ptr;
const auto dexAt313 = tl_dexAt[2][2].ptr;
const auto dexAt322 = tl_dexAt[3][2].ptr;
const auto dexAt323 = tl_dexAt[4][2].ptr;
const auto dexAt333 = tl_dexAt[5][2].ptr;
const auto dtrGt11 = tl_dtrGt[0][0].ptr;
const auto dtrGt12 = tl_dtrGt[1][0].ptr;
const auto dtrGt13 = tl_dtrGt[2][0].ptr;
const auto dtrGt21 = tl_dtrGt[0][1].ptr;
const auto dtrGt22 = tl_dtrGt[1][1].ptr;
const auto dtrGt23 = tl_dtrGt[2][1].ptr;
const auto dtrGt31 = tl_dtrGt[0][2].ptr;
const auto dtrGt32 = tl_dtrGt[1][2].ptr;
const auto dtrGt33 = tl_dtrGt[2][2].ptr;
const auto dTheta1 = tl_dTheta[0].ptr;
const auto dTheta2 = tl_dTheta[1].ptr;
const auto dTheta3 = tl_dTheta[2].ptr;
const auto dalpha1 = tl_dalpha[0].ptr;
const auto dalpha2 = tl_dalpha[1].ptr;
const auto dalpha3 = tl_dalpha[2].ptr;
const auto dbeta11 = tl_dbeta[0][0].ptr;
const auto dbeta12 = tl_dbeta[1][0].ptr;
const auto dbeta13 = tl_dbeta[2][0].ptr;
const auto dbeta21 = tl_dbeta[0][1].ptr;
const auto dbeta22 = tl_dbeta[1][1].ptr;
const auto dbeta23 = tl_dbeta[2][1].ptr;
const auto dbeta31 = tl_dbeta[0][2].ptr;
const auto dbeta32 = tl_dbeta[1][2].ptr;
const auto dbeta33 = tl_dbeta[2][2].ptr;
const auto ddW11 = tl_ddW[0].ptr;
const auto ddW12 = tl_ddW[1].ptr;
const auto ddW13 = tl_ddW[2].ptr;
const auto ddW22 = tl_ddW[3].ptr;
const auto ddW23 = tl_ddW[4].ptr;
const auto ddW33 = tl_ddW[5].ptr;
const auto ddgamt1111 = tl_ddgamt[0][0].ptr;
const auto ddgamt1112 = tl_ddgamt[1][0].ptr;
const auto ddgamt1113 = tl_ddgamt[2][0].ptr;
const auto ddgamt1122 = tl_ddgamt[3][0].ptr;
const auto ddgamt1123 = tl_ddgamt[4][0].ptr;
const auto ddgamt1133 = tl_ddgamt[5][0].ptr;
const auto ddgamt1211 = tl_ddgamt[0][1].ptr;
const auto ddgamt1212 = tl_ddgamt[1][1].ptr;
const auto ddgamt1213 = tl_ddgamt[2][1].ptr;
const auto ddgamt1222 = tl_ddgamt[3][1].ptr;
const auto ddgamt1223 = tl_ddgamt[4][1].ptr;
const auto ddgamt1233 = tl_ddgamt[5][1].ptr;
const auto ddgamt1311 = tl_ddgamt[0][2].ptr;
const auto ddgamt1312 = tl_ddgamt[1][2].ptr;
const auto ddgamt1313 = tl_ddgamt[2][2].ptr;
const auto ddgamt1322 = tl_ddgamt[3][2].ptr;
const auto ddgamt1323 = tl_ddgamt[4][2].ptr;
const auto ddgamt1333 = tl_ddgamt[5][2].ptr;
const auto ddgamt2211 = tl_ddgamt[0][3].ptr;
const auto ddgamt2212 = tl_ddgamt[1][3].ptr;
const auto ddgamt2213 = tl_ddgamt[2][3].ptr;
const auto ddgamt2222 = tl_ddgamt[3][3].ptr;
const auto ddgamt2223 = tl_ddgamt[4][3].ptr;
const auto ddgamt2233 = tl_ddgamt[5][3].ptr;
const auto ddgamt2311 = tl_ddgamt[0][4].ptr;
const auto ddgamt2312 = tl_ddgamt[1][4].ptr;
const auto ddgamt2313 = tl_ddgamt[2][4].ptr;
const auto ddgamt2322 = tl_ddgamt[3][4].ptr;
const auto ddgamt2323 = tl_ddgamt[4][4].ptr;
const auto ddgamt2333 = tl_ddgamt[5][4].ptr;
const auto ddgamt3311 = tl_ddgamt[0][5].ptr;
const auto ddgamt3312 = tl_ddgamt[1][5].ptr;
const auto ddgamt3313 = tl_ddgamt[2][5].ptr;
const auto ddgamt3322 = tl_ddgamt[3][5].ptr;
const auto ddgamt3323 = tl_ddgamt[4][5].ptr;
const auto ddgamt3333 = tl_ddgamt[5][5].ptr;
const auto ddalpha11 = tl_ddalpha[0].ptr;
const auto ddalpha12 = tl_ddalpha[1].ptr;
const auto ddalpha13 = tl_ddalpha[2].ptr;
const auto ddalpha22 = tl_ddalpha[3].ptr;
const auto ddalpha23 = tl_ddalpha[4].ptr;
const auto ddalpha33 = tl_ddalpha[5].ptr;
const auto ddbeta111 = tl_ddbeta[0][0].ptr;
const auto ddbeta112 = tl_ddbeta[1][0].ptr;
const auto ddbeta113 = tl_ddbeta[2][0].ptr;
const auto ddbeta121 = tl_ddbeta[0][1].ptr;
const auto ddbeta122 = tl_ddbeta[1][1].ptr;
const auto ddbeta123 = tl_ddbeta[2][1].ptr;
const auto ddbeta131 = tl_ddbeta[0][2].ptr;
const auto ddbeta132 = tl_ddbeta[1][2].ptr;
const auto ddbeta133 = tl_ddbeta[2][2].ptr;
const auto ddbeta221 = tl_ddbeta[0][3].ptr;
const auto ddbeta222 = tl_ddbeta[1][3].ptr;
const auto ddbeta223 = tl_ddbeta[2][3].ptr;
const auto ddbeta231 = tl_ddbeta[0][4].ptr;
const auto ddbeta232 = tl_ddbeta[1][4].ptr;
const auto ddbeta233 = tl_ddbeta[2][4].ptr;
const auto ddbeta331 = tl_ddbeta[0][5].ptr;
const auto ddbeta332 = tl_ddbeta[1][5].ptr;
const auto ddbeta333 = tl_ddbeta[2][5].ptr;

noinline([&]() __attribute__((__flatten__, __hot__)) {
grid.loop_int_device<0, 0, 0>(
  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
const int ijk = layout2.linear(p.i, p.j, p.k);
const int ijk5 = layout5.linear(p.i, p.j, p.k);

const CCTK_REAL ceta = calceta(p.x, p.y, p.z);

const auto
dlnW1
=
dW1[ijk5]/W[ijk]
;

const auto
dlnW2
=
dW2[ijk5]/W[ijk]
;

const auto
dlnW3
=
dW3[ijk5]/W[ijk]
;

const auto
invgamt11
=
-Power(gamt23[ijk],2) + gamt22[ijk]*gamt33[ijk]
;

const auto
invgamt12
=
gamt13[ijk]*gamt23[ijk] - gamt12[ijk]*gamt33[ijk]
;

const auto
invgamt13
=
-(gamt13[ijk]*gamt22[ijk]) + gamt12[ijk]*gamt23[ijk]
;

const auto
invgamt22
=
-Power(gamt13[ijk],2) + gamt11[ijk]*gamt33[ijk]
;

const auto
invgamt23
=
gamt12[ijk]*gamt13[ijk] - gamt11[ijk]*gamt23[ijk]
;

const auto
invgamt33
=
-Power(gamt12[ijk],2) + gamt11[ijk]*gamt22[ijk]
;

const auto
invgam11
=
invgamt11*Power(W[ijk],2)
;

const auto
invgam12
=
invgamt12*Power(W[ijk],2)
;

const auto
invgam13
=
invgamt13*Power(W[ijk],2)
;

const auto
invgam22
=
invgamt22*Power(W[ijk],2)
;

const auto
invgam23
=
invgamt23*Power(W[ijk],2)
;

const auto
invgam33
=
invgamt33*Power(W[ijk],2)
;

const auto
gam11
=
gamt11[ijk]/Power(W[ijk],2)
;

const auto
gam12
=
gamt12[ijk]/Power(W[ijk],2)
;

const auto
gam13
=
gamt13[ijk]/Power(W[ijk],2)
;

const auto
gam22
=
gamt22[ijk]/Power(W[ijk],2)
;

const auto
gam23
=
gamt23[ijk]/Power(W[ijk],2)
;

const auto
gam33
=
gamt33[ijk]/Power(W[ijk],2)
;

const auto
GtDDD111
=
dgamt111[ijk5]/2.
;

const auto
GtDDD112
=
dgamt211[ijk5]/2.
;

const auto
GtDDD113
=
dgamt311[ijk5]/2.
;

const auto
GtDDD122
=
-0.5*dgamt122[ijk5] + dgamt212[ijk5]
;

const auto
GtDDD123
=
(-dgamt123[ijk5] + dgamt213[ijk5] + dgamt312[ijk5])/2.
;

const auto
GtDDD133
=
-0.5*dgamt133[ijk5] + dgamt313[ijk5]
;

const auto
GtDDD211
=
dgamt112[ijk5] - dgamt211[ijk5]/2.
;

const auto
GtDDD212
=
dgamt122[ijk5]/2.
;

const auto
GtDDD213
=
(dgamt123[ijk5] - dgamt213[ijk5] + dgamt312[ijk5])/2.
;

const auto
GtDDD222
=
dgamt222[ijk5]/2.
;

const auto
GtDDD223
=
dgamt322[ijk5]/2.
;

const auto
GtDDD233
=
-0.5*dgamt233[ijk5] + dgamt323[ijk5]
;

const auto
GtDDD311
=
dgamt113[ijk5] - dgamt311[ijk5]/2.
;

const auto
GtDDD312
=
(dgamt123[ijk5] + dgamt213[ijk5] - dgamt312[ijk5])/2.
;

const auto
GtDDD313
=
dgamt133[ijk5]/2.
;

const auto
GtDDD322
=
dgamt223[ijk5] - dgamt322[ijk5]/2.
;

const auto
GtDDD323
=
dgamt233[ijk5]/2.
;

const auto
GtDDD333
=
dgamt333[ijk5]/2.
;

const auto
GtDDU111
=
GtDDD111*invgamt11 + GtDDD112*invgamt12 + GtDDD113*invgamt13
;

const auto
GtDDU112
=
GtDDD111*invgamt12 + GtDDD112*invgamt22 + GtDDD113*invgamt23
;

const auto
GtDDU113
=
GtDDD111*invgamt13 + GtDDD112*invgamt23 + GtDDD113*invgamt33
;

const auto
GtDDU121
=
GtDDD112*invgamt11 + GtDDD122*invgamt12 + GtDDD123*invgamt13
;

const auto
GtDDU122
=
GtDDD112*invgamt12 + GtDDD122*invgamt22 + GtDDD123*invgamt23
;

const auto
GtDDU123
=
GtDDD112*invgamt13 + GtDDD122*invgamt23 + GtDDD123*invgamt33
;

const auto
GtDDU131
=
GtDDD113*invgamt11 + GtDDD123*invgamt12 + GtDDD133*invgamt13
;

const auto
GtDDU132
=
GtDDD113*invgamt12 + GtDDD123*invgamt22 + GtDDD133*invgamt23
;

const auto
GtDDU133
=
GtDDD113*invgamt13 + GtDDD123*invgamt23 + GtDDD133*invgamt33
;

const auto
GtDDU211
=
GtDDD211*invgamt11 + GtDDD212*invgamt12 + GtDDD213*invgamt13
;

const auto
GtDDU212
=
GtDDD211*invgamt12 + GtDDD212*invgamt22 + GtDDD213*invgamt23
;

const auto
GtDDU213
=
GtDDD211*invgamt13 + GtDDD212*invgamt23 + GtDDD213*invgamt33
;

const auto
GtDDU221
=
GtDDD212*invgamt11 + GtDDD222*invgamt12 + GtDDD223*invgamt13
;

const auto
GtDDU222
=
GtDDD212*invgamt12 + GtDDD222*invgamt22 + GtDDD223*invgamt23
;

const auto
GtDDU223
=
GtDDD212*invgamt13 + GtDDD222*invgamt23 + GtDDD223*invgamt33
;

const auto
GtDDU231
=
GtDDD213*invgamt11 + GtDDD223*invgamt12 + GtDDD233*invgamt13
;

const auto
GtDDU232
=
GtDDD213*invgamt12 + GtDDD223*invgamt22 + GtDDD233*invgamt23
;

const auto
GtDDU233
=
GtDDD213*invgamt13 + GtDDD223*invgamt23 + GtDDD233*invgamt33
;

const auto
GtDDU311
=
GtDDD311*invgamt11 + GtDDD312*invgamt12 + GtDDD313*invgamt13
;

const auto
GtDDU312
=
GtDDD311*invgamt12 + GtDDD312*invgamt22 + GtDDD313*invgamt23
;

const auto
GtDDU313
=
GtDDD311*invgamt13 + GtDDD312*invgamt23 + GtDDD313*invgamt33
;

const auto
GtDDU321
=
GtDDD312*invgamt11 + GtDDD322*invgamt12 + GtDDD323*invgamt13
;

const auto
GtDDU322
=
GtDDD312*invgamt12 + GtDDD322*invgamt22 + GtDDD323*invgamt23
;

const auto
GtDDU323
=
GtDDD312*invgamt13 + GtDDD322*invgamt23 + GtDDD323*invgamt33
;

const auto
GtDDU331
=
GtDDD313*invgamt11 + GtDDD323*invgamt12 + GtDDD333*invgamt13
;

const auto
GtDDU332
=
GtDDD313*invgamt12 + GtDDD323*invgamt22 + GtDDD333*invgamt23
;

const auto
GtDDU333
=
GtDDD313*invgamt13 + GtDDD323*invgamt23 + GtDDD333*invgamt33
;

const auto
Gt111
=
GtDDD111*invgamt11 + GtDDD211*invgamt12 + GtDDD311*invgamt13
;

const auto
Gt112
=
GtDDD112*invgamt11 + GtDDD212*invgamt12 + GtDDD312*invgamt13
;

const auto
Gt113
=
GtDDD113*invgamt11 + GtDDD213*invgamt12 + GtDDD313*invgamt13
;

const auto
Gt122
=
GtDDD122*invgamt11 + GtDDD222*invgamt12 + GtDDD322*invgamt13
;

const auto
Gt123
=
GtDDD123*invgamt11 + GtDDD223*invgamt12 + GtDDD323*invgamt13
;

const auto
Gt133
=
GtDDD133*invgamt11 + GtDDD233*invgamt12 + GtDDD333*invgamt13
;

const auto
Gt211
=
GtDDD111*invgamt12 + GtDDD211*invgamt22 + GtDDD311*invgamt23
;

const auto
Gt212
=
GtDDD112*invgamt12 + GtDDD212*invgamt22 + GtDDD312*invgamt23
;

const auto
Gt213
=
GtDDD113*invgamt12 + GtDDD213*invgamt22 + GtDDD313*invgamt23
;

const auto
Gt222
=
GtDDD122*invgamt12 + GtDDD222*invgamt22 + GtDDD322*invgamt23
;

const auto
Gt223
=
GtDDD123*invgamt12 + GtDDD223*invgamt22 + GtDDD323*invgamt23
;

const auto
Gt233
=
GtDDD133*invgamt12 + GtDDD233*invgamt22 + GtDDD333*invgamt23
;

const auto
Gt311
=
GtDDD111*invgamt13 + GtDDD211*invgamt23 + GtDDD311*invgamt33
;

const auto
Gt312
=
GtDDD112*invgamt13 + GtDDD212*invgamt23 + GtDDD312*invgamt33
;

const auto
Gt313
=
GtDDD113*invgamt13 + GtDDD213*invgamt23 + GtDDD313*invgamt33
;

const auto
Gt322
=
GtDDD122*invgamt13 + GtDDD222*invgamt23 + GtDDD322*invgamt33
;

const auto
Gt323
=
GtDDD123*invgamt13 + GtDDD223*invgamt23 + GtDDD323*invgamt33
;

const auto
Gt333
=
GtDDD133*invgamt13 + GtDDD233*invgamt23 + GtDDD333*invgamt33
;

const auto
trGtd1
=
Gt111*invgamt11 + 2*Gt112*invgamt12 + 2*Gt113*invgamt13 + Gt122*invgamt22 +
  2*Gt123*invgamt23 + Gt133*invgamt33
;

const auto
trGtd2
=
Gt211*invgamt11 + 2*Gt212*invgamt12 + 2*Gt213*invgamt13 + Gt222*invgamt22 +
  2*Gt223*invgamt23 + Gt233*invgamt33
;

const auto
trGtd3
=
Gt311*invgamt11 + 2*Gt312*invgamt12 + 2*Gt313*invgamt13 + Gt322*invgamt22 +
  2*Gt323*invgamt23 + Gt333*invgamt33
;

const auto
dgam111
=
(dgamt111[ijk5] - 2*dlnW1*gamt11[ijk])/Power(W[ijk],2)
;

const auto
dgam112
=
(dgamt112[ijk5] - 2*dlnW1*gamt12[ijk])/Power(W[ijk],2)
;

const auto
dgam113
=
(dgamt113[ijk5] - 2*dlnW1*gamt13[ijk])/Power(W[ijk],2)
;

const auto
dgam122
=
(dgamt122[ijk5] - 2*dlnW1*gamt22[ijk])/Power(W[ijk],2)
;

const auto
dgam123
=
(dgamt123[ijk5] - 2*dlnW1*gamt23[ijk])/Power(W[ijk],2)
;

const auto
dgam133
=
(dgamt133[ijk5] - 2*dlnW1*gamt33[ijk])/Power(W[ijk],2)
;

const auto
dgam211
=
(dgamt211[ijk5] - 2*dlnW2*gamt11[ijk])/Power(W[ijk],2)
;

const auto
dgam212
=
(dgamt212[ijk5] - 2*dlnW2*gamt12[ijk])/Power(W[ijk],2)
;

const auto
dgam213
=
(dgamt213[ijk5] - 2*dlnW2*gamt13[ijk])/Power(W[ijk],2)
;

const auto
dgam222
=
(dgamt222[ijk5] - 2*dlnW2*gamt22[ijk])/Power(W[ijk],2)
;

const auto
dgam223
=
(dgamt223[ijk5] - 2*dlnW2*gamt23[ijk])/Power(W[ijk],2)
;

const auto
dgam233
=
(dgamt233[ijk5] - 2*dlnW2*gamt33[ijk])/Power(W[ijk],2)
;

const auto
dgam311
=
(dgamt311[ijk5] - 2*dlnW3*gamt11[ijk])/Power(W[ijk],2)
;

const auto
dgam312
=
(dgamt312[ijk5] - 2*dlnW3*gamt12[ijk])/Power(W[ijk],2)
;

const auto
dgam313
=
(dgamt313[ijk5] - 2*dlnW3*gamt13[ijk])/Power(W[ijk],2)
;

const auto
dgam322
=
(dgamt322[ijk5] - 2*dlnW3*gamt22[ijk])/Power(W[ijk],2)
;

const auto
dgam323
=
(dgamt323[ijk5] - 2*dlnW3*gamt23[ijk])/Power(W[ijk],2)
;

const auto
dgam333
=
(dgamt333[ijk5] - 2*dlnW3*gamt33[ijk])/Power(W[ijk],2)
;

const auto
GamDDD111
=
dgam111/2.
;

const auto
GamDDD112
=
dgam211/2.
;

const auto
GamDDD113
=
dgam311/2.
;

const auto
GamDDD122
=
-0.5*dgam122 + dgam212
;

const auto
GamDDD123
=
(-dgam123 + dgam213 + dgam312)/2.
;

const auto
GamDDD133
=
-0.5*dgam133 + dgam313
;

const auto
GamDDD211
=
dgam112 - dgam211/2.
;

const auto
GamDDD212
=
dgam122/2.
;

const auto
GamDDD213
=
(dgam123 - dgam213 + dgam312)/2.
;

const auto
GamDDD222
=
dgam222/2.
;

const auto
GamDDD223
=
dgam322/2.
;

const auto
GamDDD233
=
-0.5*dgam233 + dgam323
;

const auto
GamDDD311
=
dgam113 - dgam311/2.
;

const auto
GamDDD312
=
(dgam123 + dgam213 - dgam312)/2.
;

const auto
GamDDD313
=
dgam133/2.
;

const auto
GamDDD322
=
dgam223 - dgam322/2.
;

const auto
GamDDD323
=
dgam233/2.
;

const auto
GamDDD333
=
dgam333/2.
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
exAtUU11
=
Power(invgamt11,2)*exAt11[ijk] + 2*invgamt11*invgamt12*exAt12[ijk] +
  2*invgamt11*invgamt13*exAt13[ijk] + Power(invgamt12,2)*exAt22[ijk] +
  2*invgamt12*invgamt13*exAt23[ijk] + Power(invgamt13,2)*exAt33[ijk]
;

const auto
exAtUU12
=
invgamt11*invgamt12*exAt11[ijk] +
  (Power(invgamt12,2) + invgamt11*invgamt22)*exAt12[ijk] +
  invgamt12*invgamt13*exAt13[ijk] + invgamt11*invgamt23*exAt13[ijk] +
  invgamt12*invgamt22*exAt22[ijk] + invgamt13*invgamt22*exAt23[ijk] +
  invgamt12*invgamt23*exAt23[ijk] + invgamt13*invgamt23*exAt33[ijk]
;

const auto
exAtUU13
=
invgamt11*invgamt13*exAt11[ijk] +
  (invgamt12*invgamt13 + invgamt11*invgamt23)*exAt12[ijk] +
  Power(invgamt13,2)*exAt13[ijk] + invgamt11*invgamt33*exAt13[ijk] +
  invgamt12*invgamt23*exAt22[ijk] + invgamt13*invgamt23*exAt23[ijk] +
  invgamt12*invgamt33*exAt23[ijk] + invgamt13*invgamt33*exAt33[ijk]
;

const auto
exAtUU22
=
Power(invgamt12,2)*exAt11[ijk] + 2*invgamt12*invgamt22*exAt12[ijk] +
  2*invgamt12*invgamt23*exAt13[ijk] + Power(invgamt22,2)*exAt22[ijk] +
  2*invgamt22*invgamt23*exAt23[ijk] + Power(invgamt23,2)*exAt33[ijk]
;

const auto
exAtUU23
=
invgamt12*invgamt13*exAt11[ijk] +
  (invgamt13*invgamt22 + invgamt12*invgamt23)*exAt12[ijk] +
  invgamt13*invgamt23*exAt13[ijk] + invgamt12*invgamt33*exAt13[ijk] +
  invgamt22*invgamt23*exAt22[ijk] + Power(invgamt23,2)*exAt23[ijk] +
  invgamt22*invgamt33*exAt23[ijk] + invgamt23*invgamt33*exAt33[ijk]
;

const auto
exAtUU33
=
Power(invgamt13,2)*exAt11[ijk] + 2*invgamt13*invgamt23*exAt12[ijk] +
  2*invgamt13*invgamt33*exAt13[ijk] + Power(invgamt23,2)*exAt22[ijk] +
  2*invgamt23*invgamt33*exAt23[ijk] + Power(invgamt33,2)*exAt33[ijk]
;

const auto
tDtDW11
=
ddW11[ijk5] - Gt111*dW1[ijk5] - Gt211*dW2[ijk5] - Gt311*dW3[ijk5]
;

const auto
tDtDW12
=
ddW12[ijk5] - Gt112*dW1[ijk5] - Gt212*dW2[ijk5] - Gt312*dW3[ijk5]
;

const auto
tDtDW13
=
ddW13[ijk5] - Gt113*dW1[ijk5] - Gt213*dW2[ijk5] - Gt313*dW3[ijk5]
;

const auto
tDtDW22
=
ddW22[ijk5] - Gt122*dW1[ijk5] - Gt222*dW2[ijk5] - Gt322*dW3[ijk5]
;

const auto
tDtDW23
=
ddW23[ijk5] - Gt123*dW1[ijk5] - Gt223*dW2[ijk5] - Gt323*dW3[ijk5]
;

const auto
tDtDW33
=
ddW33[ijk5] - Gt133*dW1[ijk5] - Gt233*dW2[ijk5] - Gt333*dW3[ijk5]
;

const auto
DDalpha11
=
-(Gam111*dalpha1[ijk5]) - Gam211*dalpha2[ijk5] - Gam311*dalpha3[ijk5] +
  ddalpha11[ijk5]
;

const auto
DDalpha12
=
-(Gam112*dalpha1[ijk5]) - Gam212*dalpha2[ijk5] - Gam312*dalpha3[ijk5] +
  ddalpha12[ijk5]
;

const auto
DDalpha13
=
-(Gam113*dalpha1[ijk5]) - Gam213*dalpha2[ijk5] - Gam313*dalpha3[ijk5] +
  ddalpha13[ijk5]
;

const auto
DDalpha22
=
-(Gam122*dalpha1[ijk5]) - Gam222*dalpha2[ijk5] - Gam322*dalpha3[ijk5] +
  ddalpha22[ijk5]
;

const auto
DDalpha23
=
-(Gam123*dalpha1[ijk5]) - Gam223*dalpha2[ijk5] - Gam323*dalpha3[ijk5] +
  ddalpha23[ijk5]
;

const auto
DDalpha33
=
-(Gam133*dalpha1[ijk5]) - Gam233*dalpha2[ijk5] - Gam333*dalpha3[ijk5] +
  ddalpha33[ijk5]
;

const auto
RtW11
=
(tDtDW11 + gamt11[ijk]*(invgamt11*tDtDW11 + 2*invgamt12*tDtDW12 +
       2*invgamt13*tDtDW13 + invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 +
       invgamt33*tDtDW33 - 2*(Power(dlnW1,2)*invgamt11 +
          2*dlnW1*dlnW2*invgamt12 + 2*dlnW1*dlnW3*invgamt13 +
          Power(dlnW2,2)*invgamt22 + 2*dlnW2*dlnW3*invgamt23 +
          Power(dlnW3,2)*invgamt33)*W[ijk]))/W[ijk]
;

const auto
RtW12
=
(tDtDW12 + gamt12[ijk]*(invgamt11*tDtDW11 + 2*invgamt12*tDtDW12 +
       2*invgamt13*tDtDW13 + invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 +
       invgamt33*tDtDW33 - 2*(Power(dlnW1,2)*invgamt11 +
          2*dlnW1*dlnW2*invgamt12 + 2*dlnW1*dlnW3*invgamt13 +
          Power(dlnW2,2)*invgamt22 + 2*dlnW2*dlnW3*invgamt23 +
          Power(dlnW3,2)*invgamt33)*W[ijk]))/W[ijk]
;

const auto
RtW13
=
(tDtDW13 + gamt13[ijk]*(invgamt11*tDtDW11 + 2*invgamt12*tDtDW12 +
       2*invgamt13*tDtDW13 + invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 +
       invgamt33*tDtDW33 - 2*(Power(dlnW1,2)*invgamt11 +
          2*dlnW1*dlnW2*invgamt12 + 2*dlnW1*dlnW3*invgamt13 +
          Power(dlnW2,2)*invgamt22 + 2*dlnW2*dlnW3*invgamt23 +
          Power(dlnW3,2)*invgamt33)*W[ijk]))/W[ijk]
;

const auto
RtW22
=
(tDtDW22 + gamt22[ijk]*(invgamt11*tDtDW11 + 2*invgamt12*tDtDW12 +
       2*invgamt13*tDtDW13 + invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 +
       invgamt33*tDtDW33 - 2*(Power(dlnW1,2)*invgamt11 +
          2*dlnW1*dlnW2*invgamt12 + 2*dlnW1*dlnW3*invgamt13 +
          Power(dlnW2,2)*invgamt22 + 2*dlnW2*dlnW3*invgamt23 +
          Power(dlnW3,2)*invgamt33)*W[ijk]))/W[ijk]
;

const auto
RtW23
=
(tDtDW23 + gamt23[ijk]*(invgamt11*tDtDW11 + 2*invgamt12*tDtDW12 +
       2*invgamt13*tDtDW13 + invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 +
       invgamt33*tDtDW33 - 2*(Power(dlnW1,2)*invgamt11 +
          2*dlnW1*dlnW2*invgamt12 + 2*dlnW1*dlnW3*invgamt13 +
          Power(dlnW2,2)*invgamt22 + 2*dlnW2*dlnW3*invgamt23 +
          Power(dlnW3,2)*invgamt33)*W[ijk]))/W[ijk]
;

const auto
RtW33
=
(tDtDW33 + gamt33[ijk]*(invgamt11*tDtDW11 + 2*invgamt12*tDtDW12 +
       2*invgamt13*tDtDW13 + invgamt22*tDtDW22 + 2*invgamt23*tDtDW23 +
       invgamt33*tDtDW33 - 2*(Power(dlnW1,2)*invgamt11 +
          2*dlnW1*dlnW2*invgamt12 + 2*dlnW1*dlnW3*invgamt13 +
          Power(dlnW2,2)*invgamt22 + 2*dlnW2*dlnW3*invgamt23 +
          Power(dlnW3,2)*invgamt33)*W[ijk]))/W[ijk]
;

const auto
Rt11
=
3*Gt111*GtDDU111 + 3*Gt112*GtDDU112 + 3*Gt113*GtDDU113 + 2*Gt211*GtDDU121 +
  2*Gt212*GtDDU122 + 2*Gt213*GtDDU123 + 2*Gt311*GtDDU131 +
  2*Gt312*GtDDU132 + 2*Gt313*GtDDU133 + Gt211*GtDDU211 + Gt212*GtDDU212 +
  Gt213*GtDDU213 + Gt311*GtDDU311 + Gt312*GtDDU312 + Gt313*GtDDU313 +
  GtDDD111*trGtd1 + GtDDD112*trGtd2 + GtDDD113*trGtd3 -
  (invgamt11*ddgamt1111[ijk5])/2. - invgamt12*ddgamt1211[ijk5] -
  invgamt13*ddgamt1311[ijk5] - (invgamt22*ddgamt2211[ijk5])/2. -
  invgamt23*ddgamt2311[ijk5] - (invgamt33*ddgamt3311[ijk5])/2. +
  dtrGt11[ijk5]*gamt11[ijk] + dtrGt12[ijk5]*gamt12[ijk] +
  dtrGt13[ijk5]*gamt13[ijk]
;

const auto
Rt12
=
(2*Gt112*GtDDU111 + 2*Gt122*GtDDU112 + 2*Gt123*GtDDU113 +
    2*Gt111*GtDDU121 + 2*Gt212*GtDDU121 + 2*Gt112*GtDDU122 +
    2*Gt222*GtDDU122 + 2*Gt113*GtDDU123 + 2*Gt223*GtDDU123 +
    2*Gt312*GtDDU131 + 2*Gt322*GtDDU132 + 2*Gt323*GtDDU133 +
    2*Gt111*GtDDU211 + 2*Gt112*GtDDU212 + 2*Gt113*GtDDU213 +
    4*Gt211*GtDDU221 + 4*Gt212*GtDDU222 + 4*Gt213*GtDDU223 +
    2*Gt311*GtDDU231 + 2*Gt312*GtDDU232 + 2*Gt313*GtDDU233 +
    2*Gt311*GtDDU321 + 2*Gt312*GtDDU322 + 2*Gt313*GtDDU323 +
    GtDDD112*trGtd1 + GtDDD211*trGtd1 + GtDDD122*trGtd2 + GtDDD212*trGtd2 +
    GtDDD123*trGtd3 + GtDDD213*trGtd3 - invgamt11*ddgamt1112[ijk5] -
    2*invgamt12*ddgamt1212[ijk5] - 2*invgamt13*ddgamt1312[ijk5] -
    invgamt22*ddgamt2212[ijk5] - 2*invgamt23*ddgamt2312[ijk5] -
    invgamt33*ddgamt3312[ijk5] + dtrGt21[ijk5]*gamt11[ijk] +
    dtrGt11[ijk5]*gamt12[ijk] + dtrGt22[ijk5]*gamt12[ijk] +
    dtrGt23[ijk5]*gamt13[ijk] + dtrGt12[ijk5]*gamt22[ijk] +
    dtrGt13[ijk5]*gamt23[ijk])/2.
;

const auto
Rt13
=
(2*Gt113*GtDDU111 + 2*Gt123*GtDDU112 + 2*Gt133*GtDDU113 +
    2*Gt213*GtDDU121 + 2*Gt223*GtDDU122 + 2*Gt233*GtDDU123 +
    2*Gt111*GtDDU131 + 2*Gt313*GtDDU131 + 2*Gt112*GtDDU132 +
    2*Gt323*GtDDU132 + 2*Gt113*GtDDU133 + 2*Gt333*GtDDU133 +
    2*Gt211*GtDDU231 + 2*Gt212*GtDDU232 + 2*Gt213*GtDDU233 +
    2*Gt111*GtDDU311 + 2*Gt112*GtDDU312 + 2*Gt113*GtDDU313 +
    2*Gt211*GtDDU321 + 2*Gt212*GtDDU322 + 2*Gt213*GtDDU323 +
    4*Gt311*GtDDU331 + 4*Gt312*GtDDU332 + 4*Gt313*GtDDU333 +
    GtDDD113*trGtd1 + GtDDD311*trGtd1 + GtDDD123*trGtd2 + GtDDD312*trGtd2 +
    GtDDD133*trGtd3 + GtDDD313*trGtd3 - invgamt11*ddgamt1113[ijk5] -
    2*invgamt12*ddgamt1213[ijk5] - 2*invgamt13*ddgamt1313[ijk5] -
    invgamt22*ddgamt2213[ijk5] - 2*invgamt23*ddgamt2313[ijk5] -
    invgamt33*ddgamt3313[ijk5] + dtrGt31[ijk5]*gamt11[ijk] +
    dtrGt32[ijk5]*gamt12[ijk] + dtrGt11[ijk5]*gamt13[ijk] +
    dtrGt33[ijk5]*gamt13[ijk] + dtrGt12[ijk5]*gamt23[ijk] +
    dtrGt13[ijk5]*gamt33[ijk])/2.
;

const auto
Rt22
=
Gt112*GtDDU121 + Gt122*GtDDU122 + Gt123*GtDDU123 + 2*Gt112*GtDDU211 +
  2*Gt122*GtDDU212 + 2*Gt123*GtDDU213 + 3*Gt212*GtDDU221 +
  3*Gt222*GtDDU222 + 3*Gt223*GtDDU223 + 2*Gt312*GtDDU231 +
  2*Gt322*GtDDU232 + 2*Gt323*GtDDU233 + Gt312*GtDDU321 + Gt322*GtDDU322 +
  Gt323*GtDDU323 + GtDDD212*trGtd1 + GtDDD222*trGtd2 + GtDDD223*trGtd3 -
  (invgamt11*ddgamt1122[ijk5])/2. - invgamt12*ddgamt1222[ijk5] -
  invgamt13*ddgamt1322[ijk5] - (invgamt22*ddgamt2222[ijk5])/2. -
  invgamt23*ddgamt2322[ijk5] - (invgamt33*ddgamt3322[ijk5])/2. +
  dtrGt21[ijk5]*gamt12[ijk] + dtrGt22[ijk5]*gamt22[ijk] +
  dtrGt23[ijk5]*gamt23[ijk]
;

const auto
Rt23
=
(2*Gt112*GtDDU131 + 2*Gt122*GtDDU132 + 2*Gt123*GtDDU133 +
    2*Gt113*GtDDU211 + 2*Gt123*GtDDU212 + 2*Gt133*GtDDU213 +
    2*Gt213*GtDDU221 + 2*Gt223*GtDDU222 + 2*Gt233*GtDDU223 +
    2*Gt212*GtDDU231 + 2*Gt313*GtDDU231 + 2*Gt222*GtDDU232 +
    2*Gt323*GtDDU232 + 2*Gt223*GtDDU233 + 2*Gt333*GtDDU233 +
    2*Gt112*GtDDU311 + 2*Gt122*GtDDU312 + 2*Gt123*GtDDU313 +
    2*Gt212*GtDDU321 + 2*Gt222*GtDDU322 + 2*Gt223*GtDDU323 +
    4*Gt312*GtDDU331 + 4*Gt322*GtDDU332 + 4*Gt323*GtDDU333 +
    GtDDD213*trGtd1 + GtDDD312*trGtd1 + GtDDD223*trGtd2 + GtDDD322*trGtd2 +
    GtDDD233*trGtd3 + GtDDD323*trGtd3 - invgamt11*ddgamt1123[ijk5] -
    2*invgamt12*ddgamt1223[ijk5] - 2*invgamt13*ddgamt1323[ijk5] -
    invgamt22*ddgamt2223[ijk5] - 2*invgamt23*ddgamt2323[ijk5] -
    invgamt33*ddgamt3323[ijk5] + dtrGt31[ijk5]*gamt12[ijk] +
    dtrGt21[ijk5]*gamt13[ijk] + dtrGt32[ijk5]*gamt22[ijk] +
    dtrGt22[ijk5]*gamt23[ijk] + dtrGt33[ijk5]*gamt23[ijk] +
    dtrGt23[ijk5]*gamt33[ijk])/2.
;

const auto
Rt33
=
Gt113*GtDDU131 + Gt123*GtDDU132 + Gt133*GtDDU133 + Gt213*GtDDU231 +
  Gt223*GtDDU232 + Gt233*GtDDU233 + 2*Gt113*GtDDU311 + 2*Gt123*GtDDU312 +
  2*Gt133*GtDDU313 + 2*Gt213*GtDDU321 + 2*Gt223*GtDDU322 +
  2*Gt233*GtDDU323 + 3*Gt313*GtDDU331 + 3*Gt323*GtDDU332 +
  3*Gt333*GtDDU333 + GtDDD313*trGtd1 + GtDDD323*trGtd2 + GtDDD333*trGtd3 -
  (invgamt11*ddgamt1133[ijk5])/2. - invgamt12*ddgamt1233[ijk5] -
  invgamt13*ddgamt1333[ijk5] - (invgamt22*ddgamt2233[ijk5])/2. -
  invgamt23*ddgamt2333[ijk5] - (invgamt33*ddgamt3333[ijk5])/2. +
  dtrGt31[ijk5]*gamt13[ijk] + dtrGt32[ijk5]*gamt23[ijk] +
  dtrGt33[ijk5]*gamt33[ijk]
;

const auto
R11
=
Rt11 + RtW11
;

const auto
R12
=
Rt12 + RtW12
;

const auto
R13
=
Rt13 + RtW13
;

const auto
R22
=
Rt22 + RtW22
;

const auto
R23
=
Rt23 + RtW23
;

const auto
R33
=
Rt33 + RtW33
;

const auto
trR
=
invgam11*R11 + 2*invgam12*R12 + 2*invgam13*R13 + invgam22*R22 +
  2*invgam23*R23 + invgam33*R33
;

const auto
rho
=
(Power(beta1[ijk],2)*eT11[ijk] + Power(beta2[ijk],2)*eT22[ijk] +
    Power(beta3[ijk],2)*eT33[ijk] +
    2*beta1[ijk]*(beta2[ijk]*eT12[ijk] + beta3[ijk]*eT13[ijk] -
       eTt1[ijk]) + 2*beta2[ijk]*(beta3[ijk]*eT23[ijk] - eTt2[ijk]) -
    2*beta3[ijk]*eTt3[ijk] + eTtt[ijk])/Power(alpha[ijk],2)
;

const auto
Sm1
=
(beta1[ijk]*eT11[ijk] + beta2[ijk]*eT12[ijk] + beta3[ijk]*eT13[ijk] -
    eTt1[ijk])/alpha[ijk]
;

const auto
Sm2
=
(beta1[ijk]*eT12[ijk] + beta2[ijk]*eT22[ijk] + beta3[ijk]*eT23[ijk] -
    eTt2[ijk])/alpha[ijk]
;

const auto
Sm3
=
(beta1[ijk]*eT13[ijk] + beta2[ijk]*eT23[ijk] + beta3[ijk]*eT33[ijk] -
    eTt3[ijk])/alpha[ijk]
;

const auto
Ss11
=
eT11[ijk]
;

const auto
Ss12
=
eT12[ijk]
;

const auto
Ss13
=
eT13[ijk]
;

const auto
Ss22
=
eT22[ijk]
;

const auto
Ss23
=
eT23[ijk]
;

const auto
Ss33
=
eT33[ijk]
;

const auto
trSs
=
invgam11*Ss11 + 2*invgam12*Ss12 + 2*invgam13*Ss13 + invgam22*Ss22 +
  2*invgam23*Ss23 + invgam33*Ss33
;


dtW[ijk]
=
beta1[ijk]*dW1[ijk5] + beta2[ijk]*dW2[ijk5] + beta3[ijk]*dW3[ijk5] -
  ((dbeta11[ijk5] + dbeta22[ijk5] + dbeta33[ijk5])*W[ijk])/3. +
  (alpha[ijk]*exKh[ijk]*W[ijk])/3. + (2*alpha[ijk]*Theta[ijk]*W[ijk])/3.
;

dtgamt11[ijk]
=
beta1[ijk]*dgamt111[ijk5] + beta2[ijk]*dgamt211[ijk5] +
  beta3[ijk]*dgamt311[ijk5] - 2*alpha[ijk]*exAt11[ijk] +
  2*dbeta11[ijk5]*gamt11[ijk] -
  (2*(dbeta11[ijk5] + dbeta22[ijk5] + dbeta33[ijk5])*gamt11[ijk])/3. +
  2*dbeta12[ijk5]*gamt12[ijk] + 2*dbeta13[ijk5]*gamt13[ijk]
;

dtgamt12[ijk]
=
beta1[ijk]*dgamt112[ijk5] + beta2[ijk]*dgamt212[ijk5] +
  beta3[ijk]*dgamt312[ijk5] - 2*alpha[ijk]*exAt12[ijk] +
  dbeta21[ijk5]*gamt11[ijk] + (dbeta11[ijk5]*gamt12[ijk])/3. +
  (dbeta22[ijk5]*gamt12[ijk])/3. - (2*dbeta33[ijk5]*gamt12[ijk])/3. +
  dbeta23[ijk5]*gamt13[ijk] + dbeta12[ijk5]*gamt22[ijk] +
  dbeta13[ijk5]*gamt23[ijk]
;

dtgamt13[ijk]
=
beta1[ijk]*dgamt113[ijk5] + beta2[ijk]*dgamt213[ijk5] +
  beta3[ijk]*dgamt313[ijk5] - 2*alpha[ijk]*exAt13[ijk] +
  dbeta31[ijk5]*gamt11[ijk] + dbeta32[ijk5]*gamt12[ijk] +
  (dbeta11[ijk5]*gamt13[ijk])/3. - (2*dbeta22[ijk5]*gamt13[ijk])/3. +
  (dbeta33[ijk5]*gamt13[ijk])/3. + dbeta12[ijk5]*gamt23[ijk] +
  dbeta13[ijk5]*gamt33[ijk]
;

dtgamt22[ijk]
=
beta1[ijk]*dgamt122[ijk5] + beta2[ijk]*dgamt222[ijk5] +
  beta3[ijk]*dgamt322[ijk5] - 2*alpha[ijk]*exAt22[ijk] +
  2*dbeta21[ijk5]*gamt12[ijk] + 2*dbeta22[ijk5]*gamt22[ijk] -
  (2*(dbeta11[ijk5] + dbeta22[ijk5] + dbeta33[ijk5])*gamt22[ijk])/3. +
  2*dbeta23[ijk5]*gamt23[ijk]
;

dtgamt23[ijk]
=
beta1[ijk]*dgamt123[ijk5] + beta2[ijk]*dgamt223[ijk5] +
  beta3[ijk]*dgamt323[ijk5] - 2*alpha[ijk]*exAt23[ijk] +
  dbeta31[ijk5]*gamt12[ijk] + dbeta21[ijk5]*gamt13[ijk] +
  dbeta32[ijk5]*gamt22[ijk] - (2*dbeta11[ijk5]*gamt23[ijk])/3. +
  (dbeta22[ijk5]*gamt23[ijk])/3. + (dbeta33[ijk5]*gamt23[ijk])/3. +
  dbeta23[ijk5]*gamt33[ijk]
;

dtgamt33[ijk]
=
beta1[ijk]*dgamt133[ijk5] + beta2[ijk]*dgamt233[ijk5] +
  beta3[ijk]*dgamt333[ijk5] - 2*alpha[ijk]*exAt33[ijk] +
  2*dbeta31[ijk5]*gamt13[ijk] + 2*dbeta32[ijk5]*gamt23[ijk] +
  2*dbeta33[ijk5]*gamt33[ijk] -
  (2*(dbeta11[ijk5] + dbeta22[ijk5] + dbeta33[ijk5])*gamt33[ijk])/3.
;

dtexKh[ijk]
=
-(DDalpha11*invgam11) - 2*DDalpha12*invgam12 - 2*DDalpha13*invgam13 -
  DDalpha22*invgam22 - 2*DDalpha23*invgam23 - DDalpha33*invgam33 +
  beta1[ijk]*dexKh1[ijk5] + beta2[ijk]*dexKh2[ijk5] +
  beta3[ijk]*dexKh3[ijk5] + alpha[ijk]*
   (4*cpi*rho + 4*cpi*trSs + exAtUU11*exAt11[ijk] + 2*exAtUU12*exAt12[ijk] +
     2*exAtUU13*exAt13[ijk] + exAtUU22*exAt22[ijk] +
     2*exAtUU23*exAt23[ijk] + exAtUU33*exAt33[ijk] + Power(exKh[ijk],2)/3. +
     ckappa1*Theta[ijk] - ckappa1*ckappa2*Theta[ijk] +
     (4*exKh[ijk]*Theta[ijk])/3. + (4*Power(Theta[ijk],2))/3.)
;

dtexAt11[ijk]
=
(3*beta1[ijk]*dexAt111[ijk5] + 3*beta2[ijk]*dexAt211[ijk5] +
    3*beta3[ijk]*dexAt311[ijk5] + 4*dbeta11[ijk5]*exAt11[ijk] -
    2*dbeta22[ijk5]*exAt11[ijk] - 2*dbeta33[ijk5]*exAt11[ijk] -
    6*invgamt11*alpha[ijk]*Power(exAt11[ijk],2) +
    6*dbeta12[ijk5]*exAt12[ijk] -
    12*invgamt12*alpha[ijk]*exAt11[ijk]*exAt12[ijk] -
    6*invgamt22*alpha[ijk]*Power(exAt12[ijk],2) +
    6*dbeta13[ijk5]*exAt13[ijk] -
    12*invgamt13*alpha[ijk]*exAt11[ijk]*exAt13[ijk] -
    12*invgamt23*alpha[ijk]*exAt12[ijk]*exAt13[ijk] -
    6*invgamt33*alpha[ijk]*Power(exAt13[ijk],2) +
    3*alpha[ijk]*exAt11[ijk]*exKh[ijk] +
    6*alpha[ijk]*exAt11[ijk]*Theta[ijk] - 3*DDalpha11*Power(W[ijk],2) +
    DDalpha11*gam11*invgam11*Power(W[ijk],2) +
    2*DDalpha12*gam11*invgam12*Power(W[ijk],2) +
    2*DDalpha13*gam11*invgam13*Power(W[ijk],2) +
    DDalpha22*gam11*invgam22*Power(W[ijk],2) +
    2*DDalpha23*gam11*invgam23*Power(W[ijk],2) +
    DDalpha33*gam11*invgam33*Power(W[ijk],2) +
    3*R11*alpha[ijk]*Power(W[ijk],2) -
    gam11*invgam11*R11*alpha[ijk]*Power(W[ijk],2) -
    2*gam11*invgam12*R12*alpha[ijk]*Power(W[ijk],2) -
    2*gam11*invgam13*R13*alpha[ijk]*Power(W[ijk],2) -
    gam11*invgam22*R22*alpha[ijk]*Power(W[ijk],2) -
    2*gam11*invgam23*R23*alpha[ijk]*Power(W[ijk],2) -
    gam11*invgam33*R33*alpha[ijk]*Power(W[ijk],2) -
    24*cpi*Ss11*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam11*invgam11*Ss11*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam11*invgam12*Ss12*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam11*invgam13*Ss13*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam11*invgam22*Ss22*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam11*invgam23*Ss23*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam11*invgam33*Ss33*alpha[ijk]*Power(W[ijk],2))/3.
;

dtexAt12[ijk]
=
(3*beta1[ijk]*dexAt112[ijk5] + 3*beta2[ijk]*dexAt212[ijk5] +
    3*beta3[ijk]*dexAt312[ijk5] + 3*dbeta21[ijk5]*exAt11[ijk] +
    dbeta11[ijk5]*exAt12[ijk] + dbeta22[ijk5]*exAt12[ijk] -
    2*dbeta33[ijk5]*exAt12[ijk] -
    6*invgamt11*alpha[ijk]*exAt11[ijk]*exAt12[ijk] -
    6*invgamt12*alpha[ijk]*Power(exAt12[ijk],2) +
    3*dbeta23[ijk5]*exAt13[ijk] -
    6*invgamt13*alpha[ijk]*exAt12[ijk]*exAt13[ijk] +
    3*dbeta12[ijk5]*exAt22[ijk] -
    6*invgamt12*alpha[ijk]*exAt11[ijk]*exAt22[ijk] -
    6*invgamt22*alpha[ijk]*exAt12[ijk]*exAt22[ijk] -
    6*invgamt23*alpha[ijk]*exAt13[ijk]*exAt22[ijk] +
    3*dbeta13[ijk5]*exAt23[ijk] -
    6*invgamt13*alpha[ijk]*exAt11[ijk]*exAt23[ijk] -
    6*invgamt23*alpha[ijk]*exAt12[ijk]*exAt23[ijk] -
    6*invgamt33*alpha[ijk]*exAt13[ijk]*exAt23[ijk] +
    3*alpha[ijk]*exAt12[ijk]*exKh[ijk] +
    6*alpha[ijk]*exAt12[ijk]*Theta[ijk] - 3*DDalpha12*Power(W[ijk],2) +
    DDalpha11*gam12*invgam11*Power(W[ijk],2) +
    2*DDalpha12*gam12*invgam12*Power(W[ijk],2) +
    2*DDalpha13*gam12*invgam13*Power(W[ijk],2) +
    DDalpha22*gam12*invgam22*Power(W[ijk],2) +
    2*DDalpha23*gam12*invgam23*Power(W[ijk],2) +
    DDalpha33*gam12*invgam33*Power(W[ijk],2) -
    gam12*invgam11*R11*alpha[ijk]*Power(W[ijk],2) +
    3*R12*alpha[ijk]*Power(W[ijk],2) -
    2*gam12*invgam12*R12*alpha[ijk]*Power(W[ijk],2) -
    2*gam12*invgam13*R13*alpha[ijk]*Power(W[ijk],2) -
    gam12*invgam22*R22*alpha[ijk]*Power(W[ijk],2) -
    2*gam12*invgam23*R23*alpha[ijk]*Power(W[ijk],2) -
    gam12*invgam33*R33*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam12*invgam11*Ss11*alpha[ijk]*Power(W[ijk],2) -
    24*cpi*Ss12*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam12*invgam12*Ss12*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam12*invgam13*Ss13*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam12*invgam22*Ss22*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam12*invgam23*Ss23*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam12*invgam33*Ss33*alpha[ijk]*Power(W[ijk],2))/3.
;

dtexAt13[ijk]
=
(3*beta1[ijk]*dexAt113[ijk5] + 3*beta2[ijk]*dexAt213[ijk5] +
    3*beta3[ijk]*dexAt313[ijk5] + 3*dbeta31[ijk5]*exAt11[ijk] +
    3*dbeta32[ijk5]*exAt12[ijk] + dbeta11[ijk5]*exAt13[ijk] -
    2*dbeta22[ijk5]*exAt13[ijk] + dbeta33[ijk5]*exAt13[ijk] -
    6*invgamt11*alpha[ijk]*exAt11[ijk]*exAt13[ijk] -
    6*invgamt12*alpha[ijk]*exAt12[ijk]*exAt13[ijk] -
    6*invgamt13*alpha[ijk]*Power(exAt13[ijk],2) +
    3*dbeta12[ijk5]*exAt23[ijk] -
    6*invgamt12*alpha[ijk]*exAt11[ijk]*exAt23[ijk] -
    6*invgamt22*alpha[ijk]*exAt12[ijk]*exAt23[ijk] -
    6*invgamt23*alpha[ijk]*exAt13[ijk]*exAt23[ijk] +
    3*dbeta13[ijk5]*exAt33[ijk] -
    6*invgamt13*alpha[ijk]*exAt11[ijk]*exAt33[ijk] -
    6*invgamt23*alpha[ijk]*exAt12[ijk]*exAt33[ijk] -
    6*invgamt33*alpha[ijk]*exAt13[ijk]*exAt33[ijk] +
    3*alpha[ijk]*exAt13[ijk]*exKh[ijk] +
    6*alpha[ijk]*exAt13[ijk]*Theta[ijk] - 3*DDalpha13*Power(W[ijk],2) +
    DDalpha11*gam13*invgam11*Power(W[ijk],2) +
    2*DDalpha12*gam13*invgam12*Power(W[ijk],2) +
    2*DDalpha13*gam13*invgam13*Power(W[ijk],2) +
    DDalpha22*gam13*invgam22*Power(W[ijk],2) +
    2*DDalpha23*gam13*invgam23*Power(W[ijk],2) +
    DDalpha33*gam13*invgam33*Power(W[ijk],2) -
    gam13*invgam11*R11*alpha[ijk]*Power(W[ijk],2) -
    2*gam13*invgam12*R12*alpha[ijk]*Power(W[ijk],2) +
    3*R13*alpha[ijk]*Power(W[ijk],2) -
    2*gam13*invgam13*R13*alpha[ijk]*Power(W[ijk],2) -
    gam13*invgam22*R22*alpha[ijk]*Power(W[ijk],2) -
    2*gam13*invgam23*R23*alpha[ijk]*Power(W[ijk],2) -
    gam13*invgam33*R33*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam13*invgam11*Ss11*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam13*invgam12*Ss12*alpha[ijk]*Power(W[ijk],2) -
    24*cpi*Ss13*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam13*invgam13*Ss13*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam13*invgam22*Ss22*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam13*invgam23*Ss23*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam13*invgam33*Ss33*alpha[ijk]*Power(W[ijk],2))/3.
;

dtexAt22[ijk]
=
(3*beta1[ijk]*dexAt122[ijk5] + 3*beta2[ijk]*dexAt222[ijk5] +
    3*beta3[ijk]*dexAt322[ijk5] + 6*dbeta21[ijk5]*exAt12[ijk] -
    6*invgamt11*alpha[ijk]*Power(exAt12[ijk],2) -
    2*dbeta11[ijk5]*exAt22[ijk] + 4*dbeta22[ijk5]*exAt22[ijk] -
    2*dbeta33[ijk5]*exAt22[ijk] -
    12*invgamt12*alpha[ijk]*exAt12[ijk]*exAt22[ijk] -
    6*invgamt22*alpha[ijk]*Power(exAt22[ijk],2) +
    6*dbeta23[ijk5]*exAt23[ijk] -
    12*invgamt13*alpha[ijk]*exAt12[ijk]*exAt23[ijk] -
    12*invgamt23*alpha[ijk]*exAt22[ijk]*exAt23[ijk] -
    6*invgamt33*alpha[ijk]*Power(exAt23[ijk],2) +
    3*alpha[ijk]*exAt22[ijk]*exKh[ijk] +
    6*alpha[ijk]*exAt22[ijk]*Theta[ijk] - 3*DDalpha22*Power(W[ijk],2) +
    DDalpha11*gam22*invgam11*Power(W[ijk],2) +
    2*DDalpha12*gam22*invgam12*Power(W[ijk],2) +
    2*DDalpha13*gam22*invgam13*Power(W[ijk],2) +
    DDalpha22*gam22*invgam22*Power(W[ijk],2) +
    2*DDalpha23*gam22*invgam23*Power(W[ijk],2) +
    DDalpha33*gam22*invgam33*Power(W[ijk],2) -
    gam22*invgam11*R11*alpha[ijk]*Power(W[ijk],2) -
    2*gam22*invgam12*R12*alpha[ijk]*Power(W[ijk],2) -
    2*gam22*invgam13*R13*alpha[ijk]*Power(W[ijk],2) +
    3*R22*alpha[ijk]*Power(W[ijk],2) -
    gam22*invgam22*R22*alpha[ijk]*Power(W[ijk],2) -
    2*gam22*invgam23*R23*alpha[ijk]*Power(W[ijk],2) -
    gam22*invgam33*R33*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam22*invgam11*Ss11*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam22*invgam12*Ss12*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam22*invgam13*Ss13*alpha[ijk]*Power(W[ijk],2) -
    24*cpi*Ss22*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam22*invgam22*Ss22*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam22*invgam23*Ss23*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam22*invgam33*Ss33*alpha[ijk]*Power(W[ijk],2))/3.
;

dtexAt23[ijk]
=
(3*beta1[ijk]*dexAt123[ijk5] + 3*beta2[ijk]*dexAt223[ijk5] +
    3*beta3[ijk]*dexAt323[ijk5] + 3*dbeta31[ijk5]*exAt12[ijk] +
    3*dbeta21[ijk5]*exAt13[ijk] -
    6*invgamt11*alpha[ijk]*exAt12[ijk]*exAt13[ijk] +
    3*dbeta32[ijk5]*exAt22[ijk] -
    6*invgamt12*alpha[ijk]*exAt13[ijk]*exAt22[ijk] -
    2*dbeta11[ijk5]*exAt23[ijk] + dbeta22[ijk5]*exAt23[ijk] +
    dbeta33[ijk5]*exAt23[ijk] -
    6*invgamt12*alpha[ijk]*exAt12[ijk]*exAt23[ijk] -
    6*invgamt13*alpha[ijk]*exAt13[ijk]*exAt23[ijk] -
    6*invgamt22*alpha[ijk]*exAt22[ijk]*exAt23[ijk] -
    6*invgamt23*alpha[ijk]*Power(exAt23[ijk],2) +
    3*dbeta23[ijk5]*exAt33[ijk] -
    6*invgamt13*alpha[ijk]*exAt12[ijk]*exAt33[ijk] -
    6*invgamt23*alpha[ijk]*exAt22[ijk]*exAt33[ijk] -
    6*invgamt33*alpha[ijk]*exAt23[ijk]*exAt33[ijk] +
    3*alpha[ijk]*exAt23[ijk]*exKh[ijk] +
    6*alpha[ijk]*exAt23[ijk]*Theta[ijk] - 3*DDalpha23*Power(W[ijk],2) +
    DDalpha11*gam23*invgam11*Power(W[ijk],2) +
    2*DDalpha12*gam23*invgam12*Power(W[ijk],2) +
    2*DDalpha13*gam23*invgam13*Power(W[ijk],2) +
    DDalpha22*gam23*invgam22*Power(W[ijk],2) +
    2*DDalpha23*gam23*invgam23*Power(W[ijk],2) +
    DDalpha33*gam23*invgam33*Power(W[ijk],2) -
    gam23*invgam11*R11*alpha[ijk]*Power(W[ijk],2) -
    2*gam23*invgam12*R12*alpha[ijk]*Power(W[ijk],2) -
    2*gam23*invgam13*R13*alpha[ijk]*Power(W[ijk],2) -
    gam23*invgam22*R22*alpha[ijk]*Power(W[ijk],2) +
    3*R23*alpha[ijk]*Power(W[ijk],2) -
    2*gam23*invgam23*R23*alpha[ijk]*Power(W[ijk],2) -
    gam23*invgam33*R33*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam23*invgam11*Ss11*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam23*invgam12*Ss12*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam23*invgam13*Ss13*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam23*invgam22*Ss22*alpha[ijk]*Power(W[ijk],2) -
    24*cpi*Ss23*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam23*invgam23*Ss23*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam23*invgam33*Ss33*alpha[ijk]*Power(W[ijk],2))/3.
;

dtexAt33[ijk]
=
(3*beta1[ijk]*dexAt133[ijk5] + 3*beta2[ijk]*dexAt233[ijk5] +
    3*beta3[ijk]*dexAt333[ijk5] + 6*dbeta31[ijk5]*exAt13[ijk] -
    6*invgamt11*alpha[ijk]*Power(exAt13[ijk],2) +
    6*dbeta32[ijk5]*exAt23[ijk] -
    12*invgamt12*alpha[ijk]*exAt13[ijk]*exAt23[ijk] -
    6*invgamt22*alpha[ijk]*Power(exAt23[ijk],2) -
    2*dbeta11[ijk5]*exAt33[ijk] - 2*dbeta22[ijk5]*exAt33[ijk] +
    4*dbeta33[ijk5]*exAt33[ijk] -
    12*invgamt13*alpha[ijk]*exAt13[ijk]*exAt33[ijk] -
    12*invgamt23*alpha[ijk]*exAt23[ijk]*exAt33[ijk] -
    6*invgamt33*alpha[ijk]*Power(exAt33[ijk],2) +
    3*alpha[ijk]*exAt33[ijk]*exKh[ijk] +
    6*alpha[ijk]*exAt33[ijk]*Theta[ijk] - 3*DDalpha33*Power(W[ijk],2) +
    DDalpha11*gam33*invgam11*Power(W[ijk],2) +
    2*DDalpha12*gam33*invgam12*Power(W[ijk],2) +
    2*DDalpha13*gam33*invgam13*Power(W[ijk],2) +
    DDalpha22*gam33*invgam22*Power(W[ijk],2) +
    2*DDalpha23*gam33*invgam23*Power(W[ijk],2) +
    DDalpha33*gam33*invgam33*Power(W[ijk],2) -
    gam33*invgam11*R11*alpha[ijk]*Power(W[ijk],2) -
    2*gam33*invgam12*R12*alpha[ijk]*Power(W[ijk],2) -
    2*gam33*invgam13*R13*alpha[ijk]*Power(W[ijk],2) -
    gam33*invgam22*R22*alpha[ijk]*Power(W[ijk],2) -
    2*gam33*invgam23*R23*alpha[ijk]*Power(W[ijk],2) +
    3*R33*alpha[ijk]*Power(W[ijk],2) -
    gam33*invgam33*R33*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam33*invgam11*Ss11*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam33*invgam12*Ss12*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam33*invgam13*Ss13*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam33*invgam22*Ss22*alpha[ijk]*Power(W[ijk],2) +
    16*cpi*gam33*invgam23*Ss23*alpha[ijk]*Power(W[ijk],2) -
    24*cpi*Ss33*alpha[ijk]*Power(W[ijk],2) +
    8*cpi*gam33*invgam33*Ss33*alpha[ijk]*Power(W[ijk],2))/3.
;

dttrGt1[ijk]
=
(-6*exAtUU11*dalpha1[ijk5] - 6*exAtUU12*dalpha2[ijk5] -
    6*exAtUU13*dalpha3[ijk5] - trGtd1*dbeta11[ijk5] -
    3*trGtd2*dbeta21[ijk5] + 2*trGtd1*dbeta22[ijk5] -
    3*trGtd3*dbeta31[ijk5] + 2*trGtd1*dbeta33[ijk5] +
    4*invgamt11*ddbeta111[ijk5] + 7*invgamt12*ddbeta121[ijk5] +
    invgamt11*ddbeta122[ijk5] + 7*invgamt13*ddbeta131[ijk5] +
    invgamt11*ddbeta133[ijk5] + 3*invgamt22*ddbeta221[ijk5] +
    invgamt12*ddbeta222[ijk5] + 6*invgamt23*ddbeta231[ijk5] +
    invgamt13*ddbeta232[ijk5] + invgamt12*ddbeta233[ijk5] +
    3*invgamt33*ddbeta331[ijk5] + invgamt13*ddbeta333[ijk5] +
    3*beta1[ijk]*dtrGt11[ijk5] + 3*beta2[ijk]*dtrGt21[ijk5] +
    3*beta3[ijk]*dtrGt31[ijk5] -
    2*alpha[ijk]*(9*dlnW1*exAtUU11 + 9*dlnW2*exAtUU12 + 9*dlnW3*exAtUU13 -
       3*exAtUU11*Gt111 - 6*exAtUU12*Gt112 - 6*exAtUU13*Gt113 -
       3*exAtUU22*Gt122 - 6*exAtUU23*Gt123 - 3*exAtUU33*Gt133 +
       24*cpi*invgamt11*Sm1 + 24*cpi*invgamt12*Sm2 + 24*cpi*invgamt13*Sm3 -
       3*ckappa1*trGtd1 + 2*invgamt11*dexKh1[ijk5] +
       2*invgamt12*dexKh2[ijk5] + 2*invgamt13*dexKh3[ijk5] +
       invgamt11*dTheta1[ijk5] + invgamt12*dTheta2[ijk5] +
       invgamt13*dTheta3[ijk5] + 3*ckappa1*trGt1[ijk]))/3.
;

dttrGt2[ijk]
=
(-6*exAtUU12*dalpha1[ijk5] - 6*exAtUU22*dalpha2[ijk5] -
    6*exAtUU23*dalpha3[ijk5] + 2*trGtd2*dbeta11[ijk5] -
    3*trGtd1*dbeta12[ijk5] - trGtd2*dbeta22[ijk5] -
    3*trGtd3*dbeta32[ijk5] + 2*trGtd2*dbeta33[ijk5] +
    invgamt12*ddbeta111[ijk5] + 3*invgamt11*ddbeta112[ijk5] +
    invgamt22*ddbeta121[ijk5] + 7*invgamt12*ddbeta122[ijk5] +
    invgamt23*ddbeta131[ijk5] + 6*invgamt13*ddbeta132[ijk5] +
    invgamt12*ddbeta133[ijk5] + 4*invgamt22*ddbeta222[ijk5] +
    7*invgamt23*ddbeta232[ijk5] + invgamt22*ddbeta233[ijk5] +
    3*invgamt33*ddbeta332[ijk5] + invgamt23*ddbeta333[ijk5] +
    3*beta1[ijk]*dtrGt12[ijk5] + 3*beta2[ijk]*dtrGt22[ijk5] +
    3*beta3[ijk]*dtrGt32[ijk5] -
    2*alpha[ijk]*(9*dlnW1*exAtUU12 + 9*dlnW2*exAtUU22 + 9*dlnW3*exAtUU23 -
       3*exAtUU11*Gt211 - 6*exAtUU12*Gt212 - 6*exAtUU13*Gt213 -
       3*exAtUU22*Gt222 - 6*exAtUU23*Gt223 - 3*exAtUU33*Gt233 +
       24*cpi*invgamt12*Sm1 + 24*cpi*invgamt22*Sm2 + 24*cpi*invgamt23*Sm3 -
       3*ckappa1*trGtd2 + 2*invgamt12*dexKh1[ijk5] +
       2*invgamt22*dexKh2[ijk5] + 2*invgamt23*dexKh3[ijk5] +
       invgamt12*dTheta1[ijk5] + invgamt22*dTheta2[ijk5] +
       invgamt23*dTheta3[ijk5] + 3*ckappa1*trGt2[ijk]))/3.
;

dttrGt3[ijk]
=
(-6*exAtUU13*dalpha1[ijk5] - 6*exAtUU23*dalpha2[ijk5] -
    6*exAtUU33*dalpha3[ijk5] + 2*trGtd3*dbeta11[ijk5] -
    3*trGtd1*dbeta13[ijk5] + 2*trGtd3*dbeta22[ijk5] -
    3*trGtd2*dbeta23[ijk5] - trGtd3*dbeta33[ijk5] +
    invgamt13*ddbeta111[ijk5] + 3*invgamt11*ddbeta113[ijk5] +
    invgamt23*ddbeta121[ijk5] + invgamt13*ddbeta122[ijk5] +
    6*invgamt12*ddbeta123[ijk5] + invgamt33*ddbeta131[ijk5] +
    7*invgamt13*ddbeta133[ijk5] + invgamt23*ddbeta222[ijk5] +
    3*invgamt22*ddbeta223[ijk5] + invgamt33*ddbeta232[ijk5] +
    7*invgamt23*ddbeta233[ijk5] + 4*invgamt33*ddbeta333[ijk5] +
    3*beta1[ijk]*dtrGt13[ijk5] + 3*beta2[ijk]*dtrGt23[ijk5] +
    3*beta3[ijk]*dtrGt33[ijk5] -
    2*alpha[ijk]*(9*dlnW1*exAtUU13 + 9*dlnW2*exAtUU23 + 9*dlnW3*exAtUU33 -
       3*exAtUU11*Gt311 - 6*exAtUU12*Gt312 - 6*exAtUU13*Gt313 -
       3*exAtUU22*Gt322 - 6*exAtUU23*Gt323 - 3*exAtUU33*Gt333 +
       24*cpi*invgamt13*Sm1 + 24*cpi*invgamt23*Sm2 + 24*cpi*invgamt33*Sm3 -
       3*ckappa1*trGtd3 + 2*invgamt13*dexKh1[ijk5] +
       2*invgamt23*dexKh2[ijk5] + 2*invgamt33*dexKh3[ijk5] +
       invgamt13*dTheta1[ijk5] + invgamt23*dTheta2[ijk5] +
       invgamt33*dTheta3[ijk5] + 3*ckappa1*trGt3[ijk]))/3.
;

dtTheta[ijk]
=
beta1[ijk]*dTheta1[ijk5] + beta2[ijk]*dTheta2[ijk5] +
  beta3[ijk]*dTheta3[ijk5] - (alpha[ijk]*
     (48*cpi*rho - 3*trR + 3*exAtUU11*exAt11[ijk] +
       6*exAtUU12*exAt12[ijk] + 6*exAtUU13*exAt13[ijk] +
       3*exAtUU22*exAt22[ijk] + 6*exAtUU23*exAt23[ijk] +
       3*exAtUU33*exAt33[ijk] - 2*Power(exKh[ijk],2) +
       12*ckappa1*Theta[ijk] + 6*ckappa1*ckappa2*Theta[ijk] -
       8*exKh[ijk]*Theta[ijk] - 8*Power(Theta[ijk],2)))/6.
;

dtalpha[ijk]
=
beta1[ijk]*dalpha1[ijk5] + beta2[ijk]*dalpha2[ijk5] +
  beta3[ijk]*dalpha3[ijk5] - cmuL*alpha[ijk]*exKh[ijk]
;

dtbeta1[ijk]
=
beta1[ijk]*(-ceta + dbeta11[ijk5]) + beta2[ijk]*dbeta21[ijk5] +
  beta3[ijk]*dbeta31[ijk5] + cmuS*trGt1[ijk]
;

dtbeta2[ijk]
=
beta1[ijk]*dbeta12[ijk5] + beta2[ijk]*(-ceta + dbeta22[ijk5]) +
  beta3[ijk]*dbeta32[ijk5] + cmuS*trGt2[ijk]
;

dtbeta3[ijk]
=
beta1[ijk]*dbeta13[ijk5] + beta2[ijk]*dbeta23[ijk5] +
  beta3[ijk]*(-ceta + dbeta33[ijk5]) + cmuS*trGt3[ijk]
;


});
});

#endif // #ifndef Z4COWGPU_SET_RHS_GF3D5_HXX

/* Z4cowGPU_set_rhs_GF3D5.hxx */
