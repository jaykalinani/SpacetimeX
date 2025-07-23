/* Z4cowGPU_set_constraint.hxx */
/* Produced with Generato */

#ifndef Z4COWGPU_SET_CONSTRAINT_HXX
#define Z4COWGPU_SET_CONSTRAINT_HXX

const auto ZtC1 = gf_ZtC[0];
const auto ZtC2 = gf_ZtC[1];
const auto ZtC3 = gf_ZtC[2];
const auto MtC1 = gf_MtC[0];
const auto MtC2 = gf_MtC[1];
const auto MtC3 = gf_MtC[2];
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

noinline([&]() __attribute__((__flatten__, __hot__)) {
grid.loop_int_device<0, 0, 0>(
  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
const int ijk = layout2.linear(p.i, p.j, p.k);

const auto dW1 = calcderivs1_1(W, p.i, p.j, p.k);
const auto dW2 = calcderivs1_2(W, p.i, p.j, p.k);
const auto dW3 = calcderivs1_3(W, p.i, p.j, p.k);
const auto dgamt111 = calcderivs1_1(gamt11, p.i, p.j, p.k);
const auto dgamt112 = calcderivs1_1(gamt12, p.i, p.j, p.k);
const auto dgamt113 = calcderivs1_1(gamt13, p.i, p.j, p.k);
const auto dgamt122 = calcderivs1_1(gamt22, p.i, p.j, p.k);
const auto dgamt123 = calcderivs1_1(gamt23, p.i, p.j, p.k);
const auto dgamt133 = calcderivs1_1(gamt33, p.i, p.j, p.k);
const auto dgamt211 = calcderivs1_2(gamt11, p.i, p.j, p.k);
const auto dgamt212 = calcderivs1_2(gamt12, p.i, p.j, p.k);
const auto dgamt213 = calcderivs1_2(gamt13, p.i, p.j, p.k);
const auto dgamt222 = calcderivs1_2(gamt22, p.i, p.j, p.k);
const auto dgamt223 = calcderivs1_2(gamt23, p.i, p.j, p.k);
const auto dgamt233 = calcderivs1_2(gamt33, p.i, p.j, p.k);
const auto dgamt311 = calcderivs1_3(gamt11, p.i, p.j, p.k);
const auto dgamt312 = calcderivs1_3(gamt12, p.i, p.j, p.k);
const auto dgamt313 = calcderivs1_3(gamt13, p.i, p.j, p.k);
const auto dgamt322 = calcderivs1_3(gamt22, p.i, p.j, p.k);
const auto dgamt323 = calcderivs1_3(gamt23, p.i, p.j, p.k);
const auto dgamt333 = calcderivs1_3(gamt33, p.i, p.j, p.k);
const auto dexKh1 = calcderivs1_1(exKh, p.i, p.j, p.k);
const auto dexKh2 = calcderivs1_2(exKh, p.i, p.j, p.k);
const auto dexKh3 = calcderivs1_3(exKh, p.i, p.j, p.k);
const auto dexAt111 = calcderivs1_1(exAt11, p.i, p.j, p.k);
const auto dexAt112 = calcderivs1_1(exAt12, p.i, p.j, p.k);
const auto dexAt113 = calcderivs1_1(exAt13, p.i, p.j, p.k);
const auto dexAt122 = calcderivs1_1(exAt22, p.i, p.j, p.k);
const auto dexAt123 = calcderivs1_1(exAt23, p.i, p.j, p.k);
const auto dexAt133 = calcderivs1_1(exAt33, p.i, p.j, p.k);
const auto dexAt211 = calcderivs1_2(exAt11, p.i, p.j, p.k);
const auto dexAt212 = calcderivs1_2(exAt12, p.i, p.j, p.k);
const auto dexAt213 = calcderivs1_2(exAt13, p.i, p.j, p.k);
const auto dexAt222 = calcderivs1_2(exAt22, p.i, p.j, p.k);
const auto dexAt223 = calcderivs1_2(exAt23, p.i, p.j, p.k);
const auto dexAt233 = calcderivs1_2(exAt33, p.i, p.j, p.k);
const auto dexAt311 = calcderivs1_3(exAt11, p.i, p.j, p.k);
const auto dexAt312 = calcderivs1_3(exAt12, p.i, p.j, p.k);
const auto dexAt313 = calcderivs1_3(exAt13, p.i, p.j, p.k);
const auto dexAt322 = calcderivs1_3(exAt22, p.i, p.j, p.k);
const auto dexAt323 = calcderivs1_3(exAt23, p.i, p.j, p.k);
const auto dexAt333 = calcderivs1_3(exAt33, p.i, p.j, p.k);
const auto dtrGt11 = calcderivs1_1(trGt1, p.i, p.j, p.k);
const auto dtrGt12 = calcderivs1_1(trGt2, p.i, p.j, p.k);
const auto dtrGt13 = calcderivs1_1(trGt3, p.i, p.j, p.k);
const auto dtrGt21 = calcderivs1_2(trGt1, p.i, p.j, p.k);
const auto dtrGt22 = calcderivs1_2(trGt2, p.i, p.j, p.k);
const auto dtrGt23 = calcderivs1_2(trGt3, p.i, p.j, p.k);
const auto dtrGt31 = calcderivs1_3(trGt1, p.i, p.j, p.k);
const auto dtrGt32 = calcderivs1_3(trGt2, p.i, p.j, p.k);
const auto dtrGt33 = calcderivs1_3(trGt3, p.i, p.j, p.k);
const auto dTheta1 = calcderivs1_1(Theta, p.i, p.j, p.k);
const auto dTheta2 = calcderivs1_2(Theta, p.i, p.j, p.k);
const auto dTheta3 = calcderivs1_3(Theta, p.i, p.j, p.k);
const auto ddW11 = calcderivs2_11(W, p.i, p.j, p.k);
const auto ddW12 = calcderivs2_12(W, p.i, p.j, p.k);
const auto ddW13 = calcderivs2_13(W, p.i, p.j, p.k);
const auto ddW22 = calcderivs2_22(W, p.i, p.j, p.k);
const auto ddW23 = calcderivs2_23(W, p.i, p.j, p.k);
const auto ddW33 = calcderivs2_33(W, p.i, p.j, p.k);
const auto ddgamt1111 = calcderivs2_11(gamt11, p.i, p.j, p.k);
const auto ddgamt1112 = calcderivs2_11(gamt12, p.i, p.j, p.k);
const auto ddgamt1113 = calcderivs2_11(gamt13, p.i, p.j, p.k);
const auto ddgamt1122 = calcderivs2_11(gamt22, p.i, p.j, p.k);
const auto ddgamt1123 = calcderivs2_11(gamt23, p.i, p.j, p.k);
const auto ddgamt1133 = calcderivs2_11(gamt33, p.i, p.j, p.k);
const auto ddgamt1211 = calcderivs2_12(gamt11, p.i, p.j, p.k);
const auto ddgamt1212 = calcderivs2_12(gamt12, p.i, p.j, p.k);
const auto ddgamt1213 = calcderivs2_12(gamt13, p.i, p.j, p.k);
const auto ddgamt1222 = calcderivs2_12(gamt22, p.i, p.j, p.k);
const auto ddgamt1223 = calcderivs2_12(gamt23, p.i, p.j, p.k);
const auto ddgamt1233 = calcderivs2_12(gamt33, p.i, p.j, p.k);
const auto ddgamt1311 = calcderivs2_13(gamt11, p.i, p.j, p.k);
const auto ddgamt1312 = calcderivs2_13(gamt12, p.i, p.j, p.k);
const auto ddgamt1313 = calcderivs2_13(gamt13, p.i, p.j, p.k);
const auto ddgamt1322 = calcderivs2_13(gamt22, p.i, p.j, p.k);
const auto ddgamt1323 = calcderivs2_13(gamt23, p.i, p.j, p.k);
const auto ddgamt1333 = calcderivs2_13(gamt33, p.i, p.j, p.k);
const auto ddgamt2211 = calcderivs2_22(gamt11, p.i, p.j, p.k);
const auto ddgamt2212 = calcderivs2_22(gamt12, p.i, p.j, p.k);
const auto ddgamt2213 = calcderivs2_22(gamt13, p.i, p.j, p.k);
const auto ddgamt2222 = calcderivs2_22(gamt22, p.i, p.j, p.k);
const auto ddgamt2223 = calcderivs2_22(gamt23, p.i, p.j, p.k);
const auto ddgamt2233 = calcderivs2_22(gamt33, p.i, p.j, p.k);
const auto ddgamt2311 = calcderivs2_23(gamt11, p.i, p.j, p.k);
const auto ddgamt2312 = calcderivs2_23(gamt12, p.i, p.j, p.k);
const auto ddgamt2313 = calcderivs2_23(gamt13, p.i, p.j, p.k);
const auto ddgamt2322 = calcderivs2_23(gamt22, p.i, p.j, p.k);
const auto ddgamt2323 = calcderivs2_23(gamt23, p.i, p.j, p.k);
const auto ddgamt2333 = calcderivs2_23(gamt33, p.i, p.j, p.k);
const auto ddgamt3311 = calcderivs2_33(gamt11, p.i, p.j, p.k);
const auto ddgamt3312 = calcderivs2_33(gamt12, p.i, p.j, p.k);
const auto ddgamt3313 = calcderivs2_33(gamt13, p.i, p.j, p.k);
const auto ddgamt3322 = calcderivs2_33(gamt22, p.i, p.j, p.k);
const auto ddgamt3323 = calcderivs2_33(gamt23, p.i, p.j, p.k);
const auto ddgamt3333 = calcderivs2_33(gamt33, p.i, p.j, p.k);

const auto
dlnW1
=
dW1/W[ijk]
;

const auto
dlnW2
=
dW2/W[ijk]
;

const auto
dlnW3
=
dW3/W[ijk]
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
GtDDD111
=
dgamt111/2.
;

const auto
GtDDD112
=
dgamt211/2.
;

const auto
GtDDD113
=
dgamt311/2.
;

const auto
GtDDD122
=
-0.5*dgamt122 + dgamt212
;

const auto
GtDDD123
=
(-dgamt123 + dgamt213 + dgamt312)/2.
;

const auto
GtDDD133
=
-0.5*dgamt133 + dgamt313
;

const auto
GtDDD211
=
dgamt112 - dgamt211/2.
;

const auto
GtDDD212
=
dgamt122/2.
;

const auto
GtDDD213
=
(dgamt123 - dgamt213 + dgamt312)/2.
;

const auto
GtDDD222
=
dgamt222/2.
;

const auto
GtDDD223
=
dgamt322/2.
;

const auto
GtDDD233
=
-0.5*dgamt233 + dgamt323
;

const auto
GtDDD311
=
dgamt113 - dgamt311/2.
;

const auto
GtDDD312
=
(dgamt123 + dgamt213 - dgamt312)/2.
;

const auto
GtDDD313
=
dgamt133/2.
;

const auto
GtDDD322
=
dgamt223 - dgamt322/2.
;

const auto
GtDDD323
=
dgamt233/2.
;

const auto
GtDDD333
=
dgamt333/2.
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
ddW11 - dW1*Gt111 - dW2*Gt211 - dW3*Gt311
;

const auto
tDtDW12
=
ddW12 - dW1*Gt112 - dW2*Gt212 - dW3*Gt312
;

const auto
tDtDW13
=
ddW13 - dW1*Gt113 - dW2*Gt213 - dW3*Gt313
;

const auto
tDtDW22
=
ddW22 - dW1*Gt122 - dW2*Gt222 - dW3*Gt322
;

const auto
tDtDW23
=
ddW23 - dW1*Gt123 - dW2*Gt223 - dW3*Gt323
;

const auto
tDtDW33
=
ddW33 - dW1*Gt133 - dW2*Gt233 - dW3*Gt333
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
  Gt213*GtDDU213 + Gt311*GtDDU311 + Gt312*GtDDU312 + Gt313*GtDDU313 -
  (ddgamt1111*invgamt11)/2. - ddgamt1211*invgamt12 - ddgamt1311*invgamt13 -
  (ddgamt2211*invgamt22)/2. - ddgamt2311*invgamt23 -
  (ddgamt3311*invgamt33)/2. + GtDDD111*trGtd1 + GtDDD112*trGtd2 +
  GtDDD113*trGtd3 + dtrGt11*gamt11[ijk] + dtrGt12*gamt12[ijk] +
  dtrGt13*gamt13[ijk]
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
    2*Gt311*GtDDU321 + 2*Gt312*GtDDU322 + 2*Gt313*GtDDU323 -
    ddgamt1112*invgamt11 - 2*ddgamt1212*invgamt12 -
    2*ddgamt1312*invgamt13 - ddgamt2212*invgamt22 -
    2*ddgamt2312*invgamt23 - ddgamt3312*invgamt33 + GtDDD112*trGtd1 +
    GtDDD211*trGtd1 + GtDDD122*trGtd2 + GtDDD212*trGtd2 + GtDDD123*trGtd3 +
    GtDDD213*trGtd3 + dtrGt21*gamt11[ijk] +
    (dtrGt11 + dtrGt22)*gamt12[ijk] + dtrGt23*gamt13[ijk] +
    dtrGt12*gamt22[ijk] + dtrGt13*gamt23[ijk])/2.
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
    4*Gt311*GtDDU331 + 4*Gt312*GtDDU332 + 4*Gt313*GtDDU333 -
    ddgamt1113*invgamt11 - 2*ddgamt1213*invgamt12 -
    2*ddgamt1313*invgamt13 - ddgamt2213*invgamt22 -
    2*ddgamt2313*invgamt23 - ddgamt3313*invgamt33 + GtDDD113*trGtd1 +
    GtDDD311*trGtd1 + GtDDD123*trGtd2 + GtDDD312*trGtd2 + GtDDD133*trGtd3 +
    GtDDD313*trGtd3 + dtrGt31*gamt11[ijk] + dtrGt32*gamt12[ijk] +
    dtrGt11*gamt13[ijk] + dtrGt33*gamt13[ijk] + dtrGt12*gamt23[ijk] +
    dtrGt13*gamt33[ijk])/2.
;

const auto
Rt22
=
Gt112*GtDDU121 + Gt122*GtDDU122 + Gt123*GtDDU123 + 2*Gt112*GtDDU211 +
  2*Gt122*GtDDU212 + 2*Gt123*GtDDU213 + 3*Gt212*GtDDU221 +
  3*Gt222*GtDDU222 + 3*Gt223*GtDDU223 + 2*Gt312*GtDDU231 +
  2*Gt322*GtDDU232 + 2*Gt323*GtDDU233 + Gt312*GtDDU321 + Gt322*GtDDU322 +
  Gt323*GtDDU323 - (ddgamt1122*invgamt11)/2. - ddgamt1222*invgamt12 -
  ddgamt1322*invgamt13 - (ddgamt2222*invgamt22)/2. - ddgamt2322*invgamt23 -
  (ddgamt3322*invgamt33)/2. + GtDDD212*trGtd1 + GtDDD222*trGtd2 +
  GtDDD223*trGtd3 + dtrGt21*gamt12[ijk] + dtrGt22*gamt22[ijk] +
  dtrGt23*gamt23[ijk]
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
    4*Gt312*GtDDU331 + 4*Gt322*GtDDU332 + 4*Gt323*GtDDU333 -
    ddgamt1123*invgamt11 - 2*ddgamt1223*invgamt12 -
    2*ddgamt1323*invgamt13 - ddgamt2223*invgamt22 -
    2*ddgamt2323*invgamt23 - ddgamt3323*invgamt33 + GtDDD213*trGtd1 +
    GtDDD312*trGtd1 + GtDDD223*trGtd2 + GtDDD322*trGtd2 + GtDDD233*trGtd3 +
    GtDDD323*trGtd3 + dtrGt31*gamt12[ijk] + dtrGt21*gamt13[ijk] +
    dtrGt32*gamt22[ijk] + dtrGt22*gamt23[ijk] + dtrGt33*gamt23[ijk] +
    dtrGt23*gamt33[ijk])/2.
;

const auto
Rt33
=
Gt113*GtDDU131 + Gt123*GtDDU132 + Gt133*GtDDU133 + Gt213*GtDDU231 +
  Gt223*GtDDU232 + Gt233*GtDDU233 + 2*Gt113*GtDDU311 + 2*Gt123*GtDDU312 +
  2*Gt133*GtDDU313 + 2*Gt213*GtDDU321 + 2*Gt223*GtDDU322 +
  2*Gt233*GtDDU323 + 3*Gt313*GtDDU331 + 3*Gt323*GtDDU332 +
  3*Gt333*GtDDU333 - (ddgamt1133*invgamt11)/2. - ddgamt1233*invgamt12 -
  ddgamt1333*invgamt13 - (ddgamt2233*invgamt22)/2. - ddgamt2333*invgamt23 -
  (ddgamt3333*invgamt33)/2. + GtDDD313*trGtd1 + GtDDD323*trGtd2 +
  GtDDD333*trGtd3 + dtrGt31*gamt13[ijk] + dtrGt32*gamt23[ijk] +
  dtrGt33*gamt33[ijk]
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
trdexAtUU1
=
-2*dgamt111*exAtUU11*invgamt11 - dgamt211*exAtUU12*invgamt11 -
  2*dgamt113*exAtUU13*invgamt11 - dgamt311*exAtUU13*invgamt11 -
  dgamt212*exAtUU22*invgamt11 - dgamt213*exAtUU23*invgamt11 -
  dgamt312*exAtUU23*invgamt11 - dgamt313*exAtUU33*invgamt11 +
  dexAt111*Power(invgamt11,2) - dgamt211*exAtUU11*invgamt12 -
  2*dgamt122*exAtUU12*invgamt12 - 2*dgamt212*exAtUU12*invgamt12 -
  2*dgamt123*exAtUU13*invgamt12 - dgamt213*exAtUU13*invgamt12 -
  dgamt312*exAtUU13*invgamt12 - dgamt222*exAtUU22*invgamt12 -
  dgamt223*exAtUU23*invgamt12 - dgamt322*exAtUU23*invgamt12 -
  dgamt323*exAtUU33*invgamt12 + 2*dexAt112*invgamt11*invgamt12 +
  dexAt211*invgamt11*invgamt12 + dexAt122*Power(invgamt12,2) +
  dexAt212*Power(invgamt12,2) -
  2*dgamt112*(exAtUU12*invgamt11 + exAtUU11*invgamt12) -
  2*dgamt113*exAtUU11*invgamt13 - dgamt311*exAtUU11*invgamt13 -
  2*dgamt123*exAtUU12*invgamt13 - dgamt213*exAtUU12*invgamt13 -
  dgamt312*exAtUU12*invgamt13 - 2*dgamt133*exAtUU13*invgamt13 -
  2*dgamt313*exAtUU13*invgamt13 - dgamt223*exAtUU22*invgamt13 -
  dgamt233*exAtUU23*invgamt13 - dgamt323*exAtUU23*invgamt13 -
  dgamt333*exAtUU33*invgamt13 + 2*dexAt113*invgamt11*invgamt13 +
  dexAt311*invgamt11*invgamt13 + 2*dexAt123*invgamt12*invgamt13 +
  dexAt213*invgamt12*invgamt13 + dexAt312*invgamt12*invgamt13 +
  dexAt133*Power(invgamt13,2) + dexAt313*Power(invgamt13,2) -
  dgamt212*exAtUU11*invgamt22 - dgamt222*exAtUU12*invgamt22 -
  dgamt223*exAtUU13*invgamt22 + dexAt212*invgamt11*invgamt22 +
  dexAt222*invgamt12*invgamt22 + dexAt223*invgamt13*invgamt22 -
  dgamt213*exAtUU11*invgamt23 - dgamt312*exAtUU11*invgamt23 -
  dgamt223*exAtUU12*invgamt23 - dgamt322*exAtUU12*invgamt23 -
  dgamt233*exAtUU13*invgamt23 - dgamt323*exAtUU13*invgamt23 +
  dexAt213*invgamt11*invgamt23 + dexAt312*invgamt11*invgamt23 +
  dexAt223*invgamt12*invgamt23 + dexAt322*invgamt12*invgamt23 +
  dexAt233*invgamt13*invgamt23 + dexAt323*invgamt13*invgamt23 -
  dgamt313*exAtUU11*invgamt33 - dgamt323*exAtUU12*invgamt33 -
  dgamt333*exAtUU13*invgamt33 + dexAt313*invgamt11*invgamt33 +
  dexAt323*invgamt12*invgamt33 + dexAt333*invgamt13*invgamt33
;

const auto
trdexAtUU2
=
-(dgamt113*exAtUU23*invgamt11) - 2*dgamt211*exAtUU12*invgamt12 -
  dgamt113*exAtUU13*invgamt12 - dgamt311*exAtUU13*invgamt12 -
  dgamt122*exAtUU22*invgamt12 - 2*dgamt212*exAtUU22*invgamt12 -
  dgamt123*exAtUU23*invgamt12 - 2*dgamt213*exAtUU23*invgamt12 -
  dgamt312*exAtUU23*invgamt12 - dgamt313*exAtUU33*invgamt12 +
  dexAt111*invgamt11*invgamt12 + dexAt112*Power(invgamt12,2) +
  dexAt211*Power(invgamt12,2) -
  dgamt111*(exAtUU12*invgamt11 + exAtUU11*invgamt12) -
  dgamt113*exAtUU12*invgamt13 - dgamt311*exAtUU12*invgamt13 -
  dgamt123*exAtUU22*invgamt13 - dgamt312*exAtUU22*invgamt13 -
  dgamt133*exAtUU23*invgamt13 - dgamt313*exAtUU23*invgamt13 +
  dexAt113*invgamt12*invgamt13 + dexAt311*invgamt12*invgamt13 -
  dgamt122*exAtUU12*invgamt22 - 2*dgamt212*exAtUU12*invgamt22 -
  dgamt123*exAtUU13*invgamt22 - dgamt312*exAtUU13*invgamt22 -
  2*dgamt222*exAtUU22*invgamt22 - 2*dgamt223*exAtUU23*invgamt22 -
  dgamt322*exAtUU23*invgamt22 - dgamt323*exAtUU33*invgamt22 +
  dexAt112*invgamt11*invgamt22 + dexAt122*invgamt12*invgamt22 +
  2*dexAt212*invgamt12*invgamt22 + dexAt123*invgamt13*invgamt22 +
  dexAt312*invgamt13*invgamt22 + dexAt222*Power(invgamt22,2) -
  dgamt112*(exAtUU22*invgamt11 + 2*exAtUU12*invgamt12 +
     exAtUU11*invgamt22) - dgamt113*exAtUU11*invgamt23 -
  dgamt123*exAtUU12*invgamt23 - 2*dgamt213*exAtUU12*invgamt23 -
  dgamt312*exAtUU12*invgamt23 - dgamt133*exAtUU13*invgamt23 -
  dgamt313*exAtUU13*invgamt23 - 2*dgamt223*exAtUU22*invgamt23 -
  dgamt322*exAtUU22*invgamt23 - 2*dgamt233*exAtUU23*invgamt23 -
  2*dgamt323*exAtUU23*invgamt23 - dgamt333*exAtUU33*invgamt23 +
  dexAt113*invgamt11*invgamt23 + dexAt123*invgamt12*invgamt23 +
  2*dexAt213*invgamt12*invgamt23 + dexAt312*invgamt12*invgamt23 +
  dexAt133*invgamt13*invgamt23 + dexAt313*invgamt13*invgamt23 +
  2*dexAt223*invgamt22*invgamt23 + dexAt322*invgamt22*invgamt23 +
  dexAt233*Power(invgamt23,2) + dexAt323*Power(invgamt23,2) -
  dgamt313*exAtUU12*invgamt33 - dgamt323*exAtUU22*invgamt33 -
  dgamt333*exAtUU23*invgamt33 + dexAt313*invgamt12*invgamt33 +
  dexAt323*invgamt22*invgamt33 + dexAt333*invgamt23*invgamt33
;

const auto
trdexAtUU3
=
-(dgamt113*exAtUU33*invgamt11) - dgamt211*exAtUU13*invgamt12 -
  dgamt122*exAtUU23*invgamt12 - dgamt212*exAtUU23*invgamt12 -
  dgamt123*exAtUU33*invgamt12 - dgamt213*exAtUU33*invgamt12 -
  dgamt211*exAtUU12*invgamt13 - 2*dgamt113*exAtUU13*invgamt13 -
  2*dgamt311*exAtUU13*invgamt13 - dgamt212*exAtUU22*invgamt13 -
  dgamt123*exAtUU23*invgamt13 - dgamt213*exAtUU23*invgamt13 -
  2*dgamt312*exAtUU23*invgamt13 - dgamt133*exAtUU33*invgamt13 -
  2*dgamt313*exAtUU33*invgamt13 + dexAt111*invgamt11*invgamt13 +
  dexAt112*invgamt12*invgamt13 + dexAt211*invgamt12*invgamt13 +
  dexAt113*Power(invgamt13,2) + dexAt311*Power(invgamt13,2) -
  dgamt111*(exAtUU13*invgamt11 + exAtUU11*invgamt13) -
  dgamt212*exAtUU13*invgamt22 - dgamt222*exAtUU23*invgamt22 -
  dgamt223*exAtUU33*invgamt22 + dexAt212*invgamt13*invgamt22 -
  dgamt122*exAtUU12*invgamt23 - dgamt212*exAtUU12*invgamt23 -
  dgamt123*exAtUU13*invgamt23 - dgamt213*exAtUU13*invgamt23 -
  2*dgamt312*exAtUU13*invgamt23 - dgamt222*exAtUU22*invgamt23 -
  2*dgamt223*exAtUU23*invgamt23 - 2*dgamt322*exAtUU23*invgamt23 -
  dgamt233*exAtUU33*invgamt23 - 2*dgamt323*exAtUU33*invgamt23 +
  dexAt112*invgamt11*invgamt23 + dexAt122*invgamt12*invgamt23 +
  dexAt212*invgamt12*invgamt23 + dexAt123*invgamt13*invgamt23 +
  dexAt213*invgamt13*invgamt23 + 2*dexAt312*invgamt13*invgamt23 +
  dexAt222*invgamt22*invgamt23 + dexAt223*Power(invgamt23,2) +
  dexAt322*Power(invgamt23,2) -
  dgamt112*(exAtUU23*invgamt11 + exAtUU13*invgamt12 + exAtUU12*invgamt13 +
     exAtUU11*invgamt23) - dgamt113*exAtUU11*invgamt33 -
  dgamt123*exAtUU12*invgamt33 - dgamt213*exAtUU12*invgamt33 -
  dgamt133*exAtUU13*invgamt33 - 2*dgamt313*exAtUU13*invgamt33 -
  dgamt223*exAtUU22*invgamt33 - dgamt233*exAtUU23*invgamt33 -
  2*dgamt323*exAtUU23*invgamt33 - 2*dgamt333*exAtUU33*invgamt33 +
  dexAt113*invgamt11*invgamt33 + dexAt123*invgamt12*invgamt33 +
  dexAt213*invgamt12*invgamt33 + dexAt133*invgamt13*invgamt33 +
  2*dexAt313*invgamt13*invgamt33 + dexAt223*invgamt22*invgamt33 +
  dexAt233*invgamt23*invgamt33 + 2*dexAt323*invgamt23*invgamt33 +
  dexAt333*Power(invgamt33,2)
;


ZtC1[ijk]
=
(-trGtd1 + trGt1[ijk])/2.
;

ZtC2[ijk]
=
(-trGtd2 + trGt2[ijk])/2.
;

ZtC3[ijk]
=
(-trGtd3 + trGt3[ijk])/2.
;

HC[ijk]
=
(-48*cpi*rho + 3*trR - 3*exAtUU11*exAt11[ijk] - 6*exAtUU12*exAt12[ijk] -
    6*exAtUU13*exAt13[ijk] - 3*exAtUU22*exAt22[ijk] -
    6*exAtUU23*exAt23[ijk] - 3*exAtUU33*exAt33[ijk] +
    2*Power(exKh[ijk],2) + 8*exKh[ijk]*Theta[ijk] + 8*Power(Theta[ijk],2))/3.
;

MtC1[ijk]
=
(-9*dlnW1*exAtUU11 - 9*dlnW2*exAtUU12 - 9*dlnW3*exAtUU13 +
    3*exAtUU11*Gt111 + 6*exAtUU12*Gt112 + 6*exAtUU13*Gt113 +
    3*exAtUU22*Gt122 + 6*exAtUU23*Gt123 + 3*exAtUU33*Gt133 -
    2*dexKh1*invgamt11 - 4*dTheta1*invgamt11 - 2*dexKh2*invgamt12 -
    4*dTheta2*invgamt12 - 2*dexKh3*invgamt13 - 4*dTheta3*invgamt13 -
    24*cpi*invgamt11*Sm1 - 24*cpi*invgamt12*Sm2 - 24*cpi*invgamt13*Sm3 +
    3*trdexAtUU1)/3.
;

MtC2[ijk]
=
(-9*dlnW1*exAtUU12 - 9*dlnW2*exAtUU22 - 9*dlnW3*exAtUU23 +
    3*exAtUU11*Gt211 + 6*exAtUU12*Gt212 + 6*exAtUU13*Gt213 +
    3*exAtUU22*Gt222 + 6*exAtUU23*Gt223 + 3*exAtUU33*Gt233 -
    2*dexKh1*invgamt12 - 4*dTheta1*invgamt12 - 2*dexKh2*invgamt22 -
    4*dTheta2*invgamt22 - 2*dexKh3*invgamt23 - 4*dTheta3*invgamt23 -
    24*cpi*invgamt12*Sm1 - 24*cpi*invgamt22*Sm2 - 24*cpi*invgamt23*Sm3 +
    3*trdexAtUU2)/3.
;

MtC3[ijk]
=
(-9*dlnW1*exAtUU13 - 9*dlnW2*exAtUU23 - 9*dlnW3*exAtUU33 +
    3*exAtUU11*Gt311 + 6*exAtUU12*Gt312 + 6*exAtUU13*Gt313 +
    3*exAtUU22*Gt322 + 6*exAtUU23*Gt323 + 3*exAtUU33*Gt333 -
    2*dexKh1*invgamt13 - 4*dTheta1*invgamt13 - 2*dexKh2*invgamt23 -
    4*dTheta2*invgamt23 - 2*dexKh3*invgamt33 - 4*dTheta3*invgamt33 -
    24*cpi*invgamt13*Sm1 - 24*cpi*invgamt23*Sm2 - 24*cpi*invgamt33*Sm3 +
    3*trdexAtUU3)/3.
;


});
});

#endif // #ifndef Z4COWGPU_SET_CONSTRAINT_HXX

/* Z4cowGPU_set_constraint.hxx */
