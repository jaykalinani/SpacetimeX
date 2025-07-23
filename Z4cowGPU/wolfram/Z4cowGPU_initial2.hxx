/* Z4cowGPU_initial2.hxx */
/* Produced with Generato */

#ifndef Z4COWGPU_INITIAL2_HXX
#define Z4COWGPU_INITIAL2_HXX

const auto trGt1 = gf_trGt[0];
const auto trGt2 = gf_trGt[1];
const auto trGt3 = gf_trGt[2];
const auto gamt11 = gf_gamt[0];
const auto gamt12 = gf_gamt[1];
const auto gamt13 = gf_gamt[2];
const auto gamt22 = gf_gamt[3];
const auto gamt23 = gf_gamt[4];
const auto gamt33 = gf_gamt[5];

noinline([&]() __attribute__((__flatten__, __hot__)) {
grid.loop_int_device<0, 0, 0>(
  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
const int ijk = layout2.linear(p.i, p.j, p.k);

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


trGt1[ijk]
=
GtDDD111*Power(invgamt11,2) + 2*GtDDD112*invgamt11*invgamt12 +
  GtDDD211*invgamt11*invgamt12 + 2*GtDDD212*Power(invgamt12,2) +
  2*GtDDD113*invgamt11*invgamt13 + GtDDD311*invgamt11*invgamt13 +
  2*GtDDD213*invgamt12*invgamt13 + 2*GtDDD312*invgamt12*invgamt13 +
  2*GtDDD313*Power(invgamt13,2) + GtDDD122*invgamt11*invgamt22 +
  GtDDD222*invgamt12*invgamt22 + GtDDD322*invgamt13*invgamt22 +
  2*GtDDD123*invgamt11*invgamt23 + 2*GtDDD223*invgamt12*invgamt23 +
  2*GtDDD323*invgamt13*invgamt23 + GtDDD133*invgamt11*invgamt33 +
  GtDDD233*invgamt12*invgamt33 + GtDDD333*invgamt13*invgamt33
;

trGt2[ijk]
=
GtDDD111*invgamt11*invgamt12 + 2*GtDDD112*Power(invgamt12,2) +
  2*GtDDD113*invgamt12*invgamt13 + GtDDD211*invgamt11*invgamt22 +
  GtDDD122*invgamt12*invgamt22 + 2*GtDDD212*invgamt12*invgamt22 +
  2*GtDDD213*invgamt13*invgamt22 + GtDDD222*Power(invgamt22,2) +
  GtDDD311*invgamt11*invgamt23 + 2*GtDDD123*invgamt12*invgamt23 +
  2*GtDDD312*invgamt12*invgamt23 + 2*GtDDD313*invgamt13*invgamt23 +
  2*GtDDD223*invgamt22*invgamt23 + GtDDD322*invgamt22*invgamt23 +
  2*GtDDD323*Power(invgamt23,2) + GtDDD133*invgamt12*invgamt33 +
  GtDDD233*invgamt22*invgamt33 + GtDDD333*invgamt23*invgamt33
;

trGt3[ijk]
=
GtDDD111*invgamt11*invgamt13 + 2*GtDDD112*invgamt12*invgamt13 +
  2*GtDDD113*Power(invgamt13,2) + GtDDD122*invgamt13*invgamt22 +
  GtDDD211*invgamt11*invgamt23 + 2*GtDDD212*invgamt12*invgamt23 +
  2*GtDDD123*invgamt13*invgamt23 + 2*GtDDD213*invgamt13*invgamt23 +
  GtDDD222*invgamt22*invgamt23 + 2*GtDDD223*Power(invgamt23,2) +
  GtDDD311*invgamt11*invgamt33 + 2*GtDDD312*invgamt12*invgamt33 +
  GtDDD133*invgamt13*invgamt33 + 2*GtDDD313*invgamt13*invgamt33 +
  GtDDD322*invgamt22*invgamt33 + GtDDD233*invgamt23*invgamt33 +
  2*GtDDD323*invgamt23*invgamt33 + GtDDD333*Power(invgamt33,2)
;


});
});

#endif // #ifndef Z4COWGPU_INITIAL2_HXX

/* Z4cowGPU_initial2.hxx */
