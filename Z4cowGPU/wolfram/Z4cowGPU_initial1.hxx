/* Z4cowGPU_initial1.hxx */
/* Produced with Generato */

#ifndef Z4COWGPU_INITIAL1_HXX
#define Z4COWGPU_INITIAL1_HXX

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
const auto alpha = gf_alpha;
const auto beta1 = gf_beta[0];
const auto beta2 = gf_beta[1];
const auto beta3 = gf_beta[2];
const auto ADMgam11 = gf_ADMgam[0];
const auto ADMgam12 = gf_ADMgam[1];
const auto ADMgam13 = gf_ADMgam[2];
const auto ADMgam22 = gf_ADMgam[3];
const auto ADMgam23 = gf_ADMgam[4];
const auto ADMgam33 = gf_ADMgam[5];
const auto ADMK11 = gf_ADMK[0];
const auto ADMK12 = gf_ADMK[1];
const auto ADMK13 = gf_ADMK[2];
const auto ADMK22 = gf_ADMK[3];
const auto ADMK23 = gf_ADMK[4];
const auto ADMK33 = gf_ADMK[5];
const auto ADMalpha = gf_ADMalpha;
const auto ADMbeta1 = gf_ADMbeta[0];
const auto ADMbeta2 = gf_ADMbeta[1];
const auto ADMbeta3 = gf_ADMbeta[2];

noinline([&]() __attribute__((__flatten__, __hot__)) {
grid.loop_int_device<0, 0, 0>(
  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
const int ijk = layout2.linear(p.i, p.j, p.k);

const auto
invdetgamma
=
-(1/(Power(ADMgam13[ijk],2)*ADMgam22[ijk] -
      2*ADMgam12[ijk]*ADMgam13[ijk]*ADMgam23[ijk] +
      Power(ADMgam12[ijk],2)*ADMgam33[ijk] +
      ADMgam11[ijk]*(Power(ADMgam23[ijk],2) - ADMgam22[ijk]*ADMgam33[ijk])))
;

const auto
invgamma11
=
invdetgamma*(-Power(ADMgam23[ijk],2) + ADMgam22[ijk]*ADMgam33[ijk])
;

const auto
invgamma12
=
invdetgamma*(ADMgam13[ijk]*ADMgam23[ijk] - ADMgam12[ijk]*ADMgam33[ijk])
;

const auto
invgamma13
=
invdetgamma*(-(ADMgam13[ijk]*ADMgam22[ijk]) + ADMgam12[ijk]*ADMgam23[ijk])
;

const auto
invgamma22
=
invdetgamma*(-Power(ADMgam13[ijk],2) + ADMgam11[ijk]*ADMgam33[ijk])
;

const auto
invgamma23
=
invdetgamma*(ADMgam12[ijk]*ADMgam13[ijk] - ADMgam11[ijk]*ADMgam23[ijk])
;

const auto
invgamma33
=
invdetgamma*(-Power(ADMgam12[ijk],2) + ADMgam11[ijk]*ADMgam22[ijk])
;

const auto
trexK
=
invgamma11*ADMK11[ijk] + 2*invgamma12*ADMK12[ijk] +
  2*invgamma13*ADMK13[ijk] + invgamma22*ADMK22[ijk] +
  2*invgamma23*ADMK23[ijk] + invgamma33*ADMK33[ijk]
;

const auto
tempTheta
=
0
;

const auto tempW = cbrt(sqrt(invdetgamma));

W[ijk]
=
tempW
;

gamt11[ijk]
=
Power(tempW,2)*ADMgam11[ijk]
;

gamt12[ijk]
=
Power(tempW,2)*ADMgam12[ijk]
;

gamt13[ijk]
=
Power(tempW,2)*ADMgam13[ijk]
;

gamt22[ijk]
=
Power(tempW,2)*ADMgam22[ijk]
;

gamt23[ijk]
=
Power(tempW,2)*ADMgam23[ijk]
;

gamt33[ijk]
=
Power(tempW,2)*ADMgam33[ijk]
;

exKh[ijk]
=
-2*tempTheta + trexK
;

exAt11[ijk]
=
Power(tempW,2)*(-0.3333333333333333*(trexK*ADMgam11[ijk]) + ADMK11[ijk])
;

exAt12[ijk]
=
Power(tempW,2)*(-0.3333333333333333*(trexK*ADMgam12[ijk]) + ADMK12[ijk])
;

exAt13[ijk]
=
Power(tempW,2)*(-0.3333333333333333*(trexK*ADMgam13[ijk]) + ADMK13[ijk])
;

exAt22[ijk]
=
Power(tempW,2)*(-0.3333333333333333*(trexK*ADMgam22[ijk]) + ADMK22[ijk])
;

exAt23[ijk]
=
Power(tempW,2)*(-0.3333333333333333*(trexK*ADMgam23[ijk]) + ADMK23[ijk])
;

exAt33[ijk]
=
Power(tempW,2)*(-0.3333333333333333*(trexK*ADMgam33[ijk]) + ADMK33[ijk])
;

Theta[ijk]
=
tempTheta
;

alpha[ijk]
=
ADMalpha[ijk]
;

beta1[ijk]
=
ADMbeta1[ijk]
;

beta2[ijk]
=
ADMbeta2[ijk]
;

beta3[ijk]
=
ADMbeta3[ijk]
;


});
});

#endif // #ifndef Z4COWGPU_INITIAL1_HXX

/* Z4cowGPU_initial1.hxx */
