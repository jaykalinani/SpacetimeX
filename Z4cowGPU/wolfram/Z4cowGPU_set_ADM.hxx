/* Z4cowGPU_set_ADM.hxx */
/* Produced with Generato */

#ifndef Z4COWGPU_SET_ADM_HXX
#define Z4COWGPU_SET_ADM_HXX

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

noinline([&]() __attribute__((__flatten__, __hot__)) {
grid.loop_all_device<0, 0, 0>(
  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
const int ijk = layout2.linear(p.i, p.j, p.k);

ADMgam11[ijk]
=
gamt11[ijk]/Power(W[ijk],2)
;

ADMgam12[ijk]
=
gamt12[ijk]/Power(W[ijk],2)
;

ADMgam13[ijk]
=
gamt13[ijk]/Power(W[ijk],2)
;

ADMgam22[ijk]
=
gamt22[ijk]/Power(W[ijk],2)
;

ADMgam23[ijk]
=
gamt23[ijk]/Power(W[ijk],2)
;

ADMgam33[ijk]
=
gamt33[ijk]/Power(W[ijk],2)
;

ADMK11[ijk]
=
(exAt11[ijk] + (gamt11[ijk]*(exKh[ijk] + 2*Theta[ijk]))/3.)/Power(W[ijk],2)
;

ADMK12[ijk]
=
(exAt12[ijk] + (gamt12[ijk]*(exKh[ijk] + 2*Theta[ijk]))/3.)/Power(W[ijk],2)
;

ADMK13[ijk]
=
(exAt13[ijk] + (gamt13[ijk]*(exKh[ijk] + 2*Theta[ijk]))/3.)/Power(W[ijk],2)
;

ADMK22[ijk]
=
(exAt22[ijk] + (gamt22[ijk]*(exKh[ijk] + 2*Theta[ijk]))/3.)/Power(W[ijk],2)
;

ADMK23[ijk]
=
(exAt23[ijk] + (gamt23[ijk]*(exKh[ijk] + 2*Theta[ijk]))/3.)/Power(W[ijk],2)
;

ADMK33[ijk]
=
(exAt33[ijk] + (gamt33[ijk]*(exKh[ijk] + 2*Theta[ijk]))/3.)/Power(W[ijk],2)
;

ADMalpha[ijk]
=
alpha[ijk]
;

ADMbeta1[ijk]
=
beta1[ijk]
;

ADMbeta2[ijk]
=
beta2[ijk]
;

ADMbeta3[ijk]
=
beta3[ijk]
;


});
});

#endif // #ifndef Z4COWGPU_SET_ADM_HXX

/* Z4cowGPU_set_ADM.hxx */
