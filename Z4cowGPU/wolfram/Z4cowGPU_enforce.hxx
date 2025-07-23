/* Z4cowGPU_enforce.hxx */
/* Produced with Generato */

#ifndef Z4COWGPU_ENFORCE_HXX
#define Z4COWGPU_ENFORCE_HXX

const auto gamt11 = gf_gamt[0];
const auto gamt12 = gf_gamt[1];
const auto gamt13 = gf_gamt[2];
const auto gamt22 = gf_gamt[3];
const auto gamt23 = gf_gamt[4];
const auto gamt33 = gf_gamt[5];
const auto exAt11 = gf_exAt[0];
const auto exAt12 = gf_exAt[1];
const auto exAt13 = gf_exAt[2];
const auto exAt22 = gf_exAt[3];
const auto exAt23 = gf_exAt[4];
const auto exAt33 = gf_exAt[5];
const auto alpha = gf_alpha;

noinline([&]() __attribute__((__flatten__, __hot__)) {
grid.loop_all_device<0, 0, 0>(
  grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
const int ijk = layout2.linear(p.i, p.j, p.k);

const auto
Wold
=
W[ijk]
;

const auto
gamtold11
=
gamt11[ijk]
;

const auto
gamtold12
=
gamt12[ijk]
;

const auto
gamtold13
=
gamt13[ijk]
;

const auto
gamtold22
=
gamt22[ijk]
;

const auto
gamtold23
=
gamt23[ijk]
;

const auto
gamtold33
=
gamt33[ijk]
;

const auto
exAtold11
=
exAt11[ijk]
;

const auto
exAtold12
=
exAt12[ijk]
;

const auto
exAtold13
=
exAt13[ijk]
;

const auto
exAtold22
=
exAt22[ijk]
;

const auto
exAtold23
=
exAt23[ijk]
;

const auto
exAtold33
=
exAt33[ijk]
;

const auto
alphaold
=
alpha[ijk]
;

W[ijk] = fmax(cWfloor, Wold);
alpha[ijk] = fmax(calphafloor, alphaold);

const auto
invdetgamtold
=
1/(-(Power(gamtold13,2)*gamtold22) + 2*gamtold12*gamtold13*gamtold23 -
    gamtold11*Power(gamtold23,2) - Power(gamtold12,2)*gamtold33 +
    gamtold11*gamtold22*gamtold33)
;

const auto Wfac = cbrt(sqrt(invdetgamtold));

gamt11[ijk]
=
gamtold11*Power(Wfac,2)
;

gamt12[ijk]
=
gamtold12*Power(Wfac,2)
;

gamt13[ijk]
=
gamtold13*Power(Wfac,2)
;

gamt22[ijk]
=
gamtold22*Power(Wfac,2)
;

gamt23[ijk]
=
gamtold23*Power(Wfac,2)
;

gamt33[ijk]
=
gamtold33*Power(Wfac,2)
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
trexAtold
=
exAtold11*invgamt11 + 2*exAtold12*invgamt12 + 2*exAtold13*invgamt13 +
  exAtold22*invgamt22 + 2*exAtold23*invgamt23 + exAtold33*invgamt33
;


exAt11[ijk]
=
exAtold11 - (trexAtold*gamt11[ijk])/3.
;

exAt12[ijk]
=
exAtold12 - (trexAtold*gamt12[ijk])/3.
;

exAt13[ijk]
=
exAtold13 - (trexAtold*gamt13[ijk])/3.
;

exAt22[ijk]
=
exAtold22 - (trexAtold*gamt22[ijk])/3.
;

exAt23[ijk]
=
exAtold23 - (trexAtold*gamt23[ijk])/3.
;

exAt33[ijk]
=
exAtold33 - (trexAtold*gamt33[ijk])/3.
;


});
});

#endif // #ifndef Z4COWGPU_ENFORCE_HXX

/* Z4cowGPU_enforce.hxx */
