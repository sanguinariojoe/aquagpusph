
#ifndef _ARRHENIUS_DETONATION_H_INCLUDED_
#define _ARRHENIUS_DETONATION_H_INCLUDED_

#ifndef SPECIES_HEADER
#error "species.xml module requires to load a backend module"
#endif
#include SPECIES_HEADER
//#include "resources/Scripts/cfd/species/species_auxiliary.hcl"
#include "resources/Scripts/cfd/reactive/reaction_generic.hcl"

#define E_ch 7000.0f
#define K_ch 1.0e7f

inline float
arrhenius_detonation(float zeta, float T)
{

	return K_ch * (1.0f - zeta) * exp(-E_ch / T);
}

#define LOW_CONC 4.0e-2f
inline float
zeta_dot_calc_arrhenius(float z, float T, float y_0, float MMix)
{

	// float MMix;
	// MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);
	// const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};

	// calculatate molar fraction
	// of mixture

    const species_t Mis = MIS;
	float zx = MMix / Mis.SPECIES_COMPONENT0 * z;

	float zeta;

	// check that there is mroe than 4%
	if (zx > LOW_CONC) {
		zeta = give_zeta(z, y_0);
		return arrhenius_detonation(zeta, T);
	} else {
		//        zeta = 1.0;
		return 0.0f;
	}
}

void
w_rhos_arrhenius_det(
    float z,
    float T, 
    species_t ys,
    __global species_t* w_rhos,
    __global float* deintdt,
    __global float* zeta_dot)
{

    float MMix = molar_mass_mixture(ys);

	// float inv_MMix = 1.0f / MMix;

	*zeta_dot = zeta_dot_calc_arrhenius(z, T, ys.SPECIES_COMPONENT0, MMix);

	w_from_zeta_dot(w_rhos, deintdt, zeta_dot, MMix);

	return;
}

#endif // _ARRHENIUS_DETONATION_H_INCLUDED_