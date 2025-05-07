
#ifndef _ARRHENIUS_DETONATION_H_INCLUDED_
#define _ARRHENIUS_DETONATION_H_INCLUDED_

#include "resources/Scripts/cfd/species/species_auxiliary.hcl"
#include "resources/Scripts/cfd/reaction/reaction_generic.hcl"

#define E_ch 7000.0f
#define K_ch 1000000.0f

inline float arrhenius_detonation(float zeta, float T) {

  return K_ch * (1.0f - zeta) * exp(-E_ch / T);
  
}

float zeta_dot_calc_arrhenius(float z, float T, float y_H2, float MMix){
    
    //float MMix;
    //MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);
    //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
    
    // calculatate molar fraction
    // of mixture
    float zx = MMix / Mis[0] * z;
    
    float zeta;

    //check that there is mroe than 4%
    if(zx>4.e-2f){    
        zeta = give_zeta(z, y_H2);
        return arrhenius_detonation(zeta, T);
    }else{
//        zeta = 1.0;
        return 0.0f;
    }

}


void w_rhos_arrhenius_det(float z, float T, float y_H2, float y_O2, float y_N2, float y_H2O, 
            __global float* w_rho_H2, __global float* w_rho_O2, 
            __global float* w_rho_N2, __global float* w_rho_H2O, 
            __global float* deintdt, __global float* zeta_dot){

    float MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);

    //float inv_MMix = 1.0f / MMix;

    *zeta_dot = zeta_dot_calc_arrhenius(z, T, y_H2, MMix);
    
    w_from_zeta_dot(w_rho_H2, w_rho_O2, w_rho_N2, w_rho_H2O, deintdt, zeta_dot, MMix);

    return;
}


#endif    // _ARRHENIUS_DETONATION_H_INCLUDED_