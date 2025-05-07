
#ifndef _ARRHENIUS_DETONATION_H_INCLUDED_
#define _ARRHENIUS_DETONATION_H_INCLUDED_

#include "resources/Scripts/cfd/species/species_auxiliary.hcl"
#include "resources/Scripts/cfd/reaction/reaction_generic.hcl"


__constant float S_L = 3.0f;
__constant float sigma = 5.0f;

inline float zimont_deflagration(float mod_grad_zeta) {

  return sigma * S_L * mod_grad_zeta;
  
}

float zeta_dot_calc_zimont(float z, float MMix, float mod_grad_zeta){
    
    //float MMix;
    //MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);
    //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
    
    // calculatate molar fraction
    // of mixture
    float zx = MMix / Mis[0] * z;
    
    //check that there is mroe than 4%
    if(zx>4.e-2f){            
        return zimont_deflagration(mod_grad_zeta)
    }else{
//        zeta = 1.0;
        return 0.0f;
    }

}


void w_rhos_arrhenius_det(float z, float T, 
    float y_H2, float y_O2, float y_N2, float y_H2O, 
    __global float* w_rho_H2, __global float* w_rho_O2, 
    __global float* w_rho_N2, __global float* w_rho_H2O, 
    __global float* deintdt, __global float* zeta_dot,
    vec* grad_zeta){

    float MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);

    float mod_grad_zeta = sqrt(dot(grad_zeta, grad_zeta));
    //float inv_MMix = 1.0f / MMix;

    *zeta_dot = zeta_dot_calc(z, T, y_H2, MMix, mod_grad_zeta);
    
    w_from_zeta_dot(w_rho_H2, w_rho_O2, w_rho_N2, w_rho_H2O, deintdt, zeta_dot, MMix);

    //*zeta_dot = zeta_dot_val;

/*    float inv_M_times_zeta_dot = *zeta_dot / MMix;
 
    //const float nus[4] = {-1.0f, -0.5f, 0.0f, 1.0f};
    //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
    //const float h_h20 = -285.83e3f/0.018f;

    *w_rho_H2 += nus[0] * Mis[0] * inv_M_times_zeta_dot;
    *w_rho_O2 += nus[1] * Mis[1] * inv_M_times_zeta_dot;
    *w_rho_N2 += nus[2] * Mis[2] * inv_M_times_zeta_dot;
    *w_rho_H2O += nus[3] * Mis[3] * inv_M_times_zeta_dot;
    *deintdt -= h_h20 * *w_rho_H2O;
    
    //printf("%g\n",*w_rho_H2O);
*/

    return;
}


#endif    // _ARRHENIUS_DETONATION_H_INCLUDED_