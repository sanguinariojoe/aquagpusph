
#ifndef _REACTION_GENERIC_H_INCLUDED_
#define _REACTION_GENERIC_H_INCLUDED_

#ifndef SPECIES_HEADER
#error "working with species requires to load a backend module"
#endif
#include SPECIES_HEADER
//#include "resources/Scripts/cfd/species/species_auxiliary.hcl"
//#include "resources/Scripts/cfd/reaction/arrhenius_detonation.hcl"

//__constant float h_h20 = -285.83e3f/0.018f;

/*
float arrhenius_detonation(float zeta, float T);
float calc_cp_mix(float y_H2, float y_O2, float y_N2, float y_H2O);
void calc_gamma_cp_cv(float y_H2, float y_O2, float y_N2, float y_H2O, float* gamma, float* cp, float* cv);
float molar_mass_mixture(float y_H2, float y_O2, float y_N2, float y_H2O);
float M_zeta_dot(float z, float T, float y_H2, float y_O2, float y_N2, float y_H2O);
void w_rhos(float z, float T, float y_H2, float y_O2, float y_N2, float y_H2O, 
            float* w_rho_H2, float* w_rho_O2, float* w_rho_N2, float* w_rho_H2O, float* deintdt);
*/

inline float give_zeta(float z, float y_0){

    // zeta prograss variable of combustion
    return (z - y_0) / (z+1.0e-12f);
}

void w_from_zeta_dot(__global species_t* w_rhos,
    /*__global float* w_rho_H2, 
    __global float* w_rho_O2, 
    __global float* w_rho_N2, 
    __global float* w_rho_H2O, */
    __global float* deintdt, 
    __global float* zeta_dot, 
    float MMix){

        float inv_M_times_zeta_dot = *zeta_dot / MMix;
 
        //const float nus[4] = {-1.0f, -0.5f, 0.0f, 1.0f};
        //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
        //const float h_h20 = -285.83e3f/0.018f;
        
        const species_t Mis = MIS;
        const species_t nus = NUS;

        *w_rhos += nus * Mis * inv_M_times_zeta_dot;

        /**w_rho_H2 += nus[0] * Mis[0] * inv_M_times_zeta_dot;
        *w_rho_O2 += nus[1] * Mis[1] * inv_M_times_zeta_dot;
        *w_rho_N2 += nus[2] * Mis[2] * inv_M_times_zeta_dot;
        *w_rho_H2O += nus[3] * Mis[3] * inv_M_times_zeta_dot;*/
        
        const species_t hplus_mass = HPLUS_MASS; 

        *deintdt -= dot(hplus_mass, *w_rhos);
        //*deintdt -= h_h20 * (*w_rho).H2O;
        
        //printf("%g\n",*w_rho_H2O);
        
        return;
    }

#endif    // _REACTION_GENERIC_H_INCLUDED_