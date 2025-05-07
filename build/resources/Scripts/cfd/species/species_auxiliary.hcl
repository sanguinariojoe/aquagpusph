
#ifndef _SPECIES_AUXILIARY_H_INCLUDED_
#define _SPECIES_AUXILIARY_H_INCLUDED_

#define R_gas 8.31f

__constant float cps[4] = {14200.0f, 918.0f, 1040.0f, 2050.0f};
__constant float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
__constant float nus[4] = {-1.0f, -0.5f, 0.0f, 1.0f};

/*
float arrhenius_detonation(float zeta, float T);
float calc_cp_mix(float y_H2, float y_O2, float y_N2, float y_H2O);
void calc_gamma_cp_cv(float y_H2, float y_O2, float y_N2, float y_H2O, float* gamma, float* cp, float* cv);
float molar_mass_mixture(float y_H2, float y_O2, float y_N2, float y_H2O);
float M_zeta_dot(float z, float T, float y_H2, float y_O2, float y_N2, float y_H2O);
void w_rhos(float z, float T, float y_H2, float y_O2, float y_N2, float y_H2O, 
            float* w_rho_H2, float* w_rho_O2, float* w_rho_N2, float* w_rho_H2O, float* deintdt);
*/

float molar_mass_mixture(float y_H2, float y_O2, float y_N2, float y_H2O){
    //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
    return 1.0f /(y_H2 / Mis[0] + y_O2 / Mis[1] + y_N2 / Mis[2] + y_H2O / Mis[3]);

}

float calc_cp_mix(float y_H2, float y_O2, float y_N2, float y_H2O){

    //const float cps[4] = {14200.0f, 913.0f, 1040.0f, 2050.0f};
    return y_H2 * cps[0] + y_O2 * cps[1] + y_N2 * cps[2] + y_H2O * cps[3];

}

void X_from_Y(float y_H2, float y_O2, float y_N2, float y_H2O, __global float* x_H2, __global float* x_O2, __global float* x_N2, __global float* x_H2O){
    float MMix;
    //const float Mis[4] = {0.002f, 0.032f, 0.028f, 0.018f};
    MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);

    *x_H2 = MMix / Mis[0] * y_H2;
    *x_O2 = MMix / Mis[1] * y_O2;
    *x_N2 = MMix / Mis[2] * y_N2;
    *x_H2O = MMix / Mis[3] * y_H2O;
    
    //printf("x_O2 = %f   x_N2 = %f\n", *x_O2, *x_N2);
    return;
}

void calc_gamma_cv(float y_H2, float y_O2, float y_N2, float y_H2O, __global float* gamma, __global float* cv){

    float MMix, R_mix, cp_local, cv_local;
    
    MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);
    cp_local = calc_cp_mix(y_H2, y_O2, y_N2, y_H2O);
    
    R_mix = R_gas / MMix;

    //cv_local = cp_local / (cp_local-R_mix);
    cv_local = cp_local - R_mix;

    *gamma = cp_local / cv_local;
    
    *cv = cv_local;
    //printf("%f, %f, %f\n", MMix, cp_local, cv_local);

    return;
}

void calc_gamma_cp_cv(float y_H2, float y_O2, float y_N2, float y_H2O, __global float* gamma, __global float* cv, __global float* cp){

    float MMix, R_mix, cp_local, cv_local;
    
    MMix = molar_mass_mixture(y_H2, y_O2, y_N2, y_H2O);
    cp_local = calc_cp_mix(y_H2, y_O2, y_N2, y_H2O);
    
    R_mix = R_gas / MMix;

    //cv_local = cp_local / (cp_local-R_mix);
    cv_local = cp_local - R_mix;

    *gamma = cp_local / cv_local;
    
    *cv = cv_local;
    //printf("%f, %f, %f\n", MMix, cp_local, cv_local);

    *cp = cp_local;
    
    return;
}


#endif    // _SPECIES_AUXILIARY_H_INCLUDED_