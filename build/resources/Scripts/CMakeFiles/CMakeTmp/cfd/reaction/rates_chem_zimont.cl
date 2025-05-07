#include "resources/Scripts/types/types.h"
#include "resources/Scripts/cfd/reaction/reaction_generic.hcl"
#include "resources/Scripts/cfd/reaction/arrhenius_detonation.hcl"

__kernel void entry(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global float* rho,
                    const __global float* eint,
                    const __global float* p,
                    const __global float* T,
                    const __global float* z,
                    const __global float* y_H2,
                    const __global float* y_O2,
                    const __global float* y_N2,
                    const __global float* y_H2O,
                    const __global vec* grad_zeta,
                    //const __global float* div_u,
                    //const __global float* work_density,
                    //__global vec* dudt,
                    //__global float* drhodt,
                    __global float* deintdt,
                    __global float* dy_H2dt,
                    __global float* dy_O2dt,
                    __global float* dy_N2dt,
                    __global float* dy_H2Odt, 
                    __global float* zeta_dot,                    
                    usize N,
                    LINKLIST_LOCAL_PARAMS)
{
    const usize i = get_global_id(0);
    
    if(i >= N)
        return;
    if(imove[i] != 1){
        return;
    }
    /*drhodt[i] = -div_u[i];
    dudt[i] = -grad_p[i] + g;
    deintdt[i] = -work_density[i];*/
    
    w_rhos_zimont_def(z[i], T[i], y_H2[i], y_O2[i], y_N2[i], y_H2O[i], 
                        dy_H2dt+i, dy_O2dt+i, dy_N2dt+i, dy_H2Odt+i, 
                        deintdt+i, 
                        zeta_dot+i,
                        grad_zeta+i);
}