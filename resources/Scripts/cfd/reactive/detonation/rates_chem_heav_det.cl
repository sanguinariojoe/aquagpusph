#include "resources/Scripts/types/types.h"
#include "resources/Scripts/cfd/reactive/reaction_generic.hcl"
#include "resources/Scripts/cfd/reactive/detonation/heav_detonation.hcl"

__kernel void entry(const __global unsigned int* iset,
                    const __global int* imove,
                    const __global float* rho,
                    const __global float* eint,
                    const __global float* p,
                    const __global float* T,
                    const __global float* z,
                    const __global species_t* ys,
                    const __global int* trigger,
                    __global float* deintdt,
                    __global species_t* dysdt,
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

    w_rhos_heav_det(z[i], T[i], ys[i], trigger[i], dysdt+i, deintdt+i, zeta_dot+i);
}