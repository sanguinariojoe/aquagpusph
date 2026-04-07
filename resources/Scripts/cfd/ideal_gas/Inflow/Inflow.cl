/** @addtogroup cfd
 * @{
 */

/** @file
 * @brief Vanish the velocity and desnity rates of variation of the velocity
 * and density for the dummy particles of the inflow.
 */

#include EOS_MODEL
#include "resources/Scripts/types/types.h"

/** @brief Particles generation at the inflow boundary condition.
 *
 * Particles are generated just when the inflow is starving, i.e. the
 * previously generated layer of particles have moved more than dr. To do that
 * /outlet is extracting the particles from the "buffer", which are the last
 * particles in the sorted list.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param iset Set of particles index.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param rho Density \f$ \rho \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param m Mass \f$ m \f$.
 * @param p Pressure \f$ p \f$.
 * @param refd Density of reference of the fluid \f$ \rho_0 \f$.
 * @param N Number of particles.
 * @param nbuffer Number of buffer particles.
 * @param dt Time step \f$ \Delta t \f$.
 * @param cs Speed of sound \f$ c_s \f$.
 * @param p0 Background pressure \f$ p_0 \f$.
 * @param g Gravity acceleration \f$ \mathbf{g} \f$.
 * @param dr Distance between particles \f$ \Delta r \f$.
 * @param inflow_r Lower corner of the inflow square.
 * @param inflow_ru Square U vector.
 * @param inflow_rv Square V vector.
 * @param inflow_N Number of particles to be generated in each direction.
 * @param inflow_n = Velocity direction of the generated particles.
 * @param inflow_U = Constant inflow velocity magnitude
 * @param inflow_rFS The point where the pressure is the reference one (0 Pa).
 * @param inflow_R Accumulated displacement (to be added to the generation point)
 * @param inflow_starving Is the inflow starving, so we need to feed it?
 */
__kernel void feed(svec2 inflow_N,
                   int inflow_starving,
                   float inflow_rho,
                   float inflow_e,
                   float inflow_gamma,   
                   usize N,
                   usize nbuffer,
                   __global float* restrict p,
                   __global float* restrict rho,
                   __global float* restrict drhodt,
                   __global float* restrict eint,
                   __global float* restrict gamma,
                   __global float* restrict deintdt
                   )
{
    // find position in global arrays
    const usize i = get_global_id(0);
    if(inflow_starving == 0)
        return;
    if((i >= nbuffer) || (i >= (inflow_N.x * inflow_N.y))){
        // Either the thread has not a buffer particle to consume or such buffer
        // particle is not required
        return;
    }
    const usize ii = i + N - nbuffer;
 
    rho[ii] = inflow_rho;
    drhodt[ii] = 0.f;    
    eint[ii] = inflow_e;
    deintdt[ii] = 0.f;
    gamma[ii] = inflow_gamma;
    p[ii] = p_from_rho_eint(inflow_gamma, inflow_rho, inflow_e);
    
    //p[ii] = (inflow_gamma - 1.0f) * inflow_rho * inflow_e;
}

/** @brief Vanish the velocity and desnity rates of variation of the velocity
 * and density for the dummy particles of the inflow.
 *
 * @param imove Moving flags.
 *   - imove > 0 for regular fluid particles.
 *   - imove = 0 for sensors.
 *   - imove < 0 for boundary elements/particles.
 * @param r Position \f$ \mathbf{r} \f$.
 * @param u Velocity \f$ \mathbf{u} \f$.
 * @param dudt Velocity rate of change \f$ \frac{d \mathbf{u}}{d t} \f$.
 * @param drhodt Density rate of change \f$ \frac{d \rho}{d t} \f$.
 * @param N Number of particles.
 * @param inflow_r Lower corner of the inflow square.
 * @param inflow_U Velocity magnitude of the generated particles.
 * @param inflow_n Velocity direction of the generated particles.
 */
__kernel void rates(__global int* restrict imove,
                    __global vec* restrict r,
                    __global float* restrict deintdt,
                    usize N,
                    vec inflow_r,
                    vec inflow_n)
{
    // find position in global arrays
    const usize i = get_global_id(0);
    if(i >= N)
        return;
    if(imove[i] != 1)
        return;

    // Discard the particles already passed through the inflow
    if(dot(r[i] - inflow_r, inflow_n) > 0.f)
        return;

    deintdt[i] = 0.f;
}