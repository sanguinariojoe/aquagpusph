
#ifndef _SOUND_SPEED_PERFECT_GAS_H_INCLUDED_
#define _SOUND_SPEED_PERFECT_GAS_H_INCLUDED_

float SoundSpeedPerfectGas(float gamma, float p, float rho)
{
    return sqrt(gamma * p / rho);
}

#endif    // _SOUND_SPEED_PERFECT_GAS_H_INCLUDED_