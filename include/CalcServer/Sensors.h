/*
 *  This file is part of AQUAgpusph, a free CFD program based on SPH.
 *  Copyright (C) 2012  Jose Luis Cercos Pita <jl.cercos@upm.es>
 *
 *  AQUAgpusph is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AQUAgpusph is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AQUAgpusph.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SENSORS_H_INCLUDED
#define SENSORS_H_INCLUDED

#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Sensors Sensors.h CalcServer/Sensors.h
 * @brief Sensors fields computation.
 */
class Sensors : public Aqua::CalcServer::Kernel
{
public:
    /** Constructor.
     */
    Sensors();

    /** Destructor.
     */
    ~Sensors();

    /** Sensors calculation.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

    /** Get the sensors positions.
     * @return Sensors positions.
     */
    vec* positions();

protected:
    /** Retrieve data form the device, printing it in the output file.
     * @return false if all gone right, true otherwise.
     */
    bool printOutput();

private:
    /** Setup the OpenCL stuff
     * @return false if all gone right, true otherwise.
     */
    bool setupOpenCL();

    /** Start the output file
     * @return false if all gone right, true otherwise.
     */
    bool initOutput();

    /// Number of sensors
    unsigned int _n;

    /// OpenCL script path
    char* _path;

    /// Output file
    FILE* _output;
    /// Last time when a file was printed
    float _output_time;

    /// Device stored pressure variance
    cl_mem _dev_dens_var;

    /// Positions
    vec *_pos;
    /// Pressure
    cl_float *_press;
    /// Density
    cl_float *_dens;
    /// Pressure variance
    cl_float *_dens_var;
    /// Shepard term
    cl_float *_sum_W;

    /// OpenCL program
    cl_program _program;
    /// OpenCL kernel
    cl_kernel _kernel;
    /// Global work size
    size_t _global_work_size;
    /// Local work size
    size_t _local_work_size;
};

}}  // namespace

#endif // SENSORS_H_INCLUDED
