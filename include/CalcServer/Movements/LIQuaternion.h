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

/** @file
 * @brief Linearly interpolated quaternion based motion.
 * (See Aqua::CalcServer::Movement::LIQuaternion for details)
 */

#ifndef LIQUATERNION_H_INCLUDED
#define LIQUATERNION_H_INCLUDED

#include <CalcServer/Movements/Quaternion.h>
#include <CalcServer/Movements/LinearInterpolation.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class LIQuaternion LIQuaternion.h CalcServer/Movements/LIQuaternion.h
 * @brief Linearly interpolated quaternion based motion.
 *
 * A prescribed motion will be imposed by a data file usage, which its data will
 * be interpolated with Aqua::CalcServer::Movement::LinearInterpolation tool.
 *
 * @warning It is strongly recommended to use
 * Aqua::CalcServer::Movement::C1Quaternion instead of this motion model.
 * @see Aqua::CalcServer::Movement::LinearInterpolation
 * @see Aqua::CalcServer::Movement::Quaternion
 * @see Aqua::CalcServer::Movement::C1Quaternion
 */
class LIQuaternion : public Aqua::CalcServer::Movement::Quaternion
{
public:
    /// Constructor.
    LIQuaternion();

    /// Destructor.
    ~LIQuaternion();

    /** @brief Execute the motion.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

protected:
    /** @brief Get the prescribed motion data file from the XML definition file.
     * @param root Input node of the parser.
     * @return false if all gone right, true otherwise.
     */
    bool _parse(xercesc::DOMElement *root);

private:
    /** @brief Read and set the initial quaternion position.
     * @return false if all gone right, true otherwise.
     */
    bool setInitialPositions();

    /// Data file linear interpolation tool
    LinearInterpolation *_data;
};

}}} // namespace

#endif // LIQUATERNION_H_INCLUDED
