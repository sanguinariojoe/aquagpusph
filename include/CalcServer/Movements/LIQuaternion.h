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

#ifndef LIQUATERNION_H_INCLUDED
#define LIQUATERNION_H_INCLUDED

// ----------------------------------------------------------------------------
// Include the Movement
// ----------------------------------------------------------------------------
#include <CalcServer/Movements/Quaternion.h>

// ----------------------------------------------------------------------------
// Include the Linear interpolator
// ----------------------------------------------------------------------------
#include <CalcServer/Movements/LinearInterpolation.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class LIQuaternion LIQuaternion.h CalcServer/Movements/LIQuaternion.h
 * @brief Linear interpolation solid quaternion based motion.
 */
class LIQuaternion : public Aqua::CalcServer::Movement::Quaternion
{
public:
    /** Constructor.
     */
    LIQuaternion();

    /** Destructor.
     */
    ~LIQuaternion();

    /** Execute the motion.
     * @return false if all gone right, true otherwise.
     */
    bool execute();

protected:
    /** Parse the input definition file.
     * @param root Input node of the parser.
     * @return false if all gone right, true otherwise.
     */
    bool _parse(xercesc::DOMElement *root);

private:
    /** Read the initial quaternion position.
     * @return false if all gone right, true otherwise.
     */
    bool setInitialPositions();

    /// Data file linear interpolator
    LinearInterpolation *_data;
};

}}} // namespace

#endif // LIQUATERNION_H_INCLUDED
