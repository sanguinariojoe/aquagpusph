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

#ifndef C1Quaternion_H_INCLUDED
#define C1Quaternion_H_INCLUDED

// ----------------------------------------------------------------------------
// Include the Movement
// ----------------------------------------------------------------------------
#include <CalcServer/Movements/Quaternion.h>

// ----------------------------------------------------------------------------
// Include the C1 interpolator
// ----------------------------------------------------------------------------
#include <CalcServer/Movements/C1Interpolation.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class C1Quaternion C1Quaternion.h CalcServer/Movements/C1Quaternion.h
 * @brief C1 interpolated solid quaternion based motion.
 */
class C1Quaternion : public Aqua::CalcServer::Movement::Quaternion
{
public:
    /** Constructor.
     */
    C1Quaternion();

    /** Destructor.
     */
    ~C1Quaternion();

    /** Execute the movement.
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

    /// Data file Contiguous C1 interpolator
    C1Interpolation *_data;
};

}}} // namespace

#endif // C1Quaternion_H_INCLUDED
