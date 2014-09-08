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
 * @brief Continuous C1 quaternion based motion.
 * (See Aqua::CalcServer::Movement::C1Quaternion for details)
 */

#ifndef C1Quaternion_H_INCLUDED
#define C1Quaternion_H_INCLUDED

#include <CalcServer/Movements/Quaternion.h>
#include <CalcServer/Movements/C1Interpolation.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class C1Quaternion C1Quaternion.h CalcServer/Movements/C1Quaternion.h
 * @brief C1 interpolated quaternion based motion.
 *
 * A prescribed motion will be imposed by a data file usage, which its data will
 * be interpolated with Aqua::CalcServer::Movement::C1Interpolation tool.
 *
 * @see Aqua::CalcServer::Movement::C1Interpolation
 * @see Aqua::CalcServer::Movement::Quaternion
 */
class C1Quaternion : public Aqua::CalcServer::Movement::Quaternion
{
public:
    /// Constructor.
    C1Quaternion();

    /// Destructor.
    ~C1Quaternion();

    /** @brief Execute the movement.
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

    /// Data file Contiguous C1 interpolation tool
    C1Interpolation *_data;
};

}}} // namespace

#endif // C1Quaternion_H_INCLUDED
