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
 * @brief Quaternion based motions base class.
 * (See Aqua::CalcServer::Movement::Quaternion for details)
 */

#ifndef QUATERNION_H_INCLUDED
#define QUATERNION_H_INCLUDED

#include <ProblemSetup.h>
#include <CalcServer/Movements/Movement.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class Quaternion Quaternion.h CalcServer/Movements/Quaternion.h
 * @brief Solid quaternion based motion.
 *
 * The quaternion is defined by the center COR, and the axis, referenced to an
 * earth fixed reference system.
 *
 * This class is not designed to be used directly by the AQUAgpusph binary but
 * as an utility for the overloaded quaternion classes.
 *
 * @see Quaternion.cl
 */
class Quaternion : public Aqua::CalcServer::Movement::Movement
{
public:
    /// Constructor.
    Quaternion();

    /// Destructor.
    virtual ~Quaternion();

    /** @brief Set the quaternion (by the center point and the axis vectors).
     * @param cor Center of the quaternion.
     * @param axis Axis matrix. The axis matrix must contains one axis per row.
     * @param initial true if it is the quaternion at \f$ t = 0 \f$ seconds,
     * false otherwise.
     * For the initial quaternion the relative positions of the particles to the
     * local reference system (controlled by the quaternion) will be computed.
     * @return false if all gone right, true otherwise.
     */
    bool set(vec cor, mat axis, bool initial=false);

    /** @brief Get the updated COR.
     * @return Center of the quaternion.
     */
    vec getCOR(){return _cor;}

    /** @brief Get the updated quaternion axis.
     * @return Quaternion axis matrix (distributed by rows).
     */
    mat getAxis(){return _axis;}

    /** @brief Perform the work.
     * @return false if all gone right, true otherwise.
     */
    virtual bool execute();

protected:
    /** @brief Parse the input definition file.
     * @param root Input node of the parser.
     * @return false if all gone right, true otherwise.
     */
    virtual bool _parse(xercesc::DOMElement *root);

    /** @brief Executes the motion of ghost particles walls.
     * @return false if all gone right, true otherwise.
     */
    virtual bool executeWalls();

    /** @brief Executes the motion of the computational domain.
     *
     * Domain will follow solid translations but no rotations.
     *
     * @return false if all gone right, true otherwise.
     */
    virtual bool executeDomain();

    /// Quaternion center.
    vec _cor;
    /// Quaternion axis vectors.
    mat _axis;

    /// Old quaternion center.
    vec _old_cor;
    /// Old quaternion axis vectors.
    mat _old_axis;
private:
    /** @brief Compute the particles positions respect to the Quaternion.
     * @return false if all gone right, true otherwise.
     */
    bool computePos();

    /** @brief Compute the wall vertexes (for the ghost particles) respect to
     * the Quaternion.
     * @return false if all gone right, true otherwise.
     */
    bool computeWalls();

    /** @brief Compute the domain bounds respect to the COR.
     * @return false if all gone right, true otherwise.
     */
    bool computeDomain();

    /// Relative positions to the quaternion.
    cl_mem _pos;
    /// Relative normals to the quaternion.
    cl_mem _normal;

    /// Stored walls (relative to axes)
    std::deque<InputOutput::ProblemSetup::sphGhostParticles::Wall*> _walls;

    /// Minimum domain coordinates relative to COR
    vec _domain_min;
    /// Maximum domain coordinates relative to COR
    vec _domain_max;
};

}}} // namespace

#endif // QUATERNION_H_INCLUDED
