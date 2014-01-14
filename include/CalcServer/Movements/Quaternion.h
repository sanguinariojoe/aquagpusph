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

#ifndef QUATERNION_H_INCLUDED
#define QUATERNION_H_INCLUDED

#include <ProblemSetup.h>
#include <CalcServer/Movements/Movement.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class Quaternion Quaternion.h CalcServer/Movements/Quaternion.h
 * @brief Solid quaternion based motion. Quaternion, that must be manually
 * provided, defines the solid position for each time instant.
 * This class is not designed to be used by the AQUAgpusph binary but as an
 * utility for the overloaded quaternion classes.
 */
class Quaternion : public Aqua::CalcServer::Movement::Movement
{
public:
	/** Constructor.
	 * @param data_file Data file path.
	 * @note Data file can be omissed at construction, but be ensure
	 * to provide it later.
	 */
	Quaternion();

	/** Destructor.
	 */
	~Quaternion();

	/** Set the quaternion (by the center point and the axis vectors).
	 * @param cor Center of the quaternion.
	 * @param axis Axis matrix. The axis matrix must
	 * contains one axis per row.
	 * @param initial true if it is the initial set of the quaternion.
	 * Relative positions will be recomputed, and backup values asssigned.
	 * @return false if all gone right, true otherwise.
	 */
	bool set(vec cor, mat axis, bool initial=false);

	/** Get the updated COR.
	 * @return Center of the quaternion.
	 */
	vec getCOR(){return _cor;}

	/** Get the updated quaternion axis.
	 * @return Quaternion axis matrix (distributed by rows).
	 */
	mat getAxis(){return _axis;}

	/** Execute the motion.
	 * @return false if all gone right, true otherwise.
	 */
	virtual bool execute();

protected:
	/** Parse the input definition file.
	 * @param root Input node of the parser.
	 * @return false if all gone right, true otherwise.
	 */
	virtual bool _parse(xercesc::DOMElement *root);

	/** Executes the movement over walls (ghost particles).
	 * @return false if all gone right, true otherwise.
	 */
	virtual bool executeWalls();

	/** Executes the domain motion if requested. Domain
	 * will follow solid translations but not rotations.
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
	/** Compute the particles positions respect to the Quaternion.
	 * @return false if all gone right, true otherwise.
	 */
	bool computePos();

	/** Compute the wall vertexes (for the ghost particles) respect to the
	 * Quaternion.
	 * @return false if all gone right, true otherwise.
	 */
	bool computeWalls();

	/** Compute the domain bounds respect to the COR. (the domain bounds are
     * not affected by the rotations)
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
