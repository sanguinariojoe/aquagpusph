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

// ----------------------------------------------------------------------------
// Include Problem setup
// ----------------------------------------------------------------------------
#include <ProblemSetup.h>

// ----------------------------------------------------------------------------
// Include the Movement
// ----------------------------------------------------------------------------
#include <CalcServer/Movements/Movement.h>

namespace Aqua{ namespace CalcServer{ namespace Movement{

/** @class Quaternion Quaternion.h CalcServer/Movements/Quaternion.h
 * @brief Solid quaternion movement. Quaternion, that must be procided
 * manually defines the solid position at any time. This class can be
 * used if AQUAgpusph is used as library, or as base for other more
 * complex Quaternion based movements.
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

	/** Set quaternion (as point and axis vectors).
	 * @param cor Center of rotation.
	 * @param axis Axis matrix. The axis matrix must
	 * contains each axis at each row.
	 * @param initial true if initial set of quaternion
	 * being performed. Relative positions will be
	 * recomputed, and backup values asssigned.
	 * @return false if all gone right, true otherwise.
	 */
	bool set(vec cor, mat axis, bool initial=false);

	/** Get updated COR.
	 * @return Center of rotation.
	 */
	vec getCOR(){return _cor;}

	/** Get updated quaternion axis.
	 * @return Quaternion axis matrix.
	 */
	mat getAxis(){return mAxis;}

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
	mat mAxis;

	/// Old quaternion center.
	vec mOldCOR;
	/// Old quaternion axis vectors.
	mat mOldAxis;
private:
	/** Computes particles positions relative to Quaternion.
	 * @return false if all gone right, true otherwise.
	 */
	bool computePos();

	/** Computes wall vertexes (ghost particles) relative
	 * to Quaternion.
	 * @return false if all gone right, true otherwise.
	 */
	bool computeWalls();

	/** Computes domain boundsrelative to COR.
	 * @return false if all gone right, true otherwise.
	 */
	bool computeDomain();

	/// Relative positions to the quaternion.
	cl_mem mRelPos;
	/// Relative normals to the quaternion.
	cl_mem mRelNormal;

	/// Stored walls (relative to axes)
	std::deque<InputOutput::ProblemSetup::sphGhostParticles::Wall*> walls;

	/// Minimum domain coordinates relative to COR
	vec domain_min;
	/// Maximum domain coordinates relative to COR
	vec domain_max;
};

}}} // namespace

#endif // QUATERNION_H_INCLUDED
