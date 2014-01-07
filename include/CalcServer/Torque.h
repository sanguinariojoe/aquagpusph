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

#ifndef TORQUE_H_INCLUDED
#define TORQUE_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

// ----------------------------------------------------------------------------
// Include Reduction tool
// ----------------------------------------------------------------------------
#include <CalcServer/Reduction.h>

namespace Aqua{ namespace CalcServer{

/** @class Torque Torque.h CalcServer/Torque.h
 * @brief Computes fluid torque relative to a
 * provided point called COR (Center of rotation).
 */
class Torque : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Torque();

	/** Destructor.
	 */
	~Torque();

	/** Set Center of rotation.
	 * @param COR Center of rotation.
	 */
	void cor(vec COR){mCOR = COR;}

	/** Get Center of rotation.
	 * @return Center of rotation.
	 */
	vec cor(){return mCOR;}

	/** Get resultant torque.
	 * @return Fluid torque.
	 * @warning Remember call execute before use this method.
	 */
	vec torque(){return mTorque;}

	/** Get resultant force.
	 * @return Fluid force.
	 * @warning Remember call execute before use this method.
	 */
	vec force(){return mForce;}

	/** Performs space discretization.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
	/** Setup Torque OpenCL stuff.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupTorque();

	/** Setup Reduction
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupReduction();

	/// Center of rotation
	vec mCOR;

	/// Server allocated torque.
	cl_mem mDevTorque;
	/// Host allocated torque
	vec mTorque;
	/// Server allocated torque.
	cl_mem mDevForce;
	/// Host allocated torque
	vec mForce;
	/// Kernel path
	char *_path;
	/// OpenCL program
	cl_program _program;
	/// OpenCL kernel
	cl_kernel _kernel;
	/// Global work size
	size_t _global_work_size;
	/// Local work size
	size_t _local_work_size;
    /// Torque value reduction tool
    Reduction *torqueReduction;
    /// Force value reduction tool
    Reduction *forceReduction;
};

}}  // namespace

#endif // TORQUE_H_INCLUDED
