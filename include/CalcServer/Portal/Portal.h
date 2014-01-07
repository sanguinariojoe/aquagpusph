/*
 * This source file is part of AQUA-gpusph.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor Boston, MA 02110-1301,  USA
 */
/*
	Authors:
	- Cercos Pita, Jose Luis
	- Miguel Gonzalez, Leo
	- Rey Villaverde, Anton
	- Saelices, Jaime
	- Souto Iglesias, Antonio
*/

#ifndef PORTAL_H_INCLUDED
#define PORTAL_H_INCLUDED

// ----------------------------------------------------------------------------
// Include the Problem setup manager header
// ----------------------------------------------------------------------------
#include <ProblemSetup.h>

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{ namespace Portal{

/** @class Portal Portal.h CalcServer/Portal/Portal.h
 * @brief Base portal class, that only teleport particles that
 * pass throught outlet portal to the inlet portal.
 */
class Portal : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Portal(InputOutput::ProblemSetup::sphPortal *portal);

	/** Destructor.
	 */
	~Portal();

	/** Teleport particles from outlet to inlet.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
	/** Setup OpenCL kernel
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupOpenCL();

	/// OpenCL script path
	char* _path;
	/// Main portal
	InputOutput::ProblemSetup::sphPortal *mPortal;

	/// OpenCL program
	cl_program program;
	/// OpenCL vertices set kernel
	cl_kernel kernel;
	/// Global work size (calculated with local_work_size).
	size_t _global_work_size;
	/// Local work size (default value = 256)
	size_t _local_work_size;
};

}}}  // namespace

#endif // PORTAL_H_INCLUDED
