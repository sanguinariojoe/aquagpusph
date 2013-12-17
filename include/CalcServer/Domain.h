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

#ifndef DOMAIN_H_INCLUDED
#define DOMAIN_H_INCLUDED

// ----------------------------------------------------------------------------
// Include Generic kernel
// ----------------------------------------------------------------------------
#include <CalcServer/Kernel.h>

namespace Aqua{ namespace CalcServer{

/** @class Domain Domain.h CalcServer/Domain.h
 * @brief Kernel designed to test if any particle
 * is out of domain bounds. If any particle is
 * detected out of domain bounds must be converted
 * into a fixed zero mass particle.
 */
class Domain : public Aqua::CalcServer::Kernel
{
public:
	/** Constructor.
	 */
	Domain();

	/** Destructor.
	 */
	~Domain();

	/** Executes time integration Domain stage.
	 * @return false if all gone right. \n true otherwise.
	 */
	bool execute();

private:
	/** Setup OpenCL kernel
	 * @return false if all gone right. \n true otherwise.
	 */
	bool setupOpenCL();

	/// OpenCL script path
	char* mPath;

	/// OpenCL program
	cl_program clProgram;
	/// OpenCL kernel
	cl_kernel clKernel;
	/// Global work size
	size_t clGlobalWorkSize;
	/// Local work size
	size_t clLocalWorkSize;
};

}}  // namespace

#endif // DOMAIN_H_INCLUDED
