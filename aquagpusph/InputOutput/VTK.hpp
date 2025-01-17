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
 * @brief Particles VTK data files loader/saver.
 * (See Aqua::InputOutput::VTK for details)
 */

#ifndef VTK_H_INCLUDED
#define VTK_H_INCLUDED

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkTypeInt32Array.h>
#include <vtkTypeInt64Array.h>
#include <vtkTypeUInt32Array.h>
#include <vtkTypeUInt64Array.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkVertex.h>
#include <vtkCellArray.h>

#ifdef VTK_NLOHMANN_JSON
#include <vtk_nlohmannjson.h>
#include VTK_NLOHMANN_JSON(json.hpp)
#else
#include <nlohmann/json.hpp>
#endif // VTK_NLOHMANN JSON

#include "aquagpusph/sphPrerequisites.hpp"
#include "Particles.hpp"

namespace Aqua {
namespace InputOutput {

/** \class VTK VTK.h InputOutput/VTK.h
 * @brief VTK particles data files loader/saver.
 *
 * VTK is a visualization format, to learn more about it please visit the
 * following web page:
 *
 * http://www.vtk.org
 *
 * This type of files can be easily post-processed with Paraview:
 *
 * http://www.paraview.org
 *
 * The fields to be saved/loaded are:
 *   -# \f$ \mathbf{r} \f$ . \f$ x \f$
 *   -# \f$ \mathbf{r} \f$ . \f$ y \f$
 *   -# \f$ \mathbf{r} \f$ . \f$ z \f$ (For 3D cases only)
 *   -# \f$ \mathbf{n} \f$ . \f$ x \f$
 *   -# \f$ \mathbf{n} \f$ . \f$ y \f$
 *   -# \f$ \mathbf{n} \f$ . \f$ z \f$ (For 3D cases only)
 *   -# \f$ \mathbf{u} \f$ . \f$ x \f$
 *   -# \f$ \mathbf{u} \f$ . \f$ y \f$
 *   -# \f$ \mathbf{u} \f$ . \f$ z \f$ (For 3D cases only)
 *   -# \f$ \frac{d \mathbf{u}}{dt} \f$ . \f$ x \f$
 *   -# \f$ \frac{d \mathbf{u}}{dt} \f$ . \f$ y \f$
 *   -# \f$ \frac{d \mathbf{u}}{dt} \f$ . \f$ z \f$ (For 3D cases only)
 *   -# \f$ \rho \f$
 *   -# \f$ \frac{d \rho}{dt} \f$
 *   -# \f$ m \f$
 *   -# moving flag (see Aqua::InputOutput::Fluid::imove)
 */
class VTK : public Particles
{
  public:
	/** @brief Constructor
	 * @param sim_data Simulation data
	 * @param iset Particles set index.
	 * @param offset First particle managed by this saver/loader.
	 * @param n Number of particles managed by this saver/loader. If 0,
	 * the number of particles will be obtained from the input file (thus only
	 * valid for loaders)
	 */
	VTK(ProblemSetup& sim_data,
	    unsigned int iset,
	    size_t offset,
	    size_t n = 0);

	/// Destructor
	~VTK();

	/** @brief Save the data.
	 *
	 * @param t Simulation time
	 */
	void save(float t);

	/** @brief Load the data.
	 */
	void load();

	/** @brief Print the data to a file
	 * @note This method is public to work with the OpenCL callbacks, but it is
	 * not meant to be called by the users
	 */
	void print_file() final;

  private:
	/** @brief Build the vtkUnstructuredGrid object
	 * @note This method is quite slow, so it should be called just once.
	 * Afterwards the returned object can be edited and saved for optimum
	 * performance
	 * @see ::editVTK()
	 */
	vtkUnstructuredGrid* makeVTK();

	/** @brief Edit the vtkUnstructuredGrid object
	 * @see ::makeVTK()
	 */
	void editVTK();

	/** @brief Compute the number of particles handled by this instance
	 * @return Number of particles
	 */
	size_t compute_n();

	/** @brief Create a new file to write.
	 * @return The file handler, NULL if errors happened.
	 * @see Aqua::InputOutput::Particles::file(const char* basename,
	 *                                         unsigned int start_index,
	 *                                         unsigned int digits=5)
	 */
	vtkXMLUnstructuredGridWriter* create();

	/** @brief Create/Update the VTK series json file.
	 *
	 * Such file provides info about the simulation time associated to each
	 * snapshot.
	 * @param t Simulation time
	 */
	void updateSeries(float t);

	/** @brief VTK series file name
	 * @return the VTK series file name
	 */
	const std::string filenameSeries();

	/// Next output file index
	unsigned int _next_file_index;

	/// VTK object to be edited and saved
	vtkUnstructuredGrid* _vtk;

	/// VTK series file name
	std::string _name_series;

	/// VTK series data
	nlohmann::json _data_series;

}; // class InputOutput

}
} // namespaces

#endif // VTK_H_INCLUDED
