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

#include <pthread.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkVertex.h>
#include <vtkCellArray.h>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMNodeIterator.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/util/XMLUni.hpp>

#include <sphPrerequisites.h>
#include <InputOutput/Particles.h>

namespace Aqua{
namespace InputOutput{

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
     * @param first First particle managed by this saver/loader.
     * @param n Number of particles managed by this saver/loader.
     * @param iset Particles set index.
     */
    VTK(ProblemSetup sim_data,
        unsigned int first,
        unsigned int n,
        unsigned int iset);

    /// Destructor
    ~VTK();

    /** @brief Save the data.
     */
    void save();

    /** @brief Load the data.
     */
    void load();

private:
    /** @brief Create a new file to write.
     * @return The file handler, NULL if errors happened.
     * @see Aqua::InputOutput::Particles::file(const char* basename,
     *                                         unsigned int start_index,
     *                                         unsigned int digits=5)
     */
    vtkXMLUnstructuredGridWriter* create();

    /** @brief Create/Update the Paraview Data File.
     *
     * Such file is used to indicates Paraview the list of files which compose
     * an animation, and the time instant of each one.
     * @return false if all gone right, true otherwise.
     */
    bool updatePVD();

    /** @brief Check if the Paraview Data File exist, creating it otherwise.
     *
     * Such file is used to indicates Paraview the list of files which compose
     * an animation, and the time instant of each one.
     * @param generate true if the file should be generated in case it does not
     * exist, false if the document should be associated to an existing file.
     * @return The document object. NULL if the file cannot be open/generated
     */
    xercesc::DOMDocument* getPVD(bool generate=true);

    /** @brief PVD file name
     * @return the PVD file name
     */
    const char* filenamePVD();

    /// PVD file name
    char* _namePVD;

    /// Launched threads ids
    std::deque<pthread_t> _tids;

    /// Next output file index
    unsigned int _next_file_index;
};  // class InputOutput

}}  // namespaces

#endif // VTK_H_INCLUDED
