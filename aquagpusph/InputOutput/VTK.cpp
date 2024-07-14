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

#include "aquagpusph/sphPrerequisites.hpp"

#ifdef HAVE_VTK

#include <unistd.h>
#include <signal.h>

#include "VTK.hpp"
#include "Logger.hpp"
#include "aquagpusph/ProblemSetup.hpp"
#include "aquagpusph/CalcServer/CalcServer.hpp"
#include "aquagpusph/AuxiliarMethods.hpp"

#include <vector>

using json = nlohmann::json;

namespace Aqua {
namespace InputOutput {

VTK::VTK(ProblemSetup& sim_data,
         unsigned int iset,
         size_t first,
         size_t n_in)
  : Particles(sim_data, iset, first, n_in)
  , _next_file_index(0)
  , _vtk(NULL)
  , _name_series("")
  , _data_series({ { "file-series-version", "1.0" },
                   { "files", json::array() } })
{
	if (n() == 0) {
		n(compute_n());
	}
}

VTK::~VTK()
{
	waitForSavers();
	if (_vtk)
		_vtk->Delete();
}

void
VTK::load()
{
	unsigned int n, N, progress;
	int aux;
	cl_int err_code;
	CalcServer::CalcServer* C = CalcServer::CalcServer::singleton();

	loadDefault();

	std::ostringstream msg;
	msg << "Loading particles from VTK file \""
	    << simData().sets.at(setId())->inputPath() << "\"..." << std::endl;
	LOG(L_INFO, msg.str());

	vtkSmartPointer<vtkXMLUnstructuredGridReader> f =
	    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

	if (!f->CanReadFile(simData().sets.at(setId())->inputPath().c_str())) {
		LOG(L_ERROR, "The file cannot be read.\n");
		throw std::runtime_error("Failure reading file");
	}

	f->SetFileName(simData().sets.at(setId())->inputPath().c_str());
	f->Update();

	vtkSmartPointer<vtkUnstructuredGrid> grid = f->GetOutput();

	// Assert that the number of particles is right
	n = bounds().y - bounds().x;
	N = (size_t)grid->GetNumberOfPoints();
	if (n != N) {
		std::ostringstream msg;
		msg << "Expected " << n << " particles, but the file contains just "
		    << N << " ones." << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Invalid number of particles in file");
	}

	// Check the fields to read
	std::vector<std::string> fields = simData().sets.at(setId())->inputFields();
	if (!fields.size()) {
		LOG0(L_ERROR, "0 fields were set to be read from the file.\n");
		throw std::runtime_error("No fields have been marked to read");
	}
	bool have_r = false;
	for (auto field : fields) {
		if (!field.compare("r")) {
			have_r = true;
			break;
		}
	}
	if (!have_r) {
		LOG0(L_ERROR, "\"r\" field was not set to be read from the file.\n");
		throw std::runtime_error("\"r\" field is mandatory");
	}

	// Setup an storage
	std::vector<void*> data;
	Variables* vars = C->variables();
	for (auto field : fields) {
		if (!vars->get(field)) {
			std::ostringstream msg;
			msg << "Undeclared variable \"" << field << "\" set to be read."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable");
		}
		if (vars->get(field)->type().find('*') == std::string::npos) {
			std::ostringstream msg;
			msg << "Can't read scalar variable \"" << field << "\"."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable type");
		}
		ArrayVariable* var = (ArrayVariable*)vars->get(field);
		size_t typesize = vars->typeToBytes(var->type());
		size_t len = var->size() / typesize;
		if (len < bounds().y) {
			std::ostringstream msg;
			msg << "Array variable \"" << field << "\" is not long enough."
			    << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::runtime_error("Invalid variable length");
		}
		void* store = malloc(typesize * n);
		if (!store) {
			std::ostringstream msg;
			msg << "Failure allocating " << typesize * n
			    << "bytes for variable \"" << field << "\"." << std::endl;
			LOG(L_ERROR, msg.str());
			throw std::bad_alloc();
		}
		data.push_back(store);
	}

	progress = -1;
	vtkSmartPointer<vtkPoints> vtk_points = grid->GetPoints();
	vtkSmartPointer<vtkPointData> vtk_data = grid->GetPointData();
	for (size_t i = 0; i < n; i++) {
		for (unsigned int j = 0; j < fields.size(); j++) {
			if (!fields.at(j).compare("r")) {
				double* vect = vtk_points->GetPoint(i);
				if (C->have_3d()) {
					vec4* ptr = (vec4*)data.at(j);
					ptr[i].x = vect[0];
					ptr[i].y = vect[1];
					ptr[i].z = vect[2];
					ptr[i].w = 0.f;
				} else {
					vec2* ptr = (vec2*)data.at(j);
					ptr[i].x = vect[0];
					ptr[i].y = vect[1];
				}
				continue;
			}
			ArrayVariable* var = (ArrayVariable*)vars->get(fields.at(j));
			size_t type_size = vars->typeToBytes(var->type());
			unsigned int n_components = vars->typeToN(var->type());
			if (startswith(var->type(), "unsigned int") ||
			    startswith(var->type(), "uivec")) {
				vtkSmartPointer<vtkTypeUInt32Array> vtk_array =
				    (vtkTypeUInt32Array*)(vtk_data->GetArray(
				        fields.at(j).c_str(), aux));
				for (unsigned int k = 0; k < n_components; k++) {
					uicl component = vtk_array->GetComponent(i, k);
					size_t offset = type_size * i + sizeof(uicl) * k;
					memcpy((char*)data.at(j) + offset,
					       &component,
					       sizeof(uicl));
				}
			} else if (startswith(var->type(), "unsigned long") ||
			           startswith(var->type(), "ulvec")) {
				vtkSmartPointer<vtkTypeUInt64Array> vtk_array =
				    (vtkTypeUInt64Array*)(vtk_data->GetArray(
				        fields.at(j).c_str(), aux));
				for (unsigned int k = 0; k < n_components; k++) {
					ulcl component = vtk_array->GetComponent(i, k);
					size_t offset = type_size * i + sizeof(ulcl) * k;
					memcpy((char*)data.at(j) + offset,
					       &component,
					       sizeof(ulcl));
				}
			} else if (startswith(var->type(), "int") ||
			    startswith(var->type(), "ivec")) {
				vtkSmartPointer<vtkTypeInt32Array> vtk_array =
				    (vtkTypeInt32Array*)(vtk_data->GetArray(
				        fields.at(j).c_str(), aux));
				for (unsigned int k = 0; k < n_components; k++) {
					icl component = vtk_array->GetComponent(i, k);
					size_t offset = type_size * i + sizeof(icl) * k;
					memcpy((char*)data.at(j) + offset,
					       &component,
					       sizeof(icl));
				}
			} else if (startswith(var->type(), "long") ||
			           startswith(var->type(), "lvec")) {
				vtkSmartPointer<vtkTypeInt64Array> vtk_array =
				    (vtkTypeInt64Array*)(vtk_data->GetArray(
				        fields.at(j).c_str(), aux));
				for (unsigned int k = 0; k < n_components; k++) {
					lcl component = vtk_array->GetComponent(i, k);
					size_t offset = type_size * i + sizeof(lcl) * k;
					memcpy((char*)data.at(j) + offset,
					       &component,
					       sizeof(lcl));
				}
			} else if (var->type().find("float") != std::string::npos ||
			           var->type().find("vec") != std::string::npos ||
			           var->type().find("matrix") != std::string::npos) {
				vtkSmartPointer<vtkFloatArray> vtk_array =
				    (vtkFloatArray*)(vtk_data->GetArray(fields.at(j).c_str(),
				                                        aux));
				for (unsigned int k = 0; k < n_components; k++) {
					fcl component = vtk_array->GetComponent(i, k);
					size_t offset = type_size * i + sizeof(fcl) * k;
					memcpy((char*)data.at(j) + offset,
					       &component,
					       sizeof(fcl));
				}
			} else if (var->type().find("double") != std::string::npos ||
			           var->type().find("dvec") != std::string::npos) {
				vtkSmartPointer<vtkFloatArray> vtk_array =
				    (vtkFloatArray*)(vtk_data->GetArray(fields.at(j).c_str(),
				                                        aux));
				for (unsigned int k = 0; k < n_components; k++) {
					dcl component = vtk_array->GetComponent(i, k);
					size_t offset = type_size * i + sizeof(dcl) * k;
					memcpy((char*)data.at(j) + offset,
					       &component,
					       sizeof(dcl));
				}
			}
		}
		if (progress != i * 100 / n) {
			progress = i * 100 / n;
			if (!(progress % 10)) {
				std::ostringstream msg;
				msg << "\t\t" << progress << "%" << std::endl;
				LOG(L_DEBUG, msg.str());
			}
		}
	}

	// Send the data to the server and release it
	for (unsigned int i = 0; i < fields.size(); i++) {
		ArrayVariable* var = (ArrayVariable*)vars->get(fields.at(i));
		size_t typesize = vars->typeToBytes(var->type());
		cl_mem mem = *(cl_mem*)var->get();
		err_code = clEnqueueWriteBuffer(C->command_queue(),
		                                mem,
		                                CL_TRUE,
		                                typesize * bounds().x,
		                                typesize * n,
		                                data.at(i),
		                                0,
		                                NULL,
		                                NULL);
		free(data.at(i));
		data.at(i) = NULL;
		if (err_code != CL_SUCCESS) {
			std::ostringstream msg;
			msg << "Failure sending variable \"" << fields.at(i)
			    << "\" to the computational device." << std::endl;
			LOG(L_ERROR, msg.str());
			Logger::singleton()->printOpenCLError(err_code);
			throw std::runtime_error("OpenCL error");
		}
	}
	data.clear();
}

/** @brief Data structure to send the data to a parallel writer thread
 */
typedef struct
{
	/// The field names
	std::vector<std::string> fields;
	/// Bounds of the particles index managed by this writer
	uivec2 bounds;
	/// Screen manager
	Logger* S;
	/// VTK arrays
	CalcServer::CalcServer* C;
	/// The data associated to each field
	std::vector<void*> data;
	/// The VTK file decriptor
	vtkXMLUnstructuredGridWriter* f;
} data_pthread;

void
VTK::print_file()
{
	// Do a fast check about VTK and 64 bits
	const size_t n = bounds().y - bounds().x;
	if (n > VTK_ID_MAX) {
		LOG(L_ERROR,
		    std::string("Failure writing the VTK file for set ") +
		    std::to_string(setId()) + "\n");
		LOG0(L_DEBUG,
		     std::string("VTK can handle ") + std::to_string(VTK_ID_MAX) +
		     " particles, but " + std::to_string(n) + " ones were found\n");
		LOG0(L_DEBUG, "Some possible solutions are:\n");
		LOG0(L_DEBUG, " - Reducing the number of particles\n");
		LOG0(L_DEBUG, " - Splitting the particles set on smaller ones\n");
		LOG0(L_DEBUG, " - Recompiling VTK with -DVTK_USE_64BIT_IDS=ON\n");
	}

	// Call create ASAP so we can log that we are working on printing the file
	auto f = create();

	if (!_vtk)
		_vtk = makeVTK();
	else {
		editVTK();
	}

	// Write file
#if VTK_MAJOR_VERSION <= 5
	f->SetInput(_vtk);
#else  // VTK_MAJOR_VERSION
	f->SetInputData(_vtk);
#endif // VTK_MAJOR_VERSION

	if (!f->Write()) {
		LOG(L_ERROR,
		    std::string("Failure writing \"") + f->GetFileName() +
		        "\" VTK file.\n");
	}

	LOG(L_INFO, std::string("Wrote \"") + f->GetFileName() + "\" VTK file.\n");
	f->Delete();

	updateSeries(time());

	Particles::print_file();
}

vtkUnstructuredGrid*
VTK::makeVTK()
{
	auto C = CalcServer::CalcServer::singleton();
	auto fields = simData().sets.at(setId())->outputFields();

	// Create storage arrays
	const size_t n = bounds().y - bounds().x;
	std::vector<vtkSmartPointer<vtkDataArray>> vtk_arrays;
	Variables* vars = C->variables();
	for (auto field : fields) {
		// If we are here is because the field exists and it has been correctly
		// downloaded
		ArrayVariable* var = (ArrayVariable*)vars->get(field);

		unsigned int n_components = vars->typeToN(var->type());
		if (startswith(var->type(), "unsigned int") ||
		    startswith(var->type(), "uivec")) {
			vtkSmartPointer<vtkTypeUInt32Array> vtk_array =
			    vtkSmartPointer<vtkTypeUInt32Array>::New();
			vtk_array->SetNumberOfComponents(n_components);
			vtk_array->Allocate(n);
			vtk_array->SetNumberOfTuples(n);
			vtk_array->SetName(field.c_str());
			vtk_arrays.push_back(vtk_array);
		} else if (startswith(var->type(), "unsigned long") ||
		           startswith(var->type(), "ulvec")) {
			vtkSmartPointer<vtkTypeUInt64Array> vtk_array =
			    vtkSmartPointer<vtkTypeUInt64Array>::New();
			vtk_array->SetNumberOfComponents(n_components);
			vtk_array->Allocate(n);
			vtk_array->SetNumberOfTuples(n);
			vtk_array->SetName(field.c_str());
			vtk_arrays.push_back(vtk_array);
		} else if (startswith(var->type(), "int") ||
		           startswith(var->type(), "ivec")) {
			vtkSmartPointer<vtkTypeInt32Array> vtk_array =
			    vtkSmartPointer<vtkTypeInt32Array>::New();
			vtk_array->SetNumberOfComponents(n_components);
			vtk_array->Allocate(n);
			vtk_array->SetNumberOfTuples(n);
			vtk_array->SetName(field.c_str());
			vtk_arrays.push_back(vtk_array);
		} else if (startswith(var->type(), "long") ||
		           startswith(var->type(), "lvec")) {
			vtkSmartPointer<vtkTypeInt64Array> vtk_array =
			    vtkSmartPointer<vtkTypeInt64Array>::New();
			vtk_array->SetNumberOfComponents(n_components);
			vtk_array->Allocate(n);
			vtk_array->SetNumberOfTuples(n);
			vtk_array->SetName(field.c_str());
			vtk_arrays.push_back(vtk_array);
		} else if (startswith(var->type(), "float") ||
		           startswith(var->type(), "vec") ||
		           startswith(var->type(), "matrix")) {
			vtkSmartPointer<vtkFloatArray> vtk_array =
			    vtkSmartPointer<vtkFloatArray>::New();
			vtk_array->SetNumberOfComponents(n_components);
			vtk_array->Allocate(n);
			vtk_array->SetNumberOfTuples(n);
			vtk_array->SetName(field.c_str());
			vtk_arrays.push_back(vtk_array);
		} else if (startswith(var->type(), "double") ||
		           startswith(var->type(), "dvec")) {
			vtkSmartPointer<vtkDoubleArray> vtk_array =
			    vtkSmartPointer<vtkDoubleArray>::New();
			vtk_array->SetNumberOfComponents(n_components);
			vtk_array->Allocate(n);
			vtk_array->SetNumberOfTuples(n);
			vtk_array->SetName(field.c_str());
			vtk_arrays.push_back(vtk_array);
		}
	}

	vtkSmartPointer<vtkVertex> vtk_vertex;
	vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
	vtk_points->Allocate(n);
	vtk_points->SetNumberOfPoints(n);
	vtkSmartPointer<vtkCellArray> vtk_cells =
	    vtkSmartPointer<vtkCellArray>::New();
	vtk_cells->Allocate(n);
	vtk_cells->SetNumberOfCells(n);

	for (unsigned int i = 0; i < fields.size(); i++) {
		void* ptr = data()[fields[i]];

		if (!fields[i].compare("r")) {
			// The mesh points are a bit of a special case
			for (size_t j = 0; j < n; j++) {
				if (C->have_3d()) {
					vec4* coords = (vec4*)ptr;
					vtk_points->SetPoint(j, coords[j].x,
					                        coords[j].y,
					                        coords[j].z);
				} else {
					vec2* coords = (vec2*)ptr;
					vtk_points->SetPoint(j, coords[j].x,
					                        coords[j].y,
					                        0.0);
				}
				vtk_vertex = vtkSmartPointer<vtkVertex>::New();
				vtk_vertex->GetPointIds()->SetId(0, j);
				vtk_cells->InsertNextCell(vtk_vertex);
			}
			continue;
		}

		ArrayVariable* var = (ArrayVariable*)(vars->get(fields[i]));
		const size_t typesize = vars->typeToBytes(var->type());
		const unsigned int n_components = vars->typeToN(var->type());
		void *vect = malloc(typesize);
		if (!vect) {
			LOG(L_ERROR, std::string("Failure allocating ") +
				std::to_string(typesize) + "bytes.\n");
			throw std::bad_alloc();
		}
		for (size_t j = 0; j < n; j++) {
			const size_t offset = typesize * j;
			if (startswith(var->type(), "unsigned int") ||
			    startswith(var->type(), "uivec")) {
				memcpy(vect,
				       (char*)ptr + offset,
				       n_components * sizeof(vtkTypeUInt32));
				vtkSmartPointer<vtkTypeUInt32Array> vtk_array =
				    (vtkTypeUInt32Array*)(vtk_arrays[i].GetPointer());
				vtk_array->SetTypedTuple(j, (vtkTypeUInt32*)vect);
			} else if (startswith(var->type(), "unsigned long") ||
			           startswith(var->type(), "ulvec")) {
				memcpy(vect,
				       (char*)ptr + offset,
				       n_components * sizeof(vtkTypeUInt64));
				vtkSmartPointer<vtkTypeUInt64Array> vtk_array =
				    (vtkTypeUInt64Array*)(vtk_arrays[i].GetPointer());
				vtk_array->SetTypedTuple(j, (vtkTypeUInt64*)vect);
			} else if (startswith(var->type(), "int") ||
			           startswith(var->type(), "ivec")) {
				memcpy(vect,
				       (char*)ptr + offset,
				       n_components * sizeof(vtkTypeInt32));
				vtkSmartPointer<vtkTypeInt32Array> vtk_array =
				    (vtkTypeInt32Array*)(vtk_arrays[i].GetPointer());
				vtk_array->SetTypedTuple(j, (vtkTypeInt32*)vect);
			} else if (startswith(var->type(), "long") ||
			           startswith(var->type(), "lvec")) {
				memcpy(vect,
				       (char*)ptr + offset,
				       n_components * sizeof(vtkTypeInt64));
				vtkSmartPointer<vtkTypeInt64Array> vtk_array =
				    (vtkTypeInt64Array*)(vtk_arrays[i].GetPointer());
				vtk_array->SetTypedTuple(j, (vtkTypeInt64*)vect);
			} else if (startswith(var->type(), "float") ||
			           startswith(var->type(), "vec") ||
			           startswith(var->type(), "matrix")) {
				memcpy(vect,
				       (char*)ptr + offset,
				       n_components * sizeof(fcl));
				vtkSmartPointer<vtkFloatArray> vtk_array =
				    (vtkFloatArray*)(vtk_arrays[i].GetPointer());
				vtk_array->SetTypedTuple(j, (fcl*)vect);
			} else if (startswith(var->type(), "double") ||
			           startswith(var->type(), "dvec")) {
				memcpy(vect,
				       (char*)ptr + offset,
				       n_components * sizeof(dcl));
				vtkSmartPointer<vtkDoubleArray> vtk_array =
				    (vtkDoubleArray*)(vtk_arrays[i].GetPointer());
				vtk_array->SetTypedTuple(j, (dcl*)vect);
			}
		}
		free(vect);
	}

	// Setup the unstructured grid
	vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
	grid->SetPoints(vtk_points);
	grid->SetCells(vtk_vertex->GetCellType(), vtk_cells);
	for (unsigned int i = 0; i < fields.size(); i++) {
		if (!fields.at(i).compare("r")) {
			continue;
		}

		ArrayVariable* var = (ArrayVariable*)(vars->get(fields.at(i)));
		if (startswith(var->type(), "unsigned int") ||
		    startswith(var->type(), "uivec")) {
			vtkSmartPointer<vtkTypeUInt32Array> vtk_array =
			    (vtkTypeUInt32Array*)(vtk_arrays.at(i).GetPointer());
			grid->GetPointData()->AddArray(vtk_array);
		} else if (startswith(var->type(), "unsigned long") ||
		           startswith(var->type(), "ulvec")) {
			vtkSmartPointer<vtkTypeUInt64Array> vtk_array =
			    (vtkTypeUInt64Array*)(vtk_arrays.at(i).GetPointer());
			grid->GetPointData()->AddArray(vtk_array);
		} else if (startswith(var->type(), "int") ||
		           startswith(var->type(), "ivec")) {
			vtkSmartPointer<vtkTypeInt32Array> vtk_array =
			    (vtkTypeInt32Array*)(vtk_arrays.at(i).GetPointer());
			grid->GetPointData()->AddArray(vtk_array);
		} else if (startswith(var->type(), "long") ||
		           startswith(var->type(), "lvec")) {
			vtkSmartPointer<vtkTypeInt64Array> vtk_array =
			    (vtkTypeInt64Array*)(vtk_arrays.at(i).GetPointer());
			grid->GetPointData()->AddArray(vtk_array);
		} else if (startswith(var->type(), "float") ||
		           startswith(var->type(), "vec") ||
		           startswith(var->type(), "matrix")) {
			vtkSmartPointer<vtkFloatArray> vtk_array =
			    (vtkFloatArray*)(vtk_arrays.at(i).GetPointer());
			grid->GetPointData()->AddArray(vtk_array);
		} else if (startswith(var->type(), "double") ||
		           startswith(var->type(), "dvec")) {
			vtkSmartPointer<vtkDoubleArray> vtk_array =
			    (vtkDoubleArray*)(vtk_arrays.at(i).GetPointer());
			grid->GetPointData()->AddArray(vtk_array);
		}
	}

	return grid;
}

void
VTK::editVTK()
{
	unsigned int i, j;
	auto C = CalcServer::CalcServer::singleton();
	auto vars = C->variables();
	auto fields = simData().sets.at(setId())->outputFields();
	const unsigned int n = bounds().y - bounds().x;

	for (i = 0; i < fields.size(); i++) {
		void* ptr = data()[fields[i]];

		if (!fields[i].compare("r")) {
			// The mesh points are a bit of a special case
			auto points = _vtk->GetPoints();
			for (j = 0; j < n; j++) {
				if (C->have_3d()) {
					vec4* coords = (vec4*)ptr;
					points->SetPoint(j, coords[j].x, coords[j].y, coords[j].z);
				} else {
					vec2* coords = (vec2*)ptr;
					points->SetPoint(j, coords[j].x, coords[j].y, 0.0);
				}
			}
			// _vtk->SetPoints(points);
			points->Modified();
			continue;
		}

		auto array = _vtk->GetPointData()->GetAbstractArray(fields[i].c_str());
		auto var = (ArrayVariable*)(vars->get(fields[i]));
		const size_t typesize = vars->typeToBytes(var->type());
		for (j = 0; j < n; j++) {
			const size_t offset = typesize * j;
			auto shifted = (char*)ptr + offset;
			if (var->type().find("unsigned int") != std::string::npos ||
			    var->type().find("uivec") != std::string::npos) {
				auto vtk_array = (vtkUnsignedIntArray*)array;
				vtk_array->SetTypedTuple(j, (unsigned int*)shifted);
			} else if (var->type().find("int") != std::string::npos ||
			           var->type().find("ivec") != std::string::npos) {
				auto vtk_array = (vtkIntArray*)array;
				vtk_array->SetTypedTuple(j, (int*)shifted);
			} else if (var->type().find("float") != std::string::npos ||
			           var->type().find("vec") != std::string::npos ||
			           var->type().find("matrix") != std::string::npos) {
				auto vtk_array = (vtkFloatArray*)array;
				vtk_array->SetTypedTuple(j, (float*)shifted);
			}
		}
		array->Modified();
	}

	_vtk->Modified();
}

void
VTK::save(float t)
{
	// Check the fields to write
	std::vector<std::string> fields =
	    simData().sets.at(setId())->outputFields();
	bool have_r = false;
	for (auto field : fields) {
		if (!field.compare("r")) {
			have_r = true;
			break;
		}
	}
	if (!have_r) {
		LOG(L_ERROR, "\"r\" field was not set to be saved into the file.\n");
		throw std::runtime_error("\"r\" field is mandatory");
	}

	Particles::save(t);
}

size_t
VTK::compute_n()
{
	vtkSmartPointer<vtkXMLUnstructuredGridReader> f =
	    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

	if (!f->CanReadFile(simData().sets.at(setId())->inputPath().c_str())) {
		std::ostringstream msg;
		msg << "Cannot load VTK file \""
		    << simData().sets.at(setId())->inputPath() << "\"!" << std::endl;
		LOG(L_ERROR, msg.str());
		throw std::runtime_error("Failure reading file");
	}

	f->SetFileName(simData().sets.at(setId())->inputPath().c_str());
	f->Update();

	vtkSmartPointer<vtkUnstructuredGrid> grid = f->GetOutput();

	return (size_t)grid->GetNumberOfPoints();
}

vtkXMLUnstructuredGridWriter*
VTK::create()
{
	vtkXMLUnstructuredGridWriter* f = NULL;

	std::string basename = simData().sets.at(setId())->outputPath();
	// Check that {index} scape string is present, for backward compatibility
	if (basename.find("{index}") == std::string::npos) {
		basename += ".{index}.vtu";
	}
	_next_file_index = file(basename, _next_file_index);

	std::ostringstream msg;
	msg << "Writing \"" << file() << "\" VTK file..." << std::endl;
	LOG(L_INFO, msg.str());

	f = vtkXMLUnstructuredGridWriter::New();
	basename = file();
	f->SetFileName(basename.c_str());
	_next_file_index++;

	return f;
}

void
VTK::updateSeries(float t)
{
	std::ostringstream msg;
	msg << "Writing \"" << filenameSeries() << "\" Paraview data file..."
	    << std::endl;
	LOG(L_INFO, msg.str());

	auto jfile = json::object();
	jfile["name"] = file();
	jfile["time"] = t;
	_data_series["files"].push_back(jfile);

	std::ofstream f(filenameSeries(), std::ios::out | std::ios::trunc);
	if (!f.is_open()) {
		std::ostringstream msg;
		msg << "Failure writing on '" << filenameSeries() << "'" << std::endl;
		LOG(L_ERROR, msg.str());
	}
	f << _data_series.dump(4) << std::endl;
	f.close();
}

const std::string
VTK::filenameSeries()
{
	const std::string ext(".vtu.series");
	if (_name_series == "") {
		try {
			unsigned int i = 0;
			_name_series = newFilePath(
			    simData().sets.at(setId())->outputPath() + ext, i, 1);
		} catch (std::invalid_argument& e) {
			std::ostringstream msg;
			_name_series =
			    setStrConstantsCopy(simData().sets.at(setId())->outputPath()) +
			    ext;
			msg << "Overwriting '" << _name_series << "'" << std::endl;
			LOG(L_WARNING, msg.str());
		}
	}
	return _name_series;
}

}
} // namespace

#endif // HAVE_VTK
