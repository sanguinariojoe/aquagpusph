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

#include <sphPrerequisites.h>

#ifdef HAVE_VTK

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>

#include <InputOutput/VTK.h>
#include <ScreenManager.h>
#include <ProblemSetup.h>
#include <TimeManager.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>

#include <vector>
#include <deque>
static std::deque<char*> cpp_str;
static std::deque<XMLCh*> xml_str;
static std::deque<xercesc::XercesDOMParser*> parsers;

static char *xmlTranscode(const XMLCh *txt)
{
    char *str = xercesc::XMLString::transcode(txt);
    cpp_str.push_back(str);
    return str;
}

static XMLCh *xmlTranscode(const char *txt)
{
    XMLCh *str = xercesc::XMLString::transcode(txt);
    xml_str.push_back(str);
    return str;
}

static void xmlClear()
{
    unsigned int i;
    for(i = 0; i < cpp_str.size(); i++){
        xercesc::XMLString::release(&cpp_str.at(i));
    }
    cpp_str.clear();
    for(i = 0; i < xml_str.size(); i++){
        xercesc::XMLString::release(&xml_str.at(i));
    }
    xml_str.clear();
    for(i = 0; i < parsers.size(); i++){
        delete parsers.at(i);
    }
    parsers.clear();
}

#ifdef xmlS
    #undef xmlS
#endif // xmlS
#define xmlS(txt) xmlTranscode(txt)

#ifdef xmlAttribute
    #undef xmlAttribute
#endif
#define xmlAttribute(elem, att) xmlS( elem->getAttribute(xmlS(att)) )

#ifdef xmlHasAttribute
    #undef xmlHasAttribute
#endif
#define xmlHasAttribute(elem, att) elem->hasAttribute(xmlS(att))

#ifndef REQUESTED_FIELDS
    #ifdef HAVE_3D
        #define REQUESTED_FIELDS 17
    #else
        #define REQUESTED_FIELDS 13
    #endif
#endif // REQUESTED_FIELDS

using namespace xercesc;

namespace Aqua{ namespace InputOutput{

VTK::VTK(ProblemSetup sim_data,
         unsigned int first,
         unsigned int n,
         unsigned int iset)
    : Particles(sim_data, first, n, iset)
    , _namePVD(NULL)
    , _next_file_index(0)
{
}

VTK::~VTK()
{
    unsigned int i;
    ScreenManager *S = ScreenManager::singleton();

    // Wait for the writers working
    S->addMessageF(L_INFO, "Waiting for the writers...\n");
    for(i = 0; i < _tids.size(); i++){
        pthread_join(_tids.at(i), NULL);
    }
    _tids.clear();

    if(_namePVD) delete[] _namePVD; _namePVD=NULL;
}

void VTK::load()
{
    unsigned int i, j, k, n, N, progress;
    int aux;
    cl_int err_code;
    char msg[1024];
    ScreenManager *S = ScreenManager::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    loadDefault();

    sprintf(msg,
            "Loading particles from VTK file \"%s\"\n",
            simData().sets.at(setId())->inputPath().c_str());
    S->addMessageF(L_INFO, msg);

    vtkSmartPointer<vtkXMLUnstructuredGridReader> f =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    if(!f->CanReadFile(simData().sets.at(setId())->inputPath().c_str())){
        S->addMessageF(L_ERROR, "The file cannot be read.\n");
        throw std::runtime_error("Failure reading file");
    }

    f->SetFileName(simData().sets.at(setId())->inputPath().c_str());
    f->Update();

    vtkSmartPointer<vtkUnstructuredGrid> grid = f->GetOutput();

    // Assert that the number of particles is right
    n = bounds().y - bounds().x;
    N = (unsigned int)grid->GetNumberOfPoints();
    if( n != N){
        sprintf(msg,
                "Expected %u particles, but the file contains %u ones.\n",
                n,
                N);
        S->addMessage(L_ERROR, msg);
        throw std::runtime_error("Invalid number of particles");
    }

    // Check the fields to read
    std::vector<std::string> fields = simData().sets.at(setId())->inputFields();
    if(!fields.size()){
        S->addMessage(L_ERROR, "0 fields were set to be read from the file.\n");
        throw std::runtime_error("No fields have been marked to read");
    }
    bool have_r = false;
    for(i = 0; i < fields.size(); i++){
        if(!strcmp(fields.at(i).c_str(), "r")){
            have_r = true;
            break;
        }
    }
    if(!have_r){
        S->addMessage(L_ERROR, "\"r\" field was not set to be read from the file.\n");
        throw std::runtime_error("\"r\" field is mandatory");
    }

    // Setup an storage
    std::deque<void*> data;
    Variables* vars = C->variables();
    for(i = 0; i < fields.size(); i++){
        if(!vars->get(fields.at(i).c_str())){
            sprintf(msg,
                    "\"%s\" field has been set to be read, but it was not declared.\n",
                    fields.at(i).c_str());
            S->addMessage(L_ERROR, msg);
            throw std::runtime_error("Invalid field");
        }
        if(!strchr(vars->get(fields.at(i).c_str())->type(), '*')){
            sprintf(msg,
                    "\"%s\" field has been set to be read, but it was declared as an scalar.\n",
                    fields.at(i).c_str());
            S->addMessage(L_ERROR, msg);
            throw std::runtime_error("Invalid field type");
        }
        ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(i).c_str());
        size_t typesize = vars->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < bounds().y){
            sprintf(msg,
                    "Failure reading \"%s\" field, which has not length enough.\n",
                    fields.at(i).c_str());
            S->addMessage(L_ERROR, msg);
            throw std::runtime_error("Invalid field length");
        }
        void *store = malloc(typesize * n);
        if(!store){
            sprintf(msg,
                    "Failure allocating memory for \"%s\" field.\n",
                    fields.at(i).c_str());
            S->addMessage(L_ERROR, msg);
            throw std::bad_alloc();
        }
        data.push_back(store);
    }

    progress = -1;
    vtkSmartPointer<vtkPoints> vtk_points = grid->GetPoints();
    vtkSmartPointer<vtkPointData> vtk_data = grid->GetPointData();
    for(i = 0; i < n; i++){
        for(j = 0; j < fields.size(); j++){
            if(!fields.at(j).compare("r")){
                double *vect = vtk_points->GetPoint(i);
                vec *ptr = (vec*)data.at(j);
                ptr[i].x = vect[0];
                ptr[i].y = vect[1];
                #ifdef HAVE_3D
                    ptr[i].z = vect[2];
                    ptr[i].w = 0.f;
                #endif
                continue;
            }
            ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(j).c_str());
            size_t type_size = vars->typeToBytes(var->type());
            unsigned int n_components = vars->typeToN(var->type());
            if(strstr(var->type(), "unsigned int") ||
               strstr(var->type(), "uivec")){
                vtkSmartPointer<vtkUnsignedIntArray> vtk_array =
                    (vtkUnsignedIntArray*)(vtk_data->GetArray(fields.at(j).c_str(), aux));
                for(k = 0; k < n_components; k++){
                    unsigned int component = vtk_array->GetComponent(i, k);
                    size_t offset = type_size * i + sizeof(unsigned int) * k;
                    memcpy((char*)data.at(j) + offset,
                           &component,
                           sizeof(unsigned int));
                }
            }
            else if(strstr(var->type(), "int") ||
                    strstr(var->type(), "ivec")){
                vtkSmartPointer<vtkIntArray> vtk_array =
                    (vtkIntArray*)(vtk_data->GetArray(fields.at(j).c_str(), aux));
                for(k = 0; k < n_components; k++){
                    int component = vtk_array->GetComponent(i, k);
                    size_t offset = type_size * i + sizeof(int) * k;
                    memcpy((char*)data.at(j) + offset,
                           &component,
                           sizeof(int));
                }
            }
            else if(strstr(var->type(), "float") ||
                    strstr(var->type(), "vec") ||
                    strstr(var->type(), "matrix")){
                vtkSmartPointer<vtkFloatArray> vtk_array =
                    (vtkFloatArray*)(vtk_data->GetArray(fields.at(j).c_str(), aux));
                for(k = 0; k < n_components; k++){
                    float component = vtk_array->GetComponent(i, k);
                    size_t offset = type_size * i + sizeof(float) * k;
                    memcpy((char*)data.at(j) + offset,
                           &component,
                           sizeof(float));
                }
            }
        }
        if(progress != i * 100 / n){
            progress = i * 100 / n;
            if(!(progress % 10)){
                sprintf(msg, "\t\t%u%%\n", progress);
                S->addMessage(L_DEBUG, msg);
            }
        }
    }

    // Send the data to the server and release it
    for(i = 0; i < fields.size(); i++){
        ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(i).c_str());
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
        free(data.at(i)); data.at(i) = NULL;
        if(err_code != CL_SUCCESS){
            sprintf(msg,
                    "Failure sending variable \"%s\" to the server.\n",
                    fields.at(i).c_str());
            S->addMessageF(L_ERROR, msg);
            S->printOpenCLError(err_code);
        }
    }
    data.clear();
}

/** @brief Data structure to send the data to a parallel writer thread
 */
typedef struct{
    /// The field names
    std::vector<std::string> fields;
    /// Bounds of the particles index managed by this writer
    uivec2 bounds;
    /// Screen manager
    ScreenManager *S;
    /// VTK arrays
    CalcServer::CalcServer *C;
    /// The data associated to each field
    std::deque<void*> data;
    /// The VTK file decriptor
    vtkXMLUnstructuredGridWriter *f;
}data_pthread;

/** @brief Parallel thread to write the data
 * @param data_void Input data of type data_pthread* (dynamically casted as
 * void*)
 */
void* save_pthread(void *data_void)
{
    unsigned int i, j;
    char msg[1024];
    data_pthread *data = (data_pthread*)data_void;

    // Create storage arrays
    std::deque< vtkSmartPointer<vtkDataArray> > vtk_arrays;
    Variables* vars = data->C->variables();
    for(i = 0; i < data->fields.size(); i++){
        if(!vars->get(data->fields.at(i).c_str())){
            sprintf(msg,
                    "\"%s\" field has been set to be saved, but it was not declared.\n",
                    data->fields.at(i).c_str());
            data->S->addMessage(L_ERROR, msg);
            for(i = 0; i < data->fields.size(); i++){
                free(data->data.at(i)); data->data.at(i) = NULL;
            }
            data->data.clear();
            data->f->Delete();
            delete data; data=NULL;
            return NULL;
        }
        if(!strchr(vars->get(data->fields.at(i).c_str())->type(), '*')){
            sprintf(msg,
                    "\"%s\" field has been set to be saved, but it was declared as a scalar.\n",
                    data->fields.at(i).c_str());
            data->S->addMessage(L_ERROR, msg);
            for(i = 0; i < data->fields.size(); i++){
                free(data->data.at(i)); data->data.at(i) = NULL;
            }
            data->data.clear();
            data->f->Delete();
            delete data; data=NULL;
            return NULL;
        }
        ArrayVariable *var = (ArrayVariable*)vars->get(data->fields.at(i).c_str());
        size_t typesize = vars->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < data->bounds.y){
            sprintf(msg,
                    "Failure saving \"%s\" field, which has not length enough.\n",
                    data->fields.at(i).c_str());
            data->S->addMessage(L_ERROR, msg);
            for(i = 0; i < data->fields.size(); i++){
                free(data->data.at(i)); data->data.at(i) = NULL;
            }
            data->data.clear();
            data->f->Delete();
            delete data; data=NULL;
            return NULL;
        }

        unsigned int n_components = vars->typeToN(var->type());
        if(strstr(var->type(), "unsigned int") ||
           strstr(var->type(), "uivec")){
            vtkSmartPointer<vtkUnsignedIntArray> vtk_array =
                vtkSmartPointer<vtkUnsignedIntArray>::New();
            vtk_array->SetNumberOfComponents(n_components);
            vtk_array->SetName(data->fields.at(i).c_str());
            vtk_arrays.push_back(vtk_array);
        }
        else if(strstr(var->type(), "int") ||
                strstr(var->type(), "ivec")){
            vtkSmartPointer<vtkIntArray> vtk_array =
                vtkSmartPointer<vtkIntArray>::New();
            vtk_array->SetNumberOfComponents(n_components);
            vtk_array->SetName(data->fields.at(i).c_str());
            vtk_arrays.push_back(vtk_array);
        }
        else if(strstr(var->type(), "float") ||
                strstr(var->type(), "vec") ||
                strstr(var->type(), "matrix")){
            vtkSmartPointer<vtkFloatArray> vtk_array =
                vtkSmartPointer<vtkFloatArray>::New();
            vtk_array->SetNumberOfComponents(n_components);
            vtk_array->SetName(data->fields.at(i).c_str());
            vtk_arrays.push_back(vtk_array);
        }
    }

    vtkSmartPointer<vtkVertex> vtk_vertex;
    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> vtk_cells = vtkSmartPointer<vtkCellArray>::New();

    for(i = 0; i < data->bounds.y - data->bounds.x; i++){
        for(j = 0; j < data->fields.size(); j++){
            if(!strcmp(data->fields.at(j).c_str(), "r")){
                vec *ptr = (vec*)(data->data.at(j));
                #ifdef HAVE_3D
                    vtk_points->InsertNextPoint(ptr[i].x, ptr[i].y, ptr[i].z);
                #else
                    vtk_points->InsertNextPoint(ptr[i].x, ptr[i].y, 0.f);
                #endif
                continue;
            }
            ArrayVariable *var = (ArrayVariable*)(
                vars->get(data->fields.at(j).c_str()));
            size_t typesize = vars->typeToBytes(var->type());
            unsigned int n_components = vars->typeToN(var->type());
            if(strstr(var->type(), "unsigned int") ||
               strstr(var->type(), "uivec")){
                unsigned int vect[n_components];
                size_t offset = typesize * i;
                memcpy(vect,
                       (char*)(data->data.at(j)) + offset,
                       n_components * sizeof(unsigned int));
                vtkSmartPointer<vtkUnsignedIntArray> vtk_array =
                    (vtkUnsignedIntArray*)(vtk_arrays.at(j).GetPointer());
                #if VTK_MAJOR_VERSION < 7
                    vtk_array->InsertNextTupleValue(vect);
                #else
                    vtk_array->InsertNextTypedTuple(vect);
                #endif // VTK_MAJOR_VERSION
            }
            else if(strstr(var->type(), "int") ||
                    strstr(var->type(), "ivec")){
                int vect[n_components];
                size_t offset = typesize * i;
                memcpy(vect,
                       (char*)(data->data.at(j)) + offset,
                       n_components * sizeof(int));
                vtkSmartPointer<vtkIntArray> vtk_array =
                    (vtkIntArray*)(vtk_arrays.at(j).GetPointer());
                #if VTK_MAJOR_VERSION < 7
                    vtk_array->InsertNextTupleValue(vect);
                #else
                    vtk_array->InsertNextTypedTuple(vect);
                #endif // VTK_MAJOR_VERSION
            }
            else if(strstr(var->type(), "float") ||
                    strstr(var->type(), "vec") ||
                    strstr(var->type(), "matrix")){
                float vect[n_components];
                size_t offset = typesize * i;
                memcpy(vect,
                       (char*)(data->data.at(j)) + offset,
                       n_components * sizeof(float));
                vtkSmartPointer<vtkFloatArray> vtk_array =
                    (vtkFloatArray*)(vtk_arrays.at(j).GetPointer());
                #if VTK_MAJOR_VERSION < 7
                    vtk_array->InsertNextTupleValue(vect);
                #else
                    vtk_array->InsertNextTypedTuple(vect);
                #endif // VTK_MAJOR_VERSION
            }
        }
        vtk_vertex = vtkSmartPointer<vtkVertex>::New();
        vtk_vertex->GetPointIds()->SetId(0, i);
        vtk_cells->InsertNextCell(vtk_vertex);
    }

    // Setup the unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(vtk_points);
    grid->SetCells(vtk_vertex->GetCellType(), vtk_cells);
    for(i = 0; i < data->fields.size(); i++){
        if(!data->fields.at(i).compare("r")){
            continue;
        }

        ArrayVariable *var = (ArrayVariable*)(
            vars->get(data->fields.at(i).c_str()));
        if(strstr(var->type(), "unsigned int") ||
           strstr(var->type(), "uivec")){
            vtkSmartPointer<vtkUnsignedIntArray> vtk_array =
                (vtkUnsignedIntArray*)(vtk_arrays.at(i).GetPointer());
            grid->GetPointData()->AddArray(vtk_array);
        }
        else if(strstr(var->type(), "int") ||
                strstr(var->type(), "ivec")){
            vtkSmartPointer<vtkIntArray> vtk_array =
                (vtkIntArray*)(vtk_arrays.at(i).GetPointer());
            grid->GetPointData()->AddArray(vtk_array);
        }
        else if(strstr(var->type(), "float") ||
                strstr(var->type(), "vec") ||
                strstr(var->type(), "matrix")){
            vtkSmartPointer<vtkFloatArray> vtk_array =
                (vtkFloatArray*)(vtk_arrays.at(i).GetPointer());
            grid->GetPointData()->AddArray(vtk_array);
        }
    }

    // Write file
    #if VTK_MAJOR_VERSION <= 5
        data->f->SetInput(grid);
    #else // VTK_MAJOR_VERSION
        data->f->SetInputData(grid);
    #endif // VTK_MAJOR_VERSION

    if(!data->f->Write()){
        data->S->addMessageF(L_ERROR, "Failure writing the VTK file.\n");
    }

    // Clean up
    for(i = 0; i < data->fields.size(); i++){
        free(data->data.at(i)); data->data.at(i) = NULL;
    }
    data->data.clear();
    data->f->Delete();
    delete data; data=NULL;
    return NULL;
}

void VTK::save()
{
    unsigned int i;
    ScreenManager *S = ScreenManager::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    // Check the fields to write
    std::vector<std::string> fields = simData().sets.at(setId())->outputFields();
    if(!fields.size()){
        S->addMessage(L_ERROR, "0 fields were set to be saved into the file.\n");
        throw std::runtime_error("No fields have been marked to be saved");
    }
    bool have_r = false;
    for(i = 0; i < fields.size(); i++){
        if(!fields.at(i).compare("r")){
            have_r = true;
            break;
        }
    }
    if(!have_r){
        S->addMessage(L_ERROR, "\"r\" field was not set to be saved into the file.\n");
        throw std::runtime_error("\"r\" field is mandatory");
    }

    // Setup the data struct for the parallel thread
    data_pthread *data = new data_pthread;
    data->fields = fields;
    data->bounds = bounds();
    data->C = C;
    data->S = S;
    data->data = download(fields);
    if(!data->data.size()){
        throw std::runtime_error("Failure downloading data");
    }
    data->f = create();
    if(!data->f){
        throw std::runtime_error("Failure creating file");
    }

    // Launch the thread
    pthread_t tid;
    int err;
    err = pthread_create(&tid, NULL, &save_pthread, (void*)data);
    if(err){
        S->addMessageF(L_ERROR, "Failure launching the parallel thread.\n");
        char err_str[strlen(strerror(err)) + 2];
        strcpy(err_str, strerror(err));
        strcat(err_str, "\n");
        S->addMessage(L_DEBUG, err_str);
        throw std::runtime_error("Failure launching saving thread");
    }
    _tids.push_back(tid);

    // Clear the already finished threads
    if(_tids.size() > 0){
        i = _tids.size() - 1;
        while(true){
            if(pthread_kill(_tids.at(i), 0)){
                pthread_join(_tids.at(i), NULL);
                _tids.erase(_tids.begin() + i);
            }
            if(i == 0) break;
            i--;
        }
    }

    // Check and limit the number of active writing processes
    if(_tids.size() > 2){
        S->addMessageF(L_WARNING, "More than 2 active writing tasks\n");
        S->addMessageF(L_DEBUG, "This may result in heavy performance penalties, and hard disk failures\n");
        S->addMessageF(L_DEBUG, "Please, consider a reduction of the output printing rate\n");
        while(_tids.size() > 2){
            pthread_join(_tids.at(0), NULL);
            _tids.erase(_tids.begin());
        }
    }

    // Update the PVD file
    if(updatePVD()){
        throw std::runtime_error("Failure updating PVD description file");
    }
}

vtkXMLUnstructuredGridWriter* VTK::create(){
    char *basename, msg[1024];
    size_t len;
    vtkXMLUnstructuredGridWriter *f = NULL;
    ScreenManager *S = ScreenManager::singleton();

    // Create the file base name
    len = strlen(simData().sets.at(setId())->outputPath().c_str()) + 8;
    basename = new char[len];
    strcpy(basename, simData().sets.at(setId())->outputPath().c_str());
    strcat(basename, ".%d.vtu");

    _next_file_index = file(basename, _next_file_index);
    if(!_next_file_index){
        delete[] basename;
        S->addMessageF(L_ERROR, "Failure getting a valid filename.\n");
        S->addMessage(L_DEBUG, "\tHow do you received this message?.\n");
        return NULL;
    }
    delete[] basename;

    sprintf(msg, "Writing \"%s\" VTK output...\n", file());
    S->addMessageF(L_INFO, msg);

    f = vtkXMLUnstructuredGridWriter::New();
    f->SetFileName(file());

    return f;
}

bool VTK::updatePVD(){
    unsigned int n;
    char msg[1024];
    ScreenManager *S = ScreenManager::singleton();
    TimeManager *T = TimeManager::singleton();

    sprintf(msg, "Writing \"%s\" Paraview data file...\n", filenamePVD());
    S->addMessageF(L_INFO, msg);

    bool should_release_doc = false;
    DOMDocument* doc = getPVD(false);
    if(!doc){
        should_release_doc = true;
        doc = getPVD(true);
    }
    DOMElement* root = doc->getDocumentElement();
    if(!root){
        S->addMessageF(L_ERROR, "Empty XML file found!\n");
        return true;
    }
    n = doc->getElementsByTagName(xmlS("VTKFile"))->getLength();
    if(n != 1){
        sprintf(msg,
                "Expected 1 VTKFile root section, but %u has been found\n",
                n);
        S->addMessageF(L_ERROR, msg);
        return true;
    }

    DOMNodeList* nodes = root->getElementsByTagName(xmlS("Collection"));
    if(nodes->getLength() != 1){
        sprintf(msg,
                "Expected 1 collection, but %u has been found\n",
                nodes->getLength());
        S->addMessageF(L_ERROR, msg);
    }
    DOMNode* node = nodes->item(0);
    DOMElement* elem = dynamic_cast<xercesc::DOMElement*>(node);

    DOMElement *s_elem;
    s_elem = doc->createElement(xmlS("DataSet"));
    sprintf(msg, "%g", T->time());
    s_elem->setAttribute(xmlS("timestep"), xmlS(msg));
    s_elem->setAttribute(xmlS("group"), xmlS(""));
    s_elem->setAttribute(xmlS("part"), xmlS("0"));
    s_elem->setAttribute(xmlS("file"), xmlS(file()));
    elem->appendChild(s_elem);

    // Save the XML document to a file
    DOMImplementation* impl;
    DOMLSSerializer* saver;
    impl = DOMImplementationRegistry::getDOMImplementation(xmlS("LS"));
    saver = ((DOMImplementationLS*)impl)->createLSSerializer();

    if(saver->getDomConfig()->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
        saver->getDomConfig()->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
    saver->setNewLine(xmlS("\r\n"));

    XMLFormatTarget *target = new LocalFileFormatTarget(filenamePVD());
    // XMLFormatTarget *target = new StdOutFormatTarget();
    DOMLSOutput *output = ((DOMImplementationLS*)impl)->createLSOutput();
    output->setByteStream(target);
    output->setEncoding(xmlS("UTF-8"));

    try {
        saver->write(doc, output);
    }
    catch( XMLException& e ){
        char* message = xmlS(e.getMessage());
        S->addMessageF(L_ERROR, "XML toolkit writing error.\n");
        sprintf(msg, "\t%s\n", message);
        S->addMessage(L_DEBUG, msg);
        xmlClear();
        return true;
    }
    catch( DOMException& e ){
        char* message = xmlS(e.getMessage());
        S->addMessageF(L_ERROR, "XML DOM writing error.\n");
        sprintf(msg, "\t%s\n", message);
        S->addMessage(L_DEBUG, msg);
        xmlClear();
        return true;
    }
    catch( ... ){
        S->addMessageF(L_ERROR, "Writing error.\n");
        S->addMessage(L_DEBUG, "\tUnhandled exception\n");
        xmlClear();
        return true;
    }

    target->flush();

    delete target;
    saver->release();
    output->release();
    if(should_release_doc)
        doc->release();
    xmlClear();

    return false;
}

DOMDocument* VTK::getPVD(bool generate)
{
    DOMDocument* doc = NULL;
    FILE *dummy=NULL;

    // Try to open as ascii file, just to know if the file already exist
    dummy = fopen(filenamePVD(), "r");
    if(!dummy){
        if(!generate){
            return NULL;
        }
        DOMImplementation* impl;
        impl = DOMImplementationRegistry::getDOMImplementation(xmlS("Range"));
        DOMDocument* doc = impl->createDocument(
            NULL,
            xmlS("VTKFile"),
            NULL);
        DOMElement* root = doc->getDocumentElement();
        root->setAttribute(xmlS("type"), xmlS("Collection"));
        root->setAttribute(xmlS("version"), xmlS("0.1"));
        DOMElement *elem;
        elem = doc->createElement(xmlS("Collection"));
        root->appendChild(elem);
        return doc;
    }
    fclose(dummy);
    XercesDOMParser *parser = new XercesDOMParser();
    parser->setValidationScheme(XercesDOMParser::Val_Never);
    parser->setDoNamespaces(false);
    parser->setDoSchema(false);
    parser->setLoadExternalDTD(false);
    parser->parse(filenamePVD());
    doc = parser->getDocument();
    parsers.push_back(parser);
    return doc;
}

const char* VTK::filenamePVD()
{
    if(!_namePVD){
        size_t len = strlen(simData().sets.at(setId())->outputPath().c_str()) + 5;
        _namePVD = new char[len];
        strcpy(_namePVD, simData().sets.at(setId())->outputPath().c_str());
        strcat(_namePVD, ".pvd");
    }
    return (const char*)_namePVD;
}

}}  // namespace

#endif // HAVE_VTK
