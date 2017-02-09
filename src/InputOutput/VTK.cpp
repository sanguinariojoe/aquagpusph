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

VTK::VTK(unsigned int first, unsigned int n, unsigned int iset)
    : Particles(first, n, iset)
    , _namePVD(NULL)
    , _next_file_index(0)
{
}

VTK::~VTK()
{
    unsigned int i;
    ScreenManager *S = ScreenManager::singleton();

    // Wait for the writers working
    S->addMessageF(1, "Waiting for the writers...\n");
    for(i = 0; i < _tids.size(); i++){
        if(!pthread_kill(_tids.at(i), 0))
        {
            pthread_join(_tids.at(i), NULL);
        }
    }
    _tids.clear();

    if(_namePVD) delete[] _namePVD; _namePVD=NULL;
}

bool VTK::load()
{
    unsigned int i, j, k, n, N, progress;
    int aux;
    cl_int err_code;
    char msg[1024];
    ScreenManager *S = ScreenManager::singleton();
    ProblemSetup *P = ProblemSetup::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    loadDefault();

    sprintf(msg,
            "Loading particles from VTK file \"%s\"\n",
            P->sets.at(setId())->inputPath());
    S->addMessageF(1, msg);

    vtkSmartPointer<vtkXMLUnstructuredGridReader> f =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    if(!f->CanReadFile(P->sets.at(setId())->inputPath())){
        S->addMessageF(3, "The file cannot be read.\n");
        return true;
    }

    f->SetFileName(P->sets.at(setId())->inputPath());
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
        S->addMessage(3, msg);
        return true;
    }

    // Check the fields to read
    std::deque<char*> fields = P->sets.at(setId())->inputFields();
    if(!fields.size()){
        S->addMessage(3, "0 fields were set to be read from the file.\n");
        return true;
    }
    bool have_r = false;
    for(i = 0; i < fields.size(); i++){
        if(!strcmp(fields.at(i), "r")){
            have_r = true;
            break;
        }
    }
    if(!have_r){
        S->addMessage(3, "\"r\" field was not set to be read from the file.\n");
        return true;
    }

    // Setup an storage
    std::deque<void*> data;
    Variables* vars = C->variables();
    for(i = 0; i < fields.size(); i++){
        if(!vars->get(fields.at(i))){
            sprintf(msg,
                    "\"%s\" field has been set to be read, but it was not declared.\n",
                    fields.at(i));
            S->addMessage(3, msg);
            return true;
        }
        if(!strchr(vars->get(fields.at(i))->type(), '*')){
            sprintf(msg,
                    "\"%s\" field has been set to be read, but it was declared as an scalar.\n",
                    fields.at(i));
            S->addMessage(3, msg);
            return true;
        }
        ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(i));
        size_t typesize = vars->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < bounds().y){
            sprintf(msg,
                    "Failure reading \"%s\" field, which has not length enough.\n",
                    fields.at(i));
            S->addMessage(3, msg);
            return true;
        }
        void *store = malloc(typesize * n);
        if(!store){
            sprintf(msg,
                    "Failure allocating memory for \"%s\" field.\n",
                    fields.at(i));
            S->addMessage(3, msg);
            return true;
        }
        data.push_back(store);
    }

    progress = -1;
    vtkSmartPointer<vtkPoints> vtk_points = grid->GetPoints();
    vtkSmartPointer<vtkPointData> vtk_data = grid->GetPointData();
    for(i = 0; i < n; i++){
        for(j = 0; j < fields.size(); j++){
            if(!strcmp(fields.at(j), "r")){
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
            ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(j));
            size_t type_size = vars->typeToBytes(var->type());
            unsigned int n_components = vars->typeToN(var->type());
            if(strstr(var->type(), "unsigned int") ||
               strstr(var->type(), "uivec")){
                vtkSmartPointer<vtkUnsignedIntArray> vtk_array =
                    (vtkUnsignedIntArray*)(vtk_data->GetArray(fields.at(j), aux));
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
                    (vtkIntArray*)(vtk_data->GetArray(fields.at(j), aux));
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
                    (vtkFloatArray*)(vtk_data->GetArray(fields.at(j), aux));
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
                S->addMessage(0, msg);
            }
        }
    }

    // Send the data to the server and release it
    for(i = 0; i < fields.size(); i++){
        ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(i));
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
                    fields.at(i));
            S->addMessageF(3, msg);
            S->printOpenCLError(err_code);
        }
    }
    data.clear();

    return false;
}

/** @brief Data structure to send the data to a parallel writer thread
 */
typedef struct{
    /// The field names
    std::deque<char*> fields;
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
        if(!vars->get(data->fields.at(i))){
            sprintf(msg,
                    "\"%s\" field has been set to be saved, but it was not declared.\n",
                    data->fields.at(i));
            data->S->addMessage(3, msg);
            for(i = 0; i < data->fields.size(); i++){
                free(data->data.at(i)); data->data.at(i) = NULL;
            }
            data->data.clear();
            data->f->Delete();
            delete data; data=NULL;
            return NULL;
        }
        if(!strchr(vars->get(data->fields.at(i))->type(), '*')){
            sprintf(msg,
                    "\"%s\" field has been set to be saved, but it was declared as a scalar.\n",
                    data->fields.at(i));
            data->S->addMessage(3, msg);
            for(i = 0; i < data->fields.size(); i++){
                free(data->data.at(i)); data->data.at(i) = NULL;
            }
            data->data.clear();
            data->f->Delete();
            delete data; data=NULL;
            return NULL;
        }
        ArrayVariable *var = (ArrayVariable*)vars->get(data->fields.at(i));
        size_t typesize = vars->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < data->bounds.y){
            sprintf(msg,
                    "Failure saving \"%s\" field, which has not length enough.\n",
                    data->fields.at(i));
            data->S->addMessage(3, msg);
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
            vtk_array->SetName(data->fields.at(i));
            vtk_arrays.push_back(vtk_array);
        }
        else if(strstr(var->type(), "int") ||
                strstr(var->type(), "ivec")){
            vtkSmartPointer<vtkIntArray> vtk_array =
                vtkSmartPointer<vtkIntArray>::New();
            vtk_array->SetNumberOfComponents(n_components);
            vtk_array->SetName(data->fields.at(i));
            vtk_arrays.push_back(vtk_array);
        }
        else if(strstr(var->type(), "float") ||
                strstr(var->type(), "vec") ||
                strstr(var->type(), "matrix")){
            vtkSmartPointer<vtkFloatArray> vtk_array =
                vtkSmartPointer<vtkFloatArray>::New();
            vtk_array->SetNumberOfComponents(n_components);
            vtk_array->SetName(data->fields.at(i));
            vtk_arrays.push_back(vtk_array);
        }
    }

    vtkSmartPointer<vtkVertex> vtk_vertex;
    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> vtk_cells = vtkSmartPointer<vtkCellArray>::New();

    for(i = 0; i < data->bounds.y - data->bounds.x; i++){
        for(j = 0; j < data->fields.size(); j++){
            if(!strcmp(data->fields.at(j), "r")){
                vec *ptr = (vec*)(data->data.at(j));
                #ifdef HAVE_3D
                    vtk_points->InsertNextPoint(ptr[i].x, ptr[i].y, ptr[i].z);
                #else
                    vtk_points->InsertNextPoint(ptr[i].x, ptr[i].y, 0.f);
                #endif
                continue;
            }
            ArrayVariable *var = (ArrayVariable*)(
                vars->get(data->fields.at(j)));
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
                vtk_array->InsertNextTupleValue(vect);
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
                vtk_array->InsertNextTupleValue(vect);
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
                vtk_array->InsertNextTupleValue(vect);
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
        if(!strcmp(data->fields.at(i), "r")){
            continue;
        }

        ArrayVariable *var = (ArrayVariable*)(
            vars->get(data->fields.at(i)));
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
        data->S->addMessageF(3, "Failure writing the VTK file.\n");
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

bool VTK::save()
{
    unsigned int i;
    ScreenManager *S = ScreenManager::singleton();
    ProblemSetup *P = ProblemSetup::singleton();
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    // Check the fields to write
    std::deque<char*> fields = P->sets.at(setId())->outputFields();
    if(!fields.size()){
        S->addMessage(3, "0 fields were set to be saved into the file.\n");
        return true;
    }
    bool have_r = false;
    for(i = 0; i < fields.size(); i++){
        if(!strcmp(fields.at(i), "r")){
            have_r = true;
            break;
        }
    }
    if(!have_r){
        S->addMessage(3, "\"r\" field was not set to be saved into the file.\n");
        return true;
    }

    // Setup the data struct for the parallel thread
    data_pthread *data = new data_pthread;
    data->fields = fields;
    data->bounds = bounds();
    data->C = C;
    data->S = S;
    data->data = download(fields);
    if(!data->data.size()){
        return true;
    }
    data->f = create();
    if(!data->f){
        return true;
    }

    // Launch the thread
    pthread_t tid;
    int err;
    err = pthread_create(&tid, NULL, &save_pthread, (void*)data);
    if(err){
        S->addMessageF(3, "Failure launching the parallel thread.\n");
        char err_str[strlen(strerror(err)) + 2];
        strcpy(err_str, strerror(err));
        strcat(err_str, "\n");
        S->addMessage(0, err_str);
        return true;
    }
    _tids.push_back(tid);

    // Clear the already finished threads
    for(i = 0; i < _tids.size(); i++){
        if(pthread_kill(_tids.at(i), 0)){
            _tids.erase(_tids.begin() + i);
        }
    }

    // Check and limit the number of active writing processes
    if(_tids.size() > 2){
        S->addMessageF(2, "More than 2 active writing tasks\n");
        S->addMessageF(0, "This may result in heavy performance penalties, and hard disk failures\n");
        S->addMessageF(0, "Please, consider a reduction of the output printing rate\n");
        while(_tids.size() > 2){
            if(!pthread_kill(_tids.at(0), 0))
            {
                pthread_join(_tids.at(0), NULL);
            }
            _tids.erase(_tids.begin());
        }
    }

    // Update the PVD file
    if(updatePVD()){
        return true;
    }

    return false;
}

vtkXMLUnstructuredGridWriter* VTK::create(){
    char *basename, msg[1024];
    size_t len;
    vtkXMLUnstructuredGridWriter *f = NULL;
    ScreenManager *S = ScreenManager::singleton();
    ProblemSetup *P = ProblemSetup::singleton();

    // Create the file base name
    len = strlen(P->sets.at(setId())->outputPath()) + 8;
    basename = new char[len];
    strcpy(basename, P->sets.at(setId())->outputPath());
    strcat(basename, ".%d.vtu");

    _next_file_index = file(basename, _next_file_index);
    if(!_next_file_index){
        delete[] basename;
        S->addMessageF(3, "Failure getting a valid filename.\n");
        S->addMessage(0, "\tHow do you received this message?.\n");
        return NULL;
    }
    delete[] basename;

    sprintf(msg, "Writing \"%s\" VTK output...\n", file());
    S->addMessageF(1, msg);

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
    S->addMessageF(1, msg);

    bool should_release_doc = false;
    DOMDocument* doc = getPVD(false);
    if(!doc){
        should_release_doc = true;
        doc = getPVD(true);
    }
    DOMElement* root = doc->getDocumentElement();
    if(!root){
        S->addMessageF(3, "Empty XML file found!\n");
        return true;
    }
    n = doc->getElementsByTagName(xmlS("VTKFile"))->getLength();
    if(n != 1){
        sprintf(msg,
                "Expected 1 VTKFile root section, but %u has been found\n",
                n);
        S->addMessageF(3, msg);
        return true;
    }

    DOMNodeList* nodes = root->getElementsByTagName(xmlS("Collection"));
    if(nodes->getLength() != 1){
        sprintf(msg,
                "Expected 1 collection, but %u has been found\n",
                nodes->getLength());
        S->addMessageF(3, msg);
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
        S->addMessageF(3, "XML toolkit writing error.\n");
        sprintf(msg, "\t%s\n", message);
        S->addMessage(0, msg);
        xmlClear();
        return true;
    }
    catch( DOMException& e ){
        char* message = xmlS(e.getMessage());
        S->addMessageF(3, "XML DOM writing error.\n");
        sprintf(msg, "\t%s\n", message);
        S->addMessage(0, msg);
        xmlClear();
        return true;
    }
    catch( ... ){
        S->addMessageF(3, "Writing error.\n");
        S->addMessage(0, "\tUnhandled exception\n");
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
        ProblemSetup *P = ProblemSetup::singleton();
        size_t len = strlen(P->sets.at(setId())->outputPath()) + 5;
        _namePVD = new char[len];
        strcpy(_namePVD, P->sets.at(setId())->outputPath());
        strcat(_namePVD, ".pvd");
    }
    return (const char*)_namePVD;
}

}}  // namespace

#endif // HAVE_VTK
