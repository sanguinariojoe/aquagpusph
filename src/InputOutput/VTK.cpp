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

#include <InputOutput/VTK.h>
#include <ScreenManager.h>
#include <ProblemSetup.h>
#include <TimeManager.h>
#include <Fluid.h>
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
{
}

VTK::~VTK()
{
}

bool VTK::load()
{
    unsigned int i, n, N, cell=0, progress;
    int aux;
    char msg[1024];
    ScreenManager *S = ScreenManager::singleton();
    ProblemSetup *P = ProblemSetup::singleton();
    Fluid *F = Fluid::singleton();

    sprintf(msg,
            "Loading fluid from VTK file \"%s\"\n",
            P->sets.at(setId())->inputPath());
    S->addMessageF(1, msg);

    vtkSmartPointer<vtkXMLUnstructuredGridReader> f =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    if(!f->CanReadFile(P->sets.at(setId())->inputPath())){
        S->addMessageF(3, "Teh file cannot be readed.\n");
        return true;
    }

    f->SetFileName(P->sets.at(setId())->inputPath());
    f->Update();

    vtkSmartPointer<vtkUnstructuredGrid> grid = f->GetOutput();

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

    vtkSmartPointer<vtkPoints> points = grid->GetPoints();
    vtkSmartPointer<vtkPointData> data = grid->GetPointData();
    vtkSmartPointer<vtkFloatArray> normal = (vtkFloatArray*)(data->GetArray("n", aux));
    vtkSmartPointer<vtkFloatArray> v = (vtkFloatArray*)(data->GetArray("v", aux));
    vtkSmartPointer<vtkFloatArray> dv = (vtkFloatArray*)(data->GetArray("dv/dt", aux));
    vtkSmartPointer<vtkFloatArray> dens = (vtkFloatArray*)(data->GetArray("dens", aux));
    vtkSmartPointer<vtkFloatArray> drdt = (vtkFloatArray*)(data->GetArray("ddens/dt", aux));
    vtkSmartPointer<vtkFloatArray> m = (vtkFloatArray*)(data->GetArray("mass", aux));
    vtkSmartPointer<vtkIntArray> imove = (vtkIntArray*)(data->GetArray("imove", aux));

    progress = -1;
    for(i=bounds().x; i<bounds().y; i++){
        F->ifluid[i] = setId();

        double *vect;
        vect = points->GetPoint(cell);
        F->pos[i].x = vect[0];
        F->pos[i].y = vect[1];
        #ifdef HAVE_3D
            F->pos[i].z = vect[2];
        #endif
        vect = normal->GetTuple(cell);
        F->normal[i].x = vect[0];
        F->normal[i].y = vect[1];
        #ifdef HAVE_3D
            F->normal[i].z = vect[2];
        #endif
        vect = v->GetTuple(cell);
        F->v[i].x = vect[0];
        F->v[i].y = vect[1];
        #ifdef HAVE_3D
            F->v[i].z = vect[2];
        #endif
        vect = dv->GetTuple(cell);
        F->f[i].x = vect[0];
        F->f[i].y = vect[1];
        #ifdef HAVE_3D
            F->f[i].z = vect[2];
        #endif
        F->dens[i] = dens->GetComponent(cell, 0);
        F->drdt[i] = drdt->GetComponent(cell, 0);
        F->mass[i] = m->GetComponent(cell, 0);
        F->imove[i] = imove->GetComponent(cell, 0);

        if(progress != cell * 100 / n){
            progress = cell * 100 / n;
            if(!(progress % 10)){
                sprintf(msg, "\t\t%u%%\n", progress);
                S->addMessage(0, msg);
            }
        }
        cell++;
    }

    return false;
}

bool VTK::save()
{
    unsigned int i;
    vtkSmartPointer<vtkVertex> vertex;
    float vect[3] = {0.f, 0.f, 0.f};
    ScreenManager *S = ScreenManager::singleton();
    TimeManager *T = TimeManager::singleton();
    Fluid *F = Fluid::singleton();

    vtkXMLUnstructuredGridWriter *f = create();
    if(!f)
        return true;

    // Setup the arrays to write
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> normal = vtkSmartPointer<vtkFloatArray>::New();
    normal->SetNumberOfComponents(3);
    normal->SetName("n");
    vtkSmartPointer<vtkFloatArray> v = vtkSmartPointer<vtkFloatArray>::New();
    v->SetNumberOfComponents(3);
    v->SetName("v");
    vtkSmartPointer<vtkFloatArray> dv = vtkSmartPointer<vtkFloatArray>::New();
    dv->SetNumberOfComponents(3);
    dv->SetName("dv/dt");
    vtkSmartPointer<vtkFloatArray> p = vtkSmartPointer<vtkFloatArray>::New();
    p->SetNumberOfComponents(1);
    p->SetName("press");
    vtkSmartPointer<vtkFloatArray> dens = vtkSmartPointer<vtkFloatArray>::New();
    dens->SetNumberOfComponents(1);
    dens->SetName("dens");
    vtkSmartPointer<vtkFloatArray> drdt = vtkSmartPointer<vtkFloatArray>::New();
    drdt->SetNumberOfComponents(1);
    drdt->SetName("ddens/dt");
    vtkSmartPointer<vtkFloatArray> m = vtkSmartPointer<vtkFloatArray>::New();
    m->SetNumberOfComponents(1);
    m->SetName("mass");
    vtkSmartPointer<vtkFloatArray> W = vtkSmartPointer<vtkFloatArray>::New();
    W->SetNumberOfComponents(1);
    W->SetName("sumW");
    vtkSmartPointer<vtkIntArray> imove = vtkSmartPointer<vtkIntArray>::New();
    imove->SetNumberOfComponents(1);
    imove->SetName("imove");
    vtkSmartPointer<vtkUnsignedIntArray> id = vtkSmartPointer<vtkUnsignedIntArray>::New();
    id->SetNumberOfComponents(1);
    id->SetName("id");
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    // Fill them with the data
    for(i=bounds().x; i<bounds().y; i++){
        #ifdef HAVE_3D
            points->InsertNextPoint(F->pos[i].x, F->pos[i].y, F->pos[i].z);
            vect[0] = F->normal[i].x;
            vect[1] = F->normal[i].y;
            vect[2] = F->normal[i].z;
            normal->InsertNextTupleValue(vect);
            vect[0] = F->v[i].x;
            vect[1] = F->v[i].y;
            vect[2] = F->v[i].z;
            v->InsertNextTupleValue(vect);
            vect[0] = F->f[i].x;
            vect[1] = F->f[i].y;
            vect[2] = F->f[i].z;
            dv->InsertNextTupleValue(vect);
        #else // HAVE_3D
            points->InsertNextPoint(F->pos[i].x, F->pos[i].y, 0.f);
            vect[0] = F->normal[i].x;
            vect[1] = F->normal[i].y;
            normal->InsertNextTupleValue(vect);
            vect[0] = F->v[i].x;
            vect[1] = F->v[i].y;
            v->InsertNextTupleValue(vect);
            vect[0] = F->f[i].x;
            vect[1] = F->f[i].y;
            dv->InsertNextTupleValue(vect);
        #endif // HAVE_3D
        p->InsertNextValue(F->press[i]);
        dens->InsertNextValue(F->dens[i]);
        drdt->InsertNextValue(F->drdt[i]);
        m->InsertNextValue(F->mass[i]);
        W->InsertNextValue(F->shepard[i]);
        imove->InsertNextValue(F->imove[i]);
        id->InsertNextValue(i);

        vertex = vtkSmartPointer<vtkVertex>::New();
        vertex->GetPointIds()->SetId(0, i - bounds().x);
        cells->InsertNextCell(vertex);
    }

    // Setup the unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    grid->SetCells(vertex->GetCellType(), cells);
    grid->GetPointData()->AddArray(normal);
    grid->GetPointData()->AddArray(v);
    grid->GetPointData()->AddArray(dv);
    grid->GetPointData()->AddArray(p);
    grid->GetPointData()->AddArray(dens);
    grid->GetPointData()->AddArray(drdt);
    grid->GetPointData()->AddArray(m);
    grid->GetPointData()->AddArray(W);
    grid->GetPointData()->AddArray(imove);
    grid->GetPointData()->AddArray(id);

    // Write file
    #if VTK_MAJOR_VERSION <= 5
        f->SetInput(grid);
    #else // VTK_MAJOR_VERSION
        f->SetInputData(grid);
    #endif // VTK_MAJOR_VERSION

    if(!f->Write()){
        S->addMessageF(3, "Failure writing the VTK file.\n");
        return true;
    }
    f->Delete();

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

    if(file(basename, 0)){
        delete[] basename;
        S->addMessageF(3, "Failure getting a valid filename.\n");
        S->addMessageF(0, "\tHow do you received this message?.\n");
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

    DOMDocument* doc = getPVD();
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
    // doc->release();
    xmlClear();

    return false;
}

DOMDocument* VTK::getPVD()
{
    DOMDocument* doc = NULL;
    FILE *dummy=NULL;

    // Try to open as ascii file, just to know if the file already exist
    dummy = fopen(filenamePVD(), "r");
    if(!dummy){
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

static char* namePVD = NULL;

const char* VTK::filenamePVD()
{
    if(!namePVD){
        ProblemSetup *P = ProblemSetup::singleton();
        size_t len = strlen(P->sets.at(setId())->outputPath()) + 5;
        namePVD = new char[len];
        strcpy(namePVD, P->sets.at(setId())->outputPath());
        strcat(namePVD, ".pvd");
    }
    return (const char*)namePVD;
}

}}  // namespace

#endif // HAVE_VTK
