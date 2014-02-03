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

#ifndef REQUESTED_FIELDS
    #ifdef HAVE_3D
        #define REQUESTED_FIELDS 17
    #else
        #define REQUESTED_FIELDS 13
    #endif
#endif // REQUESTED_FIELDS

namespace Aqua{ namespace InputOutput{

VTK::VTK(unsigned int first, unsigned int n, unsigned int ifluid)
    : Particles(first, n, ifluid)
{
}

VTK::~VTK()
{
}

bool VTK::load()
{
    unsigned int i, n, N, cell=0;
    int aux;
    char msg[1024];
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();
	Fluid *F = Fluid::singleton();

	sprintf(msg,
            "Loading fluid from VTK file \"%s\"\n",
            P->fluids[fluidId()].in_path);
	S->addMessageF(1, msg);

    vtkSmartPointer<vtkXMLUnstructuredGridReader> f =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    if(!f->CanReadFile(P->fluids[fluidId()].in_path)){
        S->addMessageF(3, "Teh file cannot be readed.\n");
        return true;
    }

    f->SetFileName(P->fluids[fluidId()].in_path);
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

    for(i=bounds().x; i<bounds().y; i++){
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

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> f = create();
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
    vtkSmartPointer<vtkFloatArray> h = vtkSmartPointer<vtkFloatArray>::New();
    h->SetNumberOfComponents(1);
    h->SetName("h");
    vtkSmartPointer<vtkFloatArray> m = vtkSmartPointer<vtkFloatArray>::New();
    m->SetNumberOfComponents(1);
    m->SetName("mass");
    vtkSmartPointer<vtkFloatArray> W = vtkSmartPointer<vtkFloatArray>::New();
    W->SetNumberOfComponents(1);
    W->SetName("sumW");
    vtkSmartPointer<vtkFloatArray> dW = vtkSmartPointer<vtkFloatArray>::New();
    dW->SetNumberOfComponents(3);
    dW->SetName("gradW");
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
            vect[0] = F->shepard_gradient[i].x;
            vect[1] = F->shepard_gradient[i].y;
            vect[2] = F->shepard_gradient[i].z;
            dW->InsertNextTupleValue(vect);
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
            vect[0] = F->shepard_gradient[i].x;
            vect[1] = F->shepard_gradient[i].y;
            dW->InsertNextTupleValue(vect);
        #endif // HAVE_3D
        p->InsertNextValue(F->press[i]);
        dens->InsertNextValue(F->dens[i]);
        drdt->InsertNextValue(F->drdt[i]);
        h->InsertNextValue(F->hp[i]);
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
    grid->GetPointData()->AddArray(h);
    grid->GetPointData()->AddArray(m);
    grid->GetPointData()->AddArray(W);
    grid->GetPointData()->AddArray(dW);
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

    return false;
}

vtkSmartPointer<vtkXMLUnstructuredGridWriter> VTK::create(){
    char *basename, msg[1024];
    size_t len;
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> f = NULL;
	ScreenManager *S = ScreenManager::singleton();
	ProblemSetup *P = ProblemSetup::singleton();

    // Create the file base name
    len = strlen(P->fluids[fluidId()].out_path) + 8;
    basename = new char[len];
    strcpy(basename, P->fluids[fluidId()].out_path);
    strcat(basename, ".%d.vtu");

    if(file(basename, 0)){
        delete[] basename;
        S->addMessageF(3, "Failure getting a valid filename.\n");
        S->addMessageF(0, "\tHow do you received this message?.\n");
        return NULL;
    }
    delete[] basename;

    f = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    f->SetFileName(file());

	return f;
}

}}  // namespace

#endif // HAVE_VTK
