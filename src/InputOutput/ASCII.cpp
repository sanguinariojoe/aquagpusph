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
 * @brief Particles plain text data files loader/saver (with math expressions
 * evaluator).
 * (See Aqua::InputOutput::ASCII for details)
 */

#include <InputOutput/ASCII.h>
#include <InputOutput/Logger.h>
#include <ProblemSetup.h>
#include <CalcServer.h>
#include <AuxiliarMethods.h>

#ifndef MAX_LINE_LEN
    #define MAX_LINE_LEN 1024
#endif // MAX_LINE_LEN

namespace Aqua{ namespace InputOutput{

ASCII::ASCII(ProblemSetup& sim_data,
             unsigned int iset,
             unsigned int first,
             unsigned int n_in)
    : Particles(sim_data, iset, first, n_in)
{
    if(n() == 0) {
        n(compute_n());
    }
}

ASCII::~ASCII()
{
}

void ASCII::load()
{
    std::ifstream f;
    cl_int err_code;
    char *pos = NULL;
    unsigned int n, N, n_fields;
    CalcServer::CalcServer *C = CalcServer::CalcServer::singleton();

    loadDefault();

    std::ostringstream msg;
    msg << "Loading particles from ASCII file \""
             <<  simData().sets.at(setId())->inputPath()
             << "\"..." << std::endl;
    LOG(L_INFO, msg.str());

    f.open(simData().sets.at(setId())->inputPath());
    if(!f) {
        std::ostringstream msg;
        msg << "Failure reading the file \"" <<
               simData().sets.at(setId())->inputPath() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::ifstream::failure(msg.str());
    }

    // Assert that the number of particles is right
    n = bounds().y - bounds().x;
    N = readNParticles(f);
    if(n != N){
        std::ostringstream msg;
        msg << "Expected " << n << " particles, but the file contains just "
            << N << " ones." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::runtime_error("Invalid number of particles in file");
    }

    // Check the fields to read
    std::vector<std::string> fields = simData().sets.at(setId())->inputFields();
    if(!fields.size()){
        LOG(L_ERROR, "0 fields were set to be read from the file.\n");
        throw std::runtime_error("No fields have to be read");
    }
    bool have_r = false;
    for(auto field : fields){
        if(!field.compare("r")){
            have_r = true;
            break;
        }
    }
    if(!have_r){
        LOG(L_ERROR, "\"r\" field was not set to be read from the file.\n");
        throw std::runtime_error("Reading \"r\" field is mandatory");
    }
    // Setup an storage
    std::vector<void*> data;
    Variables *vars = C->variables();
    n_fields = 0;
    for(auto field : fields){
        if(!vars->get(field)){
            std::ostringstream msg;
            msg << "Undeclared variable \"" << field
                << "\" set to be read." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable");
        }
        if(vars->get(field)->type().find('*') == std::string::npos){
            std::ostringstream msg;
            msg << "Can't read scalar variable \"" << field
                << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable type");
        }
        ArrayVariable *var = (ArrayVariable*)vars->get(field);
        n_fields += vars->typeToN(var->type());
        size_t typesize = vars->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < bounds().y) {
            std::ostringstream msg;
            msg << "Array variable \"" << field
                << "\" is not long enough." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable length");
        }
        void *store = malloc(typesize * n);
        if(!store){
            std::ostringstream msg;
            msg << "Failure allocating " << typesize * n
                << "bytes for variable \"" << field
                << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::bad_alloc();
        }
        data.push_back(store);
    }

    // Read the particles
    f.clear();
    f.seekg(0);
    unsigned int i=0, i_line=0, progress=-1;
    std::string line;
    while(getline(f, line))
    {
        i_line++;

        formatLine(line);
        if(line == "")
            continue;

        unsigned int n_available_fields = readNFields(line);
        if(n_available_fields != n_fields){
            std::ostringstream msg;
            msg << "Line " << i_line << " has " << n_available_fields
                << " fields, but " << n_fields << " are required." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Bad formatted file");
        }

        unsigned int j = 0;
        for(auto field : fields){
            line = readField(field, line, i, data.at(j));
            j++;
        }

        i++;

        if(progress != i * 100 / n){
            progress = i * 100 / n;
            if(!(progress % 10)){
                msg.str("");
                msg << "\t\t" << progress << "%" << std::endl;
                LOG0(L_DEBUG, msg.str());
            }
        }
    }

    // Send the data to the server and release it
    i = 0;
    for(auto field : fields){
        ArrayVariable *var = (ArrayVariable*)vars->get(field);
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
            std::ostringstream msg;
            msg << "Failure sending variable \"" << field
                << "\" to the server." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("OpenCL error");
        }
        i++;
    }
    data.clear();

    f.close();
}

void ASCII::save(float t)
{
    unsigned int i, j;
    cl_int err_code;
    Variables *vars = CalcServer::CalcServer::singleton()->variables();

    std::vector<std::string> fields = simData().sets.at(setId())->outputFields();
    if(!fields.size()){
        LOG(L_ERROR, "0 fields were set to be saved into the file.\n");
        throw std::runtime_error("No fields have been set to be saved");
    }

    std::ofstream f;
    create(f);

    // Write a head
    f << "#########################################################" << std::endl;
    f << "#                                                       #" << std::endl;
    f << "#    #    ##   #  #   #                           #     #" << std::endl;
    f << "#   # #  #  #  #  #  # #                          #     #" << std::endl;
    f << "#  ##### #  #  #  # #####  ##  ###  #  #  ## ###  ###   #" << std::endl;
    f << "#  #   # #  #  #  # #   # #  # #  # #  # #   #  # #  #  #" << std::endl;
    f << "#  #   # #  #  #  # #   # #  # #  # #  #   # #  # #  #  #" << std::endl;
    f << "#  #   #  ## #  ##  #   #  ### ###   ### ##  ###  #  #  #" << std::endl;
    f << "#                            # #             #          #" << std::endl;
    f << "#                          ##  #             #          #" << std::endl;
    f << "#                                                       #" << std::endl;
    f << "#########################################################" << std::endl;
    f << "#" << std::endl;
    f << "#    File autogenerated by AQUAgpusph" << std::endl;
    f << "#    t = " << t << " s" << std::endl;
    f << "#" << std::endl;
    f << "#########################################################" << std::endl;
    f << std::endl;
    f.flush();

    for(auto field : fields){
        if(!vars->get(field)){
            std::ostringstream msg;
            msg << "Can't save undeclared variable \"" << field
                << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable");
        }
        if(vars->get(field)->type().find('*') == std::string::npos){
            std::ostringstream msg;
            msg << "Can't save scalar variable \"" << field
                << "\"." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable type");
        }
        ArrayVariable *var = (ArrayVariable*)vars->get(field);
        size_t typesize = vars->typeToBytes(var->type());
        size_t len = var->size() / typesize;
        if(len < bounds().y){
            std::ostringstream msg;
            msg << "Variable \"" << field
                << "\" is not long enough." << std::endl;
            LOG(L_ERROR, msg.str());
            throw std::runtime_error("Invalid variable length");
        }
    }
    std::vector<void*> data = download(fields);

    for(i = 0; i < bounds().y - bounds().x; i++){
        for(j = 0; j < fields.size(); j++){
            ArrayVariable *var = (ArrayVariable*)vars->get(fields.at(j).c_str());
            std::string type_name = var->type();
            if(!type_name.compare("int*")){
                int* v = (int*)data.at(j);
                f << v[i] << ",";
            }
            else if(!type_name.compare("unsigned int*")){
                unsigned int* v = (unsigned int*)data.at(j);
                f << v[i] << ",";
            }
            else if(!type_name.compare("float*")){
                float* v = (float*)data.at(j);
                f << v[i] << ",";
            }
            else if(!type_name.compare("ivec*")){
                #ifdef HAVE_3D
                    ivec* v = (ivec*)data.at(j);
                    f << v[i].x << " "
                      << v[i].y << " "
                      << v[i].z << " "
                      << v[i].w << ",";
                #else
                    ivec* v = (ivec*)data.at(j);
                    f << v[i].x << " "
                      << v[i].y << ",";
                #endif // HAVE_3D
            }
            else if(!type_name.compare("ivec2*")){
                ivec2* v = (ivec2*)data.at(j);
                f << v[i].x << " "
                  << v[i].y << ",";
            }
            else if(!type_name.compare("ivec3*")){
                ivec3* v = (ivec3*)data.at(j);
                f << v[i].x << " "
                  << v[i].y << " "
                  << v[i].z << ",";
            }
            else if(!type_name.compare("ivec4*")){
                ivec4* v = (ivec4*)data.at(j);
                f << v[i].x << " "
                  << v[i].y << " "
                  << v[i].z << " "
                  << v[i].w << ",";
            }
            else if(!type_name.compare("uivec*")){
                #ifdef HAVE_3D
                    uivec* v = (uivec*)data.at(j);
                    f << v[i].x << " "
                      << v[i].y << " "
                      << v[i].z << " "
                      << v[i].w << ",";
                #else
                    uivec* v = (uivec*)data.at(j);
                    f << v[i].x << " "
                      << v[i].y << ",";
                #endif // HAVE_3D
            }
            else if(!type_name.compare("uivec2*")){
                uivec2* v = (uivec2*)data.at(j);
                f << v[i].x << " "
                  << v[i].y << ",";
            }
            else if(!type_name.compare("uivec3*")){
                uivec3* v = (uivec3*)data.at(j);
                f << v[i].x << " "
                  << v[i].y << " "
                  << v[i].z << ",";
            }
            else if(!type_name.compare("uivec4*")){
                uivec4* v = (uivec4*)data.at(j);
                f << v[i].x << " "
                  << v[i].y << " "
                  << v[i].z << " "
                  << v[i].w << ",";
            }
            else if(!type_name.compare("vec*")){
                #ifdef HAVE_3D
                    vec* v = (vec*)data.at(j);
                    f << v[i].x << " "
                      << v[i].y << " "
                      << v[i].z << " "
                      << v[i].w << ",";
                #else
                    vec* v = (vec*)data.at(j);
                    f << v[i].x << " "
                      << v[i].y << ",";
                #endif // HAVE_3D
            }
            else if(!type_name.compare("vec2*")){
                vec2* v = (vec2*)data.at(j);
                f << v[i].x << " "
                  << v[i].y << ",";
            }
            else if(!type_name.compare("vec3*")){
                vec3* v = (vec3*)data.at(j);
                f << v[i].x << " "
                  << v[i].y << " "
                  << v[i].z << ",";
            }
            else if(!type_name.compare("vec4*")){
                vec4* v = (vec4*)data.at(j);
                f << v[i].x << " "
                  << v[i].y << " "
                  << v[i].z << " "
                  << v[i].w << ",";
            }
            else if(!type_name.compare("matrix*")){
                #ifdef HAVE_3D
                    matrix* m = (matrix*)data.at(j);
                    f << m[i].s0 << " " << m[i].s1 << " " << m[i].s2 << " " << m[i].s3 << " "
                      << m[i].s4 << " " << m[i].s5 << " " << m[i].s6 << " " << m[i].s7 << " "
                      << m[i].s8 << " " << m[i].s9 << " " << m[i].sA << " " << m[i].sB << " "
                      << m[i].sC << " " << m[i].sD << " " << m[i].sE << " " << m[i].sF << ",";
                #else
                    matrix* m = (matrix*)data.at(j);
                    f << m[i].s0 << " " << m[i].s1 << " "
                      << m[i].s2 << " " << m[i].s3 << ",";
                #endif // HAVE_3D
            }
        }
        f << std::endl;
        f.flush();
    }

    for(auto d : data){
        free(d);
    }
    data.clear();

    f.close();
}

const unsigned int ASCII::compute_n()
{
    std::ifstream f;
    f.open(simData().sets.at(setId())->inputPath());
    if(!f) {
        std::ostringstream msg;
        msg << "Failure reading the file \"" <<
               simData().sets.at(setId())->inputPath() << "\"." << std::endl;
        LOG(L_ERROR, msg.str());
        throw std::ifstream::failure(msg.str());
    }

    unsigned int n = readNParticles(f);
    f.close();
    return (const unsigned int)n;
}

unsigned int ASCII::readNParticles(std::ifstream& f)
{
    if(!f.is_open())
        return 0;
    f.clear();
    f.seekg(0);

    std::string line;
    unsigned int n=0;
    while(getline(f, line))
    {
        formatLine(line);
        if(line == "")
            continue;

        n++;
    }

    return n;
}

void ASCII::formatLine(std::string& l)
{
    if(l == "")
        return;

    unsigned int i;

    // Look for comments and discard them
    if(l.find('#') != std::string::npos)
        l.erase(l.begin() + l.find('#'), l.end());

    trim(l);
    // Replace all the separators by commas
    const char *separators = " ;()[]{}\t";
    for(i=0; i<strlen(separators); i++){
        replaceAll(l, std::string(1, separators[i]), ",");
    }

    // Remove all the concatenated separators
    while(l.find(",,") != std::string::npos){
        replaceAll(l, ",,", ",");
    }

    // Remove the preceding and trailing separators
    while ((l.size() > 0) && (l.front() == ',')) {
        l.erase(l.begin());
    }
    while ((l.size() > 0) && (l.back() == ',')) {
        l.pop_back();
    }
}

unsigned int ASCII::readNFields(std::string l)
{
    if (l == "") {
        return 0;
    }

    unsigned int n = 0;
    std::istringstream fields(l);
    std::string s;
    // The user may ask to break lines by semicolon usage. In such a case,
    // we are splitting the input and recursively calling this function with
    // each piece
    while (getline(fields, s, ',')) {
        n++;
    }

    return n;
}

static std::string _remaining;

std::string ASCII::readField(const std::string field,
                             const std::string line,
                             unsigned int index,
                             void* data)
{
    unsigned int i;
    Variables *vars = CalcServer::CalcServer::singleton()->variables();
    ArrayVariable *var = (ArrayVariable*)vars->get(field);

    unsigned int n = vars->typeToN(var->type());
    size_t type_size = vars->typeToBytes(var->type());

    void* ptr = (void*)((char*)data + type_size * index);
    try {
        vars->solve(var->type(), line, ptr);
    } catch (...) {
        return NULL;
    }

    std::string _remaining = line;
    for(i = 0; i < n; i++){
        std::size_t sep = _remaining.find(',');
        if (sep == std::string::npos) {
            _remaining = "";
            break;
        }
        _remaining = _remaining.substr(sep + 1);
    }
    return _remaining;
}

void ASCII::create(std::ofstream& f){
    std::string basename = simData().sets.at(setId())->outputPath();
    // Check that {index} scape string is present, for backward compatibility
    if(basename.find("{index}") == std::string::npos){
        basename += ".{index}.dat";
    }
    _next_file_index = file(basename, _next_file_index);

    std::ostringstream msg;
    msg << "Writing \"" << file() << "\" ASCII file..." << std::endl;
    LOG(L_INFO, msg.str());
    _next_file_index++;

    f.open(file());
}

}}  // namespace
