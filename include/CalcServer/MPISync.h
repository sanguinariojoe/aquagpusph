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
 * @brief Synchronize arrays between processes, sending information by
 * the network.
 * (See Aqua::CalcServer::MPISync for details)
 * @note Hardcoded versions of the files CalcServer/MPISync.cl.in and
 * CalcServer/MPISync.hcl.in are internally included as a text array.
 */

#ifndef MPISYNC_H_INCLUDED
#define MPISYNC_H_INCLUDED

#ifndef HAVE_MPI
#error MPI not available
#endif

#include <mpi.h>

#include <CalcServer.h>
#include <CalcServer/Kernel.h>
#include <CalcServer/RadixSort.h>
#include <CalcServer/Reduction.h>
#include <CalcServer/UnSort.h>
#include <CalcServer/SetScalar.h>
#include <CalcServer/Set.h>

namespace Aqua{ namespace CalcServer{

/** @class MPISync MPISync.h CalcServer/MPISync.h
 * @brief Synchronize arrays between processes.
 * 
 * When MPI is enabled, several instances/processes of AQUAgpusph can be
 * simultaneously launched, letting each process compute a subset of the whole
 * simulation, in such a way the global simulation computation can be
 * accelerated.
 *
 * The data synchronization is a quite expensive operation, both in
 * computational and physical time, so it shall be thoroughly used.
 *
 * First, the tool should dispose the data in a convenient way,
 * which implies a sorting algorithm, as well as reduction operations to compute
 * the amount of data to be sent to each process.
 * After that, the data is sent by the network to the rest of processes, which
 * would take some time, depending on the amount of data to send, and the
 * network speed.
 *
 * In parallel, the tool will prepare everything to download the data incoming
 * from the other processes, INSIDE THE SAME ARRAYS.
 *
 * To reduce the computational cost and avoid aside effects, it is strongly
 * recommended to copy the actual data into helper arrays before synchronizing.
 *
 * @note Since the mask array shall be sorted, power of 2 arrays are required.
 */
class MPISync : public Aqua::CalcServer::Tool
{
public:
    /** Constructor
     * @param name Tool name
     * @param mask Mask of the data to be sent to each process. Numbers out
     * of bounds (i.e. bigger or equal to the number of processes) will be
     * ignored, and therefore not sent anywhere
     * @param fields Fields to be synchronized between processes
     * @param procs Processes to be considered to send information. If an empty
     * list is provided, all the processes will be considered. Providing a list
     * of processes is reducing the number of reductions to be carried out, and
     * therefore the computational cost
     * @param once Run this tool just once. Useful to make initializations
     *
     * @warning The used mask will be overwritten
     */
    MPISync(const std::string name,
            const std::string mask,
            const std::vector<std::string> fields,
            const std::vector<unsigned int> procs,
            bool once=false);

    /** Destructor.
     */
    ~MPISync();

    /** Initialize the tool.
     */
    void setup();

protected:
    /** Execute the tool
     * @param events List of events that shall be waited before safe execution
     * @return OpenCL event to be waited before accessing the dependencies
     */
    cl_event _execute(const std::vector<cl_event> events);

private:
    /** Get the input variables
     */
    void variables();

    /** @brief Create the mask sort tool
     *
     * Creating the sort tool requires creating also a set of additional
     * variables to store the permutations. Those variables will be named alike
     * the mask, with a '__' prefix and a suffix to describe them, i.e.
     * '_sorted' for the unsorted to sorted permutations, and '_unsorted' for
     * the sorted to unsorted permutations.
     */
    void setupSort();

    /** @brief Create a field sort tool
     *
     * The sorted field will be stored in a new field named like the original
     * field, with a '__' prefix and a '_sorted' suffix.
     */
    void setupFieldSort(InputOutput::ArrayVariable* field);

    /** @brief Create the senders to each process
     */
    void setupSenders();

    /** @brief Create the receivers from each process
     */
    void setupReceivers();

    /// Mask name
    std::string _mask_name;
    /// Mask variable
    InputOutput::ArrayVariable *_mask;

    /// List of field names
    std::vector<std::string> _field_names;
    /// List of fields
    std::vector<InputOutput::ArrayVariable*> _fields;
    /// List of sorted fields
    std::vector<InputOutput::ArrayVariable*> _fields_sorted;

    /// List of processes to be considered at the time of sending data
    std::vector<unsigned int> _procs;

    /** Auxiliar variable to store the original index of each sorted component
     * of the mask.
     */
    InputOutput::ArrayVariable *_unsorted_id;

    /** Auxiliar variable to store the sorted index of each unsorted component
     * of the mask.
     */
    InputOutput::ArrayVariable *_sorted_id;

    /// Sorting by cells computation tool
    RadixSort *_sort;

    /// List of field sorters
    std::vector<UnSort*> _field_sorters;

    /// Total number of elements
    unsigned int _n;

public:
    /** @class Exchanger MPISync.h CalcServer/MPISync.h
     * @brief Interprocess array synchronization base class.
     * 
     * This class is used to setup host storages for the send & receive
     * operations, to provide some helper functions to translate
     */
    class Exchanger
    {
    public:
        /** Constructor
         * @param name The same name that the owner tool (See
         * Aqua::CalcServer::MPISync)
         * @param mask Already sorted mask
         * @param fields Already sorted fields
         * @param field_hosts Allocated host memory to temporary copy the data,
         * while it is sent to another process or it is uploaded to the
         * computational device either
         * @param proc Process to which the data shall be sent
         */
        Exchanger(const std::string name,
                  InputOutput::ArrayVariable *mask,
                  const std::vector<InputOutput::ArrayVariable*> fields,
                  const std::vector<void*> field_hosts,
                  const unsigned int proc);

        /** Destructor.
         */
        ~Exchanger();

        /** @brief Parent tool name
         * @return Parent tool name
         */
        const std::string name(){return _name;}

        /** @brief Data structure to store the type information required by MPI
         *
         * MPI requires its own data type enumerate to let the interface know
         * how the data shall be packed/unpacked.
         * To this end, we should store the underlying type, and the number of
         * components
         */
        typedef struct {
            /// Number of components
            unsigned int n;
            /// Underlying type, in MPI format
            MPI_Datatype t;
        } MPIType;

        /** @brief MPI type descriptor
         * @param t Type string
         * @return MPI type descriptor. In case type cannot be handled, a
         * MPI::DATATYPE_NULL type will be returned
         */
        static const MPIType typeToMPI(std::string t);

    protected:
        /// Mask
        InputOutput::ArrayVariable *_mask;

        /// Field
        std::vector<InputOutput::ArrayVariable*> _fields;

        /// Processor
        unsigned int _proc;

        /// Total number of elements
        unsigned int _n;

        /// Host memory arrays to download, send, receive and upload the data
        std::vector<void*> _fields_host;
    private:
        /// Owner tool name
        std::string _name;
    };

    /** @class Sender MPISync.h CalcServer/MPISync.h
     * @brief Synchronize arrays between processes.
     * 
     * 
     */
    class Sender : public Exchanger
    {
    public:
        /** Constructor
         * @param name The same name that the owner tool (See
         * Aqua::CalcServer::MPISync)
         * @param mask Already sorted mask
         * @param fields Already sorted fields
         * @param field_hosts Allocated host memory to temporary copy the data,
         * while it is sent to another process or it is uploaded to the
         * computational device either
         * @param proc Process to which the data shall be sent
         * @param n_offset Variable where the number of already sent particles
         * should be stored. Sending instances will use this variable to both,
         * know when they can proceed, and which is the offset to read from
         * fields
         */
        Sender(const std::string name,
               InputOutput::ArrayVariable *mask,
               const std::vector<InputOutput::ArrayVariable*> fields,
               const std::vector<void*> field_hosts,
               const unsigned int proc,
               InputOutput::UIntVariable *n_offset);

        /** Destructor.
         */
        ~Sender();

        /** @brief Send the information
         */
        void execute(void);
    private:
        /** Create the submask array
        */
        void setupSubMaskMems();

        /** Setup the OpenCL stuff
         * @param kernel_name Name of the kernel, either "n_offset_mask" or
         * "n_send_mask"
         */
        void setupOpenCL(const std::string kernel_name);

        /** Register the "number of elements to send" variable
         * @param var_name Variable name, either "n_offset" or "n_send"
         */
        void setupReduction(const std::string var_name);

        /// Accumulated number of elements sent
        InputOutput::UIntVariable *_n_offset;

        /// Submask memory object
        InputOutput::ArrayVariable* _n_offset_mask;

        /// OpenCL kernel to compute the offset mask
        cl_kernel _n_offset_kernel;

        /// Submaks of elements already processed
        Reduction *_n_offset_reduction;

        /// Number of elements to be sent
        InputOutput::UIntVariable *_n_send;

        /// Submask of elements to be sent
        InputOutput::ArrayVariable* _n_send_mask;

        /// OpenCL kernel to compute the sending mask
        cl_kernel _n_send_kernel;

        /// Reduction to compute the number of elements to send
        Reduction *_n_send_reduction;

        /// Global work sizes in each step
        size_t _global_work_size;
        /// Local work sizes in each step
        size_t _local_work_size;
    };

    /** @class Receiver MPISync.h CalcServer/MPISync.h
     * @brief Synchronize arrays between processes.
     * 
     * 
     */
    class Receiver : public Exchanger
    {
    public:
        /** Constructor
         * @param name The same name that the owner tool (See
         * Aqua::CalcServer::MPISync)
         * @param mask Incoming data process mask
         * @param fields Fields to store the incoming data
         * @param field_hosts Allocated host memory to temporary copy the
         * incoming data, which will be uploaded to the computational device
         * afterwards
         * @param proc Process from which the data shall be received
         * @param n_offset Variable where the number of already received
         * particles should be stored.
         */
        Receiver(const std::string name,
                 InputOutput::ArrayVariable *mask,
                 const std::vector<InputOutput::ArrayVariable*> fields,
                 const std::vector<void*> field_hosts,
                 const unsigned int proc,
                 InputOutput::UIntVariable *n_offset);

        /** Destructor.
         */
        ~Receiver();

        /** @brief Receive the information
         */
        void execute(void);
    private:
        /** Setup the OpenCL stuff
        */
        void setupOpenCL();

        /// OpenCL kernel
        cl_kernel _kernel;

        /// Accumulated number of received elements
        InputOutput::UIntVariable *_n_offset;

        /// Local work sizes in each step
        size_t _local_work_size;
    };


private:
    /// Cumulative number of particles sent
    InputOutput::UIntVariable *_n_offset_send;

    /// Offset reinitialization tool
    SetScalar *_n_offset_send_reinit;

    /// Host memory arrays to download and send data to other processes
    std::vector<void*> _fields_send;

    /// Set of information senders
    std::vector<Sender*> _senders;

    /// Mask reinitialization tool
    Set *_mask_reinit;

    /// Cumulative number of particles sent
    InputOutput::UIntVariable *_n_offset_recv;

    /// Offset reinitialization tool
    SetScalar *_n_offset_recv_reinit;

    /// Host memory arrays to download and send data to other processes
    std::vector<void*> _fields_recv;

    /// Set of information senders
    std::vector<Receiver*> _receivers;
};

}}  // namespace

#endif // MPISYNC_H_INCLUDED
