/*!
 * \file floe/io/mpi_utilities.hpp
 * \brief mpi send and receive serialized objects
 * \author Quentin Jouet
 */

#ifndef IO_MPI_UTILS_HPP
#define IO_MPI_UTILS_HPP

#include <mpi.h>
#include <iostream> // debug
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/list.hpp>
#include "floe/io/inter_process_message.hpp"

namespace floe { namespace io
{

class MPIWrapper{
public:
    MPIWrapper(int buffer_size=1e6) : m_buffer_size{buffer_size}, m_msg_buffer{new char[buffer_size]} { 
        MPI_Comm_rank( MPI_COMM_WORLD, &m_mpi_rank );
    }

    template<typename TObject>
    void send_serial(TObject const& msg, int process_id, int tag){
        std::stringstream ss; // any stream can be used
        cereal::BinaryOutputArchive oarchive(ss); // Create an this->output archive
        oarchive(msg); // Write the data to the archive
        MPI_Send((void*)ss.str().c_str(), ss.str().length(), MPI_BYTE, process_id, tag, MPI_COMM_WORLD);
    }

    template<typename TObject>
    TObject receive_serial(int source, int tag){
        MPI_Status status;
        MPI_Recv(msg_buffer(), buffer_size(), MPI_BYTE, source, tag, MPI_COMM_WORLD, &status);
        int len;
        MPI_Get_count(&status, MPI_BYTE, &len);
        // int source = status.MPI_SOURCE;
        // auto tag = status.MPI_TAG;
        std::stringstream ss; 
        ss.write(msg_buffer(), len);
        cereal::BinaryInputArchive iarchive(ss); // Create an input archive
        
        // if (tag==0){ // todo when different messages class are needed
            TObject response;
            iarchive(response);
        // }
        response.mpi_source(status.MPI_SOURCE);
        return response;
    }

    inline int process_rank(){ return m_mpi_rank; }
private:
    int m_mpi_rank;
    int m_buffer_size;
    char* m_msg_buffer;
    inline char* msg_buffer() { return m_msg_buffer; }
    inline int buffer_size() { return m_buffer_size; }
};

}} // namespace floe::io


#endif // IO_MPI_UTILS_HPP