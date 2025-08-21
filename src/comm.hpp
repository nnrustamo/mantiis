#include <mpi.h>

namespace mantiis_parallel
{
    void launchMPI(int &proc_id, int &num_of_proc, int &num_thread)
    {
        omp_set_num_threads(num_thread);            
        MPI_Init(NULL, NULL);                       
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);    
        MPI_Comm_size(MPI_COMM_WORLD, &num_of_proc);
    }
}
