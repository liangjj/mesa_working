#include <vector>
#include "mpi.h"

int
main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );

    int commSize, commRank;
    MPI_Comm_size( MPI_COMM_WORLD, &commSize );
    MPI_Comm_rank( MPI_COMM_WORLD, &commRank );

    std::vector<int> recvCounts( commSize, 1 );
    std::vector<double> sendBuf( commSize );
    std::vector<double> recvBuf( 1 );

    MPI_Reduce_scatter
    ( &sendBuf[0], &recvBuf[0], &recvCounts[0], 
      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    MPI_Finalize();
    return 0;
}
