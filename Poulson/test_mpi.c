#include "mpi.h"

int
main( int argc, char* argv[] )
{
    int i;
    int commSize;
    int *recvCounts;
    double *sendBuf, *recvBuf;
    MPI_Init( &argc, &argv );

    MPI_Comm_size( MPI_COMM_WORLD, &commSize );

    recvCounts = malloc( commSize*sizeof(int) ); 
    for( i=0; i<commSize; ++i )
        recvCounts[i] = 1;
    sendBuf = malloc( commSize*sizeof(double) );
    recvBuf = malloc( sizeof(double) );

    MPI_Reduce_scatter
    ( sendBuf, recvBuf, recvCounts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    MPI_Finalize();
    return 0;
}
