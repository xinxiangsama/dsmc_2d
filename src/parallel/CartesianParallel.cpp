#include "CartesianParallel.h"

void CartesianParallel::ZoneDcomposition()
{
    int dims[2];
    int periods[2] = {0, 0}; // non-periodic
    MPI_Dims_create(numprocs, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &m_cartesian_comm);

    int coords[2];
    MPI_Cart_coords(m_cartesian_comm, myid, 2, coords);
}