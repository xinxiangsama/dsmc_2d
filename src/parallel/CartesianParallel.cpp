#include "CartesianParallel.h"
#include <iostream>
#include <fmtmsg.h>
void CartesianParallel::ZoneDcomposition()
{
    int dims[2] {};
    int periods[2] = {0, 0}; // non-periodic
    MPI_Dims_create(numprocs, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &m_cartesian_comm);

    m_mesh->setnumberCpuX(dims[0]);
    m_mesh->setnumberCpuY(dims[1]);

    int coords[2];
    MPI_Cart_coords(m_cartesian_comm, myid, 2, coords);

    m_mesh->setCpuCoordX(coords[0]);
    m_mesh->setCpuCoordY(coords[1]);

    auto N1 = m_mesh->getnumberCellsXGlobal();
    auto N2 = m_mesh->getnumberCellsYGlobal();
    m_mesh->setnumberCellsX(N1 / dims[0]);
    m_mesh->setnumberCellsY(N2 / dims[1]);
    m_mesh->setoffsetX(coords[0] * m_mesh->getnumberCellsX());
    m_mesh->setoffsetY(coords[1] * m_mesh->getnumberCellsY());

    auto L1 = m_mesh->getGlobalLengthX();
    auto L2 = m_mesh->getGlobalLengthY();
    m_mesh->setLocalLengthX(L1 / dims[0]);
    m_mesh->setLocalLengthY(L2 / dims[1]);

}

void CartesianParallel::setNeibours()
{
    int left, right, bottom, top;
    MPI_Cart_shift(m_cartesian_comm, 0, 1, &left, &right);
    MPI_Cart_shift(m_cartesian_comm, 1, 1, &bottom, &top);
    m_neighbours[0] = left;
    m_neighbours[1] = right;
    m_neighbours[2] = bottom;
    m_neighbours[3] = top;

}

const int &CartesianParallel::getLeftNeibour()
{
    return m_neighbours[0];
}
const int &CartesianParallel::getRightNeibour()
{
    return m_neighbours[1];
}
const int &CartesianParallel::getTopNeibour()
{
    return m_neighbours[3];
}
const int &CartesianParallel::getBottomNeibour()
{
    return m_neighbours[2];
}

void CartesianParallel::info()
{   
    if (myid == 0)
    {
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "|               Parallel Information        |" << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "| Number of processors            | " << numprocs  << std::endl;
        std::cout << "| Number of processors in X       | " << m_mesh->getnumberCpuX() << std::endl;
        std::cout << "| Number of processors in Y       | " << m_mesh->getnumberCpuY() << std::endl;
        std::cout << "| Number of cells in X            | " << m_mesh->getnumberCellsXGlobal() << std::endl;
        std::cout << "| Number of cells in Y            | " << m_mesh->getnumberCellsYGlobal()<< std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    std::cout << "Processor " << myid << " Neighbours:" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "| Property                        | Value" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "| Left neighbour                  | " << m_neighbours[0] << std::endl;
    std::cout << "| Right neighbour                 | " << m_neighbours[1] <<std::endl;
    std::cout << "| Bottom neighbour                | " << m_neighbours[2] << std::endl;
    std::cout << "| Top neighbour                   | " << m_neighbours[3] << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    std::cout << "Processor " << myid << " Information:" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "| Property                        | Value  " << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "| Coordinates                     | (" << m_mesh->getCpuCoordX() << ", " << m_mesh->getCpuCoordY() << ")" << std::endl;
    std::cout << "| Local number of cells in X      | " << m_mesh->getnumberCellsX() << std::endl;
    std::cout << "| Local number of cells in Y      | " << m_mesh->getnumberCellsY()  << std::endl;
    std::cout << "| Local domain length in X        | " << m_mesh->getLocalLengthX()  << std::endl;
    std::cout << "| Local domain length in Y        | " << m_mesh->getLocalLengthY()  << std::endl;
    std::cout << "| Offset in X                     | " << m_mesh->getoffsetX()  << std::endl;
    std::cout << "| Offset in Y                     | " << m_mesh->getoffsetY() << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
}
