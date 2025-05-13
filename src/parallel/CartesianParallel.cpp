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

    // only decompose in 2D !!!
    auto N1 = m_mesh->getnumberCellsXGlobal();
    auto N2 = m_mesh->getnumberCellsYGlobal();
    auto N3 = m_mesh->getnumberCellsZGlobal();
    m_mesh->setnumberCellsX(N1 / dims[0]);
    m_mesh->setnumberCellsY(N2 / dims[1]);
    m_mesh->setnumberCellsZ(N3);
    m_mesh->setoffsetX(coords[0] * m_mesh->getnumberCellsX());
    m_mesh->setoffsetY(coords[1] * m_mesh->getnumberCellsY());

    auto L1 = m_mesh->getGlobalLengthX();
    auto L2 = m_mesh->getGlobalLengthY();
    auto L3 = m_mesh->getGlobalLengthZ();
    m_mesh->setLocalLengthX(L1 / dims[0]);
    m_mesh->setLocalLengthY(L2 / dims[1]);
    m_mesh->setLocalLengthZ(L3);

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

void CartesianParallel::exchangedata()
{
    MPI_Datatype MPI_Particle;
    MPI_Type_contiguous(8, MPI_DOUBLE, &MPI_Particle);
    MPI_Type_commit(&MPI_Particle);

    std::vector<ParticleExchangeType> sendbuffer(m_sendbuffer.size());
    auto sendnum = static_cast<int>(m_sendbuffer.size());
    for(size_t i = 0; i < m_sendbuffer.size(); ++i){
        sendbuffer[i].mass = m_sendbuffer[i].getmass();
        sendbuffer[i].erot = m_sendbuffer[i].getRotationalEnergy();
        sendbuffer[i].x = m_sendbuffer[i].getposition()(0);
        sendbuffer[i].y = m_sendbuffer[i].getposition()(1);
        sendbuffer[i].z = m_sendbuffer[i].getposition()(2);
        sendbuffer[i].u = m_sendbuffer[i].getvelocity()(0);
        sendbuffer[i].v = m_sendbuffer[i].getvelocity()(1);
        sendbuffer[i].w = m_sendbuffer[i].getvelocity()(2);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(int step = 1; step < numprocs; ++step){
        int sendid = myid + step;
        if(sendid > numprocs - 1){
            sendid -= numprocs;
        }
        int recvid = myid - step;
        if(recvid < 0){
            recvid += numprocs;
        }

        // exchange send num and store it
        int recvnum;
        MPI_Sendrecv(&sendnum, 1, MPI_INT, sendid, 0, &recvnum, 1, MPI_INT, recvid, 0, m_cartesian_comm, MPI_STATUS_IGNORE);
        std::vector<ParticleExchangeType> recvbuffer(recvnum);
        // exchange data
        MPI_Sendrecv(sendbuffer.data(), sendnum, MPI_Particle, sendid, 1, recvbuffer.data(), recvnum, MPI_Particle, recvid, 1, m_cartesian_comm, MPI_STATUS_IGNORE);
        // store the data
        m_recvbuffer.reserve(m_recvbuffer.size() + recvnum);
        for(int i = 0; i < recvnum; ++i){
            // check if the particle is in the local domain
            if(recvbuffer[i].x > m_mesh->getoffsetX() * m_mesh->getUnidX() &&
            recvbuffer[i].x < (m_mesh->getnumberCellsX() + m_mesh->getoffsetX()) * m_mesh->getUnidX() &&
            recvbuffer[i].y > m_mesh->getoffsetY() * m_mesh->getUnidY() &&
            recvbuffer[i].y < (m_mesh->getnumberCellsY() + m_mesh->getoffsetY()) * m_mesh->getUnidY()){
                Eigen::Vector3d position(recvbuffer[i].x, recvbuffer[i].y, recvbuffer[i].z);
                Eigen::Vector3d velocity(recvbuffer[i].u, recvbuffer[i].v, recvbuffer[i].w);
                auto Erot {recvbuffer[i].erot};
                m_recvbuffer.emplace_back(mass, position, velocity, Erot);
            }
        }
    }
}

void CartesianParallel::setsendbuffer(std::vector<Particle>& sendbuffer)
{
    m_sendbuffer.insert(m_sendbuffer.end(),sendbuffer.begin(), sendbuffer.end());
    sendbuffer.clear();
}

void CartesianParallel::writerecvbuffer(std::vector<Particle>& m_particles)
{   
    // std::cout <<"我出去了 : " << m_sendbuffer.size()<<" 个粒子"<<
    //             "我进来了 : " <<m_recvbuffer.size() <<" 个粒子"<<std::endl;
    m_particles.insert(m_particles.end(), m_recvbuffer.begin(), m_recvbuffer.end());
    m_recvbuffer.clear();
    m_sendbuffer.clear();
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
