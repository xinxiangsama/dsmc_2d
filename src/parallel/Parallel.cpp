#include "Parallel.h"

void Parallel::setMesh(Mesh *mesh)
{
    m_mesh = mesh;
}

void Parallel::setIdAndNumprocs(const int &_myid, const int &_numprocs)
{
    myid = _myid;
    numprocs = _numprocs;
}