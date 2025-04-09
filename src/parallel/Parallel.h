#pragma once
#include <mpi.h>
#include "../Param.h"
#include "../meshes/Mesh.h"
#include <vector>

class Parallel
{
public:
    Parallel() = default;
    ~Parallel() = default;
    void setMesh(Mesh* mesh);
    void setIdAndNumprocs(const int& _myid, const int& _numprocs);
    virtual void ZoneDcomposition() = 0; // subclass must have one!
    virtual void setNeibours() = 0;
    virtual void exchangedata() {};

    // Access
    virtual const int& getLeftNeibour() {};
    virtual const int& getRightNeibour() {};
    virtual const int& getTopNeibour() {};
    virtual const int& getBottomNeibour() {};
    virtual void info() {};
protected:
    Mesh* m_mesh = nullptr;
    int myid = 0;
    int numprocs = 0;
};