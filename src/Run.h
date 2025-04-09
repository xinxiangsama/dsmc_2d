#pragma once
#include <mpi.h>
#include <memory>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "parallel/CartesianParallel.h"
#include "meshes/CartesianMesh.h"
#include "cell/Cell.h"

class Run
{
public:
    Run() = default;
    ~Run() = default;

    /*============Initialize==========*/
    void initialize(int argc, char** argv);
    void finalize();

protected:
    std::vector<std::unique_ptr<Cell>> m_cells;
    std::vector<std::unique_ptr<Particle>> m_particles;

    /*parallel*/
    int myid;
    int numprocs;
    std::unique_ptr<Parallel> m_parallel;
    /*mesh*/
    std::unique_ptr<Mesh> m_mesh;
};
