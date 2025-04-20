#pragma once
#include <mpi.h>
#include <memory>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "parallel/CartesianParallel.h"
#include "meshes/CartesianMesh.h"
#include "boundary/WallBoundary.h"
#include "boundary/OutletBoundary.h"
#include "boundary/PeriodicBoundary.h"
#include "boundary/InletBoundary.h"
#include "io/Output.h"
#include "cell/Cell.h"
#include "object/Circle.h"
#include "object/Square.h"

class Output;
class Run
{
public:
    Run() = default;
    ~Run() = default;

    /*============Initialize==========*/
    void initialize(int argc, char** argv);
    /*============Main calculate loop===========*/
    void solver();
    /*============Finalize===========*/
    void finalize();

    void assignParticle(const double& coef);
    void particlemove();
    void ressignParticle();
    void collision();
    void assignParticle2cell();

protected:
    std::vector<Cell> m_cells;
    std::vector<Particle> m_particles;
    size_t numparticlelocal;

    /*parallel*/
    int myid;
    int numprocs;
    std::unique_ptr<Parallel> m_parallel;
    /*mesh*/
    std::unique_ptr<Mesh> m_mesh;
    std::unique_ptr<Geom> m_geom;
    std::unique_ptr<Boundary> inlet;
    std::unique_ptr<Boundary> outlet;
    std::unique_ptr<Boundary> wall1; // bottom wall
    std::unique_ptr<Boundary> wall2; // top wall
    std::unique_ptr<Boundary> wall3; // front wall
    std::unique_ptr<Boundary> wall4; // back wall

    /*output*/
    std::unique_ptr<Output> m_output;

    friend class Output;
};
