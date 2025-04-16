#include "Run.h"
#include <chrono>
double Vstd;
double Vmax;
double VHS_coe;
std::unique_ptr<Random> randomgenerator;

auto GammaFun = [](double xlen){
    double A, ylen, GAM;
    
    A = 1.0;
    ylen = xlen;
    
    if (ylen < 1.0) {
        A = A / ylen;
    } else {
        while (ylen >= 1.0) {
            ylen = ylen - 1;
            A = A * ylen;
        }
    }
    
    GAM = A * (1.0 - 0.5748 * ylen + 0.9512 * ylen * ylen - 0.6998 * ylen * ylen * ylen + 0.4245 * ylen * ylen * ylen * ylen - 0.1010 * ylen * ylen * ylen * ylen * ylen);
    
    return GAM;
};

void Run::initialize(int argc, char **argv)
{   
    /*cal base var*/
    Vstd = sqrt(2 * boltz * T / mass);
    Vmax = 2 * sqrt(8 / M_PI) * Vstd;
    VHS_coe = GammaFun(2.5 - Vtl);
    /*mesh part*/
    m_mesh = std::make_unique<CartesianMesh>();
    m_mesh->setGlobalLengthX(L1);
    m_mesh->setGlobalLengthY(L2);
    m_mesh->setGlobalLengthZ(L3);
    m_mesh->setnumberCellsXGlobal(N1);
    m_mesh->setnumberCellsYGlobal(N2);
    m_mesh->setnumberCellsZGlobal(N3);

    m_geom = std::make_unique<Circle>(1024, LargrangianPoint::Coord{Center_x, Center_y, 0.0}, Radius);
    // m_geom = std::make_unique<Square>(4, LargrangianPoint::Coord{Center_x, Center_y, 0.0}, Radius);
    m_geom->Initialize();

    /*parallel part*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    m_parallel = std::make_unique<CartesianParallel>();
    m_parallel->setMesh(m_mesh.get());
    m_parallel->setIdAndNumprocs(myid, numprocs);
    m_parallel->ZoneDcomposition();
    m_parallel->setNeibours();
    m_parallel->info();

    /*second mesh part*/
    m_mesh->allocateCells(m_cells);
    m_mesh->setelement();
    m_mesh->BindElementwithFace();
    m_mesh->BindCellwithElement(m_cells);
    m_mesh->cutcell(m_geom.get());
    /*boundary part*/
    bool ifdisfuse = false;
    inlet = std::make_unique<PeriodicBoundary>(Eigen::Vector3d(0.0, 0.5 * L2, 0.5 * L3), Eigen::Vector3d(1.0, 0.0, 0.0), 0, L1);
    outlet = std::make_unique<PeriodicBoundary>(Eigen::Vector3d(L1, 0.5 * L2, 0.5 * L3), Eigen::Vector3d(-1.0, 0.0, 0.0), 0, L1);
    // inlet = std::make_unique<WallBoundary>(Eigen::Vector3d(0.0, 0.5 * L2, 0.5 * L3), Eigen::Vector3d(1.0, 0.0, 0.0), ifdisfuse);
    // outlet = std::make_unique<WallBoundary>(Eigen::Vector3d(L1, 0.5 * L2, 0.5 * L3), Eigen::Vector3d(-1.0, 0.0, 0.0), ifdisfuse);
    wall1 = std::make_unique<WallBoundary>(Eigen::Vector3d(0.5 * L1, 0.0, 0.5 * L3), Eigen::Vector3d(0.0, 1.0, 0.0), ifdisfuse);
    wall2 = std::make_unique<WallBoundary>(Eigen::Vector3d(0.5 * L1, L2, 0.5 * L3), Eigen::Vector3d(0.0, -1.0, 0.0), ifdisfuse);
    // wall3 = std::make_unique<WallBoundary>(Eigen::Vector3d(0.5 * L1, 0.5 * L2, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0), ifdisfuse);
    // wall4 = std::make_unique<WallBoundary>(Eigen::Vector3d(0.5 * L1, 0.5 * L2, L3), Eigen::Vector3d(0.0, 0.0, -1.0), ifdisfuse);
    wall3 = std::make_unique<PeriodicBoundary>(Eigen::Vector3d(0.5 * L1, 0.5 * L2, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0), 2, L3);
    wall4 = std::make_unique<PeriodicBoundary>(Eigen::Vector3d(0.5 * L1, 0.5 * L2, L3), Eigen::Vector3d(0.0, 0.0, -1.0), 2, L3);
    
    /*random part*/
    randomgenerator = std::make_unique<Random>();
    /*output part*/
    m_output = std::make_unique<Output>(this);
    /*Initial particle phase*/
    assignParticle();
    for(auto& cell : m_cells){
        cell->allocatevar();
    }

    for(auto& cell : m_cells){
        cell->sample();
    }
    m_output->Write2HDF5("./res/init.h5");
    if(myid == 0){
        std::cout << "MPI Initialized" << std::endl;
        std::cout << "Mesh Initialized" << std::endl;
        std::cout << "Parallel Initialized" << std::endl;
        std::cout << "Boundary Initialized" << std::endl;
        std::cout << "Random Initialized" << std::endl;
        std::cout << "Output Initialized" << std::endl;
    }
}

void Run::assignParticle()
{
    numparticlelocal = N_Particle / numprocs;
    m_particles.reserve(numparticlelocal);
    if(N_Particle % (N1 * N2 * N3) != 0){
        std::cerr <<"particle can't be devided by mesh" << std::endl;
    }
    std::random_device rd;
    auto numparticlepercell = numparticlelocal / (m_mesh->getnumberCellsX() * m_mesh->getnumberCellsY() * m_mesh->getnumberCellsZ());
    for(auto& cell: m_cells){
        std::mt19937 gen(rd() + cell->getindex()[0] + cell->getindex()[1] + cell->getindex()[2]);
        for(int i = 0; i < numparticlepercell; ++i){
            auto particle = std::make_shared<Particle>();
            particle->setmass(mass);
            auto rx = randomgenerator->getrandom01();
            auto ry = randomgenerator->getrandom01();
            auto rz = randomgenerator->getrandom01();
            // double rx = {std::uniform_real_distribution<double>(0, 1)(gen)};
            // double ry = {std::uniform_real_distribution<double>(0, 1)(gen)};
            // double rz = {std::uniform_real_distribution<double>(0, 1)(gen)};
            double x = cell->getposition()(0) + (rx - 0.5) * m_mesh->getUnidX();
            double y = cell->getposition()(1) + (ry - 0.5) * m_mesh->getUnidY();
            double z = cell->getposition()(2) + (rz - 0.5) * m_mesh->getUnidZ();
            if((x - Center_x) * (x - Center_x) + (y - Center_y) * (y - Center_y) <= 1.5*(Radius * Radius)){
                continue;
            }
            particle->setposition(Eigen::Vector3d(x, y, z));
            auto velocity = randomgenerator->MaxwellDistribution(Vstd);
            
            velocity(0) += V_jet;
            // if(x <= L1 * 0.5){
            //     velocity(0) += V_jet;
            // }else{
            //     velocity(0) -= V_jet;
            // }
            particle->setvelocity(velocity);
            m_particles.push_back(particle);
            cell->insertparticle(particle);
        }
    }
}

void Run::particlemove()
{   
    // for(auto& particle : m_particles){
    //     particle->Move(tau);

    //     if(inlet->isHit(particle->getposition())){
    //         inlet->Reflect(particle.get(), tau);
    //     }
    //     if(outlet->isHit(particle->getposition())){
    //         outlet->Reflect(particle.get(), tau);
    //     }

    //     if(wall1->isHit(particle->getposition())){
    //         wall1->Reflect(particle.get(), tau);
    //     }
    //     if(wall2->isHit(particle->getposition())){
    //         wall2->Reflect(particle.get(), tau);
    //     }

    //     if(wall3->isHit(particle->getposition())){
    //         wall3->Reflect(particle.get(), tau);
    //     }
    //     if(wall4->isHit(particle->getposition())){
    //         wall4->Reflect(particle.get(), tau);
    //     }

    // }

    for(auto& cell : m_cells){
        for(auto& particle : cell->getparticles()){
            particle->Move(tau);

            if(cell->ifcut()){
                for(auto& segment : cell->getelement()->getsegments()){
                    if(segment->isHit(particle->getposition())){
                        // segment->Reflect(particle, tau);
                        // break;
                        particle->setposition({999.0, 999.0, 999.0});
                    }
                    // particle->setposition({999.0, 999.0, 999.0});
                }

                // particle->setposition({999.0, 999.0, 999.0});
            }

            if(inlet->isHit(particle->getposition())){
                inlet->Reflect(particle.get(), tau);
            }
            if(outlet->isHit(particle->getposition())){
                outlet->Reflect(particle.get(), tau);
            }

            if(wall1->isHit(particle->getposition())){
                wall1->Reflect(particle.get(), tau);
            }
            if(wall2->isHit(particle->getposition())){
                wall2->Reflect(particle.get(), tau);
            }

            if(wall3->isHit(particle->getposition())){
                wall3->Reflect(particle.get(), tau);
            }
            if(wall4->isHit(particle->getposition())){
                wall4->Reflect(particle.get(), tau);
            }
        }
    }
}

void Run::ressignParticle()
{   
    for(auto& cell : m_cells){
        cell->removeallparticles();
    }
    std::vector<std::shared_ptr<Particle>> particle_out;
    std::vector<std::shared_ptr<Particle>> particle_in;
    particle_out.reserve(m_particles.size());
    particle_in.reserve(m_particles.size());
    for(auto& particle : m_particles){
        auto position = particle->getposition();
        if(position(0) < m_mesh->getoffsetX() * m_mesh->getUnidX() || 
        position(0) > (m_mesh->getnumberCellsX() + m_mesh->getoffsetX()) * m_mesh->getUnidX()|| 
        position(1) < m_mesh->getoffsetY() * m_mesh->getUnidY() || 
        position(1) > (m_mesh->getnumberCellsY() + m_mesh->getoffsetY()) * m_mesh->getUnidY()){
            particle_out.push_back(particle);
        }else{
            particle_in.push_back(particle);
        }
    }


    m_particles.clear();
    m_particles.insert(m_particles.end(), particle_in.begin(), particle_in.end());;
    m_parallel->setsendbuffer(particle_out);

    m_parallel->exchangedata();

    m_parallel->writerecvbuffer(m_particles);

    assignParticle2cell();

    int N_particle_local = m_particles.size();
    int N_particle_global {};
    MPI_Reduce(&N_particle_local, &N_particle_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid == 0) {
        std::cout << "Total number of particles: " << N_particle_global << std::endl;
    }
}

void Run::collision()
{
    int Ncollision_local {};
    for(auto& cell : m_cells){
        cell->collision();
        Ncollision_local += cell->getCollisionNum();
    }

    int Ncollision_global {};
    MPI_Reduce(&Ncollision_local, &Ncollision_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid == 0) {
        std::cout << "Total number of collisions: " << Ncollision_global << std::endl;
    }
}

Cell* Run::locatecell(const Particle::Coord& position)
{
    auto i = static_cast<int>(position(0) / m_mesh->getUnidX()) - m_mesh->getoffsetX();
    auto j = static_cast<int>(position(1) / m_mesh->getUnidY()) - m_mesh->getoffsetY();
    auto k = static_cast<int>(position(2) / m_mesh->getUnidZ());
    // k = std::max(0, std::min(k, m_mesh->getnumberCellsZ() - 1));
    if(i < 0 || i >= m_mesh->getnumberCellsX() || j < 0 || j >= m_mesh->getnumberCellsY()){
        std::cerr << "particle out of range at position: (" 
              << position(0) << ", " 
              << position(1) << ", " 
              << position(2) << ")" << std::endl;
        return nullptr;
    }

    if(Eigen::Vector3d{i,j,k} != m_cells[k + j * m_mesh->getnumberCellsZ() + i * m_mesh->getnumberCellsY() * m_mesh->getnumberCellsZ()]->getindex()){
        std::cerr <<"error!"<<std::endl;
    }

    return m_cells[k + j * m_mesh->getnumberCellsZ() + i * m_mesh->getnumberCellsY() * m_mesh->getnumberCellsZ()].get();
}

void Run::assignParticle2cell()
{
    for(auto& particle : m_particles){
        auto position = particle->getposition();
        auto cell = locatecell(position);
        if(cell != nullptr){
            cell->insertparticle(particle);
        }
    }
}

void Run::solver()
{
    for(size_t iter = 0; iter < 10000; ++iter){
        if(myid == 0){
            std::cout << "iter: " << iter << std::endl;
        }
        particlemove();
        collision();
        ressignParticle();
        if (iter % 50 == 0) {
            for(auto& cell : m_cells){
                cell->sample();
            }
            m_output->Write2HDF5("./res/step" + std::to_string(iter) + ".h5");
        }
    }
    if(myid == 0){
        std::cout << "Simulation Finished" << std::endl;
    }
    m_output->Write2HDF5("./res/final.h5");
    if(myid == 0){
        std::cout << "Output Finished" << std::endl;
    }
}

void Run::finalize()
{
    m_cells.clear();
    MPI_Finalize();
    if(myid == 0){
        std::cout << "MPI Finalized" << std::endl;
    }
}
