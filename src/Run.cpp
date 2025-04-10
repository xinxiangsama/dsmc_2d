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
    m_mesh->setnumberCellsXGlobal(N1);
    m_mesh->setnumberCellsYGlobal(N2);


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

    /*boundary part*/
    bool ifdisfuse = false;
    inlet = std::make_unique<WallBoundary>(Eigen::Vector2d(0.0, L2 * 0.5), Eigen::Vector2d(1.0, 0.0), ifdisfuse);
    outlet = std::make_unique<WallBoundary>(Eigen::Vector2d(L1, L2 * 0.5), Eigen::Vector2d(-1.0, 0.0), ifdisfuse);
    wall1 = std::make_unique<WallBoundary>(Eigen::Vector2d(L1 * 0.5, 0.0), Eigen::Vector2d(0.0, 1.0), ifdisfuse);
    wall2 = std::make_unique<WallBoundary>(Eigen::Vector2d(L1 * 0.5, L2), Eigen::Vector2d(0.0, -1.0), ifdisfuse);

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
    if(N_Particle % N1 * N2 != 0){
        std::cerr <<"particle can't be devided by mesh" << std::endl;
    }
    
    auto numparticlepercell = numparticlelocal / (m_mesh->getnumberCellsX() * m_mesh->getnumberCellsY());
    for(auto& cell: m_cells){
        for(int i = 0; i < numparticlepercell; ++i){
            auto particle = std::make_unique<Particle>();
            particle->setmass(mass);
            auto rx = randomgenerator->getrandom01();
            auto ry = randomgenerator->getrandom01();
            double x = cell->getposition()(0) + (rx - 0.5) * m_mesh->getUnidX();
            double y = cell->getposition()(1) + (ry - 0.5) * m_mesh->getUnidY();
            particle->setposition(Eigen::Vector2d(x, y));
            auto velocity = randomgenerator->MaxwellDistribution(Vstd);
            // velocity(0) += V_jet;
            particle->setvelocity(velocity);
            m_particles.emplace_back(std::move(particle));
            cell->insertparticle(m_particles.back().get());
        }
    }
}

void Run::particlemove()
{   

    for(auto& particle : m_particles){
        particle->Move(tau);

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

    }
}

void Run::ressignParticle()
{   
    for(auto& cell : m_cells){
        cell->removeallparticles();
    }
    std::vector<std::unique_ptr<Particle>> particle_out;
    std::vector<std::unique_ptr<Particle>> particle_in;
    particle_out.reserve(m_particles.size());
    particle_in.reserve(m_particles.size());
    for(auto& particle : m_particles){
        auto position = particle->getposition();
        if(position(0) < m_mesh->getoffsetX() * m_mesh->getUnidX() || 
        position(0) > (m_mesh->getnumberCellsX() + m_mesh->getoffsetX()) * m_mesh->getUnidX()|| 
        position(1) < m_mesh->getoffsetY() * m_mesh->getUnidY() || 
        position(1) > (m_mesh->getnumberCellsY() + m_mesh->getoffsetY()) * m_mesh->getUnidY()){
            particle_out.emplace_back(std::move(particle));
        }else{
            particle_in.emplace_back(std::move(particle));
        }
    }


    m_particles.clear();
    m_particles.insert(m_particles.end(), std::make_move_iterator(particle_in.begin()), std::make_move_iterator(particle_in.end()));
    auto start_sendbuffer = std::chrono::high_resolution_clock::now();
    m_parallel->setsendbuffer(particle_out);
    auto end_sendbuffer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sendbuffer = end_sendbuffer - start_sendbuffer;

    auto start_exchangedata = std::chrono::high_resolution_clock::now();
    m_parallel->exchangedata();
    auto end_exchangedata = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_exchangedata = end_exchangedata - start_exchangedata;

    auto start_writerecvbuffer = std::chrono::high_resolution_clock::now();
    m_parallel->writerecvbuffer(m_particles);
    auto end_writerecvbuffer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_writerecvbuffer = end_writerecvbuffer - start_writerecvbuffer;

    auto start_assignParticle2cell = std::chrono::high_resolution_clock::now();
    assignParticle2cell();
    auto end_assignParticle2cell = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_assignParticle2cell = end_assignParticle2cell - start_assignParticle2cell;

    // if(myid == 0) {
    //     std::cout << "Time taken for setsendbuffer: " << elapsed_sendbuffer.count() << " seconds" << std::endl;
    //     std::cout << "Time taken for exchangedata: " << elapsed_exchangedata.count() << " seconds" << std::endl;
    //     std::cout << "Time taken for writerecvbuffer: " << elapsed_writerecvbuffer.count() << " seconds" << std::endl;
    //     std::cout << "Time taken for assignParticle2cell: " << elapsed_assignParticle2cell.count() << " seconds" << std::endl;
    // }
}

void Run::collision()
{
    for(auto& cell : m_cells){
        cell->collision();
    }
}

Cell* Run::locatecell(const Particle::Coord& position)
{
    auto i = static_cast<int>(position(0) / m_mesh->getUnidX()) - m_mesh->getoffsetX();
    auto j = static_cast<int>(position(1) / m_mesh->getUnidY()) - m_mesh->getoffsetY();
    if(i < 0 || i >= m_mesh->getnumberCellsX() || j < 0 || j >= m_mesh->getnumberCellsY()){
        std::cerr << "particle out of range" << std::endl;
        return nullptr;
    }

    return m_cells[j + i * m_mesh->getnumberCellsY()].get();
}

void Run::assignParticle2cell()
{
    for(auto& particle : m_particles){
        auto position = particle->getposition();
        auto cell = locatecell(position);
        if(cell != nullptr){
            cell->insertparticle(particle.get());
        }
    }
}

void Run::solver()
{
    for(size_t iter = 0; iter < 1000; ++iter){
        if(myid == 0){
            std::cout << "iter: " << iter << std::endl;
        }
        particlemove();
        // collision();
        ressignParticle();
        if (iter % 100 == 0) {
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
