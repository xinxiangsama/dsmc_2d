#include "Run.h"
double Vstd;
double Vmax;
std::unique_ptr<Random> randomgenerator;


void Run::initialize(int argc, char **argv)
{   
    /*cal base var*/
    Vstd = sqrt(2 * boltz * T / mass);
    Vmax = 2 * sqrt(8 / M_PI) * Vstd;

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
    /*Initial particle phase*/
    assignParticle();
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
            velocity(0) += V_jet;
            particle->setvelocity(velocity);
            m_particles.emplace_back(std::move(particle));
            cell->insertparticle(m_particles.back().get());
        }
    }
}

void Run::particlemove()
{   
    if(myid == 0){
        std::cout <<"particle position before move" << m_particles[1]->getposition()<<std::endl;
    }

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

    if(myid == 0){
        std::cout <<"particle position after move" << m_particles[1]->getposition()<<std::endl;
    }
}

void Run::ressignParticle()
{   
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

    std::cout <<"particle_out size" << particle_out.size() << std::endl;
    std::cout <<"particle_in size" << particle_in.size() << std::endl;

    m_particles.clear();
    m_particles.insert(m_particles.end(), std::make_move_iterator(particle_in.begin()), std::make_move_iterator(particle_in.end()));
    m_parallel->setsendbuffer(particle_out);
    m_parallel->exchangedata();
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

void Run::solver()
{
    particlemove();
    ressignParticle();
}

void Run::finalize()
{
    m_cells.clear();
    MPI_Finalize();
    if(myid == 0){
        std::cout << "MPI Finalized" << std::endl;
    }
}
