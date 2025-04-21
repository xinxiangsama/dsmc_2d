#include "Cell.h"
#include <iostream>
extern double Vstd;
extern double Vmax;
extern double VHS_coe;
extern std::unique_ptr<Random> randomgenerator;

const Cell::Coord &Cell::getposition() const
{
    return m_position;
}
const std::vector<Particle*> &Cell::getparticles() const
{
    return m_particles;
}
const Phase *Cell::getphase() const
{
    return m_phase.get();
}
Element *Cell::getelement() const
{
    return m_element;
}
const Cell::Coord &Cell::getindex() const
{
    return m_index;
}
const std::vector<std::shared_ptr<Cell>> &Cell::getchildren() const
{
    return m_children;
}
const int Cell::getCollisionNum()
{
    return N_collision;
}
const double &Cell::getdt()
{
    return m_dt;
}
bool Cell::ifcut()
{
    return m_element->ifcut();
}
void Cell::setposition(const Coord &position)
{
    m_position = position;
}
const double &Cell::getmfp()
{
    return m_mfp;
}
const AMRlevel &Cell::getAMRlevel()
{
    return m_level;
}
void Cell::setelement(Element *element)
{
    m_element = element;
}
void Cell::setindex(const Coord &index)
{
    m_index = index;
}
void Cell::setAMRlevel(const AMRlevel &level)
{
    m_level = level;
}
void Cell::allocatevar()
{
    m_phase = std::make_shared<Phase>();
}
void Cell::insertparticle(Particle*  particle)
{
    m_particles.push_back(particle);
}
void Cell::removeparticle(Particle*  particle)
{
    auto it = std::remove(m_particles.begin(), m_particles.end(), particle);
    m_particles.erase(it);
}

void Cell::removeallparticles()
{
    m_particles.clear();
    N_particles = 0;
}
void Cell::VTS()
{   
    auto cellvolume = m_element->getvolume();
    auto temperature = m_phase->gettemperature();
    auto vnorm = m_phase->getvelocity().norm();
    m_mfp = 1 / (sqrt(2) * (N_particles / cellvolume) * M_PI * diam * diam);
    m_mps = sqrt(2 * boltz * temperature / mass);
    m_dt = 0.5 * (m_mfp / (vnorm + m_mps));

    m_weight = Fn * m_dt / tau;
}
void Cell::comtimetokenleaving(Particle *particle)
{
    auto m_L1 = m_element->getL1();
    auto m_L2 = m_element->getL2();
    double xmin = m_position(0) - 0.5 * m_L1;
    double xmax = m_position(0) + 0.5 * m_L1;
    double ymin = m_position(1) - 0.5 * m_L2;
    double ymax = m_position(1) + 0.5 * m_L2;

    auto u = particle->getvelocity()[0];
    auto x = particle->getposition()[0];
    auto v = particle->getvelocity()[1];
    auto y = particle->getposition()[1];

    auto tx = u > 0 ? (xmax - x) / u : (xmin - x) / u;
    auto ty = v > 0 ? (ymax - y) / v : (ymin - y) / v;

    particle->settmove(std::min<double>(tx, ty));
}

// only Level1 cell perform this
void Cell::genAMRmesh()
{   
    // save old AMR mesh
    std::vector<std::shared_ptr<Cell>> m_oldchildren;
    m_oldchildren.insert(m_oldchildren.end(), 
                         std::make_move_iterator(m_children.begin()), 
                         std::make_move_iterator(m_children.end()));
    m_children.clear();
    
    double L1x = m_element->getL1();
    double L1y = m_element->getL2();
    // step 1 : Loop over all L3 cells within each L1 cell, find the maximum mean free path λmax
    double Maxmfp {};
    Maxmfp = findMaxmfpOverAllchild();
    // step 2 : Set the new L2 cell size to Chλmax  (here, Ch is a constant to satisfy the restriction on cell size,  usually Ch = 0.5 as discussed before).
    double Chx = 0.5;
    double Chy = 0.5;
    double L2x = Chx * Maxmfp;
    double L2y = Chy * Maxmfp;
    // Check if L1 can be evenly divided by L2
    if (L1x / L2x != static_cast<int>(L1x / L2x)) {
        // Adjust Ch so that L1 can be evenly divided by L2
        Chx = L1x / std::ceil(L1x / (Maxmfp * Chx));
        L2x = Chx * Maxmfp; // Recalculate L2 with the new Ch value
    }
    // same operation for other dimention
    if (L1y / L2y != static_cast<int>(L1y / L2y)) {
        Chy = L1y / std::ceil(L1y / (Maxmfp * Chy));
        L2y = Chy * Maxmfp; 
    }
    auto NxLv2 = static_cast<int>(L1x / L2x);
    auto NyLv2 = static_cast<int>(L1y / L2y);
    // step 3 : generate a uniform new L2 Cartesian grid  within each L1 cell
    m_element->genAMRmesh(NxLv2, NyLv2, L2x, L2y);
    auto& Lv2Elements = m_element->getchildren();
    for(int i = 0; i < NxLv2; ++i){
        for(int j = 0; j < NyLv2; ++j){
            int index = j + i * NyLv2;
            auto childcell = std::make_shared<Cell>();
            auto& element = Lv2Elements[index];
            childcell->setindex({i, j});
            childcell->setelement(element.get());
            childcell->setposition(element->getposition());
            childcell->setAMRlevel(AMRlevel::Lv2);
            insertchildern(childcell);
        }
    }

    // step 4 : Loop over each new L2 cell and find the minimum mean free path λmin within the  new L2 cell. λmin within the new L2 cell will be obtained from λ of all the old L3  cells which intersect with the new L2 cell.
    for(auto& childcell : m_children){
        double Minmfp {};
        Minmfp = findMinmfpOverL2child(childcell, m_oldchildren);
        double L2Chx = 0.5;
        double L2Chy = 0.5;
        double L3x = L2Chx * Minmfp;
        double L3y = L2Chy * Minmfp;
        // Check if L2 can be evenly divided by L3
        if (L2x / L3x != static_cast<int>(L2x / L3x)) {
            // Adjust Ch so that L2 can be evenly divided by L3
            L2Chx = L2x / std::ceil(L2x / (Minmfp * L2Chx));
            L3x = L2Chx * Minmfp; // Recalculate L3 with the new Ch value
        }
        // same operation for other dimention
        if (L2y / L3y != static_cast<int>(L2y / L3y)) {
            L2Chy = L2y / std::ceil(L2y / (Minmfp * L2Chy));
            L3y = L2Chy * Minmfp; 
        }
        auto NxLv3 = static_cast<int>(L2x / L3x);
        auto NyLv3 = static_cast<int>(L2y / L3y);
        // step 5 : Set the new L3 cell size within each new L2 cell to Chλmin.
        childcell->getelement()->genAMRmesh(NxLv3, NyLv3, L3x, L3y);
        auto Lv3Elements = childcell->getelement()->getchildren();
        for(int i = 0; i < NxLv3; ++i){
            for(int j = 0; j < NyLv3; ++j){
                int index = j + i * NyLv3;
                auto L3cell = std::make_shared<Cell>();
                auto& element = Lv3Elements[index];
                L3cell->setindex({i, j});
                L3cell->setelement(element.get());
                L3cell->setposition(element->getposition());
                L3cell->setAMRlevel(AMRlevel::Lv3);
                childcell->insertchildern(L3cell);
            }
        }
    }
    
}
void Cell::insertchildern(std::shared_ptr<Cell> child)
{
    m_children.emplace_back(std::move(child));
}
double Cell::findMaxmfpOverAllchild()
{   
    double res{-999.0};
    if(!m_children.empty()){
        for(auto& child : m_children){
            auto Maxchild {child->findMaxmfpOverAllchild()};
            res = Maxchild > res ? Maxchild : res;
        }
    }else{
        return m_mfp;
    }

    return res;
}
double Cell::findMinmfpOverL2child(std::shared_ptr<Cell> childcell, std::vector<std::shared_ptr<Cell>>& oldchildcells)
{   
    double res {999.0};
    for(auto& oldchildcell : oldchildcells){
        auto oldLv3cells = oldchildcell->getchildren();
        for(auto& Lv3cell : oldLv3cells){
            if(childcell->getelement()->isIntersecting(Lv3cell->getelement())){
                res = std::min(res, Lv3cell->getmfp());
            }
        }
    }

    return res;
}
void Cell::collision()
{   
    N_particles = m_particles.size();
    N_collision = 0;
    if (N_particles < 2) {
        // std::cerr << "[Warning] Cell @ " << m_index[0] << ", " << m_index[1]
        //           << " has fewer than 2 particles, skipping collision." << std::endl;
        return;
    }
    // until the smallest cell is reached
    if(m_children.size() != 0)
    {
        for(auto &child : m_children)
        {
            child->collision();
        }
    }else{
        auto cell_volume = m_element->getvolume();
        auto Srcmax = Vmax * M_PI * diam * diam;

        double expected_Mcand = static_cast<double>(N_particles * N_particles) * Fn * Srcmax * tau / (2 * cell_volume);
        int M_candidate = static_cast<int>(expected_Mcand + Mcand_r);
        Mcand_r = expected_Mcand - M_candidate;

        for(int i = 0; i < M_candidate; ++i){
            auto index1 = randomgenerator->getrandomint(0, N_particles - 1);
            auto index2 = randomgenerator->getrandomint(0, N_particles - 1);

            while(index1 == index2){
                index2 = randomgenerator->getrandomint(0, N_particles - 1);
            }
            if(index1 != index2){
                auto particle1 = m_particles[index1];
                auto particle2 = m_particles[index2];

                Particle::Coord v_rel = static_cast<Particle::Coord>(particle1->getvelocity() - particle2->getvelocity());

                auto v_rel_mag = v_rel.norm();
                auto Src = v_rel_mag * diam * diam * M_PI * (pow((2.0 * boltz * T / (0.5 * mass * pow(v_rel_mag, 2))),(Vtl-0.5))) / VHS_coe;
                if(Src > Srcmax){
                    Srcmax = Src;
                }
    
                auto rand01 = randomgenerator->getrandom01();
                if(rand01 < (Src / Srcmax)){
                    particle1->Collision(particle2);
                    // auto Vmean = Particle::Coord(
                    //     0.5 * (particle1->getvelocity()[0] + particle2->getvelocity()[0]),
                    //     0.5 * (particle1->getvelocity()[1] + particle2->getvelocity()[1]),
                    //     0.5 * (particle1->getvelocity()[2] + particle2->getvelocity()[2])
                    // );
                    // Particle::Coord Vmean = 0.5 * static_cast<Particle::Coord>(particle1->getvelocity() + particle2->getvelocity());
                    // Particle::Coord Vmean = 0.5 * (particle1->getvelocity() + particle2->getvelocity());
                    // double cosr = 2.0 * randomgenerator->getrandom01() - 1.0;
                    // double sinr = sqrt(1.0 - cosr * cosr);
                    // double phi = 2.0 * M_PI * randomgenerator->getrandom01();
                    // Particle::Coord Vrel_new(
                    //     cosr,
                    //     sinr * cos(phi),
                    //     sinr * sin(phi)
                    // );
        
                    // Vrel_new = Vrel_new.normalized() * v_rel_mag;

                    // particle1->setvelocity(Vmean + 0.5 * Vrel_new);
                    // particle2->setvelocity(Vmean - 0.5 * Vrel_new);

                    N_collision++;
                }
            }
        }

    }
}
void Cell::sample()
{
    N_particles = m_particles.size();
    // auto m_volume = m_element->getvolume();
    auto m_volume = (Volume) / (N1 * N2 * N3);

    auto ParticleNumdensity = static_cast<double>(N_particles * Fn) / m_volume;
    m_phase->setdensity(ParticleNumdensity * mass);

    // accumulate the mean and mean square of the velocity
    double U {}, V {}, W {}, U_squared {}, V_squared{}, W_squared{};
    for (auto &particle : m_particles){
        U += particle->getvelocity()[0];
        V += particle->getvelocity()[1];
        W += particle->getvelocity()[2];
        U_squared += particle->getvelocity()[0] * particle->getvelocity()[0];
        V_squared += particle->getvelocity()[1] * particle->getvelocity()[1];
        W_squared += particle->getvelocity()[2] * particle->getvelocity()[2];
    }
    U /= N_particles;
    V /= N_particles;
    W /= N_particles;
    U_squared /= N_particles;
    V_squared /= N_particles;
    W_squared /= N_particles;

    double Temperature = mass * ((U_squared + V_squared + W_squared) - (U * U + V * V + W * W)) / (3.0 * boltz);
    double Pressure = m_phase->getdensity() * ((U_squared + V_squared + W_squared) - (U * U + V * V + W * W)) / 3.0;

    m_phase->setvelocity({U, V, W});
    m_phase->settemperature(Temperature);
    m_phase->setpressure(Pressure);

    if(m_phase->getdensity() > 1.0e-2){
        std::cerr <<"density is error, particle num in this cell is : "<< N_particles<<
        " And cell volume is : "<<m_element->getvolume()<<std::endl;
    }
}
Cell::~Cell()
{
    removeallparticles();
}