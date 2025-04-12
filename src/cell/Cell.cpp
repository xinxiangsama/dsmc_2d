#include "Cell.h"
#include "../Param.h"
#include <iostream>
extern double Vstd;
extern double Vmax;
extern double VHS_coe;
extern std::unique_ptr<Random> randomgenerator;

const Cell::Coord &Cell::getposition() const
{
    return m_position;
}
const std::vector<std::shared_ptr<Particle>> &Cell::getparticles() const
{
    return m_particles;
}
const Phase *Cell::getphase() const
{
    return m_phase.get();
}
const Element *Cell::getelement() const
{
    return m_element;
}
const Cell::Coord &Cell::getindex() const
{
    return m_index;
}
const std::vector<std::unique_ptr<Cell>> &Cell::getchildren() const
{
    return m_children;
}
const int Cell::getCollisionNum()
{
    return N_collision;
}
void Cell::setposition(const Coord &position)
{
    m_position = position;
}
void Cell::setelement(Element* element)
{
    m_element = element;
}
void Cell::setindex(const Coord &index)
{
    m_index = index;
}
void Cell::allocatevar()
{
    m_phase = std::make_unique<Phase>();
}
void Cell::insertparticle(std::shared_ptr<Particle> particle)
{
    m_particles.push_back(particle);
}
void Cell::removeparticle(std::shared_ptr<Particle> particle)
{
    auto it = std::remove(m_particles.begin(), m_particles.end(), particle);
    m_particles.erase(it);
}

void Cell::removeallparticles()
{
    m_particles.clear();
    N_particles = 0;
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
                    particle1->Collision(particle2.get());
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
    auto m_volume = m_element->getvolume();

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
}
Cell::~Cell()
{
    removeallparticles();
}