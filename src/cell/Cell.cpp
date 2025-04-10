#include "Cell.h"
#include "../Param.h"
extern double Vstd;
extern double Vmax;
extern double VHS_coe;
extern std::unique_ptr<Random> randomgenerator;

const Cell::Coord &Cell::getposition() const
{
    return m_position;
}
const std::vector<Particle *> &Cell::getparticles() const
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
void Cell::insertparticle(Particle *particle)
{
    m_particles.push_back(particle);
}
void Cell::removeparticle(Particle *particle)
{
    auto it = std::remove(m_particles.begin(), m_particles.end(), particle);
    m_particles.erase(it);
}

void Cell::removeallparticles()
{
    m_particles.clear();
}
void Cell::collision()
{   
    N_particles = m_particles.size();
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

        auto M_candidate = static_cast<int>(N_particles * N_particles * Fn * Srcmax * tau / (2 * cell_volume));

        for(size_t i = 0; i < M_candidate; ++i){
            auto index1 = randomgenerator->getrandomint(0, N_particles - 1);
            auto index2 = randomgenerator->getrandomint(0, N_particles - 1);
            if(index1 == index2) continue;

            auto particle1 = m_particles[index1];
            auto particle2 = m_particles[index2];

            auto v_rel_mag = (particle1->getvelocity() - particle2->getvelocity()).norm();
            auto Src = v_rel_mag * diam * diam * M_PI * (pow((2.0 * boltz * T / (0.5 * mass * pow(v_rel_mag, 2))),(Vtl-0.5))) / VHS_coe;
            if(Src > Srcmax){
                Srcmax = Src;
            }

            auto rand01 = randomgenerator->getrandom01();
            if(rand01 < Src / Srcmax){
                particle1->Collision(particle2);
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
    double U {}, V {}, Usquared {}, V_squared{};
    for (auto &particle : m_particles){
        U += particle->getvelocity().x();
        V += particle->getvelocity().y();
        Usquared += particle->getvelocity().x() * particle->getvelocity().x();
        V_squared += particle->getvelocity().y() * particle->getvelocity().y();
    }
    U /= N_particles;
    V /= N_particles;
    Usquared /= N_particles;
    V_squared /= N_particles;

    double Temperature = mass * ((Usquared + V_squared) - (U * U + V * V)) / (3.0 * boltz);
    double Pressure = m_phase->getdensity() * ((Usquared + V_squared) - (U * U + V * V)) / 3.0;

    m_phase->setvelocity({U, V});
    m_phase->settemperature(Temperature);
    m_phase->setpressure(Pressure);
}
Cell::~Cell()
{
    removeallparticles();
}