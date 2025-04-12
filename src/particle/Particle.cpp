#include "Particle.h"
#include <memory>

extern std::unique_ptr<Random> randomgenerator;

const double &Particle::getmass()
{
    return m_mass;
}

const Particle::Coord &Particle::getposition()
{
    return m_position;
}
const Particle::Coord &Particle::getvelocity()
{
    return m_velocity;
}
void Particle::setmass(const double &mass)
{
    m_mass = mass;
}
void Particle::setposition(const Coord &position)
{
    m_position = position;
}
void Particle::setvelocity(const Coord &velocity)
{
    m_velocity = velocity;
}

void Particle::Move(const double &dt)
{
    m_position += m_velocity * dt;
}

void Particle::Collision(Particle *other)
{
    auto v_mean = (m_velocity + other->getvelocity()) * 0.5;
    auto v_rel_mag = (m_velocity - other->getvelocity()).norm();
    auto rand01 = randomgenerator->getrandom01();
    auto cos_theta = 2.0 * rand01 - 1.0;
    auto sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    rand01 = randomgenerator->getrandom01();
    auto phi = 2.0 * M_PI * rand01;

    auto v_rel = v_rel_mag * Eigen::Vector3d(cos_theta, sin_theta * cos(phi), sin_theta * sin(phi));

    m_velocity = v_mean + v_rel * 0.5;
    other->setvelocity(v_mean - v_rel * 0.5);
}
