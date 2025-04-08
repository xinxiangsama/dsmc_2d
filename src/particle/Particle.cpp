#include "Particle.h"

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
