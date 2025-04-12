#include "Segment.h"
#include "random/Random.h"
#include <memory>
extern std::unique_ptr<Random> randomgenerator;

bool Segment::isHit(const Particle::Coord &position) const
{   
    Particle::Coord m_point = m_leftpoint->getPosition();
    // Check if the particle is within a certain distance from the wall
    double distance = (position - m_point).dot(m_normal);
    return distance < 0.0;
}

void Segment::Reflect(Particle *particle, const double &dt) const
{
    Particle::Coord m_point = m_leftpoint->getPosition();
    auto velocity = particle->getvelocity();
    auto position = particle->getposition() - velocity * dt; //back to previous position
    auto distance = (position - m_point).dot(m_normal);
    auto v_normal = velocity.dot(m_normal);
    if(v_normal >= 0.0){
        return;
    }
    auto t_hit = - distance / v_normal;
    if(t_hit > dt || t_hit < 0.0){
        return;
    }

    auto hit_position = position + velocity * t_hit;
    // Specular reflection
    auto new_velocity = velocity  + 2 * (abs(velocity.dot(m_normal))) * m_normal;
    particle->setvelocity(new_velocity);
    // Update the position of the particle
    auto new_position = hit_position + particle->getvelocity() * (dt - t_hit);
    particle->setposition(new_position);
}

const double &Segment::getlength() const
{
    return m_length;
}

const LargrangianPoint::Coord &Segment::getsloop() const
{
    return m_sloop;
}

const LargrangianPoint::Coord &Segment::getnormal() const
{
    return m_normal;
}

LargrangianPoint *Segment::getleftpoint() const
{
    return m_leftpoint;
}

LargrangianPoint *Segment::getrightpoint() const
{
    return m_rightpoint;
}

void Segment::setlength(const double &length)
{
    m_length = length;
}

void Segment::setsloop(const Coord &sloop)
{
    m_sloop = sloop;
}

void Segment::setnormal(const Coord &normal)
{
    m_normal = normal;
}

void Segment::setleftpoint(LargrangianPoint *leftpoint)
{
    m_leftpoint = leftpoint;
}

void Segment::setrightpoint(LargrangianPoint *rightpoint)
{
    m_rightpoint = rightpoint;
}