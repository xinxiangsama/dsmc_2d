#include "OutflowBoundary.h"

OutflowBoundary::OutflowBoundary(const Particle::Coord &point, const Particle::Coord &normal) : m_point(point), m_normal(normal)
{
    m_normal.normalized();
}

bool OutflowBoundary::isHit(const Particle::Coord &position) const
{
    // Check if the particle is within a certain distance from the wall
    double distance = (position - m_point).dot(m_normal);
    return distance < 0.0;
}

void OutflowBoundary::Reflect(Particle *particle, const double& dt) const
{
    particle->setvelocity(Eigen::Vector2d(0.0, 0.0));
}
