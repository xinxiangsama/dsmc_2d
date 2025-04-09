#pragma once
#include "Boundary.h"

class OutflowBoundary : public Boundary
{
public:
    OutflowBoundary(const Particle::Coord& point, const Particle::Coord& normal);
    virtual ~OutflowBoundary() = default;

    virtual bool isHit(const Particle::Coord& position) const override;
    virtual void Reflect(Particle* particle, const double& dt) const override;
    virtual const Particle::Coord& getnormal() const override { return m_normal; }

protected:
    Particle::Coord m_point;   // a point in boundary
    Particle::Coord m_normal;  // normal
};