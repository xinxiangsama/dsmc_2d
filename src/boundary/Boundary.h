#pragma once
#include "../particle/Particle.h"

class Boundary
{
public:
    Boundary() = default;
    virtual ~Boundary() = default;

    // judge the particle whether arrive the boundary
    virtual bool isHit(const Particle::Coord& position) const = 0;

    // define the action of particle when it hit the boundary
    virtual void Reflect(Particle* particle, const double& dt) const = 0;

    virtual bool isPeriodic() const { return false; }

    virtual const Particle::Coord& getnormal() const {}
};
