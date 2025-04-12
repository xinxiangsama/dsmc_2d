#pragma once
#include "Segment.h"
#include "../particle/Particle.h"

class Geom
{
public:
    Geom() = default;
    Geom(const size_t& num) : numLagrangianPoints(num) {};
    virtual ~Geom() = default;

    virtual void Initialize() = 0;
    // judge the particle whether arrive the boundary
    virtual bool isHit(const Particle::Coord& position) const = 0;

    // define the action of particle when it hit the boundary
    virtual void Reflect(Particle* particle, const double& dt) const = 0;


protected:
    size_t numLagrangianPoints;
    std::vector<LargrangianPoint> m_points;
    std::vector<Segment> m_segments;
};