#pragma once
#include "Geom.h"

class Circle : public Geom
{
public:
    Circle() = default;
    Circle(const int& num, const LargrangianPoint::Coord& center, const double& radius);
    virtual ~Circle() = default;

    virtual void Initialize();
    // judge the particle whether arrive the boundary
    virtual bool isHit(const Particle::Coord& position) const {};

    // define the action of particle when it hit the boundary
    virtual void Reflect(Particle* particle, const double& dt) const {};
protected:
    LargrangianPoint::Coord m_center;
    double m_radius;
};