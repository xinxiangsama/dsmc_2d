#pragma once
#include <Eigen/Dense>
#include "../random/Random.h"

class Particle 
{
public:
    using Coord = Eigen::Vector3d;
    Particle() = default;
    Particle(const double& mass, const Coord& position, const Coord& velocity)
        : m_mass(mass), m_position(position), m_velocity(velocity) {}
    virtual ~Particle() = default;

    // Accessors
    const double& getmass();
    const Coord& getposition();
    const Coord& getvelocity();
    // Modifiers
    void setmass(const double& mass);
    void setposition(const Coord& position);
    void setvelocity(const Coord& velocity);
    // Functions
    virtual void Move(const double& dt);
    virtual void Collision(Particle* other);
protected:
    double m_mass;
    Coord m_position;
    Coord m_velocity;
};

struct ParticleExchangeType
{
public:
    double mass;
    double x;
    double y;
    double z;
    double u;
    double v;
    double w;
};