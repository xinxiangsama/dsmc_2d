#pragma once
#include <Eigen/Dense>

class Particle 
{
public:
    using Coord = Eigen::Vector2d;
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
protected:
    double m_mass;
    Coord m_position;
    Coord m_velocity;
};