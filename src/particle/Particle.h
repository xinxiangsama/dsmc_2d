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
    const bool& ifvalid() {return isvalid;}
    const int& getcellID();
    const double& gettmove();
    // Modifiers
    void setmass(const double& mass);
    void setposition(const Coord& position);
    void setvelocity(const Coord& velocity);
    void settmove(const double& t);
    void setInvalid() {isvalid = false;}
    void setcellID(const int& id);
     // Functions
    virtual void Move(const double& dt);
    virtual void Collision(Particle* other);
protected:
    double m_mass;
    Coord m_position;
    Coord m_velocity;
    bool isvalid {true};

    double t_move; // time token to leave current cell
    int m_cellID;
    int m_globalID;
    int m_localID;
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