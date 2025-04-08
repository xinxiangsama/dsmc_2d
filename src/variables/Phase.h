#pragma once
#include <Eigen/Dense>

class Phase
{
public:
    using Coord = Eigen::Vector2d;
    Phase() = default;
    virtual ~Phase() = default;

protected:
    double m_density;
    Coord m_velocity;
    double m_pressure;
    double m_temperature;
};