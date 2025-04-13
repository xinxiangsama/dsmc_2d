#pragma once
#include "Segment.h"
#include "../particle/Particle.h"
#include "meshes/Element.h"

class Geom
{
public:
    Geom() = default;
    Geom(const size_t& num) : numLagrangianPoints(num) {};
    virtual ~Geom() = default;

    virtual void Initialize() = 0;

    void SortSegment2Element(Element* element);
    std::unique_ptr<Eigen::Vector2d> cn_PnPolyX(const Eigen::Vector2d& P);
    std::unique_ptr<Eigen::Vector2d> cn_PnPolyY(const Eigen::Vector2d& P);

protected:
    size_t numLagrangianPoints = 0;
    std::vector<LargrangianPoint> m_points;
    std::vector<Segment> m_segments;
};