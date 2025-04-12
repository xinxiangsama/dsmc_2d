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
    bool segmentIntersectsBox(const std::array<double, 2>& p1, const std::array<double, 2>& p2,
        double xmin, double xmax, double ymin, double ymax);


protected:
    size_t numLagrangianPoints = 0;
    std::vector<LargrangianPoint> m_points;
    std::vector<Segment> m_segments;
};