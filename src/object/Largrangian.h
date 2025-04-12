#pragma once
#include <Eigen/Dense>
#include "Segment.h"

class Segment;
class LargrangianPoint
{
public:
    using Coord = Eigen::Vector3d;
    LargrangianPoint() = default;

    // Modify
    void setPosition(const Coord& position) { m_position = position; }
    void setLeftSegment(Segment* leftSegment) { m_leftsegment = leftSegment; }
    void setRightSegment(Segment* rightSegment) { m_rightsegment = rightSegment; }

    const Coord& getPosition() const { return m_position; }
    Segment* getLeftSegment() const { return m_leftsegment; }
    Segment* getRightSegment() const { return m_rightsegment; }
    // Access
protected:
    Coord m_position;
    Segment* m_leftsegment;
    Segment* m_rightsegment;
};