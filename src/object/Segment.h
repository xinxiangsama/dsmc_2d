#pragma once
#include "./Largrangian.h"
class LargrangianPoint;
class Segment
{
public:
using Coord = LargrangianPoint::Coord;
    Segment() = default;

    // Access
    const double& getlength() const;
    const Coord& getsloop() const;
    const Coord& getnormal() const;
    LargrangianPoint*  getleftpoint() const;
    LargrangianPoint*  getrightpoint() const;
    // Modify
    void setlength(const double& length);
    void setsloop(const Coord& sloop);
    void setnormal(const Coord& normal);
    void setleftpoint(LargrangianPoint* leftpoint);
    void setrightpoint(LargrangianPoint* rightpoint);

protected:
    double m_length;
    LargrangianPoint* m_leftpoint;
    LargrangianPoint* m_rightpoint;
    Coord m_sloop;
    Coord m_normal;
};