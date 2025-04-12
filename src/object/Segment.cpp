#include "Segment.h"

const double &Segment::getlength() const
{
    return m_length;
}

const LargrangianPoint::Coord &Segment::getsloop() const
{
    return m_sloop;
}

const LargrangianPoint::Coord &Segment::getnormal() const
{
    return m_normal;
}

LargrangianPoint *Segment::getleftpoint() const
{
    return m_leftpoint;
}

LargrangianPoint *Segment::getrightpoint() const
{
    return m_rightpoint;
}

void Segment::setlength(const double &length)
{
    m_length = length;
}

void Segment::setsloop(const Coord &sloop)
{
    m_sloop = sloop;
}

void Segment::setnormal(const Coord &normal)
{
    m_normal = normal;
}

void Segment::setleftpoint(LargrangianPoint *leftpoint)
{
    m_leftpoint = leftpoint;
}

void Segment::setrightpoint(LargrangianPoint *rightpoint)
{
    m_rightpoint = rightpoint;
}