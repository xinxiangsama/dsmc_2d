#include "Element.h"
#include <iostream>
#include <stdexcept>

const Element::Coord &Element::getposition() const
{
    return m_position;
}

const double &Element::getL1() const
{
    return m_L1;
}

const double &Element::getL2() const
{
    return m_L2;
}
const double &Element::getL3() const
{
    return m_L3;
}
const double &Element::getvolume() const
{
    return m_volume;
}
const std::array<std::unique_ptr<Face>, 6> &Element::getfaces() const
{
    return m_faces;
}
Face *Element::getface(size_t index) const
{
    if (index < m_faces.size())
    {
        return m_faces[index].get();
    }
    std::cerr << "Index out of bounds in getface: " << index << std::endl;
    return nullptr; // or throw an exception
}
const std::vector<Segment *> &Element::getsegments() const
{
    return m_segments;
}
bool Element::ifcut()
{
    return m_segments.size() != 0; // if not zero, means it has been cut;
}
void Element::setposition(const Coord &position)
{
    m_position = position;
}
void Element::setL1(const double &L1)
{
    m_L1 = L1;
}
void Element::setL2(const double &L2)
{
    m_L2 = L2;
}
void Element::setL3(const double &L3)
{
    m_L3 = L3;
}
void Element::setvolume(const double &volume)
{
    m_volume = volume;
}
void Element::setface(size_t index, std::unique_ptr<Face>&& face)
{
    if (index < m_faces.size())
    {
        m_faces[index] = std::move(face);
    }else
    {
        std::cerr << "Index out of bounds in setface: " << index << std::endl;
    }
}

void Element::setface(size_t index, std::unique_ptr<Face> &face)
{
    if (index < m_faces.size())
    {
        m_faces[index] = std::move(face);
    }else
    {
        std::cerr << "Index out of bounds in setface: " << index << std::endl;
    }
}

void Element::insertsegment(Segment *segment)
{
    m_segments.push_back(segment);
}
