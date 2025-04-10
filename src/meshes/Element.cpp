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

const double &Element::getvolume() const
{
    return m_volume;
}
const std::array<std::unique_ptr<Face>, 4> &Element::getfaces() const
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
