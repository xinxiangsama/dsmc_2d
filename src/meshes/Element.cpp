#include "Element.h"
#include <iostream>
#include <stdexcept>

const Element::Coord &Element::getposition() const
{
    return m_position;
}

const double &Element::getvolume() const
{
    return m_volume;
}
const std::array<Face *, 4> &Element::getfaces() const
{
    return m_faces;
}
Face *Element::getface(size_t index) const
{
    if (index < m_faces.size())
    {
        return m_faces[index];
    }
    std::cerr << "Index out of bounds in getface: " << index << std::endl;
    return nullptr; // or throw an exception
}
void Element::setposition(const Coord &position)
{
    m_position = position;
}
void Element::setvolume(const double &volume)
{
    m_volume = volume;
}
void Element::setface(size_t index, Face *face)
{
    if (index < m_faces.size())
    {
        m_faces[index] = face;
    }
}
