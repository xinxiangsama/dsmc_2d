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
const std::vector<std::unique_ptr<Segment>> &Element::getsegments() const
{
    return m_segments;
}
const std::vector<std::unique_ptr<LargrangianPoint>> &Element::getIntersectionPs() const
{
    return m_intersectionPs;
}
std::array<std::unique_ptr<Vertice>, 4> &Element::getvertices()
{
    return m_vertices;
}
bool Element::ifcut()
{
    return !((!m_vertices[0]->isWall()) && (!m_vertices[1]->isWall()) && (!m_vertices[2]->isWall()) && (!m_vertices[3]->isWall()));
}
void Element::setposition(const Coord &position)
{
    m_position = position;
    m_vertices[0] = std::make_unique<Vertice>(Eigen::Vector3d{m_position(0) - 0.5 * m_L1, m_position(1) - 0.5 * m_L2, 0.0});
    m_vertices[1] = std::make_unique<Vertice>(Eigen::Vector3d{m_position(0) + 0.5 * m_L1, m_position(1) - 0.5 * m_L2, 0.0});
    m_vertices[2] = std::make_unique<Vertice>(Eigen::Vector3d{m_position(0) + 0.5 * m_L1, m_position(1) + 0.5 * m_L2, 0.0});
    m_vertices[3] = std::make_unique<Vertice>(Eigen::Vector3d{m_position(0) - 0.5 * m_L1, m_position(1) + 0.5 * m_L2, 0.0});
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

void Element::insertsegment(std::unique_ptr<Segment>& segment)
{
    m_segments.emplace_back(std::move(segment));
}

void Element::insertIntersectionP(std::unique_ptr<Eigen::Vector2d>& P)
{
    m_intersectionPs.emplace_back(std::make_unique<LargrangianPoint>(LargrangianPoint::Coord{P->x(), P->y(), 0.0}));
}

bool Element::ifContain2d(const Eigen::Vector2d &P)
{
    auto xmin = m_position(0) - 0.5 * m_L1;
    auto xmax = m_position(0) + 0.5 * m_L1;
    auto ymin = m_position(1) - 0.5 * m_L2;
    auto ymax = m_position(1) + 0.5 * m_L2;

    if(P.x() == xmin && P.y() > ymin && P.y() < ymax){
        return true;
    }else if (P.x() == xmax && P.y() > ymin && P.y() < ymax){
        return true;
    }else if (P.y() == ymin && P.x() > xmin && P.x() < xmax){
        return true;
    }else if (P.y() == ymax && P.x() > xmin && P.x() < xmax){
        return true;
    }
    return false;
}
