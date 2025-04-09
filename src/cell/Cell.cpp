#include "Cell.h"

const Cell::Coord &Cell::getposition() const
{
    return m_position;
}
const std::vector<Particle *> &Cell::getparticles() const
{
    return m_particles;
}
const Phase *Cell::getphase() const
{
    return m_phase.get();
}
const Element *Cell::getelement() const
{
    return m_element;
}
const Cell::Coord &Cell::getindex() const
{
    return m_index;
}
const std::vector<std::unique_ptr<Cell>> &Cell::getchildren() const
{
    return m_children;
}
void Cell::setposition(const Coord &position)
{
    m_position = position;
}
void Cell::setelement(Element* element)
{
    m_element = element;
}
void Cell::setindex(const Coord &index)
{
    m_index = index;
}
void Cell::allocatevar()
{
    m_phase = std::make_unique<Phase>();
}
void Cell::insertparticle(Particle *particle)
{
    m_particles.push_back(particle);
}
void Cell::removeparticle(Particle *particle)
{
    auto it = std::remove(m_particles.begin(), m_particles.end(), particle);
    m_particles.erase(it);
}

void Cell::removeallparticles()
{
    m_particles.clear();
}
void Cell::collision()
{   
    // until the smallest cell is reached
    if(m_children.size() != 0)
    {
        for(auto &child : m_children)
        {
            child->collision();
        }
    }else{
        // collision is this cell
    }
}
void Cell::sample()
{
}
Cell::~Cell()
{
    removeallparticles();
}