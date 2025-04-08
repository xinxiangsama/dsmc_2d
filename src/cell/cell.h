#pragma once
#include "../particle/Particle.h"
#include "../meshes/Element.h"
#include <memory>
#include <vector>

class Cell
{
public:
    using Coord = Particle::Coord;
    Cell() = default;
    virtual ~Cell() = default;
protected:
    Coord m_position;
    std::vector<Particle*> m_particles;
    std::unique_ptr<Element> m_element;
    std::vector<std::unique_ptr<Cell>> m_children;
};