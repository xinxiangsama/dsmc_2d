#pragma once
#include "../particle/Particle.h"
#include "../meshes/Element.h"
#include "../variables/Phase.h"
#include <memory>
#include <vector>

class Cell
{
public:
    using Coord = Particle::Coord;
    Cell() = default;
    virtual ~Cell();

    // Accessors
    const Coord& getposition() const;
    const std::vector<Particle*>& getparticles() const;
    const Phase* getphase() const;
    const Element* getelement() const;
    const std::vector<std::unique_ptr<Cell>>& getchildren() const;
    // Modifiers
    void setposition(const Coord& position);
    void setelement(std::unique_ptr<Element> element);

    // Functions
    void allocatevar();
    void insertparticle(Particle* particle);
    void removeparticle(Particle* particle);
    void removeallparticles();

    virtual void collision();
protected:
    Coord m_position;
    std::vector<Particle*> m_particles;
    std::unique_ptr<Phase> m_phase;
    std::unique_ptr<Element> m_element;
    std::vector<std::unique_ptr<Cell>> m_children;
};