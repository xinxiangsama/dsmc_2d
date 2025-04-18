#pragma once
#include "../particle/Particle.h"
#include "../meshes/Element.h"
#include "../variables/Phase.h"
#include "../random/Random.h"
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
    const std::vector<Particle* >& getparticles() const;
    const Phase* getphase() const;
    const Element* getelement() const;
    const Coord& getindex() const;
    const std::vector<std::shared_ptr<Cell>>& getchildren() const;
    const int getCollisionNum();
    bool ifcut();
    // Modifiers
    void setposition(const Coord& position);
    void setelement(Element* element);
    void setindex(const Coord& index);

    // Functions
    void allocatevar();
    void insertparticle(Particle* particle);
    void removeparticle(Particle*  particle);
    void removeallparticles();

    // collision
    virtual void collision();
    // sample micro physics
    virtual void sample();
protected:
    Coord m_position;
    Coord m_index;
    std::vector<Particle*> m_particles;
    std::shared_ptr<Phase> m_phase;
    Element* m_element;
    std::vector<std::shared_ptr<Cell>> m_children;
    size_t N_particles {};
    size_t N_collision {};
    double Mcand_r {}; // the resident caused buy turn double to int
};