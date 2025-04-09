#pragma once
#include "Parallel.h"
#include <array>

class CartesianParallel : public Parallel
{
public:
    CartesianParallel() = default;
    ~CartesianParallel() = default;

    void ZoneDcomposition() override;
    void setNeibours() override;

    const int& getLeftNeibour() override;
    const int& getRightNeibour() override;
    const int& getTopNeibour() override;
    const int& getBottomNeibour() override;
    void info() override;
protected:
    std::array<int, 4> m_neighbours; // left, right, bottom, top
    MPI_Comm m_cartesian_comm;
};