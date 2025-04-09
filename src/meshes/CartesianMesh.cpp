#include "CartesianMesh.h"
#include <iostream>

void CartesianMesh::allocateCells(std::vector<std::unique_ptr<Cell>> &cells)
{
    cells.reserve(m_numberCellsX * m_numberCellsY);
    for (int i = 0; i < m_numberCellsX; ++i)
    {
        for (int j = 0; j < m_numberCellsY; ++j)
        {
            cells.emplace_back(std::make_unique<Cell>());
        }
    }
}

void CartesianMesh::setelement()
{
    m_UnidX = m_GlobalLengthX / static_cast<double>(m_numberCellsXGlobal);
    m_UnidY = m_GlobalLengthY / static_cast<double>(m_numberCellsYGlobal);

    m_elements.reserve(m_numberCellsX * m_numberCellsY);
    for (int i = 0; i < m_numberCellsX; ++i)
    {
        for (int j = 0; j < m_numberCellsY; ++j)
        {
            auto element = std::make_unique<Element>();
            element->setL1(m_UnidX);
            element->setL2(m_UnidY);
            element->setposition(Eigen::Vector2d((i + m_offsetX + 0.5) * m_UnidX, (j + m_offsetY + 0.5) * m_UnidY));
            element->setvolume(m_UnidX * m_UnidY);
            m_elements.push_back(std::move(element));
        }
    }

}

void CartesianMesh::BindCellwithElement(std::vector<std::unique_ptr<Cell>> &cells)
{
    for(int i = 0; i < m_numberCellsX; ++i)
    {
        for (int j = 0; j < m_numberCellsY; ++j)
        {
            auto &cell = cells[i + j * m_numberCellsX];
            auto &element = m_elements[i + j * m_numberCellsX];
            cell->setposition(element->getposition());
            cell->setelement(element.get());
        }
    }
}

void CartesianMesh::BindElementwithFace()
{
    for (int i = 0; i < m_numberCellsX; ++i)
    {
        for (int j = 0; j < m_numberCellsY; ++j)
        {
            auto &element = m_elements[i + j * m_numberCellsX];
            auto leftface = std::make_unique<Face>();
            auto rightface = std::make_unique<Face>();
            auto bottomface = std::make_unique<Face>();
            auto topface = std::make_unique<Face>();
            leftface->setposition(element->getposition() - Eigen::Vector2d(0.5 * element->getL1(), 0));
            leftface->setnormal(Eigen::Vector2d(-1, 0));
            rightface->setposition(element->getposition() + Eigen::Vector2d(0.5 * element->getL1(), 0));
            rightface->setnormal(Eigen::Vector2d(1, 0));
            bottomface->setposition(element->getposition() - Eigen::Vector2d(0, 0.5 * element->getL2()));
            bottomface->setnormal(Eigen::Vector2d(0, -1));
            topface->setposition(element->getposition() + Eigen::Vector2d(0, 0.5 * element->getL2()));
            topface->setnormal(Eigen::Vector2d(0, 1));
            element->setface(0, std::move(leftface));
            element->setface(1, std::move(rightface));
            element->setface(2, std::move(bottomface));
            element->setface(3, std::move(topface));

        }
    }
}

void CartesianMesh::setnumberCellsX(const int &N)
{
    m_numberCellsX = N;
}

void CartesianMesh::setnumberCellsY(const int &N)
{
    m_numberCellsY = N;
}
void CartesianMesh::setnumberCellsXGlobal(const int &N)
{
    m_numberCellsXGlobal = N;
}
void CartesianMesh::setnumberCellsYGlobal(const int &N)
{
    m_numberCellsYGlobal = N;
}
void CartesianMesh::setnumberCpuX(const int &N)
{
    m_numberCpuX = N;
}
void CartesianMesh::setnumberCpuY(const int &N)
{
    m_numberCpuY = N;
}
void CartesianMesh::setCpuCoordX(const int &N)
{
    m_CpuCoordX = N;
}
void CartesianMesh::setCpuCoordY(const int &N)
{
    m_CpuCoordY = N;
}
void CartesianMesh::setoffsetX(const int &N)
{
    m_offsetX = N;
}
void CartesianMesh::setoffsetY(const int &N)
{
    m_offsetY = N;
}
void CartesianMesh::setGlobalLengthX(const double &L)
{
    m_GlobalLengthX = L;
}
void CartesianMesh::setGlobalLengthY(const double &L)
{
    m_GlobalLengthY = L;
}
void CartesianMesh::setLocalLengthX(const double &L)
{
    m_LocalLengthX = L;
}
void CartesianMesh::setLocalLengthY(const double &L)
{
    m_LocalLengthY = L;
}
const int &CartesianMesh::getnumberCellsXGlobal()
{
    return m_numberCellsXGlobal;
}
const int &CartesianMesh::getnumberCellsYGlobal()
{
    return m_numberCellsYGlobal;
}
const int &CartesianMesh::getnumberCellsX()
{
    return m_numberCellsX;
}
const int &CartesianMesh::getnumberCellsY()
{
    return m_numberCellsY;
}
const double &CartesianMesh::getUnidX()
{
    return m_UnidX;
}
const double &CartesianMesh::getUnidY()
{
    return m_UnidY;
}
const int &CartesianMesh::getnumberCpuX()
{
    return m_numberCpuX;
}
const int &CartesianMesh::getnumberCpuY()
{
    return m_numberCpuY;
}
const int &CartesianMesh::getCpuCoordX()
{
    return m_CpuCoordX;
}
const int &CartesianMesh::getCpuCoordY()
{
    return m_CpuCoordY;
}
const int &CartesianMesh::getoffsetX()
{
    return m_offsetX;
}
const int &CartesianMesh::getoffsetY()
{
    return m_offsetY;
}
const double &CartesianMesh::getGlobalLengthX()
{
    return m_GlobalLengthX;
}
const double &CartesianMesh::getGlobalLengthY()
{
    return m_GlobalLengthY;
}
const double &CartesianMesh::getLocalLengthX()
{
    return m_LocalLengthX;
}
const double &CartesianMesh::getLocalLengthY()
{
    return m_LocalLengthY;
}