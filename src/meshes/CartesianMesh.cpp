#include "CartesianMesh.h"
#include <iostream>

void CartesianMesh::allocateCells(std::vector<std::unique_ptr<Cell>> &cells)
{
    cells.reserve(m_numberCellsX * m_numberCellsY * m_numberCellsZ);
    for (int i = 0; i < m_numberCellsX; ++i)
    {
        for (int j = 0; j < m_numberCellsY; ++j)
        {
            for(int k = 0; k < m_numberCellsZ; ++k)
            {
                auto cell = std::make_unique<Cell>();
                cell->setindex(Cell::Coord(i, j, k));
                cells.push_back(std::move(cell));
            }
        }
    }
}

void CartesianMesh::setelement()
{
    m_UnidX = m_GlobalLengthX / static_cast<double>(m_numberCellsXGlobal);
    m_UnidY = m_GlobalLengthY / static_cast<double>(m_numberCellsYGlobal);
    m_UnidZ = m_GlobalLengthZ / static_cast<double>(m_numberCellsZGlobal);

    m_elements.reserve(m_numberCellsX * m_numberCellsY * m_numberCellsZ);
    for (int i = 0; i < m_numberCellsX; ++i)
    {
        for (int j = 0; j < m_numberCellsY; ++j)
        {
            for(int k = 0; k < m_numberCellsZ; ++k)
            {
                auto element = std::make_unique<Element>();
                element->setL1(m_UnidX);
                element->setL2(m_UnidY);
                element->setL3(m_UnidZ);
                element->setvolume(m_UnidX * m_UnidY * m_UnidZ);
                element->setposition(Eigen::Vector3d((i + m_offsetX + 0.5) * m_UnidX, (j + m_offsetY + 0.5) * m_UnidY, (k + 0.5) * m_UnidZ));
                m_elements.push_back(std::move(element));
            }
        }
    }

}

void CartesianMesh::BindCellwithElement(std::vector<std::unique_ptr<Cell>> &cells)
{
    for(int i = 0; i < m_numberCellsX; ++i)
    {
        for (int j = 0; j < m_numberCellsY; ++j)
        {
            for(int k = 0; k < m_numberCellsZ; ++k)
            {
                auto &cell = cells[k + j * m_numberCellsZ + i * m_numberCellsY * m_numberCellsZ];
                auto &element = m_elements[k + j * m_numberCellsZ + i * m_numberCellsY * m_numberCellsZ];
                cell->setelement(element.get());
                cell->setposition(element->getposition());
            }
        }
    }
}

void CartesianMesh::BindElementwithFace()
{
    for (int i = 0; i < m_numberCellsX; ++i)
    {
        for (int j = 0; j < m_numberCellsY; ++j)
        {
            for(int k = 0; k < m_numberCellsZ; ++k)
            {
                auto &element = m_elements[k + j * m_numberCellsZ + i * m_numberCellsY * m_numberCellsZ];
                auto leftface = std::make_unique<Face>();
                leftface->setposition(element->getposition() - Eigen::Vector3d(0.5 * element->getL1(), 0.0, 0.0));
                leftface->setnormal(Eigen::Vector3d(-1.0, 0.0, 0.0));
                auto rightface = std::make_unique<Face>();
                rightface->setposition(element->getposition() + Eigen::Vector3d(0.5 * element->getL1(), 0.0, 0.0));
                rightface->setnormal(Eigen::Vector3d(1.0, 0.0, 0.0));
                auto bottomface = std::make_unique<Face>();
                bottomface->setposition(element->getposition() - Eigen::Vector3d(0.0, 0.5 * element->getL2(), 0.0));
                bottomface->setnormal(Eigen::Vector3d(0.0, -1.0, 0.0));
                auto topface = std::make_unique<Face>();
                topface->setposition(element->getposition() + Eigen::Vector3d(0.0, 0.5 * element->getL2(), 0.0));
                topface->setnormal(Eigen::Vector3d(0.0, 1.0, 0.0));
                auto frontface = std::make_unique<Face>();
                frontface->setposition(element->getposition() - Eigen::Vector3d(0.0, 0.0, 0.5 * element->getL3()));
                frontface->setnormal(Eigen::Vector3d(0.0, 0.0, -1.0));
                auto backface = std::make_unique<Face>();
                backface->setposition(element->getposition() + Eigen::Vector3d(0.0, 0.0, 0.5 * element->getL3()));
                backface->setnormal(Eigen::Vector3d(0.0, 0.0, 1.0));
                element->setface(0, std::move(leftface));
                element->setface(1, std::move(rightface));
                element->setface(2, std::move(bottomface));
                element->setface(3, std::move(topface));
                element->setface(4, std::move(frontface));
                element->setface(5, std::move(backface));
            }

        }
    }
}

void CartesianMesh::cutcell(Geom *geom)
{   
    int cutcell_num {};
    for(auto& element:m_elements){
        geom->SortSegment2Element(element.get());
        if(element->getsegments().size()){
            cutcell_num++;
        }
    }

    std::cout <<"cut cell num is : "<<cutcell_num<<std::endl;
}

void CartesianMesh::setnumberCellsX(const int &N)
{
    m_numberCellsX = N;
}

void CartesianMesh::setnumberCellsY(const int &N)
{
    m_numberCellsY = N;
}
void CartesianMesh::setnumberCellsZ(const int &N)
{
    m_numberCellsZ = N;
}
void CartesianMesh::setnumberCellsXGlobal(const int &N)
{
    m_numberCellsXGlobal = N;
}
void CartesianMesh::setnumberCellsYGlobal(const int &N)
{
    m_numberCellsYGlobal = N;
}
void CartesianMesh::setnumberCellsZGlobal(const int &N)
{
    m_numberCellsZGlobal = N;
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
void CartesianMesh::setGlobalLengthZ(const double &L)
{
    m_GlobalLengthZ = L;
}
void CartesianMesh::setLocalLengthX(const double &L)
{
    m_LocalLengthX = L;
}
void CartesianMesh::setLocalLengthY(const double &L)
{
    m_LocalLengthY = L;
}
void CartesianMesh::setLocalLengthZ(const double &L)
{
    m_LocalLengthZ = L;
}
const int &CartesianMesh::getnumberCellsXGlobal()
{
    return m_numberCellsXGlobal;
}
const int &CartesianMesh::getnumberCellsYGlobal()
{
    return m_numberCellsYGlobal;
}
const int &CartesianMesh::getnumberCellsZGlobal()
{
    return m_numberCellsZGlobal;
}
const int &CartesianMesh::getnumberCellsX()
{
    return m_numberCellsX;
}
const int &CartesianMesh::getnumberCellsY()
{
    return m_numberCellsY;
}
const int &CartesianMesh::getnumberCellsZ()
{
    return m_numberCellsZ;
}
const double &CartesianMesh::getUnidX()
{
    return m_UnidX;
}
const double &CartesianMesh::getUnidY()
{
    return m_UnidY;
}
const double &CartesianMesh::getUnidZ()
{
    return m_UnidZ;
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
const double &CartesianMesh::getGlobalLengthZ()
{
    return m_GlobalLengthZ;
}
const double &CartesianMesh::getLocalLengthX()
{
    return m_LocalLengthX;
}
const double &CartesianMesh::getLocalLengthY()
{
    return m_LocalLengthY;
}
const double &CartesianMesh::getLocalLengthZ()
{
    return m_LocalLengthZ;
}