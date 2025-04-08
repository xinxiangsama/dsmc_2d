#pragma once
#include "Mesh.h"

class CartesianMesh : public Mesh
{
public:
    CartesianMesh() = default;
    ~CartesianMesh() = default;
    void allocateCells(std::vector<std::unique_ptr<Cell>>& cells) override;
    void setelement() override;
    void setface() override;
    void BindCellwithElement(std::vector<std::unique_ptr<Cell>>& cells);
    void BindElementwithFace();
    // Modifiers
    void setnumberCellsX(const int& N) override;
    void setnumberCellsY(const int& N) override;
    void setnumberCellsXGlobal(const int& N) override;
    void setnumberCellsYGlobal(const int& N) override;
    void setnumberCpuX(const int& N) override;
    void setnumberCpuY(const int& N) override;
    void setCpuCoordX(const int& N) override;
    void setCpuCoordY(const int& N) override;
    void setGlobalLengthX(const double& L) override;
    void setGlobalLengthY(const double& L) override;
    void setLocalLengthX(const double& L) override;
    void setLocalLengthY(const double& L) override;
    // Accessers
    const int& getnumberCellsXGlobal() override;
    const int& getnumberCellsYGlobal() override;
    const int& getnumberCellsX() override;
    const int& getnumberCellsY() override;
    const double& getUnidX() override;
    const double& getUnidY() override;
    const int& getoffsetX() override;
    const int& getoffsetY() override;
    const double& getGlobalLengthX() override;
    const double& getGlobalLengthY() override;
    const double& getLocalLengthX() override;
    const double& getLocalLengthY() override;
protected:
    int m_numberCellsXGlobal = 0; // number of cells in the global domain
    int m_numberCellsYGlobal = 0; // number of cells in the global domain
    int m_numberCellsX = 0; // number of cells in the local domain
    int m_numberCellsY = 0; // number of cells in the local domain
    double m_GlobalLengthX = 0.0; // length of the global domain in x direction
    double m_GlobalLengthY = 0.0; // length of the global domain in y direction
    double m_LocalLengthX = 0.0; // length of the local domain in x direction
    double m_LocalLengthY = 0.0; // length of the local domain in y direction
    int m_offsetX = 0; // offset in x direction for the local domain
    int m_offsetY = 0; // offset in y direction for the local domain
    int m_numberCpuX = 1; // number of processors in x direction
    int m_numberCpuY = 1; // number of processors in y direction
    int m_CpuCoordX = 0; // coordinate of the processor in x direction
    int m_CpuCoordY = 0; // coordinate of the processor in y direction
    double m_UnidX = 0.0; // length of the cell in x direction
    double m_UnidY = 0.0; // length of the cell in y direction
    std::vector<double> m_dXi; // vector of cell lengths in x direction
    std::vector<double> m_dYj; // vector of cell lengths in y direction

    std::vector<std::unique_ptr<Element>> m_elements; // vector of elements
    std::vector<std::unique_ptr<Face>> m_faces; // vector of faces
};