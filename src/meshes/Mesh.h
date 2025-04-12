#pragma once
#include "../cell/Cell.h"
#include "object/Geom.h"

class Mesh
{
public:
    Mesh() = default;
    virtual ~Mesh() = default;

    virtual void allocateCells(std::vector<std::unique_ptr<Cell>>& cells) {};
    virtual void setelement() {};
    virtual void BindCellwithElement(std::vector<std::unique_ptr<Cell>>& cells) {};
    virtual void BindElementwithFace() {};
    virtual void cutcell(Geom* geom) {};

    // Modify
    virtual void setnumberCellsX(const int&) {};
    virtual void setnumberCellsY(const int&) {};
    virtual void setnumberCellsZ(const int&) {};
    virtual void setnumberCellsXGlobal(const int& N) {};
    virtual void setnumberCellsYGlobal(const int& N) {};
    virtual void setnumberCellsZGlobal(const int& N) {};
    virtual void setnumberCpuX(const int&) {};
    virtual void setnumberCpuY(const int&) {};
    virtual void setCpuCoordX(const int&) {};
    virtual void setCpuCoordY(const int&) {};
    virtual void setoffsetX(const int&) {};
    virtual void setoffsetY(const int&) {};
    virtual void setGlobalLengthX(const double&) {};
    virtual void setGlobalLengthY(const double&) {};
    virtual void setGlobalLengthZ(const double&) {};
    virtual void setLocalLengthX(const double& ) {};
    virtual void setLocalLengthY(const double& ) {};
    virtual void setLocalLengthZ(const double& ) {};
    // Access
    virtual const int& getnumberCellsXGlobal() {};
    virtual const int& getnumberCellsYGlobal() {};
    virtual const int& getnumberCellsZGlobal() {};
    virtual const int& getnumberCellsX() {};
    virtual const int& getnumberCellsY() {};
    virtual const int& getnumberCellsZ() {};
    virtual const double& getUnidX() {};
    virtual const double& getUnidY() {};
    virtual const double& getUnidZ() {};
    virtual const int& getnumberCpuX() {};
    virtual const int& getnumberCpuY() {};
    virtual const int& getCpuCoordX() {};
    virtual const int& getCpuCoordY() {};
    virtual const int& getoffsetX() {};
    virtual const int& getoffsetY() {};
    virtual const double& getGlobalLengthX() {};
    virtual const double& getGlobalLengthY() {};
    virtual const double& getGlobalLengthZ() {};
    virtual const double& getLocalLengthX() {};
    virtual const double& getLocalLengthY() {};
    virtual const double& getLocalLengthZ() {};
};