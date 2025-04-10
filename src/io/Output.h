#pragma once
#include <H5Cpp.h>
#include <string>
#include <ostream>
#include "../Run.h"

class Run;
class Output
{
public:
    Output() = default;
    Output(Run* run);

    void Write2HDF5(const std::string& filename);
protected:
    Run* m_run {nullptr};
};