#include "./Run.h"
#include <vtkVersion.h>

int main(int argc, char** argv)
{
    Run run;
    run.initialize(argc, argv);
    run.solver();
    // std::cout << "VTK Version: "
    // << VTK_MAJOR_VERSION << "."
    // << VTK_MINOR_VERSION << std::endl;
    run.finalize();
}