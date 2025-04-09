#include "./Run.h"

int main(int argc, char** argv)
{
    Run run;
    run.initialize(argc, argv);
    run.solver();
    run.finalize();
}