#include "./Run.h"

int main(int argc, char** argv)
{
    Run run;
    run.initialize(argc, argv);
    run.finalize();
    return 0;
}