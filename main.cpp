#include "MultigridSolver.h"

int main(int argc, char* argv[])
{
    std::vector<double> v(17,0);
    std::vector<double> f(17,0);
    int n = 4;
    MultigridSolver Solver(n, f, v);
    Solver.PrintInfo();
