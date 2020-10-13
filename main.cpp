#include "MultigridSolver.h"

int main(int argc, char* argv[])
{
    std::vector<double> v(17,0);
    std::vector<double> f(17,0);
    int n = 4;
    MultigridSolver Solver(n, f, v);
    std::cout << Solver._n << std::endl;
    std::cout << Solver._f.size() << std::endl;
    std::cout << Solver._u0 << std::endl;
    std::cout << Solver._maxstep << std::endl;
    MultigridSolver Solver1();
    std::cout << Solver1._nowstep;
    return 0;
}
