#include "MultigridSolver.h"

int main(int argc, char* argv[])
{
    std::vector<double> v(17,0);
    std::vector<double> f(17,0);
    int n = 4;
    MultigridSolver Solver(n, f, v);
    Solver.PrintInfo();
    /*
    double a[9] = {0, 0.7071, -1.0000, 0.7071, 0.0000,-0.7071, 1.0000, -0.7071, -0.0000};
    std::vector<double> A(a,a+9);
    Solver.pRestrictOP()->SetInput(A);
    Solver.pRestrictOP()->restrict();
    Solver.pRestrictOP()->PrintType();
    Solver.pRestrictOP()->PrintInput();
    Solver.pRestrictOP()->PrintOutput();
    Solver.pInterpolateOP()->SetInput(A);
    Solver.pInterpolateOP()->interpolate();
    Solver.pInterpolateOP()->PrintType();
    Solver.pInterpolateOP()->PrintInput();
    Solver.pInterpolateOP()->PrintOutput();
    */
    Solver.PrintIdx();
    Solver.SetNowLevel(1);
    Solver.UpdateIndex();
    Solver.PrintIdx();
    Solver.SetNowLevel(2);
    Solver.UpdateIndex();
    Solver.PrintIdx();
    Solver.SetNowLevel(3);
    Solver.UpdateIndex();
    Solver.PrintIdx();
    Solver.SetNowLevel(4);
    Solver.UpdateIndex();
    Solver.PrintIdx();
}
