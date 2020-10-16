#include "MultigridSolver.h"

#define PI (4.0 * atan(1.0))

int main(int argc, char* argv[])
{
    int n = 3;
    int a = (int)(pow(2,n));
    std::vector<double> v(a+1,0);
    std::vector<double> f(a+1,0);
    for (int i = 0; i < a+1; i++)
	f[i] = PI * PI * sin(PI * i / a);
    MultigridSolver Solver(n, f, v);
    /*
    Solver.PrintInfo();
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
    std::vector<double> AS;
    /*
    for (int j = 0; j < 5; j++)
    {
	Solver.VCycle();
	AS = Solver.ReturnSolution();
	double e = 0;
	for (int i = 0; i < a; i++)
	    e = e + pow((AS[i] - sin(PI * i / a)), 2);
	std::cout << sqrt(e) << std::endl;
	}*/
    Solver.Solve();
    AS = Solver.ReturnSolution();
    double e = 0;
    for (int i = 0; i < a; i++)
	e = e + pow((AS[i] - sin(PI * i / a)), 2);
    std::cout << sqrt(e) << std::endl;
    return 0;
}
