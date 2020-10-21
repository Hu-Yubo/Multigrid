#include "MultigridSolver.h"

double _f(double x)
{
    return PI*PI*sin(PI*x);
}

double _u(double x)
{
    return sin(PI*x);
}

int main(int argc, char* argv[])
{
    int n = 9;
    int a = (int)(pow(2,n));
    std::vector<double> v(a+1,0);
    std::vector<double> f(a+1,0);
    std::vector<double> RS(a+1,0);
    std::vector<double> x;
    for (double i = 0; i < a+1; i++)
    {
	x.push_back(double(i/a));
    }
    for (int i = 0; i < a+1; i++)
    {
	if (i == 0 || i == a)
	    v[i] = _u(x[i]);
	f[i] = _f(x[i]);
	RS[i] = _u(x[i]);
    }
    MultigridSolver Solver(n, f, v, "FullWeighting", "Linear", "VC");
    Solver.SetRS(RS);
    Solver.Solve();
    return 0;
}
