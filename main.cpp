#include "MultigridSolver.h"

int main(int argc, char* argv[])
{
    int n = 10;
    int a = (int)(pow(2,n));
    std::vector<double> v(a+1,0);
    std::vector<double> f(a+1,0);
    std::vector<double> x;
    for (double i = 0; i < a+1; i++)
    {
	x.push_back(double(i/a));
    }
    for (int i = 0; i < a+1; i++)
	/// f[i] = PI * PI * sin(PI * i / a);
	f[i] = (sin(x[i])-cos(x[i])*cos(x[i]))*exp(sin(x[i]));
    /// MultigridSolver Solver(n, f, v);
    MultigridSolver Solver(n, f, v, 1, exp(sin(1.0)));
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
    /*
    for (int i = 0; i < AS.size();i++)
    {
	std::cout << AS[i] << " ";
    }
    /*
    AS = Solver.ReturnSolution();
    double e = 0;
    for (int i = 0; i < a; i++)
	e = e + pow((AS[i] - sin(PI * i / a)), 2);
    std::cout << sqrt(e) << std::endl;
    */
    return 0;
}
