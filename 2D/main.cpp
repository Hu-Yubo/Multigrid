/**
 * @file   main.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 20 23:46:49 2020
 * 
 * @brief  calculate for test
 * 
 * 
 */

#include "MultigridSolver.h"

int CtoI(int i, int j, int SdLen)
{
    return (i*SdLen + j);
}

double _u(double x, double y)
{
    return cos(PI * x) * cos(PI * y);
}

double _f(double x, double y)
{
    return 2 * PI * PI *  cos(PI * x) * cos(PI * y);
}

int main(int argc, char* argv[])
{
    int n = 10;
    int SdLen = (int)(pow(2, n)) + 1;
    double h = 1.0 / (SdLen-1);
    int Size = SdLen * SdLen;
    std::vector<double> v(Size, 0);
    std::vector<double> f(Size, 0);
    std::vector<double> RS(Size, 0);
    for (int i = 0; i < SdLen; i++)
	for (int j = 0; j < SdLen; j++)
	{
	    if (i == 0 || i == SdLen-1 || j == 0 || j == SdLen-1)
	    {
		v[CtoI(i,j,SdLen)] = _u(i*h, j*h);
	    }
	    f[CtoI(i,j,SdLen)] = _f(i*h, j*h);
	    RS[CtoI(i,j,SdLen)] = _u(i*h, j*h);
	}
    MultigridSolver Solver(n, f, v, "FullWeighting", "Linear", "FMG");
    Solver.SetRS(RS);
    Solver.Solve();
    return 0;
}
