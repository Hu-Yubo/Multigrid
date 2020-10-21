/**
 * @file   main.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 20 23:46:49 2020
 * 
 * @brief  
 * 
 * 
 */

#include "MultigridSolver.h"

int CtoI(int i, int j, int SdLen)
{
    return (i*SdLen + j);
}

int main(int argc, char* argv[])
{
    int n =2;
    int SdLen = (int)(pow(2, n)) + 1;
    double h = 1.0 / (SdLen-1);
    int Size = SdLen * SdLen;
    std::vector<double> v(Size, 0);
    std::vector<double> f(Size, 0);
    for (int i = 0; i < SdLen; i++)
	for (int j = 0; j < SdLen; j++)
	    f[CtoI(i,j,SdLen)] = 2 * PI * PI * sin(PI * i * h) * sin(PI * j * h);
    MultigridSolver Solver(n, f, v);
    Solver.Solve();
    return 0;
}
