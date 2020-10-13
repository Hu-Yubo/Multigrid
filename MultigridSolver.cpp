/**
 * @file   MultigridSolver.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 13 21:55:00 2020
 * 
 * @brief  solve 1D possion equation by multigrid method on interval [0,1].
 * 
 * 
 */

#include "MultigridSolver.h"

MultigridSolver::MultigridSolver(int n, std::vector<double> f, std::vector<double> v, double u0, double u1, double tol, int maxstep)
{
    _n = n;
    int total_length = (int)(pow(2, n + 1)) + n - 2;
    int solution_length = (int)(pow(2, n)) + 1;
    if (f.size() != solution_length)
    {
	std::cerr << "Error!" << std::endl;
	exit(-1);
    }
    _f = std::vector<double>(total_length - solution_length, 0);
    _f.insert(_f.end(), _f.begin(), _f.end());
    if (v.size() != solution_length)
    {
	std::cerr << "Error!" << std::endl;
	exit(-1);
    }
    _v = std::vector<double>(total_length - solution_length, 0);
    _v.insert(_v.end(), _v.begin(), _v.end());
    _u0 = u0;
    _u1 = u1;
    _tol = tol;
    _maxstep = maxstep;
    _nowstep = 0;
}
