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

MultigridSolver::MultigridSolver(int n, std::vector<double> f, std::vector<double> v, double u0, double u1, double tol, int maxstep, std::string S1, std::string S2)
{
    _n = n;
    int total_length = (int)(pow(2, _n + 1)) + _n - 2;
    int solution_length = (int)(pow(2, _n)) + 1;
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
    _nowlevel = 0;
    if (S1 == "FullWeighting")
    {
        _pRestrictOP = new FullWeightingRestriction();
    }
    else if (S2 == "Injection")
    {
	_pRestrictOP = new InjectionRestriction();
    }
    if (S2 == "Linear")
    {
	_pInterpolateOP = new LinearInterpolation();
    }
    _Idx = std::vector<int>(2,0);
}

void MultigridSolver::SetGridLevel(int n)
{
    _n = n;
}

void MultigridSolver::SetRightSide(std::vector<double> f)
{
    int total_length = (int)(pow(2, _n + 1)) + _n - 2;
    int solution_length = (int)(pow(2, _n)) + 1;
    if (f.size() != solution_length)
    {
	std::cerr << "Error!" << std::endl;
	exit(-1);
    }
    _f = std::vector<double>(total_length - solution_length, 0);
    _f.insert(_f.end(), _f.begin(), _f.end());
}

void MultigridSolver::SetInitialGuess(std::vector<double> v)
{
    int total_length = (int)(pow(2, _n + 1)) + _n - 2;
    int solution_length = (int)(pow(2, _n)) + 1;
    if (v.size() != solution_length)
    {
	std::cerr << "Error!" << std::endl;
	exit(-1);
    }
    _v = std::vector<double>(total_length - solution_length, 0);
    _v.insert(_v.end(), _v.begin(), _v.end());
}

void MultigridSolver::SetBoundaryCond(double u0, double u1)
{
    _u0 = u0;
    _u1 = u1;
}

void MultigridSolver::SetTolerance(double tol)
{
    _tol = tol;
}

void MultigridSolver::SetMaxStep(int maxstep)
{
    _maxstep = maxstep;
}

void MultigridSolver::SetRestrictionType(std::string S)
{
    if (S == "FullWeighting")
    {
        _pRestrictOP = new FullWeightingRestriction();
    }
    else if (S == "Injection")
    {
	_pRestrictOP = new InjectionRestriction();
    }
}

void MultigridSolver::SetNowLevel(int nowlevel)
{
    _nowlevel = nowlevel;
}

void MultigridSolver::PrintInfo()
{
    std::cout << "The total level: " << _n << std::endl;
    std::cout << "The length of f, v:" << _f.size() << " " << _v.size() << std::endl;
    std::cout << "The boundary condition: " << _u0 << " " << _u1 << std::endl;
    std::cout << "The tolerance: " << _tol << std::endl;
    std::cout << "The upper limit of iteration steps: " << _maxstep << std::endl;
    _pRestrictOP -> PrintType();
    _pInterpolateOP->PrintType();
}

RestrictionOperator* MultigridSolver::pRestrictOP()
{
    return _pRestrictOP;
}

InterpolationOperator* MultigridSolver::pInterpolateOP()
{
    return _pInterpolateOP;
}

void MultigridSolver::UpdateIndex()
{
    _Idx[0] = (int)(pow(2, _n - _nowlevel + 1)) - 2 + _n - _nowlevel;
    _Idx[1] = _Idx[0] + (int)(pow(2, _n - _nowlevel + 1)); 
}
void MultigridSolver::PrintIdx()
{
    std::cout << _Idx[0] << " " << _Idx[1] << std::endl;
}

void MultigridSolver::WeightedJacobi()
{
    std::vector<double> v_star(_Idx[1] - _Idx[0] + 1);
    
}
