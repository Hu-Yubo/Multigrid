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
    _h = 1.0/pow(2, _n);
    int total_length = (int)(pow(2, _n + 1)) + _n - 2;
    int solution_length = (int)(pow(2, _n)) + 1;
    if (f.size() != solution_length)
    {
	std::cerr << "Error!" << std::endl;
	exit(-1);
    }
    _f = std::vector<double>(total_length - solution_length, 0);
    _f.insert(_f.end(), f.begin(), f.end());
    if (v.size() != solution_length)
    {
	std::cerr << "Error!" << std::endl;
	exit(-1);
    }
    _v = std::vector<double>(total_length - solution_length, 0);
    _v.insert(_v.end(), v.begin(), v.end());
    _u0 = u0;
    _u1 = u1;
    _tol = tol;
    _maxstep = maxstep;
    _nowlevel = 1;
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
    _w = 2.0/3;
}

void MultigridSolver::SetGridLevel(int n)
{
    _n = n;
    _h = 1/pow(2, n);
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
    _h = 1 / pow(2, _n - _nowlevel + 1);
}
void MultigridSolver::PrintIdx()
{
    std::cout << _Idx[0] << " " << _Idx[1] << std::endl;
}

void MultigridSolver::WeightedJacobi()
{
    std::vector<double> v_star(_Idx[1] - _Idx[0] + 1, 0);
    for (int i = 1; i < v_star.size() - 1; i++)
	v_star[i] = (_v[i-1+_Idx[0]] + _v[i+1+_Idx[0]] + _h*_h*_f[i+_Idx[0]]) / 2;
    for (int i = 1; i < v_star.size() - 1; i++)
	_v[i+_Idx[0]] = _v[i+_Idx[0]] + _w * (v_star[i] - _v[i+_Idx[0]]);
}

void MultigridSolver::VCycle()
{
    UpdateIndex();
    if (_nowlevel == _n)
	BottomSolve();
    else
    {
	for (int i = 0; i < _RlxTimes; i++)
	    WeightedJacobi();
	std::vector<double> r_h;
	std::vector<double> r_2h;
	for (int i = _Idx[0]; i <= _Idx[1]; i++)
	{
	    if (i == _Idx[0])
		///r_h.push_back(_f[i] - (2*_v[i] - _v[i+1])/_h/_h);
		r_h.push_back(0);
	    else if (i == _Idx[1])
	       	///r_h.push_back(_f[i] - (2*_v[i] - _v[i-1])/_h/_h);
		r_h.push_back(0);
	    else
		r_h.push_back(_f[i] - (2*_v[i] - _v[i-1]- _v[i+1])/_h/_h);
	}
	_pRestrictOP->SetInput(r_h);
	_pRestrictOP->restrict();
	r_2h = _pRestrictOP->ReturnOutput();
	_nowlevel++;
	UpdateIndex();
	/// Update _f
	_f.erase(_f.begin()+_Idx[0], _f.begin()+_Idx[1]+1);
	_f.insert(_f.begin() + _Idx[0], r_2h.begin(), r_2h.end());
	/// Does r_2h/r_h need to be destructed here?
	r_2h.assign(_Idx[1]-_Idx[0]+1, 0);
	_v.erase(_v.begin()+_Idx[0], _v.begin()+_Idx[1]+1);
	_v.insert(_v.begin() + _Idx[0], r_2h.begin(), r_2h.end());
	VCycle();
	r_2h.assign(_v.begin()+_Idx[0], _v.begin()+_Idx[1]+1);
	_pInterpolateOP->SetInput(r_2h);
	_pInterpolateOP->interpolate();
	r_h = _pInterpolateOP->ReturnOutput();
	_nowlevel--;
	UpdateIndex();
	for (int i = _Idx[0]; i <= _Idx[1]; i++)
	    _v[i] = _v[i] + r_h[i-_Idx[0]];
	for (int i = 0; i < _RlxTimes; i++)
	    WeightedJacobi();
    }
}

void MultigridSolver::BottomSolve()
{
    _v[1] = (_f[0] + _f[2] + 2 * _f[1]) / 2 * _h * _h;
    // _v[0] = (_h * _h * _f[0] + _v[1]) / 2;
    // _v[2] = (_h * _h * _f[2] + _v[1]) / 2;
    _v[0] = 0;
    _v[2] = 0;
}

void MultigridSolver::Solve()
{
    UpdateIndex();
    std::vector<double> all0(_Idx[0], 0);
    for (int i = 0; i < 10; i++)
    {
	_v.erase(_v.begin(), _v.begin()+_Idx[0]);
	_v.insert(_v.begin(), all0.begin(), all0.end());
	_f.erase(_f.begin(), _f.begin()+_Idx[0]);
	_f.insert(_f.begin(), all0.begin(), all0.end());
	VCycle();
    }
}

std::vector<double> MultigridSolver::ReturnSolution()
{
    std::vector<double> result;
    result.insert(result.begin(), _v.begin()+_Idx[0], _v.begin()+_Idx[1]+1);
    return result;
}
