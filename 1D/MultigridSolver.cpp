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

MultigridSolver::MultigridSolver(int n, std::vector<double> f, std::vector<double> v, std::string S1, std::string S2, std::string S3)
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
    _u0 = *v.begin();
    _u1 = *(v.end()-1);
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
    _TypeofCycle = S3;
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

void MultigridSolver::SetRS(std::vector<double> RS)
{
    _RS = RS;
}

void MultigridSolver::PrintInfo()
{
    std::cout << "The total level: " << _n << std::endl;
    std::cout << "The length of f, v:" << _f.size() << " " << _v.size() << std::endl;
    std::cout << "The boundary condition: " << _u0 << " " << _u1 << std::endl;
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

double MultigridSolver::RE_2Norm()
{
    double RSNorm = 0;
    double eNorm = 0;
    for (int i = 0; i < _RS.size(); i++)
    {
	RSNorm = RSNorm + pow(_RS[i], 2);
	eNorm = eNorm + pow(_RS[i] - _v[_Idx[0]+i], 2);
    }
    return sqrt(eNorm)/sqrt(RSNorm);
}

void MultigridSolver::WeightedJacobi()
{
    std::vector<double> v_star(_Idx[1] - _Idx[0] + 1, 0);
    for (int i = 1; i < v_star.size() - 1; i++)
	v_star[i] = (_v[i-1+_Idx[0]] + _v[i+1+_Idx[0]] + _h*_h*_f[i+_Idx[0]]) / 2;
    for (int i = 1; i < v_star.size() - 1; i++)
	_v[i+_Idx[0]] = _v[i+_Idx[0]] + _w * (v_star[i] - _v[i+_Idx[0]]);
}

void MultigridSolver::VCycle(int StartLevel)
{
    _nowlevel = StartLevel;
    UpdateIndex();
    if (_nowlevel == _n)
	BottomSolve(0, 0);
    else
    {
	for (int i = 0; i < _RlxTimes; i++)
	    WeightedJacobi();
	std::vector<double> r_h;
	std::vector<double> r_2h;
	for (int i = _Idx[0]; i <= _Idx[1]; i++)
	{
	    if (i == _Idx[0])
		r_h.push_back(0);
	    else if (i == _Idx[1])
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
	VCycle(_nowlevel);
	r_2h.assign(_v.begin()+_Idx[0], _v.begin()+_Idx[1]+1);
	_pInterpolateOP->SetInput(r_2h);
	_pInterpolateOP->interpolate();
	r_h = _pInterpolateOP->ReturnOutput();
	_nowlevel--;
	UpdateIndex();
	for (int i = _Idx[0]; i <= _Idx[1]; i++)
	    _v[i] = _v[i] + r_h[i-_Idx[0]];
	for (int i = 0; i < 10; i++)
	    WeightedJacobi();
    }
}

void MultigridSolver::BottomSolve(double u0, double u1)
{
    _v[0] = u0;
    _v[2] = u1;
    _v[1] = (_h*_h*_f[1]+_v[0]+_v[2])/2;
}

void MultigridSolver::Solve()
{
    UpdateIndex();
    if (_TypeofCycle == "VC")
    {
	std::vector<double> all0(_Idx[0], 0);
	for (int i = 0; i < 10; i++)
	{
	    _v.erase(_v.begin(), _v.begin()+_Idx[0]);
	    _v.insert(_v.begin(), all0.begin(), all0.end());
	    _v[_Idx[0]] = _u0;
	    _v[_Idx[1]] = _u1;
	    _f.erase(_f.begin(), _f.begin()+_Idx[0]);
	    _f.insert(_f.begin(), all0.begin(), all0.end());
	    VCycle(1);
	    std::cout << "The iteration step: " << i+1 << ", the relative error: " << RE_2Norm() << std::endl;
	}
    }
    else
    {
	FMG();
	std::cout << "The relative error: " << RE_2Norm() << std::endl;
    }
}

void MultigridSolver::FMG()
{
    UpdateIndex();
    if (_nowlevel == _n)
    {
	BottomSolve(_u0, _u1);
    }
    else
    {
	std::vector<double> f_h;
	std::vector<double> f_2h;
	f_h.insert(f_h.begin(), _f.begin()+_Idx[0], _f.begin()+_Idx[1]+1);
	_pRestrictOP->SetInput(f_h);
	_pRestrictOP->restrict();
	f_2h = _pRestrictOP->ReturnOutput();
	_nowlevel++;
	UpdateIndex();
	_f.erase(_f.begin()+_Idx[0], _f.begin()+_Idx[1]+1);
	_f.insert(_f.begin() + _Idx[0], f_2h.begin(), f_2h.end());
	FMG();
	f_2h.clear();
	f_2h.insert(f_2h.begin(), _v.begin()+_Idx[0], _v.begin()+_Idx[1]+1);
	_pInterpolateOP->SetInput(f_2h);
	_pInterpolateOP->interpolate();
	f_h = _pInterpolateOP->ReturnOutput();
	_nowlevel--;
	UpdateIndex();
	_v.erase(_v.begin()+_Idx[0], _v.begin()+_Idx[1]+1);
	_v.insert(_v.begin() + _Idx[0], f_h.begin(), f_h.end());
	VCycle(_nowlevel);
    }
}

std::vector<double> MultigridSolver::ReturnSolution()
{
    std::vector<double> result;
    result.insert(result.begin(), _v.begin()+_Idx[0], _v.begin()+_Idx[1]+1);
    return result;
}

