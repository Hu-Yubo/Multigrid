/**
 * @file   MultigridSolver.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 20 09:02:35 2020
 * 
 * @brief  solve 2D possion equation by multigrid method on domain [0,1]x[0,1].
 * 
 * 
 */

#include "MultigridSolver.h"

MultigridSolver::MultigridSolver(int n, std::vector<double> f, std::vector<double> v, std::string S1, std::string S2, std::string S3)
{
    _n = n;
    _h = 1.0/pow(2, _n);
    _SdLen = (int)(pow(2, _n)) + 1;
    for (int i = 0; i < _n; i++)
    {
	int length = (int)(pow(2, _n - i)) + 1;
	_f.push_back(std::vector<double>(length * length, 0));
	_v.push_back(std::vector<double>(length * length, 0));
    }
    _f[0] = f;
    _v[0] = v;
    _BtmBnd.push_back(_v[0][0]);
    _BtmBnd.push_back(_v[0][(_SdLen-1)/2]);
    _BtmBnd.push_back(_v[0][_SdLen-1]);
    _BtmBnd.push_back(_v[0][(_SdLen-1)/2*_SdLen]);
    _BtmBnd.push_back(_v[0][(_SdLen-1)/2*_SdLen+_SdLen-1]);
    _BtmBnd.push_back(_v[0][_SdLen*(_SdLen-1)]);
    _BtmBnd.push_back(_v[0][_SdLen*_SdLen-1-(_SdLen-1)/2]);
    _BtmBnd.push_back(_v[0][_SdLen*_SdLen-1]);
    _nowlevel = 1;
    if (S1 == "FullWeighting")
    {
        _pRestrictOP = new FullWeightingRestriction();
    }
    else if (S1 == "Injection")
    {
	_pRestrictOP = new InjectionRestriction();
    }
    if (S2 == "Linear")
    {
	_pInterpolateOP = new LinearInterpolation();
    }
    _w = 2.0/3;
    _TypeofCycle = S3;
}

void MultigridSolver::SetRS(std::vector<double> RS)
{
    _RS = RS;
}

int MultigridSolver::CortoIdx(int i, int j)
{
    return (i * _SdLen + j);
}

std::vector<int> MultigridSolver::IsBoundary()
{
    int size = _SdLen * _SdLen;
    std::vector<int> BndFlag(size, 0);
    for (int i = 0; i < _SdLen; i++)
    {
	BndFlag[i] = 1;
	BndFlag[size - 1 - i] = 1;
	BndFlag[i * _SdLen] = 1;
	BndFlag[size - i * _SdLen - 1] = 1;
    }
    return BndFlag;
}

double MultigridSolver::RE_2norm()
{
    double RSNorm = 0;
    double eNorm = 0;
    for (int i = 0; i < _RS.size(); i++)
    {
	RSNorm = RSNorm + pow(_RS[i],2);
	eNorm = eNorm + pow(_RS[i] - _v[0][i], 2);
    }
    return sqrt(eNorm)/sqrt(RSNorm);
}

void MultigridSolver::UpdateData()
{
    _SdLen = (int)(pow(2, _n - _nowlevel + 1)) + 1;
    _h = 1.0 / (_SdLen - 1);
    _BndMark = IsBoundary();
}

void MultigridSolver::WeightedJacobi()
{
    std::vector<double> v_star(_SdLen * _SdLen, 0);
    for (int i = 0; i < _v[_nowlevel-1].size(); i++)
    {
	if (_BndMark[i] == 0)
	    v_star[i] = (_v[_nowlevel-1][i-1] + _v[_nowlevel-1][i+1]
			 + _v[_nowlevel-1][i-_SdLen] + _v[_nowlevel-1][i+_SdLen]
			 + _h * _h * _f[_nowlevel-1][i]) / 4.0;
	else
	    v_star[i] = _v[_nowlevel-1][i];
    }
    for (int i = 0; i < _v[_nowlevel-1].size(); i++)
	_v[_nowlevel-1][i] = (1.0 - _w) * _v[_nowlevel-1][i] + _w * v_star[i];	  
}

void MultigridSolver::BottomSolve(std::vector<double> BtmBnd)
{
    _v[_nowlevel-1][0] = BtmBnd[0];
    _v[_nowlevel-1][1] = BtmBnd[1];
    _v[_nowlevel-1][2] = BtmBnd[2];
    _v[_nowlevel-1][3] = BtmBnd[3];
    _v[_nowlevel-1][5] = BtmBnd[4];
    _v[_nowlevel-1][6] = BtmBnd[5];
    _v[_nowlevel-1][7] = BtmBnd[6];
    _v[_nowlevel-1][8] = BtmBnd[7];
    _v[_nowlevel-1][4] = (_h*_h*_f[_nowlevel-1][4] + _v[_nowlevel-1][1] + _v[_nowlevel-1][3] + _v[_nowlevel-1][5] + _v[_nowlevel-1][7])/4.0; 
}

void MultigridSolver::VCycle(int StartLevel)
{
    _nowlevel = StartLevel;
    UpdateData();
    if (_nowlevel == _n)
    {
	std::vector<double> a(8,0);
	BottomSolve(a);
    }
    else
    {
	for (int i = 0; i < _RlxTimes; i++)
	    WeightedJacobi();
	std::vector<double> r_h;
	std::vector<double> r_2h;
	for (int i = 0; i < _v[_nowlevel-1].size(); i++)
	{
	    if (_BndMark[i] == 1)
		r_h.push_back(0);
	    else
		r_h.push_back(_f[_nowlevel-1][i] - (4*_v[_nowlevel-1][i] - _v[_nowlevel-1][i-1] - _v[_nowlevel-1][i+1]
						    - _v[_nowlevel-1][i-_SdLen] - _v[_nowlevel-1][i+_SdLen])/(_h*_h));
	}
	_pRestrictOP->SetInput(r_h);
	_pRestrictOP->SetSdLen(_SdLen);
	_pRestrictOP->restrict();
	r_2h = _pRestrictOP->ReturnOutput();
	_nowlevel++;
	UpdateData();
	_f[_nowlevel-1] = r_2h;
	VCycle(_nowlevel);
	r_2h = _v[_nowlevel-1];
	_pInterpolateOP->SetInput(r_2h);
	_pInterpolateOP->SetSdLen(_SdLen);
	_pInterpolateOP->interpolate();
	r_h = _pInterpolateOP->ReturnOutput();
	_nowlevel--;
	UpdateData();
	for (int i = 0; i < _v[_nowlevel-1].size(); i++)
	    _v[_nowlevel-1][i] = _v[_nowlevel-1][i] + r_h[i];
	for (int i = 0; i < _RlxTimes; i++)
	    WeightedJacobi();
    }
}

void MultigridSolver::Solve()
{
    UpdateData();
    if (_TypeofCycle == "VC")
    {
	for(int i = 0; i < 10; i++)
	{
	    for (int j = 1; j < _n; j++)
	    {
		int length = (int)(pow(2, _n - j)) + 1;
		_f[j] = std::vector<double>(length * length, 0);
		_v[j] = std::vector<double>(length * length, 0);
	    }
	    VCycle(1);
	    if (!_RS.empty())
		std::cout << "The iteration step: " << i+1 << ", the relative error: " << RE_2norm() << std::endl;
	}
    }
    else
    {
	FMG();
	if (!_RS.empty())
	    std::cout << "the relative error: " << RE_2norm() << std::endl;
    }
}

void MultigridSolver::FMG()
{
    UpdateData();
    if (_nowlevel == _n)
	BottomSolve(_BtmBnd);
    else
    {
	std::vector<double> f_h = _f[_nowlevel-1];
	std::vector<double> f_2h;
	_pRestrictOP->SetInput(f_h);
	_pRestrictOP->SetSdLen(_SdLen);
	_pRestrictOP->restrict();
	f_2h = _pRestrictOP->ReturnOutput();
	_nowlevel++;
	UpdateData();
	_f[_nowlevel-1] = f_2h;
	FMG();
	f_2h = _v[_nowlevel-1];
	_pInterpolateOP->SetInput(f_2h);
	_pInterpolateOP->SetSdLen(_SdLen);
	_pInterpolateOP->interpolate();
	f_h = _pInterpolateOP->ReturnOutput();
	_nowlevel--;
	UpdateData();
	_v[_nowlevel-1] = f_h;
	VCycle(_nowlevel);
    }
}


std::vector<double> MultigridSolver::ReturnSolution()
{
    return _v[0];
}
    
