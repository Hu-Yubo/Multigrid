/**
 * @file   FEMSpace.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 27 00:03:04 2020
 * 
 * @brief  A space for generating and storing some components of FEM.
 *         2D possion equation
 * 
 */

#include "FEMSpace.h"

FEMSpace::FEMSpace(int SdLen, double(*func)(double, double))
{
    _SdLen = SdLen;
    _DIM = _SdLen * _SdLen;
    _NEle = (_SdLen - 1) * (_SdLen - 1);
    _A = StiffMat(_DIM);
    _rhs = std::vector<double>(_DIM);
    _func = func;
    _BaseEle = Element();
}

FEMSpace::FEMSpace(int SdLen)
{
    _SdLen = SdLen;
    _DIM = _SdLen * _SdLen;
    _NEle = (_SdLen - 1) * (_SdLen - 1);
    _A = std::vector<std::map<int, double>>(_DIM);
    _rhs = std::vector<double>(_DIM);
    _BaseEle = Element();
}

std::vector<int> FEMSpace::NodeofEle(int i)
{
    std::vector<int> Idx;
    int a = i / (_SdLen-1);
    int b = i % (_SdLen-1);
    Idx.push_back(a*_SdLen + b);
    Idx.push_back(Idx[0] + 1);
    Idx.push_back(Idx[1] + _SdLen);
    Idx.push_back(Idx[2] - 1);
    return Idx;
}

void FEMSpace::GenerateA()
{
    /// traverse each element
    for (int k = 0; k < _NEle; k++)
    {
	std::vector<int> Idx = NodeofEle(k);
	_BaseEle = Element(Node(Idx[0], _SdLen), Node(Idx[1], _SdLen), Node(Idx[2], _SdLen), Node(Idx[3], _SdLen));
	for (int i = 1; i <= 4; i++)
	    for (int j = 1; j <= 4; j++)
	    {
		_A[_BaseEle.NdIdx(i)][_BaseEle.NdIdx(j)] += _BaseEle.a(i, j);
	    }
    }
}

void FEMSpace::PrintA()
{
    for (int i = 0; i < _A.size(); i++)
    {
	std::map<int, double>::iterator it;
	std::cout << i << ": ";
	for (it = _A[i].begin(); it != _A[i].end(); it++)
	{
	    std::cout << it->first << " " << it->second << " ";
	}
	std::cout << std::endl;
    }
}
