/**
 * @file   RestrictionOperator.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 20 21:31:36 2020
 * 
 * @brief  2 kinds of restriction operator for MultigridSolver, 2D case
 * 
 * 
 */
#include "RestrictionOperator.h"

RestrictionOperator::RestrictionOperator()
{
}

void RestrictionOperator::SetInput(std::vector<double> a)
{
    _Input = a;
}

void RestrictionOperator::SetSdLen(int a)
{
    _OldSdLen = a;
    _NewSdLen = (a + 1) / 2;
    _Output = std::vector<double>(_NewSdLen * _NewSdLen, 0); 
}

int RestrictionOperator::CortoIdx(int i, int j, int SdLen)
{
    return (i * SdLen + j);
}

std::vector<double> RestrictionOperator::ReturnOutput()
{
    return _Output;
}

void RestrictionOperator::restrict()
{
    std::cerr << "Restriction Error!" << std::endl;
    exit(-1);
}

void InjectionRestriction::restrict()
{
    for (int i = 0; i < _NewSdLen; i++)
	for (int j = 0; j < _NewSdLen; j++)
		_Output[CortoIdx(i, j, _NewSdLen)] = _Input[CortoIdx(2*i, 2*j, _OldSdLen)];
}

void FullWeightingRestriction::restrict()
{
    for (int i = 0; i < _NewSdLen; i++)
	for (int j = 0; j < _NewSdLen; j++)
	{
	    if (i == 0 || i == _NewSdLen-1 || j == 0 || j == _NewSdLen-1)
		_Output[CortoIdx(i, j, _NewSdLen)] = _Input[CortoIdx(2*i, 2*j, _OldSdLen)];
	    else
	    {
		_Output[CortoIdx(i, j, _NewSdLen)] = (_Input[CortoIdx(2*i-1, 2*j-1, _OldSdLen)] + _Input[CortoIdx(2*i-1, 2*j+1, _OldSdLen)]
						      + _Input[CortoIdx(2*i+1, 2*j-1, _OldSdLen)] + _Input[CortoIdx(2*i+1, 2*j+1, _OldSdLen)]
						      + 2*(_Input[CortoIdx(2*i, 2*j-1, _OldSdLen)] + _Input[CortoIdx(2*i, 2*j+1, _OldSdLen)]
							   + _Input[CortoIdx(2*i-1, 2*j, _OldSdLen)] + _Input[CortoIdx(2*i+1, 2*j, _OldSdLen)])
						      + 4*_Input[CortoIdx(2*i, 2*j, _OldSdLen)])/16.0;
	    }
	}
}
