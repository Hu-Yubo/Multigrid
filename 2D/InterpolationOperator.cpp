/**
 * @file   InterpolationOperator.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 20 22:42:30 2020
 * 
 * @brief  2 kinds of interpolation operator for MultigridSolver, 2D case
 * 
 * 
 */

#include "InterpolationOperator.h"

InterpolationOperator::InterpolationOperator()
{
}

void InterpolationOperator::SetInput(std::vector<double> a)
{
    _Input = a;
}

void InterpolationOperator::SetSdLen(int a)
{
    _OldSdLen = a;
    _NewSdLen = a * 2 - 1;
    _Output = std::vector<double>(_NewSdLen * _NewSdLen);
}

int InterpolationOperator::CortoIdx(int i, int j, int SdLen)
{
    return (i * SdLen + j);
}

std::vector<double> InterpolationOperator::ReturnOutput()
{
    return _Output;
}

void InterpolationOperator::interpolate()
{
    std::cerr << "Interpolation Error!" << std::endl;
    exit(-1);
}

void LinearInterpolation::interpolate()
{
    for (int i = 0; i < _NewSdLen; i++)
	for (int j = 0; j < _NewSdLen; j++)
	{
	    if (i%2 == 0 && j%2 == 0)
		_Output[CortoIdx(i, j, _NewSdLen)] = _Input[CortoIdx(i/2, j/2, _OldSdLen)];
	    else if (i%2 != 0 && j%2 == 0)
		_Output[CortoIdx(i, j, _NewSdLen)] = (_Input[CortoIdx(i/2, j/2, _OldSdLen)] + _Input[CortoIdx(i/2+1, j/2, _OldSdLen)]) / 2.0;
	    else if (i%2 == 0 && j%2 != 0)
		_Output[CortoIdx(i, j, _NewSdLen)] = (_Input[CortoIdx(i/2, j/2, _OldSdLen)] + _Input[CortoIdx(i/2, j/2+1, _OldSdLen)]) / 2.0;
	    else
		_Output[CortoIdx(i, j, _NewSdLen)] = (_Input[CortoIdx(i/2, j/2, _OldSdLen)] + _Input[CortoIdx(i/2+1, j/2, _OldSdLen)]
						   + _Input[CortoIdx(i/2, j/2+1, _OldSdLen)] + _Input[CortoIdx(i/2+1, j/2+1, _OldSdLen)])/4.0;
	}
}
