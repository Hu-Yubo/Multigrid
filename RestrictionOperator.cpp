/**
 * @file   RestrictionOperator.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Wed Oct 14 12:54:52 2020
 * 
 * @brief  2 kinds of restriction operator for MultigridSolver
 * 
 * 
 */

#include "RestrictionOperator.h"

RestrictionOperator::RestrictionOperator(std::vector<double> a)
{
    _Input = a;
}

void RestrictionOperator::PrintInput()
{
    std::vector<double>::iterator idx;
    for (idx = _Input.begin(); idx != _Input.end(); ++i)
	std::cout << *idx << " "; 
}

std::vector<double> InjectionRestriction::restrict()
{
    int n = _Input.size();
    std::vector<double> Output((n+1)/2);
    for (int i = 0; i < Output.size(); i++)
    {
	Output[i] = _Input[2*i];
    }
    return Output;
}

std::vector<double> FullWeightingRestriction::restrict()
{
    int n = _Input.size();
    std::vector<double> Output((n+1)/2);
    Output[0] = _Input[0];
    Output[(n+1)/2-1] = _Input[n-1];
    for (int i = 1; i < (n-1)/2; i++)
    {
	Output[i] = (_Input[2*i-1] + 2 * _Input[2*i] + _Input[2*i+1]);
    }
    return Output;
}
