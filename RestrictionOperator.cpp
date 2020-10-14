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

RestrictionOperator::RestrictionOperator()
{
}

RestrictionOperator::RestrictionOperator(std::vector<double> a)
{
    _Input = a;
    _Output = std::vector<double> ((_Input.size()+1)/2);
}

void RestrictionOperator::PrintInput()
{
    std::vector<double>::iterator idx;
    for (idx = _Input.begin(); idx != _Input.end(); ++idx)
	std::cout << *idx << " ";
    std::cout << std::endl;
}

void RestrictionOperator::PrintOutput()
{
    std::vector<double>::iterator idx;
    for (idx = _Output.begin(); idx != _Output.end(); ++idx)
	std::cout << *idx << " ";
    std::cout << std::endl;
}

void RestrictionOperator::SetInput(std::vector<double> a)
{
    _Input = a;
    _Output = std::vector<double> ((_Input.size()+1)/2);
}

std::vector<double> RestrictionOperator::ReturnOutput()
{
    return _Output;
}

void RestrictionOperator::PrintType()
{
    std::cout << "No type" << std::endl;
}

void RestrictionOperator::restrict()
{
    std::cerr << "Restriction Error!" << std::endl;
    exit(-1);
}

void InjectionRestriction::restrict()
{
    for (int i = 0; i < _Output.size(); i++)
	_Output[i] = _Input[2*i];
}

void InjectionRestriction::PrintType()
{
    std::cout << "The type of restriction operator: Injection" << std::endl;
}

void FullWeightingRestriction::restrict()
{
    int n = _Input.size();
    _Output[0] = _Input[0];
    _Output[(n+1)/2-1] = _Input[n-1];
    for (int i = 1; i < (n-1)/2; i++)
    {
	_Output[i] = (_Input[2*i-1] + 2 * _Input[2*i] + _Input[2*i+1]);
    }
}

void FullWeightingRestriction::PrintType()
{
    std::cout << "The type of restriction operator: Full-weighting" << std::endl;
}
