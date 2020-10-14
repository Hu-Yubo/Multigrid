/**
 * @file   InterpolationOperator.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Wed Oct 14 22:09:19 2020
 * 
 * @brief  2 kinds of interpolation operator for MultigridSolver
 * 
 * 
 */

#include "InterpolationOperator.h"

InterpolationOperator::InterpolationOperator()
{
}

InterpolationOperator::InterpolationOperator(std::vector<double> a)
{
    _Input = a;
    _Output = std::vector<double> (_Input.size()*2-1);
}

void InterpolationOperator::PrintInput()
{
    std::vector<double>::iterator idx;
    for (idx = _Input.begin(); idx != _Input.end(); ++idx)
	std::cout << *idx << " ";
    std::cout << std::endl;
}

void InterpolationOperator::PrintOutput()
{
    std::vector<double>::iterator idx;
    for (idx = _Output.begin(); idx != _Output.end(); ++idx)
	std::cout << *idx << " ";
    std::cout << std::endl;
}

void InterpolationOperator::SetInput(std::vector<double> a)
{
    _Input = a;
    _Output = std::vector<double> (_Input.size()*2-1);
}

std::vector<double> InterpolationOperator::ReturnOutput()
{
    return _Output;
}

void InterpolationOperator::PrintType()
{
    std::cout << "No type" << std::endl;
}

void InterpolationOperator::interpolate()
{
    std::cerr << "Restriction Error!" << std::endl;
    exit(-1);
}

void LinearInterpolation::PrintType()
{
    std::cout << "The type of interpolation operator: Linear" << std::endl;
}

void LinearInterpolation::interpolate()
{
    int n = _Input.size();
    for (int i = 0; i < _Output.size(); i++)
    {
	if (i % 2 == 0)
	    _Output[i] = _Input[i/2];
	else
	    _Output[i] = (_Input[i/2] + _Input[i/2 + 1]) / 2;
    }
}
