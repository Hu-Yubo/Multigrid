/**
 * @file   InterpolationOperator.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Wed Oct 14 21:37:59 2020
 * 
 * @brief  2 kinds of interpolation operator for MultigridSolver
 * 
 * 
 */

#include <iostream>
#include <vector>
#include <math.h>

class InterpolationOperator
{
public:
    std::vector<double> _Input;
    std::vector<double> _Output;
    InterpolationOperator();
    InterpolationOperator(std::vector<double> a);
    void PrintInput();
    void PrintOutput();
    void SetInput(std::vector<double> a);
    std::vector<double> ReturnOutput();
    virtual void PrintType();
    virtual void interpolate();
};

class LinearInterpolation : public InterpolationOperator
{
public:
    void interpolate();
    void PrintType();
};
