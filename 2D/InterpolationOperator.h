/**
 * @file   InterpolationOperator.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 20 22:31:55 2020
 * 
 * @brief  2 kinds of interpolation operator for MultigridSolver, 2D case
 * 
 * 
 */

#include <iostream>
#include <vector>
#include <math.h>

class InterpolationOperator
{
protected:
    std::vector<double> _Input;
    std::vector<double> _Output;
    int _OldSdLen;
    int _NewSdLen;
public:
    InterpolationOperator();
    void SetInput(std::vector<double> a);
    void SetSdLen(int a);
    std::vector<double> ReturnOutput();
    int CortoIdx(int i, int j, int SdLen);
    virtual void interpolate();
};

class LinearInterpolation : public InterpolationOperator
{
public:
    void interpolate();
};
