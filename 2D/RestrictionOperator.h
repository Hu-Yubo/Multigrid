/**
 * @file   RestrictionOperator.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 20 21:17:53 2020
 * 
 * @brief  2 kinds of restriction operator for MultigridSolver, 2D case
 * 
 * 
 */

#include <iostream>
#include <vector>
#include <math.h>

class RestrictionOperator
{
protected:
    std::vector<double> _Input;
    std::vector<double> _Output;
    int _OldSdLen;
    int _NewSdLen;
public:
    RestrictionOperator();
    void SetInput(std::vector<double> a);
    void SetSdLen(int a);
    std::vector<double> ReturnOutput();
    int CortoIdx(int i, int j, int SdLen);
    virtual void restrict();
};


class InjectionRestriction : public RestrictionOperator
{
public:
    void restrict();
};

class FullWeightingRestriction : public RestrictionOperator
{
public:
    void restrict();
};


