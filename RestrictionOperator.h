/**
 * @file   RestrictionOperator.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Wed Oct 14 12:30:33 2020
 * 
 * @brief  2 kinds of restriction operator for MultigridSolver
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
public:
    RestrictionOperator();
    RestrictionOperator(std::vector<double> a);
    void PrintInput();
    void PrintOutput();
    void SetInput(std::vector<double> a);
    std::vector<double> ReturnOutput();
    virtual void PrintType();
    virtual void restrict();
};

class InjectionRestriction : public RestrictionOperator
{
public:
    void restrict();
    void PrintType();
};

class FullWeightingRestriction : public RestrictionOperator
{
public:
    void restrict();
    void PrintType();
};
    

