/**
 * @file   RestrictionOperator.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Wed Oct 14 12:30:33 2020
 * 
 * @brief  2 kinds of restriction operator for MultigridSolver
 * 
 * 
 */

#include <vector>
#include <math.h>

class RestrictionOperator
{
protected:
    std::vector<double> _Input;
public:
    RestrictionOperator(std::vector<double> a);
    void PrintInput();
};

class InjectionRestriction : public RestrictionOperator
{
public:
    std::vector<double> restrict();
};

class FullWeightingRestriction : public RestrictionOperator
{
public:
    std::vector<double> restrict();
};
    

