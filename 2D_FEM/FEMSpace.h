/**
 * @file   FEMSpace.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Mon Oct 26 23:20:09 2020
 * 
 * @brief  A space for generating and storing some components of FEM.
 *         2D possion equation 
 * 
 */

#include "Element.h"
#include <map>

/// Use std::map to store stiff matrix, idea by LSJ
typedef std::vector<std::map<int, double>> StiffMat;

class FEMSpace
{
private:
    Element _BaseEle;
    /// The count of dofs in each row
    int _SdLen;
    int _DIM;
    /// The total count of elements
    int _NEle;
    StiffMat _A;
    std::vector<double> _rhs;
    double(*_func) (double, double);

public:
    FEMSpace(int SdLen, double(*func)(double, double));
    FEMSpace(int SdLen);
    /** 
     * @brief Output the global indexes of all dofs in ith element, in order.
     * 
     * @param i 
     * 
     * @return 
     */
    std::vector<int> NodeofEle(int i);
    void GenerateA();
    void PrintA();
    void GenerateRhs();
    void PrintRhs();
};
