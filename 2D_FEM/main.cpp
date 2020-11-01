/**
 * @file   main.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Mon Oct 26 14:11:59 2020
 * 
 * @brief  test file
 * 
 * 
 */

#include "FEMSpace.h"

double f(double x, double y)
{
    return 2*PI*PI*sin(PI*x)*sin(PI*y);
}

int main(int argc, char* argv[])
{
    FEMSpace FE(2, f);
    FE.GenerateA();
    FE.PrintA();
    FE.GenerateRhs();
    FE.PrintRhs();
    return 0;
}
