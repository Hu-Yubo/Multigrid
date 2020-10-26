/**
 * @file   main.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Mon Oct 26 14:11:59 2020
 * 
 * @brief  test file
 * 
 * 
 */

#include "Element.h"

int main(int argc, char* argv[])
{
    Node N1(0,1);
    Node N2(1,1);
    Node N3(2,1);
    Node N4(3,1);
    Element E(N1, N2, N4, N3);
    for (int i = 0; i < 4; i++)
    {
	for (int j = 0; j < 4; j++)
	    std::cout << E.a_ij(i+1, j+1) << " ";
	std::cout << std::endl;
    }
    for (int i = 0; i < 4; i++)
	std::cout << E.rhs(i+1) << std::endl;
    return 0;
}
