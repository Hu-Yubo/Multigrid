/**
 * @file   Element.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Mon Oct 26 10:00:01 2020
 * 
 * @brief  
 * 
 * 
 */

#include "Element.h"

Element::Element(Node N1,  Node N2, Node N3, Node N4)
{
    _Node[0] = N1;
    _Node[1] = N2;
    _Node[2] = N3;
    _Node[3] = N4;
}

double Element::phi(double xi, double eta, int i)
{
    switch(i)
    {
    case 1:
	return (xi-1)*(eta-1)/4.0;
    case 2:
	return -(xi+1)*(eta-1)/4.0;
    case 3:
	return (xi+1)*(eta+1)/4.0;
    case 4:
	return -(xi-1)*(eta+1)/4.0;
    }
}

