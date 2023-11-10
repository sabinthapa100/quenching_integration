/*
 
 outputroutines.h
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <vector>

#ifndef __outputroutines_h__
#define __outputroutines_h__

using namespace std;

void print_line();
void create_output_directory();
void printResult(string label, double y, double pt, double result, double error);

#endif /* __outputroutines_h__ */
