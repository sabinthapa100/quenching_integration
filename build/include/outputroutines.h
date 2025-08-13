/*
 
 outputroutines.h
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */
#ifndef __outputroutines_h__
#define __outputroutines_h__

#include <vector>
#pragma once
#include <string>

using namespace std;

void print_line();
void create_output_directory();                        // creates "output"
void create_output_directory(const std::string&);      // creates "output/<particle>"
void printResult(string label, double y, double pt, double result, double error);

#endif /* __outputroutines_h__ */
