/*
 
 paramreader.h
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */
#ifndef __paramreader_h__
#define __paramreader_h__

#include <string>
#include <vector>

void readParametersFromFile(std::string filename, int echo);
void readParametersFromCommandLine(int argc, char** argv, int echo);
void processParameters();

// Extern globals parsed by paramreader.cpp
extern int    A;                 // 208 Pb, 197 Au, etc.
extern double rho0;              // 0.17 fm^-3 default
extern int    useCentrality;     // 0=min-bias, 1=centrality
extern double cmin, cmax;        // single-bin fractions (0..1)
extern std::vector<double> centralityEdges;  // percent edges (0..100)

#endif /* __paramreader_h__ */
