/*
 
 paramreader.cpp
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <stdlib.h>
#include <chrono>
#include <string.h>
#include <algorithm>
#include "paramreader.h"
#include "outputroutines.h"

#include "main.h"

using namespace std;


double massQQ, p0, m, n; //m,n are pp params

// --- define storage for globals declared extern in header ---
int    A        = 208;
double rho0     = 0.17;       // fm^-3
int    useCentrality = 0;     // 0: min-bias
double cmin     = 0.0;        // fractions (0..1)
double cmax     = 0.0;
std::vector<double> centralityEdges;  // percent (0..100), optional

// this workhorse examines a key to see if it corresponds to a var we are setting
// and then attempts to set the var corresponding to key by converting value to the
// appropriate type.  lots of hardcoding here
void setParameter(const char *key, const char *value) {
    if (strcmp(key,"collisionType")==0) collisionType = atoi(value);
    if (strcmp(key,"A")==0)             A             = atoi(value);
    if (strcmp(key,"useCentrality")==0) useCentrality = atoi(value);   // was atof
    if (strcmp(key, "particleType")==0) particleType  = atoi(value);
    if (strcmp(key,"alphas")==0)        alphas        = atof(value);
    if (strcmp(key,"massQQ")==0)        massQQ        = atof(value);
    if (strcmp(key,"lambdaQCD")==0)     lambdaQCD     = atof(value);   // removed stray space
    if (strcmp(key,"qhat0")==0)         qhat0         = atof(value);
    if (strcmp(key,"lp")==0)            lp            = atof(value);
    if (strcmp(key,"massp")==0)         massp         = atof(value);
    if (strcmp(key,"rootsnn")==0)       rootsnn       = atof(value);
    if (strcmp(key,"p0")==0)            p0            = atof(value);
    if (strcmp(key,"m")==0)             m             = atof(value);
    if (strcmp(key,"n")==0)             n             = atof(value);
    if (strcmp(key,"Ny")==0)            Ny            = atoi(value);   // ints
    if (strcmp(key,"Npt")==0)           Npt           = atoi(value);
    if (strcmp(key,"y_min")==0)         y_min         = atof(value);
    if (strcmp(key,"y_max")==0)         y_max         = atof(value);
    if (strcmp(key,"ptmin")==0)         ptmin         = atof(value);
    if (strcmp(key,"ptmax")==0)         ptmax         = atof(value);

    // NEW: centralityEdges CSV parser (0..100 percent)
    if (strcmp(key,"centralityEdges")==0) {
        centralityEdges.clear();
        std::string s(value);
        for (char &ch : s) if (ch==';') ch=',';                 // allow semicolons
        std::stringstream ss(s);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            if (!tok.empty()) centralityEdges.push_back(atof(tok.c_str()));
        }
        // sanitize
        for (double &v : centralityEdges) v = std::min(100.0, std::max(0.0, v));
        std::sort(centralityEdges.begin(), centralityEdges.end());
        centralityEdges.erase(std::unique(centralityEdges.begin(), centralityEdges.end()), centralityEdges.end());
    }
}


//
// This routine assumes that parameters are in text file with
// each parameter on a new line in the format
//
// PARAMKEY	PARAMVALUE
//
// The PARAMKEY must begin the line and only tabs and spaces
// can appear between the PARAMKEY and PARAMVALUE.
//
// Lines which begin with 'commentmarker' defined below are ignored
//
void readParametersFromFile(string filename, int echo) {
    
    string commentmarker = "//";
    char space = ' ';
    char tab = '\t';
    
    int maxline = 128; // maximum line length used in the buffer for reading
    char buffer[maxline];
    ifstream paramFile(filename.c_str());
    
    while(!paramFile.eof()) {
        paramFile.getline(buffer,maxline,'\n');
        string line = buffer; int length = strlen(buffer);
        if (line.substr(0,commentmarker.length())!=commentmarker && line.length()>0) {
            char key[64]="",value[64]="";
            int founddelim=0;
            for (int i=0;i<length;i++) {
                if (buffer[i]==space || buffer[i]==tab) founddelim=1;
                else {
                    if (founddelim==0) key[strlen(key)] = buffer[i];
                    else value[strlen(value)] = buffer[i];
                }
            }
            if (strlen(key)>0 && strlen(value)>0) {
                setParameter(key,value);
                if (echo) cout << key << " = " << value << endl;
            }
        }
    }
    
    return;
}

//
// Read parameters from commandline
//
void readParametersFromCommandLine(int argc, char** argv, int echo) {
    int optind = 1;
    while (optind < argc)
    {
        if (argv[optind][0]=='-') {
            string key = argv[optind];
            key = key.substr(1,key.length()); // remove '-'
            string value = argv[optind+1]; // load value
            if (echo) cout << key << " = " << value << endl;
            setParameter(key.c_str(),value.c_str());
            optind++;
        }
        optind++;
    }
    return;
}


void processParameters() {
    beamRap = acosh(rootsnn / (2.0 * massp));
    xA0 = 1.0 / (2.0 * massp * lA / HBARC);
    xB0 = 1.0 / (2.0 * massp * lB / HBARC);
    dy = (y_max - y_min) / (Ny - 1);
    dpt = (ptmax - ptmin) / (Npt - 1);
    
    
    if (particleType == 0) { // Upsilon
       switch (upsilonState) {
          case 0:
            massQQ = 9.95; // Average mass [DEFAULT]
            cout << "<<<<< Upsilon (Average Mass): " << massQQ << " GeV >>>>>>>" << endl;
          break;
         case 1:
            massQQ = 9.46; //1S state
            cout << "f'<<<<< Upsilon (1S): " << massQQ << " GeV >>>>>>>" << endl;
          break;
        case 2:
            massQQ = 10.02326; //2S state
            cout << "<<<<< Upsilon (2S): " << massQQ << " GeV >>>>>>>" << endl;
          break;
        case 3:
            massQQ = 10.3552; //3S state
            cout << "<<<<< Upsilon (3S): " << massQQ << " GeV >>>>>>>" << endl;
            break;
        default:
          massQQ = 9.95;// Average  Bottomonia Mass
            // Handle unexpected state
          cerr << "Invalid Upsilon state specified." << endl;
        
          cout << "<<<<< Taking Upsilon (Average Mass): " << massQQ << " GeV >>>>>>>" << endl;
        
          //DEFAULT-- Upsilon-- parameters useful in the pp cross-section parametrization
      
          cout << "<<<<< Upsilon State with Average Mass: " << massQQ << " GeV >>>>>>>" << endl;
          p0 = 6.6; // Ref: https://arxiv.org/abs/1304.0901
          m = 2.8;
          n = 13.8;
          cout << "<<<<< pp differential cross-section parameters for Upsilon (m, n): " << m << "," << n << endl;
        break;
        }
        
    } 
    
    else if (particleType == 1) { // J/Psi
        double massQQ = 3.43; //average mass of J/psi, Chi(1P), Psi(2S)
        // J/Psi -- parameters useful in the pp cross-section parametrization
	//double massQQ = 3.0969; //J/Psi mass only
        p0 = 4.2;
        m = 3.5;
        n = 19.2;
        cout << "<<<<< J/Psi State with Mass = " << massQQ << " GeV >>>>>>>" << endl;
        cout << "<<<<< pp differential cross-section parameters for J/psi (m, n): " << m << "," << n << endl;
    } 
    
    else {
        // Handle unexpected particle type
        cerr << "Invalid particle type specified." << endl;
    }

    print_line();

    if (alphas == 0) {
        cout << "==> Using 4-loop running coupling" << endl;
    }
}

