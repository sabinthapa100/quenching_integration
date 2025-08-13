/*
 
 processCollision.cpp
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include "crosssections.h"
#include "processCollision.h"
#include "outputroutines.h"

#include "main.h"
#include <fstream>
#include <string>
using namespace std;

// 0 for Upsilon, 1 for J/Psi
int particleType = 0;  
// Upsilon state: 0 for average [0 Default], 1 for 1S, 2 for 2S, 3 for 3S
int upsilonState = 0; // Default

// --- NEW: compact particle tag for filenames/labels ---
static std::string particleTag()
{
    if (particleType == 1) return "JPsi";
    // Upsilon family
    switch (upsilonState) {
        case 0: default: return "Upsilon";
        case 1: return "Upsilon1S";
        case 2: return "Upsilon2S";
        case 3: return "Upsilon3S";
    }
}

void processCollision(int collisionType) {
    double resultpA, resultRpA, resultAB, resultRAB;
    double errorpA, errorAB, errorRpA, errorRAB;
    auto sanitize = [](double& x){ if (!std::isfinite(x) || x < 0.0) x = 0.0; };
    auto safe_div = [](double num, double den){
          const double tiny = 1e-300;
          return (std::isfinite(num) && std::isfinite(den) && std::fabs(den) > tiny) ? num/den : 0.0; };

    if (collisionType == 0) {cout << "DEFAULT: pA and AB Collision!" << endl;}
    if (collisionType == 1) {cout << "DEFAULT: pA Collision Only!" << endl;}
    if (collisionType == 2) {cout << "DEFAULT: AB Collision Only!" << endl;}
    
    
    ofstream output_file_1, output_file_2, output_file_3, output_file_4, output_file_5;
    
    // --- NEW: prefix includes outTag + particle tag ---
    const std::string ptag   = particleTag();
    create_output_directory(ptag);
    const std::string prefix = "output/" + ptag + "/" + outTag + ptag + "_";
    
    // pp
    output_file_1.open( (prefix + "pp-cross-section.tsv").c_str() );

    // pA
    if (collisionType == 0 || collisionType == 1) {
        output_file_2.open( (prefix + "pA-cross-section.tsv").c_str() );
        output_file_3.open( (prefix + "RpA.tsv").c_str() );
    }
    // AB
    if (collisionType == 0 || collisionType == 2) {
        output_file_4.open( (prefix + "AB-cross-section.tsv").c_str() );
        output_file_5.open( (prefix + "RAB.tsv").c_str() );
    }
    
    for (int i = 0; i < Ny; i++) {
        double y = y_min + i * dy;
        for (int j = 0; j < Npt; j++) {
            double pt = ptmin + j * dpt;

            // output pp cross section
            double resultPP = dsigdyd2pt(y, pt);
            output_file_1 << y << "\t" << pt << "\t" << resultPP << endl;
            	
            if (collisionType == 0 || collisionType == 1) {
                // compute pA cross section
                pACrossSection(y, pt, &resultpA, &errorpA);
                sanitize(resultpA);  sanitize(errorpA);
                output_file_2 << y << "\t" << pt << "\t" << resultpA << "\t" << errorpA << endl;

                // output RpA
                resultRpA = safe_div(resultpA, resultPP);
                errorRpA  = (resultpA > 1e-300) ? safe_div(errorpA, resultpA) : 0.0;

                // --- NEW: include particle name in printed label ---
                std::string lblRpA = "RpA: " + ptag;
                printResult(lblRpA.c_str(), y, pt, resultRpA, errorRpA);

                output_file_3 << y << "\t" << pt << "\t" << resultRpA << "\t" << errorRpA << endl;
            }

            if (collisionType == 0 || collisionType == 2) {
                // compute AB cross section
                ABCrossSection(y, pt, &resultAB, &errorAB);
                sanitize(resultAB);  sanitize(errorAB);
                output_file_4 << y << "\t" << pt << "\t" << resultAB << "\t" << errorAB << endl;

                // output RAB
                resultRAB = safe_div(resultAB, resultPP);
                errorRAB  = (resultAB > 1e-300) ? safe_div(errorAB, resultAB) : 0.0;

                // --- NEW: include particle name in printed label ---
                std::string lblRAB = "RAB: " + ptag;
                printResult(lblRAB.c_str(), y, pt, resultRAB, errorRAB);

                output_file_5 << y << "\t" << pt << "\t" << resultRAB << "\t" << errorRAB << endl;
            }
        }
    }
		
    if (collisionType == 0 || collisionType == 1) {
    	output_file_1.close();
        output_file_2.close();
        output_file_3.close();
    }
    if (collisionType == 0 || collisionType == 2) {
    	output_file_1.close();
        output_file_4.close();
        output_file_5.close();
    }
}
