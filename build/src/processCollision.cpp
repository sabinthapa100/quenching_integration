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
using namespace std;

// 0 for Upsilon, 1 for J/Psi
int particleType = 0;  
// Upsilon state: 0 for 1S, 1 for 2S, 2 for 3S, 3 for average [3 Default]
int upsilonState = 3;

void processCollision(int collisionType) {
    double resultpA, resultRpA, resultAB, resultRAB;
    double errorpA, errorAB, errorRpA, errorRAB;

    if (collisionType == 0) {cout << "DEFAULT: pA and AB Collision!" << endl;}
    if (collisionType == 1) {cout << "DEFAULT: pA Collision Only!" << endl;}
    if (collisionType == 2) {cout << "DEFAULT: AB Collision Only!" << endl;}
    
    
    ofstream output_file_1, output_file_2, output_file_3, output_file_4, output_file_5;

    
    output_file_1.open("output/pp-cross-section.tsv");

    if (collisionType == 0 || collisionType == 1) {
        output_file_2.open("output/pA-cross-section.tsv");
        output_file_3.open("output/RpA.tsv");
    }
    if (collisionType == 0 || collisionType == 2) {
        output_file_4.open("output/AB-cross-section.tsv");
        output_file_5.open("output/RAB.tsv");
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
                output_file_2 << y << "\t" << pt << "\t" << resultpA << "\t" << errorpA << endl;

                // output RpA
                resultRpA = resultpA / resultPP;
                errorRpA = errorpA / resultpA;
                printResult("RpA", y, pt, resultRpA, errorRpA);
                output_file_3 << y << "\t" << pt << "\t" << resultRpA << "\t" << errorRpA << endl;
            }

            if (collisionType == 0 || collisionType == 2) {
                // compute AB cross section
                ABCrossSection(y, pt, &resultAB, &errorAB);
                output_file_4 << y << "\t" << pt << "\t" << resultAB << "\t" << errorAB << endl;

                // output RAB
                resultRAB = resultAB / resultPP;
                errorRAB = errorAB / resultAB;
                printResult("RAB", y, pt, resultRAB, errorRAB);
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

