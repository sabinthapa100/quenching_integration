/*
 
 main.cpp
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <iostream>
#include <chrono>

// https://github.com/JohannesBuchner/cuba
#include "cuba.h"

#include "main.h"
#include "paramreader.h"
#include "outputroutines.h"
#include "runningcoupling.h"
#include "processCollision.h"
#include "crosssections.h"

using namespace std;

int main(int argc, char *argv[]) {
    
    print_line(); // cosmetic
    auto starttime = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Quenching calculation started: " << ctime(&starttime);
    print_line(); // cosmetic
    
    // read parameters from file and command line
    readParametersFromFile("input/params.txt", 1);
    if (argc > 1) {
        print_line();
        cout << "Parameters from command line" << endl;
        print_line();
        readParametersFromCommandLine(argc, argv, 1);
    }
    processParameters();
    
    print_line(); // cosmetic
    
    create_output_directory();
    
    print_line(); // cosmetic

    switch (collisionType) {
        case 0: // pA and AB results (DEFAULT - BOTH)
        case 1: // pA Only Result
        case 2: // AB Only Result
            processCollision(collisionType);
            break;
        default:
            cerr << "Invalid collision type. Exiting..." << endl;
            return 1;
    }

    // print done!
    print_line();

    auto endtime = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Done: " << ctime(&endtime);
    print_line(); // cosmetic

    return 0;
}

