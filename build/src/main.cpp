/*
 
 main.cpp
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <string>

// https://github.com/JohannesBuchner/cuba
#include "cuba.h"

#include "main.h"
#include "paramreader.h"
#include "outputroutines.h"
#include "runningcoupling.h"
#include "processCollision.h"
#include "crosssections.h"
#include "glauber.h"

using namespace std;

std::string outTag = "";  // prefix for output files (set per run/bin)

// choose sigmaNN based on sqrt(s_NN) (mb)
static double select_sigmaNN_mb(double roots) {
    if (fabs(roots - 200.0)  < 1e-6) return 42.0;   // RHIC
    if (fabs(roots - 2760.0) < 1e-6) return 62.0;   // early LHC
    if (fabs(roots - 5023.0) < 1e-6) return 67.6;   // LHC 5.02/5.03
    if (fabs(roots - 8160.0) < 1e-6) return 71.0;   // LHC 8.16
    return 67.6;                                    // fallback
}

// default centrality edges in percent (0..100)
static std::vector<double> default_edges(double roots, int A) {
    std::vector<double> e;
    bool isRHIC = (fabs(roots - 200.0) < 1e-6);
    if (isRHIC && A==197) {       // d+Au at RHIC
        e = {0,20,40,60,80,100};
    } else {                      // p+Pb at LHC
        for (int i=0;i<=20;++i) e.push_back(5.0*i); // 0,5,...,100
    }
    return e;
}

int main(int argc, char *argv[]) {

    print_line(); // cosmetic
    auto starttime = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Reading Parameters and Setting Up Calculation: " << ctime(&starttime);
    print_line(); // cosmetic

    // 1) Read parameters from file and command line
    readParametersFromFile("input/params.txt", 1);
    if (argc > 1) {
        print_line();
        cout << "Parameters from command line" << endl;
        print_line();
        readParametersFromCommandLine(argc, argv, 1);
    }
    processParameters();

    print_line(); // cosmetic

    // 2) sigmaNN for Glauber
    double sigmaNN_mb = select_sigmaNN_mb(rootsnn);
    set_sigmaNN_mb(sigmaNN_mb);

    cout << "Computing Effective Path Lengths (sigmaNN = " << sigmaNN_mb << " mb) ...\n";

    // 3) Branch: min-bias vs centrality
    if (useCentrality == 0) {
        // ---- MIN-BIAS ----
        double LA = compute_LA_minbias_pA(A, rho0, lp);
        cout << std::fixed << setprecision(6)
             << "[MIN-BIAS] rootsNN=" << rootsnn << " GeV, A =" << A
             << " => L_eff = " << LA << " fm\n";

        // Use it automatically
        if (collisionType == 2) {               // AA
            lA = LA; lB = LA;
        } else {                               
            lA = LA; lB = lp; // lB is ignored for pA
        }

        outTag = "minbias_";
        print_line();
        cout << "Quenching calculation started: " << ctime(&starttime);
        create_output_directory();
        print_line();

        processCollision(collisionType);
    } 
    else {
        // ---- CENTRALITY ----
        std::vector<double> edges;
        if (!centralityEdges.empty()) {
            edges = centralityEdges;                    // already in percent
        } else if (cmax > cmin && cmax > 0.0) {
            edges = {100.0*cmin, 100.0*cmax};          // one bin
        } else {
            edges = default_edges(rootsnn, A);
        }

        if (edges.size() < 2) {
            cerr << "centralityEdges has fewer than 2 entries; falling back to min-bias.\n";
            double LA = compute_LA_minbias_pA(A, rho0, lp);
            cout << "[MIN-BIAS fallback] L_eff = " << LA << " fm\n";
        } else {
            cout << "Centrality edges (percent):";
            for (size_t i=0;i<edges.size();++i) cout << (i?", ":" ") << edges[i];
            cout << " %\n";

            // Compute and summarize L_eff for all bins
            std::vector<double> leff_values;
            std::vector<double> npart_values;
            std::vector<std::string> leff_labels;
            for (size_t i=0; i+1<edges.size(); ++i) {
                double c0p = edges[i];
                double c1p = edges[i+1];
                double c0  = c0p/100.0;
                double c1  = c1p/100.0;
                double LA = compute_LA_centrality_pA(c0, c1, A, sigmaNN_mb, rho0, lp);
                leff_values.push_back(LA);
                leff_labels.push_back(std::to_string((int)c0p) + "-" + std::to_string((int)c1p));
                double npart = compute_Npart_centrality_pA(c0, c1, A, sigmaNN_mb, rho0, lp);
                npart_values.push_back(npart);
            }
            print_line();
            std::cout << "Summary of L_eff for all centrality bins:\n";
            for (size_t i=0; i<leff_values.size(); ++i) {
                std::cout << "  Bin " << leff_labels[i] << "%: L_eff = " << leff_values[i] << " fm, " << "<N_part> = " << npart_values[i] << "\n";
            }
            print_line();

            // Now process collisions for each bin
            for (size_t i=0; i<leff_values.size(); ++i) {
                double LA = leff_values[i];
                if (collisionType == 2) { lA = LA; lB = LA; }
                else                   { lA = LA; lB = lp; }
                outTag = "cent_" + leff_labels[i] + "_";
                print_line();
                cout << "Quenching calculation started (bin " << leff_labels[i] << "%): " << ctime(&starttime);
                create_output_directory();
                print_line();
                processCollision(collisionType);
            }
        }
    }
    // done
    print_line();
    auto endtime = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Done: " << ctime(&endtime);
    print_line(); // cosmetic

    return 0;
}

