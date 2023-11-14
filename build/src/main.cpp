/*
 
 main.cpp
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <iostream>
#include <fstream>

#include "main.h"
#include "cuba.h"
#include "paramreader.h"
#include "outputroutines.h"

using namespace std;

//
//  main.cpp
//
//  Implements energy loss and momentum broadening in AB collisions
//

//
// Global params; these default values are overridden by the params file
//
//Default Collision Type is pA, For AB, collisionType = 1 (overridden by params)
int collisionType = 0; 
int nc = 3;
double alphas = 0.5; // QCD coupling constant
double lambdaQCD = 0.25;
double qhat0 = 0.075; // GeV^2/fm
double lp = 1.5; // in fm
double lA = 10.11; // in fm
double lB = 10.11; // in fm
double massp = 0.938; // mass of proton in GeV, it can also be 1 GeV
double rootsnn = 5023; // collision energy, sqrt(s_NN)

// Computed from params file
double beamRap, xA0, xB0;

//DEFAULT-- Upsilon-- parameters useful in the pp cross-section parametrization
double massQQ = 9.46; //Upsilon mass
double p0 = 6.6; // in GeV for Upsilon State at 7 TeV
double m = 2.8;
double n = 13.8;

// J/Psi -- parameters useful in the pp cross-section parametrization
//double massQQ = 3.0969; //J/Psi mass
//double p0 = 4.2; // in GeV for J/Psi State at 7 TeV
//double m = 3.5;
//double n = 19.0;

// Parameters for loops in y and pt
int Ny = 10*2+1;
int Npt = 40*2+1;
double y_min = -5.0;
double y_max = 5.0;
double ptmin = 0.1;
double ptmax = 40.1;
double dy = (y_max - y_min) / (Ny-1);
double dpt = (ptmax - ptmin) / (Npt-1);

inline double Mperp2(double pt) { return pt * pt + massQQ * massQQ; }
inline double Mperp(double pt) { return sqrt(Mperp2(pt));}

inline double ymax(double pt) { return log(rootsnn / Mperp(pt));}

// transport momentum coefficient
inline double qhat(double x) { return qhat0 * pow(pow(10, -2) / x, 0.3); }

//x-values
inline double xA2(double y, double pt) { return Mperp(pt) / rootsnn * exp(-y); }
inline double xB2(double y, double pt) { return Mperp(pt) / rootsnn * exp(-y); }
inline double myXA(double y, double pt) { return min(xA0, xA2(y, pt)); }
inline double myXB(double y, double pt) { return min(xB0, xB2(y, pt));}

//transverse momentum broadening values
inline double lA2(double y, double pt) { return qhat(myXA(y, pt)) * lA;}
inline double lB2(double y, double pt) { return qhat(myXB(y, pt)) * lB;}
inline double lAp2(double y, double pt) { return qhat(myXA(y, pt)) * lp;}
inline double lBp2(double y, double pt) { return qhat(myXB(y, pt)) * lp;}

inline double dptA(double y, double pt) { return sqrt(lA2(y, pt) - lAp2(y, pt));}
inline double dptB(double y, double pt) { return sqrt(lB2(y, pt) - lBp2(y, pt));}

inline double LambdaAp2(double y, double pt) { return max(lambdaQCD * lambdaQCD, lAp2(y, pt));}
inline double LambdaBp2(double y, double pt) { return max(lambdaQCD * lambdaQCD, lBp2(y, pt));}

inline double dymax(double y, double pt) {
    double log2 = log(2.0);
    double ymax_val = ymax(pt);
    double result = min(log2, ymax_val - y);
    return result;
}

inline double calculatePolyLog1(double z, double y, double pt) {
    double coeff;
    coeff = - lA2(y, pt)/(z*z*Mperp2(pt));
    return gsl_sf_dilog(coeff);
}

inline double calculatePolyLog2(double z, double y, double pt) {
    double coeff;
    coeff = - LambdaAp2(y, pt)/(z*z*Mperp2(pt));
    return gsl_sf_dilog(coeff);
}

inline double calculatePolyLog3(double z, double y, double pt) {
    double coeff;
    coeff = - lB2(y, pt)/(z*z*Mperp2(pt));
    return gsl_sf_dilog(coeff);
}

inline double calculatePolyLog4(double z, double y, double pt) {
    double coeff;
    coeff = - LambdaBp2(y, pt)/(z*z*Mperp2(pt));
    return gsl_sf_dilog(coeff);
}

//calculating quenching weights
double PhatA(double z, double y, double pt) {
    double polylog1_result = calculatePolyLog1(z, y, pt);
    double polylog2_result = calculatePolyLog2(z, y, pt);
    double exponent = alphas * nc * (polylog1_result - polylog2_result) / (2 * M_PI);
    double logterms = 2 * log(1 + lA2(y, pt) / (z * z * Mperp2(pt))) / z - 2 * log(1 + LambdaAp2(y, pt) / (z * z * Mperp2(pt))) / z;
    double result = alphas * exp(exponent) * nc * logterms / (2 * M_PI);
    if (lA == lp){return 1;}
    else {return result;}
}

double PhatB(double z, double y, double pt) {
    double polylog3_result = calculatePolyLog3(z, y, pt);
    double polylog4_result = calculatePolyLog4(z, y, pt);
    double exponent = alphas * nc * (polylog3_result - polylog4_result) / (2 * M_PI);
    double logterms = 2 * log(1 + lB2(y, pt) / (z * z * Mperp2(pt))) / z - 2 * log(1 + LambdaBp2(y, pt) / (z * z * Mperp2(pt))) / z;
    double result = alphas * exp(exponent) * nc * logterms / (2 * M_PI);
    if (lB == lp){return 1;}
    else {return result;}
}

//pp-cross section parametrization
inline double f1(double pt) { return pow(p0 * p0 / (p0 * p0 + pt * pt), m); }

inline double f2(double y, double pt) {
    double arg;
    if (fabs(y) < ymax(pt)) arg = (1 - 2 * Mperp(pt) / rootsnn * cosh(y));
    else arg = 1e-30;
    return pow(arg, n);
}

inline double dsigdyd2pt(double y, double pt) {
    if (fabs(y) > ymax(pt)) return 1e-30;
    double f1Val = f1(pt);
    double f2Val = f2(y,pt);
    return f1Val * f2Val;
}

//calculation of momentum shift (for pA)
inline double shiftedPTpA(double pt, double dpta, double phiA) {
    return sqrt(pow(-dpta + cos(phiA) * pt, 2) + pow(sin(phiA) * pt, 2));
}

//calculation of momentum shift (for AB)
inline double shiftedPTAB(double pt, double dptb, double dpta, double phiB, double phiA) {
    double component1 = pt - dpta * cos(phiA) - dptb * cos(phiB);
    double component2 = dpta * sin(phiA) + dptb * sin(phiB);
    return sqrt(component1 * component1 + component2 * component2);
}

int scaledpAIntegrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata)
{
    Parameters* P = reinterpret_cast<Parameters*>(userdata);
    
    // Rescale ua and phiA to [0,1]
    // scaled values
    double ua = uMin + xx[0]*(P->uaMax-uMin);
    double phiA = xx[1]*2*M_PI;
    
    double dpta = dptA(P->y,P->pt);
    double shiftedPt = shiftedPTpA(P->pt, dpta, phiA);
    double phatAVal = PhatA(exp(exp(ua)) - 1, P->y, P->pt);
    double dsigVal = dsigdyd2pt(P->y + exp(ua), shiftedPt);
    
    ff[0] = exp(ua) * phatAVal * dsigVal;
    
    return 0;
}

int scaledABIntegrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata)
{
    Parameters* p = reinterpret_cast<Parameters*>(userdata);
    
    // Rescale ua and ub to [0, 1]
    // scaled values
    double ub = uMin + xx[0]*(p->ubMax-uMin);
    double ua = uMin + xx[1]*(p->uaMax-uMin);
    double phiB = xx[2]*2*M_PI;
    double phiA = xx[3]*2*M_PI;
    
    double dptb = dptB(p->y,p->pt);
    double dpta = dptA(-p->y,p->pt);
    double shiftedPt = shiftedPTAB(p->pt, dptb, dpta, phiB, phiA);
    double phatAVal = PhatA(exp(exp(ua)) - 1, -p->y, p->pt);
    double phatBVal = PhatB(exp(exp(ub)) - 1, p->y, p->pt);
    double dsigVal = dsigdyd2pt(p->y + exp(ub) - exp(ua), shiftedPt);
    ff[0] = exp(ub) * exp(ua) * phatBVal * phatAVal * dsigVal;
    
    return 0;
}

void pACrossSection(double y, double pt, double* res, double* err) {
    
    if (fabs(y) > ymax(pt)) {
        *res = 1e-30;
        *err = 1e-30;
        return;
    }
    
    Parameters P;
    P.y = y;
    P.pt = pt;
    P.uaMax = log(dymax(P.y, P.pt));
    
    // Integration parameters
    cubareal integral_result, error, prob;
    int nregions, neval, fail;
    
    Cuhre(
          NDIM2, NCOMP, scaledpAIntegrand, &P, NVEC, epsrel, epsabs, VERBOSE | LAST,
          MINEVAL, maxeval, KEY, nullptr, nullptr, &nregions, &neval, &fail,
          &integral_result, &error, &prob
          );
    
    if (fail==1) cout << "Error during integration." << endl;
    *res = (P.uaMax-uMin)*integral_result;
    *err = (P.uaMax-uMin)*error;
}

void ABCrossSection(double y, double pt, double* res, double* err) {
    
    if (fabs(y) > ymax(pt)) {
        *res = 1e-30;
        *err = 1e-30;
        return;
    }
    
    Parameters p;
    p.y = y;
    p.pt = pt;
    p.uaMax = log(dymax(-p.y, p.pt));
    p.ubMax = log(dymax(p.y, p.pt));
    
    // Integration parameters
    cubareal integral_result, error, prob;
    int nregions, neval, fail;
    
    Cuhre(
          NDIM4, NCOMP, scaledABIntegrand, &p, NVEC, epsrel, epsabs, VERBOSE | LAST,
          MINEVAL, maxeval, KEY, nullptr, nullptr, &nregions, &neval, &fail,
          &integral_result, &error, &prob
          );
    
    if (fail==1) cout << "Error during integration." << endl;
    *res = (p.ubMax-uMin)*(p.uaMax-uMin)*integral_result;
    *err = (p.ubMax-uMin)*(p.uaMax-uMin)*error;
}

int main()
{
    
    print_line(); // cosmetic
    auto starttime = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Quenching calculation started: " << ctime(&starttime);
    print_line(); // cosmetic
    
    readParametersFromFile("input/params.txt",1);
    processParameters();
    
    print_line(); // cosmetic
    
    create_output_directory();
    
    print_line(); // cosmetic
    string filename_1 = "output/pp-cross-section.tsv";
    ofstream output_file_1(filename_1);
    

    double resultpA, resultRpA, resultAB, resultRAB;
    double errorpA, errorAB,  errorRAB;
    
    switch (collisionType){
    case 0:
        {
            string filename_2 = "output/pA-cross-section.tsv";
            string filename_3 = "output/RpA.tsv";
            ofstream output_file_2(filename_2);
            ofstream output_file_3(filename_3);

            double resultpA, resultRpA;
            double errorpA, errorRpA;

            for (int i = 0; i < Ny; i++) {
                double y = y_min + i * dy;
                for (int j = 0; j < Npt; j++) {
                    double pt = ptmin + j * dpt;

                    // output pp cross section
                    double resultPP = dsigdyd2pt(y, pt);
                    output_file_1 << y << "\t" << pt << "\t" << resultPP << endl;

                    // compute pA cross section
                    pACrossSection(y, pt, &resultpA, &errorpA);
                    output_file_2 << y << "\t" << pt << "\t" << resultpA << "\t" << errorpA << endl;

                    // output RpA
                    resultRpA = resultpA / resultPP;
                    errorRpA = errorpA / resultpA;
                    printResult("RpA", y, pt, resultRpA, errorRpA);
                    output_file_3 << y << "\t" << pt << "\t" << resultRpA << "\t" << errorRpA << endl;
                }
            }

            output_file_2.close();
            output_file_3.close();
            break;
        }
        
        case 1:
        {
            string filename_4 = "output/AB-cross-section.tsv";
            string filename_5 = "output/RAB.tsv";
            ofstream output_file_4(filename_4);
            ofstream output_file_5(filename_5);

            double resultAB, resultRAB;
            double errorAB, errorRAB;

            for (int i = 0; i < Ny; i++) {
                double y = y_min + i * dy;
                for (int j = 0; j < Npt; j++) {
                    double pt = ptmin + j * dpt;

                    // output pp cross section
                    double resultPP = dsigdyd2pt(y, pt);
                    output_file_1 << y << "\t" << pt << "\t" << resultPP << endl;

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

            output_file_4.close();
            output_file_5.close();
            break;
        }
        default:
            cerr << "Invalid collision type. Exiting." << endl;
            return 1;
    }

    output_file_1.close();

    // print done!
    print_line();

    auto endtime = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << "Done: " << ctime(&endtime);
    print_line(); // cosmetic

    return 0;
}
