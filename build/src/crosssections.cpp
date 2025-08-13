/*
 
 crosssections.cpp
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <numeric> 

#include "cuba.h"
#include "paramreader.h"
#include "runningcoupling.h"
#include "crosssections.h"

#include <omp.h>

#include "main.h"
#include "glauber.h"

using namespace std;


//
// Global params; these default values are overridden by the params file
//
//Default Collision Type is pA, For AB, collisionType = 1 (overridden by params)

int collisionType = 0; 
int nc = 3;
double alphas = 0.5; // QCD coupling constant
double qhat0 = 0.075; // GeV^2/fm
double lp = 1.5; // in fm
double lA = 10.11; // in fm for Pb, https://arxiv.org/pdf/1304.0901
double lB = 10.11; // in fm for Pb
double massp = 0.938; // mass of proton in GeV, it can also be 1 GeV
double rootsnn = 5023; // collision energy, sqrt(s_NN)

//Default value of lambdaQCD 
// QCD scale (Lambda QCD) in GeV
// this is used in running coupling and phat calculation
double lambdaQCD = 0.308;
// double lambdaQCD = 0.25; 

// Computed from params file
double beamRap, xA0, xB0;

// lower limit for u integrations
const double uMin = -30.0;

// Define the desired absolute and relative error limits (global)
const double epsabs = 1e-12;
const double epsrel = 1e-12;
const int maxeval = 1e7;

// Parameters for loops in y and pt
int Ny = 10*2+1;
int Npt = 40*2+1;
double y_min = -5.0;
double y_max = 5.0;
double ptmin = 0.1;
double ptmax = 40.1;
double dy = (y_max - y_min) / (Ny-1);
double dpt = (ptmax - ptmin) / (Npt-1);

// numerical safety helpers
static inline double safeSqrt(double x) { return std::sqrt(x < 0.0 ? 0.0 : x); }
static inline double clampExpArg(double x, double lo, double hi){ return (x<lo?lo:(x>hi?hi:x)); }
constexpr double zMinFloor = 1e-8; // avoid z->0 where cancellations get noisy

// Physics Functions
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

// inline double dptA(double y, double pt) { return sqrt(lA2(y, pt) - lAp2(y, pt));}
// inline double dptB(double y, double pt) { return sqrt(lB2(y, pt) - lBp2(y, pt));}
inline double dptA(double y, double pt) { return safeSqrt(lA2(y, pt) - lAp2(y, pt)); }
inline double dptB(double y, double pt) { return safeSqrt(lB2(y, pt) - lBp2(y, pt)); }

inline double LambdaAp2(double y, double pt) { return max(lambdaQCD * lambdaQCD, lAp2(y, pt));}
inline double LambdaBp2(double y, double pt) { return max(lambdaQCD * lambdaQCD, lBp2(y, pt));}

inline double dymax(double y, double pt) {
    double log2 = log(2.0);
    double ymax_val = ymax(pt);
    double result = min(log2, ymax_val - y);
    return (result > 1e-12 ? result : 1e-12);
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
double PhatA(double z, double y, double pt, double alphas) {
    if (z < 1e-12) z = 1e-12;
    double polylog1_result = calculatePolyLog1(z, y, pt);
    double polylog2_result = calculatePolyLog2(z, y, pt);
    double exponent = alphas * nc * (polylog1_result - polylog2_result) / (2 * M_PI);
    exponent = clampExpArg(exponent, -700.0, 700.0);
    //stable logs
    const double M2 = Mperp2(pt);
    const double invz2M2 = 1.0 / (z*z*M2);
   //  double logterms = 2 * log(1 + lA2(y, pt) / (z * z * Mperp2(pt))) / z - 2 * log(1 + LambdaAp2(y, pt) / (z * z * Mperp2(pt))) / z;
    double logterms = 2.0 * std::log1p(lA2(y, pt) * invz2M2) / z
                    - 2.0 * std::log1p(LambdaAp2(y, pt) * invz2M2) / z;
    double result = alphas * std::exp(exponent) * nc * logterms / (2 * M_PI);
    
    // if (lA == lp){return 1;}
    if (!std::isfinite(result) || result < 0.0) result = 0.0;
    return result;
}

double PhatB(double z, double y, double pt, double alphas) {
    if (!(z > zMinFloor)) z = zMinFloor;
    if (std::fabs(lB - lp) < 1e-12) return 1.0;

    const double M2 = Mperp2(pt);
    const double invz2M2 = 1.0 / (z*z*M2);

    double poly3 = calculatePolyLog3(z, y, pt);
    double poly4 = calculatePolyLog4(z, y, pt);
    double expo = alphas * nc * (poly3 - poly4) / (2.0 * M_PI);
    expo = clampExpArg(expo, -700.0, 700.0);

    double logterms = 2.0 * std::log1p(lB2(y, pt) * invz2M2) / z
                    - 2.0 * std::log1p(LambdaBp2(y, pt) * invz2M2) / z;

    double res = (alphas * std::exp(expo) * nc * logterms) / (2.0 * M_PI);
    // if (lB == lp){return 1;}
    if (!std::isfinite(res) || res < 0.0) res = 0.0;
    return res;
}


//pp-cross section parametrization
inline double f1(double pt) { return pow(p0 * p0 / (p0 * p0 + pt * pt), m); }

inline double f2(double y, double pt) {
    double arg;
    if (fabs(y) < ymax(pt)) arg = (1 - 2 * Mperp(pt) / rootsnn * cosh(y));
    else arg = 1e-30;
    return pow(arg, n);
}

double dsigdyd2pt(double y, double pt) {
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
    double z = std::exp(std::exp(ua)) - 1.0;
    if (!(z > zMinFloor)) z = zMinFloor;                // avoid z→0 singularity
    double phatAVal = PhatA(z, P->y, P->pt, P->alphas_a);
    double dsigVal = dsigdyd2pt(P->y + exp(ua), shiftedPt);
    double val = std::exp(ua) * phatAVal * dsigVal;
    if (!std::isfinite(val) || val <= 0.0) val = 0.0;
    ff[0] = val;
    
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
    double zA = std::exp(std::exp(ua)) - 1.0;
    double zB = std::exp(std::exp(ub)) - 1.0;
    if (!(zA > zMinFloor)) zA = zMinFloor;
    if (!(zB > zMinFloor)) zB = zMinFloor;
    double phatAVal = PhatA(zA, -p->y, p->pt, p->alphas_a);
    double phatBVal = PhatB(zB, p->y, p->pt, p->alphas_b);
    double dsigVal = dsigdyd2pt(p->y + exp(ub) - exp(ua), shiftedPt);
    
    double val = std::exp(ua) * std::exp(ub) * phatAVal * phatBVal * dsigVal;
    if (!std::isfinite(val) || val <= 0.0) val = 0.0;
    ff[0] = val;
    
    return 0;
}

// Function to calculate the pA cross section
void pACrossSection(double y, double pt, double* res, double* err) {
    
    if (fabs(y) > ymax(pt)) {
        *res = 1e-30;
        *err = 1e-30;
        return;
    }
    
    Parameters P;
    P.y = y;
    P.pt = pt;
    P.uaMax = log(dymax(y, pt));
    P.alphas_a = (alphas != 0 ? alphas : runningCoupling(dptA(y, pt)));

    // Integration parameters
    cubareal integral_result, error, prob;
    int nregions, neval, fail;
    //INTEGRATION ROUTINE
    Cuhre(
          NDIM2, NCOMP, scaledpAIntegrand, &P, NVEC, epsrel, epsabs, VERBOSE | LAST,
          MINEVAL, maxeval, KEY, nullptr, nullptr, &nregions, &neval, &fail,
          &integral_result, &error, &prob
          );

    //ERROR HANDLING 
   if (fail!=0 || !std::isfinite(integral_result)) cout << ">>>> Error (pA) or NaN! <<<<< " << endl;
   *res = (P.uaMax-uMin)*integral_result;
   *err = (P.uaMax-uMin)*error;
   if(!std::isfinite(*res) || *res<0.0) *res=0.0;
}

// Function to calculate the AB cross section
void ABCrossSection(double y, double pt, double* res, double* err) {
    
    if (fabs(y) > ymax(pt)) {
        *res = 1e-30;
        *err = 1e-30;
        return;
    }
    
    Parameters p;
    p.y = y;
    p.pt = pt;
    p.uaMax = log(dymax(-y, pt));
    p.ubMax = log(dymax(y, pt));
    p.alphas_a = (alphas != 0 ? alphas : runningCoupling(dptA(-y, pt)));
    p.alphas_b = (alphas != 0 ? alphas : runningCoupling(dptB(y, pt)));

    // Integration parameters
    cubareal integral_result, error, prob;
    int nregions, neval, fail;
    //INTEGRATION ROUTINE
    Cuhre(
          NDIM4, NCOMP, scaledABIntegrand, &p, NVEC, epsrel, epsabs, VERBOSE | LAST,
          MINEVAL, maxeval, KEY, nullptr, nullptr, &nregions, &neval, &fail,
          &integral_result, &error, &prob
          );
      
 
     if (fail!=0 || !std::isfinite(integral_result)) { cout << ">>>> Error (AB) or NaN!" << endl; }

   *res = (p.ubMax-uMin)*(p.uaMax-uMin)*integral_result;
   *err = (p.ubMax-uMin)*(p.uaMax-uMin)*error;
   if(!std::isfinite(*res) || *res<0.0) *res=0.0;
}
