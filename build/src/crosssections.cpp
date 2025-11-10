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
#include <fstream>
#include "cuba.h"
#include "paramreader.h"
#include "runningcoupling.h"
#include "crosssections.h"
#include <iomanip>
#include <omp.h>

#include "main.h"
#include "glauber.h"

using namespace std;

//
// Global params; these default values are overridden by the params file
//
int collisionType = 0; 
int nc = 3;
double alphas  = 0.5;    // QCD coupling constant
double qhat0   = 0.075;  // GeV^2/fm
double lp      = 1.5;    // fm
double lA      = 10.11;  // fm (Pb default), from https://arxiv.org/pdf/1304.0901
double lB      = 10.11;  // fm (Pb)
double massp   = 0.938;  // GeV
double rootsnn = 5023;   // GeV

// QCD scale (Lambda_QCD) in GeV (used in running coupling and Phat)
double lambdaQCD = 0.308;

// Computed from params file
double beamRap, xA0, xB0;

// lower limit for u integrations
const double uMin = -30.0;

// CUBA error limits
const double epsabs  = 1e-12;
const double epsrel  = 1e-12;
const int    maxeval = 1e7;

// (grids kept for completeness)
int Ny = 10*2+1;
int Npt = 40*2+1;
double y_min = -5.0;
double y_max =  5.0;
double ptmin = 0.1;
double ptmax = 40.1;
double dy  = (y_max - y_min) / (Ny-1);
double dpt = (ptmax - ptmin) / (Npt-1);

// ------------------------------------------------------------------
// numeric safety & Helper
// ------------------------------------------------------------------
static inline double safeSqrt(double x) { return std::sqrt(x < 0.0 ? 0.0 : x); }
static inline double clampExpArg(double x, double lo, double hi){ return (x<lo?lo:(x>hi?hi:x)); }
constexpr double zMinFloor = 1e-8; // avoid z→0 singularity

// ------------------------------------------------------------------
// Debugger
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// Single-CSV diagnostics (compact): ΔpT scales, running αs, 1 quench sample
// Writes: output/<outTag>diagnostics.csv
// Columns:
//   y,pt,mu_A,alphaA,muB,alphaB,dyAmax,dyBmax,zA_med,PhA_med,zB_med,PhB_med
// Logs only every LOG_STRIDE_Y and LOG_STRIDE_PT grid point to keep files small.
// ------------------------------------------------------------------
void debug_dump_scales(double y, double pt) {
    static std::ofstream f;
    static std::string last_prefix;
    static bool header_written = false;

    // --- reduce rows: log every Kth grid point
    static const int LOG_STRIDE_Y  = 2;  // change to 3/4 if you want even fewer rows
    static const int LOG_STRIDE_PT = 2;
    const int iy = (int) llround((y  - y_min )/dy );
    const int ip = (int) llround((pt - ptmin)/dpt);
    if (iy < 0 || ip < 0) return;
    if ((iy % LOG_STRIDE_Y) != 0 || (ip % LOG_STRIDE_PT) != 0) return;

    const std::string filepath = "output/" + outTag + "diagnostics.csv";

    // (Re)open file if first time or outTag changed (new centrality/bin)
    if (!f.is_open() || last_prefix != outTag) {
        if (f.is_open()) f.close();
        last_prefix = outTag;
        header_written = false;
        f.open(filepath.c_str(), std::ios::out | std::ios::trunc);
        if (!f) { std::cerr << "Cannot open " << filepath << " for writing.\n"; return; }
    }

    if (!header_written) {
        f << "y,pt,"
          << "muA,dpT_A,dpT_B,alphaA,muB,alphaB,"
          << "dyAmax,dyBmax,"
          << "zA_med,PhA_med,zB_med,PhB_med\n";
        header_written = true;
    }

    // Scales & couplings actually used
    // const double mu_A = dptA(y, pt);
    // const double mu_B = dptB(y, pt);
    double mu_A = Mperp (pt);
    double mu_B = Mperp (pt);
    const double aA  = (alphas != 0 ? alphas : runningCoupling(mu_A));
    const double aB  = (alphas != 0 ? alphas : runningCoupling(mu_B));

    // Rapidity windows
    const double dyAmax = dymax(+y, pt);
    const double dyBmax = dymax(-y, pt);

    // Median quenching sample (δy = 0.5 * δy_max)
    const double zA = std::exp(0.5*dyAmax) - 1.0;
    const double zB = std::exp(0.5*dyBmax) - 1.0;
    const double PhA = PhatA(zA, y, pt, aA);
    const double PhB = PhatB(zB, y, pt, aB);

    #pragma omp critical
    {
        f << std::setprecision(10)
          << y  << ',' << pt << ','
          << mu_A << ',' << dptA(y, pt) << ',' << dptB(y, pt)<< ',' << aA << ',' << mu_B << ',' << aB << ','
          << dyAmax << ',' << dyBmax << ','
          << zA << ',' << PhA << ',' << zB << ',' << PhB << '\n';
    }
}

// ------------------------------------------------------------------
// physics helpers
// ------------------------------------------------------------------
inline double Mperp2(double pt) { return pt*pt + massQQ*massQQ; }
inline double Mperp (double pt) { return std::sqrt(Mperp2(pt)); }

inline double ymax(double pt) { return std::log(rootsnn / Mperp(pt)); }

// transport coefficient qhat(x)
double qhat(double x) { return qhat0 * std::pow(1.0e-2 / x, 0.3); }

// Bjorken-x on A and B sides (COM rapidity y)
// *** IMPORTANT SIGN CONVENTION ***
//   x_A ~ e^{-y}, x_B ~ e^{-y}.  Sign flip is whenever needed;
inline double xA2(double y, double pt) { return Mperp(pt) / rootsnn * std::exp(-y); }
inline double xB2(double y, double pt) { return Mperp(pt) / rootsnn * std::exp(-y); }

inline double myXA(double y, double pt) { return std::min(xA0, xA2(y, pt)); }
inline double myXB(double y, double pt) { return std::min(xB0, xB2(y, pt)); }

// broadening scales
inline double lA2 (double y, double pt) { return qhat(myXA(y, pt)) * lA; }
inline double lB2 (double y, double pt) { return qhat(myXB(y, pt)) * lB; }
inline double lAp2(double y, double pt) { return qhat(myXA(y, pt)) * lp; }
inline double lBp2(double y, double pt) { return qhat(myXB(y, pt)) * lp; }

// ΔpT (safe)
inline double dptA(double y, double pt) { return safeSqrt(lA2(y, pt) - lAp2(y, pt)); }
inline double dptB(double y, double pt) { return safeSqrt(lB2(y, pt) - lBp2(y, pt)); }

// Λ_p^2 (GeV^2)
inline double LambdaAp2(double y, double pt) { return std::max(lambdaQCD*lambdaQCD, lAp2(y, pt)); }
inline double LambdaBp2(double y, double pt) { return std::max(lambdaQCD*lambdaQCD, lBp2(y, pt)); }

// rapidity window for δy (floor to avoid log(≤0))
// For AB: uaMax = log(dymax(-y,pt)), ubMax = log(dymax(+y,pt))  (see below).
inline double dymax(double y, double pt) {
    const double r = std::min(std::log(2.0), ymax(pt) - y);
    return (r > 1e-12 ? r : 1e-12);           // original guard
}

// ------------------------------------------------------------------
// Quenching weights (Arleo–Peigné form) with stable numerics
// Must be non-static to match crosssections.h prototypes.
// ------------------------------------------------------------------
double PhatA(double z, double y, double pt, double a) {
    if (!(z > zMinFloor)) z = zMinFloor;

    const double M2      = Mperp2(pt);
    const double invz2M2 = 1.0 / (z*z*M2);
    const double l2      = lA2(y, pt);
    const double Lp2     = LambdaAp2(y, pt);

    // Exact AP step function: if no extra path length on A, no induced radiation
    if (l2 <= Lp2) return 0.0;

    double expo = a*nc*(gsl_sf_dilog(-l2*invz2M2) - gsl_sf_dilog(-Lp2*invz2M2)) / (2.0*M_PI);
    expo = clampExpArg(expo, -700.0, 700.0);

    const double logterms = 2.0 * ( std::log1p(l2*invz2M2) - std::log1p(Lp2*invz2M2) ) / z;

    const double res = (a * std::exp(expo) * nc * logterms) / (2.0*M_PI);
    return (std::isfinite(res) && res > 0.0) ? res : 0.0;
}

double PhatB(double z, double y, double pt, double a) {
    if (!(z > zMinFloor)) z = zMinFloor;

    const double M2      = Mperp2(pt);
    const double invz2M2 = 1.0 / (z*z*M2);
    const double l2      = lB2(y, pt);
    const double Lp2     = LambdaBp2(y, pt);

    if (l2 <= Lp2) return 0.0;

    double expo = a*nc*(gsl_sf_dilog(-l2*invz2M2) - gsl_sf_dilog(-Lp2*invz2M2)) / (2.0*M_PI);
    expo = clampExpArg(expo, -700.0, 700.0);

    const double logterms = 2.0 * ( std::log1p(l2*invz2M2) - std::log1p(Lp2*invz2M2) ) / z;

    const double res = (a * std::exp(expo) * nc * logterms) / (2.0*M_PI);
    return (std::isfinite(res) && res > 0.0) ? res : 0.0;
}

// ------------------------------------------------------------------
// pp cross section shape
// ------------------------------------------------------------------
inline double f1(double pt) {
    const double denom = (p0*p0 + pt*pt);
    return std::pow(p0*p0 / denom, m);
}
inline double f2(double y, double pt) {
    if (std::fabs(y) >= ymax(pt)) return 1e-30;
    double arg = 1.0 - 2.0 * Mperp(pt) / rootsnn * std::cosh(y);
    if (arg <= 0.0) arg = 1e-30;              // original guard
    return std::pow(arg, n);
}

double dsigdyd2pt(double y, double pt) {
    if (std::fabs(y) > ymax(pt)) return 1e-30;
    return f1(pt) * f2(y, pt);
}

// ------------------------------------------------------------------
// momentum shifts
// ------------------------------------------------------------------
inline double shiftedPTpA(double pt, double dpta, double phiA) {
    const double c = std::cos(phiA), s = std::sin(phiA);
    const double x = -dpta + c*pt;
    return std::sqrt(x*x + (s*pt)*(s*pt));
}
inline double shiftedPTAB(double pt, double dptb, double dpta, double phiB, double phiA) {
    const double cA = std::cos(phiA), sA = std::sin(phiA);
    const double cB = std::cos(phiB), sB = std::sin(phiB);
    const double comp1 = pt - dpta*cA - dptb*cB;
    const double comp2 =       dpta*sA + dptb*sB;
    return std::sqrt(comp1*comp1 + comp2*comp2);
}

// ------------------------------------------------------------------
// pA integrand (2D): x = (ua, phiA)
// ------------------------------------------------------------------
int scaledpAIntegrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata)
{
    Parameters* P = reinterpret_cast<Parameters*>(userdata);
    
    // Rescale ua and phiA to [0,1]
    // scaled values
    double ua = uMin + xx[0]*(P->uaMax-uMin);
    double phiA = xx[1]*2*M_PI;
    
    double dpta = dptA(P->y,P->pt);
    double shiftedPt = shiftedPTpA(P->pt, dpta, phiA);
    double phatAVal = PhatA(exp(exp(ua)) - 1, P->y, P->pt, P->alphas_a);
    double dsigVal = dsigdyd2pt(P->y + exp(ua), shiftedPt);
    
    ff[0] = exp(ua) * phatAVal * dsigVal;
    
    return 0;
}

// ------------------------------------------------------------------
// AB integrand (4D): x = (ub, ua, phiB, phiA)
// y + δy_B - δy_A, with δy_{A,B} = exp(u_{a,b})
// Limits: uaMax=log(dymax(-y,pt)), ubMax=log(dymax(+y,pt))
// ------------------------------------------------------------------
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
    double phatAVal = PhatA(exp(exp(ua)) - 1, -p->y, p->pt, p->alphas_a); // + in A, - in B
    double phatBVal = PhatB(exp(exp(ub)) - 1, p->y, p->pt, p->alphas_b);
    double dsigVal = dsigdyd2pt(p->y + exp(ub) - exp(ua), shiftedPt);
    ff[0] = exp(ub) * exp(ua) * phatBVal * phatAVal * dsigVal;
    
    return 0;
}

// ------------------------------------------------------------------
// pA and AB drivers
// ------------------------------------------------------------------
void pACrossSection(double y, double pt, double* res, double* err) {
    if (std::fabs(y) > ymax(pt)) { *res = 1e-30; *err = 1e-30; return; }
    
    // δ–limit: if no extra medium (l_A^2 ≤ Λ_p^2), pA → pp exactly
    if (lA2(y, pt) <= LambdaAp2(y, pt)) {
        *res = dsigdyd2pt(y, pt); *err = 0.0; return;
    }

    Parameters P;
    P.y = y;
    P.pt = pt;
    P.uaMax    = std::log(dymax( y, pt));                      
    // P.alphas_a = (alphas != 0 ? alphas : runningCoupling(dptA(y, pt)));
    // double mu_A = dptA(y, pt);
    double mu_A = std::max( Mperp(pt), 1.0 );
    P.alphas_a = (alphas != 0 ? alphas : runningCoupling(mu_A));
    if (P.uaMax < uMin) P.uaMax = uMin;
    
    debug_dump_scales(y, pt);
    
    cubareal integral_result, error, prob;
    int nregions, neval, fail;

    Cuhre(
        NDIM2, NCOMP, scaledpAIntegrand, &P, NVEC, epsrel, epsabs, VERBOSE | LAST,
        MINEVAL, maxeval, KEY, nullptr, nullptr, &nregions, &neval, &fail,
        &integral_result, &error, &prob
    );
    if (fail!=0 || !std::isfinite(integral_result)) std::cout << ">>>> Error (pA) or NaN! <<<<<\n";
    *res = (P.uaMax - uMin) * integral_result;
    *err = (P.uaMax - uMin) * error;
    if (!std::isfinite(*res) || *res < 0.0) *res = 0.0;
}

void ABCrossSection(double y, double pt, double* res, double* err) {
    if (std::fabs(y) > ymax(pt)) { *res = 1e-30; *err = 1e-30; return; }

    // δ–limit: if no extra medium (l_A^2 ≤ Λ_p^2), pA → pp exactly
    if (lA2(y, pt) <= LambdaAp2(y, pt)) {
        *res = dsigdyd2pt(y, pt); *err = 0.0; return;
    }

    Parameters p;
    p.y  = y;
    p.pt = pt;
    p.uaMax    = std::log(dymax(-y, pt));                      // δy_A^max(-y)
    p.ubMax    = std::log(dymax( y, pt));                      // δy_B^max(+y)
    // p.alphas_a = (alphas != 0 ? alphas : runningCoupling(dptA(y, pt)));
    // p.alphas_b = (alphas != 0 ? alphas : runningCoupling(dptB(y, pt)));
    double mu_A = Mperp (pt);
    double mu_B = Mperp (pt);
    p.alphas_a = (alphas != 0 ? alphas : runningCoupling(mu_A));
    p.alphas_b = (alphas != 0 ? alphas : runningCoupling(mu_B));
    
    if (p.uaMax < uMin) p.uaMax = uMin;
    if (p.ubMax < uMin) p.ubMax = uMin;
    
    cubareal integral_result, error, prob;
    int nregions, neval, fail;

    Cuhre(
        NDIM4, NCOMP, scaledABIntegrand, &p, NVEC, epsrel, epsabs, VERBOSE | LAST,
        MINEVAL, maxeval, KEY, nullptr, nullptr, &nregions, &neval, &fail,
        &integral_result, &error, &prob
    );

    if (fail!=0 || !std::isfinite(integral_result)) std::cout << ">>>> Error (AB) or NaN!\n";

    *res = (p.ubMax - uMin) * (p.uaMax - uMin) * integral_result;
    *err = (p.ubMax - uMin) * (p.uaMax - uMin) * error;

    if (!std::isfinite(*res) || *res < 0.0) *res = 0.0;
}


