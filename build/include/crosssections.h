#ifndef __crosssections_h__
#define __crosssections_h__

#define NDIM2 2 //2-D integral
#define NDIM4 4 //4-D integral
#define NCOMP 1 
#define NVEC 1
#define VERBOSE 0
#define LAST 0
#define MINEVAL 0
#define KEY 9 // sets the integration rule
#define HBARC 0.197326938

#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "cuba.h"

using namespace std;

extern int collisionType, nc, Ny, Npt;
extern double y_min, y_max, ptmin, ptmax, dy, dpt;
extern double alphas, qhat0,lp,lA,lB;
extern double massp,rootsnn, beamRap,xA0,xB0;

struct Parameters
{
    double y, pt, uaMax, ubMax, alphas_a, alphas_b;
};

extern double alphas,lambdaQCD, qhat0, lp, lA, lB, massp, rootsnn;

extern double beamRap, xA0, xB0;

extern int Ny, Npt;
extern double y_min, y_max;
extern double ptmin, ptmax;
extern double dy, dpt;

inline double Mperp2(double pt);
inline double Mperp(double pt);
inline double ymax(double pt);
inline double qhat(double x);
inline double xA2(double y, double pt);
inline double xB2(double y, double pt);
inline double myXA(double y, double pt);
inline double myXB(double y, double pt);
inline double lA2(double y, double pt);
inline double lB2(double y, double pt);
inline double lAp2(double y, double pt);
inline double lBp2(double y, double pt);
inline double dptA(double y, double pt);
inline double dptB(double y, double pt);
inline double LambdaAp2(double y, double pt);
inline double LambdaBp2(double y, double pt);
inline double dymax(double y, double pt);
inline double calculatePolyLog1(double z, double y, double pt);
inline double calculatePolyLog2(double z, double y, double pt);
inline double calculatePolyLog3(double z, double y, double pt);
inline double calculatePolyLog4(double z, double y, double pt);
double PhatA(double z, double y, double pt, double alphas);
double PhatB(double z, double y, double pt, double alphas);
inline double f1(double pt);
inline double f2(double y, double pt);
double dsigdyd2pt(double y, double pt);
inline double shiftedPTpA(double pt, double dpta, double phiA);
inline double shiftedPTAB(double pt, double dptb, double dpta, double phiB, double phiA);
int scaledpAIntegrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata);
int scaledABIntegrand(const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata);
void averageNeighbors(double y, double pt, double* res, double* err, int num_neighbors, void (*crossSectionFunc)(double, double, double*, double*));
void pACrossSection(double y, double pt, double* res, double* err);
void ABCrossSection(double y, double pt, double* res, double* err);

#endif /**__crosssections_h__**/

