/*
 
 main.h
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#define NDIM2 2 //2-D integral
#define NDIM4 4 //4-D integral
#define NCOMP 1 
#define NVEC 1
#define VERBOSE 0
#define LAST 0
#define MINEVAL 0
#define KEY 9 // sets the integration rule
#define HBARC 0.197326938

using namespace std;

// lower limit for u integrations
const double uMin = -30.0;

// Define the desired absolute and relative error limits (global)
const double epsabs = 1e-12;
const double epsrel = 1e-12;
const int maxeval = 1e7;

struct Parameters
{
    double y, pt, uaMax, ubMax;
};

extern int collisionType, nc, Ny, Npt;
extern double y_min, y_max, ptmin, ptmax, dy, dpt;
extern double alphas,massQQ,lambdaQCD,qhat0,lp,lA,lB;
extern double massp,rootsnn,p0,m,n,beamRap,xA0,xB0;

