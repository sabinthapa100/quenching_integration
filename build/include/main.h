#define NDIM 4
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

extern int nc,Ny,Npt;
extern double y_min,y_max,ptmin,ptmax, dy, dpt;
extern double alphas,massQQ,lambdaQCD,qhat0,lp,lA,lB;
extern double massp,rootsnn,p0,m,n,beamRap,xA0,xB0;
