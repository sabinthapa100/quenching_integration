/*

   glauber.cpp

   Copyright (c) Michael Strickland & Sabin Thapa

   GNU General Public License (GPLv3)
*/

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <vector>
using namespace std;

#include "glauber.h"

double SIGMANN = 62.0; // default; will be overridden by set_sigmaNN_mb()
void set_sigmaNN_mb(double sNN_mb) { SIGMANN = sNN_mb; }

double SMALL = 1e-15;

struct ta_params { double x; double y; double A;};

double woodsSaxonDist(double r, double A) {
	double n0 = 0.17; // central density in fm^(-3)
	double Rn = 1.12*pow(A,1./3.) - 0.86*pow(A,-1./3.); // radius in fm 
	double d = 0.54; // thickness in fm
	return n0/(1+exp((r-Rn)/d));
}

// generates A samples corresponding to Woods-Saxon nucleus with A nucleons
// results are loaded into x and y arrays (z is not used in current context)
// note: assumes that the the random number generator has already been seeded
void sampleWoodsSaxon(int A, double *x, double *y) {
	const double M = 1.01*4.5310551374155095;
	const double rmax = 20;
	int m = 0;
	double r,v,f,u;
	while (m<A) {
	  r = rmax*((double) rand())/((double)RAND_MAX);
	  u = ((double) rand())/((double)RAND_MAX);
	  if (u < r*r*woodsSaxonDist(r,A)/M) {
	    v = 2*(((double) rand())/((double)RAND_MAX)-0.5);
	    f = 2*M_PI*((double) rand())/((double)RAND_MAX);
	    x[m] = r*sqrt(1-v*v)*cos(f);
	    y[m] = r*sqrt(1-v*v)*sin(f);
	    m++;
	  }
	}
}

// finds colliding nucleons using hard sphere model
// returns number of wounded nucleons found and loads their
// positions into the x and y arrays passed down to it
int MCcollisions(int A, double b, double *x, double *y) {
	double *x1,*y1,*x2,*y2,dn,dist;	
	int *l1,*l2;
	x1 = new double[A];
	y1 = new double[A];
	l1 = new int[A];
	x2 = new double[A];
	y2 = new double[A];
	l2 = new int[A];
	sampleWoodsSaxon(A,x1,y1);
	sampleWoodsSaxon(A,x2,y2);
	for (int i=0; i<A; i++) {
		x1[i] -= b/2;
		x2[i] += b/2;
	}
	dn = sqrt(0.1*SIGMANN/M_PI);
	for (int i=0; i<A; i++) {
		l1[i] = 0;
		l2[i] = 0;
	}
	for (int i=0; i<A; i++)
	  for (int j=0; j<A; j++) {
	    dist = pow(x1[i]-x2[j],2);
	    dist += pow(y1[i]-y2[j],2);
	    dist = sqrt(dist);
	    if (dist<dn) {
		l1[i] += 1;
		l2[j] += 1;
	    }
	  }
	int n = 0;
	for (int i=0; i<A; i++) {
	  if (l1[i]>0) {
	    x[n] = x1[i];
	    y[n] = y1[i];
	    n++;
	  }
	  if (l2[i]>0) {
	    x[n] = x2[i];
	    y[n] = y2[i];
	    n++;
	  }
	}
	delete[] x1;
	delete[] y1;
	delete[] l1;
	delete[] x2;
	delete[] y2;
	delete[] l2;
	return n;
}

double TAintegrand(double z, void * params) {
	struct ta_params * my_params = (struct ta_params *)params;
        double x = my_params->x;
        double y = my_params->y;
        double A = my_params->A;
	return woodsSaxonDist(sqrt(x*x+y*y+z*z),A);
}


// thickness function
double TA(double x, double y, double A) {
	double result, error;
	double n = 100;
        struct ta_params int_params = {x,y,A};
	gsl_function F;
        F.function = &TAintegrand;
        F.params = &int_params;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(n);
	gsl_integration_qags(&F, 0., 50., 0, 1.0e-6, n, w, &result, &error);
	gsl_integration_workspace_free(w);
	return 2*result;
}

double Tp(double x, double y) {
 	//return 0.479598*exp(-1.51107*pow(x*x+y*y,0.925));
	return 0.400905*exp(-1.28022*pow(x*x+y*y,0.925));
}

double TAB(double x, double y, double A, double B, double b) {
	return TA(x+b/2,y,A)*TA(x-b/2,y,B);
}

double TpA(double x, double y, double A, double b) {
	return TA(x+b,y,A)*Tp(x,y);
}

double nPartAB(double x, double y, double A, double B, double b) {
	// 1 fm^(-2) = 0.1 mb^(-1)
	double t,rval;
	t = SIGMANN*TA(x-b/2,y,B)*0.1/B;
	if (t > SMALL) 
		rval = TA(x+b/2,y,A)*(1. - pow(1. - t,(double)B));
	else
		rval = TA(x+b/2,y,A)*(B*t*(6 + (-1 + B)*t*(-3 + (-2 + B)*t)))/6.;
	t = SIGMANN*TA(x+b/2,y,A)*0.1/A;
	if (t > SMALL)
		rval += TA(x-b/2,y,B)*(1. - pow(1. - t,(double)A));
	else
		rval += TA(x-b/2,y,B)*(A*t*(6 + (-1 + A)*t*(-3 + (-2 + A)*t)))/6.;
	return rval;
}

// wounded nucleon profile
double tProfileABpart(double x, double y, double A, double B, double b) {
        return pow(nPartAB(x,y,A,B,b)/nPartAB(0.,0.,A,B,0.),0.25);
}

// binary collision profile
double tProfileABbin(double x, double y, double A, double B, double b) {
        return pow(TAB(x,y,A,B,b)/TAB(0.,0.,A,B,0.),0.25);
}
// Min Bias Collision Case L_eff = L_p + 
double compute_LA_minbias_pA(int A, double rho0, double Lp)
{
    // ∫ d^2b T_A^2(b) = 2π ∫_0^{bMax} b db [T_A(b,0,A)]^2
    const double bMax = 20.0;   // fm (safe for Pb/Au tails)
    const int    Nb   = 2000;
    const double db   = bMax / Nb;

    auto TA_b = [&](double b){ return TA(b, 0.0, A); };

    double I = 0.0;
    for (int i=0; i<=Nb; ++i) {
        double b = i * db;
        double w = (i==0 || i==Nb) ? 0.5 : 1.0; // trapezoid
        double t = TA_b(b);
        I += w * (2.0*M_PI) * b * (t*t);
    }
    I *= db;

    double corr = (A>0 ? (A-1.0)/(A*1.0*A*1.0) : 0.0);
    double LA = Lp + corr * I / std::max(1e-16, rho0);
    if (!std::isfinite(LA) || LA <= Lp) return Lp;
    return LA;
}

// Different centrality classes 
// p+A centrality slicing per Eq. (B.9) of arXiv:1304.0901
// Centrality is provided as percent fractions cmin,cmax in [0,1]:
// e.g. 0-0.20 (most-central 0–20%) etc.
double compute_LA_centrality_pA(double cmin, double cmax,
                                int A, double sigmaNN_mb,
                                double rho0, double Lp)
{
    // units: 1 mb = 0.1 fm^2
    const double sigma = sigmaNN_mb * 0.1;   // fm^2

    // --- Step 1: compute σ_N = ∫ d^2b Binom(A,N) p(b)^N (1-p)^(A-N), N>=1
    const double bMax = 20.0;   // fm
    const int    Nb   = 4000;
    const double db   = bMax / Nb;

    // precompute log binomial coeffs
    std::vector<double> logC(A+1, 0.0);
    const double lA1 = std::lgammal(A+1);
    for (int N=0; N<=A; ++N)
        logC[N] = lA1 - std::lgammal(N+1) - std::lgammal(A-N+1);

    std::vector<double> sigmaN(A+1, 0.0); // σ_N for N=0..A, we use N>=1
    for (int i=0; i<=Nb; ++i) {
        const double b   = i*db;
        const double w   = (i==0 || i==Nb) ? 0.5 : 1.0;      // trapezoid
        const double jac = w * (2.0*M_PI) * b * db;          // d^2b
        const double TA_b = TA(b, 0.0, A);                   // fm^-2
        double p = sigma * TA_b / (double)A;                 // dimensionless
        if (p <= 0.0) continue;
        if (p >= 1.0) p = 1.0 - 1e-12;
        const double lp = std::log(p);
        const double lq = std::log(1.0 - p);

        for (int N=1; N<=A; ++N) {
            const double lnw = logC[N] + N*lp + (A-N)*lq;
            const double wN  = std::exp(lnw);
            if (wN>0.0) sigmaN[N] += jac * wN;
        }
    }

    // total inelastic σ = sum_{N>=1} σ_N
    double sig_inel = 0.0;
    for (int N=1; N<=A; ++N) sig_inel += sigmaN[N];
    if (sig_inel <= 0.0) return Lp;

    // --- Step 2: map percent edges to integer N ranges using the CENTRAL cumulative
    // F(N) = fraction of σ with multiplicity >= N
    std::vector<double> F(A+2, 0.0);
    double tail = 0.0;
    for (int N=A; N>=1; --N) { tail += sigmaN[N]; F[N] = tail / sig_inel; }
    F[A+1] = 0.0;

    auto N_from_cent = [&](double c)->int {
        // find N with F(N) >= c and F(N+1) < c
        if (c <= 0.0) return A+1; // sentinel (so the first bin gets N_high=A)
        if (c >= 1.0) return 1;
        for (int N=A; N>=1; --N) if (F[N] >= c && F[N+1] < c) return N;
        return 1;
    };

    // bin [cmin,cmax]: N ∈ [N_low..N_high]
    int N_low  = N_from_cent(cmax);                 // more peripheral edge
    int N_high = (cmin <= 0.0 ? A : N_from_cent(cmin)-1); // more central edge
    if (N_low < 1) N_low = 1;
    if (N_high > A) N_high = A;
    if (N_low > N_high) return Lp; // empty bin

    // --- Step 3: evaluate Eq. (B.9)
    long double Num = 0.0L, Den = 0.0L;
    for (int N=N_low; N<=N_high; ++N) {
        Num += (long double)N*(N-1) * sigmaN[N];
        Den += (long double)N         * sigmaN[N];
    }
    if (Den <= 0.0L) return Lp;

    const long double denom = (long double)sigma * (long double)rho0 * Den;
    double LA = Lp + (double)(Num / std::max(1e-30L, denom));
    if (!std::isfinite(LA) || LA <= Lp) return Lp;
    return LA;
}

