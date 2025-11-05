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
#include <algorithm>
#include <string>
#include "paramreader.h"
using namespace std;

#include "glauber.h"

double SIGMANN = 71; // default; can be overridden via set_sigmaNN_mb()
void set_sigmaNN_mb(double sNN_mb) { SIGMANN = sNN_mb; }

double SMALL = 1e-15;

struct ta_params { double x; double y; double A; };

// ---- Woods–Saxon -----------------------------------------------------------

double woodsSaxonDist(double r, double A) {
    // central density rho0 is defined in glauber.h (extern double rho0)
    const double Rn = 1.12*pow(A,1./3.) - 0.86*pow(A,-1./3.); // fm
    const double d  = 0.549;                                  // fm
    return rho0/(1.0 + exp((r - Rn)/d));
}

// ---------------------------------------------------------------------------
// MC helper used elsewhere (unchanged)
// ---------------------------------------------------------------------------

void sampleWoodsSaxon(int A, double *x, double *y) {
    const double M = 1.01*4.5310551374155095;
    const double rmax = 20.0;
    int m = 0;
    double r,v,f,u;
    while (m < A) {
        r = rmax * ((double) rand())/((double)RAND_MAX);
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

int MCcollisions(int A, double b, double *x, double *y) {
    double *x1,*y1,*x2,*y2,dn,dist;
    int *l1,*l2;
    x1 = new double[A]; y1 = new double[A]; l1 = new int[A];
    x2 = new double[A]; y2 = new double[A]; l2 = new int[A];
    sampleWoodsSaxon(A,x1,y1);
    sampleWoodsSaxon(A,x2,y2);
    for (int i=0; i<A; i++) { x1[i] -= b/2; x2[i] += b/2; }
    dn = sqrt(0.1*SIGMANN/M_PI);
    for (int i=0; i<A; i++) { l1[i]=0; l2[i]=0; }
    for (int i=0; i<A; i++)
      for (int j=0; j<A; j++) {
        dist = (x1[i]-x2[j])*(x1[i]-x2[j]) + (y1[i]-y2[j])*(y1[i]-y2[j]);
        dist = sqrt(dist);
        if (dist < dn) { l1[i] += 1; l2[j] += 1; }
      }
    int n = 0;
    for (int i=0; i<A; i++) {
      if (l1[i]>0) { x[n]=x1[i]; y[n]=y1[i]; n++; }
      if (l2[i]>0) { x[n]=x2[i]; y[n]=y2[i]; n++; }
    }
    delete[] x1; delete[] y1; delete[] l1;
    delete[] x2; delete[] y2; delete[] l2;
    return n;
}

// ---------------------------------------------------------------------------
// Thickness TA(x,y) = ∫dz ρ(√(x²+y²+z²))  (LUT + linear interpolation)
// ---------------------------------------------------------------------------

static double TAintegrand(double z, void * params) {
    ta_params *p = (ta_params*)params;
    const double x = p->x, y = p->y, A = p->A;
    return woodsSaxonDist( sqrt(x*x + y*y + z*z), A );
}

static double TA_line_integral(double r, double A) {
    double result = 0.0, error = 0.0;
    const size_t n = 200;
    ta_params int_params = {r, 0.0, A};
    gsl_function F; F.function = &TAintegrand; F.params = &int_params;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(n);
    gsl_integration_qags(&F, 0., 50., 0, 1.0e-6, n, w, &result, &error);
    gsl_integration_workspace_free(w);
    return 2.0*result; // by symmetry
}

namespace {
    struct TA_LUT {
        double A_cached   = -1.0;
        double rmax       = 50.0;
        double dr         = 0.02;
        std::vector<double> vals;        // TA(r_i)
    };
    static TA_LUT g_ta;
}

static inline void ensure_TA_table(double A) {
    if (A == g_ta.A_cached && !g_ta.vals.empty()) return;
    g_ta.A_cached = A;
    const int Nr = (int)std::floor(g_ta.rmax / g_ta.dr) + 1;
    g_ta.vals.assign(Nr, 0.0);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < Nr; ++i) {
        const double r = i * g_ta.dr;
        g_ta.vals[i] = TA_line_integral(r, A);
    }
}

double TA(double x, double y, double A) {
    ensure_TA_table(A);
    const double r = std::sqrt(x*x + y*y);
    const int Nr = (int)std::floor(g_ta.rmax / g_ta.dr) + 1;
    // hard clamp to zero beyond LUT domain to avoid spurious peripheral weight
    if (r >= g_ta.rmax) return 0.0;
    const double idx = r / g_ta.dr;
    int i = (int)std::floor(idx);
    if (i >= Nr-1) return 0.0;
    const double t = idx - i;
    return (1.0 - t)*g_ta.vals[i] + t*g_ta.vals[i+1];
}

// ---------------------------------------------------------------------------
// Proton profile (normalized ≈1 over d^2s; used in optical Tp⊗TA path)
// ---------------------------------------------------------------------------

double Tp(double x, double y) {
    // m≈1.85, Rp≈0.975 fm; numeric prefactor chosen to normalize ∫d^2s Tp ≈ 1
    return 0.400905*exp(-1.28022*pow(x*x+y*y,0.925));
}

// Convolution kernel: T_pA(b) = ∫ d^2s TA(s+b/2) Tp(s-b/2)
static double T_pA_conv(double b, int A) {
    const double xmin = b - 5.0, xmax = b + 5.0;
    const double ymin = -15.0   , ymax = 15.0;
    const int Nx = 160, Ny = 160;
    const double dx = (xmax - xmin) / Nx;
    const double dy = (ymax - ymin) / Ny;

    double sum = 0.0;
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2) reduction(+:sum) schedule(static)
    #endif
    for (int ix = 0; ix < Nx; ++ix)
        for (int iy = 0; iy < Ny; ++iy) {
            const double x = xmin + (ix + 0.5) * dx;
            const double y = ymin + (iy + 0.5) * dy;
            sum += TA(x + b/2.0, y, A) * Tp(x - b/2.0, y);
        }
    return sum * dx * dy; // fm^-2
}

double TAB(double x, double y, double A, double B, double b) {
    return TA(x+b/2.,y,A)*TA(x-b/2.,y,B);
}

double TpA(double x, double y, double A, double b) {
    return TA(x+b/2., y, A) * Tp(x-b/2., y);
}

// ---------------------------------------------------------------------------
// AA wounded nucleon profile (unchanged)
// ---------------------------------------------------------------------------

double nPartAB(double x, double y, double A, double B, double b) {
    // 1 fm^(-2) = 0.1 mb^(-1)
    double t, rval;
    t = SIGMANN*TA(x-b/2.,y,B)*0.1/B;
    if (t > SMALL)
        rval = TA(x+b/2.,y,A)*(1. - pow(1. - t,(double)B));
    else
        rval = TA(x+b/2.,y,A)*(B*t*(6 + (-1 + B)*t*(-3 + (-2 + B)*t)))/6.;
    t = SIGMANN*TA(x+b/2.,y,A)*0.1/A;
    if (t > SMALL)
        rval += TA(x-b/2.,y,B)*(1. - pow(1. - t,(double)A));
    else
        rval += TA(x-b/2.,y,B)*(A*t*(6 + (-1 + A)*t*(-3 + (-2 + A)*t)))/6.;
    return rval;
}

double tProfileABpart(double x, double y, double A, double B, double b) {
    return pow(nPartAB(x,y,A,B,b)/nPartAB(0.,0.,A,B,0.),0.25);
}

double tProfileABbin(double x, double y, double A, double B, double b) {
    return pow(TAB(x,y,A,B,b)/TAB(0.,0.,A,B,0.),0.25);
}

// ---------------------------------------------------------------------------
// pA building blocks
//   * Binomial/point-like (Arleo–Peigné) uses TA(b) only.
//   * Optical/smeared uses Tp⊗TA.
// ---------------------------------------------------------------------------

// --- Point-like proton (1304.0901, Appendix B) ---
static inline double lambda_pointlike(double b, int A, double sigmaNN_mb) {
    const double sigma = sigmaNN_mb * 0.1; // fm^2
    return sigma * TA(b, 0.0, A);
}
static inline double Pinel_pointlike(double b, int A, double sigmaNN_mb) {
    const double lam = lambda_pointlike(b,A,sigmaNN_mb);
    return (lam<=0.0 ? 0.0 : 1.0 - exp(std::max(-700.0, -lam)));
}

// --- Smeared proton (optical with Tp ⊗ TA) ---
static inline double lambda_smeared(double b, int A, double sigmaNN_mb) {
    const double sigma = sigmaNN_mb * 0.1; // fm^2
    return sigma * T_pA_conv(b, A);
}
static inline double Pinel_smeared(double b, int A, double sigmaNN_mb) {
    const double lam = lambda_smeared(b,A,sigmaNN_mb);
    return (lam<=0.0 ? 0.0 : 1.0 - exp(std::max(-700.0, -lam)));
}

// ---------------------------------------------------------------------------
// Centrality helpers (two flavors to keep models independent)
// ---------------------------------------------------------------------------

static void pA_bin_edges_pointlike(double c0, double c1, int A, double sigmaNN_mb,
                                   double &bmin, double &bmax)
{
    const double bMax = 20.0, db = bMax / 600.0;
    std::vector<double> cum(601, 0.0);
    double total = 0.0;
    for (int i = 0; i <= 600; ++i) {
        const double b = i * db;
        const double prob = Pinel_pointlike(b, A, sigmaNN_mb);
        total += prob * (2.0 * M_PI * b * db);
        cum[i] = total;
    }
    const double tmin = c0*total, tmax = c1*total;
    bmin = 0.0; bmax = bMax;
    for (int i=0;i<=600;++i){ if(cum[i]>=tmin){ bmin=i*db; break; } }
    for (int i=0;i<=600;++i){ if(cum[i]>=tmax){ bmax=i*db; break; } }
}

static void pA_bin_edges_smeared(double c0, double c1, int A, double sigmaNN_mb,
                                 double &bmin, double &bmax)
{
    const double bMax = 20.0, db = bMax / 600.0;
    std::vector<double> cum(601, 0.0);
    double total = 0.0;
    for (int i = 0; i <= 600; ++i) {
        const double b = i * db;
        const double prob = Pinel_smeared(b, A, sigmaNN_mb);
        total += prob * (2.0 * M_PI * b * db);
        cum[i] = total;
    }
    const double tmin = c0*total, tmax = c1*total;
    bmin = 0.0; bmax = bMax;
    for (int i=0;i<=600;++i){ if(cum[i]>=tmin){ bmin=i*db; break; } }
    for (int i=0;i<=600;++i){ if(cum[i]>=tmax){ bmax=i*db; break; } }
}

// ---------------------------------------------------------------------------
// Public: <N_part> in centrality bin — BINOMIAL (point-like) as in 1304.0901
// ---------------------------------------------------------------------------

double compute_Npart_centrality_pA(double c0, double c1, int A,
                                   double sigmaNN_mb, double /*rho0*/, double /*lp*/)
{
    double bmin=0.0, bmax=0.0;
    pA_bin_edges_pointlike(c0, c1, A, sigmaNN_mb, bmin, bmax);

    const int Nb = 600;
    const double db = (bmax - bmin) / Nb;

    long double num = 0.0L, den = 0.0L;
    for (int i = 0; i < Nb; ++i) {
        const double b = bmin + (i + 0.5) * db;
        const double lam   = lambda_pointlike(b, A, sigmaNN_mb);
        const double pinel = Pinel_pointlike(b, A, sigmaNN_mb);
        if (pinel<=0.0) continue;
        const double Npart_cond = 1.0 + lam/pinel; // 1 (proj) + <N_coll|inel>
        const double w = 2.0 * M_PI * b * db * pinel;
        num += (long double)Npart_cond * w;
        den += (long double)w;
    }
    if (den <= 0.0L) return 1.0;
    return (double)(num / den);
}

// ---------------------------------------------------------------------------
// Minimum-bias L_eff  (unchanged; used for checks)
// ---------------------------------------------------------------------------

double compute_LA_minbias_pA(int A, double rho0, double Lp)
{
    const double bMax = 20.0, db = bMax / 1000.0;
    auto TA_b = [&](double b){ return TA(b, 0.0, A); };

    double I = 0.0;
    for (int i=0; i<=1000; ++i) {
        const double b = i * db;
        const double w = (i==0 || i==1000) ? 0.5 : 1.0;
        const double t = TA_b(b);
        I += w * (2.0*M_PI) * b * (t*t);
    }
    I *= db;

    const double corr = (A>0 ? (A-1.0)/(A*1.0*A*1.0) : 0.0);
    const double LA = Lp + corr * I / std::max(1e-16, rho0);
    if (!std::isfinite(LA) || LA <= Lp) return Lp;
    return LA;
}

// ---------------------------------------------------------------------------
// Optical (Poisson) L_eff in a centrality bin — SMEARED (Tp ⊗ TA)
//   L_A = L_p + [⟨N(N−1)⟩ / (σ ρ0 ⟨N⟩)], with class weights ∝ P_inel(b).
//   Implemented as: Num = ∫ λ^2 d^2b, Den = ∫ λ d^2b (P cancels).
// ---------------------------------------------------------------------------

static double compute_LA_centrality_pA_optical(double cmin, double cmax,
                                               int A, double sigmaNN_mb,
                                               double rho0, double Lp)
{
    double bmin=0.0, bmax=0.0;
    pA_bin_edges_smeared(cmin, cmax, A, sigmaNN_mb, bmin, bmax);

    const int Nb = 600;
    const double db = (bmax - bmin) / Nb;
    const double sigma = sigmaNN_mb * 0.1; // fm^2

    long double Num = 0.0L, Den = 0.0L;
    for (int i = 0; i < Nb; ++i) {
        const double b = bmin + (i + 0.5) * db;
        const double lam = lambda_smeared(b, A, sigmaNN_mb);
        const double wA  = (2.0 * M_PI * b * db); // area-only weight (P cancels)
        Num += (long double)(lam*lam) * wA;
        Den += (long double)(lam)     * wA;
    }
    if (Den <= 0.0L) return Lp;
    const long double denom = (long double)(sigma) * (long double)(rho0) * Den;
    double LA = Lp + (double)(Num / std::max(1.0e-18L, denom));
    if (!std::isfinite(LA)) LA = Lp;
    return std::max(Lp, LA);
}

// ---------------------------------------------------------------------------
// Binomial L_eff (Eq. B.9 of 1304.0901) — POINT-LIKE proton (exact A&P)
// ---------------------------------------------------------------------------

static double compute_LA_centrality_pA_binomial(double cmin, double cmax,
                                                int A, double sigmaNN_mb,
                                                double rho0, double Lp,
                                                bool allow_fallback_optical)
{
    const double sigma = sigmaNN_mb * 0.1;   // fm^2
    const double bMax = 20.0;
    const int    Nb   = 1000;
    const double db   = bMax / Nb;

    // --- Step 1: compute σ_N integrated over b (point-like p → p(b)=σTA/A)
    std::vector<double> logC(A+1, 0.0);
    const double lA1 = std::lgammal(A+1);
    for (int N=0; N<=A; ++N)
        logC[N] = lA1 - std::lgammal(N+1) - std::lgammal(A-N+1);

    std::vector<double> sigmaN(A+1, 0.0);
    for (int i=0; i<=Nb; ++i) {
        const double b   = i*db;
        const double w   = (i==0 || i==Nb) ? 0.5 : 1.0; // trapezoid
        const double jac = w * (2.0*M_PI) * b * db;     // d^2b
        const double TA_b = TA(b, 0.0, A);              // fm^-2
        double p = sigma * TA_b / (double)A;            // dimensionless
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

    double sig_inel = 0.0;
    for (int N=1; N<=A; ++N) sig_inel += sigmaN[N];
    if (sig_inel <= 0.0) {
        const double LAopt = allow_fallback_optical
            ? compute_LA_centrality_pA_optical(cmin,cmax,A,sigmaNN_mb,rho0,Lp)
            : Lp;
        return std::max(Lp, LAopt);
    }

    // --- Step 2: map centrality to [N_low, N_high]
    std::vector<double> F(A+2, 0.0);
    double tail = 0.0;
    for (int N=A; N>=1; --N) { tail += sigmaN[N]; F[N] = tail / sig_inel; }
    F[A+1] = 0.0;

    auto N_from_cent = [&](double c)->int {
        if (c <= 0.0) return A+1;
        if (c >= 1.0) return 1;
        for (int N=A; N>=1; --N) if (F[N] >= c && F[N+1] < c) return N;
        return 1;
    };

    int N_low  = N_from_cent(cmax);
    int N_high = (cmin <= 0.0 ? A : N_from_cent(cmin)-1);
    if (N_low < 1) N_low = 1;
    if (N_high > A) N_high = A;
    if (N_low > N_high) {
        const double LAopt = allow_fallback_optical
            ? compute_LA_centrality_pA_optical(cmin,cmax,A,sigmaNN_mb,rho0,Lp)
            : Lp;
        return std::max(Lp, LAopt);
    }

    // --- Step 3: Eq. (B.9)
    long double Num = 0.0L, Den = 0.0L;
    for (int N=N_low; N<=N_high; ++N) {
        Num += (long double)N*(N-1) * sigmaN[N];
        Den += (long double)N         * sigmaN[N];
    }

    double LA_binom = Lp;
    if (Den > 0.0L) {
        const long double denom = (long double)sigma * (long double)rho0 * Den;
        LA_binom = Lp + (double)(Num / std::max(1.0e-18L, denom));
        if (!std::isfinite(LA_binom)) LA_binom = Lp;
    }

    if (allow_fallback_optical && (LA_binom <= Lp + 1.0e-6)) {
        const double LA_opt = compute_LA_centrality_pA_optical(cmin,cmax,A,sigmaNN_mb,rho0,Lp);
        return std::max(Lp, LA_opt);
    }
    return std::max(Lp, LA_binom);
}

// ---------------------------------------------------------------------------
// Unified selector (binomial = A&P point-like; optical = Tp ⊗ TA)
// ---------------------------------------------------------------------------

double compute_LA_centrality_pA_flex(double cmin, double cmax,
                                     int A, double sigmaNN_mb,
                                     double rho0, double Lp,
                                     const std::string& method,
                                     bool allow_fallback_optical)
{
    if (method == "optical") {
        return compute_LA_centrality_pA_optical(cmin,cmax,A,sigmaNN_mb,rho0,Lp);
    }
    return compute_LA_centrality_pA_binomial(cmin,cmax,A,sigmaNN_mb,rho0,Lp,
                                             allow_fallback_optical);
}

// Backward-compatible wrapper (keeps existing callsites working)
double compute_LA_centrality_pA(double cmin, double cmax,
                                int A, double sigmaNN_mb,
                                double rho0, double Lp)
{
    //"optical", "binomial"
    return compute_LA_centrality_pA_flex(cmin, cmax, A, sigmaNN_mb, rho0, Lp,
                                         "binomial", /*allow_fallback_optical=*/true);
}

// ===========================================================================
// OPTIONAL HELPERS (unchanged public surface)
// ===========================================================================

// ⟨N_coll⟩ in a centrality bin (point-like by design; conditional-on-inelastic)
double compute_Ncoll_centrality_pA(double c0, double c1, int A, double sigmaNN_mb)
{
    double bmin=0.0, bmax=0.0;
    pA_bin_edges_pointlike(c0, c1, A, sigmaNN_mb, bmin, bmax);

    const int Nb = 600;
    const double db = (bmax - bmin) / Nb;

    long double num = 0.0L, den = 0.0L;
    for (int i = 0; i < Nb; ++i) {
        const double b = bmin + (i + 0.5) * db;
        const double lam   = lambda_pointlike(b, A, sigmaNN_mb);
        const double pinel = Pinel_pointlike(b, A, sigmaNN_mb);
        if (pinel <= 0.0) continue;
        const double Kcond = lam / pinel;                // ⟨N_coll | inel, b⟩
        const double w     = 2.0 * M_PI * b * db * pinel; // class weight ∝ P(b)
        num += (long double)Kcond * w;                   // = lam * area
        den += (long double)w;                           // = P * area
    }
    if (den <= 0.0L) return 1.0;
    return (double)(num/den);
}

// Map centrality bin to [N_min, N_max] using the same binomial spectrum.
static void pA_bin_class_from_centrality(double cmin, double cmax,
                                         int A, double sigmaNN_mb,
                                         int &Nmin, int &Nmax, double &Nmean)
{
    const double sigma = sigmaNN_mb * 0.1;   // fm^2
    const double bMax = 20.0;
    const int    Nb   = 1000;
    const double db   = bMax / Nb;

    std::vector<double> logC(A+1, 0.0);
    const double lA1 = std::lgammal(A+1);
    for (int N=0; N<=A; ++N)
        logC[N] = lA1 - std::lgammal(N+1) - std::lgammal(A-N+1);

    std::vector<double> sigmaN(A+1, 0.0);
    for (int i=0; i<=Nb; ++i) {
        const double b   = i*db;
        const double w   = (i==0 || i==Nb) ? 0.5 : 1.0; // trapezoid
        const double jac = w * (2.0*M_PI) * b * db;     // d^2b
        const double TA_b = TA(b, 0.0, A);              // fm^-2
        double p = sigma * TA_b / (double)A;
        if (p <= 0.0) continue;
        if (p >= 1.0) p = 1.0 - 1e-12;
        const double lp = std::log(p), lq = std::log(1.0 - p);
        for (int N=1; N<=A; ++N) {
            const double lnw = logC[N] + N*lp + (A-N)*lq;
            const double wN  = std::exp(lnw);
            if (wN>0.0) sigmaN[N] += jac * wN;
        }
    }
    double sig_inel = 0.0;
    for (int N=1; N<=A; ++N) sig_inel += sigmaN[N];
    if (sig_inel <= 0.0) { Nmin=1; Nmax=1; Nmean=1.0; return; }

    std::vector<double> F(A+2, 0.0);
    double tail = 0.0;
    for (int N=A; N>=1; --N) { tail += sigmaN[N]; F[N] = tail / sig_inel; }
    F[A+1] = 0.0;

    auto N_from_cent = [&](double c)->int {
        if (c <= 0.0) return A+1;
        if (c >= 1.0) return 1;
        for (int N=A; N>=1; --N) if (F[N] >= c && F[N+1] < c) return N;
        return 1;
    };

    Nmin  = N_from_cent(cmax);
    Nmax  = (cmin <= 0.0 ? A : N_from_cent(cmin)-1);
    if (Nmin < 1) Nmin = 1;
    if (Nmax > A) Nmax = A;
    if (Nmin > Nmax) { Nmin=1; Nmax=1; Nmean=1.0; return; }

    long double num = 0.0L, den = 0.0L;
    for (int N=Nmin; N<=Nmax; ++N) { num += (long double)N * sigmaN[N]; den += (long double)sigmaN[N]; }
    Nmean = (den>0.0L) ? (double)(num/den) : 1.0;
}

// Pretty printer (unchanged API)
void print_pA_centrality_table(const std::vector<double>& edges_percent,
                               int A, double sigmaNN_mb, double rho0, double Lp,
                               const std::string& Leff_method /*"binomial" or "optical"*/)
{
    if (edges_percent.size() < 2) return;
    cout << "---------------------------------------------------------------\n";
    cout << "Centrality  |  N_part(min..max)  |  N_coll(min..max)  |  <N_coll>  |  L_eff [fm]\n";
    cout << "---------------------------------------------------------------\n";
    for (size_t i=0; i+1<edges_percent.size(); ++i) {
        const double c0 = edges_percent[i  ]/100.0;
        const double c1 = edges_percent[i+1]/100.0;

        int Nmin=1, Nmax=1; double Nmean=1.0;
        pA_bin_class_from_centrality(c0, c1, A, sigmaNN_mb, Nmin, Nmax, Nmean);

        const double Ncoll_mean = compute_Ncoll_centrality_pA(c0, c1, A, sigmaNN_mb);
        const double Leff = compute_LA_centrality_pA_flex(c0, c1, A, sigmaNN_mb, rho0, Lp,
                                                          Leff_method, /*fallback=*/true);

        // In pA: N_coll == N_part^A; total N_part = N_part^A + 1 (projectile proton)
        const int Npart_min_total = Nmin + 1;
        const int Npart_max_total = Nmax + 1;

        cout.setf(std::ios::fixed); cout.precision(2);
        cout << " " << edges_percent[i] << "-" << edges_percent[i+1] << "%  |  "
             << Npart_min_total << " , " << Npart_max_total << "        |  "
             << Nmin << " , " << Nmax << "            |  "
             << Ncoll_mean << "    |  "
             << Leff << "\n";
    }
    cout << "---------------------------------------------------------------\n";
}

