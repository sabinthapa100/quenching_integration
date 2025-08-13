/*
 
 runningcoupling.cpp
 
 Implements four-loop running coupling as a function of mu in GeV
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_zeta.h>

#include <main.h>

using namespace std;

#define Nf  3.
#define Zeta3   gsl_sf_zeta(3.)
#define b0  (11. - 2.*Nf/3.)
#define b1  (102. - 38.*Nf/3.)
#define b2  (2857./2. - (5033./18.)*Nf + (325./54.)*Nf*Nf)
#define b3  ((149753./6. + 3564.*Zeta3) - (1078361./162. + 6508.*Zeta3/27.)*Nf + (50065./162. + 6472.*Zeta3/81.)*Nf*Nf + 1093.*Nf*Nf*Nf/729.)
#define b4  0
#define alphasRef (0.326/(4*M_PI)) // taken from lattice calculation
#define muRef 1.5 // GeV
#define EPS 1e-4
#define STEP 1e-3
//#define lambdaQCD 0.308 // GeV
double runningCouplingApprox(double mu);
int func(double t, const double y[], double f[], void *params)
{
    f[0] = -2*(b0*y[0]*y[0] + b1*y[0]*y[0]*y[0] + b2*y[0]*y[0]*y[0]*y[0] + b3*y[0]*y[0]*y[0]*y[0]*y[0])/t;
    return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 1, 1);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, -2*(2*b0*y[0] + 3*b1*y[0]*y[0] + 4*b2*y[0]*y[0]*y[0] + 5*b3*y[0]*y[0]*y[0]*y[0])/t);
    dfdt[0] = 2*(b0*y[0]*y[0] + b1*y[0]*y[0]*y[0] + b2*y[0]*y[0]*y[0]*y[0] + b3*y[0]*y[0]*y[0]*y[0]*y[0])/t/t;
    return GSL_SUCCESS;
}

// approximate form
double runningCouplingApprox(double mu) {
    // Make sure mu is not below Lambda_QCD
    const double mu_floor = 1.01*lambdaQCD;
    mu = std::max(mu, mu_floor);

    // Asymptotic expansion variable
    double t = 2.0*std::log(mu/lambdaQCD);

    // Outside the domain of the asymptotic series? -> use 1-loop with a t floor
    const double T_MIN = 0.5;         // minimal "safe" t for using higher-loop series
    const double ALPHA_MAX = 1.0;     // freeze alpha_s to avoid unphysical blow-ups

    if (t < T_MIN || !std::isfinite(t)) {
        double a1 = 4.0*M_PI/(b0*std::max(t, 1e-12));   // 1-loop fallback
        return std::min(std::max(a1, 0.0), ALPHA_MAX);
    }

    // 4-loop series (your original code)
    const double l  = std::log(t);

    const double b02 = b0*b0,  b04 = b02*b02,  b06 = b02*b04,  b08 = b04*b04;
    const double b12 = b1*b1,  b13 = b1*b12,  b14 = b12*b12;
    const double t2 = t*t, t3 = t*t2, t4 = t2*t2;
    const double l2 = l*l, l3 = l*l2, l4 = l2*l2;

    double res = 1.0;
    res += -b1*l/b02/t;
    res += (b12*(l2 - l - 1) + b0*b2)/b04/t2;
    res += (b13*(-2*l3 + 5*l2 + 4*l - 1) - 6*b0*b1*b2*l + b02*b3)/(2.0*b06*t3);
    res += (18*b0*b2*b12*(2*l2 - l - 1) + b14*(6*l4 - 26*l3 - 9*l2 + 24*l + 7)
          - b02*b3*b1*(12*l + 1) + 2*b02*(5*b2*b2 + b0*b4))/(6.0*b08*t4);

    double alpha = 4.0*M_PI*res/(b0*t);

    // If the 4-loop series misbehaves, fall back to 1-loop and cap
    if (!std::isfinite(alpha) || alpha <= 0.0 || alpha > ALPHA_MAX) {
        double a1 = 4.0*M_PI/(b0*std::max(t, T_MIN));
        alpha = std::min(std::max(a1, 0.0), ALPHA_MAX);
    }
    return alpha;
}

// mu should be in GeV
double runningCoupling(double mu) {
    // Disable GSL's default abort-on-error once
    static bool handlerOff = (gsl_set_error_handler_off(), true);
    (void)handlerOff;

    // Cache last known-good alpha_s (4π * alphasRef ≈ 0.326 from lattice)
    static double lastAlpha = 4.0*M_PI*alphasRef;

    // Keep mu away from the Landau pole
    const double muMin = 1.01*lambdaQCD;
    if (!(mu > muMin) || !std::isfinite(mu)) mu = muMin;

    double t1 = muRef;
    double y[1] = { alphasRef };
    gsl_odeiv2_system sys = {func, jac, 1, NULL};

    // Ensure step sign matches target; GSL dislikes mismatched directions
    double step = STEP * ((mu >= muRef) ? +1.0 : -1.0);
    gsl_odeiv2_driver *d =
        gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                      std::fabs(step), EPS, 0.0);

    int status = gsl_odeiv2_driver_apply(d, &t1, mu, y);
    gsl_odeiv2_driver_free(d);

    if (status == GSL_SUCCESS && std::isfinite(y[0]) && y[0] > 0.0) {
        double alpha = 4.0*M_PI*y[0];
        if (std::isfinite(alpha) && alpha > 0.0) lastAlpha = alpha;
        std::cout << " Using 4-loop Running Coupling: " << alpha << std::endl;
        return alpha;
    }

    // Failure -> use previous good value (or the constant fallback if needed)
    if (!std::isfinite(lastAlpha) || lastAlpha <= 0.0)
        lastAlpha = 4.0*M_PI*alphasRef;

    std::cerr << "==> running coupling error, return value = "
              << status << " ; using cached/constant alpha_s = "
              << lastAlpha << std::endl;
    return lastAlpha;
}

