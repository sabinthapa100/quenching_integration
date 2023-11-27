/*
 
 runningcoupling.cpp
 
 Implements four-loop running coupling as a function of mu in GeV
 
 Copyright (c) Michael Strickland and Sabin Thapa
 
 GNU General Public License (GPLv3)
 See detailed text in license directory
 
 */

#include <iostream>
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

// mu should be in GeV
double runningCoupling(double mu) {
    double t1 = muRef;
    double y[1] = { alphasRef };
    gsl_odeiv2_system sys = {func, jac, 1, NULL};
    double step = STEP*(mu >= muRef ? +1 : -1);
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, step, EPS, 0.0);
    int status = gsl_odeiv2_driver_apply(d, &t1, mu, y);
    if (status != GSL_SUCCESS) cout << "==> running coupling error, return value = " << status << endl;
    gsl_odeiv2_driver_free (d);
    return 4*M_PI*y[0];
}

// approximate form
double runningCouplingApprox(double mu) {
    
    const double t = 2*log(mu/lambdaQCD);
    const double l = log(t);
    
    const double b02 = b0*b0;
    const double b04 = b02*b02;
    const double b06 = b02*b04;
    const double b08 = b04*b04;
    const double b12 = b1*b1;
    const double b13 = b1*b12;
    const double b14 = b12*b12;
    const double t2 = t*t;
    const double t3 = t*t2;
    const double t4 = t2*t2;
    const double l2 = l*l;
    const double l3 = l*l2;
    const double l4 = l2*l2;
    
    double res = 1;
    res += -b1*l/b02/t;
    res += (b12*(l2 - l -1)+b0*b2)/b04/t2;
    res += (b13*(-2*l3 + 5*l2 + 4*l -1) - 6*b0*b1*b2*l + b02*b3)/2/b06/t3;
    res += (18*b0*b2*b12*(2*l2 - l - 1) + b14*(6*l4 - 26*l3 - 9*l2 + 24*l + 7) - b02*b3*b1*(12*l+1)+2*b02*(5*b2*b2 + b0*b4))/6/b08/t4;
    return 4*M_PI*res/b0/t;
}
