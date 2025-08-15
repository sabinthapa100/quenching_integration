// build/include/glauber.h  (or your include path)

#ifndef __glauber_h__
#define __glauber_h__

#pragma once

// Geometry / sampling
double woodsSaxonDist(double r, double A);
void   sampleWoodsSaxon(int A, double *x, double *y);
int    MCcollisions(int A, double b, double *x, double *y);

// Thickness / profiles
double TA(double x, double y, double A);
double TAB(double x, double y, double A, double B, double b);
double TpA(double x, double y, double A, double b);
double nPartAB(double x, double y, double A, double B, double b);
double tProfileABpart(double x, double y, double A, double B, double b);
double tProfileABbin (double x, double y, double A, double B, double b);

// Sigma_NN control (mb)
void   set_sigmaNN_mb(double sNN_mb);

// Effective path lengths
double compute_LA_minbias_pA(int A, double rho0, double Lp);
double compute_LA_centrality_pA(double cmin, double cmax,
                                int A, double sigmaNN_mb,
                                double rho0, double Lp);

double compute_Npart_centrality_pA(double c0, double c1, int A, double sigmaNN_mb, double rho0, double lp);

#endif // __glauber_h__

