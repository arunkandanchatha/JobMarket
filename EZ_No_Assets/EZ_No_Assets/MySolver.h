#pragma once

#include <cmath>
#include "stdafx.h"
#include "OLGModel.h"

#define TINY 1.0e-20
#define GLIMIT 100.0
#define GOLD  ((sqrt(5.) + 1.) / 2.)
#define R (GOLD-1)
#define C (1-R) 
#define MAX_ITERATIONS (1000)

typedef double (OLGModel::*OLGModelMemFn)(double x);

class MySolver
{
public:
	MySolver(OLGModel &p_model, OLGModelMemFn p_fn, double p_tol);
	~MySolver();

	double solve(double lowerLimit, double upperLimit);
	double value(double x);

private:
	const double TOLERANCE;

	OLGModelMemFn m_functionToSolve;
	OLGModel m_model;
	double goldenSearch(double a, double b, double c);
	const double EPSILON = 1.0e-10;

	/* From Numerical Recipes*/
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);
};

