#include "MySolver.h"

MySolver::MySolver(FunctorBase &p_functionToSolve, double p_tol)
:TOLERANCE(p_tol), m_functionToSolve(&p_functionToSolve)
{
}

MySolver::~MySolver()
{
}

double MySolver::solve(double lowerLimit, double upperLimit) {
	
	m_min = lowerLimit;
	m_max = upperLimit;
	double ax = lowerLimit + (upperLimit - lowerLimit) / 3;
	double bx = lowerLimit + (upperLimit - lowerLimit) / 2;
	double cx, fa, fb, fc;

	mnbrak(&ax, &bx, &cx, &fa, &fb, &fc);
	return goldenSearch(ax, bx, cx);
}

double MySolver::value(double x) {
	return (*m_functionToSolve)(x);
}

double MySolver::goldenSearch(double ax, double bx, double cx) {
	double f1, f2, x0, x1, x2, x3;
	x0 = ax;
	x3 = cx;
	if (ABS(cx - bx) > ABS(bx - ax)) {
		x1 = bx;
		x2 = bx + C*(cx - bx);
	}
	else {
		x2 = bx;
		x1 = bx - C*(bx - ax);
	}
	f1 = (*m_functionToSolve)(x1);
	f2 = (*m_functionToSolve)(x2);
	int iterationCount = 0;
	double absVal = ABS(x3 - x0);
	double tolVal = TOLERANCE*(ABS(x1) + ABS(x2));
	while ((iterationCount++ < MAX_ITERATIONS) && (absVal > tolVal)) {
		//std::cout << iterationCount << ": f(" << x1 << ")=" << f1 << ",   f(" << x2 << ")=" << f2 << std::endl;
		if (f2 < f1) {
			SHFT3(x0, x1, x2, R*x1 + C*x3);
			SHFT2(f1, f2, (*m_functionToSolve)(x2));
		}
		else {
			SHFT3(x3, x2, x1, R*x2 + C*x0);
			SHFT2(f2, f1, (*m_functionToSolve)(x1));
		}
		absVal = ABS(x3 - x0);
		tolVal = TOLERANCE*(ABS(x1) + ABS(x2));
	}
	//std::cout << tolVal << std::endl;
	if (f1 < f2) {
		return x1;
	}
	else {
		return x2;
	}
}

void MySolver::mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc)
{
	double ulim, u, r, q, fu, dum;
	*fa = (*m_functionToSolve)(*ax);
	*fb = (*m_functionToSolve)(*bx);
	if (*fb > *fa) {
		SHFT3(dum, *ax, *bx, dum);
		SHFT3(dum, *fb, *fa, dum);
	}
	*cx = (*bx) + GOLD*(*bx - *ax);
	if (*cx < m_min) {
		*cx = m_min;
	}
	*fc = (*m_functionToSolve)(*cx);
	while (*fb > *fc) {
		r = (*bx - *ax)*(*fb - *fc);
		q = (*bx - *cx)*(*fb - *fa);
		u = (*bx) - ((*bx - *cx)*q - (*bx - *ax)*r) /
			(2.0*SIGN(MAX(fabs(q - r), TINY), q - r));
		ulim = (*bx) + GLIMIT*(*cx - *bx);
		if ((*bx - u)*(u - *cx) > 0.0) {
			fu = (*m_functionToSolve)(u);
			if (fu < *fc) {
				*ax = (*bx);
				*bx = u;
				*fa = (*fb);
				*fb = fu;
				return;
			}
			else if (fu > *fb) {
				*cx = u;
				*fc = fu;
				return;
			}
			u = (*cx) + GOLD*(*cx - *bx);
			fu = (*m_functionToSolve)(u);
		}
		else if ((*cx - u)*(u - ulim) > 0.0) {
			fu = (*m_functionToSolve)(u);
			if (fu < *fc) {
				SHFT3(*bx, *cx, u, *cx + GOLD*(*cx - *bx));
				SHFT3(*fb, *fc, fu, (*m_functionToSolve)(u));
			}
		}
		else if ((u - ulim)*(ulim - *cx) >= 0.0) {
			u = ulim;
			fu = (*m_functionToSolve)(u);
		}
		else {
			u = (*cx) + GOLD*(*cx - *bx);
			fu = (*m_functionToSolve)(u);
		}
		SHFT3(*ax, *bx, *cx, u);
		SHFT3(*fa, *fb, *fc, fu);
	}
}