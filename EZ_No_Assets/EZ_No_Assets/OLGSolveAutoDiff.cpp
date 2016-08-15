#include "OLGSolveAutoDiff.h"



OLGSolveAutoDiff::OLGSolveAutoDiff(int gens, std::vector<double> ys):m_gens(gens)
{
	m_Y.resize(ys.size());
	m_Y = ys;
}


OLGSolveAutoDiff::~OLGSolveAutoDiff()
{
}

double OLGSolveAutoDiff::solveProblem(std::vector<double> x) {

	//setup input
	int numStates = x.size();
	std::vector<adouble> m_thetas(numStates);
	for (unsigned int i = 0; i < x.size(); i++) {
		m_thetas[i] = x[i];
	}

	//now, solve for wages
	{
		std::vector<adouble> E_vals(numStates);
		std::vector<adouble> U_vals(numStates);
		std::vector<adouble> W_vals(numStates);
		std::vector<adouble> wages(numStates);
		for (unsigned int i = m_gens - 1; i < m_gens; i--) {
			std::vector<adouble> latestE = E_vals;
			std::vector<adouble> latestU = U_vals;
			std::vector<adouble> latestW = W_vals;
			std::vector<adouble> latestWages = wages;
			for (int j = 0; j < numStates; j++) {
#if 0
				std::cout << i << ":" << j << std::endl;
#endif
				if (i == m_gens - 1) {
					E_vals[j] = (1 - D_BETA)*pow(D_b, D_RHO);
					U_vals[j] = (1 - D_BETA)*pow(D_b, D_RHO);
					W_vals[j] = 0;
					wages[j] = D_b;
					continue;
				}
				double prodInState = m_Y[j];

				//find solution
				{
					adouble lwbnd = latestWages[j] - D_b;
					double upbnd = prodInState - D_b;

					if (lwbnd == upbnd) {
						wages[j] = latestWages[j];
						U_vals[j] = latestU[j];
						E_vals[j] = latestE[j];
						W_vals[j] = latestW[j];
						continue;
					}

					//otherwise do brent to solve
					adouble arg;
					adouble c;
					adouble d;
					adouble e;
					adouble eps;
					adouble fu;
					adouble fv;
					adouble fw;
					adouble fx;
					adouble midpoint;
					adouble p;
					adouble q;
					adouble r;
					double tol;
					adouble tol1;
					adouble tol2;
					adouble u;
					adouble v;
					adouble w;
					adouble x;

					adouble a=lwbnd;
					adouble b=upbnd;
					int status;
					adouble value;

					for (int convCounter = 0; convCounter < 50; convCounter++) {
						if (convCounter > 0) {
							value = f(arg);
						}
						status = convCounter;

						//double local_min_rc(double &a, double &b, int &status, double value);
						//code from brent.cpp (Jburkardt)
						//****************************************************************************80
						//
						//  Purpose:
						//
						//    LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
						//
						//  Discussion:
						//
						//    This routine seeks an approximation to the point where a function
						//    F attains a minimum on the interval (A,B).
						//
						//    The method used is a combination of golden section search and
						//    successive parabolic interpolation.  Convergence is never much
						//    slower than that for a Fibonacci search.  If F has a continuous
						//    second derivative which is positive at the minimum (which is not
						//    at A or B), then convergence is superlinear, and usually of the
						//    order of about 1.324...
						//
						//    The routine is a revised version of the Brent local minimization
						//    algorithm, using reverse communication.
						//
						//    It is worth stating explicitly that this routine will NOT be
						//    able to detect a minimizer that occurs at either initial endpoint
						//    A or B.  If this is a concern to the user, then the user must
						//    either ensure that the initial interval is larger, or to check
						//    the function value at the returned minimizer against the values
						//    at either endpoint.
						//
						//  Licensing:
						//
						//    This code is distributed under the GNU LGPL license.
						//
						//  Modified:
						//
						//    17 July 2011
						//
						//  Author:
						//
						//    John Burkardt
						//
						//  Reference:
						//
						//    Richard Brent,
						//    Algorithms for Minimization Without Derivatives,
						//    Dover, 2002,
						//    ISBN: 0-486-41998-3,
						//    LC: QA402.5.B74.
						//
						//    David Kahaner, Cleve Moler, Steven Nash,
						//    Numerical Methods and Software,
						//    Prentice Hall, 1989,
						//    ISBN: 0-13-627258-4,
						//    LC: TA345.K34.
						//
						//  Parameters
						//
						//    Input/output, double &A, &B.  On input, the left and right
						//    endpoints of the initial interval.  On output, the lower and upper
						//    bounds for an interval containing the minimizer.  It is required
						//    that A < B.
						//
						//    Input/output, int &STATUS, used to communicate between
						//    the user and the routine.  The user only sets STATUS to zero on the first
						//    call, to indicate that this is a startup call.  The routine returns STATUS
						//    positive to request that the function be evaluated at ARG, or returns
						//    STATUS as 0, to indicate that the iteration is complete and that
						//    ARG is the estimated minimizer.
						//
						//    Input, double VALUE, the function value at ARG, as requested
						//    by the routine on the previous call.
						//
						//    Output, double LOCAL_MIN_RC, the currently considered point.
						//    On return with STATUS positive, the user is requested to evaluate the
						//    function at this point, and return the value in VALUE.  On return with
						//    STATUS zero, this is the routine's estimate for the function minimizer.
						//
						//  Local parameters:
						//
						//    C is the squared inverse of the golden ratio.
						//
						//    EPS is the square root of the relative machine precision.
						//
						{

							//
							//  STATUS (INPUT) = 0, startup.
							//
							if (status == 0)
							{
								if (b <= a)
								{
									std::cout << "\n";
									std::cout << "LOCAL_MIN_RC - Fatal error!\n";
									std::cout << "  A < B is required, but\n";
									std::cout << "  A = " << a << "\n";
									std::cout << "  B = " << b << "\n";
									status = -1;
									exit(-1);
								}
								c = 0.5 * (3.0 - sqrt(5.0));

								eps = 1e-15;
								tol = 1e-9;

								v = a + c * (b - a);
								w = v;
								x = v;
								e = 0.0;

								status = 1;
								arg = x;

								continue;
							}
							//
							//  STATUS (INPUT) = 1, return with initial function value of FX.
							//
							else if (status == 1)
							{
								fx = value;
								fv = fx;
								fw = fx;
							}
							//
							//  STATUS (INPUT) = 2 or more, update the data.
							//
							else if (2 <= status)
							{
								fu = value;

								if (fu <= fx)
								{
									if (x <= u)
									{
										a = x;
									}
									else
									{
										b = x;
									}
									v = w;
									fv = fw;
									w = x;
									fw = fx;
									x = u;
									fx = fu;
								}
								else
								{
									if (u < x)
									{
										a = u;
									}
									else
									{
										b = u;
									}

									if (fu <= fw || w == x)
									{
										v = w;
										fv = fw;
										w = u;
										fw = fu;
									}
									else if (fu <= fv || v == x || v == w)
									{
										v = u;
										fv = fu;
									}
								}
							}
							//
							//  Take the next step.
							//
							midpoint = 0.5 * (a + b);
							tol1 = eps * ABS(x) + tol / 3.0;
							tol2 = 2.0 * tol1;
							//
							//  If the stopping criterion is satisfied, we can exit.
							//
							if (ABS(x - midpoint) <= (tol2 - 0.5 * (b - a)))
							{
								status = 0;
								break;
							}
							//
							//  Is golden-section necessary?
							//
							if (ABS(e) <= tol1)
							{
								if (midpoint <= x)
								{
									e = a - x;
								}
								else
								{
									e = b - x;
								}
								d = c * e;
							}
							//
							//  Consider fitting a parabola.
							//
							else
							{
								r = (x - w) * (fx - fv);
								q = (x - v) * (fx - fw);
								p = (x - v) * q - (x - w) * r;
								q = 2.0 * (q - r);
								if (0.0 < q)
								{
									p = -p;
								}
								q = ABS(q);
								r = e;
								e = d;
								//
								//  Choose a golden-section step if the parabola is not advised.
								//
								if (
									(ABS(0.5 * q * r) <= ABS(p)) ||
									(p <= q * (a - x)) ||
									(q * (b - x) <= p))
								{
									if (midpoint <= x)
									{
										e = a - x;
									}
									else
									{
										e = b - x;
									}
									d = c * e;
								}
								//
								//  Choose a parabolic interpolation step.
								//
								else
								{
									d = p / q;
									u = x + d;

									if ((u - a) < tol2)
									{
										d = tol1 * (((midpoint - x)<0)?-1:1);
									}

									if ((b - u) < tol2)
									{
										d = tol1 * (((midpoint - x)<0) ? -1 : 1);
									}
								}
							}
							//
							//  F must not be evaluated too close to X.
							//
							if (tol1 <= ABS(d))
							{
								u = x + d;
							}
							if (ABS(d) < tol1)
							{
								u = x + tol1 * ((d<0) ? -1 : 1);
							}
							//
							//  Request value of F(U).
							//
							arg = u;
							status = status + 1;

							continue;
						}
					}
				}
			}
		}
	}

	//finally, return distance
	double retVal = 0;
	for (int i = 0; i < x.size(); i++) {
		retVal += pow(D_C / D_BETA - m_f->calculatedF(m_thetas[i]) / m_thetas[i] * expectedW(i), 2);
	}
	return sqrt(retVal);
}
