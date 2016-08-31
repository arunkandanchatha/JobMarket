#pragma once
#include "stdafx.h"
#include "adept.h"
using adept::adouble;

#include <vector>

class OLGSolveAutoDiff
{
public:
	OLGSolveAutoDiff(int gens, std::vector<double> ys, MatrixXd conditionalProbs, double parameter,
		/*std::vector<double> bargaining,*/ double s, MatrixXd &wages);
	OLGSolveAutoDiff(OLGSolveAutoDiff &orig);
	~OLGSolveAutoDiff();

	double solveProblem(std::vector<double>& input, std::vector<double>& bargaining);
	double solveProblem(std::vector<double>& input, std::vector<double>&grad, std::vector<double>& bargaining);

	void setBargaining(std::vector<double>& b);

private:
	const int m_gens;
	std::vector<double> m_Y;
	std::vector<double> m_bargaining;
	MatrixXd* m_wgs;
	MatrixXd m_conditionalProbs;

	adept::Stack stack_;

	//FIX THESE TWO PARAMETERS (HOW???)
	adouble m_parameter;
	adouble m_s;
#if 0
public:
	//this is ugly but has to be here because of silly template/header conflicts
	template <typename T> bool operator()(T const * const* parameters, T* residual) const {
		//setup input
		int numStates = m_Y.size();
		std::vector<T> x(numStates);
		for (int i = 0; i < numStates; i++) {
			x[i] = parameters[0][i];
		}
		std::vector<std::vector<T>> W_vals;
		std::vector<std::vector<T>> wages;
		W_vals.resize(m_gens);
		wages.resize(m_gens);
		for (int i = 0; i < m_gens; i++) {
			W_vals[i].resize(numStates);
			wages[i].resize(numStates);
		}

		//now, solve for wages
		{
			std::vector<T> E_vals(numStates);
			std::vector<T> U_vals(numStates);
			for (unsigned int i = m_gens - 1; i < m_gens; i--) {
				std::vector<T> latestE(numStates);
				latestE = E_vals;
				std::vector<T> latestU(numStates);
				latestU = U_vals;
				std::vector<T> latestW(numStates);
				std::vector<T> latestWages(numStates);
				if (i < m_gens - 1) {
					latestW = W_vals[i + 1];
					latestWages = wages[i + 1];
				}
				for (int j = 0; j < numStates; j++) {
#if 0
					std::cout << i << ":" << j << std::endl;
#endif
					if (i == m_gens - 1) {
						E_vals[j] = T(1 - D_BETA)*T(pow(D_b, D_RHO));
						U_vals[j] = T(1 - D_BETA)*T(pow(D_b, D_RHO));
						W_vals[i][j] = T(0);
						wages[i][j] = T(D_b);
						continue;
					}
					T prodInState = T(m_Y[j]);

					//find solution
					{
						T lwbnd = latestWages[j] - T(D_b);
						T upbnd = prodInState - T(D_b);

						if (lwbnd == upbnd) {
							wages[i][j] = latestWages[j];
							U_vals[j] = latestU[j];
							E_vals[j] = latestE[j];
							W_vals[i][j] = latestW[j];
							continue;
						}

						//otherwise do brent to solve
						T arg;
						T c;
						T d;
						T e;
						T eps;
						T fu;
						T fv;
						T fw;
						T fx;
						T midpoint;
						T p;
						T q;
						T r;
						double tol;
						T tol1;
						T tol2;
						T u;
						T v;
						T w;
						T x;

						T a = lwbnd;
						T b = upbnd;
						int status;
						T value;

						for (int convCounter = 0; convCounter < 100; convCounter++) {
							if (convCounter > 0) {
								T calcE;
								{
									T total = 0;

									VectorXd nextPDF = m_conditionalProbs.row(j);
									for (int ii = MAX(0, j - MAX_SHOCKS_PER_MONTH); ii < MIN(numStates, j + MAX_SHOCKS_PER_MONTH + 1); ii++) {
										total += T(nextPDF(ii))*T(pow(latestE[ii] - T(m_s)*(latestE[ii] - latestU[ii]), T(D_RHO)));
									}
									total *= T(D_BETA);
									calcE = T(pow(T(1 - D_BETA)*T(pow(T(D_b) + arg, T(D_RHO))) + total, T(1.0) / T(D_RHO)));
									E_vals[j] = calcE;
								}
								T calcU;
								{
									T total = 0;

									VectorXd nextPDF = m_conditionalProbs.row(j);
									for (int ii = MAX(0, j - MAX_SHOCKS_PER_MONTH); ii < MIN(numStates, j + MAX_SHOCKS_PER_MONTH + 1); ii++) {
										T nextProb = T(nextPDF(ii));
										T calcF = T(D_MY_ALPHA*D_MY_MU)*x[ii] / T(pow(T(1) + T(pow(T(D_MY_MU)*x[ii], T(m_parameter))), T(1.0 / m_parameter)));
										T nextVal = T(pow(latestU[ii] + calcF*(latestE[ii] - latestU[ii]), T(D_RHO)));
										total += nextProb*nextVal;
									}
									total *= T(D_BETA);
									calcU = T(pow(T(1 - D_BETA)*T(pow(D_b, D_RHO)) + total, T(1.0 / D_RHO)));
									U_vals[j] = calcU;
								}
								T calcW;
								{
									T total = 0;

									VectorXd nextPDF = m_conditionalProbs.row(j);
									for (int ii = MAX(0, j - MAX_SHOCKS_PER_MONTH); ii < MIN(numStates, j + MAX_SHOCKS_PER_MONTH + 1); ii++) {
										total += T(nextPDF(ii))*T(1 - m_s)*latestW[ii];
									}
									calcW = prodInState - T(D_b) - arg + T(D_BETA)*total;
									W_vals[i][j] = calcW;
								}
								T partialE_partialDel;
								{
									T total = 0;

									VectorXd nextPDF = m_conditionalProbs.row(j);
									for (int ii = MAX(0, j - MAX_SHOCKS_PER_MONTH); ii < MIN(numStates, j + MAX_SHOCKS_PER_MONTH + 1); ii++) {
										T nextProb = nextPDF(ii);
										total += nextProb*T(pow(latestE[ii] - T(m_s)*(latestE[ii] - latestU[ii]), T(D_RHO)));
									}
									total *= T(D_BETA);


									T insideBracket = T(1 - D_BETA)*T(pow(T(D_b) + arg, T(D_RHO))) + total;
									partialE_partialDel = T(pow(insideBracket, T(1.0 / D_RHO - 1))*T(1 - D_BETA)*T(pow(T(D_b) + arg, T(D_RHO - 1))));
								}

								T origRetVal = calcE - calcU - m_bargaining[j] / (1 - m_bargaining[j])*calcW*partialE_partialDel;
								T retVal = ABS(origRetVal);

								if (retVal < 0) {
									std::cout << "OLGSolveAutoDiff.cpp-operator(): return value < 0. How is this possible?" << std::endl;
									exit(-1);
								}

								value = retVal;
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
								midpoint = T(0.5) * (a + b);
								tol1 = T(eps) * ABS(x) + T(tol / 3.0);
								tol2 = T(2.0) * tol1;
								//
								//  If the stopping criterion is satisfied, we can exit.
								//
								if (ABS(x - midpoint) <= (tol2 - T(0.5) * (b - a)))
								{
									status = 0;
									arg = x;
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
									q = T(2.0) * (q - r);
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
										(ABS(T(0.5) * q * r) <= ABS(p)) ||
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
											d = tol1 * (((midpoint - x)<0) ? T(-1) : T(1));
										}

										if ((b - u) < tol2)
										{
											d = tol1 * (((midpoint - x)<0) ? T(-1) : T(1));
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
									u = x + tol1 * ((d<0) ? T(-1) : T(1));
								}
								//
								//  Request value of F(U).
								//
								arg = u;
								status = status + 1;

								continue;
							}
						}//loop through 50
						 //at this point, arg is the correct wage
						wages[i][j] = T(D_b) + arg;
					}
				}
			}
		}

		//finally, return distance
		for (int i = 0; i < numStates; i++) {
			int state = i;
			T expectedW = 0;
			{
				T total = 0;
				VectorXd temp = m_conditionalProbs.row(state);
				for (int ii = 0; ii < numStates; ii++) {
					T tempVal = 0;
					for (unsigned int j = 0; j < m_gens; j++) {
						tempVal += W_vals[j][ii];
					}
					tempVal /= T(m_gens);
					total += tempVal * T(temp(ii));
				}
				expectedW = total;
			}

			T calculatedF = T(D_MY_ALPHA*D_MY_MU)*x[ii] / T(pow(T(1) + T(pow(T(D_MY_MU)*x[ii], T(m_parameter))), T(1.0 / m_parameter)));
			residual[state]=T(pow(T(D_C / D_BETA) - calculatedF / x[i] * expectedW, T(2)));
		}
		for (int i = 0; i < m_gens; i++) {
			for (int j = 0; j < numStates; j++) {
				(*m_wgs)(i, j) = double(wages[i][j]);
			}
		}
		return true;
	}
#endif
};

