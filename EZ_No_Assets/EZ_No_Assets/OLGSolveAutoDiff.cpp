#include "OLGSolveAutoDiff.h"

OLGSolveAutoDiff::OLGSolveAutoDiff(int gens, std::vector<double> ys, MatrixXd conditionalProbs, double parameter,
	/*std::vector<double> bargaining,*/ double s, MatrixXd &wages)
	:m_gens(gens),m_conditionalProbs(conditionalProbs), m_parameter(parameter), /*m_bargaining(bargaining),*/m_s(s),
	m_wgs(&wages)
{
	m_Y.resize(ys.size());
	m_Y = ys;
}

OLGSolveAutoDiff::OLGSolveAutoDiff(OLGSolveAutoDiff &orig) 
	: m_gens(orig.m_gens), m_conditionalProbs(orig.m_conditionalProbs), m_parameter(orig.m_parameter), /*m_bargaining(bargaining),*/m_s(orig.m_s),
m_wgs(orig.m_wgs) {
	m_Y.resize(orig.m_Y.size());
	m_Y = orig.m_Y;
}

OLGSolveAutoDiff::~OLGSolveAutoDiff()
{
}

double OLGSolveAutoDiff::solveProblem(std::vector<double>& x, std::vector<double>& bargaining) {
	std::vector<double> grad(x.size());
	return solveProblem(x, grad, bargaining);
}

double OLGSolveAutoDiff::solveProblem(std::vector<double>& xx, std::vector<double>& grad, std::vector<double>& bargaining) {

	//setup input
	int numStates = xx.size();
	std::vector<adouble> m_thetas(numStates);
	for (unsigned int i = 0; i < xx.size(); i++) {
		m_thetas[i] = xx[i];
	}

	stack_.new_recording();

	std::vector<adouble> m_bargaining(xx.size());
	for (unsigned int i = 0; i < xx.size(); i++) {
		m_bargaining[i] = bargaining[i];
	}

	std::vector<std::vector<adouble>> W_vals;
	std::vector<std::vector<adouble>> wages;
	W_vals.resize(m_gens);
	wages.resize(m_gens);
	for (int i = 0; i < m_gens; i++) {
		W_vals[i].resize(numStates);
		wages[i].resize(numStates);
	}

	//now, solve for wages
	{
		std::vector<adouble> E_vals(numStates);
		std::vector<adouble> U_vals(numStates);
		for (unsigned int i = m_gens - 1; i < m_gens; i--) {
			std::vector<adouble> latestE(numStates);
			latestE = E_vals;
			std::vector<adouble> latestU(numStates);
			latestU = U_vals;
			std::vector<adouble> latestW(numStates);
			std::vector<adouble> latestWages(numStates);
			if (i < m_gens - 1) {
				latestW = W_vals[i + 1];
				latestWages = wages[i + 1];
			}
			for (int j = 0; j < numStates; j++) {
#if 0
				std::cout << i << ":" << j << std::endl;
#endif
				if (i == m_gens - 1) {
					E_vals[j] = (1 - D_BETA)*pow(D_b, D_RHO);
					U_vals[j] = (1 - D_BETA)*pow(D_b, D_RHO);
					W_vals[i][j] = 0;
					wages[i][j] = D_b;
					continue;
				}
				double prodInState = m_Y[j];

				//find solution
				{
					adouble lwbnd = latestWages[j] - D_b;
					double upbnd = prodInState - D_b;

					if (lwbnd == upbnd) {
						wages[i][j] = latestWages[j];
						U_vals[j] = latestU[j];
						E_vals[j] = latestE[j];
						W_vals[i][j] = latestW[j];
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

					for (int convCounter = 0; convCounter < 100; convCounter++) {
						if (convCounter > 0) {
							adouble calcE;
							{
								adouble total = 0;

								VectorXd nextPDF = m_conditionalProbs.row(j);
								for (int ii = MAX(0, j - MAX_SHOCKS_PER_MONTH); ii < MIN(numStates, j + MAX_SHOCKS_PER_MONTH + 1); ii++) {
									total += nextPDF(ii)*pow(latestE[ii] - m_s*(latestE[ii] - latestU[ii]), D_RHO);
								}
								total *= D_BETA;
								calcE = pow((1 - D_BETA)*pow(D_b + arg, D_RHO) + total, 1.0 / D_RHO);
								E_vals[j] = calcE;
							}
							adouble calcU;
							{
								adouble total = 0;

								VectorXd nextPDF = m_conditionalProbs.row(j);
								for (int ii = MAX(0, j - MAX_SHOCKS_PER_MONTH); ii < MIN(numStates, j + MAX_SHOCKS_PER_MONTH + 1); ii++) {
									double nextProb = nextPDF(ii);
									adouble calcF = D_MY_ALPHA*D_MY_MU*m_thetas[j] / pow(1 + pow(D_MY_MU*m_thetas[j], m_parameter), 1.0 / m_parameter);
									adouble nextVal = pow(latestU[ii] + calcF*(latestE[ii] - latestU[ii]), D_RHO);
									total += nextProb*nextVal;
								}
								total *= D_BETA;
								calcU = pow((1 - D_BETA)*pow(D_b, D_RHO) + total, 1.0 / D_RHO);
								U_vals[j] = calcU;
							}
							adouble calcW;
							{
								adouble total = 0;

								VectorXd nextPDF = m_conditionalProbs.row(j);
								for (int ii = MAX(0, j - MAX_SHOCKS_PER_MONTH); ii < MIN(numStates, j + MAX_SHOCKS_PER_MONTH + 1); ii++) {
									total += nextPDF(ii)*(1 - m_s)*latestW[ii];
								}
								calcW = prodInState - D_b - arg + D_BETA*total;
								W_vals[i][j] = calcW;
							}
							adouble partialE_partialDel;
							{
								adouble total = 0;

								VectorXd nextPDF = m_conditionalProbs.row(j);
								for (int ii = MAX(0, j - MAX_SHOCKS_PER_MONTH); ii < MIN(numStates, j + MAX_SHOCKS_PER_MONTH + 1); ii++) {
									double nextProb = nextPDF(ii);
									total += nextProb*pow(latestE[ii] - m_s*(latestE[ii] - latestU[ii]), D_RHO);
								}
								total *= D_BETA;


								adouble insideBracket = (1 - D_BETA)*pow(D_b + arg, D_RHO) + total;
								partialE_partialDel = pow(insideBracket, 1.0 / D_RHO - 1)*(1 - D_BETA)*pow(D_b + arg, D_RHO - 1);
							}

							adouble origRetVal = calcE - calcU - m_bargaining[j] / (1 - m_bargaining[j])*calcW*partialE_partialDel;
							adouble retVal = ABS(origRetVal);

							if (retVal < 0) {
								std::cout << "OLGSolveAutoDiff.cpp-solve(): return value < 0. How is this possible?" << std::endl;
								std::cout << ABS(origRetVal) << std::endl;
								std::cout << adept::value(origRetVal) << ":" << (adept::value(origRetVal) < 0) << ":" << adept::value(retVal) << std::endl;
								exit(-1);
							}

							value = retVal;
#if 0
							std::cout.precision(15);
							std::cout << "OLGSolveAutoDiff: f(" << i << "," << j << "," << adept::value(arg) << ")=" << adept::value(retVal) << std::endl;
							if (i == 0 && j == 0) {
								std::cout << "E: ";
								for (int ii = 0; ii < numStates; ii++) {
									std::cout << latestE[ii] << ",";
								}
								std::cout << std::endl;
								std::cout << "U: ";
								for (int ii = 0; ii < numStates; ii++) {
									std::cout << latestU[ii] << ",";
								}
								std::cout << std::endl;
								std::cout << "W: ";
								for (int ii = 0; ii < numStates; ii++) {
									std::cout << W_vals[i+1][ii] << ",";
								}
								std::cout << std::endl;
								std::cout << calcE << ":" << calcU << ":" << calcW << ":" << partialE_partialDel << ":" << m_bargaining << std::endl;
								exit(-1);
							}
#endif
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
					}//loop through 50
					//at this point, arg is the correct wage
					wages[i][j] = D_b + arg;
					if (i < m_gens - 1) {
						if (wages[i][j] < wages[i + 1][j]) {
							std::cout << "ERROR! OLGSolveAutoDiff.solveProblem() - wages not non-decreasing. " << std::endl;
							std::cout << "Wage(" << i << ", " << j << ") = " << wages[i][j] << std::endl;
							std::cout << "Wage(" << i+1 << ", " << j << ") = " << wages[i+1][j] << std::endl;
							exit(-1);
						}
					}
				}
			}
		}
	}

	//finally, return distance
	adouble retVal = 0;
	for (int i = 0; i < xx.size(); i++) {
		int state = i;
		adouble expectedW = 0;
		{
			adouble total = 0;
			VectorXd temp = m_conditionalProbs.row(state);
			for (int ii = 0; ii < numStates; ii++) {
				adouble tempVal = 0;
				for (unsigned int j = 0; j < m_gens; j++) {
					tempVal += W_vals[j][ii];
				}
				tempVal /= m_gens;
				total += tempVal * temp(ii);
			}
			expectedW = total;
		}

		adouble calculatedF = D_MY_ALPHA*D_MY_MU*m_thetas[i] / pow(1 + pow(D_MY_MU*m_thetas[i], m_parameter), 1.0 / m_parameter);
		retVal += pow(D_C / D_BETA - calculatedF / m_thetas[i] * expectedW, 2);
		if (value(retVal) != value(retVal)) {
			std::cout << "OLGSolveAUtoDiff.solveProblem(): retval is NaN" << std::endl;
			std::cout << "Theta(" << i << ")=" << m_thetas[i] << std::endl;
			std::cout << "ExpectedW=" << expectedW << std::endl;
			exit(-1);
		}
	}
	retVal.set_gradient(1.0);
	stack_.compute_adjoint();
	for (int i = 0; i < xx.size(); i++) {
		grad[i] = m_thetas[i].get_gradient();
	}

	for (int i = m_gens-1; i >= 0; i--) {
		for (int j = 0; j < numStates; j++) {
			(*m_wgs)(i,j) = adept::value(wages[i][j]);
			if (i < m_gens - 1) {
				if ((*m_wgs)(i, j) < (*m_wgs)(i + 1, j)) {
					std::cout << "ERROR! OLGSolveAutoDiff.solveProblem() - wages not non-decreasing. " << std::endl;
					std::cout << "Wage(" << i << ", " << j << ") = " << wages[i][j] << std::endl;
					std::cout << "Wage(" << i + 1 << ", " << j << ") = " << wages[i + 1][j] << std::endl;
					exit(-1);
				}
			}
		}
	}
	return value(retVal);
}

void OLGSolveAutoDiff::setBargaining(std::vector<double>& b) {
	m_bargaining = b;
}
