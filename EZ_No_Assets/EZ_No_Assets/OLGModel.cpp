#include "OLGModel.h"
#include "nlopt.hpp"
#include <exception>

OLGModel::OLGModel(unsigned int generations, double y, double s, MatchingFunction &f, ShockProcess &sp, double bargaining, bool autoDiff)
	: m_gens(generations), m_f(&f), /*m_bargaining(bargaining),*/ m_Es(s), m_sp(&sp), m_Y(sp.numStates()),
	E_vals(generations, sp.numStates()), U_vals(generations, sp.numStates()),
	W_vals(generations, sp.numStates()), wages(generations, sp.numStates()),
	m_thetas(sp.numStates()), m_autodiff(autoDiff)
{
	for (unsigned int i = 0; i < sp.numStates(); i++) {
		m_thetas[i] = f.getTheta();
	}

	const VectorXd* x = sp.states();
	for (unsigned int i = 0; i < sp.numStates(); i++) {
		m_Y(i) = D_b + exp((*x)(i))*(y - D_b);
	}

	//error check
	bool foundError = false;
	for (unsigned int i = 0; i < sp.numStates(); i++) {
		if (m_Y(i) < D_b) {
			foundError = true;
			std::cout << "Error! OLGModel.constructor(): Y(" << i << ")=" << m_Y(i) << " which is less than outside option, b=" << D_b << std::endl;
		}
	}
	if (foundError) {
		exit(-1);
	}

}

OLGModel::~OLGModel()
{
}

void OLGModel::solveWages()
{
	for (unsigned int i = m_gens - 1; i < m_gens; i--) {
		lastSolveGen = i + 1;

		#pragma omp parallel for num_threads(3)
		for (int j = 0; j < m_sp->numStates(); j++) {
#if 0
			std::cout << i << ":" << j << std::endl;
#endif
			if (i == m_gens - 1) {
				E_vals(i, j) = (1 - D_BETA)*pow(D_b, D_RHO);
				U_vals(i, j) = (1 - D_BETA)*pow(D_b, D_RHO);
				W_vals(i, j) = 0;
				wages(i, j) = D_b;
			}
			else {
				VectorXd latestU = U_vals.row(i + 1);
				VectorXd latestE = E_vals.row(i + 1);
				VectorXd latestW = W_vals.row(i + 1);
				double latestWages = wages(i + 1, j);
				double prodInState = m_Y(j);

				Generation toSolve(*this, &OLGModel::nonLinearWageEquation, j, latestU, latestE, latestW);

#if 0
				if (i == 0 && j == 0) {
					std::cout << "f value = " << toSolve(0.440577330755482) << std::endl;
					int numStates = m_sp->numStates();
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
						std::cout << latestW[ii] << ",";
					}
					std::cout << std::endl;
					exit(-1);
				}
		
#endif
#if 0
				MySolver solver(toSolve, pow(10,-30));

				double del = solver.solve(latestWages, prodInState);
#else
				nlopt::opt opt(nlopt::LN_BOBYQA, 1);
				double lwbnd = latestWages - D_b;
				double upbnd = prodInState - D_b;

				if (lwbnd == upbnd) {
#if 0
					std::cout << "Error! OLGModel.cpp-solveWages() : can only set wage equal to previous one. Really don't want this as it is wrong." << std::endl;
					std::cout << "w=" << lwbnd + D_b << std::endl;
					exit(-1);
#endif
					wages(i, j) = latestWages;
					U_vals(i, j) = U_vals(i + 1, j);
					E_vals(i, j) = E_vals(i + 1, j);
					W_vals(i, j) = W_vals(i + 1, j);

				}
				else {

					opt.set_lower_bounds(lwbnd);
					opt.set_upper_bounds(upbnd);
					opt.set_min_objective(Generation::wrap, &toSolve);

					opt.set_xtol_rel(1e-10);

					std::vector<double> x(1);
					x[0] = (lwbnd + upbnd) / 2;

					double minf;
					nlopt::result result;
					try {
						result = opt.optimize(x, minf);
					}
					catch (const nlopt::roundoff_limited& e) {
						std::cout << e.what() << std::endl;
					}
					catch (const nlopt::forced_stop& e) {
						std::cout << e.what() << std::endl;
						exit(-1);
					}
					catch (const std::runtime_error& e) {
						std::cout << e.what() << std::endl;
						exit(-1);
					}
					catch (const std::exception &e) {
						std::cout << e.what() << std::endl;
						exit(-1);
					}
					double del = x[0];
#endif
					double newWages = D_b + del;

					if (newWages < latestWages) {
						std::cout << "OLGModel-solveWages(): cohort " << i << " shock " << j << " has wageT < wage(T+1). HOW?" << std::endl;
						std::cout << "w" << i << "=" << newWages << ", w" << i + 1 << "=" << latestWages << std::endl;
						std::cout << "diff: " << newWages - latestWages << std::endl;
						std::cout << "y: " << m_Y(j) << "   y-w(T+1): " << m_Y(j) - latestWages << std::endl;
						for (int lll = m_gens - 1; lll > i; lll--) {
							std::cout << "y-w(" << lll << "): " << m_Y(j) - wages(lll, j) << std::endl;
						}
						std::cout << "value: " << nonLinearWageEquation(j, del, latestU, latestE, latestW) << std::endl;
						std::cout << "value0: " << nonLinearWageEquation(j, m_Y(j) - D_b, latestU, latestE, latestW) << std::endl;
						exit(-1);
					}
					if (m_Y(j) - newWages < 0) {
						std::cout << "OLGModel-solveWages(): cohort " << i << " shock " << j << " has y-w<0. HOW?" << std::endl;
						std::cout << "y: " << m_Y(j) << "   w: " << newWages << "   y-w(T+1): " << m_Y(j) - newWages << std::endl;
						for (int lll = m_gens - 1; lll > i; lll--) {
							std::cout << "y-w(" << lll << "): " << m_Y(j) - wages(lll, j) << std::endl;
						}
						std::cout << "value: " << nonLinearWageEquation(j, del, latestU, latestE, latestW) << std::endl;
						std::cout << "value0: " << nonLinearWageEquation(j, m_Y(j) - D_b, latestU, latestE, latestW) << std::endl;
						exit(-1);
					}

					wages(i, j) = newWages;
					U_vals(i, j) = calcU(j, latestU, latestE);
					E_vals(i, j) = calcE(j, del, latestU, latestE);
					W_vals(i, j) = calcW(j, del, latestW);
				}
			}
		}
	}
	return;
}

double OLGModel::calcU(int state, VectorXd &Up1, VectorXd &Ep1) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	for (int i = MAX(0,state-MAX_SHOCKS_PER_MONTH); i < MIN(m_sp->numStates(),state+MAX_SHOCKS_PER_MONTH+1); i++) {
//	for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		double nextProb = nextPDF(i, 0);
		double calcF = m_f->calculatedF(m_thetas[i]);
		double nextVal = pow(Up1(i) + calcF*(Ep1(i) - Up1(i)), D_RHO);
		total += nextProb*nextVal;
	}
	total *= D_BETA;
	double retVal = pow((1 - D_BETA)*pow(D_b, D_RHO) + total, 1.0 / D_RHO);
	return retVal;
}

double OLGModel::calcE(int state, double delta, VectorXd &Up1, VectorXd &Ep1) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1); i++) {
//		for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		total += nextPDF(i, 0)*pow(Ep1(i) - m_Es*(Ep1(i) - Up1(i)), D_RHO);
	}
	total *= D_BETA;
	double retVal = pow((1 - D_BETA)*pow(D_b + delta, D_RHO) + total, 1.0 / D_RHO);
	return retVal;
}

double OLGModel::calcW(int state, double delta, VectorXd &Wp1) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1); i++) {
		//for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		total += nextPDF(i, 0)*(1 - m_Es)*Wp1(i);
	}
	double retVal = m_Y(state) - D_b - delta + D_BETA*total;
	return retVal;
}

double OLGModel::nonLinearWageEquation(int state, double x, VectorXd& Up1, VectorXd& Ep1, VectorXd& Wp1) {
	double t_calcE = calcE(state, x, Up1, Ep1);
	double t_calcU = calcU(state, Up1, Ep1);
	double t_calcW = calcW(state, x, Wp1);
	double t_partialE_partialDel = partialE_partialDel(state, x, Up1, Ep1);


#if 0
	std::cout.precision(15);
	std::cout << t_calcE << ":" << t_calcU << ":" << t_calcW << ":" << t_partialE_partialDel << ":" << m_bargaining << std::endl;
#endif
	double bargaining = 1 - m_f->getElasticity(m_thetas(state));
	double retVal = ABS(
//			t_calcE - t_calcU - m_bargaining / (1 - m_bargaining)*t_calcW*t_partialE_partialDel
			t_calcE - t_calcU - bargaining / (1 - bargaining)*t_calcW*t_partialE_partialDel
		);

	if (retVal < 0) {
		std::cout << "OLGModel.cpp-nonLinearWageEquation(): return value < 0. How is this possible?" << std::endl;
		exit(-1);
	}
	return retVal;
}

double OLGModel::partialE_partialDel(int state, double x, VectorXd& Up1, VectorXd& Ep1) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1); i++) {
		//for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		total += nextPDF(i, 0)*pow(Ep1(i) - m_Es*(Ep1(i) - Up1(i)), D_RHO);
	}
	total *= D_BETA;


	double insideBracket = (1 - D_BETA)*pow(D_b + x, D_RHO) + total;
	double dE_dDel = pow(insideBracket, 1.0 / D_RHO - 1)*(1-D_BETA)*pow(D_b+x,D_RHO-1);
	return dE_dDel;
}

double OLGModel::expectedW(int state, bool forceNoShocks) {
	double total = 0;
	if (forceNoShocks || (m_sp->numStates() == 1)) {
		for (unsigned int i = 0; i < m_gens; i++) {
			total += W_vals(i, state);
		}
		total /= m_gens;
		return total;
	}

	pdfMatrix temp = m_sp->nextPeriodPDF(state);
	for (int i = 0; i < m_sp->numStates(); i++) {
		double tempVal = 0;
		for (unsigned int j = 0; j < m_gens; j++) {
			tempVal += W_vals(j, i);
		}
		tempVal /= m_gens;
		total += tempVal * temp(i, 0);
	}
	return total;
}

double OLGModel::elasticityWRTymb() {
	unsigned int expectedState = (m_sp->numStates() - 1) / 2;
	double Ey = m_Y(expectedState);
	double bargaining = 1 - m_f->getElasticity(m_thetas(expectedState));
	OLGModel thetaChange(m_gens, Ey, m_Es, *(m_f->dTheta()), *m_sp, bargaining, false);
	OLGModel yChange(m_gens, 1.0001*Ey, m_Es, *m_f, *m_sp, bargaining, false);

	thetaChange.solveWages();
	yChange.solveWages();

	double EW = expectedW(expectedState, true);
	double num = (Ey - D_b)*(yChange.expectedW(expectedState, true) - EW) / (1.0001*Ey - Ey);
	double denom = bargaining*EW - m_f->getTheta()*
		(thetaChange.expectedW(expectedState, true) - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

	delete thetaChange.m_f;
	return num/denom;
}

double OLGModel::elasticityWRTs() {
	OLGModel thetaChange(*this);
	thetaChange.m_f = m_f->dTheta();

	unsigned int expectedState = (m_sp->numStates() - 1) / 2;
	double bargaining = 1 - m_f->getElasticity(m_thetas(expectedState));

	OLGModel sChange(m_gens, m_Y(expectedState), 1.0001*m_Es, *m_f, *m_sp, bargaining, false);

	thetaChange.solveWages();
	sChange.solveWages();

	double EW = expectedW(expectedState, true);
	double num = (m_Es)*(sChange.expectedW(expectedState, true) - EW) / (sChange.m_Es - m_Es);
	double denom = bargaining*EW - m_f->getTheta()*
		(thetaChange.expectedW(expectedState, true) - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

	delete thetaChange.m_f;
	return num / denom;
}

void OLGModel::printWages() {
	std::cout.precision(15);
	for (unsigned int i = 0; i < m_gens; i++) {
		for (unsigned int j = 0; j < m_sp->numStates(); j++) {
			std::cout << "Cohort (" << i << "," << j << "): y=" << m_Y(j) << " b=" << D_b << " w=" << wages(i, j)<< std::endl;
		}
	}
	return;
}

double OLGModel::operator()(const std::vector<double> &x, std::vector<double> &grad)
{
	static int counter = 0;

	if (x.size() != m_sp->numStates()) {
		std::cout << "Error! OLGModel.cpp::operator() - numStates and xsize (size of theta being solved) are not equal."
			<< std::endl;
		std::cout << "x:" << x.size() << ", states:" << m_sp->numStates() << std::endl;
		exit(-1);
	}
	double retVal = 0;
	for (unsigned int i = 0; i < x.size(); i++) {
		m_thetas(i) = x[i];
	}
	if (!m_autodiff) {
		solveWages();
		for (int i = 0; i < x.size(); i++) {
			retVal += pow(D_C / D_BETA - m_f->calculatedF(m_thetas[i]) / m_thetas[i] * expectedW(i), 2);
		}
	}
	else {
		std::vector<double> myYs(m_Y.size());
		VectorXd::Map(&myYs[0], m_Y.size()) = m_Y;
		std::vector<double> bargaining(m_Y.size());
		for (int i = 0; i < x.size(); i++) {
			bargaining[i] = 1 - m_f->getElasticity(m_thetas[i]);
		}
		OLGSolveAutoDiff soln(m_gens, myYs, m_sp->getProbMatrix(), m_f->getParameter()/*, bargaining*/, m_Es, wages);

		std::vector<double> myThetas(m_thetas.size());
		VectorXd::Map(&myThetas[0], m_thetas.size()) = m_thetas;
		retVal = soln.solveProblem(myThetas,grad,bargaining);
	}
	return sqrt(retVal);

}

double OLGModel::wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {

	static double bestSoFar = 10000;
	static unsigned int m_counter=0;
	static int lastStateSize = 3;

	if (lastStateSize != x.size()) {
		lastStateSize = x.size();
		m_counter = 0;
		bestSoFar = 10000;
	}

	OLGModel modelToSolve = *(OLGModel *)data;
	std::cout << "OLGModel.cpp-wrap(): Starting evaluation " << ++m_counter << std::endl << std::flush;
	double value = (*reinterpret_cast<OLGModel*>(data))(x, grad);

	std::cout << "OLGModel.cpp-wrap(): Evaluation " << m_counter << ": " << value << std::endl << std::flush;

	if (value < bestSoFar) {
		bestSoFar = value;
		printStatus(x, m_counter, value);
	}

	return value;
}

void OLGModel::printStatus(const std::vector<double>& solution, int numCalls, double distance) {
	using namespace std;

	std::cout << "OLGModel.printStatus(): At evals=" << numCalls << " new best has distance=" << distance
		<< " where current solution is: " << std::endl;
	for (int i = 0; i < solution.size(); i++) {
		printf("%4i %10lf\r\n", i, solution[i]);
	}
	std::cout << std::flush;

	ostringstream os;
	ofstream out_stream;
	os << "intermediateResults.dat";
	out_stream.precision(15);
	out_stream << std::scientific;
	out_stream.open(os.str(), std::ofstream::out/* | std::ofstream::app*/);

	for (int i = 0; i < solution.size(); i++) {
		out_stream << solution[i] << std::endl;
	}
	out_stream.close();

	return;
}
