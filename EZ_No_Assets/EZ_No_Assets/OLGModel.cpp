#include "OLGModel.h"
#include "MySolver.h"
#include "nlopt.hpp"
#include <exception>

OLGModel::OLGModel(int generations, double y, double s, MatchingFunction &f, ShockProcess &sp, double bargaining, bool autoDiff)
	: m_gens(generations), m_f(&f), m_bargaining(bargaining), m_Es(s), m_sp(&sp), m_Y(generations),
	E_vals(generations), U_vals(generations, sp.numStates()),
	W_vals(generations), wages(generations),
	m_thetas(sp.numStates()), m_autodiff(autoDiff)
{
	if (generations < 1) {
		std::cout << "OLGModel.constructor() - cannot pass fewer than 1 generation. gens=" << generations << std::endl;
		exit(-1);
	}

	for (unsigned int i = 0; i < sp.numStates(); i++) {
		m_thetas[i] = f.getTheta();
	}

	const VectorXd* x = sp.states();
	for (int genIndex = 0; genIndex < generations; genIndex++) {
		m_Y[genIndex] = MatrixXd(generations, sp.numStates());
		E_vals[genIndex] = MatrixXd(generations, sp.numStates());
		W_vals[genIndex] = MatrixXd(generations, sp.numStates());
		wages[genIndex] = MatrixXd(generations, sp.numStates());
		for (int genIndex2 = 0; genIndex2 < generations; genIndex2++) {
			for (unsigned int i = 0; i < sp.numStates(); i++) {
				m_Y[genIndex](genIndex2, i) = D_b + exp((*x)(i))*(pow(D_TENURE_INCREASE, genIndex)*pow(D_PROD_INCREASE, genIndex2)*y - D_b);
			}
		}
	}

	//error check
	bool foundError = false;
	for (int genIndex = 0; genIndex < generations; genIndex++) {
		for (int genIndex2 = 0; genIndex2 < generations; genIndex2++) {
			for (unsigned int i = 0; i < sp.numStates(); i++) {
				if (m_Y[genIndex2](genIndex, i) < D_b) {
					foundError = true;
					std::cout << "Error! OLGModel.constructor(): Y_" << genIndex2 << "(" << genIndex << "," << i << ")="
						<< m_Y[genIndex2](genIndex, i) << " which is less than outside option, b=" << D_b << std::endl;
				}
			}
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
	for (int i = m_gens - 1; i >= 0; i--) {
//#pragma omp parallel for num_threads(3)
		for (int tenureIndex = 0; tenureIndex <= i; tenureIndex++) {
			for (int j = 0; j < m_sp->numStates(); j++) {
#if 0
				std::cout << i << ":" << j << std::endl;
#endif
				if (i == m_gens - 1) {
					E_vals[tenureIndex](i, j) = pow(1 - D_BETA, 1.0 / D_RHO)*D_b;
					U_vals(i, j) = pow(1 - D_BETA, 1.0 / D_RHO)*D_b;
					W_vals[tenureIndex](i, j) = 0;
					wages[tenureIndex](i, j) = D_b;
				}
				else {
					/*
					VectorXd nextU = U_vals.row(i + 1);
					VectorXd nextE = E_vals[tenureIndex + 1].row(i + 1);
					VectorXd nextW = W_vals[tenureIndex + 1].row(i + 1);
					*/
					double prodInState = m_Y[tenureIndex](i, j);

					Generation toSolve(*this, &OLGModel::nonLinearWageEquation, j, i, tenureIndex/*, nextU, nextE, nextW */);
					nlopt::opt opt(nlopt::LN_BOBYQA, 1);
					//				double lwbnd = latestWages - D_b;
					double lwbnd = 0;
					double upbnd = prodInState - D_b;

#if 0
					//don't think we can do this anymore
					if (lwbnd == upbnd) {
						wages(i, j) = latestWages;
						U_vals(i, j) = U_vals(i + 1, j);
						E_vals(i, j) = E_vals(i + 1, j);
						W_vals(i, j) = W_vals(i + 1, j);
					}
					else {
#endif

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
						double newWages = D_b + del;

#if 0
						if (newWages < latestWages) {
							std::cout << "OLGModel-solveWages(): cohort " << i << " shock " << j << " has wageT < wage(T+1). HOW?" << std::endl;
							std::cout << "w" << i << "=" << newWages << ", w" << i + 1 << "=" << latestWages << std::endl;
							std::cout << "diff: " << newWages - latestWages << std::endl;
							std::cout << "y: " << m_Y(i, j) << "   y-w(T+1): " << m_Y(i, j) - latestWages << std::endl;
							for (int lll = m_gens - 1; lll > i; lll--) {
								std::cout << "y-w(" << lll << "): " << m_Y(i, j) - wages(lll, j) << std::endl;
							}
							std::cout << "value: " << nonLinearWageEquation(j, del, i, latestU, latestE, latestW) << std::endl;
							std::cout << "value0: " << nonLinearWageEquation(j, m_Y(i, j) - D_b, i, latestU, latestE, latestW) << std::endl;
							exit(-1);
						}
#endif
						if (m_Y[tenureIndex](i, j) - newWages < 0) {
							std::cout << "OLGModel-solveWages(): cohort " << i << " shock " << j << " tenure " << tenureIndex << " has y-w<0. HOW?" << std::endl;
							std::cout << "y: " << m_Y[tenureIndex](i, j) << "   w: " << newWages << "   y-w(T+1): " << m_Y[tenureIndex](i, j) - newWages << std::endl;
							for (int lll = m_gens - 1; lll > i; lll--) {
								std::cout << "y-w(" << lll << "): " << m_Y[tenureIndex](lll, j) - wages[tenureIndex](lll, j) << std::endl;
							}
							double val1 = nonLinearWageEquation(j, i, tenureIndex, del/*, nextU, nextE, nextW*/);
							double val2 = nonLinearWageEquation(j, i, tenureIndex, m_Y[tenureIndex](i, j) - D_b/*, nextU, nextE, nextW*/);
							std::cout << "value: " << val1 << std::endl;
							std::cout << "value0: " << val2 << std::endl;
							std::cout << "diff: " << val2 - val1 << std::endl;
							exit(-1);
						}

						wages[tenureIndex](i, j) = newWages;
						if (tenureIndex == 0) {
							//VectorXd tempNextE = E_vals[0].row(i + 1);
							U_vals(i, j) = calcU(j, i/*, nextU, tempNextE*/);
						}
						E_vals[tenureIndex](i, j) = calcE(j, i, tenureIndex, del/*, nextU, nextE*/);
						W_vals[tenureIndex](i, j) = calcW(j, i, tenureIndex, del/*, nextW*/ );
				}
			}
		}
	}
#if 0
	}
#endif
	return;
}

double OLGModel::calcU(int state, int whichGen/*, VectorXd &Up1, VectorXd &Ep1*/) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	int maxIters = MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1);
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
		double nextProb = nextPDF(i, 0);
		double calcF = m_f->calculatedF(m_thetas[state]);
		double nextVal = D_DEATH*(U_vals(m_gens-1,i))+(1-D_DEATH)*(pow(U_vals(whichGen+1, i) + calcF*(E_vals[0](whichGen+1,i) - U_vals(whichGen + 1, i)), D_RHO));
		total += nextProb*nextVal;
	}
	total *= D_BETA;
	double retVal = pow((1 - D_BETA)*pow(D_b, D_RHO) + total, 1.0 / D_RHO);
	return retVal;
}

double OLGModel::calcE(int state, int whichGen, int tenure, double delta/*, VectorXd &Up1, VectorXd &Ep1*/) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	int maxIters = MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1);
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
		double nextVal = D_DEATH*pow(1 - D_BETA, 1.0 / D_RHO)*D_b 
			+ (1 - D_DEATH)*(pow(E_vals[tenure+1](whichGen+1,i) - m_Es*(E_vals[tenure + 1](whichGen + 1, i) - U_vals(whichGen+1,i)), D_RHO));
		total += nextPDF(i, 0)*nextVal;
	}
	total *= D_BETA;
	double retVal = pow((1 - D_BETA)*pow(D_b + delta, D_RHO) + total, 1.0 / D_RHO);
	return retVal;
}

double OLGModel::calcW(int state, int whichGen, int whichTenure, double delta/*, VectorXd &Wp1*/) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	int maxIters = MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1);
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
		total += nextPDF(i, 0)*(1 - m_Es - D_DEATH)*W_vals[whichTenure + 1](whichGen + 1, i);
	}
	double retVal = m_Y[whichTenure](whichGen, state) - D_b - delta + D_BETA*total;
	return retVal;
}

double OLGModel::nonLinearWageEquation(int state, int whichGen, int tenure, double x/*, VectorXd& Up1, VectorXd& Ep1, VectorXd& Wp1*/) {
	double t_calcE = calcE(state, whichGen, tenure, x/*, Up1, Ep1*/);
	double t_calcU = calcU(state, whichGen/*, Up1, Ep1*/);
	double t_calcW = calcW(state, whichGen, tenure, x/*, Wp1*/);
	double t_partialE_partialDel = partialE_partialDel(state, whichGen, tenure, x/*, Up1, Ep1*/);

	if (m_bargaining == 1) {
		return -(t_calcE - t_calcU);
	}
	//double bargaining = 1 - m_f->getElasticity(m_thetas(state));
	double retVal = ABS(
		t_calcE - t_calcU - m_bargaining / (1 - m_bargaining)*t_calcW*t_partialE_partialDel
		//			t_calcE - t_calcU - bargaining / (1 - bargaining)*t_calcW*t_partialE_partialDel
		);

#if 0
	std::cout << "trial=" << x << std::endl;
	std::cout << "E=" << t_calcE << std::endl;
	std::cout << "U=" << t_calcU << std::endl;
	std::cout << "W=" << t_calcW << std::endl;
	std::cout << "p=" << t_partialE_partialDel << std::endl;
	std::cout << "r=" << retVal << std::endl;
	std::cout << "============================" << std::endl;
#endif
	if (retVal < 0) {
		std::cout << "OLGModel.cpp-nonLinearWageEquation(): return value < 0. How is this possible?" << std::endl;
		exit(-1);
	}
	return retVal;
}

double OLGModel::partialE_partialDel(int state, int whichGen, int tenure, double x/*, VectorXd& Up1, VectorXd& Ep1*/) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1); i++) {
		double nextVal = D_DEATH*pow(1 - D_BETA, 1.0 / D_RHO)*D_b 
			+ (1 - D_DEATH)*(pow(E_vals[tenure+1](whichGen+1,i) - m_Es*(E_vals[tenure + 1](whichGen + 1, i) - U_vals(whichGen+1,i)), D_RHO));
		total += nextPDF(i, 0)*nextVal;
	}
	total *= D_BETA;


	double insideBracket = (1 - D_BETA)*pow(D_b + x, D_RHO) + total;
	double dE_dDel = pow(insideBracket, 1.0 / D_RHO - 1)*(1 - D_BETA)*pow(D_b + x, D_RHO - 1);
	return dE_dDel;
}

double OLGModel::expectedW0(int state, bool forceNoShocks) {
	double total = 0;
	if (forceNoShocks || (m_sp->numStates() == 1)) {
		for (int i = 0; i < m_gens; i++) {
			total += W_vals[0](i, state);
		}
		total /= m_gens;
		return total;
	}

	pdfMatrix temp = m_sp->nextPeriodPDF(state);
	for (int i = 0; i < m_sp->numStates(); i++) {
		double tempVal = 0;
		for (int j = 0; j < m_gens; j++) {
			tempVal += W_vals[0](j, i);
		}
		tempVal /= m_gens;
		total += tempVal * temp(i, 0);
	}
	return total;
}

double OLGModel::elasticityWRTymb() {
	unsigned int expectedState = (m_sp->numStates() - 1) / 2;
	double Ey = m_Y[0](0, expectedState);
	OLGModel thetaChange(m_gens, Ey, m_Es, *(m_f->dTheta()), *m_sp, m_bargaining, false);
	OLGModel yChange(m_gens, 1.0001*Ey, m_Es, *m_f, *m_sp, m_bargaining, false);

	Ey = 0;
	double Eb = 0;
	for (int i = 0; i < m_gens-1; i++) {
		Ey += m_Y[0](i, expectedState);
		Eb += D_b;
	}
	Ey /= m_gens-1;
	Eb /= m_gens-1;
	thetaChange.solveWages();
	yChange.solveWages();

	double EW = expectedW0(expectedState, true);
	double tEW = thetaChange.expectedW0(expectedState, true);
	double yEW = yChange.expectedW0(expectedState, true);
	double num = (Ey - Eb)*(yEW - EW) / (1.0001*Ey - Ey);
	double denom = (1 - m_f->getElasticity(m_thetas(expectedState)))*EW - m_f->getTheta()*
		(tEW - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

	if (denom == 0) {
		std::cout << "OLGModel-elasticityWRTymb() - denominator = 0. Should not be possible." << std::endl;
		std::cout << "EW=" << EW << std::endl;
		std::cout << "tEW=" << tEW << std::endl;
		std::cout << "yEW=" << yEW << std::endl;
		exit(-1);
	}
	delete thetaChange.m_f;
	return num / denom;
}

std::vector<double> OLGModel::wageElasticityWRTymb() {

	std::cout << "ERROR! OLGModel.wageElasticityWRTymb() is not working!" << std::endl;
	exit(-1);

#if 0
	unsigned int expectedState = (m_sp->numStates() - 1) / 2;
	double Ey = m_Y(0, expectedState);
	double bargaining = 1 - m_f->getElasticity(m_thetas(expectedState));

	//FIX THIS. AUTODIFF IF WE HAVE TO?
	OLGModel yChange(m_gens, 1.0001*Ey, m_Es, *m_f, *m_sp, bargaining, false);

	solveWages();
	yChange.solveWages();

	printWages();
	yChange.printWages();

	std::vector<double> expectedWage(m_sp->numStates());
	std::vector<double> newExpectedWage(m_sp->numStates());
	for (int j = 0; j < m_sp->numStates(); j++) {
		for (int i = 0; i < m_gens; i++) {
			if (i == 0) {
				expectedWage[j] = 0;
				newExpectedWage[j] = 0;
			}
			expectedWage[j] += (wages(i, j) / m_gens);
			newExpectedWage[j] += (yChange.wages(i, j) / m_gens);
		}
	}
	std::vector<double> retVal(m_sp->numStates());
	for (int i = 0; i < m_sp->numStates(); i++) {
		if (newExpectedWage[i] < expectedWage[i]) {
			std::cout << "ERROR! New wage should always be higher than old wage. wage " << i << std::endl;
			std::cout << "Average New: " << newExpectedWage[i] << ", Old: " << expectedWage[i] << std::endl;
			for (int j = 0; j < m_gens; j++) {
				std::cout << "New: " << yChange.wages(j, i) << ", Old: " << wages(j, i) << std::endl;
			}
			exit(-1);
		}
		retVal[i] = (Ey - D_b) / (expectedWage[i] - D_b)*(newExpectedWage[i] - expectedWage[i]) / (1.0001*Ey - Ey);
	}

	return retVal;
#else
	std::vector<double> retVal(m_sp->numStates());
	return retVal;
#endif
}

double OLGModel::elasticityWRTs() {
	std::cout << "OLGModel.cpp-elasticityWRTs(): Probably not correct either. " << std::endl;
	exit(-1);
#if 0
	OLGModel thetaChange(*this);
	thetaChange.m_f = m_f->dTheta();

	unsigned int expectedState = (m_sp->numStates() - 1) / 2;
	//double bargaining = 1 - m_f->getElasticity(m_thetas(expectedState));

	OLGModel sChange(m_gens, m_Y(0, expectedState), 1.0001*m_Es, *m_f, *m_sp, m_bargaining, false);

	thetaChange.solveWages();
	sChange.solveWages();

	double EW = expectedW(expectedState, true);
	double num = (m_Es)*(sChange.expectedW(expectedState, true) - EW) / (sChange.m_Es - m_Es);
	double denom = (1 - m_f->getElasticity(m_thetas(expectedState)))*EW - m_f->getTheta()*
		(thetaChange.expectedW(expectedState, true) - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

	delete thetaChange.m_f;
	return num / denom;
#else
	return 0;
#endif
}

void OLGModel::printWages() {
	std::cout.precision(15);
	for (int i = 0; i < m_gens; i++) {
		for (int k = 0; k <= ((D_TENURE_INCREASE==1)?0:i); k++) {
			for (unsigned int j = 0; j < m_sp->numStates(); j++) {
				std::cout << "Cohort_" << i << " (" << j << "," << k << "): y=" << m_Y[k](i, j) << " b=" << D_b
					<< " w=" << wages[k](i, j) << std::endl;
			}
		}
	}
	return;
}

double OLGModel::operator()(const std::vector<double> &x, std::vector<double> &grad)
{
	std::cout << "OLGModel::operator() not implemented yet." << std::endl;
	exit(-1);
	return 0;
#if 0
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
			retVal += pow(D_C / D_BETA - m_f->calculatedF(m_thetas[i]) / m_thetas[i] * expectedW0(i), 2);
		}
	}
	else {
		std::cout << "ERROR! OLGModel.cpp-operator() - need to fix autodiff to account for different productivity." << std::endl;
		exit(-1);
		std::vector<double> myYs(m_Y.size());
		VectorXd::Map(&myYs[0], m_Y.size()) = m_Y;
		std::vector<double> bargaining(m_Y.size());
		for (int i = 0; i < x.size(); i++) {
			bargaining[i] = m_bargaining;// 1 - m_f->getElasticity(m_thetas[i]);
		}
		OLGSolveAutoDiff soln(m_gens, myYs, m_sp->getProbMatrix(), m_f->getParameter()/*, bargaining*/, m_Es, wages);

		std::vector<double> myThetas(m_thetas.size());
		VectorXd::Map(&myThetas[0], m_thetas.size()) = m_thetas;
		retVal = soln.solveProblem(myThetas, grad, bargaining);
	}
	return sqrt(retVal);
#endif
}

OLGSolveAutoDiff OLGModel::getSolver() {
	std::cout << "OLGModel::getSolver() not implemented yet." << std::endl;
	exit(-1);
	std::vector<double> temp(2);
	MatrixXd temp2;
	MatrixXd temp3;
	OLGSolveAutoDiff soln(1, temp, temp2, 0.0, 0.0, temp3);

	return soln;

#if 0
	std::vector<double> myYs(m_Y.size());
	VectorXd::Map(&myYs[0], m_Y.size()) = m_Y;
	std::vector<double> bargaining(m_Y.size());
	OLGSolveAutoDiff soln(m_gens, myYs, m_sp->getProbMatrix(), m_f->getParameter(), m_Es, wages);

	for (int i = 0; i < myYs.size(); i++) {
		bargaining[i] = m_bargaining;// 1 - m_f->getElasticity(1);
	}

	soln.setBargaining(bargaining);
	return soln;
#endif
}

double OLGModel::wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {

	static double bestSoFar = 10000;
	static unsigned int m_counter = 0;
	static int lastStateSize = 3;

	if (lastStateSize != x.size()) {
		lastStateSize = x.size();
		m_counter = 0;
		bestSoFar = 10000;
	}

	OLGModel modelToSolve = *(OLGModel *)data;
	std::cout << "OLGModel.cpp-wrap(): Starting evaluation " << ++m_counter << std::endl << std::flush;
	double value = (*reinterpret_cast<OLGModel*>(data))(x, grad);

	if (value != value) {
		std::cout << "OLGModel.wrap() - solution is NaN" << std::endl;
		exit(-1);
	}
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

