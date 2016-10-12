#include "OLGModel.h"
#include "MySolver.h"
#include "nlopt.hpp"
#include <exception>

OLGModel::OLGModel(int generations, double y, double s, MatchingFunction &f, ShockProcess &sp, double bargaining, bool autoDiff) :
#ifdef DO_TENURE_SOLVE
	E_vals(generations), U_vals(generations),
	W_vals(generations), wages(generations),
	m_Y(generations),
#else
	E_vals(1), U_vals(1),
	W_vals(1), wages(1),
	m_Y(1),
#endif
	m_gens(generations), m_bargaining(bargaining), 
	m_thetas(sp.numStates()), m_autodiff(autoDiff), m_Es(s), m_sp(&sp), m_f(&f),
	habits(generations,2*generations+1), habitProb(generations, 2 * generations + 1),
	oldWages(generations)
{
	if (generations < 3) {
		std::cout << "OLGModel.constructor() - cannot pass fewer than " << 3 << " generations. gens=" << generations << std::endl;
		exit(-1);
	}

	for (int i = 0; i < generations; i++) {
		oldWages[i].resize(WAGE_GRID_SIZE);
	}

	for (int i = 0; i < sp.numStates(); i++) {
		m_thetas[i] = f.getTheta();
	}

	for (int i = 0; i < habits.rows(); i++) {
		for (int j = 0; j < habits.cols(); j++) {
			habits(i, j) = 0;
			habitProb(i, j) = 0;
		}
	}

	double probSuccess = 1 - D_S;
	double habitIncrement = D_b / 4 / generations;
	for (int trialIndex = 0; trialIndex < generations; trialIndex++) {
		for (int successIndex = 0; successIndex <= trialIndex; successIndex++) {
			int numFailures = trialIndex - successIndex;
			habits(trialIndex, generations + successIndex - numFailures) = D_b / 2 + (successIndex - numFailures)*habitIncrement;
			habitProb(trialIndex, generations + successIndex - numFailures) = nCr(trialIndex, successIndex)*pow(probSuccess, successIndex)*pow(1 - probSuccess, numFailures);
		}
	}

	int solveGenerations = 1;

#ifdef DO_TENURE_SOLVE
	solveGenerations = generations;
#endif

	const VectorXd* x = sp.states();
	for (int genIndex = 0; genIndex < solveGenerations; genIndex++) {
		m_Y[genIndex] = MatrixXd(generations, sp.numStates());
		U_vals[genIndex].resize(2 * generations + 1);
		E_vals[genIndex].resize(2 * generations + 1);
		W_vals[genIndex].resize(2 * generations + 1);
		wages[genIndex].resize(2 * generations + 1);
		for (int habitIndex = 0; habitIndex < 2 * generations + 1; habitIndex++) {
			U_vals[genIndex][habitIndex].resize(WAGE_GRID_SIZE);
			E_vals[genIndex][habitIndex].resize(WAGE_GRID_SIZE);
			W_vals[genIndex][habitIndex].resize(WAGE_GRID_SIZE);
			wages[genIndex][habitIndex].resize(WAGE_GRID_SIZE);
			for (int wageIndex2 = 0; wageIndex2 < WAGE_GRID_SIZE; wageIndex2++) {
				U_vals[genIndex][habitIndex][wageIndex2] = MatrixXd(generations, sp.numStates());
				E_vals[genIndex][habitIndex][wageIndex2] = MatrixXd(generations, sp.numStates());
				W_vals[genIndex][habitIndex][wageIndex2] = MatrixXd(generations, sp.numStates());
				wages[genIndex][habitIndex][wageIndex2] = MatrixXd(generations, sp.numStates());
			}
		}
		for (int genIndex2 = 0; genIndex2 < generations; genIndex2++) {
			for (int i = 0; i < sp.numStates(); i++) {
				m_Y[genIndex](genIndex2, i) = D_b + exp((*x)(i))*((pow(D_TENURE_INCREASE, genIndex))*pow(D_PROD_INCREASE, genIndex2)*y - D_b);
			}
			for (int i = 0; i < WAGE_GRID_SIZE; i++) {
				oldWages[genIndex2][i] = ((double)i) / (WAGE_GRID_SIZE - 1)*(m_Y[genIndex](genIndex2, sp.numStates() - 1) - D_b);
			}
		}
	}

	//error check
	bool foundError = false;
	for (int genIndex = 0; genIndex < generations; genIndex++) {
		for (int genIndex2 = 0; genIndex2 < solveGenerations; genIndex2++) {
			for (int i = 0; i < sp.numStates(); i++) {
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

double OLGModel::adjustmentCost(const double original, const double updateVal) {
	if (original == D_b) {
		return 0;
	}
	return 0;// ABS(updateVal - original);
}

void OLGModel::solveWages()
{
	for (int tenureIndex = 0; tenureIndex < m_Y.size(); tenureIndex++) {
		for (int i = m_gens - 1; i >= 0; i--) {
#pragma omp parallel for num_threads(3)
			for (int j = 0; j < m_sp->numStates(); j++) {
				for (int habitIndex = 0; habitIndex < habits.cols(); habitIndex++) {
					for (int wageIndex = 0; wageIndex < WAGE_GRID_SIZE; wageIndex++) {
						if (i == m_gens - 1) {
							E_vals[tenureIndex][habitIndex][wageIndex](i, j) = pow(1 - D_BETA, 1.0 / D_RHO)*(D_b - (1.0*habitIndex) / WAGE_GRID_SIZE * D_b / 2);
							U_vals[tenureIndex][habitIndex][wageIndex](i, j) = pow(1 - D_BETA, 1.0 / D_RHO)*(D_b - (1.0*habitIndex) / WAGE_GRID_SIZE * D_b / 2);
							W_vals[tenureIndex][habitIndex][wageIndex](i, j) = 0;
							wages[tenureIndex][habitIndex][wageIndex](i, j) = D_b;
						}
						else {
							double prodInState = m_Y[tenureIndex](i, j);

							Generation toSolve(*this, &OLGModel::nonLinearWageEquation, i, j, habitIndex, tenureIndex, wageIndex);
							nlopt::opt opt(nlopt::LN_BOBYQA, 1);
							double lwbnd = 0;
							double upbnd = prodInState - D_b;

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

							if (m_Y[tenureIndex](i, j) - newWages < 0) {
								std::cout << "OLGModel-solveWages(): cohort " << i << " shock " << j << " tenure " << tenureIndex << " has y-w<0. HOW?" << std::endl;
								std::cout << "y: " << m_Y[tenureIndex](i, j) << "   w: " << newWages << "   y-w(T+1): " << m_Y[tenureIndex](i, j) - newWages << std::endl;
								for (int lll = m_gens - 1; lll > i; lll--) {
									std::cout << "y-w(" << lll << "): " << m_Y[tenureIndex](lll, j) - wages[tenureIndex][habitIndex][wageIndex](lll, j) << std::endl;
								}
								double val1 = nonLinearWageEquation(i, j, habitIndex, tenureIndex, wageIndex, del);
								double val2 = nonLinearWageEquation(i, j, habitIndex, tenureIndex, wageIndex, m_Y[tenureIndex](i, j) - D_b);
								std::cout << "value: " << val1 << std::endl;
								std::cout << "value0: " << val2 << std::endl;
								std::cout << "diff: " << val2 - val1 << std::endl;
								exit(-1);
							}

							wages[tenureIndex][habitIndex][wageIndex](i, j) = newWages;
							U_vals[tenureIndex][habitIndex][wageIndex](i, j) = calcU(i, j, habitIndex, tenureIndex, wageIndex);
							E_vals[tenureIndex][habitIndex][wageIndex](i, j) = calcE(i, j, habitIndex, tenureIndex, wageIndex, del);
							W_vals[tenureIndex][habitIndex][wageIndex](i, j) = calcW(i, j, habitIndex, tenureIndex, wageIndex, del);
						}
					}
				}
			}
		}
	}
	return;
}

double OLGModel::calcU(int generation, int state, int habit, int tenure, int wageLastPeriod) {

#ifndef DO_TENURE_SOLVE
	if (tenure != 0) {
		std::cout << "Error. OLGModel.cpp-calcU(): Not solving for tenure. Should not get tenure input != 0. tenureInput="
			<< tenure << std::endl;
		exit(-1);
	}
#endif

	double total = 0;
	int whichGen = generation;
	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	int maxIters = MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1);
	int nextHabit = MAX(habit - 1, 0);
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
		double nextProb = nextPDF(i, 0);
		double calcF = m_f->calculatedF(m_thetas[state]);
		double deathPart = D_DEATH*(U_vals[TENURE_INCREASE_UU(tenure)][nextHabit][0](m_gens - 1, i));
		double nextVal = (1-D_DEATH)*(pow(U_vals[TENURE_INCREASE_UU(tenure)][nextHabit][0](whichGen+1, i)
			+ calcF*(E_vals[TENURE_INCREASE_UE(tenure)][nextHabit][0](whichGen+1,i)
				- U_vals[TENURE_INCREASE_UU(tenure)][nextHabit][0](whichGen + 1, i)), D_RHO));
		total += nextProb*(deathPart+nextVal);
	}
	total *= D_BETA;
	double currentConsUtil = (1 - D_BETA)*pow(D_b-habits(generation,habit), D_RHO);
	double retVal = pow(currentConsUtil + total, 1.0 / D_RHO);
	return retVal;
}

double OLGModel::calcE(int generation, int state, int habit, int tenure, int wageLastPeriod, double delta) {
#ifndef DO_TENURE_SOLVE
	if (tenure != 0) {
		std::cout << "Error. OLGModel.cpp-calcE(): Not solving for tenure. Should not get tenure input != 0. tenureInput="
			<< tenure << std::endl;
		exit(-1);
	}
#endif

	double total = 0;
	int whichGen = generation;
	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	int maxIters = MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1);

	std::vector<std::vector<double>> eVec(maxIters);
	for (int i = 0; i < maxIters; i++) {
		eVec[i].resize(WAGE_GRID_SIZE);
	}

	int nextHabit = MIN(habit + 1, 2*m_gens);
	int zIndex = 0;
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
		for (int j = 0; j < WAGE_GRID_SIZE; j++) {
			eVec[zIndex][j] = E_vals[TENURE_INCREASE_EE(tenure)][nextHabit][j](whichGen + 1, i);
		}
		zIndex++;
	}

	zIndex = 0;
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
		double deathPart = D_DEATH*pow(1 - D_BETA, 1.0 / D_RHO)*(D_b - habits(whichGen, nextHabit));
		//double wageNow = D_b + delta;

		//interpolate E value using wageNow
		double nextE = utilities::interpolate(oldWages[whichGen], eVec[zIndex], delta);
		zIndex++;
		double nextU = U_vals[TENURE_INCREASE_EU(tenure)][nextHabit][0](whichGen + 1, i);

		double nextVal = (1 - D_DEATH)*(pow(nextE - m_Es*(nextE - nextU), D_RHO));
		total += nextPDF(i, 0)*(deathPart+nextVal);
	}
	total *= D_BETA;
	double retVal = pow((1 - D_BETA)*pow(D_b + delta - habits(generation, habit), D_RHO) + total, 1.0 / D_RHO);
	return retVal;
}

double OLGModel::calcW(int generation, int state, int habit, int tenure, int wageLastPeriod, double delta) {
#ifndef DO_TENURE_SOLVE
	if (tenure != 0) {
		std::cout << "Error. OLGModel.cpp-calcW(): Not solving for tenure. Should not get tenure input != 0. tenureInput="
			<< tenure << std::endl;
		exit(-1);
	}
#endif

	double total = 0;
	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	int maxIters = MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1);
	int whichGen = generation;

	std::vector<std::vector<double>> wVec(maxIters);
	for (int i = 0; i < maxIters; i++) {
		wVec[i].resize(WAGE_GRID_SIZE);
	}

	int zIndex = 0;
	int nextHabit = MIN(habit + 1, 2 * m_gens);
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
		for (int j = 0; j < WAGE_GRID_SIZE; j++) {
			wVec[zIndex][j] = W_vals[TENURE_INCREASE_EE(tenure)][nextHabit][j](whichGen + 1, i);
		}
		zIndex++;
	}

	zIndex = 0;
	//double wageThisPeriod = D_b + delta;
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
		//interpolate W values using wageNow
		double nextW = utilities::interpolate(oldWages[whichGen], wVec[zIndex],	delta);
		zIndex++;
		total += nextPDF(i, 0)*(1 - m_Es - D_DEATH)*nextW;
	}
	double retVal = m_Y[tenure](whichGen, state) - D_b - delta 
		- adjustmentCost(oldWages[whichGen][wageLastPeriod],delta) + D_BETA*total;
	return retVal;
}

double OLGModel::nonLinearWageEquation(int generation, int state, int habit, int tenure, int wageLastPeriod, double x) {

#ifndef DO_TENURE_SOLVE
	if (tenure != 0) {
		std::cout << "Error. OLGModel.cpp-nonLinearWageEquation(): Not solving for tenure. "
			<<"Should not get tenure input != 0. tenureInput="
			<< tenure << std::endl;
		exit(-1);
	}
#endif

	int whichGen = generation;
	double t_calcE = calcE(generation, state, habit, tenure, wageLastPeriod, x);
	double t_calcU = calcU(generation, state, habit, tenure, wageLastPeriod);
	double t_calcW = calcW(generation, state, habit, tenure, wageLastPeriod, x);
	double t_partialE_partialDel = partialE_partialDel(generation, state, habit, tenure, wageLastPeriod, x);

	if (m_bargaining == 1) {
		return -(t_calcE - t_calcU);
	}
	double retVal = ABS(
		t_calcE - t_calcU - m_bargaining / (1 - m_bargaining)*t_calcW*t_partialE_partialDel
		);

	if (retVal < 0) {
		std::cout << "OLGModel.cpp-nonLinearWageEquation(): return value < 0. How is this possible?" << std::endl;
		exit(-1);
	}
	return retVal;
}

double OLGModel::partialE_partialDel(int generation, int state, int habit, int tenure, int wageLastPeriod, double delta) {
#ifndef DO_TENURE_SOLVE
	if (tenure != 0) {
		std::cout << "Error. OLGModel.cpp-partialE_partialDel(): Not solving for tenure."
			<<" Should not get tenure input != 0. tenureInput="
			<< tenure << std::endl;
		exit(-1);
	}
#endif

	double total = 0;
	int whichGen = generation;
	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	int maxIters = MIN(m_sp->numStates(), state + MAX_SHOCKS_PER_MONTH + 1);

	std::vector<std::vector<double>> eVec(maxIters);
	for (int i = 0; i < maxIters; i++) {
		eVec[i].resize(WAGE_GRID_SIZE);
	}

	int nextHabit = MIN(habit + 1, 2 * m_gens);
	int zIndex = 0;
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
		for (int j = 0; j < WAGE_GRID_SIZE; j++) {
			eVec[zIndex][j] = E_vals[TENURE_INCREASE_EE(tenure)][nextHabit][j](whichGen + 1, i);
		}
		zIndex++;
	}

	zIndex = 0;
	for (int i = MAX(0, state - MAX_SHOCKS_PER_MONTH); i < maxIters; i++) {
//		double wageNow = D_b + delta;

		//interpolate E value using wageNow
		double nextE = utilities::interpolate(oldWages[whichGen], eVec[zIndex], delta);
		zIndex++;
		double nextU = U_vals[TENURE_INCREASE_EU(tenure)][nextHabit][0](whichGen + 1, i);

		double deathPart = D_DEATH*pow(1 - D_BETA, 1.0 / D_RHO)*(D_b - habits(whichGen, nextHabit));
		double nextVal = (1 - D_DEATH)*(pow(nextE - m_Es*(nextE - nextU), D_RHO));
		total += nextPDF(i, 0)*(deathPart+nextVal);
	}
	total *= D_BETA;

	double insideBracket = (1 - D_BETA)*pow(D_b + delta - habits(generation, habit), D_RHO) + total;
	double dE_dDel = pow(insideBracket, 1.0 / D_RHO - 1)*(1 - D_BETA)*pow(D_b + delta - habits(generation, habit), D_RHO - 1);
	return dE_dDel;
}

VectorXd OLGModel::getSteadyStateDistrib(double jobFind) {
	if (jobFind > 1) {
		std::cout << "OLGModel.cpp-getSteadStateDistrib() - probability of finding a job must be <= 1" << std::endl;
		exit(-1);
	}
	if (jobFind < 0) {
		std::cout << "OLGModel.cpp-getSteadStateDistrib() - probability of finding a job must be >= 0" << std::endl;
		exit(-1);
	}
	MatrixXd trans(m_gens, m_gens);
	VectorXd startDistrib(m_gens);
	for (int i = 0; i < m_gens; i++) {
		for (int j = 0; j < m_gens; j++) {
			trans(i, j) = 0;
		}
		if (i == 0) {
			trans(i, i) = 1-jobFind;
			trans(i, i + 1) = jobFind;
		}
		else if (i == (m_gens - 1)) {
			trans(i, 0) = (1-jobFind)*D_S;
			trans(i, 1) = jobFind*D_S;
			trans(i, i) = 1 - D_S;
		}
		else {
			trans(i, 0) = (1 - jobFind)*D_S;
			trans(i, 1) = jobFind*D_S;
			trans(i, i + 1) = 1 - D_S;
		}
		startDistrib(i) = 0;
	}
	startDistrib(0) = 1;
	MatrixXd transT = trans.transpose();

	VectorXd mystat(m_gens);
	mystat = transT*startDistrib;
	for (int i = 0; i < 1000 * m_gens; i++) {
		mystat = transT*mystat;
	}
	return mystat;
}

double OLGModel::expectedW0(int state, bool forceNoShocks) {
	double total = 0;
	double denom = (1 - pow(D_DEATH, m_gens)) / (1 - D_DEATH);
	if (forceNoShocks || (m_sp->numStates() == 1)) {
		for (int i = 0; i < m_gens-1; i++) {
			for (int habitIndex = m_gens - 1 - i; habitIndex <= m_gens - 1 + i; habitIndex++) {
				double wval = W_vals[0][habitIndex][0](i, state);
				double myChange = (pow(1 - D_DEATH, i) / denom)*habitProb(i, habitIndex)*wval;
				total += myChange;
			}
		}
		return total;
	}

	pdfMatrix temp = m_sp->nextPeriodPDF(state);
	for (int j = 0; j < m_sp->numStates(); j++) {
		double tempVal = 0;
		for (int i = 0; i < m_gens-1; i++) {
			for (int habitIndex = m_gens - 1 - i; habitIndex <= m_gens - 1 + i; habitIndex++) {
				double wval = W_vals[0][habitIndex][0](i, j);
				double myChange = (pow(1 - D_DEATH, i) / denom)*habitProb(i, habitIndex)*wval;
				tempVal += myChange;
			}
		}
		total += tempVal * temp(j, 0);
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
	double dEy = 0;
	double denom1 = (1 - pow(D_DEATH, m_gens)) / (1 - D_DEATH);
	for (int i = 0; i < m_gens-1; i++) {
		Ey += (pow(1 - D_DEATH, i) / denom1)*m_Y[0](i, expectedState);
		dEy += (pow(1 - D_DEATH, i) / denom1)*
			(yChange.m_Y[0](i, expectedState) - m_Y[0](i, expectedState));
	}
	Eb = D_b;

	thetaChange.solveWages();
	yChange.solveWages();
	//yChange.solveWithWages(wages);

	double EW = expectedW0(expectedState, true);
	double tEW = thetaChange.expectedW0(expectedState, true);
	double yEW = yChange.expectedW0(expectedState, true);

	double num = (Ey - Eb)*(yEW - EW) / (dEy);

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
		for (int j = 0; j < m_sp->numStates(); j++) {
			for (int habitIndex = 0; habitIndex < habits.row(0).size(); habitIndex++) {
#ifdef DO_TENURE_SOLVE
				for (int k = 0; k <= i; k++)
#else
				for (int k = 0; k < 1; k++)
#endif
				{
					for (int wageIndex = 0; wageIndex < WAGE_GRID_SIZE; wageIndex++) {
						std::cout << "Cohort_" << i << " (" << j << "," << habitIndex << ","
							<< k << "," << wageIndex <<"): y=" << m_Y[k](i, j) << " b=" << D_b
							<< " w=" << wages[k][habitIndex][wageIndex](i, j) << std::endl;
					}
				}
			}
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
			retVal += pow(D_C / D_BETA - m_f->calculatedF(m_thetas[i]) / m_thetas[i] * expectedW0(i), 2);
		}
	}
	else {
		std::cout << "ERROR! OLGModel.cpp-operator() - need to fix autodiff to account for different productivity." << std::endl;
		exit(-1);
#if 0
		std::vector<double> bargaining(m_Y.size());
		for (int i = 0; i < x.size(); i++) {
			bargaining[i] = m_bargaining;// 1 - m_f->getElasticity(m_thetas[i]);
		}
		OLGSolveAutoDiff soln(m_gens, m_Y, m_sp->getProbMatrix(), m_f->getParameter()/*, bargaining*/, m_Es, wages);

		std::vector<double> myThetas(m_thetas.size());
		VectorXd::Map(&myThetas[0], m_thetas.size()) = m_thetas;
		retVal = soln.solveProblem(myThetas, grad, bargaining);
#endif
	}
	return sqrt(retVal);
}

OLGSolveAutoDiff OLGModel::getSolver() {
	std::cout << "OLGModel::getSolver() not implemented yet." << std::endl;
	exit(-1);
	std::vector<MatrixXd> temp(2);
	MatrixXd temp2;
	std::vector<MatrixXd> temp3(2);
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

#if 0
void OLGModel::printWageElasticity() {
	unsigned int expectedState = (m_sp->numStates() - 1) / 2;
	double Ey = m_Y[0](0, expectedState);
	OLGModel yChange(m_gens, 1.0001*(Ey-D_b)+D_b, m_Es, *m_f, *m_sp, m_bargaining, false);

	yChange.solveWages();

	std::cout << "generation,tenure,z,y,y-b,w,w-b,y',y'-b,w',w'-b,elast" << std::endl;
	for (int gen_i = 0; gen_i < m_gens-1; gen_i++) {
		for (int tenure_i = 0; tenure_i <= ((D_TENURE_INCREASE == 1) ? 0 : gen_i); tenure_i++) {
			for (int shock_i = 0; shock_i < m_sp->numStates(); shock_i++) {
				std::cout << gen_i << "," << tenure_i << "," << shock_i << ",";
				std::cout << m_Y[tenure_i](gen_i, shock_i);
				std::cout << "," << (m_Y[tenure_i](gen_i, shock_i) - D_b) << "," << wages[tenure_i](gen_i, shock_i);
				std::cout << "," << (wages[tenure_i](gen_i, shock_i) - D_b) << ",";
				std::cout << yChange.m_Y[tenure_i](gen_i, shock_i);
				std::cout << "," << (yChange.m_Y[tenure_i](gen_i, shock_i) - D_b) << "," << yChange.wages[tenure_i](gen_i, shock_i);
				std::cout << "," << (yChange.wages[tenure_i](gen_i, shock_i) - D_b) << ",";
				double elast = (yChange.wages[tenure_i](gen_i, shock_i) - D_b) - (wages[tenure_i](gen_i, shock_i) - D_b);
				elast /= (yChange.m_Y[tenure_i](gen_i, shock_i) - D_b) - (m_Y[tenure_i](gen_i, shock_i) - D_b);
				elast *= ((m_Y[tenure_i](gen_i, shock_i) - D_b)/ (wages[tenure_i](gen_i, shock_i) - D_b));
				std::cout << elast << std::endl;
			}
		}
	}
}
#endif