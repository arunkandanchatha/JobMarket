#include "OLGModel.h"

OLGModel::OLGModel(unsigned int generations, double y, double s, MatchingFunction &f, ShockProcess &sp)
	: m_gens(generations),m_f(&f),m_gamma(1-f.getEta()),m_Es(s),m_sp(&sp),m_Y(sp.numStates()),
	E_vals(generations,sp.numStates()), U_vals(generations,sp.numStates()),
	W_vals(generations,sp.numStates()),wages(generations,sp.numStates()),
	m_thetas(sp.numStates())
{
	for (unsigned int i = 0; i < sp.numStates(); i++) {
		m_thetas(i) = f.getTheta();
	}

	const VectorXd* x = sp.states();
	for (unsigned int i = 0; i < sp.numStates(); i++) {
		m_Y(i) = D_b + exp((*x)(i))*(y - D_b);
	}
}

OLGModel::~OLGModel()
{
}

void OLGModel::solveWages()
{
	for (unsigned int i = m_gens - 1; i < m_gens; i--) {
		lastSolveGen = i + 1;
		for (unsigned int j = 0; j < m_sp->numStates(); j++) {
			if (i == m_gens-1) {
				E_vals(i,j) = (1 - D_BETA)*pow(D_b, D_RHO);
				U_vals(i, j) = (1 - D_BETA)*pow(D_b, D_RHO);
				W_vals(i, j) = 0;
				wages(i, j) = D_b;
			}
			else {
				VectorXd latestU = U_vals.row(i + 1);
				VectorXd latestE = E_vals.row(i + 1);
				VectorXd latestW = W_vals.row(i + 1);
				double latestWages = wages(i + 1, j);

				Generation toSolve(*this, &OLGModel::nonLinearWageEquation, j, latestU, latestE, latestW);
				MySolver solver(toSolve, 1.0E-20);

				double prodInState = m_Y(j);
				double del = solver.solve(latestWages - D_b, prodInState);
				double newWages = D_b + del;
				if (newWages <= latestWages) {
					std::cout << "OLGModel-solveWages(): cohort " << i << " shock " << j << " has wageT <= wage(T+1). HOW?" << std::endl;
					std::cout << "w" << i << "=" << newWages << ", w" << i + 1 << "=" << latestWages << std::endl;
					std::cout << "diff: " << newWages - latestWages << std::endl;
					exit(-1);
				}
				wages(i, j) = newWages;
				U_vals(i,j) = calcU(j, latestU, latestE);
				E_vals(i, j) = calcE(j, del, latestU, latestE);
				W_vals(i, j) = calcW(j, del, latestW);
			}
		}
	}
}

double OLGModel::calcU(int state, VectorXd &Up1, VectorXd &Ep1) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
//	for (unsigned int i = 0; i < m_sp->numStates(); i++) {
	for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		double nextProb = nextPDF(i, 0);
		double calcF = m_f->calculatedF(m_thetas(i));
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
	//	for (unsigned int i = 0; i < m_sp->numStates(); i++) {
	for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		total += nextPDF(i, 0)*pow(Ep1(i) - m_Es*(Ep1(i) - Up1(i)), D_RHO);
	}
	total *= D_BETA;
	return pow((1 - D_BETA)*pow(D_b+delta, D_RHO) + total, 1.0 / D_RHO);
}

double OLGModel::calcW(int state, double delta, VectorXd &Wp1) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	//	for (unsigned int i = 0; i < m_sp->numStates(); i++) {
	for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		total += nextPDF(i, 0)*(1 - m_Es)*Wp1(i);
	}
	return m_Y(state) - D_b - delta + D_BETA*total;
	return 1 - D_b - delta + D_BETA*total;
}

double OLGModel::nonLinearWageEquation(int state, double x, VectorXd& Up1, VectorXd& Ep1, VectorXd& Wp1) {
	double penalty = 0;
	if (D_b + x > m_Y(state)) {
		penalty += 100;
	}
	else if (x < 0) {
		penalty += 100;
	}
	if (lastSolveGen < (m_gens - 1)) {
		if (D_b + x <= wages(lastSolveGen, state)) {
			penalty += 100;
		}
	}

	return penalty + ABS(calcE(state, x, Up1, Ep1) - calcU(state, Up1, Ep1) 
		- m_gamma / (1 - m_gamma)*calcW(state, x, Wp1)*partialE_partialDel(state, x, Up1, Ep1));
}

double OLGModel::partialE_partialDel(int state, double x, VectorXd& Up1, VectorXd& Ep1) {
	double total = 0;

	pdfMatrix& nextPDF = m_sp->nextPeriodPDF(state);
	//	for (unsigned int i = 0; i < m_sp->numStates(); i++) {
	for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
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
	OLGModel thetaChange(m_gens, Ey, m_Es, *(m_f->dTheta()), *m_sp);
	OLGModel yChange(m_gens, 1.0001*Ey, m_Es, *m_f, *m_sp);

	thetaChange.solveWages();
	yChange.solveWages();

	double EW = expectedW(expectedState, true);
	double num = (Ey - D_b)*(yChange.expectedW(expectedState, true) - EW) / (1.0001*Ey - Ey);
	double denom = (1 - m_f->getEta())*EW - m_f->getTheta()*
		(thetaChange.expectedW(expectedState, true) - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

	delete thetaChange.m_f;
	return num/denom;
}

double OLGModel::elasticityWRTs() {
	OLGModel thetaChange(*this);
	thetaChange.m_f = m_f->dTheta();

	unsigned int expectedState = (m_sp->numStates() - 1) / 2;
	OLGModel sChange(m_gens, m_Y(expectedState), 1.0001*m_Es, *m_f, *m_sp);

	thetaChange.solveWages();
	sChange.solveWages();

	double EW = expectedW(expectedState, true);
	double num = (m_Es)*(sChange.expectedW(expectedState, true) - EW) / (sChange.m_Es - m_Es);
	double denom = (1 - m_f->getEta())*EW - m_f->getTheta()*
		(thetaChange.expectedW(expectedState, true) - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

	delete thetaChange.m_f;
	return num / denom;
}

void OLGModel::printWages() {
	for (unsigned int i = 0; i < m_gens; i++) {
		for (unsigned int j = 0; j < m_sp->numStates(); j++) {
			std::cout << "Cohort (" << i << "," << j << "): y=" << m_Y(j) << " b=" << D_b << " w=" << wages(i, j)
				<< "      W: " << W_vals(i,j) << std::endl;
		}
	}
	return;
}

double OLGModel::operator()(const std::vector<double> &x, std::vector<double> &grad)
{
	if (x.size() != m_sp->numStates()) {
		std::cout << "Error! OLGModel.cpp::operator() - numStates and xsize (size of theta being solved) are not equal."
			<< std::endl;
		exit(-1);
	}
	for (unsigned int i = 0; i < x.size(); i++) {
		m_thetas(i) = x[i];
	}
	solveWages();

	double retVal = 0;
	for (int i = 0; i < x.size(); i++) {
//		std::cout << m_thetas[i] << ":" 
//			<< abs(D_C / D_BETA - m_f->calculatedF(m_thetas[i]) / m_thetas[i] * expectedW(i)) << std::endl;
		retVal += abs(D_C / D_BETA - m_f->calculatedF(m_thetas[i]) / m_thetas[i] * expectedW(i));
	}
	std::cout << retVal << std::endl;
	std::cout << "===================" << std::endl;
	return retVal;
}

double OLGModel::wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
	OLGModel modelToSolve = *(OLGModel *)data;
	return (*reinterpret_cast<OLGModel*>(data))(x, grad);
}