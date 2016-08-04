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

				pdfMatrix nextPDF = m_sp->nextPeriodPDF(j);

				Generation toSolve(*this, &OLGModel::nonLinearWageEquation, j, latestU, latestE, latestW, nextPDF);
				MySolver solver(toSolve, 1.0E-20);

				double prodInState = m_Y(j);
				double del = solver.solve(latestWages - D_b, prodInState);
				double newWages = D_b + del;
				if (newWages <= latestWages) {
					std::cout << "OLGModel-solveWages(): cohort " << i << " shock " << j << " has wageT <= wage(T+1). HOW?" << std::endl;
					std::cout << "w" << i << "=" << newWages << ", w" << i + 1 << "=" << latestWages << std::endl;
					std::cout << newWages - latestWages << std::endl;
					exit(-1);
				}
				wages(i, j) = newWages;
				U_vals(i,j) = calcU(j, latestU, latestE, nextPDF);
				E_vals(i, j) = calcE(j, del, latestU, latestE, nextPDF);
				W_vals(i, j) = calcW(j, del, latestW, nextPDF);
			}
		}
	}
}

double OLGModel::calcU(int state, VectorXd &Up1, VectorXd &Ep1, pdfMatrix& nextPDF) {
	double total = 0;

	for (unsigned int i = 0; i < m_sp->numStates(); i++) {
//	for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
	total += nextPDF(i, 0)*pow(Up1(i) + m_f->calculatedF(m_thetas(i))*(Ep1(i) - Up1(i)), D_RHO);
	}
	total *= D_BETA;
	return pow((1 - D_BETA)*pow(D_b, D_RHO) + total, 1.0/D_RHO);
}

double OLGModel::calcE(int state, double delta, VectorXd &Up1, VectorXd &Ep1, pdfMatrix& nextPDF) {
	double total = 0;

	for (unsigned int i = 0; i < m_sp->numStates(); i++) {
		//	for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		total += nextPDF(i, 0)*pow(Ep1(i) - m_Es*(Ep1(i) - Up1(i)), D_RHO);
	}
	total *= D_BETA;
	return pow((1 - D_BETA)*pow(D_b+delta, D_RHO) + total, 1.0 / D_RHO);
}

double OLGModel::calcW(int state, double delta, VectorXd &Wp1, pdfMatrix& nextPDF) {
	double total = 0;

	for (unsigned int i = 0; i < m_sp->numStates(); i++) {
		//	for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		total += nextPDF(i, 0)*(1 - m_Es)*Wp1(i);
	}
	return m_Y(state) - D_b - delta + D_BETA*total;
	return 1 - D_b - delta + D_BETA*total;
}

double OLGModel::nonLinearWageEquation(int state, double x, VectorXd& Up1, VectorXd& Ep1, VectorXd& Wp1, pdfMatrix& nextPDF) {
	double penalty = 0;
	if (D_b + x > m_Y(state)) {
		penalty += 100;
	}
	else if (x < 0) {
		penalty += 100;
	}

	return penalty + ABS(calcE(state, x, Up1, Ep1, nextPDF) - calcU(state, Up1, Ep1, nextPDF) 
		- m_gamma / (1 - m_gamma)*calcW(state, x, Wp1, nextPDF)*partialE_partialDel(state, x, Up1, Ep1, nextPDF));
}

double OLGModel::partialE_partialDel(int state, double x, VectorXd& Up1, VectorXd& Ep1, pdfMatrix& nextPDF) {
	double total = 0;

	for (unsigned int i = 0; i < m_sp->numStates(); i++) {
		//	for (unsigned int i = ((state==0)?0:(state-1)); i < ((state == (m_sp->numStates()-1)) ? m_sp->numStates() : (state + 2)); i++) {
		total += nextPDF(i, 0)*pow(Ep1(i) - m_Es*(Ep1(i) - Up1(i)), D_RHO);
	}
	total *= D_BETA;


	double insideBracket = (1 - D_BETA)*pow(D_b + x, D_RHO) + total;
	double dE_dDel = pow(insideBracket, 1.0 / D_RHO - 1)*(1-D_BETA)*pow(D_b+x,D_RHO-1);
	return dE_dDel;
}

double OLGModel::expectedW(int state) {
	double total = 0;
	for (unsigned int i = 0; i < m_gens; i++) {
		total += W_vals(i,state);
	}
	total /= m_gens;
	return total;
}

double OLGModel::elasticityWRTymb() {
	unsigned int expectedState = (m_sp->numStates() - 1) / 2;
	double Ey = m_Y(expectedState);
	OLGModel thetaChange(m_gens, Ey, m_Es, *(m_f->dTheta()), *m_sp);
	OLGModel yChange(m_gens, 1.0001*Ey, m_Es, *m_f, *m_sp);

	thetaChange.solveWages();
	yChange.solveWages();

	double EW = expectedW(expectedState);
	double num = (Ey - D_b)*(yChange.expectedW(expectedState) - EW) / (1.0001*Ey - Ey);
	double denom = (1 - m_f->getEta())*EW - m_f->getTheta()*
		(thetaChange.expectedW(expectedState) - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

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

	double EW = expectedW(expectedState);
	double num = (m_Es)*(sChange.expectedW(expectedState) - EW) / (sChange.m_Es - m_Es);
	double denom = (1 - m_f->getEta())*EW - m_f->getTheta()*
		(thetaChange.expectedW(expectedState) - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

	delete thetaChange.m_f;
	return num / denom;
}

void OLGModel::printWages() {
	for (unsigned int i = 0; i < m_gens; i++) {
		for (unsigned int j = 0; j < m_sp->numStates(); j++) {
			std::cout << "Cohort (" << i << "," << j << "): " << wages(m_gens - 1 - i, j) << std::endl;
		}
	}
	return;
}