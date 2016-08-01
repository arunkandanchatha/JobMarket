#include "OLGModel.h"

OLGModel::OLGModel(unsigned int generations, double y, double s, MatchingFunction &f)
	: m_gens(generations),m_f(&f),m_y(y),m_gamma(1-f.getEta()),m_s(s)
{
}

OLGModel::OLGModel(OLGModel &orig) 
	: m_gens(orig.m_gens), m_f(orig.m_f), m_y(orig.m_y),m_gamma(orig.m_gamma),m_s(orig.m_s)
{
}

OLGModel::~OLGModel()
{
}

void OLGModel::solveWages()
{
	return solveWages(m_gens);
}

void OLGModel::solveWages(unsigned int generations)
{
	if (generations > m_gens) {
		std::cout << "OLGModel.solveWages(generations): Error! found a generation higher than max." << std::endl;
		exit(-1);
	}
	if (generations < 1) {
		std::cout << "OLGModel.solveWages(generations): Error! found a generation < 1." << std::endl;
		exit(-1);
	}

	for (unsigned int i = generations; i > 0; i--) {
		if (i == m_gens) {
			E_vals.push_back((1 - D_BETA)*pow(D_b, D_RHO));
			U_vals.push_back((1 - D_BETA)*pow(D_b, D_RHO));
			W_vals.push_back(0);
			wages.push_back(D_b);
		}
		else {
			double latestU = U_vals.back();
			double latestE = E_vals.back();
			double latestW = W_vals.back();
			double latestWages = wages.back();
			Generation toSolve(*this, &OLGModel::nonLinearWageEquation, latestU, latestE, latestW);
			MySolver solver(toSolve, 1.0E-20);
			double del = solver.solve(latestWages - D_b, m_y);
			double newWages = D_b + del;
			if (newWages <= latestWages) {
				std::cout << "OLGModel-solveWages(): cohort " << i << " has wageT < wage(T+1). HOW?" << std::endl;
				std::cout << "w" << i << "=" << newWages << ", w" << i + 1 << "=" << latestWages << std::endl;
				exit(-1);
			}
			wages.push_back(newWages);
			U_vals.push_back(calcU(latestU, latestE));
			E_vals.push_back(calcE(del, latestU, latestE));
			W_vals.push_back(calcW(del, latestW));
		}
	}
}

double OLGModel::calcU(double Up1, double Ep1) {
	return pow((1 - D_BETA)*pow(D_b, D_RHO) + D_BETA*pow(Up1 + m_f->f()*(Ep1 - Up1), D_RHO), 1.0/D_RHO);
}

double OLGModel::calcE(double delta, double Up1, double Ep1) {
	return pow((1 - D_BETA)*pow(D_b+delta, D_RHO) + D_BETA*pow(Ep1 -m_s*(Ep1 - Up1), D_RHO), 1.0 / D_RHO);
}

double OLGModel::calcW(double delta, double Wp1) {
	return m_y - D_b - delta + D_BETA*((1 - m_s)*Wp1);
}

double OLGModel::nonLinearWageEquation(double x, double Up1, double Ep1, double Wp1) {
	double penalty = 0;
	if (D_b + x > m_y) {
		penalty += 100;
	}
	else if (x < 0) {
		penalty += 100;
	}
	return penalty + ABS(calcE(x, Up1, Ep1) - calcU(Up1, Ep1) 
		- m_gamma / (1 - m_gamma)*calcW(x, Wp1)*partialE_partialDel(x, Up1, Ep1));
}

double OLGModel::partialE_partialDel(double x, double Up1, double Ep1) {
	double insideBracket = (1 - D_BETA)*pow(D_b + x, D_RHO) + D_BETA*(pow(Up1+(1-m_s)*(Ep1-Up1),D_RHO));
	double dE_dDel = pow(insideBracket, 1.0 / D_RHO - 1)*(1-D_BETA)*pow(D_b+x,D_RHO-1);
	return dE_dDel;
}

double OLGModel::expectedW() {
	double total = 0;
	for (unsigned int i = 0; i < m_gens; i++) {
		total += W_vals[i];
	}
	total /= m_gens;
	return total;
}

double OLGModel::elasticityWRTymb() {
	OLGModel thetaChange(*this);
	thetaChange.m_f = m_f->dTheta();
	OLGModel yChange(m_gens,1.0001*m_y,m_s,*m_f);

	thetaChange.solveWages();
	yChange.solveWages();

	double EW = expectedW();
	double num = (m_y - D_b)*(yChange.expectedW() - EW) / (yChange.m_y - m_y);
	double denom = (1 - m_f->getEta())*EW - m_f->getTheta()*
		(thetaChange.expectedW() - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

	delete thetaChange.m_f;
	return num/denom;
}

double OLGModel::elasticityWRTs() {
	OLGModel thetaChange(*this);
	thetaChange.m_f = m_f->dTheta();
	OLGModel sChange(m_gens, m_y, 1.0001*m_s, *m_f);

	thetaChange.solveWages();
	sChange.solveWages();

	double EW = expectedW();
	double num = (m_s)*(sChange.expectedW() - EW) / (sChange.m_s - m_s);
	double denom = (1 - m_f->getEta())*EW - m_f->getTheta()*
		(thetaChange.expectedW() - EW) / (thetaChange.m_f->getTheta() - m_f->getTheta());

	delete thetaChange.m_f;
	return num / denom;
}

void OLGModel::printWages() {
	for (unsigned int i = 0; i < m_gens; i++) {
		std::cout << "Cohort " << i << ": " << wages[m_gens - 1 - i] << std::endl;
	}
	return;
}