#include "OLGModel.h"

OLGModel::OLGModel(unsigned int generations, double y, MatchingFunction &f) 
	: m_gens(generations),m_f(&f),m_y(y)
{
}

OLGModel::OLGModel(OLGModel &orig) 
	: m_gens(orig.m_gens), m_f(orig.m_f), m_y(orig.m_y)
{
}

OLGModel::~OLGModel()
{
}

void OLGModel::solveWages()
{
	return solveWages(1);
}

void OLGModel::solveWages(unsigned int generations)
{
	if (generations > m_gens) {
		std::cout << "OLGModel.solveWages(generations): Error! found a generation higher than max." << std::endl;
		exit(-1);
	}
	if (generations == m_gens) {
		E_vals.push_back((1 - D_BETA)*pow(D_b, D_RHO));
		U_vals.push_back((1 - D_BETA)*pow(D_b, D_RHO));
		W_vals.push_back(0);
		wages.push_back(D_b);
		return;
	}
	else {
		solveWages(generations + 1);

		//solve for delGen using E_vals, U_vals, and W_vals
		double latestU = U_vals.back();
		double latestE = E_vals.back();
		double latestW = W_vals.back();
		double latestWages = wages.back();
		Generation toSolve(*this, &OLGModel::nonLinearWageEquation, latestU, latestE, latestW);
		MySolver solver(toSolve, 10E-15);
		double del = solver.solve(latestWages-D_b,m_y);
		wages.push_back(D_b + del);
		U_vals.push_back(calcU(latestU, latestE));
		E_vals.push_back(calcE(del, latestU, latestE));
		W_vals.push_back(calcW(del, latestW));
	}
	return;
}

double OLGModel::calcU(double Up1, double Ep1) {
	return pow((1 - D_BETA)*pow(D_b, D_RHO) + D_BETA*pow(Up1 + m_f->f()*(Ep1 - Up1), D_RHO), 1.0/D_RHO);
}

double OLGModel::calcE(double delta, double Up1, double Ep1) {
	return pow((1 - D_BETA)*pow(D_b+delta, D_RHO) + D_BETA*pow(Ep1 -D_S*(Ep1 - Up1), D_RHO), 1.0 / D_RHO);
}

double OLGModel::calcW(double delta, double Wp1) {
	return m_y - D_b - delta + D_BETA*((1 - D_S)*Wp1);
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
		- D_GAMMA / (1 - D_GAMMA)*calcW(x, Wp1)*partialE_partialDel(x, Up1, Ep1));
}

double OLGModel::partialE_partialDel(double x, double Up1, double Ep1) {
	double insideBracket = (1 - D_BETA)*pow(D_b + x, D_RHO) + D_BETA*(pow(Up1+(1-D_S)*(Ep1-Up1),D_RHO));
	double dE_dDel = pow(insideBracket, 1.0 / D_RHO - 1)*(1-D_BETA)*pow(D_b+x,D_RHO-1);
	return dE_dDel;
}


double OLGModel::elasticity(OLGModel &thetaChange, OLGModel &yChange) {
#if 0
	double w11 = D_b + del1;
	double w21 = D_b + del2;
	double w12 = D_b + yChange.del1;
	double w22 = D_b + yChange.del2;
	double w1theta = D_b + thetaChange.del1;
	double w2theta = D_b + thetaChange.del2;

	double y = D_Y;
	double y2 = yChange.D_Y;

	double theta = D_THETA;
	double theta2 = thetaChange.D_THETA;

	double num1 = (2 * y - w11 - w21) / 2 + D_BETA / 2 * (1 - D_S)*(y - w21);
	double num2 = (2 * y2 - w12 - w22) / 2 + D_BETA / 2 * (1 - D_S)*(y2 - w22);
	double partialWRTymb = (num2 - num1) / (y2 - y);

	double numTheta = (2 * y - w1theta - w2theta) / 2 + D_BETA / 2 * (1 - D_S)*(y - w2theta);
	double partialWRTTheta = (numTheta - num1) / (theta2 - theta);

	double numDel2 = 0;// 0.5*(-1 - D_BETA*(1 - D_S)) * (yChange.del2 - del2) / (y2 - y);
	double numDel1 = 0;// -0.5 * (yChange.del1 - del1) / (y2 - y);

	double denom = (1 - D_ETA)*num1 - D_THETA*partialWRTTheta;

	return (y-D_b)*partialWRTymb/denom;
#endif
	return 0;
}

void OLGModel::printWages() {
	for (int i = 0; i < m_gens; i++) {
		std::cout << "Cohort " << i << ": " << wages[m_gens - 1 - i] << std::endl;
	}
	return;
}