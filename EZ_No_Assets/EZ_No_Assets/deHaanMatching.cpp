#include "deHaanMatching.h"
deHaanMatching::deHaanMatching(double fTarget, double etaTarget) : MatchingFunction(etaTarget), m_mu(fTarget)
{
	m_theta = 2;
}

deHaanMatching::deHaanMatching(deHaanMatching &orig) : MatchingFunction(orig.m_bargaining), m_mu(orig.m_mu)
{
	m_theta = orig.m_theta;
}


deHaanMatching::~deHaanMatching()
{
}

adouble deHaanMatching::calculatedF(adouble x)
{
	return x/pow(1+pow(x, m_mu),1.0/ m_mu);
}

#if 0
double deHaanMatching::calculatedF(double x)
{
	return x / pow(1 + pow(x, m_mu), 1.0 / m_mu);
}
#endif

double deHaanMatching::f()
{
	return value(calculatedF(m_theta));
}

MatchingFunction *deHaanMatching::dTheta()
{
	deHaanMatching *newMatch = new deHaanMatching(*this);
	newMatch->m_theta = 1.0001*m_theta;
	return newMatch;
}
