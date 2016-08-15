#include "deHaanMatching.h"
deHaanMatching::deHaanMatching(double parameter) : MatchingFunction(parameter), m_mu(19.89107045)
{
	m_theta = 0.5;
}

deHaanMatching::deHaanMatching(deHaanMatching &orig) : MatchingFunction(orig.m_parameter), m_mu(orig.m_mu)
{
	m_theta = orig.m_theta;
}


deHaanMatching::~deHaanMatching()
{
}

double deHaanMatching::calculatedF(double x)
{
	return m_mu*x / pow(1 + pow(m_mu*x, m_parameter), 1.0 / m_parameter);
}

#if 0
double deHaanMatching::calculatedF(double x)
{
	return x / pow(1 + pow(x, m_mu), 1.0 / m_mu);
}
#endif

double deHaanMatching::f()
{
	return calculatedF(m_theta);
}

MatchingFunction *deHaanMatching::dTheta()
{
	deHaanMatching *newMatch = new deHaanMatching(*this);
	newMatch->m_theta = 1.0001*m_theta;
	return newMatch;
}
