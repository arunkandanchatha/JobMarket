#include "deHaanMatching.h"
deHaanMatching::deHaanMatching(double parameter) : MatchingFunction(D_MY_PARAMETER), m_mu(D_MY_MU), m_alpha(D_MY_ALPHA)
{
	m_theta = 1;
}

deHaanMatching::deHaanMatching(deHaanMatching &orig) : MatchingFunction(orig.m_parameter), m_mu(orig.m_mu), m_alpha(orig.m_alpha)
{
	m_theta = orig.m_theta;
}


deHaanMatching::~deHaanMatching()
{
}

double deHaanMatching::calculatedF(double x)
{
	return m_alpha*m_mu*x / pow(1 + pow(m_mu*x, m_parameter), 1.0 / m_parameter);
}

double deHaanMatching::getElasticity(double x)
{
//	return 1 / (1 + pow(m_mu*x, m_parameter));
	return 0.28;
}

#if 0
double deHaanMatching::calculatedF(double x)
{
	return x / pow(1 + pow(x, m_mu), 1.0 / m_mu);
}
#endif

#if 0
double deHaanMatching::f()
{
	return calculatedF(m_theta);
}
#endif

MatchingFunction *deHaanMatching::dTheta()
{
	deHaanMatching *newMatch = new deHaanMatching(*this);
	newMatch->m_theta = 1.0001*m_theta;
	return newMatch;
}
