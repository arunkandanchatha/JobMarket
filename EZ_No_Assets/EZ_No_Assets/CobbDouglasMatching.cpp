#include "CobbDouglasMatching.h"

CobbDouglasMatching::CobbDouglasMatching(double fTarget, double parameter) : MatchingFunction(parameter), m_mu(fTarget)
{
	m_theta = 1; 
}

CobbDouglasMatching::CobbDouglasMatching(CobbDouglasMatching &orig) : MatchingFunction(orig.m_parameter), m_mu(orig.m_mu)
{
	m_theta = orig.m_theta;
}


CobbDouglasMatching::~CobbDouglasMatching()
{
}

double CobbDouglasMatching::getElasticity(double x)
{
	return m_parameter;
}


double CobbDouglasMatching::calculatedF(double x)
{
	return m_mu*pow(x, m_parameter);
}

double CobbDouglasMatching::f()
{
	return calculatedF(m_theta);
}

MatchingFunction *CobbDouglasMatching::dTheta()
{
	CobbDouglasMatching *newMatch = new CobbDouglasMatching(*this);
	newMatch->m_theta = 1.0001*m_theta;
	return newMatch;
}
