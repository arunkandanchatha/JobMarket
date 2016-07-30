#include "CobbDouglasMatching.h"

CobbDouglasMatching::CobbDouglasMatching(double fTarget, double etaTarget) : MatchingFunction(etaTarget), m_mu(fTarget)
{
	m_theta = 1; 
}

CobbDouglasMatching::CobbDouglasMatching(CobbDouglasMatching &orig) : MatchingFunction(orig.m_eta), m_mu(orig.m_mu)
{
	m_theta = orig.m_theta;
}


CobbDouglasMatching::~CobbDouglasMatching()
{
}

double CobbDouglasMatching::calculatedF(double x)
{
	return m_mu*pow(x, m_eta);
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
