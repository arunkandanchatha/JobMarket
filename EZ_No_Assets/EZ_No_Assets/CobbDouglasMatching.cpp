#include "CobbDouglasMatching.h"

CobbDouglasMatching::CobbDouglasMatching(double fTarget, double etaTarget) : MatchingFunction(1-etaTarget), m_mu(fTarget)
{
	m_theta = 1; 
}

CobbDouglasMatching::CobbDouglasMatching(CobbDouglasMatching &orig) : MatchingFunction(orig.m_bargaining), m_mu(orig.m_mu)
{
	m_theta = orig.m_theta;
}


CobbDouglasMatching::~CobbDouglasMatching()
{
}

adouble CobbDouglasMatching::calculatedF(adouble x)
{
	return m_mu*pow(x, 1- m_bargaining);
}

double CobbDouglasMatching::f()
{
	return value(calculatedF(m_theta));
}

MatchingFunction *CobbDouglasMatching::dTheta()
{
	CobbDouglasMatching *newMatch = new CobbDouglasMatching(*this);
	newMatch->m_theta = 1.0001*m_theta;
	return newMatch;
}
