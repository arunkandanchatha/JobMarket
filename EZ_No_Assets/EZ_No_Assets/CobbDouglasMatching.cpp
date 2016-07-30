#include "CobbDouglasMatching.h"

CobbDouglasMatching::CobbDouglasMatching(double fTarget, double etaTarget) : MatchingFunction(etaTarget)
{
	m_theta = exp(log(fTarget) / etaTarget);
}

CobbDouglasMatching::CobbDouglasMatching(CobbDouglasMatching &orig) : MatchingFunction(orig.m_eta)
{
	m_theta = orig.m_theta;
}


CobbDouglasMatching::~CobbDouglasMatching()
{
}

double CobbDouglasMatching::calculatedF(double x)
{
	return pow(x, m_eta);
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
