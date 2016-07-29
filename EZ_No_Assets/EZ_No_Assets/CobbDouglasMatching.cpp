#include "CobbDouglasMatching.h"

CobbDouglasMatching::CobbDouglasMatching(double fTarget, double etaTarget) : MatchingFunction(etaTarget)
{
	m_theta = exp(log(fTarget) / etaTarget);
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