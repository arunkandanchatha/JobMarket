#include "MatchingFunction.h"



MatchingFunction::MatchingFunction(double etaTarget) : m_eta(etaTarget)
{
}


MatchingFunction::~MatchingFunction()
{
}

double MatchingFunction::getTheta()
{
	return m_theta;
}
