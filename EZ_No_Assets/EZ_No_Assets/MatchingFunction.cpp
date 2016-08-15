#include "MatchingFunction.h"



MatchingFunction::MatchingFunction(double bargainingPower) : m_bargaining(bargainingPower)
{
}


MatchingFunction::~MatchingFunction()
{
}

double MatchingFunction::getTheta()
{
	return m_theta;
}

double MatchingFunction::getBargaining()
{
	return m_bargaining;
}
