#include "MatchingFunction.h"



MatchingFunction::MatchingFunction(double parameter) : m_parameter(parameter)
{
}


MatchingFunction::~MatchingFunction()
{
}

double MatchingFunction::getTheta()
{
	return m_theta;
}

double MatchingFunction::getParameter()
{
	return m_parameter;
}
