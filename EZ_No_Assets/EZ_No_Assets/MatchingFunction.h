#pragma once
#include "stdafx.h"
class MatchingFunction
{
public:
	MatchingFunction(double parameter);
	~MatchingFunction();

	virtual double calculatedF(double newTheta) = 0;
	virtual double f() = 0;
	virtual MatchingFunction *dTheta() = 0;

	double getTheta();
	double getParameter();

protected:
	double m_theta;
	const double m_parameter;
};

