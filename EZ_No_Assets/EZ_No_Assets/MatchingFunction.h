#pragma once
#include "stdafx.h"

class MatchingFunction
{
	friend class OLGModel;
public:
	MatchingFunction(double parameter);
	~MatchingFunction();

	virtual double f() = 0;
	virtual MatchingFunction *dTheta() = 0;

	double getTheta();
	double getParameter();

protected:
	virtual double calculatedF(double newTheta) = 0;
	double m_theta;
	const double m_parameter;
};

