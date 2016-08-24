#pragma once
#include "stdafx.h"

#define D_MY_MU (15.14)
#define D_MY_ALPHA (1)
#define D_MY_PARAMETER (0.38)

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
	virtual double calculatedF(double newTheta) = 0;
	virtual double getElasticity(double newTheta) = 0;

protected:
	double m_theta;
	const double m_parameter;
};

