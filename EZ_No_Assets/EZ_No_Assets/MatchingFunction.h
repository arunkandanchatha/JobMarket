#pragma once
#include "stdafx.h"

#define D_MY_MU (2.57)
#define D_MY_ALPHA (0.54)
#define D_MY_PARAMETER (1.31)

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

