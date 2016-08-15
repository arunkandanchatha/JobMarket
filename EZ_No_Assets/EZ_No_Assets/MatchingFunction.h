#pragma once
#include "stdafx.h"

class MatchingFunction
{
	friend class OLGModel;
public:
	MatchingFunction(double etaTarget);
	~MatchingFunction();

	virtual double f() = 0;
	virtual MatchingFunction *dTheta() = 0;

	double getTheta();
	double getBargaining();

protected:
	virtual double calculatedF(double newTheta) = 0;
	double m_theta;
	const double m_bargaining;
};

