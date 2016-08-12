#pragma once
#include "stdafx.h"
#include "adept.h"
using adept::adouble;

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
	virtual adouble calculatedF(adouble newTheta) = 0;
	adouble m_theta;
	const double m_bargaining;
};

