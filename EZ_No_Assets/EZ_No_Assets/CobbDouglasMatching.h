#pragma once
#include "stdafx.h"
#include "MatchingFunction.h"
#include "adept.h"
using adept::adouble;

class CobbDouglasMatching :
	public MatchingFunction
{
private:
	CobbDouglasMatching(CobbDouglasMatching &orig);
	const double m_mu;
	adouble calculatedF(adouble newTheta) final;

public:
	CobbDouglasMatching(double fTarget, double etaTarget);
	~CobbDouglasMatching();

	double f() final;
	MatchingFunction *dTheta() final;
};

