#pragma once
#include "stdafx.h"
#include "MatchingFunction.h"
class CobbDouglasMatching :
	public MatchingFunction
{
private:
	CobbDouglasMatching(CobbDouglasMatching &orig);
	const double m_mu;

public:
	CobbDouglasMatching(double fTarget, double etaTarget);
	~CobbDouglasMatching();

	double calculatedF(double newTheta) final;
	double f() final;
	MatchingFunction *dTheta() final;
};

