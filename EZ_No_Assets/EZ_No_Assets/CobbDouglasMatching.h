#pragma once
#include "stdafx.h"
#include "MatchingFunction.h"

class CobbDouglasMatching :
	public MatchingFunction
{
private:
	CobbDouglasMatching(CobbDouglasMatching &orig);
	const double m_mu;
	double calculatedF(double newTheta) final;
	double getElasticity(double newTheta) final;

public:
	CobbDouglasMatching(double fTarget, double etaTarget);
	~CobbDouglasMatching();

	double f() final;
	MatchingFunction *dTheta() final;
};

