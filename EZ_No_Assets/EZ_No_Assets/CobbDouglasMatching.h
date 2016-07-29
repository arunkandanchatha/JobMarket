#pragma once
#include "stdafx.h"
#include "MatchingFunction.h"
class CobbDouglasMatching :
	public MatchingFunction
{
public:
	CobbDouglasMatching(double fTarget, double etaTarget);
	~CobbDouglasMatching();

	double calculatedF(double newTheta) final;
	double f() final;

};

