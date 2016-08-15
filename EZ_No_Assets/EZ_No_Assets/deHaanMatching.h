#pragma once
#include "MatchingFunction.h"
class deHaanMatching :
	public MatchingFunction
{
private:
	deHaanMatching(deHaanMatching &orig);
	const double m_mu;

public:
	deHaanMatching(double fTarget);
	~deHaanMatching();

	double calculatedF(double newTheta) final;
	double f() final;
	MatchingFunction *dTheta() final;
};

