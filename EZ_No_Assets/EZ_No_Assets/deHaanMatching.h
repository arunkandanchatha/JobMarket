#pragma once
#include "MatchingFunction.h"

class deHaanMatching :
	public MatchingFunction
{
	friend class OLGModel;
private:
	deHaanMatching(deHaanMatching &orig);
	const double m_mu;
	double calculatedF(double newTheta) final;

public:
	deHaanMatching(double fTarget, double etaTarget);
	~deHaanMatching();

	double f() final;
	MatchingFunction *dTheta() final;
};

