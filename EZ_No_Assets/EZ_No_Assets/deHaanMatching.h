#pragma once
#include "MatchingFunction.h"
#include "adept.h"
using adept::adouble;

class deHaanMatching :
	public MatchingFunction
{
	friend class OLGModel;
private:
	deHaanMatching(deHaanMatching &orig);
	const double m_mu;
	adouble calculatedF(adouble newTheta) final;

public:
	deHaanMatching(double fTarget, double etaTarget);
	~deHaanMatching();

	double f() final;
	MatchingFunction *dTheta() final;
};

