#pragma once
#include "MatchingFunction.h"

class deHaanMatching :
	public MatchingFunction
{
	friend class OLGModel;
private:
	deHaanMatching(deHaanMatching &orig);
	const double m_mu;
	const double m_alpha;

public:
	deHaanMatching(double fTarget);
	~deHaanMatching();
	double calculatedF(double newTheta) final;
	double getElasticity(double newTheta) final;

#if 0
	double f() final;
#endif
	MatchingFunction *dTheta() final;
};

