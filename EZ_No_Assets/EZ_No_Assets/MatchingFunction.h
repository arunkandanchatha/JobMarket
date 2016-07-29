#pragma once
class MatchingFunction
{
public:
	MatchingFunction(double etaTarget);
	~MatchingFunction();

	virtual double calculatedF(double newTheta) = 0;
	virtual double f() = 0;

	double getTheta();

protected:
	double m_theta;
	const double m_eta;
};

