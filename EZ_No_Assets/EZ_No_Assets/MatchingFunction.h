#pragma once
class MatchingFunction
{
public:
	MatchingFunction(double etaTarget);
	~MatchingFunction();

	virtual double calculatedF(double newTheta) = 0;
	virtual double f() = 0;
	virtual MatchingFunction *dTheta() = 0;

	double getTheta();
	double getEta();

protected:
	double m_theta;
	const double m_eta;
};

