#pragma once
#include "ShockProcess.h"
class ShimerProcess :
	public ShockProcess
{
public:
	ShimerProcess(int states, double sigma, double gamma);
	~ShimerProcess();

	pdfMatrix nextPeriodPDF(int state);

private:
	const double m_sigma;
	const double m_gamma;
	VectorXd m_values;

	void calculateConditionalProbabilities();
};

