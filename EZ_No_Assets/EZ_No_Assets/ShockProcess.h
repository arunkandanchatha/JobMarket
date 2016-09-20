#pragma once

#include "stdafx.h"

typedef Matrix<double, Dynamic, 2> pdfMatrix;

class ShockProcess
{
public:
	ShockProcess(unsigned int states);
	~ShockProcess();

	MatrixXd getProbMatrix();

	pdfMatrix nextPeriodPDF(unsigned int state);
	const VectorXd* states();
	int numStates();
	void printStates();

private:
	virtual void calculateConditionalProbabilities() = 0;

protected:
	MatrixXd m_conditionalProbs;
	VectorXd m_values;
};