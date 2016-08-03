#pragma once

#include <Eigen/Dense>
using namespace Eigen;

typedef Matrix<double, Dynamic, 2> pdfMatrix;

class ShockProcess
{
public:
	ShockProcess(int states);
	~ShockProcess();

	virtual pdfMatrix nextPeriodPDF(int state) = 0;

private:
	virtual void calculateConditionalProbabilities() = 0;

protected:
	const int m_states;
	MatrixXd m_conditionalProbs;
};