#include "ShimerProcess.h"



ShimerProcess::ShimerProcess(int states, double sigma, double gamma)
	: ShockProcess(2*states+1), m_sigma(sigma), m_gamma(gamma), m_values(m_states)
{
	calculateConditionalProbabilities();
}


ShimerProcess::~ShimerProcess()
{
}

pdfMatrix ShimerProcess::nextPeriodPDF(int state)
{
	pdfMatrix retVal;
	retVal.col(1) = m_conditionalProbs.col(state).transpose();
	retVal.col(2) = m_values;
	return retVal;
}

void ShimerProcess::calculateConditionalProbabilities()
{
	int states = (m_states - 1) / 2;
	double lambda = states*m_gamma;
	double delta = m_sigma / sqrt(lambda);

	for (int i = 0; i < m_states; i++) {
		m_values(i) = (i-states)*delta;
	}

	for (int i = 0; i < m_states; i++) {
		if (i == 0) {
			m_conditionalProbs(i, 1) = 1;
			continue;
		}
		if (i == (m_states - 1)) {
			m_conditionalProbs(i, m_states - 2) = 1;
			continue;
		}
		m_conditionalProbs(i, i - 1) = 0.5*(1 - m_values(i) / (states*delta));
		m_conditionalProbs(i, i + 1) = 0.5*(1 + m_values(i) / (states*delta));
	}
	return;
}
