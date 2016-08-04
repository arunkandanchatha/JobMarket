#include "ShimerProcess.h"



ShimerProcess::ShimerProcess(int states, double sigma, double gamma)
	: ShockProcess(2*states+1), m_sigma(sigma), m_gamma(gamma)
{
	calculateConditionalProbabilities();
}


ShimerProcess::~ShimerProcess()
{
}

void ShimerProcess::calculateConditionalProbabilities()
{
	int states = (numStates() - 1) / 2;
	double lambda = states*m_gamma;
	double delta = m_sigma / sqrt(lambda);

	for (unsigned int i = 0; i < numStates(); i++) {
		m_values(i) = ((int)i-states)*delta;
	}

	for (unsigned int i = 0; i < numStates(); i++) {
		for (unsigned int j = 0; j < numStates(); j++) {
			m_conditionalProbs(i, j) = 0;
		}
	}
	for (unsigned int i = 0; i < numStates(); i++) {
		if (i == 0) {
			m_conditionalProbs(i, 1) = 1;
			continue;
		}
		if (i == (numStates() - 1)) {
			m_conditionalProbs(i, numStates() - 2) = 1;
			continue;
		}
		m_conditionalProbs(i, i - 1) = 0.5*(1 - m_values(i) / (states*delta));
		m_conditionalProbs(i, i + 1) = 0.5*(1 + m_values(i) / (states*delta));
	}

	return;
}
