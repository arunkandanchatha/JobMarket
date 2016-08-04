#include "NoShocksProcess.h"



NoShocksProcess::NoShocksProcess():ShockProcess(1)
{
	calculateConditionalProbabilities();
}


NoShocksProcess::~NoShocksProcess()
{
}

void NoShocksProcess::calculateConditionalProbabilities()
{
	m_values(0) = 0;
	m_conditionalProbs(0, 0) = 1;
	return;
}
