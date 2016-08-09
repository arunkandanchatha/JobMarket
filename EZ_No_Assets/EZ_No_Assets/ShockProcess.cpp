#include "ShockProcess.h"



ShockProcess::ShockProcess(unsigned int states) : m_conditionalProbs(states, states), m_values(states)
{
}


ShockProcess::~ShockProcess()
{
}

pdfMatrix ShockProcess::nextPeriodPDF(unsigned int state)
{
	pdfMatrix retVal(m_values.size(),2);
	retVal.col(0) = m_conditionalProbs.row(state).transpose();
	retVal.col(1) = m_values;
	return retVal;
}

unsigned int ShockProcess::numStates()
{
	return m_values.size();
}

const VectorXd* ShockProcess::states()
{
	return &m_values;
}

void ShockProcess::printStates()
{
	for (int i = 0; i < m_values.size(); i++) {
		std::cout << "Shimer state " << i << "=" << m_values[i] << std::endl;
	}
	return;
}
