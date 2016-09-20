#include "ShimerProcess.h"



ShimerProcess::ShimerProcess(int states, double sigma, double gamma)
	: ShockProcess(2*states+1), m_sigma(sigma), m_gamma(gamma)
{
	if (states < 1) {
		std::cout << "must have at least two states for shimer process" << std::endl;
		exit(-1);
	}
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

	for (int i = 0; i < numStates(); i++) {
		m_values(i) = ((int)i-states)*delta;
	}

#if 0
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
		m_conditionalProbs(i, i + 1) = 0.5*(1 - m_values(i) / (states*delta));
		m_conditionalProbs(i, i - 1) = 0.5*(1 + m_values(i) / (states*delta));
	}
#else
	//let's try something different
	//one row: probability of transitioning from index i to index i+1,
	//conditional of getting a shock. Obviously, the probability of 
	//transitioning to i-1 is 1-p(i+1)
	VectorXd posProb(numStates());
	for (int i = 0; i < numStates(); i++) {
//		posProb(i) = 0.5*(1 - m_values(i) / (states*delta));
		posProb(i) = 0.5*(1 - m_values(i) / (1000*delta));
	}
	posProb(0) = 1;
	posProb(numStates() - 1) = 0;
	//we have a poisson process with (on average) one shock per quarter, or
	//alternatively, 1/3 of a shock per month. We can find the probability
	//of having x shocks using the formula:
	//         P(X=x)=exp(-1/3)(1/3)^(x)/(x!)
	if ((MAX_SHOCKS_PER_MONTH + 1) != 11) {
		std::cout << "ERROR! ShimerProcess.cpp:calculateConditionalProbabilities() - expect only 10 shocks,"
			<< " as that is as many entries in factorial array as we have entered." << std::endl;
		exit(-1);
	}
	std::vector<double> shockProb(MAX_SHOCKS_PER_MONTH+1);
	int fact[MAX_SHOCKS_PER_MONTH+1] = { 1,1,2,6,24,120,720,5040,40320,362880,3628800};
	double cumeVal = 0;
	for (int i = 0; i < MAX_SHOCKS_PER_MONTH; i++) {
		shockProb[i] = exp(-4.0 / 3)*pow(4.0 / 3, i) / fact[i];
		cumeVal += shockProb[i];
	}
	shockProb[MAX_SHOCKS_PER_MONTH] = 1 - cumeVal;


	MatrixXd probMatrix(numStates(), numStates());
	//the starting probability for not moving is the probability of no shocks;
	for (int i = 0; i < numStates(); i++) {
		for (int j = 0; j < numStates(); j++) {
			probMatrix(i, j) = 0;
		}
		probMatrix(i, i) = shockProb[0];
	}
	//now, given x shocks, we are at most x steps away from our starting state.
	//but maybe we got a mix of good and bad shocks. So, given x shocks,
	//where g are good, and b are bad (with g+b=x), 
	for (int i = 0; i < numStates(); i++) {
		VectorXd temp = probMatrix.row(i);
		for (int j = 1; j <= MAX_SHOCKS_PER_MONTH; j++) {
			setProbMatrix(i, j, shockProb[j], 1, posProb, temp);
		}
		probMatrix.row(i) = temp;
		double mySum = 0;
		for (int j = 0; j < numStates(); j++) {
			mySum += probMatrix(i, j);
		}
		//dealing with floating point, so let's enforce sum
		//of row to be 1
		for (int j = 0; j < numStates(); j++) {
			probMatrix(i, j) /= mySum;
		}
	}
	m_conditionalProbs = probMatrix;
#endif
	return;
}

void ShimerProcess::setProbMatrix(int currentState, int remainingShocks, double overallMult, double cumeProb, VectorXd& transProb, VectorXd& values) {
	if (remainingShocks == 0) {
		values[currentState] += (overallMult*cumeProb);
		return;
	}
	if (currentState < (numStates() - 1)) {
		setProbMatrix(currentState + 1, remainingShocks - 1, overallMult, cumeProb * transProb[currentState], transProb, values);
	}
	if (currentState > 0) {
		setProbMatrix(currentState - 1, remainingShocks - 1, overallMult, cumeProb * (1-transProb[currentState]), transProb, values);
	}
	return;
}

void ShimerProcess::printMatrix() {
	using namespace std;

	std::cout << "state,";
	int maxVal = (int)m_conditionalProbs.rows();
	for (int i = 0; i < maxVal; i++) {
		std::cout << i;
		if (i != m_conditionalProbs.size() - 1) {
			std::cout << ",";
		}
	}
	std::cout << std::endl;

	for (int i = 0; i < maxVal; i++) {
		std::cout << i << ",";
		for (int j = 0; j < maxVal-1; j++) {
			std::cout << m_conditionalProbs(i, j) << ",";
		}
		std::cout << m_conditionalProbs(i, maxVal-1) << std::endl;
	}

	ostringstream os;
	ofstream out_stream;
	os << "transition.csv";
	//out_stream.precision(15);
	out_stream << std::scientific;
	out_stream.open(os.str(), std::ofstream::out/* | std::ofstream::app*/);

	for (int i = 0; i < maxVal; i++) {
		out_stream << i;
		if (i != m_conditionalProbs.size() - 1) {
			out_stream << ",";
		}
	}
	out_stream << std::endl;

	for (int i = 0; i < maxVal; i++) {
		out_stream << i << ",";
		for (int j = 0; j < maxVal - 1; j++) {
			out_stream << m_conditionalProbs(i, j) << ",";
		}
		out_stream << m_conditionalProbs(i, maxVal - 1) << std::endl;
	}

	out_stream.close();

}
