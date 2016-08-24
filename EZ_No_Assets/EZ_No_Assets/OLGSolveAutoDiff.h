#pragma once
#include "adept.h"
using adept::adouble;

#include <vector>
#include "macros.h"
#include "OLGModel.h"

class OLGSolveAutoDiff
{
public:
	OLGSolveAutoDiff(int gens, std::vector<double> ys, MatrixXd conditionalProbs, double parameter,
		/*std::vector<double> bargaining,*/ double s, MatrixXd &wages);
	~OLGSolveAutoDiff();

	double solveProblem(std::vector<double>& input, std::vector<double>& bargaining);
	double solveProblem(std::vector<double>& input, std::vector<double>&grad, std::vector<double>& bargaining);

private:
	const int m_gens;
	std::vector<double> m_Y;
	MatrixXd* m_wgs;
	MatrixXd m_conditionalProbs;

	adept::Stack stack_;

	adouble m_parameter;
	adouble m_s;
};

