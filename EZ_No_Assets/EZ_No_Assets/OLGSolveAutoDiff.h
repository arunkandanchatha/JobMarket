#pragma once
#include "adept.h"
using adept::adouble;

#include <vector>
#include "OLGModel.h"

class OLGSolveAutoDiff
{
public:
	OLGSolveAutoDiff(int gens, std::vector<double> ys);
	~OLGSolveAutoDiff();

	double solveProblem(std::vector<double> input);

private:
	const int m_gens;
	std::vector<double> m_Y;
};

