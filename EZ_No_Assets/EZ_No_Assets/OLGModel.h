#pragma once

#include "stdafx.h"
#include "Generation.h"
#include "MatchingFunction.h"
#include "ShockProcess.h"
#include "OLGSolveAutoDiff.h"
#include "utilities.h"
//#include "MySolver.h"

#include <vector>

class OLGModel
{
	friend class Generation;
public:
	OLGModel(int generations, double y, double s, MatchingFunction &f, ShockProcess &sp, double bargaining, bool autoDiff);
	~OLGModel();
	void solveWages();
	double elasticityWRTymb();
	double elasticityWRTs();
	std::vector<double> wageElasticityWRTymb();
	void printWages();
	static void printStatus(const std::vector<double>& solution, int numCalls, double distance);
	double operator()(const std::vector<double> &x, std::vector<double> &grad);
	static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data);
	OLGSolveAutoDiff getSolver();
//	void printWageElasticity();

protected:
	double nonLinearWageEquation(int generation, int state, int habit, int tenure, int wageLastPeriod, double x);

private:
	VectorXd OLGModel::getSteadyStateDistrib(double jobFind);
	double calcU(int generation, int state, int habit, int tenure, int wageLastPeriod);
	double calcE(int generation, int state, int habit, int tenure, int wageLastPeriod, double x);
	double calcW(int generation, int state, int habit, int tenure, int wageLastPeriod, double x);
	double partialE_partialDel(int generation, int state, int habit, int tenure, int wageLastPeriod, double x);
	double expectedW0(int state, bool forceNoShocks=false);
	double adjustmentCost(const double original, const double updateVal);

	VectorXd m_thetas;

	std::vector<std::vector<std::vector<MatrixXd>>> E_vals;
	std::vector<std::vector<std::vector<MatrixXd>>> U_vals;
	std::vector<std::vector<std::vector<MatrixXd>>> W_vals;
	std::vector<std::vector<std::vector<MatrixXd>>> wages;
	std::vector<MatrixXd> m_Y;
	MatrixXd habits;
	MatrixXd habitProb;
	std::vector<std::vector<double>> oldWages;
	const double m_Es;
	const int m_gens;
	MatchingFunction *m_f;
	ShockProcess *m_sp;
	const double m_bargaining;
	const bool m_autodiff;
};