#pragma once

#include "stdafx.h"
#include "Generation.h"
#include "MatchingFunction.h"
#include "ShockProcess.h"
#include "OLGSolveAutoDiff.h"
//#include "MySolver.h"

#include <vector>

class OLGModel
{
	friend class Generation;
public:
	OLGModel(unsigned int generations, double y, double s, MatchingFunction &f, ShockProcess &sp, double bargaining, bool autoDiff);
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

protected:
	double nonLinearWageEquation(int state, int whichGen, int tenure, double x/*, VectorXd &Up1, VectorXd &Ep1, VectorXd &Wp1*/);

private:
	double calcU(int state, int whichGen/*, VectorXd &Up1, VectorXd &Ep1*/);
	double calcE(int state, int whichGen, int tenure, double delta/*, VectorXd &Up1, VectorXd &Ep1*/);
	double calcW(int state, int whichGen, int tenure, double delta/*, VectorXd &Wp1*/);
	double partialE_partialDel(int state, int whichGen, int tenure, double x/*, VectorXd &Up1, VectorXd &Ep1*/);
	double expectedW0(int state, bool forceNoShocks=false);

	VectorXd m_thetas;

	std::vector<MatrixXd> E_vals;
	MatrixXd U_vals;
	std::vector<MatrixXd> W_vals;
	std::vector<MatrixXd> wages;
	std::vector<MatrixXd> m_Y;
	const double m_Es;
	const unsigned int m_gens;
	MatchingFunction *m_f;
	ShockProcess *m_sp;
	const double m_bargaining;
	const bool m_autodiff;
};