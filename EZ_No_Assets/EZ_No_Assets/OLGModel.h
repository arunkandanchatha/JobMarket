#pragma once

#include "stdafx.h"
#include "Generation.h"
#include "MatchingFunction.h"
#include "MySolver.h"
#include "ShockProcess.h"
#include <vector>

#define D_b (0.4)
#define D_BETA (1.0 / 1.004)
#define D_RHO (-1.0)
#define D_C (1)
#define D_S (0.034)
#define D_ETA (0.28)

class OLGModel
{
	friend class Generation;
public:
	OLGModel(unsigned int generations, double y, double s, MatchingFunction &f, ShockProcess &sp, double bargaining);
	~OLGModel();
	void solveWages();
	double elasticityWRTymb();
	double elasticityWRTs();
	void printWages();
	static void printStatus(const std::vector<double>& solution, int numCalls, double distance);
	double operator()(const std::vector<double> &x, std::vector<double> &grad);
	static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data);

protected:
	double nonLinearWageEquation(int state, double x, VectorXd &Up1, VectorXd &Ep1, VectorXd &Wp1);

private:
	double calcU(int state, VectorXd &Up1, VectorXd &Ep1);
	double calcE(int state, double delta, VectorXd &Up1, VectorXd &Ep1);
	double calcW(int state, double delta, VectorXd &Wp1);
	double partialE_partialDel(int state, double x, VectorXd &Up1, VectorXd &Ep1);
	double expectedW(int state, bool forceNoShocks=false);

	int lastSolveGen;

	MatrixXd E_vals;
	MatrixXd U_vals;
	MatrixXd W_vals;
	MatrixXd wages;
	VectorXd m_thetas;

	VectorXd m_Y;
	const double m_Es;
	const unsigned int m_gens;
	MatchingFunction *m_f;
	ShockProcess *m_sp;
	const double m_bargaining;
};