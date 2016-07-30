#pragma once

#include "stdafx.h"
#include "Generation.h"
#include "MatchingFunction.h"
#include "MySolver.h"

#define D_b (0.4)
#define D_BETA (1.0 / 1.004)
#define D_RHO (-1.0)
#define D_GAMMA (0.72)
#define D_S (0.034)
#define D_ETA (1-D_GAMMA)
#define D_l (log(1-D_ETA)/log(0.45)) 

class OLGModel
{
	friend class Generation;
public:
	OLGModel(unsigned int generations, double y, MatchingFunction &f);
	OLGModel(OLGModel &orig);
	~OLGModel();
	void solveWages();
	double elasticity();
	void printWages();

protected:
	double nonLinearWageEquation(double x, double Up1, double Ep1, double Wp1);

private:
	void solveWages(unsigned int generation);
	double calcU(double Up1, double Ep1);
	double calcE(double delta, double Up1, double Ep1);
	double calcW(double delta, double Wp1);
	double partialE_partialDel(double x, double Up1, double Ep1);
	double expectedW();

	std::vector<double> E_vals;
	std::vector<double> U_vals;
	std::vector<double> W_vals;
	std::vector<double> wages;

	const double m_y;
	const unsigned int m_gens;
	MatchingFunction *m_f;
};