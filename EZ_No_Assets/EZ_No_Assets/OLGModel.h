#pragma once

#include "stdafx.h"
#include "Generation.h"
#include "MatchingFunction.h"
#include "MySolver.h"
#include "ShockProcess.h"

#define D_b (0.4)
#define D_BETA (1.0 / 1.004)
#define D_RHO (-1.0)

class OLGModel
{
	friend class Generation;
public:
	OLGModel(unsigned int generations, double y, double s, MatchingFunction &f, ShockProcess &sp);
	~OLGModel();
	void solveWages();
	double elasticityWRTymb();
	double elasticityWRTs();
	void printWages();
//	VectorXd operator()(VectorXd& p_thetaGuess);

protected:
	double nonLinearWageEquation(int state, double x, VectorXd &Up1, VectorXd &Ep1, VectorXd &Wp1, pdfMatrix& nextPDF);

private:
	double calcU(int state, VectorXd &Up1, VectorXd &Ep1, pdfMatrix& nextPDF);
	double calcE(int state, double delta, VectorXd &Up1, VectorXd &Ep1, pdfMatrix& nextPDF);
	double calcW(int state, double delta, VectorXd &Wp1, pdfMatrix& nextPDF);
	double partialE_partialDel(int state, double x, VectorXd &Up1, VectorXd &Ep1, pdfMatrix& nextPDF);
	double expectedW(int state);

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
	const double m_gamma;
};