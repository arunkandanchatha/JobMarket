#pragma once

#include "stdafx.h"
#include "FunctorBase.h"
#include "ShockProcess.h"
#include <vector>
class OLGModel;

//typedef double (OLGModel::*OLGModelMemFn)(int state, int whichGen, double x, VectorXd &Up1, VectorXd &Ep1, VectorXd &Wp1);
typedef double (OLGModel::*OLGModelMemFn)(int state, double x, VectorXd &Up1, VectorXd &Ep1, VectorXd &Wp1);

class Generation : public FunctorBase
{
public:
//	Generation(OLGModel &p_model, OLGModelMemFn p_fn, int state, int whichGen, VectorXd &p_Up1, VectorXd &p_Ep1, VectorXd &p_Wp1);
	Generation(OLGModel &p_model, OLGModelMemFn p_fn, int state, VectorXd &p_Up1, VectorXd &p_Ep1, VectorXd &p_Wp1);
	~Generation();

	double operator()(double x);
	double operator()(const std::vector<double> &x, std::vector<double> &grad);
	static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data);

private:
	VectorXd *Ep1,*Up1,*Wp1;
	pdfMatrix *m_nextPDF;
	int m_state;
//	int m_whichGen;
	OLGModel *m_model;
	OLGModelMemFn m_fn;

};