#pragma once

#include "stdafx.h"
#include "FunctorBase.h"
#include "ShockProcess.h"
#include <vector>
class OLGModel;

typedef double (OLGModel::*OLGModelMemFn)(int generation, int state, int habit, int tenure, int wageLastPeriod, double x);

class Generation : public FunctorBase
{
public:
	Generation(OLGModel &p_model, OLGModelMemFn p_fn, int p_generation, int p_state, int p_habit, int p_tenure, int p_wageLastPeriod);
	~Generation();

	double operator()(double x);
	double operator()(const std::vector<double> &x, std::vector<double> &grad);
	static double wrap(const std::vector<double> &x, std::vector<double> &grad, void *data);

private:
	VectorXd *Ep1,*Up1,*Wp1;
	pdfMatrix *m_nextPDF;
	int m_state;
	int m_whichGen;
	int tenure;
	int m_habit;
	int m_wageLastPeriod;
	OLGModel *m_model;
	OLGModelMemFn m_fn;

};