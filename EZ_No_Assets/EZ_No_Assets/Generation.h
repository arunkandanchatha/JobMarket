#pragma once

#include <vector>
#include "macros.h"
#include "FunctorBase.h"
class OLGModel;

typedef double (OLGModel::*OLGModelMemFn)(double x, double Up1, double Ep1, double Wp1);

class Generation : public FunctorBase
{
public:
	Generation(OLGModel &p_model, OLGModelMemFn p_fn, double p_Up1, double p_Ep1, double p_Wp1);
	~Generation();

	double operator()(double x);

private:
	const double Ep1,Up1,Wp1;
	OLGModel *m_model;
	OLGModelMemFn m_fn;

};

