#include "Generation.h"

Generation::Generation(OLGModel &p_model, OLGModelMemFn p_fn, double p_Up1, double p_Ep1, double p_Wp1)
	: m_model(&p_model), m_fn(p_fn), Ep1(p_Ep1), Up1(p_Up1), Wp1(p_Wp1)
{
}


Generation::~Generation()
{
}

double Generation::operator()(double x) {
	return CALL_MEMBER_FN(*m_model, m_fn)(x,Up1,Ep1,Wp1);
}