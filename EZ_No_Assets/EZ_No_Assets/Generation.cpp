#include "Generation.h"

Generation::Generation(OLGModel &p_model, OLGModelMemFn p_fn, int state, VectorXd &p_Up1, VectorXd &p_Ep1, VectorXd &p_Wp1)
	: m_model(&p_model), m_fn(p_fn), Ep1(&p_Ep1), Up1(&p_Up1), Wp1(&p_Wp1), m_state(state)
{
}


Generation::~Generation()
{
}

double Generation::operator()(double x) {
	return CALL_MEMBER_FN(*m_model, m_fn)(m_state, x, *Up1, *Ep1, *Wp1);
}

double Generation::operator()(const std::vector<double> &x, std::vector<double> &grad) {
	return CALL_MEMBER_FN(*m_model, m_fn)(m_state, x[0], *Up1, *Ep1, *Wp1);
}

double Generation::wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
	double retVal = (*reinterpret_cast<Generation*>(data))(x, grad);
	return retVal;
}