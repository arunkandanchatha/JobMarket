#include "Generation.h"

#if 1
Generation::Generation(OLGModel &p_model, OLGModelMemFn p_fn, int state, int whichGen, int p_tenure/*, VectorXd &p_Up1, 
	VectorXd &p_Ep1, VectorXd &p_Wp1*/)
	: m_model(&p_model), m_fn(p_fn), /*Ep1(&p_Ep1), Up1(&p_Up1), Wp1(&p_Wp1),*/ m_state(state), m_whichGen(whichGen),
	tenure(p_tenure)
{
}
#else
Generation::Generation(OLGModel &p_model, OLGModelMemFn p_fn, int state, VectorXd &p_Up1, VectorXd &p_Ep1, VectorXd &p_Wp1)
	: m_model(&p_model), m_fn(p_fn), Ep1(&p_Ep1), Up1(&p_Up1), Wp1(&p_Wp1), m_state(state)
{
}
#endif

Generation::~Generation()
{
}

double Generation::operator()(double x) {
	return CALL_MEMBER_FN(*m_model, m_fn)(m_state, m_whichGen, tenure, x/*, *Up1, *Ep1, *Wp1*/);
//	return CALL_MEMBER_FN(*m_model, m_fn)(m_state, x, *Up1, *Ep1, *Wp1);
}

double Generation::operator()(const std::vector<double> &x, std::vector<double> &grad) {
	return CALL_MEMBER_FN(*m_model, m_fn)(m_state, m_whichGen, tenure, x[0]/*, *Up1, *Ep1, *Wp1*/);
//	return CALL_MEMBER_FN(*m_model, m_fn)(m_state, x[0], *Up1, *Ep1, *Wp1);
}

double Generation::wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
	double retVal = (*reinterpret_cast<Generation*>(data))(x, grad);
	return retVal;
}