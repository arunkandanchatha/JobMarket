#include "Generation.h"

#if 1
Generation::Generation(OLGModel &p_model, OLGModelMemFn p_fn, int p_generation, int p_state, int p_habit, int p_tenure, int p_wageLastPeriod)
	: m_model(&p_model), m_fn(p_fn), m_state(p_state), m_whichGen(p_generation), tenure(p_tenure), m_habit(p_habit), m_wageLastPeriod(p_wageLastPeriod)
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
	return CALL_MEMBER_FN(*m_model, m_fn)(m_whichGen,m_state,m_habit,tenure,m_wageLastPeriod,x);
}

double Generation::operator()(const std::vector<double> &x, std::vector<double> &grad) {
	return CALL_MEMBER_FN(*m_model, m_fn)(m_whichGen, m_state, m_habit, tenure, m_wageLastPeriod, x[0]);
}

double Generation::wrap(const std::vector<double> &x, std::vector<double> &grad, void *data) {
	double retVal = (*reinterpret_cast<Generation*>(data))(x, grad);
	return retVal;
}