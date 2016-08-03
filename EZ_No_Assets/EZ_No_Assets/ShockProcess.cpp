#include "ShockProcess.h"



ShockProcess::ShockProcess(int states) : m_states(states), m_conditionalProbs(states, states)
{
}


ShockProcess::~ShockProcess()
{
}
