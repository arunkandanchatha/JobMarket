#pragma once
#include "ShockProcess.h"
class NoShocksProcess :
	public ShockProcess
{
public:
	NoShocksProcess();
	~NoShocksProcess();

private:
	void calculateConditionalProbabilities();
};

