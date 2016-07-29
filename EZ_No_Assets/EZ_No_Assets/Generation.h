#pragma once

#include <vector>

class Generation
{
public:
	Generation(std::vector<double> &p_E, std::vector<double> &p_U, std::vector<double> &p_dE);
	~Generation();

	double operator()(double x);

private:
	std::vector<double> E;
	std::vector<double> U;
	std::vector<double> dE;

};

