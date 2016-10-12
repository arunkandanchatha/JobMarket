#pragma once

#include "stdafx.h"

class utilities
{
public:
	utilities();
	~utilities();

	static double interpolate(const std::vector<double> &x, const std::vector<double> &y, double a_prime);
};

