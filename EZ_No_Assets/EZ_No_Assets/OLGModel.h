#pragma once

#include "stdafx.h"

#define D_b (0.4)
#define D_BETA (1.0 / 1.004)
#define D_RHO (1.0)
#define D_Y (1.0)
#define D_GAMMA (0.72)
#define D_F (0.45)
#define D_S (0.034)

class OLGModel
{
public:
	OLGModel();
	~OLGModel();
	double solveDel2(double x);
	double solveDel1(double x);
	void setDel1(double x);
	void setDel2(double x);

private:
	double E2(double x);
	double E1(double x);
	double U2(double x);
	double U1(double x);
	double dE2_dDel(double x);
	double dE1_dDel(double x);

	double del1, del2;
};