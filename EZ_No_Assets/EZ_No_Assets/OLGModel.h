#pragma once

#include "stdafx.h"

#define D_b (0.4)
#define D_BETA (1.0 / 1.004)
#define D_RHO (-1.0)
#define D_GAMMA (0.72)
#define D_S (0.034)
#define D_ETA (1-D_GAMMA)
#define D_l (log(1-D_ETA)/log(0.45)) 

class OLGModel
{
public:
	OLGModel(double y, double theta);
	OLGModel(OLGModel &orig);
	~OLGModel();
	double solveDel2(double x);
	double solveDel1(double x);
	void setDel1(double x);
	void setDel2(double x);
	double elasticity(OLGModel &thetaChange, OLGModel &yChange);
	void setTheta(double x);
	//double D_F() { return D_THETA / pow(1 + pow(D_THETA, D_l), 1 / D_l); };
	double D_F() { return pow(D_THETA,D_ETA); };

private:
	double E2(double x);
	double E1(double x);
	double U1(double x);
	double U2(double x);
	double partialE2_partialDel2(double x);
	double partialE1_partialDel1(double x);
#if 0
	double totalDel1_totalymb(double x);
	double totalDel2_totalymb(double x);
	double dsqE1_dDelsq(double x);
	double dsqE2_dDelsq(double x);
	double dsqE1_dDeldDel(double x);
	double dsqE1_dDeldYmB(double x);
	double dsqE1_dDeldTheta(double x);

	double dU1_dTheta(double x);
#endif
	double del1, del2;
	double D_Y;
	double D_THETA;
};