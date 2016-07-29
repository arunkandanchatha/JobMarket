// EZ_No_Assets.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "OLGModel.h"
#include "MySolver.h"

int main()
{

	double theta = 0.05774;
	double theta2 = theta + 0.0001;

	double y = 1.0;
	OLGModel myModel(y,theta);

	MySolver delta2(myModel, &OLGModel::solveDel2, 1E-15);
	double del2 = delta2.solve(0, 1);
	std::cout << "y2(" << del2 << ")=" << delta2.value(del2) << std::endl;
	myModel.setDel2(del2);

	MySolver delta1(myModel, &OLGModel::solveDel1, 1E-15);
	double del1 = delta1.solve(del2, 1);
	std::cout << "y1(" << del1 << ")=" << delta1.value(del1) << std::endl;
	myModel.setDel1(del1);

	OLGModel myModelTheta(myModel);
	myModelTheta.setTheta(theta2);
	MySolver deltaTheta(myModelTheta, &OLGModel::solveDel2, 1E-15);
	double del2Theta = deltaTheta.solve(0, 1);
	std::cout << "y2(" << del2Theta << ")=" << deltaTheta.value(del2Theta) << std::endl;
	myModelTheta.setDel2(del2Theta);

	MySolver deltaTheta1(myModelTheta, &OLGModel::solveDel1, 1E-15);
	double del1Theta = deltaTheta1.solve(del2Theta, 1);
	std::cout << "y1(" << del1Theta << ")=" << deltaTheta1.value(del1Theta) << std::endl;
	myModelTheta.setDel1(del1Theta);

	double y2 = 1.0001;
	OLGModel myModel2(y2, theta);
	MySolver delta22(myModel2, &OLGModel::solveDel2, 1E-15);
	double del22 = delta22.solve(0, 1);
	std::cout << "y2(" << del22 << ")=" << delta22.value(del22) << std::endl;
	myModel2.setDel2(del22);

	MySolver delta12(myModel2, &OLGModel::solveDel1, 1E-15);
	double del12 = delta12.solve(del22, 1);
	std::cout << "y1(" << del12 << ")=" << delta12.value(del12) << std::endl;
	myModel2.setDel1(del12);

	std::cout << "elasticity(" << D_ETA << ")=" << myModel.elasticity(myModelTheta, myModel2) << std::endl;

	std::cout << "No curvature" << std::endl;
	std::cout << "del1=" << D_GAMMA*(y - D_b) << "," << D_GAMMA*(y2 - D_b) << std::endl;
	std::cout << "del2=" << D_GAMMA*(y - D_b)+D_BETA*D_GAMMA*(1-D_GAMMA)*myModel.D_F()*(y-D_b) << "," << D_GAMMA*(y2 - D_b) + D_BETA*D_GAMMA*(1 - D_GAMMA)*myModel.D_F()*(y2 - D_b) << std::endl;


	std::cout << (1 + 0.5*D_BETA*(1 - D_S - D_GAMMA*myModel.D_F())) / ((1 + 0.5*D_BETA*(1 - D_S - myModel.D_F()*D_GAMMA))*(1 - D_ETA) + 0.5*D_BETA*D_GAMMA*myModel.D_F()*D_ETA) << std::endl;
	return 0;
}