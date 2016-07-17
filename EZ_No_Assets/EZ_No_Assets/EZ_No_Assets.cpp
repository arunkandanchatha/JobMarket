// EZ_No_Assets.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "OLGModel.h"
#include "MySolver.h"

int main()
{
	OLGModel myModel;
	MySolver delta2(myModel, &OLGModel::solveDel2, 1E-10);
	double del2 = delta2.solve(0, 1);
	std::cout << "y2(" << del2 << ")=" << delta2.value(del2) << std::endl;
	myModel.setDel2(del2);

	MySolver delta1(myModel, &OLGModel::solveDel1, 1E-15);
	double del1 = delta1.solve(del2, 1);
	myModel.setDel1(del1);
	std::cout << "y1(" << del1 << ")=" << delta1.value(del1) << std::endl;
	return 0;
}