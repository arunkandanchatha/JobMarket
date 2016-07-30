// EZ_No_Assets.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "OLGModel.h"
#include "CobbDouglasMatching.h"

int main()
{

	//create matching function targetting f=X and eta=Y
	CobbDouglasMatching myF(0.45, D_ETA);

	//create model with N generations, F matching function
	OLGModel model(300, 1.0, myF);

	//solve for all wages
	model.solveWages();
	//model.printWages();

	//solve for elasticity
	std::cout << model.elasticity();
	return 0;
}