// EZ_No_Assets.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "OLGModel.h"
#include "CobbDouglasMatching.h"

int main()
{

	//create matching function targetting f=X and eta=Y
	CobbDouglasMatching myF(0.45, 0.28);

	//create model with N generations, F matching function
	OLGModel model(3, 1.0, myF);

	//solve for all wages
	model.solveWages();

	//solve for elasticity
	return 0;
}