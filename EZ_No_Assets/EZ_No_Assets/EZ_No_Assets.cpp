// EZ_No_Assets.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "OLGModel.h"
#include "CobbDouglasMatching.h"

int main(int argc, char *argv[])
{

	int numGens = 50;

	if (argc == 2) {
		numGens = atoi(argv[1]);
	}

	std::cout << "Number of Generations = " << numGens << std::endl;

	for (int i = 280; i < 281; i++) {
		//create matching function targetting f=X and eta=Y
		CobbDouglasMatching myF(0.45, i*0.001);

		//create model with N generations, F matching function
		OLGModel model(numGens, 1.0, 0.034, myF);

		//solve for all wages
		model.solveWages();

		//solve for elasticity
		//std::cout << myF.getEta() << "," << model.elasticityWRTymb() << std::endl;
		std::cout << myF.getEta() << "," << model.elasticityWRTs() << std::endl;
	}
	return 0;
}