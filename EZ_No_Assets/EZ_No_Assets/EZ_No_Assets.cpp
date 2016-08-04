// EZ_No_Assets.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "OLGModel.h"
#include "CobbDouglasMatching.h"
#include "NoShocksProcess.h"
#include "ShimerProcess.h"

int main(int argc, char *argv[])
{

	int numGens = 50;

	if (argc == 2) {
		numGens = atoi(argv[1]);
	}

	std::cout << "Number of Generations = " << numGens << std::endl;

//	for (int i = 1; i < 100; i++) {

	int i = 28;

		//create matching function targetting f=X and eta=Y
		CobbDouglasMatching myF(0.45, i*0.01);

		//create shock process
		//NoShocksProcess p;
		ShimerProcess p(10, 0.0165, 0.004);

		//create model with N generations, F matching function
		OLGModel model(numGens, 1.0, 0.034, myF, p);

		//solve for all wages
		model.solveWages();

		//solve for elasticity
		std::cout << myF.getEta() << "," << model.elasticityWRTymb() << std::endl;
		std::cout << myF.getEta() << "," << model.elasticityWRTs() << std::endl;
//	}
	return 0;
}