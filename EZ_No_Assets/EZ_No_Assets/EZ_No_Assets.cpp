// EZ_No_Assets.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "OLGModel.h"
#include "CobbDouglasMatching.h"
#include "NoShocksProcess.h"
#include "ShimerProcess.h"
#include "nlopt.hpp"

double myConstraint(const std::vector<double> &x, std::vector<double> &grad, void*data);

int main(int argc, char *argv[])
{

	int numGens;
	int numStates;

	switch (argc) {
	case 3:
		numGens = atoi(argv[1]);
		numStates = atoi(argv[2]);
		if (numStates % 2 == 0) {
			numStates++;
		}
		break;
	default:
		std::cout << "invoke as follows:" << std::endl
			<< "./executable gens states" << std::endl;
		return 1;
	}

	std::cout << "Number of Generations = " << numGens << std::endl;
	std::cout << "Number of States (always odd) = " << numStates << std::endl;

	//create matching function targetting f=X and eta=Y
	CobbDouglasMatching myF(0.45, 0.28);

	//create shock process
	//NoShocksProcess p;
	ShimerProcess p((numStates-1)/2, 0.0165, 0.004);

	//create model with N generations, F matching function
	OLGModel model(numGens, 1.0, 0.034, myF, p);

#if 0
	//solve for all wages
	model.solveWages();

	//solve for elasticity
	std::cout << myF.getEta() << "," << model.elasticityWRTymb() << std::endl;
	std::cout << myF.getEta() << "," << model.elasticityWRTs() << std::endl;
#elif 0
	std::vector<double> x(numStates);
	for (int i = 0; i < numStates; i++) {
		x[i] = 1;
	}
	std::vector<double> temp;
	std::cout << "Distance: " << model(x, temp) << std::endl;
	model.printWages();
	p.printStates();

#else
//	nlopt::opt opt(nlopt::LN_COBYLA, numStates);
//	nlopt::opt opt(nlopt::GN_ORIG_DIRECT, numStates);
	nlopt::opt opt(nlopt::LN_BOBYQA, numStates);
	opt.set_lower_bounds(0.01);
	opt.set_min_objective(OLGModel::wrap, &model);

	std::vector<int> data(numStates);
	for (int i = 1; i < numStates; i++) {
		data[i] = i;
//		opt.add_inequality_constraint(myConstraint, &data[i], 1e-2);
	}
	opt.set_xtol_rel(1e-4);

#if 1
	std::vector<double> x(numStates);
	for (int i = 0; i < numStates; i++) {
		x[i] = 1;
	}
#else
	std::vector<double> x(11);
	x[0] = 0.63;
	x[1] = 0.66;
	x[2] = 0.72;
	x[3] = 0.89;
	x[4] = 1.15;
	x[5] = 1.51;
	x[6] = 2.04;
	x[7] = 2.57;
	x[8] = 3.16;
	x[9] = 3.37;
	x[10] = 3.4;

#endif

	double minf;
	nlopt::result result = opt.optimize(x, minf);
	std::cout << "found minimum value " << minf << " at " << std::endl;
	for (int i = 0; i < numStates; i++) {
		std::cout << "theta(" << i << ")=" << x[i] << std::endl;
	}
	model.printWages();
#endif

	return 0;
}

double myConstraint(const std::vector<double> &x, std::vector<double> &grad, void*data)
{
	int topX = *(int *)(data);
	if (topX < 1) {
		std::cout << "for constraints, must have topX > 0. " << topX << " is not valid" << std::endl;
		exit(-1);
	}
	return x[topX - 1] - x[topX];
}