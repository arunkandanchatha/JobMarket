// EZ_No_Assets.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "OLGModel.h"
#include "CobbDouglasMatching.h"
#include "deHaanMatching.h"
#include "NoShocksProcess.h"
#include "ShimerProcess.h"
#include "nlopt.hpp"
#include "SimAnnealForOLGModel.h"
//#include "vld.h"

double myConstraint(const std::vector<double> &x, std::vector<double> &grad, void*data);
void initialize(std::vector<double> &x, char *filename);

int main(int argc, char *argv[])
{

	int numGens;
	int numStates;
	char which;
	double fTarget;
	bool readFile = false;
	switch (argc) {
	case 6:
		readFile = true;
	case 5:
		which = argv[1][0];
		numGens = atoi(argv[2]);
		numStates = atoi(argv[3]);
		fTarget = atof(argv[4]);
		if (numStates % 2 == 0) {
			numStates++;
		}
		break;
	default:
		std::cout << "invoke as follows:" << std::endl
			<< "./executable <s|c|b|e> gens states fTarget [thetaGuess.file]" << std::endl
			<< "where s-solve (simulated annealing)" << std::endl
			<< "      c - solve(COBYLA)"<<std::endl
			<< "      b - bobyqa" <<std::endl
			<< "      e - elasticity" << std::endl;
		return 1;
	}

	std::cout << "Number of Generations = " << numGens << std::endl;
	std::cout << "Number of States (always odd) = " << numStates << std::endl;

	if(which == 'e')
	{
		//solve for all wages
		//create matching function targetting f=X and eta=Y
		CobbDouglasMatching myF(fTarget, 0.28);

		//create shock process
		NoShocksProcess p;

		//create model with N generations, F matching function
		OLGModel model(numGens, 1.0, 0.034, myF, p);
		model.solveWages();

		//solve for elasticity
		std::cout << myF.getBargaining() << "," << model.elasticityWRTymb() << std::endl;
		std::cout << myF.getBargaining() << "," << model.elasticityWRTs() << std::endl;

	}else{
		std::vector<double> x(numStates);
		if (readFile) {
			initialize(x, argv[5]);
		}
		else {
			x.resize(3);
			double midX = 8.15;
			for (int i = 0; i < 3; i++) {
				x[i] = midX + 0.01*(i - 1);
			}
		}
		for (int solveIndex = x.size(); solveIndex <= numStates; solveIndex += 2) {
			//create matching function targetting f=X and eta=Y
			deHaanMatching myF(fTarget, .28);

			//create shock process
			ShimerProcess p((solveIndex - 1) / 2, 0.0165, 4.0 / ((solveIndex - 1) / 2));

			//create model with N generations, F matching function
			OLGModel model(numGens, 1.0, 0.1, myF, p);

			if (which == 's') {
				const double targ = 0;
				SimAnnealForOLGModel worldSolve = SimAnnealForOLGModel(x, targ, 1, 1.0E-8, model, (x[solveIndex - 1] - x[0]) / (solveIndex - 1) / 5, solveIndex*1000);
				std::vector<double> *soln = worldSolve.solve();
				std::vector<double> myTempVector;
				OLGModel::printStatus(*soln, -1, model(*soln, myTempVector));
				for (int i = 0; i < solveIndex; i++) {
					std::cout << "theta(" << i << ")=" << (*soln)[i] << std::endl;
				}
				delete soln;
			}
			else {
				nlopt::algorithm algoChoice;
				switch (which) {
				case 'c':
					algoChoice = nlopt::LN_COBYLA;
					break;
				case 'b':
					algoChoice = nlopt::LN_BOBYQA;
					break;
				default:
					std::cout << "Unknown algorithm type " << which << std::endl;
					exit(-1);
				}
				
				nlopt::opt opt(algoChoice, solveIndex);
				//nlopt::opt opt(nlopt::GN_ORIG_DIRECT, numStates);
				//nlopt::opt opt(nlopt::GN_ISRES, numStates);
				//nlopt::opt opt(nlopt::LN_BOBYQA, numStates);
				//nlopt::opt opt2(nlopt::LN_BOBYQA, numStates);
				//opt2.set_lower_bounds(0.01);
				//opt2.set_xtol_rel(1e-4);
				//opt2.set_maxeval(10000);
				//opt.set_local_optimizer(opt2);
				opt.set_maxeval(solveIndex * 500);
				opt.set_lower_bounds(x[0] / 2.0);
				opt.set_upper_bounds(2 * x[solveIndex - 1]);
				opt.set_min_objective(OLGModel::wrap, &model);
				opt.set_population(5 * solveIndex);

				std::vector<int> data(solveIndex);
				for (int i = 0; i < solveIndex - 1; i++) {
					data[i] = i + 1;
					opt.add_inequality_constraint(myConstraint, &data[i], 0);
				}
				//			data[0] = numStates - 2;
				//			opt.add_inequality_constraint(myConstraint2, &data[0], 1e-4);
				opt.set_xtol_rel(1e-8);

				double minf = 200;
				nlopt::result result = opt.optimize(x, minf);
				std::cout << "found minimum value " << minf << " at " << std::endl;
				for (int i = 0; i < solveIndex; i++) {
					std::cout << "theta(" << i << ")=" << x[i] << std::endl;
				}
			}
			std::vector<double> newX(solveIndex + 2);
			newX[0] = MAX(x[0] - 0.01,0);
			for (int myIndex = 0; myIndex < solveIndex; myIndex++) {
				newX[myIndex + 1] = x[myIndex];
			}
			newX[solveIndex + 1] = newX[solveIndex] + 0.01;

			x.clear();
			x.resize(solveIndex + 2);
			x = newX;
			//	model.printWages();
		}
	}

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

void initialize(std::vector<double> &x, char *filename) {

	using namespace std;

	ifstream file(filename);
	std::vector<string> value;
	string line;

	unsigned int counter = 0;

	while (getline(file, line)) {
		value.push_back(line);
	}

	if (value.size() > x.size()) {
		cout << "ERROR! EZ_No_Assets.cpp-initialize(file) : file has unexpected number of lines. " << value.size() << " instead of " << x.size() << endl;
		exit(-1);
	}
	else {
		x.resize(value.size());
	}

	for (int i = 0; i < x.size(); i++) {
		stringstream ss(value[i].c_str());
		string substring;
		//value
		getline(ss, substring);
		x[i] = atof(substring.c_str());
	}

	return;
}
