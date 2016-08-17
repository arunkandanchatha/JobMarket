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
#include "adept_source.h"
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
			<< "./executable <s|c|e|a> gens states fTarget [thetaGuess.file]" << std::endl
			<< "where s-solve (simulated annealing), c-solve (COBYLA) and e-elasticity" << std::endl;
		return 1;
	}

	std::cout << "Number of Generations = " << numGens << std::endl;
	std::cout << "Number of States (always odd) = " << numStates << std::endl;

	if(which == 'e')
	{
		//solve for all wages
		//create matching function targetting f=X and eta=Y
		CobbDouglasMatching myF(fTarget, 1-D_ETA);

		//create shock process
		NoShocksProcess p;

		//create model with N generations, F matching function
		OLGModel model(numGens, 1.0, 0.034, myF, p, 1 - D_ETA, false);
		model.solveWages();
		model.printWages();
		//solve for elasticity
		std::cout << myF.getParameter() << "," << model.elasticityWRTymb() << std::endl;
		std::cout << myF.getParameter() << "," << model.elasticityWRTs() << std::endl;

	}else{
		bool autoDiff = false;
		if (which == 'a') {
			autoDiff = true;
		}
		std::vector<double> x(numStates);
		if (readFile) {
			initialize(x, argv[5]);
		}
		else {
			x.resize(3);
			double midX = 0.50;
			for (int i = 0; i < x.size(); i++) {
				x[i] = midX + 0.01*(i - 1);
			}
		}
		for (int solveIndex = x.size(); solveIndex <= numStates; solveIndex += 2) {
			//create matching function targetting f=X and eta=Y
			deHaanMatching myF(fTarget);

			//create shock process
			ShimerProcess p((solveIndex - 1) / 2, 0.0165/sqrt(3), (4.0 / 3)/((solveIndex - 1) / 2));

			//create model with N generations, F matching function
			OLGModel model(numGens, 1.0, D_S, myF, p, 1 - D_ETA,autoDiff);

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
				case 'a':
					algoChoice = nlopt::LD_MMA;
					break;
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
				opt.set_maxeval(solveIndex * 250);
				opt.set_lower_bounds(readFile ? (x[0] / 2.0) : 0.00001);
				opt.set_upper_bounds(readFile ? (2 * x[solveIndex - 1]) : 5);
				opt.set_min_objective(OLGModel::wrap, &model);
				opt.set_population(5 * solveIndex);

				std::vector<int> data(solveIndex);
				if (which == 'c') {
					for (int i = 0; i < solveIndex - 1; i++) {
						data[i] = i + 1;
						opt.add_inequality_constraint(myConstraint, &data[i], 0);
					}
				}
				opt.set_stopval(5e-3);

				double minf = 200;
				nlopt::result result;
				try {
					result = opt.optimize(x, minf);
				}
				catch (const nlopt::roundoff_limited& e) {
					std::cout << e.what() << std::endl;
				}
				catch (const nlopt::forced_stop& e) {
					std::cout << e.what() << std::endl;
					exit(-1);
				}
				catch (const std::runtime_error& e) {
					std::cout << e.what() << std::endl;
					exit(-1);
				}
				catch (const std::exception &e) {
					std::cout << e.what() << std::endl;
					exit(-1);
				}
				std::cout << "found minimum value " << minf << " at " << std::endl;
				for (int i = 0; i < solveIndex; i++) {
					std::cout << "theta(" << i << ")=" << x[i] << std::endl;
				}
				std::vector<double> myTempVector;
				OLGModel::printStatus(x, -1, minf);
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
			//model.printWages();
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
