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
			<< "./executable <s|c|e|a|x> gens states fTarget [thetaGuess.file]" << std::endl
			<< "where " << std::endl
			<< "    s-solve (simulated annealing)" << std::endl
			<< "    c-solve (COBYLA)" << std::endl
			<< "    e-elasticity" << std::endl
			<< "    x-simulate only" << std::endl;
		return 1;
	}

	std::cout << "Number of Generations = " << numGens << std::endl;
	std::cout << "Number of States (always odd) = " << numStates << std::endl;

	std::vector<double> x(numStates);

	if (which == 'e')
	{
		for (int i = 1; i < 81; i++) {
			//solve for all wages
			//create matching function targetting f=X and eta=Y
			CobbDouglasMatching myF(fTarget, 0.01*i/*D_ETA*/);

			//create shock process
			NoShocksProcess p;

			//create model with N generations, F matching function
			OLGModel model(numGens, 1.0, D_S, myF, p, 1 - D_ETA, false);
			model.solveWages();
			//model.printWages();
			//solve for elasticity
			std::cout << myF.getParameter() << "," << model.elasticityWRTymb() << std::endl;
			//std::cout << myF.getParameter() << "," << model.elasticityWRTs() << std::endl;
		}
		return 1;
	}
	bool autoDiff = false;
	if (which == 'a') {
		autoDiff = true;
	}
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
	if (which != 'x') {
		for (int solveIndex = x.size(); solveIndex <= numStates; solveIndex += 2) {
			//create matching function targetting f=X and eta=Y
			deHaanMatching myF(fTarget);

			//create shock process
			ShimerProcess p((solveIndex - 1) / 2, 0.0165, 4.0 / ((solveIndex - 1) / 2));

			//create model with N generations, F matching function
			OLGModel model(numGens, 1.0, D_S, myF, p, 1 - D_ETA, autoDiff);

			if (which == 's') {
				const double targ = 0;
				SimAnnealForOLGModel worldSolve = SimAnnealForOLGModel(x, targ, 1, 1.0E-8, model, (x[solveIndex - 1] - x[0]) / (solveIndex - 1) / 5, solveIndex * 1000);
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
				opt.set_stopval(1e-4);

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
			newX[0] = MAX(x[0] - 0.01, 0);
			for (int myIndex = 0; myIndex < solveIndex; myIndex++) {
				newX[myIndex + 1] = x[myIndex];
			}
			newX[solveIndex + 1] = newX[solveIndex] + 0.01;

			x.clear();
			x.resize(solveIndex + 2);
			x = newX;
			//model.printWages();
		}

		return 1;
	}
	//now we simulate
	srand(1234);

	//create shock process
	ShimerProcess p((x.size() - 1) / 2, 0.0165, 4.0 / ((x.size() - 1) / 2));
	p.printMatrix();
	exit(-1);

	//create matching
	deHaanMatching myF(fTarget);

	//simulate 100 times
	std::cout << "simulation      u     theta (state)" << std::endl;
	for (int i = 0; i < 10000; i++) {
		//choose random starting unemployment between 0 and 1
		double u = rand() / double(RAND_MAX);
		//choose random starting "shock" state
		double randNum = rand() / double(RAND_MAX);
		int currState = floor((x.size()-1) * randNum);

		double currTheta = x[currState];

		//std::cout << i << "    " << u << "     " << currTheta << " " << currState << std::endl;
		pdfMatrix nextPDF = p.nextPeriodPDF(currState);
		//each time, simulate 5000 periods
		for (int j = 0; j < 5000; j++) {
			//update unemployment rate
			double loseJobs = D_S*(1 - u);
			double gainJobs = myF.calculatedF(currTheta)*u;
			double changeU = loseJobs - gainJobs;
			u += changeU;
			if (j > 4400) {
				std::cout << i << "    " << u << "     " << currTheta << "  " << currState << "    " << std::endl;
			}

			double randNumForShocks = rand() / double(RAND_MAX+1);
			double cumeTotal = 0;
			bool found = false;
			for (int k = MAX(currState - MAX_SHOCKS_PER_MONTH, 0); 
				k < MIN(x.size(),currState + MAX_SHOCKS_PER_MONTH + 1); k++) {
				cumeTotal += nextPDF(k, 0);
				if (randNumForShocks <= cumeTotal) {
					found = true;
					currState = k;
					currTheta = x[currState];
					nextPDF = p.nextPeriodPDF(currState);
					break;
				}
			}
			if (!found) {
				std::cout << "ERROR! EZ_No_Assets.cpp-main(): Couldn't find new state for randNum=" << randNumForShocks << std::endl;
				std::cout << "     CumeTotal = " << cumeTotal << std::endl;
				std::cout << "     State = " << currState << std::endl;
				exit(-1);
			}
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
