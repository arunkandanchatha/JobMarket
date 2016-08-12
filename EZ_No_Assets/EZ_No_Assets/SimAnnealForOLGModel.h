#pragma once
#include "SimulatedAnnealingInterface.h"
#include <vector>
#include "OLGModel.h"
#include "macros.h"

class SimAnnealForOLGModel :
	public SimulatedAnnealingInterface<std::vector<double>,double>
{
public:

	void printStatus(const std::vector<double>& solution, double temp);

	SimAnnealForOLGModel(const std::vector<double>& startSolution, const double& target,
		double starttemp, double precision, OLGModel& f, double p_maxstep, int max_mult);

	~SimAnnealForOLGModel();

	std::vector<double>* giveRandomNeighbour(const std::vector<double>& lastSolution) const;

	double calcDistanceToTarget(const std::vector<double>& solution) const;

	bool shouldStopHook(const std::vector<double>& solution);


	double getEnergy(const std::vector<double>& lastSolution) const {
		std::vector<double> temp;
		return (*func)(lastSolution, temp);
	};

private:
	double unifRand() const {
		return rand() / double(RAND_MAX);
	};

	double unifRand(double a, double b) const {
		double x;
		x = (b - a) * unifRand() + a;
		return x;
	};

	OLGModel *func;
	const double maxstep;
	long int numCalls = 0;
	int maxMult;

};

