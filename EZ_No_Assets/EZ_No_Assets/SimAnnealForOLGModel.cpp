#include "SimAnnealForOLGModel.h"


SimAnnealForOLGModel::SimAnnealForOLGModel(const std::vector<double>& startSolution, const double& target,
	double starttemp, double precision, OLGModel& f, double p_maxstep, int max_mult) :
	SimulatedAnnealingInterface<std::vector<double>, double>(startSolution, target, starttemp, precision), maxstep(p_maxstep),
	maxMult(max_mult) {
	static long int seed = 1;
	srand(seed++);
	func = &f;
}

SimAnnealForOLGModel::~SimAnnealForOLGModel()
{
}

std::vector<double>* SimAnnealForOLGModel::giveRandomNeighbour(
	const std::vector<double>& lastPoint)  const {
	std::vector<double> *newPoint = new std::vector<double>(lastPoint.size());
	*newPoint = lastPoint;

	for (int i = 0; i < lastPoint.size(); i++) {
		double diffFromPrev = lastPoint[i] - lastPoint[MAX(i - 1, 0)];
		double diffFromNext = lastPoint[MIN(i + 1, lastPoint.size() - 1)]- lastPoint[i];
		double dist;
		if (i == 0) {
			//don't go above 1/2 way to the next point
			double modifiedMaxStep = MIN(diffFromNext / 2, maxstep);
			//don't go below 1/2 way to zero
			dist = (unifRand(-MIN(maxstep,lastPoint[i]/2), modifiedMaxStep));
		}
		else if (i == lastPoint.size() - 1) {
			//don't go below 1/2 way to previous point
			double modifiedMaxStep2 = MIN(diffFromPrev / 2, maxstep);
			dist = (unifRand(-modifiedMaxStep2, maxstep));
		}
		else {
			//don't go above 1/2 way to the next point
			double modifiedMaxStep = MIN(diffFromNext / 2, maxstep);
			//don't go below 1/2 way to previous point
			double modifiedMaxStep2 = MIN(diffFromPrev/2, maxstep);
			dist = (unifRand(-modifiedMaxStep2, modifiedMaxStep));
		}
		(*newPoint)[i] = lastPoint[i] + dist;

	}
	return newPoint;
}

double SimAnnealForOLGModel::calcDistanceToTarget(const std::vector<double>& angle) const {
	double distance = getEnergy(angle) - *TARGET;

	if (distance != distance) {
		std::cout << "Error! SimAnnealForOLGModel.calcDistanceToTarget: distance is NaN" << std::endl;
		for (int i = 0; i < angle.size(); i++) {
			std::cout << i << ":" << angle[i] << std::endl;
		}
		exit(-1);
	}

	if (distance > 0) {
		return distance;
	}
	else {
		return (-distance);
	}
}

bool SimAnnealForOLGModel::shouldStopHook(const std::vector<double>& solution) {
	numCalls++;
	return numCalls > maxMult;
}

void SimAnnealForOLGModel::printStatus(const std::vector<double>& solution, double temp) {
	using namespace std;

	std::cout << "SimAnnealForOLGModel.printStatus(): At temp=" << temp << ", evals="<< numCalls << " new best has distance=" << calcDistanceToTarget(*best)
		<< " where current solution is: " << std::endl;
	for (int i = 0; i < solution.size(); i++) {
		printf("%4i %10lf\r\n", i, solution[i]);
	}
	std::cout << std::flush;

	ostringstream os;
	ofstream out_stream;
	os << "intermediateResults.dat";
	out_stream.precision(15);
	out_stream << std::scientific;
	out_stream.open(os.str(), std::ofstream::out/* | std::ofstream::app*/);

	for (int i = 0; i < solution.size(); i++) {
		out_stream << solution[i] << std::endl;
	}
	out_stream.close();

	return;
}
