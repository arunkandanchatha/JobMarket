#include "SimulatedAnnealingInterface.h"
//because of templated class, all this stuff has to be in the header file
#if 0
template<class Solution, class Target>
bool SimulatedAnnealingInterface<Solution, Target>::shouldStop(
	const Solution& solution) {

	if (shouldStopHook(solution)) {
		//		std::cout << std::endl << "STOP REASON: ShouldStopHook" << std::endl;
		return true;
	}
	else if (calcDistanceToTarget(solution) < PRECISION) {
		std::cout << std::endl
			<< "STOP REASON: Distance to target is smaller than the required precision. Solution found."
			<< std::endl;
		return true;
	}
	else {
		return false;
	}

}

template<class Solution, class Target>
double SimulatedAnnealingInterface<Solution, Target>::calcDistanceChange(
	const Solution& oldSolution, const Solution& newSolution) const {

	double oldDistance = calcDistanceToTarget(oldSolution);
	double newDistance = calcDistanceToTarget(newSolution);
	assert(oldDistance >= 0 && newDistance >= 0); // distances are always positive
	return (newDistance - oldDistance);

}

template<class Solution, class Target>
double SimulatedAnnealingInterface<Solution, Target>::calcProbability(double change,
	double temp) const { // should be overwritten

	return exp(-1.0 * change / temp);

}

template<class Solution, class Target>
double SimulatedAnnealingInterface<Solution, Target>::calcNewTemp(
	double lastTemp) const {
	return lastTemp * ALPHA;
}

template<class Solution, class Target>
bool SimulatedAnnealingInterface<Solution, Target>::accept(const Solution& oldSolution,
	const Solution& newSolution, double temp) const {

	//std::cout << "Testing acceptance: Old Solution: " << oldSolution << " New Solution: " << newSolution << std::endl;
	double change = calcDistanceChange(oldSolution, newSolution);
	//std::cout << "\t ->Change: " << change << std::endl;
	if (change < 0) {
		//		std::cout << "Result is better so: ";
		return true;
	}
	else {
		double probability = calcProbability(change, temp);
		//DEBUG:BEGIN
		if (!(probability >= 0 && probability <= 1)) {
			std::cout << "Probability wrong: " << probability << " temp: "
				<< temp << std::endl;
			std::cin.get();
		}
		//DEBUG:END
		assert(probability >= 0 && probability <= 1); // probability has to be checked
													  //		std::cout << "Probability we're gonna accept: " << probability << " result: ";
		return (rand() % ((int)(1.0 / probability)) == 0);
	}

}

template<class Solution, class Target>
Solution* SimulatedAnnealingInterface<Solution, Target>::solve() {

	//printStatus(*solution, temp);
	long int counter = 0;
	while (!shouldStop(*solution)) {
		//std::cout<<counter++<<std::endl<<std::flush;
		Solution *solPtr = solution;
		//std::cout << "GET NEIGHBOR" << std::endl << std::flush;
		Solution* newSolution = giveRandomNeighbour(*solution);
		//std::cout << "CHECK SOLUTION" << std::endl << std::flush;
		if (accept(*solution, *newSolution, temp)) {
			solution = newSolution;
			//std::cout << "ACCEPTED" << std::endl<<std::flush;
		}
		else {
			//std::cout << "REJECTED" << std::endl<<std::flush;
			delete newSolution;
			continue;
		}
		//std::cout << "energy(solution)" << std::endl<<std::flush;
		double newEn = calcDistanceToTarget(*solution);
		//std::cout << "energy(best)" << std::endl<<std::flush;
		double bestEn = calcDistanceToTarget(*best);
		if (newEn < bestEn) {
			if (solPtr == best) {
				//std::cout << "delete solPtr" << std::endl<<std::flush;
				delete solPtr;
			}
			else {
				//std::cout << "delete solPtr2" << std::endl<<std::flush;
				delete solPtr;
				//std::cout << "delete best" << std::endl<<std::flush;
				delete best;
			}
			best = solution;
			std::cout << "simAnneal.hpp-solve(): " << counter++ << ":" << (*solution)[0] << ":" << newEn << std::endl;
		}
		else {
			if ((solPtr != best) && (solPtr != solution)) {
				//std::cout << "delete solPtr3" << std::endl<<std::flush;
				delete solPtr;
			}
			else {
				//std::cout << "don't delete solPtr3" << std::endl<<std::flush;
			}
		}
		//std::cout << "pointer deletion done?" << std::endl<<std::flush;

		//std::cout << "\n*****************\n" << std::endl << std::flush;
		if (solution == NULL) {
			std::cout << "ERROR! simAnneal.hpp-solve(): Solution disappeared!" << std::endl << std::flush;
			exit(-1);
		}
		//printStatus(*solution, temp);
		temp = calcNewTemp(temp);
		//std::cout << "NEW TEMP " << temp << std::endl << std::flush;
		assert(temp >= 0); // only positive temperatures are allowed!
						   //std::cout << "ASSERT COMPLETE." << std::endl << std::flush;
	}
	//std::cout << "\n\n\n*************************************************************\n\n"
	//<< "We're done, Solution: \n\n" << (*solution)[0] << "\n\n*************************************************************\n" << std::endl;

	//std::cout<<"in solve"<<(best==solution)<<std::endl<<std::flush;

	if (solution != best) {
		delete solution;
	}
	delete TARGET;

	//	std::cout << "\n\n\nPress enter to continue..." << std::endl;

	//	std::cin.get();
	return best;
}

template<class Solution, class Target>
bool SimulatedAnnealingInterface<Solution, Target>::shouldStopHook(
	const Solution& solution) {
	return false;
}

template<class Solution, class Target>
void SimulatedAnnealingInterface<Solution, Target>::printStatus(const Solution& solution,
	double temp) {
	std::cout << "Current solution: " << solution[0] << " at Temp: " << temp << std::endl << std::flush;
}

template<class Solution, class Target>
SimulatedAnnealingInterface<Solution, Target>::SimulatedAnnealingInterface(
	const Solution& startSolution, const Target& target, double starttemp,
	double precision) :
	TARGET(new Target(target)), solution(new Solution(startSolution)), temp(
		starttemp), PRECISION(precision) {
	assert(starttemp >= 0); // only positive temperatures are allowed!
	best = solution;
}

template<class Solution, class Target>
SimulatedAnnealingInterface<Solution, Target>::~SimulatedAnnealingInterface() {
}

#endif