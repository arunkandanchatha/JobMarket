#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <GenEigsSolver.h>
#include <MatOp/SparseGenMatProd.h>
#include <iostream>

#define ABS(x) ((x<0)?-x:x)
#define MAX(x,y) (((x)>(y))?(x):(y))

const int numJobStates = 2;
const int genCount = 1000;
//const int education = genCount + 1;
//const int oneSet = genCount * education;
const int oneSet = genCount * (genCount + 1) / 2;
const int n = numJobStates * oneSet;

#if 0
#define AGG_JOB_LOSS_PROB (.1)
inline double PROB_JOB_LOSS(int gen) {
	double total = (genCount * (genCount + 1)) / 2;
	return (1.0 - ((double)(gen+1)) / total) * genCount * AGG_JOB_LOSS_PROB;
}
#else
#define PROB_JOB_LOSS(x) (0.5)
#endif

#define PROB_FIND_JOB (.5)
#define GEN_OFFSET(gen) (((gen) * ((gen) + 1)) / 2)
#define TO_INDEX(job,gen,edu) ((job)*oneSet+GEN_OFFSET((gen))+(edu))

using namespace Spectra;

int main()
{
#if 1
	Eigen::SparseMatrix<double> M(n, n);
	M.reserve(Eigen::VectorXi::Constant(n, 3));
	for (int i = 0; i < genCount; i++)
	{
		if (i == (genCount - 1)) {
			for (int j = 0; j <=  i; j++) {
				M.insert(TO_INDEX(0, i, j), TO_INDEX(0, 0, 0)) = 1 - PROB_FIND_JOB;
				M.insert(TO_INDEX(0, i, j), TO_INDEX(1, 0, 0)) = PROB_FIND_JOB;

				double pjl = PROB_JOB_LOSS(i);
				M.insert(TO_INDEX(1, i, j), TO_INDEX(0, 0, 0)) = pjl;
				M.insert(TO_INDEX(1, i, j), TO_INDEX(1, 0, 0)) = 1 - PROB_JOB_LOSS(i);
			}
			continue;
		}
		for (int j = 0; j <= i; j++) {
			int row = TO_INDEX(0, i, j);
			int col = TO_INDEX(0, i + 1, j);
			M.insert(row, col) = 1 - PROB_FIND_JOB;
			M.insert(TO_INDEX(0, i, j), TO_INDEX(1, i + 1, j)) = PROB_FIND_JOB;

			M.insert(TO_INDEX(1, i, j), TO_INDEX(0, i + 1, j + 1)) = PROB_JOB_LOSS(i);
			M.insert(TO_INDEX(1, i, j), TO_INDEX(1, i + 1, j + 1)) = 1 - PROB_JOB_LOSS(i);
		}
	}
#endif
#if 0
	std::cout << M;
	exit(-1);
	M = M.transpose();

	// Construct matrix operation object using the wrapper class SparseGenMatProd
	SparseGenMatProd<double> op(M);

	// Construct eigen solver object, requesting the largest three eigenvalues
	GenEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op, 1, 2*genCount);

	// Initialize and compute
	eigs.init();
	int nconv = eigs.compute();

	// Retrieve results
	Eigen::VectorXcd evalues;
	if (eigs.info() == SUCCESSFUL) {
		evalues = eigs.eigenvalues();
		std::cout << "Eigenvalues found:\n" << evalues << std::endl;

		evalues = eigs.eigenvectors(1);
		double total = 0;
		for (int j = 0; j < evalues.size(); j++) {
			total += ABS(evalues(j).real());
		}

#if 1
		std::cout << "Unemployed" << std::endl;
		for (int i = 0; i < genCount; i++) {
			for (int j = 0; j <= i; j++) {
				std::cout << i << "," << j << "," << ABS(evalues(TO_INDEX(0, i, j)).real()) / total * genCount << std::endl;
			}
		}

		std::cout << "Employed" << std::endl;
		for (int i = 0; i < genCount; i++) {
			for (int j = 0; j <= i; j++) {
				std::cout << i << "," << j << "," << ABS(evalues(TO_INDEX(1, i, j)).real()) / total * genCount << std::endl;
			}
		}
#endif
	}
	else {
		std::cout << "Eigenvalues not found:\n" << evalues << std::endl;
		std::cout << eigs.info() << std::endl;
		exit(-1);
	}
#else

#if 0
	Eigen::SparseMatrix<double> M(3, 3);
	M.reserve(Eigen::VectorXi::Constant(3, 3));
	int row = 0;
	int col = 0;
	M.insert(row, col++) = .5;
	M.insert(row, col++) = .25;
	M.insert(row++, col++) = .25;
	col = 0;
	M.insert(row, col++) = .5;
	M.insert(row, col++) = 0;
	M.insert(row++, col++) = .5;
	col = 0;
	M.insert(row, col++) = .25;
	M.insert(row, col++) = .25;
	M.insert(row++, col++) = .5;

	Eigen::SparseMatrix<double> M2(3, 3);
#else
	Eigen::SparseMatrix<double> M2(n, n);
	Eigen::SparseMatrix<double> myTempSparse(n, n);
#endif
	M2 = M;

	Eigen::VectorXd temp(n);
	Eigen::VectorXd oldtemp(n);
	Eigen::VectorXd start(n);

	for (int i = 0; i < genCount; i++) {
		start(TO_INDEX(0, i, 0)) = 1.0 / genCount;
	}

	temp = start;
	oldtemp = temp;

	double oldSize = 100;

	for (int i = 0; i < 100; i++) {
		std::cout << i;
		myTempSparse = (M * M2).pruned(1e-10);

		temp = start.transpose()*myTempSparse;
		double size = MAX(ABS((temp - oldtemp).maxCoeff()), ABS((temp - oldtemp).minCoeff()));
		oldSize = size;
		std::cout << "," << size << std::endl;
		M2 = myTempSparse;

		if (size > oldSize) {
			std::cout << size << " > " << oldSize << std::endl;
//			std::cout << "Temp: " << temp << std::endl;
//			std::cout << "Old:  " << oldtemp << std::endl;
			exit(-1);
		}

		if (size < 1e-6) {
			break;
		}
		oldtemp = temp;
	}

	//std::cout << M2;
	//std::cout << "======================" << std::endl;


	//std::cout << std::endl << start.transpose() * M2;
#endif
	return 0;
}