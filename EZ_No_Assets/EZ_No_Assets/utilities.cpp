#include "utilities.h"

utilities::utilities()
{
}


utilities::~utilities()
{
}

//****************************************************************************80

int utilities::bracket(std::vector<double>xd, double xi)

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
//
//  Discussion:
//
//    We assume XD is sorted.
//
//    If XI is contained in the interval [XD(1),XD(N)], then the returned 
//    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
//
//    If XI < XD(1) then B=0
//    IF XI > XD(N) then B=N-2
//
//    This code implements a version of binary search which is perhaps more
//    understandable than the usual ones.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ND, the number of data values.
//
//    Input, double XD[N], the sorted data.
//
//    Input, double XD, the query value.
//
//    Output, int R8VEC_BRACKET5, the bracket information.
//
{
	int b;
	int l;
	int m;
	int r;
	int nd = (int)xd.size();

	if (xi < xd[0])
	{
		b = 0;
	}
	else if (xd[nd - 1] < xi)
	{
		b = nd-2;
	}
	else
	{
		l = 0;
		r = nd - 1;

		while (l + 1 < r)
		{
			m = (l + r) / 2;
			if (xi < xd[m])
			{
				r = m;
			}
			else
			{
				l = m;
			}
		}
		b = l;
	}

	return b;
}

double utilities::interpolate(const std::vector<double> &x, const std::vector<double> &y, double x_prime) {
	using namespace std;

	double t;
	double yi;

	if (x.size() == 1)
	{
		return y[0];
	}

	int left = bracket(x, x_prime);

	if (x[left] > x_prime && left != 0) {
		std::cout << "Error! utilities.cpp-interpolate() had wrong bracket." << std::endl;
		exit(-1);
	}
	t = (x_prime - x[left]) / (x[left + 1] - x[left]);
	yi = (1.0 - t) * y[left] + t * y[left + 1];

	if (yi != yi) {
		std::cout << "ERROR! utilities.cpp-interpolate() - return value is NaN. "
			<< y[left + 1] << "-" << y[left] << "/(" << x[left + 1] << "-"
			<< x[left] << "),  slope=" << t
			<< " x_prime=" << x_prime << std::endl << std::flush;
		exit(-1);
	}

	return yi;
}
