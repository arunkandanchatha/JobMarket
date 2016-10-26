#pragma once
#define D_SCALE (1e10)
#define D_y (1*D_SCALE)

class myDBClass {
public:
	static double myDB;
};

#define D_b (myDBClass::myDB*D_SCALE)

#define D_BETA (1.0 / 1.004)
#define D_RHO (-4)
#define D_FIRM_RHO (1)
#define D_C (0*D_SCALE)
#define D_S(x,y) (0.034)
#define D_ETA (0.28)
#define MAX_SHOCKS_PER_MONTH (10)
#define D_DEATH (0) 
#define D_PROD_INCREASE (1.001)
#define D_TENURE_INCREASE (1)
#if 0
#define DO_TENURE_SOLVE (1)
#endif

#define WAGE_GRID_SIZE (50)
#if WAGE_GRID_SIZE>1
#define DO_ADJUSTMENT_COSTS (1)
#endif
#define D_MAX_ADJ_COST (0.25*D_SCALE)

#if 0
#define DO_HABIT_FORMATION (1)
#endif

#if 0
#define D_MY_MU (15.14)
#define D_MY_ALPHA (1)
#define D_MY_PARAMETER (0.38)
#else
#define D_MY_MU (1)
#define D_MY_ALPHA (1)
#define D_MY_PARAMETER (1.6)
#endif

#define D_SD_CONF_FOR_BIN_EST (3)

#define SHFT3(a,b,c,d) ((a)=(b),(b)=(c),(c)=(d))
#define SHFT2(a,b,c) ((a)=(b),(b)=(c))
#include "adept.h"
using adept::adouble;

inline adouble SIGN(const adouble &a, const adouble &b) {
	return (value(b) >= 0) ? (value(a) >= 0 ? a : (-a)) : (value(a) >= 0 ? (-a) : a);
}
inline double SIGN(const double &a, const double &b)
{
	return (b >= 0) ? (a >= 0 ? a : (-a)) : (a >= 0 ? (-a) : a);
}
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)>(b))?(b):(a))
inline adouble ABS(const adouble &a)
{
	double tempVal = adept::value(a);
	if (tempVal < 0) {
		return -a;
	}
	else {
		return a;
	}
}
inline double ABS(const double &a)
{
	return (a < 0) ? (-a) : a;
}
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

inline int nCr(int n, int r) {
	if (r > n - r) {
		r = n - r;
	}
	double ans = 1.0;
	for (int j = 1; j <= r; ++j) {
		ans = ans*n / j;
		n--;
	}
	return (int)ans;
}

inline int TENURE_INCREASE_EE(int t) {
#ifdef DO_TENURE_SOLVE
	if (t < 0) {
		std::cout << "Error! macros.h - TENURE_INCREASE_E: Must pass a positive tenure level. t=" << t << std::endl;
		exit(-1);
	}
	return t + 1;
#else
	return 0;
#endif
}

inline int TENURE_INCREASE_EU(int t) {
#ifdef DO_TENURE_SOLVE
	if (t < 0) {
		std::cout << "Error! macros.h - TENURE_INCREASE_U: Must pass a positive tenure level. t=" << t << std::endl;
		exit(-1);
	}
	return 0;
#else
	return 0;
#endif
}

inline int TENURE_INCREASE_UE(int t) {
#ifdef DO_TENURE_SOLVE
	if (t < 0) {
		std::cout << "Error! macros.h - TENURE_INCREASE_E: Must pass a positive tenure level. t=" << t << std::endl;
		exit(-1);
	}
	return t;
#else
	return 0;
#endif
}

inline int TENURE_INCREASE_UU(int t) {
#ifdef DO_TENURE_SOLVE
	if (t < 0) {
		std::cout << "Error! macros.h - TENURE_INCREASE_U: Must pass a positive tenure level. t=" << t << std::endl;
		exit(-1);
	}
	return 0;
#else
	return 0;
#endif
}

