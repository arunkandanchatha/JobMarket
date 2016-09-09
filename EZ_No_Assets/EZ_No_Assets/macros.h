#pragma once
#define D_y (1.0)
#define D_b (0.4)
#define D_BETA (1.0 / 1.004)
#define D_RHO (1)
#define D_C (0.21)
#define D_S (0.034)
#define D_ETA (0.28)
#define MAX_SHOCKS_PER_MONTH (10)
#define D_DEATH (0) 
//#define D_PROD_INCREASE (1)//(1.0014)

#if 0
#define D_MY_MU (15.14)
#define D_MY_ALPHA (1)
#define D_MY_PARAMETER (0.38)
#else
#define D_MY_MU (1)
#define D_MY_ALPHA (1)
#define D_MY_PARAMETER (1.6)
#endif

#define SHFT3(a,b,c,d) ((a)=(b),(b)=(c),(c)=(d))
#define SHFT2(a,b,c) ((a)=(b),(b)=(c))
#include "adept.h"
using adept::adouble;

//#define SIGN(a) (((a)<0)?-1:1)
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
