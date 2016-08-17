#pragma once
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
