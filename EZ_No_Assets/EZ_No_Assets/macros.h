#pragma once
#define SHFT3(a,b,c,d) ((a)=(b),(b)=(c),(c)=(d))
#define SHFT2(a,b,c) ((a)=(b),(b)=(c))
//#define SIGN(a) (((a)<0)?-1:1)
template<class T>
inline T SIGN(const T &a, const T &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)>(b))?(b):(a))
template<class T>
inline T ABS(const T &a)
{
	return a < 0 ? -a : a;
}
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
