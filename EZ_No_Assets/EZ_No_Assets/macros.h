#pragma once
#define SHFT3(a,b,c,d) ((a)=(b),(b)=(c),(c)=(d))
#define SHFT2(a,b,c) ((a)=(b),(b)=(c))
#define SIGN(a) (((a)<0)?-1:1)
#define MAX(a,b) (((a)>(b))?(a):(b))
#define ABS(a) (((a)<0)?(-(a)):(a))
#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))
