#pragma once
class FunctorBase
{
public:
	FunctorBase();
	~FunctorBase();
	virtual double operator()(double x) = 0;
};

