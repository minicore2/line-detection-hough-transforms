#pragma once
#include "cv.h"
#include "highgui.h"

class nttUtil
{
public:
	nttUtil(void);
public:
	~nttUtil(void);

	bool isArrayContainAnInteger(CArray<int, int> *cArray, int n);
	int generateRandomNumber(int max);
	int generateRandomNumberOpenCV(int max);
	double roundDouble(double doValue, int nPrecision);
	double factorial(int n);
	double mean(CvSeq *seq);
	double sumOfSquare(CvSeq *seq);
	double variance(CvSeq *seq, double mean);

	// optimization with ASM code
	float addFloat(float f1, float f2);
	void floatToInt(int *int_pointer, float f);
};
