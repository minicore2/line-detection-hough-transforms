#pragma once
#include "cv.h"
#include "highgui.h"
#include "nttUtil.h"

#define DOUBLE_PRECISION (7)
#define PI (3.1415927) 
#define HALF_PI (1.5707963) 

class nttStandardHoughNew
{
public:
	nttStandardHoughNew(void);
public:
	~nttStandardHoughNew(void);

	void run(IplImage* src, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax);

private:
	double* sinThetaList;	
	double* cosThetaList;

	double *accData;
	int accCols;

	nttUtil util;

	//
	int generateSinThetaAndCosThetaList(double thetaResolution);
	void accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho);
	
	static int cmp_func( const void* _a, const void* _b, void* userdata ) {
		double accValueA = ((CvPoint3D64f*)_a)->x;
		double accValueB = ((CvPoint3D64f*)_b)->x;
		return (int)(accValueB - accValueA);
	}
};
