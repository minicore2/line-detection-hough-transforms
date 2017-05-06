#pragma once
#include "cv.h"
#include "highgui.h"
#include "nttUtil.h"

#define DOUBLE_PRECISION (2)
#define PI (3.1415927) 
#define HALF_PI (1.5707963) 

class nttStandardHough
{
public:
	nttStandardHough(void);
public:
	~nttStandardHough(void);

	// main method
	void run(IplImage* src, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax);

private:
	double* sinThetaListStd;	
	double* cosThetaListStd;

	int *accDataStd;
	int accColsStd;

	nttUtil utilStd;
	
	//
	int generateSinThetaAndCosThetaList(double thetaResolution);
	void accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho);
	bool checkLocalMaximumWithWindow(int winHeight, int winWidth, int numOfTheta, int numOfRho, int rowOfCenter, int columnOfCenter, int threshold);

	static int cmp_func( const void* _a, const void* _b, void* userdata ) {
		CvRect* a = (CvRect*)_a;
		CvRect* b = (CvRect*)_b;
		int accValueA = a->x;
		int accValueB = b->x;
		return accValueB - accValueA;
	}
};
