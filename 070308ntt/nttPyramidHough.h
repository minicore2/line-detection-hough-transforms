#pragma once
#include "cv.h"
#include "highgui.h"

class nttPyramidHough
{
public:
	nttPyramidHough(void);
public:
	~nttPyramidHough(void);

	//
	void run(IplImage* img, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax);

	int generateSinThetaAndCosThetaList(double thetaResolution);
	void accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho);
};
