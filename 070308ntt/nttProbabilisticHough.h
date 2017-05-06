#pragma once
#include "cv.h"
#include "highgui.h"
#include "nttUtil.h"

#define DOUBLE_PRECISION (7)
#define PI (3.1415927) 
#define HALF_PI (1.5707963) 

class nttProbabilisticHough
{
public:
	nttProbabilisticHough(void);
public:
	~nttProbabilisticHough(void);
	
	// main method
	void run(IplImage* src, CvSeq *linesFound, double rhoResolution, 
		double thetaResolution, int threshold, int minLineLength, int maxGap);

private:
	double* sinThetaList;
	double* cosThetaList;

	CvMemStorage* storage;
	CvSeq* onPixels;

	int *accData;
	int accCols;

	uchar *copiedData;

	nttUtil util;

	//
	int generateSinThetaAndCosThetaList(double thetaResolution);
	void getAllOnPixels(int height, int width, int step);
	void accumulateOrDeaccumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho, int accOrDeacc);
	void removeOnPixel(int x, int y);
	bool isOnPixelsContain(int x, int y);
	void findPixelsLeftAndRight(int x, int y, int thetaIndex, double rho, double rhoResolution, int numOfTheta, 
								int halfOfNumOfRho, int step, int width, int height, int maxGap, int leftOrRight);
	void findPixelsUpAndDown(int x, int y, int thetaIndex, double rho, double rhoResolution, int numOfTheta, 
		 					 int halfOfNumOfRho, int step, int width, int height, int maxGap, int upOrDown);
};
