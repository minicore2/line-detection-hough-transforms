#pragma once
#include "cv.h"
#include "highgui.h"

class nttProgressiveProbabilisticHough
{
public:
	nttProgressiveProbabilisticHough(void);
public:
	~nttProgressiveProbabilisticHough(void);

	//
	void run(IplImage* src, CvSeq *linesFound, double rhoResolution, 
		double thetaResolution, int threshold, int minLineLength, int maxGap);

	// === private ===
	int generateSinThetaAndCosThetaList(double thetaResolution);
	void getAllOnPixels(int height, int width, int step);
	void accumulateOrDeaccumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho, int accOrDeacc);
	void removeOnPixel(int x, int y);
	bool isOnPixelsContain(int x, int y);
	void findPixelsLeftAndRight(int x, int y, int thetaIndex, double rho, double rhoResolution, int numOfTheta, 
								int halfOfNumOfRho, int step, int width, int height, int maxGap, int leftOrRight);

	double calculateGaussianDistributionFunction(double x, double mean, double standardDeviation);
	double calculateIntegralOfGaussianDistributionFunction(int lowerThreshold, int upperThreshold, 
		double mean, double standardDeviation);
};
