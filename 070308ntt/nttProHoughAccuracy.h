#pragma once
#include "cv.h"
#include "highgui.h"

class nttProHoughAccuracy
{
public:
	nttProHoughAccuracy(void);
public:
	~nttProHoughAccuracy(void);

	//
	void run(IplImage* src, CvSeq *linesFound, double rhoResolution, double thetaResolution, double threshold, int minLineLength, int maxGap);

	void init();

	// === private ===
	int generateSinThetaAndCosThetaList(double thetaResolution);
	void getAllOnPixels(int height, int width, int step);
	void accumulateOrDeaccumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho, int accOrDeacc);
	void removeOnPixel(int x, int y);
	bool isOnPixelsContain(int x, int y);
	void findPixelsLeftAndRight(int x, int y, double sinThetaPeak, double cosThetaPeak, double rho, double rhoResolution, double thetaResolution,
			int minLineLength, int numOfTheta, int numOfRho, int halfOfNumOfRho, int step, int width, int height, int maxGap, int leftOrRight);

	void spreadOrUnspread(int numOfTheta, int numOfRho, int halfOfNumOfRho, int minLineLength, 
			double thetaResolution, double rhoResolution, int sprOrUnspr);
	double calculateStandardDeviation(int L, double thetaRes, double theta);
	double calculateConvolutionValue(double standardDeviation, double rhoI);

	double calculateAlpha(CvSeq *eightH, CvSeq *eightRho);
	double calculateBeta(CvSeq *eightH, CvSeq *eightRho);
	double calculateGamma(CvSeq *eightRho);
	double calculateLogDeltaH(CvSeq *eightH, CvSeq *eightRho);
};
