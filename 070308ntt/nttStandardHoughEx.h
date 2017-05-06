#pragma once

#include "cv.h"
#include "highgui.h"
#include "Util.h"

#define DOUBLE_PRECISION (7)
#define PI (3.1415927) 
#define HALF_PI (1.5707963) 

class nttStandardHoughEx
{
public:
	nttStandardHoughEx(void);
public:
	~nttStandardHoughEx(void);

	// main method
	void run(IplImage* src, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax, 
		int lineLength, int maxGap, int windowWidthEdgeThick, int windowHeightEdgeThick);
	
	void getHoughSpaceInIplImage(IplImage* hs);
	void releaseAccumulator();

private:
	double* sinThetaListStd;	
	double* cosThetaListStd;

	int *accDataStd;
	int accColsStd;	
	CvMat* accumulatorStd; // luu ra ngoai day de co the get IplImage cua acccumulator space
	
	CvMat* startPointCorrespondingToAccumulator;	int *startPointCorrespondingToAccumulatorData;
	CvMat* endPointCorrespondingToAccumulator;		int *endPointCorrespondingToAccumulatorData;	

	//
	int generateSinThetaAndCosThetaList(double thetaResolution);
	void accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho, int width, int height);	
	bool checkLocalMaximumWithWindow(int winHeight, int winWidth, int numOfTheta, int numOfRho, int rowOfCenter, int columnOfCenter, int threshold);
	bool isOnPixelsContain(CvSeq* visitedPoints, CvPoint p);

	void findBetweenPos(int startX, int startY, int endX, int endY, CvSeq* result/*, int beforeIndex*/);

	void findInHorizontalDirection(int x, int y, IplImage* img, CvSeq *linesFound, double rho, double rhoResolution, 
		int maxGap, int thetaIndex, int lineLength);
	void findInVerticalDirection(int x, int y, IplImage* img, CvSeq *linesFound, double rho, double rhoResolution,
		int maxGap, int thetaIndex, int lineLength);

	void find(int x, int y, IplImage* img, CvSeq *linesFound, double rho, double rhoResolution,	int maxGap, int thetaIndex, int lineLength,
		double cosTheta, double sinTheta);

	void findFiniteLines(IplImage* img, CvSeq *linesFound, CvSeq* betPos, int maxGap, int lineLength);

	static int cmp_func( const void* _a, const void* _b, void* userdata ) {
		CvRect* a = (CvRect*)_a;
		CvRect* b = (CvRect*)_b;
		int accValueA = a->x;
		int accValueB = b->x;
		return accValueB - accValueA;
	}
};
