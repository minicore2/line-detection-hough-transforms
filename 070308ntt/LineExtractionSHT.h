#pragma once
#include "cv.h"
#include "highgui.h"
#include "nttUtil.h"

#define DOUBLE_PRECISION (7)
#define PI (3.1415927) 
#define HALF_PI (1.5707963) 

typedef struct LineExtractedSHT {
	int startX;
	int startY;
	int endX;
	int endY;
	double theta;
	double rho;
} LineExtractedSHT;

// ---------

class LineExtractionSHT
{
public:
	LineExtractionSHT(void);
public:
	~LineExtractionSHT(void);

	// main method
	void run(IplImage* src, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, 
			 int linesMax, int maxGap, int lineLength);

private:
	double* sinThetaList;	
	double* cosThetaList;

	double *accData;
	int accCols;

	CvRect* startAndEndPoints;

	nttUtil util;

	//
	int generateSinThetaAndCosThetaList(double thetaResolution);
	void accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho);
	void updateStartAndEndPoints(int position, int x, int y);
	void findPointsBetween(CvSeq* result, int startX, int startY, int endX, int endY);

	void filterByMaxNumberOfLines(CvSeq *linesFound, int linesMax);
	bool checkLocalMaximumWithWindow(int winHeight, int winWidth, int numOfTheta, int numOfRho, int rowOfCenter, int columnOfCenter, int threshold);
	
	// static method
	static int cmp_func( const void* _a, const void* _b, void* userdata ) {
		LineExtractedSHT* lineA = (LineExtractedSHT*)_a;
		LineExtractedSHT* lineB = (LineExtractedSHT*)_b;

		int lengthA = max( abs(lineA->startX - lineA->endX), abs(lineA->startY - lineA->endY) );
		int lengthB = max( abs(lineB->startX - lineB->endX), abs(lineB->startY - lineB->endY) );

		return lengthB - lengthA;
	}
};
