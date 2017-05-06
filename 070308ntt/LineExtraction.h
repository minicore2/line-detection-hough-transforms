#pragma once
#include "cv.h"
#include "highgui.h"
#include "nttUtil.h"

#define DOUBLE_PRECISION (2)
#define PI (3.1415927) 
#define HALF_PI (1.5707963) 

typedef struct LineExtracted {
	int startX;
	int startY;
	int endX;
	int endY;
	double theta;
	double rho;
} LineExtracted;

// ---------

class LineExtraction
{
public:
	LineExtraction(void);
public:
	~LineExtraction(void);

	// main method
	void run(IplImage* src, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, 
			 int linesMax, int maxGap, int lineLength, int lineInLengthFrom, int lineInLengthTo, double thetaFrom, double thetaTo);

private:
	double* sinThetaList;	
	double* cosThetaList;

	CvMat* accumulator;
	double *accData;
	int accCols;

	nttUtil util;

	//
	int generateSinThetaAndCosThetaList(double thetaResolution);
	void accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho);
	void filterByMaxNumberOfLines(CvSeq *linesFound, CvSeq *tmpLinesFound, int linesMax);
	void filterByLineLength(CvSeq *linesFound, int lowerThrehold, int upperThreshold);
	void filterByTheta(CvSeq *linesFound, double thetaFrom, double thetaTo);

	bool checkLocalMaximumWithWindow(int winHeight, int winWidth, int numOfTheta, int numOfRho, 
												   int rowOfCenter, int columnOfCenter, int threshold);
	
	// static method
	static int cmp_func( const void* _a, const void* _b, void* userdata ) {
		LineExtracted* lineA = (LineExtracted*)_a;
		LineExtracted* lineB = (LineExtracted*)_b;

		int lengthA = max( abs(lineA->startX - lineA->endX), abs(lineA->startY - lineA->endY) );
		int lengthB = max( abs(lineB->startX - lineB->endX), abs(lineB->startY - lineB->endY) );

		return lengthB - lengthA;
	}
};
