#pragma once
#include "cv.h"
#include "highgui.h"

class nttHoughGreen
{
public:
	nttHoughGreen(void);
public:
	~nttHoughGreen(void);

	//
	void run(IplImage* src, CvSeq *linesFound, double thetaResolution, double threshold);

	int generateSinThetaAndCosThetaList(double thetaResolution);
	void getPixelsOnContour(CvSeq* contour);
	void rotatePixelsOnContour(int thetaIndex);
	bool resampleRho(double deltaRho, int thetaIndex);
	bool isContain(CvSeq *linesFound, int startX, int startY, int endX, int endY);
};
