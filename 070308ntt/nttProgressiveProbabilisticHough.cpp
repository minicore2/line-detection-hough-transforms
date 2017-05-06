#include "StdAfx.h"
#include "nttProgressiveProbabilisticHough.h"
#include "cv.h"
#include "highgui.h"
#include "nttUtil.h"

// ----- for computation -----
#define DOUBLE_PRECISION (7)
#define PI (3.1415927) // 3.1415926535897932384626433832795
#define HALF_PI (1.5707963) // 1.5707963267948966192313216916398
// ---------------------------

double* sinThetaListThres;	
double* cosThetaListThres;

CvMemStorage* storageThres;
CvSeq* onPixelsThres;

int *accDataThres;
int accColsThres;

uchar *copiedDataThres;

CvPoint line1Thres[2] = {{0,0}, {0,0}};
CvPoint line2Thres[2] = {{0,0}, {0,0}};

nttUtil utilThres;

// ------------

nttProgressiveProbabilisticHough::nttProgressiveProbabilisticHough(void) { }
nttProgressiveProbabilisticHough::~nttProgressiveProbabilisticHough(void) { }

/**********************************************************************************
 * 1. Check input image, if empty => finish
 * 2. Update accumulator with a pixel randomly selected from image
 * 3. Remove selected pixel from image
 * 4. Check if the highest peak in accumulator 
 *	  (that was modified by new pixel) > threshold ? If not => step 1
 * 5. Look along a corridor specified by the peak in accumulator & find longest 
 *	  segment that either is continuous or exhibits a gap not exceeding a given 
 *    threshold
 * 6. Remove pixels in segment from image
 * 7. "Unvote" from accumulator all pixels from the line that previously voted
 * 8. If line segment > minimum length, add it into output list
 * 9. Go to step 1
 **********************************************************************************/
void nttProgressiveProbabilisticHough::run(IplImage* src, CvSeq *linesFound, double rhoResolution, double thetaResolution, 
										int threshold, int minLineLength, int maxGap) {
	IplImage* copiedImageThres = cvCloneImage(src);
	int height = copiedImageThres->height, width = copiedImageThres->width, step = copiedImageThres->widthStep;
	copiedDataThres = (uchar *)copiedImageThres->imageData;
	
	int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);
	int numOfRho = cvRound( ((width + height)*2 + 1) / rhoResolution );
	int halfOfNumOfRho = cvRound(numOfRho / 2);

	// init accumulator
	CvMat* accumulatorThres = cvCreateMat(numOfTheta, numOfRho, CV_32SC1);
	cvZero(accumulatorThres);
	accColsThres = accumulatorThres->cols;
	accDataThres = accumulatorThres->data.i;

	getAllOnPixels(height, width, step); // collect ON pixels (in src image) to onPixelsThres

	// check input image, if empty => finish
	double oneSubtractP = 1 - (double)1 / (double)numOfRho;
	while (onPixelsThres->total > 0) {
		// select a random ON pixel from onPixelsThres
		int randomNum = utilThres.generateRandomNumberOpenCV(onPixelsThres->total);
		CvPoint* point = (CvPoint*)cvGetSeqElem(onPixelsThres, randomNum);
		int x = point->x, y = point->y;

		accumulateOrDeaccumulate(x, y, rhoResolution, numOfTheta, halfOfNumOfRho, 1); // update accumulator with the selected random pixel
		cvSeqRemove(onPixelsThres, randomNum); // 3. remove selected pixel from image

		// find highest peak in accumulator		
		int highestPeakValue, thetaIndex, rhoIndex, numOfAccumulated, numOfBinAccumulated;
		numOfBinAccumulated = numOfAccumulated = highestPeakValue = thetaIndex = rhoIndex = 0;
		for(int i = 0; i < numOfTheta; i++) 
			for(int j = 0; j < numOfRho; j++) {
				int accValue = accDataThres[i*accColsThres + j];
				if (accValue > 0) {
					numOfBinAccumulated++;
					numOfAccumulated += accValue;
					if (highestPeakValue < accValue) {
						highestPeakValue = accValue; // accValue is a double
						thetaIndex = i; // thetaIndex = row
						rhoIndex = j;
					}
				}
			}

		double mean = (double)numOfAccumulated / (double)numOfBinAccumulated;
		double standardDeviation = sqrt(mean * oneSubtractP);
		if (calculateIntegralOfGaussianDistributionFunction(0, highestPeakValue, mean, standardDeviation) < 0.92) continue;

		//if (highestPeakValue < threshold) continue; // 4.2 check if highest peak > threshold 
				
		// look along a corridor specified by the peak in the accumulator
		//    (from the current point, walk in each direction along the found line and extract the line segment)
		double rho = (rhoIndex - halfOfNumOfRho) * rhoResolution;
		findPixelsLeftAndRight(x, y, thetaIndex, rho, rhoResolution, numOfTheta, halfOfNumOfRho, step, width, height, maxGap, 1);
		findPixelsLeftAndRight(x, y, thetaIndex, rho, rhoResolution, numOfTheta, halfOfNumOfRho, step, width, height, maxGap, -1);

		// if line segment > min length, add it into output list (connect 2 lines,which run in 2 directions, together)
		int startX = line1Thres[1].x;	int startY = line1Thres[1].y;
		int endX = line2Thres[1].x;		int endY = line2Thres[1].y;
		if (abs(startX - endX) >= minLineLength || abs(startY - endY) >= minLineLength) {
			CvRect lr = { startX, startY, endX, endY }; // <-- CvRect = line
			cvSeqPush(linesFound, &lr);
		}
	}
	cvReleaseMat(&accumulatorThres);
	cvReleaseImage(&copiedImageThres);	
	cvReleaseMemStorage(&storageThres);
	delete sinThetaListThres;
	delete cosThetaListThres;
}

// ================= private =================

int nttProgressiveProbabilisticHough::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// 180 degree [ -90 ; +90 ] ~ CV_PI radian [ -PI/2 ; PI/2 ]
	int numOfTheta = cvRound(PI / thetaResolution); // cvRound: convert floating-point number to integer
	sinThetaListThres = new double[numOfTheta];
	cosThetaListThres = new double[numOfTheta];
	for (int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI;
		sinThetaListThres[i] = utilThres.roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaListThres[i] = utilThres.roundDouble(cos(theta), DOUBLE_PRECISION);
	}
	return numOfTheta;
}

void nttProgressiveProbabilisticHough::getAllOnPixels(int height, int width, int step) {
	storageThres = cvCreateMemStorage(0);
	onPixelsThres = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), storageThres);
	for(int i = 0; i < height; i++) 
		for(int j = 0; j < width; j++) { 
			int position = i*step + j;
			if (copiedDataThres[position] > 0) {
				CvPoint p = {j, i}; // (x, y)
				cvSeqPush(onPixelsThres, &p);
			}
		}
}

void nttProgressiveProbabilisticHough::accumulateOrDeaccumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho, int accOrDeacc) {
	for (int i = 0; i < numOfTheta; i++) {
		double rho = x*cosThetaListThres[i] + y*sinThetaListThres[i];
		int indexRho = cvRound( rho / rhoResolution ) + halfOfNumOfRho;
		double remainder = rho - (indexRho - halfOfNumOfRho) * rhoResolution;
		if (remainder > rhoResolution/2) indexRho++;
		accDataThres[i*accColsThres + indexRho] += accOrDeacc;
	}
}

void nttProgressiveProbabilisticHough::removeOnPixel(int x, int y) {
	for (int i = 0; i < onPixelsThres->total; i++) {
		CvPoint* point = (CvPoint*)cvGetSeqElem(onPixelsThres, i);
		if (point->x == x && point->y == y) {
			cvSeqRemove(onPixelsThres, i);
			break;
		}
	}
}

bool nttProgressiveProbabilisticHough::isOnPixelsContain(int x, int y) {
	for (int i = 0; i < onPixelsThres->total; i++) {
		CvPoint* point = (CvPoint*)cvGetSeqElem(onPixelsThres, i);
		if (point->x == x && point->y == y) return true;
	}
	return false;	
}

void nttProgressiveProbabilisticHough::findPixelsLeftAndRight(int x, int y, int thetaIndex, double rho, double rhoResolution, int numOfTheta, 
												   int halfOfNumOfRho, int step, int width, int height, int maxGap, int leftOrRight) {
	int currentX = x;	int currentY = y;
	bool foundOnPixelInCorridor = false;
	int maxX;

	while (true) {
		int upperY = -min(maxGap/2, currentY);
		int lowerY = min(maxGap/2, height - currentY);		
		if (leftOrRight == 1) maxX = min(maxGap, width - currentX); // leftOrRight = -1 means left, otherwise, right
		else maxX = min(maxGap, currentX);

		for (int i = upperY; i <= lowerY; i++) {
			for (int j = 0; j < maxX; j++) {
				if (i == 0 && j == 0) continue;
				int nextX = currentX + leftOrRight * j; //
				int nextY = currentY + i; 		
				int nextPosition = nextY*step + nextX; 
				if (copiedDataThres[nextPosition] > 0 
							&& fabs(rho - (nextX*cosThetaListThres[thetaIndex] + nextY*sinThetaListThres[thetaIndex])) < rhoResolution) {
					// 7. "Unvote" from accumulator all pixels from the line that previously voted
					if (!isOnPixelsContain(nextX, nextY)) 
						accumulateOrDeaccumulate(nextX, nextY, rhoResolution, numOfTheta, halfOfNumOfRho, -1);

					// 6. remove pixels in segment from image
					removeOnPixel(nextX, nextY); 
					copiedDataThres[nextPosition] = 0;

					currentX = nextX;	currentY = nextY;
					foundOnPixelInCorridor = true;
					break;
				} 
			}
			if (foundOnPixelInCorridor) break;			
		}
		if (!foundOnPixelInCorridor) break;	
		foundOnPixelInCorridor = false; 
	}

	if (leftOrRight == 1) {
		line1Thres[0].x = x;		line1Thres[0].y = y;
		line1Thres[1].x = currentX;	line1Thres[1].y = currentY;
	} else {
		line2Thres[0].x = x;		line2Thres[0].y = y;
		line2Thres[1].x = currentX;	line2Thres[1].y = currentY;
	}
}

double nttProgressiveProbabilisticHough::calculateGaussianDistributionFunction(double x, double mean, double standardDeviation) {
	double xMinusMean = x - mean;
	double variance = standardDeviation*standardDeviation;
	return exp( -xMinusMean*xMinusMean / (2*variance) ) / ( sqrt(2*PI) * standardDeviation );
	//return exp( -xMinusMean*xMinusMean / (2*variance) );
}

double nttProgressiveProbabilisticHough::calculateIntegralOfGaussianDistributionFunction(int lowerThreshold, int upperThreshold, 
																			  double mean, double standardDeviation) {
	const double interval = 0.01;
	double result = 0;
	int numOfInterval = (int) ((upperThreshold - lowerThreshold) / interval);
	for (int i = 0; i < numOfInterval; i++) { 
		double xi = lowerThreshold + i * interval;
		double xi_1 = xi + interval;
		double yi = calculateGaussianDistributionFunction(xi, mean, standardDeviation);
		double yi_1 = calculateGaussianDistributionFunction(xi_1, mean, standardDeviation);
		double area = (yi + yi_1) * interval / 2;
		result += area;
	}
	return result;
}
