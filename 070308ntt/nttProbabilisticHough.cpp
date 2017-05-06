#include "StdAfx.h"
#include "nttProbabilisticHough.h"
#include "cv.h"
#include "highgui.h"

CvPoint line1[2] = {{0,0}, {0,0}};
CvPoint line2[2] = {{0,0}, {0,0}};

nttProbabilisticHough::nttProbabilisticHough(void) { }
nttProbabilisticHough::~nttProbabilisticHough(void) { }

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
void nttProbabilisticHough::run(IplImage* src, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int minLineLength, int maxGap) {
	IplImage* copiedImage = cvCloneImage(src);
	int height = copiedImage->height, width = copiedImage->width, step = copiedImage->widthStep;
	copiedData = (uchar *)copiedImage->imageData;
	
	int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);
	int numOfRho = cvRound( ((width + height)*2 + 1) / rhoResolution );
	int halfOfNumOfRho = cvRound(numOfRho / 2);

	// init accumulator
	CvMat* accumulator = cvCreateMat(numOfTheta, numOfRho, CV_32SC1);
	cvZero(accumulator);
	accCols = accumulator->cols;
	accData = accumulator->data.i;

	getAllOnPixels(height, width, step); // collect ON pixels (in src image) to onPixels

	// check input image, if empty => finish
	while (onPixels->total > 0) {
		// select a random ON pixel from onPixels
		int randomNum = util.generateRandomNumberOpenCV(onPixels->total);
		//randomNum = 0;
		CvPoint* point = (CvPoint*)cvGetSeqElem(onPixels, randomNum);		
		
		int x = point->x, y = point->y;

		accumulateOrDeaccumulate(x, y, rhoResolution, numOfTheta, halfOfNumOfRho, 1); // update accumulator with the selected random pixel
		cvSeqRemove(onPixels, randomNum); // remove selected pixel from image

		// find highest peak in accumulator		
		int highestPeakValue = 0, thetaIndex = 0, rhoIndex = 0;
		for(int i = 0; i < numOfTheta; i++) 
			for(int j = 0; j < numOfRho; j++) {
				int accValue = accData[i*accCols + j];
				if (accValue > 0) {
					if (highestPeakValue < accValue) {
						highestPeakValue = accValue; 
						thetaIndex = i; // thetaIndex = row
						rhoIndex = j;
					}
				}
			}
		if (highestPeakValue < threshold) continue; // check if highest peak > threshold 
				
		// look along a corridor specified by the peak in the accumulator
		//    (from the current point, walk in each direction along the found line and extract the line segment)
		double rho = (rhoIndex - halfOfNumOfRho) * rhoResolution; // CACH TINH RHO TU index
		double theta = thetaIndex * thetaResolution - HALF_PI; // CACH TINH THETA TU index

		if (fabs(theta) < 5*HALF_PI/90) {
			findPixelsUpAndDown(x, y, thetaIndex, rho, rhoResolution, numOfTheta, halfOfNumOfRho, step, width, height, maxGap, 1);
			findPixelsUpAndDown(x, y, thetaIndex, rho, rhoResolution, numOfTheta, halfOfNumOfRho, step, width, height, maxGap, -1);
		} else {
			findPixelsLeftAndRight(x, y, thetaIndex, rho, rhoResolution, numOfTheta, halfOfNumOfRho, step, width, height, maxGap, 1);
			findPixelsLeftAndRight(x, y, thetaIndex, rho, rhoResolution, numOfTheta, halfOfNumOfRho, step, width, height, maxGap, -1);
		}

		// if line segment > min length, add it into output list (connect 2 lines,which run in 2 directions, together)
		int startX = line1[1].x;	int startY = line1[1].y;
		int endX = line2[1].x;		int endY = line2[1].y;
		if (abs(startX - endX) >= minLineLength || abs(startY - endY) >= minLineLength) {
			CvRect lr = { startX, startY, endX, endY }; // <-- CvRect = line
			cvSeqPush(linesFound, &lr);
		}
	}

	cvReleaseMat(&accumulator);
	cvReleaseImage(&copiedImage);
	cvReleaseMemStorage(&storage);
	delete sinThetaList;
	delete cosThetaList;
}

// ================= private =================

int nttProbabilisticHough::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// 180 degree [ -90 ; +90 ] ~ CV_PI radian [ -PI/2 ; PI/2 ]
	int numOfTheta = cvRound(PI / thetaResolution) + 1; // cvRound: convert floating-point number to integer
	sinThetaList = new double[numOfTheta];
	cosThetaList = new double[numOfTheta];
	for (int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI; // CACH TINH THETA TU index
		sinThetaList[i] = util.roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaList[i] = util.roundDouble(cos(theta), DOUBLE_PRECISION);
	}
	return numOfTheta;
}

void nttProbabilisticHough::getAllOnPixels(int height, int width, int step) {
	storage = cvCreateMemStorage(0);
	onPixels = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), storage);
	for(int i = 0; i < height; i++) 
		for(int j = 0; j < width; j++) { 
			int position = i*step + j;
			if (copiedData[position] > 0) {
				CvPoint p = {j, i}; // (x, y)
				cvSeqPush(onPixels, &p);
			}
		}
}

void nttProbabilisticHough::accumulateOrDeaccumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho, int accOrDeacc) {
	for (int i = 0; i < numOfTheta; i++) {
		double rho = x*cosThetaList[i] + y*sinThetaList[i];
		int indexRho = cvRound( rho / rhoResolution ) + halfOfNumOfRho;
		double remainder = rho - (indexRho - halfOfNumOfRho) * rhoResolution;
		if (remainder > rhoResolution/2) indexRho++;
		accData[i*accCols + indexRho] += accOrDeacc;
	}
}

void nttProbabilisticHough::removeOnPixel(int x, int y) {
	for (int i = 0; i < onPixels->total; i++) {
		CvPoint* point = (CvPoint*)cvGetSeqElem(onPixels, i);
		if (point->x == x && point->y == y) {
			cvSeqRemove(onPixels, i);
			break;
		}
	}
}

bool nttProbabilisticHough::isOnPixelsContain(int x, int y) {
	for (int i = 0; i < onPixels->total; i++) {
		CvPoint* point = (CvPoint*)cvGetSeqElem(onPixels, i);
		if (point->x == x && point->y == y) return true;
	}
	return false;	
}

void nttProbabilisticHough::findPixelsLeftAndRight(int x, int y, int thetaIndex, double rho, double rhoResolution, int numOfTheta, 
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
			for (int j = 1; j < maxX; j++) {
			//for (int j = 0; j < maxX; j++) {
				if (i == 0 && j == 0) continue;
				int nextX = currentX + leftOrRight * j; //
				int nextY = currentY + i; 		
				int nextPosition = nextY*step + nextX; 
				if ( copiedData[nextPosition] > 0 && fabs(rho - (nextX*cosThetaList[thetaIndex] + nextY*sinThetaList[thetaIndex])) < rhoResolution ) {
					// "Unvote" from accumulator all pixels from the line that previously voted
					if (!isOnPixelsContain(nextX, nextY)) 
						accumulateOrDeaccumulate(nextX, nextY, rhoResolution, numOfTheta, halfOfNumOfRho, -1); // -1 means deaccumulate

					// remove pixels in segment from image
					removeOnPixel(nextX, nextY); 
					//copiedData[nextPosition] = 0; // neu ko co dong nay se loop bat tan, vi ko xoa diem vua moi xet thi lan sau se xet no nua

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
		line1[0].x = x;			line1[0].y = y;
		line1[1].x = currentX;	line1[1].y = currentY;
	} else {
		line2[0].x = x;			line2[0].y = y;
		line2[1].x = currentX;	line2[1].y = currentY;
	}
}

void nttProbabilisticHough::findPixelsUpAndDown(int x, int y, int thetaIndex, double rho, double rhoResolution, int numOfTheta, 
												int halfOfNumOfRho, int step, int width, int height, int maxGap, int upOrDown) {
	int currentX = x;	int currentY = y;
	bool foundOnPixelInCorridor = false;
	int maxY;

	while (true) {
		int leftX = -min(maxGap/2, currentX);
		int rightX = min(maxGap/2, width - currentX - 1);
		if (upOrDown == 1) maxY = min(maxGap, height - currentY - 1); // upOrDown = -1 means up, otherwise, down
		else maxY = min(maxGap, currentY);

		for (int i = 1; i < maxY; i++) {
			for (int j = leftX; j <= rightX; j++) {
				if (i == 0 && j == 0) continue;
				int nextX = currentX + j; 
				int nextY = currentY + upOrDown * i; 		
				int nextPosition = nextY*step + nextX; 
				if ( copiedData[nextPosition] > 0 && fabs(rho - (nextX*cosThetaList[thetaIndex] + nextY*sinThetaList[thetaIndex])) < rhoResolution ) {
					// "Unvote" from accumulator all pixels from the line that previously voted
					if (!isOnPixelsContain(nextX, nextY)) 
						accumulateOrDeaccumulate(nextX, nextY, rhoResolution, numOfTheta, halfOfNumOfRho, -1); // -1 means deaccumulate

					// remove pixels in segment from image
					removeOnPixel(nextX, nextY); 
					//copiedData[nextPosition] = 0; // neu ko co dong nay se loop bat tan, vi ko xoa diem vua moi xet thi lan sau se xet no nua

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

	if (upOrDown == 1) {
		line1[0].x = x;			line1[0].y = y;
		line1[1].x = currentX;	line1[1].y = currentY;
	} else {
		line2[0].x = x;			line2[0].y = y;
		line2[1].x = currentX;	line2[1].y = currentY;
	}
}
