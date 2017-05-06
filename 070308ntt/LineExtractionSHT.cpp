#include "StdAfx.h"
#include "LineExtractionSHT.h"

LineExtractionSHT::LineExtractionSHT(void) { }
LineExtractionSHT::~LineExtractionSHT(void) { }

void LineExtractionSHT::run(IplImage* img, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax, int maxGap, int lineLength) {
	int height = img->height, width = img->width, step = img->widthStep;
	const uchar* imgData = (uchar *)img->imageData;

	int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);
	int numOfRho = cvRound( ((width + height)*2 + 1) / rhoResolution );
	int halfOfNumOfRho = cvRound(numOfRho / 2);

	// init startAndEndPoints
	int areaOfParamSpace = numOfTheta * numOfRho;
	startAndEndPoints = new CvRect[areaOfParamSpace];
	for (int i = 0; i < areaOfParamSpace; i++) startAndEndPoints[i] = cvRect( 0, 0, 0, 0 );

	// init satisfyList(accValue, thetaIndex, rhoIndex)
	CvMemStorage* satisfyListStorage = cvCreateMemStorage(0);
	CvSeq* satisfyList = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint3D64f), satisfyListStorage);

	// init accumulator
	CvMat* accumulator = cvCreateMat(numOfTheta, numOfRho, CV_32FC2);
	cvZero(accumulator);
	accCols = accumulator->cols;
	accData = accumulator->data.db;
	
	// select pixels and accumulate to accumulator
	for(int i = 0; i < height; i++) 
		for(int j = 0; j < width; j++) 
			if (imgData[i*step + j] > 0) accumulate(j, i, rhoResolution, numOfTheta, halfOfNumOfRho);

	// find in accumulator the peak is greater than threshold (find local maximum)
	for(int i = 0; i < numOfTheta; i++) 
		for(int j = 0; j < numOfRho; j++) {
			if (checkLocalMaximumWithWindow(2, 2, numOfTheta, numOfRho, i, j, threshold)) {
				int accValue = accData[i*accCols + j];
				CvPoint3D64f s = { accValue, i, j }; 
				cvSeqPush(satisfyList, &s);
			}
		}

	// operate only with startPoint & endPoint
	for (int i = 0; i < satisfyList->total; i++) {
		CvPoint3D64f* s = (CvPoint3D64f*)cvGetSeqElem(satisfyList, i);
		int thetaIndex = (int)(s->y), rhoIndex = (int)(s->z);
		double cosTheta = cosThetaList[thetaIndex], sinTheta = sinThetaList[thetaIndex];
		double rho = (rhoIndex - halfOfNumOfRho) * rhoResolution;
		double theta = thetaIndex * thetaResolution - HALF_PI;

		CvRect points = startAndEndPoints[thetaIndex*accCols + rhoIndex];

		LineExtractedSHT le = { points.x, points.y, points.width, points.height, theta, rho };
		cvSeqPush(linesFound, &le);
	}

	// take in the satisfyList a number of line(theta, rho) to return
	/*for (int i = 0; i < satisfyList->total; i++) {
		CvPoint3D64f* s = (CvPoint3D64f*)cvGetSeqElem(satisfyList, i);
		int thetaIndex = (int)(s->y), rhoIndex = (int)(s->z);
		double cosTheta = cosThetaList[thetaIndex], sinTheta = sinThetaList[thetaIndex];
		double rho = (rhoIndex - halfOfNumOfRho) * rhoResolution;
		double theta = thetaIndex * thetaResolution - HALF_PI;
		
		// find the finite lines
		CvRect startAndEnd = startAndEndPoints[thetaIndex*accCols + rhoIndex];

		// search points "between" startPoint & endPoint
		// look around each between-point to find points satisfying the straight-line condition
		int startX = startAndEnd.x, startY = startAndEnd.y, endX = startAndEnd.width, endY = startAndEnd.height;

		CvMemStorage* betPointsStorage = cvCreateMemStorage(0);
		CvSeq* betPoints = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), betPointsStorage);
		findPointsBetween(betPoints, startX, startY, endX, endY);

		int x0 = -1, y0 = -1, xt = -1, yt = -1; // beginX, beginY, lastX, lastY
		for (int j = 0; j < betPoints->total; j++) { // run along the corridor
			CvPoint* betP = (CvPoint*)cvGetSeqElem(betPoints, j);
			int betX = betP->x, betY = betP->y; 
			
			int scanRangeBegin = 0, scanRangeEnd = 0;
			BYTE cX, cY;
			if (startX - endX == 0) { // vertical
				scanRangeBegin = -min(maxGap, betX);	scanRangeEnd = min(maxGap, width - betX - 1);
				cX = 1;		cY = 0;
			} else { // diagonal & horizontal
				scanRangeBegin = -min(maxGap, betY);	scanRangeEnd = min(maxGap, height - betY - 1);
				cX = 0;		cY = 1;
			}

			for (int k = scanRangeBegin; k <= scanRangeEnd; k++) {				
				int corridorX = betX + k * cX;
				int corridorY = betY + k * cY;
				if ( imgData[corridorY*step + corridorX] > 0 && ( fabs(rho - (corridorX*cosTheta + corridorY*sinTheta)) <= rhoResolution) ) {
					if (xt == -1 && yt == -1) {
						x0 = corridorX;		y0 = corridorY;
						xt = corridorX;		yt = corridorY;
					} else {
						if (max(abs(corridorX - xt), abs(corridorY - yt)) <= maxGap) {
							xt = corridorX;
							yt = corridorY;
						} else {
							if ( max(abs(xt - x0), abs(yt - y0)) >= lineLength ) {
								LineExtractedSHT le = { x0, y0, xt, yt, theta, rho };
								cvSeqPush(linesFound, &le);
								j--; // xet lai diem-tren-truc (betX, betY)
							}
							xt = -1; yt = -1; // reset
						}
					}
				}
			}
		} 
	}*/

	// sort linesFound in the descending order of length
	cvSeqSort( linesFound, cmp_func, 0 /* userdata is not used here */ );

	// filter by conditions
	filterByMaxNumberOfLines(linesFound, linesMax); // get the maximum number of lines

	// release
	cvReleaseMemStorage(&satisfyListStorage);
	cvReleaseMat(&accumulator);
	delete sinThetaList;
	delete cosThetaList;
	delete startAndEndPoints;
}

// ================= private =================

int LineExtractionSHT::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// 180 degree [ -90 ; +90 ] ~ CV_PI radian [ -PI/2 ; PI/2 ]
	int numOfTheta = cvRound(PI / thetaResolution) + 1; // cvRound: convert floating-point number to integer
	sinThetaList = new double[numOfTheta + 1];
	cosThetaList = new double[numOfTheta + 1];
	for (int i = 0; i <= numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI;
		sinThetaList[i] = util.roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaList[i] = util.roundDouble(cos(theta), DOUBLE_PRECISION);
	}
	return numOfTheta;
}

void LineExtractionSHT::accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho) {
	for (int i = 0; i < numOfTheta; i++) {
		
		/* approximation in accumulating
		double rhoExact = x*cosThetaList[i] + y*sinThetaList[i];
		int indexRhoLow = cvRound( rhoExact / rhoResolution ) + halfOfNumOfRho;
		double rhoLow = (indexRhoLow - halfOfNumOfRho) * rhoResolution;

		int position = i*accCols + indexRhoLow;

		// acccumulate to parameter space
		accData[position] += fabs((rhoLow + rhoResolution) - rhoExact); // rhoHigh - rhoExact
		accData[position + 1] += fabs(rhoExact - rhoLow);

		// log start point & end point
		updateStartAndEndPoints(position, x, y);
		updateStartAndEndPoints(position + 1, x, y);
		*/

		// --- +1 ---
		double rho = x*cosThetaList[i] + y*sinThetaList[i];
		int indexRho = cvRound( rho / rhoResolution ) + halfOfNumOfRho;
		double remainder = rho - (indexRho - halfOfNumOfRho) * rhoResolution;
		if (remainder > rhoResolution/2) indexRho++;

		int position = i*accCols + indexRho;
		accData[position] += 1;
		updateStartAndEndPoints(position, x, y);		
	}
}

void LineExtractionSHT::updateStartAndEndPoints(int position, int x, int y) {
	CvRect r = startAndEndPoints[position];
	int startX = r.x, startY = r.y, endX = r.width, endY = r.height;

	if (startX == 0 && startY == 0 && endX == 0 && endY == 0) {
		startAndEndPoints[position] = cvRect(x, y, x, y);
	} else {
		if (startX == endX && startX == x) {
			if (startY < endY) {
				if (y < startY) startAndEndPoints[position] = cvRect(x, y, x, endY);
				if (endY < y) startAndEndPoints[position] = cvRect(x, startY, x, y);
			} else { // startY >= endY
				if (y < endY) startAndEndPoints[position] = cvRect(x, startY, x, y);
				if (startY < y) startAndEndPoints[position] = cvRect(x, y, x, endY);
			}
		} else {
			if (startX < endX) {
				if (x < startX) startAndEndPoints[position] = cvRect(x, y, endX, endY);
				if (endX < x) startAndEndPoints[position] = cvRect(startX, startY, x, y);
			} else { // startX >= endX
				if (x < endX) startAndEndPoints[position] = cvRect(startX, startY, x, y);
				if (startX < x) startAndEndPoints[position] = cvRect(x, y, endX, endY);
			}
		}
	}
}

void LineExtractionSHT::findPointsBetween(CvSeq* result, int startX, int startY, int endX, int endY) {
	int x0 = startX, y0 = startY, xt = endX, yt = endY;
	int meanX = (startX + endX)/2, meanY = (startY + endY)/2;

	if (abs(meanX - x0) > 1) { // || abs(meanY - y0) > 1) {		
		findPointsBetween(result, startX, startY, meanX, meanY);
		CvPoint p = { meanX, meanY };
		cvSeqPush(result, &p);	
		findPointsBetween(result, meanX, meanY, endX, endY);
	} 
}

void LineExtractionSHT::filterByMaxNumberOfLines(CvSeq *linesFound, int linesMax) {
	int linesTotal = linesFound->total;
	for (int i = linesMax; i < linesTotal; i++)
		cvSeqRemove(linesFound, linesFound->total - 1);
}

bool LineExtractionSHT::checkLocalMaximumWithWindow(int winHeight, int winWidth, int numOfTheta, int numOfRho, 
												   int rowOfCenter, int columnOfCenter, int threshold) { // 8-connectivity
	int accValue = accData[rowOfCenter*accCols + columnOfCenter];
	if ( (accValue > threshold) == false) return false;
	
	int startColumn = -min(winWidth, columnOfCenter);	int endColumn = min(winWidth, numOfRho - columnOfCenter - 1);
	int startRow = -min(winHeight, rowOfCenter);		int endRow = min(winHeight, numOfTheta - rowOfCenter - 1);

	for (int i = startRow; i <= endRow; i++) {
		for (int j = startColumn; j <= endColumn; j++) {
			if (i == 0 && j == 0) continue;

			int newRow = rowOfCenter + i;
			int newColumn = columnOfCenter + j;
			int acc = accData[newRow*accCols + newColumn];
			if ( newRow < rowOfCenter || (newRow == rowOfCenter && startColumn < columnOfCenter) )
				if ((accValue > acc) == false) return false;
			else
				if ((accValue >= acc) == false) return false;
		}
	}

	return true;
}