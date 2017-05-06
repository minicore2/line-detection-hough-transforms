#include "StdAfx.h"
#include "LineExtraction.h"

LineExtraction::LineExtraction(void) { }
LineExtraction::~LineExtraction(void) { }

void LineExtraction::run(IplImage* img, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax, int maxGap, int lineLength, 
						 int lineInLengthFrom, int lineInLengthTo, double thetaFrom, double thetaTo) {
	int height = img->height, width = img->width, step = img->widthStep;
	const uchar* imgData = (uchar *)img->imageData;

	int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);
	int numOfRho = cvRound( ((width + height)*2 + 1) / rhoResolution );
	int halfOfNumOfRho = cvRound(numOfRho / 2);

	CvMemStorage* tmpLinesFoundStorage = cvCreateMemStorage(0);
	CvSeq* tmpLinesFound = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(LineExtracted), tmpLinesFoundStorage);

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
			int deltaI = 1, deltaJ = 1;
			if (i == 0 || i == numOfTheta - 1) deltaI = 0;
			if (j == 0 || j == numOfRho - 1) deltaJ = 0;

			double accValue = accData[i*accCols + j];
			double accValueUp = accData[(i - deltaI)*accCols + j];
			double accValueLeft = accData[i*accCols + j - deltaJ];
			double accValueRight = accData[i*accCols + j + deltaJ];			
			double accValueDown = accData[(i + deltaI)*accCols + j];
			
			//if (accValue > threshold) {
			/*if ( (accValue > threshold) && (accValue > accValueUp) && (accValue >= accValueDown)
										&& (accValue > accValueLeft) && (accValue >= accValueRight)) { // find local maximum
				CvPoint3D64f s = { accValue, i, j }; 
				cvSeqPush(satisfyList, &s);
			}*/

			if (checkLocalMaximumWithWindow(5, 5, numOfTheta, numOfRho, i, j, threshold)) {
				CvPoint3D64f s = { accValue, i, j }; 
				cvSeqPush(satisfyList, &s);
			}

		}

	// take in the satisfyList a number of line(theta, rho) to return
	for (int i = 0; i < satisfyList->total; i++) {
		CvPoint3D64f* s = (CvPoint3D64f*)cvGetSeqElem(satisfyList, i);
		int thetaIndex = (int)(s->y), rhoIndex = (int)(s->z);

		double cosTheta = cosThetaList[thetaIndex], sinTheta = sinThetaList[thetaIndex];
		double rho = (rhoIndex - halfOfNumOfRho) * rhoResolution;

		//double theta = thetaIndex * thetaResolution - HALF_PI;
		double theta = util.roundDouble(thetaIndex * thetaResolution - HALF_PI, DOUBLE_PRECISION);
		
		// find the finite lines
		bool cont = false;
		CvPoint lastPixel = { -1, -1 };
		int beginX = -1, beginY = -1;
		for(int m = 0; m < height; m++) { // y
			for(int n = 0; n < width; n++) { // x
				int lastX = lastPixel.x, lastY = lastPixel.y;
				if (imgData[m*step + n] > 0 && ( fabs(rho - (n*cosTheta + m*sinTheta)) <= rhoResolution) ) {					
					if (lastX == -1 && lastY == -1) {
						beginX = n;
						beginY = m;
						lastPixel = cvPoint(n, m);
					} else {
						if (max(abs(n - lastX), abs(m - lastY)) <= maxGap) {
							lastPixel = cvPoint(n, m);	
							cont = true;
						} else {
							/*if ( max(abs(lastX - beginX), abs(lastY - beginY)) >= lineLength ) {
								LineExtracted lr = { beginX, beginY, lastX, lastY, theta, rho };
								cvSeqPush(tmpLinesFound, &lr);
								n--;
							}*/
							lastPixel = cvPoint(-1, -1); // reset
							cont = false;
						}

						if (cont) {
							if ( max(abs(lastX - beginX), abs(lastY - beginY)) >= lineLength ) {
								int size = tmpLinesFound->total;
								if (size > 0) {
									LineExtracted* l = (LineExtracted*)cvGetSeqElem(tmpLinesFound, size - 1);
									if (max(abs(l->endX - lastX), abs(l->endY - lastY)) <= maxGap)
										cvSeqRemove(tmpLinesFound, size - 1);
								}
								LineExtracted lr = { beginX, beginY, lastX, lastY, theta, rho };
								cvSeqPush(tmpLinesFound, &lr);
							}
						}
					}
				}
			}
		}
	}


	// sort linesFound in the descending order of length
	cvSeqSort( tmpLinesFound, cmp_func, 0 /* userdata is not used here */ );

	// filter by conditions
	filterByMaxNumberOfLines(linesFound, tmpLinesFound, linesMax); // get the maximum number of lines

/*	if (lineInLengthFrom != -1 && lineInLengthTo != -1 ) 
		filterByLineLength(linesFound, lineInLengthFrom, lineInLengthTo);
	
	if (thetaFrom != -1 && thetaTo != -1 ) 
		filterByTheta(linesFound, thetaFrom, thetaTo);
*/
	// release
	cvReleaseMemStorage(&satisfyListStorage);
	cvReleaseMat(&accumulator);
	delete sinThetaList;
	delete cosThetaList;
}

// ================= private =================

int LineExtraction::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// 180 degree [ -90 ; +90 ] ~ CV_PI radian [ -PI/2 ; PI/2 ]
	//int numOfTheta = cvRound(PI / thetaResolution) + 1; // cvRound: convert floating-point number to integer
	int numOfTheta = cvRound(PI / thetaResolution) - 1;
	//int numOfTheta = cvRound(PI / thetaResolution);

	sinThetaList = new double[numOfTheta];
	cosThetaList = new double[numOfTheta];
	for (int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI;
		sinThetaList[i] = util.roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaList[i] = util.roundDouble(cos(theta), DOUBLE_PRECISION);
		//sinThetaList[i] = sin(theta);
		//cosThetaList[i] = cos(theta);
	}
	return numOfTheta;
}

void LineExtraction::accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho) {
	for (int i = 0; i < numOfTheta; i++) {
		double rhoExact = x*cosThetaList[i] + y*sinThetaList[i];
		int indexRhoLow = cvRound( rhoExact / rhoResolution ) + halfOfNumOfRho;
		double rhoLow = (indexRhoLow - halfOfNumOfRho) * rhoResolution;

		int position = i*accCols + indexRhoLow;

		// acccumulate to parameter space
		accData[position] += fabs((rhoLow + rhoResolution) - rhoExact); // rhoHigh - rhoExact
		accData[position + 1] += fabs(rhoExact - rhoLow);

		// --- +1 ---
		/*double rho = x*cosThetaList[i] + y*sinThetaList[i];
		int indexRho = cvRound( rho / rhoResolution ) + halfOfNumOfRho;
		double remainder = rho - (indexRho - halfOfNumOfRho) * rhoResolution;
		if (remainder > rhoResolution/2) indexRho++;
		accData[i*accCols + indexRho] += 1;*/
	}
}

void LineExtraction::filterByMaxNumberOfLines(CvSeq *linesFound, CvSeq *tmpLinesFound, int linesMax) {
	for (int i = 0; i < min(linesMax, tmpLinesFound->total); i++)
		cvSeqPush(linesFound, (LineExtracted*)cvGetSeqElem(tmpLinesFound, i) );
}

void LineExtraction::filterByLineLength(CvSeq *linesFound, int lowerThrehold, int upperThreshold) {
	for (int i = 0; i < linesFound->total; i++) {
		LineExtracted* lr = (LineExtracted*)cvGetSeqElem(linesFound, i);
		int lineLength = max( abs(lr->startX - lr->endX), abs(lr->startY - lr->endY) );
		if (lineLength < lowerThrehold || upperThreshold < lineLength) {
			cvSeqRemove(linesFound, i);
			i--;
		}
	}
}

void LineExtraction::filterByTheta(CvSeq *linesFound, double thetaFrom, double thetaTo) {
	for (int i = 0; i < linesFound->total; i++) {
		LineExtracted* lr = (LineExtracted*)cvGetSeqElem(linesFound, i);
		if (lr->theta < thetaFrom || thetaTo < lr->theta) {
			cvSeqRemove(linesFound, i);
			i--;
		}
	}
}

bool LineExtraction::checkLocalMaximumWithWindow(int winHeight, int winWidth, int numOfTheta, int numOfRho, 
												   int rowOfCenter, int columnOfCenter, int threshold) { // 8-connectivity
	double accValue = accData[rowOfCenter*accCols + columnOfCenter];
	if (accValue <= threshold) return false;
	
	int startColumn = -min(winWidth, columnOfCenter);	int endColumn = min(winWidth, numOfRho - columnOfCenter - 1);
	int startRow = -min(winHeight, rowOfCenter);		int endRow = min(winHeight, numOfTheta - rowOfCenter - 1);

	for (int i = startRow; i <= endRow; i++) {
		for (int j = startColumn; j <= endColumn; j++) {
			if (i == 0 && j == 0) continue;

			int newRow = rowOfCenter + i;
			int newColumn = columnOfCenter + j;
			double acc = accData[newRow*accCols + newColumn];
			if ( (newRow < rowOfCenter) || (newRow == rowOfCenter && newColumn < columnOfCenter) ) {
				if (accValue <= acc) 
					return false;
			} else {
				if (accValue < acc) 
					return false;
			}
		}
	}

	// xet truong hop 0 radian << theta && theta << 3.14 radian
	if (winHeight > rowOfCenter) {
		// (1) xet neu theta nam gan canh tren (gan 0 radian)
		int outOfUpBorder = winHeight - (rowOfCenter + 1);
		for (int i = 0; i < outOfUpBorder; i++) {
			for (int j = startColumn; j <= endColumn; j++) {
				int newRow = (numOfTheta - 1) - i;
				int newColumn = columnOfCenter + j;
				double a = accData[newRow*accCols + newColumn];
				if (accValue <= a) 
					return false;
				
			}
		}

		// (2) xet neu theta nam gan canh tren (gan 3.14 radian)
		int outOfBottomBorder = winHeight - (numOfTheta - rowOfCenter);
		for (int i = 0; i < outOfBottomBorder; i++) {
			for (int j = startColumn; j <= endColumn; j++) {
				int newRow = i;
				int newColumn = columnOfCenter + j;
				if (accValue <= accData[newRow*accCols + newColumn]) 
					return false;
			}
		}
	}



	return true;
}
