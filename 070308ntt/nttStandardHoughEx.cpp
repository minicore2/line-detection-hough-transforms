#include "StdAfx.h"
#include "nttStandardHoughEx.h"

nttStandardHoughEx::nttStandardHoughEx(void) {}
nttStandardHoughEx::~nttStandardHoughEx(void) {}

void nttStandardHoughEx::run(IplImage* img, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax,
							 int lineLength, int maxGap, int windowWidthEdgeThick, int windowHeightEdgeThick) {
	int height = img->height, width = img->width, step = img->widthStep;
	const uchar* imgData = (uchar *)img->imageData;

	int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);
	int numOfRho = cvRound( ((width + height)*sqrt(2.0)/2 + 1) / rhoResolution);  // <-- giam kich thuoc param space
	int halfOfNumOfRho = numOfRho / 2;

	// init satisfyListStd(accValue, thetaIndex, rhoIndex, 0 /* not used */ )
	CvMemStorage* satisfyListStorageStd = cvCreateMemStorage(0);
	CvSeq* satisfyListStd = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvRect), satisfyListStorageStd);

	// init accumulator
	accumulatorStd = cvCreateMat(numOfTheta, numOfRho, CV_32SC1);	
	cvZero(accumulatorStd);
	accColsStd = accumulatorStd->cols;
	accDataStd = accumulatorStd->data.i;

	startPointCorrespondingToAccumulator = cvCreateMat(numOfTheta, numOfRho, CV_32SC1);
	cvZero(startPointCorrespondingToAccumulator);	
	startPointCorrespondingToAccumulatorData = startPointCorrespondingToAccumulator->data.i;

	endPointCorrespondingToAccumulator = cvCreateMat(numOfTheta, numOfRho, CV_32SC1);
	cvZero(endPointCorrespondingToAccumulator);	
	endPointCorrespondingToAccumulatorData = endPointCorrespondingToAccumulator->data.i;
	
	// select pixels and accumulate to accumulator
	for(int i=0; i<height; i++) for(int j=0; j<width; j++) if (imgData[i*step+j]>0) accumulate(j, i, rhoResolution, numOfTheta, halfOfNumOfRho, width, height);

	// lay ra cac gia tri nao tich luy trong accumulator ma > threshold (co tim theo local maximum) de tao thanh duong thang
	for(int i = 0; i < numOfTheta; i++) 
		for(int j = 0; j < numOfRho; j++) {
/*
			int deltaI = 1, deltaJ = 1;
			if (i == 0 || i == numOfTheta - 1) deltaI = 0;
			if (j == 0 || j == numOfRho - 1) deltaJ = 0;

			int accValue = accDataStd[i*accColsStd + j];
			int accValueUp = accDataStd[(i - deltaI)*accColsStd + j];			
			int accValueLeft = accDataStd[i*accColsStd + j - deltaJ];
			int accValueRight = accDataStd[i*accColsStd + j + deltaJ];			
			int accValueDown = accDataStd[(i + deltaI)*accColsStd + j];
			//int accValueUpLeft = accDataStd[(i - deltaI)*accColsStd + j - deltaJ];
			//int accValueUpRight = accDataStd[(i - deltaI)*accColsStd + j + deltaJ];
			//int accValueDownLeft = accDataStd[(i + deltaI)*accColsStd + j - deltaJ];			
			//int accValueDownRight = accDataStd[(i + deltaI)*accColsStd + j + deltaJ];

			if ( (accValue > threshold) && (accValue > accValueUp) && (accValue >= accValueDown)
										&& (accValue > accValueLeft) && (accValue >= accValueRight)) { // find local maximum
				CvRect s = { accValue, i, j, 0 }; // 0: not used
				cvSeqPush(satisfyListStd, &s);
			}*/
			
			if (checkLocalMaximumWithWindow(windowWidthEdgeThick, windowHeightEdgeThick, numOfTheta, numOfRho, i, j, threshold)) {
				int accValue = accDataStd[i*accColsStd + j];
				CvRect s = { accValue, i, j, 0 }; // 0: not used
				cvSeqPush(satisfyListStd, &s);
			}
		}

	// sap xep satisfyListStd theo thu tu giam dan cua accValue
	cvSeqSort( satisfyListStd, cmp_func, 0 /* userdata is not used here */ );

	// lay ra n = linesMax phan tu de dua vao linesFound
	for (int i = 0; i < min(linesMax, satisfyListStd->total); i++) {
		CvRect* s = (CvRect*)cvGetSeqElem(satisfyListStd, i);
		int thetaIndex = s->y, rhoIndex = s->width;

		double cosTheta = cosThetaListStd[thetaIndex], sinTheta = sinThetaListStd[thetaIndex];
		double rho = (rhoIndex - halfOfNumOfRho) * rhoResolution;
		
		// tim infiniteLine
		/*double x0 = cosTheta*rho + width, y0 = sinTheta*rho + height; // <-- giam kich thuoc param space		
		CvRect lr = { cvRound(x0 + 1000*(-sinTheta)), cvRound(y0 + 1000*(cosTheta)), 
					  cvRound(x0 - 1000*(-sinTheta)), cvRound(y0 - 1000*(cosTheta)) }; 
		cvSeqPush(linesFound, &lr);*/

		// tim finiteLine
		// find(x, y, img, linesFound, rho, rhoResolution, maxGap, thetaIndex, lineLength, cosTheta, sinTheta);

		int startXY = startPointCorrespondingToAccumulatorData[thetaIndex*accColsStd + rhoIndex];
		CvPoint start = cvPoint(startXY / 1000, startXY % 1000);
		int endXY = endPointCorrespondingToAccumulatorData[thetaIndex*accColsStd + rhoIndex];
		CvPoint end = cvPoint(endXY / 1000, endXY % 1000);				

		CvMemStorage* betPosStorage = cvCreateMemStorage(0);
		CvSeq* betPos = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), betPosStorage);
		//cvSeqPush(betPos, &start);
		//cvSeqPush(betPos, &end);

		findBetweenPos(start.x, start.y, end.x, end.y, betPos);
		findFiniteLines(img, linesFound, betPos, maxGap, lineLength);

		cvReleaseMemStorage(&betPosStorage);
	}

	// release
	cvReleaseMemStorage(&satisfyListStorageStd);
	cvReleaseMat(&accumulatorStd);
	cvReleaseMat(&startPointCorrespondingToAccumulator);
	cvReleaseMat(&endPointCorrespondingToAccumulator);
	delete sinThetaListStd;
	delete cosThetaListStd;
}

// ================= private =================

int nttStandardHoughEx::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// 180 degree [ -90 ; +90 ] ~ CV_PI radian [ -PI/2 ; PI/2 ]
	int numOfTheta = cvRound(PI / thetaResolution) + 1; 
	sinThetaListStd = new double[numOfTheta];
	cosThetaListStd = new double[numOfTheta];
	for (int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI;

		sinThetaListStd[i] = Util::roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaListStd[i] = Util::roundDouble(cos(theta), DOUBLE_PRECISION);

		//sinThetaListStd[i] = Util::roundDouble(sin(theta)*pow(2.0,31.0), 0); // doan ma tinh cos & sin theo cach cua DongKyun
		//cosThetaListStd[i] = Util::roundDouble(cos(theta)*pow(2.0,31.0), 0);
	}
	return numOfTheta;
}

void nttStandardHoughEx::accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho, int width, int height) {
	for (int i = 0; i < numOfTheta; i++) {
		double rho = (x-width)*cosThetaListStd[i] + (y-height)*sinThetaListStd[i]; // <-- giam kich thuoc param space
		int indexRho = cvRound( rho / rhoResolution ) + halfOfNumOfRho;
		double remainder = rho - (indexRho - halfOfNumOfRho) * rhoResolution;
		if (remainder > rhoResolution/2) indexRho++;
		accDataStd[i*accColsStd + indexRho]++;

		int start = startPointCorrespondingToAccumulatorData[i*accColsStd + indexRho];
		int end = endPointCorrespondingToAccumulatorData[i*accColsStd + indexRho];
		if (start == 0 && end == 0) // luc ban dau
			startPointCorrespondingToAccumulatorData[i*accColsStd + indexRho] = x*1000 + y; // encode x,y
		else
			endPointCorrespondingToAccumulatorData[i*accColsStd + indexRho] = x*1000 + y; // encode x,y
	}	

	// doan ma tinh cos & sin theo cach cua DongKyun
/*	double r = Util::roundDouble(rhoResolution*pow(2.0,31.0), 0); 
	for (int i = 0; i < numOfTheta; i++) {
		double rho = (x-width)*cosThetaListStd[i] + (y-height)*sinThetaListStd[i]; // <-- giam kich thuoc param space
		int indexRho = cvRound( rho / r ) + halfOfNumOfRho;
		double remainder = rho - (indexRho - halfOfNumOfRho) * r;
		if (remainder > r/2) indexRho++;
		accDataStd[i*accColsStd + indexRho]++;

		endPointCorrespondingToAccumulatorData[i*accColsStd + indexRho] = x*1000 + y; // encode x,y
	}*/
}

bool nttStandardHoughEx::checkLocalMaximumWithWindow(int winHeight, int winWidth, int numOfTheta, int numOfRho, 
												   int rowOfCenter, int columnOfCenter, int threshold) { // 8-connectivity
	int accValue = accDataStd[rowOfCenter*accColsStd + columnOfCenter];
	if (accValue <= threshold) return false;
	
	int startColumn = -min(winWidth, columnOfCenter);	int endColumn = min(winWidth, numOfRho - columnOfCenter - 1);
	int startRow = -min(winHeight, rowOfCenter);		int endRow = min(winHeight, numOfTheta - rowOfCenter - 1);

	for (int i = startRow; i <= endRow; i++) {
		for (int j = startColumn; j <= endColumn; j++) {
			if (i == 0 && j == 0) continue;

			int newRow = rowOfCenter + i;
			int newColumn = columnOfCenter + j;
			int acc = accDataStd[newRow*accColsStd + newColumn];
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

	// (1) xet neu theta nam gan canh tren (gan 0 radian)
	/*int outOfUpBorder = winHeight - (rowOfCenter + 1);
	for (int i = 0; i < outOfUpBorder; i++) {
		for (int j = startColumn; j <= endColumn; j++) {
			int newRow = (numOfTheta - 1) - i;
			int newColumn = columnOfCenter + j;
			if (accValue <= accDataStd[newRow*accColsStd + newColumn]) return false;
		}
	}

	// (2) xet neu theta nam gan canh tren (gan 3.14 radian)
	int outOfBottomBorder = winHeight - (numOfTheta - rowOfCenter);
	for (int i = 0; i < outOfBottomBorder; i++) {
		for (int j = startColumn; j <= endColumn; j++) {
			int newRow = i;
			int newColumn = columnOfCenter + j;
			if (accValue <= accDataStd[newRow*accColsStd + newColumn]) return false;
		}
	}*/

	return true;
}

bool nttStandardHoughEx::isOnPixelsContain(CvSeq* visitedPoints, CvPoint p) { // not used
	for (int i = 0; i < visitedPoints->total; i++) {
		CvPoint* point = (CvPoint*)cvGetSeqElem(visitedPoints, i);
		if (point->x == p.x && point->y == p.y) return true;
	}
	return false;	
}

void nttStandardHoughEx::findInHorizontalDirection(int x, int y, IplImage* img, CvSeq *linesFound, double rho, double rhoResolution,
												   int maxGap, int thetaIndex, int lineLength) { // not used
	int width = img->width;	int height = img->height;	int step = img->widthStep;
	uchar* imgData = (uchar *)img->imageData;

	int currentX = int(x*1.0);	int currentY = int(y*1.0);
	bool found = false;	

	CvMemStorage* visitedPointsStorage = cvCreateMemStorage(0);
	CvSeq* visitedPoints = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), visitedPointsStorage);

	//while(true) {

		// tim finite line
		while(true) {
			int upperY = -min(maxGap/2, currentY);
			int lowerY = min(maxGap/2, height - currentY);
			int maxLeftX = min(maxGap, currentX);

			for (int i = 0; i < maxLeftX; i++) {
				for (int j = upperY; j <= lowerY; j++) {
					if (i == 0 && j == 0) continue;
					int nextX = currentX - i; //
					int nextY = currentY + j; 		
					int nextPosition = nextY*step + nextX; 

					if (isOnPixelsContain(visitedPoints, cvPoint(nextX, nextY))) break;
					CvPoint p = cvPoint(currentX, currentY);
					cvSeqPush(visitedPoints, &p);

					if (imgData[nextPosition] > 0 
						&& fabs(rho - ((nextX-width)*cosThetaListStd[thetaIndex] + (nextY-height)*sinThetaListStd[thetaIndex])) < rhoResolution ) { // <-- giam kich thuoc param space
					
						currentX = nextX;
						currentY = nextY;
						found = true;

						break;
					}
				}
				if(found) break;
			}
			if(!found) break;
			found = false;
		}

		// neu tim duoc line thoa lineLength thi add vao linesFound
		if (abs(x - currentX) > lineLength || abs(y - currentY) > lineLength) {
			CvRect line = { x, y, currentX, currentY }; 
			cvSeqPush(linesFound, &line);
		}

		// dich chuyen ve be trai maxGap pixels de tim xem con pixels nao nua ko
	//	currentX -= maxGap;
	//	if (currentX < 0) break;
	//}

	cvReleaseMemStorage(&visitedPointsStorage);
}

void nttStandardHoughEx::findInVerticalDirection(int x, int y, IplImage* img, CvSeq *linesFound, double rho, double rhoResolution,
												   int maxGap, int thetaIndex, int lineLength) { // not used
	int width = img->width;	int height = img->height;	int step = img->widthStep;
	uchar* imgData = (uchar *)img->imageData;

	int currentX = int(x*1.0);	int currentY = int(y*1.0);
	bool found = false;	

	CvMemStorage* visitedPointsStorage = cvCreateMemStorage(0);
	CvSeq* visitedPoints = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), visitedPointsStorage);

	//while(true) {

		// tim finite line
		while(true) {
			int leftX = -min(maxGap/2, currentX);
			int rightX = min(maxGap/2, width - currentX);
			int maxUpperY = min(maxGap, currentY);
			
			for (int i = 0; i < maxUpperY; i++) {
				for (int j = leftX; j <= rightX; j++) {
					if (i == 0 && j == 0) continue;
					int nextX = currentX + j; //
					int nextY = currentY - i; 		
					int nextPosition = nextY*step + nextX; 

					if (isOnPixelsContain(visitedPoints, cvPoint(nextX, nextY))) break;
					CvPoint p = cvPoint(currentX, currentY);
					cvSeqPush(visitedPoints, &p);

					if (imgData[nextPosition] > 0 
						&& fabs(rho - ((nextX-width)*cosThetaListStd[thetaIndex] + (nextY-height)*sinThetaListStd[thetaIndex])) < rhoResolution ) { // <-- giam kich thuoc param space
					
						currentX = nextX;
						currentY = nextY;
						found = true;

						break;
					}
				}
				if(found) break;
			}
			if(!found) break;
			found = false;
		}

		// neu tim duoc line thoa lineLength thi add vao linesFound
		if (abs(x - currentX) > lineLength || abs(y - currentY) > lineLength) {
			CvRect line = { x, y, currentX, currentY }; 
			cvSeqPush(linesFound, &line);
		}

		// dich chuyen de tim xem con pixels nao nua ko
	//	currentY -= maxGap;
	//	if (currentY < 0) break;
	//}

	cvReleaseMemStorage(&visitedPointsStorage);
}

void nttStandardHoughEx::find(int x, int y, IplImage* img, CvSeq *linesFound, double rho, double rhoResolution,
							  int maxGap, int thetaIndex, int lineLength, double cosTheta, double sinTheta) { // not used
	int width = img->width;	int height = img->height;	int step = img->widthStep;
	uchar* imgData = (uchar *)img->imageData;

	int currentX = int(x*1.0);	int currentY = int(y*1.0);
	bool found = false;	

	CvMemStorage* visitedPointsStorage = cvCreateMemStorage(0);
	CvSeq* visitedPoints = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), visitedPointsStorage);

	while(true) {

		// tim finite line
		while(true) {
			int leftX = -min(maxGap, currentX);
			int rightX = min(maxGap, width - currentX);
			int upY = -min(maxGap, currentY);
			int downY = min(maxGap, height - currentY);
			
			for (int i = upY; i < downY; i++) {
				for (int j = leftX; j <= rightX; j++) {
					if (i == 0 && j == 0) continue;
					int nextX = currentX + j; //
					int nextY = currentY + i; 		
					int nextPosition = nextY*step + nextX; 

					if (isOnPixelsContain(visitedPoints, cvPoint(nextX, nextY))) break;
					CvPoint p = cvPoint(currentX, currentY);
					cvSeqPush(visitedPoints, &p);

					if (imgData[nextPosition] > 0 
						&& fabs(rho - ((nextX-width)*cosThetaListStd[thetaIndex] + (nextY-height)*sinThetaListStd[thetaIndex])) < rhoResolution ) { // <-- giam kich thuoc param space
					
						currentX = nextX;
						currentY = nextY;
						found = true;

						break;
					}
				}
				if(found) break;
			}
			if(!found) break;
			found = false;
		}

		// neu tim duoc line thoa lineLength thi add vao linesFound
		if (abs(x - currentX) > lineLength || abs(y - currentY) > lineLength) {
			CvRect line = { x, y, currentX, currentY }; 
			cvSeqPush(linesFound, &line);
		}

		// dich chuyen de tim xem con pixels nao nua ko
		if (fabs(sinTheta) < 0.08) { // horizontal
			currentX -= maxGap;
			if (currentX < 0) break;
		} else {
			double cs = -Util::roundDouble(cosTheta, 2) / Util::roundDouble(sinTheta, 2);
			if (cs == 0) { // vertical
				currentY -= maxGap;
				if (currentY < 0) break;
			} else if (cs < 0) { // high to low direction
				currentX -= maxGap;
				currentY -= maxGap;
				if (currentX < 0 || currentY < 0) break;
			} else { // low to high direction
				currentX += maxGap;
				currentY -= maxGap;
				if (currentX >= width || currentY < 0) break;
			}
		}
		
	}

	cvReleaseMemStorage(&visitedPointsStorage);
}

void nttStandardHoughEx::findBetweenPos(int startX, int startY, int endX, int endY, CvSeq* result/*, int beforeIndex*/) {
	// dinh dung de quy, nhung ket qua sai bet ve thu tu cua cac pixel
	/*int betX = (startX + endX)/2;
	int betY = (startY + endY)/2;
	if ( ! ((betX == startX && betY == startY) || (betX == endX && betY == endY)) ) {
		CvPoint bet = cvPoint(betX, betY);
		//cvSeqPush(result, &bet);
		cvSeqInsert(result, beforeIndex, &bet);

		int newBeforeIndex = beforeIndex + 2*int(beforeIndex*1.0/2);
		findBetweenPos(startX, startY, betX, betY, result, newBeforeIndex);
		findBetweenPos(betX, betY, endX, endY, result, newBeforeIndex + 2);
	}*/

	// http://www.cs.unc.edu/~mcmillan/comp136/Lecture6/Lines.html

    int dy = endY - startY;
    int dx = endX - startX;
    int stepx, stepy;

    if (dy < 0) { dy = -dy;  stepy = -1; } else { stepy = 1; }
    if (dx < 0) { dx = -dx;  stepx = -1; } else { stepx = 1; }
    dy <<= 1; // dy is now 2*dy
    dx <<= 1; // dx is now 2*dx

    // raster.setPixel(pix, x0, y0);
	CvPoint bet = cvPoint(startX, startY);
	cvSeqPush(result, &bet);

    if (dx > dy) {
        int fraction = dy - (dx >> 1);  // same as 2*dy - dx
        while (startX != endX) {
            if (fraction >= 0) {
                startY += stepy;
                fraction -= dx;        // same as fraction -= 2*dx
            }
            startX += stepx;
            fraction += dy;            // same as fraction -= 2*dy
            
			//raster.setPixel(pix, startX, startY);
			CvPoint bet = cvPoint(startX, startY);
			cvSeqPush(result, &bet);
        }
    } else {
        int fraction = dx - (dy >> 1);
        while (startY != endY) {
            if (fraction >= 0) {
                startX += stepx;
                fraction -= dy;
            }
            startY += stepy;
            fraction += dx;
            
			//raster.setPixel(pix, startX, startY);
			CvPoint bet = cvPoint(startX, startY);
			cvSeqPush(result, &bet);
        }
    }
}

void nttStandardHoughEx::findFiniteLines(IplImage* img, CvSeq *linesFound, CvSeq* betPos, int maxGap, int lineLength) {
	int width = img->width;	int height = img->height;	int step = img->widthStep;
	uchar* imgData = (uchar *)img->imageData;

	bool found = false;
	int countToMaxGap = 0;
	CvPoint* startPoint = (CvPoint*)cvGetSeqElem(betPos, 0); // startPoint > 0 (co gia tri, chi ko phai la position)
	CvPoint* endPoint = startPoint;

	int i = 1;
	while (i < betPos->total) {
		CvPoint* pos = (CvPoint*)cvGetSeqElem(betPos, i);
		
		for (int j = -min(maxGap/2, pos->x); j <= min(maxGap/2, width - pos->x); j++) {
			int xChecked = pos->x + j; //
			int yChecked = pos->y; 	
			if (imgData[yChecked*step + xChecked] > 0) {
				found = true;
				countToMaxGap = 0;
				//endPoint->x = xChecked;
				//endPoint->y = yChecked;
				endPoint = pos;
				break;
			}
		}

		if (!found) countToMaxGap++; // dem so lan ko tim thay lien tuc

		if (countToMaxGap == maxGap || i == betPos->total - 1) {
			if (abs(startPoint->x - endPoint->x) > lineLength || abs(startPoint->y - endPoint->y) > lineLength) {
				CvRect line = { startPoint->x, startPoint->y, endPoint->x, endPoint->y }; 
				cvSeqPush(linesFound, &line);
			}

			i++;
			startPoint = (CvPoint*)cvGetSeqElem(betPos, i); // luc nay start point co the ko co gia tri, co the chi la position thoi
			endPoint = startPoint;

			countToMaxGap = 0;
		}

		i++;
		found = false;
	}
}