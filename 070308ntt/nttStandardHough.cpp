#include "StdAfx.h"
#include "nttStandardHough.h"
#include "cv.h"
#include "highgui.h"

nttStandardHough::nttStandardHough(void) { }
nttStandardHough::~nttStandardHough(void) { }

void nttStandardHough::run(IplImage* img, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax) {
	int height = img->height, width = img->width, step = img->widthStep;
	const uchar* imgData = (uchar *)img->imageData;

	int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);
	int numOfRho = cvRound( ((width + height)*2 + 1) / rhoResolution );
	int halfOfNumOfRho = cvRound(numOfRho / 2);

	// init satisfyListStd(accValue, thetaIndex, rhoIndex, 0 /* not used */ )
	CvMemStorage* satisfyListStorageStd = cvCreateMemStorage(0);
	CvSeq* satisfyListStd = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvRect), satisfyListStorageStd);

	// init accumulator
	CvMat* accumulatorStd = cvCreateMat(numOfTheta, numOfRho, CV_32SC1);
	cvZero(accumulatorStd);
	accColsStd = accumulatorStd->cols;
	accDataStd = accumulatorStd->data.i;
	
	// select pixels and accumulate to accumulator
	for(int i = 0; i < height; i++) 
		for(int j = 0; j < width; j++) 
			if (imgData[i*step + j] > 0) accumulate(j, i, rhoResolution, numOfTheta, halfOfNumOfRho);

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

			if (checkLocalMaximumWithWindow(3, 3, numOfTheta, numOfRho, i, j, threshold)) {
				int accValue = accDataStd[i*accColsStd + j];
				CvRect s = { accValue, i, j, 0 }; // 0: not used
				cvSeqPush(satisfyListStd, &s);
			}
		}

	// sap xep satisfyListStd theo thu tu giam dan cua accValue
	cvSeqSort( satisfyListStd, cmp_func, 0 /* userdata is not used here */ );

	// lay ra n = linesMax phan tu de dua va linesFound
	for (int i = 0; i < min(linesMax, satisfyListStd->total); i++) {
		CvRect* s = (CvRect*)cvGetSeqElem(satisfyListStd, i);
		int thetaIndex = s->y, rhoIndex = s->width;

		double cosTheta = cosThetaListStd[thetaIndex], sinTheta = sinThetaListStd[thetaIndex];
		double rho = (rhoIndex - halfOfNumOfRho) * rhoResolution;
		double x0 = cosTheta*rho, y0 = sinTheta*rho;
		CvRect lr = { cvRound(x0 + 1000*(-sinTheta)), cvRound(y0 + 1000*(cosTheta)), 
					  cvRound(x0 - 1000*(-sinTheta)), cvRound(y0 - 1000*(cosTheta)) }; 
		cvSeqPush(linesFound, &lr);
	}

	// release
	cvReleaseMemStorage(&satisfyListStorageStd);
	cvReleaseMat(&accumulatorStd);	
	delete sinThetaListStd;
	delete cosThetaListStd;
}

// ================= private =================

int nttStandardHough::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// 180 degree [ -90 ; +90 ] ~ CV_PI radian [ -PI/2 ; PI/2 ]
	int numOfTheta = cvRound(PI / thetaResolution) + 1; // cvRound: convert floating-point number to integer
	sinThetaListStd = new double[numOfTheta];
	cosThetaListStd = new double[numOfTheta];
	for (int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI;
		sinThetaListStd[i] = utilStd.roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaListStd[i] = utilStd.roundDouble(cos(theta), DOUBLE_PRECISION);
		//sinThetaListStd[i] = sin(theta);
		//cosThetaListStd[i] = cos(theta);
	}
	return numOfTheta;
}

void nttStandardHough::accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho) {
	for (int i = 0; i < numOfTheta; i++) {
		double rho = x*cosThetaListStd[i] + y*sinThetaListStd[i];
		int indexRho = cvRound( rho / rhoResolution ) + halfOfNumOfRho;
		double remainder = rho - (indexRho - halfOfNumOfRho) * rhoResolution;
		if (remainder > rhoResolution/2) indexRho++;
		accDataStd[i*accColsStd + indexRho] += 1;
	}
}

bool nttStandardHough::checkLocalMaximumWithWindow(int winHeight, int winWidth, int numOfTheta, int numOfRho, 
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

	return true;
}

