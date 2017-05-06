#include "StdAfx.h"
#include "nttStandardHoughNew.h"

nttStandardHoughNew::nttStandardHoughNew(void) { }
nttStandardHoughNew::~nttStandardHoughNew(void) { }

void nttStandardHoughNew::run(IplImage* img, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax) {
	int height = img->height, width = img->width, step = img->widthStep;
	const uchar* imgData = (uchar *)img->imageData;

	int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);
	int numOfRho = cvRound( ((width + height)*2 + 1) / rhoResolution );
	int halfOfNumOfRho = cvRound(numOfRho / 2);

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
			
			double accValueUpLeft = accData[(i - deltaI)*accCols + j - deltaJ];
			double accValueUpRight = accData[(i - deltaI)*accCols + j + deltaJ];
			double accValueDownLeft = accData[(i + deltaI)*accCols + j - deltaJ];			
			double accValueDownRight = accData[(i + deltaI)*accCols + j + deltaJ];

			//if ( (accValue > threshold) && (accValue > accValueUp) && (accValue >= accValueDown)
			//							&& (accValue > accValueLeft) && (accValue >= accValueRight)) { // find local maximum
			if ( (accValue > threshold) 
					&& (accValue > accValueUpLeft) && (accValue > accValueUp) && (accValue > accValueUpRight) && (accValue > accValueLeft)
					&& (accValue >= accValueRight) && (accValue >= accValueDownLeft) && (accValue >= accValueDown) && (accValue >= accValueDownRight) ) { // find local maximum
				CvPoint3D64f s = { accValue, i, j }; 
				cvSeqPush(satisfyList, &s);
			}
		}
		
	// sort satisfyList in the descending order of accValue
	cvSeqSort( satisfyList, cmp_func, 0 /* userdata is not used here */ );

	// take in the satisfyList a number of line(theta, rho) to return
	for (int i = 0; i < min(linesMax, satisfyList->total); i++) {
		CvPoint3D64f* s = (CvPoint3D64f*)cvGetSeqElem(satisfyList, i);
		int thetaIndex = (int)(s->y), rhoIndex = (int)(s->z);

		double cosTheta = cosThetaList[thetaIndex], sinTheta = sinThetaList[thetaIndex];
		double rho = (rhoIndex - halfOfNumOfRho) * rhoResolution;
		double x0 = cosTheta*rho, y0 = sinTheta*rho;
		CvRect lr = { cvRound(x0 + 1000*(-sinTheta)), cvRound(y0 + 1000*(cosTheta)), 
					  cvRound(x0 - 1000*(-sinTheta)), cvRound(y0 - 1000*(cosTheta)) }; 
		cvSeqPush(linesFound, &lr);
	}

	// release
	cvReleaseMemStorage(&satisfyListStorage);
	cvReleaseMat(&accumulator);		
	delete sinThetaList;
	delete cosThetaList;
}

// ================= private =================

int nttStandardHoughNew::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// 180 degree [ -90 ; +90 ] ~ CV_PI radian [ -PI/2 ; PI/2 ]
	int numOfTheta = cvRound(PI / thetaResolution) + 1; // cvRound: convert floating-point number to integer
	sinThetaList = new double[numOfTheta];
	cosThetaList = new double[numOfTheta];
	for (int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI;
		sinThetaList[i] = util.roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaList[i] = util.roundDouble(cos(theta), DOUBLE_PRECISION);
	}
	return numOfTheta;
}

void nttStandardHoughNew::accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho) {
	for (int i = 0; i < numOfTheta; i++) {
		double rhoExact = x*cosThetaList[i] + y*sinThetaList[i];
		int indexRhoLow = cvRound( rhoExact / rhoResolution ) + halfOfNumOfRho;
		double rhoLow = (indexRhoLow - halfOfNumOfRho) * rhoResolution;
		accData[i*accCols + indexRhoLow] += fabs((rhoLow + rhoResolution) - rhoExact); // rhoHigh - rhoExact
		accData[i*accCols + (indexRhoLow + 1)] += fabs(rhoExact - rhoLow);
	}
}
