#include "StdAfx.h"
#include "nttPyramidHough.h"
#include "cv.h"
#include "highgui.h"
#include "nttUtil.h"

// ----- for computation -----
#define DOUBLE_PRECISION (7)
#define PI (3.1415927) // 3.1415926535897932384626433832795
#define HALF_PI (1.5707963) // 1.5707963267948966192313216916398
// ---------------------------

#define SUBIMAGE_SIZE (64)
#define NUM_CHILDREN (4)

double* sinThetaListPy;	
double* cosThetaListPy;

int *accDataPy;
int accColsPy;

nttUtil utilPy;

// -------

nttPyramidHough::nttPyramidHough(void) { }
nttPyramidHough::~nttPyramidHough(void) { }

void nttPyramidHough::run(IplImage* img, CvSeq *linesFound, double rhoResolution, double thetaResolution, int threshold, int linesMax) {
	int height = img->height, width = img->width, step = img->widthStep;
	uchar* imgData = (uchar *)img->imageData;

	int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);
	int numOfRho = cvRound( ((width + height)*2 + 1) / rhoResolution );
	int halfOfNumOfRho = cvRound(numOfRho / 2);

	// init accumulator
	CvMat* accumulatorPy = cvCreateMat(numOfTheta, numOfRho, CV_32SC1);
	accColsPy = accumulatorPy->cols;
	accDataPy = accumulatorPy->data.i;

	// pyramid divide
	int partWidth = width, numOfDivide = 0;
	while (partWidth > SUBIMAGE_SIZE) { 
		partWidth = partWidth / 2;
		numOfDivide++;
	}

	int numOfPart = pow(2., (double)numOfDivide);
	int partHeight = height / numOfPart;

	// init thetaAndRhoList
	int numOfPartSqr = numOfPart * numOfPart;
	CvMemStorage** thetaAndRhoListStorages = new CvMemStorage*[numOfPartSqr];
	CvSeq** thetaAndRhoLists = new CvSeq*[numOfPartSqr];
	for (int i = 0; i < numOfPartSqr; i++) {
		thetaAndRhoListStorages[i] = cvCreateMemStorage(0);
		thetaAndRhoLists[i] = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint2D64f), thetaAndRhoListStorages[i]);
	}

	// tim ra line(theta, rho) trong moi part
	int partCount = 0;
	for (int i = 0; i < (numOfPart == 1 ? 1 : numOfPart / 2); i++) { // height
		for (int j = 0; j < (numOfPart == 1 ? 1 : numOfPart / 2); j++) { // width

			for (int k = 0; k <= 1; k++) {
				for (int l = 0; l <= 1; l++) {
					int h = (2*i + k) * partHeight;
					int w = (2*j + l) * partWidth;

					cvZero(accumulatorPy);
					for (int m = h; m < min(h + partHeight, height); m++) { // accumulate from the subimage
						for (int n = w; n < min(w + partWidth, width); n++) {
							if (imgData[m*step + n] > 0) 
								accumulate(n, m, rhoResolution, numOfTheta, halfOfNumOfRho);
						}
					}

					// lay ra cac gia tri nao tich luy trong accumulator ma > threshold 
					for(int m = 0; m < numOfTheta; m++) {
						for(int n = 0; n < numOfRho; n++) {
							int deltaI = 1, deltaJ = 1;
							if (m == 0 || m == numOfTheta - 1) deltaI = 0;
							if (n == 0 || n == numOfRho - 1) deltaJ = 0;

							int accValue = accDataPy[m*accColsPy + n];			
							int accValueUp = accDataPy[(m - deltaI)*accColsPy + n];			
							int accValueLeft = accDataPy[m*accColsPy + n - deltaJ];
							int accValueRight = accDataPy[m*accColsPy + n + deltaJ];			
							int accValueDown = accDataPy[(m + deltaI)*accColsPy + n];

							if ( (accValue > threshold) && (accValue > accValueUp) && (accValue >= accValueDown)
										&& (accValue > accValueLeft) && (accValue >= accValueRight)) { // find local maximum
								double theta = m * thetaResolution - HALF_PI;
								double rho = (n - halfOfNumOfRho) * rhoResolution;
								CvPoint2D64f thetaAndRho = { theta, rho }; 
								cvSeqPush(thetaAndRhoLists[partCount], &thetaAndRho);
							}
						}
					}
					if (numOfPart == 1) break;
					partCount++;
				}
				if (numOfPart == 1) break;
			}
			if (numOfPart == 1) break;
		}
		if (numOfPart == 1) break;
	}

	// ghep cac line nay lai
	CvMemStorage* tmpStorage = cvCreateMemStorage(0);
	CvSeq* tmp = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint2D64f), tmpStorage);

	int remainNumOfPart = numOfPart;
	while (numOfPart > 1 && remainNumOfPart > 0) {

		for (int i = 0; i < pow(numOfPart, 2.) / NUM_CHILDREN; i++) {				
			tmp = cvCloneSeq(thetaAndRhoLists[i * NUM_CHILDREN], tmpStorage);
			int tmpTotal = tmp->total;

			for (int j = 1; j < NUM_CHILDREN; j++) {
				
				// merging
				int thetaAndRhoListIndex = i * NUM_CHILDREN + j;				
				for (int k = 0; k < tmpTotal; k++) {
					CvPoint2D64f* line1 = (CvPoint2D64f*)cvGetSeqElem(tmp, k);
					double theta1 = line1->x, rho1 = line1->y;

					bool canMergeFlag = false;
					for (int l = 0; l < thetaAndRhoLists[thetaAndRhoListIndex]->total; l++) {
						CvPoint2D64f* line2 = (CvPoint2D64f*)cvGetSeqElem(thetaAndRhoLists[thetaAndRhoListIndex], l);
						double theta2 = line2->x, rho2 = line2->y;

						if (fabs(theta1 - theta2) < thetaResolution && fabs(rho1 - rho2) < rhoResolution) {
							double newTheta = utilPy.roundDouble( (theta1 + theta2) / 2, DOUBLE_PRECISION );

							double c = utilPy.roundDouble( cos( (theta1 - theta2) / 2 ) , DOUBLE_PRECISION );
							double newRho = utilPy.roundDouble( (rho1 + rho2) / (2 * c), DOUBLE_PRECISION );

							CvPoint2D64f thetaAndRho = { newTheta, newRho };
							cvSeqPush(tmp, &thetaAndRho); 

							canMergeFlag = true;
						} else {
							cvSeqPush(tmp, line2); 
						}
					}

					if (canMergeFlag) cvSeqRemove(tmp, k);
				}				
			}

			thetaAndRhoLists[i] = cvCloneSeq(tmp, thetaAndRhoListStorages[i]);
			cvClearSeq(tmp);
		}

		remainNumOfPart /= 4;
	}

	// ve line
	for (int i = 0; i < min(linesMax, thetaAndRhoLists[0]->total); i++) {
		CvPoint2D64f* thetaAndRho = (CvPoint2D64f*)cvGetSeqElem(thetaAndRhoLists[0], i);
		double theta = thetaAndRho->x, rho = thetaAndRho->y;

		double cosTheta = utilPy.roundDouble(cos(theta), DOUBLE_PRECISION);
		double sinTheta = utilPy.roundDouble(sin(theta), DOUBLE_PRECISION);
		double x0 = cosTheta*rho, y0 = sinTheta*rho;
		CvRect lr = { cvRound(x0 + 1000*(-sinTheta)), cvRound(y0 + 1000*(cosTheta)), 
					  cvRound(x0 - 1000*(-sinTheta)), cvRound(y0 - 1000*(cosTheta)) }; 
		cvSeqPush(linesFound, &lr);
	}

	// release
	delete thetaAndRhoListStorages;
	cvReleaseMemStorage(&tmpStorage);
	cvReleaseMat(&accumulatorPy);		
	delete sinThetaListPy;
	delete cosThetaListPy;
}

// ================= private =================

int nttPyramidHough::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// 180 degree [ -90 ; +90 ] ~ CV_PI radian [ -PI/2 ; PI/2 ]
	int numOfTheta = cvRound(PI / thetaResolution); // cvRound: convert floating-point number to integer
	sinThetaListPy = new double[numOfTheta];
	cosThetaListPy = new double[numOfTheta];
	for (int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI;
		sinThetaListPy[i] = utilPy.roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaListPy[i] = utilPy.roundDouble(cos(theta), DOUBLE_PRECISION);
	}
	return numOfTheta;
}

void nttPyramidHough::accumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho) {
	for (int i = 0; i < numOfTheta; i++) {
		double rho = x*cosThetaListPy[i] + y*sinThetaListPy[i];
		int indexRho = cvRound( rho / rhoResolution ) + halfOfNumOfRho;
		double remainder = rho - (indexRho - halfOfNumOfRho) * rhoResolution;
		if (remainder > rhoResolution/2) indexRho++;
		accDataPy[i*accColsPy + indexRho] += 1;
	}
}