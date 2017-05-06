#include "StdAfx.h"
#include "nttHoughGreen.h"
#include "cv.h"
#include "highgui.h"
#include "nttUtil.h"

// ----- for computation -----
#define DOUBLE_PRECISION (7)
#define PI (3.1415927) // 3.1415926535897932384626433832795
#define HALF_PI (1.5707963) // 1.5707963267948966192313216916398
// ---------------------------

double* sinThetaListHG;
double* cosThetaListHG;

CvMemStorage* pixelsOnContourStorageHG;			CvSeq* pixelsOnContourHG;
CvMemStorage* pixelsOnContourRotatedStorageHG;	CvSeq* pixelsOnContourRotatedHG;

nttUtil utilHG;

// -----------

nttHoughGreen::nttHoughGreen(void) { }
nttHoughGreen::~nttHoughGreen(void) { }

/********************************************************
 * Trace obj's R contour C, store all (x, y) or (r, phi)
 *	for all theta
 *		for all L of C {
 *			rotation by theta
 *			resample rhoL
 *			update HGT(theta, rhoL)
 *		}
 ********************************************************/
void nttHoughGreen::run(IplImage* src, CvSeq *linesFound, double thetaResolution, double threshold) {
    int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);

	CvMemStorage* contourStorage = cvCreateMemStorage(0);
    CvSeq* contour = 0;
    cvFindContours( src, contourStorage, &contour, sizeof(CvContour), CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE );	
	
    for( ; contour != 0; contour = contour->h_next ) {
		double hgtMax = 0, rhoMax = 0; 
		int thetaIndex = -1;
		getPixelsOnContour(contour);

		for (int i = 0; i < numOfTheta; i++) {
			double hgt = 0, rho = 0;
			rotatePixelsOnContour(i);

			for (int j = 0; j < pixelsOnContourHG->total; j++) {
				// rotation by theta
				CvPoint2D32f* p = (CvPoint2D32f*)cvGetSeqElem(pixelsOnContourRotatedHG, j);
				double xRotated = p->x, yRotated = p->y;

				p = (CvPoint2D32f*)cvGetSeqElem(pixelsOnContourRotatedHG, j == 0 ? pixelsOnContourRotatedHG->total - 1: j - 1);
				double deltaRhoRotated = xRotated - p->x;
				
				// resample rhoL
				if ( !resampleRho(deltaRhoRotated, i) ) continue;

				// chieu (xRotated, yRotated) xuong truc X
				int length = 0;
				for (int k = 0; k < pixelsOnContourHG->total; k++) {
					if (k == j) continue;

					p = (CvPoint2D32f*)cvGetSeqElem(pixelsOnContourRotatedHG, k);
					double xRotatedProj = p->x, yRotatedProj = p->y;

					p = (CvPoint2D32f*)cvGetSeqElem(pixelsOnContourRotatedHG, k == 0 ? pixelsOnContourRotatedHG->total - 1: k - 1);
					double deltaRhoRotatedProj = xRotatedProj - p->x;
					
					// resample rhoL
					/*if ( resampleRho(deltaRhoRotatedProj, i) ) {
						//double length = 0; 
						//if ( (yRotated > 0 && yRotatedProj > 0) || (yRotated < 0 && yRotatedProj < 0) ) length = fabs(yRotated - yRotatedProj);
						//else length = fabs(yRotated) + fabs(yRotatedProj);

						//double length = yRotated * deltaRhoRotated + yRotatedProj * deltaRhoRotatedProj; 

						if (hgt < length) hgt = length; // update HGT
					}*/

					if ( resampleRho(deltaRhoRotatedProj, i) )
						if (fabs(xRotated - xRotatedProj) <= 0.5) 
							if (fabs(yRotatedProj) > fabs(yRotated)) length++;
				}

				//if (hgt < length) hgt = length; // update HGT

				//rho = xRotated;

				if (length >= threshold) {
					double cosTheta = cosThetaListHG[i], sinTheta = sinThetaListHG[i];
					double x0 = cosTheta*xRotated, y0 = sinTheta*xRotated;
					int startX = cvRound(x0 + 1000*(-sinTheta)), startY = cvRound(y0 + 1000*(cosTheta));
					int endX = cvRound(x0 - 1000*(-sinTheta)), endY = cvRound(y0 - 1000*(cosTheta));
					
					if (!isContain(linesFound, startX, startY, endX, endY))	{
						CvRect lr = { startX, startY, endX, endY }; 
						cvSeqPush(linesFound, &lr);
					}
				}
			}

			/*if (hgtMax < hgt) {
				hgtMax = hgt;
				thetaIndex = i;
				rhoMax = rho;
			}*/
		}

		/*if (hgtMax >= threshold) {
			double cosTheta = cosThetaListHG[thetaIndex], sinTheta = sinThetaListHG[thetaIndex];
			double x0 = cosTheta*rhoMax, y0 = sinTheta*rhoMax;
			CvRect lr = { cvRound(x0 + 1000*(-sinTheta)), cvRound(y0 + 1000*(cosTheta)), 
						  cvRound(x0 - 1000*(-sinTheta)), cvRound(y0 - 1000*(cosTheta)) }; 
			cvSeqPush(linesFound, &lr);
		}*/
    }	
	cvReleaseMemStorage(&contourStorage);
	cvReleaseMemStorage(&pixelsOnContourStorageHG);
	cvReleaseMemStorage(&pixelsOnContourRotatedStorageHG);
	delete sinThetaListHG;
	delete cosThetaListHG;	
}

// ================= private =================

int nttHoughGreen::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// CV_PI radian [ 0 ; PI )
	int numOfTheta = cvRound(PI / thetaResolution); 
	sinThetaListHG = new double[numOfTheta];
	cosThetaListHG = new double[numOfTheta];
	for (int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution;
		sinThetaListHG[i] = utilHG.roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaListHG[i] = utilHG.roundDouble(cos(theta), DOUBLE_PRECISION);
	}
	return numOfTheta;
}

void nttHoughGreen::getPixelsOnContour(CvSeq* contour) {
	pixelsOnContourStorageHG = cvCreateMemStorage(0);
	pixelsOnContourHG = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), pixelsOnContourStorageHG);

	CvSeqReader reader;
	cvStartReadSeq( contour, &reader, 0 );
	CvPoint p;
	CV_READ_SEQ_ELEM( p, reader );

	int x = p.x, y = p.y, startX = p.x, startY = p.y; // (x, y)s are just SAMPLE points	
	while(true) {
		cvSeqPush(pixelsOnContourHG, &p); //
		CV_READ_SEQ_ELEM( p, reader );
		if (p.x == startX && p.y == startY) break;

		int rangeX = x - p.x, rangeY = y - p.y;
		int coeffX = 0, coeffY = 0;
		if (rangeX > 0) coeffX = -1; else if (rangeX < 0) coeffX = 1;
		if (rangeY > 0) coeffY = -1; else if (rangeY < 0) coeffY = 1;

		for (int i = 1; i < max(abs(rangeX), abs(rangeY)); i++) {
			CvPoint ep = { x + i*coeffX, y + i*coeffY };
			cvSeqPush(pixelsOnContourHG, &ep); //
		}
		x = p.x;
		y = p.y;			
	}
}

void nttHoughGreen::rotatePixelsOnContour(int thetaIndex) {	
	pixelsOnContourRotatedStorageHG = cvCreateMemStorage(0);
	pixelsOnContourRotatedHG = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint2D32f), pixelsOnContourRotatedStorageHG);
	
	for (int j = 0; j < pixelsOnContourHG->total; j++) {
		CvPoint* p = (CvPoint*)cvGetSeqElem(pixelsOnContourHG, j);
		int x = p->x, y = p->y;

		CvPoint2D32f pf = { (float)(x*cosThetaListHG[thetaIndex] + y*sinThetaListHG[thetaIndex]), // xRotated = rhoRotated
							(float)(-x*sinThetaListHG[thetaIndex] + y*cosThetaListHG[thetaIndex]) };
		cvSeqPush(pixelsOnContourRotatedHG, &pf); 
	}
}

bool nttHoughGreen::resampleRho(double deltaRho, int thetaIndex) {
	double c1 = (cosThetaListHG[thetaIndex] > 0 ? -cosThetaListHG[thetaIndex] : cosThetaListHG[thetaIndex]);
	double s1 = (sinThetaListHG[thetaIndex] > 0 ? -sinThetaListHG[thetaIndex] : sinThetaListHG[thetaIndex]);

	double c2 = (cosThetaListHG[thetaIndex] < 0 ? -cosThetaListHG[thetaIndex] : cosThetaListHG[thetaIndex]);
	double s2 = (sinThetaListHG[thetaIndex] < 0 ? -sinThetaListHG[thetaIndex] : sinThetaListHG[thetaIndex]);

	return min(c1, s1) <= deltaRho && deltaRho <= max(c2, s2);

	//return min(cosThetaListHG[thetaIndex], sinThetaListHG[thetaIndex]) <= deltaRho 
	//	 && deltaRho <= max(cosThetaListHG[thetaIndex], sinThetaListHG[thetaIndex]);
}

bool nttHoughGreen::isContain(CvSeq *linesFound, int startX, int startY, int endX, int endY) {
	for (int i = 0; i < linesFound->total; i++) {
		CvRect* ln = (CvRect*)cvGetSeqElem(linesFound, i);
		if (ln->x == startX && ln->y == startY && ln->width == endX && ln->height == endY) return true;
	}
	return false;	
}