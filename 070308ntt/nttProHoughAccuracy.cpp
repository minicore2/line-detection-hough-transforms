#include "StdAfx.h"
#include "nttProHoughAccuracy.h"
#include "cv.h"
#include "highgui.h"
#include "nttUtil.h"

// ----- for computation -----
#define DOUBLE_PRECISION (7)
#define PI (3.1415927) // 3.1415926535897932384626433832795
#define HALF_PI (1.5707963) // 1.5707963267948966192313216916398
#define SQRT_2PI (2.5066283) // 2.5066282931459941533857660217468
// ---------------------------

#define NOISE (1)
#define APPROXIMATE_ACCVALUE (0.00001)

double* sinThetaListPa;
double* cosThetaListPa;

CvMemStorage* storagePa;				CvSeq* onPixelsPa;
CvMemStorage* removedPixelsStoragePa;	CvSeq* removedPixelsPa; // can be disabled if not use "spreading"

double *accDataPa;
int accColsPa;

uchar *copiedDataPa;

CvPoint line1Pa[2] = {{0,0}, {0,0}};
CvPoint line2Pa[2] = {{0,0}, {0,0}};

nttUtil utilPa;
// 
CvMemStorage* eightHStoragePa;			CvSeq* eightHPa;
CvMemStorage* eightRhoStoragePa;		CvSeq* eightRhoPa;
CvMemStorage* listOfHmaxStoragePa;		CvSeq* listOfHmaxPa;
CvMemStorage* listOfRhoHmaxStoragePa;	CvSeq* listOfRhoHmaxPa;
CvMemStorage* listOfThetaHmaxStoragePa;	CvSeq* listOfThetaHmaxPa;

// ------------

nttProHoughAccuracy::nttProHoughAccuracy(void) { }
nttProHoughAccuracy::~nttProHoughAccuracy(void) { }

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
void nttProHoughAccuracy::run(IplImage* src, CvSeq *linesFound, double rhoResolution, double thetaResolution, 
								double threshold, int minLineLength, int maxGap) {
	IplImage* copiedImagePa = cvCloneImage(src);
	int height = copiedImagePa->height, width = copiedImagePa->width, step = copiedImagePa->widthStep;
	copiedDataPa = (uchar *)copiedImagePa->imageData;
	
	int numOfTheta = generateSinThetaAndCosThetaList(thetaResolution);
	int numOfRho = cvRound( ((width + height)*2 + 1) / rhoResolution );
	int halfOfNumOfRho = cvRound(numOfRho / 2);

	init();
	
	// init accumulator
	CvMat* accumulatorPa = cvCreateMat(numOfTheta, numOfRho, CV_32FC2);
	cvZero(accumulatorPa);
	accColsPa = accumulatorPa->cols;
	accDataPa = accumulatorPa->data.db;

	getAllOnPixels(height, width, step); // collect ON pixels (in src image) to onPixelsPa

	// check input image, if empty => finish
	while (onPixelsPa->total > 0) {		
		int randomNum = utilPa.generateRandomNumberOpenCV(onPixelsPa->total);
		CvPoint* point = (CvPoint*)cvGetSeqElem(onPixelsPa, randomNum); // select a random ON pixel from onPixelsPa
		int x = point->x, y = point->y;

		accumulateOrDeaccumulate(x, y, rhoResolution, numOfTheta, halfOfNumOfRho, 1); // update accumulator with the selected random pixel
		cvSeqPush(removedPixelsPa, point); // can be disabled if not use "spreading"
		cvSeqRemove(onPixelsPa, randomNum); // remove selected pixel from image

		//spreadOrUnspread(numOfTheta, numOfRho, halfOfNumOfRho, minLineLength, thetaResolution, rhoResolution, 1); // can be disabled if not use "spreading"

		// find the highest peak in accumulator	
		double highestPeak = 0;
		for(int i = 0; i < numOfTheta; i++) 
			for(int j = 0; j < numOfRho; j++) {
				double accValue = accDataPa[i*accColsPa + j];
				if (highestPeak < accValue) highestPeak = accValue; 
			}

		// fit a Gaussian to 8-points surrounding the max in the column (in rho direction)
		double alpha, beta, gamma, delta, A;
		cvClearSeq(listOfHmaxPa);
		cvClearSeq(listOfRhoHmaxPa);
		cvClearSeq(listOfThetaHmaxPa);		
		for(int i = 0; i < numOfTheta; i++) 
			for(int j = 0; j < numOfRho; j++) {
				double accValue = accDataPa[i*accColsPa + j];
				if (accValue > 0 && fabs(highestPeak - accValue) <= APPROXIMATE_ACCVALUE) // accValue is also the highest peak
					for (int m = -min(1, i); m <= min(1, numOfTheta - i - 1); m++) { // m in [-1, 1]						
						for (int mr = -min(4, j); mr <= min(4, numOfRho - j - 1); mr++) // mr in [-4, 0) and (0, 4]
							if (mr != 0) { // tru "diem" co highest peak ra, chi xe 4 diem tren & 4 diem duoi cua highest peak
								double h = accDataPa[(i+m)*accColsPa + (j+mr)];
								if (h != 0) {
									double rho = ((j+mr) - halfOfNumOfRho) * rhoResolution;
									cvSeqPush(eightRhoPa, &rho);
									cvSeqPush(eightHPa, &h);
								}
							}
						alpha = calculateAlpha(eightHPa, eightRhoPa); // c
						beta = calculateBeta(eightHPa, eightRhoPa); // b
						gamma = calculateGamma(eightRhoPa); // a													
						delta = beta*beta - 4*gamma*(-alpha); // delta = b^2 - 4ac, in the equation: gamma*A^2 + beta*A - alpha = 0	
						A = (-beta + sqrt(delta)) / (2*gamma); // A is the constant we need to find out
						// A = (-beta - sqrt(delta)) / (2*gamma);

						double hMax = utilPa.roundDouble(A, DOUBLE_PRECISION); // hMax = A * e*(..rhoMax-rhoMean..); we have rhoMax = rhoMean
						double rhoHmax = utilPa.mean(eightRhoPa); // rhoHmax = rho_highestPeak						
						double thetaHmax = (i+m) * thetaResolution - HALF_PI;
						cvSeqPush(listOfHmaxPa, &hMax);
						cvSeqPush(listOfRhoHmaxPa, &rhoHmax);
						cvSeqPush(listOfThetaHmaxPa, &thetaHmax);
						cvClearSeq(eightHPa);
						cvClearSeq(eightRhoPa);
					}
			}

		// fit a Gaussian is theta direction, using the above hMax(s) and rhoMax(s)
		alpha = calculateAlpha(listOfHmaxPa, listOfThetaHmaxPa); // c
		beta = calculateBeta(listOfHmaxPa, listOfThetaHmaxPa); // b
		gamma = calculateGamma(listOfThetaHmaxPa); // a
		delta = beta*beta - 4*gamma*(-alpha);
		A = (-beta + sqrt(delta)) / (2*gamma);
		double hPeak = utilPa.roundDouble(A, DOUBLE_PRECISION);

		// check vote threshold
		bool isNotANumber = (hPeak != hPeak);
		if (isNotANumber || hPeak < threshold) continue; 

		// linear interpolation to find rhoPeak
		double thetaPeak = utilPa.mean(listOfThetaHmaxPa);

		int indexLow = 0;
		for (int i = 0; i < listOfThetaHmaxPa->total; i++) {
			double* t = (double*)cvGetSeqElem(listOfThetaHmaxPa, i);
			if ( fabs(thetaPeak - (*t)) <= thetaResolution ) {
				indexLow = i;	
				break;
			}
		}
		double* thetaLow = (double*)cvGetSeqElem(listOfThetaHmaxPa, indexLow);
		double* thetaHigh = (double*)cvGetSeqElem(listOfThetaHmaxPa, indexLow + 1);
		double* rhoLow = (double*)cvGetSeqElem(listOfRhoHmaxPa, indexLow);
		double* rhoHigh = (double*)cvGetSeqElem(listOfRhoHmaxPa, indexLow + 1);

		double rhoPeak = ( (thetaPeak - (*thetaLow)) / ((*thetaHigh) - (*thetaLow)) ) * ((*rhoHigh) - (*rhoLow)) + (*rhoLow);
		
		// find x,y satisfying the condition rhoPeak = x*sin(thetaPeak) + y*cos(thetaPeak)
		double sinThetaPeak = utilPa.roundDouble(sin(thetaPeak), DOUBLE_PRECISION);
		double cosThetaPeak = utilPa.roundDouble(cos(thetaPeak), DOUBLE_PRECISION);

		for (int i = 0; i < removedPixelsPa->total; i++) {
			CvPoint* p = (CvPoint*)cvGetSeqElem(removedPixelsPa, i);
			x = p->x;
			y = p->y;
			if (fabs(rhoPeak - (x*cosThetaPeak + y*sinThetaPeak)) <= rhoResolution) break;
		}

		// look along a corridor specified by the peak in the accumulator
		//    (from the current point, walk in each direction along the found line and extract the line segment)
		findPixelsLeftAndRight(x, y, sinThetaPeak, cosThetaPeak, rhoPeak, rhoResolution, thetaResolution, minLineLength, 
			numOfTheta, numOfRho, halfOfNumOfRho, step, width, height, maxGap, 1);
		findPixelsLeftAndRight(x, y, sinThetaPeak, cosThetaPeak, rhoPeak, rhoResolution, thetaResolution, minLineLength, 
			numOfTheta, numOfRho, halfOfNumOfRho, step, width, height, maxGap, -1);

		// if line segment > min length, add it into output list (connect 2 lines,which run in 2 directions, together)
		int startX = line1Pa[1].x;	int startY = line1Pa[1].y;
		int endX = line2Pa[1].x;	int endY = line2Pa[1].y;
		if (abs(startX - endX) >= minLineLength || abs(startY - endY) >= minLineLength) {
			CvRect lr = { startX, startY, endX, endY }; // <-- CvRect = line
			cvSeqPush(linesFound, &lr);
		}
	}
	cvReleaseMat(&accumulatorPa);
	cvReleaseImage(&copiedImagePa);	
	cvReleaseMemStorage(&storagePa);
	cvReleaseMemStorage(&removedPixelsStoragePa);
	cvReleaseMemStorage(&listOfHmaxStoragePa);
	cvReleaseMemStorage(&listOfRhoHmaxStoragePa);
	cvReleaseMemStorage(&listOfThetaHmaxStoragePa);
	cvReleaseMemStorage(&eightHStoragePa);
	cvReleaseMemStorage(&eightRhoStoragePa);
	delete sinThetaListPa;
	delete cosThetaListPa;
}

// ================= private =================

int nttProHoughAccuracy::generateSinThetaAndCosThetaList(double thetaResolution) { 
	// 180 degree [ -90 ; +90 ] ~ CV_PI radian [ -PI/2 ; PI/2 ]
	int numOfTheta = cvRound(PI / thetaResolution); // cvRound: convert floating-point number to integer
	sinThetaListPa = new double[numOfTheta];
	cosThetaListPa = new double[numOfTheta];
	for (int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI;
		sinThetaListPa[i] = utilPa.roundDouble(sin(theta), DOUBLE_PRECISION);
		cosThetaListPa[i] = utilPa.roundDouble(cos(theta), DOUBLE_PRECISION);
	}
	return numOfTheta;
}

void nttProHoughAccuracy::getAllOnPixels(int height, int width, int step) {
	storagePa = cvCreateMemStorage(0);
	onPixelsPa = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), storagePa);
	for(int i = 0; i < height; i++) 
		for(int j = 0; j < width; j++) { 
			int position = i*step + j;
			if (copiedDataPa[position] > 0) {
				CvPoint p = {j, i}; // (x, y)
				cvSeqPush(onPixelsPa, &p);
			}
		}
}

void nttProHoughAccuracy::accumulateOrDeaccumulate(int x, int y, double rhoResolution, int numOfTheta, int halfOfNumOfRho, int accOrDeacc) {
	for (int i = 0; i < numOfTheta; i++) {
		double rhoExact = x*cosThetaListPa[i] + y*sinThetaListPa[i];
		int indexRhoLow = cvRound( rhoExact / rhoResolution ) + halfOfNumOfRho; 
		double rhoLow = (indexRhoLow - halfOfNumOfRho) * rhoResolution;
		accDataPa[i*accColsPa + indexRhoLow] += accOrDeacc * fabs((rhoLow + rhoResolution) - rhoExact); // rhoHigh - rhoExact
		accDataPa[i*accColsPa + (indexRhoLow + 1)] += accOrDeacc * fabs(rhoExact - rhoLow);
	}
}

void nttProHoughAccuracy::removeOnPixel(int x, int y) {
	for (int i = 0; i < onPixelsPa->total; i++) {
		CvPoint* point = (CvPoint*)cvGetSeqElem(onPixelsPa, i);
		if (point->x == x && point->y == y) {
			cvSeqRemove(onPixelsPa, i);
			break;
		}
	}
}

bool nttProHoughAccuracy::isOnPixelsContain(int x, int y) {
	for (int i = 0; i < onPixelsPa->total; i++) {
		CvPoint* point = (CvPoint*)cvGetSeqElem(onPixelsPa, i);
		if (point->x == x && point->y == y) return true;
	}
	return false;	
}

void nttProHoughAccuracy::findPixelsLeftAndRight(int x, int y, double sinThetaPeak, double cosThetaPeak, double rho, double rhoResolution, 
												double thetaResolution, int minLineLength, int numOfTheta, int numOfRho, int halfOfNumOfRho, 
												int step, int width, int height, int maxGap, int leftOrRight) {
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
				if (copiedDataPa[nextPosition] > 0 && fabs(rho - (nextX*cosThetaPeak + nextY*sinThetaPeak)) <= rhoResolution) {
					// "Unvote" from accumulator all pixels from the line that previously voted
					if (!isOnPixelsContain(nextX, nextY)) {
						//spreadOrUnspread(numOfTheta, numOfRho, halfOfNumOfRho, minLineLength, thetaResolution, rhoResolution, 1);
						accumulateOrDeaccumulate(x, y, rhoResolution, numOfTheta, halfOfNumOfRho, -1);
					}

					// remove pixels in segment from image
					CvPoint p = {nextX, nextY}; 
					cvSeqPush(removedPixelsPa, &p); // can be disabled if not use "spreading"
					removeOnPixel(nextX, nextY); 
					copiedDataPa[nextPosition] = 0;

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
		line1Pa[0].x = x;			line1Pa[0].y = y;
		line1Pa[1].x = currentX;	line1Pa[1].y = currentY;
	} else {
		line2Pa[0].x = x;			line2Pa[0].y = y;
		line2Pa[1].x = currentX;	line2Pa[1].y = currentY;
	}
}

// ================= methods for accuracy =================

void nttProHoughAccuracy::spreadOrUnspread(int numOfTheta, int numOfRho, int halfOfNumOfRho, int minLineLength,  
										   double thetaResolution, double rhoResolution, int sprOrUnspr) {
	for(int i = 0; i < numOfTheta; i++) {
		double theta = i * thetaResolution - HALF_PI;
		double standardDeviation = calculateStandardDeviation(minLineLength, thetaResolution, theta);
		for(int j = 0; j < numOfRho; j++) {	
			double rho = (j - halfOfNumOfRho) * rhoResolution;
			double convolutionValue = calculateConvolutionValue(standardDeviation, rho);
			accDataPa[i*accColsPa + j] = utilPa.roundDouble(accDataPa[i*accColsPa + j]*pow(convolutionValue, sprOrUnspr), DOUBLE_PRECISION);
		}
	}
}

double nttProHoughAccuracy::calculateStandardDeviation(int L, double thetaRes, double theta) {
	double thetaResDiv2 = utilPa.roundDouble(thetaRes / 2, DOUBLE_PRECISION);
	double sinThetaResDiv2 = utilPa.roundDouble(sin(thetaResDiv2), DOUBLE_PRECISION);
	double cosThetaResDiv2 = utilPa.roundDouble(cos(thetaResDiv2), DOUBLE_PRECISION);
	double sinTheta = utilPa.roundDouble(sin(theta), DOUBLE_PRECISION);
	return utilPa.roundDouble( (L*sinThetaResDiv2 + 6*NOISE*cosThetaResDiv2) * sinTheta / 2 , DOUBLE_PRECISION);
}

double nttProHoughAccuracy::calculateConvolutionValue(double standardDeviation, double rhoI) {
	double variance = standardDeviation*standardDeviation;
	double ofE = utilPa.roundDouble( -rhoI*rhoI / (2*variance), DOUBLE_PRECISION);
	return utilPa.roundDouble( exp(ofE) / (SQRT_2PI * standardDeviation), DOUBLE_PRECISION);
}

// alpha = 2*variance*numOfH*log(deltaH) - sumOfSquare(H)
double nttProHoughAccuracy::calculateAlpha(CvSeq *eightHPa, CvSeq *eightRhoPa) {
	double rhoMean = utilPa.mean(eightRhoPa);
	double variance = utilPa.variance(eightRhoPa, rhoMean);
	int numOfH = eightHPa->total;
	double logDeltaH = calculateLogDeltaH(eightHPa, eightRhoPa);
	double sumOfSquareOfH = utilPa.sumOfSquare(eightHPa);
	return 2*variance*numOfH*logDeltaH - sumOfSquareOfH;
}

// beta = - 2*H1*(e^..rho1..) - 2*H2*(e^..rho2..) - ....... - 2*Hn*(e^..rhoN..)
double nttProHoughAccuracy::calculateBeta(CvSeq *eightHPa, CvSeq *eightRhoPa) {
	double result = 0;
	double rhoMean = utilPa.mean(eightRhoPa);
	double variance = utilPa.variance(eightRhoPa, rhoMean);
	for (int i = 0; i < eightHPa->total; i++) {
		double* hI = (double*)cvGetSeqElem(eightHPa, i);
		double* rhoI = (double*)cvGetSeqElem(eightRhoPa, i);
		double ex = utilPa.roundDouble( exp( -(((*rhoI) - rhoMean) * ((*rhoI) - rhoMean) / (2 * variance)) ), DOUBLE_PRECISION);
		result += ( -2 * (*hI) * ex );
	}
	return result;
}

// gamma = (e^..rho1..)^2 + (e^..rho2..)^2 + ....... + (e^..rhoN..)^2
double nttProHoughAccuracy::calculateGamma(CvSeq *eightRhoPa) {
	double result = 0;
	double rhoMean = utilPa.mean(eightRhoPa);
	double variance = utilPa.variance(eightRhoPa, rhoMean);
	for (int i = 0; i < eightRhoPa->total; i++) {
		double* rhoI = (double*)cvGetSeqElem(eightRhoPa, i);
		double ex = utilPa.roundDouble( exp( -(((*rhoI) - rhoMean) * ((*rhoI) - rhoMean) / (2 * variance)) ), DOUBLE_PRECISION);
		result += ex*ex;
	}
	return result;
}

// log(deltaH) >= ((beta^2) / (4*gamma) - sumOfSquare(H)) / (2*variance*numOfH)
double nttProHoughAccuracy::calculateLogDeltaH(CvSeq *eightHPa, CvSeq *eightRhoPa) {
	double beta = calculateBeta(eightHPa, eightRhoPa);
	double gamma = calculateGamma(eightRhoPa);
	double sumOfSquareOfH = utilPa.sumOfSquare(eightHPa);

	double rhoMean = utilPa.mean(eightRhoPa);
	double variance = utilPa.variance(eightRhoPa, rhoMean);
	int numOfH = eightHPa->total;

	double a = utilPa.roundDouble( (beta*beta) / (4*gamma) , DOUBLE_PRECISION);
	double b = utilPa.roundDouble( 1 / (2*variance*numOfH) , DOUBLE_PRECISION);

	double more = 0.0001; // the smaller number, the more precision; but maybe cannot come to result
	return more + fabs( (a - sumOfSquareOfH) * b ); // take notice !!!
}

void nttProHoughAccuracy::init() {
	listOfHmaxStoragePa = cvCreateMemStorage(0);
	listOfHmaxPa = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(double), listOfHmaxStoragePa);

	listOfRhoHmaxStoragePa = cvCreateMemStorage(0);
	listOfRhoHmaxPa = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(double), listOfRhoHmaxStoragePa);

	listOfThetaHmaxStoragePa = cvCreateMemStorage(0);
	listOfThetaHmaxPa = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(double), listOfThetaHmaxStoragePa);

	eightHStoragePa = cvCreateMemStorage(0);
	eightHPa = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(double), eightHStoragePa);

	eightRhoStoragePa = cvCreateMemStorage(0);
	eightRhoPa = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(double), eightRhoStoragePa);

	removedPixelsStoragePa = cvCreateMemStorage(0);
	removedPixelsPa = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvPoint), removedPixelsStoragePa);
}