#include "StdAfx.h"
#include "nttUtil.h"
#include "math.h"
#include "cv.h"

nttUtil::nttUtil(void) { }
nttUtil::~nttUtil(void) { }

bool nttUtil::isArrayContainAnInteger(CArray<int, int> *cArray, int n) {
	for (int i = 0; i <= cArray->GetUpperBound(); i++)
		if (cArray->GetAt(i) == n) return true;
	return false;
}

int nttUtil::generateRandomNumber(int max) {
	// http://members.cox.net/srice1/random/crandom.html
	//double maxOfRandom = (double)(RAND_MAX) + (double)(1);
	//double fractionOfRandom = (double)rand() / ((double)(RAND_MAX) + (double)(1));
	//return (int) ( fractionOfRandom * (max + 1) );
	return (int) ( (double)rand() / ((double)(RAND_MAX) + (double)(1)) * (max + 1) );
}

int nttUtil::generateRandomNumberOpenCV(int max) {
	CvRNG rng_state = cvRNG(0xffffffff);
	return cvRandInt(&rng_state) % max;
}

// http://www.codeproject.com/cpp/floatutils.asp?df=100&forumid=208&exp=0&select=14154
double nttUtil::roundDouble(double doValue, int nPrecision) {
    static const double doBase = 10.0;
    double doComplete5, doComplete5i;
    
    doComplete5 = doValue * pow(doBase, (double) (nPrecision + 1));
    
    if(doValue < 0.0) doComplete5 -= 5.0;
    else doComplete5 += 5.0;
    
    doComplete5 /= doBase;
    modf(doComplete5, &doComplete5i);
    
    return doComplete5i / pow(doBase, (double) nPrecision);
}

double nttUtil::factorial(int n) {
	double result = 1;
	for (int i = 2; i <= n; i++) result *= i;
	return result;
}

double nttUtil::mean(CvSeq *seq) {
	double total = 0;
	for(int i = 0; i < seq->total; i++) {
		double* v = (double*)cvGetSeqElem(seq, i);
		total += (*v);
	}
	return roundDouble(total / seq->total, 7);
}

double nttUtil::sumOfSquare(CvSeq *seq) {
	double total = 0;
	for(int i = 0; i < seq->total; i++) {
		double* v = (double*)cvGetSeqElem(seq, i);
		total += (*v) * (*v);
	}
	return total;
}

double nttUtil::variance(CvSeq *seq, double mean) {
	double total = 0;
	for(int i = 0; i < seq->total; i++) {
		double* v = (double*)cvGetSeqElem(seq, i);
		total += ((*v) - mean) * ((*v) - mean);
	}
	return roundDouble(total / seq->total, 7);
}

/**************************** optimization with ASM code ****************************/
float nttUtil::addFloat(float f1, float f2) {
	float fsum;
	_asm { 
		finit 
		fld    DWORD PTR f1 
		fadd   DWORD PTR f2 
		fstp   DWORD PTR fsum 
		fwait 
    }; 
	return fsum;
}

// http://www.flipcode.com/files/code/CosSin.cpp
// usage:
//		float f = 1.23;
//		int i;
//		floatToInt(&i, f);
//
void nttUtil::floatToInt(int *int_pointer, float f) {
	__asm  fld  f
	__asm  mov  edx,int_pointer
	__asm  FRNDINT
	__asm  fistp dword ptr [edx];
}