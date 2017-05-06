#pragma once

#include "cv.h"
#include "highgui.h"

class Util
{
public:
	Util(void);
public:
	~Util(void);

	// -----------------------------------------------------------

	static int generateRandomNumberOpenCV(int max) {
		CvRNG rng_state = cvRNG(0xffffffff);
		return cvRandInt(&rng_state) % max;
	}

	// http://www.codeproject.com/cpp/floatutils.asp?df=100&forumid=208&exp=0&select=14154
	static double roundDouble(double doValue, int nPrecision) {
		static const double doBase = 10.0;
		double doComplete5, doComplete5i;
	    
		doComplete5 = doValue * pow(doBase, (double) (nPrecision + 1));
	    
		if(doValue < 0.0) doComplete5 -= 5.0;
		else doComplete5 += 5.0;
	    
		doComplete5 /= doBase;
		modf(doComplete5, &doComplete5i);
	    
		return doComplete5i / pow(doBase, (double) nPrecision);
	}

	static char* ConvertToChar(const CString &s) {
		int nSize = s.GetLength();
		char *pAnsiString = new char[nSize+1];
		memset(pAnsiString,0,nSize+1);
		wcstombs(pAnsiString, s, nSize+1); // declare _CRT_SECURE_NO_DEPRECATE in stdafx.h
		return pAnsiString;
	} 

	static IplImage* subtractImage(IplImage* src1, IplImage* src2) {
		int step = src1->widthStep, width = src1->width, height = src1->height; 
		uchar *src1Data = (uchar *)src1->imageData;
		uchar *src2Data = (uchar *)src2->imageData;
		
		IplImage* result = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
		uchar *resultData = (uchar *)result->imageData;
		for (int i = 0; i < height; i++) 
			for (int j = 0; j < width; j++) {
				int position = i*step + j;
				if (src1Data[position] > src2Data[position]) resultData[position] = src1Data[position] - src2Data[position];
				else resultData[position] = 0;				
			}
	}

	// -----------------------------------------------------------

	static BYTE * CreateBMPInfo(int width, int height, int bpp) { 
		int total, sizeofbits;
		sizeofbits = width * bpp * height;
		total = sizeof(BITMAPINFOHEADER);
		if(bpp==1) total += 256*sizeof(RGBQUAD);
		total += sizeofbits;
		BYTE * pBuff = new BYTE[total];	
		
		BITMAPINFOHEADER * pBmih;
		RGBQUAD * pRgb4;

		//fill BITMAPINFOHEADER
		pBmih = (BITMAPINFOHEADER *)pBuff;
		pBmih->biSize=sizeof(BITMAPINFOHEADER);
		pBmih->biWidth=width;
		pBmih->biHeight=-height;		//chu y
		pBmih->biPlanes=1;
		pBmih->biBitCount=bpp*8;
		pBmih->biCompression=BI_RGB;
		pBmih->biSizeImage=width*height;
		pBmih->biXPelsPerMeter=0;
		pBmih->biYPelsPerMeter=0;
		pBmih->biClrUsed=0;
		pBmih->biClrImportant=0;
		
		if(bpp==1){
			pRgb4 = (RGBQUAD *)(pBuff + sizeof(BITMAPINFOHEADER));
			for(int i=0;i<256;i++){
				pRgb4->rgbRed=(BYTE)i;
				pRgb4->rgbGreen=(BYTE)i;
				pRgb4->rgbBlue=(BYTE)i;
				pRgb4->rgbReserved=(BYTE)0;		
				pRgb4++;
			}		
		}	
		return pBuff;
	}

};