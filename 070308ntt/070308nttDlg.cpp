// 070308nttDlg.cpp : implementation file

#include "stdafx.h"
#include "070308ntt.h"
#include "070308nttDlg.h"
#include "cv.h"
#include "highgui.h"
#include <math.h>
#include "nttUtil.h"

#include "nttStandardHough.h"
#include "nttStandardHoughNew.h"
#include "nttProbabilisticHough.h"
#include "nttProgressiveProbabilisticHough.h"
#include "nttProHoughAccuracy.h"
#include "nttHoughGreen.h"
#include "nttPyramidHough.h"
#include "LineExtraction.h"
#include "LineExtractionSHT.h"
#include "nttStandardHoughEx.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// ntt ------------
nttStandardHough myStandardHough;
nttStandardHoughNew myStandardHoughNew;
nttProbabilisticHough myProbabilisticHough;
nttProgressiveProbabilisticHough myProgressiveProbabilisticHough;
nttProHoughAccuracy myProHoughAccuracy;
nttHoughGreen myHoughGreen;
nttPyramidHough myPyramidHough;
nttStandardHoughEx myStandardHoughEx;
LineExtraction myLineExtraction;
LineExtractionSHT myLineExtractionSHT;

IplImage *srcImage, *grayImage, *bwImage, *boundaryImage, *lineExtractImage, *printImage;
int idcPicture_with, idcPicture_height, numOfSavedImage;
BYTE *bitmapInfo;

FILE *outfile; // dung de luu tick count

HWND m_hOut1Wnd, m_hOut2Wnd, m_hOut3Wnd, m_hOut4Wnd;
// ----------------

// CAboutDlg dialog used for App About
class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CMy070308nttDlg dialog
CMy070308nttDlg::CMy070308nttDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CMy070308nttDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMy070308nttDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CMy070308nttDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_BUTTON_LOAD_AN_IMAGE, &CMy070308nttDlg::OnBnClickedButtonLoadAnImage)
	ON_BN_CLICKED(IDC_BUTTON_QUIT, &CMy070308nttDlg::OnBnClickedButtonQuit)
	ON_BN_CLICKED(IDC_BUTTON_SAVE_RESULT_IMG, &CMy070308nttDlg::OnBnClickedButtonSaveResultImage)

	ON_BN_CLICKED(IDC_BUTTON_OPENCV_HOUGHLINE2, &CMy070308nttDlg::OnBnClickedButtonOpencvHoughline2)
	ON_BN_CLICKED(IDC_BUTTON_NTT_STANDARD_HOUGH, &CMy070308nttDlg::OnBnClickedButtonNttStandardHough)
	ON_BN_CLICKED(IDC_BUTTON_NTT_STANDARD_HOUGH_NEW, &CMy070308nttDlg::OnBnClickedButtonNttStandardHoughNew)	
	ON_BN_CLICKED(IDC_BUTTON_NTT_PROBABILISTIC_HOUGH, &CMy070308nttDlg::OnBnClickedButtonNttProbabilisticHough)	
	ON_BN_CLICKED(IDC_BUTTON_NTT_PROGRESSIVE_PROBABILISTIC_HOUGH, &CMy070308nttDlg::OnBnClickedButtonNttProgressiveProbabilisticHough)
	ON_BN_CLICKED(IDC_BUTTON_NTT_HOUGH_ACCURACY, &CMy070308nttDlg::OnBnClickedButtonNttHoughAccuracy)
	ON_BN_CLICKED(IDC_BUTTON_NTT_HOUGH_GREEN, &CMy070308nttDlg::OnBnClickedButtonNttHoughGreen)
	ON_BN_CLICKED(IDC_BUTTON_NTT_PYRAMID_HOUGH, &CMy070308nttDlg::OnBnClickedButtonNttPyramidHough)
	ON_BN_CLICKED(IDC_BUTTON_LAST_DECISION, &CMy070308nttDlg::OnBnClickedButtonLastDecision)
END_MESSAGE_MAP()


// CMy070308nttDlg message handlers

BOOL CMy070308nttDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here
	
	Init(); // ntt

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CMy070308nttDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CMy070308nttDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CMy070308nttDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

// ===================================================================

void CMy070308nttDlg::OnBnClickedButtonLoadAnImage() {
	CFileDialog dlg(TRUE, _T("*.avi"), _T(""), OFN_FILEMUSTEXIST|OFN_PATHMUSTEXIST|OFN_HIDEREADONLY,
						_T("Image (*.bmp)|*.bmp|All Files (*.*)|*.*||"),NULL);	
	if (dlg.DoModal() != IDOK) return;
	
	// load image
	char filenameInMultiByte[256];
	WideCharToMultiByte(CP_ACP, 0, dlg.GetPathName(), -1, filenameInMultiByte, 256, NULL, NULL);	
	srcImage = cvvLoadImage(filenameInMultiByte);

	UpdateScreen();
}

void CMy070308nttDlg::OnBnClickedButtonQuit() {	
	Sleep(500); //Wait all thread to finish his current jobs
	cvDestroyAllWindows();
	ReleaseResources();   
	OnOK();
}

void CMy070308nttDlg::OnBnClickedButtonSaveResultImage() {
	CString fn;
	fn.Format(L".\\%04d.bmp", numOfSavedImage++);
	CStringA filename (fn); // convert CString to char *
	//cvSaveImage(filename, lineExtractImage);
	cvSaveImage(filename, printImage);

	fn.Format(L".\\edge%04d.bmp", 1);
	CStringA filename2 (fn); // convert CString to char *

	for(int i = 0; i < boundaryImage->height; i++) 
		for(int j = 0; j < boundaryImage->width; j++) 
			boundaryImage->imageData[i*boundaryImage->widthStep + j] = 255 - boundaryImage->imageData[i*boundaryImage->widthStep + j];
	cvSaveImage(filename2, boundaryImage);
}

void CMy070308nttDlg::OnBnClickedButtonOpencvHoughline2() {
	// ----------------------------------------------------
	/*printImage = cvCreateImage(cvGetSize(srcImage), 8, 3); 
	int pw = printImage->width, ph = printImage->height, ps = printImage->widthStep, pc = printImage->nChannels;
	for(int i = 0; i < ph; i++) 
		for(int j = 0; j < pw; j++) 
			for (int k = 0; k < pc; k++)
				printImage->imageData[i*ps + j*pc + k] = 255;*/
	// ----------------------------------------------------



	// ----------------- tick count START -----------------
	__int64 freq, tStart, tStop;
	unsigned long TimeDiff;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq); // Get the frequency of the hi-res timer
	QueryPerformanceCounter((LARGE_INTEGER*)&tStart); // Assuming that has worked you can then use hi-res timer 
	// ----------------------------------------------------



	ConvertSrcImageToGrayscale();
	DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY); // convert gray image to binary image

	lineExtractImage = cvCloneImage(srcImage); // dung de ve len ddo' cac line extract
	boundaryImage = cvCreateImage(cvGetSize(bwImage), 8, 1); 
	
	cvCanny(grayImage, boundaryImage, 50, 200, 3); 
	//boundaryImage = cvCloneImage(grayImage); // KHONG DUNG CANNY
	
	CvSeq* lines = 0;
	
	CvMemStorage* storage = cvCreateMemStorage(0); // block_size = 0
		// Size of the storage blocks in bytes. If it is 0, the block size is set to default value (~64K)
		// The function cvCreateMemStorage creates a memory storage and returns pointer to it. 
		// Initially the storage is empty. All fields of the header, except the block_size, are set to 0.

	// -----------
	// rho, theta edit control
	CWnd *editControl;	
	CString rhoString, thetaString, thresholdString;
	editControl = this->GetDlgItem(IDC_EDIT_OPENCV_HOUGH_RHO);		editControl->GetWindowTextW(rhoString);
	editControl = this->GetDlgItem(IDC_EDIT_OPENCV_HOUGH_THETA);	editControl->GetWindowTextW(thetaString);
	editControl = this->GetDlgItem(IDC_EDIT_OPENCV_HOUGH_THRESHOLD);editControl->GetWindowTextW(thresholdString);

	// IDC_EDIT_OPENCV_HOUGH_METHOD combo box
	CComboBox *methodCombo;
	methodCombo = (CComboBox *) this->GetDlgItem(IDC_COMBO_OPENCV_HOUGH_METHOD);	

	double rho = 1; // Distance resolution in pixel-related units
	//double theta = CV_PI/180; // Angle resolution measured in radians
	double theta = 0.01;
	int threshold = 120;	

	int minLineLength = 20;
	int maxGapBetweenLineSegments = 6;
	
	if (!rhoString.IsEmpty() && !thetaString.IsEmpty() && !thresholdString.IsEmpty()) {
		rho = wcstod(rhoString, NULL);
		theta = wcstod(thetaString, NULL);
		threshold = cvRound(wcstod(thresholdString, NULL)); // convert to (int) threshold

		switch (methodCombo->GetCurSel()) {

			case 0: 
				MessageBox(L"CV_HOUGH_MULTI_SCALE: not yet implement");
				break;

			case 1 :
				//MessageBox(L"CV_HOUGH_PROBABILISTIC");		

				lines = cvHoughLines2(boundaryImage, storage, CV_HOUGH_PROBABILISTIC, rho, theta, threshold, minLineLength, maxGapBetweenLineSegments); 
				for(int i = 0; i < lines->total; i++) {
					CvPoint* line = (CvPoint*)cvGetSeqElem(lines,i);					
					//cvLine(printImage, line[0], line[1], CV_RGB(255,0,0), 1, 8); // ve len lineExtractImage, cac line extract tu boundaryImage
					cvLine(lineExtractImage, line[0], line[1], CV_RGB(255,0,0), 1, 8); // ve len lineExtractImage, cac line extract tu boundaryImage
				} 

				break;

			case 2:  
				//MessageBox(L"CV_HOUGH_STANDARD");
				lines = cvHoughLines2(boundaryImage, storage, CV_HOUGH_STANDARD, rho, theta, threshold, 0, 0 );
				int linesMax = 20;
				for(int i = 0; i < MIN(lines->total, linesMax); i++) {
					float* line = (float*)cvGetSeqElem(lines,i); // char* cvGetSeqElem(const CvSeq* seq, int index);						
					float rho = line[0], theta = line[1];
					CvPoint pt1, pt2;
					double a = cos(theta), b = sin(theta);
					double x0 = a*rho, y0 = b*rho;
					pt1.x = cvRound(x0 + 1000*(-b));	pt2.x = cvRound(x0 - 1000*(-b));
					pt1.y = cvRound(y0 + 1000*(a));     pt2.y = cvRound(y0 - 1000*(a));
					cvLine( lineExtractImage, pt1, pt2, CV_RGB(255,0,0), 1, 8 );
				} 
				break;
		}

	} else {
		//MessageBox(L"Default: CV_HOUGH_PROBABILISTIC");
		lines = cvHoughLines2(boundaryImage, storage, CV_HOUGH_PROBABILISTIC, rho, theta, threshold, minLineLength, maxGapBetweenLineSegments); 
				// CvSeq* cvHoughLines2(CvArr* image, void* line_storage, int method, double rho, 
				//						double theta, int threshold, double param1=0, double param2=0);	
		for(int i = 0; i < lines->total; i++) {
			CvPoint* line = (CvPoint*)cvGetSeqElem(lines,i);
			cvLine(lineExtractImage, line[0], line[1], CV_RGB(255,0,0), 1); // ve len lineExtractImage, cac line extract tu boundaryImage
		}
	}

	cvNamedWindow("prob", CV_WINDOW_AUTOSIZE);
	cvShowImage("prob", lineExtractImage);
	//cvShowImage("prob", boundaryImage);	
	//cvShowImage("prob", printImage);

	// hien thi ket qua bao nhieu lines detect duoc
	CString outTextString;		outTextString.Format(L"%3d", lines->total);
	GetDlgItem(IDC_STATIC_NUM_DETECTED_LINES)->SetWindowTextW((LPCTSTR)outTextString); 

	// --- save result to BMP file
/*	CString fn;
	fn.Format(L".\\%04d.bmp", numOfSavedImage);
	CStringA filename (fn); // convert CString to char *
	cvSaveImage(filename, lineExtractImage);
*/


	// ----------------- tick count STOP -----------------	
	QueryPerformanceCounter((LARGE_INTEGER*)&tStop); // Perform operations that require timing
	TimeDiff = (unsigned long)(((tStop - tStart) * 1000000) / freq); // Calculate time difference

	// write TimeDiff to file
	//FILE *outfile; 
	outfile = fopen("TimeDiff.txt", "a"); 
	fprintf(outfile, "%10d \n", TimeDiff);
	fclose(outfile);
	// ----------------------------------------------------




	// =======================
	UpdateScreen();
}

void CMy070308nttDlg::OnBnClickedButtonNttStandardHough() { 
	ConvertSrcImageToGrayscale();
	DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY); // convert gray image to binary image

	lineExtractImage = cvCloneImage(srcImage); // dung de ve len ddo' cac line extract
	boundaryImage = cvCreateImage(cvGetSize(bwImage), 8, 1); 
	cvCanny(grayImage, boundaryImage, 50, 200, 3);

	// -------------



	// ----------------- tick count START -----------------
	__int64 freq, tStart, tStop;
	unsigned long TimeDiff;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq); // Get the frequency of the hi-res timer
	QueryPerformanceCounter((LARGE_INTEGER*)&tStart); // Assuming that has worked you can then use hi-res timer 
	// ----------------------------------------------------



	double rhoResolution = 1; 
	double thetaResolution = 0.01; // in radian
	int voteThreshold = 100; 
	int linesMax = 100;

	// -------- SHT -------- 
	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* result = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvRect), storage);
	myStandardHough.run(boundaryImage, result, rhoResolution, thetaResolution, voteThreshold, linesMax);
 
	for(int i = 0; i < result->total; i++) {	
		CvRect* lr = (CvRect*)cvGetSeqElem(result,i);
		CvPoint startPoint = { lr->x, lr->y };
		CvPoint endPoint = { lr->width, lr->height };
		cvLine(lineExtractImage, startPoint, endPoint, CV_RGB(255,0,0), 1); // ve len lineExtractImage, cac line extract tu boundaryImage
	}



	// ----------------- tick count STOP -----------------	
	QueryPerformanceCounter((LARGE_INTEGER*)&tStop); // Perform operations that require timing
	TimeDiff = (unsigned long)(((tStop - tStart) * 1000000) / freq); // Calculate time difference

	// write TimeDiff to file
	//FILE *outfile; 
	outfile = fopen("TimeDiff.txt", "a"); 
	fprintf(outfile, "%10d \n", TimeDiff);
	fclose(outfile);
	// ----------------------------------------------------



	// -------- SHT-more --------
	/*int maxGap = 4;
	int lineLength = 10;

	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* result = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(LineExtracted), storage);
	myLineExtractionSHT.run(boundaryImage, result, rhoResolution, thetaResolution, voteThreshold, linesMax, maxGap, lineLength);
 
	for(int i = 0; i < result->total; i++) {	
		LineExtracted* lr = (LineExtracted*)cvGetSeqElem(result,i);
		CvPoint startPoint = { lr->startX, lr->startY };
		CvPoint endPoint = { lr->endX, lr->endY };
		double theta = lr->theta;
		double rho = lr->rho;
		cvLine(lineExtractImage, startPoint, endPoint, CV_RGB(255,0,0), 3);
	}*/

	// --- save edge-image to BMP file
/*	CString fnEdge;
	fnEdge.Format(L".\\_edge%04d.bmp", ++numOfSavedImage);
	CStringA filenameEdge (fnEdge); // convert CString to char *
	cvSaveImage(filenameEdge, boundaryImage);

	// --- save result to BMP file
	CString fn;
	fn.Format(L".\\%04d.bmp", numOfSavedImage);
	CStringA filename (fn); // convert CString to char *
	cvSaveImage(filename, lineExtractImage);
*/
	cvNamedWindow("prob", CV_WINDOW_AUTOSIZE);
	cvShowImage("prob", lineExtractImage);

	//
	cvReleaseMemStorage( &storage );
	UpdateScreen();
}

void CMy070308nttDlg::OnBnClickedButtonNttStandardHoughNew() { // tich luy co tro.ng so^' 
	ConvertSrcImageToGrayscale();
	DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY); // convert gray image to binary image

	lineExtractImage = cvCloneImage(srcImage); // dung de ve len ddo' cac line extract
	boundaryImage = cvCreateImage(cvGetSize(bwImage), 8, 1); 
	cvCanny(grayImage, boundaryImage, 50, 200, 3);

	// -------------

	CWnd *editControl;	
	CString rhoString, thetaString, thresholdString, linesMaxString;
	editControl = this->GetDlgItem(IDC_EDIT_SHT_NEW_RHO);		editControl->GetWindowTextW(rhoString);
	editControl = this->GetDlgItem(IDC_EDIT_SHT_NEW_THETA);		editControl->GetWindowTextW(thetaString);
	editControl = this->GetDlgItem(IDC_EDIT_SHT_NEW_THRESHOLD);	editControl->GetWindowTextW(thresholdString);
	editControl = this->GetDlgItem(IDC_EDIT_SHT_NEW_LINESMAX);	editControl->GetWindowTextW(linesMaxString);

	double rhoResolution = 1; 
	double thetaResolution = 0.01; // in radian
	int voteThreshold = 70; 
	int linesMax = 20;

	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* result = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvRect), storage);

	if (!rhoString.IsEmpty() && !thetaString.IsEmpty() && !thresholdString.IsEmpty() && !linesMaxString.IsEmpty()) {
		rhoResolution = wcstod(rhoString, NULL);
		thetaResolution = wcstod(thetaString, NULL);
		voteThreshold = cvRound(wcstod(thresholdString, NULL)); // convert to (int) threshold
		linesMax = cvRound(wcstod(linesMaxString, NULL)); 		
	} 
	myStandardHoughNew.run(boundaryImage, result, rhoResolution, thetaResolution, voteThreshold, linesMax);
 
	// draw lines
	for(int i = 0; i < result->total; i++) {	
		CvRect* lr = (CvRect*)cvGetSeqElem(result,i);
		CvPoint startPoint = { lr->x, lr->y };
		CvPoint endPoint = { lr->width, lr->height };
		cvLine(lineExtractImage, startPoint, endPoint, CV_RGB(255,0,0), 3); // ve len lineExtractImage, cac line extract tu boundaryImage
	}

	cvReleaseMemStorage( &storage );
	UpdateScreen();
}

void CMy070308nttDlg::OnBnClickedButtonLastDecision() { 
	// ----------------------------------------------------
	printImage = cvCreateImage(cvGetSize(srcImage), 8, 3); 
	int pw = printImage->width, ph = printImage->height, ps = printImage->widthStep, pc = printImage->nChannels;
	for(int i = 0; i < ph; i++) 
		for(int j = 0; j < pw; j++) 
			for (int k = 0; k < pc; k++)
				printImage->imageData[i*ps + j*pc + k] = 255;
	// ----------------------------------------------------

	ConvertSrcImageToGrayscale();
	DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY); // convert gray image to binary image

	lineExtractImage = cvCloneImage(srcImage); // dung de ve len ddo' cac line extract
	boundaryImage = cvCreateImage(cvGetSize(bwImage), 8, 1); 
	cvCanny(grayImage, boundaryImage, 50, 200, 3);
	//cvCanny(grayImage, boundaryImage, 100, 250, 3);

	// -------------

	CWnd *editControl;	

	CString rhoString, thetaString, thresholdString, linesMaxString, maxGapString, lineLengthString;	
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_RHO);			editControl->GetWindowTextW(rhoString);
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_THETA);		editControl->GetWindowTextW(thetaString);
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_THRESHOLD);	editControl->GetWindowTextW(thresholdString);
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_LINESMAX);	editControl->GetWindowTextW(linesMaxString);
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_MAXGAP);		editControl->GetWindowTextW(maxGapString);
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_LINELEN);		editControl->GetWindowTextW(lineLengthString);

	CString lineInLengthFromString, lineInLengthToString, declinationThetaFromString, declinationThetaToString;
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_LINEINLENGTH_FROM);	editControl->GetWindowTextW(lineInLengthFromString);
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_LINEINLENGTH_TO);		editControl->GetWindowTextW(lineInLengthToString);
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_DECLINATION_THETAFROM);	editControl->GetWindowTextW(declinationThetaFromString);
	editControl = this->GetDlgItem(IDC_EDIT_LAST_DECISION_DECLINATION_THETATO);	editControl->GetWindowTextW(declinationThetaToString);

	double rhoResolution = 1; 
	double thetaResolution = 0.01; // in radian
	int voteThreshold = 120; 
	int linesMax = 100;
	int maxGap = 6; 
	int lineLength = 50;

	int lineInLengthFrom = -1;
	int lineInLengthTo = -1;
	double thetaFrom = -1; // in radian
	double thetaTo = -1; // in radian

	if (!rhoString.IsEmpty() && !thetaString.IsEmpty() && !thresholdString.IsEmpty() && !linesMaxString.IsEmpty() 
							 && !maxGapString.IsEmpty() && !lineLengthString.IsEmpty()) {
		rhoResolution = wcstod(rhoString, NULL);
		thetaResolution = wcstod(thetaString, NULL);
		voteThreshold = cvRound(wcstod(thresholdString, NULL)); // convert to (int) threshold
		linesMax = cvRound(wcstod(linesMaxString, NULL));
		maxGap = cvRound(wcstod(maxGapString, NULL));
		lineLength = cvRound(wcstod(lineLengthString, NULL));
	} 

	if (!lineInLengthFromString.IsEmpty() && !lineInLengthToString.IsEmpty()) {
		lineInLengthFrom = cvRound(wcstod(lineInLengthFromString, NULL));
		lineInLengthTo = cvRound(wcstod(lineInLengthToString, NULL));
	} 
	
	if (!declinationThetaFromString.IsEmpty() && !declinationThetaToString.IsEmpty()) {
		thetaFrom = wcstod(declinationThetaFromString, NULL);
		thetaTo = wcstod(declinationThetaToString, NULL);
	}

	//
	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* result = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(LineExtracted), storage);
	//myLineExtraction.run(boundaryImage, result, rhoResolution, thetaResolution, voteThreshold, linesMax, maxGap, lineLength,
	//					 lineInLengthFrom, lineInLengthTo, thetaFrom, thetaTo);

	myStandardHoughEx.run(boundaryImage, result, rhoResolution, thetaResolution, voteThreshold, linesMax, lineLength, maxGap, 3, 3);
 
	// draw lines
	for(int i = 0; i < result->total; i++) {	
		LineExtracted* lr = (LineExtracted*)cvGetSeqElem(result,i);
		CvPoint startPoint = { lr->startX, lr->startY };
		CvPoint endPoint = { lr->endX, lr->endY };
		double theta = lr->theta;
		double rho = lr->rho;
		//cvLine(printImage, startPoint, endPoint, CV_RGB(255,0,0), 1); 
		cvLine(lineExtractImage, startPoint, endPoint, CV_RGB(255,0,0), 1); // ve len lineExtractImage, cac line extract tu boundaryImage
	}

	cvNamedWindow("ntt", CV_WINDOW_AUTOSIZE);
	cvShowImage("ntt", lineExtractImage);
	//cvShowImage("ntt", printImage);

	//
	cvReleaseMemStorage( &storage );
	UpdateScreen();
}

void CMy070308nttDlg::OnBnClickedButtonNttProbabilisticHough() { 
	ConvertSrcImageToGrayscale();
	//DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY_INV); // convert gray image to binary image
	DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY); // convert gray image to binary image

	lineExtractImage = cvCloneImage(srcImage); // dung de ve len ddo' cac line extract
	boundaryImage = cvCreateImage(cvGetSize(bwImage), 8, 1); 
	cvCanny(grayImage, boundaryImage, 50, 300, 3);
	//boundaryImage = bwImage;
	//cvSobel(grayImage, boundaryImage, 2, 1, 5);

	// -------------
	


	// ----------------- tick count START -----------------
	__int64 freq, tStart, tStop;
	unsigned long TimeDiff;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq); // Get the frequency of the hi-res timer
	QueryPerformanceCounter((LARGE_INTEGER*)&tStart); // Assuming that has worked you can then use hi-res timer 
	// ----------------------------------------------------



	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* result = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvRect), storage);

	double rhoResolution = 1; // in pixel 
	double thetaResolution = 0.1; // in radian
	int voteThreshold = 100; 
	int minLineLength = 20;
	int maxGap = 6;	
	myProbabilisticHough.run(boundaryImage, result, rhoResolution, thetaResolution, voteThreshold, minLineLength, maxGap);
 
	for(int i = 0; i < result->total; i++) {	
		CvRect* lr = (CvRect*)cvGetSeqElem(result,i);
		CvPoint startPoint = { lr->x, lr->y };
		CvPoint endPoint = { lr->width, lr->height };
		cvLine(lineExtractImage, startPoint, endPoint, CV_RGB(255,0,0), 2); // ve len lineExtractImage, cac line extract tu boundaryImage
	}



	// ----------------- tick count STOP -----------------	
	QueryPerformanceCounter((LARGE_INTEGER*)&tStop); // Perform operations that require timing
	TimeDiff = (unsigned long)(((tStop - tStart) * 1000000) / freq); // Calculate time difference

	// write TimeDiff to file
	//FILE *outfile; 
	outfile = fopen("TimeDiff.txt", "a"); 
	fprintf(outfile, "%10d \n", TimeDiff);
	fclose(outfile);
	// ----------------------------------------------------




	// --- save result to BMP file	
/*	CString fn;
	fn.Format(L".\\_result%04d.bmp", ++numOfSavedImage);
	CStringA filename (fn); // convert CString to char *
	cvSaveImage(filename, lineExtractImage);

	// --- save edge-image to BMP file
	CString fnEdge;
	fnEdge.Format(L".\\_edge%04d.bmp", numOfSavedImage);
	CStringA filenameEdge (fnEdge); // convert CString to char *
	cvSaveImage(filenameEdge, boundaryImage);
*/
	//
	cvReleaseMemStorage( &storage );
	UpdateScreen();
}

void CMy070308nttDlg::OnBnClickedButtonNttProgressiveProbabilisticHough() { 
	ConvertSrcImageToGrayscale();
	DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY); // convert gray image to binary image

	lineExtractImage = cvCloneImage(srcImage); // dung de ve len ddo' cac line extract
	boundaryImage = cvCreateImage(cvGetSize(bwImage), 8, 1); 
	cvCanny(grayImage, boundaryImage, 50, 200, 3);

	// -------------
	


	// ----------------- tick count START -----------------
	__int64 freq, tStart, tStop;
	unsigned long TimeDiff;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq); // Get the frequency of the hi-res timer
	QueryPerformanceCounter((LARGE_INTEGER*)&tStart); // Assuming that has worked you can then use hi-res timer 
	// ----------------------------------------------------



	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* result = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvRect), storage);

	double rhoResolution = 1; // in pixel 
	double thetaResolution = 0.01; // in radian
	int voteThreshold = 0; // not used
	int minLineLength = 20;
	int maxGap = 6;	
	myProgressiveProbabilisticHough.run(boundaryImage, result, rhoResolution, thetaResolution, voteThreshold, minLineLength, maxGap);
 
	for(int i = 0; i < result->total; i++) {	
		CvRect* lr = (CvRect*)cvGetSeqElem(result,i);
		CvPoint startPoint = { lr->x, lr->y };
		CvPoint endPoint = { lr->width, lr->height };
		cvLine(lineExtractImage, startPoint, endPoint, CV_RGB(255,0,0), 3); // ve len lineExtractImage, cac line extract tu boundaryImage
	}
	cvReleaseMemStorage( &storage );


	// ----------------- tick count STOP -----------------	
	QueryPerformanceCounter((LARGE_INTEGER*)&tStop); // Perform operations that require timing
	TimeDiff = (unsigned long)(((tStop - tStart) * 1000000) / freq); // Calculate time difference

	// write TimeDiff to file
	//FILE *outfile; 
	outfile = fopen("TimeDiff.txt", "a"); 
	fprintf(outfile, "%10d \n", TimeDiff);
	fclose(outfile);
	// ----------------------------------------------------



	UpdateScreen();
}

void CMy070308nttDlg::OnBnClickedButtonNttHoughAccuracy() {
	ConvertSrcImageToGrayscale();
	DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY); // convert gray image to binary image

	lineExtractImage = cvCloneImage(srcImage); // dung de ve len ddo' cac line extract
	boundaryImage = cvCreateImage(cvGetSize(bwImage), 8, 1); 
	cvCanny(grayImage, boundaryImage, 50, 200, 3);

	// -------------
	
	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* result = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvRect), storage);

	double rhoResolution = 1; 
	double thetaResolution = 0.01; // in radian
	double voteThreshold = 0.05; 
	int minLineLength = 20;
	int maxGap = 6;	
	myProHoughAccuracy.run(boundaryImage, result, rhoResolution, thetaResolution, voteThreshold, minLineLength, maxGap);

	for(int i = 0; i < result->total; i++) {	
		CvRect* lr = (CvRect*)cvGetSeqElem(result,i);
		CvPoint startPoint = { lr->x, lr->y };
		CvPoint endPoint = { lr->width, lr->height };
		cvLine(lineExtractImage, startPoint, endPoint, CV_RGB(255,0,0), 3); // ve len lineExtractImage, cac line extract tu boundaryImage
	}
	cvReleaseMemStorage( &storage );	

	UpdateScreen();
}

void CMy070308nttDlg::OnBnClickedButtonNttHoughGreen() {
	ConvertSrcImageToGrayscale();
	DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY); // convert gray image to binary image

	lineExtractImage = cvCloneImage(srcImage); // dung de ve len ddo' cac line extract
	IplImage* edgeImage = cvCreateImage(cvGetSize(bwImage), 8, 1); 
	cvCanny(grayImage, edgeImage, 50, 200, 3);

	/* ------ cvDrawContours ------ */
	boundaryImage = cvCreateImage( cvGetSize(grayImage), 8, 3 );
    CvMemStorage* storage = cvCreateMemStorage(0);
    CvSeq* contour = 0;

    cvFindContours( edgeImage, storage, &contour, sizeof(CvContour), CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE );
    cvZero( boundaryImage );
    for( ; contour != 0; contour = contour->h_next ) {        
		CvScalar color = CV_RGB( rand()&255, rand()&255, rand()&255 );
		cvDrawContours( boundaryImage, contour, color, color, -1, CV_FILLED, 8 ); // replace CV_FILLED with 1 to see the outlines
    }

    //cvNamedWindow( "Components", 1 );
    //cvShowImage( "Components", boundaryImage );
	//cvWaitKey(0);

	cvReleaseMemStorage(&storage);
	/* --------------------------------- */

	CvMemStorage* resultStorage = cvCreateMemStorage(0);
	CvSeq* result = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvRect), resultStorage);

	double thetaResolution = CV_PI/180;
	double threshold = 50;
	myHoughGreen.run(edgeImage, result, thetaResolution, threshold);

	for(int i = 0; i < result->total; i++) {	
		CvRect* lr = (CvRect*)cvGetSeqElem(result,i);
		CvPoint startPoint = { lr->x, lr->y };
		CvPoint endPoint = { lr->width, lr->height };
		cvLine(lineExtractImage, startPoint, endPoint, CV_RGB(255,0,0), 3); // ve len lineExtractImage, cac line extract tu boundaryImage
	}
	cvReleaseMemStorage( &resultStorage );
	cvReleaseImage(&edgeImage);

	UpdateScreen();
}

void CMy070308nttDlg::OnBnClickedButtonNttPyramidHough() {
	ConvertSrcImageToGrayscale();
	DoThresholdOnGrayImage(127, 255, CV_THRESH_BINARY); // convert gray image to binary image

	lineExtractImage = cvCloneImage(srcImage); // dung de ve len ddo' cac line extract
	boundaryImage = cvCreateImage(cvGetSize(bwImage), 8, 1); 
	cvCanny(grayImage, boundaryImage, 50, 200, 3);

	// -------------
	
	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* result = cvCreateSeq(CV_SEQ_ELTYPE_GENERIC, sizeof(CvSeq), sizeof(CvRect), storage);

	double rhoResolution = 1; 
	double thetaResolution = 0.01; // in radian
	int voteThreshold = 25; 
	int linesMax = 50;
	myPyramidHough.run(boundaryImage, result, rhoResolution, thetaResolution, voteThreshold, linesMax);

	for(int i = 0; i < result->total; i++) {	
		CvRect* lr = (CvRect*)cvGetSeqElem(result,i);
		CvPoint startPoint = { lr->x, lr->y };
		CvPoint endPoint = { lr->width, lr->height };
		cvLine(lineExtractImage, startPoint, endPoint, CV_RGB(255,0,0), 3); // ve len lineExtractImage, cac line extract tu boundaryImage
	}
	cvReleaseMemStorage( &storage );

	UpdateScreen();
}

// ===================================================================

void CMy070308nttDlg::Init() {
	// init N picture controls
	CRect rect;
	CWnd * picWindow;

	picWindow = this->GetDlgItem(IDC_PICTURE1);
	picWindow->GetWindowRect(&rect); // > GetClientRect(&rect);
	this->ScreenToClient(rect); // map to client co-ords
	idcPicture_with = rect.Width();
	idcPicture_height = rect.Height();
	m_hOut1Wnd = picWindow->GetSafeHwnd();	

	picWindow = this->GetDlgItem(IDC_PICTURE2);	m_hOut2Wnd = picWindow->GetSafeHwnd();
	picWindow = this->GetDlgItem(IDC_PICTURE3);	m_hOut3Wnd = picWindow->GetSafeHwnd();	
	picWindow = this->GetDlgItem(IDC_PICTURE4);	m_hOut4Wnd = picWindow->GetSafeHwnd();	

	// init IDC_COMBO_OPENCV_HOUGH_METHOD
	CComboBox *openCVHoughMethodCombo;
	openCVHoughMethodCombo = (CComboBox *) this->GetDlgItem(IDC_COMBO_OPENCV_HOUGH_METHOD);
	openCVHoughMethodCombo->AddString(L"MULTI_SCALE"); // 0
	openCVHoughMethodCombo->AddString(L"PROBABILISTIC"); // 1	
	openCVHoughMethodCombo->AddString(L"STANDARD"); // 2
	openCVHoughMethodCombo->SetCurSel(2);

	//
	CString s;
	s.Format(L"%1d", 1);			this->GetDlgItem(IDC_EDIT_OPENCV_HOUGH_RHO)->SetWindowTextW((LPCTSTR)s);
	s.Format(L"%1.2f", 0.01);		this->GetDlgItem(IDC_EDIT_OPENCV_HOUGH_THETA)->SetWindowTextW((LPCTSTR)s);
	s.Format(L"%3d", 50);			this->GetDlgItem(IDC_EDIT_OPENCV_HOUGH_THRESHOLD)->SetWindowTextW((LPCTSTR)s);

	// 
	numOfSavedImage = 0;
}

void CMy070308nttDlg::ReleaseResources() {	
	// images
	if(srcImage)	cvReleaseImage(&srcImage);	srcImage = NULL;
	if(grayImage)	cvReleaseImage(&grayImage);	grayImage = NULL;
	if(bwImage)		cvReleaseImage(&bwImage);	bwImage = NULL;
	if(boundaryImage)			cvReleaseImage(&boundaryImage);	boundaryImage = NULL;
	if(lineExtractImage)		cvReleaseImage(&lineExtractImage);	lineExtractImage = NULL;

	// bitmapInfo
	if(bitmapInfo) delete[] bitmapInfo;
	bitmapInfo = NULL;
}

// ----- GUI ----- 

void CMy070308nttDlg::ShowImageToPictureControl(HWND picControl, BYTE *imageData, int width, int height, int typeOfImage) {
	if (bitmapInfo) delete[] bitmapInfo;
	bitmapInfo = CreateBMPInfo(width, height, typeOfImage); // 1 = black-white; 3 = RGB

	HDC hDC=::GetDC(picControl);
	::SetStretchBltMode(hDC, COLORONCOLOR);
	::StretchDIBits(hDC, 0, 0, idcPicture_with, idcPicture_height, 0, 0, width, height,
			imageData,					// lpBits
			(LPBITMAPINFO) bitmapInfo,	// lpBitsInfo
			DIB_RGB_COLORS,				// iUsage
			SRCCOPY);

	::ReleaseDC(picControl, hDC);
}

BYTE * CMy070308nttDlg::CreateBMPInfo(int width, int height, int bpp) { // khong hieu
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

void CMy070308nttDlg::UpdateScreen() {
	const int GRAY_CHANNEL = 1;
	const int RGB_CHANNEL = 3;
	if(srcImage) ShowImageToPictureControl(m_hOut1Wnd, (BYTE*)srcImage->imageData, srcImage->width, srcImage->height, RGB_CHANNEL);
	if(bwImage) ShowImageToPictureControl(m_hOut2Wnd, (BYTE*)bwImage->imageData, bwImage->width, bwImage->height, GRAY_CHANNEL);
	if(boundaryImage) ShowImageToPictureControl(m_hOut3Wnd, (BYTE*)boundaryImage->imageData, boundaryImage->width, boundaryImage->height, GRAY_CHANNEL);
	//if(boundaryImage) ShowImageToPictureControl(m_hOut3Wnd, (BYTE*)boundaryImage->imageData, boundaryImage->width, boundaryImage->height, RGB_CHANNEL); // neu muon hien thi colored contours
	if(lineExtractImage) ShowImageToPictureControl(m_hOut4Wnd, (BYTE*)lineExtractImage->imageData, lineExtractImage->width, lineExtractImage->height, RGB_CHANNEL);
}

// ----- image processing -----

void CMy070308nttDlg::ConvertSrcImageToGrayscale() {
	if (srcImage) {
		grayImage = cvCreateImage( cvSize(srcImage->width, srcImage->height), IPL_DEPTH_8U, 1 ); // create an "empty" image
		cvConvertImage(srcImage, grayImage);
	}
}

void CMy070308nttDlg::DoThresholdOnGrayImage(double threshold, double maxValue, int thresholdType) {
	//int threshold = 127, max_value = 255; // param to convert to Black-White image
	if (grayImage) {
		bwImage = cvCreateImage( cvSize(grayImage->width, grayImage->height), IPL_DEPTH_8U, 1 ); // create an "empty" image
		cvThreshold(grayImage, bwImage, threshold, maxValue, thresholdType); 
	}
}
