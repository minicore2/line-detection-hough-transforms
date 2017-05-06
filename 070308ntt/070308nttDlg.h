// 070308nttDlg.h : header file
#pragma once
#include "cv.h"
#include "highgui.h"
#include "afxwin.h"
#include <math.h>

// CMy070308nttDlg dialog
class CMy070308nttDlg : public CDialog
{
// Construction
public:
	CMy070308nttDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_MY070308NTT_DIALOG };

	// --- ntt ---

	void Init(); 
	void ReleaseResources();

	void ShowImageToPictureControl(HWND picControl, BYTE *imageData, int width, int height, int typeOfImage); 
	BYTE * CreateBMPInfo(int width, int height, int bpp);
	void UpdateScreen(); 

	void ConvertSrcImageToGrayscale(); 	
	void DoThresholdOnGrayImage(double threshold, double maxValue, int thresholdType);

	// ----------

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support	

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonLoadAnImage();	
public:
	afx_msg void OnBnClickedButtonQuit();
	afx_msg void OnBnClickedButtonSaveResultImage();
public:
	afx_msg void OnBnClickedButtonOpencvHoughline2();
public:
	afx_msg void OnBnClickedButtonNttStandardHough();
public:
	afx_msg void OnBnClickedButtonNttStandardHoughNew();
public:
	afx_msg void OnBnClickedButtonNttProbabilisticHough();
public:
	afx_msg void OnBnClickedButtonNttProgressiveProbabilisticHough();
public:
	afx_msg void OnBnClickedButtonNttHoughAccuracy();
public:
	afx_msg void OnBnClickedButtonNttHoughGreen();
public:
	afx_msg void OnBnClickedButtonNttPyramidHough();	
public:
	afx_msg void OnBnClickedButtonLastDecision();	
};
