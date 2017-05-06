// 070308ntt.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols


// CMy070308nttApp:
// See 070308ntt.cpp for the implementation of this class
//

class CMy070308nttApp : public CWinApp
{
public:
	CMy070308nttApp();

// Overrides
	public:
	virtual BOOL InitInstance();

// Implementation

	DECLARE_MESSAGE_MAP()
};

extern CMy070308nttApp theApp;