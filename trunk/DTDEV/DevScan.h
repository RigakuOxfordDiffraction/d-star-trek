//DevScan.h:	header file for CDevScan class
//		Modified MSC det_rot.h for use at the SBC
//		The CDevScan class facilitates a series of independent data measurements 
//		with or without shutter control and with or without goniometer control
//		such as, multiple "dark" images and data scans
//		11/14/95  Need to add "dark" info
//
//		Regarding image filenames:  Present in d*TREK image header are
//		SCAN_TEMPLATE, SCAN_SEQ_INFO, and USERNAME.  
//		See Section “How an SBC image filename is formed” for more
//		information.
/*
 *
 * Copyright (C) 2014 Rigaku Americas Corporation
 *                    9009 New Trails Drive
 *                    The Woodlands, TX, USA  77381
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice(s), this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice(s), this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of the Rigaku Americas Corporation nor the
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL RIGAKU AMERICAS CORPORATION BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA OR PROFITS; OR BUSINESS INTERUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
 */

#ifndef _DEVSCAN_H_
#define _DEVSCAN_H_

#include <time.h>
#include "Cstring.h"

#include "DevDetector.h"
#include "DevGoniom.h"
#include "DevShutter.h"
#include "DevCounter.h"
#include "Cscan.h"

enum eDevScan_State
{
    eDevScanState_Unknown,
    eDevScanState_NotReady,
    eDevScanState_Ready,
    eDevScanState_InProgress,
    eDevScanState_InProgressPaused,    
    eDevScanState_ScanDone,	        //3 ways to determine ImageDone:
					//1)  Ask for last filename
					//2)  Get last sequence number
					//3)  Check disk to see if file exists
    eDevScanState_Aborted
};

enum eDevScan_Mode
{
    eDevScanMode_Unknown,			//An unknown mode
    eDevScanMode_StillClosed,			//For multiple DARK images
    eDevScanMode_StillOpen,			//For Multiple STILLS
    eDevScanMode_ScanClosed,			//For Testing Purposes only
    eDevScanMode_ScanOpen			//For X-ray Diffraction Data 
};

// CDevScan class
class CDevScan
{
// Member variables
private:  
    CDevDetector  *m_poDetector;
    CDevShutter   *m_poShutter;
    CDevGoniom    *m_poSampleGoniometer;
    CDevGoniom    *m_poDetectorGoniometer;
    CDevCounter   *m_poCounter;
    double         m_dStartPosition;		//Scan Sample Start Position
    int            m_nSeqStart;
    int            m_nSeqIncr;
    int            m_nNumImages;		//Number of Images in the Scan
    int            m_nSeqCurr;
    int            m_nScanAxis;
    double         m_dIncrement;		//Scan Image Increment
    double         m_dExposureTime;		//Scan Image Exposure Time
    Cstring        m_sTemplate;			//Image Filename Template
    eDevScan_State m_eState;
    eDevScan_Mode  m_eMode;
    int            m_nDetDatumAxes; 
    int            m_nSampleDatumAxes; 
    double        *m_pdDetDatum;
    double        *m_pdSampleDatum;
    Cstring        m_sHeader;
    double         m_dStartTime;
    double         m_dElapsedTime;
    double         m_dCurrentPosition;
    double         m_dEndPosition;
    double         m_dTotalTimeNeeded;

// Functions & methods
public:
    CDevScan (CDevGoniom   *poDetectorGoniometer, 
	      CDevGoniom   *poSampleGoniometer, 
	      CDevDetector *poDetector, 
	      CDevShutter  *poShutter, 
	      CDevCounter  *poCounter=NULL);
   ~CDevScan ();
    int    nGetScanAxis (void);		      //Get Scan Axis
    double dGetDatumDetPosition (const int nAxis);
    double dGetDatumDetPosition (const Cstring& sAxisName);
    double dGetDatumSamplePosition (const int nAxis);
    double dGetDatumSamplePosition (const Cstring& sAxisName);
    double dSetStartPosition(const double dStartPosition);	
                                              //Set Relative Scan Start Position
					      //Returns Scan StartPosition
    int    nSetNumImages (const int nNumImages); //Set Scan Number of Images
					      //Returns Scan Number of Images
inline int nGetNumImages (void) { return (m_nNumImages);} 
					      //Returns Scan Number of Images
    double dSetExposureTime (const double dExposureTime);	
					      //Set Scan Image Exposure Time
					      //Returns Scan Image Exposure Time
    double dSetIncrement (const double dIncrement);
    double dGetIncrement (void);
    int    nSetHeader(const Cstring& sImageHeader);
    double dSetScanAxis (const int nAxis, const double dPosition);
    double dSetScanAxis (const Cstring& sAxisName, const double dPosition);
    double dSetDatumDetPosition (const int nAxis, const double dPosition);
    double dSetDatumDetPosition(const Cstring& sAxisName, const double dPosition);
    double dSetDatumSamplePosition (const int nAxis, const double dPosition);
    double dSetDatumSamplePosition (const Cstring& sAxisName, const double dPosition);
    int    nCheckTemplate(const Cstring& sTemplate);
    int    nCheckDiskSpace(const Cstring& sTemplate);

    int    nScan (const int nAxis, const double dStartPosition,
		  const int nNumImages, const double dIncrement, 
		  const double dExposureTime, const eDevScan_Mode eMode);
    int    nScan (const Cstring& sAxisName, const double dStartPosition,
		  const int nNumImages, const double dIncrement,
		  const double dExposureTime, const eDevScan_Mode eMode);
    Cstring sGetLastImageFilename(void);
    int    nGetLastSeqNum(void);
    eDevScan_Mode eGetScanMode (void);		//Get Scan Mode
    int    nSetScanMode (const eDevScan_Mode eMode=eDevScanMode_ScanOpen);
						//Default Scan Mode=Diffraction Data
    int nSetup(void);					
    int nSetup(Cscan *poScan);
    int nStart(void);
    eDevScan_State ePoll(double *pdPosition=NULL, double *pdElapsed=NULL,
			 int *pnLastSeqNum=NULL); //Returns Scan State and axial
    eDevScan_State eWait(void);
    int nAbort(void);
    int nDone(void);
    int nInit(void);
};
#endif 	//_DEVSCAN_H_
