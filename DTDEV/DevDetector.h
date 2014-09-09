//DevDetector.h:	header file for CDevDetector class
//			Modified MSC detector.h for use at the SBC
//			The CDevDetector class will facilitate a STILL, a single data
//			measurement with or without shutter control and without goniometer 
//			control.  For multiple data measurements use the CDevScan class.
//			Also, supports a scan of 1 image.
//
//			Regarding Image Filename:  The image filename template and 
//			username must be present in d*TREK image header.
//			If a single image, the sequence number is 0.  If a scan, the sequence 
//			number is determined by hardware.  
//			See Section How an SBC image filename is formed.

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

#ifndef __DEVDETECTOR_H__
#define __DEVDETECTOR_H__
                      
#include "Cstring.h"
#include <time.h>

enum eDevDetector_Mode
{
    eDevDetMode_StillClosed,
    eDevDetMode_StillOpen,
    eDevDetMode_ScanClosed,					
    eDevDetMode_ScanOpen					
};

enum eDevDetector_State
{ 
    eDevDetState_Unknown,
    eDevDetState_NotReady,
    eDevDetState_Ready,
    eDevDetState_Exposing, 
    eDevDetState_Reading,
    eDevDetState_Storing,
    eDevDetState_Done,
    eDevDetState_Failed
 };

// CDevDetector class
class CDevDetector
{

// Member variables
private:                   

  Cstring	     m_sName;					// Detector ID
  unsigned int       m_unDim0;
  unsigned int       m_unDim1;
  double             m_dSetTime;
  eDevDetector_Mode  m_eMode;
  eDevDetector_State m_eState;
  Cstring            m_sHeaderString;
  int                m_nAxis;
  double             m_dInc;
  double             m_dElapsedTime;
  double             m_dStartTime;
  
// Functions & methods
public:       
    CDevDetector(const Cstring sName ="");
    ~CDevDetector();
    unsigned int unGetDim0(void);		//Return size of Detector
    unsigned int unGetDim1(void);		//Return size of Detector
    double dSetTime (double dTime);		//Return Actual Exposure Time Set
    int nSetDetectorMode(const eDevDetector_Mode eMode = eDevDetMode_StillClosed);
    eDevDetector_Mode eGetDetectorMode (void);	//Get Detector Mode
    int nExpose(double dTime);			//Does an nSetTime,
						//nSetup,nStart,nWait
    int nSetHeader(Cstring& sImageHeader);	//Transfer header to Control System
    double dSetScanAxis (int nAxis, double dIncrement);	//Returns actual set increment
    double dSetScanAxis (Cstring& sAxisName, double dIncrement);
						//Returns actual set increment
    int nInit(void);
    int nSetup(void);
    int nStart(void);
    eDevDetector_State ePoll(double *pdElapsedTime=NULL);
						//Returns detector state
						//and elapsed time if not NULL
    eDevDetector_State eWait(void);
    int nAbort(void);
    int nDone(void);
};                  

#endif  // __DEVDETECTOR_H__
