//
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// ADevScan.h      Initial author: MW, RD, JWP           10-Dec-1995
//  This file contains the abstract base class ADevScan definition which
//     implements a dummy Scan object of d*TREK.
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
#ifndef _ADEVSCAN_H_
#define _ADEVSCAN_H_

#ifdef USE_MSVC_CSTRING
#include "stdafx.h"
#define Cstring CString
#else
#include "Cstring.h"
#endif

#include <time.h>
#include "DevErrCode.h" 
#include "DevEnumeration.h"

#include "CDevDetector.h"
#include "CDevGoniom.h"
#include "CDevShutter.h"
#include "CDevCounter.h"

#ifdef MSVC
#include "CscanMS.h"
#else 
#include "Cscan.h"
#endif

#ifdef MSVC
#define G_nNumScanStates 12
#else
const int G_nNumScanStates = (int)eDevScanState_NumStates;
#endif

#ifdef MSVC
#define G_nNumScanModes 7
#else
const int G_nNumScanModes = (int)eDevScanMode_NumModes;
#endif

// ADevScan class

class ADevScan
{
// Member variables

protected:

    CDevDetector  *m_poDevDetector;
    CDevShutter   *m_poDevShutter;
    CDevGoniom    *m_poDevGoniomSample;
    CDevGoniom    *m_poDevGoniomDetector;
    CDevCounter   *m_poDevCounter;
    double         m_dStartPosition;            //Scan Sample Start Position
    int            m_nLastSeqNumber;
    int            m_nSeqStart;
    int            m_nSeqIncr;
    int            m_nNumImages;               //Number of Images in the Scan
    int            m_nSeqCurr;
    int            m_nScanAxis;
    double         m_dIncrement;               //Scan Image Increment SET
    double         m_dRequestedIncrement;      //Scan Image Increment requested
    double         m_dExposureTime;            //Scan Image Exposure Time set
    double         m_dRequestedExposureTime;   //Scan Image Exposure Time req.
    Cstring        m_sTemplate;                //Image Filename Template
    eDevScan_State m_eState;
    eDevScan_Mode  m_eMode;
    int            m_nDetDatumAxes; 
    int            m_nSampleDatumAxes; 
    int            m_nPauseSeqNumber;
    double        *m_pdDetDatum;
    double        *m_pdSampleDatum;
    Cstring        m_sHeader;
    double         m_dStartTime;
    double         m_dStartPause;
    double         m_dElapsedTime;
    double         m_dCurrentPosition;
    double         m_dEndPosition;
    double         m_dTotalTimeNeeded;
    Cstring        m_sStatusMsg;
    Cstring        m_sDetectorOptions;

    static Cstring ms_asScanStates[G_nNumScanStates];
    static Cstring ms_asScanModes[G_nNumScanModes];

// Functions & methods

public:

    ADevScan (CDevGoniom   *poDevGoniomDetector, 
              CDevGoniom   *poDevGoniomSample, 
              CDevDetector *poDevDetector, 
              CDevShutter  *poDevShutter, 
              CDevCounter  *poDevCounter=NULL);

    virtual ~ADevScan ();

    virtual void vConstruct(void) = 0;
    virtual void vDestruct(void)  = 0;

    inline void vGetStatusMsg(Cstring *psStatus) { *psStatus = m_sStatusMsg; }

    // Get and Set Sample Goniometer Scan Axis & Set Start Position

    virtual int nGetScanAxis (int *pnAxis);
    virtual int nSetScanAxis (const int nAxis, const double dPosition = 0.0);
    virtual int nSetScanAxis (const Cstring& sAxisName, const double dPosition);


    // Get and Set Sample Goniometer Reference Positions

    virtual int nGetDatumSamplePosition (const int nAxis, double *pdPosition);
    virtual int nGetDatumSamplePosition (const Cstring& sAxisName, double *pdPosition);
    virtual int nSetDatumSamplePosition (const int nAxis, const double dPosition);
    virtual int nSetDatumSamplePosition (const Cstring& sAxisName, 
	                                 const double dPosition);


    // Get and Set Detector Goniometer Reference Positions
 
    virtual int        nGetDatumDetPosition (const int nAxis, double *pdPosition);
    virtual int nGetDatumDetPosition (const Cstring& sAxisName, double *pdPosition);
    virtual int nSetDatumDetPosition(const Cstring& sAxisName, const double dPosition);
    virtual int nSetDatumDetPosition (const int nAxis, const double dPosition);


    // Get and Set Sample Goniometer Scan Axis Start Position
    //   relative to Reference Positions

    virtual int nSetStartPosition(const double dStartPosition);        
    virtual int nGetStartPosition(double* pdStartPosition);        


    // Get and Set Number of Images in a Scan

    virtual int nSetNumImages (const int nNumImages);
    inline int nGetNumImages (void) { return (m_nNumImages);}


    // Get and Set Image Exposure Time in a Scan

    virtual int nGetExposureTime (double* pdExposureTime);
    virtual int nSetExposureTime (const double dExposureTime);


    // Get and Set Sample Goniometer Increment between Images in a Scan

    virtual int nGetIncrement (double *pdIncrement);
    virtual int nSetIncrement (const double dIncrement);
    virtual int nSetHeader (const Cstring& sImageHeader);


   //Get and Set Scan Mode select a still/scan with/without shutter open

    virtual int nGetScanMode (eDevScan_Mode* peMode);                
    virtual int nSetScanMode (const eDevScan_Mode eMode=eDevScanMode_ScanOpen);


   //System Related Checks

    virtual int nCheckTemplate(const Cstring& sTemplate);
    virtual int nCheckDiskSpace(const Cstring& sTemplate);


   //Before Starting a SCAN either
   // .....Do a Setup with above parameters or

    virtual int nSetup(void);


   // .....Do a Setup with parameters from  D*TREK Scan Object

    virtual int nSetup(Cscan *poScan);


    // Start a SCAN 

    virtual int nStart(void);


    // Abort a SCAN 

    virtual int nAbort(void);


    // Set & Get Parameters, Setup, and Start All-In-One 

    virtual int nScan (const int nAxis, const double dStartPosition,
                       const int nNumImages, const double dIncrement, 
                       const double dExposureTime, const eDevScan_Mode eMode);

    virtual int nScan (const Cstring& sAxisName, const double dStartPosition,
                       const int nNumImages, const double dIncrement,
                       const double dExposureTime, const eDevScan_Mode eMode);


    // Updates the Scan State & optionally provides current Scan info.

    virtual int nPoll(eDevScan_State *peState=NULL, double *pdPosition=NULL, 
                      double *pdElapsed=NULL, int *pnLastSeqNum=NULL); 


    // Obtains state stored in object as opposed from Control System

    inline eDevScan_State    ePoll(void) {return(m_eState);} 
 

    // Application can keep track of SCAN progress via 
    //....... LastImageFilename written to Disk
    //....... Sequence # of LastImageFilename  written to Disk
    //............where  LastSeqNumber ranges form 0 to nNumImages
    //............and 0 means no image written to disk yet

    virtual int nGetLastImageFilename(Cstring *psLastImageFilename);
    virtual int nGetLastSeqNum(int *pnLastSeqNumber);

    
    // Wait for image to become available

    virtual int nWait(const int nImageNumber=0);

    // Misc.

    virtual int nPauseResume(const int nMode = 0);
    virtual int nDone(void);
    virtual int nInit(void);
};
#endif  //_ADEVSCAN_H_
