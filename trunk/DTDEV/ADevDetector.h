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
// ADevDetector.h      Initial author: MW, RD JWP           10-Dec-1995
//  This file contains the ABSTRACT base class BDevDetector definition which
//     implements a dummy Detector device of d*TREK.  Derived subclasses
//     MUST supply nContruct and nDestruct methods!
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
//+Description
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#ifndef _DEVDETECTOR_H_
#define _DEVDETECTOR_H_
                      
#include "Dtrek.h"

#include <time.h>
#ifdef USE_MSVC_CSTRING
#include "stdafx.h"
#define Cstring CString
#else
#include "Cstring.h"
#endif

#include "DevErrCode.h"
#include "DevEnumeration.h"

// ADevDetector class
class ADevDetector
{

// Member variables

protected:                   

  Cstring            m_sName;             // Detector name
  unsigned int       m_unDim0;            // First dimension of detector in pix
  unsigned int       m_unDim1;            // 2nd   dimension of detector in pix
  double             m_dSetTime;          // Time in seconds actually set
  double             m_dRequestedTime;    // Time in seconds requested
  eDevDetector_Mode  m_eMode;             // Mode of detector
  eDevDetector_State m_eState;            // State of detector
  Cstring            m_sHeaderString;     // A header string for ?
  double             m_dElapsedTime;      // Elapsed time since last start
  double             m_dMinTime;          // Min allowed exposure time
  double             m_dMaxTime;          // Max allowed exposure time

  double             m_dStartTime;        // Internal clock time of start
  Cstring            m_sStatusMsg;        // Retrievable status message 
  Cstring            m_sOptions;          // Detector specific options

// Functions & methods

public:       

  ADevDetector(const Cstring &sName);
  ADevDetector(const char *pcName="");
  virtual ~ADevDetector();

  inline void vGetHeader(Cstring *psHeader) { *psHeader = m_sHeaderString; }
  inline void vGetStatusMsg(Cstring *psStatus) { *psStatus = m_sStatusMsg; }

  virtual void vConstruct(void) = 0;          // Pure virtual function
  virtual void vDestruct(void)  = 0;          // Pure virtual function

  virtual int nGetDim0(int *pnDim0);          // First dimension of detector
  virtual int nGetDim1(int *pnDim1);          // Second dimension of detector
  virtual int nSetTime (const double dTime);  // Set requested exposure time
  virtual int nGetTime (double *pdTime);      // Get actual set exposure time

  virtual int nPoll(eDevDetector_State *peState = NULL,// Returns detector state
		    double *pdElapsedTime = NULL);     // elapsed time 

  // Issue single commands to the detector, all return DEV_SUCCESS
  // on successful completion

  virtual int nInit(void);
  virtual int nSetup(const double dTime);
  virtual int nWait(void);
  virtual int nStart(void);
  virtual int nStop(void);
  virtual int nDone(void);
  virtual int nAbort(void);

  virtual int nExpose(const double dTime);    // Does an nSetTime,
                                              // nSetup,nStart,nWait

// Some things that might be better in a derived class:

  virtual int nSetDetectorMode(const eDevDetector_Mode eMode 
			 = eDevDetMode_StillClosed);

  virtual eDevDetector_Mode eGetDetectorMode (void); // Get Detector Mode

  virtual int nSetHeader(Cstring& sImageHeader); // Transfer header to Control System

  virtual int nSetOptions(Cstring &sOptions) // Detector specific options
    { m_sOptions = sOptions; return DEV_SUCCESS; }
  virtual Cstring sGetOptions(void)
    { return m_sOptions; }
};                  

#endif  // _DEVDETECTOR_H_
