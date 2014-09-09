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
// ADevSource.h      Initial author: MW, RD, JWP           10-Dec-1995
//  This file contains the base class ADevSource definition which
//     implements a dummy Source device of d*TREK.
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


#ifndef _ADEVSOURCE_H_
#define _ADevSOURCE_H_

#ifdef USE_MSVC_CSTRING
#include "stdafx.h"
#define Cstring CString
#else
#include "Cstring.h"
#endif

#include <time.h>
#include "DevErrCode.h"
#include "DevEnumeration.h"

// ADevSource class
class ADevSource
{
// Member variables

 protected:

    Cstring m_sName;                       // Source ID
    double m_dEnergy;                      // Beamline Energy
    double m_dWavelength;                  // Beamline Wavelength
    double m_dStart;
    double m_dSet;
    double m_dRequested;
    double m_dMin;
    double m_dMax;
    eDevSource_State m_eState;
    int    m_nOptimizeMode;
    double m_dElapsedTime;
    double m_dStartTime;
    Cstring m_sStatusMsg;

// Functions & methods

public:

    ADevSource (const Cstring &sName);
    ADevSource (const char *pcName = "");
    virtual ~ADevSource ();
    static double ms_dFactor;

    virtual void vConstruct(void) = 0;    // Pure virtual function
    virtual void vDestruct(void)  = 0;    // Pure virtual function

    inline void vGetStatusMsg(Cstring *psStatus) { *psStatus = m_sStatusMsg; }

    virtual int nGetEnergy (double *pdEnergy);          // Get Current Beamline Energy
    virtual int nGetWavelength (double *pdWavelength);  // Get Current Beamline Wavelength
    virtual int nGetSetEnergy (double *pdEnergy);       // Get set Beamline Energy
    virtual int nGetSetWavelength (double *pdWavelength);  // Get set Beamline Wavelength
    virtual int nSetRequestedEnergy (const double dEnergy);      // Set Requested Beamline Energy

    virtual int nSetRequestedWavelength (const double dWavelength); // Set Requested Beamline Wavelength

    virtual int nSetOptimize(const int nMode);          // For future use
    virtual int nStart(void); 
    virtual int nPoll(eDevSource_State *peState=NULL, double *pdEnergy=NULL);
    virtual int nWait(void);
    virtual int nAbort(void);
    virtual int nInit(void);
    virtual int nSetup(void);
    virtual int nDone(void);
};
#endif  //_ADEVSOURCE_H_


