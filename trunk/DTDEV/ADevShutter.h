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
// ADevShutter.h      Initial author: MW, RD, JWP           10-Dec-1995
//  This file contains the BASE class ADevShutter definition which
//     implements a dummy Shutter device of d*TREK.
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

#ifndef _ADEVSHUTTER_H_
#define _ADEVSHUTTER_H_

#ifdef USE_MSVC_SCSTRING
#include "stdafx.h"
#define   Cstring CString
#else
#include "Cstring.h"
#endif
#include <time.h>
#include "DevErrCode.h"
#include "DevEnumeration.h"

// ADevShutter class
class ADevShutter
{
// Member variables

protected:

    Cstring  m_sName;                           //Shutter ID
    double   m_dOpenTime;                       //Time that shutter is opened

    eDevShutter_State m_eState;
    double   m_dStartTime;
    double   m_dElapsedTime;
    double   m_dSetTime;
    double   m_dRequestedTime;
    double   m_dMaxTime;
    double   m_dMinTime;
    double   m_dTimeToOpen;
    double   m_dTimeToClose;
    Cstring  m_sStatusMsg;

// Functions & methods

public:

    ADevShutter (const Cstring& sName);
    ADevShutter (const char *pcName = "");
    virtual ~ADevShutter ();

    inline void vGetStatusMsg(Cstring *psStatus) { *psStatus = m_sStatusMsg; }

    virtual void vConstruct(void) = 0;   // Pure virtual function
    virtual void vDestruct(void)  = 0;   // Pure virtual function

    virtual int nGetTimeToOpen (double *pdTimeToOpen);
    virtual int nGetTimeToClose (double *pdTimeToClose);
    virtual int nGetElapsedTime (double *pdElapsedTime); // Use UNIX timing, not from hardware
                                                 // Elapsed time=Present-time -Open-time
    virtual int nGetMaxTime(double *pdMaxTime);  // Get value from file, later from db
    virtual int nGetMinTime(double *pdMinTime);  // Get value from file, later from db
    virtual int nSetRequestedTime(const double dTime);   // dTime=0.0 for untimed
    virtual int nGetSetTime( double *pdSetTime);
    virtual int nOpen (const int IgnoreState = 0); // Open shutter
                                                   // If IgnoreState=0, If already in requested
                                                   //     state, generate an error
                                                   // If IgnoreState=1, do not generate an error
    virtual int nClose (const int IgnoreState = 0); // Close shutter-see IgnoreState
                                                    //    comments above
    virtual int nInit (void);                       // Get values from file, later from db
    virtual int nStart(void);
    virtual int nSetup(const double dTime);
    virtual int nPoll (eDevShutter_State *peState = NULL, double *pdElapsedTime=NULL);
    virtual int nWait(void); 
    virtual int nAbort(void);
    virtual int nDone(void);
};
#endif  //_ADEVSHUTTER_H_


