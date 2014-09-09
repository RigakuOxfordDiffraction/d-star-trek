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
// ADevCounter.h      Initial author: MW, RD, JWP           10-Dec-1995
//  This file contains the abstract base class ADevCounter definition which
//     implements a dummy Counter device of d*TREK.
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

#ifndef _ADEVCOUNTER_H_
#define _ADEVCOUNTER_H_

#ifdef USE_MSVC_CSTRING
#include "stdafx.h"
#define Cstring CString
#else
#include "Cstring.h"
#endif

#include <time.h>
#include "DevErrCode.h"
#include "DevEnumeration.h"

const int g_nNumCounterStates = (int) eDevCounterState_NumStates;  // Global!

//ADevCounter class
class ADevCounter
{

//Member variables

protected:

    Cstring *m_psChannelName;   // Pointer to an Array of Cstring Objects
    Cstring  m_sName;
    double  *m_pdCounts;        // Pointer to an Array of  Counts
    int     m_nNumChannels;     // Number of Channels on the Counter
    eDevCounter_State m_eState;
    double   m_dStartTime;
    double   m_dElapsedTime;
    double   m_dSetTime;
    double   m_dRequestedTime;
    double   m_dMaxTime;
    double   m_dMinTime;
    Cstring  m_sStatusMsg;

//Functions & methods
public:
    ADevCounter(const Cstring& sName);
    ADevCounter(const char *pcName="");
    virtual ~ADevCounter();

    inline void vGetStatusMsg(Cstring *psStatus) { *psStatus = m_sStatusMsg; }

    virtual void vConstruct(void) = 0;          // Pure virtual function
    virtual void vDestruct(void)  = 0;          // Pure virtual function

    // Obtain the limits of scaler counting times

    inline double dGetMaxTime (void) { return (m_dMaxTime); }
    inline double dGetMinTime (void) { return (m_dMinTime); }


    // Number of counting channels in the scaler..16 at SBC

    inline int nGetNumChannels (void) { return (m_nNumChannels); }


    // Set and Get scaler channel names

    virtual int nSetChannelName (const int nChannelNum, const Cstring& sChannelName);
    virtual int nGetChannelName (const int nChannelNum, Cstring* psName);


    // Get channel counts from one or many channels

    virtual int nGetCount(const Cstring& sChannelName, double *pdCounts);
    virtual int nGetCount(const int nChannelNum, double *pdCounts);
    virtual int nGetCounts(const int nNumChannels, double *pdCounts);


    // Set the total count time (0 if forever), 
    //         i.e support for timed gating only

    virtual int nSetRequestedTime(const double dTime);
    virtual int nGetRequestedTime(double* pdTime);
    virtual int nGetSetTime(double* pdTime);


    // Get the time elapsed since counting started 

    virtual int nGetElapsedTime (double* pdElapsedTime);


    // Obtain the state from object data, not from the hardware device

    inline eDevCounter_State ePoll(void) { return (m_eState); }


    //Update state & obtain other data also

    virtual int nPoll(eDevCounter_State* peState=NULL, const int nDimIn = 0,
                      double *pdElapsedTime=NULL, double *pdCounts=NULL);


    //Issue single commands to scaler

    virtual int nInit(void);
    virtual int nSetup(const double dTime);
    virtual int nWait(void);
    virtual int nStart(void);
    virtual int nStop(void);
    virtual int nDone(void);
    virtual int nAbort(void);
};

#endif   //_ADEVCOUNTER_H_
