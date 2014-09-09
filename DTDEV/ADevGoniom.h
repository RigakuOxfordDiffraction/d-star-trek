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
// ADevGoniom.h      Initial author: MW, RD, JWP           10-Dec-1995
//  This file contains the base class ADevGoniom definition which
//     implements a goniometer device of d*TREK.
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

#ifndef _ADEVGONIOM_H_
#define _ADEVGONIOM_H_

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

static const int g_nNumStates = (int)eDevGoniomState_NumStates;  // Global!

// ADevGoniom class

class ADevGoniom
{

// Member variables

protected:                   

    Cstring *m_psAxisName;
    Cstring  m_sName;                                  // Goniometer Name
    int      m_nNumAxes;
    int      m_nScanAxis;
    double  *m_pdCurrent;
    double  *m_pdSet;
    double  *m_pdRequested;
    double  *m_pdCurrentSpeed;
    double  *m_pdSetSpeed;
    double  *m_pdRequestedSpeed;
    double   m_dStartTime;
    double   m_dElapsedTime;
    double  *m_pdStart;

    double  *m_pdMin;
    double  *m_pdMax;

    double  *m_pdRequestedAccel;
    double  *m_pdSetAccel;
    double  *m_pdCurrentAccel;

    double  *m_pdMinSpeed;
    double  *m_pdMaxSpeed;
    double  *m_pdMinAccel;
    double  *m_pdMaxAccel;

    eDevGoniom_State  m_eState;
    eDevGoniom_State *m_peStateAxis;
    eDevGoniom_Mode   m_eMode;

    Cstring           m_sStatusMsg;

    static Cstring     ms_asGoniomStates[g_nNumStates];

// Functions & methods

public:

     ADevGoniom(const Cstring &sName);
     ADevGoniom(const char *pcName="");
     virtual ~ADevGoniom();

     virtual void vConstruct(void)  = 0;  // Pure virtual function
     virtual void vDestruct(void)   = 0;  // Pure virtual function

     inline void vGetStatusMsg(Cstring *psStatus) { *psStatus = m_sStatusMsg; }

    // Returns goniometer name from data member

    inline Cstring sGetName(void) { return(m_sName); }


    // Gets goniometer axis name given number from member data

    Cstring sGetAxisName (const int nAxis);


    // Returns number of axis from data member

    inline int nGetNumAxes() { return (m_nNumAxes); }


    // Set & Get mode

    virtual int nSetMode(const eDevGoniom_Mode eMode);
    virtual int nGetMode(eDevGoniom_Mode *peMode);
      
    // Get actual (current) axis position

    virtual int nGetPosition (const int nAxis, double* pdPosition);	
    virtual int nGetPosition (const Cstring& sAxisName, double* pdPosition);


    // Set & Get position which axis is REQUESTED to move to	

    virtual int nSetRequestedPosition(const int nAxis, const double dPosition);
    virtual int nSetRequestedPosition(const Cstring& sAxisName, const double dPosition);
    virtual int nGetRequestedPosition(const int nAxis, double* pdPosition);
    virtual int nGetRequestedPosition(const Cstring& sAxisName, double* pdPosition);


    // Get position that the axis is actual SET to move to

    virtual int nGetSetPosition(const int nAxis, double* pdPosition);
    virtual int nGetSetPosition(const Cstring& sAxisName, double* pdPosition);


    // Set & Get axis speed

    virtual int nSetRequestedSpeed(const int nAxis, const double dSpeed);
    virtual int nSetRequestedSpeed(const Cstring& sAxisName, const double dSpeed);
    virtual int nGetSpeed(const int nAxis, double *pdSpeed);
    virtual int nGetSpeed(const Cstring& sAxisName, double *pdSpeed);


    // Setup the move in hardware interface

    virtual int nSetup(void);


    // Start the axis moving						 	

    virtual int nStart(void);


    // Waits for move to end or a timeout 

    virtual int nWait(void);


    // Abort all moving axes

    virtual int nAbort(void);


    // MOVE one axis in slew mode
    // ..Set Position, Setup, Start, & Wait ALL_IN_ONE

    virtual int nMove(const int nAxis, const double dPosition);
    virtual int nMove(const Cstring& sAxisName, const double dPosition);


    //Obtains the current hardware state & stores internally
    //..used only in object construction

    virtual int nInit(void);

    // Obtains hardware state & stores internally
    // ... optionally Gets positions of all axes

    virtual int nPoll(eDevGoniom_State *peState=NULL,
		      const int nDim=0, double *pdPositions=NULL);

    // Returns the stored state in object

    virtual eDevGoniom_State ePoll(void);


    // Does cleanup and disconnect
    //  ...used only in object destruction			

    virtual int nDone(void);				
};

#endif  // _ADEVGONIOM_H_
