// DevGoniom.h: header file for CDevGoniom class
//		The CDevGoniom class facilitates single goniometer motor "slew"
//		movements

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

#ifndef __DEVGONIOM_H__
#define __DEVGONIOM_H__

#include "Cstring.h"
#include <time.h>

//Only goniometer mode is SLEWing, although we have added support for user selectable
//motor speed.

enum eDevGoniom_State
{
    eDevGoniomState_Unknown,
    eDevGoniomState_Collided,
    eDevGoniomState_Jammed,
    eDevGoniomState_Moving,
    eDevGoniomState_NotReady,
    eDevGoniomState_Ready,
    eDevGoniomState_Aborted,
    eDevGoniomState_SuccessMove,
    eDevGoniomState_TimedOut
};


// CGoniometer class
class CDevGoniom
{

// Member variables
private:                   
    Cstring *m_psAxisName;
    Cstring  m_sName;				//Goniometer ID
    int      m_nNumAxes;
    double  *m_pdCurrent;
    double  *m_pdSet;
    double  *m_pdRequested;
    double  *m_pdCurrentSpeed;
    double  *m_pdSetSpeed;
    double  *m_pdRequestedSpeed;
    double   m_dStartTime;
    double   m_dElapsedTime;
    double  *m_pdStart;

    eDevGoniom_State m_eState;
    eDevGoniom_State *m_peStateAxis;

// Functions & methods

public:
     CDevGoniom(const Cstring sName ="");
     ~CDevGoniom();
    Cstring sGetAxisName (const int nAxis);		//Returns axis name
    double dGetPosition(const int nAxis);		//Get current position of nAxis motor

    double dGetPosition (const Cstring& sAxisName);	//Returns current position of AxisName motor
						//Using a reference rather than passing 
						//whole object
    double dSetPosition(const int nAxis, const double dPosition);
						//Returns the actual set position
    double dSetPosition(const Cstring& sAxisName, const double dPosition);
						//Returns the actual set position
    double dSetSpeed(const int nAxis, const double dSpeed);
    double dGetSpeed(const int nAxis);
						//Returns the actual set speed
						//Speed in units of deg/s or mm/s
						//Speed=0.0 means move at slew speed
    double dSetSpeed(const Cstring& sAxisName, const double dSpeed);
						//Returns the actual set speed
    int nGetNumAxes(void);			//Returns the number of axes
    Cstring sGetNumAxes(void);			//Returns the number of axes
    int nMove(const int nAxis,const double dPosition);
						//Move 1 axis=nAxis
						//Pass which axis by number and position 
						//to move to
						//Does dSetRequestedPosition, 	
						//nSetup, nStart, and nWait
    int nMove(const Cstring& sAxisName, const double dPosition);
						//Move 1 axis=srtAxisName	
						//Pass which axis by name and position
						//to move to
						//Does dSetRequestedPosition, 	
						//nSetup, nStart, and nWait
    int nInit(void);				//Reads init file
    int nSetup(void);						 	
    int nStart(void);
    eDevGoniom_State ePoll(const int nDim=0, double *pdPositions=NULL);			
						//Returns goniometer state and axes 	
						//current positions
    eDevGoniom_State eWait(void);
    int nAbort(void);
    int nDone(void);				//Does cleanup and disconnect
};

#endif  // __DEVGONIOM_H__
