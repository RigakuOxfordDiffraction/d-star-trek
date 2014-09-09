//DevShutter.h:	header file for CDevShutter class
//			Modified from MSC shutter.h for use at the SBC
//			The CDevShutter class will facilitate shutter control

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

#ifndef _DEVSHUTTER_H_
#define _DEVSHUTTER_H_

#include "Cstring.h"
#include <time.h>

enum eDevShutter_State
{
    eDevShutterState_Unknown,
    eDevShutterState_Open,
    eDevShutterState_Closed,
    eDevShutterState_Opening,
    eDevShutterState_Closing,
    eDevShutterState_TimedOut
};
// CDevShutter class
class CDevShutter
{
// Member variables
private:                   
    Cstring  m_sName;				//Shutter ID
    double   m_dOpenTime;			//Time that shutter is opened

    eDevShutter_State m_eState;
    double   m_dStartTime;
    double   m_dElapsedTime;
    double   m_dSetTime;
    double   m_dRequestedTime;
    double   m_dMaxTime;
    double   m_dMinTime;
    double   m_dTimeToOpen;
    double   m_dTimeToClose;

// Functions & methods
public:
    CDevShutter (const Cstring& sName = "");
    ~CDevShutter ();
    double dGetTimeToOpen (void);
    double dGetTimeToClose (void);
    double dGetElapsedTime (void); 		//Use UNIX timing, not from hardware
						//Elapsed time=Present-time -Open-time
    double dGetMaxTime(void);			//Get value from file, later from db
    double dGetMinTime(void);			//Get value from file, later from db
    double dSetTime(const double dTime);	//Returns actual set time	
						//Time=0.0 for untimed
    int nOpen (const int IgnoreState = 0);	//Open shutter
						//If IgnoreState=0, If already in requested
						//state, generate an error
						//If IgnoreState=1, do not generate an error
    int nClose (const int IgnoreState = 0);	//Close shutter-see IgnoreState 	
						//comments above
    int nInit (void);				//Get values from file, later from db
    int nStart(void);
    double dSetup(const double dTime);
    eDevShutter_State ePoll (double *pdElapsedTime=NULL); //Return Shutter state
    eDevShutter_State eWait(void); 
    int nAbort(void);
    int nDone(void);
};
#endif	//_DEVSHUTTER_H_


