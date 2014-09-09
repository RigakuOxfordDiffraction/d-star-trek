//DevCounter.h:	header file for CDevCounter class
//		The CDevCounter class will facilitate counter control

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
#ifndef _DEVCOUNTER_H_
#define _DEVCOUNTER_H_

#include "Cstring.h"
#include <time.h>

enum eDevCounter_State
{
    eDevCounterState_Unknown,
    eDevCounterState_Counting,
    eDevCounterState_Idle
};

//CDevCounter class
class CDevCounter
{

//Member variables
private:
    Cstring *m_psChannelName;				//Pointer to an Array of Cstring Objects
    Cstring  m_sName;
    double  *m_pdCounts;				//Pointer to an Array of  Counts
    int     m_nNumChannels;				//Number of Channels on the Counter
    eDevCounter_State m_eState;
    double   m_dStartTime;
    double   m_dElapsedTime;
    double   m_dSetTime;
    double   m_dRequestedTime;
    double   m_dMaxTime;
    double   m_dMinTime;

//Functions & methods
public:
    CDevCounter(const Cstring& sName="");
    ~CDevCounter();
    int nSetChannelName (const int nChannelNum, const Cstring& sChannelName);
    Cstring sGetChannelName (const int nChannelNum);
    int nGetNumChannels (void);
    int nGetCount(const Cstring& sChannelName, double *pdCounts);	//Get Counts from a single channel
    int nGetCount(const int nChannelNum, double *pdCounts);		//Get Counts from a single channel
 

    int nGetCounts(const int nNumChannels, double *pdCounts); 
							//Get counts from all channels on a 	
							//single counter
    double dGetElapsedTime (void);
    double dGetMaxTime (void);
    double dGetMinTime (void);
    double dSetTime(const double dTime);		//Returns actual set time
							//Time=0.0 for Untimed
    int nInit (void);
    double dSetup(const double dTime);
    int nStart (void);
    eDevCounter_State ePoll(const int nDim = 0, double *pdElapsedTime=NULL,
			double *pdCounts=NULL);
    eDevCounter_State eWait(void);
    int nStop (void);
    int nDone(void);
    int nAbort(void);
};
#endif	//_DEVCOUNTER_H_
