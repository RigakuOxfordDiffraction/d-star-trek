#1//DevSource.h:	header file for CDevSource class
//		The CDevSource class will facilitate control of beamline 
//		energy/wavelength.

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

#ifndef _DEVSOURCE_H_
#define _DEVSOURCE_H_

#include "Cstring.h"
#include <time.h>

enum eDevSource_State
{
    eDevSourceState_Unknown,
    eDevSourceState_Available,				
    eDevSourceState_NotAvailable			
};

// CDevSource class
class CDevSource
{
// Member variables
private:                   
    Cstring m_sName;				//Source ID
    double m_dEnergy;				//Beamline Energy
    double m_dWavelength;			//Beamline Wavelength
    double m_dStart;
    double m_dSet;
    double m_dRequested;
    eDevSource_State m_eState;
    int    m_nOptimizeMode;
    double m_dElapsedTime;
    double m_dStartTime;

// Functions & methods
public:
    static double m_sdFactor;
    CDevSource (const Cstring sName = "");
    ~CDevSource ();
    double dGetEnergy (void);			//Get Current Beamline Energy
    double dGetWavelength (void);		//Get Current Beamline Wavelength
    double dSetEnergy (const double dEnergy);	//Set Requested Beamline Energy
						//Returns actual set energy
    double dSetWavelength (const double dWavelength); // Set Requested Beamline Wavelength
						//Returns actual set wavelength
    int nSetOptimize(const int nMode);		//For future use
    int nStart(void); 
    eDevSource_State ePoll(double *pdEnergy=NULL);
    eDevSource_State eWait(double *pdEnergy=NULL);
    int nAbort(void);
    int nInit(void);
    int nSetup(void);
    int nDone(void);
};
#endif	//_DEVSOURCE_H_


