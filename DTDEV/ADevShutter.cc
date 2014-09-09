//
// Copyright (c) 1996 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// ADevShutter.cc      Initial author: J.W. Pflugrath           01-Mar-1996
//  This file contains the member functions of BASE class ADevShutter which
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

#include "ADevShutter.h"

//+Code begin

//+Public functions

// Constructors, destructors and assignments

ADevShutter::ADevShutter(const Cstring& sName)
{
//  cout << "ADevShutter:: constructor called!" << endl;
}

ADevShutter::ADevShutter(const char *pcName)
{
//  cout << "ADevShutter:: constructor called!" << endl;
}

ADevShutter::~ADevShutter()
{
//  cout << "ADevShutter:: destructor called!" << endl;
}

int
ADevShutter::nGetTimeToOpen(double *pdTimeToOpen)
{
  *pdTimeToOpen = m_dTimeToOpen;
  return (DEV_SUCCESS);
}

int 
ADevShutter::nGetTimeToClose(double *pdTimeToClose)
{
  *pdTimeToClose = m_dTimeToClose;
  return (DEV_SUCCESS);
}

int
ADevShutter::nGetElapsedTime(double *pdElapsedTime)
{
  return (nPoll(NULL, pdElapsedTime));
}

int
ADevShutter::nGetMinTime(double *pdMinTime)
{
  *pdMinTime = m_dMinTime;
  return (DEV_SUCCESS);
}

int
ADevShutter::nGetMaxTime(double *pdMaxTime)
{
  *pdMaxTime = m_dMaxTime;
  return (DEV_SUCCESS);
}

int
ADevShutter::nSetRequestedTime(const double dTime)
{
  m_dRequestedTime = dTime;
  m_dSetTime       = m_dRequestedTime;
  return (DEV_SUCCESS);
}

int 
ADevShutter::nGetSetTime(double *pdSetTime)
{
  *pdSetTime = m_dSetTime;
  return (DEV_SUCCESS);
}

int
ADevShutter::nOpen(const int nIgnore)
{
  int nStat;
  nStat = 0;
  if (m_eState != eDevShutterState_Unknown)
    {
      if ( (0 == nIgnore) && (m_eState == eDevShutterState_Open) )
	nStat = 1;
      nStat = nSetup(m_dMaxTime);
      if (DEV_SUCCESS == nStat)
	nStat = nStart();
    }
  else
    {
      m_sStatusMsg = "ERROR in ADevShutter::nOpen, invalid state.\n";
      nStat = DEV_INVALIDSTATE;
    }
  return (nStat);
}

int
ADevShutter::nClose(const int nIgnore)
{
  // If nIgnore == 0, report error if shutter already closed
  // Else nIgnore != 0, report error only if shutter did not close

  int nStat;
  nStat = 0;
  if (m_eState != eDevShutterState_Unknown)
    {
      if ( (0 == nIgnore) && (m_eState == eDevShutterState_Closed) )
	{
	  m_sStatusMsg = "ERROR in ADevShutter::nClose, shutter already closed.\n";
	  nStat = DEV_INVALIDSTATE;
	}
      m_eState = eDevShutterState_Closed;
    }
  else
    {
      m_sStatusMsg = "ERROR in ADevShutter::nClose, invalid state.\n";
      nStat = DEV_INVALIDSTATE;
    }
  return (nStat);
}

int
ADevShutter::nInit(void)
{
  m_eState         = eDevShutterState_Closed;
  m_dElapsedTime   = 0.0;
  m_dSetTime       = 60000.0;
  m_dRequestedTime = 0.0;
  return (DEV_SUCCESS);
}

int
ADevShutter::nSetup(const double dTime)
{
  if (m_eState == eDevShutterState_Closed)
    {
      m_dRequestedTime = dTime;
      m_dSetTime       = m_dRequestedTime;
      m_dElapsedTime   = (double) 0.0;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}

int
ADevShutter::nStart(void)
{
  time_t tStartTime;
  if (m_eState == eDevShutterState_Closed)
    {
      (void) time(&tStartTime);
      m_dStartTime   = (double) tStartTime;
      m_eState       = eDevShutterState_Open;
      return (DEV_SUCCESS);
    }
  else
    {
      m_sStatusMsg = "ERROR in ADevShutter::nStart, invalid state.\n";
      return (DEV_INVALIDSTATE);
    }
}

int
ADevShutter::nPoll(eDevShutter_State *peState, double *pdElapsedTime)
{
  time_t tCurrentTime;  

  if (m_eState == eDevShutterState_Open)
    {
      (void) time(&tCurrentTime);
      m_dElapsedTime = (double) tCurrentTime - m_dStartTime;

      if (m_dElapsedTime > m_dSetTime)
	{
	  m_dElapsedTime = m_dSetTime;
	  m_eState       = eDevShutterState_Closed;
	}
    }
  if (NULL != pdElapsedTime)
    {
      *pdElapsedTime = m_dElapsedTime;
    }
  if (NULL != peState)
    {
      *peState = m_eState;
    }
  if (   (eDevShutterState_Unknown == m_eState)
      || (eDevShutterState_TimedOut == m_eState) )
    {
      m_sStatusMsg = "ERROR in ADevShutter::nPoll, state unknown or timeout.\n";
      return (DEV_INVALIDSTATE);
    }
  else
    return (DEV_SUCCESS);
}

int
ADevShutter::nWait(void)
{
  
  int nStat = DEV_SUCCESS;
  while (   (DEV_SUCCESS == nStat)
	 && (   (eDevShutterState_Open    == m_eState) 
	     || (eDevShutterState_Opening == m_eState)
	     || (eDevShutterState_Closing == m_eState) ) )
    nStat = nPoll();             // Poll until closed or error
  return (nStat);
}

int
ADevShutter::nAbort(void)
{
  return (nClose());
}

int
ADevShutter::nDone(void)
{
  m_eState = eDevShutterState_Unknown;
  return (DEV_SUCCESS);
}
