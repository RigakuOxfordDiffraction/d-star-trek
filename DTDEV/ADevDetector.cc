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
// ADevDetector.cc      Initial author: J.W. Pflugrath           01-Mar-1995
//  This file contains the member functions of abstract class ADevDetector which
//     implements a dummy Detector device of d*TREK.
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

#include "ADevDetector.h"

//+Code begin

//+Public functions

// Constructors, destructors and assignments

ADevDetector::ADevDetector(const Cstring &sName)
{
//  cout << "ADevDetector constructor called!\n";
}

ADevDetector::ADevDetector(const char *pcName)
{
//  cout << "ADevDetector constructor called!\n";
}

ADevDetector::~ADevDetector()
{
//  cout << "ADevDetector destructor called!\n";
}

int
ADevDetector::nGetDim0(int *pnDim0)
{
  *pnDim0 = (int) m_unDim0;
  return (DEV_SUCCESS);
}

int
ADevDetector::nGetDim1(int *pnDim1)
{
  *pnDim1 = (int) m_unDim1;
  return (DEV_SUCCESS);
}

int
ADevDetector::nSetTime(const double dTime)
{
  m_dRequestedTime = dTime;
  if ( (m_dMinTime <= dTime) && (m_dMaxTime >= dTime) )
    {
      // Within time limits

      m_dSetTime       = dTime;
      return (DEV_SUCCESS);
      
      // Note: There is a possibility that the detector cannot make an
      //       exposure for this time.  For example, some detectors can
      //       only expose for times that are integral seconds.  This
      //       routine should check for that and set the time to the
      //       nearest valid time.  Of course in this case, the calling
      //       routine needs to know that the requested time was not 
      //       actually set.
    }
  else
    {
      // Outside of time limits

      m_dSetTime = m_dMinTime;
      return (DEV_INVALIDARG);
    }
}

int
ADevDetector::nGetTime (double *pdTime)
{
  *pdTime = m_dSetTime;
  return (DEV_SUCCESS);
}

int
ADevDetector::nSetDetectorMode(const eDevDetector_Mode eMode)
{
  m_eMode = eMode;
  return (DEV_SUCCESS);
}

eDevDetector_Mode
ADevDetector::eGetDetectorMode(void)
{
  return (m_eMode);
}

int
ADevDetector::nExpose(const double dTime)
{
  int nStat = nSetup(dTime);
  if (DEV_SUCCESS == nStat)
    {
      nStat = nStart();
      if (DEV_SUCCESS == nStat)
	{
	  nStat = nWait();
	}
    }
  return (nStat);
}

int
ADevDetector::nSetHeader(Cstring& sHeader)
{
  m_sHeaderString = sHeader;

  return (DEV_SUCCESS);
}

int
ADevDetector::nInit(void)
{
  m_eState = eDevDetState_NotReady;
  return (DEV_SUCCESS);
}

int 
ADevDetector::nSetup(const double dTime)
{
  if (   (eDevDetState_NotReady == m_eState)
      || (eDevDetState_Ready == m_eState) 
      || (eDevDetState_Done == m_eState) )
    {
      int nStat = nSetTime(dTime);
      if (DEV_SUCCESS == nStat)
	{
	  m_eState       = eDevDetState_Ready;
	  m_dElapsedTime = 0.0;
	}
      return (nStat);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}

int
ADevDetector::nStart(void)
{
  time_t tStartTime;
  if (m_eState == eDevDetState_Ready)
    {
      (void) time(&tStartTime);
      m_dStartTime   = (double) tStartTime;
      m_dElapsedTime = (double) 0.0;
      m_eState       = eDevDetState_Exposing;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}

int
ADevDetector::nPoll(eDevDetector_State *peState, double *pdElapsedTime)
{

  time_t tCurrentTime;  

  if (m_eState == eDevDetState_Exposing)
    {
      (void) time(&tCurrentTime);
      m_dElapsedTime = (double) tCurrentTime - m_dStartTime;
      if (m_dElapsedTime >=  m_dSetTime)
	{
	  m_dElapsedTime = m_dSetTime;
	  m_eState       = eDevDetState_Done;
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

  if (   (eDevDetState_Failed  != m_eState)
      || (eDevDetState_Unknown != m_eState) )
    return (DEV_SUCCESS);
  else
    return (DEV_INVALIDMODE);
}

int
ADevDetector::nWait(void)
{
  
  int nStat = DEV_SUCCESS;
  while ( (DEV_SUCCESS == nStat) && (eDevDetState_Exposing == m_eState) )
    nStat = nPoll();                      // Poll until finished
  return (nStat);
}

int
ADevDetector::nAbort(void)
{
  m_eState = eDevDetState_Failed;
  return (DEV_SUCCESS);
}

int
ADevDetector::nDone(void)
{
  m_eState = eDevDetState_Unknown;
  return (DEV_SUCCESS);
}

int
ADevDetector::nStop(void)
{
  m_eState = eDevDetState_Unknown;
  return (DEV_SUCCESS);
}


