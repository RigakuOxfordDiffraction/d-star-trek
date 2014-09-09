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
// ADevCounter.cc      Initial author: J.W. Pflugrath           10-Dec-1995
//  This file contains the member functions of abstract base class ADevCounter
//      which implements a dummy Counter device of d*TREK.
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

#include "ADevCounter.h"

//+Code begin

//+Public functions

// Constructors, destructors and assignments

ADevCounter::ADevCounter(const Cstring& sName)
{
  // All construction is done in the derived class
}

ADevCounter::ADevCounter(const char *pcName)
{
  // All construction is done in the derived class
}

ADevCounter::~ADevCounter()
{
  // All destruction is done in the derived class
}

int
ADevCounter::nSetChannelName(const int nChannelNum, const Cstring& sChannelName)
{
  if ( (0 <= nChannelNum) &&  (nChannelNum < m_nNumChannels) )
    {
      m_psChannelName[nChannelNum] = sChannelName;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevCounter::nGetChannelName(const int nChannelNum, Cstring *psName)
{
  if ( (0 <= nChannelNum) && (nChannelNum < m_nNumChannels) )
    {
      *psName = m_psChannelName[nChannelNum];
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevCounter::nGetCount(const Cstring& sChannelName, double *pdCounts)
{
  int nChannel;
  for (nChannel = 0; (nChannel < m_nNumChannels) 
                     && (sChannelName != m_psChannelName[nChannel]); nChannel++)
    ; // Loop to find matching channel name

  return (nGetCount(nChannel, pdCounts));
}

int
ADevCounter::nGetCount(const int nChannelNum, double *pdCounts)
{
  if ( (0 <= nChannelNum) &&  (nChannelNum < m_nNumChannels) )
    {
      *pdCounts = m_pdCounts[nChannelNum];
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevCounter::nGetCounts(const int nNumChannels, double *pdCounts)
{
  if (   (0 >= nNumChannels) || (nNumChannels > m_nNumChannels) 
      || (NULL == pdCounts))
    return (DEV_INVALIDARG);

  for (int nChan = 0; nChan < nNumChannels; nChan++)
    {
      (void) nGetCount(nChan, &pdCounts[nChan]);
    }
  return (DEV_SUCCESS);
}

int
ADevCounter::nSetRequestedTime(const double dTime)
{
  m_dRequestedTime = dTime;
  if (m_dMaxTime > m_dRequestedTime)
    {
      m_dSetTime = m_dMaxTime;
      return (DEV_INVALIDSETPOINT);
    }
  else if (m_dMinTime > m_dRequestedTime)
    {
      m_dSetTime = m_dMinTime;
      return (DEV_INVALIDSETPOINT);
    }
  m_dSetTime       = m_dRequestedTime;  // Not always true on some hardware!
  return (DEV_SUCCESS);
}

int
ADevCounter::nGetRequestedTime(double *pdTime)
{
  *pdTime = m_dRequestedTime;
  return (DEV_SUCCESS);
}

int
ADevCounter::nGetSetTime(double *pdTime)
{
  *pdTime = m_dSetTime;
  return (DEV_SUCCESS);
}

int
ADevCounter::nGetElapsedTime(double *pdElapsedTime)
{
  if (NULL == pdElapsedTime) return (DEV_INVALIDARG);
  return (nPoll(NULL, 0, pdElapsedTime, NULL));
}

int
ADevCounter::nInit(void)
{
  int i;
  for (i = 0; i < m_nNumChannels; i++)
    {
      m_pdCounts[i] = 0.0;
    }
  m_eState = eDevCounterState_Idle;
  return (DEV_SUCCESS);
}

int
ADevCounter::nSetup(const double dTime)
{
  int i;
  if (eDevCounterState_Idle == m_eState)
    {
      m_dRequestedTime = dTime;
      m_dSetTime       = m_dRequestedTime;
      m_dElapsedTime   = (double) 0.0;
      for (i = 0; i < m_nNumChannels; i++)
	{
	  m_pdCounts[i] = 0.0;
	}
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}

int
ADevCounter::nStart(void)
{
  time_t tStartTime;
  if (eDevCounterState_Idle == m_eState)
    {
      (void) time(&tStartTime);
      m_dStartTime   = (double) tStartTime;
      m_eState       = eDevCounterState_Counting;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}

int
ADevCounter::nPoll(eDevCounter_State* peState, const int nDimIn,
		   double *pdElapsedTime, double *pdCounts)

{
  time_t tCurrentTime;  
  int    nDim;
  int    i;
  int    nCounting;

  nDim = nDimIn;
  if (0 > nDim) nDim = 0;
  if (nDim > m_nNumChannels) nDim = m_nNumChannels;

  if (eDevCounterState_Counting == m_eState)
    {
      (void) time(&tCurrentTime);
      m_dElapsedTime = (double) tCurrentTime - m_dStartTime;

      // Calculate current position based on elapsed time and speed

      nCounting = 0;
      if (m_dElapsedTime > m_dSetTime)
	{
	  m_dElapsedTime = m_dSetTime;
	  m_eState       = eDevCounterState_Idle;
	}
      else
	{
	  nCounting++;
	}
      for (i = 0; i < m_nNumChannels; i++)
	{
	  m_pdCounts[i] = m_dElapsedTime * (double) (i + 1);
	}
    }
  if (NULL != pdElapsedTime)
    {
      *pdElapsedTime = m_dElapsedTime;
    }
  if (NULL != pdCounts)
    {
      for (i = 0; i < nDim; i++)
	{
	  pdCounts[i] = m_pdCounts[i];
	}
    }
  if (NULL != peState)
    *peState = m_eState;
  
  if (   (eDevCounterState_Unknown == m_eState)
      || (eDevCounterState_NotAvailable == m_eState) )
    return (DEV_INVALIDSTATE);
  else
    return (DEV_SUCCESS);
}

int
ADevCounter::nWait(void)
{
  
  int nError = DEV_SUCCESS;
  while (   (DEV_SUCCESS == nError) 
	 && (eDevCounterState_Counting == m_eState) )
    nError = nPoll();             // Poll until not counting or error
  return (nError);
}

int
ADevCounter::nAbort(void)
{
  m_eState = eDevCounterState_Idle;  // Maybe unknown?
  return (DEV_SUCCESS);
}

int
ADevCounter::nStop(void)
{
  return (nAbort());
}

int
ADevCounter::nDone(void)
{
  m_eState = eDevCounterState_Unknown;
  return (DEV_SUCCESS);
}

