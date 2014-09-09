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
// ADevSource.cc      Initial author: J.W. Pflugrath           10-Dec-1995
//  This file contains the member functions of base class ADevSource which
//     implements a dummy Source device of d*TREK.
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

#include "ADevSource.h"

//+Code begin

//+Public functions

// Constructors, destructors and assignments

double ADevSource::ms_dFactor = 12398.541;

ADevSource::ADevSource(const Cstring &sName)
{
//  cout << "ADevSource constructor called!\n";
}

ADevSource::ADevSource(const char *pcName)
{
//  cout << "ADevSource constructor called!\n";
}

ADevSource::~ADevSource()
{
//  cout << "ADevSource destructor called!\n";
}

int
ADevSource::nGetEnergy(double *pdEnergy)
{
  *pdEnergy = m_dWavelength * ms_dFactor;
  return (DEV_SUCCESS);
}

int
ADevSource::nGetWavelength(double *pdWavelength)
{
  *pdWavelength = m_dWavelength;
  return (DEV_SUCCESS);
}

int
ADevSource::nGetSetEnergy(double *pdEnergy)
{
  *pdEnergy = m_dSet * ms_dFactor;
  return (DEV_SUCCESS);
}

int
ADevSource::nGetSetWavelength(double *pdWavelength)
{
  *pdWavelength = m_dSet;
  return (DEV_SUCCESS);
}

int
ADevSource::nSetRequestedEnergy(const double dEnergy)
{
  m_dSet       = dEnergy / ms_dFactor;
  m_dRequested = m_dSet;
  return (DEV_SUCCESS);
}

int 
ADevSource::nSetRequestedWavelength(const double dWavelength)
{
  if ( (m_dMin > dWavelength) ||  (m_dMax < dWavelength) )
    return (DEV_INVALIDARG);
  m_dSet       = dWavelength;
  m_dRequested = m_dSet;

  return (DEV_SUCCESS);
}

int
ADevSource::nSetOptimize(const int nMode)
{
  m_nOptimizeMode = nMode;
  return (DEV_SUCCESS);
}

int
ADevSource::nInit(void)
{
  m_eState      = eDevSourceState_Available;
  if (m_dWavelength <= 0.0)
    {
      m_dWavelength = 1.54178;
    }
  return (DEV_SUCCESS);
}

int 
ADevSource::nSetup(void)
{
  if (m_eState != eDevSourceState_Unknown)
    {
      m_eState       = eDevSourceState_Available;
      m_dElapsedTime = 0.0;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}

int
ADevSource::nStart(void)
{
  time_t tStartTime;
  if (m_eState == eDevSourceState_Available)
    {
      (void) time(&tStartTime);
      m_dStartTime   = (double) tStartTime;
      m_dElapsedTime = (double) 0.0;
      m_eState       = eDevSourceState_NotAvailable;
      m_dStart       = m_dWavelength;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDMODE);
    }
}

int
ADevSource::nPoll(eDevSource_State *peState, double *pdWavelength)
{

  time_t tCurrentTime;  

  if (m_eState == eDevSourceState_NotAvailable)
    {
      if (m_dWavelength != m_dSet)
	{
	  double dDelta;
	  (void) time(&tCurrentTime);
	  m_dElapsedTime = (double) tCurrentTime - m_dStartTime;
	  dDelta = m_dElapsedTime * 0.5;           // 0.5 Angstrom per second
	  if (m_dSet < m_dStart) dDelta = -dDelta;
	  m_dWavelength = m_dStart + dDelta;
	  if (    ( (dDelta >  0.0) && (m_dWavelength >= m_dSet) )
	      || ( (dDelta <  0.0) && (m_dWavelength <= m_dSet) ) )
	    {
	      m_dWavelength = m_dSet;
	      m_eState       = eDevSourceState_Available;
	    }
	}
      else
	{
	  m_eState = eDevSourceState_Available;
	}
    }
  if (NULL != pdWavelength)
    {
      *pdWavelength = m_dWavelength;
    }
  if (NULL != peState)
    {
      *peState = m_eState;
    }

  return (DEV_SUCCESS);
}

int
ADevSource::nWait(void)
{
  
  int nStat = DEV_SUCCESS;
  while (   (DEV_SUCCESS == nStat)
	 && (eDevSourceState_Available == m_eState) )
    nStat = nPoll();             // Poll until not moving or error
  return (nStat);
}

int
ADevSource::nAbort(void)
{
  m_eState = eDevSourceState_Unknown;
  return (DEV_SUCCESS);
}

int
ADevSource::nDone(void)
{
  m_eState = eDevSourceState_Unknown;
//  cout << "ADevSource nDone called!\n";
  return (DEV_SUCCESS);
}
