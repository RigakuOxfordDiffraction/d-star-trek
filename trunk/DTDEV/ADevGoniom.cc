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
// ADevGoniom.cc      Initial author: J.W. Pflugrath           10-Dec-1995
//  This file contains the member functions of base class ADevGoniom which
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

#include "ADevGoniom.h"

//+Code begin

//+Definitions, constants, and initialization of static member variables

Cstring	ADevGoniom::ms_asGoniomStates[] = {
	"Unknown",
	"NotReady",
	"Ready",
	"Moving",
	"Done",
	"Aborted",
	"TimedOut",
	"Failed",
	"Collided",
	"Jammed"
	};

//+Public functions

// Constructors, destructors and assignments

ADevGoniom::ADevGoniom(const Cstring &sName)
{
//  cout << "ADevGoniom constructor called!\n";
}

ADevGoniom::ADevGoniom(const char *pcName)
{
//  cout << "ADevGoniom constructor called!\n";
}

ADevGoniom::~ADevGoniom()
{
//  cout << "ADevGoniom destructor called!\n";
}

Cstring
ADevGoniom::sGetAxisName(const int nAxis)
{
  if ( (0 <= nAxis) && (nAxis < m_nNumAxes) )
    {
      return (m_psAxisName[nAxis]);
    }
  else
    {
      return (Cstring("Error"));
    }
}

int
ADevGoniom::nSetMode(const eDevGoniom_Mode eMode)
{
  m_eMode = eMode;
  return (DEV_SUCCESS);
}

int
ADevGoniom::nGetMode(eDevGoniom_Mode *peMode)
{
  *peMode = m_eMode;
  return (DEV_SUCCESS);  
}

int 
ADevGoniom::nGetSetPosition(const int nAxis, double* pdPosition)
{
  if ( (0 <= nAxis) && (nAxis < m_nNumAxes) )
    {
      *pdPosition = m_pdSet[nAxis];
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevGoniom::nGetSetPosition(const Cstring& sAxisName, double* pdPosition)
{
  int nAxis;
  for (nAxis = 0; (nAxis < m_nNumAxes) && (sAxisName != m_psAxisName[nAxis]); 
       nAxis++) ; // Find the axis number
  return (nGetSetPosition(nAxis, pdPosition));
}


int
ADevGoniom::nGetPosition(const int nAxis, double *pdPosition)
{
  if ( (0 <= nAxis) && (nAxis < m_nNumAxes) )
    {
      *pdPosition = m_pdCurrent[nAxis];
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevGoniom::nGetPosition(const Cstring& sAxisName, double* pdPosition)
{
  int nAxis;
  for (nAxis = 0; (nAxis < m_nNumAxes) && (sAxisName != m_psAxisName[nAxis]); 
       nAxis++) ; // Find the axis number
  return (nGetPosition(nAxis, pdPosition));
}

int
ADevGoniom::nSetRequestedPosition(const int nAxis, const double dPosition)
{
  if ( (0 <= nAxis) && (nAxis < m_nNumAxes) )
    {
      if ( (dPosition >= m_pdMin[nAxis]) && (dPosition <= m_pdMax[nAxis]) )
	{
	  m_pdRequested[nAxis] = dPosition;

	  // BE CAREFUL, not all goniometers can reach all set positions!
	  // In particular, those with stepping motors cannot normally go to
	  // non-integral steps!

	  m_pdSet[nAxis]       = m_pdRequested[nAxis];
	  return (DEV_SUCCESS);
	}
      else
	{
	  return (DEV_INVALIDARG);
	}
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevGoniom::nSetRequestedPosition(const Cstring& sAxisName,
				  const double dPosition)
{
  int nAxis;
  for (nAxis = 0; (nAxis < m_nNumAxes) && (sAxisName != m_psAxisName[nAxis]); 
       nAxis++) ; // Find the axis number
  return (nSetRequestedPosition(nAxis, dPosition));
}

int
ADevGoniom::nGetRequestedPosition(const int nAxis, double* pdPosition)
{
  if ( (0 <= nAxis) && (nAxis < m_nNumAxes) )
    {
      *pdPosition = m_pdSet[nAxis];
      return (DEV_SUCCESS);
    }
  else
    {
     return (DEV_INVALIDARG);
    }
}

int 
ADevGoniom::nGetRequestedPosition(const Cstring& sAxisName, double* pdPosition)
{
  int nAxis;
  for (nAxis = 0; (nAxis < m_nNumAxes) && (sAxisName != m_psAxisName[nAxis]); 
       nAxis++) ; // Find the axis number
  return (nGetSetPosition(nAxis, pdPosition));
}


int
ADevGoniom::nGetSpeed(const int nAxis, double *pdSpeed)
{
  if ( (0 <= nAxis) && (nAxis < m_nNumAxes) )
    {
      *pdSpeed = m_pdSetSpeed[nAxis];
      return (DEV_SUCCESS);
    }
  else
    {
     return (DEV_INVALIDARG);
    }
}

int
ADevGoniom::nGetSpeed(const Cstring& sAxisName, double *pdSpeed)
{
  int nAxis;
  for (nAxis = 0; (nAxis < m_nNumAxes) && (sAxisName != m_psAxisName[nAxis]); 
       nAxis++) ; // Find the axis number
  return (nGetSpeed(nAxis, pdSpeed));
}

int
ADevGoniom::nSetRequestedSpeed(const int nAxis, const double dSpeed)
{
  double dTemp;
  if ( (0 > nAxis) || (nAxis >= m_nNumAxes) )
    return (DEV_INVALIDARG);

  if (0.0 == dSpeed)
    {
      // Move at max speed if input speed is 0.0
      dTemp = m_pdMaxSpeed[nAxis];
    }
  else
    {
      dTemp = dSpeed;
    }
  if ( (dTemp < m_pdMinSpeed[nAxis]) || (dTemp > m_pdMaxSpeed[nAxis]) )
    return (DEV_INVALIDARG);

  m_pdRequestedSpeed[nAxis] = dTemp;
  m_pdSetSpeed[nAxis]       = m_pdRequestedSpeed[nAxis];
  return (DEV_SUCCESS);
}

int
ADevGoniom::nSetRequestedSpeed(const Cstring& sAxisName, const double dSpeed)
{
  int nAxis = 0;
  for (nAxis = 0; 
        (nAxis < m_nNumAxes) && (sAxisName != m_psAxisName[nAxis]); nAxis++) ;

  return (nSetRequestedSpeed(nAxis, dSpeed));
}

int
ADevGoniom::nMove(const int nAxis, const double dPosition)
{
  int nStat;
  nStat = nSetRequestedPosition(nAxis, dPosition);
  if (DEV_SUCCESS == nStat)
    {
      nStat = nSetup();
      if (DEV_SUCCESS == nStat)
	{
	  nStat = nStart();
	  if (DEV_SUCCESS == nStat)
	    {
	      nStat = nWait();
	    }
	}
    }
  return (nStat);
}

int
ADevGoniom::nMove(const Cstring& sAxisName, const double dPosition)
{
  int nAxis = 0;
  for (nAxis = 0; 
       (nAxis < m_nNumAxes) && (sAxisName != m_psAxisName[nAxis]); nAxis++) ; //

  return (nMove(nAxis, dPosition));
}

int
ADevGoniom::nInit(void)
{
  m_eState = eDevGoniomState_NotReady;
  return (DEV_SUCCESS);
}

int 
ADevGoniom::nSetup(void)
{
  if (   (eDevGoniomState_NotReady == m_eState)
      || (eDevGoniomState_Ready    == m_eState)
      || (eDevGoniomState_Done     == m_eState) )
    {
      m_eState       = eDevGoniomState_Ready;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}

int
ADevGoniom::nStart(void)
{
  int i;
  time_t tStartTime;
  if (eDevGoniomState_Ready == m_eState)
    {
      (void) time(&tStartTime);
      m_dStartTime   = (double) tStartTime;
      m_dElapsedTime = (double) 0.0;
      m_eState       = eDevGoniomState_Moving;
      for (i = 0; i < m_nNumAxes; i++)
	{
	  m_pdStart[i] = m_pdCurrent[i];    // Save start for simulation
	}
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}

int
ADevGoniom::nPoll(eDevGoniom_State *peState, const int nDimIn, double *pdPos)
{
  time_t tCurrentTime;  
  int    nDim;
  int    i;
  int    nNotMoving;
  double dDelta;

  nDim = nDimIn;
  if (nDim < 0) nDim = 0;
  if (nDim > m_nNumAxes) nDim = m_nNumAxes;

  if (eDevGoniomState_Moving == m_eState)
    {
      (void) time(&tCurrentTime);
      m_dElapsedTime = (double) tCurrentTime - m_dStartTime;

      // Calculate current position based on elapsed time and speed

      nNotMoving = 0;                 // Number of motors not moving
      for (i = 0; i < m_nNumAxes; i++)
	{
	  if (m_pdSet[i] != m_pdCurrent[i])
	    {
	      dDelta = m_dElapsedTime * m_pdSetSpeed[i];  // how far it should have moved
     	      if (m_pdSet[i] < m_pdCurrent[i]) dDelta = -dDelta;
	      m_pdCurrent[i] = m_pdStart[i] + dDelta;

	      if (    ( (dDelta < 0.0) && (m_pdCurrent[i] <= m_pdSet[i]) )
		   || ( (dDelta > 0.0) && (m_pdCurrent[i] >= m_pdSet[i]) ) )
		{
		  nNotMoving++;                // Got there, so not moving
		  m_pdCurrent[i] = m_pdSet[i];
		}
	    }
	  else
	    {
	      nNotMoving++;
	    }
	}

      if (nNotMoving == m_nNumAxes)
	{
	  m_eState       = eDevGoniomState_Done;
	}
    }
  if (NULL != pdPos)
    {
      for (i = 0; i < nDim; i++)
	{
	  pdPos[i] = m_pdCurrent[i];
	}
    }
  if (NULL != peState)
    {
      *peState = m_eState;
    }
  if (   (eDevGoniomState_NotReady  == m_eState)
      || (eDevGoniomState_Ready     == m_eState) 
      || (eDevGoniomState_Moving    == m_eState) 
      || (eDevGoniomState_Done      == m_eState) 
      )
    return (DEV_SUCCESS);
  else
    return (DEV_INVALIDSTATE);
}

eDevGoniom_State
ADevGoniom::ePoll(void)
{
  return (m_eState);
}

int
ADevGoniom::nWait(void)
{
  int nStat = DEV_SUCCESS;
  while ( (DEV_SUCCESS == nStat) && (eDevGoniomState_Moving == m_eState) )
    nStat = nPoll();                      // Poll until not moving
  return (nStat);
}

int
ADevGoniom::nAbort(void)
{
  m_eState = eDevGoniomState_Aborted;
  return (DEV_SUCCESS);
}

int
ADevGoniom::nDone(void)
{
  m_eState = eDevGoniomState_Unknown;
  return (DEV_SUCCESS);
}

