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
// ADevScan.cc      Initial author: J.W. Pflugrath           10-Dec-1995
//  This file contains the member functions of abstract base class ADevScan
//      which implements a dummy Scan device of d*TREK.
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

#include "ADevScan.h"

//+Code begin

//+Public functions

// Constructors, destructors and assignments

ADevScan::ADevScan(CDevGoniom   *poDevGoniomDetector, 
                   CDevGoniom   *poDevGoniomSample, 
                   CDevDetector *poDevDetector, 
                   CDevShutter  *poDevShutter, 
                   CDevCounter  *poDevCounter)
{
  // Derived classes do all the construction
}

ADevScan::~ADevScan()
{
  // Derived classes do all the destruction
}

int
ADevScan::nGetScanAxis(int *pnAxis)
{
  *pnAxis = m_nScanAxis;
  return (DEV_SUCCESS);
}

int
ADevScan::nGetDatumDetPosition (const int nAxis, double *pdPosition)
{
  if ( (0 <= nAxis) && (nAxis < m_nDetDatumAxes) )
    {
      *pdPosition = m_pdDetDatum[nAxis];
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevScan::nGetDatumDetPosition (const Cstring& sAxisName, double *pdPosition)
{
  if (NULL != m_poDevGoniomDetector)
    {
      int nAxis;
      for (nAxis = 0;    (nAxis < m_nDetDatumAxes)
                      && (sAxisName 
                          != m_poDevGoniomDetector->sGetAxisName(nAxis));
           nAxis++) ;
      
      return (nGetDatumDetPosition(nAxis, pdPosition));
    }
  else
    {
      return (DEV_NOTCONNECTED);
    }
}

int
ADevScan::nGetDatumSamplePosition (const int nAxis, double *pdPosition)
{
  if ( (0 <= nAxis) && (nAxis < m_nSampleDatumAxes) )
    {
      *pdPosition = m_pdSampleDatum[nAxis];
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevScan::nGetDatumSamplePosition (const Cstring& sAxisName, double *pdPosition)
{
  if (NULL != m_poDevGoniomSample)
    {
      int nAxis;
      for (nAxis = 0;    (nAxis < m_nSampleDatumAxes)
                      && (sAxisName 
                          != m_poDevGoniomSample->sGetAxisName(nAxis));
           nAxis++) ;
      
      return (nGetDatumSamplePosition(nAxis, pdPosition));
    }
  else
    {
      return (DEV_NOTCONNECTED);
    }
}

int 
ADevScan::nGetStartPosition(double *pdStartPosition)
{
  *pdStartPosition = m_dStartPosition;
  return (DEV_SUCCESS);
}

int
ADevScan::nSetStartPosition(const double dStartPosition)
{
  // Need error checking that we can actually get to this start position!
  // Which axis is it?  It is m_nScanAxis!

  m_dStartPosition = dStartPosition;
  return (DEV_SUCCESS);
}

int
ADevScan::nSetNumImages (const int nNumImages)
{
  m_nNumImages = nNumImages;
  return (DEV_SUCCESS);
}

int 
ADevScan::nGetExposureTime(double *pdExposureTime)
{
  *pdExposureTime = m_dStartPosition;
  return (DEV_SUCCESS);
}

int
ADevScan::nSetExposureTime (const double dExposureTime)
{
  // Need check that exposure time is within detector time limits!
  // Some detectors cannot integrate for all possible exposure times, so
  // we need to maintain the time actually set, too!

  m_dExposureTime          = dExposureTime;
  m_dRequestedExposureTime = dExposureTime;
  
  return (DEV_SUCCESS);
}

int
ADevScan::nSetIncrement (const double dIncrement)
{
  // Need check that increment is something goniometer can do.
  // Some goniometers cannot move for all possible increments, so
  // we need to maintain the time actually set, too!

  m_dIncrement          = dIncrement;
  m_dRequestedIncrement = dIncrement;
  return (DEV_SUCCESS);
}

int
ADevScan::nGetIncrement (double *pdIncrement)
{
  *pdIncrement = m_dIncrement;
  return (DEV_SUCCESS);
}

int
ADevScan::nSetHeader(const Cstring& sImageHeader)
{
  m_sHeader = sImageHeader;
  return (DEV_SUCCESS);
}

int
ADevScan::nSetScanAxis (const int nAxis, const double dPosition)
{
  if ( (0 <= nAxis) && (nAxis < m_poDevGoniomSample->nGetNumAxes()))
    {
      m_nScanAxis = nAxis;
      return (nSetStartPosition(dPosition));
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevScan::nSetScanAxis (const Cstring& sAxisName, const double dPosition)
{
  int nAxis;
  for (nAxis = 0;    (nAxis < m_nSampleDatumAxes)
       && (sAxisName 
           != m_poDevGoniomSample->sGetAxisName(nAxis));
       nAxis++) ;
  return (nSetScanAxis(nAxis, dPosition));
}

int
ADevScan::nSetDatumDetPosition (const int nAxis, const double dPosition)
{
  // Need more error checking that the Datum position is actually reachable

  if ( (0 <= nAxis) && (nAxis < m_nDetDatumAxes) )
    {
      m_pdDetDatum[nAxis] = dPosition;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevScan::nSetDatumDetPosition(const Cstring& sAxisName, const double dPosition)
{
  int nAxis;
  for (nAxis = 0;    (nAxis < m_nDetDatumAxes)
       && (sAxisName 
           != m_poDevGoniomDetector->sGetAxisName(nAxis));
       nAxis++) ;
  return (nSetDatumDetPosition(nAxis, dPosition));
}

int
ADevScan::nSetDatumSamplePosition (const int nAxis, const double dPosition)
{
  if ( (0 <= nAxis) && (nAxis < m_nSampleDatumAxes) )
    {
      m_pdSampleDatum[nAxis] = dPosition;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDARG);
    }
}

int
ADevScan::nSetDatumSamplePosition(const Cstring& sAxisName,
                                  const double dPosition)
{
  int nAxis;
  for (nAxis = 0;    (nAxis < m_nSampleDatumAxes)
       && (sAxisName 
           != m_poDevGoniomSample->sGetAxisName(nAxis));
       nAxis++) ;
  return (nSetDatumSamplePosition(nAxis, dPosition));
}

int
ADevScan::nCheckTemplate(const Cstring& sTemplate)
{
  // m_sTemplate is a member variable
  return (DEV_SUCCESS);
}

int
ADevScan::nCheckDiskSpace(const Cstring& sTemplate)
{
  return (DEV_SUCCESS);
}

int
ADevScan::nScan (const int nAxis, const double dStartPosition,
                 const int nNumImages, const double dIncrement,
                 const double dExposureTime, const eDevScan_Mode eMode)
{
  // Not implemented!

  return (DEV_INVALIDMODE);
}

int
ADevScan::nScan (const Cstring& sAxisName, const double dStartPosition,
                 const int nNumImages, const double dIncrement,
                 const double dExposureTime, const eDevScan_Mode eMode)
{
  int nAxis;
  for (nAxis = 0;    (nAxis < m_nSampleDatumAxes)
       && (sAxisName 
           != m_poDevGoniomSample->sGetAxisName(nAxis));
       nAxis++) ;
  return (nScan(nAxis, dStartPosition, nNumImages, dIncrement, dExposureTime,
                eMode));
}

int
ADevScan::nGetLastImageFilename(Cstring *psLastImageFilename)
{
  int nStat;
  int nLastSeqNum;

  nStat = nGetLastSeqNum(&nLastSeqNum);
  if (DEV_SUCCESS == nStat)
    {
      int i, nFirstHash, nLen;
      Cstring sName;
      char *pcTemp;

      sName  = m_sTemplate;
      pcTemp = new char [sName.GetLength()+1];  // This is max length we'll need
      nLen   = sprintf(pcTemp, "%d", abs(nLastSeqNum));
      nFirstHash = -1;
      for (i = sName.GetLength()-1; i >= 0; i--) 
        {
          // Work backwards through sName
          if (   ('#' == sName.GetAt(i))
              || ('?' == sName.GetAt(i)) ) 
            {
              nFirstHash = i;
              if (0 < nLen) 
                sName.SetAt(i, pcTemp[--nLen]);
              else if (0 == nLen)                 // No more digits in pcTemp
                sName.SetAt(i, '0');              //   so insert leading zeroes
            }
        }
      if (0 > nLastSeqNum) 
        {
          if ( (0 <= nFirstHash) && ('0' == sName.GetAt(nFirstHash)) )
            sName.SetAt(nFirstHash, '-');// Have negative sign and have space
                                        // for it in sName
          else
            nLen = -1;
        }

      if ( (0 == nLen) || (-1 == nFirstHash) )
        {
          // All the hash marks were used up or there were no hash marks

          *psLastImageFilename = sName;
          nStat = DEV_SUCCESS;
        }
      else 
        {
          nStat = DEV_INVALIDMODE;
        }
      delete [] pcTemp;
      pcTemp = NULL;
    }

  return (nStat);
}

int
ADevScan::nGetLastSeqNum(int *pnLastSeqNum)
{
  int nLast;
  if (eDevScanState_Done == m_eState)
    {
      nLast = m_nSeqStart + (m_nNumImages - 1) * m_nSeqIncr;
    }
  else if (   (eDevScanState_InProgress == m_eState) 
           || (eDevScanState_Paused     == m_eState)
           || (eDevScanState_Pausing    == m_eState) )
    {
      nLast = m_nSeqCurr - m_nSeqIncr;
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
  *pnLastSeqNum = nLast;
  return (DEV_SUCCESS);
}

int
ADevScan::nGetScanMode (eDevScan_Mode *peMode)
{
  *peMode = m_eMode;
  return (DEV_SUCCESS);
}

int
ADevScan::nSetScanMode (const eDevScan_Mode eMode)
{
  m_eMode = eMode;
  return (DEV_SUCCESS);
}

int
ADevScan::nSetup(Cscan *poScan)
{
  // Setup a hardware scan from values in a software scan object
  int    nStat, nStat1;
  int    nDetAxes, nCrysAxes;
  float  *pfDet, *pfCrys;

  nStat = 0;

  // Set mode

  if (poScan->eGetMode() == eScanMode_StillClosed)
    {
      m_eMode = eDevScanMode_StillClosed;
    }
  else if (poScan->eGetMode() == eScanMode_StillOpen)
    {
      m_eMode = eDevScanMode_StillOpen;
    }
  else if (poScan->eGetMode() == eScanMode_ScanOpen)
    {
      m_eMode = eDevScanMode_ScanOpen;
    }
  else if (poScan->eGetMode() == eScanMode_ScanClosed)
    {
      m_eMode = eDevScanMode_ScanClosed;
    }
  else
    {
      m_eMode = eDevScanMode_Unknown;
      return (DEV_INVALIDMODE);
    }

  // Assume failure in setup until proven otherwise

  m_eState = eDevScanState_Failed;

  // Set scan axis

  nStat = nSetScanAxis(poScan->m_poRotation->sGetName(), (double)0.0);

  if (DEV_SUCCESS != nStat) return (nStat);

  // Set start position

  nStat = nSetStartPosition((double)poScan->m_poRotation->fGetRotStart());
  if (DEV_SUCCESS != nStat) return (nStat);

  // Should check that the start position matches exactly or change the
  //  scan object to reflect the nearby start position.

  if (m_dStartPosition != (double)poScan->m_poRotation->fGetRotStart())
    {
      poScan->m_poRotation->vSetRotStart((float)m_dStartPosition);
    }

  // Set rot increment

  if (   (eDevScanMode_StillOpen   == m_eMode)
      || (eDevScanMode_StillClosed == m_eMode) )
    {
      m_dIncrement = 0.0;
      m_nNumImages = poScan->nGetNumImages();
    }
  else
    {
      nStat = nSetIncrement((double)poScan->m_poRotation->fGetIncrement());
      if (m_dIncrement != (double)poScan->m_poRotation->fGetIncrement())
        {
          poScan->m_poRotation->vSetIncrement((float)m_dIncrement);
        }
      // Compute number of images

      if (0.0 < m_dIncrement)
        {
          m_nNumImages = (int) ( 
                           ((double) poScan->m_poRotation->fGetRotEnd() 
                                - m_dStartPosition)
                            / m_dIncrement + 0.5); // Round up!
        }
      else
        {
          m_nNumImages = poScan->nGetNumImages();
        }
    }

  // Set exposure time

  nStat = nSetExposureTime((double)poScan->m_poRotation->fGetExposureTime());
  if (DEV_SUCCESS != nStat) return (nStat);
  if (m_dExposureTime != (double)poScan->m_poRotation->fGetExposureTime())
    {
      poScan->m_poRotation->vSetExposureTime((float) m_dExposureTime);
    }

  // Deal with datum position
  
  nCrysAxes = m_nSampleDatumAxes;
  nDetAxes  = m_nDetDatumAxes;
  pfCrys    = new float [nCrysAxes];
  pfDet     = new float [nDetAxes];

  nStat1 = poScan->nGetDatum(&nCrysAxes, pfCrys, &nDetAxes, pfDet);
  if (0 != nStat1)
    {
      return (DEV_INVALIDARG);
    }
  else
    {
      int i;
      for (i = 0; i < nDetAxes; i++)
        {
          m_pdDetDatum[i] = (double)pfDet[i];
        }
      for (i = 0; i < nCrysAxes; i++)
        {
          m_pdSampleDatum[i] = (double)pfCrys[i];
        }
    }

  delete [] pfDet;
  delete [] pfCrys;

  // Sequence start and increment

  m_nSeqStart = poScan->nGetSeqNum(0);  // Get first sequence
  m_nSeqIncr  = poScan->nGetSeqInc();

  // Check template and diskspace

  m_sTemplate = poScan->sGetTemplate();
  nStat       = nCheckTemplate(m_sTemplate);
  if (DEV_SUCCESS != nStat) return (nStat);

  nStat = nCheckDiskSpace(m_sTemplate);
  if (DEV_SUCCESS != nStat) return (nStat);
  
//  m_nNumOsc = poScan->m_poRotation->nNumOscillations;
// nDarkIntvl??

  if (DEV_SUCCESS == nStat)
    {
      // All OK, up to now, so do regular setup
      m_eState = eDevScanState_NotReady;
      nStat = nSetup();
    }
  return (nStat);
}

int
ADevScan::nSetup(void)
{
  int nStat;

  // Can only do setup if in certain states, so check

  if (   (m_eState != eDevScanState_NotReady)
      && (m_eState != eDevScanState_Ready)
      && (m_eState != eDevScanState_Done) )
    {
      m_sStatusMsg = "ERROR in ADevScan::nSetup - invalid mode!\n";
      return (DEV_INVALIDMODE);
    }
  else
    {
      // Ok, current state OK, so proceed.  But be cautious!
      // Assume failure, until proven otherwise!

      m_eState = eDevScanState_Failed;

      // Now remember, this is just a simulation!
      // so, do one big exposure and sort things out in ePoll

      double dSpeed;
      m_dTotalTimeNeeded = (double) m_nNumImages * m_dExposureTime;
      m_nSeqCurr = m_nSeqStart;
      nStat = m_poDevDetector->nSetup(m_dTotalTimeNeeded);
      if (DEV_SUCCESS != nStat) 
        {
          m_poDevDetector->vGetStatusMsg(&m_sStatusMsg);
          m_sStatusMsg 
            = "ERROR in ADevScan::nSetup - Problem with detector nSetup.\n"
              + m_sStatusMsg;
          return (nStat);
        }

      // Setup other devices as required: Detector goniometer,
      //                                  sample goniometer, shutter

      // Move to requested detector position

      double dPos;
      int i;

      for (i = 0; i < m_nDetDatumAxes; i++)
	{
	  dPos = m_pdDetDatum[i];
	  nStat = m_poDevGoniomDetector->nSetRequestedPosition(i, dPos);
	  if (DEV_SUCCESS != nStat) return (nStat);
	  nStat = m_poDevGoniomDetector->nSetRequestedSpeed(i, (double)0.0);
	  if (DEV_SUCCESS != nStat) return (nStat);
	}
      nStat = m_poDevGoniomDetector->nSetup();
      if (DEV_SUCCESS != nStat)
	{
	  m_poDevGoniomDetector->vGetStatusMsg(&m_sStatusMsg);
	  m_sStatusMsg 
	    = "ERROR in ADevScan::nSetup - Problem with det goniom setup.\n " 
	      + m_sStatusMsg;
	  return (nStat);
	}
      nStat = m_poDevGoniomDetector->nStart();  
      if (DEV_SUCCESS != nStat)
	{
	  m_poDevGoniomDetector->vGetStatusMsg(&m_sStatusMsg);
	  m_sStatusMsg 
	    = "ERROR in ADevScan::nSetup - Problem with det goniom start.\n " 
	      + m_sStatusMsg;
	  return (nStat);
	}
      nStat = m_poDevGoniomDetector->nWait();
      if (DEV_SUCCESS != nStat)
	{
	  m_poDevGoniomDetector->vGetStatusMsg(&m_sStatusMsg);
	  m_sStatusMsg 
	    = "ERROR in ADevScan::nSetup - Problem with det goniom wait.\n " 
	      + m_sStatusMsg;
	  return (nStat);
	}

      // Now if a scan move the sample goniometer to the starting position

      if (   (m_eMode == eDevScanMode_ScanClosed)
          || (m_eMode == eDevScanMode_ScanOpen) )
        {
          // Move to starting datum position plus scanstart

          for (i = 0; i < m_nSampleDatumAxes; i++)
            {
              dPos = m_pdSampleDatum[i];
              if (i == m_nScanAxis)
                {
                  dPos = dPos + m_dStartPosition;
                }
              nStat = m_poDevGoniomSample->nSetRequestedPosition(i, dPos);
              if (DEV_SUCCESS != nStat) return (nStat);
              nStat = m_poDevGoniomSample->nSetRequestedSpeed(i, (double)0.0); // Slew
              if (DEV_SUCCESS != nStat) return (nStat);
            }
          nStat = m_poDevGoniomSample->nSetup();
          if (DEV_SUCCESS != nStat)
            {
              m_poDevGoniomSample->vGetStatusMsg(&m_sStatusMsg);
              m_sStatusMsg 
                = "ERROR in ADevScan::nSetup - Problem with goniom setup.\n " 
                  + m_sStatusMsg;
              return (nStat);
            }
          nStat = m_poDevGoniomSample->nStart();  
          if (DEV_SUCCESS != nStat)
            {
              m_poDevGoniomSample->vGetStatusMsg(&m_sStatusMsg);
              m_sStatusMsg 
                = "ERROR in ADevScan::nSetup - Problem with goniom start.\n " 
                  + m_sStatusMsg;
              return (nStat);
            }
          nStat = m_poDevGoniomSample->nWait();
          if (DEV_SUCCESS != nStat)
            {
              m_poDevGoniomSample->vGetStatusMsg(&m_sStatusMsg);
              m_sStatusMsg 
                = "ERROR in ADevScan::nSetup - Problem with goniom wait.\n " 
                  + m_sStatusMsg;
              return (nStat);
            }

          // Setup scan speed
      
          if (0.0 < m_dExposureTime)
            {
              dSpeed = m_dIncrement / m_dExposureTime;
            }
          else
            {
              dSpeed = 0.0;
            }
          nStat = m_poDevGoniomSample->nSetRequestedSpeed(m_nScanAxis, dSpeed);
          if (DEV_SUCCESS != nStat) 
            {
              m_poDevGoniomSample->vGetStatusMsg(&m_sStatusMsg);
              m_sStatusMsg 
                = "ERROR in ADevScan::nSetup - Problem with goniom setspeed.\n"
                  + m_sStatusMsg;
              return (nStat);
            }
          m_dEndPosition = m_dStartPosition + m_pdSampleDatum[m_nScanAxis]
                            + ( (double)m_nNumImages * m_dIncrement);
          nStat = m_poDevGoniomSample->nSetRequestedPosition(m_nScanAxis, 
                                                              m_dEndPosition);
          if (DEV_SUCCESS != nStat)
            {
              m_poDevGoniomSample->vGetStatusMsg(&m_sStatusMsg);
              m_sStatusMsg 
                = "ERROR in ADevScan::nSetup - Problem with goniom setpos.\n"
                  + m_sStatusMsg;
              return (nStat);
            }

          nStat = m_poDevGoniomSample->nSetup();
          if (DEV_SUCCESS != nStat)
            {
              m_poDevGoniomSample->vGetStatusMsg(&m_sStatusMsg);
              m_sStatusMsg 
                = "ERROR in ADevScan::nSetup - Problem with goniom setup.\n"
                  + m_sStatusMsg;
              return (nStat);
            }
        }

      if (   (eDevScanMode_StillOpen == m_eMode)
          || (eDevScanMode_ScanOpen   == m_eMode) )
        {
          nStat = m_poDevShutter->nSetup(m_dTotalTimeNeeded);
          nStat = DEV_SUCCESS;
          if (DEV_SUCCESS != nStat) return (nStat);
        }

      // If we get here, we are probably free of errors
      m_eState       = eDevScanState_Ready;
      m_dElapsedTime = 0.0;

      return (nStat);
    }
}

int
ADevScan::nStart(void)
{
  // This simulation of a scan, does a scan as a single long exposure.
  // In reality, some systems will do scans as a series of shorter exposures.

  int nStat;
  time_t tStartTime;
  if (m_eState == eDevScanState_Ready)
    {
      (void) time(&tStartTime);
      m_dStartTime      = (double) tStartTime;
      m_dElapsedTime    = (double) 0.0;
      m_dStartPause     = (double) 0.0;
      m_nPauseSeqNumber = 32767;             // Should be > last seq number
      m_eState          = eDevScanState_InProgress;
      nStat = m_poDevDetector->nStart();
      if (DEV_SUCCESS != nStat)
        {
          m_poDevDetector->vGetStatusMsg(&m_sStatusMsg);
          m_sStatusMsg 
            = "ERROR in ADevScan::nStart - Problem with detector nStart.\n"
              + m_sStatusMsg;
          return (nStat);
        }
      else if (   (m_eMode == eDevScanMode_ScanClosed)
            || (m_eMode == eDevScanMode_ScanOpen)   )
        {
          nStat = m_poDevGoniomSample->nStart();
        }
      if ( (DEV_SUCCESS == nStat) 
          && ( (m_eMode == eDevScanMode_StillOpen)
          || (m_eMode == eDevScanMode_ScanOpen) ) )
        {
          nStat = m_poDevShutter->nStart();
        }
      return (nStat);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}

int
ADevScan::nPoll(eDevScan_State *peState, double *pdPosition,
                double *pdElapsed, int *pnLastSeqNum)
{
  time_t tCurrentTime;  
  int    nExposing;
  double dDelta;
  double dSpeed;
  int    nStat;

  if (   (m_eState == eDevScanState_InProgress)
      || (m_eState == eDevScanState_Pausing) )
    {
      (void) time(&tCurrentTime);
      m_dElapsedTime = (double) tCurrentTime - m_dStartTime;

      // Calculate current position based on elapsed time and speed

      nExposing = 0;
      if (m_dElapsedTime <= m_dTotalTimeNeeded)
        {
          // Still exposing

          nExposing++;
          if (   (eDevScanMode_ScanOpen == m_eMode)
              || (eDevScanMode_ScanClosed == m_eMode))
            {
              // A moving axis, so compute how far it should have moved

              nStat = m_poDevGoniomSample->nGetSpeed(m_nScanAxis, &dSpeed);
              if (DEV_SUCCESS != nStat) return (nStat);

              dDelta = m_dElapsedTime * dSpeed;
               m_dCurrentPosition = m_dStartPosition + dDelta;
               if (m_dCurrentPosition >=  m_dEndPosition) 
                 {
                   // Got there, so not moving anymore
                   m_dCurrentPosition = m_dEndPosition;
                 }
             }
        }
      else
        {
          // A still mode

          m_dElapsedTime = m_dTotalTimeNeeded;
        }

      // Now compute current sequence number based on
      // current elapsed time (it's simulation!!)

      if (0.0 < m_dExposureTime)
        {
          m_nSeqCurr = m_nSeqStart 
                        +  m_nSeqIncr
                        * (int) (m_dElapsedTime / m_dExposureTime);
        }

      if (0 == nExposing)
        {
          m_eState       = eDevScanState_Done;

          // We think scan is finished, so poll our device friends for clean-up

          nStat = m_poDevDetector->nPoll();
          nStat = m_poDevGoniomDetector->nPoll();
          nStat = m_poDevGoniomSample->nPoll();
        }
    }
  else if (eDevScanState_Paused == m_eState)
    {
      (void) time(&tCurrentTime);

      // Add paused time to start time, so math still works
      m_dStartTime = m_dStartTime + ((double) tCurrentTime - m_dStartPause);
    }


  if (NULL != pdPosition)
    {
      *pdPosition = m_dCurrentPosition;
    }
  if (NULL != pdElapsed)
    {
      *pdElapsed = m_dElapsedTime; 
    }
  if (NULL != pnLastSeqNum)
    {
      *pnLastSeqNum = m_nSeqCurr - m_nSeqIncr;
    }
  if (NULL != peState)
    {
      *peState = m_eState;
    }

  if (   (eDevScanState_NotReady   == m_eState)
      || (eDevScanState_Ready      == m_eState) 
      || (eDevScanState_Starting   == m_eState)
      || (eDevScanState_InProgress == m_eState)
      || (eDevScanState_Paused     == m_eState)
      || (eDevScanState_Resuming   == m_eState)
      || (eDevScanState_Done       == m_eState)
      )
    return (DEV_SUCCESS);
  else
    return (DEV_INVALIDMODE);
}

int
ADevScan::nWait(const int nImageSeqNumber)
{
  // Wait for next image .......if nImageSeqNumber = 0
  // Wait for indicated image.. if nImageSeqNumber > 0
  // Sequence numbers start at 1 so that a 0 means none have been done

  int nStat = DEV_SUCCESS;
  int nTest, nLast;
  
  if (0 == nImageSeqNumber)
    {
      nStat = nGetLastSeqNum(&nTest);
      if (DEV_SUCCESS != nStat) return (nStat);
    }
  else
    {
      nTest = nImageSeqNumber;   // Watch out, may never find it, nTest > max
    }                           // possible!

  nStat = nGetLastSeqNum(&nLast);
  if (DEV_SUCCESS != nStat) return (nStat);
  while (nLast >= nTest)
    {
      while (   (DEV_SUCCESS == nStat) 
             && (   (eDevScanState_InProgress == m_eState)
                 || (eDevScanState_Pausing    == m_eState)
                 || (eDevScanState_Paused     == m_eState)
                 || (eDevScanState_Resuming   == m_eState)
                 || (eDevScanState_Starting   == m_eState) )
             )
        {
          nStat = nPoll();                      // Poll until not moving
        }
    }
         
  return (nStat);
}

int
ADevScan::nAbort(void)
{
  if (   (m_eState == eDevScanState_InProgress) 
      || (m_eState == eDevScanState_Paused)
      || (m_eState == eDevScanState_Pausing) )

    {
      if (   (m_eMode == eDevScanMode_StillOpen)
          || (m_eMode == eDevScanMode_ScanOpen) )
        m_poDevShutter->nAbort();
      if (   (eDevScanMode_ScanOpen == m_eMode)
          || (eDevScanMode_ScanClosed == m_eMode))
        m_poDevGoniomSample->nAbort();

//      m_poDevGoniomDetector->nAbort();
      m_poDevDetector->nAbort();
      m_poDevShutter->nAbort();
      m_eState = eDevScanState_Aborted;
      
      nInit();
    }
  return (DEV_SUCCESS);
}

int
ADevScan::nDone(void)
{
  m_eState = eDevScanState_Unknown;
  return (DEV_SUCCESS);
}

int
ADevScan::nInit(void)
{
  // Initialize devices used by the scan pseudo-device if needed

  int nStat = DEV_SUCCESS;
  m_eState = eDevScanState_Failed;   // Assume failure until proven otherwise
  if (NULL != m_poDevDetector)
    {
      nStat = m_poDevDetector->nInit();
      if (DEV_SUCCESS != nStat) 
        {
          m_poDevDetector->vGetStatusMsg(&m_sStatusMsg);
          m_sStatusMsg 
            = "ERROR in ADevScan::nInit - Problem with detector nInit.\n"
              + m_sStatusMsg;
          return (nStat);
        }
    }
  if (NULL != m_poDevGoniomSample)
    {
      nStat = m_poDevGoniomSample->nPoll();
      if (DEV_SUCCESS != nStat)
        {
          nStat = m_poDevGoniomSample->nInit();
          if (DEV_SUCCESS != nStat) return (nStat);
        }
    }
  if (NULL != m_poDevShutter)
    {
      nStat = m_poDevShutter->nPoll();
      if (DEV_SUCCESS != nStat)
        {
          nStat = m_poDevShutter->nInit();
          if (DEV_SUCCESS != nStat) return (nStat);
        }
    }
  m_eState = eDevScanState_NotReady;  // Success!
  return (DEV_SUCCESS);
}

int
ADevScan::nPauseResume(const int nMode)
{
  // nMode = 0 for pause
  //      != 0 for resume

  if (0 == nMode)
    {
      if (   (m_eState == eDevScanState_InProgress) 
          || (m_eState == eDevScanState_Starting)
          || (m_eState == eDevScanState_Resuming) )
        {
          time_t tCurrentTime;  
          (void) time(&tCurrentTime);
          m_dStartPause = (double) tCurrentTime;
          m_eState = eDevScanState_Pausing;
          return (DEV_SUCCESS);
        }
      else
        {
          return (DEV_INVALIDSTATE);
        }
    }
  else if (   (m_eState == eDevScanState_Paused)
           || (m_eState == eDevScanState_Pausing) )
    {
//      m_eState = eDevScanState_Resuming;
      m_eState = eDevScanState_InProgress;
      return (DEV_SUCCESS);
    }
  else
    {
      return (DEV_INVALIDSTATE);
    }
}
