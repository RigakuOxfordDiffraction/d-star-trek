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
// CDevScan.cc      Initial author: J.W. Pflugrath           10-Dec-1995
//  This file contains the member functions of class CDevScan which
//     is derived from base class ADevScan
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

#include "CDevScan.h"

//+Code begin

//+Public functions

// Constructors, destructors and assignments

CDevScan::CDevScan (CDevGoniom   *poDevGoniomDetector, 
		    CDevGoniom   *poDevGoniomSample, 
		    CDevDetector *poDevDetector, 
		    CDevShutter  *poDevShutter, 
		    CDevCounter  *poDevCounter) : ADevScan(poDevGoniomDetector,
							   poDevGoniomSample,
							   poDevDetector,
							   poDevShutter,
							   poDevCounter)
{
  m_poDevDetector        = poDevDetector;
  m_poDevGoniomDetector  = poDevGoniomDetector;
  m_poDevGoniomSample    = poDevGoniomSample;
  m_poDevShutter         = poDevShutter;
  m_poDevCounter         = poDevCounter;

  vConstruct();
}

CDevScan::~CDevScan()
{
  (void) nDone();
  vDestruct();
}

void
CDevScan::vConstruct(void)
{
//  cout << "CDevScan::vConstruct called!\n";

  m_sTemplate            = "./image.????";
  m_nSeqStart            = 1;
  m_nSeqIncr             = 1;
  m_nSeqCurr             = 1;
  m_nNumImages           = 1;
  m_nScanAxis            = 0;
  m_dIncrement           = (double) 1.0;
  m_dStartPosition       = (double) 0.0;
  m_dExposureTime        = (double) 180.0;
  m_dElapsedTime         = (double) 0.0;
  m_dCurrentPosition     = (double) 0.0;
  m_eState               = eDevScanState_Unknown;
  m_eMode                = eDevScanMode_StillClosed;
  m_nSampleDatumAxes     = 1;
  m_nDetDatumAxes        = 1;
  if (NULL != m_poDevGoniomDetector)
    m_nDetDatumAxes = m_poDevGoniomDetector->nGetNumAxes();
  if (NULL != m_poDevGoniomSample)
    m_nSampleDatumAxes = m_poDevGoniomSample->nGetNumAxes();
  m_pdDetDatum    = new double [m_nDetDatumAxes];
  m_pdSampleDatum = new double [m_nSampleDatumAxes];
  m_sStatusMsg    = "No scan status message";
}

void
CDevScan::vDestruct(void)
{
//  cout << "CDevScan::vDestruct called!\n";
  nDone();
  delete [] m_pdDetDatum;
  delete [] m_pdSampleDatum;
}

int
CDevScan::nSetHeader(const Cimage_header& oHeader)
{
  return (DEV_SUCCESS);
}
