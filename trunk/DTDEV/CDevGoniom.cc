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
// CDevGoniom.cc      Initial author: J.W. Pflugrath           10-Dec-1995
//  This file contains the member functions of class CDevGoniom which
//     is derived from base class ADevGoniom
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

#include "CDevGoniom.h"

//+Code begin

//+Public functions

// Constructors, destructors and assignments

CDevGoniom::CDevGoniom(const Cstring sName)
{
  m_sName = sName;
  vConstruct();
}

CDevGoniom::~CDevGoniom()
{
  (void) nDone();
  vDestruct();
}

void
CDevGoniom::vConstruct(void)
{
  int i;
//  m_sName         = sName;
  if ("Det" != m_sName)
    {
      m_nNumAxes = 3;             // Crystal has Omega, Chi, & Phi
      m_psAxisName    = new Cstring [m_nNumAxes];
      m_psAxisName[0] = "Omega";
      m_psAxisName[1] = "Chi";
      m_psAxisName[2] = "Phi";
    }
  else 
    {
      m_nNumAxes = 2;             // Detector has dist and swing
      m_psAxisName    = new Cstring [m_nNumAxes];
      m_psAxisName[0] = "Dist";
      m_psAxisName[1] = "Theta";
    }
  m_pdCurrent        = new double [m_nNumAxes];
  m_pdStart          = new double [m_nNumAxes];
  m_pdSet            = new double [m_nNumAxes];
  m_pdRequested      = new double [m_nNumAxes];
  m_pdCurrentSpeed   = new double [m_nNumAxes];
  m_pdSetSpeed       = new double [m_nNumAxes];
  m_pdRequestedSpeed = new double [m_nNumAxes];

  m_pdMin            = new double [m_nNumAxes];
  m_pdMax            = new double [m_nNumAxes];

  m_pdRequestedAccel = new double [m_nNumAxes];
  m_pdSetAccel       = new double [m_nNumAxes];
  m_pdCurrentAccel   = new double [m_nNumAxes];

  m_pdMinSpeed       = new double [m_nNumAxes];
  m_pdMaxSpeed       = new double [m_nNumAxes];
  m_pdMinAccel       = new double [m_nNumAxes];
  m_pdMaxAccel       = new double [m_nNumAxes];

  m_peStateAxis      = new eDevGoniom_State [m_nNumAxes];
  for (i = 0; i < m_nNumAxes; i++)
    {
      m_pdCurrent[i]        = 0.0;
      m_pdSet[i]            = 0.0;
      m_pdRequested[i]      = 0.0;
      m_pdCurrentSpeed[i]   = 50.0;
      m_pdSetSpeed[i]       = 50.0;
      m_pdRequestedSpeed[i] = 50.0;
      m_peStateAxis[i]      = eDevGoniomState_Unknown;
      m_pdMinSpeed[i]       = 0.0;
      m_pdMaxSpeed[i]       = 1000.0;
      m_pdMinAccel[i]       = 0.0;
      m_pdMaxAccel[i]       = 1000.0;
      m_pdMin[i]            = -1000.0;
      m_pdMax[i]            = +1000.0;
    }

  if (m_sName == "Det")
    {
      m_pdCurrent[0]   = 100.0;
      m_pdSet[0]       = 100.0;
      m_pdRequested[0] = 100.0;
    }

  m_eState   = eDevGoniomState_Unknown;
  m_sStatusMsg     = "No goniometer status message";
}

void
CDevGoniom::vDestruct(void)
{
  delete [] m_psAxisName;
  delete [] m_pdCurrent;
  delete [] m_pdStart;
  delete [] m_pdSet;
  delete [] m_pdRequested;
  delete [] m_pdCurrentSpeed;
  delete [] m_pdSetSpeed;
  delete [] m_pdRequestedSpeed;
  delete [] m_pdMin;
  delete [] m_pdMax;

  delete [] m_pdRequestedAccel;
  delete [] m_pdSetAccel;
  delete [] m_pdCurrentAccel;

  delete [] m_pdMinSpeed;
  delete [] m_pdMaxSpeed;
  delete [] m_pdMinAccel;
  delete [] m_pdMaxAccel;
}

