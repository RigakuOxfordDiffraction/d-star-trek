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
// CscanMS.cc            Initial author: J.W. Pflugrath           24-Mar-1995
//  This file contains the member functions of class CscanMS which implements
//    the scan encapsulation of d*TREK for Microsoft platforms.
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

#include "Dtrek.h"
#include "Cscan.h"         // Class definition and prototypes

//+Definitions, constants, and initialization of static member variables

Cstring Cscan::m_ssScanRotPrefix = "SCAN_";
Cstring Cscan::m_ssScanTemplate  = "SCAN_TEMPLATE";
Cstring Cscan::m_ssScanKey       = "SCAN_KEY";
Cstring Cscan::m_ssScanSeqInfo   = "SCAN_SEQ_INFO";
Cstring Cscan::m_ssScanWavelength= "SCAN_WAVELENGTH";
Cstring Cscan::m_ssScanCrysDatum = "SCAN_CRYS_RELZERO";
Cstring Cscan::m_ssScanDetDatum  = "SCAN_DET_RELZERO";
Cstring Cscan::m_ssScanMode      = "SCAN_MODE";
Cstring Cscan::m_ssScanModeStillO= "Still_Open";
Cstring Cscan::m_ssScanModeStillC= "Still_Closed";
Cstring Cscan::m_ssScanModeScanO = "Scan_Open";
Cstring Cscan::m_ssScanModeScanC = "Scan_Closed";


//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cscan::Cscan()
{
 (void) nInitValues();
}

Cscan::Cscan(const Cstring sTemplate, const int nSeqNum = 0, 
	     const int nStepNum = 1)
{
 (void) nInitValues();
 sFileTemplate    = sTemplate;
 nFileStartSeqNum = nSeqNum;
 nFileStepNum     = nStepNum;
 nFileSeqNum      = nFileStartSeqNum;
}

Cscan::~Cscan()
{
  if (poRotation != NULL)
    {
      delete poRotation;
    }
  if (pfCrysDatum != NULL)
    {
      delete [] pfCrysDatum;
      pfCrysDatum = NULL;
      nNumCrysDatum = 0;
    }
  if (pfDetDatum != NULL)
    {
      delete [] pfDetDatum;
      pfDetDatum = NULL;
      nNumDetDatum = 0;
    }
}

int Cscan::nInitValues(void)
{
  eThe_State              = eScan_unknown_state;
  sScan_key               = "unknown";
  sFileTemplate           = "";
  nFileStartSeqNum        = 0;
  nFileStepNum            = 1;
  nFileSeqNum             = nFileStartSeqNum;
  nNumImages              = 0;
  poTheDCoffsetImages     = NULL;
  poRotation              = new Crotation (); // Default rotation
  nNumDetectors           = 0;
  nNumCounters            = 0;
  nNumDetDatum            = 0;
  pfDetDatum              = NULL;
  nNumCrysDatum           = 0;
  pfCrysDatum             = NULL;
  fWavelength             = 0.0;

  return (0);
}

int Cscan::nList(void)
{
  cout << "Scan key:  " << sScan_key << endl;
  if (poRotation != NULL) 
    (void) poRotation->nList();  
  else
    cout << "Rotation value unknown.\n";
  cout << "Template: " << sFileTemplate
       << "\n   Start: " << nFileStartSeqNum
       << "\n    Step: " << nFileStepNum << endl;
  return (0);
}

int Cscan::nSetRotation(const Crotation& oRotation)
{
  if (poRotation != NULL)
    {
      delete poRotation;  // Delete current Crotation object and 
    }
  poRotation = new Crotation(oRotation);
  return (0);
}

int Cscan::nGetImageName(Cstring* sName)
{

  // Create a filename from the sFileTemplate and nFileSeqNum

  int  i, j, nLen, nFirstHash;
  int  nStat;
  bool bNegative;

  char *cTemp;
  cTemp = new char [sFileTemplate.length()+1];  // This is max length we'll need

  *sName = sFileTemplate;

  bNegative = (nFileSeqNum < 0);            // Take care of negative sign later

  nLen = sprintf(cTemp, "%d", abs(nFileSeqNum));

  if (nLen <= 0) {
    cerr << "ERROR: Cscan: writing nFileSeqNum: " << nFileSeqNum << endl;
    nStat = nLen;
  }
  else {
    nFirstHash = -1;
    for (i = sName->length()-1; i >= 0; i--) {  // Work backwards through sName
      if (   (sName->at(i, 1) == "#")
	  || (sName->at(i, 1) == "?") ) {
	nFirstHash = i;
	if (nLen > 0) 
	  sName->at(i,1) = cTemp[--nLen];
	else if (nLen == 0)                 // No more digits in cTemp
	  sName->at(i,1) = "0";             //   so insert leading zeroes
      }
    }
    if (bNegative) {
      if ( (nFirstHash >= 0) && (sName->at(nFirstHash, 1) == "0") )
	sName->at(nFirstHash, 1) = "-";    // Have negative sign and have space
                                           // for it in sName
      else
	nLen = -1;
    }

    if ( (nLen == 0) || (nFirstHash == -1) )

      // All the hash marks were used up or there were no hash marks

      nStat = 0;
    else {
      cerr << "ERROR: Cscan not enough ?'s in sFileTemplate!" << endl;
      nStat = 1;
    }
  }

  delete [] cTemp;
  cTemp = NULL;

  return (nStat);
}
int
Cscan::nSetDatum(const int nCrys, float *pfCrys, const int nDet, float *pfDet)
{
  int i;
  int nStat = 0;
  if (nCrys > 0)
    {
      nNumCrysDatum = nCrys;
      if (pfCrysDatum != NULL)
	{
	  delete [] pfCrysDatum;
	}
	pfCrysDatum = new float [nNumCrysDatum];

      for (i = 0; i < nNumCrysDatum; i++)
	{
	  pfCrysDatum[i] = pfCrys[i];
	}
    }
  else
    {
      nStat = 1;
    }
  if (nDet > 0)
    {
      nNumDetDatum = nDet;
      if (pfDetDatum != NULL)
	{
	  delete [] pfDetDatum;
	}
	pfDetDatum = new float [nNumDetDatum];

      for (i = 0; i < nNumDetDatum; i++)
	{
	  pfDetDatum[i] = pfDet[i];
	}
    }
  else
    {
      nStat = nStat + 1;
    }
  return (nStat);
}

int
Cscan::nGetDatum(int *pnCrys, float *pfCrys, int *pnDet, float *pfDet)
{
  // Get the datum values in the scan object and return them in the
  // passed arguments.
  // On entry, *pnCrys has the dimensions of the *pfCrys array
  //           *pnDet  has the dimensions of the *pfDet  array
  // On exit   *pnCrys has the number of datum values transferred to *pfCrys
  //           *pnDet  has the number of datum values transferred to *pfDet
  // If there was enough space in the input arguments, return 0
  // otherwise return 1 if not enough dimensions in *pfCrys
  //                  2 if not enough dimensions in *pfDet
  //                  3 if not enough dimensions in both 

  int i;
  int nStat = 0;
  int nDim;
  nDim = nNumCrysDatum;
  if (*pnCrys < nDim)
    {
      nStat = 1; // Not enough space
      nDim = *pnCrys;
    }
  *pnCrys = nDim;
  for (i = 0; i < nDim; i++)
    {
      pfCrys[i] = pfCrysDatum[i];
    }

  nDim = nNumDetDatum;
  if (*pnDet < nDim)
    {
      nStat = nStat + 2;
      nDim = *pnDet;
    }
  *pnDet = nDim;
  for (i = 0; i < nNumDetDatum; i++)
    {
      pfDet[i] = pfDetDatum[i];
    }
  return (nStat);
}
int
Cscan::nGetNumImages(void)
{
  if (poRotation->fGetIncrement() > 0.0)
    {
      return ((int)(  (poRotation->fGetRotEnd() - poRotation->fGetRotStart())
		    / poRotation->fGetIncrement()  + 0.5)); // Round-up
    }
  else
    {
      return (1);  // Assume always at least one image in a scan?
    }
}
