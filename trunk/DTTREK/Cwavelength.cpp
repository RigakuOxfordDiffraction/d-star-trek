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
// Cwavelength.cc         Initial author: J.W. Pflugrath           01-May-1995
//    This file contains the member functions of class Cwavelength.
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

//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Dtrek.h"
#include "Cwavelength.h"         // Class definition and prototypes
                                 //  Cwavelength.h includes others
#include "dtrekdefs.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

//+Definitions, constants, and initialization of static member variables

Cstring Cwavelength::ms_sWavelength         = D_K_Wavelength;
Cstring Cwavelength::ms_sWavelengthKey      = D_K_WavelengthKey;
float Cwavelength::ms_afKnownKa1Ka2Values[] = {
    1.5405f,1.5443f,0.70926f,0.71354f,1.9360f,1.9399f,2.2896f,2.2935f,0.0f
};


//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cwavelength::Cwavelength()
{
  (void) nInitValues();
  m_pfWavelength  = new float [1];
  m_nAlloc        = 1;
  *m_pfWavelength = 0.0;
}

Cwavelength::Cwavelength(const float fLambda)
{
  (void) nInitValues();

  m_eThe_Type     = eWavelength_single_type;
  m_pfWavelength  = new float [1];
  m_nAlloc        = 1;
  *m_pfWavelength = fLambda;
  nVerifyKa1Ka2Information();
}

Cwavelength::Cwavelength(const int nNum, const float* pfLambda)
{
  (void) nInitValues();


  float  *pfTemp;

  m_eThe_Type    = ((nNum==1)?eWavelength_single_type:eWavelength_multiple_type);
  m_pfWavelength = new float [nNum];
  m_nAlloc       = nNum;

  pfTemp       = m_pfWavelength;

  for (int i = 0; i < nNum; i++)
    {
      *pfTemp++ = *(pfLambda+i);
      if (0.0 >= *(m_pfWavelength+i))
	{
        cout << "\n\n\nERROR - wavelength <= 0.0! Reset to to 1 Angstrom!\n\n\n" 
	       << flush;
      
      *(m_pfWavelength+i) = 1.0;
	}
    }    
  nVerifyKa1Ka2Information();
}

Cwavelength::Cwavelength(Cimage_header& oHeader, const Cstring sPrefix, bool bA1A2Check)
{
  int     nStat;
  Cstring sTemp;

  nStat = nInitValues();
  
  m_bA1A2Check = bA1A2Check;
  
  if (!oHeader.bIsAvailable())
    {
      cout << "Cwavelength ERROR: image header is not valid.  "
           << "Cannot construct wavelength!\n";
      nStat = 1;
    }
  else if (0 == oHeader.nGetValue(sPrefix + ms_sWavelengthKey, &sTemp))
    {
      cout << "Cscan ERROR: cannot use database key yet!\n";
      nStat = 2;
    }
  else
    {
      // Try to get required information from the header

      float fValues[2];
      int   nTemp;

      nTemp = oHeader.nGetValue(sPrefix + ms_sWavelength, 2, fValues);
      if (0 == nTemp)
	{
	  m_nAlloc = (int) fValues[0];
	  if (1 == m_nAlloc)
	    {
	      m_eThe_Type = eWavelength_single_type;
	      nStat = nSetWavelength(fValues[1]);
	    }
	  else
	    {
	      float *pfTemp;
	      pfTemp = new float [m_nAlloc+1];
	      nTemp = oHeader.nGetValue(sPrefix + ms_sWavelength, 
					m_nAlloc+1, pfTemp);
	      if (0 == nTemp)
		{
		  m_eThe_Type = eWavelength_multiple_type;
		  nStat     = nSetWavelengths(m_nAlloc, pfTemp+1);
		}
	      delete [] pfTemp;
	    }
	}
    }
  
    nVerifyKa1Ka2Information();
}

int
Cwavelength::nVerifyKa1Ka2Information()
{
    if( !m_bA1A2Check )
        return 0; // nothing to do
    
    int nWave;
    double fNominal;
    if (m_nAlloc==1) {
        fNominal = *m_pfWavelength;
        // Only one wavelength?  
        // Should have two.  Search our database for the "real" values.
        for (nWave = 0; ms_afKnownKa1Ka2Values[nWave]!=0.0; nWave+=2) {
            if ((ms_afKnownKa1Ka2Values[nWave]-0.0001<=fNominal) && (ms_afKnownKa1Ka2Values[nWave+1]+0.0001>=fNominal)) {
                // We found a match.  Use it.
                delete[] m_pfWavelength;
                m_pfWavelength = new float[3];
				m_pfWavelength[0] = fNominal;
                m_pfWavelength[1] = ms_afKnownKa1Ka2Values[nWave];
                m_pfWavelength[2] = ms_afKnownKa1Ka2Values[nWave+1];
                m_nAlloc = 3;
                break;
            };
        };
    };
    return 0;
};

Cwavelength::~Cwavelength()
{
  if (NULL != m_pfWavelength) 
    {
      delete [] m_pfWavelength;
      m_pfWavelength = NULL;
      m_nAlloc = 0;
    }
}

int Cwavelength::nInitValues(void) 
{
  m_eThe_Type    = eWavelength_unknown_type;
  m_pfWavelength = NULL;
  m_nAlloc       = 0;
  m_sInputUnits  = "unknown";

  m_bA1A2Check = TRUE;

  return (0);
}  

int Cwavelength::nSetWavelength(const float fLambda)
{
  if (NULL == m_pfWavelength)
    {
      m_pfWavelength = new float [1];
      m_nAlloc = 1;
    }
  *m_pfWavelength = fLambda;
  if (0.0 >= *m_pfWavelength)
    {
      cout << "\n\n\nERROR - wavelength <= 0.0! Reset to to 1 Angstrom!\n\n\n" << flush;
      *m_pfWavelength = 1.0;
    }
  if (m_eThe_Type == eWavelength_single_type) return (0);
  return (1);
}

int Cwavelength::nSetWavelengths(const int nNum, const float* pfLambda)
{
  if ( (NULL == m_pfWavelength) ||  (m_nAlloc <= nNum) )
    {
      if (NULL != m_pfWavelength) delete [] m_pfWavelength;
      m_nAlloc = nNum;
      m_pfWavelength = new float [m_nAlloc];
    }
  for (int i = 0; i < nNum; i++)
    {
      *(m_pfWavelength+i) = *(pfLambda+i);
      if (0.0 >= *(m_pfWavelength+i))
	{
	  cout << "\n\n\nERROR - wavelength <= 0.0! Reset to to 1 Angstrom!\n\n\n" 
	       << flush;
	  *(m_pfWavelength+i) = 1.0;
	}
    }
  return (0);
}

int Cwavelength::nList(const int nFlag)
{
 if (m_eThe_Type == eWavelength_single_type) 
   {
     cout << "  Single wavelength: " << *m_pfWavelength << endl; 
   }
 else if (m_eThe_Type == eWavelength_multiple_type) 
   {
     cout << " Multiple wavelengths: " << endl;
   }
 else
   {
     cout << " Wavelength:      " << "UNKNOWN!" << endl; 
   }
 return (0);
}

int
Cwavelength::nUpdateHeader(Cimage_header* poHeader, const Cstring &sPre)
{
  int     nStat;
  Cstring sPrefix;
  sPrefix = sPre;
  if ("" == sPrefix)
    {
      // If the input prefix is blank, then use the member variable as the
      // prefix.

//      sPrefix = m_sPrefix;
    }
  
  float *pfTemp;
  pfTemp = new float [m_nAlloc + 1];
  pfTemp[0] = (float) m_nAlloc;
  for (int i = 0; i < m_nAlloc; i++)
    {
      pfTemp[i+1] = m_pfWavelength[i];
    }
  nStat = poHeader->nReplaceValue(sPre + ms_sWavelength, m_nAlloc+1, pfTemp, 6);
  delete pfTemp;
  return (nStat);
}

