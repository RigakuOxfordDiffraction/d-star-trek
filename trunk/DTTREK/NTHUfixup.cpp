//
// Copyright (c) 2002 Rigaku/MSC, Inc.
//
// NTHUfixup.cc     Initial author: J.W. Pflugrath        04-June-2002
//    Fixup an image header from the NTHU system for use with dtorient
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
//   Check if the CRYSTAL_GONIO_* keywords in the header indicate this
//   header comes from the NTHU system.  If not, simply return.
//   If the header comes from the NTHU system, then
//   1. Determine if the theta value is 0, if not, then error!
//   2. Change the CRYSTAL_ Cgoniometer object so that is does 
//        NOT have 4 axes by removing all information for the FIRST axis.
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include <stdio.h>        // for sprintf, until we figure out sstream

#include "Cstring.h"
#include "Dtrek.h"
#include "Cimage_header.h"
#include "Cgoniometer.h"
#include "Ccrystal.h"
#include "Csource.h"
#include "dtrekvec.h"

using std::cout;
using std::endl;
using std::flush;

//+Public functions

int DTREK_EXPORT nNTHUfixup(Cimage_header &roHeader)
{
  int nStat = 1;
  bool bNTHUGonio = FALSE;
  Cgoniometer *poCrysGonio;
  poCrysGonio   = new Cgoniometer(roHeader, Ccrystal::ms_sCrystalPrefix);
  if (poCrysGonio->bIsAvailable())
    if (4 == poCrysGonio->nGetNumValues())
      if ("Psi" == poCrysGonio->sGetName(1))
	    bNTHUGonio = TRUE;

  if (bNTHUGonio)
    {
      Cstring sTemp;
      float a12fVecs[12];

      nStat = roHeader.nGetValue(Ccrystal::ms_sCrystalPrefix
				      + Cgoniometer::ms_sGonioValues,
				  4, a12fVecs);

      if (0.0 == a12fVecs[0])
	{ 
	  nStat = roHeader.nReplaceValue(Ccrystal::ms_sCrystalPrefix
					    + Cgoniometer::ms_sGonioValues,
					    3, &a12fVecs[1]);
	  nStat = roHeader.nReplaceValue(Ccrystal::ms_sCrystalPrefix
					  + Cgoniometer::ms_sGonioNumValues,
					  3);
	  nStat = roHeader.nReplaceValue(Ccrystal::ms_sCrystalPrefix
					  + Cgoniometer::ms_sGonioNames,
					  "Psi Kappa Phi");
	  nStat = roHeader.nReplaceValue(Ccrystal::ms_sCrystalPrefix
					  + Cgoniometer::ms_sGonioUnits,
					  "deg deg deg");
     
	  nStat = roHeader.nGetValue(Ccrystal::ms_sCrystalPrefix
				      + Cgoniometer::ms_sGonioVectors,
				      12, a12fVecs);
	  nStat = roHeader.nReplaceValue(Ccrystal::ms_sCrystalPrefix
					  + Cgoniometer::ms_sGonioVectors,
					  9, &a12fVecs[3]);
	  roHeader.nWrite("NTHU1.head");
	}
      else
	{
	  cout << "ERROR in NTHUfixup, Omega value is NOT 0! ABORT!\n" << flush;
	}
    }
  delete poCrysGonio;
  return (nStat);
}

int DTREK_EXPORT nNTHUfixup(Cstring& rsOrientWish, Cimage_header& roHeaderOrig,
	   Cimage_header& roHeader, 
	   Ccrystal *poCrystal)
{
  float a3fHKL[3];
  float a3fX[3];
  float fDstar;
  Cstring sTemp;
  int nStat = 0;
  float a3x3fBMat[3][3];
  sTemp = rsOrientWish.after('(');
  sTemp = sTemp.before(')');
  //cout << "HKL: " << sTemp << endl;
  nStat = sscanf(sTemp.string(), "%f %f %f", 
		 &a3fHKL[0], &a3fHKL[1], &a3fHKL[2]);
  if (3 != nStat)
    {
      return (1);
    }

  //   d* = Bh    |d*| = 1/|d|
  //   --   =-     --       -
  //
  //   Theta depends on wavelength
  //   lambda = 2d sin( theta )
  //   theta  = asin(lambda / 2d)

  Csource oSource(roHeaderOrig);

  poCrystal->nCalcBMatrix();  // Wavelength is 1 Angstrom;
  poCrystal->vGetBMatrix(&a3x3fBMat[0][0]);
  vMulMat3DVec3D(a3x3fBMat, a3fHKL, a3fX);
  fDstar = fLenVec3D(a3fX);

  float f2Theta;
  float fWave;
  fWave = oSource.fGetWavelength();

  f2Theta = 2.0 * asin( (double) (fWave * fDstar * 0.5) );
  f2Theta = f2Theta / Gs_dRADIANS_PER_DEGREE;

  cout << "\nNTHU info"
       << "\n========="
       << "\nWavelength is: " << fWave
       << "\n2theta for refln (" 
       << a3fHKL[0] << ' ' << a3fHKL[1] << ' ' << a3fHKL[2] << ") is "
       << f2Theta << endl;
	  

  float a4fCrysGonValues[4];

  nStat = roHeader.nGetValue(Ccrystal::ms_sCrystalPrefix
			     + Cgoniometer::ms_sGonioValues,
			     3, &a4fCrysGonValues[1]);
  a4fCrysGonValues[0] = f2Theta;
  nStat = roHeaderOrig.nReplaceValue(Ccrystal::ms_sCrystalPrefix
			     + Cgoniometer::ms_sGonioValues,
			     4, a4fCrysGonValues);

  float a5fDatum[5];

  nStat = roHeaderOrig.nGetValue(Cstring(D_K_ScanCrysDatum),
				 5, &a5fDatum[0]);
  int i;
  for (i = 0; i < 4; i++)
    a5fDatum[i+1] = a4fCrysGonValues[i];
  a5fDatum[0] = 4;
  nStat = roHeaderOrig.nReplaceValue(Cstring(D_K_ScanCrysDatum),
				     5, &a5fDatum[0]);

  
  Cgoniometer *poCrysGonio;
  poCrysGonio   = new Cgoniometer(roHeaderOrig, Ccrystal::ms_sCrystalPrefix);
  poCrysGonio->nSetDatum(4, &a4fCrysGonValues[0]);
  //poCrysGonio->nList();
/*
//For debugging:    
  for (i = 3; i >= 0; i--)
    {
      poCrysGonio->nGetRotVector(i, (float*)&a5fDatum[0]);
      cout << "For axis " << i << " rot vector is "
	   << a5fDatum[0] << ", " << a5fDatum[1] << ", " << a5fDatum[2] << endl;
    }
*/
  // It is the Psi axis!

  poCrysGonio->nGetRotVector(1, (float*)&a5fDatum[0]);
  roHeaderOrig.nReplaceValue(Cstring(D_K_RotVector), 3, a5fDatum);
  roHeaderOrig.nReplaceValue(Cstring(D_K_ScanPrefix D_K_RotVector), 3, a5fDatum);
  roHeader = roHeaderOrig;

  delete poCrysGonio;
  return (nStat);
}
