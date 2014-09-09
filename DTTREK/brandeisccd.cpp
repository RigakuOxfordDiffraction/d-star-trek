//
// Copyright (c) 1998 Molecular Structure Corporation
//
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
//+External prototypes

#include "dtrekdefs.h"
#include "Cimage_header.h"

int nEditBrandeisCCDHeader(const Cstring &rsFilename, Cimage_header *poHeader)
{
  float fTemp;
  float a10fTemp[10];
  int   a10nTemp[10];
  int   nTemp;
  Cstring sDet = "B1p2_";
  Cstring sTemp;

//  cout << "nEditBrandeisCCDHeader called!\n" << endl;
  poHeader->nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_BRANDEIS));

  poHeader->nGetValue("DETNAME", &sDet);

  poHeader->nReplaceValue(D_K_DetectorNames, sDet);
  poHeader->nReplaceValue(D_K_DetectorNumber, "1");

  poHeader->nReplaceValue(sDet + Cstring(D_K_GonioNames), 
			  "RotZ RotX/Swing RotY TransX TransY TransZ/Dist");
  poHeader->nReplaceValue(sDet + Cstring(D_K_GonioNumValues), "6"); 
  poHeader->nReplaceValue(sDet + Cstring(D_K_GonioUnits),
			  "deg deg deg mm mm mm");
  poHeader->nReplaceValue(sDet + Cstring(D_K_GonioVectors), 
//                        "0 0 1 -1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1");
  // Brandeis detectors at X25 and X12C have the theta angle going the other
  // direction 17-Apr-2000
			  "0 0 1  1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1");

  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;
  a10fTemp[3] = 0.0;
  a10fTemp[4] = 0.0;
  a10fTemp[5] = 100.0;
  poHeader->nGetValue("DIST", &a10fTemp[5]);
  poHeader->nGetValue("THETA", &a10fTemp[1]);

  poHeader->nReplaceValue(sDet + Cstring(D_K_GonioValues), 6, a10fTemp, 3);

  poHeader->nGetValue(D_K_Size1, &a10fTemp[6]);
  poHeader->nGetValue(D_K_Size2, &a10fTemp[7]);
  a10fTemp[8] = 0.0;
  a10fTemp[9] = 0.0;
  poHeader->nGetValue("PIXEL_SIZE", 2, &a10fTemp[8]);

  if (0.0 >= a10fTemp[8]) 
      a10fTemp[8] = 0.100f;
  if (0.0 >= a10fTemp[9]) 
      a10fTemp[9] = 0.100f;
  a10fTemp[0] = a10fTemp[7] * 0.5;
  a10fTemp[1] = a10fTemp[6] * 0.5;
  poHeader->nGetValue("BEAM_POSITION", 2, &a10fTemp[0]);
//+jwp 3-Nov-1999  
//swap beam x and y for the brandeis detector
  fTemp       = a10fTemp[0];
  a10fTemp[0] = a10fTemp[1];
  a10fTemp[1] = fTemp;
//-jwp 

  a10fTemp[2] = a10fTemp[8];
  a10fTemp[3] = a10fTemp[9];

  poHeader->nReplaceValue(sDet + Cstring(D_K_SpatialDistortionInfo),
			  4, a10fTemp, 5);
  poHeader->nReplaceValue(sDet + Cstring(D_K_SpatialBeamPosition), 
			  2, a10fTemp, 1);

  poHeader->nReplaceValue(sDet + Cstring(D_K_DetectorDescription), 
			  sDet + " conversion");
  a10nTemp[0] = 2048;
  a10nTemp[1] = 2048;
  poHeader->nGetValue(D_K_Size1, &a10nTemp[0]);
  poHeader->nGetValue(D_K_Size2, &a10nTemp[1]);
  poHeader->nReplaceValue(sDet + D_K_DetectorDimensions, 2, a10nTemp, 0);  
  a10fTemp[0] = a10fTemp[9] * (float)a10nTemp[0];
  a10fTemp[1] = a10fTemp[9] * (float)a10nTemp[1];
  poHeader->nReplaceValue(sDet + Cstring(D_K_DetectorSize), 2, a10fTemp, 2);
  poHeader->nReplaceValue(sDet + Cstring(D_K_DetectorVectors), "1 0 0 0 1 0");
  poHeader->nReplaceValue(sDet + D_K_NonunfType,
				 Cstring(D_K_NonunfStateSimpleMask));
  poHeader->nReplaceValue(sDet + D_K_NonunfInfo, 
				 "$(FirstScanImage)");
  poHeader->nReplaceValue(sDet + Cstring(D_K_SpatialDistortionType), 
			  D_K_SpatialTypeSimple);
  poHeader->nReplaceValue(sDet + Cstring(D_K_SpatialDistortionVectors),
			  "1 0 0 -1");

  nTemp = 60000;                                        // lowest precedence
  poHeader->nGetValue("SATURATED_VALUE", &nTemp);       // higher precedence
  poHeader->nReplaceValue(D_K_SaturatedValue, nTemp);

  poHeader->nReplaceValue(D_K_DtdisplayOrientation, "-X+Y");
  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;

  // Try to read crystal gonio values from existing keyword.  If keyword does
  // not exist, then the values remain unchanged so this is safe

  poHeader->nGetValue("DATUM", 3, a10fTemp);

  poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioValues, 3, a10fTemp, 3);
  poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioNumValues, 3);
  // Should use keyword AXIS to figure out rest
  poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioNames, "Omega Kappa Phi");
  poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioUnits, "deg deg deg");
  poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioDescription, 
			  "Kappa 3-circle");
  poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
			  "1 0 0 0.6427876 0 -0.7660444 1 0 0");
//                                sin(40)    -cos(40)
//                                cos(50)    -sin(50)
  // Construct scan

  a10fTemp[0] = 0.0;   // Rot Start default
  a10fTemp[1] = 0.5;   // Rot End   default
  a10fTemp[2] = 0.5;   // Rot Incr  default
  a10fTemp[3] = 5.0;   // Exp time  default
  a10fTemp[4] = 1.0;   // Num oscillations
  a10fTemp[5] = 0.0;   // ?
  a10fTemp[6] = 0.0;   // ?
  a10fTemp[7] = 0.0;   // ?
  a10fTemp[8] = 0.0;   // ?
  a10fTemp[9] = 0.0;   // ?

  poHeader->nGetValue("PHIST", &a10fTemp[0]);
  poHeader->nGetValue("PHIINC", &a10fTemp[2]);
  a10fTemp[1] = a10fTemp[0] + a10fTemp[2];
  poHeader->nGetValue("TIME", &a10fTemp[3]);
  poHeader->nReplaceValue(D_K_ScanPrefix D_K_Rotation, 10, a10fTemp, 4);
  poHeader->nReplaceValue(               D_K_Rotation, 10, a10fTemp, 4);

  //  Default axis is Phi ... in case ROT_AXIS does not exit.
  sTemp = "Phi";
  poHeader->nGetValue("ROT_AXIS", &sTemp);
  sTemp.upcase();

  a10fTemp[0] = 1.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;
  if ("OMEGA" != sTemp)
    {
      // Need to build rot_vector
    }

  poHeader->nReplaceValue(D_K_ScanPrefix D_K_RotAxisName, sTemp);
  poHeader->nReplaceValue(               D_K_RotAxisName, sTemp);
  poHeader->nReplaceValue(D_K_ScanPrefix D_K_RotVector, 3, a10fTemp, 3);
  poHeader->nReplaceValue(               D_K_RotVector, 3, a10fTemp, 3);
  poHeader->nReplaceValue("SCAN_FIELDS", "RotStart RotEnd RotInc RotTime");
  poHeader->nReplaceValue("ROTATION_FIELDS", "RotStart RotEnd RotInc RotTime");

/////////
  // Determine scan template and seq info from the filename 
 
  // Change last 3 digits in rsFilename to ?s 

  int     a3nSeqInfo[3] = {0, 1, 0}; 
  sTemp = rsFilename; 
  sTemp = sBuildScanTemplate(rsFilename, 3, &a3nSeqInfo[0]); 
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanTemplate), sTemp); 
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanSeqInfo), 3, a3nSeqInfo); 
/////////
 
  poHeader->nReplaceValue(D_K_SourceIntensity, "1.0");
  a10fTemp[1] = 1.00;
  poHeader->nGetValue("WAVE", &a10fTemp[1]);
  a10fTemp[0] = 1.0;

  poHeader->nReplaceValue(D_K_SourceWavelength, 2, a10fTemp, 5);
  poHeader->nReplaceValue(D_K_SourceCrossfire, "0.0002 0.0002 0.0 0.0");
  poHeader->nReplaceValue(D_K_SourcePolarz, "0.95 0 1 0");
  poHeader->nReplaceValue(D_K_SourceSize, "0.0 0.0 0.0 0.0");
  poHeader->nReplaceValue(D_K_SourceSpectralDispersion, "0.0002 0.0002");
  poHeader->nReplaceValue(D_K_SourceVectors, "0 0 1 0 1 0 1 0 0");
  poHeader->nReplaceValue(Cstring (D_K_SourceValues),
			  Cstring ("0 0"));

  return (0);
}
