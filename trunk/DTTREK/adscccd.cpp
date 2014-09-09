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

int nEditADSCCCDHeader(const Cstring &rsFilename, Cimage_header *poHeader)
{
  float a10fTemp[10];
  int   a10nTemp[10];
  int   nTemp;
  int   nSerialNumber;
  Cstring sDet = "A4_";

  //  cout << "nEditADSCCDHeader called!\n" << endl;

  poHeader->nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_ADSC));

  poHeader->nReplaceValue(D_K_DetectorNames, sDet);
  poHeader->nReplaceValue(D_K_DetectorNumber, "1");
  poHeader->nReplaceValue(sDet + D_K_GonioNames, 
			  "RotZ RotX/Swing RotY TransX TransY TransZ/Dist");
  poHeader->nReplaceValue(sDet + D_K_GonioNumValues, "6"); 
  poHeader->nReplaceValue(sDet + D_K_GonioUnits,
			  "deg deg deg mm mm mm");
  poHeader->nReplaceValue(sDet + D_K_GonioVectors, 
			  "0 0 1 -1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1");

  nSerialNumber = 0;
  (void) poHeader->nGetValue("DETECTOR_SN", &nSerialNumber);
  if ( (416 == nSerialNumber) || (443 == nSerialNumber) )
    {
      // It's the X26C detector or the IMCA-CAT ID17 one
      poHeader->nReplaceValue(sDet + D_K_GonioVectors, 
			      "0 0 1 1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1");
    }
  
  a10nTemp[0] = 2304;
  a10nTemp[1] = 2304;
  poHeader->nGetValue(D_K_Size1, &a10nTemp[0]);
  poHeader->nGetValue(D_K_Size2, &a10nTemp[1]);
  poHeader->nReplaceValue(sDet + D_K_DetectorDimensions, 2, a10nTemp, 0);  

  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;
  a10fTemp[3] = 0.0;
  a10fTemp[4] = 0.0;
  a10fTemp[5] = 100.0;
  poHeader->nGetValue("DISTANCE", &a10fTemp[5]);
  poHeader->nGetValue("TWOTHETA", &a10fTemp[1]);
  poHeader->nReplaceValue(sDet + D_K_GonioValues, 6, a10fTemp, 3);

  poHeader->nGetValue("BEAM_CENTER_X", &a10fTemp[7]);
  poHeader->nGetValue("BEAM_CENTER_Y", &a10fTemp[8]);
  poHeader->nGetValue("PIXEL_SIZE", &a10fTemp[9]);

  if (0.0 >= a10fTemp[9]) 
      a10fTemp[9] = 0.0816f;

  // Convert beam x, y into pixel values.

  a10fTemp[0] = a10fTemp[8] / a10fTemp[9];
  a10fTemp[1] = a10fTemp[7] / a10fTemp[9];
  if (908 == nSerialNumber)
    {
      // Detector with serial number 908 has X and Y in the image header
      // switched from other ADSC detectors ...

      // Switch beam x and y
      a10fTemp[0] = a10fTemp[7] / a10fTemp[9];
      a10fTemp[1] = a10nTemp[1] - (a10fTemp[8] / a10fTemp[9]);
    }

  a10fTemp[2] = a10fTemp[9];
  a10fTemp[3] = a10fTemp[9];
  poHeader->nReplaceValue(sDet + D_K_SpatialBeamPosition, 2, a10fTemp, 5);
  poHeader->nReplaceValue(sDet + D_K_SpatialDistortionInfo,4, a10fTemp, 5);

  poHeader->nReplaceValue(sDet + D_K_DetectorDescription, sDet + " conversion");

  a10fTemp[0] = a10fTemp[9] * (float)a10nTemp[0];
  a10fTemp[1] = a10fTemp[9] * (float)a10nTemp[1];
  poHeader->nReplaceValue(sDet + D_K_DetectorSize, 2, a10fTemp, 2);
  poHeader->nReplaceValue(sDet + D_K_DetectorVectors, "1 0 0 0 1 0");
  poHeader->nReplaceValue(sDet + D_K_NonunfType,
				 Cstring(D_K_NonunfStateSimpleMask));
  poHeader->nReplaceValue(sDet + D_K_NonunfInfo, 
				 "$(FirstScanImage)");

  poHeader->nReplaceValue(sDet + D_K_SpatialDistortionType, 
			  D_K_SpatialTypeSimple);
  poHeader->nReplaceValue(sDet + D_K_SpatialDistortionVectors, "1 0 0 -1");

  nTemp = 65535;
  poHeader->nGetValue("CCD_IMAGE_SATURATION", &nTemp);
  poHeader->nReplaceValue(D_K_SaturatedValue, nTemp);

  poHeader->nReplaceValue(D_K_DtdisplayOrientation, "-X+Y");

  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;

  nTemp = poHeader->nGetValue(D_K_CrystalPrefix D_K_GonioValues, 3, &a10fTemp[3]);

//+27-Apr-2010 jwp

  if (0 == nTemp)
    {
      // Found a CRYSTAL_GONIO_VALUES keyword to use

      a10fTemp[0] = a10fTemp[3];
      a10fTemp[1] = a10fTemp[4];
      a10fTemp[2] = a10fTemp[5];
    }
  else
    {
      // If the keyword DATUM exists and has 3 values,

      nTemp = poHeader->nGetValue("DATUM", 3, &a10fTemp[3]);
      if (0 == nTemp)
	{
	  a10fTemp[0] = a10fTemp[3];
	  a10fTemp[1] = a10fTemp[4];
	  a10fTemp[2] = a10fTemp[5];
	}
    }
//-27-Apr-2010 jwp

//+12-Jul-2000 jwp
  // Make crystal_gonio_value for omega==0, not the starting phi value!
  //  poHeader->nGetValue("PHI", &a10fTemp[0]);
//-12-Jul-2000 jwp

  poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioValues, 3, a10fTemp, 3);
  poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioNumValues, 3);
  // Should use keyword AXIS to figure out rest

  if (903 != nSerialNumber)
  {
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioNames, "Omega Chi Phi");
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioUnits, "deg deg deg");
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioDescription, 
			      "Eulerian 3-circle");

      // WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!
      // See lines way below where these are really set again
      // WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!

      if (   (414 == nSerialNumber)
          || (406 == nSerialNumber) 
          || (457 == nSerialNumber) // Australia Light Source
          || (928 == nSerialNumber) // Australia Light Source
          ) 
          {
            poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
                                    "-1 0 0 0 1 0 1 0 0");
          }
      else if (415 == nSerialNumber)
          {
            poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
                                    "-1 0 0 0 1 0 -1 0 0");
          }
      else if ( (914 == nSerialNumber) || (909 == nSerialNumber) )
          {
            poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
                                    "-1 0 0 0 1 0 -1 0 0");
          }
          //mrp begin 4992 bug fix.  I am guessing about omega.
      else if (471 == nSerialNumber)
        {
          poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors,
                                  "-1 0 0 0 1 0 -1 0 0");
        }
          //mrp end
      else
      {
         poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
				  "1 0 0 0 1 0 1 0 0");
      }
  }
  else
  {
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioNames, "Omega Kappa Phi");
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioUnits, "deg deg deg");
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioDescription, 
			      "Crystal Logic Kappa orienter for Q315");
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
			      "1.0 0.0 0.0 0.64279 0.0 0.76604 1.0 0.0 0.0");
  }

  // Construct scan

  a10fTemp[0] = 0.0;   // Rot Start default
  a10fTemp[1] = 1.0;   // Rot End   default
  a10fTemp[2] = 1.0;   // Rot Incr  default
  a10fTemp[3] = 10.0;  // Exp time  default
  a10fTemp[4] = 1.0;   // Num oscillations
  a10fTemp[5] = 0.0;   // ?
  a10fTemp[6] = 0.0;   // ?
  a10fTemp[7] = 0.0;   // ?
  a10fTemp[8] = 0.0;   // ?
  a10fTemp[9] = 0.0;   // ?

  poHeader->nGetValue("OSC_START", &a10fTemp[0]);
  poHeader->nGetValue("OSC_RANGE", &a10fTemp[2]);
  if  (   
           (413 == nSerialNumber)
        || (420 == nSerialNumber)
        || (428 == nSerialNumber) )
    {
      // ESRF detectors put the oscillation start into the PHI keyword.
      // 2004-Aug-04 OK, there is a problem here there can be 
      // OSC_START, OSC_RANGE, OMEGA and PHI keywords.
      // If the rotation axis is phi and OSC_START equals PHI, then do NOT
      // add OSC_START and PHI to get the rotation_start.
      // However, if they are different, then add OSC_START and PHI

      a10fTemp[9] = 0.0;
      poHeader->nGetValue("PHI", &a10fTemp[9]);
      if (a10fTemp[0] != a10fTemp[9])
	a10fTemp[0] += a10fTemp[9];
    }

  poHeader->nReplaceValue(sDet + "OSC_START", a10fTemp[0]);
  poHeader->nReplaceValue(sDet + "OSC_RANGE", a10fTemp[2]);
  poHeader->nDelete("OSC_START");
  poHeader->nDelete("OSC_RANGE");

  a10fTemp[1] = a10fTemp[0] + a10fTemp[2];
  poHeader->nGetValue("TIME", &a10fTemp[3]);
  poHeader->nReplaceValue(D_K_ScanPrefix D_K_Rotation, 10, a10fTemp, 4);
  poHeader->nReplaceValue(               D_K_Rotation, 10, a10fTemp, 4);
  if  (    (443 == nSerialNumber)
        || (413 == nSerialNumber)
        || (420 == nSerialNumber)
        || (428 == nSerialNumber) )
    {
      poHeader->nReplaceValue(D_K_ScanPrefix D_K_RotAxisName, "Phi");
      poHeader->nReplaceValue(               D_K_RotAxisName, "Phi");
    }
  else
    {
      poHeader->nReplaceValue(D_K_ScanPrefix D_K_RotAxisName, "Omega");
      poHeader->nReplaceValue(               D_K_RotAxisName, "Omega");
    }
  if (418 == nSerialNumber)
    {
      //+jwp 16-April-2008  The X4A-C goniometers changed
      // WARNING: I don't know about 2theta swing angle yet

      poHeader->nReplaceValue(D_K_ScanPrefix D_K_RotVector, "1 0 0");
      poHeader->nReplaceValue(               D_K_RotVector, "1 0 0");
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
			      "1 0 0 0 0 1 0 1 0"); // Just a guess
    }
  else if (914 == nSerialNumber)
    {
      poHeader->nReplaceValue(D_K_ScanPrefix D_K_RotVector, "-1 0 0");
      poHeader->nReplaceValue(               D_K_RotVector, "-1 0 0");
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
			      "-1 0 0 0 1 0 -1 0 0");
    }
  else if (   (414 == nSerialNumber) || (415 == nSerialNumber)
           || (471 == nSerialNumber)
           || (406 == nSerialNumber)
           || (457 == nSerialNumber) // Australia Light Source
           || (928 == nSerialNumber) // Australia Light Source
          ) 
  {
      // CHESS detector rotates the other way 
      //  and so does detector 415 at Bio-CARS

      poHeader->nReplaceValue(D_K_ScanPrefix D_K_RotVector, "-1 0 0");
      poHeader->nReplaceValue(               D_K_RotVector, "-1 0 0");
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
			      "-1 0 0 0 1 0 1 0 0");
    }
  else
    {
      poHeader->nReplaceValue(D_K_ScanPrefix D_K_RotVector, "1 0 0");
      poHeader->nReplaceValue(               D_K_RotVector, "1 0 0");
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
			      "1 0 0 0 1 0 1 0 0");
    }
  poHeader->nReplaceValue("SCAN_FIELDS", "RotStart RotEnd RotInc RotTime");
  poHeader->nReplaceValue("ROTATION_FIELDS", "RotStart RotEnd RotInc RotTime");

/////////
  // Determine scan template and seq info from the filename 
 
  // Change last 4 digits in rsFilename to ?s 

  int     a3nSeqInfo[3] = {0, 1, 0}; 
  Cstring sTemp;

  sTemp = sBuildScanTemplate(rsFilename, 4, &a3nSeqInfo[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanTemplate), sTemp); 
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanSeqInfo), 3, a3nSeqInfo); 
/////////
 
  poHeader->nReplaceValue(D_K_SourceIntensity, "1.0");
  a10fTemp[1] = 1.00;
  poHeader->nGetValue("WAVELENGTH", &a10fTemp[1]);
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
