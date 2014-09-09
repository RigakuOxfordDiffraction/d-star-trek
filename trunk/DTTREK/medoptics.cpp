//
// Copyright (c) 1998 Molecular Structure Corporation
//
#include <string.h>
#include "dtreksys.h"
#include "medoptics.h"
#include "Cimage_header.h"
#include "dskio.h"
#include "dtrekdefs.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

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

//
// WARNING // WARNING // WARNING // WARNING // WARNING // WARNING // WARNING //
// The following code has NOT been test on a 64-bit compilation
//

int 
nReadMedOpticsHeader(const int* pnFile, const Cstring& rsFilename,
		   char *pcBuffer, Cimage_header *poHeader)
{
  // Read a no-header MedOptics image that may or may not have been 
  // partially read already.
  // Convert the header to a d*TREK header.
  //
  // Use the dskio routines for i/o.
  // The file to read must have been opened with the dskbor routine already.
  // Arguments:
  //   pnFile    pointer to open file tag for dskio routines.
  //   pcBuffer  pointer to first 512-bytes of header if header partially read
  //                  or NULL if header not partially read.
  //   poHeader  pointer to image header to modify.
  //   rsFilename image filename; the filename contains some sequence info
  //
  // Returns
  //     0 success
  // not 0 failure
  //

  int   nStat;       // Local status flag
  int   nLen;        // Temp variable for lengths
  float a10fTemp[10];
  int   a10nTemp[10];
  int   nDim0, nDim1;

  static Cstring sMedO = "MEDOPT_";
  
//  cout << "in nReadMedOpticsHeader" << endl;

  nStat = 0;
  nLen  = lFileGetSize(rsFilename);
  if (10 + 512 * 512 * 2 == nLen)
    {
      nDim0 = 512;
      nDim1 = 512;
    }
  else if (10 + 1024 * 1024 * 2 == nLen)
    {
      nDim0 = 1024;
      nDim1 = 1024;
    }
  else if (45785088 == nLen)
    {
      // For SPring8 images (they're big!)

      nDim0 = 4968;
      nDim1 = 4608;
    }

  a10nTemp[0] = nDim0;
  a10nTemp[1] = nDim1;
  
  poHeader->nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_MEDOPTICS));

  (void) poHeader->nReplaceValue(Cstring (D_K_Type), Cstring("mad"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Dim), Cstring("2"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Size1), a10nTemp[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_Size2), a10nTemp[1]);
  (void) poHeader->nReplaceValue(Cstring (D_K_DataType), 
				 Cstring("unsigned short int"));
  (void) poHeader->nReplaceValue(Cstring (D_K_DtdisplayOrientation), 
				 Cstring("-X+Y"));


  poHeader->nReplaceValue(Cstring (D_K_ByteOrder), 
			  Cstring (D_K_LittleEndian));

  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNumber), (int)1);
  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNames), sMedO);
  (void) poHeader->nReplaceValue(sMedO + D_K_DetectorDimensions, 2, a10nTemp);

  a10fTemp[0] = 0.100 * a10nTemp[0];
  a10fTemp[1] = 0.100 * a10nTemp[1];
  (void) poHeader->nReplaceValue(sMedO + D_K_DetectorSize, 2, a10fTemp, 5);
  (void) poHeader->nReplaceValue(sMedO + D_K_DetectorDescription,
				 Cstring ("MedOptics conversion"));

  // The following may not be correct:

  (void) poHeader->nReplaceValue(sMedO + D_K_DetectorVectors,
				 Cstring ("1 0 0 0 1 0"));

  (void) poHeader->nReplaceValue(sMedO + D_K_GonioNumValues,
				 (int) 6);
  (void) poHeader->nReplaceValue(sMedO + D_K_GonioUnits,
				 Cstring ("deg deg deg mm mm mm"));
  (void) poHeader->nReplaceValue(sMedO + D_K_GonioNames,
                 Cstring ("RotZ RotX/Swing RotY TransX TransY TransZ/Dist"));

  (void) poHeader->nReplaceValue(sMedO + D_K_GonioVectors,
		 Cstring ("0 0 1 -1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1"));

  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;   // Two theta
  a10fTemp[2] = 0.0;
  a10fTemp[3] = 0.0;
  a10fTemp[4] = 0.0;
  a10fTemp[5] = 100.0; // xtal_to_detector
      
  (void) poHeader->nReplaceValue(sMedO + D_K_GonioValues, 6, a10fTemp, 3);

  (void) poHeader->nReplaceValue(sMedO + D_K_NonunfType, 
				 Cstring(D_K_NonunfStateNone));

      
  (void) poHeader->nReplaceValue(sMedO + D_K_SpatialDistortionType,
				 Cstring (D_K_SpatialTypeSimple));
  a10fTemp[0] = 0.5 * (float)nDim0;
  a10fTemp[1] = 0.5 * (float)nDim1;
  a10fTemp[2] = 0.100f;
  a10fTemp[3] = 0.100f;
  (void) poHeader->nReplaceValue(sMedO + D_K_SpatialDistortionInfo,
				 4, a10fTemp, 5);
      
  (void) poHeader->nReplaceValue(sMedO + D_K_SpatialDistortionVectors,
				 Cstring("1 0 0 -1"));

  (void) poHeader->nReplaceValue(Cstring (D_K_SaturatedValue), 65535);
  a10fTemp[0] = 1.0;
  a10fTemp[1] = 1.54178f;

  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePrefix) + 
				 D_K_Wavelength, 2, a10fTemp);
  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
				 Cstring ("0.50 1 0 0"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceSpectralDispersion),
				 Cstring ("0.0002 0.0002"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceCrossfire),
				 Cstring ("0.0002 0.0002 0.0 0.0"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceSize),
				 Cstring ("0.0 0.0 0.0 0.0"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceVectors),
				 Cstring ("0 0 1 0 1 0 1 0 0"));
  (void) poHeader->nReplaceValue(Cstring (D_K_SourceValues),
				 Cstring ("0 0"));

  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + 
				 D_K_GonioNumValues, (int) 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioNames,
				 Cstring ("Omega Chi Phi"));
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioUnits,
				 Cstring ("deg deg deg"));
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
				 Cstring ("1 0 0 0 0 1 1 0 0"));
  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
				 3, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + 
				 D_K_GonioDescription,
				 Cstring ("Single w/2-arc goniometer head"));
  a10fTemp[0] = 0.0;   // (float)tDIP2030_header.start_phi * 0.001;
  a10fTemp[1] = 1.0;   // (float)tDIP2030_header.end_phi   * 0.001;
  a10fTemp[2] = a10fTemp[1] - a10fTemp[0];
  a10fTemp[3] = 10.0;  // (float)tDIP2030_header.exposure_time * 0.001;
  a10fTemp[4] = 1;
  a10fTemp[5] = 0.0;
  a10fTemp[6] = 0.0;
  a10fTemp[7] = 0.0;
  a10fTemp[8] = 0.0;
  a10fTemp[9] = 0.0;
  (void) poHeader->nReplaceValue(Cstring (D_K_Rotation), 10, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_Rotation,
				 10, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_RotAxisName),
				 Cstring ("Phi"));
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + 
				 D_K_RotAxisName,
				 Cstring ("Phi"));

  a10fTemp[0] = 1.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;

  (void) poHeader->nReplaceValue(Cstring (D_K_RotVector),
				 3, &a10fTemp[0], 4);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_RotVector,
				 3, &a10fTemp[0], 4);

  // Determine scan template and seq info from the filename

  // Change last 3 digits in rsFilename to ?s

  int a3nSeqInfo[3] = {0, 1, 0};
  Cstring sTemp;

  sTemp = sBuildScanTemplate(rsFilename, 3, &a3nSeqInfo[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanTemplate), sTemp);

  (void) poHeader->nReplaceValue(Cstring (D_K_ScanSeqInfo), 3, a3nSeqInfo);

  if ( (0 != nStat) && ("" != sGetEnv("DTREK_RAXIS_IGNORE")))
    {
      cout << "Warning in medoptics, error reading header, but ignored!\n";
      nStat = 0;
    }
  return (nStat);
}
