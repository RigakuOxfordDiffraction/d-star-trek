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

//
// WARNING // WARNING // WARNING // WARNING // WARNING // WARNING // WARNING //
// The following code has NOT been test on a 64-bit compilation
//

#include <string.h>
#include "dtreksys.h"
#include "dip2030.h"
#include "raxis.h"
#include "Cimage_header.h"
#include "dskio.h"
#include "dtrekdefs.h"
#include "dtrekvec.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

int 
nReadDIP2030Header(const int* pnFile, const long lFileSize, 
		   const Cstring& rsFilename,
		   char *pcBuffer, Cimage_header *poHeader)
{
  // Build a temporary DIP2030 header.  For DIP2030 images, the image info is
  // found in a trailer that is read in nReadDIP2030Trailer() after the
  // binary pixel info is read in.
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
  // The DIP2030 header is always 1024 bytes long.

  int   nStat;       // Local status flag
  int   nLen;        // Temp variable for lengths
  float a10fTemp[10];
  int   a10nTemp[10];

  static Cstring sDIP2030 = "DIP2030_";
  
  //  cout << "about to build nReadDIP2030Header" << endl;
  nStat = 0;
  nLen  = sizeof(tagDIP2030_header);
  if (NULL != pcBuffer)
    {
      nLen = nLen - 512;
    }

  if (18001024 == lFileSize)
    {
      a10nTemp[0] = 3000;
      a10nTemp[1] = 3000;
    }
  else if (32001024 == lFileSize)
    {
      a10nTemp[0] = 4000;
      a10nTemp[1] = 4000;
    }
  else
    {
      cout << "ERROR with filesize (" << lFileSize 
           << ") in nReadDIP20??Header!\n\n";
      return (-1);
    }
    
   poHeader->nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_DIP));

  (void) poHeader->nReplaceValue(Cstring (D_K_Type), Cstring("mad"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Dim), Cstring("2"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Size1), a10nTemp[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_Size2), a10nTemp[1]);
  (void) poHeader->nReplaceValue(Cstring (D_K_DataType), 
				 Cstring("unsigned short int"));
  (void) poHeader->nReplaceValue(Cstring (D_K_DtdisplayOrientation), 
				 Cstring("-X+Y"));


  poHeader->nReplaceValue(Cstring (D_K_ByteOrder), 
			  Cstring (D_K_BigEndian));

  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNumber), (int)1);
  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNames), sDIP2030);
  (void) poHeader->nReplaceValue(sDIP2030 + D_K_DetectorDimensions, 2, a10nTemp);

  a10fTemp[0] = 0.100 * (float)a10nTemp[0];
  a10fTemp[1] = 0.100 * (float)a10nTemp[1];
  (void) poHeader->nReplaceValue(sDIP2030 + D_K_DetectorSize, 2, a10fTemp, 5);
  (void) poHeader->nReplaceValue(sDIP2030 + D_K_DetectorDescription,
				 Cstring ("DIP2030 conversion"));

  // The following may not be correct:

  (void) poHeader->nReplaceValue(sDIP2030 + D_K_DetectorVectors,
				 Cstring ("1 0 0 0 1 0"));

  (void) poHeader->nReplaceValue(sDIP2030 + D_K_GonioNumValues,
				 (int) 6);
  (void) poHeader->nReplaceValue(sDIP2030 + D_K_GonioUnits,
				 Cstring ("deg deg deg mm mm mm"));
  (void) poHeader->nReplaceValue(sDIP2030 + D_K_GonioNames,
                 Cstring ("RotZ RotX/Swing RotY TransX TransY TransZ/Dist"));

  (void) poHeader->nReplaceValue(sDIP2030 + D_K_GonioVectors,
		 Cstring ("0 0 1 -1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1"));

  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;   // Two theta
  a10fTemp[2] = 0.0;
  a10fTemp[3] = 0.0;
  a10fTemp[4] = 0.0;
  a10fTemp[5] = 100.0; // xtal_to_detector
      
  (void) poHeader->nReplaceValue(sDIP2030 + D_K_GonioValues, 6, a10fTemp, 3);

  (void) poHeader->nReplaceValue(sDIP2030 + D_K_NonunfType, 
				 Cstring(D_K_NonunfStateNone));

      
  (void) poHeader->nReplaceValue(sDIP2030 + D_K_SpatialDistortionType,
				 Cstring (D_K_SpatialTypeSimple));
  a10fTemp[0] = 0.5 * (float)a10nTemp[0];
  a10fTemp[1] = 0.5 * (float)a10nTemp[1];
  a10fTemp[2] = 0.100f;
  a10fTemp[3] = 0.100f;
  (void) poHeader->nReplaceValue(sDIP2030 + D_K_SpatialDistortionInfo,
				 4, a10fTemp, 5);
      
  (void) poHeader->nReplaceValue(sDIP2030 + D_K_SpatialDistortionVectors,
				 Cstring("1 0 0 -1"));

  (void) poHeader->nReplaceValue(Cstring (D_K_SaturatedValue), 500000);
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
				 Cstring ("Omega Kappa Phi"));
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioUnits,
				 Cstring ("deg deg deg"));
  a10fTemp[0] =  1.0f;
  a10fTemp[1] =  0.0f;
  a10fTemp[2] =  0.0f;

  a10fTemp[3] = -cos(50.0 * Gs_dRADIANS_PER_DEGREE);
  a10fTemp[4] =  0.0f;
  a10fTemp[5] =  sin(50.0 * Gs_dRADIANS_PER_DEGREE);

  a10fTemp[6] =  1.0f;
  a10fTemp[7] =  0.0f;
  a10fTemp[8] =  0.0f;

  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
				 9, a10fTemp, 4);

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
  a10fTemp[3] = 300.0; // (float)tDIP2030_header.exposure_time * 0.001;
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

  Cstring sTemp;

  int a3nSeqInfo[3] = {0, 1, 0};

  sTemp = sBuildScanTemplate(rsFilename, 3, &a3nSeqInfo[0]);

  (void) poHeader->nReplaceValue(Cstring (D_K_ScanTemplate), sTemp);

  (void) poHeader->nReplaceValue(Cstring (D_K_ScanSeqInfo), 3, a3nSeqInfo);

  if ( (0 != nStat) && ("" != sGetEnv("DTREK_RAXIS_IGNORE")))
    {
      cout << "Warning in dip2030, error reading header, but ignored!\n";
      nStat = 0;
    }
  return (nStat);
}

int nReadDIP2030Trailer(const int nFile, Cimage_header *poHeader)
{
  static Cstring sDIP2030 = "DIP2030_";
  int nStat;
  int nFileL;
  int nBytes;
  tagDIP2030_header tDIP2030_header;
  float a10fTemp[10];

  nBytes = sizeof(tDIP2030_header);
  nFileL = nFile;
  (void) dskbr(&nFileL, (char *)&tDIP2030_header, &nBytes, &nStat);

  // Mark trailer as read

  poHeader->nReplaceValue(sDIP2030 + "TRAILER", (int)1);

  if (0 != nStat)
    {
      cout << "ERROR reading DIP2030 trailer." << endl;
    }
  else
    {
      // Parse header and modify input header keywords

      int   *pnHead;
      float *pfHead;
      char  *pcHead;

      pnHead = (int *)  &tDIP2030_header;
      pcHead = (char *) &tDIP2030_header;
      pfHead = (float *)&tDIP2030_header;
//      printf("id: %s    monochromator:  %s  comment:\n%s\n",
//	     &pcHead[0], &pcHead[128], &pcHead[40]);
      //      cout << "dtype is: " << tDIP2030_header.tPart1.dtype << endl;
      if ( (0 >= pnHead[4]) || (500 <= pnHead[4]) )  // Could look at byte-order
	{
	  // Swap some binary values


          tDIP2030_header.tPart1.m_type1 = nSwapLong(tDIP2030_header.tPart1.m_type1);
          tDIP2030_header.tPart1.m_type2 = nSwapLong(tDIP2030_header.tPart1.m_type2);
          tDIP2030_header.tPart1.dtype = nSwapLong(tDIP2030_header.tPart1.dtype);
          tDIP2030_header.tPart1.pixelsize = nSwapLong(tDIP2030_header.tPart1.pixelsize);
          tDIP2030_header.tPart1.pixelsize2 = nSwapLong(tDIP2030_header.tPart1.pixelsize2);
          tDIP2030_header.tPart1.xsize = nSwapLong(tDIP2030_header.tPart1.xsize);
          tDIP2030_header.tPart1.ysize = nSwapLong(tDIP2030_header.tPart1.ysize);
          tDIP2030_header.tPart1.ipno = nSwapLong(tDIP2030_header.tPart1.ipno);
          tDIP2030_header.tPart1.ipx = nSwapLong(tDIP2030_header.tPart1.ipx);
          tDIP2030_header.tPart1.ipy = nSwapLong(tDIP2030_header.tPart1.ipy);
          tDIP2030_header.tPart1.repet = nSwapLong(tDIP2030_header.tPart1.repet);
          tDIP2030_header.tPart1.osc_axis = nSwapLong(tDIP2030_header.tPart1.osc_axis);
          tDIP2030_header.tPart1.xstart = nSwapLong(tDIP2030_header.tPart1.xstart);
          tDIP2030_header.tPart1.ystart = nSwapLong(tDIP2030_header.tPart1.ystart);

          tDIP2030_header.tPart1.radius = fSwapFloat(tDIP2030_header.tPart1.radius);
          tDIP2030_header.tPart1.x_lamda = fSwapFloat(tDIP2030_header.tPart1.x_lamda);
          tDIP2030_header.tPart1.cdist = fSwapFloat(tDIP2030_header.tPart1.cdist);
          tDIP2030_header.tPart1.pttheta = fSwapFloat(tDIP2030_header.tPart1.pttheta);
          tDIP2030_header.tPart1.exposure = fSwapFloat(tDIP2030_header.tPart1.exposure);
          tDIP2030_header.tPart1.kv = fSwapFloat(tDIP2030_header.tPart1.kv);
          tDIP2030_header.tPart1.ma = fSwapFloat(tDIP2030_header.tPart1.ma);
          tDIP2030_header.tPart1.collimator = fSwapFloat(tDIP2030_header.tPart1.collimator);
          tDIP2030_header.tPart1.coupling = fSwapFloat(tDIP2030_header.tPart1.coupling);
          tDIP2030_header.tPart1.phi1 = fSwapFloat(tDIP2030_header.tPart1.phi1);
          tDIP2030_header.tPart1.phi2 = fSwapFloat(tDIP2030_header.tPart1.phi2);
          tDIP2030_header.tPart1.phispeed = fSwapFloat(tDIP2030_header.tPart1.phispeed);
          tDIP2030_header.tPart1.g_omega = fSwapFloat(tDIP2030_header.tPart1.g_omega);
          tDIP2030_header.tPart1.g_kappa = fSwapFloat(tDIP2030_header.tPart1.g_kappa);
          tDIP2030_header.tPart1.g_phi = fSwapFloat(tDIP2030_header.tPart1.g_phi);



          /*
	  pnHead[1]  = nSwapLong(pnHead[1]);
	  pnHead[4]  = nSwapLong(pnHead[4]);
	  pnHead[5]  = nSwapLong(pnHead[5]);
	  pnHead[7]  = nSwapLong(pnHead[7]);
	  pnHead[8]  = nSwapLong(pnHead[8]);
	  pnHead[9]  = nSwapLong(pnHead[9]);
	  pnHead[41] = nSwapLong(pnHead[41]);
	  pnHead[42] = nSwapLong(pnHead[42]);
	  pnHead[51] = nSwapLong(pnHead[51]);


	  pfHead[30] = fSwapFloat(pfHead[30]);
	  pfHead[31] = fSwapFloat(pfHead[31]);
	  pfHead[43] = fSwapFloat(pfHead[43]);
	  pfHead[46] = fSwapFloat(pfHead[46]);
	  pfHead[48] = fSwapFloat(pfHead[48]);
	  pfHead[49] = fSwapFloat(pfHead[49]);
	  pfHead[50] = fSwapFloat(pfHead[50]);
	  pfHead[53] = fSwapFloat(pfHead[53]);
	  pfHead[60] = fSwapFloat(pfHead[60]);
	  tDIP2030_header.tPart1.g_omega = fSwapFloat(tDIP2030_header.tPart1.g_omega);
	  tDIP2030_header.tPart1.g_kappa = fSwapFloat(tDIP2030_header.tPart1.g_kappa);
	  tDIP2030_header.tPart1.g_phi   = fSwapFloat(tDIP2030_header.tPart1.g_phi);
	  tDIP2030_header.tPart1.dtype   = nSwapLong(tDIP2030_header.tPart1.dtype);
	  tDIP2030_header.tPart1.dtype   = nSwapLong(tDIP2030_header.tPart1.dtype);
	  tDIP2030_header.tPart1.osc_axis= nSwapLong(tDIP2030_header.tPart1.osc_axis);
      */
	}

      if (1 == tDIP2030_header.tPart1.osc_axis)
	{
	  // Rotation axis vector is along -X

	  a10fTemp[0] = -1.0;
	  a10fTemp[1] = 0.0;
	  a10fTemp[2] = 0.0;

	  (void) poHeader->nReplaceValue(Cstring (D_K_RotVector),
					 3, &a10fTemp[0], 4);
	  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_RotVector,
					 3, &a10fTemp[0], 4);

	  a10fTemp[0] =  -1.0f;
	  a10fTemp[1] =  0.0f;
	  a10fTemp[2] =  0.0f;

	  a10fTemp[3] = -cos(50.0 * Gs_dRADIANS_PER_DEGREE);
	  a10fTemp[4] =  0.0f;
	  a10fTemp[5] =  sin(50.0 * Gs_dRADIANS_PER_DEGREE);

	  a10fTemp[6] =  1.0f;
	  a10fTemp[7] =  0.0f;
	  a10fTemp[8] =  0.0f;

	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
					 9, a10fTemp, 4);

	}
      else if (2 == tDIP2030_header.tPart1.osc_axis)
	{
	  cout << "WARNING! oscillation axis appears to be unsupported KAPPA!\n" << flush;
	}

      a10fTemp[0] = 5.0f;
      if (   (0 == tDIP2030_header.tPart1.dtype)
	  || (21 == tDIP2030_header.tPart1.dtype) )
	{
	  a10fTemp[0] = 8.0f;
	}   
      else if (   (300 != tDIP2030_header.tPart1.dtype) 
	       || (2   != tDIP2030_header.tPart1.dtype) )
	{
	  cout << "WARNING!!  DIP2030 dtype is: " 
               << tDIP2030_header.tPart1.dtype << ".  Please contact MSC!\n"
               << endl;
	}
      (void) poHeader->nReplaceValue(sDIP2030 + "Compression", a10fTemp[0]);

      a10fTemp[0] = (float) pnHead[41];
      a10fTemp[1] = (float) pnHead[42];
      a10fTemp[2] = (float) pnHead[4] * 0.001;
      a10fTemp[3] = (float) pnHead[5] * 0.001;
      (void) poHeader->nReplaceValue(sDIP2030 + D_K_SpatialDistortionInfo,
				     4, a10fTemp, 5);
      a10fTemp[0] = 1.0;
      a10fTemp[1] = pfHead[30];
      (void) poHeader->nReplaceValue(Cstring (D_K_SourcePrefix) + 
				     D_K_Wavelength, 2, a10fTemp);

      a10fTemp[0] = 0.0;
      a10fTemp[1] = 0.0;   // Two theta
      a10fTemp[2] = 0.0;
      a10fTemp[3] = 0.0;
      a10fTemp[4] = 0.0;
      a10fTemp[5] = pfHead[31]; // xtal_to_detector
      (void) poHeader->nReplaceValue(sDIP2030 + D_K_GonioValues, 6,
				     a10fTemp, 3);

      a10fTemp[0] = pfHead[48];
      a10fTemp[1] = pfHead[49];
      a10fTemp[2] = a10fTemp[1] - a10fTemp[0];
      a10fTemp[3] = pfHead[43];
      a10fTemp[4] = (float) pnHead[51];
      a10fTemp[5] = 0.0;
      a10fTemp[6] = 0.0;
      a10fTemp[7] = 0.0;
      a10fTemp[8] = 0.0;
      a10fTemp[9] = 0.0;
      (void) poHeader->nReplaceValue(Cstring (D_K_Rotation), 10, a10fTemp, 3);
      (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_Rotation,
				     10, a10fTemp, 3);

      (void) poHeader->nReplaceValue(Cstring(D_K_Comment), 
				     Cstring(&pcHead[40]));

      a10fTemp[0] = tDIP2030_header.tPart1.g_omega;
      a10fTemp[1] = tDIP2030_header.tPart1.g_kappa;
      a10fTemp[2] = tDIP2030_header.tPart1.g_phi;
      (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
				 3, a10fTemp, 3);

/*
TODO: Set polarization based on Monochromator type found in pcHead[0].

      printf("pixel %d %d microns, pixels in frame %d x %d   which ip %d\n",
	     pnHead[4],pnHead[5],pnHead[7],pnHead[8],pnHead[9]);
      printf(" nominal center?  %d %d  xray wavelength %f  distance %6.1f mm \n",
	     pnHead[41],pnHead[42], pfHead[30],pfHead[31]);
      printf("collimator  %4.1f mm,  ip offset  %8.2f degrees,  time %6.0f sec.\n",
	     pfHead[46], pfHead[60], pfHead[43]);
      printf("initial angle %8.2f deg., final  %8.2f deg. scan rate%6.3f degr/min\n",
	     pfHead[48],pfHead[49],pfHead[50]);
      printf("scan repetitions  %d, uninterpreted values %d  %d  %d  %f\n",
	     pnHead[51],pnHead[1],pnHead[3],pnHead[52],pfHead[53]);
*/
    }

  return (nStat);
}
