//
// Copyright (c) 1999 Molecular Structure Corporation
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

#include <string.h>
#include "dtreksys.h"
#include "marip.h"
#include "Cimage_header.h"
#include "dskio.h"
#include "dtrekdefs.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

//
// WARNING // WARNING // WARNING // WARNING // WARNING // WARNING // WARNING //
// The following code has NOT been test on a 64-bit compilation
//

int 
nReadMARIPHeader(const int* pnFile, const Cstring& rsFilename,
		 char *pcBuffer, Cimage_header *poHeader)
{
  // Read a MARIP header that may or may not have been partially read already.
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
  // Ooops, the MARIP header is not always 4096 bytes long.  I found one
  // with 4000 byte headers.  Oh well.  Some headers have ASCII keywords,
  // some do not.
  //
  // The MARIP header is always 4096 bytes long.  It contains alot of binary
  // data that may be big or little endian.  So we have to figure this out,
  // maybe swap bytes.

  //MAR345_HEADER tMAR345_header;
  //MAR300_HEADER tMAR300_header;

  int   i;           // Loop counter
  int   nStat;       // Local status flag
  int   nSwap;       // Flag whether to swap (=1) header bytes or not (=0) 
  int   nLen;        // Temp variable for lengths
  bool  bSwap = FALSE;  // Flag whether to swap byteorder of image data or not
  float a10fTemp[10];
  int   a10nTemp[10];
  int   a32nTemp[32];
  Cstring a10sTemp[10];

  static Cstring sMARIP = "MARIP_";

  char  a4096cMAR[4096];
  Cstring sMARtext;

//  cout << "Inside MARIP stuff\n" << flush;
  nLen = sizeof(a4096cMAR);
  bool bUGAformat = FALSE;

  if (NULL != pcBuffer)
    {
      if (0 == strncmp("MAR", (char *)(pcBuffer+124), 3))
	{
	  // Header is probably 4000 bytes, not 4096 bytes.

	  nLen       = 4000;
	  bUGAformat = TRUE;
	}

      memcpy((void*)a4096cMAR, (void*)pcBuffer, 512);
      nLen = nLen - 512;
      dskbr((int *)pnFile, (char *)&a4096cMAR[512], &nLen, &nStat);
    }
  else
    {
      dskbr((int *)pnFile, (char *)&a4096cMAR, &nLen, &nStat);
    }

  a4096cMAR[4095] = '\0';
  sMARtext = Cstring(&a4096cMAR[128]);

  // Look at first 4 characters as a 4-byte integer

  nSwap = 0;
  memcpy((void*)a32nTemp, (void*)a4096cMAR, 32 * 4);
  nLen = a32nTemp[0];
//  cout << "nLen before is: " << nLen << endl;
  if ( (1200 > nLen) || (3450 < nLen) )
    {
      // Need to swap nLen

      nLen  = nSwapLong( (UINT4) nLen);
      nSwap = 1;
    }
  if (1 == nSwap)
    {
//      cout << "Swap first 32 longwords.\n";
      for (i = 0; i < 32; i++)
	{
	  a32nTemp[i] = nSwapLong( (UINT4) a32nTemp[i]);
	}
    }

  // Is it a MAR345 or MAR300 header?  Decide and adjust a32nTemp
  // accordingly

//  cout << "nLen after is: " << nLen << endl;

  if (1234 == nLen)
    {
      // It is a MAR345

//      cout << "MAR 345 image found.\n" << flush;
    }
  else
    {
      // It is not a 345 header, so assume a different header
      // but with keywords, use the keyword to get results 

      if (sMARtext.contains("FORMAT")) 
	{
//	  cout << "using keywords\n";

	  nStat = split(sMARtext.after("FORMAT"), a10sTemp, 3, " \n");
	  if (0 < nStat)
	    a32nTemp[1] = atoi(a10sTemp[0].string());  //Number of pixels
	  else
	    a32nTemp[1] = 2000;

	  nStat = split(sMARtext.after("HIGH"), a10sTemp, 1, " \n");
	  if (0 < nStat)
	    a32nTemp[2] = atoi(a10sTemp[0].string());  //Number of pixels
	  else
	    a32nTemp[2] = 0;

	  nStat = split(sMARtext.after("LENGTH"), a10sTemp, 1, " \n");      
	  if (0 < nStat)
	    a32nTemp[6] = atoi(a10sTemp[0].string());

	  nStat = split(sMARtext.after("HEIGHT"), a10sTemp, 1, " \n");      
	  if (0 < nStat)
	    a32nTemp[7] = atoi(a10sTemp[0].string());

	  nStat = split(sMARtext.after("TWOTHETA"), a10sTemp, 1, " \n");      
	  if (0 < nStat)
	    a32nTemp[15] = int(atof(a10sTemp[0].string()) * 1000.0f);

	  nStat = split(sMARtext.after("DISTANCE"), a10sTemp, 1, " \n");      
	  if (0 < nStat)
	    a32nTemp[8] = int(atof(a10sTemp[0].string()) * 1000.0f);

	  nStat = split(sMARtext.after("OMEGA"), a10sTemp, 3, " \n");      
	  if (0 < nStat)
	    a32nTemp[12] = int(atof(a10sTemp[1].string()));

	  nStat = split(sMARtext.after("CHI"), a10sTemp, 1, " \n");      
	  if (0 < nStat)
	    a32nTemp[14] = int(atof(a10sTemp[0].string()) * 1000.0f);
	  
	  nStat = split(sMARtext.after("PHI"), a10sTemp, 3, " \n");      
	  if (1 < nStat)
	    {
	      a32nTemp[10] = int(atof(a10sTemp[1].string()) * 1000.0f);
	      a32nTemp[11] = int(atof(a10sTemp[2].string()) * 1000.0f);
	    }
	  a32nTemp[30]  = (int) 0;
	  a32nTemp[31]  = (int) 0;
	  a32nTemp[29]  = (int) 0;
	  nStat = 0;
	}
      else
	{
	  // No keywords, but a32nTemp is loaded up with mostly proper values
	  
	  /* First 32 longs contain: */

//	  cout << "NOT using keywords\n";
	  nStat = 0;

	  MAR300_HEADER t300head;
	  memcpy((void*)&t300head, (void*)a4096cMAR, sizeof(MAR300_HEADER));

	  if (1 == nSwap)
	    {
	      // Must swap bytes in t300head.

	      int *pnTemp;

	      pnTemp = (int*) &t300head.pixels_x;

	      for (i = 0; i < 10; i++, pnTemp++)
		{
		  *pnTemp = nSwapLong( (UINT4) *pnTemp);
		}
	      float *pfTemp;
	      pfTemp = (float*)&t300head.prog_time;
	      for (i = 0; i < 15; i++, pfTemp++)
		{
		  *pfTemp = fSwapFloat(*pfTemp);
		}
	    }

	  a32nTemp[2] = t300head.high_pixels;

	  if (t300head.pixels_x != t300head.pixels_y)
	    cout << "WARNING MARIP image is NOT square!\n";
	  a32nTemp[1] = t300head.pixels_x;
	  a32nTemp[6]  = (int) 150;
	  a32nTemp[7]  = (int) 150;

	  // Are other pixel sizes possible with a MAR300?

	  a32nTemp[15] = (int) 0; // (t300head.theta        * 1000.0f);
	  a32nTemp[8]  = (int) (t300head.distance     * 1000.0f);
	  a32nTemp[9]  = (int) (t300head.lambda       * 1000000.0f);

	  if ( ((2 * t300head.pixels_x * t300head.pixels_y) + t300head.lrecl)
	       <= lFileGetSize(rsFilename))
	    {
	      // Not PCK format since not compressed

	      // (void) poHeader->nReplaceValue(Cstring (D_K_Compression),
	      //   "MAR");
	      // Set byte order here

	      bSwap = (1 == nSwap);
	    }
	  else
	    {
	      // Force PCK format later on
	      (void) poHeader->nReplaceValue(Cstring (D_K_Compression),
					     "PCK");
	    }

	  if ( (t300head.omega != 0) && (t300head.omega < 1000) )
	    {
	      a32nTemp[12]  = (int) (t300head.omega         * 1000.0f);
	    }
	  else
	    {
	      a32nTemp[12]  = (int) (t300head.omega);
	    }

	  a32nTemp[14]  = (int) 0; // (t300head.chi           * 1000.0f);
	  a32nTemp[10]  = (int) (t300head.phi_start     * 1000.0f);
	  a32nTemp[11]  = (int) (t300head.phi_end       * 1000.0f);

	  a32nTemp[30]  = (int) (t300head.centre_x / 0.150f * 100.0f + 0.5f);
	  a32nTemp[31]  = (int) (t300head.centre_y / 0.150f * 100.0f + 0.5f);
	  a32nTemp[29]  = (int) (t300head.exptime_sec    * 1000.0f);
	  // Missing from header: polarization, exposure time

/*
	  h.byteorder 	= (int  )head[ 0];
	  h.size 		= (short)head[ 1];
	  h.high 		= (int  )head[ 2];
	  h.format	= (char )head[ 3];
	  h.mode		= (char )head[ 4];
	  h.pixels	= (int  )head[ 5];
	  h.pixel_length 	= (float)head[ 6];
	  h.pixel_height 	= (float)head[ 7];
	  h.wave      	= (float)head[ 9]/1000000.;
	  h.dist    	= (float)head[ 8]/1000.;
	  h.phibeg	= (float)head[10]/1000.;
	  h.phiend	= (float)head[11]/1000.;
	  h.omebeg	= (float)head[12]/1000.;
	  h.omeend	= (float)head[13]/1000.;
	  h.chi   	= (float)head[14]/1000.;
	  h.theta 	= (float)head[15]/1000.;
*/	  
	}
    }

  a10nTemp[0] = a32nTemp[1];
  a10nTemp[1] = a32nTemp[1];
    
  poHeader->nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_MAR_IP));

  (void) poHeader->nReplaceValue(Cstring (D_K_Type), Cstring("mad"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Dim), Cstring("2"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Size1), a10nTemp[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_Size2), a10nTemp[1]);
  (void) poHeader->nReplaceValue(Cstring (D_K_DataType), 
				 Cstring("unsigned short int"));
  if (bUGAformat)
    {
      (void) poHeader->nReplaceValue(Cstring (D_K_DtdisplayOrientation), 
				     Cstring("-Y-X"));
    }
  else
    {
      (void) poHeader->nReplaceValue(Cstring (D_K_DtdisplayOrientation), 
				     Cstring("-X+Y"));
    }

  if (bSwap)
    {
      if (eCPU_big_endian == nGetByteOrderCPU())
	(void) poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
				       (Cstring)D_K_LittleEndian);
      else
	(void) poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
				       (Cstring)D_K_BigEndian);

    }
  else
    {
      (void) poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
				     sGetByteOrderCPU());
    }
  

  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNumber), (int)1);
  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNames), sMARIP);
  (void) poHeader->nReplaceValue(sMARIP + D_K_DetectorDimensions, 2, a10nTemp);

  a10fTemp[0] = (float) a10nTemp[0] * (float)a32nTemp[6] * 0.001;
  a10fTemp[1] = (float) a10nTemp[1] * (float)a32nTemp[7] * 0.001;
  (void) poHeader->nReplaceValue(sMARIP + D_K_DetectorSize, 2, a10fTemp, 5);
  (void) poHeader->nReplaceValue(sMARIP + D_K_DetectorDescription,
				 Cstring ("MARIP conversion"));

  // The following may not be correct:

  (void) poHeader->nReplaceValue(sMARIP + D_K_DetectorVectors,
				 Cstring ("1 0 0 0 1 0"));

  (void) poHeader->nReplaceValue(sMARIP + D_K_GonioNumValues,
				 (int) 6);
  (void) poHeader->nReplaceValue(sMARIP + D_K_GonioUnits,
				 Cstring ("deg deg deg mm mm mm"));
  (void) poHeader->nReplaceValue(sMARIP + D_K_GonioNames,
                 Cstring ("RotZ RotX/Swing RotY TransX TransY TransZ/Dist"));

  (void) poHeader->nReplaceValue(sMARIP + D_K_GonioVectors,
		 Cstring ("0 0 1 1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1"));

  a10fTemp[0] = 0.0;
  a10fTemp[1] = (float)a32nTemp[15] * 0.001;  // theta
  a10fTemp[2] = 0.0;
  a10fTemp[3] = 0.0;
  a10fTemp[4] = 0.0;

  nStat = split(sMARtext.after("DISTANCE"), a10sTemp, 1, " \n");
  if (1 == nStat)
    {
      a10fTemp[5] = (float)atof(a10sTemp[0].string());
    }
  else
    {
      a10fTemp[5] = (float)a32nTemp[8] *  0.001; // distance
    }
  nStat = 0;
      
  if (0.0 >= a10fTemp[5])
    a10fTemp[5] = 100.0;

  (void) poHeader->nReplaceValue(sMARIP + D_K_GonioValues, 6, a10fTemp, 3);

  (void) poHeader->nReplaceValue(sMARIP + D_K_NonunfType, 
				 Cstring(D_K_NonunfStateSimpleMask));
  (void) poHeader->nReplaceValue(sMARIP + D_K_NonunfInfo, 
				 "$(FirstScanImage)");
  (void) poHeader->nReplaceValue(sMARIP + D_K_SpatialDistortionType,
				 Cstring (D_K_SpatialTypeSimple));

  // Try to find direct beam center 

  nStat = split(sMARtext.after("CENTER"), a10sTemp, 4, " \n");
  if (4 == nStat)
    {
      a10fTemp[0] = (float)atof(a10sTemp[1].string());
      a10fTemp[1] = (float)atof(a10sTemp[3].string());
    }
  else if (0 < a32nTemp[30])
    {
      a10fTemp[0] = (float)a32nTemp[30] * 0.01f;
      a10fTemp[1] = (float)a32nTemp[31] * 0.01f;
    }
  else
    {
      a10fTemp[0] = (float) a32nTemp[1] * 0.5;
      a10fTemp[1] = (float) a32nTemp[1] * 0.5;
    }
  
  a10fTemp[2] = (float)a32nTemp[6] * 0.001;
  a10fTemp[3] = (float)a32nTemp[7] * 0.001;

  (void) poHeader->nReplaceValue(sMARIP + D_K_SpatialDistortionInfo,
				 4, a10fTemp, 5);
      
  if (bUGAformat)
    {
      (void) poHeader->nReplaceValue(sMARIP + D_K_SpatialDistortionVectors,
				     Cstring("0 1 1 0"));
    }
  else
    {
      (void) poHeader->nReplaceValue(sMARIP + D_K_SpatialDistortionVectors,
				     Cstring("1 0 0 -1"));
    }

  (void) poHeader->nReplaceValue(Cstring (D_K_SaturatedValue), 
				 65535);

  nStat = split(sMARtext.after("WAVELENGTH"), a10sTemp, 1, " \n");
  if (1 == nStat)
    {
      a10fTemp[1] = (float)atof(a10sTemp[0].string());
    }
  else
    {
      a10fTemp[1] = (float)a32nTemp[9] * 0.000001;
    }
  if (0.0 >= a10fTemp[1]) 
    a10fTemp[1] = 1.0000;
  a10fTemp[0] = 1.0;
  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePrefix) + 
				 D_K_Wavelength, 2, a10fTemp);

  nStat = split(sMARtext.after("FORMAT"), a10sTemp, 3, " \n");
  if (3 == nStat)
    {
      (void) poHeader->nReplaceValue(sMARIP + "FORMAT", a10sTemp[1]);
      if ("PCK" == a10sTemp[1])
	{
	  (void) poHeader->nReplaceValue(Cstring (D_K_Compression),
					 a10sTemp[1]);
	}
    }
  nStat = 0;

  nStat = split(sMARtext.after("POLAR"), a10sTemp, 1, " \n");

  if (1 == nStat)
    a10fTemp[0] = (float)atof(a10sTemp[0].string());
  if (0.0 >= a10fTemp[0])
    a10fTemp[0] = 0.5;

  nStat = 0;

  a10fTemp[1] = 0.0;
  a10fTemp[2] = 1.0;
  a10fTemp[3] = 0.0;
  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
				 4, a10fTemp, 5);

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
				 Cstring ("Omega LgArc SmArc"));
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioUnits,
				 Cstring ("deg deg deg"));
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
				 Cstring ("1 0 0 0 0 1 0 1 0"));
  a10fTemp[0] = (float)a32nTemp[12] * 0.001; // omega
  a10fTemp[1] = (float)a32nTemp[14] * 0.001; // chi
//  a10fTemp[2] = (float)a32nTemp[10] * 0.001; // phi
  a10fTemp[2] = 0.0f;
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
				 3, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + 
				 D_K_GonioDescription,
				 Cstring ("Single w/2-arc goniometer head"));

  // Determine scan template and seq info from the filename
  // We need the sequence number in order to check special circumstances
  // in the rotation start, end, range

  // Change last 3 digits in rsFilename to ?s

  Cstring sTemp = rsFilename;
  int a3nSeqInfo[3] = {0, 1, 0};

  sTemp = sBuildScanTemplate(rsFilename, 4, &a3nSeqInfo[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanTemplate), sTemp);

  (void) poHeader->nReplaceValue(Cstring (D_K_ScanSeqInfo), 3, a3nSeqInfo);

  // TODO: Watch out for cross from 359 to 0 to 1 degree.

  a10fTemp[0] = (float)a32nTemp[10] * 0.001;
  a10fTemp[1] = (float)a32nTemp[11] * 0.001;

  // BEWARE: rotation angle increment must be positive!

  a10fTemp[2] = a10fTemp[1] - a10fTemp[0];
  if (180.0 < ABS(a10fTemp[2]))
    {
      // Image probably crosses the boundary at 360.  Is it at the beginning?
      // Or at the end?

      cout << "WARNING, inconsistence rotation start and end adjusted, likely cause:\n"
           << "    crosses from less than 360 to more than 0 degrees.\n" << endl;
      if (5 > a3nSeqInfo[0])
	{
	  // Probably near the beginning of the scan, so adjust a10fTemp[0]

	  a10fTemp[0] = a10fTemp[0] - 360.0;
	  a10fTemp[2] = a10fTemp[1] - a10fTemp[0];
	}
      else
	{
	  // Probably near the end of the scan, so adjust a10fTemp[1]

	  a10fTemp[1] = a10fTemp[1] + 360.0;
	  a10fTemp[2] = a10fTemp[1] - a10fTemp[0];
	}
    }


  nStat = split(sMARtext.after("TIME"), a10sTemp, 1, " \n");
  if (1 == nStat)
    a10fTemp[3] = (float)atof(a10sTemp[0].string());
  else if (0 < a32nTemp[29])
    a10fTemp[3] = (float)a32nTemp[29] * 0.001f;
  else
    a10fTemp[3] = 20.0;    
  nStat = 0;

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
				 Cstring ("Omega"));
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + 
				 D_K_RotAxisName,
				 Cstring ("Omega"));

  a10fTemp[0] = 1.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;

  (void) poHeader->nReplaceValue(Cstring (D_K_RotVector),
				 3, &a10fTemp[0], 4);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_RotVector,
				 3, &a10fTemp[0], 4);

  if ( (0 != nStat) && ("" != sGetEnv("DTREK_RAXIS_IGNORE")))
    {
      cout << "Warning in MARIP, error reading header, but ignored!\n";
      nStat = 0;
    }

  // Read the overflow table (if it exists)

  int   nNumOverflows;

  nNumOverflows = a32nTemp[2];

  poHeader->nReplaceValue("OVERFLOW_COUNT", nNumOverflows);
  int   *pnOffset;
  int   *pnValues;
  int   nToRead;
//  cout << "Num overflows is " << nNumOverflows << endl;
  if (0 < nNumOverflows)
    {
      pnOffset = new int   [nNumOverflows];
      pnValues = new int   [nNumOverflows];
      for (i = 0; (0 == nStat) && (i < nNumOverflows); i++)
	{
	  nToRead = 4; // sizeof(int);
	  (void) dskbr((int *)pnFile, (char *)&pnOffset[i], &nToRead, &nStat);
	  if (1 == nSwap) 
	    pnOffset[i] = nSwapLong( (UINT4) pnOffset[i]);
	  pnOffset[i] -= 1;
	  nToRead = 4; // sizeof(int);
	  (void) dskbr((int *)pnFile, (char *)&pnValues[i], &nToRead, &nStat);
	  if (1 == nSwap) 
	    pnValues[i] = nSwapLong( (UINT4) pnValues[i]);
//	  cout << "offset, value: " << pnOffset[i] << ", " << pnValues[i] << endl;
	}

      poHeader->nReplaceValue("OVERFLOW_OFFSET", nNumOverflows, pnOffset);
      poHeader->nReplaceValue("OVERFLOW_VALUES", nNumOverflows, pnValues);

      // Read any dummy overflows in the overflow table, so the file pointer
      // is at the beginning of pixel data (overflows are written in groups
      // of 8

      nNumOverflows = 8 - (nNumOverflows % 8);
      for (i = 0; (0 == nStat) && (i < nNumOverflows); i++)
	{
	  nToRead = 4;
	  (void) dskbr((int *)pnFile, (char *)&pnOffset[0], &nToRead, &nStat);
	  nToRead = 4;
	  (void) dskbr((int *)pnFile, (char *)&pnValues[0], &nToRead, &nStat);
	}

      delete [] pnOffset;
      delete [] pnValues;
    }

  return (nStat);
}
