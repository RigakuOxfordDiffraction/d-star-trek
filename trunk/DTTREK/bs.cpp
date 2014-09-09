//
// Copyright (c) 1997 Molecular Structure Corporation
//
#include <string.h>
#include "dtreksys.h"
#include "raxis.h"
#include "bs.h"
#include "Cimage.h"
#include "Cimage_header.h"
#include "dskio.h"
#include "dtrekdefs.h"
#include "dtrekvec.h"
#include "Cgoniometer.h"

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

int 
nReadSiemensHeader(const int* pnFile, const Cstring& rsFilename,
		 char *pcBuffer, Cimage_header *poHeader)
{
  // Read a Siemens header that may or may not have been partially read already.
  // Convert the header to a d*TREK header.

  // =======                    ===
  // WARNING: This routine does NOT read ALL Bruker image file formats!
  // =======                    ===

  //
  // Use the dskio routines for i/o.
  // The file to read must have been opened with the dskbor routine already.
  // Arguments:
  //   pnFile    pointer to open file tag for dskio routines.
  //   pcBuffer  pointer to first 512-bytes of header if header partially read
  //                  or NULL if header not partially read.
  //   poHeader  pointer to image header to create.
  //   rsFilename image filename; the filename contains some sequence info
  //
  // Returns
  //     0 success
  // not 0 failure
  //
  // The Siemens header is sometimes multiple of 512 bytes long.

  int i,j,k;  // Loop counter
  int nStat;  // Local status flag
  int nLen;   // Pixels per line
  int nSwap;  // Flag whether to swap (=1) bytes in integer*4's or not (=0)
  int nCount; 
  Cimage_header oTempHeader;

  char       a512cTemp[513];
  Cstring    sSHeader;
  Cstring    sTemp;
  Cstring    sBSPrefix = "BS_";
  bool       bUseLinearInterp = false;
  double     fDistanceOffset = 0.0;       // The distance in the image "DISTANC" is sometimes missleading.
  double     a2fImageDim[2];
  int        a2nImageDim[2];

//  cout << "nReadSiemensHeader called!\n";
  if (NULL != pcBuffer)
    {
      memcpy((void*)&a512cTemp, (void*)pcBuffer, 512);
      nStat = 0;
    }
  else
    {
      nLen = 512;
      dskbr((int *)pnFile, (char *)&a512cTemp[0], &nLen, &nStat);
    }

  // Have 512 of header, now compute size of rest of header and read it in

  if (0 != nStat)
    {
//      cout << "nstat not 0\n";
      return (nStat);
    }
  else
    {
      a512cTemp[512] = '\0';
      sSHeader = a512cTemp;
      nCount   = sSHeader.find("HDRBLKS:");
      sTemp    = sSHeader.substr(nCount + 8, 72);
      nCount   = 0;
      nCount   = atoi(sTemp.string());
      if (0 < nCount)
	{
	  char *pcTemp;
	  
	  nLen = (nCount - 1) * 512;
	  pcTemp = new char [1 + nLen];
	  pcTemp[nLen] = '\0';

	  // Read the rest of the header in 512 byte increments;

//	  cout << "nLen is: " << nLen << endl;
	  dskbr((int *)pnFile, pcTemp, &nLen, &nStat);
	  if (0 == nStat)
	    {
	      sSHeader += pcTemp;
	    }
	  else
	    {
	      cout << "ERROR reading header!\n";
	    }
	  delete pcTemp;
	}
    }
  if (0 == nStat)
    {
      // Put all the Siemens keywords into the header with BS_ in front
      // of them

      Cstring sTemp2, sPrev;

      sPrev = "";
      for (i = 0; i < sSHeader.length(); i=i+80)
	{
	  sTemp = sSHeader.substr(i, 7);
	  if ('.' < sTemp.GetAt(0))
	    {
	      // No binary so legit keyword.  Remove trailing spaces

          for (k=sTemp.length()-1;(k>0) && (sTemp.GetAt(k)==' ');k--);
          sTemp = sTemp.substr(0,k+1);

	      // Get value, strip off leading and trailing spaces, non-prntables

	      sTemp2 = sSHeader.substr(i+8,72);
          for (j=0;((j < sTemp2.length()) && ((char) sTemp2.GetAt(j) == ' '));j++);
          for (k=sTemp2.length()-1;(k>=0) && ((char) sTemp2.GetAt(k) ==' ');k--);
          if ((k>=0) && (j<=k)) 
              sTemp2 = sTemp2.substr(j,k-j+1);
          else
              sTemp2 = "";

	      // Look for previous info of same keyword in header
          
          nStat = oTempHeader.nGetValue(sBSPrefix + sTemp, &sPrev);
          if (0 == nStat)
          {
              // Previous value detected, so sTemp2 is a continuation
              if (sTemp2.length())		  
                  sTemp2 = sPrev + ' ' + sTemp2;
              else 
                  sTemp2 = sPrev;
          }
          //	      cout << "keyword,value: " << sBSPrefix + sTemp << '\n' << sTemp2 << endl;
	      oTempHeader.nReplaceValue(sBSPrefix + sTemp, sTemp2);
	    }
	}
      nStat = 0;
    }

  // Now convert Siemens items to d*TREK items

  oTempHeader.nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_BRUKER_SIEMENS));

  float afTemp[10];
  int   anTemp[10];
  (void) oTempHeader.nReplaceValue(Cstring (D_K_DetectorNumber), (int)1);
  (void) oTempHeader.nReplaceValue(Cstring (D_K_DetectorNames), 
				   sBSPrefix);

  (void) oTempHeader.nReplaceValue(Cstring (D_K_Type), Cstring("mad"));
  (void) oTempHeader.nReplaceValue(Cstring (D_K_Dim), Cstring("2"));

  nStat = oTempHeader.nGetValue(sBSPrefix + "NCOLS", &a2nImageDim[0]);
  if (0 != nStat) 
    {
      a2nImageDim[0] = 512;
      cout << "nstat cols non-zero" << endl;
    }
  nStat = oTempHeader.nGetValue(sBSPrefix + "NROWS", &a2nImageDim[1]);
  if (0 != nStat) 
    {
      a2nImageDim[1] = 512;
      cout << "nstat rows non-zero" << endl;
    }

  nStat = 0;

//  cout << "cols, rows: " << a2nImageDim[0] << ' ' << a2nImageDim[1] << endl;
  (void) oTempHeader.nReplaceValue(sBSPrefix + 
				   D_K_DetectorDimensions, 2, 
				   a2nImageDim);
  (void) oTempHeader.nReplaceValue(Cstring (D_K_Size1), a2nImageDim[0]);
  (void) oTempHeader.nReplaceValue(Cstring (D_K_Size2), a2nImageDim[1]);
  (void) oTempHeader.nReplaceValue(Cstring (D_K_DataType), 
                                 Cstring("unsigned short int"));
  (void) oTempHeader.nReplaceValue(Cstring (D_K_Compression),
                                 Cstring("BS1"));
  nSwap = 0;
  (void) oTempHeader.nGetValue(sBSPrefix + "WORDORD", &nSwap);
  
/*
  if (0 == nSwap)
    {
      // Image is little_indian
      
      oTempHeader.nReplaceValue(Cstring (D_K_ByteOrder),
				Cstring (D_K_LittleEndian));
    }
  else
    {
      oTempHeader.nReplaceValue(Cstring (D_K_ByteOrder),
				Cstring (D_K_BigEndian));
    }
    // Have byte order of image match byte order of CPU.  The read routines
    // will swap bytes so image does have this byte order in the end. 
    // Normally this would be done by the Cimage:: class, but we want to
    // do it here within the bs.cc routines.
*/
  oTempHeader.nReplaceValue(Cstring (D_K_ByteOrder),
			    sGetByteOrderCPU());

  // We don't know the pixel size since that is in separate spatial distortion
  // tables

  (void) oTempHeader.nGetValue(sBSPrefix + "DETTYPE", &sTemp);
  if (sTemp.contains("CCD-PXL-ARR"))
    {
      fDistanceOffset = 0.0;
      a2fImageDim[0] = 160.0; 
      a2fImageDim[1] = 160.0;
      (void) oTempHeader.nReplaceValue(D_K_DtdisplayOrientation, 
				       Cstring("-Y-X"));
      (void) oTempHeader.nReplaceValue(Cstring (D_K_SourcePolarz),
				       Cstring ("0.99 0 1 0"));

      // Ion chamber readings are sometimes available in the header
      // as the second value to the NCOUNTS keyword

      anTemp[1] = 0;
      oTempHeader.nGetValue(sBSPrefix + "NCOUNTS", 2, &anTemp[0]);
      if (0 < anTemp[1])
	oTempHeader.nReplaceValue( Cstring(D_K_SourceIntensity), anTemp[1]);

      // Bruker-Siemens detector has saturated value of 65535 or 131071
      // depending on the NEXP keyword, if 1, then correlated sum NOT performed.

      anTemp[0] = 1;
      anTemp[1] = 65535;
      nStat = oTempHeader.nGetValue( sBSPrefix + "NEXP", &anTemp[0]);
      if (1 != anTemp[0]) anTemp[1] = 131071;
      (void) oTempHeader.nReplaceValue(Cstring (D_K_SaturatedValue),
				       anTemp[1]);
      
    }
  else if (sTemp.contains("CCD-PXL-2K"))
    {
      fDistanceOffset = 0.0;
      a2fImageDim[0] = 80.4;   // 512 * 0.157
      a2fImageDim[1] = 80.4;
      (void) oTempHeader.nReplaceValue(D_K_DtdisplayOrientation,
				       Cstring ("-X+Y"));
      (void) oTempHeader.nReplaceValue(Cstring (D_K_SourcePolarz),
				       Cstring ("0.5 1 0 0"));
    }
  else if (sTemp.contains("CCD-PXL-L6000"))
    {
      // It is a Smart6000

      (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_NonunfType, 
				       "$(FirstScanImage)");

      fDistanceOffset = 0.0;
      a2fImageDim[0] = 160.0;
      a2fImageDim[1] = 160.0;
      (void) oTempHeader.nReplaceValue(D_K_DtdisplayOrientation, 
				       Cstring("-X+Y"));
      (void) oTempHeader.nReplaceValue(Cstring (D_K_SourcePolarz),
				       Cstring ("0.50 1 0 0"));

      // Ion chamber readings are sometimes available in the header
      // as the second value to the NCOUNTS keyword

      anTemp[1] = 0;
      oTempHeader.nGetValue(sBSPrefix + "NCOUNTS", 2, &anTemp[0]);
      if (0 < anTemp[1])
	oTempHeader.nReplaceValue( Cstring(D_K_SourceIntensity), anTemp[1]);

      anTemp[1] = 262136;
      (void) oTempHeader.nReplaceValue(Cstring (D_K_SaturatedValue),
				       anTemp[1]);
    }
  else if (sTemp.contains("Multiwire")) {
      fDistanceOffset = 8.0;
      a2fImageDim[0] = 92.16;   // 1024 * 0.090
      a2fImageDim[1] = 92.16;
      (void) oTempHeader.nReplaceValue(D_K_DtdisplayOrientation,
				       Cstring ("+X+Y"));
      (void) oTempHeader.nReplaceValue(Cstring (D_K_SourcePolarz),
				       Cstring ("0.5 1 0 0"));
      bUseLinearInterp = true;

  } else {
      fDistanceOffset = 0.0;

      printf("WARNING:  d*TREK does not know the exact type of detector used!\n");

      a2fImageDim[0] = 80.4;   // 512 * 0.157
      a2fImageDim[1] = 80.4;
      (void) oTempHeader.nReplaceValue(D_K_DtdisplayOrientation,
				       Cstring ("-X+Y"));
      (void) oTempHeader.nReplaceValue(Cstring (D_K_SourcePolarz),
				       Cstring ("0.5 1 0 0"));
    }

  (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_DetectorSize,
				 2, a2fImageDim, 5);
  (void) oTempHeader.nReplaceValue(sBSPrefix + 
					      D_K_DetectorDescription,
				     Cstring ("BS conversion"));
// The following may not be correct:
  (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_DetectorVectors,
				 Cstring ("1 0 0 0 1 0"));
  (void) oTempHeader.nReplaceValue(sBSPrefix + 
				 D_K_DetectorRefineFlags,
				 Cstring ("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"));
  (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_GonioNumValues,
				 (int) 6);
  (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_GonioUnits,
				 Cstring ("deg deg deg mm mm mm"));
  (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_GonioNames,
				 Cstring ("RotZ RotX/Swing RotY TransX TransY TransZ/Dist"));
  
  (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_GonioVectors,
				 Cstring ("0 0 1 -1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1"));

  (void) oTempHeader.nGetValue(sBSPrefix + "ANGLES", 4, afTemp);
  if (afTemp[0] > 180.0f)
    afTemp[0] = afTemp[0] - 360.0;
  afTemp[1] = afTemp[0];  // Detector 2theta
  afTemp[0] = 0.0;
  afTemp[2] = 0.0;
  afTemp[3] = 0.0;
  afTemp[4] = 0.0;
  (void) oTempHeader.nGetValue(sBSPrefix + "DISTANC", &afTemp[5]);
  afTemp[5] *= 10.0;       // Detector distance
  afTemp[5] += fDistanceOffset;
  (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_GonioValues, 
				 6, afTemp, 3);

  (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_NonunfType, 
				 Cstring(D_K_NonunfStateNone));

      
  // Beam position is stored as (x, SizeY-y)

  (void) oTempHeader.nGetValue(sBSPrefix + "CENTER", 2, afTemp);
  (void) oTempHeader.nGetValue(Cstring (D_K_Size2), &afTemp[4]);
  afTemp[1] = afTemp[4] - afTemp[1];

  (void) oTempHeader.nReplaceValue(sBSPrefix + 
				   D_K_SpatialBeamPosition,
				   2, afTemp, 3);

  sTemp = "";
  (void)oTempHeader.nGetValue(sBSPrefix + "WARPFIL", &sTemp);
  if (("LINEAR" == sTemp) || (bUseLinearInterp))
  {
      // It is an APEX detector
      
      (void) oTempHeader.nReplaceValue(sBSPrefix + 
          D_K_SpatialDistortionType,
          D_K_SpatialTypeSimple);
      afTemp[2] = a2fImageDim[0]/a2nImageDim[0];
      afTemp[3] = a2fImageDim[1]/a2nImageDim[1];
      
      (void) oTempHeader.nReplaceValue(sBSPrefix + 
          D_K_SpatialDistortionInfo, 4,
          afTemp, 4);
  }
  else
  {
      // Need a spatial distortion file.
      // 1.  Look for a .p4p or .P4P file with same basename in the image
      //     directory, then in the current working directory
      
      // 2. Else use the d*TREK distor fileset convention
      
      (void) oTempHeader.nReplaceValue(sBSPrefix + 
          D_K_SpatialDistortionType,
          D_K_SpatialTypeInterp);
      sTemp = rsFilename;
      if (!sTemp.contains(".sfrm"))
      {
          sTemp = sTemp.before((int)sTemp.length()-3) + "p4p";
          if (!bFileExists(sTemp))
          {
              sTemp = sTemp.before((int)sTemp.length()-3) + "P4P";
              if (!bFileExists(sTemp))      
              {
                  sTemp = sGetCWD() + sFileGetBasename(rsFilename);
                  sTemp = sTemp.before((int)sTemp.length()-3) + "p4p";
                  if (!bFileExists(sTemp))
                      sTemp = sTemp.before((int)sTemp.length()-3) + "P4P";
                  if (!bFileExists(sTemp))
                      sTemp = "$(distor)";    // The default
              }
          }
      }
      else if (sTemp.contains(".sfrm"))
      {
          sTemp = sTemp.before((int)sTemp.length()-4) + "spin";
          if (!bFileExists(sTemp))
          {
              sTemp = sGetCWD() + sFileGetBasename(rsFilename);
              sTemp = sTemp.before((int)sTemp.length()-4) + "spin";
              if (!bFileExists(sTemp))
                  sTemp = "$(distor)";
          }
      }
      else
          sTemp = "$(distor)";	
      
      (void) oTempHeader.nReplaceValue(sBSPrefix + 
          D_K_SpatialDistortionInfo,
          sTemp);
  }

  (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_SpatialDistortionVectors,
				 Cstring("0 1 1 0"));
  
      
  afTemp[0] = 1.0;
  afTemp[1] = 0.0;
  (void) oTempHeader.nGetValue(sBSPrefix + "WAVELEN", &afTemp[1]);
  if (0.0 >= afTemp[1]) 
      afTemp[1] = 1.54178f;
  (void) oTempHeader.nReplaceValue(Cstring (D_K_SourcePrefix) + 
				 D_K_Wavelength, 2, afTemp);

  (void) oTempHeader.nReplaceValue(Cstring (D_K_SourceSpectralDispersion),
				   Cstring ("0.0002 0.0002"));
  (void) oTempHeader.nReplaceValue(Cstring (D_K_SourceCrossfire),
				   Cstring ("0.0002 0.0002 0.0 0.0"));
  (void) oTempHeader.nReplaceValue(Cstring (D_K_SourceSize),
				   Cstring ("0.0 0.0 0.0 0.0"));
  (void) oTempHeader.nReplaceValue(Cstring (D_K_SourceVectors),
				   Cstring ("0 0 1 0 1 0 1 0 0"));
  (void) oTempHeader.nReplaceValue(Cstring (D_K_SourceValues),
				   Cstring ("0 0"));
  (void) oTempHeader.nReplaceValue(Cstring (D_K_SourceRefineFlags),
				 Cstring ("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"));

  if (0 == oTempHeader.nGetValue(sBSPrefix + "CELL", 6, afTemp))
    {
      (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + (D_K_CrystalXUnitCell),
				     6, afTemp, 3);
    }

  (void) oTempHeader.nGetValue(sBSPrefix + "START", &afTemp[0]);
  (void) oTempHeader.nGetValue(sBSPrefix + "RANGE", &afTemp[2]);
  (void) oTempHeader.nGetValue(sBSPrefix + "CUMULAT", &afTemp[3]);
  (void) oTempHeader.nGetValue(sBSPrefix + "NSTEPS", &afTemp[4]);
  (void) oTempHeader.nGetValue(sBSPrefix + "INCREME", &afTemp[5]);

  // If the crystal rotation was negative, then try some tricks
  // to make it look like it is running forwards.

  // We may need to invert direction of rotation axis, too.

  bool bNegRotate = FALSE;
  if (0.0 > afTemp[5])
    {
      bNegRotate = TRUE;
      // Negate the start, and call it the end
//+18-June-2000
      // A little test
      //new:      afTemp[1] = afTemp[0];
      //old:
      afTemp[1] = -afTemp[0];
//-18-June-2000

      // Compute the start from the range

      afTemp[0] = afTemp[1] - afTemp[2];
    }
  else
    {
      // Compute the end from the range

      afTemp[1] = afTemp[0] + afTemp[2];
    }

  if (0.0 > afTemp[2])
    {
      cout << "BS_WARNING: RANGE is negative!" << endl;
    }

  afTemp[5] = 0.0;
  afTemp[6] = 0.0;
  afTemp[7] = 0.0;
  afTemp[8] = 0.0;
  afTemp[9] = 0.0;

  (void) oTempHeader.nReplaceValue(Cstring (D_K_Rotation),
				 10, afTemp, 3);

  (void) oTempHeader.nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_Rotation,
				   10, afTemp, 3);

  // Home lab rotation vector is usually 1,0,0; but at IMCA it is -1,0,0
  // if the rotation axis is omega

  (void) oTempHeader.nGetValue(sBSPrefix + "PROGRAM", &sTemp);

  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + 
				   D_K_GonioNumValues,
				   (int) 3);
  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + 
                   D_K_GonioUnits,
				   Cstring ("deg deg deg"));

  float fOmega, fKappaChi, fPhi, fKappaAxis;

  if (sTemp.contains("KAPPA"))
    {
      // Items here are specific for the Bruker 2x2 detector at IMCA-CAT

      // Get the Kappa angles

      (void) oTempHeader.nGetValue(sBSPrefix + "AXES2", 4, afTemp);
      
      // Siemens AXES2 are Omega, Phi, Kappa, Kappa_axis_wrt_omega

      fOmega     = afTemp[0];
      fPhi       = afTemp[1];
      fKappaChi  = afTemp[2];
      fKappaAxis = afTemp[3];

      (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + 
				       D_K_GonioDescription,
				       Cstring ("Kappa 3-circle"));

      (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioNames,
				       Cstring ("Omega Kappa Phi"));

      afTemp[0] = -1.0f;
      afTemp[1] =  0.0f;
      afTemp[2] =  0.0f;

      afTemp[3] = -cos((double)fKappaAxis * Gs_dRADIANS_PER_DEGREE);
      afTemp[4] =  0.0f;
      afTemp[5] =  sin((double)fKappaAxis * Gs_dRADIANS_PER_DEGREE);

      afTemp[6] =  1.0f;
      afTemp[7] =  0.0f;
      afTemp[8] =  0.0f;

      (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
				       9, afTemp, 4);

      (void) oTempHeader.nReplaceValue(Cstring (D_K_RotVector),
				       Cstring ("-1.0 0.0 0.0"));
      (void) oTempHeader.nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_RotVector,
				       Cstring ("-1.0 0.0 0.0"));

      (void) oTempHeader.nReplaceValue(Cstring (D_K_SourcePolarz),
				       Cstring ("0.98 1 0 0"));
      (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_NonunfType, 
				       Cstring(D_K_NonunfStateSimpleMask));
      (void) oTempHeader.nReplaceValue(sBSPrefix + D_K_NonunfInfo, 
				       "$(bruker.mask)");


    }
  else
    {
      // Items here are specific for Bruker detectors in the home lab

      // Siemens ANGLES are 2theta, omega, phi or chi
      // Rotation axis is always 0.0

      (void) oTempHeader.nGetValue(sBSPrefix + "ANGLES", 4, afTemp);

      fOmega    = afTemp[1];
      fPhi      = afTemp[2];
      fKappaChi = afTemp[3];

      (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
				       3, &afTemp[1], 3);
      (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + 
				       D_K_GonioDescription,
				       Cstring ("Eulerian 3-circle"));
      (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioNames,
				       Cstring ("Omega Chi Phi"));

      sTemp = "";
      oTempHeader.nGetValue(sBSPrefix + "MODEL", &sTemp);
      if ("P4" == sTemp)
	{
	  // It is a P4 goniometer
	  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
					   Cstring ("-1 0 0 0 0 1 1 0 0"));
	  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + 
					   D_K_GonioDescription,
					   Cstring ("Bruker/Siemens P4 3-circle"));
	}
      else if (sTemp.contains("PLATFORM"))
	{
	  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
					   Cstring ("-1 0 0 0 0 1 1 0 0"));
	  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + 
					   D_K_GonioDescription,
					   Cstring ("Bruker/Siemens Platform 3-circle, fixed Chi"));

	}
      else if (sTemp.contains("KAPPA"))
	{
	  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + 
					   D_K_GonioDescription,
					   Cstring ("Bruker Kappa"));

	  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioNames,
				       Cstring ("Omega Kappa Phi"));
	  // Try to determine the Kappa axis alpha angle
	  sTemp = sTemp.after("KAPPA");
	  sTemp = sTemp.after("[");
	  sTemp = sTemp.before("]");
	  afTemp[0] = atof(sTemp.string());
	  if (48.0 > afTemp[0] || 52.0 < afTemp[0]) afTemp[0] = 50.0;
	  afTemp[0] = afTemp[0] * Gs_dRADIANS_PER_DEGREE;
	  sTemp = Cstring("-1 0 0 ") + Cstring(cos(afTemp[0])) + " 0 -"
	    + Cstring(sin(afTemp[0])) + Cstring(" -1 0 0");
	  //	  cout << "sTemp is: " << sTemp << endl;
	  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
					   sTemp);
//               Cstring ("-1 0 0 0.6427876 0 -0.7660444 -1 0 0"));
//                                sin(40)    -cos(40)
//                                cos(50)    -sin(50)
	}
      else
	{
	  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
					   Cstring ("-1 0 0 0 0 1 1 0 0"));
	  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + 
					   D_K_GonioDescription,
					   Cstring ("Bruker/Siemens Unknown"));

	}
      (void) oTempHeader.nReplaceValue(Cstring (D_K_SourcePolarz),
				       Cstring ("0.5 1 0 0"));
    }

  (void) oTempHeader.nGetValue(sBSPrefix + "AXIS", (int *)&anTemp[0]);

  if (2 == anTemp[0])
    {
      (void) oTempHeader.nReplaceValue(Cstring (D_K_RotAxisName),
				     Cstring ("Omega"));
      (void) oTempHeader.nReplaceValue(Cstring (D_K_ScanPrefix) + 
				     D_K_RotAxisName,
				     Cstring ("Omega"));
      // Omega datum is 0.
      afTemp[0] = 0.0f;
      afTemp[1] = fKappaChi;
      afTemp[2] = fPhi;
      anTemp[0] = 0;
    }
  else if (3 == anTemp[0])
    {
      (void) oTempHeader.nReplaceValue(Cstring (D_K_RotAxisName),
				     Cstring ("Phi"));
      (void) oTempHeader.nReplaceValue(Cstring (D_K_ScanPrefix) + 
				     D_K_RotAxisName,
				     Cstring ("Phi"));

      // Phi datum is 0.

      afTemp[0] = fOmega;
      afTemp[1] = fKappaChi;
      afTemp[2] = 0.0f;
      anTemp[0] = 2;
    }
  else
    {
      // WARNING!

      anTemp[0] = 0;
      afTemp[0] = 0.0f;
      afTemp[1] = fKappaChi;
      afTemp[2] = fPhi;

      cout << "ERROR NON-Omega and NON-Phi rotation axis!\n";
    }

  
  (void) oTempHeader.nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
				   3, afTemp, 3);

  Cgoniometer *poGoniometer;
  poGoniometer = new Cgoniometer(oTempHeader, D_K_CrystalPrefix);
  if (poGoniometer->bIsAvailable())
    {
      // anTemp[0] designates whether omega or phi axis

      poGoniometer->nGetRotVector(anTemp[0], (float *)&afTemp[5]);
      if (bNegRotate)
	{
	  vMulVec3DScalar((float *)&afTemp[5], -1.0, (float *)&afTemp[5]);
	}
    }
  else
    {
      afTemp[5] = 1.0f;
      afTemp[6] = 0.0f;
      afTemp[7] = 0.0f;
      cout << "WARNING crystal goniometer vector problem!\n";
    }

  (void) oTempHeader.nReplaceValue(Cstring (D_K_RotVector),
				   3, &afTemp[5], 4);
  (void) oTempHeader.nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_RotVector,
				   3, &afTemp[5], 4);

  delete poGoniometer;

  // Determine scan template and seq info from the filename

  // Change last 4 digits in rsFilename to ?s

  int a3nSeqInfo[3] = {0, 1, 0};

  sTemp = sBuildScanTemplate(rsFilename, 4, &a3nSeqInfo[0]);
  (void) oTempHeader.nReplaceValue(Cstring (D_K_ScanTemplate),
				 sTemp);

  (void) oTempHeader.nReplaceValue(Cstring (D_K_ScanSeqInfo),
				 3, a3nSeqInfo);

  oTempHeader.vSetFilename(rsFilename);
  *poHeader = oTempHeader;

  if ( (0 != nStat) && ("" != sGetEnv("DTREK_SIEMENS_IGNORE")))
    {
      cout << "Warning in bs, error reading header, but ignored!\n";
      nStat = 0;
    }

  return (nStat);
}

int 
nReadSiemensData(const int *pnFile, unsigned short int *puiData,
	       Cimage_header *poHeader, int *pnBSType)
{
  // Read Siemens binary data with the dskio routines.  The Siemens file is
  // already open, the header has been read, and the read position is at
  // the start of binary data.
  // This routine is normally called from the Cimage class.
  //
  // *pnFile   pointer to integer file tag for dskio routines
  // *puiData  pointer to where to put the read data
  // *pvHeader pointer to current image header
  //

  int i=0;            // Loop counter
  int nStat;           // Local status
  int a10nTemp[10];
  int nDim0, nDim1;    // Dimensions of the image.
  int nBytesPerPixel;
  int nToRead;
  unsigned short int *puiTemp;
  unsigned char *pucTemp;

//  cout << "nReadSiemensData called" << endl;
  nStat = poHeader->nGetValue(D_K_Size1, &nDim0);
  nStat = nStat + poHeader->nGetValue(D_K_Size2, &nDim1);
  nStat = nStat + poHeader->nGetValue("BS_NPIXELB", &nBytesPerPixel);

//  cout << "size is " << nDim0 << ' ' << nDim1 << endl;
//  cout << "bytes per pixel is " << nBytesPerPixel << endl;

  if (0 == nStat) 
    {
      // Number of bytes to read

      nToRead = nDim0 * nDim1 * nBytesPerPixel;

      // Read into the second half of reserved memory location

      if (1 == nBytesPerPixel)
	{
	  // Read into second half of image for "uncompression later"

	  puiTemp = puiData + (nToRead / 2);
	  (void) dskbr((int *)pnFile, (char *)puiTemp, &nToRead, &nStat);
	}
      else if (2 == nBytesPerPixel)
	{
	  // Read directly into image buffer

//	  cout << "ntoread is " << nToRead << endl;
	  (void) dskbr((int *)pnFile, (char *)puiData, &nToRead, &nStat);
//	  cout << "nstat, ntoread is: " << nStat << ' ' << nToRead << endl;
	}
      else
	{
	  cout << "ERROR, unable to read 4-bytes per pixel image!\n";
	  nStat = -2;
	  return (nStat);
	}
    }

  if (0 == nStat)
    {
      if (2 == nBytesPerPixel)
	{
	  if (eCPU_big_endian == nGetByteOrderCPU())
	    {
	      // Swap bytes since we know the image is always little_endian

	      unsigned short int *puiTemp2;
	      puiTemp2 = puiData;
	      for (i = 0; i < nToRead/2; i++)
		{
		  *puiTemp2 = ((*puiTemp2 & 0xff) << 8) | ((*puiTemp2) >> 8);
		  puiTemp2++;
		}
	    }
	}
      else if (1 == nBytesPerPixel)
	{
	  // Unpack from byte to unsigned short int, swap bytes if necessary

	  unsigned char *pucTemp2;
	  pucTemp  = (unsigned char *) puiTemp;
	  pucTemp2 = (unsigned char *) puiData;
	  if (eCPU_little_endian == nGetByteOrderCPU())
	    {
	      for (i = 0; i < nToRead; i++)
		{
		  *pucTemp2++ = *pucTemp++;
		  *pucTemp2++ = 0;
		}
	    }
	  else
	    {
	      for (i = 0; i < nToRead; i++)
		{
		  *pucTemp2++ = 0;
		  *pucTemp2++ = *pucTemp++;
		}
	    }
	}
    }

  // At this point we have the binary pixel values, but not
  // any underflow and overflow info

  union {
    short int           *pi;
    unsigned short int  *pui;
    int                 *pn;
    unsigned int        *pun;
    INT4                *pl;
    UINT4               *pul;
    char                *pc;
    unsigned char       *puc;
    float               *pf;
    void                *pv;
  } u_UnderflowTable;

  union {
    short int           *pi;
    unsigned short int  *pui;
    int                 *pn;
    unsigned int        *pun;
    INT4                *pl;
    UINT4               *pul;
    char                *pc;
    unsigned char       *puc;
    float               *pf;
    void                *pv;
  } u_OverflowTable;

  union {
    short int           *pi;
    unsigned short int  *pui;
    int                 *pn;
    unsigned int        *pun;
    INT4                *pl;
    UINT4               *pul;
    char                *pc;
    unsigned char       *puc;
    float               *pf;
    void                *pv;
  } u_BigOverflowTable;

  u_UnderflowTable.pv = NULL;
  u_OverflowTable.pv = NULL;
  u_BigOverflowTable.pv = NULL;

  int nNumUnderflows    = 0;
  int nNumOverflows     = 0;
  int nNumBigOverflows  = 0;
  int nUnderSize        = 0;
  int nOverSize         = 0;
  int nBigOverSize      = 0;
  int nBaseLineOffset   = 0;
  float fBaseLineOffset   = 0.0;

  int nFormat;
  nFormat = 0;
  nStat = poHeader->nGetValue("BS_FORMAT", &nFormat);
  if (101 == nFormat)
    {
      cout << "\nWARNING, BS FORMAT = 101 is not supported!\n"
           << "         please contact Rigaku/MSC for help!\n\n" << endl;
    }
  if (100 <= nFormat)
    {
      // New format, read number and type of the underflow and overflow tables
      // It is possible that the 3rd value does not exist in the header.

      (void) poHeader->nGetValue("BS_NOVERFL", 3, &a10nTemp[0]);
      nNumUnderflows   = a10nTemp[0];
      nNumOverflows    = a10nTemp[1];
      nNumBigOverflows = a10nTemp[2];

      (void) poHeader->nGetValue("BS_NPIXELB", 2, &a10nTemp[0]);
      nUnderSize   = a10nTemp[1];
      nOverSize    = 2;
      nBigOverSize = 4;
    }
  else
    {
      (void) poHeader->nGetValue("BS_NOVERFL", &nNumOverflows);
    }

  if (0 < nNumUnderflows)
    {
      // Read the binary underflow table.  It is a multiple of 16 bytes
      // Convert underflows to long int for use later

      nToRead = nUnderSize * nNumUnderflows;
      nToRead = ((nToRead + 15) / 16) * 16;

      // Do it like this to make sure long word boundary restrictions are
      // obeyed

      u_UnderflowTable.pul = new UINT4 [nToRead];
      (void) dskbr((int *)pnFile, (char *)u_UnderflowTable.pul,
		   &nToRead, &nStat);

      if (0 != nStat)
	{
	  cout << "ERROR reading underflow table!\n" << flush;
	}
      else if (eCPU_big_endian == nGetByteOrderCPU())
	{
	  // Need to swap bytes if 2 or 4 bytes per underflow pixel
	  
	  if (2 == nUnderSize)
	    {
	      unsigned short int *puiTemp;
	      puiTemp = u_UnderflowTable.pui;
	      i       = nNumUnderflows;
	      while (0 < i)  // Tests show this is 10-20% faster than the for loop.
		{
		  *puiTemp = ((*puiTemp & 0xff) << 8) | ((*puiTemp) >> 8);
		  puiTemp++;
		  i--;
		} 
	    }
	  else if (4 == nUnderSize)
	    {
	      for (i = 0; i < nNumUnderflows; i++)
		{
		  u_UnderflowTable.pul[i] = nSwapLong(u_UnderflowTable.pul[i]);
		}
	    }
	}

      // Now convert all underflows into 4 bytes for use later

      if (1 == nUnderSize)
	{
	  UINT4 *pulTemp;
	  unsigned char *pucTemp;
	  pucTemp = u_UnderflowTable.puc;
	  pulTemp = u_UnderflowTable.pul;

/*	  for (i = 0; i < 10; i++)
	    cout << "Char: " << i << ", " << (unsigned int)pucTemp[i] << endl;
	  for (i = nNumUnderflows-1; i >= nNumUnderflows-10; i--)
	    cout << "Char: " << i << ", " << (unsigned int)pucTemp[i] << endl;
*/
	  for (i = nNumUnderflows-1; i >= 0; i--)
	    {
	      pulTemp[i] = (UINT4) pucTemp[i];
	    }
/*
	  for (i = 0; i < 10; i++)
	    cout << "UINT4: " << i << ", " << pulTemp[i] << endl;
	  for (i = nNumUnderflows-1; i >= nNumUnderflows-10; i--)
	    cout << "UINT4: " << i << ", " << pulTemp[i] << endl;
*/
	}
      else if (2 == nUnderSize)
	{
	  UINT4 *pulTemp;
	  unsigned short int *puiTemp;
	  puiTemp = u_UnderflowTable.pui;
	  pulTemp = u_UnderflowTable.pul;

	  for (i = nNumUnderflows-1; i >= 0; i--)
	    {
	      pulTemp[i] = (UINT4) puiTemp[i];
	    }
	}
    }

  if ( (0 < nNumOverflows) && (100 > nFormat) )
    {
      // Old format overflow table

      char  a17cTemp[17];
      a17cTemp[16] = '\0';

      // We're gonna put both values and offsets into u_OverflowTables.pul

      u_OverflowTable.pul = new UINT4 [2 * nNumOverflows];
      UINT4 *pulOffset;
      pulOffset = &u_OverflowTable.pul[nNumOverflows];
      nOverSize    = 4;

      //  cout << "num overflows is " << nNumOverflows << endl;
      for (i = 0; (0 == nStat) && (i < nNumOverflows); i++)
	{
	  nToRead = 16;
	  (void) dskbr((int *)pnFile, (char *)a17cTemp, &nToRead, &nStat);
	  if (0 == nStat)
	    {
	      pulOffset[i] = atoi(&a17cTemp[9]);
	      a17cTemp[9] = '\0';
	      u_OverflowTable.pul[i] = atoi(a17cTemp);
	    }
	  else
	    {
	      cout << "WARNING in bs, error reading old-style ASCII overflow table!\n";
              cout << "   Required " << nNumOverflows << ", but error "
                      "at overflow number " << i << endl << flush;
	      nNumOverflows = i;
	    }
	}
    }

  if ( (0 < nNumOverflows) && (100 <= nFormat) )
    {
      // Read the new-style binary first overflow table (2-bytes per binary value)
      // File data is in little_endian format

      nToRead = nOverSize * nNumOverflows;
      nToRead = ( (nToRead + 15) / 16) * 16;
      u_OverflowTable.pul = new UINT4 [nToRead];
      (void) dskbr((int *)pnFile, (char *)u_OverflowTable.pul, &nToRead, &nStat);
      if (0 != nStat)
	{
	  cout << "ERROR reading 1st overflow table!\n" << flush;
	}
      else if (eCPU_big_endian == nGetByteOrderCPU())
	{
	  // Need to swap bytes (nOverSize == 2)

	  unsigned short int *puiTemp;
	  puiTemp = u_OverflowTable.pui;
	  i = nNumOverflows;
	  while (0 < i)  // Tests show this is 10-20% faster than the for loop.
	    {
	      *puiTemp = ((*puiTemp & 0xff) << 8) | ((*puiTemp) >> 8);
	      puiTemp++;
	      i--;
	    } 
	}

      // Now convert all overflows into 4 bytes for use later

      if (2 == nOverSize)
	{
	  UINT4 *pulTemp;
	  unsigned short int *puiTemp;
	  puiTemp = u_OverflowTable.pui;
	  pulTemp = u_OverflowTable.pul;

/*	  for (i = 0; i < 10; i++)
	    cout << "I2: " << i << ", " << puiTemp[i] << endl;
	  for (i = nNumOverflows-1; i >= nNumOverflows-10; i--)
	    cout << "I2: " << i << ", " << puiTemp[i] << endl;
*/
	  for (i = nNumOverflows-1; i >= 0; i--)
	    {
	      pulTemp[i] = (UINT4) puiTemp[i];
	    }
/*
	  for (i = 0; i < 10; i++)
	    cout << "I4: " << i << ", " << pulTemp[i] << endl;
	  for (i = nNumOverflows-1; i >= nNumOverflows-10; i--)
	    cout << "I4: " << i << ", " << pulTemp[i] << endl;
*/

	}
    }
  if ( (0 < nNumBigOverflows) && (100 <= nFormat) )
    {
      // Read the new-style binary second overflow table (4-bytes per binary value)
      // File data is in little_endian format

      nToRead = nBigOverSize * nNumBigOverflows;
      nToRead = ( (nToRead + 15) / 16) * 16;
      u_BigOverflowTable.pul = new UINT4 [nToRead];
      (void) dskbr((int *)pnFile, (char *)u_BigOverflowTable.pul, &nToRead, &nStat);
      if (0 != nStat)
	{
	  if ( (-1 == nStat) && (nNumBigOverflows < 30) )
	    {
	      nStat = 0;
	      cout << "INFO: 2nd overflow table missing or incomplete!\n" 
                   << "      Error ignored!\n" << flush;
	    }
	  else
	    {
	      cout << "ERROR reading 2nd overflow table!\n" << flush;
	    }
	}
      else if (eCPU_big_endian == nGetByteOrderCPU())
	{
	  // Need to swap bytes (nBigOverSize == 4)

	  for (i = 0; i < nNumBigOverflows; i++)
	    {
	      u_BigOverflowTable.pul[i] = nSwapLong(u_BigOverflowTable.pul[i]);
	    }
	}
    }

  // Get the baseline offset

  if (0 == poHeader->nGetValue("BS_NEXP", 5, &a10nTemp[0]))
    {
      nBaseLineOffset = a10nTemp[2];
      fBaseLineOffset = (float)nBaseLineOffset;
    }

  // OK, now have raw binary pixel values (as unsigned short int), underflow table, 
  // and overflow table(s)  ALL have the correct byte order
  // so convert to R-AXIS internal format unless the max value is < 65535

  float fPixel;

  if (100 > nFormat)
    {
      //cout << "Converting old format BS image to internal format.\n" << flush;
      // Old format
      // Should be no underflows
      // Convert to RAXIS2 style format
      // 1. Loop through original pixel values and pack those values above 32767
      // 2. Loop through overflow values and insert into original values if required

      puiTemp = puiData;
      for (i = 0; i < nDim0 * nDim1; i++)
	{
	  // What about the BaseLineOffset here?

	  if (32767 < *puiTemp)
	    {
	      fPixel = (float) *puiTemp + fBaseLineOffset;
	      *puiTemp = (unsigned short int) ((fPixel / 8.0) + 32768.0);
	    }
	  else
	    {
	      *puiTemp = *puiTemp + nBaseLineOffset;
	    }
	  puiTemp++;
	}

      UINT4 *pulOffset;
      pulOffset = &u_OverflowTable.pul[nNumOverflows];
      for (i = 0; i < nNumOverflows; i++)
	{
	  if (32768 > u_OverflowTable.pui[i])
	    puiData[pulOffset[i]] = (unsigned short int) u_OverflowTable.pul[i];
	  else
	    {
	      fPixel = (float) u_OverflowTable.pul[i];
	      puiData[pulOffset[i]] = (unsigned short int) 
		                        ((fPixel / 8.0) + 32768.0);
	    }
	}
    }
  else
    {
      // New format, must apply all the offsets, underflows and overflows

      puiTemp = puiData;
      UINT4 *pulOver, *pulUnder, *pulBigOver;
      pulOver    = u_OverflowTable.pul;
      pulUnder   = u_UnderflowTable.pul;
      pulBigOver = u_BigOverflowTable.pul;
      
      int nErrorUnder = 0;
      int nErrorOver = 0;
      int nErrorBigOver = 0;

      int nNumUnderflowsUsed;
      int nNumOverflowsUsed;
      int nNumBigOverflowsUsed;
      nNumBigOverflowsUsed  = nNumBigOverflows;
      nNumOverflowsUsed  = nNumOverflows;
      nNumUnderflowsUsed = nNumUnderflows;

      for (i = 0; i < nDim0 * nDim1; i++)
	{
	  if (0 == *puiTemp)
	    {
	      // Use the underflow table, and no baseline offset

	      if (0 < nNumUnderflows)
		{
		  nNumUnderflows--;
		  *puiTemp = (unsigned short int) *pulUnder++;
		}
	      else if (0 == nErrorUnder)
		{
		  //cout << "WARNING pixel with value of 0 and not enough underflows!\n" << flush;
		  nErrorUnder++;
		}
	    }
	  else if (255 == *puiTemp)
	    {
	      // Use the first overflow table

	      if (0 < nNumOverflows)
		{
		  nNumOverflows--;
		  fPixel = (float) *pulOver++; // + fBaseLineOffset;
		  if (65535  != fPixel)
		    fPixel += fBaseLineOffset;
		  else
		    {
		      if (0 < nNumBigOverflows)
			{
			  nNumBigOverflows--;
			  fPixel = (float) *pulBigOver++ + fBaseLineOffset;
			  if (262136.0 < fPixel)
			    fPixel = 262136.0;
			}
		      else
			{
			  cout << "ERROR not enough big overflows!\n" << flush;
			}
		    }
		  if (32767 < fPixel)
		    *puiTemp = (unsigned short int) ((fPixel / 8.0) + 32768.0);
		  else
		    *puiTemp = (unsigned short int) fPixel;
		}
	      else if (0 == nErrorOver)
		{
		  //cout << "WARNING pixel with value of 255 and not enough overflows!\n" << flush;
		  nErrorOver++;
		}
	    }
	  else if (65535 == *puiTemp)
	    {
	      // Use the second overflow table

	      if (0 < nNumBigOverflows)
		{
		  nNumBigOverflows--;
		  fPixel = (float) *pulBigOver++ + fBaseLineOffset;
		  if (262136.0 < fPixel)
		    fPixel = 262136.0;
		  if (32767 < fPixel)
		    *puiTemp = (unsigned short int) ((fPixel / 8.0) + 32768.0);
		  else
		    *puiTemp = (unsigned short int) fPixel;
		}
	      else
		{
		  cout << "ERROR not enough big overflows!\n" << flush;
		}
	    }
	  else if (32767 < *puiTemp)
	    {
	      fPixel = (float) *puiTemp + fBaseLineOffset;
	      *puiTemp = (unsigned short int) ((fPixel / 8.0) + 32768.0);
	    }
	  else
	    {
	      // Apply the baseline offset

	      *puiTemp = *puiTemp + nBaseLineOffset;
	    }
	  puiTemp++;
	}
    }

  if ( (0 != nNumOverflows) && (0 != nNumUnderflows) )
    {
      // We have not used up the right amount of over/under flows
      
      cout << "WARNING in bs, mis-match among over/under flows!\n";
      cout << "nNumOverflows = " << nNumOverflows << endl;
      cout << "nNumUnderflows = " << nNumUnderflows << endl << flush;
      nStat = -1;
    }
  if ( (0 != nStat) && ("" != sGetEnv("DTREK_SIEMENS_IGNORE")))
    {
      cout << "WARNING in bs, error reading file, but ignored!\n";
      nStat = 0;
    }
  if (NULL != u_UnderflowTable.pul)
    {
      delete [] u_UnderflowTable.pul;
      u_UnderflowTable.pul = NULL;
    }
  if (NULL != u_OverflowTable.pul)
    {
      delete [] u_OverflowTable.pul;
      u_OverflowTable.pul = NULL;
    }
  if (NULL != u_BigOverflowTable.pul)
    {
      delete [] u_BigOverflowTable.pul;
      u_BigOverflowTable.pul = NULL;
    }

  return (nStat);
}

int
nWriteSiemens(Cimage & oImage, const Cstring& sName)
{
  // Write an image to a Siemens style file
  // This is currently not-implemented. 
  // The code here is from the R-AXIS routine.
  
  RAXIS_header tRheader;
  int i, j, nStat;
  float afTemp[10];
  Cstring sPrefix;

  strncpy(tRheader.a10cDevice, "R-AXIS    ", 10);
  strncpy(tRheader.a10cVersion,"CCD       ", 10); 
  strncpy(tRheader.a20cCrystal,   "Unknown             ", 20);
  strncpy(tRheader.a12cCry_system,"Unknown     ", 12);
  tRheader.a6fCell[0]     = 1.0;
  tRheader.a6fCell[1]     = 1.0;
  tRheader.a6fCell[2]     = 1.0;
  tRheader.a6fCell[3]     = 90.0;
  tRheader.a6fCell[4]     = 90.0;
  tRheader.a6fCell[5]     = 90.0;

  // We should not use m_oHeader.nGetValue since it defeats encapsulation!!!

  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_CrystalPrefix) + D_K_CrystalXUnitCell, 6, tRheader.a6fCell);
  strncpy(tRheader.a12cSpace, "P1          ", 12);
  tRheader.fMosaic      = 0.20f;
  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_CrystalPrefix) + D_K_CrystalXMosaicity, &tRheader.fMosaic);

  strncpy(tRheader.a80cMemo, "memo                                                                            ", 80);
  strncpy(tRheader.a84cReserve1, "Reserve1                                                                            ", 84);
  strncpy(tRheader.a12cDate, "39-Jug-2001 ", 12);
  strncpy(tRheader.a20cOperatorname, "Unknown d*TREK      ", 20);
  strncpy(tRheader.a4cTarget, "Unk ", 4);

  tRheader.fWave        = 1.54178f;
  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_SourcePrefix) + 
				   D_K_Wavelength, 2, afTemp);
  if (0 == nStat) tRheader.fWave = afTemp[1];

  strncpy(tRheader.a20cMonochro, "Unknown             ", 20);
  tRheader.fMono_2      = 0.0;
  strncpy(tRheader.a20cCollimator, "Unknown             ", 20);
  strncpy(tRheader.a4cFilter, "Unk ", 4);
  tRheader.fCamera      = 100.0;

  nStat = oImage.m_oHeader.nGetValue(D_K_DetectorNames, &sPrefix);
  if (0 != nStat) sPrefix = "";
  nStat = oImage.m_oHeader.nGetValue(sPrefix + D_K_DetectorTranslation, 3, 
			      afTemp);
  if (0 == nStat) tRheader.fCamera = afTemp[2];

  tRheader.fKv          = 0.0;
  tRheader.fMa          = 0.0;
  tRheader.nCylinder    = 0;
  tRheader.fWeissenberg = 0.0;
  strncpy(tRheader.a12cFocus, "Unknown     ", 12);
  strncpy(tRheader.a80cXraymemo, "xraymemo                                                                        ", 80);
  strncpy(tRheader.a56cReserve2, "reserve2                                                ", 56);
  strncpy(tRheader.a4cSpindle, "Unk ", 4);
  strncpy(tRheader.a4cXray_axis, "Unk ", 4);

  tRheader.a3fPhi[0]    = 0.0;
  tRheader.a3fPhi[1]    = 0.0;
  tRheader.a3fPhi[2]    = 0.0;
  tRheader.nOsc         = 1;
  tRheader.fEx_time     = 0.0;
  tRheader.a2fXray[0]   = 100.0;  // X-ray beam position in pixels
  tRheader.a2fXray[1]   = 100.0;

  nStat = oImage.nGetDimensions(&i, &j);
  if (0 < nStat)
    {
      tRheader.a2nPix_num[0]  = i;
      tRheader.a2nPix_num[1]  = j;
      tRheader.a2fXray[0] = (float) i * 0.5;
      tRheader.a2fXray[1] = (float) j * 0.5;
    }

  tRheader.a3fCircle[0]  = 0.0;    // Crystal goniometer datum positions, omega
  tRheader.a3fCircle[1]  = 0.0;    // chi
  tRheader.a3fCircle[2]  = 0.0;    // Detector goniometer datum positions

  nStat = oImage.m_oHeader.nGetValue(D_K_Rotation, 10, afTemp);
  if (0 == nStat)
    {
      tRheader.a3fPhi[1] = afTemp[0];
      tRheader.a3fPhi[2] = afTemp[1];
      tRheader.fEx_time = afTemp[3] / 60.0;
      tRheader.nOsc     = (int) afTemp[4];    
    }

  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
				   3, afTemp);
  if (0 == nStat)
    {
      tRheader.a3fPhi[0]    = afTemp[2];
      tRheader.a3fCircle[0] = afTemp[0];
      tRheader.a3fCircle[1] = afTemp[1];
    }

  tRheader.fMu          = 0.0;

#ifndef RAXIS_HEADER_BINARY_STRESS_INFO 
  for (i = 0; i < 204; i++) 
      tRheader.a204cScanTemplate[i] = ' ';
#else
  tRheader.posn = 0;
  tRheader.side = 0;
  tRheader.slit = 0;
  
  for (i = 0; i < 64; i++) 
      tRheader.posz[i] = ' ';

  for (i = 0; i < 128; i++) 
      tRheader.psis[i] = ' ';
#endif

  tRheader.a2fPix_size[0] = 0.100f;
  tRheader.a2fPix_size[1] = 0.100f;

  
  nStat = oImage.m_oHeader.nGetValue(sPrefix + D_K_SpatialDistortionInfo,
				   4, afTemp);
  if (0 == nStat)
    {
      // Set pixel size and beam center from spatial distortion info

      tRheader.a2fXray[0]     = afTemp[0]; // X-ray beam position in pixels
      tRheader.a2fXray[1]     = afTemp[1];
      tRheader.a2fPix_size[0] = afTemp[2];
      tRheader.a2fPix_size[1] = afTemp[3];
    }
  else 
    {
      // Try to get pixel size from detector size and detector dimensions

      nStat = oImage.m_oHeader.nGetValue(sPrefix + D_K_DetectorSize, 2, afTemp);
      nStat = nStat +  oImage.m_oHeader.nGetValue(sPrefix + D_K_DetectorDimensions,
						2, &afTemp[2]);
      if (0 == nStat)
	{
	  tRheader.a2fPix_size[0] = afTemp[0] / afTemp[2];
	  tRheader.a2fPix_size[1] = afTemp[1] / afTemp[3];
	}
    }

  // First get minimum size of a record.  Then see if a larger size
  // is suggested by the header.   Images derived from R-AXIS
  // images will have a RAXIS_RECORD_SIZE keyword.  Others must have
  // records at least the size of the header.
    
  int nTemp;
  nTemp                 = 2 * tRheader.a2nPix_num[0];
  tRheader.a2nRecord[0] = max(nTemp,1024); // Record length minimum is 1024
  nStat = oImage.m_oHeader.nGetValue(sPrefix + D_K_RaxisRecordSize, &nTemp);
  if ( (0 == nStat) && (nTemp > tRheader.a2nRecord[0]) )
    {
      // Found a larger suggested record size, so use it
      tRheader.a2nRecord[0]   = nTemp;
    }
  tRheader.a2nRecord[1] = tRheader.a2nPix_num[1];
  tRheader.nRead_start  = 0;
  tRheader.nIP_num      = 1;
  tRheader.fRatio       = 1.0;
  tRheader.a2fFading[0] = 1.0;
  tRheader.a2fFading[1] = 1.0;
  strncpy(tRheader.a10cCpu, "IRIS      ", 10);
  strncpy(tRheader.a10cIp,  "Unknown   ", 10);
  tRheader.a3nDrxz[0]     = 0;
  tRheader.a3nDrxz[1]     = 0;
  tRheader.a3nDrxz[2]     = 0;

  //+ 3-Apr-2002 What to do with the new header structure elements?
  tRheader.fPixShiftOdd   = 0.0;
  tRheader.fIntRatioOdd   = 1.0;
  tRheader.nMagicNum      = 0;
  memset((void*)&tRheader.nMagicNum, 0, 180);

  //strncpy(tRheader.a180cReserve4, "Unknown                                                                                                                                                                             ", 180);
  //- 3-Apr-2002

  // Header is complete.  Now write the file.

  int nFile, nBytes, nSize, nPad, istat;

  nFile = 1;

  // Compute total size of image file in bytes.  This is needed for dskbow
  // when on VMS platforms.
  // Compute padding needed after header.

  nPad   = tRheader.a2nRecord[0] - sizeof(tRheader);
  nSize  = tRheader.a2nRecord[0] * (tRheader.a2nRecord[1] + 1);
  nBytes = sName.length();
//  (void) dskbow(&nFile, (const char *)sName, &nBytes, &nSize, &istat);
  (void) dskbow(&nFile, sName.string(), &nBytes, &nSize, &istat);

  if (0 != istat)
    {
      cout << "     Error opening file " << sName << "Error is: " 
           << istat << "\n";
      return (istat);
    }
  else
    cout << "     File " << sName << " successfully opened.\n";

  // See if data is already compressed in R-AXIS style and adjust
  // header accordingly
  
  int nAlreadyRAXISCompressed;
  nAlreadyRAXISCompressed = oImage.m_oHeader.nGetValue(
			      Cstring(D_K_RaxisCompressionRatio), 
			      &tRheader.fRatio);
  // Write the header

  nBytes = sizeof(tRheader);
  dskbw(&nFile, (char *)&tRheader, &nBytes, &istat);  

  // Write the padding required after the header, use the header for the pad

  if ( (0 == istat) && (0 < nPad) )
    {
      dskbw(&nFile, (char *)&tRheader, &nPad, &istat);  
    }

  // Figure out padding needed after every line

  nPad  = tRheader.a2nRecord[0] - sizeof(short int) * tRheader.a2nPix_num[0];

  unsigned short int *puiData, *puiTemp;
  if (0 != nAlreadyRAXISCompressed)
    {
      // Data was not already compressed into R-AXIS style, so compress it.
	
      // Get pointer to the data
      // We should be sure to transform values above 32767
      // to be (i / fRatio) - 32768.
      // unless fRatio == 1, then transform is no transform
      // We may also want to re-orient (rotate) the data.

      // Transfer image data to local array for packing in the RAXIS way

      float fPixel;
      puiData = new unsigned short int [  tRheader.a2nPix_num[0]
					* tRheader.a2nPix_num[1]
					* sizeof(unsigned short int)];
      puiTemp = puiData;
  
      for (j = 0; j < tRheader.a2nPix_num[1]; j++)
	{
	  for (i = 0; i < tRheader.a2nPix_num[0]; i++)
	    {
	      fPixel = (oImage.*oImage.prfGetPixel)(i, j);

	      if (fPixel > 32767.0) 
		{
		  if (tRheader.fRatio > 1.0)
		    {
		      *puiTemp = (unsigned short int) (fPixel / tRheader.fRatio
						       + 32767.0);
		    }
		  else
		    {
		      *puiTemp = (unsigned short int) fPixel;
		    }
		}
	      else
		{
		  *puiTemp = (unsigned short int) fPixel;
		}
	      puiTemp++;
	    }
	}
      puiTemp = puiData;
    }
  else
    {
      puiTemp = oImage.m_The_Data.pui;  // Access data directly
    }

  // Now write out image

  if (0 >= nPad)
    {
      // No padding needed, so write entire image in one shot
      nBytes = sizeof(short int) * tRheader.a2nPix_num[0] 
                                 * tRheader.a2nPix_num[1];
      dskbw(&nFile, (char *)puiTemp, &nBytes, &istat);
    }
  else
    {
      // Padding required, so write line by line and add padding
      nBytes = sizeof(short int) * tRheader.a2nPix_num[0] + nPad;
      for (i = 0; (i < tRheader.a2nPix_num[1]-1) && (0 == istat); i++)
	{
	  // Write out one line at a time, note that we pad not with 0's
	  // but with the data from the next line so that we do not need
	  // an extra call to dskbw for the padding.
	    
          dskbw(&nFile, (char *)puiTemp, &nBytes, &istat);
          puiTemp = puiTemp + tRheader.a2nPix_num[0];
	}
      // Write out last line

      if (0 == istat)
	{
	  nBytes = sizeof(short int) * tRheader.a2nPix_num[0];
	  dskbw(&nFile, (char *)puiTemp, &nBytes, &istat);
	  if (0 < nPad)
	    {
	      dskbw(&nFile, (char *)puiTemp, &nPad, &istat);
	    }
	}
    }

  if (0 != istat)
    cout << "     Error writing file " << sName << "!  Error is: " 
         << istat << "\n";
  else
    cout << "     Success writing file " << sName << "!\n";

  // Close the output file

  (void) dskbcw(&nFile, &nBytes);  // nBytes is a temp var here

  if (0 != nAlreadyRAXISCompressed) 
    delete [] puiData;  // Delete mem new'd above
  return (istat);
}

