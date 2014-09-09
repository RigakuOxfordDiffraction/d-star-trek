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

#include <string.h>
#include "dtreksys.h"
#include "marccd.h"
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
nReadMARCCDHeader(const int* pnFile, const Cstring& rsFilename,
                 char *pcBuffer, Cimage_header *poHeader)
{
  // Read a MARCCD header that may or may not have been partially read already.
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
  // The MARCCD header is always 4096 bytes long.  It contains alot of binary
  // data that may be big or little endian.  So we have to figure this out,
  // maybe swap bytes.

  tagMARCCD_header tMARCCD_header;

  int   i;           // Loop counter
  int   nStat;       // Local status flag
  int   nLen;        // Temp variable for lengths
  float a10fTemp[10];
  int   a10nTemp[10];
  int   nSerialNumber = 0;
  Cstring sTemp;
  static Cstring sMARCCD = "MARCCD_";
  
  nLen = sizeof(tMARCCD_header);
  if (NULL != pcBuffer)
    {
      memcpy((void*)&tMARCCD_header, (void*)pcBuffer, 512);
      nLen = nLen - 512;
      dskbr((int *)pnFile, (char *)&tMARCCD_header.acTIFF[512], &nLen, &nStat);
    }
  else
    {
      dskbr((int *)pnFile, (char *)&tMARCCD_header, &nLen, &nStat);
    }

  if (   ( MAR_BIG_ENDIAN != tMARCCD_header.header_byte_order)
      && (MAR_LITTLE_ENDIAN != tMARCCD_header.header_byte_order) )
    tMARCCD_header.header_byte_order
      = nSwapLong(tMARCCD_header.header_byte_order);
      
  if (   (   (MAR_BIG_ENDIAN == tMARCCD_header.header_byte_order)
          && (eCPU_little_endian == nGetByteOrderCPU()))
      || (   (   (MAR_LITTLE_ENDIAN == tMARCCD_header.header_byte_order)
              && (eCPU_big_endian == nGetByteOrderCPU()))
          )
      )
    {
      // Swap bytes in the header
  
      int *pnTemp;
      pnTemp = (int*)&tMARCCD_header;
      for (i = 0; i < sizeof(tMARCCD_header)/sizeof(int); i++)
        {
          *pnTemp = nSwapLong(*pnTemp);
          pnTemp++;
        }
    }

  a10nTemp[0] = tMARCCD_header.nfast;
  a10nTemp[1] = tMARCCD_header.nslow;
  if (   (1 > a10nTemp[0]) || (10000 < a10nTemp[0])
      || (1 > a10nTemp[1]) || (10000 < a10nTemp[1]) )
    {
      // It is not really a MAR CCD image
      nStat = -100;
      return (nStat);
    }

  // Try to figure out which detector this image came from based on the
  // comments in the image header.  Adjust vectors later on based on this.

  sTemp = Cstring(tMARCCD_header.file_comments);
  sTemp = sTemp.after("Detector Serial Number = ");
  sTemp = sTemp.before("\n");
  nSerialNumber = atoi(sTemp.string());
  //cout << "DetSerialnum: >" << sTemp << "< " << nSerialNumber << endl;
  
  poHeader->nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_MAR_CCD));

  (void) poHeader->nReplaceValue(Cstring (D_K_Type), Cstring("mad"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Dim), Cstring("2"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Size1), a10nTemp[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_Size2), a10nTemp[1]);
  (void) poHeader->nReplaceValue(Cstring (D_K_DataType), 
                                 Cstring("unsigned short int"));
  (void) poHeader->nReplaceValue(Cstring (D_K_DtdisplayOrientation), 
                                 Cstring("-X+Y"));

  if (MAR_LITTLE_ENDIAN == tMARCCD_header.data_byte_order)
    {
      poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
                              Cstring (D_K_LittleEndian));
    }
  else
    {
      poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
                              Cstring (D_K_BigEndian));
    }
  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNumber), (int)1);
  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNames), sMARCCD);
  (void) poHeader->nReplaceValue(sMARCCD + D_K_DetectorDimensions, 2, a10nTemp);

  if ( (0 >= tMARCCD_header.pixelsize_x) && (tMARCCD_header.nfast == 3072) )
    {
      cout << "WARNING: pixel size is: " << tMARCCD_header.pixelsize_x << endl;
      tMARCCD_header.pixelsize_x = 225 * 1000000 / 3072;
    }
  if ( (0 >= tMARCCD_header.pixelsize_y) && (tMARCCD_header.nfast == 3072) )
    {
      tMARCCD_header.pixelsize_y = 225 * 1000000 / 3072;
    }

  a10fTemp[0] = (float) a10nTemp[0] * (float)tMARCCD_header.pixelsize_x * 0.000001;
  a10fTemp[1] = (float) a10nTemp[1] * (float)tMARCCD_header.pixelsize_y * 0.000001;

  (void) poHeader->nReplaceValue(sMARCCD + D_K_DetectorSize, 2, a10fTemp, 5);
  (void) poHeader->nReplaceValue(sMARCCD + D_K_DetectorDescription,
                                 Cstring ("MARCCD conversion"));

  // The following may not be correct:

  (void) poHeader->nReplaceValue(sMARCCD + D_K_DetectorVectors,
                                 Cstring ("1 0 0 0 1 0"));

  (void) poHeader->nReplaceValue(sMARCCD + D_K_GonioNumValues,
                                 (int) 6);
  (void) poHeader->nReplaceValue(sMARCCD + D_K_GonioUnits,
                                 Cstring ("deg deg deg mm mm mm"));
  (void) poHeader->nReplaceValue(sMARCCD + D_K_GonioNames,
                 Cstring ("RotZ RotX/Swing RotY TransX TransY TransZ/Dist"));

  (void) poHeader->nReplaceValue(sMARCCD + D_K_GonioVectors,
                 Cstring ("0 0 1 1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1"));

  a10fTemp[0] = 0.0;
  a10fTemp[1] = (float)tMARCCD_header.start_twotheta * 0.001;
  a10fTemp[2] = 0.0;
  a10fTemp[3] = 0.0;
  a10fTemp[4] = 0.0;
  a10fTemp[5] = (float)tMARCCD_header.xtal_to_detector * (float) 0.001; 
      
  if (0.0 >= a10fTemp[5])
    {
      a10fTemp[5] = (float)tMARCCD_header.start_xtal_to_detector 
        * (float) 0.001; 
      if (0.0 >= a10fTemp[5])
        a10fTemp[5] = 100.0;
    }

  (void) poHeader->nReplaceValue(sMARCCD + D_K_GonioValues, 6, a10fTemp, 3);

  (void) poHeader->nReplaceValue(sMARCCD + D_K_NonunfType, 
                                 Cstring(D_K_NonunfStateSimpleMask));
  (void) poHeader->nReplaceValue(sMARCCD + D_K_NonunfInfo, 
                                 "$(FirstScanImage)");
      
  (void) poHeader->nReplaceValue(sMARCCD + D_K_SpatialDistortionType,
                                 Cstring (D_K_SpatialTypeSimple));
  a10fTemp[0] = (float)tMARCCD_header.beam_x / (float)0.001;
  if (0.0 >= a10fTemp[0])
    a10fTemp[0] = 0.5 * (float)tMARCCD_header.nfast;
  a10fTemp[1] = (float)tMARCCD_header.beam_y / (float)0.001;
  if (0.0 >= a10fTemp[1])
    a10fTemp[1] = 0.5 * (float)tMARCCD_header.nslow;
  a10fTemp[2] = (float)tMARCCD_header.pixelsize_x * 0.000001;
  a10fTemp[3] = (float)tMARCCD_header.pixelsize_y * 0.000001;
  if (tMARCCD_header.nfast*2 <= (int)a10fTemp[0])
    a10fTemp[0] = a10fTemp[0]* 0.000001;
  if (tMARCCD_header.nslow*2 <= (int)a10fTemp[1])
    a10fTemp[1] = a10fTemp[1]* 0.000001;
  if (3 == nSerialNumber)
    a10fTemp[0] = (float)tMARCCD_header.nfast - a10fTemp[0];

  //+2012-Jul-02 JWP
  // Look at the beam center and pixel sizes.  If the beam center appears to 
  // be in HKL-3000 format, calculate the format for d*TREK.  
  // Characteristics of HKL format:  In millimeters and X & Y switched.
  // So if conversion puts beam center closer to image center, use it.
  float a2fConvertedBeamCenter[2];
  a2fConvertedBeamCenter[0] = a10fTemp[1] / a10fTemp[3];
  a2fConvertedBeamCenter[1] = a10fTemp[0] / a10fTemp[2];
  if ( 
      (ABS(a2fConvertedBeamCenter[0] - (float) tMARCCD_header.nfast * 0.5)
       < ABS(a10fTemp[0] - (float) tMARCCD_header.nfast * 0.5))
      &&
      (ABS(a2fConvertedBeamCenter[1] - (float) tMARCCD_header.nslow * 0.5)
       < ABS(a10fTemp[1] - (float) tMARCCD_header.nslow * 0.5))
      )
    {
      a10fTemp[0] = a2fConvertedBeamCenter[0];
      a10fTemp[1] = a2fConvertedBeamCenter[1];
    }
  //-2012-Jul-02 JWP

  (void) poHeader->nReplaceValue(sMARCCD + D_K_SpatialDistortionInfo,
                                 4, a10fTemp, 5);
      

  (void) poHeader->nReplaceValue(sMARCCD + D_K_SpatialDistortionVectors,
                                 Cstring("1 0 0 -1"));

  (void) poHeader->nReplaceValue(Cstring (D_K_SaturatedValue), 
                                 (INT4)tMARCCD_header.saturated_value);
  a10fTemp[0] = 1.0;
  a10fTemp[1] = (float)tMARCCD_header.source_wavelength * 0.00001;
  if (0.0 >= a10fTemp[1]) 
    a10fTemp[1] = 1.0000;
  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePrefix) + 
                                 D_K_Wavelength, 2, a10fTemp);
  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
                                 Cstring ("0.99 0 1 0"));
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
  if (3 == nSerialNumber)
    {
      // This thought to be a MARCCD at SER-CAT BM line
      (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
//+JWP 2010-04-13
//                               Cstring ("-1 0 0 0 0 1 -1 0 0"));
                                 Cstring ("-1 0 0 0 0 1 1 0 0"));
//-JWP 2010-04-13
    }
//+JWP 2010-04-13
  else if (7 == nSerialNumber)
    {
      // This thought to be a MARCCD at SER-CAT ID line
      (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
                                 Cstring ("-1 0 0 0 0 1 -1 0 0"));
    }
//-JWP 2010-04-13
  else
    (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
                                 Cstring ("1 0 0 0 0 1 1 0 0"));


  if (   ("" != sGetEnv("DTREK_MARCCD_USEPHI") )
      || (19 == nSerialNumber)         // Serial number 19 is APS ???
      || (25 == nSerialNumber) )       // Serial number 25 is APS 21-ID-G w/MD2
    tMARCCD_header.rotation_axis = 4;  // Force use of phi axis
  if ("" != sGetEnv("DTREK_MARCCD_USEOMEGA"))
    tMARCCD_header.rotation_axis = 1;  // Force use of omega axis
  if (4 == tMARCCD_header.rotation_axis) 
    {
      a10fTemp[0] = (float)tMARCCD_header.start_omega * 0.001;
      a10fTemp[1] = (float)tMARCCD_header.start_chi   * 0.001;
      a10fTemp[2] = 0.0;
    }
  else
    {
      // Assume Omega for now
      a10fTemp[0] = 0.0;
      a10fTemp[1] = (float)tMARCCD_header.start_chi   * 0.001;
      a10fTemp[2] = (float)tMARCCD_header.start_phi   * 0.001;
    }
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
                                 3, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + 
                                 D_K_GonioDescription,
                                 Cstring ("Single w/2-arc goniometer head"));

  Cstring sAxis;
  if (4 == tMARCCD_header.rotation_axis)
    {
      a10fTemp[0] = (float)tMARCCD_header.start_phi * 0.001;
      a10fTemp[1] = (float)tMARCCD_header.end_phi   * 0.001;
      if (    (a10fTemp[0] == a10fTemp[1])
	   && (0.0 != tMARCCD_header.rotation_range))
	{
          a10fTemp[1] = (float)(  tMARCCD_header.start_phi 
                                  + tMARCCD_header.rotation_range)  * 0.001;
	}
      else if ( (0 ==   tMARCCD_header.end_phi)
          || (NULL != strstr(tMARCCD_header.file_comments, "Number = 38")) )
        {
          a10fTemp[1] = (float)(  tMARCCD_header.start_phi 
                                  + tMARCCD_header.rotation_range)  * 0.001;
        }
      sAxis = "Phi";
    }
  else
    {
      if (1 != tMARCCD_header.rotation_axis)
        {
          cout << "WARNING: rotation axis neither phi nor omega." << endl;
        }
      a10fTemp[0] = (float)tMARCCD_header.start_omega * 0.001;
      a10fTemp[1] = (float)tMARCCD_header.end_omega   * 0.001;
      if (0 ==   tMARCCD_header.end_omega)
        {
          a10fTemp[1] = (float)(  tMARCCD_header.start_omega 
                                  + tMARCCD_header.rotation_range)  * 0.001;
        }
      sAxis = "Omega";
    }
  if (a10fTemp[0] == a10fTemp[1])
    {
      cout << "WARNING: rotation start == rotation end." << endl;

      cout << "MARhead: start_phi, end_phi, start_delta, end_delta: "
           << (float)tMARCCD_header.start_phi * 0.001 << ", "
           << (float)tMARCCD_header.end_phi * 0.001 << ", "
           << (float)tMARCCD_header.start_delta * 0.001 << ", "
           << (float)tMARCCD_header.end_delta * 0.001 << ", " 
           << endl;
      cout << "MARhead: start_omega, end_omega, start_chi, end_chi: "
           << (float)tMARCCD_header.start_omega * 0.001 << ", "
           << (float)tMARCCD_header.end_omega * 0.001 << ", "
           << (float)tMARCCD_header.start_chi * 0.001 << ", "
           << (float)tMARCCD_header.end_chi * 0.001 << ", " 
           << endl;
      cout << "MARhead: rotation_axis, rotation_range: "
           << tMARCCD_header.rotation_axis << ", " 
           << tMARCCD_header.rotation_range << endl;
    }

  a10fTemp[2] = a10fTemp[1] - a10fTemp[0];
  a10fTemp[3] = (float)tMARCCD_header.exposure_time * 0.001;
  a10fTemp[4] = 1;
  a10fTemp[5] = 0.0;
  a10fTemp[6] = 0.0;
  a10fTemp[7] = 0.0;
  a10fTemp[8] = 0.0;
  a10fTemp[9] = 0.0;
  (void) poHeader->nReplaceValue(Cstring (D_K_Rotation), 10, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_Rotation,
                                 10, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_RotAxisName), sAxis);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + 
                                 D_K_RotAxisName, sAxis);
  a10fTemp[0] = 1.0;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 0.0;

  (void) poHeader->nReplaceValue(Cstring (D_K_RotVector),
                                 3, &a10fTemp[0], 4);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_RotVector,
                                 3, &a10fTemp[0], 4);

  // Determine scan template and seq info from the filename

  // Change last 3 digits in rsFilename to ?s

  sTemp = rsFilename;
  int a3nSeqInfo[3] = {0, 1, 0};

  sTemp = sBuildScanTemplate(rsFilename, 4, &a3nSeqInfo[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanTemplate), sTemp);

  (void) poHeader->nReplaceValue(Cstring (D_K_ScanSeqInfo), 3, a3nSeqInfo);

  if ( (0 != nStat) && ("" != sGetEnv("DTREK_RAXIS_IGNORE")))
    {
      cout << "Warning in marccd, error reading header, but ignored!\n";
      nStat = 0;
    }
  return (nStat);
}
