//
// Copyright (c) 2007 Rigaku Americas Corp.
//
// Please see CBFLib disclaimers and acknowledgement at the bottom of this file.
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
#include "cbfdtrek.h"
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
nReadCBFHeader(const int* pnFile, const Cstring& rsFilename,
                 char *pcBuffer, Cimage_header *poHeader)
{
  // Read a possible "pseudo" CBF header that may or may not have 
  // been partially read already.
  // Convert the header to a d*TREK header.
  //
  // With CBF images sometimes they will be read in and converted to
  // unsigned short int (for example when used as a mask file) and
  // sometimes they will be read in as long int (signed)
  //
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
  // The CBF header is of unknown size. 

  int   i;           // Loop counter
  int   nStat;       // Local status flag
  int   nLen;        // Temp variable for lengths
  float a10fTemp[10];
  int   a10nTemp[10];
  int   nSerialNumber = 0;
  Cstring sTemp;
  static Cstring sPILATUS = "PILT_";
  tagCBFDTREK_header tCBFDTREK_header;

  //cout << "cbfdtrek called\n" << flush;
  
  // The calling routine should have read the first 512 bytes of the image
  // and have it available to copy into any header structure.
  // This routine should then read any remaining bytes of the header, 
  // then parse that header and construct a d*TREK header object.
  // 
  // It turns out that the CBF format is not really used in many cases in that
  // the header is not in CBF format.  Instead you get 
  // (a) no real header, 
  // (b) a different format for the header where items are simply embedded
  //     in a comment string, or
  // (c) something else.
  //

  memset((void*)&tCBFDTREK_header, 0, sizeof(tCBFDTREK_header));
  nLen = sizeof(tCBFDTREK_header.ac2048text);

  if (NULL != pcBuffer)
    {
      memcpy((void*)&tCBFDTREK_header, (void*)pcBuffer, 512);
      nLen = nLen - 512;
      dskbr((int *)pnFile, (char *)&tCBFDTREK_header.ac2048text[512], 
            &nLen, &nStat);
    }
  else
    {
      dskbr((int *)pnFile, (char *)&tCBFDTREK_header, &nLen, &nStat);
    }

  // Populate a d*TREK cbf temporary header with
  // Some defaults that are not zero, just in case

  tCBFDTREK_header.fPolarz      = 0.99;
  tCBFDTREK_header.fExpTime     = 1.0;
  tCBFDTREK_header.nSize1      = 100;
  tCBFDTREK_header.nSize2      = 100;
  tCBFDTREK_header.fPixSize0   = 0.172;
  tCBFDTREK_header.fPixSize1   = 0.172;
  tCBFDTREK_header.fWavelength = 1.000;
  tCBFDTREK_header.fDet2theta  = 0.000;
  tCBFDTREK_header.fBeam0      = tCBFDTREK_header.nSize1 * 0.5;
  tCBFDTREK_header.fBeam1      = tCBFDTREK_header.nSize2 * 0.5;
  tCBFDTREK_header.fDetDist    = 100.0;
  tCBFDTREK_header.fOmega      = 0.0;
  tCBFDTREK_header.fChi        = 0.0;
  tCBFDTREK_header.fKappa      = 0.0;
  tCBFDTREK_header.fPhi        = 0.0;
  tCBFDTREK_header.fRotStart   = 0.0;
  tCBFDTREK_header.fRotInc     = 1.0;

  // Now have either no header or ...
  // Now have about 2048 characters of the image file in memory.
  // Parse the header into tokens and place them into asHeader[] for use later.

  bool bBeamFound = FALSE;
  // TODO: Cut the string off before 2048 if a control-Z character is found first.

  sTemp = Cstring(tCBFDTREK_header.ac2048text);
  std::vector<Cstring> asHeader;
  sTemp.nListToVector(asHeader, "# ,;\r\t\n=()");
  for (i = 0; i < asHeader.size(); i++)
    {
      //cout << asHeader[i] << endl;
      if (    ("time:" == asHeader[i])
           || ("Exposure_time" == asHeader[i]) )
        tCBFDTREK_header.fExpTime    = atof(asHeader[i+1].string());
      else if ("Pixel_size" == asHeader[i])
        {
          tCBFDTREK_header.fPixSize0   = 1000.0f * atof(asHeader[i+1].string());
          tCBFDTREK_header.fPixSize1   = 1000.0f * atof(asHeader[i+4].string());
        }
      else if ("Wavelength" == asHeader[i])
          tCBFDTREK_header.fWavelength = atof(asHeader[i+1].string());
      else if ("Polarization" == asHeader[i])
          tCBFDTREK_header.fPolarz = atof(asHeader[i+1].string());
      else if ("Detector_2theta" == asHeader[i])
        {
          tCBFDTREK_header.fDet2theta = atof(asHeader[i+1].string());
        }
      else if ("Beam_xy" == asHeader[i])
        {
          bBeamFound = TRUE;
          tCBFDTREK_header.fBeam0     = atof(asHeader[i+1].string());
          tCBFDTREK_header.fBeam1     = atof(asHeader[i+2].string());
        }
      else if ("Detector_distance" == asHeader[i])
        {
          tCBFDTREK_header.fDetDist    = atof(asHeader[i+1].string());
          if ( ("m" == asHeader[i+2]) || (2.0 > tCBFDTREK_header.fDetDist) ) 
            tCBFDTREK_header.fDetDist *= 1000.0; // Convert m to mm;        
        }
      else if ("Start_angle" == asHeader[i])      
        {
          tCBFDTREK_header.fRotStart   = atof(asHeader[i+1].string());
        }
      else if ("Angle_increment" == asHeader[i])      
        tCBFDTREK_header.fRotInc     = atof(asHeader[i+1].string());
      else if ("Phi" == asHeader[i])      
        tCBFDTREK_header.fPhi        = atof(asHeader[i+1].string());
      else if ("Chi" == asHeader[i])      
        tCBFDTREK_header.fChi        = atof(asHeader[i+1].string());
      else if ("Kappa" == asHeader[i])      
        tCBFDTREK_header.fKappa      = atof(asHeader[i+1].string());
      else if ("X-Binary-Size-Fastest-Dimension:" == asHeader[i])      
        tCBFDTREK_header.nSize1      = atoi(asHeader[i+1].string());
      else if ("X-Binary-Size-Second-Dimension:" == asHeader[i])      
        tCBFDTREK_header.nSize2      = atoi(asHeader[i+1].string());
      else if ("Silicon" == asHeader[i])
        {
          tCBFDTREK_header.fSensorThickness  = 0.0;
          tCBFDTREK_header.fSensorThickness  = atof(asHeader[i+3].string());
          if (0.0 < tCBFDTREK_header.fSensorThickness)
            {
              //cout << "Silicon3: " << asHeader[i+3] << endl;
              // Only put the keyword in if the units are "m" 
              // and the thickness is positive
              // also use the string so that all precision is kept
              if ( ("thickness" == asHeader[i+2]) && ("m" == asHeader[i+4]) )
                poHeader->nReplaceValue("PILT_SENSOR_THICKNESS", asHeader[i+3]);
              else
                cout << "WARNING in cbfdtrek, error parsing Silicon thickness!\n";
            }
        }
      else if ("Detector:" == asHeader[i])
        {
          // Look for "S/N" token, then next token is serial number text.
          // Change the way serial numbers are parsed and presented
          // so that they do not have a dash or hyphen in them:
          // Multiply before the dash by 100 and append after the dash.
          // In d*TREK serial numbers are integers

          int     nMaxToken = 5;
          int     nFoundSN = 0;
          Cstring sParseDash = "";
          for (int j = 1; j < nMaxToken; j++)
            {
              if ("S/N"  == asHeader[i+j])
                nFoundSN = j+1; // The bit with the actual S/N is the NEXT token, so add 1
            }
          //cout << "Detector 1 " << asHeader[i] << endl;
          //cout << "Detector 2:" << asHeader[i+nFoundSN] << endl;
          if (0 < nFoundSN)
            {
              int     nBefore = 0;
              sParseDash = asHeader[i+nFoundSN].before('-');
              //cout << "sParseDash 1:" << sParseDash << endl;
              if (0 < sParseDash.length())
                {
                  // Has a dash or -
                  nBefore = atoi(sParseDash.string());
                  nBefore = 100 * nBefore;
                  sParseDash = asHeader[i+nFoundSN].after('-');
                }
              if (0 < sParseDash.length())
                if (sParseDash.contains('-'))
                  sParseDash = sParseDash.before('-');
              //cout << "sParseDash 2:" << sParseDash << endl;
              if (0 >= sParseDash.length())
                // Did/Does not have a dash
                sParseDash ="0000";
              sParseDash = Cstring(nBefore) + sParseDash;
              nSerialNumber = atoi(sParseDash.string());
            }
          //cout << "New serial number: " << nSerialNumber << endl;
          poHeader->nReplaceValue("DETECTOR_SN", nSerialNumber);
        }
    }

  // At this point, whatever comment was at the top of the file has been
  // read and parsed, so we have either the defaults if no real header
  // or something from the comments

  if (!bBeamFound)
    {
      // No direct beam info found so put beam at the center

      tCBFDTREK_header.fBeam0      = tCBFDTREK_header.nSize1 * 0.5;
      tCBFDTREK_header.fBeam1      = tCBFDTREK_header.nSize2 * 0.5;
    }

  tCBFDTREK_header.fOmega      = 0.0;
  a10nTemp[0] = tCBFDTREK_header.nSize1;
  a10nTemp[1] = tCBFDTREK_header.nSize2;
  poHeader->nReplaceValue(D_K_OriginalImageFormat, 
               Cimage_header::sGetFormatTypeName(enImage_format_CBFDTREK));
  (void) poHeader->nReplaceValue(Cstring (D_K_Type), Cstring("mad"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Dim), Cstring("2"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Size1), a10nTemp[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_Size2), a10nTemp[1]);

  if (2 == atoi(sGetEnv(D_K_DtrekCBFPixDataType).string()))
    (void) poHeader->nReplaceValue(Cstring (D_K_DataType), 
                                 Cstring(D_K_UnsignedShortInt));
  else
    (void) poHeader->nReplaceValue(Cstring (D_K_DataType), 
                                 Cstring(D_K_LongInt));

  (void) poHeader->nReplaceValue(Cstring (D_K_DtdisplayOrientation), 
                                 Cstring("-X+Y"));

  // CBFLib should return image pixel values in a 'native' byte order
  (void) poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
                                 sGetByteOrderCPU());
  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNumber), (int)1);
  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNames), sPILATUS);
  (void) poHeader->nReplaceValue(sPILATUS + D_K_DetectorDimensions, 2, a10nTemp);
  a10fTemp[0] = (float) a10nTemp[0] * tCBFDTREK_header.fPixSize0;
  a10fTemp[1] = (float) a10nTemp[1] * tCBFDTREK_header.fPixSize1;

  (void) poHeader->nReplaceValue(sPILATUS + D_K_DetectorSize, 2, a10fTemp, 5);
  (void) poHeader->nReplaceValue(sPILATUS + D_K_DetectorDescription,
                                 Cstring ("PILATUS conversion"));

  // The following may not be correct:

  (void) poHeader->nReplaceValue(sPILATUS + D_K_DetectorVectors,
                                 Cstring ("1 0 0 0 1 0"));

  (void) poHeader->nReplaceValue(sPILATUS + D_K_GonioNumValues,
                                 (int) 6);
  (void) poHeader->nReplaceValue(sPILATUS + D_K_GonioUnits,
                                 Cstring ("deg deg deg mm mm mm"));
  (void) poHeader->nReplaceValue(sPILATUS + D_K_GonioNames,
                 Cstring ("RotZ RotX/Swing RotY TransX TransY TransZ/Dist"));

  (void) poHeader->nReplaceValue(sPILATUS + D_K_GonioVectors,
                 Cstring ("0 0 1 1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1"));

  a10fTemp[0] = 0.0;
  a10fTemp[1] = tCBFDTREK_header.fDet2theta; 
  a10fTemp[2] = 0.0;
  a10fTemp[3] = 0.0;
  a10fTemp[4] = 0.0;
  a10fTemp[5] = tCBFDTREK_header.fDetDist; 
      
  (void) poHeader->nReplaceValue(sPILATUS + D_K_GonioValues, 6, a10fTemp, 3);

  (void) poHeader->nReplaceValue(sPILATUS + D_K_NonunfType, 
                                 Cstring(D_K_NonunfStateSimpleMask));
  (void) poHeader->nReplaceValue(sPILATUS + D_K_NonunfInfo, 
                                 "$(FirstScanImage)");
      
  (void) poHeader->nReplaceValue(sPILATUS + D_K_SpatialDistortionType,
                                 Cstring (D_K_SpatialTypeSimple));
  a10fTemp[0] = tCBFDTREK_header.fBeam0;
  a10fTemp[1] = tCBFDTREK_header.fBeam1;
  a10fTemp[2] = tCBFDTREK_header.fPixSize0;
  a10fTemp[3] = tCBFDTREK_header.fPixSize1;
  (void) poHeader->nReplaceValue(sPILATUS + D_K_SpatialDistortionInfo,
                                 4, a10fTemp, 5);
      
  (void) poHeader->nReplaceValue(sPILATUS + D_K_SpatialDistortionVectors,
                                 Cstring("1 0 0 -1"));

  (void) poHeader->nReplaceValue(Cstring (D_K_SaturatedValue), 
                                 (INT4)1048512);
  a10fTemp[0] = 1.0;
  a10fTemp[1] = tCBFDTREK_header.fWavelength;

  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePrefix) + 
                                 D_K_Wavelength, 2, a10fTemp);
  a10fTemp[0] = tCBFDTREK_header.fPolarz;
  a10fTemp[1] = 0.0;
  a10fTemp[2] = 1.0;
  a10fTemp[3] = 0.0;
  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
                                 4, a10fTemp);
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

  // Assume rotation is around Omega for now

  a10fTemp[0] = tCBFDTREK_header.fOmega;
  a10fTemp[1] = tCBFDTREK_header.fChi;
  a10fTemp[2] = tCBFDTREK_header.fPhi;

  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
                                 3, a10fTemp, 3);
  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + 
                                 D_K_GonioDescription,
                                 Cstring ("Pilatus 3-circle"));

  a10fTemp[0] = tCBFDTREK_header.fRotStart;
  a10fTemp[1] = tCBFDTREK_header.fRotStart + tCBFDTREK_header.fRotInc;
  Cstring sAxis = "Omega";
  if (a10fTemp[0] == a10fTemp[1])
    {
      cout << "WARNING: rotation start == rotation end." << endl;
    }

  a10fTemp[2] = tCBFDTREK_header.fRotInc;
  a10fTemp[3] = tCBFDTREK_header.fExpTime;
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

  // Change last 4 digits in rsFilename to ?s

  sTemp = rsFilename;
  int a3nSeqInfo[3] = {0, 1, 0};

  sTemp = sBuildScanTemplate(rsFilename, 4, &a3nSeqInfo[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_ScanTemplate), sTemp);

  (void) poHeader->nReplaceValue(Cstring (D_K_ScanSeqInfo), 3, a3nSeqInfo);

  (void) poHeader->nReplaceValue(Cstring(D_K_Compression), Cstring("CBF"));


  // Detector serial number specific changes to the header values

  if (3000104 == nSerialNumber)
    {
      poHeader->nReplaceValue(D_K_CrystalPrefix D_K_GonioVectors, 
                              "1 0 0 0.642788 -0.766044 0  1 0 0");
    }
  if ( (0 != nStat) && ("" != sGetEnv("DTREK_CBFDTREK_IGNORE")))
    {
      cout << "Warning in cbfdtrek, error reading header, but ignored!\n";
      nStat = 0;
    }
  //+JWP 2008-07-24
  // Back door to read and apply a header from another image to this image
  // The other image MUST have the same name EXCEPT instead of ending in .cbf
  // its name ends in .img

  sTemp = sGetEnv("DTREK_CBF_HEADER_OTHER");
  if ("" != sTemp)
    {
      Cstring sHeaderName;

      // Get the current scan sequence number and use it to create the name
      // of the file that contains the header you want to use.
      
      sHeaderName = rsFilename;
      sHeaderName = sHeaderName.substr(0, rsFilename.length()-4);
      sHeaderName += Cstring(".img");
      sTemp = sGetEnv("DTREK_CBF_HEADER_OTHERDIR");
      if ("" != sTemp)
        {
          sHeaderName = sTemp + sFileGetBasename(sHeaderName);
        }
      // Now that we have the name, get the header

      // ERROR: One cannot read another image header because of interference
      //        in the dskio routines (dskbr, dskbcr) without taking special
      //        precautions
      // Solution: read header as a string.

      char acHeader[2*4096];
      char acSize[16];
      memset(acHeader, 0, sizeof(acHeader));
      memset(acSize, 0, sizeof(acSize));
      
      FILE *pFile;
      pFile = fopen(sHeaderName, "r");
      if (!pFile)
        {
          cout << "ERROR opening other header file: " << sHeaderName << endl;
        }
      else
        {
          int nOtherRead;
          int nOtherSize;
          int jj;
          nOtherRead = fread(acHeader, 1, sizeof(acHeader), pFile);
          if (sizeof(acHeader) == nOtherRead)
            {
              // HEADER_BYTES should start at 15 and end at 18-22.
              // ii should be 15
              for (jj = 15; jj < 37; jj++)
                {
                  if (';' == acHeader[jj])
                    break;
                }
            }
          if (jj < 22)
            {
              strncpy(acSize, &acHeader[15], jj-15);
              nOtherSize = atoi(acSize);
            }
          else
            {
              nOtherSize = 0;
            }
          fclose(pFile);
          if (0 < nOtherSize)
            {
              cout << "INFO other header: " << sHeaderName << endl;
              acHeader[nOtherSize] = '\0';
              nStat = poHeader->nCopyMask(Cstring(acHeader), "CRYS*");
              nStat = poHeader->nCopyMask(Cstring(acHeader), "PLT*");
              nStat = poHeader->nCopyMask(Cstring(acHeader), "DET*");
              sPILATUS = "PLT_";
              (void) poHeader->nReplaceValue(sPILATUS + D_K_SpatialDistortionVectors,
                                             Cstring("0 1 1 0"));
              // Beam center needs changing nDim0 - beam1
              a10fTemp[0] = tCBFDTREK_header.fBeam0;
              a10fTemp[1] = tCBFDTREK_header.fBeam1;
              a10fTemp[2] = tCBFDTREK_header.fPixSize0;
              a10fTemp[3] = tCBFDTREK_header.fPixSize1;
              (void) poHeader->nGetValue(sPILATUS + D_K_SpatialDistortionInfo,
                                         4, a10fTemp);
              a10fTemp[0] = tCBFDTREK_header.nSize1 - a10fTemp[0];
              (void) poHeader->nReplaceValue(sPILATUS + D_K_SpatialDistortionInfo,
                                             4, a10fTemp);
              (void) poHeader->nReplaceValue(sPILATUS + D_K_SpatialBeamPosition,
                                             2, a10fTemp);

              nStat = poHeader->nCopyMask(Cstring(acHeader), "SOURCE*");
              nStat = poHeader->nCopyMask(Cstring(acHeader), "SCAN_ROT*");
              nStat = poHeader->nCopyMask(Cstring(acHeader), "ROT*");
              nStat = poHeader->nCopyMask(Cstring(acHeader), "OPT*");
              nStat = poHeader->nCopyMask(Cstring(acHeader), "COM*");
              (void) poHeader->nReplaceValue(Cstring(D_K_Compression), Cstring("CBF"));
            }
          else
            {
              cout << "WARNING in nReadCBFHeader(), other header not available!\n";
            }
        }
    }
  //-JWP 2008-07-24

  return (nStat);
}

int
nReadCBFData(const Cstring &rsFilename, void *pvData,
             Cimage_header *poHeader)
{
  // This is supposed to be called from an Cimage object.
  // Read a CBF-style image data.  
  // The image file header should already have been read.
  // TODO: is file open already?

  int nStat = 0;
  FILE *fp;

  cbf_handle hCBFile;

  int nBinary_id;
  unsigned int compression;
  int binary_id;
  size_t elsize;
  int elsigned;
  int elunsigned;
  size_t elements, elements_read;
  int minelement;
  int maxelement; 
  char *byteorder ="little_endian";
  size_t dim1;
  size_t dim2;
  size_t dim3; 
  size_t padding;

  int nPixSize = 2;
  Cstring sTemp;

  sTemp = Cstring(D_K_UnsignedShortInt);
  (void) poHeader->nGetValue(Cstring (D_K_DataType), &sTemp);
  if (Cstring(D_K_LongInt) == sTemp)
    nPixSize = 4;

  // Re-open file
  // Read the data

  if ( ( fp = fopen( rsFilename.string(), "rb" ) ) == NULL )
    {
      cout << "     Error opening file " << rsFilename << endl;
      nStat = 1;
      return (nStat);
    }
  // cbf_* routines return an integer

  nStat = cbf_make_handle(&hCBFile);
  if (0 != nStat) cout << "ERROR in cbf_make_handle\n" << flush;

  nStat = cbf_read_widefile(hCBFile, fp, MSG_DIGEST);
  if (0 != nStat) cout << "ERROR in cbf_read_widefile\n" << flush;

  nStat = cbf_find_tag(hCBFile, "_array_data.data");
  if (0 != nStat) cout << "ERROR in cbf_find_tag\n" << flush;

  nStat = cbf_find_column   (hCBFile, "data");
  if (0 != nStat) cout << "ERROR in cbf_find_column data\n" << flush;

  nStat = cbf_get_integerarrayparameters_wdims (hCBFile, &compression, 
                                                &nBinary_id, 
                                                &elsize, 
                                                &elsigned, 
                                                &elunsigned, 
                                                &elements, 
                                                &minelement, 
                                                &maxelement,
                                                (const char **) &byteorder, 
                                                &dim1, &dim2, &dim3, &padding);
  size_t nNumRead = 0;
  int *pnData = NULL;
  int nDim0, nDim1;
  poHeader->nGetDimensions(&nDim0, &nDim1);
  //dim1 = nDim0; dim2 = nDim1;
  if ( (dim1 != nDim0) || (dim2 != nDim1) )
    {
      cout << "ERROR in cbfdtrek: image dimensions do NOT match!\n" << flush;
      nStat = 1;
      
    }
  else
    {
      if (4 == nPixSize)
// If image is of type 'long int', then read directly into the image location
        pnData = (int *) pvData;
      else
// If image is of type 'unsigned short int', then do something else
        pnData = new int [nDim0 * nDim1];

      nStat =  cbf_get_integerarray (hCBFile, 
                                     &binary_id, (void *)pnData, 
                                     elsize, elsigned, 
                                     elements, &elements_read);
      nStat = 0;
      if (elements != elements_read)
        { 
          cout << "ERROR in cbfdtrek: cbf_get_integerarray failure!\n" << flush;
          nStat = 1;
        }
    }
  // cbf_free_handle hopefully closes the file?

  nStat = cbf_free_handle (hCBFile);

  // Close the file and free up memory that was used.
  // Not sure about should this be closed
  // call fclose() causes a core dump
  //  fclose (fp);

  // Transfer integer data to our image data
  // Note that PILATUS images have bad pixels, both cold and hot.  
  // One should use an active mask file to help flag these pixels.
  // The nonunf file could have these bad pixels flagged as well.

  int i;
  unsigned short int *pui;
  unsigned short int *puiData;
  puiData = (unsigned short int*) pvData;
  pui     = (unsigned short int*) pvData;

  int *pnTemp;
  int *pn;
  pnTemp = (int*)pnData;
  pn     = (int*)pnData;

  unsigned short int uiTemp;
  int nPixVal;
  int nPixValMin, nPixValMax;
  nPixValMin = nPixValMax = *pnTemp;

  // A comment about the ADD ONE mechanism:
  // Pilatus images are apparently signed long int.  Pixels values that
  // are negative are considered to be bad pixels.  Pixel values of 0 are 
  // normally considered to be good pixels, except in the known gaps of the
  // detector.  Some versions of the images do have 0 for the 'gap' pixels
  // which would normally be considered good.  
  //
  // In a trick to deal with the images there are 2 internal methods of 
  // storing them.  
  //
  // The first method is to store them as unsigned short int.
  // The negative pixel values present a problem, so all pixels values are made
  // positive by changing negative pixels to 0 and adding 1 to all pixel values
  // that are greater than or equal to 0.  The gaps are still set to 0.
  // Furthermore, pixel values above 
  // 32768 are compressed using the RAXIS compression schem.
  //
  // The second method is to keep the pixel values in their native format of
  // signed long int.  Negative pixels values are still set to 0 and 1 is
  // added to all other pixel values (except the gaps are set to 0).
  // Internally, a member variable of the Cimage class, namely
  // m_nPixValOffset, is set to 1.  When a pixel value is served up, the
  // m_nPixValOffset is subtracted giving original pixel value.
  //
  // In the future, the m_nPixValOffset may be subtracted from the 
  // unsigned short int internal format as well.
  //

  bool bAddOne = FALSE;
  bool bDivTwo = TRUE;
  int nAddNum = 0;
  bAddOne = ("" != sGetEnv("PILT_DO_ADD_1"));

  //int nDivNum = 1;
  //bDivTwo = ("" == sGetEnv("PILT_DO_NOT_DIV_2"));
  //if (bDivTwo) nDivNum = 2;

  if (bAddOne) nAddNum = 1; // * nDivNum;
  //std::cout << "BADDONE is " << bAddOne << std::endl <<  std::flush;
  //std::cout << "nAddNum is " << nAddNum << std::endl <<  std::flush;

  int nSatVal = 1048512;

  // NOTE: No need to transfer if image type is 'long int'
  //       but we may still need to manipulate the negative pixels
  
  if (1 == 1)
    {
      //std::cout << "NO EXTRA WORK\n" << std::endl <<  std::flush;
    }
  else
  { // Start of extra work
    //std::cout << "EXTRA WORK\n" << std::endl <<  std::flush;

    for (i = 0; i < nDim0 * nDim1; i++)
      {
        nPixVal = *pnTemp; 
        if (nPixVal > nPixValMax) nPixValMax = nPixVal;
        if (nPixVal < nPixValMin) nPixValMin = nPixVal;

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        nPixVal += nAddNum;   // setting nAddNum removes an 'if' from the loop
        //nPixVal /= nDivNum;
        *pnTemp++ = nPixVal;  // Get the adjusted added to number back in there

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        if (nSatVal < nPixVal) nPixVal = nSatVal;

        if (4 != nPixSize)
          {
            //2011-01-12 JWP+ Moved this check inside the 4 != nPixSize block
            if (0 > nPixVal)
              {
                // For negative pixel values, this will set them to 0
                //cout << "WARNING CBF pixval at I " << i << ": " << nPixVal << endl;
                nPixVal = -nAddNum; // Converts to a 0 when nAddNum is added below
              }

            // If needed, copy the long int data over to the
            // unsigned short int memory location

            if (32767 < nPixVal)
              {
                // Use an R-AXIS compression scheme for things larger than 32767
                // JWP: the +32 is to do implicit rounding!
                uiTemp = (unsigned short int)( ((nPixVal+32) / 32) + 32767);
              }
            else
              {
                uiTemp = (unsigned short int) nPixVal;
              }
            *pui++ = uiTemp;
          } // endif 4 != nPixSize
      }

    // I am not sure if the (bAddOne) is necessary.  We probably want to
    // make these gap pixels bad no matter what the stored image says they are.

    // Set the gap or seams to a value of nGapValue which is -1 or 0?

    long int nGapValue = -2;

    int nNumVertSections = 0; // 4
    int anColStart[] = {487, 981, 1475, 1969};
    //[i =487->493, 981->987, 1475->1481, 1969->1975, j=0->(nDim1-1)]
    //  490, 984, 1478, 1972
    for (i = 0; i < sizeof(anColStart)/sizeof(int); i++)
      if (nDim0 > anColStart[i])
        nNumVertSections++;

    int anRowStart[] = {195, 407, 619, 831, 1043, 1255, 1467, 1679, 1891, 2103, 2315 }; 
    int nNumHorzSections = 0; // 11
    for (i = 0; i < sizeof(anRowStart)/sizeof(int); i++)
      if (nDim1 > anRowStart[i])
        nNumHorzSections++;
    
    if ( (bAddOne) && (0 < nNumVertSections) || (0 < nNumHorzSections) )
      {
        // Now restore 0's or nGapValue in the seams.
        // There are nNumVertSections columns of 0s and nNumHorzSections rows of 0s

        int icol, irow, j;

        if (4 != nPixSize)
          {
            for (irow = 0; irow < nDim1; irow++) // For each row of the detector
              { 
                for (j = 0; j < nNumVertSections; j++)          // There are ? vertical sections
                  {
                    // Position pui to start of the section
                    pui = &puiData[anColStart[j] + irow * nDim0];
                    for (i = 0; i < 7; i++)      // Each section is 7 pixels wide
                      {
                        *pui++ = 0;
                      }
                  }
              }
          }
        else  // Already 4-byte integer data
          {
            for (irow = 0; irow < nDim1; irow++) // For each row of the detector
              { 
                for (j = 0; j < nNumVertSections; j++)  // There are ? vertical sections
                  {
                    // Position pn to start of the section
                    pn = &pnData[anColStart[j] + irow * nDim0];
                    for (i = 0; i < 7; i++)      // Each section is 7 pixels wide
                      {
                        *pn++ = nGapValue;
                      }
                  }
              }
          }
        // Now work on nNumHorzSections row-sections of 0s
        if (4 != nPixSize)
          {
            for (irow = 0; irow < nNumHorzSections; irow++) // Instead of 11, maybe use sizeof()?
              {
                pui = &puiData[anRowStart[irow]*nDim0];
                for (i = 0; i < 17 * nDim0; i++)  // Each section is 17 pixels wide
                  *pui++ = 0;
              }
          }
        else
          {
            for (irow = 0; irow < nNumHorzSections; irow++) // Instead of 11, maybe use sizeof()?
              {
                pn = &pnData[anRowStart[irow]*nDim0];
                for (i = 0; i < 17 * nDim0; i++)  // Each section is 17 pixels wide
                  *pn++ = nGapValue;
              }
          }
      }

    cout << "INFO: Min, Max pixel values in pseudo-CBF image: " 
         << nPixValMin << ", " << nPixValMax << endl;
  } // End of extra work

  if (4 != nPixSize)
    {
      // Free up any allocated memory that will not be used again
      delete [] pnData;
      pnData = NULL;
    }
  return (nStat);
}

/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE THE CBFLIB PACKAGE UNDER THE TERMS OF THE GPL *
 * WHILE YOU MAY ALTERNATIVE DISTRIBUTE THE API UNDER THE LGPL        *
 * YOU MAY ***NOT*** DISTRBUTE THIS PROGRAM UNDER THE LGPL            *
 *                                                                    *                                                                    *
 **********************************************************************/

/*************************** GPL NOTICES ******************************
 *                                                                    *
 * This program is free software; you can redistribute it and/or      *
 * modify it under the terms of the GNU General Public License as     *
 * published by the Free Software Foundation; either version 2 of     *
 * (the License, or (at your option) any later version.               *
 *                                                                    *
 * This program is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      *
 * GNU General Public License for more details.                       *
 *                                                                    *
 * You should have received a copy of the GNU General Public License  *
 * along with this program; if not, write to the Free Software        *
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA           *
 * 02111-1307  USA                                                    *
 *                                                                    *
 **********************************************************************/


/**********************************************************************
 *                                                                    *
 *                    Stanford University Notices                     *
 *  for the CBFlib software package that incorporates SLAC software   *
 *                 on which copyright is disclaimed                   *
 *                                                                    *
 * This software                                                      *
 * -------------                                                      *
 * The term 'this software', as used in these Notices, refers to      *
 * those portions of the software package CBFlib that were created by *
 * employees of the Stanford Linear Accelerator Center, Stanford      *
 * University.                                                        *
 *                                                                    *
 * Stanford disclaimer of copyright                                   *
 * --------------------------------                                   *
 * Stanford University, owner of the copyright, hereby disclaims its  *
 * copyright and all other rights in this software.  Hence, anyone    *
 * may freely use it for any purpose without restriction.             *
 *                                                                    *
 * Acknowledgement of sponsorship                                     *
 * ------------------------------                                     *
 * This software was produced by the Stanford Linear Accelerator      *
 * Center, Stanford University, under Contract DE-AC03-76SFO0515 with *
 * the Department of Energy.                                          *
 *                                                                    *
 * Government disclaimer of liability                                 *
 * ----------------------------------                                 *
 * Neither the United States nor the United States Department of      *
 * Energy, nor any of their employees, makes any warranty, express or *
 * implied, or assumes any legal liability or responsibility for the  *
 * accuracy, completeness, or usefulness of any data, apparatus,      *
 * product, or process disclosed, or represents that its use would    *
 * not infringe privately owned rights.                               *
 *                                                                    *
 * Stanford disclaimer of liability                                   *
 * --------------------------------                                   *
 * Stanford University makes no representations or warranties,        *
 * express or implied, nor assumes any liability for the use of this  *
 * software.                                                          *
 *                                                                    *
 * Maintenance of notices                                             *
 * ----------------------                                             *
 * In the interest of clarity regarding the origin and status of this *
 * software, this and all the preceding Stanford University notices   *
 * are to remain affixed to any copy or derivative of this software   *
 * made or distributed by the recipient and are to be affixed to any  *
 * copy of software made or distributed by the recipient that         *
 * contains a copy or derivative of this software.                    *
 *                                                                    *
 * Based on SLAC Software Notices, Set 4                              *
 * OTT.002a, 2004 FEB 03                                              *
 **********************************************************************/


/**********************************************************************
 *                                 NOTICE                             *
 * Creative endeavors depend on the lively exchange of ideas. There   *
 * are laws and customs which establish rights and responsibilities   *
 * for authors and the users of what authors create.  This notice     *
 * is not intended to prevent you from using the software and         *
 * documents in this package, but to ensure that there are no         *
 * misunderstandings about terms and conditions of such use.          *
 *                                                                    *
 * Please read the following notice carefully.  If you do not         *
 * understand any portion of this notice, please seek appropriate     *
 * professional legal advice before making use of the software and    *
 * documents included in this software package.  In addition to       *
 * whatever other steps you may be obliged to take to respect the     *
 * intellectual property rights of the various parties involved, if   *
 * you do make use of the software and documents in this package,     *
 * please give credit where credit is due by citing this package,     *
 * its authors and the URL or other source from which you obtained    *
 * it, or equivalent primary references in the literature with the    *
 * same authors.                                                      *
 *                                                                    *
 * Some of the software and documents included within this software   *
 * package are the intellectual property of various parties, and      *
 * placement in this package does not in any way imply that any       *
 * such rights have in any way been waived or diminished.             *
 *                                                                    *
 * With respect to any software or documents for which a copyright    *
 * exists, ALL RIGHTS ARE RESERVED TO THE OWNERS OF SUCH COPYRIGHT.   *
 *                                                                    *
 * Even though the authors of the various documents and software      *
 * found here have made a good faith effort to ensure that the        *
 * documents are correct and that the software performs according     *
 * to its documentation, and we would greatly appreciate hearing of   *
 * any problems you may encounter, the programs and documents any     *
 * files created by the programs are provided **AS IS** without any   *
 * warranty as to correctness, merchantability or fitness for any     *
 * particular or general use.                                         *
 *                                                                    *
 * THE RESPONSIBILITY FOR ANY ADVERSE CONSEQUENCES FROM THE USE OF    *
 * PROGRAMS OR DOCUMENTS OR ANY FILE OR FILES CREATED BY USE OF THE   *
 * PROGRAMS OR DOCUMENTS LIES SOLELY WITH THE USERS OF THE PROGRAMS   *
 * OR DOCUMENTS OR FILE OR FILES AND NOT WITH AUTHORS OF THE          *
 * PROGRAMS OR DOCUMENTS.                                             *
 **********************************************************************/

/**********************************************************************
 *                                                                    *
 *                           The IUCr Policy                          *
 *      for the Protection and the Promotion of the STAR File and     *
 *     CIF Standards for Exchanging and Archiving Electronic Data     *
 *                                                                    *
 * Overview                                                           *
 *                                                                    *
 * The Crystallographic Information File (CIF)[1] is a standard for   *
 * information interchange promulgated by the International Union of  *
 * Crystallography (IUCr). CIF (Hall, Allen & Brown, 1991) is the     *
 * recommended method for submitting publications to Acta             *
 * Crystallographica Section C and reports of crystal structure       *
 * determinations to other sections of Acta Crystallographica         *
 * and many other journals. The syntax of a CIF is a subset of the    *
 * more general STAR File[2] format. The CIF and STAR File approaches *
 * are used increasingly in the structural sciences for data exchange *
 * and archiving, and are having a significant influence on these     *
 * activities in other fields.                                        *
 *                                                                    *
 * Statement of intent                                                *
 *                                                                    *
 * The IUCr's interest in the STAR File is as a general data          *
 * interchange standard for science, and its interest in the CIF,     *
 * a conformant derivative of the STAR File, is as a concise data     *
 * exchange and archival standard for crystallography and structural  *
 * science.                                                           *
 *                                                                    *
 * Protection of the standards                                        *
 *                                                                    *
 * To protect the STAR File and the CIF as standards for              *
 * interchanging and archiving electronic data, the IUCr, on behalf   *
 * of the scientific community,                                       *
 *                                                                    *
 * * holds the copyrights on the standards themselves,                *
 *                                                                    *
 * * owns the associated trademarks and service marks, and            *
 *                                                                    *
 * * holds a patent on the STAR File.                                 *
 *                                                                    *
 * These intellectual property rights relate solely to the            *
 * interchange formats, not to the data contained therein, nor to     *
 * the software used in the generation, access or manipulation of     *
 * the data.                                                          *
 *                                                                    *
 * Promotion of the standards                                         *
 *                                                                    *
 * The sole requirement that the IUCr, in its protective role,        *
 * imposes on software purporting to process STAR File or CIF data    *
 * is that the following conditions be met prior to sale or           *
 * distribution.                                                      *
 *                                                                    *
 * * Software claiming to read files written to either the STAR       *
 * File or the CIF standard must be able to extract the pertinent     *
 * data from a file conformant to the STAR File syntax, or the CIF    *
 * syntax, respectively.                                              *
 *                                                                    *
 * * Software claiming to write files in either the STAR File, or     *
 * the CIF, standard must produce files that are conformant to the    *
 * STAR File syntax, or the CIF syntax, respectively.                 *
 *                                                                    *
 * * Software claiming to read definitions from a specific data       *
 * dictionary approved by the IUCr must be able to extract any        *
 * pertinent definition which is conformant to the dictionary         *
 * definition language (DDL)[3] associated with that dictionary.      *
 *                                                                    *
 * The IUCr, through its Committee on CIF Standards, will assist      *
 * any developer to verify that software meets these conformance      *
 * conditions.                                                        *
 *                                                                    *
 * Glossary of terms                                                  *
 *                                                                    *
 * [1] CIF:  is a data file conformant to the file syntax defined     *
 * at http://www.iucr.org/iucr-top/cif/spec/index.html                *
 *                                                                    *
 * [2] STAR File:  is a data file conformant to the file syntax       *
 * defined at http://www.iucr.org/iucr-top/cif/spec/star/index.html   *
 *                                                                    *
 * [3] DDL:  is a language used in a data dictionary to define data   *
 * items in terms of "attributes". Dictionaries currently approved    *
 * by the IUCr, and the DDL versions used to construct these          *
 * dictionaries, are listed at                                        *
 * http://www.iucr.org/iucr-top/cif/spec/ddl/index.html               *
 *                                                                    *
 * Last modified: 30 September 2000                                   *
 *                                                                    *
 * IUCr Policy Copyright (C) 2000 International Union of              *
 * Crystallography                                                    *
 **********************************************************************/
