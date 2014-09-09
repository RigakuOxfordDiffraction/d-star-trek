//
// Copyright (c) 2002 Rigaku/MSC
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
#include "raxis.h"
#include "winbmp.h"
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
nReadWindowsBMPHeader(const int* pnFile, const Cstring& rsFilename,
		   char *pcBuffer, Cimage_header *poHeader)
{
  // Build a temporary header, this doesn't actually read anything
  // things are read in nReadWindowsBMPData()

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
  float a10fTemp[10];
  int   a10nTemp[10];

  static Cstring sWinBMP = "BMP_";
  
  BITMAPINFOHEADER  bmpih;
  BITMAPFILEHEADER  bmpfh;
  int nDataBytes, nHeaderBytes;
  int nColors;

  memcpy(&bmpfh, pcBuffer, sizeof(BITMAPFILEHEADER));
  memcpy(&bmpih, pcBuffer+sizeof(BITMAPFILEHEADER), sizeof(BITMAPINFOHEADER));

  if (bmpih.biClrUsed != 0)
	  nColors = bmpih.biClrUsed;
  else
	  nColors = 1 << bmpih.biBitCount;

  if ( ( bmpih.biBitCount != 24 ) && ( bmpih.biBitCount != 8 ) )
  {
	return 1;
  }

  nDataBytes = bmpih.biBitCount / 8;
  nHeaderBytes = bmpfh.bfOffBits;

  //cout << "about to build nReadWinBMPHeader" << endl;
  nStat = 0;

  a10nTemp[0] = bmpih.biWidth;
  a10nTemp[1] = bmpih.biHeight;
  
  poHeader->nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_WINBMP));

  (void) poHeader->nReplaceValue(Cstring (D_K_Type), Cstring("mad"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Dim), Cstring("2"));
  (void) poHeader->nReplaceValue(Cstring (D_K_Size1), a10nTemp[0]);
  (void) poHeader->nReplaceValue(Cstring (D_K_Size2), a10nTemp[1]);
  (void) poHeader->nReplaceValue(Cstring (D_K_DataType), 
				 Cstring(D_K_Compressed));
  (void) poHeader->nReplaceValue(Cstring (D_K_Compression), 
				 Cstring(D_K_WinBMP));

  (void) poHeader->nReplaceValue(Cstring (D_K_DtdisplayOrientation), 
				 Cstring("-X-Y"));


  poHeader->nReplaceValue(Cstring (D_K_ByteOrder), 
			  Cstring (D_K_BigEndian));

  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNumber), (int)1);
  (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNames), sWinBMP);
  (void) poHeader->nReplaceValue(sWinBMP + D_K_DetectorDimensions, 2, a10nTemp);

  a10fTemp[0] = 0.0105f * bmpih.biWidth;
  a10fTemp[1] = 0.0105f * bmpih.biHeight;
  (void) poHeader->nReplaceValue(sWinBMP + D_K_DetectorSize, 2, a10fTemp, 5);
  (void) poHeader->nReplaceValue(sWinBMP + D_K_DetectorDescription,
				 Cstring ("WinBMP conversion"));

  // The following may not be correct:

  (void) poHeader->nReplaceValue(sWinBMP + D_K_DetectorVectors,
				 Cstring ("1 0 0 0 1 0"));

  (void) poHeader->nReplaceValue(sWinBMP + D_K_GonioNumValues,
				 (int) 6);
  (void) poHeader->nReplaceValue(sWinBMP + D_K_GonioUnits,
				 Cstring ("deg deg deg mm mm mm"));
  (void) poHeader->nReplaceValue(sWinBMP + D_K_GonioNames,
                 Cstring ("RotZ RotX/Swing RotY TransX TransY TransZ/Dist"));

  (void) poHeader->nReplaceValue(sWinBMP + D_K_GonioVectors,
		 Cstring ("0 0 1 -1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1"));

  a10fTemp[0] = 0.0;
  a10fTemp[1] = 0.0;   // Two theta
  a10fTemp[2] = 0.0;
  a10fTemp[3] = 0.0;
  a10fTemp[4] = 0.0;
  a10fTemp[5] = 100.0; // xtal_to_detector
      
  (void) poHeader->nReplaceValue(sWinBMP + D_K_GonioValues, 6, a10fTemp, 3);

  (void) poHeader->nReplaceValue(sWinBMP + D_K_NonunfType, 
				 Cstring(D_K_NonunfStateNone));

      
  (void) poHeader->nReplaceValue(sWinBMP + D_K_SpatialDistortionType,
				 Cstring (D_K_SpatialTypeSimple));
  a10fTemp[0] = 0.5 * bmpih.biWidth;
  a10fTemp[1] = 0.5 * bmpih.biHeight;
  a10fTemp[2] = 0.0105f;
  a10fTemp[3] = 0.0105f;
  (void) poHeader->nReplaceValue(sWinBMP + D_K_SpatialDistortionInfo,
				 4, a10fTemp, 5);
      
  (void) poHeader->nReplaceValue(sWinBMP + D_K_SpatialDistortionVectors,
				 Cstring("1 0 0 -1"));

  (void) poHeader->nReplaceValue(Cstring (D_K_SaturatedValue), 500000);
  a10fTemp[0] = 1.0;
  a10fTemp[1] = 1.54178f;

  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePrefix) + 
				 D_K_Wavelength, 2, a10fTemp, 6);
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
  a10fTemp[0] = 0.0;   // Start_phi
  a10fTemp[1] = 1.0;   // End_phi
  a10fTemp[2] = a10fTemp[1] - a10fTemp[0];
  a10fTemp[3] = 300.0; // Exposure_time
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

  //Add
  (void) poHeader->nReplaceValue(sWinBMP + "DATA_BYTE", nDataBytes);//int

  (void) poHeader->nReplaceValue(sWinBMP + "HEADER_BYTE", nHeaderBytes);//int

  (void) poHeader->nReplaceValue(sWinBMP + "COLORS", nColors);//int

  if ( (0 != nStat) && ("" != sGetEnv("DTREK_RAXIS_IGNORE")))
    {
      cout << "Warning in WinBMP, error reading header, but ignored!\n";
      nStat = 0;
    }
  return (nStat);
}

int nReadWindowsBMPData(const int *pnFile, unsigned short int *puiData,
		   Cimage_header *poHeader)
{
  int i,j;
  int nStat = 0;
  int nDim0, nDim1;    // Dimensions of the image.
  int nToRead;
  int nDataBytes, nHeaderBytes, nColors;
  unsigned short int *puiTemp;
  unsigned char *pucTemp;

  //cout << "nReadWindowsBMPData called" << endl;
  nStat = poHeader->nGetValue(D_K_Size1, &nDim0);
  nStat = nStat + poHeader->nGetValue(D_K_Size2, &nDim1);

  //cout << "size is " << nDim0 << ' ' << nDim1 << endl;
  //cout << "bytes per pixel is " << nBytesPerPixel << endl;

  if (0 == nStat) 
    nStat = poHeader->nGetValue(Cstring("BMP_DATA_BYTE"), &nDataBytes);
  if (0 == nStat) 
    nStat = poHeader->nGetValue(Cstring("BMP_HEADER_BYTE"), &nHeaderBytes);
  if (0 == nStat) 
    nStat = poHeader->nGetValue(Cstring("BMP_COLORS"), &nColors);

  if (0 == nStat) 
    {
      // Number of bytes to read

      nToRead = nDim0 * nDim1 * nDataBytes + (nDim0%4) * nDim1;
      
      pucTemp = new unsigned char [nToRead];
	
	  RGBQUAD* pRGB = NULL;

	  if (nColors > 256)
	  {
        // Skip the header bytes
	    (void) dskbr((int *)pnFile, (char *)pucTemp, &nHeaderBytes, &nStat);
	  }
	  else
	  {
	    // 54 is the size of the BMP header that we expect to have at the top
	    // of the file
	    nHeaderBytes = 54;
	    (void) dskbr((int *)pnFile, (char *)pucTemp, &nHeaderBytes, &nStat);

		pRGB = new RGBQUAD[nColors];
	    int nRGBTableBytes = 4 * nColors;
        if (0 == nStat)
		{
	      (void) dskbr((int *)pnFile, (char *)pRGB, &nRGBTableBytes, &nStat);
		}
	  }

      if (0 == nStat)
	  {
	    // Maybe analyze header info at this point

	    (void) dskbr((int *)pnFile, (char *)pucTemp, &nToRead, &nStat);
	    if (0 == nStat)
	    {
	      // In the expected format, there are 3 bytes per pixel (RGB)
	      // Convert to an unsigned short int
	      
	      unsigned char *pucTemp2;
	      puiTemp  = puiData;
	      pucTemp2 = pucTemp;
		  int r,g,b;
	      for (j = 0; j < nDim1; j++)
		  {
			  for (i = 0; i < nDim0; i++)
			  {
				  if(3 == nDataBytes)
				  {
				    r = pucTemp2[2];
				    g = pucTemp2[1];
				    b = pucTemp2[0];
				  }
				  else
				  {
				    r = pRGB[pucTemp2[0]].rgbRed;
				    g = pRGB[pucTemp2[0]].rgbGreen;
				    b = pRGB[pucTemp2[0]].rgbBlue;
				  }

				  *puiTemp++ = (unsigned short int)255 - (unsigned short int)
						   (2 + r + g + b) / 3;
				  pucTemp2 += nDataBytes;
			  }
			  if(nDim0 % 4)
			  {
			    pucTemp2 += (nDim0 % 4);
			  }
		  }
	    }
	}
      delete [] pucTemp;

	  if(pRGB)
		  delete [] pRGB;
    }


  if (0 != nStat)
    cout << "ERROR in nReadWindowsBMP!\n" << flush;
  return (nStat);
}

