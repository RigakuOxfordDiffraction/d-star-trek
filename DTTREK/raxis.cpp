//
// Copyright (c) 1995 Molecular Structure Corporation
//
#include <string.h>
#include "dtreksys.h"
#include "raxis.h"
#include "Cimage.h"
#include "Cgoniometer.h"
#include "dskio.h"
#include "dtrekdefs.h"

#ifdef RAXIS_HEADER_BINARY_STRESS_INFO
#include "Cstress.h"
#endif

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
/*
inline float
fSwapFloat(float fFloat)
{
  return ( (fFloat << 24) | ((fFloat << 8) & 0x00ff0000) |
                  ((fFloat >> 8) & 0x0000ff00) | (fFloat >> 24));
}
*/

static void vTrimTrailingSpaces(Cstring& s);
static void vSetHeaderRapidShadowInfo(Cimage_header* poHeader, bool bIsDMaxRapid);

#ifdef RAXIS_HEADER_BINARY_STRESS_INFO
static bool bGetSetStressMeasurementInformation(RAXIS_header& stRAXIS_Header, Cimage_header& poHeader, bool bSwapBytes);
#endif

static const char*   s_pcc_SinglePhiGonio = "SinglePhi Goniometer";
static const char*   s_pcc_QuarterChiGonio = "QuaterChi Goniometer";


// Use this as a buffer for storing RAXIS header.
RAXIS_header g_oRaxisHeader;

int
nReadRAXISheader(const int* pnFile, const Cstring& rsFilename,
                 char *pcBuffer, Cimage_header *poHeader)
{
  // Read an RAXIS header that may or may not have been partially read already.
  // Convert the header to a d*TREK header.
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
  // The RAXIS header is always 1024 bytes long.  It contains alot of binary
  // data that may be VMS native or SGI native.  So we have to figure this out,
  // maybe swap bytes and maybe convert binary floating point stuff.
  // Also the header may be padded.  This routine reads and discards the
  // padding, so that subsequent routines can read the binary pixel data.

  RAXIS_header tRAXIS_header;

  int i;        // Loop counter
  int nStat;    // Local status flag
  int nLen;     // Pixels per line
  int nSwap;    // Flag whether to swap (=1) bytes in integer*4's or not (=0)
  int nRotAxis; // 1,2,3 for omega, chi or phi rot. axes.
  
  bool      bIsRapid = false;// Is this a rapid detector.  If so, we must do some special stuff with the goniometer data.
  bool      bSinglePhiGoniometer = false; // special consideration for the crystal goniometer type.
  
  float *pfTemp;

  Cstring s;
  char ac[256];

  memset(ac,0,sizeof(ac));
  //cout << "sizeof raxis_heder: " << sizeof(tRAXIS_header) << endl;
  if (NULL != pcBuffer)
    {
      memcpy((void*)&tRAXIS_header, (void*)pcBuffer, 512);
      nLen = sizeof(tRAXIS_header) - 512;
      dskbr((int *)pnFile, &tRAXIS_header.a4cSpindle[0], &nLen, &nStat);
    }
  else
    {
      nLen = sizeof(tRAXIS_header);
      dskbr((int *)pnFile, (char *)&tRAXIS_header, &nLen, &nStat);
    }


  nLen  = tRAXIS_header.a2nPix_num[0];
//  cout << "pixnum[0] is " << nLen << endl;

  nSwap = 0;                          // Not swapped;
  nStat = 1;
  if ((nLen < 1) || (nLen > 10000) )
    {
      // Probably need to swap, so try it

       nSwap =   nSwapLong(nLen);

       if ((nSwap < 1) || (nSwap > 10000) )
         {
           // Error
         }
       else
         {
           nLen  = nSwap;
           nSwap = 1;   // Swapping need!
           nStat = 0;
         }
     }
  else
    {
      nStat = 0;
    }


//  cout << "pixnum[0] now is " << nLen << endl;
  // Here we know if need to swap, but don't know if VAX or IEEE float

  if (0 == nStat)
    {
      if (1 == nSwap)
        {
          tRAXIS_header.nOsc        = nSwapLong(tRAXIS_header.nOsc);
          //for (i = 0; i < 51; i++)
            //{
              //tRAXIS_header.a51nReserve3[i] = nSwapLong(tRAXIS_header.a51nReserve3[i]);
            //}
          tRAXIS_header.a2nPix_num[0] = nLen;
          tRAXIS_header.a2nPix_num[1] = nSwapLong(tRAXIS_header.a2nPix_num[1]);
          tRAXIS_header.a2nRecord[0]  = nSwapLong(tRAXIS_header.a2nRecord[0]);
          tRAXIS_header.a2nRecord[1]  = nSwapLong(tRAXIS_header.a2nRecord[1]);
          tRAXIS_header.nRead_start   = nSwapLong(tRAXIS_header.nRead_start);
          tRAXIS_header.nIP_num       = nSwapLong(tRAXIS_header.nIP_num);
          tRAXIS_header.a3nDrxz[0]    = nSwapLong(tRAXIS_header.a3nDrxz[0]);
          tRAXIS_header.a3nDrxz[1]    = nSwapLong(tRAXIS_header.a3nDrxz[1]);
          tRAXIS_header.a3nDrxz[2]    = nSwapLong(tRAXIS_header.a3nDrxz[2]);

// WARNING!  Linux may not necessarily be running on an IntelPC
#if defined(WIN32) || defined(LINUX) || defined(OSF1)
// WARNING!  Linux may not necessarily be running on an IntelPC
/*
if we convert with VAXtoIEEE no need to swapfloat when going from VAX to IRIX,
but what about VAX to PC? 
2005-09-01 Yes, we need to swapfloat when going from VAX to PC even though
           integers apparently do not need to be swapped
*/
         for (i = 0; i < 6; i++)
            {
              tRAXIS_header.a6fCell[i] = fSwapFloat(tRAXIS_header.a6fCell[i]);
            }

	 //        cout << "before swap fWave: " << tRAXIS_header.fWave << endl;
          tRAXIS_header.fMosaic       = fSwapFloat(tRAXIS_header.fMosaic);
          tRAXIS_header.fWave         = fSwapFloat(tRAXIS_header.fWave);
	  //        cout << "after swap fWave: " << tRAXIS_header.fWave << endl;
          tRAXIS_header.fMono_2       = fSwapFloat(tRAXIS_header.fMono_2);
          tRAXIS_header.fCamera       = fSwapFloat(tRAXIS_header.fCamera);
          tRAXIS_header.fKv           = fSwapFloat(tRAXIS_header.fKv);
          tRAXIS_header.fMa           = fSwapFloat(tRAXIS_header.fMa);
          tRAXIS_header.fWeissenberg  = fSwapFloat(tRAXIS_header.fWeissenberg);
          tRAXIS_header.a3fPhi[0]     = fSwapFloat(tRAXIS_header.a3fPhi[0]);
          tRAXIS_header.a3fPhi[1]     = fSwapFloat(tRAXIS_header.a3fPhi[1]);
          tRAXIS_header.a3fPhi[2]     = fSwapFloat(tRAXIS_header.a3fPhi[2]);
          tRAXIS_header.fEx_time      = fSwapFloat(tRAXIS_header.fEx_time);
          tRAXIS_header.a2fXray[0]    = fSwapFloat(tRAXIS_header.a2fXray[0]);
          tRAXIS_header.a2fXray[1]    = fSwapFloat(tRAXIS_header.a2fXray[1]);
          tRAXIS_header.a3fCircle[0]  = fSwapFloat(tRAXIS_header.a3fCircle[0]);
          tRAXIS_header.a3fCircle[1]  = fSwapFloat(tRAXIS_header.a3fCircle[1]);
          tRAXIS_header.a3fCircle[2]  = fSwapFloat(tRAXIS_header.a3fCircle[2]);
          tRAXIS_header.fMu           = fSwapFloat(tRAXIS_header.fMu);
          tRAXIS_header.a2fPix_size[0]= fSwapFloat(tRAXIS_header.a2fPix_size[0]);
          tRAXIS_header.a2fPix_size[1]= fSwapFloat(tRAXIS_header.a2fPix_size[1]);
          tRAXIS_header.fRatio        = fSwapFloat(tRAXIS_header.fRatio);
          tRAXIS_header.a2fFading[0]  = fSwapFloat(tRAXIS_header.a2fFading[0]);
          tRAXIS_header.a2fFading[1]  = fSwapFloat(tRAXIS_header.a2fFading[1]);

	  tRAXIS_header.nMagicNum     = nSwapLong(tRAXIS_header.nMagicNum);
	  tRAXIS_header.nNumGonAxes   = nSwapLong(tRAXIS_header.nNumGonAxes);
	  tRAXIS_header.nScanAxisNum  = nSwapLong(tRAXIS_header.nScanAxisNum);

	  pfTemp = &tRAXIS_header.a5x3fGonVecs[0][0];
	  for (i = 0; i < 30; i++)
	    {
	      *pfTemp = fSwapFloat(*pfTemp);
	      pfTemp++;
	    }
#endif
        }

      // Next check if floats need to be converted by looking
      // at the pix_size.
      // Assume reasonable values are >= 0.01 and <= 1.0 mm

      float fTemp1, fTemp2, fTemp3;
      fTemp1 = tRAXIS_header.a2fPix_size[0];

      //      cout << "No swap ftemp1 before: " << fTemp1 << endl;
      fTemp2 = fIEEEtoVAX(&fTemp1);  // Careful about the order here
      fTemp3 = fVAXtoIEEE(&fTemp1);
      //  cout << "after: fTemp2, fTemp3: " << fTemp2 << ' ' << fTemp3 << endl;

      if ( (fTemp1 < 0.04) || (fTemp1 > 1.0) )
        {
          // Looks bogus, try to convert 2 ways:

          fTemp3 = fSwapFloat(fTemp1);
          fTemp2 = fIEEEtoVAX(&fTemp1);  // Careful about the order here
          fTemp1 = fVAXtoIEEE(&fTemp1);
	  //	  cout << "LATE: fTemp2: " << fTemp2 << endl;
          if ( (fTemp2 >= 0.04) && (fTemp2 <= 1.0) )
            {
              // We need to switch all floats in the header from IEEE to VAX
              for (i = 0; i < 6; i++)
                {
                  tRAXIS_header.a6fCell[i] = fIEEEtoVAX(&tRAXIS_header.a6fCell[i]);
                }
              tRAXIS_header.fMosaic       = fIEEEtoVAX(&tRAXIS_header.fMosaic);
              tRAXIS_header.fWave         = fIEEEtoVAX(&tRAXIS_header.fWave);
//            cout << "after IEEEtoVAX fWave: " << tRAXIS_header.fWave << endl;
              tRAXIS_header.fMono_2       = fIEEEtoVAX(&tRAXIS_header.fMono_2);
              tRAXIS_header.fKv           = fIEEEtoVAX(&tRAXIS_header.fKv);
              tRAXIS_header.fMa           = fIEEEtoVAX(&tRAXIS_header.fMa);
              tRAXIS_header.a3fPhi[0]     = fIEEEtoVAX(&tRAXIS_header.a3fPhi[0]);
              tRAXIS_header.a3fPhi[1]     = fIEEEtoVAX(&tRAXIS_header.a3fPhi[1]);
              tRAXIS_header.a3fPhi[2]     = fIEEEtoVAX(&tRAXIS_header.a3fPhi[2]);
              tRAXIS_header.fEx_time      = fIEEEtoVAX(&tRAXIS_header.fEx_time);
              tRAXIS_header.a2fXray[0]    = fIEEEtoVAX(&tRAXIS_header.a2fXray[0]);
              tRAXIS_header.a2fXray[1]    = fIEEEtoVAX(&tRAXIS_header.a2fXray[1]);
              tRAXIS_header.a3fCircle[0]   = fIEEEtoVAX(&tRAXIS_header.a3fCircle[0]);
              tRAXIS_header.a3fCircle[1]   = fIEEEtoVAX(&tRAXIS_header.a3fCircle[1]);
              tRAXIS_header.a3fCircle[2]   = fIEEEtoVAX(&tRAXIS_header.a3fCircle[2]);
              tRAXIS_header.fMu           = fIEEEtoVAX(&tRAXIS_header.fMu);
              tRAXIS_header.a2fPix_size[0] = fIEEEtoVAX(&tRAXIS_header.a2fPix_size[0]);
              tRAXIS_header.a2fPix_size[1] = fIEEEtoVAX(&tRAXIS_header.a2fPix_size[1]);
              tRAXIS_header.fRatio        = fIEEEtoVAX(&tRAXIS_header.fRatio);
              tRAXIS_header.a2fFading[0]   = fIEEEtoVAX(&tRAXIS_header.a2fFading[0]);
              tRAXIS_header.a2fFading[1]   = fIEEEtoVAX(&tRAXIS_header.a2fFading[1]);
	      pfTemp = &tRAXIS_header.a5x3fGonVecs[0][0];
	      for (i = 0; i < 30; i++)
		{
		  *pfTemp = fIEEEtoVAX(pfTemp);
		  pfTemp++;
		}
            }
          else if ( (fTemp1 >= 0.04) && (fTemp1 <= 1.00) )
            {
	      //	      cout << "LATE: fTemp1: " << fTemp1 << endl;
              // We need to switch all floats in the header from VAX to IEEE
              for (i = 0; i < 6; i++)
                {
                  tRAXIS_header.a6fCell[i] = fVAXtoIEEE(&tRAXIS_header.a6fCell[i]);
                }
              tRAXIS_header.fMosaic       = fVAXtoIEEE(&tRAXIS_header.fMosaic);
              tRAXIS_header.fWave         = fVAXtoIEEE(&tRAXIS_header.fWave);
//            cout << "after VAXtoIEEE fWave: " << tRAXIS_header.fWave << endl;
              tRAXIS_header.fMono_2       = fVAXtoIEEE(&tRAXIS_header.fMono_2);
              tRAXIS_header.fCamera       = fVAXtoIEEE(&tRAXIS_header.fCamera);
              tRAXIS_header.fKv           = fVAXtoIEEE(&tRAXIS_header.fKv);
              tRAXIS_header.fMa           = fVAXtoIEEE(&tRAXIS_header.fMa);
              tRAXIS_header.a3fPhi[0]     = fVAXtoIEEE(&tRAXIS_header.a3fPhi[0]);
              tRAXIS_header.a3fPhi[1]     = fVAXtoIEEE(&tRAXIS_header.a3fPhi[1]);
              tRAXIS_header.a3fPhi[2]     = fVAXtoIEEE(&tRAXIS_header.a3fPhi[2]);
              tRAXIS_header.fEx_time      = fVAXtoIEEE(&tRAXIS_header.fEx_time);
              tRAXIS_header.a2fXray[0]    = fVAXtoIEEE(&tRAXIS_header.a2fXray[0]);
              tRAXIS_header.a2fXray[1]    = fVAXtoIEEE(&tRAXIS_header.a2fXray[1]);
              tRAXIS_header.a3fCircle[0]  = fVAXtoIEEE(&tRAXIS_header.a3fCircle[0]);
              tRAXIS_header.a3fCircle[1]  = fVAXtoIEEE(&tRAXIS_header.a3fCircle[1]);
              tRAXIS_header.a3fCircle[2]  = fVAXtoIEEE(&tRAXIS_header.a3fCircle[2]);
              tRAXIS_header.fMu           = fVAXtoIEEE(&tRAXIS_header.fMu);
              tRAXIS_header.a2fPix_size[0] = fVAXtoIEEE(&tRAXIS_header.a2fPix_size[0]);
              tRAXIS_header.a2fPix_size[1] = fVAXtoIEEE(&tRAXIS_header.a2fPix_size[1]);
              tRAXIS_header.fRatio        = fVAXtoIEEE(&tRAXIS_header.fRatio);
              tRAXIS_header.a2fFading[0]   = fVAXtoIEEE(&tRAXIS_header.a2fFading[0]);
              tRAXIS_header.a2fFading[1]   = fVAXtoIEEE(&tRAXIS_header.a2fFading[1]);
	      pfTemp = &tRAXIS_header.a5x3fGonVecs[0][0];
	      for (i = 0; i < 30; i++)
		{
		  *pfTemp = fVAXtoIEEE(pfTemp);
		  pfTemp++;
		}
            }
          else if ( (fTemp3 >= 0.04) && (fTemp3 <= 1.0) )
            {
	      //	      cout << "LATE: fTemp3: " << fTemp3 << endl;
              // We need to swap all floats in the header
              for (i = 0; i < 6; i++)
                {
                  tRAXIS_header.a6fCell[i] = fSwapFloat(tRAXIS_header.a6fCell[i]);
                }
              tRAXIS_header.fMosaic       = fSwapFloat(tRAXIS_header.fMosaic);
              tRAXIS_header.fWave         = fSwapFloat(tRAXIS_header.fWave);
              tRAXIS_header.fMono_2       = fSwapFloat(tRAXIS_header.fMono_2);
              tRAXIS_header.fCamera       = fSwapFloat(tRAXIS_header.fCamera);
              tRAXIS_header.fKv           = fSwapFloat(tRAXIS_header.fKv);
              tRAXIS_header.fMa           = fSwapFloat(tRAXIS_header.fMa);
              tRAXIS_header.a3fPhi[0]     = fSwapFloat(tRAXIS_header.a3fPhi[0]);
              tRAXIS_header.a3fPhi[1]     = fSwapFloat(tRAXIS_header.a3fPhi[1]);
              tRAXIS_header.a3fPhi[2]     = fSwapFloat(tRAXIS_header.a3fPhi[2]);
              tRAXIS_header.fEx_time      = fSwapFloat(tRAXIS_header.fEx_time);
              tRAXIS_header.a2fXray[0]    = fSwapFloat(tRAXIS_header.a2fXray[0]);
              tRAXIS_header.a2fXray[1]    = fSwapFloat(tRAXIS_header.a2fXray[1]);
              tRAXIS_header.a3fCircle[0]  = fSwapFloat(tRAXIS_header.a3fCircle[0]);
              tRAXIS_header.a3fCircle[1]  = fSwapFloat(tRAXIS_header.a3fCircle[1]);
              tRAXIS_header.a3fCircle[2]  = fSwapFloat(tRAXIS_header.a3fCircle[2]);
              tRAXIS_header.fMu           = fSwapFloat(tRAXIS_header.fMu);
              tRAXIS_header.a2fPix_size[0] = fSwapFloat(tRAXIS_header.a2fPix_size[0]);
              tRAXIS_header.a2fPix_size[1] = fSwapFloat(tRAXIS_header.a2fPix_size[1]);
              tRAXIS_header.fRatio        = fSwapFloat(tRAXIS_header.fRatio);
              tRAXIS_header.a2fFading[0]   = fSwapFloat(tRAXIS_header.a2fFading[0]);
              tRAXIS_header.a2fFading[1]   = fSwapFloat(tRAXIS_header.a2fFading[1]);

	      tRAXIS_header.nMagicNum     = nSwapLong(tRAXIS_header.nMagicNum);
	      tRAXIS_header.nNumGonAxes   = nSwapLong(tRAXIS_header.nNumGonAxes);
	      tRAXIS_header.nScanAxisNum  = nSwapLong(tRAXIS_header.nScanAxisNum);

	      pfTemp = &tRAXIS_header.a5x3fGonVecs[0][0];
	      for (i = 0; i < 30; i++)
		{
		  *pfTemp = fSwapFloat(*pfTemp);
		  pfTemp++;
		}
            }
          else
            {
              nStat = 2;
            }
        }

//      cout << "after head stuff fWave: " << tRAXIS_header.fWave << endl;
      // Now construct a d*TREK header from the R-AXIS info
      poHeader->nReplaceValue(D_K_OriginalImageFormat, Cimage_header::sGetFormatTypeName(enImage_format_RAXIS));

      float fTemp[10];

      strncpy(ac, tRAXIS_header.a12cDate, 12);
      s = ac;
      memset(ac, 0, 12);
      vTrimTrailingSpaces(s);
      (void)poHeader->nReplaceValue( Cstring (D_K_RxPrefix)
                                    + Cstring ("CREATE_DATE"), s);
      
      strncpy(ac, tRAXIS_header.a10cDevice, 10);
      s = ac;
      memset(ac, 0, 10);
      vTrimTrailingSpaces(s);
      bIsRapid = (bool) (s.contains("R-AXIS-CS"));

      (void)poHeader->nReplaceValue(Cstring (D_K_DetectorType), s);
      
      if (bIsRapid)
        {
          // Set some oblique incidence correction factors if Chromium wavelength

          if ( (2.28 < tRAXIS_header.fWave) && (2.30 > tRAXIS_header.fWave) )
            {
              Cstring sOblique;
              double dAirMuT;

              // mu * t for (all from JDF)
              // 1. 2.5 cm of air at front before helium cone [31.8227 * 0.0012269 * 2.5 cm]
              // 2. polycarbonate window [0.07258]
              // 3. Helium cone of about 120 mm radius [0.498 * 1.66321 x 10^-4 * 10.2]
              // 4. Mylar window [0.03777]
              // 5. Air gap between helium cone and IP phosphor
          
              // Air gap 5 is the excess distance inside the detector housing 
              // past 127.4 mm.  The 0.1 in the next line is 0.1 cm per mm to 
              // convert the t to cm since mu is in units of per cm.

              sOblique = "5 0 0.08938 0 0.07258 0 0.0008448 0 0.03777 ";
	      dAirMuT = ((double)tRAXIS_header.fCamera - 127.4);
	      if (0.0 < dAirMuT)
		{
		  dAirMuT = 31.8227 * 0.00112269
		    * 0.1 * dAirMuT;
		  sOblique = sOblique + Cstring("0 ") + Cstring(dAirMuT);
		}
	      (void) poHeader->nReplaceValue(Cstring (D_K_IntegrateOblique),
                                             sOblique);
            }
        }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      bool      bIsDMaxRapid = false;
      // In case the header does not say it IS dmax rapid, the user can define an environment variable
      if( "" != sGetEnv("DTREK_IMAGE_DMAX_RAPID") )  
      {
        bIsDMaxRapid = true;
      }
      else
      {
        strncpy(ac, tRAXIS_header.a80cMemo, 15);
        s = ac;
        memset(ac, 0, 15);
        vTrimTrailingSpaces(s);
        bIsDMaxRapid = (bool)(s.contains("RINT RAPID"));
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////

      (void)poHeader->nReplaceValue(Cstring (D_K_DtdisplayOrientation), 
                                    Cstring("-X-Y"));

      (void) poHeader->nReplaceValue(Cstring (D_K_Type), Cstring("mad"));
      (void) poHeader->nReplaceValue(Cstring (D_K_Dim), Cstring("2"));
      (void) poHeader->nReplaceValue(Cstring (D_K_Size1), tRAXIS_header.a2nPix_num[0]);
      (void) poHeader->nReplaceValue(Cstring (D_K_Size2), tRAXIS_header.a2nPix_num[1]);
      (void) poHeader->nReplaceValue(Cstring (D_K_DataType),
                                   Cstring("unsigned short int"));
      (void) poHeader->nReplaceValue(Cstring (D_K_Compression),
                                     Cstring(D_K_Raxis));
      (void) poHeader->nReplaceValue(Cstring (D_K_RaxisCompressionRatio),
                                     tRAXIS_header.fRatio);
      // Set expected saturated value
      //+jwp 21-Nov-2002      nLen = 1048544;
      nLen = 1048512;
      //-jwp 21-Nov-2002
      //      cout << "INFO fRatio: "
      //<< tRAXIS_header.fRatio << '\n';
      if (8.0 == tRAXIS_header.fRatio)
        nLen = 262136;
      else if (32.0 != tRAXIS_header.fRatio)
        cout << "WARNING, assuming R-AXIS IV image, fRatio: "
             << tRAXIS_header.fRatio << '\n';

      (void) poHeader->nReplaceValue(Cstring (D_K_SaturatedValue), nLen);

      (void) poHeader->nReplaceValue(Cstring (D_K_RaxisRecordSize),
                                     tRAXIS_header.a2nRecord[0]);
      (void) poHeader->nReplaceValue(Cstring (D_K_RaxisReadStart),
                                     tRAXIS_header.nRead_start);
      (void) poHeader->nReplaceValue(Cstring (D_K_RaxisReadLines),
                                     tRAXIS_header.a2nRecord[1]);
      (void) poHeader->nReplaceValue(Cstring ("RAXIS_DETNUM"),
                                     tRAXIS_header.nIP_num);
      int nByte_order = nGetByteOrderCPU();

/*
 * Images collected on a PC or on an SGI always have the image data written
 * out as an array of chars, and so the image data is always big-endian since
 * that is the way that it comes from the R-AXIS through the SCSI bus.
 * Note that this may be a different endian from how the header is written
 * out.
 */
      strncpy(ac,tRAXIS_header.a10cCpu,10);
      s = ac;
      memset(ac,0,10);
      vTrimTrailingSpaces(s);
      if(s.contains("PC") || s.contains("IRIS")){
         poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
                                 Cstring (D_K_BigEndian));
      }
      else if (1 == nSwap)
        {
          // The orders are different
          if (nByte_order == eCPU_big_endian)
            {
              poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
                                      Cstring (D_K_LittleEndian));
            }
          else if (nByte_order == eCPU_little_endian)
            {
              poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
                                      Cstring (D_K_BigEndian));
            }
        }
      else
        {
          // The orders are same
          if (nByte_order == eCPU_big_endian)
            {
              poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
                                      Cstring (D_K_BigEndian));
            }
          else if (nByte_order == eCPU_little_endian)
            {
              poHeader->nReplaceValue(Cstring (D_K_ByteOrder),
                                      Cstring (D_K_LittleEndian));
            }
        }

      (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNumber), (int)1);
      (void) poHeader->nReplaceValue(Cstring (D_K_DetectorNames),
                                     Cstring (D_K_RxPrefix));
      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) +
                                              D_K_DetectorDimensions, 2,
                                     tRAXIS_header.a2nPix_num);

      // Special kludge to fix R-AXIS2 pixel size info

      if (tRAXIS_header.a2fPix_size[0] == tRAXIS_header.a2fPix_size[1])
        {
          if (   (tRAXIS_header.a2fPix_size[1] >= 0.104)
              && (tRAXIS_header.a2fPix_size[1] <= 0.106) )
            {
              tRAXIS_header.a2fPix_size[0] = 0.1017f;
            }
        }

      fTemp[0] = (float) tRAXIS_header.a2nPix_num[0]
                       * tRAXIS_header.a2fPix_size[0];
      fTemp[1] = (float) tRAXIS_header.a2nPix_num[1]
                       * tRAXIS_header.a2fPix_size[1];
      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) + D_K_DetectorSize,
                                     2, fTemp, 5);
      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) +
                                              D_K_DetectorDescription,
                                     Cstring ("RAXIS conversion"));
// The following may not be correct:
      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) + D_K_DetectorVectors,
                                     Cstring ("1 0 0 0 1 0"));
      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) +
                                              D_K_DetectorRefineFlags,
                                     Cstring ("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"));

      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) + D_K_GonioNumValues,
                                     (int) 6);
      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) + D_K_GonioUnits,
                                     Cstring ("deg deg deg mm mm mm"));
      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) + D_K_GonioNames,
                                     Cstring ("RotZ RotX/Swing RotY TransX TransY TransZ/Dist"));

      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) + D_K_GonioVectors,
                                     Cstring ("0 0 1 -1 0 0 0 1 0 1 0 0 0 1 0 0 0 -1"));

      // Optics information from header

      memset(ac, 0, sizeof(ac));
      strncpy(ac, tRAXIS_header.a20cMonochro, 20);
      s = ac;
      memset(ac, 0, 20);
      vTrimTrailingSpaces(s);
      s.downcase();
      if (s.contains("none") || (0.0 == tRAXIS_header.fMono_2) )
        {
          strncpy(ac,tRAXIS_header.a80cXraymemo, 80);
          s = ac;
          memset(ac, 0, 80);
	  vTrimTrailingSpaces(s);
          poHeader->nReplaceValue(Cstring(D_K_OpticsType), s);
          poHeader->nReplaceValue(Cstring(D_K_OpticsAngle), 0.0);
        }
      else
        {
          poHeader->nReplaceValue(Cstring(D_K_OpticsType), s);
          poHeader->nReplaceValue(Cstring(D_K_OpticsAngle), tRAXIS_header.fMono_2);
      }
      strncpy(ac, tRAXIS_header.a20cCollimator, 20);
      s = ac;
      memset(ac, 0, 20);
      vTrimTrailingSpaces(s);
      poHeader->nReplaceValue(Cstring(D_K_OpticsCollimator), s);
      strncpy(ac, tRAXIS_header.a4cFilter, 4);
      s = ac;
      memset(ac, 0, 4);
      vTrimTrailingSpaces(s);
      poHeader->nReplaceValue(Cstring(D_K_OpticsFilter), s);

      fTemp[0] = 0.0;
      if (bIsRapid)
          fTemp[1] = 0.0;
      else
          fTemp[1] = tRAXIS_header.a3fCircle[2];  // 2theta
      fTemp[2] = 0.0;
      fTemp[3] = 0.0;
      fTemp[4] = 0.0;
      fTemp[5] = tRAXIS_header.fCamera;

      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) + D_K_GonioValues,
                                     6, fTemp, 3);

      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) + D_K_NonunfType,
                                     Cstring(D_K_NonunfStateNone));


      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) +
                                              D_K_SpatialDistortionType,
                                     Cstring (D_K_SpatialTypeSimple));
      fTemp[0] = tRAXIS_header.a2fXray[0];
      if (0.0 >= fTemp[0])
        fTemp[0] = (float)(tRAXIS_header.a2nPix_num[0] / 2);
      fTemp[1] = tRAXIS_header.a2fXray[1];
      if (0.0 >= fTemp[1])
        fTemp[1] = (float)(tRAXIS_header.a2nPix_num[1] / 2);
      fTemp[2] = tRAXIS_header.a2fPix_size[0];
      fTemp[3] = tRAXIS_header.a2fPix_size[1];
      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) +
                                              D_K_SpatialDistortionInfo,
                                     4, fTemp, 4);

      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) +
                                              D_K_SpatialBeamPosition,
                                     2,fTemp,4);

      (void) poHeader->nReplaceValue(Cstring (D_K_RxPrefix) +
                                              D_K_SpatialDistortionVectors,
                                   Cstring("0 1 -1 0"));
      fTemp[0] = 1.0;
      fTemp[1] = tRAXIS_header.fWave;
      if (fTemp[1] <= 0.0) fTemp[1] = 1.54178f;
      (void) poHeader->nReplaceValue(Cstring (D_K_SourcePrefix) +
                                     D_K_Wavelength, 2, fTemp, 6);
      fTemp[2] = ABS(fTemp[1] - 1.542);
      fTemp[3] = ABS(fTemp[1] - 0.717);
      fTemp[4] = ABS(fTemp[1] - 2.29);
  
      (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
                                     Cstring ("0.5 1 0 0"));
      bool bMonochrom = FALSE;
      s = "";
      poHeader->nGetValue(Cstring(D_K_OpticsType), &s);
      s.downcase();
      bMonochrom = s.contains("graphite");

      // Start with default of no polarization

     (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
                                     Cstring ("0.5 1 0 0"));
      if (0.01 > fTemp[3])
    	{
	  // Probably MoKalpha
	  if (bMonochrom)
	    (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
					   Cstring ("0.4886 1 0 0"));
	}
      else if (0.01 > fTemp[2])
	{
	  if (bMonochrom)
	    (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
					   Cstring ("0.4442 1 0 0"));
	}
      else if (0.01 > fTemp[4])
	{
	  if (bMonochrom)
	    (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
					   Cstring ("0.3698 1 0 0"));
	}
      else
        {
	  // Not CuKalpha, not MoKa, not CrKa probably synchrotron
	  // but vector is still a problem

	  (void) poHeader->nReplaceValue(Cstring (D_K_SourcePolarz),
					 Cstring ("0.95 1 0 0"));
	  // Note: Could also be "0.95 0 1 0"
	}
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
      (void) poHeader->nReplaceValue(Cstring (D_K_SourceRefineFlags),
                                     Cstring ("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"));
      poHeader->nReplaceValue(Cstring (D_K_SourceVoltage),tRAXIS_header.fKv);
      poHeader->nReplaceValue(Cstring (D_K_SourceAmperage),tRAXIS_header.fMa);
      memset(ac,0,sizeof(ac));
      strncpy(ac,tRAXIS_header.a12cFocus,12);
      s = ac;
      memset(ac,0,12);
      vTrimTrailingSpaces(s);
      poHeader->nReplaceValue(Cstring(D_K_SourceFocus),s);

      if (RAXIS_MAGIC_NUM == tRAXIS_header.nMagicNum)
	{
	  // Magic number is 1, so have
	  // NEW binary format and definition of crystal goniometer and scan vectors

	  //cout << "Magic num found in RAXIS binary header!\n" << flush;

	  // First deal with the CRYSTAL_GONIO_* items

	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) +
					 D_K_GonioNumValues,
					 tRAXIS_header.nNumGonAxes);
	  strncpy(ac, tRAXIS_header.a40cAxesNames, 40);
	  ac[41] = '\0';
	  s = ac;
	  vTrimTrailingSpaces(s);
	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioNames,
					 s);
	  s = "";
	  for (i = 0; i < tRAXIS_header.nNumGonAxes; i++)
	    s += "deg ";
	  vTrimTrailingSpaces(s);
	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioUnits,
					 s);
	  s  = "";
	  pfTemp = &tRAXIS_header.a5x3fGonVecs[0][0];
	  for (i = 0; i < tRAXIS_header.nNumGonAxes*3; i++, pfTemp++)
	    s += Cstring(*pfTemp) + ' ';
	  vTrimTrailingSpaces(s);
	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
					 s);

	  if (tRAXIS_header.nNumGonAxes > 10)
	    {
	      cout << "WARNING in raxis(), too many goniometer axes: " 
		   << tRAXIS_header.nNumGonAxes << ", reset to 10.\n" << endl;
	      tRAXIS_header.nNumGonAxes = 10;
	      tRAXIS_header.nMagicNum = 0;
	      goto MAGIC_ERROR;
	    }
	  for (i = 0; i < tRAXIS_header.nNumGonAxes; i++)
	    {
	      fTemp[i] = tRAXIS_header.a5fGonStart[i];
	    }
	  if (10 <= tRAXIS_header.nScanAxisNum)
	    {
	      cout << "WARNING in raxis(), bogus Scan Axis Num in binary header: " 
		   << tRAXIS_header.nScanAxisNum << ".\n" << endl;
	      tRAXIS_header.nMagicNum = 0;
	      goto MAGIC_ERROR;
	    }
	  fTemp[tRAXIS_header.nScanAxisNum] = 0.0;  // This might be the offset
	  for (i = 0; i < tRAXIS_header.nNumGonAxes; i++)
	    {
	      if (0.0 != tRAXIS_header.a5fGonOffset[i])
		{
		  cout << "WARNING in raxis(), unexpected non-zero offset\n"
                       << "        for axis " << i << " (range 0-4).  Offset is "
                       << tRAXIS_header.a5fGonOffset[i] << endl;
		}
	    }
	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
					 tRAXIS_header.nNumGonAxes, fTemp, 3);
	  
      ////////////////////////////////////////////////////////////////////////////
      // RB 9/27/05 Get additional info about the crystal goniometer type
      strncpy(ac, tRAXIS_header.a80cMemo, 25);
      s = ac;
      memset(ac, 0, 25);
      vTrimTrailingSpaces(s);
      
      if( s.contains(s_pcc_SinglePhiGonio) )
      {
          poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) +
					              D_K_GonioDescription,
					              Cstring(s_pcc_SinglePhiGonio));
          
          bSinglePhiGoniometer = true;

      }
      else if( s.contains(s_pcc_QuarterChiGonio) )
      {
          poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) +
					              D_K_GonioDescription,
					              Cstring(s_pcc_QuarterChiGonio));
      }
	  else  
      {
          // RB 9/27/05 If no phi/chi goniometer description is given. keep the header keyword general.
		  poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) +
					              D_K_GonioDescription,
					              Cstring ("Crystal Goniometer"));
      }
      ///////////////////////////////////////////////////////////////////////////////
      
      // Next deal with the SCAN_* items
	  fTemp[0] = tRAXIS_header.a5fGonStart[tRAXIS_header.nScanAxisNum];
	  fTemp[1] = tRAXIS_header.a5fGonEnd[tRAXIS_header.nScanAxisNum];
	  fTemp[2] = fTemp[1] - fTemp[0];
	  fTemp[3] = tRAXIS_header.fEx_time * 60.0f;
	  fTemp[4] = (float) tRAXIS_header.nOsc;
	  fTemp[5] = 0.0;
	  fTemp[6] = 0.0;
	  fTemp[7] = 0.0;
	  fTemp[8] = 0.0;
	  fTemp[9] = 0.0;
	  (void) poHeader->nReplaceValue(Cstring (D_K_Rotation),
					 10, fTemp, 3);
	  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_Rotation,
					 10, fTemp, 3);
	  Cimage_header *poTempHeader;
	  poTempHeader = new Cimage_header(poHeader->sGet().string());
	  Cgoniometer *poGoniometer;
	  poGoniometer = new Cgoniometer(*poTempHeader, D_K_CrystalPrefix);
	  delete poTempHeader;
	  if (poGoniometer->bIsAvailable())
	    {
	      poGoniometer->nGetRotVector(tRAXIS_header.nScanAxisNum, (float *)&fTemp[0]);
	      s = poGoniometer->sGetName(tRAXIS_header.nScanAxisNum);
	      (void) poHeader->nReplaceValue(Cstring (D_K_RotAxisName),
					     s);
	      (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) +
					     D_K_RotAxisName,
					     s);
	    }
	  else
	    {
	      // The X vector which is the default rotation axis
	      
	      fTemp[0] = 1.0f;
	      fTemp[1] = 0.0f;
	      fTemp[2] = 0.0f;
	      s        = "omega";
	      (void) poHeader->nReplaceValue(Cstring (D_K_RotAxisName),
					     s);
	      (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) +
					     D_K_RotAxisName,
					     s);
	      cout << "WARNING crystal goniometer vector problem!\n";
	    }
	  delete poGoniometer;
	  (void) poHeader->nReplaceValue(Cstring (D_K_RotVector),
					 3, &fTemp[0], 4);
	  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_RotVector,
					 3, &fTemp[0], 4);
	}

    if( bIsRapid && !bSinglePhiGoniometer )
        vSetHeaderRapidShadowInfo(poHeader, bIsDMaxRapid);

      // Note: the next block is executed after the above block,
      //       since if there is an error when the magic number is 1,
      //       then we can ignore it and enter this block to try to recover

MAGIC_ERROR:
    if (RAXIS_MAGIC_NUM != tRAXIS_header.nMagicNum)
    {
	  //cout << "No magic num in RAXIS binary header!\n" << flush;

	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) +
					 D_K_GonioNumValues,
					 (int) 3);
	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioNames,
					 Cstring ("Omega Chi Phi"));
	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioUnits,
					 Cstring ("deg deg deg"));
      
        ///////////////////////////////////////////////////////////////////////////////////////// 
        Cstring         strGonioVectors("");
        if( bIsRapid )
        {          
            if( bIsDMaxRapid )
                strGonioVectors = "-1 0 0 0 0 1 1 0 0";    
            else
	      {
                strGonioVectors = "-1 0 0 0 0 -1 -1 0 0";
	      }
        } 
        else 
            strGonioVectors = "1 0 0 0 0 1 1 0 0";

        poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors, strGonioVectors);
	    /////////////////////////////////////////////////////////////////////////////////////////

      fTemp[0] = tRAXIS_header.a3fCircle[0];  // Omega
	  fTemp[1] = tRAXIS_header.a3fCircle[1];  // Chi
	  fTemp[2] = tRAXIS_header.a3fPhi[0];     // Phi datum
	  
      //////////////////////////////////////////////////////////////////////////////////////////////////
      // RB 6/9/04 Need to enforce fTemp[2] = 0.0 for the following reason. CrystalClear server writes a3fPhi[0] 
      // as zero always. d*TREK data processing then adds "start" and "end" rotation axis positions to that zero value. 
      // But RINT/XRD writes a3fPhi[0] as non-zero, which breaks the d*TREK data processing.
      if( bIsDMaxRapid )
        fTemp[2] = 0.0;
///////////////////////////////////////////////////////////////////////////////

      if (bIsRapid) 
	    {
	      std::swap(fTemp[0],fTemp[2]);
	      nRotAxis = 1;
	    } 
	  else
	    nRotAxis = 3;
	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
					 3, fTemp, 3);
	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) +
					 D_K_GonioDescription,
					 Cstring ("Eulerian 3-circle"));
	  (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + 
					 D_K_CrystalXUnitCell,
					 6, tRAXIS_header.a6fCell, 3);
	  fTemp[0] = tRAXIS_header.a3fPhi[1];
	  fTemp[1] = tRAXIS_header.a3fPhi[2];
	  fTemp[2] = fTemp[1] - fTemp[0];
	  fTemp[3] = tRAXIS_header.fEx_time * 60.0f;
	  fTemp[4] = (float) tRAXIS_header.nOsc;
	  fTemp[5] = 0.0;
	  fTemp[6] = 0.0;
	  fTemp[7] = 0.0;
	  fTemp[8] = 0.0;
	  fTemp[9] = 0.0;
	  (void) poHeader->nReplaceValue(Cstring (D_K_Rotation),
					 10, fTemp, 3);
	  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_Rotation,
					 10, fTemp, 3);
	  if (nRotAxis==3)
	    {
	      (void) poHeader->nReplaceValue(Cstring (D_K_RotAxisName),
					     Cstring ("phi"));
	      (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) +
					     D_K_RotAxisName,
					     Cstring ("phi"));
	      /*
	      (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioNames,
					     Cstring ("Phi LgArc SmArc"));
	      (void) poHeader->nReplaceValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
					     Cstring ("1 0 0 0 1 0 0 0 1"));
	      */


	    }
	  else if (nRotAxis==2)
	    {
	      (void) poHeader->nReplaceValue(Cstring (D_K_RotAxisName),
					     Cstring ("chi"));
	      (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) +
					     D_K_RotAxisName,
					     Cstring ("chi"));
	    } 
	  else if (nRotAxis==1)
	    {
	      (void) poHeader->nReplaceValue(Cstring (D_K_RotAxisName),
					     Cstring ("omega"));
	      (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) +
					     D_K_RotAxisName,
					     Cstring ("omega"));
	    }

	  //+jwp 22-Sep-2000
	  // Now use Chi for the Chi angle and ignore fMu!
	  // Rotation axis is normally along X (1,0,0) when chi=0
	  // and rotate (1,0,0) accordingly by fChi around (0, 0, 1);

	  Cimage_header *poTempHeader;
	  poTempHeader = new Cimage_header(poHeader->sGet().string());
	  Cgoniometer *poGoniometer;
	  poGoniometer = new Cgoniometer(*poTempHeader, D_K_CrystalPrefix);
	  delete poTempHeader;
	  if (poGoniometer->bIsAvailable())
	    {
	      poGoniometer->nGetRotVector(nRotAxis-1, (float *)&fTemp[0]);
	    }
	  else
	    {
	      // The X vector which is the default rotation axis
	      
	      fTemp[0] = 1.0f;
	      fTemp[1] = 0.0f;
	      fTemp[2] = 0.0f;
	      cout << "WARNING crystal goniometer vector problem!\n";
	    }
	  delete poGoniometer;
	  (void) poHeader->nReplaceValue(Cstring (D_K_RotVector),
					 3, &fTemp[0], 4);
	  (void) poHeader->nReplaceValue(Cstring (D_K_ScanPrefix) + D_K_RotVector,
					 3, &fTemp[0], 4);
	}
    
    int         a3nSeqInfo[3] = {0, 1, 0};
    bool        bHaveValidTemplateFromHeader = false;
    Cstring     sFileBase = sFileGetBasename(rsFilename);
    Cstring     sTemp("");

#ifndef RAXIS_HEADER_BINARY_STRESS_INFO 
      strncpy(ac,tRAXIS_header.a204cScanTemplate,204);
      s = ac;
      memset(ac,0,204);
      vTrimTrailingSpaces(s);

        // If scan template contains a ?, then this is probably a valid template name
        // If not, this is either an image collected on an SGI or an older PC image,
        // so figure it out from the filename

      if(s.contains('?'))
        {
          bHaveValidTemplateFromHeader = true;
          sTemp = s;

          // Get the sequence number of an image from its filename and the
          // current template.
          
          Cstring       sScanBase = sFileGetBasename(s);
          int nSeq     = 0;
          int nDecimal = 1;
          int i;
          char cFile;
          char cScan;
          for(i = sScanBase.length()-1; i >= 0; i--)
            {
              cFile = sFileBase.GetAt(i);
              cScan = sScanBase.GetAt(i);
              // Filename different from scan template,
              if(cFile != cScan)
                {
                  // Not a legitimate wildcard char in the scan template
                  if(('?' != cScan) && ('#' != cScan))
                    {
                      bHaveValidTemplateFromHeader = false;
                      break;  // out of for(i >= 0) loop
                    }
                  // Only 1 minus sign allowed, so we are done
                  if('-' == cFile)
                    {
                      nSeq = -nSeq;
                      break;  // out of for(i >= 0) loop
                    }
                  // A digit between 0 and 9, inclusive
                  else if(('0' <= cFile) && ('9' >= cFile))
                    {
                      nSeq      = nSeq + nDecimal * ((int) cFile - (int) '0');
                      nDecimal *= 10;
                    }
                  else
                    {
                      bHaveValidTemplateFromHeader = false;
                      break;  // out of for(i >= 0) loop
                    }
                }
            }
          if(!bHaveValidTemplateFromHeader)
            {
              //             cout << "\nWARNING Name mismatch between filename: ";
              //             cout << sFileBase << endl;
              // cout << "and template found in header: " << sScanBase << endl;
              a3nSeqInfo[0] = nSeq;
            }
          else
            a3nSeqInfo[0] = nSeq;
      }
#else
        bGetSetStressMeasurementInformation(tRAXIS_header, *poHeader, nSwap == 1 );
#endif

      // If could not find ? in template string or couldn't get sequence info
      // from the template and filename, then use the default method which
      // gets the template and sequence info from the filename.
      if( !bHaveValidTemplateFromHeader )
      {
          // Determine scan template and seq info from the filename
          // Change last 4 digits in rsFilename to ?s
          sTemp = sBuildScanTemplate(rsFilename, 4, &a3nSeqInfo[0]);
      }

      // Prepend a directory specification if it is not already there
      if (sTemp == sFileGetBasename(sTemp))
        sTemp = sGetCWD() + sTemp;

      (void) poHeader->nReplaceValue(Cstring (D_K_ScanTemplate),
                                     sTemp);

      (void) poHeader->nReplaceValue(Cstring (D_K_ScanSeqInfo),
                                     3, a3nSeqInfo);

      // Read any padding bytes left over in the header.
      // The amount of padding is
      // calculated from the record length in the original header minus
      // the header size.

//      cout << "Record length in bytes is: "
//           << tRAXIS_header.a2nRecord[0] << endl;
      nLen = tRAXIS_header.a2nRecord[0] - sizeof(tRAXIS_header); // How much padding?
      if (0 < nLen)
        {                                  // If there is padding, skip it
          char *pcTemp;
          pcTemp = new char[nLen];
          dskbr((int *)pnFile, pcTemp, &nLen, &nStat);
          delete [] pcTemp;
        }
//Debugging:      poHeader->nList();
    }
  if ( (0 != nStat) && ("" != sGetEnv("DTREK_RAXIS_IGNORE")))
    {
      cout << "Warning in raxis, error reading header, but ignored!\n";
      nStat = 0;
    }

  if ((!nStat) && (NULL != pcBuffer))
    {
      memcpy((void*)&g_oRaxisHeader,(void*)&tRAXIS_header, sizeof(RAXIS_header));
  };

  return (nStat);
}

int
nReadRAXISdata(const int *pnFile, unsigned short int *puiData,
               Cimage_header *poHeader)
{
  // Read RAXIS binary data with the dskio routines.  The raxis file is
  // already open, the header has been read, and the read position is at
  // the start of binary data.
  // This routine is normally called from the Cimage class.
  //
  // *pnFile   pointer to integer file tag for dskio routines
  // *puiData  pointer to where to put the read data
  // *pvHeader pointer to current image header
  //

  int i, j;            // Loop counter
  int nStat;           // Local status
  int nPad;            // Number of padding bytes at end of each line;
  int nDim0, nDim1;    // Dimensions of the image.
  int nBytesPerLine;   // Can you say Bytes Per Line?

  int nReadStart;      // Starting line of readout
  int nReadLines;      // Number of lines in the file
  int nToRead;         // Number bytes in a record

  nStat = poHeader->nGetValue(D_K_Size1, &nDim0);
  nStat = nStat + poHeader->nGetValue(D_K_Size2, &nDim1);
  nStat = nStat + poHeader->nGetValue(D_K_RaxisReadStart, &nReadStart);
  nStat = nStat + poHeader->nGetValue(D_K_RaxisReadLines, &nReadLines);
//+jwp kludge around bogus image conversion software
  if ( (1 == nReadStart) && ( (nReadStart + nReadLines) > nDim1) )
    nReadStart = 0;
//-jwp
  nStat = nStat + poHeader->nGetValue(D_K_RaxisRecordSize, &nToRead);
  if (0 == nStat)
    {
      // Each line might be padded, so figure this out

      nBytesPerLine = nDim0 * 2;
      nPad          = nToRead - nBytesPerLine;

      // Use one read for every line including the padding, but increment
      // the location of where to put the lines by only the amount of data.
      // Of course, do not do this for the last line as we will go beyond
      // our bit of allocated memory for the data.

      unsigned short int *puiTemp;
      puiTemp     = puiData;
      nToRead     = nBytesPerLine + nPad;

      // Check for a partial image.
      // Fill in the lines at start of image not in file.

      for (i = 0; i < nReadStart; i++)
        {
          for (j = 0; j < nDim0; j++)
            {
              *puiTemp++ = 0;
            }
        }

      // Read the lines in the file (except the last line)

      if (0 >= nPad)
        {
          // No padding, so can use single fast read
          nToRead = nToRead * nReadLines;
          (void) dskbr((int *)pnFile, (char *)puiTemp, &nToRead, &nStat);
          puiTemp = puiTemp + (nToRead / 2);
        }
      else
        {
          for (i = 0; (i < nReadLines-1) && (0 == nStat); i++)
            {
              (void) dskbr((int *)pnFile, (char *)puiTemp, &nToRead, &nStat);
              puiTemp = puiTemp + nDim0;
            }

          // Now read the last line

            if (0 == nStat)
              {
                (void) dskbr((int *)pnFile, (char *)puiTemp, &nBytesPerLine, &nStat);
                puiTemp = puiTemp + nDim0;
              }
        }

      // Fill in rest of lines not in file if partial readout

      for (i = (nReadStart+nReadLines); i < nDim1; i++)
        {
          for (j = 0; j < nDim0; j++)
            {
              *puiTemp++ = 0;
            }
        }
    }

//+jwp 23-Oct-2003
  // Special R-AXIS2 de-jitter code
  int nDeJit = 0;
  poHeader->nGetValue("DTREK_RAXIS_DEJIT", &nDeJit);
  if (0 != nDeJit)
    {
      // Shift even lines by nDeJit; Shift odd lines by -nDeJit

      unsigned short int *puiSrc;
      unsigned short int *puiDst;
      int nBytes = 2 * nDim1;
      int nDstOffset = 0;
      int nSrcOffset = 0;
      int nStart = 0;  // Shift EVEN lines
      if (0 > nDeJit)
	{
	  nStart = 1;  // Shift ODD lines
	  nDeJit = -nDeJit;
	}
      nDstOffset = 0;
      nSrcOffset = nDeJit;
      for (i = nStart; i < nDim1; i = i + 2)
	{
	  puiSrc = puiData + (i * nDim0) + nSrcOffset;
	  puiDst = puiData + (i * nDim0) + nDstOffset;
	  for (j = 0; j < nDim0; j++)
	    {
	      *puiDst++ = *puiSrc++;
	    }
	}
    }
//-jwp 23-Oct-2003

  if ( (0 != nStat) && ("" != sGetEnv("DTREK_RAXIS_IGNORE")))
    {
      cout << "Warning in raxis, error reading file, but ignored!\n";
      nStat = 0;
    }
  return (nStat);
}


int nBuildRAXISHeader(RAXIS_header& rtHeader,Cimage& oImage) {

  int i, j, nStat;
  float fTemp[10];
  Cstring sPrefix("");
  Cstring sTemp("");
  strncpy(rtHeader.a10cDevice, "R-AXIS    ", 10);
  strncpy(rtHeader.a10cVersion,"CCD       ", 10);
  strncpy(rtHeader.a20cCrystal,   "Unknown             ", 20);
  strncpy(rtHeader.a12cCry_system,"Unknown     ", 12);
  rtHeader.a6fCell[0]     = 1.0;
  rtHeader.a6fCell[1]     = 1.0;
  rtHeader.a6fCell[2]     = 1.0;
  rtHeader.a6fCell[3]     = 90.0;
  rtHeader.a6fCell[4]     = 90.0;
  rtHeader.a6fCell[5]     = 90.0;

  // We should not use m_oHeader.nGetValue since it defeats encapsulation!!!

  nStat = oImage.m_oHeader.nGetValue(Cstring(D_K_CrystalPrefix) + D_K_CrystalXUnitCell, 6, rtHeader.a6fCell);
  strncpy(rtHeader.a12cSpace, "P1          ", 12);
  rtHeader.fMosaic      = 0.20f;
  nStat = oImage.m_oHeader.nGetValue(Cstring(D_K_CrystalPrefix) + D_K_CrystalXMosaicity, &rtHeader.fMosaic);

  strncpy(rtHeader.a80cMemo, "memo                                                                            ", 80);
  strncpy(rtHeader.a84cReserve1, "Reserve1                                                                            ", 84);
  strncpy(rtHeader.a12cDate, "39-Jug-2001 ", 12);
  strncpy(rtHeader.a20cOperatorname, "Unknown d*TREK      ", 20);
  strncpy(rtHeader.a4cTarget, "Unk ", 4);

  rtHeader.fWave        = 1.54178f;
  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_SourcePrefix) +
                                   D_K_Wavelength, 2, fTemp);
  if (0 == nStat) rtHeader.fWave = fTemp[1];

  strncpy(rtHeader.a20cMonochro, "Unknown             ", 20);
  rtHeader.fMono_2      = 0.0;
  strncpy(rtHeader.a20cCollimator, "Unknown             ", 20);
  strncpy(rtHeader.a4cFilter, "Unk ", 4);
  rtHeader.fCamera      = 100.0;

  nStat = oImage.m_oHeader.nGetValue(D_K_DetectorNames, &sPrefix);
  if (0 != nStat) sPrefix = "";
  nStat = oImage.m_oHeader.nGetValue(sPrefix + D_K_GonioValues, 6,
                              fTemp);
  if (0 == nStat) rtHeader.fCamera = fTemp[5];

  rtHeader.fKv          = 0.0;
  rtHeader.fMa          = 0.0;
  rtHeader.nCylinder    = 0;
  rtHeader.fWeissenberg = 0.0;
  strncpy(rtHeader.a12cFocus, "Unknown     ", 12);
  strncpy(rtHeader.a80cXraymemo, "xraymemo                                                                        ", 80);
  strncpy(rtHeader.a56cReserve2, "reserve2                                                ", 56);
  strncpy(rtHeader.a4cSpindle, "Unk ", 4);
  strncpy(rtHeader.a4cXray_axis, "Unk ", 4);

  rtHeader.a3fPhi[0]    = 0.0;
  rtHeader.a3fPhi[1]    = 0.0;
  rtHeader.a3fPhi[2]    = 0.0;
  rtHeader.nOsc         = 1;
  rtHeader.fEx_time     = 0.0;
  rtHeader.a2fXray[0]   = 100.0;  // X-ray beam position in pixels
  rtHeader.a2fXray[1]   = 100.0;

  nStat = oImage.nGetDimensions(&i, &j);
  if (0 < nStat)
    {
      rtHeader.a2nPix_num[0]  = i;
      rtHeader.a2nPix_num[1]  = j;
      rtHeader.a2fXray[0] = (float) i * 0.5f;
      rtHeader.a2fXray[1] = (float) j * 0.5f;
    }

  rtHeader.a3fCircle[0]  = 0.0;    // Crystal goniometer datum positions, omega
  rtHeader.a3fCircle[1]  = 0.0;    // chi
  rtHeader.a3fCircle[2]  = 0.0;    // Detector goniometer datum positions

  nStat = oImage.m_oHeader.nGetValue(D_K_Rotation, 10, fTemp);
  if (0 == nStat)
    {
      rtHeader.a3fPhi[1] = fTemp[0];
      rtHeader.a3fPhi[2] = fTemp[1];
      rtHeader.fEx_time = fTemp[3] / 60.0f;
      rtHeader.nOsc     = (int) fTemp[4];
    }

  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
                                   3, fTemp);
  if (0 == nStat)
    {
      rtHeader.a3fCircle[0] = fTemp[0]; // Omega
      rtHeader.a3fCircle[1] = fTemp[1]; // Chi
      rtHeader.a3fPhi[0]    = fTemp[2]; // Phi datum
    }

  rtHeader.fMu          = 0.0;
  //for (i = 0; i < 51; i++) rtHeader.a51nReserve3[i] = 0;


#ifndef RAXIS_HEADER_BINARY_STRESS_INFO
  Cstring       s("");
  oImage.m_oHeader.nGetValue(D_K_ScanTemplate,&s);
  memset(rtHeader.a204cScanTemplate,' ',204);
  strncpy(rtHeader.a204cScanTemplate,s.string(),s.length());
#else
  rtHeader.posn = 0;
  rtHeader.side = 0;
  rtHeader.slit = 0;
  
  for (i = 0; i < 64; i++) 
      rtHeader.posz[i] = ' ';

  for (i = 0; i < 128; i++) 
      rtHeader.psis[i] = ' ';
// todo: get stress information from d*TREK header to RAXIS header
#endif

  rtHeader.a2fPix_size[0] = 0.100f;
  rtHeader.a2fPix_size[1] = 0.100f;


  nStat = oImage.m_oHeader.nGetValue(sPrefix + D_K_SpatialDistortionInfo,
                                   4, fTemp);
  if (0 == nStat)
    {
      // Set pixel size and beam center from spatial distortion info
      // But check for spatial beam position

      nStat = oImage.m_oHeader.nGetValue(sPrefix + D_K_SpatialBeamPosition,
                                         2,fTemp);

      rtHeader.a2fXray[0]     = fTemp[0]; // X-ray beam position in pixels
      rtHeader.a2fXray[1]     = fTemp[1];
      rtHeader.a2fPix_size[0] = fTemp[2];
      rtHeader.a2fPix_size[1] = fTemp[3];
    }
  else
    {
      // Try to get pixel size from detector size and detector dimensions

      nStat = oImage.m_oHeader.nGetValue(sPrefix + D_K_DetectorSize, 2, fTemp);
      nStat = nStat +  oImage.m_oHeader.nGetValue(sPrefix + D_K_DetectorDimensions,
                                                2, &fTemp[2]);
      if (0 == nStat)
        {
          rtHeader.a2fPix_size[0] = fTemp[0] / fTemp[2];
          rtHeader.a2fPix_size[1] = fTemp[1] / fTemp[3];
        }
    }

  // First get minimum size of a record.  Then see if a larger size
  // is suggested by the header.   Images derived from R-AXIS
  // images will have a RAXIS_RECORD_SIZE keyword.  Others must have
  // records at least the size of the header.

  int nTemp;
  nTemp                 = 2 * rtHeader.a2nPix_num[0];
  rtHeader.a2nRecord[0] = max(nTemp,1024); // Record length minimum is 1024
  nStat = oImage.m_oHeader.nGetValue(sPrefix + D_K_RaxisRecordSize, &nTemp);
  if ( (0 == nStat) && (nTemp > rtHeader.a2nRecord[0]) )
    {
      // Found a larger suggested record size, so use it
      rtHeader.a2nRecord[0]   = nTemp;
    }
  rtHeader.a2nRecord[1] = rtHeader.a2nPix_num[1];
  rtHeader.nRead_start  = 0;
  rtHeader.nIP_num      = 1;
  rtHeader.fRatio       = 1.0;
  rtHeader.a2fFading[0] = 1.0;
  rtHeader.a2fFading[1] = 1.0;
#ifdef WIN32
  strncpy(rtHeader.a10cCpu, "WIN32 PC  ", 10);
#else
  strncpy(rtHeader.a10cCpu, "IRIS      ", 10);
#endif
  strncpy(rtHeader.a10cIp,  "Unknown   ", 10);
  rtHeader.a3nDrxz[0]     = 0;
  rtHeader.a3nDrxz[1]     = 0;
  rtHeader.a3nDrxz[2]     = 0;

  //+ 3-Apr-2002 What to do with the new header structure elements?

  rtHeader.fPixShiftOdd   = 0.0;
  rtHeader.fIntRatioOdd   = 1.0;
  rtHeader.nMagicNum      = 0;
  memset((void*)&rtHeader.nMagicNum, 0, 172);

  //+JWP 2009-06-08

  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_CrystalPrefix) +
				 D_K_GonioNumValues,
				 &nTemp);
  if (0 == nStat) 
    rtHeader.nNumGonAxes     = nTemp;
  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_CrystalPrefix) + D_K_GonioValues,
				     3, fTemp);
  if (0 == nStat)
    {
      rtHeader.a5fGonStart[0] = fTemp[0];
      rtHeader.a5fGonStart[1] = fTemp[1];
      rtHeader.a5fGonStart[2] = fTemp[2];
      rtHeader.a5fGonEnd[0]   = fTemp[0];
      rtHeader.a5fGonEnd[1]   = fTemp[1];
      rtHeader.a5fGonEnd[2]   = fTemp[2];
    }

  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_CrystalPrefix) + D_K_GonioVectors,
				     9, fTemp);
  
  rtHeader.a5x3fGonVecs[0][0] = fTemp[0];
  rtHeader.a5x3fGonVecs[0][1] = fTemp[1];
  rtHeader.a5x3fGonVecs[0][2] = fTemp[2];
  rtHeader.a5x3fGonVecs[1][0] = fTemp[3];
  rtHeader.a5x3fGonVecs[1][1] = fTemp[4];
  rtHeader.a5x3fGonVecs[1][2] = fTemp[5];
  rtHeader.a5x3fGonVecs[2][0] = fTemp[6];
  rtHeader.a5x3fGonVecs[2][1] = fTemp[7];
  rtHeader.a5x3fGonVecs[2][2] = fTemp[8];
  nStat = oImage.m_oHeader.nGetValue(Cstring(D_K_ScanPrefix) + Cstring (D_K_RotAxisName), &sTemp);
  sTemp.upcase();
  //std::cout << "STEMP for AXIS: >" << sTemp << "<<<<\n"; 
  if ('O' == sTemp.GetAt(0))
    rtHeader.nScanAxisNum = 0;
  else if ('C' == sTemp.GetAt(0))
    rtHeader.nScanAxisNum = 1;
  else if ('P' == sTemp.GetAt(0))
    rtHeader.nScanAxisNum = 2;
  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_Rotation), 10, fTemp);
  std::cout << "TEMP2: " << fTemp[2] << std::endl;
  // The 3rd value is the image rotation width
  if (0 == nStat)
    {
      nTemp = rtHeader.nScanAxisNum; 
      rtHeader.a5fGonStart[nTemp] = fTemp[0];
      rtHeader.a5fGonEnd[nTemp]   = fTemp[1];
      //rtHeader.a5fGonEnd[nTemp] = rtHeader.a5fGonEnd[nTemp] + fTemp[2];      
    }

  nStat = oImage.m_oHeader.nGetValue(Cstring (D_K_CrystalPrefix) + D_K_GonioNames,
				     &sTemp);
  strncpy(rtHeader.a40cAxesNames, sTemp.string(), sTemp.length());

  if (0 == nStat)
    rtHeader.nMagicNum       = RAXIS_MAGIC_NUM;

  //-JWP 2009-06-08

  //strncpy(rtHeader.a180cReserve4, "Unknown                                                                                                                                                                             ", 180);
  //- 3-Apr-2002

  // Header is complete.
  return 0;
};


int
nWriteRAXIS(Cimage & oImage, const Cstring& sName,RAXIS_header* ptKnownHeader)
{
    int i,j;
    RAXIS_header* ptHeader;
    RAXIS_header  tRAXISHeader;

  // Write an image to a RAXIS style file

  
    // Build the header.
  if (ptKnownHeader)
      // If the user provided a header, use that.
      ptHeader = ptKnownHeader;
  else {
      // Otherwise, build a header.
      ptHeader = &tRAXISHeader;
      nBuildRAXISHeader(*ptHeader,oImage);
  };
  

  int nFile, nBytes, nSize, nPad, istat;

  nFile = 1;

  // Compute total size of image file in bytes.  This is needed for dskbow
  // when on VMS platforms.
  // Compute padding needed after header.

  nPad   = ptHeader->a2nRecord[0] - sizeof(RAXIS_header);
  nSize  = ptHeader->a2nRecord[0] * (ptHeader->a2nRecord[1] + 1);
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
                              &ptHeader->fRatio);

// The file should always be written out using the same endian-ness as that
// on an SGI, which is big endian

   if(eCPU_big_endian != nGetByteOrderCPU()){
      ptHeader->nOsc        = nSwapLong(ptHeader->nOsc);
      ptHeader->a2nPix_num[0] = nSwapLong(ptHeader->a2nPix_num[0]);
      ptHeader->a2nPix_num[1] = nSwapLong(ptHeader->a2nPix_num[1]);
      ptHeader->a2nRecord[0]  = nSwapLong(ptHeader->a2nRecord[0]);
      ptHeader->a2nRecord[1]  = nSwapLong(ptHeader->a2nRecord[1]);
      ptHeader->nRead_start   = nSwapLong(ptHeader->nRead_start);
      ptHeader->nIP_num       = nSwapLong(ptHeader->nIP_num);
      ptHeader->a3nDrxz[0]    = nSwapLong(ptHeader->a3nDrxz[0]);
      ptHeader->a3nDrxz[1]    = nSwapLong(ptHeader->a3nDrxz[1]);
      ptHeader->a3nDrxz[2]    = nSwapLong(ptHeader->a3nDrxz[2]);
      for(i = 0; i < 6; i++)
         ptHeader->a6fCell[i] = fSwapFloat(ptHeader->a6fCell[i]);
      ptHeader->fMosaic       = fSwapFloat(ptHeader->fMosaic);
      ptHeader->fWave         = fSwapFloat(ptHeader->fWave);
      ptHeader->fMono_2       = fSwapFloat(ptHeader->fMono_2);
      ptHeader->fCamera       = fSwapFloat(ptHeader->fCamera);
      ptHeader->fKv           = fSwapFloat(ptHeader->fKv);
      ptHeader->fMa           = fSwapFloat(ptHeader->fMa);
      ptHeader->fWeissenberg  = fSwapFloat(ptHeader->fWeissenberg);
      ptHeader->a3fPhi[0]     = fSwapFloat(ptHeader->a3fPhi[0]);
      ptHeader->a3fPhi[1]     = fSwapFloat(ptHeader->a3fPhi[1]);
      ptHeader->a3fPhi[2]     = fSwapFloat(ptHeader->a3fPhi[2]);
      ptHeader->fEx_time      = fSwapFloat(ptHeader->fEx_time);
      ptHeader->a2fXray[0]    = fSwapFloat(ptHeader->a2fXray[0]);
      ptHeader->a2fXray[1]    = fSwapFloat(ptHeader->a2fXray[1]);
      ptHeader->a3fCircle[0]  = fSwapFloat(ptHeader->a3fCircle[0]);
      ptHeader->a3fCircle[1]  = fSwapFloat(ptHeader->a3fCircle[1]);
      ptHeader->a3fCircle[2]  = fSwapFloat(ptHeader->a3fCircle[2]);
      ptHeader->fMu           = fSwapFloat(ptHeader->fMu);
      ptHeader->a2fPix_size[0]= fSwapFloat(ptHeader->a2fPix_size[0]);
      ptHeader->a2fPix_size[1]= fSwapFloat(ptHeader->a2fPix_size[1]);
      ptHeader->fRatio        = fSwapFloat(ptHeader->fRatio);
      ptHeader->a2fFading[0]  = fSwapFloat(ptHeader->a2fFading[0]);
      ptHeader->a2fFading[1]  = fSwapFloat(ptHeader->a2fFading[1]);

      ptHeader->nMagicNum     = nSwapLong(ptHeader->nMagicNum);
      ptHeader->nNumGonAxes   = nSwapLong(ptHeader->nNumGonAxes);
      ptHeader->nScanAxisNum  = nSwapLong(ptHeader->nScanAxisNum);

      float *pfTemp;
      pfTemp = &ptHeader->a5x3fGonVecs[0][0];
      for (i = 0; i < 30; i++)
	{
	  *pfTemp = fSwapFloat(*pfTemp);
	  pfTemp++;
	}

   }

  // Write the header

  nBytes = sizeof(RAXIS_header);
  dskbw(&nFile, (char *) ptHeader, &nBytes, &istat);

  // Write the padding required after the header, use the header for the pad
  // tjn: NO  If the header buffer size is not large enough, an exception 
  // will be thrown.  Allocate memory on the spot for this purpose.

  if ( (0 == istat) && (0 < nPad) )
    {
        char* cpTemp=new char[nPad+1];
      dskbw(&nFile, cpTemp, &nPad, &istat);
        delete[] cpTemp;
    }

// Swap the header back because some of it's values are used later in this routine

   if(eCPU_big_endian != nGetByteOrderCPU()){
      ptHeader->nOsc        = nSwapLong(ptHeader->nOsc);
      ptHeader->a2nPix_num[0] = nSwapLong(ptHeader->a2nPix_num[0]);
      ptHeader->a2nPix_num[1] = nSwapLong(ptHeader->a2nPix_num[1]);
      ptHeader->a2nRecord[0]  = nSwapLong(ptHeader->a2nRecord[0]);
      ptHeader->a2nRecord[1]  = nSwapLong(ptHeader->a2nRecord[1]);
      ptHeader->nRead_start   = nSwapLong(ptHeader->nRead_start);
      ptHeader->nIP_num       = nSwapLong(ptHeader->nIP_num);
      ptHeader->a3nDrxz[0]    = nSwapLong(ptHeader->a3nDrxz[0]);
      ptHeader->a3nDrxz[1]    = nSwapLong(ptHeader->a3nDrxz[1]);
      ptHeader->a3nDrxz[2]    = nSwapLong(ptHeader->a3nDrxz[2]);
      for(i = 0; i < 6; i++)
         ptHeader->a6fCell[i] = fSwapFloat(ptHeader->a6fCell[i]);
      ptHeader->fMosaic       = fSwapFloat(ptHeader->fMosaic);
      ptHeader->fWave         = fSwapFloat(ptHeader->fWave);
      ptHeader->fMono_2       = fSwapFloat(ptHeader->fMono_2);
      ptHeader->fCamera       = fSwapFloat(ptHeader->fCamera);
      ptHeader->fKv           = fSwapFloat(ptHeader->fKv);
      ptHeader->fMa           = fSwapFloat(ptHeader->fMa);
      ptHeader->fWeissenberg  = fSwapFloat(ptHeader->fWeissenberg);
      ptHeader->a3fPhi[0]     = fSwapFloat(ptHeader->a3fPhi[0]);
      ptHeader->a3fPhi[1]     = fSwapFloat(ptHeader->a3fPhi[1]);
      ptHeader->a3fPhi[2]     = fSwapFloat(ptHeader->a3fPhi[2]);
      ptHeader->fEx_time      = fSwapFloat(ptHeader->fEx_time);
      ptHeader->a2fXray[0]    = fSwapFloat(ptHeader->a2fXray[0]);
      ptHeader->a2fXray[1]    = fSwapFloat(ptHeader->a2fXray[1]);
      ptHeader->a3fCircle[0]  = fSwapFloat(ptHeader->a3fCircle[0]);
      ptHeader->a3fCircle[1]  = fSwapFloat(ptHeader->a3fCircle[1]);
      ptHeader->a3fCircle[2]  = fSwapFloat(ptHeader->a3fCircle[2]);
      ptHeader->fMu           = fSwapFloat(ptHeader->fMu);
      ptHeader->a2fPix_size[0]= fSwapFloat(ptHeader->a2fPix_size[0]);
      ptHeader->a2fPix_size[1]= fSwapFloat(ptHeader->a2fPix_size[1]);
      ptHeader->fRatio        = fSwapFloat(ptHeader->fRatio);
      ptHeader->a2fFading[0]  = fSwapFloat(ptHeader->a2fFading[0]);
      ptHeader->a2fFading[1]  = fSwapFloat(ptHeader->a2fFading[1]);
   }

  // Figure out padding needed after every line

  nPad  = ptHeader->a2nRecord[0] - sizeof(short int) * ptHeader->a2nPix_num[0];

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
      puiData = new unsigned short int [  ptHeader->a2nPix_num[0]
                                        * ptHeader->a2nPix_num[1]
                                        * sizeof(unsigned short int)];
      puiTemp = puiData;

      for (j = 0; j < ptHeader->a2nPix_num[1]; j++)
        {
          for (i = 0; i < ptHeader->a2nPix_num[0]; i++)
            {
              fPixel = (oImage.*oImage.prfGetPixel)(i, j);

              if (fPixel > 32767.0)
                {
                  if (ptHeader->fRatio > 1.0)
                    {
                      *puiTemp = (unsigned short int) (fPixel / ptHeader->fRatio
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

// Do we need to swap the byte order of the data as well?
   if(eCPU_big_endian != nGetByteOrderCPU()){
      int nBytes = oImage.nDataSizeNeeded();
      unsigned short int *pTemp;
      pTemp = puiTemp;
      for(i = 0; i < nBytes/2; i++){
         *pTemp = ((*pTemp & 0xff) << 8) | ((*pTemp) >> 8);
         pTemp++;
      }
   }

  // Now write out image

  if (0 >= nPad)
    {
      // No padding needed, so write entire image in one shot
      nBytes = sizeof(short int) * ptHeader->a2nPix_num[0]
                                 * ptHeader->a2nPix_num[1];
      dskbw(&nFile, (char *)puiTemp, &nBytes, &istat);
    }
  else
    {
      // Padding required, so write line by line and add padding
      nBytes = sizeof(short int) * ptHeader->a2nPix_num[0] + nPad;
      for (i = 0; (i < ptHeader->a2nPix_num[1]-1) && (0 == istat); i++)
        {
          // Write out one line at a time, note that we pad not with 0's
          // but with the data from the next line so that we do not need
          // an extra call to dskbw for the padding.
        
          dskbw(&nFile, (char *)puiTemp, &nBytes, &istat);
          puiTemp = puiTemp + ptHeader->a2nPix_num[0];
        }
      // Write out last line

      if (0 == istat)
        {
          nBytes = sizeof(short int) * ptHeader->a2nPix_num[0];
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
  // Don't forget to swap back if we were using the actual data
   else if(eCPU_big_endian != nGetByteOrderCPU()){
      int nBytes = oImage.nDataSizeNeeded();
      puiTemp = oImage.m_The_Data.pui;  // Access data directly
      for(i = 0; i < nBytes/2; i++){
         *puiTemp = ((*puiTemp & 0xff) << 8) | ((*puiTemp) >> 8);
         puiTemp++;
      }
   }

  return (istat);
}

////////////////////////////////////////////////
//Remove trailing white spaces from a string
void vTrimTrailingSpaces(Cstring& s)
{
    int     ii = s.length();
    if(0 == ii)
        return;

    while( ii > 0 && ' ' == s.GetAt(ii-1) )
        ii--;

    if( 0 == ii )
        s = "";
    else
        s = s.before(ii);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Add shadow information to header, if it is not already defined.  See CGonioMask.h/.cc for more info.
void vSetHeaderRapidShadowInfo(Cimage_header* poHeader, bool bIsDMaxRapid)
{        
    int     nNumberOfShadows = 0;
    (void) poHeader->nGetValue("CRYSTAL_GONIO_NUM_SHADOWS", &nNumberOfShadows);
    if( nNumberOfShadows > 0 )
        return; // the header already has shadow information
    
    int     ii = 0; // shadow number

    if( !bIsDMaxRapid )
    {
        // CHI = -16 - -15
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 344 345 -999 999 150 22.166667 126 128 -1 2 0.1 0.1 3078 0 2587 567 913 567 706 0");
        
        // CHI = -10
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 350 350 -999 999 150 22.166667 126 128 -1 2 0.1 0.1 2993 0 2570 408 916 436 748 0");
        
        // CHI = -5
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 355 355 -999 999 150 22.166667 126 128 -1 2 0.1 0.1 2897 0 2603 255 1737 196 1717 0");
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 355 355 -999 999 150 22.166667 126 128 -1 2 0.1 0.1 1517 0 1557 125 973 301 777 157 737 0");
        
        // CHI = 0
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 0 0 -999 999 20 22.166667 126 128 -1 2 0.1 0.1 4000 0 4000 80 3500 80 3500 0");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 0 0 -999 999 20 22.166667 126 128 -1 2 0.1 0.1 1200 0 1180 1040 290 1040 260 0");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 0 0 -999 999 20 22.166667 126 128 -1 2 0.1 0.1 950 990 950 1100 530 1100 530 990");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 0 0 -999 999 80  22.166667 126 128 -1 2 0.1 0.1 3590 0 3590 90 3310 90 3310 0");
        
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 0 0 -999 999 150 22.166667 126 128 -1 2 0.1 0.1 2535 0 2535 27 2233 26 2233 0");

        // CHI = 10
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 10 10.0 -999 999 -30 22.166667 126 128 -1 2 0.1 0.1 2319 0 2191 1114 1418 1112 1405 795 1523 691 1407 0");

        // CHI = 20
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 20 20 -999 999 -30 22.166667 126 128 -1 2 0.1 0.1 2417 0 2212 1110 1420 1110 1420 800 1523 749 1282 0");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 20 20 -999 999 -30 22.166667 126 128 -1 2 0.1 0.1 1143 0 1106 66 887 70 865 0");
        
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 20 20 -999 999 -30 22.166667 126 128 -1 2 0.1 0.1 640 0 606 73 511 75 486 0");
        
        // CHI = 30
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 30 30 -999 999 -30 22.166667 126 128 -1 2 0.1 0.1 2711 0 2609 268 2478 222 2193 1104"
        " 1410 1118 1410 794 1519 748 1496 462 1412 393 1293 301 1109 403 866 153 752 284 680 236 697 0");
        
        // CHI = 40
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 40 40 -999 999 -30 22.166667 126 128 -1 2 0.1 0.1 3027 0 2801 277 2664 277 2512 532 2401 531 2210 1111"
        " 1372 1108 1372 702 1171 706 1019 492 778 465 752 0");
        
        // CHI = 45
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 45 45 -999 999 20 22.166667  126 128 -1 2 0.1 0.1 1830 0 1700 400 0 400 0 0");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 45 45 -999 999 20 22.166667 126 128 -1 2 0.1 0.1 1470 400 1400 650 0 650 0 400");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 45 45 -999 999 20 22.166667 126 128 -1 2 0.1 0.1 1240 650 1060 1100 260 1100 260 650");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 45 45 -999 999 20 22.166667 126 128 -1 2 0.1 0.1 260 650 260 740 0 740 0 650");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii), 
        "0 -999 999 45 45 -999 999 198 22.166667 126 128 -1 2 0.1 0.1 4050 530 3800 530 3720 0 4050 0");

        // CHI = 50
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 50 50 -999 999 -20 22.166667 126 128 -1 2 0.1 0.1 2710 0 2550 620 770 620 540 0");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 50 50 -999 999 -20 22.166667 126 128 -1 2 0.1 0.1 2290 620 2230 850 800 850  770 620");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 50 50 -999 999 -20 22.166667 126 128 -1 2 0.1 0.1 2090 850 2020 1080 1340 1080 920 850");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 50 50 -999 999 -30 22.166667 126 128 -1 2 0.1 0.1 970 627 929 612 916 519 977 519");

        // CHI = 54
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 54 54 -999 999 -20 22.166667 126 128 -1 2 0.1 0.1 2710 0 2550 620 770 620 540 0");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 54 54 -999 999 -20 22.166667 126 128 -1 2 0.1 0.1 2290 620 2230 850 800 850  770 620");

        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 54 54 -999 999 -20 22.166667 126 128 -1 2 0.1 0.1 2090 850 2020 1080 1340 1080 920 850");
    }
    else
    {
        // PHI = 0 - 10
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 1 9 113.5 22.166667 126 128 -1 2 0.1 0.1 4516 0 3721 1229 2318 1324 2141 1485"
        " 1824 1340 1685 1055 1337 698 1081 543 821 0");
        
        // PHI = 10 - 20
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 11 19 113.5 22.166667 126 128 -1 2 0.1 0.1 4410 0 4158 520 3638 1103 2127 1493"
        " 1821 1118 1280 744 1021 385 1030 212 1067 0");

        // PHI = 20 - 30
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 21 29 113.5 22.166667 126 128 -1 2 0.1 0.1 4484 0 3789 1085 2305 1603 1115 677"
        "  1170 120 1075 0");
        
        // PHI = 30 - 40
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 31 39 113.5 22.166667 126 128 -1 2 0.1 0.1 4328 0 3931 883 2780 1347 2736 1533 2342"
        " 1506 1249 852 1271 365 1162 280 1001 367 843 0");

        // PHI = 40 - 50
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 41 49 113.5 22.166667 126 128 -1 2 0.1 0.1 4520 0 3852 873 3098 1159 2929 1360 2906"
        "  1507 2554 1543 1859 1178 1613 1176 1394 903 1391 572 1293 486 1138 580 820 0");
        
        // PHI = 50 - 60
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 51 59 113.5 22.166667 126 128 -1 2 0.1 0.1 4429 0 3988 641 3067 1262 3123 1508 2751"
        " 1570 1505 1168 1505 747 1246 751 968 312 947 0");

        // PHI = 60 - 70
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 61 69 113.5 22.166667 126 128 -1 2 0.1 0.1 4420 0 3306 1083 3314 1434 2889 1536 1662"
        " 1272 1611 895 1521 828 1383 934 1012 361 984 0");

        // PHI = 70 - 80
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 71 79 113.5 22.166667 126 128 -1 2 0.1 0.1 4448 0 3405 1156 3426 1431 3076 1510 2985"
        " 1473 2335 1477 2165 1485 1811 1310 1756 1021 1725 956 1663 953 1528 1094 1112 576 1160 529 942 0");
        
        // PHI = 80 - 90
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 81 89 113.5 22.166667 126 128 -1 2 0.1 0.1 4246 0 3572 1031 3566 1335 3174 1441 2739"
        " 1509 1888 1440 1860 1108 1520 1144 1224 678 1010 0");
        
        // PHI = 90 - 100
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 91 99 113.5 22.166667 126 128 -1 2 0.1 0.1 4100 0 3673 975 3659 965 3754 1154 3395"
        " 1380 2177 1527 2044 1201 1940 1125 1808 1318 1341 894 1426 795 1166 431 1209 370 835 0");

        // PHI = 100 - 110
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 101 109 113.5 22.166667 126 128 -1 2 0.1 0.1 4104 0 3778 839 3908 981 3587 1265 2496"
        "  1580 2283 1468 2201 1278 2103 1275 2021 1390 1967 1390 1567 1133 960 312 841 0");

        // PHI = 110 - 120
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 111 119 113.5 22.166667 126 128 -1 2 0.1 0.1 4144 0 3939 692 4033 879 3711 1175 3448"
        "  1164 3194 1424 2531 1617 2274 1324 2105 1490 1672 1174 1390 731 1022 499 862 0");
        
        // PHI = 120 - 130
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 121 129 113.5 22.166667 126 128 -1 2 0.1 0.1 4088 0 4012 509 4163 683 3843 1084 3589"
        "  1051 3519 1219 3036 1523 2126 1517 1595 1127 1535 909 1169 776 969 356 1051 164 994 0");

        // PHI = 130 - 140
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 131 139 113.5 22.166667 126 128 -1 2 0.1 0.1 4275 0 4238 913 3828 938 3424 1293 3048"
        "  1550 2822 1520 2687 1395 2288 1594 1698 1189 1203 827 985 0");
        
        // PHI = 140 - 150
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 141 149 113.5 22.166667 126 128 -1 2 0.1 0.1 4147 0 4210 116 4404 139 4367 563 3557"
        "  1187 3555 1278 3049 1527 2835 1355 2516 1600 1857 1304 1378 971 1202 726 1167 215 948 0");

        // PHI = 150 - 160
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 151 159 113.5 22.166667 126 128 -1 2 0.1 0.1 4534 0 4420 462 3975 731 3972 844 3664"
        " 1106 3645 1255 3196 1480 2887 1337 2740 1583 1997 1437 1941 1240 1565 1167 1332 882 1260 451 1062 346 983 0");

        // PHI = 160 - 170
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 161 169 113.5 22.166667 126 128 -1 2 0.1 0.1 4530 0 4490 297 4101 608 3746 1205 3381"
        "  1388 2522 1553 1525 1203 925 0");

        // PHI = 170 - 180
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 171 179 113.5 22.166667 126 128 -1 2 0.1 0.1 4542 0 4183 365 4162 920 3300 1507 2237"
        " 1517 2237 1380 1761 1369 1513 878 1137 455 1280 187 1164 0");

        // PHI = 180 - 190
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 181 189 113.5 22.166667 126 128 -1 2 0.1 0.1 4373 0 4242 587 3778 1193 3413 1184 3480"
        " 1397 2894 1596 2890 1498 1968 1485 1195 576 1371 403 1002 0");

        // PHI = 190 - 200
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 191 199 113.5 22.166667 126 128 -1 2 0.1 0.1 4441 0 4349 643 3765 1086 3741 971 3568"
        " 1120 3630 1266 3082 1545 2097 1533 1318 757 1499 566 1176 307 1170 148 1053 0");

        // PHI = 200 - 210
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 201 209 113.5 22.166667 126 128 -1 2 0.1 0.1 4572 0 4572 491 3404 1472 2294 1629 1485"
        " 987 1562 765 1167 433 1099 0");

        // PHI = 210 - 220
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 211 219 113.5 22.166667 126 128 -1 2 0.1 0.1 4596 0 4144 907 3445 1390 2981 1515 2914"
        " 1418 2456 1665 1594 1053 1655 832 1194 567 1146 0");

        // PHI = 220 - 230
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 221 229 113.5 22.166667 126 128 -1 2 0.1 0.1 4548 0 4166 663 4020 687 4018 1028 2589"
        "  1680 1781 1210 1783 966 1247 909 1199 0");

        // PHI = 230 - 240
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 231 239 113.5 22.166667 126 128 -1 2 0.1 0.1 4570 0 4119 590 4237 792 2767 1697 2567"
        " 1473 1922 1293 1913 1024 1404 964 1279 0");

        // PHI = 240 - 250
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 241 249 113.5 22.166667 126 128 -1 2 0.1 0.1 4580 0 4210 419 4366 585 3462 1419 3466"
        " 1283 3014 1631 2608 1423 1984 1387 2077 1154 1492 1046 1328 0");

        // PHI = 250 - 260
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 251 259 113.5 22.166667 126 128 -1 2 0.1 0.1 4503 0 4323 209 4479 398 3944 1023 3945"
        " 922 3382 1435 2987 1508 2985 1462 2189 1442 2219 1197 1667 1186 1306 0");

        // PHI = 260 - 270
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 261 269 113.5 22.166667 126 128 -1 2 0.1 0.1 4586 0 4276 692 3792 1005 3560 1397 2477"
        " 1441 2467 1313 1743 1333 1371 0");

        // PHI = 270 - 280
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 271 279 113.5 22.166667 126 128 -1 2 0.1 0.1 4594 0 4313 539 3521 1367 2862 1364 2760"
        " 1460 2625 1448 2582 1358 2086 1357 1713 956 1388 274 1172 0");

        // PHI = 280 - 290
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 281 289 113.5 22.166667 126 128 -1 2 0.1 0.1 4598 0 4075 740 3981 742 3848 1207 2957"
        " 1430 2468 1352 2053 1384 2040 1309 1290 382 830 379 828 0");

        // PHI = 290 - 300
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 291 299 113.5 22.166667 126 128 -1 2 0.1 0.1 4442 0 3895 1106 3316 1154 3342 1375"
        " 2991 1424 2910 1301 2402 1422 1890 970 1378 500 1194 615 938 237 921 0");

        // PHI = 300 - 310
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 301 309 113.5 22.166667 126 128 -1 2 0.1 0.1 4427 0 4048 1008 3500 1064 3319 1435"
        " 3112 1446 3057 1306 2482 1420 1488 667 1168 717 1033 385 1002 0");

        // PHI = 310 - 320
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 311 319 113.5 22.166667 126 128 -1 2 0.1 0.1 4358 0 4215 833 3359 1289 2643 1409"
        " 1601 794 1422 902 1159 633 1082 0");
        
        // PHI = 320 - 330
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 321 329 113.5 22.166667 126 128 -1 2 0.1 0.1 4431 0 4323 684 3101 1544"
        "  1700 929 1590 1049 1291 899 1183 308 900 0");

        // PHI = 330 - 340
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 331 339 113.5 22.166667 126 128 -1 2 0.1 0.1 4503 0 4371 519 3951 745 3869 932"
        " 3550 1157 3470 1015 3317 1112 3339 1422 3043 1503 3037 1426 1885 1064 1634 1190 1387 925 1007 0");

        // PHI = 340 - 350
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 341 349 113.5 22.166667 126 128 -1 2 0.1 0.1 4550 0 4461 339 4220 339 4200 448 4000"
        "  606 3988 808 3722 977 3605 936 3373 1247 3299 1440 3082 1428 2012 1150 1762 1267 1506 1047 773 0");

        // PHI = 350 - 360
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 351 359 113.5 22.166667 126 128 -1 2 0.1 0.1 4310 0 4288 243 4092 405 4171 568 3873"
        " 923 3765 833 3141 1392 2196 1178 2066 1410 1672 1176 1205 474 989 381 881 0");

        
        // PHI = 0
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),                       // #36
        "0 -999 999 45.0 45.0 -1 1 113.5 22.166667 126 128 -1 2 0.1 0.1 4328 0 3462 1200 1708 1194 754 0");
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "1 -999 999 45.0 45.0 -1 1 113.5 22.166667 126 128 -1 2 0.1 0.1 1966 1194 112 0 0 0 0 0");
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "1 -999 999 45.0 45.0 -1 1 113.5 22.166667 126 128 -1 2 0.1 0.1 3967  515 147 0 0 0 0 0");

        // PHI = 10 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 9 11 113.5 22.166667 126 128 -1 2 0.1 0.1  4118  0 4075 104 4239  305 4138"
        " 533 3958 445 3727 1130 1942 1276 990  432 1078  0");
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "1 -999 999 45.0 45.0 9 11 113.5 22.166667 126 128 -1 2 0.1 0.1  2119 1264 112 0 0 0 0 0");
    
        // PHI = 20 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 19 21 113.5 22.166667 126 128 -1 2 0.1 0.1 4396 0 4248 316 4025  215  3733"
        " 1019 2055 1371 1821 1027 1707  1098 1124  576  1212 153 1090  0");
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "1 -999 999 45.0 45.0 19 21 113.5 22.166667 126 128 -1 2 0.1 0.1 2283 1322 112 0 0 0 0 0");

        // PHI = 30 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 29 31 113.5 22.166667 126 128 -1 2 0.1 0.1 4396   0 4248 316 4025  215  3733"
        " 1019 2590 1251 2553 1438 2416 1502 2227 1333 2055 1371 1926 1180 1265  802 1270  313  1020 0");

        // PHI = 40 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 39 41 113.5 22.166667 126 128 -1 2 0.1 0.1 4358 0  4020 795  2442 1435 2132 1190 2011"
        " 1293 1402 923 1402 632 1307  567 1306 216 1175 172 1005 327 857 0");
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_21", 
        "1 -999 999 45.0 45.0 39 41 113.5 22.166667 126 128 -1 2 0.1 0.1 2614 1375 125 0 0 0 0 0");
        
        // PHI = 50 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 49 51 113.5 22.166667 126 128 -1 2 0.1 0.1 4358 0 4020 795 2927 1241 2923 1464 2767 1546 2522 1411"
        " 2454 1435 2314 1326 2023 1294 1412 1021 1412 440 1297 420 1150 561 1063 344 969 326 887 72 948 0");

        // PHI = 60 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 59 61 113.5 22.166667 126 128 -1 2 0.1 0.1 4541 0 3100 1433 2932 1557 2930 1503 2450 1322"
        " 2301 1396 1637 1210 1626 877 1494 783 1486 680 1230 766 1170 599 984 373 973 0");
    
        // PHI = 70 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 69 71 113.5 22.166667 126 128 -1 2 0.1 0.1 4531 0 3331 1178 3329 1410 3080 1504 1836 1344 1580 819"
        " 1342 927 882 0");
    
        // PHI = 80 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 79 81 113.5 22.166667 126 128 -1 2 0.1 0.1 4308 0 3430 1137 3441 1413 2983 1429 2804 1391 2669 1458"
        " 1993 1413 1730 1042 1726 938 1506 1073 1353 857 1049 186 1139 79 1140 0");

        // PHI = 90 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 89 91 113.5 22.166667 126 128 -1 2 0.1 0.1 4091 0 3496 1480 2135 1500 1860 1086 1667 1176 1440 959"
        " 1170 425 1240 308 997 0");
    
        // PHI = 100 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 99 101 113.5 22.166667 126 128 -1 2 0.1 0.1 4101 0 3726 874 3821 1134 3591 1296 3096 1307 3106 1465"
        "  2704 1481 2323 1552 2165 1291 1659 1280 1243 586 1243 457 842 0");

        // PHI = 110 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 109 111 113.5 22.166667 126 128 -1 2 0.1 0.1 4000 0 3796 887 3856 1088 3338 1182 3165 1364 2806 1441"
        "  2746 1528 2516 1532 2352 1303 2244 1298 2174 1203 2072 1211 1994 1333 1761 1224 1739 1111 1444 819 1437 647 900 254 882 0");

        // PHI = 120 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 119 121 113.5 22.166667 126 128 -1 2 0.1 0.1 4011 0 3976 998 3415 1126 3406 1260 2707 1561 2521 1335"
        "  2336 1330 2320 1274 2229 1279 2158 1431 1866 1302 1861 1167 1553 936 1594 827 1003 462 992 0");

        // PHI = 130 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 129 131 113.5 22.166667 126 128 -1 2 0.1 0.1 4106 0 4098 838 3520 1032 3547 1135 3169 1364 3162 1445"
        "  2909 1536 2466 1311 2343 1505 2005 1360 2008 1266 1677 1064 1698 977 1170 675 1092 0");

        // PHI = 140 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 139 141 113.5 22.166667 126 128 -1 2 0.1 0.1 4120 0 4116 202 4028 204 4018 307 4203 411 4218 678 3643"
        "  951 3697 1043 3075 1519 2575 1320 2541 1495 2176 1480 2166 1374 1822 1174 1831 1058 1300 847 1169 0");

        // PHI = 150 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 149 151 113.5 22.166667 126 128 -1 2 0.1 0.1 4126 0 4334 201 4340 485 3758 805 3849 965 3243 1460 2735"
        "  1332 2702 1515 2355 1500 2351 1390 1978 1255 1988 1154 1438 986 1215 175 1041 142 999 0");

        // PHI = 160 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 159 161 113.5 22.166667 126 128 -1 2 0.1 0.1 4486 0 4494 241 3871 666 3953 794 3479 1332 2913 1298"
        "  2882 1529 2527 1521 2527 1459 2122 1328 2122 1224 1550 1125 1350 438 1108 340 1083 0");

        // PHI = 170 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 169 171 113.5 22.166667 126 128 -1 2 0.1 0.1 4526 0 3999 496 4064 623 3688 1233 3067 1268 3059 1478"
        "  2791 1531 2714 1442 2284 1396 2284 1295 1703 1245 1454 602 1224 568 1193 0");

        // PHI = 180 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 179 181 113.5 22.166667 126 128 -1 2 0.1 0.1 4350 0 4096 289 4173 414 3841 1164 3315 1199 3292 1413"
        "  2987 1510 2463 1421 2448 1354 1878 1344 1581 759 1375 761 1183 0");

        // PHI = 190 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 189 191 113.5 22.166667 126 128 -1 2 0.1 0.1 4295 0 4305 352 3947 1079 3412 1105 3422 1434 2607 1445"
        "  2607 1360 2033 1430 1715 894 1461 899 1458 550 1028 0");

        // PHI = 200 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 199 201 113.5 22.166667 126 128 -1 2 0.1 0.1 4385 0 4063 883 3537 1011 3616 1328 2832 1441 2778 1362"
        "  2239 1478 1845 1010 1564 1027 1581 700 1130 224 1105 0");

        // PHI = 210 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 209 211 113.5 22.166667 126 128 -1 2 0.1 0.1 4491 0 4432 617 3646 881 3782 1211 3000 1414 2951 1313"
        " 2434 1535 1995 1135 1697 1128 1723 818 1298 485 1197 201 1211 0");

        // PHI = 220 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 219 221 113.5 22.166667 126 128 -1 2 0.1 0.1 4448 0 4321 519 3762 808 3907 1025 3688 1178 3167 1389"
        "  3089 1280 2617 1559 2112 1183 1981 1234 1832 1156 1851 882 1412 676 1283 412 1285 0");

        // PHI = 230 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 229 231 113.5 22.166667 126 128 -1 2 0.1 0.1 4528 0 4521 228 3865 624 4022 884 3824 1087 2778 1540"
        "  2233 1232 2171 1313 1942 1241 2005 1034 1534 842 1414 655 1350 0");

        // PHI = 240 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 239 241 113.5 22.166667 126 128 -1 2 0.1 0.1 4570 0 4116 265 4156 790 3008 1506 2401 1299 2271 1370"
        "  2083 1313 2072 1077 1641 957 1515 782 1377 0");

        // PHI = 250 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 249 251 113.5 22.166667 126 128 -1 2 0.1 0.1 4366 0 4352 630 3074 1571 1674 1238 1395 0");

        // PHI = 260 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 259 261 113.5 22.166667 126 128 -1 2 0.1 0.1 4254 0 4462 488 3736 1055 3294 1468 2358"
        "  1414 2354 1271 1777 1148 1342 0");

        // PHI = 270 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 269 271 113.5 22.166667 126 128 -1 2 0.1 0.1 4588 0 4070 875 3797 875 3531 1380 1934 1337 1336 0");
        
        // PHI = 280 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 279 281 113.5 22.166667 126 128 -1 2 0.1 0.1 4591 0 4375 307 4312 307 4041 732 3933 705"
        " 3665 1251 3485 1146 3214 1338 3065 1178 2145 1365 1201 0");

        // PHI = 290 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 289 291 113.5 22.166667 126 128 -1 2 0.1 0.1 4440 0 4184 553 4033 559 3869 1141 3301 1201"
        " 3078 1431 2822 1242 2318 1425 1240 301 970 372 817 0");

        // PHI = 300 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 299 301 113.5 22.166667 126 128 -1 2 0.1 0.1 4451 0 4253 347 4166 328 4031 996 3320 1127"
        " 3388 1335 3071 1278 2375 1484 1360 482 1065 594 1003 0");

        // PHI = 310 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 309 311 113.5 22.166667 126 128 -1 2 0.1 0.1 4255 0 4081 887 3541 1044 3530 971 2512"
        " 1453 1483 669 1234 777 1142 652 1081 0");

        // PHI = 320 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 319 321 113.5 22.166667 126 128 -1 2 0.1 0.1 4331 0 4234 748 3611 978 3602 886 2823"
        " 1400 2037 1070 1260 861 1128 0");

        // PHI = 330 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 329 331 113.5 22.166667 126 128 -1 2 0.1 0.1 4428 0 4359 561 3775 842 3811 955 3617"
        " 1135 3436 965 3325 1282 2897 1417 1730 961 1466 1092 1259 318 1108 339 1010 0");

        // PHI = 340 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 339 341 113.5 22.166667 126 128 -1 2 0.1 0.1 4522 0 4490 254 3894 739 3902 971 3662"
        " 1047 3631 963 3354 1269 2942 1386 1840 1085 1580 1191 1363 509 1200 522 989 0");

        // PHI = 350 
        (void) poHeader->nReplaceValue("CRYSTAL_GONIO_SHADOW_" + Cstring(++ii),
        "0 -999 999 45.0 45.0 349 351 113.5 22.166667 126 128 -1 2 0.1 0.1 4529 0 3994 550 4027 793 3813 895 3726 819 3337 1318"
        "  1778 1282 1511 874 1489 703 1337 703 1086 214 987 173 878 0");
    }
    
    (void)poHeader->nReplaceValue("CRYSTAL_GONIO_NUM_SHADOWS", ii);
}

#ifdef RAXIS_HEADER_BINARY_STRESS_INFO
bool bGetSetStressMeasurementInformation(RAXIS_header& stRAXIS_Header, Cimage_header& oHeader, bool bSwapBytes)
{
    if( bSwapBytes )
    {
        stRAXIS_Header.posn  = nSwapLong(stRAXIS_Header.posn);
        stRAXIS_Header.side  = nSwapLong(stRAXIS_Header.side);
        stRAXIS_Header.slit  = nSwapLong(stRAXIS_Header.slit);
    }

    if( stRAXIS_Header.posn <= 0 )
        return false;

    Cstring     sTemp(""); 
    sTemp.assign(stRAXIS_Header.posz, 64);
    vTrimTrailingSpaces(sTemp);

    std::vector<int>     anPositions;
    sTemp.nListToVector(anPositions, " ");

    Cstress     oStress;
    oStress.vSetZPositions(anPositions.size(), &anPositions[0]);
    
    oStress.vSetNumOfZpositions(stRAXIS_Header.posn);

    sTemp.assign(stRAXIS_Header.psis, 128);
    vTrimTrailingSpaces(sTemp);

    std::vector<double>     afPositions;
    sTemp.nListToVector(afPositions, " ");

    oStress.vSetPsiAngles(afPositions.size(), &afPositions[0]);

    if( stRAXIS_Header.side != 0 && stRAXIS_Header.side != 1 )
        return false;

    oStress.vSetInclinationType(stRAXIS_Header.side == 0 ? Cstring("Iso") : Cstring("Side"));
 
    if( stRAXIS_Header.slit <= 0 )
        return false;

    oStress.vSetSlitWidth(stRAXIS_Header.slit);

    oStress.nUpdateHeader(oHeader);

    return true;
}
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
