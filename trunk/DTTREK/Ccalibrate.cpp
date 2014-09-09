//
// Copyright (c) 2006 Rigaku
//
// Ccalibrate.cc        Initial author: J.W. Pflugrath           02-June-1998
//  This file contains the member functions of class Ccalibrate which implements
//    the calibration encapsulation of d*TREK.  This class was built
//    from Dr. Marty Stanton's (Brandeis University) non-C++ calibrate code 
//    of 1994 that was modified by JWP.
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
//
//+ToDo
//
//   Error messages need implementing
//
//+Include files

#ifdef ANL_TIMER
#include "anlTimer.h"
#endif

#include "Ccalibrate.h"         // Class definition and prototypes
#include "Cspatial.h"
#include "Cscan.h"

using namespace std;

//+Definitions, constants, and initialization of static member variables

Cstring Ccalibrate::ms_sDtcalibFiles   = D_K_DtcalibFiles;
Cstring Ccalibrate::ms_sDtcalibOptions = D_K_DtcalibOptions;
Cstring Ccalibrate::ms_asExt[] = { ".inv_x_int",
                                   ".inv_y_int",
                                   ".x_int",
                                   ".y_int",
                                   ".calpar" };

//+Code begin


int nInterpolateMethod(int nPix0,int nPix1,int nOutDim0,int nOutDim1,int nBadValue,Cimage* poImgBadInterp,Cimage* poImageOut);
int nInterpolationFunction(int nOutDim0,int nOutDim1,int nBadValue,int nNumBadInterpOut,Cimage* poImgBadInterp,Cimage* poImageOut);

//+Public functions

// Constructors, destructors and assignments

Ccalibrate::Ccalibrate(Cimage_header& oHeader)
{
  // Careful with the following NULL's, they are also set in ::nInitValues();

  m_poNonunf         = NULL;
  m_poImgNonunf      = NULL;
  m_poImgMark        = NULL;
  m_poImgDistorXINT  = NULL;
  m_poImgDistorXINV  = NULL;
  m_poImgDistorYINT  = NULL;
  m_poImgDistorYINV  = NULL;
  m_poImgMask        = NULL;
  m_poImgRefer       = NULL;
  m_poImgDark        = NULL;
  m_poImgFlood       = NULL;
  m_poImgBadPix      = NULL;
  m_poImgBadInterp   = NULL;
  m_poHeaderSave     = NULL;
  m_poImgTransform   = NULL;
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  m_poXprop          = NULL;
#endif
  m_poReflnlistMask  = NULL;
  m_poOutDump        = NULL;

  m_poImgMergeDistorXINT  = NULL;
  m_poImgMergeDistorXINV  = NULL;
  m_poImgMergeDistorYINT  = NULL;
  m_poImgMergeDistorYINV  = NULL;
  m_poHeadMergeDistor     = NULL;

 (void) nInitValues(oHeader);  
}

Ccalibrate::~Ccalibrate() 
{
  if (NULL != m_poOutDump)
    {
      m_poOutDump->close();
      delete m_poOutDump;
      m_poOutDump = NULL;
    }
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  if (NULL != m_poXprop)
    {
      delete m_poXprop;
      m_poXprop = NULL;
    }
#endif
  if (NULL != m_poReflnlistMask)
    {
      delete m_poReflnlistMask;
      m_poReflnlistMask = NULL;
    }
  if (NULL != m_poNonunf)
    {
      delete m_poNonunf;
      m_poNonunf = NULL;
    }
  if (NULL != m_poImgNonunf)
    {
      delete m_poImgNonunf;
      m_poImgNonunf = NULL;
    }
  if (NULL != m_poImgBadPix)
    {
      delete m_poImgBadPix;
      m_poImgBadPix = NULL;
    }
  if (NULL != m_poImgMark)
    {
      delete m_poImgMark;
      m_poImgMark = NULL;
    }
  if (NULL != m_poImgDistorXINT)
    {
      delete m_poImgDistorXINT;
      m_poImgDistorXINT = NULL;
    }
  if (NULL != m_poImgDistorXINV)
    {
      delete m_poImgDistorXINV;
      m_poImgDistorXINV = NULL;
    }
  if (NULL != m_poImgDistorYINT)
    {
      delete m_poImgDistorYINT;
      m_poImgDistorYINT = NULL;
    }
  if (NULL != m_poImgDistorYINV)
    {
      delete m_poImgDistorYINV;
      m_poImgDistorYINV = NULL;
    }
  if (NULL != m_poImgMask)
    {
      delete m_poImgMask;
      m_poImgMask = NULL;
    }

  if (NULL != m_poImgRefer)
    {
      delete m_poImgRefer;
      m_poImgRefer = NULL;
    }

  if (NULL != m_poImgDark)
    {
      delete m_poImgDark;
      m_poImgDark = NULL;
    }

  if (NULL != m_poImgTransform)
    {
      delete m_poImgTransform;
      m_poImgTransform = NULL;
    }

  if (NULL != m_poImgFlood)
    {
      delete m_poImgFlood;
      m_poImgFlood = NULL;
    }

  if (NULL != m_poHeaderSave)
    {
      delete m_poHeaderSave;
      m_poHeaderSave = NULL;
    }
  if (NULL != m_poImgMergeDistorXINT)
    {
      delete m_poImgMergeDistorXINT;
      m_poImgMergeDistorXINT = NULL;
    }
  if (NULL != m_poImgMergeDistorXINV)
    {
      delete m_poImgMergeDistorXINV;
      m_poImgMergeDistorXINV = NULL;
    }
  if (NULL != m_poImgMergeDistorYINT)
    {
      delete m_poImgMergeDistorYINT;
      m_poImgMergeDistorYINT = NULL;
    }
  if (NULL != m_poImgMergeDistorYINV)
    {
      delete m_poImgMergeDistorYINV;
      m_poImgMergeDistorYINV = NULL;
    }
  if (NULL != m_poHeadMergeDistor)
    {
      delete m_poHeadMergeDistor;
      m_poHeadMergeDistor = NULL;
    }

  if (NULL != m_poImgBadInterp)
    {
      delete m_poImgBadInterp;
      m_poImgBadInterp = NULL;
      m_nNumBadInterpOut = 0;
    }
}

int Ccalibrate::nInitValues(const int nDim0, const int nDim1)
{
  // Initialize starting values.  These are usually overwritten
  // from values from the input header/save file.

  // We may want a different mechanism for getting the image dimensions

  m_a2nDim[0]        = nDim0;
  m_a2nDim[1]        = nDim1;
  m_a2nTOutDim[0]    = 0;
  m_a2nTOutDim[1]    = 0;
  m_a2nTOutOrig[0]   = 0;
  m_a2nTOutOrig[1]   = 0;
  m_a2nNumModules[0] = 1;
  m_a2nNumModules[1] = 1;

  m_nNumBadInterpOut = 0;
  m_poImgBadInterp   = NULL;
  m_nVerbose         = 2;
  m_nCorrectionMask  = CAL_DO_DARK | CAL_DO_NONUNF | CAL_DO_SPATIAL;

  m_sOutHeader       = sDtrekGetPrefix() + "dtcalibrate.head";
  m_sNonunfName      = sTransSymbol("$(NONUNF)");
  m_sDistorName      = sTransSymbol("$(DISTOR)");
  m_sMaskName        = sTransSymbol("$(MASK)");
  m_sReferName       = sTransSymbol("$(REFER)");
  m_sDarkName        = sTransSymbol("$(DARK)");
  m_sFloodName       = sTransSymbol("$(FLOOD)");
  m_sTransformName   = sTransSymbol("$(TRANSFORM)");
  m_sPostName        = "";
  m_sMergeInfoName   = sTransSymbol("$(MERGEINFO)");
  m_sSaveName        = m_sOutHeader;
  m_sBadpixName      = sTransSymbol("$(BADPIX)");
  m_sTBadpixName     = sTransSymbol("$(TBADPIX)");
  m_sDumpName        = sTransSymbol("$(DUMP)");
  m_sMarkName        = sTransSymbol("$(MARK)");
  m_sMaskPeaksName   = "";
  m_fMaskPeaksDivisor = 1.0;

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  m_poXprop          = NULL;
#endif
  m_poNonunf         = NULL;
  m_poImgNonunf      = NULL;
  m_poImgMark        = NULL;
  m_poImgDistorXINT  = NULL;
  m_poImgDistorXINV  = NULL;
  m_poImgDistorYINT  = NULL;
  m_poImgDistorYINV  = NULL;
  m_poImgMask        = NULL;
  m_poImgRefer       = NULL;
  m_poImgDark        = NULL;
  m_poImgFlood       = NULL;
  m_poImgBadPix      = NULL;
  m_poHeaderSave     = NULL;
  m_poOutDump        = NULL;

  m_poImgMergeDistorXINT  = NULL;
  m_poImgMergeDistorXINV  = NULL;
  m_poImgMergeDistorYINT  = NULL;
  m_poImgMergeDistorYINV  = NULL;
  m_poHeadMergeDistor     = NULL;

  m_nflood_radial           = 0;
  m_nflood_interpolate      = 1;
  m_nflood_geometry         = 1;
  m_nflood_film             = 0;
  m_fdarksub_const          = 0.0;
  m_fdarksub_scale          = 1.0;
  m_a2fpixel_sd[DARKIMG]    = 200.0;
  m_a2fpixel_sd[FLOODIMG]   = 200.0;
  m_a3nmax_pixel[DARKIMG]   = 65535;
  m_a3nmax_pixel[FLOODIMG]  = 65535;
  m_a3nmax_pixel[FLOOD_DARKIMG]= 65535;
  m_a3nmin_pixel[DARKIMG]   = 0;
  m_a3nmin_pixel[FLOODIMG]  = 0;
  m_a3nmin_pixel[FLOOD_DARKIMG]= 0;
  m_nradial_distortion      = 0;
  m_nbad_peaks              = 2;
  m_fmask_angle             = 0.0;
  m_nmin_peak               = 30000;
  m_fcent_to_maxval         = 2.0;
  m_fcent_to_pred           = 3.0;
  m_nbackground_pixels      = 25;
  m_npeak_size              = 9;
  m_fpeak_distance          = 14.0;
  m_a2ncenter_peak[0]       = m_a2nDim[0] / 2;
  m_a2ncenter_peak[1]       = m_a2nDim[1] / 2;
  m_fsearch_radius          = sqrtf( (float)
                                    (m_a2ncenter_peak[0] * m_a2ncenter_peak[0]
                                   + m_a2ncenter_peak[1] + m_a2ncenter_peak[1])
                                    );
  m_a2nsearch_center[0]     = m_a2nDim[0] / 2;
  m_a2nsearch_center[1]     = m_a2nDim[1] / 2;
  m_a2nvertical_limits[0]   = 0;
  m_a2nvertical_limits[1]   = m_a2nDim[1] - 1;
  m_a2nhorizontal_limits[0] = 0;
  m_a2nhorizontal_limits[1] = m_a2nDim[0] - 1;  // Fast direction is horizontal
  m_fxtod_distance          = 100.0;
  m_fmasktod_distance       = 0.0;
  m_a2fbeam_position[0]     = (float)m_a2nDim[0] * 0.5;
  m_a2fbeam_position[1]     = (float)m_a2nDim[1] * 0.5;
  m_fmask_spacing           = 1.0;
  m_fpixel_size             = 0.0;
  m_ndump                   = 0;
  m_nedges                  = 0;
  m_bCreateTransform        = FALSE;

  m_fReferenceMax           = 20000.0;
  m_fReferenceDistortMaxFactor = 1.15;

  m_a2nSeqNum[0]             = 0;
  m_a2nSeqNum[1]             = 1000;
//  m_tDistorParams;

  m_tCorrectParams.x_center      = 0.0;
  m_tCorrectParams.y_center      = 0.0;
  m_tCorrectParams.x_pt_center   = 0.0;
  m_tCorrectParams.y_pt_center   = 0.0;
  m_tCorrectParams.x_scale       = 1.0;
  m_tCorrectParams.y_scale       = 1.0;
  m_tCorrectParams.ratio         = 1.0;
  m_tCorrectParams.ver_slope     = 0.0;
  m_tCorrectParams.horz_slope    = 0.0;
  m_tCorrectParams.xint_step     = 2;
  m_tCorrectParams.yint_step     = 2;
  m_tCorrectParams.xinv_step     = 2;
  m_tCorrectParams.yinv_step     = 2;
  m_tCorrectParams.xint_start    = 1;
  m_tCorrectParams.yint_start    = 1;
  m_tCorrectParams.xinv_start    = 1;
  m_tCorrectParams.yinv_start    = 1;
  m_tCorrectParams.badflag       = INTFLAG;
  m_tCorrectParams.pscale        = 1.0 / 512.0;

  // There is some confusion about whether to use other member
  // variables or m_tMaskParams (e.g. m_fpixel_size vs m_tMaskParams.pixel_size)

  m_tMaskParams.away             = m_fcent_to_maxval;
  m_tMaskParams.back             = m_nbackground_pixels;
  m_tMaskParams.beam[0]          = m_a2fbeam_position[0];
  m_tMaskParams.beam[1]          = m_a2fbeam_position[1];
  m_tMaskParams.center[0]        = m_a2ncenter_peak[0];
  m_tMaskParams.center[1]        = m_a2ncenter_peak[1];
  m_tMaskParams.dis              = m_fpeak_distance;
  m_tMaskParams.lim              = m_nmin_peak;
  m_tMaskParams.maskd            = 0;
  m_tMaskParams.mask_angle       = m_fmask_angle;
  m_tMaskParams.dump             = m_ndump;
  m_tMaskParams.misses           = m_nbad_peaks;
  m_tMaskParams.move             = m_fcent_to_pred;
  m_tMaskParams.pkrat            = 0.0;
  m_tMaskParams.pkwid            = m_npeak_size;
  m_tMaskParams.radial           = m_nradial_distortion;
  m_tMaskParams.search_center[0] = m_a2nsearch_center[0];
  m_tMaskParams.search_center[1] = m_a2nsearch_center[1];
  m_tMaskParams.resrad           = (int) m_fsearch_radius;
  m_tMaskParams.resrad2          = (int)(m_fsearch_radius*m_fsearch_radius);
  m_tMaskParams.vertical[0]      = m_a2nvertical_limits[0];
  m_tMaskParams.vertical[1]      = m_a2nvertical_limits[1];
  m_tMaskParams.source           = m_fxtod_distance;
  m_tMaskParams.horizontal[0]    = m_a2nhorizontal_limits[0];
  m_tMaskParams.horizontal[1]    = m_a2nhorizontal_limits[1];
  m_tMaskParams.work             = 1;
  m_tMaskParams.maskparams[0]    = 90.0;       // Angle between mask holes
  m_tMaskParams.maskparams[1]    = 1.0;        // Ratio of ?
  m_tMaskParams.pixel_size       = m_fpixel_size;
  m_tMaskParams.mask_spacing     = m_fmask_spacing;

  return (0);
}

int Ccalibrate::nInitValues(Cimage_header& oHeader)
{
  int    i;       // Loop counter
  int     nStat, nTemp;
  int     nStatL;
  float   fTemp;
  float   a3fTemp[3];
  int     a3nTemp[3];
  Cstring sTemp;
  Cstring sOutString = "";
  static Cstring ssEndl = "\n";
  static Cstring ssWarning 
           = "WARNING in Ccalibrate::nInitValues(), problem reading keyword: ";

  (void) nInitValues();

  if (!oHeader.bIsAvailable())
    {
      cout << "Ccalibrate ERROR: image header is not valid."
           << "  Cannot construct refine!\n";
      nStat = 1;
      return (nStat);
    }

  // Try to get required information and parameters from the header

  nStat  = 0;

  nStatL = oHeader.nGetValue(D_K_Calib_number_modules, 2, a3nTemp);
  if (0 == nStatL)
    {
      m_a2nNumModules[0] = a3nTemp[0];
      m_a2nNumModules[1] = a3nTemp[1];
    }
  else
    {
      m_a2nNumModules[0] = 1;
      m_a2nNumModules[1] = 1;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_flood_radial, &nTemp);
  if (0 == nStatL)
    {
      m_nflood_radial = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_flood_radial + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_flood_interpolate, &nTemp);
  if (0 == nStatL)
    {
      m_nflood_interpolate = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_flood_interpolate + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_flood_geometry, &nTemp);
  if (0 == nStatL)
    {
      m_nflood_geometry = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_flood_geometry + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_flood_film, &nTemp);
  if (0 == nStatL)
    {
      m_nflood_film = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_flood_film + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_darksub_const, &fTemp);
  if (0 == nStatL)
    {
      m_fdarksub_const = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_darksub_const + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_darksub_scale, &fTemp);
  if (0 == nStatL)
    {
      m_fdarksub_scale = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_darksub_scale + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_pixel_sd, 2, a3fTemp);
  if (0 == nStatL)
    {
      m_a2fpixel_sd[0] = a3fTemp[0];
      m_a2fpixel_sd[1] = a3fTemp[1];
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_pixel_sd + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_max_pixel, 3, a3nTemp);
  if (0 == nStatL)
    {
      m_a3nmax_pixel[0] = a3nTemp[0];
      m_a3nmax_pixel[1] = a3nTemp[1];
      m_a3nmax_pixel[2] = a3nTemp[2];
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_max_pixel + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_min_pixel, 3, a3nTemp);
  if (0 == nStatL)
    {
      m_a3nmin_pixel[0] = a3nTemp[0];
      m_a3nmin_pixel[1] = a3nTemp[1];
      m_a3nmin_pixel[2] = a3nTemp[2];
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_min_pixel + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_radial_distortion, &nTemp);
  if (0 == nStatL)
    {
      m_nradial_distortion = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_radial_distortion + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_bad_peaks, &nTemp);
  if (0 == nStatL)
    {
      m_nbad_peaks = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_bad_peaks + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_mask_angle, &fTemp);
  if (0 == nStatL)
    {
      m_fmask_angle = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_mask_angle + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_min_peak, &nTemp);
  if (0 == nStatL)
    {
      m_nmin_peak = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_min_peak + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_cent_to_maxval, &fTemp);
  if (0 == nStatL)
    {
      m_fcent_to_maxval = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_cent_to_maxval + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_cent_to_pred, &fTemp);
  if (0 == nStatL)
    {
      m_fcent_to_pred = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_cent_to_pred + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_background_pixels, &nTemp);
  if (0 == nStatL)
    {
      m_nbackground_pixels = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_background_pixels + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_peak_size, &nTemp);
  if (0 == nStatL)
    {
      m_npeak_size = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_peak_size + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_peak_distance, &fTemp);
  if (0 == nStatL)
    {
      m_fpeak_distance = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_peak_distance + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_center_peak, 2, a3nTemp);
  if (0 == nStatL)
    {
      m_a2ncenter_peak[0] = a3nTemp[0];
      m_a2ncenter_peak[1] = a3nTemp[1];
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_center_peak + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_search_radius, &fTemp);
  if (0 == nStatL)
    {
      m_fsearch_radius = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_search_radius + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_search_center, 2, a3nTemp);
  if (0 == nStatL)
    {
      m_a2nsearch_center[0] = a3nTemp[0];
      m_a2nsearch_center[1] = a3nTemp[1];
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_search_center + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_vertical_limits, 2, a3nTemp);
  if (0 == nStatL)
    {
      m_a2nvertical_limits[0] = a3nTemp[0];
      m_a2nvertical_limits[1] = a3nTemp[1];
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_vertical_limits + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_horizontal_limits, 2, a3nTemp);
  if (0 == nStatL)
    {
      m_a2nhorizontal_limits[0] = a3nTemp[0];
      m_a2nhorizontal_limits[1] = a3nTemp[1];
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_horizontal_limits + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_xtod_distance, &fTemp);
  if (0 == nStatL)
    {
      m_fxtod_distance = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_xtod_distance + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_masktod_distance, &fTemp);
  if (0 == nStatL)
    {
      m_fmasktod_distance = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_masktod_distance + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_beam_position, 2, a3fTemp);
  if (0 == nStatL)
    {
      m_a2fbeam_position[0] = a3fTemp[0];
      m_a2fbeam_position[1] = a3fTemp[1];
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_beam_position + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_mask_spacing, &fTemp);
  if (0 == nStatL)
    {
      m_fmask_spacing = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_mask_spacing + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_pixel_size, &fTemp);
  if (0 == nStatL)
    {
      m_fpixel_size = fTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_pixel_size + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_edges, &nTemp);
  if (0 == nStatL)
    {
      m_nedges = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_edges + ssEndl);
      nStat++;
    }

  nStatL = oHeader.nGetValue(D_K_Calib_dump_mode, &nTemp);
  if (0 == nStatL)
    {
      m_ndump = nTemp;
    }
  else
    {
      sOutString +=  (ssWarning + D_K_Calib_dump_mode + ssEndl);
      nStat++;
    }

  if ("" != sOutString)
    {
#ifndef CORRECT_ONLY      
      cout << sOutString << flush;
#endif
    }

  // The following, if absent, do not cause a return error

  nStatL = oHeader.nGetValue(ms_sDtcalibFiles, &m_sCalibFiles);
  nStatL = oHeader.nGetValue(ms_sDtcalibOptions, &m_sCalibOptions);


  // There is some confusion about whether to use other member
  // variables or m_tMaskParams (e.g. m_fpixel_size vs m_tMaskParams.pixel_size)

  m_tMaskParams.away             = m_fcent_to_maxval;
  m_tMaskParams.back             = m_nbackground_pixels;
  m_tMaskParams.beam[0]          = m_a2fbeam_position[0];
  m_tMaskParams.beam[1]          = m_a2fbeam_position[1];
  m_tMaskParams.center[0]        = m_a2ncenter_peak[0];
  m_tMaskParams.center[1]        = m_a2ncenter_peak[1];
  m_tMaskParams.dis              = m_fpeak_distance;
  m_tMaskParams.lim              = m_nmin_peak;
  m_tMaskParams.maskd            = 0;
  m_tMaskParams.mask_angle       = m_fmask_angle;
  m_tMaskParams.dump             = m_ndump;
  m_tMaskParams.misses           = m_nbad_peaks;
  m_tMaskParams.move             = m_fcent_to_pred;
  m_tMaskParams.pkrat            = 0.0;
  m_tMaskParams.pkwid            = m_npeak_size;
  m_tMaskParams.radial           = m_nradial_distortion;
  m_tMaskParams.search_center[0] = m_a2nsearch_center[0];
  m_tMaskParams.search_center[1] = m_a2nsearch_center[1];
  m_tMaskParams.resrad           = (int) m_fsearch_radius;
  m_tMaskParams.resrad2          = (int)(m_fsearch_radius*m_fsearch_radius);
  m_tMaskParams.vertical[0]      = m_a2nvertical_limits[0];
  m_tMaskParams.vertical[1]      = m_a2nvertical_limits[1];
  m_tMaskParams.source           = m_fxtod_distance;
  m_tMaskParams.horizontal[0]    = m_a2nhorizontal_limits[0];
  m_tMaskParams.horizontal[1]    = m_a2nhorizontal_limits[1];
  m_tMaskParams.work             = 1;
  m_tMaskParams.maskparams[0]    = 90.0;       // Angle between mask holes
  m_tMaskParams.maskparams[1]    = 1.0;        // Ratio of ?
  m_tMaskParams.pixel_size       = m_fpixel_size;
  m_tMaskParams.mask_spacing     = m_fmask_spacing;

  return (nStat);
}

int Ccalibrate::nList(const int nFlag) 
{
  Cstring  a2sYesNo[2];
  a2sYesNo[0] = "No";
  a2sYesNo[1] = "Yes";
  cout << "Calibrate listing:\n";
  cout 
//       << "\nImage size (pixels):                  " << m_a2nDim[0] << ", " << m_a2nDim[1]
       << "\nNumber of modules:                    " << m_a2nNumModules[0] << ", " 
                                                     << m_a2nNumModules[1]
       << "\nFlood radial:                         " << a2sYesNo[m_nflood_radial]
       << "\nFlood interpolate:                    " << a2sYesNo[m_nflood_interpolate]
//Not used
//       << "\nFlood film:                           " << a2sYesNo[m_nflood_film]
       << "\nRadial distortion:                    " << a2sYesNo[m_nradial_distortion]
       << "\nDark subtract scale factor:           " << m_fdarksub_scale
       << "\nDark subtract constant:               " << m_fdarksub_const
       << "\nBad pixel dark, flood image min:      " << m_a3nmin_pixel[0]
                                             << "  " << m_a3nmin_pixel[1]
       << "\nBad pixel dark, flood image max:      " << m_a3nmax_pixel[0]
                                             << "  " << m_a3nmax_pixel[1]
       << "\nBad pixel dark, flood image SD:       " << m_a2fpixel_sd[0]
                                             << "  " << m_a2fpixel_sd[1]
       << "\nBad pixel flood-dark min:             " << m_a3nmin_pixel[2]
       << "\nMax bad mask peaks in search:         " << m_nbad_peaks
       << "\nMask angle:                           " << m_fmask_angle
       << "\nMin intensity of a mask peak:         " << m_nmin_peak
       << "\nMax allowed diff |cent-max|:          " << m_fcent_to_maxval
       << "\nMax allowed diff {cent-pred|:         " << m_fcent_to_pred
       << "\nMin required background pixels:       " << m_nbackground_pixels
           << " (obsolete) "
       << "\nBox size for peak search:             " << m_npeak_size
       << "\nStarting distance between peaks (px): " << m_fpeak_distance
       << "\nStarting peak position (fast, slow):  " << m_a2ncenter_peak[0]
                            << ", " << m_a2ncenter_peak[1]
       << "\nMask search limit circle center:      " << m_a2nsearch_center[0]
             << ", " << m_a2nsearch_center[1]
       << "\nMask search limit circle radius (px): " << m_fsearch_radius
       << "\nMask search limit rect   fast   (px): " << m_a2nhorizontal_limits[0]
            << ", " << m_a2nhorizontal_limits[1]
       << "\nMask search limit rect   slow   (px): " << m_a2nvertical_limits[0]
            << ", " << m_a2nvertical_limits[1]
       << "\nFlood source to detector distance:    " << m_fxtod_distance
       << "\nMask to detector distance :           " << m_fmasktod_distance
       << "\nDirect beam position (px) :           " << m_a2fbeam_position[0]
          << ", " << m_a2fbeam_position[1]
       << "\nMask spacing (mm):                    " << m_fmask_spacing
       << "\nPixel size (mm):                      " << setw(7) << m_fpixel_size
       << "\nCreate/write dump file:               " << a2sYesNo[m_ndump]
       << "\nTreat edges:                          " << a2sYesNo[m_nedges]
       << "\n"
       << "\nMask  image:                          " << m_sMaskName
       << "\nFlood image:                          " << m_sFloodName
       << "\nDark  image:                          " << m_sDarkName
       << "\nReference image:                      " << m_sReferName
       << "\nDistortion tables basename:           " << m_sDistorName
       << "\nNonunf image:                         " << m_sNonunfName
       << "\nBad pixel name:                       " << m_sBadpixName
       << "\nOutput bad pixel name:                " << m_sTBadpixName
       << "\nTransform image:                      " << m_sTransformName
       << "\nDump file:                            " << m_sDumpName
       << "\nMark file:                            " << m_sMarkName
       << '\n' << endl;

  return (0);
}

int Ccalibrate::nSetup(void) 
{
  // Setup before refinement loops

  return (0);
}
int
Ccalibrate::nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre)
{
  // Place the current values of the calibrate 
  // properties in the given header

  if ("" != sPre)
    {
      cout << "WARNING! prefix not used in Ccalibrate::nUpdateHeader()" << endl;
    }
  else if (NULL == poHeader)
    {
      cout << "ERROR, no header in Ccalibrate::nUpdateHeader()" << endl;
      return (-1);
    }

  (void) poHeader->nReplaceValue(D_K_Calib_number_modules, 2, m_a2nNumModules);
  (void) poHeader->nReplaceValue(D_K_Calib_flood_radial, m_nflood_radial);
  (void) poHeader->nReplaceValue(D_K_Calib_flood_interpolate,
                               m_nflood_interpolate);
  (void) poHeader->nReplaceValue(D_K_Calib_flood_geometry, m_nflood_geometry);
  (void) poHeader->nReplaceValue(D_K_Calib_flood_film, m_nflood_film);
  (void) poHeader->nReplaceValue(D_K_Calib_darksub_const, m_fdarksub_const);
  (void) poHeader->nReplaceValue(D_K_Calib_darksub_scale, m_fdarksub_scale);
  (void) poHeader->nReplaceValue(D_K_Calib_pixel_sd, 2, m_a2fpixel_sd);
  (void) poHeader->nReplaceValue(D_K_Calib_max_pixel, 3, m_a3nmax_pixel);
  (void) poHeader->nReplaceValue(D_K_Calib_min_pixel, 3, m_a3nmin_pixel);
  (void) poHeader->nReplaceValue(D_K_Calib_radial_distortion,
                               m_nradial_distortion);
  (void) poHeader->nReplaceValue(D_K_Calib_bad_peaks, m_nbad_peaks);
  (void) poHeader->nReplaceValue(D_K_Calib_mask_angle, m_fmask_angle);
  (void) poHeader->nReplaceValue(D_K_Calib_min_peak, m_nmin_peak);
  (void) poHeader->nReplaceValue(D_K_Calib_cent_to_maxval, m_fcent_to_maxval);
  (void) poHeader->nReplaceValue(D_K_Calib_cent_to_pred, m_fcent_to_pred);
  (void) poHeader->nReplaceValue(D_K_Calib_background_pixels,
                               m_nbackground_pixels);
  (void) poHeader->nReplaceValue(D_K_Calib_peak_size, m_npeak_size);
  (void) poHeader->nReplaceValue(D_K_Calib_peak_distance, m_fpeak_distance);
  (void) poHeader->nReplaceValue(D_K_Calib_center_peak, 2, m_a2ncenter_peak);
  (void) poHeader->nReplaceValue(D_K_Calib_search_radius, m_fsearch_radius);
  (void) poHeader->nReplaceValue(D_K_Calib_search_center, 2, m_a2nsearch_center);
  (void) poHeader->nReplaceValue(D_K_Calib_vertical_limits, 2, 
                               m_a2nvertical_limits);
  (void) poHeader->nReplaceValue(D_K_Calib_horizontal_limits, 2, 
                               m_a2nhorizontal_limits);
  (void) poHeader->nReplaceValue(D_K_Calib_xtod_distance, m_fxtod_distance);
  (void) poHeader->nReplaceValue(D_K_Calib_masktod_distance, m_fmasktod_distance);
  (void) poHeader->nReplaceValue(D_K_Calib_beam_position, 2, m_a2fbeam_position);
  (void) poHeader->nReplaceValue(D_K_Calib_mask_spacing, m_fmask_spacing);
  (void) poHeader->nReplaceValue(D_K_Calib_pixel_size, m_fpixel_size);
  (void) poHeader->nReplaceValue(D_K_Calib_edges, m_nedges);
  (void) poHeader->nReplaceValue(D_K_Calib_dump_mode, m_ndump);
  (void) poHeader->nReplaceValue(ms_sDtcalibFiles, m_sCalibFiles);
  (void) poHeader->nReplaceValue(ms_sDtcalibOptions, m_sCalibOptions);
  
  // Also write the normal .calpar stuff

  if (0.0 >= m_fpixel_size)
    {
      if (0.0 != m_tCorrectParams.x_scale)
        m_fpixel_size = m_tDistorParams.spacing / m_tCorrectParams.x_scale;
      cout << "Pixel size is 0.0, so try to calculate it... set to: " 
           << setw(9) << m_fpixel_size << endl;
    }
  (void) poHeader->nReplaceValue(D_K_CalibPixelSize, 
                                 m_fpixel_size, 9);

  (void) poHeader->nReplaceValue(D_K_CalibPscale, 
                                 m_tCorrectParams.pscale, 9);
  (void) poHeader->nGetValue(D_K_CalibPscale,
                             &m_tCorrectParams.pscale);

  (void) poHeader->nReplaceValue(D_K_CalibXintStart,
                                 m_tCorrectParams.xint_start);
  (void) poHeader->nReplaceValue(D_K_CalibYintStart,
                                 m_tCorrectParams.yint_start);
  (void) poHeader->nReplaceValue(D_K_CalibXinvStart,
                                 m_tCorrectParams.xinv_start);
  (void) poHeader->nReplaceValue(D_K_CalibYinvStart,
                                 m_tCorrectParams.yinv_start);
  (void) poHeader->nReplaceValue(D_K_CalibXintStep,
                                 m_tCorrectParams.xint_step);
  (void) poHeader->nReplaceValue(D_K_CalibYintStep,
                                 m_tCorrectParams.yint_step);
  (void) poHeader->nReplaceValue(D_K_CalibXinvStep,
                                 m_tCorrectParams.xinv_step);
  (void) poHeader->nReplaceValue(D_K_CalibYinvStep,
                                 m_tCorrectParams.yinv_step);
  (void) poHeader->nReplaceValue(D_K_CalibXSize,
                                 m_tDistorParams.x_size);
  (void) poHeader->nReplaceValue(D_K_CalibYSize,
                                 m_tDistorParams.y_size);
  (void) poHeader->nReplaceValue(D_K_CalibXBeam,
                                 m_tDistorParams.x_beam);
  (void) poHeader->nReplaceValue(D_K_CalibYBeam,
                                 m_tDistorParams.y_beam);
  (void) poHeader->nReplaceValue(D_K_CalibBadFlag,
                                 m_tCorrectParams.badflag);
  (void) poHeader->nReplaceValue(D_K_CalibXCenter,
                                 m_tCorrectParams.x_center, 9);
  (void) poHeader->nReplaceValue(D_K_CalibYCenter,
                                 m_tCorrectParams.y_center, 9);
  (void) poHeader->nReplaceValue(D_K_CalibXPtCenter,
                                 m_tCorrectParams.x_pt_center, 9);
  (void) poHeader->nReplaceValue(D_K_CalibYPtCenter,
                                 m_tCorrectParams.y_pt_center, 9);
  (void) poHeader->nReplaceValue(D_K_CalibXScale,
                                 m_tCorrectParams.x_scale, 9);
  (void) poHeader->nReplaceValue(D_K_CalibYScale,
                                 m_tCorrectParams.y_scale, 9);
  (void) poHeader->nReplaceValue(D_K_CalibRatio,
                                 m_tCorrectParams.ratio, 9);
  (void) poHeader->nReplaceValue(D_K_CalibVertSlope,
                                 m_tCorrectParams.ver_slope, 9);
  (void) poHeader->nReplaceValue(D_K_CalibHorzSlope,
                                 m_tCorrectParams.horz_slope);
  (void) poHeader->nReplaceValue(D_K_CalibRadialA1,
                                 m_tDistorParams.a1, 9);
  (void) poHeader->nReplaceValue(D_K_CalibRadialA,
                                 m_tDistorParams.a);
  (void) poHeader->nReplaceValue(D_K_CalibRadialB,
                                 m_tDistorParams.b);
  (void) poHeader->nReplaceValue(D_K_CalibRadialC,
                                 m_tDistorParams.c);
  (void) poHeader->nReplaceValue(D_K_CalibSpacing,
                                 m_tDistorParams.spacing);
  (void) poHeader->nReplaceValue(D_K_CalibXMaskPoints,
                                 m_a2nhorizontal_limits[1] -
                                 m_a2nhorizontal_limits[0]);
  (void) poHeader->nReplaceValue(D_K_CalibYMaskPoints,
                                 m_a2nvertical_limits[1] -
                                 m_a2nvertical_limits[0]);

// Next should be arrays of 4 numbers
//  (void) poHeader->nReplaceValue(D_K_CalibLeftMaskPoint, 4,
//                                 m_a2nhorizontal_limits[0]);
//  (void) poHeader->nReplaceValue(D_K_CalibCenterMaskPoint, 4, 
//  (void) poHeader->nReplaceValue(D_K_CalibRightMaskPoint, 4, 
//                                 m_a2nhorizontal_limits[1]);
//  (void) poHeader->nReplaceValue(D_K_CalibBottomMaskPoint, 4, 
//                                 m_a2nvertical_limits[1]);
//  (void) poHeader->nReplaceValue(D_K_CalibTopMaskPoint, 4, 
//                                 m_a2nvertical_limits[0]);
/*
X_CENTER=577.983459;
Y_CENTER=577.978210;
X_PT_CENTER=150.000000;
Y_PT_CENTER=150.000000;
X_SCALE=14.214501;
Y_SCALE=14.144815;
RATIO=0.995098;
VER_SLOPE=0.037659;
HORZ_SLOPE=-0.027023;
RADIAL_A1=0.995059;
RADIAL_A=1;
RADIAL_B=0;
RADIAL_C=0;
SPACING=1.000000;
X_BEAM=576.000000;
Y_BEAM=576.000000;
X_SIZE=1152;
Y_SIZE=1152;
PIXEL_SIZE=0.070700;
XINT_START=-7;
XINT_STEP=4;
YINT_START=-7;
YINT_STEP=4;
XINV_START=-42;
XINV_STEP=4;
YINV_START=-60;
YINV_STEP=4;
BAD_FLAG=100000000;
PSCALE=1.95312500e-03;
COMMENT=These fields have been added to determine module orientation;
X_MASK_POINTS=81;
Y_MASK_POINTS=81;
LEFT___MASK_POINT=      2.915     597.162      -2.932     576.978;
CENTER_MASK_POINT=    561.897     575.407     562.839     576.978;
RIGHT__MASK_POINT=   1128.406     565.124    1128.610     576.978;
BOTTOM_MASK_POINT=    543.297      12.200     562.839      11.207;
TOP____MASK_POINT=    585.228    1137.982     562.839    1142.749;
*/
  poHeader->vSetState(eImage_header_available_state);
  return (0);
}

int 
Ccalibrate::nDoCalibrate(const Cstring& sCalibOptions)
{
  // Perform the refinement with the command options in sCalibOptions

  int      i, j;
  int      nStat;
  Cstring  sTemp;
  float    fTemp;
  int      nTemp;
  int      a4nTemp[4];
  float    a3fTemp[3];

  Cstring *psWords;

  static Cstring ssErrorMissing 
    = "ERROR: Ccalibrate::nDoCalibrate - missing option argument!";
  static Cstring ssErrorInvalid
    = "ERROR: Ccalibrate::nDoCalibrate - invalid argument> "; 

  char     a2cSpace[2]= {' ', '\0'};
  int      argc;

  bool     bFirstTime = TRUE;

  // Remove leading whitespace from sCalibOptions

  sTemp = sCalibOptions;
  while (   (' '  == sTemp.GetAt(0))
         || ('\n' == sTemp.GetAt(0)) 
         || ('\t' == sTemp.GetAt(0)))
    sTemp = sTemp.after(0);

  // Find number of words in sCalibOptions

  argc   = 0;
  i      = -1;
  do
    {
      i = sTemp.find(a2cSpace[0], i+1);
      argc++;
    } while (0 <= i);

  // argc will be at least 1
  // Split sTemp into words

  psWords = new Cstring [argc];
  split(sTemp, psWords, argc, a2cSpace);

  nStat = 0;

  // Now parse each word

  m_nVerbose = 6;
  for (i = 0; (i < argc) && (0 == nStat); i++) 
    {
      sTemp = psWords[i];
      if (3 < m_nVerbose)
        cout << "***Option string: >>" << sTemp << "<<" << endl;
      if ("-verbose" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%d", &nTemp);
              if (1 == nStat) 
                {
                  m_nVerbose = nTemp;
                  nStat = 0;
                }
              else
                cout << ssErrorInvalid << psWords[i] << " < !\n";
            }
          else
            {
              cout << ssErrorMissing;
            }
        }
      else if ("-dump" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sDumpName = sTransSymbol(psWords[i]);
              m_ndump     = 1;
              if (NULL != m_poOutDump)
                delete m_poOutDump;
              (void) nFileAppendVersion(m_sDumpName, TRUE);
              m_poOutDump = new ofstream(m_sDumpName.string());
              if (m_poOutDump->rdbuf()->is_open()) 
                {
                  cout << "Dump output file " << m_sDumpName << " opened!" 
                       << endl;
                  *m_poOutDump << "! dtcalibrate:  Copyright (c) 2006 Rigaku"
                               << "\n! " << D_K_DTREKVersion 
                               << "\n! Dump file: " << m_sDumpName
                               << "\n! " << sGetDate() << ' ' << sGetTime() << '\n' << endl;

                }
              else
                {
                  cout << "Dump output file not opened!\n" << endl;
                }

                
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-dscale" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%f", &fTemp);
              if (1 == nStat) 
                {
                  m_fdarksub_scale = fTemp;
                  nStat   = 0;
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-dadd" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%f", &fTemp);
              if (1 == nStat) 
                {
                  m_fdarksub_const = fTemp;
                  nStat   = 0;
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-beamx" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%f", &fTemp);
              if (1 == nStat) 
                {
                  m_a2fbeam_position[0] = fTemp;
                  nStat   = 0;
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-beamy" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%f", &fTemp);
              if (1 == nStat) 
                {
                  m_a2fbeam_position[1] = fTemp;
                  nStat   = 0;
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-seq" == sTemp)
        {
          i++;
          for (j=0; (j < 2) && (i < argc); j++, i++)
            {
              nStat = sscanf(psWords[i].string(), "%d", &nTemp);
              if (1 == nStat)
              {
                m_a2nSeqNum[j] = nTemp;
                nStat          = 0;
              }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          i--;
          if (2 != j) 
            {
              cout << ssErrorMissing;
            }
        }
      else if ("-write" == sTemp)
        {
          i++;
          if (i < argc)
            {
              m_sDistorName = psWords[i];
              nStat = nWriteDistor();
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
#ifndef CORRECT_ONLY
#ifndef DISTOR_ONLY
      else if ("-shiftmm" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%f", &a3fTemp[0]);
              if (1 == nStat)               
                {
                  i++;
                  nStat = 2;
                  if (i < argc)
                    nStat = sscanf(psWords[i].string(), "%f", &a3fTemp[1]);
                }
              if (1 == nStat) 
                {
                  i++;
                  nStat = 2;
                  if (i < argc)
                    nStat = sscanf(psWords[i].string(), "%f", &a3fTemp[2]);
                }
              if (1 == nStat) 
                {
                  nStat = nShiftMM(a3fTemp[0], a3fTemp[1], a3fTemp[2]);
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
#endif
      else if ("-peakmin" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%f", &fTemp);
              if (1 == nStat) 
                {
                  m_nmin_peak = nint(fTemp);
                  nStat   = 0;
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-peakcent" == sTemp)
        {
          i++;
          for (j=0; (j < 2) && (i < argc); j++, i++)
            {
              nStat = sscanf(psWords[i].string(), "%d", &nTemp);
              if (1 == nStat)
              {
                m_a2ncenter_peak[j] = nTemp;
                nStat          = 0;
              }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          i--;
          if (2 != j) 
            {
              cout << ssErrorMissing;
            }
        }
#ifndef DISTOR_ONLY
      else if (   ("-maxdistor" == sTemp) 
               || ("-maxdistor" == sTemp) )
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%f", &fTemp);
              if (1 == nStat) 
                {
                  m_fReferenceDistortMaxFactor = fTemp;
                  nStat   = 0;
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-post" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sPostName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
#endif
      else if ("-mask" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sMaskName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
#ifndef DISTOR_ONLY
      else if ("-mark" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sMarkName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-badpix" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sBadpixName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-tbadpix" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sTBadpixName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-refer" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sReferName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-rej" == sTemp) 
        {
        }
      else if ("-list" == sTemp) 
        {
          nStat = nList();
        }
      else if ("-gomarks" == sTemp) 
        {
          nStat = nMakeMarks();
        }
      else if ("-gopost" == sTemp) 
        {
          if ("" == m_sPostName)
            m_sPostName = "posttransform.img";
          nStat = nPostTransform();
        }
      else if ("-gotest" == sTemp) 
        {
          nStat = nMakeT();
        }
#endif
      else if (   ("-godistor" == sTemp) 
               || ("-dodistor" == sTemp) 
               || ("-calcdistor" == sTemp) )
        {
          nStat = nGoDistor();
        }
      else if ( ("-peaksfile" == sTemp) || ("-peaksfile" == sTemp) )
        {
          i++;
          if (i < argc)
            {
              m_sMaskPeaksName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
	      nStat = -1;
            }
	  i++;
	  if (i < argc)
	    {
              nStat = sscanf(psWords[i].string(), "%f", &fTemp);
              if (1 == nStat) 
                {
                  m_fMaskPeaksDivisor = fTemp;
		  if (1.0 > m_fMaskPeaksDivisor)
		    {
		      cout << ssErrorInvalid;
		      nStat = -1;
		    }
		  else
		    {
		      nStat   = 0;
		    }
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
	    }
	  else
	    {
              cout << ssErrorInvalid;
	      nStat = -1;
	    }
        }
#ifndef DISTOR_ONLY
      else if (   ("-gononunf" == sTemp) 
               || ("-dononunf" == sTemp) 
               || ("-calcnonunf" == sTemp) )
        {
          nStat = nGoNonunf();
        }
      else if ("-gointerp" == sTemp) 
        {
          nStat = nGoInterp();
        }
      else if ("-gomerge" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sMergeInfoName = psWords[i];
              nStat = nMergeTables();
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
#endif
#endif
#ifndef DISTOR_ONLY
      else if ("-nonunf" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sNonunfName = psWords[i];

              // Delete any existing nonunf object, which will
              // force a new one to be read/created, if needed

              if (NULL != m_poImgNonunf)
                {
                  delete m_poImgNonunf;
                  m_poImgNonunf = NULL;
                }
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
#endif
      else if ("-dist" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%f", &fTemp);
              if (1 == nStat) 
                {
                  m_fxtod_distance = fTemp;
                  nStat   = 0;
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-maskangle" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              nStat = sscanf(psWords[i].string(), "%f", &fTemp);
              if (1 == nStat) 
                {
                  m_fmask_angle = fTemp;
                  nStat   = 0;
                }
              else
                cout << ssErrorInvalid << psWords[i].string() << " < !\n";
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
#ifndef DISTOR_ONLY
      else if ("-flood" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sFloodName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-dark" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sDarkName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-transform" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sTransformName = psWords[i];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
#endif
      else if ("-distor" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sDistorName = psWords[i];
              if (NULL != m_poImgDistorXINT)
                {
                  cout << "WARNING! new distor name with existing distor object!"
                       << endl;
                }
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
#ifndef DISTOR_ONLY
      else if ("-tlimits" == sTemp) 
        {
          nStat = 0;
          for (j = 0; (j < 4) && (0 == nStat); j++)
            {
              i++;
              if (i < argc)
                {
                  nStat = sscanf(psWords[i].string(), "%d", &nTemp);
                  if (1 == nStat) 
                    {
                      a4nTemp[j] = nTemp;
                      nStat   = 0;
                    }
                  else
                    cout << ssErrorInvalid << psWords[i].string() << " < !\n";
                }
            }
          if (0 == nStat)
            {
              m_a2nTOutOrig[0] = a4nTemp[0];
              m_a2nTOutOrig[1] = a4nTemp[1];
              m_a2nTOutDim[0]  = a4nTemp[2];
              m_a2nTOutDim[1]  = a4nTemp[3];
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if (   ("-gocorrect" == sTemp) 
               || ("-docorrect" == sTemp) 
               || ("-gotransform" == sTemp) 
               || ("-calccorrect" == sTemp) )
        {
          i = i + 2;
          if (i < argc)
            {
#ifndef CORRECT_ONLY
              m_bCreateTransform = ("-gotransform" == sTemp);
#else
              m_bCreateTransform = FALSE;
#endif
              if (psWords[i-1].contains('?'))
                nStat = nGoCorrectScan(psWords[i-1], psWords[i]);
              else
                nStat = nGoCorrect(psWords[i-1], psWords[i]);
            }
          else
            {
              cout << ssErrorMissing;
            }
        }
      else if ("-correct" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              sTemp = psWords[i];
              sTemp.downcase();
              nTemp = 0;
              if (sTemp.contains('d'))
                nTemp = nTemp | CAL_DO_DARK;
              if (sTemp.contains('n'))
                nTemp = nTemp | CAL_DO_NONUNF;
              if (sTemp.contains('s'))
                nTemp = nTemp | CAL_DO_SPATIAL;
              if (sTemp.contains('t'))
                nTemp = nTemp | CAL_DO_TRANSFORM;
              if (0 == nTemp)
                {
                  nTemp = CAL_DO_DARK | CAL_DO_TRANSFORM;
                  cout << ssErrorInvalid << psWords[i].string() << " < !\n";
                }
              m_nCorrectionMask = nTemp;
            }
          else
            {
              cout << ssErrorInvalid;
            }
        }
      else if ("-out" == sTemp) 
        {
          i++;
          if (i < argc)
            {
              m_sOutHeader =  psWords[i];
            }
          else
            {
              cout << ssErrorMissing;
            }
        }
      else if ("-display" == sTemp)
        {
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
           if (NULL == m_poXprop)
             {
               m_poXprop = new CXprop("DTDISPLAY_WINDOW_ID");
             }
#endif
        }
      else if ("-prompt" == sTemp)
        {
          // Ignore this
        }
#endif
      else if ("" != sTemp)
        {
          cout << ssErrorInvalid << psWords[i] << " < !\n";
        }
    }
 delete [] psWords;
 return (nStat);
}

int
Ccalibrate::nNotifyDisplay(const Cstring sImageIn, 
                           const Cstring sReflnlistIn, const Cstring sHeaderIn)
{
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  if (NULL != m_poXprop)
    {
      Cstring sUpdateString1 = "";
      Cstring sUpdateString2 = "";
      Cstring sReflnlist;

      if ("" != sImageIn)
        {
          // Tell dtdisplay to display an image

          sUpdateString1 = "DTDISPLAY_IMAGE_UPDATE ";
          sUpdateString2 =   sFileGetDirectory(sImageIn)
                           + sFileGetBasename(sImageIn);
        }
      else if ("" != sReflnlistIn)
        {
          sUpdateString1 = "DTDISPLAY_REFLN_UPDATE ";
        }
      else if ("" != sHeaderIn)
        {
          sUpdateString1 = "DTDISPLAY_REFLN_UPDATE ";
        }

      if ("" != sReflnlistIn)
        {
          if (sUpdateString1.contains("IMAGE"))
            {
              sUpdateString2 +=  " Reflnlist: " 
                               + sFileGetDirectory(sReflnlistIn)
                               + sFileGetBasename(sReflnlistIn);
              
            }
          else
            {
              sUpdateString2 =  sFileGetDirectory(sReflnlistIn)
                              + sFileGetBasename(sReflnlistIn);
            }
        }

      if ("" != sHeaderIn)
        {
          sUpdateString2 +=  " Header: " + sFileGetDirectory(sHeaderIn)
                           + sFileGetBasename(sHeaderIn);
        }
      cout << "\nDisplay update: >>>" << sUpdateString1 
           << " * " << sUpdateString2 << "<<<" << endl;
      m_poXprop->hSetProperty(sUpdateString1, sUpdateString2);
    }

#endif
  return 0;
}

int
Ccalibrate::nGoDistor(void)
{
  int nStat = 0;

  cout << "***Calculate spatial distortion tables..." << endl;

  // 1. Subtract dark image from mask image (nSubtract)
  // 2. Scan the mask for peaks (nGetMask)
  // 3. Fill in expected, but missing peaks
  // 4. Calculate interpolation tables

  if (NULL != m_poImgMask)
    delete m_poImgMask;

  cout << "***   Reading mask image... " << m_sMaskName << endl;
  m_poImgMask = new Cimage(m_sMaskName);

  if (!m_poImgMask->bIsAvailable())
    {
      cout << "ERROR, mask image unavailable!" << endl;
      nStat = -1;
      return (nStat);
    }
  
  nNotifyDisplay(m_sMaskName);

  if (0.0 != m_fdarksub_scale)
    {
      // Need dark image, re-read it for now
      
      if (NULL != m_poImgDark)
        delete m_poImgDark;

      cout << "***   Reading dark image... " << m_sDarkName << endl;
      m_poImgDark = new Cimage(m_sDarkName);

      if (!m_poImgDark->bIsAvailable())
        {
          cout << "ERROR, dark image unavailable!" << endl;
          nStat = -1;
          return (nStat);
        }

      // Subtract dark image from mask image

      cout << "\n***   Correct mask image for dark image..." << endl;
      nStat = nSubtract(m_poImgMask, m_poImgDark, m_fdarksub_scale,
                        m_fdarksub_const, (float)m_a3nmin_pixel[2]);
      if (0 != nStat)
        {
          cout << "ERROR performing MASK - DARK!" << endl;
          return (nStat);
        }
    }

  m_tMaskParams.away             = m_fcent_to_maxval;
  m_tMaskParams.back             = m_nbackground_pixels;
  m_tMaskParams.beam[0]          = m_a2fbeam_position[0];
  m_tMaskParams.beam[1]          = m_a2fbeam_position[1];
  m_tMaskParams.center[0]        = m_a2ncenter_peak[0];
  m_tMaskParams.center[1]        = m_a2ncenter_peak[1];
  m_tMaskParams.dis              = m_fpeak_distance;
  m_tMaskParams.lim              = m_nmin_peak;
  m_tMaskParams.maskd            = 0;
  m_tMaskParams.mask_angle       = m_fmask_angle;
  m_tMaskParams.dump             = m_ndump;
  m_tMaskParams.misses           = m_nbad_peaks;
  m_tMaskParams.move             = m_fcent_to_pred;
  m_tMaskParams.pkrat            = 0.0;
  m_tMaskParams.pkwid            = m_npeak_size;
  m_tMaskParams.radial           = m_nradial_distortion;
  m_tMaskParams.search_center[0] = m_a2nsearch_center[0];
  m_tMaskParams.search_center[1] = m_a2nsearch_center[1];
  m_tMaskParams.resrad           = (int) m_fsearch_radius;
  m_tMaskParams.resrad2          = (int)(m_fsearch_radius*m_fsearch_radius);
  m_tMaskParams.vertical[0]      = m_a2nvertical_limits[0];
  m_tMaskParams.vertical[1]      = m_a2nvertical_limits[1];
  m_tMaskParams.source           = m_fxtod_distance;
  m_tMaskParams.horizontal[0]    = m_a2nhorizontal_limits[0];
  m_tMaskParams.horizontal[1]    = m_a2nhorizontal_limits[1];
  m_tMaskParams.work             = 1;
  m_tMaskParams.maskparams[0]    = 90.0;       // Angle between mask holes
  m_tMaskParams.maskparams[1]    = 1.0;        // Ratio of ?
  m_tMaskParams.pixel_size       = m_fpixel_size;
  m_tMaskParams.mask_spacing     = m_fmask_spacing;

  // Find peaks in the mask image

  nStat = nGetMask(m_poImgMask, m_tMaskParams);

  if (QGOOD != nStat)
    {
      cout << "nStat from nGetMask: " << nStat << endl;
    }
  else
    {
      //***************** Calculate Interpolation Tables *********************

      nStat = nCalcInterpTable();

      if (0 == nStat)
        {
          // Write out distor files
          
          nStat = nWriteDistor();
        }
    }
  
  return (nStat);
}

int 
Ccalibrate::nSubtract(Cimage *poImage1, Cimage *poImage2,
                      const float fScale, const float fOffset,
                      const float fMin, Cimage *poImageOut)
{
  // If NULL == poImageOut:
  //    poImage1 = poImage1 - (fScale * poImage2) + fOffset (min fMin)
  // else
  //  poImageOut = poImage1 - (fScale * poImage2) + fOffset (min fMin)
  
  // Subtract poImage2 from poImage1 on a pixel by pixel basis.
  // This should be in the Cimage class.

  // Do this as floating point

  float fFloodScale = 1.0;
  if ("" != sGetEnv("DTEXPOSE_SCALE"))
    {
      fFloodScale = atof(sGetEnv("DTEXPOSE_SCALE").string());
      cout << "Dark-subtracted result image scaled by a factor of " 
           << fFloodScale << endl;
    }

  if ( (0.0 == fScale) && (0.0 == fOffset) )
    {
      // Nothing to do

      return (0);
    }

  if (  (!poImage1->bIsAvailable() || !poImage2->bIsAvailable()) )
    {
      cout << "ERROR in nSubtract, Input image(s) not available()!" << endl;
      return (-1);
    }

  if ( (NULL != poImageOut) && !poImageOut->bIsAvailable())
    {
      cout << "ERROR in nSubtract, Output image(s) not available()!" << endl;
      return (-1);
    }

  int   i, j;
  float fTemp;
  float fVal1, fVal2;
  float fSatVal;

  unsigned short int uiVal1;
  int   nDim0, nDim1, nNumPix;
  int   nMinCount = 0;
  int   nMaxCount = 0;
  Cimage *poImageResult;

  if (NULL != poImageOut)
    poImageResult = poImageOut;
  else
    poImageResult = poImage1;

  poImage1->nGetDimensions(&nDim0, &nDim1);
  nNumPix = nDim0 * nDim1;
  fSatVal = poImage1->fGetSatValue();

  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
        {
          // This is slow, but it will always work, no matter
          // what the data type of the pixels

          fVal1 = (poImage1->*poImage1->prfGetPixel)(i, j);
          fVal2 = (poImage2->*poImage2->prfGetPixel)(i, j);
          fTemp = fVal1 - (fVal2 * fScale) + fOffset;
	  //+jwp 2004-11-17
	  fTemp *= fFloodScale; // Scale result by a constant
	  //-jwp 2004-11-17
//          if (  (fVal1 == fSatVal) || (fVal2 == fSatVal))
//            {
//              cout << "Sat val found at i,j: " << i << ' ' << j << ' '
//                   << fVal1 << ' ' << fVal2 << ' ' << fTemp << endl;
//            }
          if (fMin > fTemp)
            {
              fTemp = fMin;
              nMinCount++;
            }
          else if (fSatVal <= fTemp)
            {
              // This should check versus Image->SatVal
              fTemp = fSatVal;
              nMaxCount++;
            }

          // TODO mark bad pixels (too much subtraction, etc) in a BadPix image

          (poImageResult->*poImageResult->prnSetPixel)(i, j, fTemp);
        }
    }
  cout << "      Subtract results:"
       << "\n      Subtract scale:                    " << setw(9) << fScale
       << "\n      Subtract constant:                 " << setw(9) << fOffset
       << "\n      Number of pixels truncated up   to " << setw(9) << fMin 
       << " : " << setw(9) << nMinCount
       << " (" << setw(4) << (float) nMinCount / (float) nNumPix * 100.0 << "%)"
       << "\n      Number of pixels truncated down to " << setw(9) << fSatVal 
       << " : " << setw(9) << nMaxCount 
       << " (" << setw(4) << (float) nMaxCount / (float) nNumPix * 100.0 << "%)" << endl;
  return (0);
}
 
int
Ccalibrate::nGetMask(Cimage *poMask, tagMask_parameters tMaskParams)
{
  // Search the mask for the expected positions of the peaks and
  // store then in the Creflnlist m_poReflnlistMask.

  Cstring sTemp;

  cout << "\n***   Get peaks in mask image..." << endl;

  if (NULL == poMask)
    {
      cout << "ERROR in nGetMask, MASK not available()!" << endl;
      return (-1);
    }

  if (!poMask->bIsAvailable())
    {
      cout << "ERROR in nGetMask, MASK not available()!" << endl;
      return (-1);
    }

  // Initialize the list of found peaks in the MASK

  if (NULL != m_poReflnlistMask)
    {
      delete m_poReflnlistMask;
    }

  m_poReflnlistMask = new Creflnlist();

  Creflnlist *poPeaks;

  poPeaks = m_poReflnlistMask;

  (void) poPeaks->nExpandGetField(poPeaks->ms_sfObsPx0);
  (void) poPeaks->nExpandGetField(poPeaks->ms_sfObsPx1);
  (void) poPeaks->nExpandGetField(poPeaks->ms_sfCalcPx0);
  (void) poPeaks->nExpandGetField(poPeaks->ms_sfCalcPx1);
  (void) poPeaks->nExpandGetField(poPeaks->ms_snNonunfFlag);

//jwp+
  // Modified from Marty Stanton's Fortran code found in calmask.f.

  // Finds peaks in poMask and store in m_poReflnlistMask

  Cstring sDmpstr;
      
  int i, j, k;
  int nXMin, nYMin, nXMax, nYMax;
  int nStat, nMissed;

//      real    x, y, inten,
//     $     predx, predy, delx, dely, cenx, ceny
//      integer    cx, cy
  float fX, fY;
  float fInten;
  float fXdel, fYdel;
  int   nCX, nCY;

  int nLineFit;
  int imask, jmask;
//      integer line_fit
//      integer imask, jmask
//      real xnum(5), ynum(5), s2, r2
//      integer pts
//      real xfit(MASK_SIZE), yfit(MASK_SIZE),
//     $     zfit(MASK_SIZE), wfit(MASK_SIZE)

  float afXfit[MASK_SIZE];
  float afYfit[MASK_SIZE];
  float afZfit[MASK_SIZE];
  float afWfit[MASK_SIZE];
  float aXnum[5], aYnum[5], s2, r2;

  Crefln *poRefln;

  // Initialize the structure which holds info about peaks in mask image

  for  (j = 0; j < MASK_SIZE; j++)
    {
      for (i = 0; i < MASK_SIZE; i++)
        {
          m_tMaskPoints.status[j][i] = QBAD;
          m_tMaskPoints.xo[j][i]     = 0.0;
          m_tMaskPoints.yo[j][i]     = 0.0;
          m_tMaskPoints.xc[j][i]     = 0.0;
          m_tMaskPoints.yc[j][i]     = 0.0;
          m_tMaskPoints.xu[j][i]     = 0.0;
          m_tMaskPoints.yu[j][i]     = 0.0;
          m_tMaskPoints.x_diff[j][i] = 0.0;
          m_tMaskPoints.y_diff[j][i] = 0.0;
        }
    }
  m_tMaskPoints.xb    = MASK_SIZE - 1;
  m_tMaskPoints.xe    = 0;
  m_tMaskPoints.yb    = MASK_SIZE - 1;
  m_tMaskPoints.ye    = 0;
  m_tMaskPoints.pts_x = 0;
  m_tMaskPoints.pts_y = 0;
      
  // Initial searching limits

  nCX   = tMaskParams.search_center[0];
  nCY   = tMaskParams.search_center[1];
  nXMin = max ((nCX-tMaskParams.resrad), tMaskParams.horizontal[0] );
  nXMax = min ((nCY+tMaskParams.resrad), tMaskParams.horizontal[1]-1);
  nYMin = max ((nCX-tMaskParams.resrad), tMaskParams.vertical[0] );
  nYMax = min ((nCY+tMaskParams.resrad), tMaskParams.vertical[1]-1);
      
  if ("" != m_sMaskPeaksName)
    {
      // SPECIAL KLUDGE: read the mask points from external reflnlist file

      if (NULL != m_poReflnlistMask)
	delete m_poReflnlistMask;

      cout << "INFO: reading mask peaks from file: "
	   << m_sMaskPeaksName << endl;
      cout << "INFO: peak positional divisor is: " << m_fMaskPeaksDivisor << endl << flush;
      m_poReflnlistMask = new Creflnlist(m_sMaskPeaksName);
      if (!m_poReflnlistMask->bIsAvailable())
	{
	  cout << "ERROR reading mask reflnlist file: " << m_sMaskPeaksName
	       << "\n\n" << flush;
	  exit (1);
	}
      poPeaks = m_poReflnlistMask;
      for (k = 0; k < poPeaks->nGetNumReflns(); k++)
	{
	  poRefln = poPeaks->poGetRefln(k);
	  i       = poRefln->nGetH();
	  j       = poRefln->nGetK();
	  m_tMaskPoints.inten[j][i]   = poRefln->fGetIntensity();
	  m_tMaskPoints.xo[j][i]      = poRefln->fGetField(poPeaks->m_nFI_fCalcPx0) / m_fMaskPeaksDivisor; 
	  m_tMaskPoints.yo[j][i]      = poRefln->fGetField(poPeaks->m_nFI_fCalcPx1) / m_fMaskPeaksDivisor;
	  poRefln->vSetField(poPeaks->m_nFI_fCalcPx0, m_tMaskPoints.xo[j][i]);
	  poRefln->vSetField(poPeaks->m_nFI_fCalcPx1, m_tMaskPoints.yo[j][i]);
	  m_tMaskPoints.status[j][i]  = poRefln->nGetL();
	  if (QGOOD == m_tMaskPoints.status [j][i])
	    {
	      if ( i < m_tMaskPoints.xb ) m_tMaskPoints.xb = i;
	      if ( i > m_tMaskPoints.xe ) m_tMaskPoints.xe = i;
	      if ( j < m_tMaskPoints.yb ) m_tMaskPoints.yb = j;
	      if ( j > m_tMaskPoints.ye ) m_tMaskPoints.ye = j;
	    }
	}
    }
  else 
    {
      // Find a good peak somewhere around center

      fX = tMaskParams.center[0];
      fY = tMaskParams.center[1];
      if (NULL != m_poOutDump)
	{
	  *m_poOutDump << "! Searching for center point at : " 
                       <<  fX << ' ' << fY << endl;
	}
  
      nStat = nSearchCenter(poMask, tMaskParams, &fX, &fY, &fInten);

      if (NULL != m_poOutDump)
	{
	  vDumpPoint (CTR, CTR, fX, fY, fInten, nStat);
	}

      if (QGOOD != nStat)
	{
	  if (NULL != m_poOutDump)
	    {
	      *m_poOutDump << "! Center point no good!" << endl;
	    }

	  cout << "ERROR in Ccalibrate::nGetMask - Cannot find center point at "
	       << tMaskParams.center[0] << ' '
	       << tMaskParams.center[1] << '\n'
	       << "        Try choosing another center" << endl;
	  return (-1);
	}

      m_tMaskPoints.xo[CTR][CTR]     = fX;
      m_tMaskPoints.yo[CTR][CTR]     = fY;
      m_tMaskPoints.inten[CTR][CTR]  = fInten;
      m_tMaskPoints.status[CTR][CTR] = QGOOD;

      // Now that a center peak has been found ...
      // find peaks immediately to the right
      // to provide a start to find the rest of the points
      // TODO: next peak should use a vector to find it, so that it does not
      //       have to be 90 degrees away

      fX = fX + tMaskParams.dis * cosf(tMaskParams.mask_angle 
				       * Gs_fRADIANS_PER_DEGREE);
      fY = fY + tMaskParams.dis * sinf(tMaskParams.mask_angle
				       * Gs_fRADIANS_PER_DEGREE);
      if (2 < m_nVerbose)
	cout << "Search for peak(s) 'right' of center..." << endl;

      i = CTR + 1;
      j = CTR;

      if (NULL != m_poOutDump)
	{
	  *m_poOutDump << "! Searching for right point at : " << endl;
	}

      nStat = nSearchCenter(poMask, tMaskParams, &fX, &fY, &fInten);
      if (NULL != m_poOutDump)
	{
	  vDumpPoint (i, j, fX, fY, fInten, nStat);
	}

      if (QGOOD != nStat)
	{
	  if (NULL != m_poOutDump)
	    {
	      *m_poOutDump << "! Right point no good !" << endl;
	    }
	  cout << "ERROR in nGetMask - Can not find a good peak to right of center"
	       << "\n   Try choosing another center or distance" << endl;
	  return (-1);
	}

      m_tMaskPoints.xo[j][i]     = fX;
      m_tMaskPoints.yo[j][i]     = fY;
      m_tMaskPoints.inten[j][i]  = fInten;
      m_tMaskPoints.status[j][i] = QGOOD;
      
      // Move to the right from center, filling in good peaks, 
      // until a limit is reached

      int   nNumPoints;
      float fXpred;
      float fYpred;
      float fI, fJ;

      i = CTR + 2;
      j = CTR;
      
      fXpred = 2.0 * m_tMaskPoints.xo[j][i-1] - m_tMaskPoints.xo[j][i-2];
      fYpred = 2.0 * m_tMaskPoints.yo[j][i-1] - m_tMaskPoints.yo[j][i-2];

      for (nNumPoints = 0; nNumPoints < MASK_SIZE; nNumPoints++)
	afWfit[nNumPoints] = 1.0;

      nMissed  = 0;

      while (   (fXpred < nXMax)
             && (  (  (fXpred-(float)nCX) * (fXpred-(float)nCX)
		    + (fYpred-(float)nCY) * (fYpred-(float)nCY)) 
		  <  tMaskParams.resrad2)
             && (nMissed < tMaskParams.misses) )
	{
	  fX = fXpred;
	  fY = fYpred;

	  // Find the true centroid of the spot predicted to be at fXpred, fYpred

	  nStat = nFindCenter(tMaskParams.pkwid,
			      tMaskParams.move,
			      tMaskParams.lim,
			      tMaskParams.away,
			      poMask, &fX, &fY, &fInten);
	  if (QGOOD == nStat)
	    {
	      m_tMaskPoints.xo[j][i]     = fX;
	      m_tMaskPoints.yo[j][i]     = fY;
	      m_tMaskPoints.inten[j][i]  = fInten;
	      m_tMaskPoints.status[j][i] = QGOOD;
	      nMissed            = 0;
	    }
	  else
	    {
	      m_tMaskPoints.xo[j][i]    = fXpred;
	      m_tMaskPoints.yo[j][i]    = fYpred;
	      m_tMaskPoints.inten[j][i] = 0.0;
	      nMissed           = nMissed + 1;
	    }

	  if (NULL != m_poOutDump)
	    {
	      vDumpPoint (i, j, fX, fY, fInten, nStat);
	      if (QGOOD != nStat)
		*m_poOutDump << "! Predicted position at: " 
                             << fXpred <<  ' ' << fYpred << endl;
	    }
      
	  nLineFit = 0;
	  if (QGOOD == m_tMaskPoints.status[j][i])
	    {
	      nNumPoints = 0;
	      for (imask = max (CTR, i-5); imask <= i; imask++)
		{
		  if (QGOOD == m_tMaskPoints.status[j][imask])
		    {
		      afXfit[nNumPoints] = m_tMaskPoints.xo[j][imask];
		      afYfit[nNumPoints] = m_tMaskPoints.yo[j][imask];
		      afZfit[nNumPoints] = imask;
		      nNumPoints++;
		    }
		}
	      // nLineFit = max(3, nNumPoints / 2);
	      if (6 <= nNumPoints)
		nLineFit = 3;
	      else if (4 <= nNumPoints)
		nLineFit = 2;
	      else if (2 <=  nNumPoints)
		nLineFit = 1;

	      if (2 <= nNumPoints)
		{
		  (void) nPolyFit ( afZfit-1, afXfit-1, afWfit-1, nNumPoints, nLineFit,
				    aXnum-1, &s2, &r2 );
		  (void) nPolyFit ( afZfit-1, afYfit-1, afWfit-1, nNumPoints, nLineFit,
				    aYnum-1, &s2, &r2 );
		}
	    }

	  // Try to predict where the next spot along the i'th direction
	  // should be (set fXpred, fYpred)

	  i  = i + 1;
	  fI = (float) i;
	  if (1 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fI;
	      fYpred = aYnum[0] + aYnum[1] * fI;
	    }
	  else if (2 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fI + aXnum[2] * fI * fI;
	      fYpred = aYnum[0] + aYnum[1] * fI + aYnum[2] * fI * fI;
	    }
	  else if (3 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fI + aXnum[2] * fI * fI 
		                                + aXnum[3] * fI * fI * fI;
	      fYpred = aYnum[0] + aYnum[1] * fI + aYnum[2] * fI * fI 
                                                + aYnum[3] * fI * fI * fI;
	    }
	  else
	    {
	      fXpred = 2.0 * m_tMaskPoints.xo[j][i-1]  - m_tMaskPoints.xo[j][i-2];
	      fYpred = 2.0 * m_tMaskPoints.yo[j][i-1]  - m_tMaskPoints.yo[j][i-2];
	    }
	  if (5 < m_nVerbose)
	    cout << "Next predicted X,Y " << fXpred << ' ' << fYpred
		 << " i, j: " << i << ' ' << j << endl;

	  // OK, have next predicted position, loop back to actually
	  // find the true centroid near that position

	} // end of while?
      //      end do
      m_tMaskPoints.xe = i - 1;

      // Move to the left from center, filling in good peaks, 
      // until a limit is lower limit is reached

      if (2 < m_nVerbose)
	cout << "Search for peak(s) 'left' of center..." << endl;

      i = CTR - 1;
      j = CTR;
      
      nNumPoints = 0;
      for (imask = i; imask <= i+5; imask++)
	{
	  if (QGOOD == m_tMaskPoints.status[j][imask])
	    {
	      afXfit[nNumPoints] = m_tMaskPoints.xo[j][imask];
	      afYfit[nNumPoints] = m_tMaskPoints.yo[j][imask];
	      afZfit[nNumPoints] = imask;
	      nNumPoints++;
	    }
	}
      nLineFit   = 0;
      if (6 <= nNumPoints)
	nLineFit = 3;
      else if (4 <= nNumPoints)
	nLineFit = 2;
      else if (2 <= nNumPoints)
	nLineFit = 1;

      if (2 <= nNumPoints)
	{
	  (void) nPolyFit ( afZfit-1, afXfit-1, afWfit-1, nNumPoints, nLineFit,
			    aXnum-1, &s2, &r2 );
	  (void) nPolyFit ( afZfit-1, afYfit-1, afWfit-1, nNumPoints, nLineFit,
			    aYnum-1, &s2, &r2 );
	}
      
      fI = (float) i;
      if (1 == nLineFit)
	{
	  fXpred = aXnum[0] + aXnum[1] * fI;
	  fYpred = aYnum[0] + aYnum[1] * fI;
	}
      else if (2 == nLineFit)
	{
	  fXpred = aXnum[0] + aXnum[1] * fI + aXnum[2] * fI * fI;
	  fYpred = aYnum[0] + aYnum[1] * fI + aYnum[2] * fI * fI;
	}
      else if (3 == nLineFit)
	{
	  fXpred = aXnum[0] + aXnum[1] * fI + aXnum[2] * fI * fI 
                                            + aXnum[3] * fI * fI * fI;
	  fYpred = aYnum[0] + aYnum[1] * fI + aYnum[2] * fI * fI 
	                                    + aYnum[3] * fI * fI * fI;
	}
      else
	{
	  fXpred = 2.0 * m_tMaskPoints.xo[j][i+1]  - m_tMaskPoints.xo[j][i+2];
	  fYpred = 2.0 * m_tMaskPoints.yo[j][i+1]  - m_tMaskPoints.yo[j][i+2];
	}
      if (5 < m_nVerbose)
	cout << "Next predicted X,Y " << fXpred << ' ' << fYpred
             << " i, j: " << i << ' ' << j << endl;
      nMissed = 0;

      while (   (fXpred > nXMin)
             && (  (  (fXpred-(float)nCX) * (fXpred-(float)nCX)
                    + (fYpred-(float)nCY) * (fYpred-(float)nCY)) 
                 <  tMaskParams.resrad2)
             && (nMissed < tMaskParams.misses) )
	{
	  fX = fXpred;
	  fY = fYpred;

	  nStat = nFindCenter(tMaskParams.pkwid,
			      tMaskParams.move,
			      tMaskParams.lim,
			      tMaskParams.away,
			      poMask, &fX, &fY, &fInten);
	  if (QGOOD == nStat)
	    {
	      m_tMaskPoints.xo[j][i]     = fX;
	      m_tMaskPoints.yo[j][i]     = fY;
	      m_tMaskPoints.inten[j][i]  = fInten;
	      m_tMaskPoints.status[j][i] = QGOOD;
	      nMissed            = 0;
	    }
	  else
	    {
	      m_tMaskPoints.xo[j][i]    = fXpred;
	      m_tMaskPoints.yo[j][i]    = fYpred;
	      m_tMaskPoints.inten[j][i] = 0.0;
	      nMissed           = nMissed + 1;
	    }

	  if (NULL != m_poOutDump)
	    {
	      vDumpPoint (i, j, fX, fY, fInten, nStat);
	      if (QGOOD != nStat)
		*m_poOutDump << "! Predicted position at: " 
			     << fXpred <<  ' ' << fYpred << endl;
	    }

	  nLineFit = 0;
	  if (QGOOD == m_tMaskPoints.status[j][i])
	    {
	      nNumPoints = 0;
	      for (imask = i; imask <= min(i+5, CTR); imask++)
		{
		  if (QGOOD == m_tMaskPoints.status[j][imask])
		    {
		      afXfit[nNumPoints] = m_tMaskPoints.xo[j][imask];
		      afYfit[nNumPoints] = m_tMaskPoints.yo[j][imask];
		      afZfit[nNumPoints] = imask;
		      nNumPoints++;
		    }
		}
	      if (6 <= nNumPoints)
		nLineFit = 3;
	      else if (4 <= nNumPoints)
		nLineFit = 2;
	      else if (2 <= nNumPoints)
		nLineFit = 1;

	      if (2 <= nNumPoints)
		{
		  (void) nPolyFit ( afZfit-1, afXfit-1, afWfit-1, nNumPoints, nLineFit,
                                   aXnum-1, &s2, &r2 );
		  (void) nPolyFit ( afZfit-1, afYfit-1, afWfit-1, nNumPoints, nLineFit,
                                   aYnum-1, &s2, &r2 );
		}
	    }

	  i = i - 1;

	  fI = (float) i;
	  if (1 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fI;
	      fYpred = aYnum[0] + aYnum[1] * fI;
	    }
	  else if (2 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fI + aXnum[2] * fI * fI;
	      fYpred = aYnum[0] + aYnum[1] * fI + aYnum[2] * fI * fI;
	    }
	  else if (3 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fI + aXnum[2] * fI * fI 
                                                + aXnum[3] * fI * fI * fI;
	      fYpred = aYnum[0] + aYnum[1] * fI + aYnum[2] * fI * fI 
                                                + aYnum[3] * fI * fI * fI;
	    }
	  else
	    {
	      fXpred = 2.0 * m_tMaskPoints.xo[j][i+1]  - m_tMaskPoints.xo[j][i+2];
	      fYpred = 2.0 * m_tMaskPoints.yo[j][i+1]  - m_tMaskPoints.yo[j][i+2];
	    }
	  if (5 < m_nVerbose)
	    cout << "Next predicted X,Y " << fXpred << ' ' << fYpred
                 << " i, j: " << i << ' ' << j << endl;
	} // end of while

      m_tMaskPoints.xb = i + 1;
      
      // Find a peak immediately above the peak at the center 
      // to provide a start to find the rest of the points

      if (2 < m_nVerbose)
	cout << "Search for peak(s) 'above' center..." << endl;

      i = CTR;
      j = CTR + 1;

      fX = m_tMaskPoints.xo[CTR][CTR] - tMaskParams.dis * sinf(tMaskParams.mask_angle 
                                                       * Gs_fRADIANS_PER_DEGREE);
      fY = m_tMaskPoints.yo[CTR][CTR] + tMaskParams.dis * cosf(tMaskParams.mask_angle
                                                       * Gs_fRADIANS_PER_DEGREE);
      if (NULL != m_poOutDump)
	{
	  *m_poOutDump << "! Searching for upper point at : " << endl;
	}

      nStat = nSearchCenter(poMask, tMaskParams, &fX, &fY, &fInten);
      if (NULL != m_poOutDump)
	{
	  vDumpPoint (i, j, fX, fY, fInten, nStat);
	}

      if (QGOOD != nStat)
	{
	  if (NULL != m_poOutDump)
	    {
	      *m_poOutDump << "! Upper point no good !" << endl;
	    }
	  cout << "ERROR in nGetMask - Can not find a good peak above center"
               << "\n   Try choosing another center or distance" << endl;
	  return (-1);
	}

      fXdel = fX - m_tMaskPoints.xo[CTR][CTR];
      fYdel = fY - m_tMaskPoints.yo[CTR][CTR];
      
      // Starting from the center and moving to the left,
      // fill in vertical lines up and down from the center
      // horizontal line

      if (2 < m_nVerbose)
	cout << "Search for peak(s) 'above left' of center..." << endl;

      for (i = CTR; i >= m_tMaskPoints.xb; i--)  // do i = CTR, m_tMaskPoints.xb, -1
	{
	  // Move up from center, filling in good peaks, until a
	  // upper limit is reached

	  j        = CTR + 1;
	  fXpred   = m_tMaskPoints.xo[j-1][i] + fXdel;
	  fYpred   = m_tMaskPoints.yo[j-1][i] + fYdel;
	  nMissed  = 0;
	  while (   (fYpred < nYMax)
                 && (nMissed < tMaskParams.misses)
                 && (fXpred < nXMax)
                 && (fXpred > nXMin)
                 && (  (  (fXpred-(float)nCX) * (fXpred-(float)nCX)
                        + (fYpred-(float)nCY) * (fYpred-(float)nCY)) 
                     <  tMaskParams.resrad2)
		    )
	    {
	      fX = fXpred;
	      fY = fYpred;
	      nStat = nFindCenter(tMaskParams.pkwid,
				  tMaskParams.move,
				  tMaskParams.lim,
				  tMaskParams.away,
				  poMask, &fX, &fY, &fInten);

	      if (QGOOD == nStat)
		{
		  m_tMaskPoints.xo[j][i]     = fX;
		  m_tMaskPoints.yo[j][i]     = fY;
		  m_tMaskPoints.inten[j][i]  = fInten;
		  m_tMaskPoints.status[j][i] = QGOOD;
		  nMissed            = 0;
		}
	      else
		{
		  m_tMaskPoints.xo[j][i]    = fXpred;
		  m_tMaskPoints.yo[j][i]    = fYpred;
		  m_tMaskPoints.inten[j][i] = 0.0;
		  nMissed           = nMissed + 1;
		}
	      if (NULL != m_poOutDump)
		{
		  vDumpPoint (i, j, fX, fY, fInten, nStat);
		  if (QGOOD != nStat)
		    *m_poOutDump << "! Predicted position at: " 
                                 << fXpred <<  ' ' << fYpred << endl;
		}

	      nLineFit = 0;
	      if (QGOOD == m_tMaskPoints.status[j][i])
		{
		  nNumPoints = 0;
		  for (jmask = max(CTR, j-5); jmask <= j; jmask++)
		    {
		      if (QGOOD == m_tMaskPoints.status[jmask][i])
			{
			  afXfit[nNumPoints] = m_tMaskPoints.xo[jmask][i];
			  afYfit[nNumPoints] = m_tMaskPoints.yo[jmask][i];
			  afZfit[nNumPoints] = jmask;
			  nNumPoints++;
			}
		    }

		  if (6 <= nNumPoints)
		    nLineFit = 3;
		  else if (4 <= nNumPoints)
		    nLineFit = 2;
		  else if (2 <= nNumPoints)
		    nLineFit = 1;

		  if (2 <= nNumPoints)
		    {
		      (void) nPolyFit ( afZfit-1, afXfit-1, afWfit-1, nNumPoints, nLineFit,
                                       aXnum-1, &s2, &r2 );
		      (void) nPolyFit ( afZfit-1, afYfit-1, afWfit-1, nNumPoints, nLineFit,
                                      aYnum-1, &s2, &r2 );
		    }
		}
            
	      j = j + 1;

	      fJ = (float) j;
	      if (1 == nLineFit)
		{
		  fXpred = aXnum[0] + aXnum[1] * fJ;
		  fYpred = aYnum[0] + aYnum[1] * fJ;
		}
	      else if (2 == nLineFit)
		{
		  fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ;
		  fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ;
		}
	      else if (3 == nLineFit)
		{
		  fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ 
                                                    + aXnum[3] * fJ * fJ * fJ;
		  fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ 
                                                    + aYnum[3] * fJ * fJ * fJ;
		}
	      else
		{
		  fXpred = 2.0 * m_tMaskPoints.xo[j-1][i]  - m_tMaskPoints.xo[j-2][i];
		  fYpred = 2.0 * m_tMaskPoints.yo[j-1][i]  - m_tMaskPoints.yo[j-2][i];
		}
	      if (5 < m_nVerbose)
		cout << "Next predicted X,Y " << fXpred << ' ' << fYpred
		     << " i, j: " << i << ' ' << j << endl;
	    } // end while

	  j = j - 1;
	  if (m_tMaskPoints.ye < j) m_tMaskPoints.ye = j;
         
	  // Move down from center, filling in good peaks, until a
	  // lower limit is reached

	  j = CTR - 1;
	  nNumPoints = 0;
	  for (jmask = j; jmask <= j+5; jmask++)
	    {
	      if (QGOOD == m_tMaskPoints.status[jmask][i])
		{
		  afXfit[nNumPoints] = m_tMaskPoints.xo[jmask][i];
		  afYfit[nNumPoints] = m_tMaskPoints.yo[jmask][i];
		  afZfit[nNumPoints] = jmask;
		  nNumPoints++;
		}
	    }

	  if (6 <= nNumPoints)
	    nLineFit = 3;
	  else if (4 <= nNumPoints)
	    nLineFit = 2;
	  else if (2 <= nNumPoints)
	    nLineFit = 1;
	  else
	    nLineFit = 0;

	  if (2 <= nNumPoints)
	    {
	      (void) nPolyFit ( afZfit-1, afXfit-1, afWfit-1, nNumPoints, nLineFit,
				aXnum-1, &s2, &r2 );
	      (void) nPolyFit ( afZfit-1, afYfit-1, afWfit-1, nNumPoints, nLineFit,
				aYnum-1, &s2, &r2 );
	    }
	  
	  fJ = (float) j;
	  if (1 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fJ;
	      fYpred = aYnum[0] + aYnum[1] * fJ;
	    }
	  else if (2 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ;
	      fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ;
	    }
	  else if (3 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ 
                                                + aXnum[3] * fJ * fJ * fJ;
	      fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ 
                                                + aYnum[3] * fJ * fJ * fJ;
	    }
	  else
	    {
	      fXpred = 2.0 * m_tMaskPoints.xo[CTR][i]  - m_tMaskPoints.xo[CTR+1][i];
	      fYpred = 2.0 * m_tMaskPoints.yo[CTR][i]  - m_tMaskPoints.yo[CTR+1][i];
	    }
	  if (5 < m_nVerbose)
	    cout << "Next predicted X,Y " << fXpred << ' ' << fYpred
		 << " i, j: " << i << ' ' << j << endl;

	  nMissed = 0;

	  while (   (fYpred > nYMin)
                 && (fYpred < nYMax)    // JWP Mar 11, 2009
		    && (nMissed < tMaskParams.misses)
		    && (fXpred < nXMax)
		    && (fXpred > nXMin)
		    && (  (  (fXpred-(float)nCX) * (fXpred-(float)nCX)
			     + (fYpred-(float)nCY) * (fYpred-(float)nCY)) 
			  <  tMaskParams.resrad2)
		    )
	    {
	      fX = fXpred;
	      fY = fYpred;
	      nStat = nFindCenter(tMaskParams.pkwid,
				  tMaskParams.move,
				  tMaskParams.lim,
				  tMaskParams.away,
				  poMask, &fX, &fY, &fInten);

	      if (QGOOD == nStat)
		{
		  m_tMaskPoints.xo[j][i]     = fX;
		  m_tMaskPoints.yo[j][i]     = fY;
		  m_tMaskPoints.inten[j][i]  = fInten;
		  m_tMaskPoints.status[j][i] = QGOOD;
		  nMissed            = 0;
		}
          else
            {
              m_tMaskPoints.xo[j][i]    = fXpred;
              m_tMaskPoints.yo[j][i]    = fYpred;
              m_tMaskPoints.inten[j][i] = 0.0;
              nMissed           = nMissed + 1;
            }
          if (NULL != m_poOutDump)
            {
              vDumpPoint (i, j, fX, fY, fInten, nStat);
              if (QGOOD != nStat)
                *m_poOutDump << "! Predicted position at: " 
                             << fXpred <<  ' ' << fYpred << endl;
            }

          if (QGOOD == m_tMaskPoints.status[j][i])
            {
              nNumPoints = 0;
              for (jmask = j; jmask <= min(j+5, CTR); jmask++)
                {
                  if (QGOOD == m_tMaskPoints.status[jmask][i])
                    {
                      afXfit[nNumPoints] = m_tMaskPoints.xo[jmask][i];
                      afYfit[nNumPoints] = m_tMaskPoints.yo[jmask][i];
                      afZfit[nNumPoints] = jmask;
                      nNumPoints++;
                    }
                }

              if (6 <= nNumPoints)
                nLineFit = 3;
              else if (4 <= nNumPoints)
                nLineFit = 2;
              else if (2 <= nNumPoints)
                nLineFit = 1;
              else
                nLineFit = 0;

              if (2 <= nNumPoints)
                {
                  (void) nPolyFit ( afZfit-1, afXfit-1, afWfit-1, nNumPoints, nLineFit,
                                   aXnum-1, &s2, &r2 );
                  (void) nPolyFit ( afZfit-1, afYfit-1, afWfit-1, nNumPoints, nLineFit,
                                   aYnum-1, &s2, &r2 );
                }
            }

          j = j - 1;
          fJ = (float) j;
          if (1 == nLineFit)
            {
              fXpred = aXnum[0] + aXnum[1] * fJ;
              fYpred = aYnum[0] + aYnum[1] * fJ;
            }
          else if (2 == nLineFit)
            {
              fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ;
              fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ;
            }
          else if (3 == nLineFit)
            {
              fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ 
                                                + aXnum[3] * fJ * fJ * fJ;
              fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ 
                                                + aYnum[3] * fJ * fJ * fJ;
            }
          else
            {
              fXpred = 2.0 * m_tMaskPoints.xo[j+1][i]  - m_tMaskPoints.xo[j+2][i];
              fYpred = 2.0 * m_tMaskPoints.yo[j+1][i]  - m_tMaskPoints.yo[j+2][i];
            }
          if (5 < m_nVerbose)
            cout << "Next predicted X,Y " << fXpred << ' ' << fYpred
                 << " i, j: " << i << ' ' << j << endl;
        } // end while{

      j = j + 1;
      if ( m_tMaskPoints.yb > j ) m_tMaskPoints.yb = j;

      fXdel = m_tMaskPoints.xo[CTR+1][i] - m_tMaskPoints.xo[CTR][i];
      fYdel = m_tMaskPoints.yo[CTR+1][i] - m_tMaskPoints.yo[CTR][i];

	} // end i loop
  //jwpX
      // Finally move to the right from the center
      // fill in vertical lines up and down from the center
      // horizontal line

      if (2 < m_nVerbose)
	cout << "Search for peak(s) 'above right' of center..." << endl;

      fXdel = m_tMaskPoints.xo[CTR+1][CTR] - m_tMaskPoints.xo[CTR][CTR];
      fYdel = m_tMaskPoints.yo[CTR+1][CTR] - m_tMaskPoints.yo[CTR][CTR];

      for (i = CTR+1; i <= m_tMaskPoints.xe; i++)
	{
	  // Move up from center, filling in good peaks, until a
	  // upper limit is reached

	  j        = CTR + 1;
	  fXpred   = m_tMaskPoints.xo[j-1][i] + fXdel;
	  fYpred   = m_tMaskPoints.yo[j-1][i] + fYdel;
	  nMissed  = 0;
	  nLineFit = 0;

	  while (   (fYpred < nYMax)
                 && (nMissed < tMaskParams.misses)
                 && (fXpred < nXMax)
                 && (fXpred > nXMin)
                 && (  (  (fXpred-(float)nCX) * (fXpred-(float)nCX)
                        + (fYpred-(float)nCY) * (fYpred-(float)nCY)) 
                     <  tMaskParams.resrad2)
		    )
	    {
	      fX = fXpred;
	      fY = fYpred;
	      nStat = nFindCenter(tMaskParams.pkwid,
				  tMaskParams.move,
				  tMaskParams.lim,
				  tMaskParams.away,
				  poMask, &fX, &fY, &fInten);
	      if (QGOOD == nStat)
		{
		  m_tMaskPoints.xo[j][i]     = fX;
		  m_tMaskPoints.yo[j][i]     = fY;
		  m_tMaskPoints.inten[j][i]  = fInten;
		  m_tMaskPoints.status[j][i] = QGOOD;
		  nMissed            = 0;
		}
	      else
		{
		  m_tMaskPoints.xo[j][i]    = fXpred;
		  m_tMaskPoints.yo[j][i]    = fYpred;
		  m_tMaskPoints.inten[j][i] = 0.0;
		  nMissed           = nMissed + 1;
		}
	      if (NULL != m_poOutDump)
		{
		  vDumpPoint (i, j, fX, fY, fInten, nStat);
		  if (QGOOD != nStat)
		    *m_poOutDump << "! Predicted position at: " 
                                 << fXpred <<  ' ' << fYpred << endl;
		}

	      if (QGOOD == m_tMaskPoints.status[j][i])
		{
		  nNumPoints = 0;
		  for (jmask = max(CTR, j-5); jmask <= j; jmask++)
		    {
		      if (QGOOD == m_tMaskPoints.status[jmask][i])
			{
			  afXfit[nNumPoints] = m_tMaskPoints.xo[jmask][i];
			  afYfit[nNumPoints] = m_tMaskPoints.yo[jmask][i];
			  afZfit[nNumPoints] = jmask;
			  nNumPoints++;
			}
		    }

		  if (6 <= nNumPoints)
		    nLineFit = 3;
		  else if (4 <= nNumPoints)
		    nLineFit = 2;
		  else if (2 <= nNumPoints)
		    nLineFit = 1;
		  else
		    nLineFit = 0;

		  if (2 <= nNumPoints)
		    {
		      (void) nPolyFit ( afZfit-1, afXfit-1, afWfit-1, nNumPoints, nLineFit,
					aXnum-1, &s2, &r2 );
		      (void) nPolyFit ( afZfit-1, afYfit-1, afWfit-1, nNumPoints, nLineFit,
					aYnum-1, &s2, &r2 );
		    }
		}
	      
	      j = j + 1;

	      fJ = (float) j;
	      if (1 == nLineFit)
		{
		  fXpred = aXnum[0] + aXnum[1] * fJ;
		  fYpred = aYnum[0] + aYnum[1] * fJ;
		}
	      else if (2 == nLineFit)
		{
		  fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ;
		  fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ;
		}
	      else if (3 == nLineFit)
		{
		  fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ 
                                                    + aXnum[3] * fJ * fJ * fJ;
		  fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ 
                                                    + aYnum[3] * fJ * fJ * fJ;
		}
	      else
		{
		  fXpred = 2.0 * m_tMaskPoints.xo[j-1][i]  - m_tMaskPoints.xo[j-2][i];
		  fYpred = 2.0 * m_tMaskPoints.yo[j-1][i]  - m_tMaskPoints.yo[j-2][i];
		}
	      if (5 < m_nVerbose)
		cout << "Next predicted X,Y " << fXpred << ' ' << fYpred
		     << " i, j: " << i << ' ' << j << endl;

	    } // end while

	  j = j - 1;
	  if (m_tMaskPoints.ye < j) m_tMaskPoints.ye = j;
	  
	  // Move down from center, filling in good peaks, until a
	  // lower limit is reached.

	  j          = CTR - 1;
	  nNumPoints = 0;
	  for (jmask = j; jmask <= j+5; jmask++)
	    {
	      if (QGOOD == m_tMaskPoints.status[jmask][i])
		{
		  afXfit[nNumPoints] = m_tMaskPoints.xo[jmask][i];
		  afYfit[nNumPoints] = m_tMaskPoints.yo[jmask][i];
		  afZfit[nNumPoints] = jmask;
		  nNumPoints++;
		}
	    }
	  if (6 <= nNumPoints)
	    nLineFit = 3;
	  else if (4 <= nNumPoints)
	    nLineFit = 2;
	  else if (2 <= nNumPoints)
	    nLineFit = 1;
	  else
	    nLineFit = 0;

	  if (2 <= nNumPoints)
	    {
	      (void) nPolyFit ( afZfit-1, afXfit-1, afWfit-1, nNumPoints, nLineFit,
				aXnum-1, &s2, &r2 );
	      (void) nPolyFit ( afZfit-1, afYfit-1, afWfit-1, nNumPoints, nLineFit,
				aYnum-1, &s2, &r2 );
	    }
            
	  fJ = (float) j;
	  if (1 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fJ;
	      fYpred = aYnum[0] + aYnum[1] * fJ;
	    }
	  else if (2 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ;
	      fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ;
	    }
	  else if (3 == nLineFit)
	    {
	      fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ 
                                                + aXnum[3] * fJ * fJ * fJ;
	      fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ 
                                   		+ aYnum[3] * fJ * fJ * fJ;
	    }
	  else
	    {
	      fXpred = 2.0 * m_tMaskPoints.xo[CTR][i]  - m_tMaskPoints.xo[CTR+1][i];
	      fYpred = 2.0 * m_tMaskPoints.yo[CTR][i]  - m_tMaskPoints.yo[CTR+1][i];
	    }

	  if (5 < m_nVerbose)
	    cout << "Next predicted X,Y " << fXpred << ' ' << fYpred
		 << " i, j: " << i << ' ' << j << endl;

	  nMissed = 0;

	  while (   (fYpred > nYMin)
                 && (fYpred < nYMax)    // JWP Mar 11, 2009
                 && (nMissed < tMaskParams.misses)
                 && (fXpred < nXMax)
                 && (fXpred > nXMin)
                 && (  (  (fXpred-(float)nCX) * (fXpred-(float)nCX)
                        + (fYpred-(float)nCY) * (fYpred-(float)nCY)) 
                     <  tMaskParams.resrad2)
		    )
	    {
	      fX = fXpred;
	      fY = fYpred;
	      nStat = nFindCenter(tMaskParams.pkwid,
				  tMaskParams.move,
				  tMaskParams.lim,
				  tMaskParams.away,
				  poMask, &fX, &fY, &fInten);

	      //+JWP 2009-03-16
	      if ( (0.0 > fY) && (QGOOD == nStat) )
		{
		  cout << "Y0: negative\n";
		  nStat = QBAD;
		}
	      if (QGOOD == nStat)
		{
		  m_tMaskPoints.xo[j][i]     = fX;
		  m_tMaskPoints.yo[j][i]     = fY;
		  m_tMaskPoints.inten[j][i]  = fInten;
		  m_tMaskPoints.status[j][i] = QGOOD;
		  nMissed            = 0;
		}
	      else
		{
		  m_tMaskPoints.xo[j][i]    = fXpred;
		  m_tMaskPoints.yo[j][i]    = fYpred;
		  m_tMaskPoints.inten[j][i] = 0.0;
		  nMissed           = nMissed + 1;
		}
	      if (NULL != m_poOutDump)
		{
		  vDumpPoint (i, j, fX, fY, fInten, nStat);
		  if (QGOOD != nStat)
		    *m_poOutDump << "! Predicted position at: " 
				 << fXpred <<  ' ' << fYpred << endl;
		}

	      if (QGOOD == m_tMaskPoints.status[j][i])
		{
		  nNumPoints = 0;
		  for (jmask = j; jmask <= min(CTR, j+5); jmask++)
		    {
		      if (QGOOD == m_tMaskPoints.status[jmask][i])
			{
			  afXfit[nNumPoints] = m_tMaskPoints.xo[jmask][i];
			  afYfit[nNumPoints] = m_tMaskPoints.yo[jmask][i];
			  afZfit[nNumPoints] = jmask;
			  nNumPoints++;
			}
		    }

		  if (6 <= nNumPoints)
		    nLineFit = 3;
		  else if (4 <= nNumPoints)
		    nLineFit = 2;
		  else if (2 <= nNumPoints)
		    nLineFit = 1;
		  else
		    nLineFit = 0;
              
		  if (2 <= nNumPoints)
		    {
		      (void) nPolyFit (afZfit-1, afXfit-1, afWfit-1, nNumPoints, nLineFit,
                                       aXnum-1, &s2, &r2 );
		      (void) nPolyFit ( afZfit-1, afYfit-1, afWfit-1, nNumPoints, nLineFit,
                                       aYnum-1, &s2, &r2 );
		    }
		}
            
	      j = j - 1;
	      fJ = (float) j;
	      if (1 == nLineFit)
		{
		  fXpred = aXnum[0] + aXnum[1] * fJ;
		  fYpred = aYnum[0] + aYnum[1] * fJ;
		}
	      else if (2 == nLineFit)
		{
		  fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ;
		  fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ;
		}
	      else if (3 == nLineFit)
		{
		  fXpred = aXnum[0] + aXnum[1] * fJ + aXnum[2] * fJ * fJ 
		    + aXnum[3] * fJ * fJ * fJ;
		  fYpred = aYnum[0] + aYnum[1] * fJ + aYnum[2] * fJ * fJ 
		    + aYnum[3] * fJ * fJ * fJ;
		}
	      else
		{
		  fXpred = 2.0 * m_tMaskPoints.xo[j+1][i]  - m_tMaskPoints.xo[j+2][i];
		  fYpred = 2.0 * m_tMaskPoints.yo[j+1][i]  - m_tMaskPoints.yo[j+2][i];
		}
	      if (5 < m_nVerbose)
		cout << "Next predicted X,Y for i,j " << fXpred << ' ' << fYpred 
		     << " i, j: " << i << ' ' << j << endl;
	    }
	  j = j + 1;
	  if ( m_tMaskPoints.yb > j ) m_tMaskPoints.yb = j;
  
	  fXdel = m_tMaskPoints.xo[CTR+1][i] - m_tMaskPoints.xo[CTR][i];
	  fYdel = m_tMaskPoints.yo[CTR+1][i] - m_tMaskPoints.yo[CTR][i];
	} // end while
      //      end do

      // Because of misses, the good peak limits might not be correct
      // Go thru array and make sure correct
      // Transfer tMask info to poPeaks list

      m_tMaskPoints.xb = MASK_SIZE;
      m_tMaskPoints.yb = MASK_SIZE;
      m_tMaskPoints.xe = 0;
      m_tMaskPoints.ye = 0;
      poRefln = new Crefln(poPeaks);

      float afWeight[5] = {1.0, 2.0, 5.0, 2.0, 1.0 };
      float fSumX, fCountX, fSumY, fCountY;
      for (j = 0; j < MASK_SIZE; j++)
	{
	  for (i = 0; i < MASK_SIZE; i++)
	    {
	      if (QGOOD == m_tMaskPoints.status [j][i])
		{
		  if ( i < m_tMaskPoints.xb ) m_tMaskPoints.xb = i;
		  if ( i > m_tMaskPoints.xe ) m_tMaskPoints.xe = i;
		  if ( j < m_tMaskPoints.yb ) m_tMaskPoints.yb = j;
		  if ( j > m_tMaskPoints.ye ) m_tMaskPoints.ye = j;
		  
		  // Smooth the peak positions without messing up original array
		  // just yet.

		  fCountX = 0.0;
		  fCountY = 0.0;
		  fSumX   = 0.0;
		  fSumY   = 0.0;
		  for (k = -2; k < 3; k++)
		    {
		      if ( (0 <= (i + k))  && ((i + k) < MASK_SIZE) )
			{
			  if (QGOOD == m_tMaskPoints.status [j][i+k])
			    {
			      fSumX   += (afWeight[k+2] * m_tMaskPoints.xo[j][i+k]);
			      fCountX += afWeight[k+2];
			    }
			}
		      if ( (0 <= (j + k))  && ((j + k) < MASK_SIZE) )
			{
			  if (QGOOD == m_tMaskPoints.status [j+k][i])
			    {
			      fSumY   += (afWeight[k+2] * m_tMaskPoints.yo[j+k][i]);
			      fCountY += afWeight[k+2];
			    }
			}
		    }
		  if (11.0 == fCountX)
		    fSumX = fSumX / fCountX;
		  else
		    fSumX = m_tMaskPoints.xo[j][i];
		  if (11.0 == fCountY)
		    fSumY = fSumY / fCountY;
		  else
		    fSumY = m_tMaskPoints.yo[j][i];

		  poRefln->vSetH(i);
		  poRefln->vSetK(j);
		  poRefln->vSetL(m_tMaskPoints.status[j][i]);
		  poRefln->vSetField(poPeaks->m_nFI_nNonunfFlag, (int)0);
		  poRefln->vSetIntensity(m_tMaskPoints.inten[j][i]);
		  poRefln->vSetField(poPeaks->m_nFI_fObsPx0, fSumX);
		  poRefln->vSetField(poPeaks->m_nFI_fObsPx1, fSumY);
		  poRefln->vSetField(poPeaks->m_nFI_fCalcPx0, m_tMaskPoints.xo[j][i]);
		  poRefln->vSetField(poPeaks->m_nFI_fCalcPx1, m_tMaskPoints.yo[j][i]);
		  //+JWP 2009-03-13
		  if ( (QGOOD == poRefln->nGetL())
		       && ( (0.0 > fSumX) || (0.0 > fSumY) ) )
		    {
		      cout << "Negative reset in poRefln\n";
		      m_tMaskPoints.status[j][i] = QBAD;
		      poRefln->vSetL(m_tMaskPoints.status[j][i]);
		    }
		  poPeaks->nInsert(poRefln);
		}
	    }
	}
    }

  if ( (NULL != m_poReflnlistMask) && (0 < m_poReflnlistMask->nGetNumReflns()) )
    {
      cout << "number refs: " 
	   << m_poReflnlistMask->nGetNumReflns() << endl << flush;
      m_poReflnlistMask->nWrite("dtcalmask.ref");
    }

  // Copy smoothed peak position back into working array

  for (k = 0; k < poPeaks->nGetNumReflns(); k++)
    {
      poRefln = poPeaks->poGetRefln(k);
      i       = poRefln->nGetH();
      j       = poRefln->nGetK();
      m_tMaskPoints.xo[j][i] = poRefln->fGetField(poPeaks->m_nFI_fObsPx0);
      m_tMaskPoints.yo[j][i] = poRefln->fGetField(poPeaks->m_nFI_fObsPx1);
    }

//  delete poRefln;

  m_tMaskPoints.pts_x = m_tMaskPoints.xe - m_tMaskPoints.yb + 1;
  m_tMaskPoints.pts_y = m_tMaskPoints.ye - m_tMaskPoints.yb + 1;
  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n! Points in x          : " << m_tMaskPoints.pts_x
                   << "\n! Minimum x mask point : " << m_tMaskPoints.xb
                   << "\n! Maximum x mask point : " << m_tMaskPoints.xe
                   << "\n! Points in y          : " << m_tMaskPoints.pts_y
                   << "\n! Minimum y mask point : " << m_tMaskPoints.yb
                   << "\n! Maximum y mask point : " << m_tMaskPoints.ye
                   << endl;
    }

  nStat = QGOOD;
  if (3 < m_nVerbose)
    cout << "...returning from Ccalibrate::nGetMask(...)\n" << endl;
  return (nStat);
}

int
Ccalibrate::nSearchCenter(Cimage *poMask, 
                          tagMask_parameters tMaskParams,
                          float *pfX, float *pfY, float *pfInten)
{
  // Find the maximum value in the box - 
  // start the search for a center from there
  // poMask      -- input  image
  // tMaskParams -- input  mask parameters
  // *pfX        -- modify starting X value, output final X value
  // *pfY        -- modify starting Y value, output final Y value
  // *pfInten    -- Integrated intensity of the found peak

  int   nStat;
  int   i, j, k;
  float fVal;
      
  C3Ddata *po3Ddata;

  po3Ddata = new C3Ddata((int)tMaskParams.dis, (int)tMaskParams.dis, 1, 
                         (int)(*pfX-tMaskParams.dis/2), 
                         (int)(*pfY-tMaskParams.dis/2), 0, e3Ddata_ushort);

  po3Ddata->m_bSpotChase = FALSE;
  nStat = po3Ddata->nFill2D(poMask);

  // TODO:  Any nonunf corrections?

  if (0 == nStat)
    {
      nStat = po3Ddata->nCalcGetMax(&i, &j, &k, &fVal);
      i += po3Ddata->nGetOffset(0);
      j += po3Ddata->nGetOffset(1);
    }

  if (NULL != po3Ddata)
    {
      delete po3Ddata;
      po3Ddata = NULL;
    }

  if (0 == nStat)
    {
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "!   Maximum Value at : " << i << ' ' << j << ' '
                       << fVal << endl;
        }
      if (2 < m_nVerbose)
        cout << "   Maximum Value at : " << i << ' ' << j << ' '
             << fVal << endl;

      *pfX = (float)i;
      *pfY = (float)j;
      nStat = nFindCenter(tMaskParams.pkwid,
                          tMaskParams.move,
                          tMaskParams.lim,
                          tMaskParams.away,
                          poMask, pfX, pfY, pfInten);

      if (2 < m_nVerbose)
        cout << "   Center at : " << *pfX << ' ' << *pfY << ' '
             << *pfInten << endl;

    }
  else
    {
      nStat = QBAD;
    }
  return (nStat);
}
      
int
Ccalibrate::nFindCenter( 
                        const int nSize,   // box size to find center in
                        const float fMove, // max distance point can move
                        const int nMinLim, // min peak value
                        const float fAway, // max dist centroid can be from peak
                        Cimage *poMask,    // Image to look for peaks in
                        float *pfX, float *pfY, // position to find center at
                        float *pfIntensity)    // intensity of peak
{
  // Find the center of mass within a box of "size" at "(x,y)"
  // in "poMask"
  //     
  // Note there is a problem with finding the center until the center doesn't
  // shift.  If the slope change is gradual and the pkwid is small, then it
  // will keep finding the same center, even if it is not correct
  //
  //  nStat is returned as 
  //  1 :: good peak
  //  2 :: The centroid should be within AWAY from the maximum value in the box
  //  4 :: The point should move less than MOVE from the original position
  //  8 :: The maximum pixel is greater than MIN

  int   nStat;
  int   nSearchCount;
  float fXstart, fYstart;
  float fXsave, fYsave;
  C3Ddata *po3DdataI;
  C3Ddata *po3DdataF;

  float a3fCOG[3];
  float fSigmaI;

  nSearchCount = 0;
  po3DdataI    = new C3Ddata(nSize, nSize, 1, 0, 0, 0, e3Ddata_ushort);
  po3DdataF    = new C3Ddata(nSize, nSize, 1, 0, 0, 0);
  po3DdataI->m_bSpotChase = FALSE;
  po3DdataF->m_bSpotChase = FALSE;
  fXsave       = *pfX;
  fYsave       = *pfY;

  float fAvg, fSD;

  do    // while nSearchCount is less than 10 and still need to search
    {
      fXstart      = *pfX;
      fYstart      = *pfY;
      po3DdataI->vSetOffsets((int)fXstart -  nSize / 2, 
                             (int)fYstart -  nSize / 2, 0);
      po3DdataI->vEmpty();
      nStat = po3DdataI->nFill2D(poMask); // Check that vEmpty gets filled
      if (0 == nStat)
        nStat = po3DdataI->nConvertToFloat(poMask, po3DdataF);
//+JWP 2008-09-03
      // Set any pixel value less than or equal to 0 to -999, since -999 is a special bad value in C3Ddata
      
      // po3DdataI->nResetLowValues(0.0, -999.0);

//-JWP
      po3DdataF->vSetSigmaAbove(2.0);
      //po3DdataF->vSetMinPixelsInPeak(1);
      if (0 == nStat)
        {
          // C3Ddata::nCalcGetPeakInfo needs COG in local system
          // (that is, no offsets applied)

          a3fCOG[0] = *pfX - (float)po3DdataF->nGetOffset(0);
          a3fCOG[1] = *pfY - (float)po3DdataF->nGetOffset(1);
          a3fCOG[2] = 0.5;
          nStat = po3DdataF->nCalcGetPeakInfo(&a3fCOG[0], pfIntensity, 
                                              &fSigmaI, &fAvg, &fSD,
					      NULL);
          // Ignore Mask oneedge? errors

          nStat = nStat & ~(  (1 << Crefln::ms_nErrorOnEdge0) 
                            | (1 << Crefln::ms_nErrorOffEdge0) 
                            | (1 << Crefln::ms_nErrorOnEdge1) 
                            | (1 << Crefln::ms_nErrorOffEdge1) 
                            | (1 << Crefln::ms_nErrorOnEdge2)
                            | (1 << Crefln::ms_nErrorOffEdge2)
                            | (1 << Crefln::ms_nErrorUnaccountedPix) 
                            | (1 << Crefln::ms_nErrorTooFewPix)
                           );

	  if (1024 == nStat) nStat = 0;
          if (0 == nStat)
            {
              // Transfer found centroid from local variables

              *pfX = a3fCOG[0] + (float)po3DdataF->nGetOffset(0);
              *pfY = a3fCOG[1] + (float)po3DdataF->nGetOffset(1);
            }
          if (20 < m_nVerbose)
            {
              cout << "Search result: " << nStat << "  "
                   << "Xstart, Ystart: " << fXstart << ", " << fYstart << "\n       "
                   << a3fCOG[0] << ", "
                   << a3fCOG[1] << ", "
                   << *pfIntensity << ", "
                   << fSigmaI << ", "
                   << fAvg << ", "
                   << fSD << ", "
                   << endl << flush;
            }
        }
      nSearchCount++;

      // While the centroid is more than 1/2 of a pixel 
      // from the predicted position
      // then recalculate centroid starting at the just calculated centroid

    } while (   (10 > nSearchCount) && (0 == nStat) 
             && (   (0.5 < fabsf(fXstart-*pfX)) 
                 || (0.5 < fabsf(fYstart-*pfY))) );
      
//+JWP
  if ( (0 == nStat) && ( (0.0 > *pfY) || (0.0 > *pfX) ) )
    {
      // Status is just fine, but centroid is at a negative pixel position
      // which cannot be.  What's going on?  Change the status
      nStat = QBAD;
      cout << "INFO: negative X or Y value with good position, so changed to bad!\n";
      cout << "Search result: " << nStat << "  "
           << "Xstart, Ystart: " << fXstart << ", " << fYstart << "\n       "
           << *pfX << ", "
           << *pfY << ", "
           << *pfIntensity << ", "
           << fSigmaI << ", "
           << fAvg << ", "
           << fSD << ", "
           << endl << flush;
    }
	
//-JWP
  if (0 != nStat)
    {
      // Something bad happened

      if (nSearchCount > 9) cout << "nSearchCount: " << nSearchCount << endl;
      nStat = QBAD;
    }
  else
    {
      // Re-find position of max value in the box

      int nMaxX, nMaxY, nMaxP;
      float fVal;
      nStat = po3DdataF->nCalcGetMax(&nMaxX, &nMaxY, &nMaxP, &fVal);
      if (0 != nStat)
        {
          nStat = QAWAY;
        }
      else if (   (fabsf(*pfX-(float)(nMaxX+po3DdataF->nGetOffset(0))) > fAway)
               || (fabsf(*pfY-(float)(nMaxY+po3DdataF->nGetOffset(1))) > fAway) )
        {
          nStat = QAWAY;
        }
      
      if (   (fabsf(*pfX-fXsave) > fMove)
          || (fabsf(*pfY-fYsave) > fMove) )
        nStat += QMOVE;
      
      if ( *pfIntensity < (float)nMinLim)
        nStat += QLOW;

      if (0 == nStat)
        nStat = QGOOD;   // No errors whatsoever
    }

  if (NULL != po3DdataI)
    {
      delete po3DdataI;
      po3DdataI = NULL;
    }

  if (NULL != po3DdataF)
    {
      delete po3DdataF;
      po3DdataF = NULL;
    }

  return (nStat);
}

void 
Ccalibrate::vDumpString(const Cstring& rsString)
{
  if (NULL != m_poOutDump)
    *m_poOutDump << rsString << endl;
}

void
Ccalibrate::vDumpPoint(const int i, const int j, const float fX,
                       const float fY, const float fInten,
                       const int nStatus)
{
  if (NULL != m_poOutDump)
    {
      *m_poOutDump << " Peak " << setw(5) << i << setw(5) << j
                   << "   " << setw(8) << fX << setw(8) << fY
                   << setw(11) << fInten << setw(5) << nStatus << '\n';
      if (QINTERPOLATED == (nStatus & QINTERPOLATED))
        {
          *m_poOutDump << "! Peak is an interpolated peak\n";
        }
      
      if (QAWAY == (nStatus & QAWAY))
        {
          *m_poOutDump << "! Peak too far from maximum value in box\n";
        }

      if (QMOVE == (nStatus & QMOVE))
        {
          *m_poOutDump << "! Peak moved too far from predicted position\n";
        }

      if (QLOW == (nStatus & QLOW))
        {
          *m_poOutDump << "! Peak intensity below minimum\n";
        }

      if (QREPEAT == (nStatus & QREPEAT))
        {
          *m_poOutDump << "! Duplicate Peak\n";
        }
      m_poOutDump->flush();
    }
}

int
Ccalibrate::nPolyFit(const float *pfX, const float *pfY, const float *pfW,
                     const int nNumPoints, const int nDegree, 
                     float *pfC, float *pfS2, float *pfR2)
// THIS ROUTINE NOT TESTED IN C++ VERSION!
// This routine is translated directly from the Fortran.
// In order to circumvent the array indexing differences between
// Fortran and C++, pass  (pfX-1) so that intern pfX[1] = extern afX[0]
{
  //      ====================================
  //      polyfit ( x, y, w, n, m, c, s2, r2 )
  //      ====================================
  //      This routine finds the coefficients of the nth degree polynomial
  //         y = c(1)+c(2)*x**2+...+c(n+1)*x**n
  //      with a least_square algorithm.
  //
  //      Input :
  //         x, y, w - arrays of points and weight
  //         n - number of points
  //         m - degree of polynomial to fit
  //
  //      Output :
  //         c  - array of coefficients
  //         s2 - residual variance
  //         r2 - coefficient of determination
  //
  //      Returns :
  //        completion status
  //         0 : normal
  //         $ : too many coefficients
  //-----------------------------------------------------------------------

#define PF_LD 1000
  static float d1[PF_LD];
  static float d2[PF_LD];
  static float d3[PF_LD];
  static float d4[PF_LD];
  static float d5[PF_LD];
  static float d6[PF_LD];
  float s1, s2, s3, s4, wt, vr, p, p1, p2;
  float r2;
  int   i, j, k;
      
  if (nNumPoints <= nDegree)
    {
      return (1);
    }

  s1 = 0.0;
  s2 = 0.0;
  s3 = 0.0;
  s4 = 0.0;

  for (i = 1; i <= nNumPoints; i++)
    {
      wt = pfW[i];
      s1 += wt * pfX[i];
      s2 += wt;
      s3 += wt * pfY[i];
      s4 += wt * pfX[i] * pfY[i];
    }
  d1[1] = 0.0;
  d2[1] = 1.0;
  d4[1] = s1 / s2;
  d5[1] = 0.0;
  d6[1] = s3 / s2;
  vr    = s4 -  s3 * d6[1];
      
  for (j = 1; j <= nDegree; j++)
    {
      s1 = 0.0;
      s2 = 0.0;
      s3 = 0.0;
      s4 = 0.0;
      for (i = 1; i <= nNumPoints; i++)
        {
          p1 = 0.0;
          p2 = 1.0;
          for (k = 1; k <= j; k++)
            {
              p = p2;
              p2 = (pfX[i] - d4[k]) * p2  -  d5[k] * p1;
              p1 = p;
            }
          wt = pfW[i];
          p = wt * p2 * p2;
          s1 += p * pfX[i];
          s2 += p;
          s3 += wt * p1 * p1;
          s4 += wt * pfY[i] * p2;
        }

      d4[j+1] = s1 / s2;
      d5[j+1] = s2 / s3;
      d6[j+1] = s4 / s2;
      d3[1] = -d4[j] * d2[1] -  d5[j] * d1[1];
      if (4 <= j)
        {
          for (k = 1; k <= j-2; k++)  // do k = 2, j-2
            {
              d3[k] = d2[k-1] -  d4[j] * d2[k]  -  d5[j] * d1[k];
            }
        }
      if (2 < j)
        d3[j-1] = d2[j-2] -  d4[j] * d2[j-1]  - d5[j];
      if (1 <  j)
        d3[j] = d2[j-1] - d4[j];
      for (k = 1; k <= j; k++)
        {
          d1[k] = d2[k];
          d2[k] = d3[k];
          d6[k] = d6[k]  +  d3[k] * d6[j+1];
        }
    }

  for (j = 1; j <= nDegree+1; j++)
    {
      pfC[j] = d6[j];
    }

  p2 = 0.0;
  for (i = 1; i <= nNumPoints; i++)
    {
      p = pfC[nDegree+1];
      for (j = 1; j <= nDegree; j++)
        {
          p = p * pfX[i] + pfC[nDegree+1-j];
        }
      p = p - pfY[i];
      p2 = p2 + pfW[i] * p * p;
    }

  s2 = 0.0;
  if ( nNumPoints > (nDegree + 1) )
    s2 = p2 / (nNumPoints-nDegree-1);

  r2 = 1.0;
  if (vr != 0.0)
    r2 = 1.0 -  p2 / vr;
  if (r2 < 0.0)
    r2 = 0.0;

  *pfS2 = s2;
  *pfR2 = r2;
  return (0);
}

int
Ccalibrate::nSineFit(const float *pfX, const float *pfY, float *pfW,
                     const int nNumPoints, const int nDegree, 
                     const int nNumIter, const float fCutoff,
                     float *pfC, float *pfS2, float *pfR2)
// THIS ROUTINE NOT TESTED IN C++ VERSION!
// This routine is translated directly from the Fortran.
// In order to circumvent the array indexing differences between
// Fortran and C++, pass  (pfX-1) so that intern pfX[1] = extern afX[0]

{
  //     ==================================================
  //     sinefit ( x, y, w, n, m, c, s2, r2, iter, cutoff )
  //     ==================================================
  //     This routine finds the coefficients of the nth degree polynomial
  //         y = c(1)+c(2)*x**2+...+c(n+1)*x**n
  //     with a weighted least_square algorithm.
  //
  //     Input :
  //         x, y, w - point and weight
  //         n - number of points
  //         m - degree of polynomial to fit
  //         iter - number of iterations
  //         cutoff - weighting factor
  //
  //     Output :
  //         c - coefficients
  //         s2 - residual varience
  //         r2 - coefficient of determination
  //
  //     Returns :
  //        completion status
  //          0 : normal
  //          $ : too many coefficients
  //
  // Notes
  //     1.  The polynomial is first fit with the weighting given on input.
  //         Then it is fit ITER number of iterations after assigning
  //         each point a weight equal to :
  //         w(i) = ( 1 - residual(i)/mean_residual*cutoff) **2
  //
  //     2.  On return w contains the last weighting values used

  //     integer m, n            ! polynomial order, number of points
  //     real    x(n), y(n), w(n)! data, weights
  //     real    c (m+1)         ! polynomial coef
  //     real    s2, r2          ! residual varience, coef determination      
  //     integer iter            ! number of iterations
  //     real    cutoff          ! cutoff weighting

  int   i;
  float fMean;
  int   nStat, nPolyfit;

  //  Get initial fit

  nStat = nPolyFit (pfX, pfY, pfW, nNumPoints, nDegree, pfC, pfS2, pfR2);
  if (0 == nStat)
    {
      for (i = 0; i < nNumIter; i++)
        {
          // Calculate the residuals

          vCalcResiduals(pfX, pfY, nNumPoints, nDegree, pfC, pfW);
            
          // Calculate the average residual

          fMean = fCalcAverage (nNumPoints, pfW);
            
          // Calculate new weighting factors

          vCalcWeights(nNumPoints, fMean, fCutoff, pfW);

          // Fit again

          nStat = nPolyFit ( pfX, pfY, pfW, nNumPoints, nDegree, pfC, pfS2, pfR2);
        }
    }
  return (nStat);
}

void
Ccalibrate::vCalcResiduals(const float *pfX, const float *pfY, 
                           const int nNumPoints, const int nDegree,
                           const float *pfC, float *pfW)
{
//--------------------------------------------------------------------
//
//     =================================
//     calcresiduals (x, y, n, w, c, m )
//     =================================
//     calculate residuals 
//
//          w(i) = abs( c(1) + sum ( c(m+1)*x(m)) - y(i) )
//                            m
//
//---------------------------------------------------------------------
      
  int    i, j;
  float  fT, fS;
      
  for (i = 1; i <= nNumPoints; i++)
    {
      fT = pfC[1];
      fS = 1.0;
      for (j = 1; j <= nDegree; j++)
        {
          fS = fS * pfX[i];
          fT = fT + pfC[j+1] * fS;
        }
      pfW[i] = ABS( pfY[i] - fT );
    }
}

float
Ccalibrate::fCalcAverage(const int nNumPoints, const float *pfW)
{  
//----------------------------------------------------------------------
//
//     ==============
//     average (w, n)
//     ==============
//     calculate the average of w(i)
//
// Returns : average
//
// Arguments : w - real array of values
//             n - number of values
//
//----------------------------------------------------------------------
  int   i;
  float fT;

  fT = 0.0;
  for (i = 1; i <= nNumPoints; i++)
    {
      fT += pfW[i];
    }
  fT = fT / (float) nNumPoints;
  return (fT);
}

void
Ccalibrate::vCalcWeights(const int nNumPoints, const float fAverage, 
                         const float fCutoff, float *pfW)
{     
//----------------------------------------------------------------------
//
//     =====================================
//     calcweights ( w, n, average, cutoff )
//     =====================================
//     calculate weighting factors based on residual
//     factors in w
//
//         w(i) = ( 1 - w(i)/average*cutoff) **2
//
//---------------------------------------------------------------------

  int   i;
  float fWT, fT;
      
  for (i = 1; i <= nNumPoints; i++)
    {
      fWT = pfW[i] / fCutoff / fAverage;
      if (1.0 < fWT)
        {
          pfW[i] = 0.0;
        }
      else
        {
          fT     = 1.0 - fWT * fWT;
          pfW[i] = fT * fT;
        }
    }
  }
      
int
Ccalibrate::nCalcInterpTable(void)
{
  int   nStat = QGOOD;
  float fX_ct_line, fY_ct_line;

  if (DO_RADIAL == m_nradial_distortion)
    {
      if (3 < m_nVerbose)
        {
          cout << "Looking for expected radial distortion." << endl;
        }
      fX_ct_line = 0.0;
      fY_ct_line = 0.0;

      (void) nFindMaskCenter(&m_tMaskPoints, &m_tDistorParams, &m_tMaskParams,
                             &fX_ct_line, &fY_ct_line);
      (void) nFindRadial(&m_tMaskPoints, &m_tDistorParams, &m_tMaskParams);
    }
  else
    {
      if (3 < m_nVerbose)
        {
          cout << "No radial distortion expected." << endl;
        }
      fX_ct_line = (float)CTR;
      fY_ct_line = (float)CTR;
      (void) nFindMaskCenter(&m_tMaskPoints, &m_tDistorParams, &m_tMaskParams,
                             &fX_ct_line, &fY_ct_line);
      m_tDistorParams.a = 1.0; 
      m_tDistorParams.b = 0.0; 
      m_tDistorParams.c = 0.0; 
    }
  if (3 < m_nVerbose)
    {
      cout << "Done Determining Center of Distortion\n"
           << "Mask center found at X,Y: " << fX_ct_line << ", " << fY_ct_line
           << endl;
    }

//  program_status_ ("Calculating Interpolation Tables");

  int nDim0, nDim1;
  int nDivisor0, nDivisor1;

  if (   (NULL == m_poImgMask)
      || !m_poImgMask->bIsAvailable())
    {
      cout << "ERROR in nCalcInterpTable, mask not available!\n" << endl;
      return (-1);
    }
  
  (void)m_poImgMask->nGetDimensions(&nDim0, &nDim1);

  if (600 < nDim0)
    nDivisor0 = 2;
  else 
    nDivisor0 = 2;
  while (0 != (nDim0 % nDivisor0))
    nDivisor0--;
  if (600 < nDim1)
    nDivisor1 = 2;
  else 
    nDivisor1 = 2;
  while (0 != (nDim1 % nDivisor1))
    nDivisor1--;

  cout << "      Mask image size: " << nDim0 << ' ' << nDim1 
       << " Divisors: " << nDivisor0 << ' ' << nDivisor1 << endl;
  //+JWP 2008-12-17
  //cout << "Divisors changed to 1, 1.\n";
  //nDivisor0 = 1;
  //nDivisor1 = 1;
  //-JWP 2008-12-17

  // Create the interpolation tables (as images)
  
  if (NULL != m_poImgDistorXINT)
    delete m_poImgDistorXINT;
  if (NULL != m_poImgDistorXINV)
    delete m_poImgDistorXINV;
  if (NULL != m_poImgDistorYINT)
    delete m_poImgDistorYINT;
  if (NULL != m_poImgDistorYINV)
    delete m_poImgDistorYINV;

  nDim0 = nDim0 / nDivisor0 + nDivisor0;
  nDim1 = nDim1 / nDivisor1 + nDivisor1;

  m_poImgDistorXINT = new Cimage(nDim0, nDim1, eImage_I4);
  m_poImgDistorXINV = new Cimage(nDim0, nDim1, eImage_I4);
  m_poImgDistorYINT = new Cimage(nDim0, nDim1, eImage_I4);
  m_poImgDistorYINV = new Cimage(nDim0, nDim1, eImage_I4);

/*
 * Adjust mask spacing for mask to detector distance
 */

//+JWP 2008-09-03
  //  m_tDistorParams.spacing = m_tMaskParams.mask_spacing
  //                           * (m_fxtod_distance + m_fmasktod_distance) 
  //                             / m_fxtod_distance;


  m_tDistorParams.spacing = m_tMaskParams.mask_spacing 
                            * m_fxtod_distance  / (m_fxtod_distance - m_fmasktod_distance);
//-JWP 2008-09-03

  cout << "\n***   Calculating interpolation tables..." << endl;

  nStat = nCalcInterpolate(&m_tMaskPoints,
                           &m_tDistorParams,
                           &m_tMaskParams,
                           &m_tCorrectParams,
                           m_poImgMask,
                           m_poImgDistorXINT,
                           m_poImgDistorXINV,
                           m_poImgDistorYINT,
                           m_poImgDistorYINV,
                           m_nedges);
  cout << "... Done Calculating Interpolation Tables" << endl;

  m_tDistorParams.x_beam =  m_a2fbeam_position[0];
  m_tDistorParams.y_beam =  m_a2fbeam_position[1];

  (void)m_poImgMask->nGetDimensions(&nDim0, &nDim1);

  m_tDistorParams.x_size = nDim0;
  m_tDistorParams.y_size = nDim1;

  return (nStat);
}

int
Ccalibrate::nFindMaskCenter(tagPoints *ptMask,
                            tagDistort_parameters *ptDistorParams,
                            tagMask_parameters *ptMaskParams,
                            float *pfXcent, float *pfYcent)
{
//------------------------------------------------------------------------
//
//      subroutine findcenter ( mask, dparams, mparams,
//     $     x_ct_line, y_ct_line)
//      
//------------------------------------------------------------------------
//     Determine the
//
//        Center in MASK coordinates    dparams.{x,y}_pt_center
//        Center in pixels              dparams.{x,y}_center
//        Mask angle                    dparams.{ver,horz}_slope
//        Spacing between lines         dparams.{x,y}_scale
//        X/Y ratio                     dparams.ratio
//
//     of the array of points MASK
//
//     Values returned in PARAMS
//------------------------------------------------------------------------

  int i, j;
  int nStat;
  float afX[MASK_SIZE], afY[MASK_SIZE];
  float afW[MASK_SIZE];
  float afC[5], s2, r2;
  float afCV2[5], afCV1[5], afCV0[5];
//      real cv2(5), cv1(5), cv0(5)
  float afCH2[5], afCH1[5], afCH0[5];
//      real ch2(5), ch1(5), ch0(5)
  float fXCenterLine, fYCenterLine;
//      real x_center_line, y_center_line
  float fXCL_slope, fYCL_slope;
//      real xcl_slope, ycl_slope
  float fXCL_int, fYCL_int;
//      real xcl_int, ycl_int
  float fNum;
//      real fnum

  float afHorz[MASK_SIZE][5], afVert[MASK_SIZE][5];
//      real horz (5, MASK_SIZE), VER (5, MASK_SIZE)

  
  char acHorz_line_fit[MASK_SIZE], acVert_line_fit[MASK_SIZE];
//      logical*1 horz_line_fit(MASK_SIZE), ver_line_fit(MASK_SIZE)
  int nNumPoints;

//      character*132 dmpstr


  for (i = 0; i < MASK_SIZE; i++) afW[i] = 1.0;

//  Find the position where the straightest line would fall.
//  The straightest line would have 1st and 2nd order terms equal to 0
//  and would be along a position equal to their constant term.
//  To do this plot: the 1st and 2nd order terms as a function of the
//  constant term.  The intercept of both of these plots should
//  be the center.  If the the mask is tilted the 1st order term
//  will not be 0 - use the first order term to determine the
//  angle of the mask

  for (i = 0; i < MASK_SIZE; i++)
    {
      acVert_line_fit[i] = (char)0;
      acHorz_line_fit[i] = (char)0;
    }
      
// Find equations to describe all the horizontal and vertical lines:

  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n!"
                   << "\n! Fitting horizontal lines"
                   << "\n! ========================"
                   << "\n!          Line   2nd (curve)     1st      "
                   << "      const         (slope)   Pts    varience"
                   << "\n!          ====   ===========   ==========="
                   << "   ===========    =========== ===  ============"
                   << endl;
    }

  for (j = ptMask->yb; j <= ptMask->ye; j++)
    {
      nNumPoints = 0;
      for (i = ptMask->xb; i <= ptMask->xe; i++)
        {
          if (QGOOD == ptMask->status[j][i])
            {
              afX[nNumPoints] = ptMask->xo[j][i];
              afY[nNumPoints] = ptMask->yo[j][i];
              afW[nNumPoints] = 1.0;
              nNumPoints++;
            }
        }
      if (3 <= nNumPoints)
        {
          acHorz_line_fit[j] = (char)1;
          if (DO_RADIAL == ptMaskParams->radial)
            {
              (void) nPolyFit(afX-1, afY-1, afW-1, nNumPoints, 2,
                              &(afHorz[j][0])-1, &(afHorz[j][3]), &(afHorz[j][4]));
//               call polyfit ( x, y, w, nNumPoints, 2, 
//     $              horz(1,j), horz(4,j), horz(5,j) )
            }
          else
            {
              (void) nPolyFit(afX-1, afY-1, afW-1, nNumPoints, 1,
                              &(afHorz[j][0])-1, &(afHorz[j][3]), 
                              &(afHorz[j][4]));
//               call polyfit ( x, y, w, nNumPoints, 1, 
//     $              horz(1,j), horz(4,j), horz(5,j) )
              afHorz[j][2] = 0.0;
//               horz(3,j) = 0
            }
          fNum = afHorz[j][1] +  2.0 * (float)CTR * afHorz[j][2];
          if (NULL != m_poOutDump)
            {
              *m_poOutDump << "HORZ_LINE   " << setw(3) << j 
                        << setw(14) << afHorz[j][2] 
                        << setw(14) << afHorz[j][1] 
                        << setw(14) << afHorz[j][0]
                        << setw(14) << fNum << setw(4) << nNumPoints
                        << setw(14) << afHorz[j][3] << endl;
//              write (dmpstr,602) j, horz(3,j), horz(2,j), horz(1,j),
//     $              fNum, nNumPoints, horz(4,j)
//               call dumpstring (dmpstr)
// 602           format ('HORZ_LINE   ',i4,2e14.4,f14.4,e14.4,i4,e14.4)
            }
        }
    }
      
  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n!"
                   << "\n! Fitting vertical lines"
                   << "\n!========================"
                   << "\n!          Line   2nd (curve)     1st      "
                   << "      const         (slope)   Pts    varience"
                   << "\n!          ====   ===========   ==========="
                   << "   ===========    =========== ===  ============"
                   << endl;
    }

  for (i = ptMask->xb; i <= ptMask->xe; i++)
    {
      nNumPoints = 0;
      for (j = ptMask->yb; j <= ptMask->ye; j++)
        {
          if (QGOOD == ptMask->status[j][i])
            {
              afX[nNumPoints] = ptMask->yo[j][i];
              afY[nNumPoints] = ptMask->xo[j][i];
              afW[nNumPoints] = 1.0;
              nNumPoints++;
            }
        }
      if (3 <= nNumPoints)
        {
          acVert_line_fit[i] = (char)1;
          if (DO_RADIAL == ptMaskParams->radial)
            {
              (void) nPolyFit(afX-1, afY-1, afW-1, nNumPoints, 2,
                              &(afVert[i][0])-1, &(afVert[i][3]), &(afVert[i][4]));
//               call polyfit ( x, y, w, nNumPoints, 2, 
//     $              ver(1,i), ver(4,i), ver(5,i) )
            }
          else
            {  
              (void) nPolyFit(afX-1, afY-1, afW-1, nNumPoints, 1,
                              &(afVert[i][0])-1, &(afVert[i][3]), 
                              &(afVert[i][4]));
//               call polyfit ( x, y, w, nNumPoints, 1, 
//     $              ver(1,i), ver(4,i), ver(5,i) )
              afVert[i][2] = 0.0;
            }
          fNum = afVert[i][1]  +  2.0 * CTR * afVert[i][2];
          if (NULL != m_poOutDump)
            {
              *m_poOutDump << "VERT_LINE   " << setw(3) << i 
                        << setw(14) << afVert[i][2] 
                        << setw(14) << afVert[i][1] 
                        << setw(14) << afVert[i][0]
                        << setw(14) << fNum << setw(4) << nNumPoints
                        << setw(14) << afVert[i][3] << endl;
//               write (dmpstr,603) i, ver(3,i), ver(2,i), ver(1,i), 
//     $              fNum, nNumPoints, ver(4,i)
// 603           format ('VERT_LINE   ',i4,2e14.4,f14.4,e14.4,i4,e14.4)
//               call dumpstring (dmpstr)
            }
        }
    }
      
//     Fit the 2nd, 1st and constant terms of horz and vertical lines
//     as function of line number
//-----------------------------------------------------------------------------
 
  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n! Fitting the 2nd, 1st and constant terms of horz"
                   << "\n! and vertical lines as function of line number"
                   << "\n! ================================================="
                   << endl;
    }

//    Fit second order of vertical
//    note that there aren't any second order terms if no radial distortion
//    =====================================================================

  if (DO_RADIAL == ptMaskParams->radial)
    {
      nNumPoints = 0;
      for (i = ptMask->xb; i <= ptMask->xe; i++)
        {
          if ((char)1 == acVert_line_fit[i])
            {
              afY[nNumPoints] = afVert[i][2]; // ver(3,i);
              afX[nNumPoints] = (float)i;
              afW[nNumPoints] = 1.0;
              nNumPoints++;
            }
        }
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 1, 
                        5, 3.0, afCV2-1, &s2, &r2);

      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! 2nd order vertical   coef : " 
                       << afCV2[0] << ' '
                       << afCV2[1] << endl;
        }
    }
      
//    Fit first order of vertical
//    ==============================

  nNumPoints = 0;
  for (i = ptMask->xb; i <= ptMask->xe; i++)
    {
      if ((char)1 == acVert_line_fit[i])
        {
          afY[nNumPoints] = afVert[i][1];
          afX[nNumPoints] = (float) i;
          afW[nNumPoints] = 1.0;
          nNumPoints++;
        }
    }
  if (DO_RADIAL == ptMaskParams->radial)
    {
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 1, 
                        5, 3.0, afCV1-1, &s2, &r2);
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! 1st order vertical   coef : "
                       << afCV1[0] << ' '
                       << afCV1[1] << endl;
        }
    }
  else
    {
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 4, 
                        5, 3.0, afCV1-1, &s2, &r2);
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! 1st order vertical   coef : " << afCV1[0] << endl;
        }
    }

  //  Fit constant of vertical
  //  ========================
  
  nNumPoints = 0;
  for (i = ptMask->xb; i <= ptMask->xe; i++)
    {
      if ((char)1 == acVert_line_fit[i])
        {
          afY[nNumPoints] = afVert[i][0]; // (1,i);
          afX[nNumPoints] = (float)i;
          afW[nNumPoints] = 1.0;
          nNumPoints++;
        }
    }

  if (DO_RADIAL == ptMaskParams->radial)
    {
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 1, 
                        5, 3.0, afCV0-1, &s2, &r2);
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! Constant  vertical   coef : "
                       << afCV0[0] << ' '
                       << afCV0[1] << endl;
        }
    }
  else
    {
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 4, 
                        5, 3.0, afCV0-1, &s2, &r2);
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! Constant  vertical   coef : " << afCV0[0] << endl;
        }
    }

  // Fit second order of horizontal
  // Note that there aren't any second order terms if no radial distortion
  // =====================================================================

  if (DO_RADIAL == ptMaskParams->radial)
    {
      nNumPoints = 0;
      for (i = ptMask->yb; i <= ptMask->ye; i++)
        {
          if ((char)1 == acHorz_line_fit[i])
            {
              afY[nNumPoints] = afHorz[i][2];
              afX[nNumPoints] = (float)i;
              afW[nNumPoints] = 1.0;
              nNumPoints++;
            }
        }
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 1, 
                        5, 3.0, afCH2-1, &s2, &r2);
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! 2nd order horizontal coef : "
                       << afCH2[0] << ' '
                       << afCH2[1] << endl;
        }
    }

  // Fit first order of horizontal
  // ==============================

  nNumPoints = 0;
  for (i = ptMask->yb; i <= ptMask->ye; i++)
    {
      if ((char)1 == acHorz_line_fit[i])
        {
          afY[nNumPoints] = afHorz[i][1]; // (2,i);
          afX[nNumPoints] = (float)i;
          afW[nNumPoints] = 1.0;
          nNumPoints++;
        }
    }
  if (DO_RADIAL == ptMaskParams->radial)
    {
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 1, 
                        5, 3.0, afCH1-1, &s2, &r2);
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! 1st order horizontal coef : "
                       << afCH1[0] << ' '
                       << afCH1[1] << endl;
        }
    }
  else
    {
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 4, 
                        5, 3.0, afCH1-1, &s2, &r2);
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! 1st order horizontal coef : " << afCH1[0] << endl;
        }
    }

  // Fit constant of horizontal
  // ========================

  nNumPoints = 0;
  for (i = ptMask->yb; i <= ptMask->ye; i++)
    {
      if ((char)1 == acHorz_line_fit[i])
        {
          afY[nNumPoints] = afHorz[i][0]; // (1,i)
          afX[nNumPoints] = (float)i;
          afW[nNumPoints] = 1.0;
          nNumPoints++;
        }
    }
  if (DO_RADIAL == ptMaskParams->radial)
    {
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 1, 
                        5, 3.0, afCH0-1, &s2, &r2);
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! Constant  horizontal coef : "
                       << afCH0[0] << ' '
                       << afCH0[1] << endl;
        }
    }
  else
    {
      nStat = nSineFit (afX-1, afY-1, afW-1, nNumPoints, 4, 
                        5, 3.0, afCH0-1, &s2, &r2);
      if (NULL != m_poOutDump)
        {
          *m_poOutDump << "! Constant  horizontal coef : " << afCH0[0] << endl;
        }
    }
      
//=========================================================================
//    Determine Center Line Number
//=========================================================================

  if ( (0.0 < *pfXcent) && (0.0 < *pfYcent) )
    {
      //    If given the center line, use it

      fXCenterLine = *pfXcent;
      fYCenterLine = *pfYcent;
    }
  else if (DO_RADIAL == ptMaskParams->radial)
    {
      //    If radial distortion, the center line should
      //    be the one with no second order term

      fXCenterLine =  - afCV2[0] / afCV2[1];
      fYCenterLine =  - afCH2[0] / afCH2[1];
    }
  else
    {
      //    Assume that the center line is the center line

      fXCenterLine = (float)CTR;
      fYCenterLine = (float)CTR;
    }

  ptDistorParams->x_pt_center = fXCenterLine;
  ptDistorParams->y_pt_center = fYCenterLine;

  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n! x center line " << fXCenterLine
                   << "\n! y center line " << fYCenterLine << endl;
    }
      
  //=========================================================================
  //    Determine Center Line Slope
  //=========================================================================

  if (DO_RADIAL == ptMaskParams->radial)
    {
      fXCL_slope = afCV1[0] +  afCV1[1] * fXCenterLine;
      fYCL_slope = afCH1[0] +  afCH1[1] * fYCenterLine;
    }
  else
    {
      fXCL_slope = afCV1[0]
                 + afCV1[1] * fXCenterLine
                 + afCV1[2] * fXCenterLine * fXCenterLine
                 + afCV1[3] * fXCenterLine * fXCenterLine
                            * fXCenterLine
                 + afCV1[4] * fXCenterLine * fXCenterLine
                            * fXCenterLine * fXCenterLine;

      fYCL_slope = afCH1[0]
                 + afCH1[1] * fYCenterLine
                 + afCH1[2] * fYCenterLine * fYCenterLine
                 + afCH1[3] * fYCenterLine * fYCenterLine
                            * fYCenterLine
                 + afCH1[4] * fYCenterLine * fYCenterLine
                            * fYCenterLine * fYCenterLine;
    }

  ptDistorParams->ver_slope  = fXCL_slope;
  ptDistorParams->horz_slope = fYCL_slope;
  if (NULL != m_poOutDump)
    {
      float fAngle;
      fAngle = atanf(ptDistorParams->ver_slope) / Gs_fRADIANS_PER_DEGREE;
      *m_poOutDump << "\n! Vertical  line slope, angle : " 
                   << ptDistorParams->ver_slope << ' ' << fAngle;
      fAngle = atanf(ptDistorParams->horz_slope)  / Gs_fRADIANS_PER_DEGREE;
      *m_poOutDump << "\n! Horzontal line slope, angle : " 
                   << ptDistorParams->horz_slope << ' ' << fAngle << endl;
    }
      
  //=========================================================================
  //    Determine Center Line Intercept
  //=========================================================================

  if (DO_RADIAL == ptMaskParams->radial)
    {
      //    With Radial Distortion the center line has no
      //    second order term

      fXCL_int = afCV0[0] + afCV0[1] * fXCenterLine;
      fYCL_int = afCH0[0] + afCH0[1] * fYCenterLine;
    }
  else
    {
      //    Without radial distortion, no second order terms
      //    were fit to the lines

      fXCL_int   = afCV0[0]
                 + afCV0[1] * fXCenterLine
                 + afCV0[2] * fXCenterLine * fXCenterLine
                 + afCV0[3] * fXCenterLine * fXCenterLine
                            * fXCenterLine
                 + afCV0[4] * fXCenterLine * fXCenterLine
                            * fXCenterLine * fXCenterLine;

      fYCL_int   = afCH0[0]
                 + afCH0[1] * fYCenterLine
                 + afCH0[2] * fYCenterLine * fYCenterLine
                 + afCH0[3] * fYCenterLine * fYCenterLine
                            * fYCenterLine
                 + afCH0[4] * fYCenterLine * fYCenterLine
                            * fYCenterLine * fYCenterLine;

    }
      
  //=====================================================================
  //    Determine the center of distortion as the
  //    intercept of the two center lines
  //=====================================================================

  ptDistorParams->x_center = (fXCL_slope * fYCL_int  +  fXCL_int) 
                             / (1.0 - fXCL_slope *fYCL_slope);
  ptDistorParams->y_center = fYCL_slope * ptDistorParams->x_center + fYCL_int;

  //    Note:
  //    the calculated value of x_center, y_center should be close
  //    to the coordinates of the center peak if no there is
  //    no radial distortion.  It can
  //    easily get thrown off if the spacing of the mask points 
  //    around the center is not correct.  This happened to me when 
  //    trying to understand an image taken with the gold detector where
  //    the center lines were missing due to a problem with the 
  //    electronics - the calculated center was several pixels away from
  //    the center peak.  I went thru the calculations carefully and 
  //    assured myself the the problem came from the missing lines.
  //
  //    The only thing that I might change is the order of the fit
  //    for the const and 1st order terms - a 4th order poly seems
  //    excessive for an essentially straight line.
  //=====================================================================
  //    Determine the Spacing between lines
  //=====================================================================
  //     
  //    Fit a curve through the centers of the horizontal
  //    and vertical lines.  Then solve for the center position
  //    and the slope of the curve through that point
  //     
  //=====================================================================

  float fT;

  fT         = ptDistorParams->y_center;
  nNumPoints = 0;

  for (i = ptMask->xb; i <= ptMask->xe; i++)
    {
      if ((char) 1 == acVert_line_fit[i])
        {
          afX[nNumPoints] = (float)i;
//          afY[nNumPoints] = ver(3,i)*t*t + ver(2,i)*t + ver(1,i);
          afY[nNumPoints] = afVert[i][2] * fT * fT + afVert[i][1] * fT 
                           + afVert[i][0];
          afW[nNumPoints] = 1.0;
          nNumPoints++;
        }
    }

  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n!"
                   << "\n! Fitting for x-spacing"
                   << "\n! ====================="
                   << "\n!       line      center"
                   << "\n!       ====      ======";
      for (i = 0; i < nNumPoints; i++)
        *m_poOutDump << "\nX_SPACE " << setw(4) << afX[i] << ' ' << setw(10) << afY[i];
      *m_poOutDump << endl;
    }
  (void) nPolyFit(afX-1, afY-1, afW-1, nNumPoints, 3, afC-1, &s2, &r2);
//      call polyfit ( x, y, w, nNumPoints, 3, c, s2, r2 )
  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n!"
                   << "\n! Fit third order polynomial:"
                   << "\n! ==========================="
                   << "\n! a, b, c, d = " << setw(10) << afC[3] << ' ' << setw(10) << afC[2]
                                   << ' ' << setw(10) << afC[1] << ' ' << setw(10) << afC[0]
                   << "\n!" << endl;
    }

  ptDistorParams->x_scale = 3.0 * afC[3] 
                     * ptDistorParams->x_pt_center * ptDistorParams->x_pt_center
                   + 2.0 * afC[2] * ptDistorParams->x_pt_center 
                   + afC[1];
  ptDistorParams->x_scale = ptDistorParams->x_scale 
    / sqrtf(1.0 + ptDistorParams->horz_slope * ptDistorParams->horz_slope);

  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "! X spacing      : " << ptDistorParams->x_scale << endl;
    }

  fT = ptDistorParams->x_center;
  nNumPoints = 0;
  for (i = ptMask->yb; i <= ptMask->ye; i++)
    {
      if ((char)1 == acHorz_line_fit[i])
        {
          afX[nNumPoints] = (float)i;
          afY[nNumPoints] = afHorz[i][2] * fT * fT + afHorz[i][1] * fT 
                            + afHorz[i][0];
          afW[nNumPoints] = 1.0;
          nNumPoints++;
        }
    }
  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n!"
                   << "\n! Fitting for y-spacing"
                   << "\n! ====================="
                   << "\n!       line      center"
                   << "\n!       ====      ======";
      for (i = 0; i < nNumPoints; i++)
        *m_poOutDump << "\nY_SPACE " << setw(4) << afX[i] << ' ' << setw(10) << afY[i];
      *m_poOutDump << endl;
    }

  (void) nPolyFit(afX-1, afY-1, afW-1, nNumPoints, 3, afC-1, &s2, &r2);
//      call polyfit ( x, y, w, nNumPoints, 3, c, s2, r2 )
  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n!"
                   << "\n! Fit third order polynomial:"
                   << "\n! ==========================="
                   << "\n! a, b, c, d = " << setw(10) << afC[3] << ' ' << setw(10) << afC[2]
                                   << ' ' << setw(10) << afC[1] << ' ' << setw(10) << afC[0]
                   << "\n!" << endl;
    }
  ptDistorParams->y_scale = 3.0 * afC[3]
                     * ptDistorParams->y_pt_center * ptDistorParams->y_pt_center
                   + 2.0 * afC[2] * ptDistorParams->y_pt_center 
                   + afC[1];
  ptDistorParams->y_scale = ptDistorParams->y_scale 
    / sqrtf(1.0 + ptDistorParams->ver_slope * ptDistorParams->ver_slope);

  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "! Y spacing      : " << ptDistorParams->y_scale << endl;
    }
      
  ptDistorParams->ratio = ptDistorParams->y_scale / ptDistorParams->x_scale;
  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "! Y/X Ratio              : " << ptDistorParams->ratio
                   << "\n!" << endl;
    }

  //    All done - write an exit at the end of the dumpstring
  //    for programs which might want to read it

  vDumpString ("EXIT");
  return (0);
}

int
Ccalibrate::nFindRadial(tagPoints *ptMaskPoints,
                        tagDistort_parameters *ptDistorParams,
                        tagMask_parameters *ptMaskParams)
{
/*
C------------------------------------------------------------------------
      
      subroutine findradial ( mask, dparams, mparams )
      
C------------------------------------------------------------------------
C
C     Find the radial distortion
C     Fit to the equation:     2
C             r = R(a + bR + cR )
C     
C     where   r = expected or corrected position ( = e )
C             R = observed position ( = o )
C     
C     The values for the observed positions are those determined
C     by the mask peak finding routine GetMASK
C
C     The expected positions are calculated based on the
C     average spacing dparams.{x,y}_scale, the mask angle 
C     dparams.{horz,vert}_slope, and the index of the point
C     There is no value in mask.{x,y}e at this
C     point, and there still isn't after this routine is complete.
C
C------------------------------------------------------------------------
*/
  int i, j;
  double dXcenter, dYcenter;
  double xe, ye, xo, yo;
  double e, o;
  double s, c;
  double o6, o5, o4, o3, o2;
  double o3e, o2e, oe;
  double del, del1, del2, del3;
  double xx, xy, yy, yx;
  float  fX, fY;

  float sinh, cosh, sinv, cosv;

  float fXscale, fYscale;
  float fXpt, fYpt;

  del = sqrt(1.0 + (double)ptDistorParams->horz_slope 
                 * (double)ptDistorParams->horz_slope);
  sinh = ptDistorParams->horz_slope / del;
  cosh = 1.0 / del;
  del = sqrt(1.0 + (double)ptDistorParams->ver_slope
                 * (double)ptDistorParams->ver_slope);
  sinv = ptDistorParams->ver_slope / del;
  cosv = 1.0 / del;
  dXcenter = (double)ptDistorParams->x_center;
  dYcenter = (double)ptDistorParams->y_center;
  fXscale  = ptDistorParams->x_scale;
  fYscale  = ptDistorParams->y_scale;
  fXpt     = ptDistorParams->x_pt_center;
  fYpt     = ptDistorParams->y_pt_center;
  o6       = 0.0;
  o5       = 0.0;
  o4       = 0.0;
  o3       = 0.0;
  o2       = 0.0;
  o3e      = 0.0;
  o2e      = 0.0;
  oe       = 0.0;

  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n! Calculating Radial distortion"
                   << "\n! Observed    Expected"
                   << "\n! ========    ========" << endl;
    }
      
  for (j = ptMaskPoints->yb; j <= ptMaskPoints->ye; j++)
    {
      fY = ((float)j - fYpt) * fYscale;
      xy = (double)fY * sinv;
      yy = (double)fY * cosv;
      for (i = ptMaskPoints->xb; i <= ptMaskPoints->xe; i++)
        {
          if (QGOOD == ptMaskPoints->status[j][i])
            {
              fX = ((float)i - fXpt) * fXscale;
              xx = (double)fX * cosh;
              yx = (double)fX * sinh;
              ye = yy + yx;
              xe = xy + xx;
              yo = (double)ptMaskPoints->yo[j][i] - dYcenter;
              xo = (double)ptMaskPoints->xo[j][i] - dXcenter;
              e = sqrt(xe*xe + ye*ye);
              o = sqrt(xo*xo + yo*yo);
              if (NULL != m_poOutDump)
                {
                  *m_poOutDump << o << ' ' << e << endl;
                }
              s = o*o;
              c = s*o;
               
              o6  = o6  + c*c;
              o5  = o5  + s*c;
              o4  = o4  + s*s;
              o3  = o3  + c;
              o2  = o2  + s;
              o3e = o3e + c*e;
              o2e = o2e + s*e;
              oe  = oe  + o*e;
            }
        }
    } 
      
/*
C      del = o2*(o4*o6-o5*o5) +  o3*(o5*o4-o6*o3) +   o4*(o3*o5-o4*o4)
C      del1= oe*(o4*o6-o5*o5) +  o3*(o5*o3e-o6*o2e) + o4*(o2e*o5-o3e*o4)
C      del2= o2*(o2e*o6-o3e*o5)+ oe*(o5*o4-o6*o3) +   o4*(o3*o3e-o4*o2e)
C      del3= o2*(o4*o3e-o5*o2e) +o3*(o2e*o4-o3e*o3) + oe*(o3*o5-o4*o4)
C
C     MS Jan 92 - Make sure this doesn't floating overflow
C
*/
  del = (o6/o3-o5/o4*o5/o3)
          + (o5/o2-o6/o4*o3/o2) 
          + (o5/o2-o4/o3*o4/o2);
  del1= oe/o2*(o6/o3-o5/o4*o5/o3)
          + (o5/o2*o3e/o4-o6/o4*o2e/o2)
          + (o2e/o3*o5/o2-o3e/o2*o4/o3);
  del2= (o2e/o3*o6/o4-o3e/o4*o5/o3)
          + oe/o2*(o5/o3-o6/o4)
          + (o3e/o2-o2e/o3*o4/o2);
  del3= (o3e/o3-o5/o4*o2e/o3)
          +(o2e/o2-o3e/o4*o3/o2)
          + oe/o2*(o5/o4-o4/o3);

  ptDistorParams->a = del1 / del;
  ptDistorParams->b = del2 / del;
  ptDistorParams->c = del3 / del;
  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "\n! Distortion parameter :"
                   << "\n!                      A : " << ptDistorParams->a
                   << "\n!                      B : " << ptDistorParams->b
                   << "\n!                      C : " << ptDistorParams->c
                   << "\n!" << endl;
    }
  return (0);
}

int
Ccalibrate::nCalcInterpolate(tagPoints             *ptMaskPoints,
                             tagDistort_parameters *ptDistorParams,
                             tagMask_parameters    *ptMaskParams,
                             tagCorrect_parameters *ptCorrectParams,
                             Cimage *poImgMask,
                             Cimage *poImgDistorXINT,
                             Cimage *poImgDistorXINV,
                             Cimage *poImgDistorYINT,
                             Cimage *poImgDistorYINV,
                             int nEdges)
{
/*      
C===========================================================================
C     Create an interpolation tables for final spatial distortion 
C     corrections
C     
C     Arguments :
C     mask - 
C            INPUT : mask record containing observed positions of
C            all observed mask points
C
C            OUTPUT : calculated, unpinned and interpolated
C            positions of all observed mask points.
C     
C     dparams - calculated values describing distortions
C        Center in MASK coordinates    dparams.{x,y}_pt_center
C        Center in pixels              dparams.{x,y}_center
C        Mask angle                    dparams.{ver,horz}_slope
C        Spacing between lines         dparams.{x,y}_scale
C        X/Y ratio                     dparams.ratio
C        Pincushion values             dparams.{a,b,c}
C        (Mask) Spacing                dparams.spacing
C     
C     mparams - parameters describing the mask
C        Pixel Size                    mparams.pixel_size
C        Image limits                  mparams.slit
C                                      mparams.window
C                                      mparams.search_center
C                                      mparams.resrad2
C
C     imx, imy - 
C        Size of full image.
C
C     x_int, y_int, x_inv, y_inv
C        Arrays to contain interpolation tables
C
C     axsize, aysize
C        Size of interpolation tables
C
C     edges
C        If edges = 1, worry about the edges.  This is not necessary
C        on older detectors (fast, siemens, blue, sit...) but is 
C        necessary with the modular CCD detectors
C
C     istat
C        Completion status
C     
C     
C===========================================================================
*/

  int   i,j, ii, jj, iii, jjj, k;      // counters
  int   timx, timy;                    // size of int tables
  int   xb, xe, yb, ye;                // mask start, end
  int   xbi, xei, ybi, yei;            // mask start, end
  float x, y;                          // pixel position
  float xm, ym;                        // x, y in mask coordinates
  float xu, yu;                        // x, y in unpinned coordinates
  float xi, yi;                        //  x, y interpolated value

  float afstempx[4][4];                // temp array for splines
  float afstempy[4][4];                // temp array for splines
      
  int   axsize, aysize;
      
  int   xint_step, yint_step, xinv_step, yinv_step;
  int   xint_start, yint_start, xinv_start, yinv_start;
  float pscale;

  const float  RFLAG = 1.0E30;

  int   indexi, indexj;
  float xmin, xmax, ymin, ymax;
  int   nxmin, nxmax, nymin, nymax;
  float x1, x2, x3, x4, y1, y2, y3, y4;

//  integer   inbox

  float max_radius, max_radiusi;
  int   radius;
  int   xcenter, ycenter;
  int   step;
//  float b_spline2;
  float a, b, c, m, xp, yp;
  int   nStat;

  float tr2, tri2, ty2, tyi2;
  int   badsq, badrad, ngood, goodsq, goodrad;
  int   shape;

  const int SQUARE = 1;
  const int CIRCLE = 2;

  int   anclimits[4];

//  byte      info(8192);
//  integer*2 dx(MASK_SIZE, MASK_SIZE), dy(MASK_SIZE, MASK_SIZE);

  float minxo, maxxo, minyo, maxyo;
  float minxc, maxxc, minyc, maxyc;

  nStat = 0;
        
  // Set values of all pixels in the interpolation arrays (images) to 
  // ptCorrectParams->badflag.

  int nDim;
  int nDim0, nDim1;
  LONG lIntBadFlag;

  lIntBadFlag = (LONG) ptCorrectParams->badflag;
  nDim = poImgDistorXINT->nGetDimension();
  (void) poImgDistorXINT->nGetDimensions(&axsize, &aysize);
  (void) poImgMask->nGetDimensions(&nDim0, &nDim1);
  poImgDistorXINT->nSetNextPixel(0,0);
  poImgDistorXINV->nSetNextPixel(0,0);
  poImgDistorYINT->nSetNextPixel(0,0);
  poImgDistorYINV->nSetNextPixel(0,0);
  for (i = 0; i < nDim; i++)
    {
      poImgDistorXINT->vSetNextPixel(lIntBadFlag);
      poImgDistorXINV->vSetNextPixel(lIntBadFlag);
      poImgDistorYINT->vSetNextPixel(lIntBadFlag);
      poImgDistorYINV->vSetNextPixel(lIntBadFlag);
    }
/*      
C     ========================================================
C     The first goal is to calculate the difference between
C     the calculated correct position for each mask point
C     and the 'unpinned' position.
C     ========================================================

C     Fill mask with calculated 'correct' (c) values based on pixel size,
C     mask spacing and peak index
*/

  vCalcExpected(ptMaskPoints, ptDistorParams, ptMaskParams, ptCorrectParams);

/*
C     Calculate 'unpinned' (u) peak position after correction for
C     pixel size and pin-cushion distortion
*/

  vCalcUnpinned(ptMaskPoints, ptDistorParams, ptMaskParams);

/*      
C     Fill the difference arrays with values for all the 
C     good points with the difference between unpinned (u)
C     and the calculated correct position (c)
C     
C     Mark all bad points with FLAG, and QINTERPOLATED
C     =============================================================
*/

  xb = MASK_SIZE-1;
  yb = MASK_SIZE-1;
  xe = 0;
  ye = 0;
  max_radius = 0.0;
  xcenter = (xb + xe) / 2;
  ycenter = (yb + ye) / 2;

  for (j = 0; j < MASK_SIZE; j++)
    {
      for (i = 0; i < MASK_SIZE; i++)
        {
          if (QGOOD == ptMaskPoints->status [j][i])
            {
              ptMaskPoints->x_diff[j][i] = ptMaskPoints->xu[j][i] 
                                            - ptMaskPoints->xc[j][i];
              ptMaskPoints->y_diff[j][i] = ptMaskPoints->yu[j][i]
                                            - ptMaskPoints->yc[j][i];
              radius = (i - xcenter) * (i - xcenter)
                     + (j - ycenter) * (j - ycenter);
              if ( radius > max_radius ) max_radius = radius;
              if ( i < xb ) xb = i;
              if ( i > xe ) xe = i;
              if ( j < yb ) yb = j;
              if ( j > ye ) ye = j;
            }
          else
            {
              ptMaskPoints->x_diff[j][i] = RFLAG;
              ptMaskPoints->y_diff[j][i] = RFLAG;
              ptMaskPoints->status[j][i] = QINTERPOLATED;
            }
        }
    }
/*
CDEBUG++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$      do i = 1, 8192
c$$$         info(i) = 0
c$$$      end do
c$$$      do j = 1, MASK_SIZE
c$$$         do i = 1, MASK_SIZE
c$$$            if ( ptMaskPoints->x_diff [j][i] .eq. RFLAG ) then
c$$$               dx[j][i] = 1000
c$$$               dy[j][i] = 1000
c$$$            else
c$$$               dx[j][i] = ptMaskPoints->x_diff[j][i]*10 + 1000
c$$$               dy[j][i] = ptMaskPoints->y_diff[j][i]*10 + 1000
c$$$            end if
c$$$         end do
c$$$      end do
c$$$      call clrhd (info)
c$$$      call wrmad('dx.mad', info, 
c$$$     $     dx, MASK_SIZE, MASK_SIZE, MASK_SIZE, MASK_SIZE,
c$$$     $     istat)
c$$$
c$$$      call clrhd (info)
c$$$      call wrmad('dy.mad', info, 
c$$$     $     dy, MASK_SIZE, MASK_SIZE, MASK_SIZE, MASK_SIZE,
c$$$     $     istat)
CDEBUG--------------------------------------------------
*/

/*
C     Is this a square or a round array ?
*/

  badrad  = 0;
  badsq   = 0;
  goodrad = 0;
  goodsq  = 0;
  ngood   = 0;
  for (j = 0; j < MASK_SIZE; j++)
    {
      ty2 = ((float)j - ycenter) * ((float)j - ycenter);
      for (i = 0; i < MASK_SIZE; i++)
        {
          tr2 = ((float)i - xcenter) * ((float)i - xcenter) + ty2;
          if (RFLAG == ptMaskPoints->x_diff[j][i])
            {
              if ( tr2 <=  max_radius ) badrad++;

              if ( (i >= xb) && (i <= xe) && (j >= yb) && (j <= ye) ) badsq++;
            }
          else
            {
             if ( tr2 <= max_radius ) goodrad++;
             if ( (i >= xb) && (i <= xe) && (j >= yb) && (j <= ye) ) goodsq++;
             ngood++;
           }
        }
    }
      
//c$$$      write (6,*) 'In Circle ', goodrad, ' good, ', badrad, ' bad'
//c$$$      write (6,*) 'In Square ', goodsq, ' good, ', badsq, ' bad'

  if (  (float)goodsq  / (float) (goodsq  + badsq) 
      > (float)goodrad / (float) (goodrad + badrad) )
    {
      shape = SQUARE;
    }
  else
    {
      shape = CIRCLE;
    }
  //+jwp 2006-06-01
  shape = SQUARE;
  //-jwp 2006-06-01

/*
C     Fill in any bad points in the mask.
C     extend array out by ?two? rows and columns.  This
C     will be the limit of the detector
C     ===============================================================
*/
  //+JWP 2008-12-08
  int nExtend = 3;  // Was hard-coded to 3 before
  //-JWP 2008-12-08
  ptMaskPoints->xb = ptMaskPoints->xb - nExtend;
  xb               = ptMaskPoints->xb;
  ptMaskPoints->yb = ptMaskPoints->yb - nExtend;
  yb               = ptMaskPoints->yb;
  ptMaskPoints->xe = ptMaskPoints->xe + nExtend;
  xe               = ptMaskPoints->xe;
  ptMaskPoints->ye = ptMaskPoints->ye + nExtend;
  ye               = ptMaskPoints->ye;

  radius = (int) (sqrtf(max_radius) + 0.99999 + 2.0);
      
  if (SQUARE == shape)
    {
      anclimits[0] = xb;
      anclimits[1] = yb;
      anclimits[2] = xe;
      anclimits[3] = ye;
      if (3 < m_nVerbose)
        cout << "Square " << xb << ' ' << xe << ' ' << yb << ' ' << ye << endl;
    }
  else if (CIRCLE == shape)
    {
      anclimits[0] = xcenter;
      anclimits[1] = ycenter;
      anclimits[2] = radius;
      anclimits[3] = 0;
      if (5 < m_nVerbose)
        cout << "Circle " << xcenter << ' ' << ycenter << ' ' << max_radius
             << ' ' << radius << endl;
    }

  cout << "***   Filling in mask points..." << endl;

/*
C     Any interpolated points will still be marked with 
C     QINTERPOLATED, but will have values in x_diff, y_diff
C     instead of RFLAG
C     ===============================================================
*/

  Cimage *poImgTempA;
  Cimage *poImgTempB;

  poImgTempA = new Cimage(MASK_SIZE, MASK_SIZE, eImage_realIEEE,
                         ptMaskPoints->x_diff);
  poImgTempB = new Cimage(*poImgTempA);

  if (5 < m_nVerbose)
    {
      float fPix;
      (void) poImgTempA->nSetNextPixel(0,0);
      (void) poImgTempB->nSetNextPixel(0,0);
      for (i = 0; i < poImgTempA->nGetDimension(); i++)
        {
          fPix = poImgTempA->fGetNextPixel();
          if (RFLAG == fPix)
            {
              fPix = 0.0;
            }
          poImgTempB->vSetNextPixel((float)(fPix * 10.0 + 1000.0));
        }
      poImgTempB->nWrite("x_diffin.img");
    }
  
  vCompleteArray(poImgTempA, poImgTempB, 
                 anclimits[0], anclimits[1], anclimits[2], anclimits[3],
                 RFLAG);

  // Copy filled in array back into ptMaskPoints

  (void) poImgTempB->nSetNextPixel(0,0);
  for (j = 0; j < MASK_SIZE; j++)
    {
      for (i = 0; i < MASK_SIZE; i++)
        {
          ptMaskPoints->x_diff[j][i] = poImgTempB->fGetNextPixel();
        }
    }

  if (5 < m_nVerbose)
    {
      float fPix;
      (void) poImgTempB->nSetNextPixel(0,0);
      for (i = 0; i < poImgTempB->nGetDimension(); i++)
        {
          fPix = poImgTempB->fGetNextPixelNoInc();
          if (RFLAG == fPix)
            {
              fPix = 0.0;
            }
          poImgTempB->vSetNextPixel((float)(fPix * 10.0 + 1000.0));
        }
      poImgTempB->nWrite("x_diffout.img");
    }
  
  delete poImgTempA;

  poImgTempA = new Cimage(MASK_SIZE, MASK_SIZE, eImage_realIEEE,
                          ptMaskPoints->y_diff);

  if (5 < m_nVerbose)
    {
      float fPix;
      (void) poImgTempA->nSetNextPixel(0,0);
      (void) poImgTempB->nSetNextPixel(0,0);
      for (i = 0; i < poImgTempA->nGetDimension(); i++)
        {
          fPix = poImgTempA->fGetNextPixel();
          if (RFLAG == fPix)
            {
              fPix = 0.0;
            }
          poImgTempB->vSetNextPixel((float)(fPix * 10.0 + 1000.0));
        }
      poImgTempB->nWrite("y_diffin.img");
    }
  
  vCompleteArray(poImgTempA, poImgTempB, 
                 anclimits[0], anclimits[1], anclimits[2], anclimits[3],
                 RFLAG);

  (void) poImgTempB->nSetNextPixel(0,0);
  for (j = 0; j < MASK_SIZE; j++)
    {
      for (i = 0; i < MASK_SIZE; i++)
        {
          ptMaskPoints->y_diff[j][i] = poImgTempB->fGetNextPixel();
        }
    }

  if (5 < m_nVerbose)
    {
      float fPix;
      (void) poImgTempB->nSetNextPixel(0,0);
      for (i = 0; i < poImgTempB->nGetDimension(); i++)
        {
          fPix = poImgTempB->fGetNextPixelNoInc();
          if (RFLAG == fPix)
            {
              fPix = 0.0;
            }
          poImgTempB->vSetNextPixel((float)(fPix * 10.0 + 1000.0));
        }
      poImgTempB->nWrite("y_diffout.img");
    }

  delete poImgTempA;
  delete poImgTempB;

/*
C     fill in 'unpinned' positions and
C     'observed' positions for all interpolated positions
C     =======================================================
*/
  minxo = ptMaskPoints->xo[CTR][CTR];
  maxxo = ptMaskPoints->xo[CTR][CTR];
  minyo = ptMaskPoints->yo[CTR][CTR];
  maxyo = ptMaskPoints->yo[CTR][CTR];
  minxc = ptMaskPoints->xc[CTR][CTR];
  maxxc = ptMaskPoints->xc[CTR][CTR];
  minyc = ptMaskPoints->yc[CTR][CTR];
  maxyc = ptMaskPoints->yc[CTR][CTR];
  
  for (j = yb; j <= ye; j++)
    {
      for (i = xb; i <= xe; i++)
        {
          if (   (ptMaskPoints->status[j][i] == QINTERPOLATED) 
              && (ptMaskPoints->x_diff[j][i] != RFLAG) )
            {
              ptMaskPoints->xu[j][i] = ptMaskPoints->xc[j][i] + ptMaskPoints->x_diff[j][i];
              ptMaskPoints->yu[j][i] = ptMaskPoints->yc[j][i] + ptMaskPoints->y_diff[j][i];
              
              vReverseUnpinPoint ( ptDistorParams,
                                  ptMaskPoints->xu[j][i], ptMaskPoints->yu[j][i], 
                                  &ptMaskPoints->xo[j][i], &ptMaskPoints->yo[j][i] );
            }
            
          if ( ptMaskPoints->xo[j][i] < minxo ) minxo = ptMaskPoints->xo[j][i];
          if ( ptMaskPoints->xo[j][i] > maxxo ) maxxo = ptMaskPoints->xo[j][i];
          if ( ptMaskPoints->yo[j][i] < minyo ) minyo = ptMaskPoints->yo[j][i];
          if ( ptMaskPoints->yo[j][i] > maxyo ) maxyo = ptMaskPoints->yo[j][i];
          if ( ptMaskPoints->xc[j][i] < minxc ) minxc = ptMaskPoints->xc[j][i];
          if ( ptMaskPoints->xc[j][i] > maxxc ) maxxc = ptMaskPoints->xc[j][i];
          if ( ptMaskPoints->yc[j][i] < minyc ) minyc = ptMaskPoints->yc[j][i];
          if ( ptMaskPoints->yc[j][i] > maxyc ) maxyc = ptMaskPoints->yc[j][i];
        }
    }
      
/*
C     If the workstation is on then mark the interpolated mask points
C    ===============================================================
*/

  if (NULL != m_poReflnlistMask)
    {
      Crefln *poRefln;
      Creflnlist *poPeaks;
      poPeaks = m_poReflnlistMask;
      poRefln = new Crefln(m_poReflnlistMask);

      for (j = yb; j <= ye; j++)
        {
          for (i = xb; i <= xe; i++)
            {
//              if ( ptMaskPoints->x_diff[j][i] != RFLAG )
              if ( ptMaskPoints->status[j][i] == QINTERPOLATED )
                {
                  poRefln->vSetH(i);
                  poRefln->vSetK(j);
                  poRefln->vSetL(ptMaskPoints->status[j][i]);
                  poRefln->vSetField(poPeaks->m_nFI_nNonunfFlag, 1);
                  poRefln->vSetIntensity(ptMaskPoints->inten[j][i]);
                  poRefln->vSetField(poPeaks->m_nFI_fObsPx0, ptMaskPoints->xo[j][i]);
                  poRefln->vSetField(poPeaks->m_nFI_fObsPx1, ptMaskPoints->yo[j][i]);
                  poPeaks->nInsert(poRefln);
                }
            }
        }
      poPeaks->nWrite("dtcalmask2.ref");
//      delete poRefln;
    }

/*
C     This section test the masks, and makes certain
C     that the conversion routines all work.  Should
C     be able to go from observed to corrected and back again
*/
  if (3 < m_nVerbose)
    {
      cout << "\nTesting Mask Tables/Unpin Routines\n"
           << "\nMask Pos   " << CTR-10 << ", " << CTR-10
           << "\nObserved   " << ptMaskPoints->xo[CTR-10][CTR-10] << ", "
                              << ptMaskPoints->yo[CTR-10][CTR-10]
           << "\nCalculated " << ptMaskPoints->xc[CTR-10][CTR-10] << ", "
                              << ptMaskPoints->yc[CTR-10][CTR-10]
           << "\nUnpinned   " << ptMaskPoints->xu[CTR-10][CTR-10] << ", "
                              << ptMaskPoints->yu[CTR-10][CTR-10]
           << "\nC-U        " << ptMaskPoints->x_diff[CTR-10][CTR-10] << ", "
                              << ptMaskPoints->y_diff[CTR-10][CTR-10]
           << endl;

      x = ptMaskPoints->xc[CTR-10][CTR-10];
      y = ptMaskPoints->yc[CTR-10][CTR-10];
      vCalcMaskPosition(ptCorrectParams, x, y, &xm, &ym);
      cout << "Calc Mask P " << xm << ", " << ym << endl;
      jjj = 0;
      for (jj = (int)ym-1; jj <= (int)ym+2; jj++, jjj++)
        {
          iii = 0;
          for (ii = (int)xm-1; ii <= (int)xm+2; ii++, iii++)
            {
              if ( ptMaskPoints->x_diff[ii][jj] == RFLAG )
                {
                  cout << "ERROR!\n" << endl;
                }
              else
                {
                  afstempx[jjj][iii] = ptMaskPoints->x_diff[jj][ii];
                  afstempy[jjj][iii] = ptMaskPoints->y_diff[jj][ii];
                }
            }
        }

      cout << "\n******* X delta matrix ********";
      vListMatMN(4, 4, (float *)afstempx);
      cout << "\n******* Y delta matrix ********";      
      vListMatMN(4, 4, (float *)afstempy);
      xu = xm - int(xm);
      yu = ym - int(ym);
      xm = fBSpline2(afstempx, xu, yu);
      ym = fBSpline2(afstempy, xu, yu);

      cout << "Calc C-U at " << xu << ", " << yu
           << "\nCalc C-U    " << xm << ", " << ym;
      xu = x + xm;
      yu = y + ym;
      cout << "\nCalc U      " << xu << ", " << yu;

      vReverseUnpinPoint(ptDistorParams, xu, yu, &xi, &yi );
      cout << "\nCalc O      " << xi << ", " << yi << endl;
    }

/*
C============================================================
C     For every pixel in a corrected image
C     1.  Calculate the position in MASK coordinates
C     2.  Interpolate from nearest MASK positions to find
C     the difference
C     3.  Apply the difference to find the UNPINNED position
C     4.  Reverse unpin to find the OBSERVED position
C
C
C     In the good old days of image intensifiers, I could count on 
C     the corrected image shrinking relative to the distorted
C     image, and that neither the distorted or the corrected images
C     would approach the edge of the detector.  With the newer
C     'Modular' detectors, the image goes all the way to the
C     edge of the detector.  Now, depending on the choice of the
C     corrected pixel size, it is possible that the expanded
C     image goes off the image.  This is no problem mathematically
C     (negative numbers are as good as positive numbers),
C     in MADNES (a corrected image is never formed), or when correcting 
C     multiple detectors (the image goes into a larger image).  This 
C     could present a problem when correcting a single image if the
C     array space for the corrected image is dimensioned to be the same 
C     as the uncorrected image.
C
C=============================================================
*/
  if (1 == nEdges)
    {
      xint_step  = nint ( (maxxo - minxo) / axsize + 0.5f );
      xint_start = nint (minxo);
      yint_step  = nint ( (maxyo - minyo) / aysize + 0.5f );
      yint_start = nint (minyo);
      xinv_step  = nint ( (maxxc - minxc) / axsize + 0.5f );
      xinv_start = nint (minxc);
      yinv_step  = nint ( (maxyc - minyc) / aysize + 0.5f );
      yinv_start = nint (minyc);
    }
  else
    {
      xint_step  = nint ( float(nDim0) / axsize + 0.4999f );
      xint_start = 1;
//      xint_start = -7;
      yint_step  = nint ( float(nDim1) / aysize + 0.4999f );
      yint_start = 1;
//      yint_start = -7;
      xinv_step  = nint ( float(nDim0) / axsize + 0.4999f );
      xinv_start = 1;
      yinv_step  = nint ( float(nDim1) / aysize + 0.4999f );
      yinv_start = 1;
    }
  pscale = ptCorrectParams->pscale;  // 1.0 / 512.0; // 0.01; // TODO What should this be?  1 / 512?

  if (3 < m_nVerbose)
    {
      cout << "Pscale: " << pscale << endl;
      cout << "\nEdges: " << nEdges
           << "\nminxo, maxxo, minyo, maxyo: " 
           << minxo << ", " << maxxo << ", " << minyo << ", " <<  maxyo
           << "\nminxc, maxxc, minyc, maxyc: "
           << minxc << ", " <<  maxxc << ", " <<  minyc << ", " <<  maxyc
           << "\nxint_step, xint_start, yint_step, yint_start: "
           << xint_step << ", " <<  xint_start << ", " <<  yint_step << ", " <<  yint_start
           << "\nxinv_step, xinv_start, yinv_step, yinv_start: "
           << xinv_step << ", " <<  xinv_start << ", " <<  yinv_step << ", " <<  yinv_start
           << endl;
    }

  ptCorrectParams->xint_step  = xint_step;
  ptCorrectParams->yint_step  = yint_step;
  ptCorrectParams->xinv_step  = xinv_step;
  ptCorrectParams->yinv_step  = yinv_step;
  ptCorrectParams->xint_start = xint_start;
  ptCorrectParams->yint_start = yint_start;
  ptCorrectParams->xinv_start = xinv_start;
  ptCorrectParams->yinv_start = yinv_start;

  // The next two lines really come from the input and should NOT be
  // set here as they were in Marty's original program
  //  ptCorrectParams->pscale     = pscale;
  //  ptCorrectParams->badflag    = INTFLAG;
  
  if (3 < m_nVerbose)
    {
      cout << "\nMore correction parameters"
           << "\nxint,yint_step:     " << ptCorrectParams->xint_step << ", "
                                       << ptCorrectParams->yint_step
           << "\nxinv,yinv_step:     " << ptCorrectParams->xinv_step << ", "
                                       << ptCorrectParams->yinv_step
           << "\nxint,yint_start:    " << ptCorrectParams->xint_start << ", "
                                       << ptCorrectParams->yint_start
           << "\npscale, 1/pscale:   " << ptCorrectParams->pscale << ", "
                                       << 1.0 / ptCorrectParams->pscale
           << "\nbadflag:            " << ptCorrectParams->badflag
           << endl;
    }

  cout << "***   Calculating REVERSE interpolation tables..." << endl;

  // Loop over each expected hole (peak) in the mask
  // Calculate its mm position

  timx = axsize;
  timy = aysize;
  for (j = 1; j <= timy; j++)
    {
      y = yinv_step * (j-1)  +  yinv_start;
      for (i = 1; i <= timx; i++)
        {
	  //+JWP 2006-06-01  Set pixel value to be BAD to start with!

	  poImgDistorXINV->nSetPixel(i-1, j-1, lIntBadFlag);
	  poImgDistorYINV->nSetPixel(i-1, j-1, lIntBadFlag);

	  //-JWP 2006-06-01  Set pixel value to be BAD to start with!

          x = xinv_step * (i-1)  + xinv_start;

          //   x and y are pixel positions in the correct image.
          //   determine the mask points (also in the corrected image) 
          //   which would surround it

          vCalcMaskPosition( ptCorrectParams, x, y, &xm, &ym );
//          cout << "A vCMPos: x,y,xm,ym: "
//               << x << ", " << y << ", " << xm << ", " << ym << endl;
          nStat = 0;
	  //+JWP 2008-12-17
          //if(   (xm <  1) || (xm >= MASK_SIZE-1)
	  //   || (ym <  1) || (ym >= MASK_SIZE-1))
          if(   (xm <  0) || (xm >= MASK_SIZE-1)
             || (ym <  0) || (ym >= MASK_SIZE-1))
          //-JWP
            {
              nStat = 1;
            }
          else
            {
              // Are there 16 points surrounding it?

              jjj = 0;
              for (jj = (int)ym-1; (jj <= (int)ym+2) && (0 == nStat); jj++)
                {
                  iii = 0;
                  for (ii = (int)xm-1; (ii <= (int)xm+2) && (0 == nStat); ii++)
                    {
                      if ( ptMaskPoints->x_diff[jj][ii] == RFLAG ) 
                        {
                          nStat = 2;
//C                       cout << "Pixel " << i << ' ' << j << endl;
                          break;
                        }
                      afstempx[jjj][iii] = ptMaskPoints->x_diff[jj][ii];
                      afstempy[jjj][iii] = ptMaskPoints->y_diff[jj][ii];
                      iii = iii + 1;
                    }
                  jjj = jjj + 1;
                }
            }
          if (0 == nStat)
            {
              /*            
                 C     Yes, there are 16 points. Spline fit
                 C     
                 C     This gives the difference that would be applied
                 C     to a point after unpinning to get to the corrected
                 C     position
                 C     =================================================
                 */

              xm = xm - (float)((int)xm);
              ym = ym - (float)((int)ym);

              xu = x + fBSpline2(afstempx, xm, ym);
              yu = y + fBSpline2(afstempy, xm, ym);

              vReverseUnpinPoint( ptDistorParams, xu, yu, &xi, &yi);
              if ( (xi == 0.0) && (yi == 0.0 ) )
                {
                  nStat = 3;
                }
              else
                {
                  poImgDistorXINV->nSetPixel(i-1, j-1, (ULONG)nint(xi / pscale));
                  poImgDistorYINV->nSetPixel(i-1, j-1, (ULONG)nint(yi / pscale));
                                             
//                  x_inv[j][i] = nint(xi / pscale);  // !TODO
//                  y_inv[j][i] = nint(yi / pscale);
                }
            }
          if (0 != nStat)
            {
              // Here if there was an error
              // ----                 -----

//              x_inv[j][i] = INTFLAG;
//              y_inv[j][i] = INTFLAG;

              poImgDistorXINV->nSetPixel(i-1, j-1, lIntBadFlag);
              poImgDistorYINV->nSetPixel(i-1, j-1, lIntBadFlag);
            }
        } // end i loop
    } // end j loop

  poImgDistorXINV->vSetSatValue((float)lIntBadFlag);
  poImgDistorYINV->vSetSatValue((float)lIntBadFlag);
  poImgDistorXINV->vSetState(eImage_available_state);
  poImgDistorYINV->vSetState(eImage_available_state);

  if (3 < m_nVerbose)
    {
      // This section tests the REVERSE interpolation tables
      // and demonstrates how to use them

      // It is not really necessary to use b-spline interpolation
      // but it was easy to do here, so why not.  Normally it
      //would be more efficient to use linear interpolation

      cout << "\nTesting REVERSE Interpolation Tables\n";
      x = ptMaskPoints->xc[CTR-10][CTR-10];
      y = ptMaskPoints->yc[CTR-10][CTR-10];
      cout << "\nMask Pos     " << CTR-10 << ", " << CTR-10
           << "\nCorrected    " << x << ", " << y;
      x = (x-xinv_start) / xinv_step + 1;
      y = (y-yinv_start) / yinv_step + 1;
      cout << "\nTable Index  " << x << ", " << y;
      jjj = 0;
      for (jj = (int)y-2; jj <= (int)y+1; jj++, jjj++)
        {
          iii = 0;
          for (ii = (int)x-2; ii <= (int)x+1; ii++, iii++)
            {
//              cout << "\nii, jj: " << ii << ", " << jj;
              afstempx [jjj][iii] = pscale 
                * (poImgDistorXINV->*poImgDistorXINV->prfGetPixel)(ii, jj);
              afstempy [jjj][iii] = pscale
                * (poImgDistorYINV->*poImgDistorYINV->prfGetPixel)(ii, jj);
            }
        }
      cout << "\n******* X REVERSE Interpolation matrix ********";
      vListMatMN(4, 4, (float *)afstempx);
      cout << "\n******* Y REVERSE Interpolation matrix ********";
      vListMatMN(4, 4, (float *)afstempy);
      x = x - (float)((int)x);
      y = y - (float)((int)y);
      cout << "\nCalc C-U at  " << x << ", " << y
           << "\nCalc Obs     " << fBSpline2(afstempx,x,y) << ", "
           << fBSpline2(afstempy,x,y) << endl;
    }

/*
C     Now we have the REVERSE interpolation tables.
C     From these (and only these) calculate the FORWARD 
C     interpolation tables.
C     1.  At each position in the REVERSE interpolation table
C     determine which pixels fall within box determined by
C     nearest pixels
C     2.  Check to make certain it is in box
C     3.  Solve by linear interpolation
C===============================================================
*/

  cout << "***   Calculating FORWARD interpolation tables..." << endl;

  if (3 < m_nVerbose)
    cout << "timx, timy: " << timx << ", " << timy << endl;

  LONG lPix[8];

  for (j = 1; j <= timy-1; j++)
    {
      //cout << "JJJJJJJJJJJJJJ is: " << j << endl;
      for (i = 1; i <= timx-1; i++)
        {
	  //cout << "IIIIIIIIIIIIIIIII is: " << i << endl;

          // All corners must have values in them

          poImgDistorXINV->nGetPixel(i-1, j-1, &lPix[0]);
          poImgDistorXINV->nGetPixel(i,   j-1, &lPix[1]);
          poImgDistorXINV->nGetPixel(i-1, j,   &lPix[2]);
          poImgDistorXINV->nGetPixel(i,   j,   &lPix[3]);
          poImgDistorYINV->nGetPixel(i-1, j-1, &lPix[4]);
          poImgDistorYINV->nGetPixel(i,   j-1, &lPix[5]);
          poImgDistorYINV->nGetPixel(i-1, j,   &lPix[6]);
          poImgDistorYINV->nGetPixel(i,   j,   &lPix[7]);

	  /****
jwp	  // Set them bad to start with

	  indexi = (ii - xint_start) / xint_step + 1;
	  indexj = (jj - yint_start) / yint_step + 1;
	  poImgDistorXINT->nSetPixel(indexi-1, indexj-1,
				     lIntBadFlag);
	  poImgDistorYINT->nSetPixel(indexi-1, indexj-1,
				     lIntBadFlag);
 ****/
          if (   (lIntBadFlag != lPix[0])
              && (lIntBadFlag != lPix[1])
              && (lIntBadFlag != lPix[2])
              && (lIntBadFlag != lPix[3])
              && (lIntBadFlag != lPix[4])
              && (lIntBadFlag != lPix[5])
              && (lIntBadFlag != lPix[6])
              && (lIntBadFlag != lPix[7]) )
            {
              // Determine min, max positions that might be
              // within this box

              x1 = (float)lPix[0] * pscale;
              x2 = (float)lPix[1] * pscale;
              x3 = (float)lPix[2] * pscale;
              x4 = (float)lPix[3] * pscale;
              y1 = (float)lPix[4] * pscale;
              y2 = (float)lPix[5] * pscale;
              y3 = (float)lPix[6] * pscale;
              y4 = (float)lPix[7] * pscale;
            
              nxmin = (int) (min(x1, min(x2, min(x3, x4 ))));
              nxmax = (int) (max(x1, max(x2, max(x3, x4 )))) + 1;
              nymin = (int) (min(y1, min(y2, min(y3, y4 ))));
              nymax = (int) (max(y1, max(y2, max(y3, y4 )))) + 1;

	      //cout << "nymin, nymax: " << nymin << ", " << nymax << endl;
              for (jj = nymin; jj <= nymax; jj++)
                {
		  //cout << "nxmin, nxmax: " << nxmin << ", " << nxmax << endl;
                  for (ii = nxmin; ii <= nxmax; ii++)
                    {
                      // Is it a step pixel ?

                      if (   (0 == ((ii - xint_start) % xint_step))
                          && (0 == ((jj - yint_start) % yint_step)) )
                        {
                          // Is it really within this box ?

                          if (0 == nInBox( float(ii), float(jj), 
                                           x1, y1, x2, y2, 
                                           x3, y3, x4, y4 ))
                            {
                              // Do a simple linear interpolation

                              a = ((float)ii-x3) / (x4-x3) 
                                   - ((float)ii-x1)/(x2-x1);
                              c = ((float)ii-x1) / (x2-x1);
                              m = ((float)jj-y2) / (y4-y2) 
                                    - ((float)jj-y1)/(y3-y1);
                              b = ((float)jj-y1) / (y3-y1);
                              xp = (a * b  + c) / (1.0 - a * m);
                              yp = m * xp + b;
                           
                              indexi = (ii - xint_start) / xint_step + 1;
                              indexj = (jj - yint_start) / yint_step + 1;
/*
                              if (   (indexi <= axsize)
                                  && (indexj <= aysize)
                                  && (indexi >= 1)
                                  && (indexj >= 1 ) )
                                {

                              x_int  (indexi, indexj) = 
     $                             nint((xinv_step*(i+xp-1)+xinv_start)
     $                             / pscale)
                              y_int  (indexi, indexj) = 
     $                             nint((yinv_step*(j+yp-1)+yinv_start)
     $                             / pscale)
*/
                              poImgDistorXINT->nSetPixel(indexi-1, indexj-1,
                                (LONG)nint((xinv_step*((float)i+xp-1.0)
                                              +xinv_start) / pscale));
                              poImgDistorYINT->nSetPixel(indexi-1, indexj-1,
                                (LONG)nint((yinv_step*((float)j+yp-1.0)
                                              +yinv_start) / pscale));
//                                }
                            }
                        }
                    }
                }
            }
        }
    }
  cout << "ALL DONE with FORWARD.\n" << endl;
  poImgDistorXINT->vSetSatValue((float)lIntBadFlag);
  poImgDistorYINT->vSetSatValue((float)lIntBadFlag);
  poImgDistorXINT->vSetState(eImage_available_state);
  poImgDistorYINT->vSetState(eImage_available_state);

/*
c$$$C     This section tests the forward interpolation tables
c$$$C     and demonstrates how to use them
c$$$C
c$$$C     It is not really necessary to use b-spline interpolation
c$$$C     but it was easy to do here, so why not.  Normally it
c$$$C     would be more efficient to use linear interpolation
// JWP note:  Could you nPxToPx() routine for the linear interpolation
*/

  if (3 < m_nVerbose)
    {
      x = ptMaskPoints->xo[CTR-10][CTR-10];
      y = ptMaskPoints->yo[CTR-10][CTR-10];
      cout << "Testing FORWARD Interpolation Tables\n"
           << "\nMask Pos     " << CTR-10 << ", " << CTR-10
           << "\nObserved     " << x << ", " << y << endl;
      x = (x-xint_start) / xint_step + 1;
      y = (y-yint_start) / yint_step + 1;
      cout << "Table Index  " << x << ", " <<  y << endl;
      jjj = 0;
      for (jj = (int)y-2; jj <= (int)y+1; jj++, jjj++)
        {
          iii = 0;
          for (ii = (int)x-2; ii <=(int)x+1; ii++, iii++)
            {
              afstempx [jjj][iii] = pscale 
                * (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ii, jj);
              afstempy [jjj][iii] = pscale
                * (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ii, jj);
            }
        }

      cout <<  "******* X FORWARD Interpolation matrix ********";
      vListMatMN(4, 4, (float*)afstempx);
      cout <<  "******* Y FORWARD Interpolation matrix ********";
      vListMatMN(4, 4, (float*)afstempy);
      x = x - (float)((int)x);
      y = y - (float)((int)y);
      cout << "\nCalc at      " << x << ", " << y
           << "\nCalc Correct " << fBSpline2(afstempx,x,y) << ", "
           << fBSpline2(afstempy,x,y) << endl;
    }
/*
C     ===============================================
C     Smooth out the funny edges
C     The edges might still not look that smooth, but
C     it does make it alot better!
C     ===============================================
*/
   cout << "***   Smoothing edges..." << endl;

   iii = timx;
   jjj = timy;
/*
C     Determine fill limits to smooth interpolation tables.
C     This could be included in the loops above, but who is
C     in that big a hurry
C
*/
  if ( shape == CIRCLE )
    {
/*
C        The position of the center in distorted space
C
C         x = (ptMaskPoints->xo(CTR,CTR)-xint_start) / xint_step + 1
C         y = (ptMaskPoints->yo(CTR,CTR)-yint_start) / yint_step + 1
C
C        The position of the center in undistorted space
C
C         xi = (ptMaskPoints->xc(CTR,CTR)-xinv_start) / xinv_step + 1
C         yi = (ptMaskPoints->yc(CTR,CTR)-yinv_start) / yinv_step + 1
C         write (6,*) 'Distortion Center ', x, y, xi, yi
C
C
C     Find the centroid of the interpolation tables
*/
      max_radius  = 0.0;
      max_radiusi = 0.0;
      x  = 0.0;
      y  = 0.0;
      xi = 0.0;
      yi = 0.0;
      for (j = 1; j <= timy; j++)
        {
          for (i = 1; i <= timx; i++)
            {
              poImgDistorXINT->nGetPixel(i-1, j-1,     &lPix[0]);
              if (lIntBadFlag != lPix[0])
                {
                  max_radius = max_radius + 1.0;
                  x = x + (float)i;
                  y = y + (float)j;
                }
              poImgDistorXINV->nGetPixel(i-1, j-1,     &lPix[0]);
              if (lIntBadFlag != lPix[0])
                {
                  max_radiusi = max_radiusi + 1.0;
                  xi = xi + (float)i;
                  yi = yi + (float)j;
                }
            }
        }

      x  = x / max_radius;
      y  = y / max_radius;
      xi = xi / max_radiusi;
      yi = yi / max_radiusi;
/*
C         write (6,*) 'Centeroid         ', x, y, xi, yi
C
C
*/
      max_radius  = 0;
      max_radiusi = 0;
      for (j = 1; j <= timy; j++)
        {
          ty2  = (float)j - y;
          ty2  = ty2 * ty2;
          tyi2 = (float)j - yi;
          tyi2 = tyi2 * tyi2;
          for (i = 1; i <= timx; i++)
            {
              poImgDistorXINT->nGetPixel(i-1, j-1,     &lPix[0]);
              if (lIntBadFlag != lPix[0])
                {
                  tr2 = (float)i - x;
                  tr2 = tr2 * tr2  +  ty2;
                  if ( tr2 > max_radius ) max_radius = tr2;
                }
              poImgDistorXINV->nGetPixel(i-1, j-1,     &lPix[0]);
              if (lIntBadFlag != lPix[0])
                {
                  tri2 = (float)i - xi;
                  tri2 = tri2 * tri2  +  tyi2;
                  if ( tri2 > max_radiusi ) max_radiusi = tri2;
                }
              }
          }
      max_radius  = sqrtf(max_radius);
      max_radiusi = sqrtf(max_radiusi);

      xb = nint(x);
      yb = nint(y);
      xe = nint(max_radius);
      ye = 0;

//C         write (6,*) 'Filling int tables'
//C         write (6,*) 'Circle ', xb, yb, xe
         

      xbi = nint(xi);
      ybi = nint(yi);
      xei = nint(max_radiusi);
      yei = 0;

      cout << "Filling inverse int tables\n"
           << "Circle " << xbi << ", " <<  ybi << ", " << xei << endl;
    }
  else if ( SQUARE == shape )
    {
      xb  = timx;
      yb  = timy;
      xe  = 1;
      ye  = 1;
      xbi = timx;
      ybi = timy;
      xei = 1;
      yei = 1;
      for (j = 1; j <= timy; j++)
        {
          for (i = 1; i <= timx; i++)
            {
              poImgDistorXINT->nGetPixel(i-1, j-1,     &lPix[0]);
              if (lIntBadFlag != lPix[0])
                {
                  if ( i < xb ) xb = i;
                  if ( i > xe ) xe = i;
                  if ( j < yb ) yb = j;
                  if ( j > ye ) ye = j;
                }
              poImgDistorXINV->nGetPixel(i-1, j-1,     &lPix[0]);
              if (lIntBadFlag != lPix[0])
                {
                  if ( i < xbi ) xbi = i;
                  if ( i > xei ) xei = i;
                  if ( j < ybi ) ybi = j;
                  if ( j > yei ) yei = j;
                }
            }
          }
         xb = max(1, xb-1);
         yb = max(1, yb-1);
         xe = min(timx, xe+1);
         ye = min(timy, ye+1);
//C         write (6,*) 'Filling int tables'
//C         write (6,*) 'Square ', xb, xe, yb, ye

         xbi = max(1, xbi-1);
         ybi = max(1, ybi-1);
         xei = min(timx, xei+1);
         yei = min(timy, yei+1);
//C         write (6,*) 'Filling inverse int tables'
//C         write (6,*) 'Square ', xbi, xei, ybi, yei
      }

  
  poImgTempB = new Cimage(axsize, aysize, eImage_I4);

  xb = yb = 0;
  xe = axsize - 1;
  ye = aysize - 1;
  vCompleteArray(poImgDistorXINT, poImgTempB, 
                 xb, yb, xe, ye, (float)lIntBadFlag);

  (void) poImgTempB->nSetNextPixel(0,0);
  (void) poImgDistorXINT->nSetNextPixel(0,0);
  for (i = 0; i < nDim; i++)
    {
      poImgDistorXINT->vSetNextPixel(poImgTempB->lGetNextPixel());
    }

  vCompleteArray(poImgDistorYINT, poImgTempB,
                 xb, yb, xe, ye, (float)lIntBadFlag);
  (void) poImgTempB->nSetNextPixel(0,0);
  (void) poImgDistorYINT->nSetNextPixel(0,0);
  for (i = 0; i < nDim; i++)
    {
      poImgDistorYINT->vSetNextPixel(poImgTempB->lGetNextPixel());
    }

  vCompleteArray(poImgDistorXINV, poImgTempB,
                 xbi, ybi, xei, yei, (float)lIntBadFlag);
  (void) poImgTempB->nSetNextPixel(0,0);
  (void) poImgDistorXINV->nSetNextPixel(0,0);
  for (i = 0; i < nDim; i++)
    {
      poImgDistorXINV->vSetNextPixel(poImgTempB->lGetNextPixel());
    }

  vCompleteArray(poImgDistorYINV, poImgTempB,
                 xbi, ybi, xei, yei, (float)lIntBadFlag);
  (void) poImgTempB->nSetNextPixel(0,0);
  (void) poImgDistorYINV->nSetNextPixel(0,0);
  for (i = 0; i < nDim; i++)
    {
      poImgDistorYINV->vSetNextPixel(poImgTempB->lGetNextPixel());
    }
  delete poImgTempB;

  return (0);
}

void
Ccalibrate::vCalcExpected(tagPoints             *ptMaskPoints,
                          tagDistort_parameters *ptDistorParams,
                          tagMask_parameters    *ptMaskParams,
                          tagCorrect_parameters *ptCorrectParams)
{
/*
C------------------------------------------------------------------------

      subroutine calc_expected ( mask, dparams, mparams, cparams )

C------------------------------------------------------------------------
C       for each point in the mask, calculate the expected position
C
C     The expected position is based on the spacing, mask angle
C     and peak index.
C
C     Mask Angle:
C     ===========
C     If this routine is called in workstation version, mask
C     angle will be set to the input mask angle
C
C     Else will prompt for correct values if too far off expected
C
C     Mask Spacing:
C     =============
C     If workstation mode, will use x/y ratio from maskparams
C     Else will prompt if different from 1
C     
C     The y spacing will be changed adjusted to reflect any change
C     in x/y ratio. This will have the effect of requiring larger
C     changes in the y position when the interpolation tables are
C     based on these expected values.
C
C     Could modify this routine so that it returns the correct
C     position based on the desired pixel size and the known 
C     mask hole spacing
C
C------------------------------------------------------------------------
*/
  float       xcenter, ycenter,
              xpt, ypt,
              xscale, yscale,
              sinv, cosv,
              sinh, cosh,
              xy, yx, xx, yy,
              xe, ye;
  
  int         i, j;
//  float       readdef;
  float       angle, hor_angle, ver_angle;
  float       yxratio;

  xscale    = ptDistorParams->x_scale;
  yscale    = ptDistorParams->y_scale;
  yxratio   = yscale / xscale;

  hor_angle = atanf(ptDistorParams->horz_slope) / Gs_fRADIANS_PER_DEGREE;
  ver_angle = atanf(-ptDistorParams->ver_slope) / Gs_fRADIANS_PER_DEGREE;

  if (NULL != m_poOutDump)
    {
      *m_poOutDump << "! Observed mask horizontal angle : " << hor_angle 
                   << "\n! Observed mask vertical angle :   " << ver_angle
                   << "\n! X scale, Y scale :               " << xscale << ", " << yscale
                   << "\n! Observed mask Y/X ratio :        " <<  yxratio
                   << endl;
    }

  cout << "! Observed mask Y/X ratio        : " <<  yxratio 
       << "\n! Observed mask horizontal angle : " << hor_angle
       << "\n! Observed mask vertical   angle : " << ver_angle 
       << "\n! X scale, Y scale :               " << xscale << ", " << yscale
       << endl;

  angle = ( hor_angle + ver_angle ) * 0.5; 

  // The following if assumes the mask has orthogonally placed holes,
  // that is rows and columns are 90 degrees apart.

  if (   (2.0 < fabsf( (90.0-hor_angle+ver_angle) -  ptMaskParams->maskparams[0]))
      || (5.0 < fabsf(angle-ptMaskParams->mask_angle))
      && (DO_UPDATE != ptMaskParams->work) )
    {
      cout << "Enter correct horizontal angle: ";
      cin  >> hor_angle;
      cout << "Enter correct vertical angle: ";
      cin  >> ver_angle;
//      hor_angle = readdef ('Enter correct horizontal angle',angle);
//      ver_angle = readdef ('Enter correct vertical angle',hor_angle);
    }
  else
    {
      hor_angle = ptMaskParams->mask_angle;
      ver_angle = ptMaskParams->mask_angle;
    }

  if (   (0.05 < fabsf(1.0-yxratio/ptMaskParams->maskparams[1]))
      && (DO_UPDATE != ptMaskParams->work) )
    {
//      yxratio = readdef ('Enter correct yxratio', 1.0 );
      cout << "Enter correct yxratio: ";
      cin  >> yxratio;
    }
  else
    {
      yxratio = ptMaskParams->maskparams[1];
    }

  yscale = yxratio * xscale;

  if (0.0 < ptMaskParams->pixel_size)
    {
      ptDistorParams->a1 = ptDistorParams->spacing 
                           / ptMaskParams->pixel_size / xscale;
      cout << "!Pixel size found in input\n"
           << "! Changing Old X scale, Y scale  : " << xscale << ", " << yscale
           << "\n! to\n";
      xscale = ptDistorParams->spacing / ptMaskParams->pixel_size;
      yscale = yxratio * xscale;
      cout << "!          New X scale, Y scale  : " << xscale << ", " << yscale
           << endl;
    }
  else
    {
      ptDistorParams->a1 = 1.0;
    }

  sinh    = sinf (hor_angle * Gs_fRADIANS_PER_DEGREE);
  cosh    = cosf (hor_angle * Gs_fRADIANS_PER_DEGREE);
  sinv    = sinf (ver_angle * Gs_fRADIANS_PER_DEGREE);
  cosv    = cosf (ver_angle * Gs_fRADIANS_PER_DEGREE);
  xcenter = ptDistorParams->x_center;
  ycenter = ptDistorParams->y_center;
  xpt = ptDistorParams->x_pt_center;
  ypt = ptDistorParams->y_pt_center;
      
//C      do j = ptMaskPoints->yb, ptMaskPoints->ye
//      do j = 1, MASK_SIZE
  for (j = 0; j < MASK_SIZE; j++)
    {
      ye = ((float)j - ypt) * yscale;
      xy = ye * (-sinv);
      yy = ye * cosv;
//C         do i = ptMaskPoints->xb, ptMaskPoints->xe
//         do i = 1, MASK_SIZE
      for (i = 0; i < MASK_SIZE; i++)
        {
          xe = ((float)i - xpt) * xscale;
          xx = xe * cosh;
          yx = xe * sinh;
          ptMaskPoints->xc[j][i] = xy + xx + xcenter;
          ptMaskPoints->yc[j][i] = yy + yx + ycenter;
        }
    }

  ptCorrectParams->x_center    = ptDistorParams->x_center;
  ptCorrectParams->y_center    = ptDistorParams->y_center;
  ptCorrectParams->x_pt_center = ptDistorParams->x_pt_center;
  ptCorrectParams->y_pt_center = ptDistorParams->y_pt_center;
  ptCorrectParams->x_scale     = xscale;
  ptCorrectParams->y_scale     = yscale;
  ptCorrectParams->horz_slope  = tanf(hor_angle  * Gs_fRADIANS_PER_DEGREE);
  ptCorrectParams->ver_slope   = -tanf(ver_angle * Gs_fRADIANS_PER_DEGREE);
  ptCorrectParams->ratio       = yxratio;

  if (3 < m_nVerbose)
    {
      cout << "\nCorrection parameters"
           << "\nx,y_center:         " << ptCorrectParams->x_center << ", "
                                       << ptCorrectParams->y_center 
           << "\nx,y_pt_center:      " << ptCorrectParams->x_pt_center << ", "
                                       << ptCorrectParams->y_pt_center 
           << "\nx,y_scale:          " << ptCorrectParams->x_scale << ", "
                                       << ptCorrectParams->y_scale
           << "\nhorz,ver_slope:     " << ptCorrectParams->horz_slope << ", "
                                       << ptCorrectParams->ver_slope
           << "\nyxratio:            " << ptCorrectParams->ratio
           << endl;
    }
}

void
Ccalibrate::vCalcUnpinned(tagPoints             *ptMaskPoints,
                          tagDistort_parameters *ptDistorParams,
                          tagMask_parameters    *ptMaskParams)
{
/*
C-----------------------------------------------------------------------
      subroutine calc_unpinned ( mask, dparams, mparams )
C-----------------------------------------------------------------------
C Calculate unpinned and scaled values for all mask points
C                        ======
C Given the observed position of each point, calculate
C      where unpinning, scaling in the y direction, and overall scaling
C      the point would bring it to
C
C------------------------------------------------------------------------
*/

  float     xcenter, ycenter;
  float     x, y;
  float     yxrat, rat, a, b, c;
  int       i, j;

  xcenter = ptDistorParams->x_center;
  ycenter = ptDistorParams->y_center;

  if (DO_RADIAL ==  ptMaskParams->radial)
    {
      a = ptDistorParams->a * ptDistorParams->a1;
      b = ptDistorParams->b;
      c = ptDistorParams->c;
    }
  else
    {
      a = ptDistorParams->a1;
      b = 0.0;
      c = 0.0;
    }

  yxrat = ptDistorParams->ratio;
  for (j = ptMaskPoints->yb; j <= ptMaskPoints->ye; j++)
    {
      for (i = ptMaskPoints->xb; i <= ptMaskPoints->xe; i++)
        {
          x = ptMaskPoints->xo[j][i] - xcenter;
          y = ptMaskPoints->yo[j][i] - ycenter;
          if (DO_RADIAL == ptMaskParams->radial)
            {
              rat = x * x  +  y * y;
              rat = a  +  b * sqrtf(rat)  +  c * rat;
            }
          else
            {
              rat = a;
            }
          ptMaskPoints->xu[j][i] = x * rat  + xcenter;
          ptMaskPoints->yu[j][i] = y * rat * yxrat  + ycenter;
        }
    }
}

void
Ccalibrate::vCompleteArray(Cimage *poImageIn, Cimage *poImageOut,
                           const int xstart, const int ystart, 
                           const int xend, const int yend, const float fFlag)
{
/*
C--------------------------------------------------------------
C     
C     ================================================
C     complete_array_i4 ( input, output, xsize, ysize, flag,
C     xstart, ystart, xend, yend)
C     ================================================
C     Given an incomplete array, fill in all points
C     between start and end by interpolation/extrapolation
C     
C     Returns :
C     none
C     
C     Arguments :
C     array - array to fill values with
C     type : real
C     access : modify
C     
C     xsize, ysize - dimensions of array
C     type : integer
C     access : readonly
C     
C     flag - value used to signal pixel to be replaced ( filled )
C     type : real
C     access : readonly
C
C     xstart, xend - x limits for fill
C     ystart, yend - y limits for fill
C     type : integer
C     access : readonly
C     
C     if yend = 0 then
C     xstart, xend = center of array
C     ystart = radius to fill
C     type : f floating
C     access : readonly
C     
C-----------------------------------------------------------------
      subroutine complete_array_i4 ( input, output,
     $     xsize, ysize, flag,
     $     xstart, ystart, xend, yend)
*/
      
  int nDim0, nDim1;
  int i, j, k;

  float mr2, rad2, y2;
  float xcenter, ycenter, max_radius;

  int xs, xe, ys, ye;

  const int MAXORDER = 2;
  float afC[MAXORDER+1];

  float *pfX, *pfY, *pfZ, *pfW, *pfXl, *pfYl, *pfWl;

  int   m;
  int   n, ii, jj, nl;
  int   x0, y0, xn, yn;
  int   nNumGood, nNumBad;
  float t1, t2;
  float fPixIJ;

  int MAXSIZE;
  //+JWP 2008-12-17
  // Change the integer 7 below to a variable
  int nEdgeSize = 9; // = 7
  //-JWP 2008-12-17
  // Allocate scratch arrays

  (void) poImageIn->nGetDimensions(&nDim0, &nDim1);

  MAXSIZE = max(nDim0, nDim1);
  pfX  = new float [MAXSIZE];
  pfY  = new float [MAXSIZE];
  pfW  = new float [MAXSIZE];
  pfXl = new float [MAXSIZE];
  pfYl = new float [MAXSIZE];
  pfWl = new float [MAXSIZE];

  // Find limits of usable points in the input array

  x0 = nDim0-1;
  xn = 0;
  y0 = nDim1-1;
  yn = 0;

  // Copy usable point input array to the output array while at the same time
  // finding the inner rectangle of all good points.

  for (j = 0; j < nDim1; j++)
    {
      for  (i = 0; i < nDim0; i++)
        {
          if ((poImageIn->*poImageIn->prfGetPixel)(i, j) != fFlag)
            {
              if ( i < x0 ) x0 = i;
              if ( i > xn ) xn = i;
              if ( j < y0 ) y0 = j;
              if ( j > yn ) yn = j;
            }
	  // The entire output array gets set to fFlag
          (poImageOut->*poImageOut->prnSetPixel)(i, j, fFlag);  // This is S L O O O W
        }
    }

/*
C     Find all pixels in the search area which are fFlagged
C
C     yend > 0 => fill rectangular region
*/

  if ( yend > 0 )
    {
      // It is a rectangular region

      xs = xstart;
      if ( xs < 0 ) xs = 0;
      xe = xend;
      if ( xe >= nDim0 ) xe = nDim0-1;
      ys = ystart;
      if ( ys < 0 ) ys = 0;
      ye = yend;
      if ( ye >= nDim1 ) ye = nDim1-1;
      if (3 < m_nVerbose)
        {
          cout << "Filling rectangle " << xs << ' ' << xe << ' ' << ys << ' ' << ye << endl
               << " nEdgeSize is " << nEdgeSize << " fFlag is " << fFlag << endl;
        }
      for (j = ys; j <= ye; j++)
        {
          nNumGood = 0;
          nNumBad  = 0;
          for (i = xs; i <= xe; i++)
            {
              fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(i, j);
              if (fPixIJ != fFlag)
                {
		  // Just copy it over (but wasn't this done previously?
                  (poImageOut->*poImageOut->prnSetPixel)(i, j, fPixIJ); // This is S L O O O W
                }
              else
                {
//==========================================================
                  // Flagged pixel in input array found, so interpolate/extrapolate

                  nl = 0;
                  for (jj = max(y0,j-nEdgeSize); jj <= min(yn,j+nEdgeSize); jj++)
                    {
                      k = 0;
                      for (ii = max(x0,i-nEdgeSize); ii <= min(xn,i+nEdgeSize); ii++)
                        {
                          fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(ii, jj);
                          if (fPixIJ != fFlag)
                            {
                              pfX[k] = (float)ii;
                              pfY[k] = fPixIJ;
                              if (ii != i)
//                                pfW[k] = min( 1./fabsf(ii-i), 1.);
                                pfW[k] = 1. / fabsf((float)(ii-i));
                              else
                                pfW[k] = 1.0;
                              k++;
                            }
                        } // end ii loop

                      if ( k >= 3 )
                        {
                          if ( k <= 5 )
                            {
                              vFitLine( pfX-1, pfY-1, pfW-1, k, afC-1);
                              t2  = afC[0] + afC[1] * (float)i;
                            }
                          else
                            {
                              vFitPoly2( pfX-1, pfY-1, pfW-1, k, afC-1 );
                              t2  = afC[0] + afC[1] * (float)i
                                           + afC[2] * (float)(i * i);
                            }

                          pfXl[nl] = jj;
                          pfYl[nl] = t2;
                          if (jj != j)
//                            pfWl[nl] = min( 1./fabsf(jj-j), 1.);
                            pfWl[nl] = 1. / fabsf((float)(jj-j));
                          else
                            pfWl[nl] = 1.0;
                          nl++;
                        }
                    } // end jj loop
                     
                  if ( nl >= 3 )
                    {
                      if ( nl <= 5 )
                        {
//                          call pf1 ( xl, yl, wl, nl, c );
                          vFitLine( pfXl-1, pfYl-1, pfWl-1, nl, afC-1 );
                          t2  = afC[0] + afC[1] * (float)j;
                        }
                      else
                        {
//                          call pf2 ( xl, yl, wl, nl, c );
/*
                          cout << "fitpoly2b\n";
                          cout << "n1, X, Y, W: " << nl << '\n';
                          int kjl;
                          for (kjl = 0; kjl < nl; kjl++)
                            cout << pfXl[kjl] << ' ';
                          cout << '\n';
                          for (kjl = 0; kjl < nl; kjl++)
                            cout << pfYl[kjl] << ' ';
                          cout << '\n';
                          for (kjl = 0; kjl < nl; kjl++)
                            cout << pfWl[kjl] << ' ';
                          cout << '\n';
*/
                          vFitPoly2( pfXl-1, pfYl-1, pfWl-1, nl, afC-1 );
                          t2  = afC[0] + afC[1] * (float)j
                                       + afC[2] * (float)(j * j);
//                          cout << "t2, afC: " << t2 << ' ' << afC[0] 
//                                                    << ' ' << afC[1] << ' ' <<  afC[2] << endl;
                        } 
                      if (5 < m_nVerbose)
                        {
                          cout << "Setting (i,j) " << setw(4) << i << ' '
                               << setw(4) << j << " to t2: " << t2 << endl; 
                        }
                      // Round up? TODO: If float, don't round up.
                      // If int, use nint function?

                      (poImageOut->*poImageOut->prnSetPixel)(i, j, t2);
//                      (poImageOut->*poImageOut->prnSetPixel)(i, j, t2 + 0.5);
/*
                      fPixIJ = (poImageOut->*poImageOut->prfGetPixel)(i, j);
                      if (fPixIJ != t2)
                      if (eImage_realIEEE != poImageOut->nGetDataType())
                        {
                          cout << "i,j, t2, fPixIJ: " << i << ' ' << j
                               << ", " << t2 << ", " << fPixIJ << endl;
                        }
*/
                      nNumGood++;
                    }
                  else
                    {
                      nNumBad++;
                    }
//==========================================================
                }
            } // end i loop
          if (5 < m_nVerbose)
            {
              cout << " Line " <<  j << " Fixed : " << nNumGood
                   << " Not Fixed : " << nNumBad << endl;
            }
        } // end j loop
    }
  else
    {

//     yend <= 0 => fill circular region

      xcenter    = (float)xstart;
      ycenter    = (float)ystart;
      max_radius = (float)xend;
      mr2        = max_radius * max_radius;
      if (3 < m_nVerbose)
        {
          cout << "Circle " <<  xcenter << ' ' << ycenter << ' ' << max_radius << endl;
          cout << "Bad Pixel Flag " << fFlag << endl;
        }
      for (j = 0; j < nDim1; j++)
        {
          y2 = (float)j - ycenter;
          y2 = y2 * y2;
          nNumGood = 0;
          nNumBad  = 0;
          for (i = 0; i < nDim0; i++)
            {
              rad2 = (float)i - xcenter;
              rad2 = rad2 * rad2 + y2;
              if ( rad2 < mr2 )
                {
                  fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(i, j);
                  if (fPixIJ != fFlag)
                    {
                      (poImageOut->*poImageOut->prnSetPixel)(i, j, fPixIJ); // This is S L O O O W
                    }
                  else
                    {
//==========================================================
                      // Flagged pixel in input array found, so interpolate/extrapolate

                      nl = 0;
                      for (jj = max(y0,j-nEdgeSize); jj <= min(yn,j+nEdgeSize); jj++)
                        {
                          k = 0;
                          for (ii = max(x0,i-nEdgeSize); ii <= min(xn,i+nEdgeSize); ii++)
                            {
                              fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(ii, jj);
                              if (fPixIJ != fFlag)
                                {
                                  pfX[k] = (float)ii;
                                  pfY[k] = fPixIJ;
                                  if (ii != i)
                                    // pfW[k] = min( 1./fabsf(ii-i), 1.);
                                    pfW[k] = 1. / fabsf((float)(ii-i));
                                  else
                                    pfW[k] = 1.0;
                                  k++;
                                }
                            }
                          if ( k >= 3 )
                            {
                              if ( k <= 5 )
                                {
                                  vFitLine( pfX-1, pfY-1, pfW-1, k, afC-1 );
                                  t2  = afC[0] + afC[1] * (float)i;
                                }
                              else
                                {
                                  vFitPoly2( pfX-1, pfY-1, pfW-1, k, afC-1 );
                                  t2  = afC[0] + afC[1] * (float)i
                                               + afC[2] * (float)(i * i);
                                }
                              
                              pfXl[nl] = jj;
                              pfYl[nl] = t2;
                              if (jj != j)
                                // pfWl[nl] = min( 1./fabsf(jj-j), 1.);
                                pfWl[nl] = 1. / fabsf((float)(jj-j));
                              else
                                pfWl[nl] = 1.0;
                              nl++;
                            }
                        }

                      if ( nl >= 3 )
                        {
                          if ( nl <= 5 )
                            {
//                          call pf1 ( xl, yl, wl, nl, c );
                              vFitLine( pfXl-1, pfYl-1, pfWl-1, nl, afC-1 );
                              t2  = afC[0] + afC[1] * (float)j;
                            }
                          else
                            {
//                          call pf2 ( xl, yl, wl, nl, c );
                              vFitPoly2( pfXl-1, pfYl-1, pfWl-1, nl, afC-1 );
			      cout << "FITPOLY2;;;;;\n";
                              t2  = afC[0] + afC[1] * (float)j
                                           + afC[2] * (float)(j * j);
                            }
                          if (5 < m_nVerbose)
                            {
                              cout << "Setting (i,j) " << setw(4) << i << ' '
                                   << setw(4) << j << " to t2: " << t2 << endl; 
                            }
                          (poImageOut->*poImageOut->prnSetPixel)(i, j, t2);
                          nNumGood++;
                        }
                      else
                        {
                          nNumBad++;
                        }
//C==========================================================
                    }
                }
            }
          if (5 < m_nVerbose)
            {
              cout << " Line " <<  j << " Fixed : " << nNumGood
                   << " Not Fixed : " << nNumBad << endl;
            }
        }
    }

  delete [] pfX;
  delete [] pfY;
  delete [] pfW;
  delete [] pfXl;
  delete [] pfYl;
  delete [] pfWl;
}

void
Ccalibrate::vFitPoly2(const float *pfX, const float *pfY, const float *pfW,
                      const int nNumPoints, float *pfC)
// This routine is translated directly from the Fortran.
// In order to circumvent the array indexing differences between
// Fortran and C++, pass  (pfX-1) so that intern pfX[1] = extern afX[0]
{
/*
C-----------------------------------------------------------------------
C     
C      ====================================
C      pf2 ( x, y, w, n, c )
C      ====================================
C      This routine fits 2nd order poly
C
C      Input :
C          x, y, w - point and weight
C          n - number of points
C
C      Output :
C          c - coefficients
C
C-----------------------------------------------------------------------
      subroutine pf2 ( x, y, w, n, c)
*/      
  double d, d1, d2, d3;
  int    i;
      
  double   Sy, Sxy, Sx2y,
           S, Sx, Sx2, Sx3, Sx4,
           xt, yt, wt, x2;

  Sy   = 0.0;
  Sxy  = 0.0;
  Sx2y = 0.0;
  S    = 0.0;
  Sx   = 0.0;
  Sx2  = 0.0;
  Sx3  = 0.0;
  Sx4  = 0.0;

  for (i = 1; i <= nNumPoints; i++)
    {
      xt    = pfX[i];
      yt    = pfY[i];
      wt    = pfW[i];
      x2    = xt * xt;
      Sy   += wt * yt;
      Sxy  += wt * xt * yt;
      Sx2y += wt * x2 * yt;
      S    += wt;
      Sx   += wt * xt;
      Sx2  += wt * x2;
      Sx3  += wt * xt * x2;
      Sx4  += wt * x2 * x2;
    }

  d  = S*(Sx2*Sx4-Sx3*Sx3) 
          + Sx*(Sx3*Sx2-Sx4*Sx) 
          + Sx2*(Sx*Sx3-Sx2*Sx2);
  d1  = Sy*(Sx2*Sx4-Sx3*Sx3) 
          + Sx*(Sx3*Sx2y-Sx4*Sxy) 
          + Sx2*(Sxy*Sx3-Sx2y*Sx2);
  d2  = S*(Sxy*Sx4-Sx2y*Sx3) 
          + Sy*(Sx3*Sx2-Sx4*Sx) 
          + Sx2*(Sx*Sx2y-Sx2*Sxy);
  d3  = S*(Sx2*Sx2y-Sx3*Sxy) 
          + Sx*(Sxy*Sx2-Sx2y*Sx) 
          + Sy*(Sx*Sx3-Sx2*Sx2);

  pfC[1] = (float)(d1 / d);
  pfC[2] = (float)(d2 / d);
  pfC[3] = (float)(d3 / d);
}
      
void
Ccalibrate::vFitLine(const float *pfX, const float *pfY, const float *pfW,
                     const int nNumPoints, float *pfC)
// This routine is translated directly from the Fortran.
// In order to circumvent the array indexing differences between
// Fortran and C++, pass  (pfX-1) so that intern pfX[1] = extern afX[0]
{
/*
C-----------------------------------------------------------------------
C     
C      ====================================
C      pf1 ( x, y, w, n, c )
C      ====================================
C      This routine fits line
C
C      Input :
C          x, y, w - point and weight
C          n - number of points
C
C      Output :
C          c - coefficients
C
C-----------------------------------------------------------------------
      subroutine pf1 ( x, y, w, n, c)
*/      
  int i;
  double  d, d1, d2,
          Sy, Sxy,
          S, Sx, Sx2;

  float wt;

  Sy  = 0.0;
  Sxy = 0.0;
  S   = 0.0;
  Sx  = 0.0;
  Sx2 = 0.0;

  for (i = 1; i <= nNumPoints; i++)
    {
      wt   = pfW[i];
      Sy  += wt * pfY[i];
      Sxy += wt * pfX[i] * pfY[i];
      S   += wt;
      Sx  += wt * pfX[i];
      Sx2 += wt * pfX[i] * pfX[i];
    }

  d  = S  * Sx2  -  Sx * Sx;
  d1 = Sy * Sx2  -  Sxy* Sx;
  d2 = S  * Sxy  -  Sx * Sy; 

  pfC[1] = (float)(d1 / d);
  pfC[2] = (float)(d2 / d);
}

void
Ccalibrate::vReverseUnpinPoint(tagDistort_parameters *ptDistorParams,
                               const float xu, const float yu, 
                               float *pfxo, float *pfyo)
{
/*
C--------------------------------------------------------------------------

      subroutine reverse_unpin_point ( dparams, xu, yu, xo, yo )

C--------------------------------------------------------------------------
C
C Calculate new "reverse unpinned" coodinates (xo, yo) given the 
C unpinned cordinates ( xu,yu ) and the pin cushion parameters dparams
C
C  xd = x - xc; yd = y - yc                      xd, yd : dummys
C
C  r = sqrt ( xd**2 + yd**2 )
C
C                2      3
C  ru = a*ro + b*ro + c*ro               ru = 'unpinned' radius
C                                        ro = 'observed' radius
C To solve equation for ro:
C      Solve equation of the form :
C                      2    3
C              A + Bx + Cx + Dx = 0
C
C      Using the Newton-Raphson method
C      -------------------------------
C                     f(x )
C      x    =   x  -      n                (if f(x) is deferentiable )
C       n+1      n      --------
C                     f'(x ) 
C                        n
C
C--------------------------------------------------------------------------
C
C                      2    3                  2      3
C   f(x)  = A + Bx + Cx + Dx = -r + a*re + b*re + c*re   re = estimate radius
C                        2                         2
C   f'(x) = B + 2Cx + 3Dx    =  a + 2*b*re + 3*c*re
C
C--------------------------------------------------------------------------
*/
  float a, b, c,
        re, re2, ro, ru;

  float xo, yo;

  int   tries;
  float dv;
  
  a = ptDistorParams->a * ptDistorParams->a1;
  b = ptDistorParams->b;
  c = ptDistorParams->c;
      
  // Calculate an 'unpinned' radius

  xo = xu-ptDistorParams->x_center;
  yo = (yu-ptDistorParams->y_center) / ptDistorParams->ratio;

  // If linear, can solve directly

  if ( (0.0 == b) && (0.0 == c) )
    {
      *pfxo = ptDistorParams->x_center +  xo / a;
      *pfyo = ptDistorParams->y_center +  yo / a;
      return;
    }

  ru = sqrtf(xo * xo  +  yo * yo);
      
  // Set initial guess at 'observed' radius to 'unpinned' radius

  re  = ru;
  re2 = re * re;
      
  // Calculate new 'observed' radius

  dv = a + 2.0 * b * re + 3.0 * c * re2;
  if (0.0 == dv)
    {
      // Error - set xo, yo to 0

      *pfxo = 0.0;
      *pfyo = 0.0;
      return;
    }

  ro = re - (-ru +  a * re +  b * re2 +  c * re * re2) / dv;
      
  // While there is a difference between the guess and the
  //  observed radius, keep guessing and recalculating

  tries = 1;
  while (fabsf( ro - re ) > 0.0003 )
    {
      tries++;
      re  = ro;
      re2 = re * re;
      dv  = a + 2.0 * b * re + 3.0 * c * re2;
      if ( ( 0.0 == dv ) || ( tries > 6 ) )
        {
          // Error - set xo, yo to 0

          *pfxo = 0.0;
          *pfyo = 0.0;
          return;
        }
      ro = re - (-ru +  a * re  +  b * re2  +  c * re * re2) / dv;
    } // end while

  // Success - set r to ratio between 'unpinned' and 
  // observed radius, and calulate observed x and y

  ru = ro / ru;
  *pfxo = ptDistorParams->x_center +  xo * ru;
  *pfyo = ptDistorParams->y_center +  yo * ru;
}
      

void
Ccalibrate::vCalcMaskPosition(tagCorrect_parameters *ptCorrectParams,
                              const float x, const float y, float *pfXm, float *pfYm)
{
/*
C------------------------------------------------------------------------

      subroutine calc_maskposition ( cparams, x, y, xm, ym )

C------------------------------------------------------------------------
C
C     given a calibrated pixel position, calculate
C     the position in mask coordinates
C
C------------------------------------------------------------------------
*/
  float del;
  float sinv, cosv,
        sinh, cosh,
        xt, yt;

  del  = sqrtf(1.0 + ptCorrectParams->horz_slope * ptCorrectParams->horz_slope);
  sinh = ptCorrectParams->horz_slope / del;
  cosh = 1.0 / del;
  del  = sqrtf(1.0 + ptCorrectParams->ver_slope * ptCorrectParams->ver_slope);
  sinv = ptCorrectParams->ver_slope / del;
  cosv = 1.0 / del;
      
  del  = cosh * cosv  -  sinh * sinv;
  xt   = x-ptCorrectParams->x_center;
  yt   = y-ptCorrectParams->y_center;
      
  *pfXm = (xt * cosv - yt * sinv) / del / ptCorrectParams->x_scale
          + ptCorrectParams->x_pt_center;
  *pfYm = (cosh * yt - sinh * xt) / del / ptCorrectParams->y_scale
          + ptCorrectParams->y_pt_center;
}

float
Ccalibrate::fBSpline2 ( const float lt[4][4], const float x, const float y )
{
  float afLine[4];
  int        j;
//  float t;

  for (j = 0; j < 4; j++)
    afLine[j] = fBSpline1 ( lt[j], x );
  return (fBSpline1 ( afLine, y ));
//  cout << "fBS2: ";
//  vListMatMN(4, 4, (float *)lt);
//  cout << "x, y: " << x << ' ' << y << endl;
//  cout << "Ans: " << t << endl;
//  return (t);
}

float
Ccalibrate::fBSpline1 (const float G[4], const float T)
{
  int   i, j;
  static float M[4][4] = {
    {-0.1666667,   0.5, -0.5, 0.1666667},
    { 0.5,        -1.0,  0.0, 0.6666667},
    {-0.5,         0.5,  0.5, 0.1666667},
    {0.1666667,    0.0,  0.0, 0.0} };

  float MG[4];
  float fTemp;

//  cout << "fBS1, G, T: " << G[0] << ' '
//       << G[1] << ' '<< G[2] << ' '<< G[3] << ' ' << T << '\n';

  for (j = 0; j < 4; j++)
    {
      MG[j] = 0.0;
      for (i = 0; i < 4; i++)
        {
          MG[j] += + M[i][j] * G[i];
        }
    }
  fTemp =  MG[0] * T * T * T + MG[1] * T * T + MG[2] * T + MG[3];
//  cout << "Ans: " << fTemp << endl;
  return (fTemp);
}

int
Ccalibrate::nInBox(const float xn, const float yn,
           const float x1, const float y1, const float x2, const float y2,
           const float x3, const float y3, const float x4, const float y4)
{
/*
C     ================================================
      integer function inbox (xn, yn, x1, y1, x2, y2,
     $                                x3, y3, x4, y4 )
C     ================================================
      
C     Determine if the point xn, yn is within a box determined by
C     the points (x1,y1) -> (x4,y4) where
C     (x1,y1) :: bottom left
C     (x2,y2) :: bottom right
C     (x3,y3) :: top left
C     (x4,y4) :: top right
C     
C     Returns the sides it is outside of by the flags
C     
C     bit 0      left      (x1,y1) -> (x3,j3)
C     bit 1      up        (x3,y3) -> (x4,y4)
C     bit 2      right     (x4,y4) -> (x2,y2)
C     bit 3      down      (x2,y2) -> (x1,y1)
C     bit 4      ERROR
C     
C     note that 2 bits could be set
C     
*/
  int nInBoxStat;

  float m, b;
      
  // For all of these fit a line to the two points that make up the
  // edge, then determine which side of the point it falls on
      
      
  //  Left edge
  //  =========

      m = y3 - y1;

      if ( m == 0.0 )
        {
          nInBoxStat = 16;
          return (nInBoxStat);

          //C     call errout 
          //C     $      ('ERROR in INBOX, box points have the same Y value',0,0)
          //C     write (6,*) '(x1, y1) = ',x1, y1
          //C     write (6,*) '(x3, y3) = ',x3, y3
        }
      else
        {
          m = ( x3 - x1 ) / m;
          b = x1 - m * y1;
          if ( xn < m*yn+b )
            {
              nInBoxStat = 1;
            }
          else
            {
              nInBoxStat = 0;
            }
          //C     nInBoxStat = fabsf(xn .lt. m*yn+b)
        }

  // Top edge
  // ========

  m = x4 - x3;
  if ( m == 0.0 )
    {
      nInBoxStat = 16;
      return (nInBoxStat);
      //C     call errout 
      //C     $      ('ERROR in INBOX, box points have the same x value',0,0)
      //C     write (6,*) '(x3, y3) = ',x3, y3
      //C     write (6,*) '(x4, y4) = ',x4, y4
    }
  else
    {
      m = (y4-y3) / m;
      b = y3 - m * x3;
      //C     nInBoxStat = nInBoxStat + 2*fabsf(yn .gt. m*xn+b)
      if (yn > m*xn+b) nInBoxStat = nInBoxStat + 2;
    }

  if ( nInBoxStat > 2 ) return (nInBoxStat);
      
  // right edge
  // =========

  m = y4-y2;
  if ( m == 0.0 )
    {
      nInBoxStat = 16;
      return (nInBoxStat);

      //C     call errout 
      //C     $      ('ERROR in INBOX, box points have the same Y value',0,0)
      //C     write (6,*) '(x4, y4) = ',x4, y4
      //C     write (6,*) '(x2, y2) = ',x2, y2
    }
  else
    {
      m = (x4-x2) / m;
      b = x2 - m*y2;
      //C     nInBoxStat = nInBoxStat + 4*fabsf(xn .gt. m*yn+b)
      if (xn > m*yn+b) nInBoxStat = nInBoxStat + 4;

      // if it is both to the left of the left edge and to the
      // right of the right line it must either be above or below
      // the box.  Keep above flag if set.  If it is below it will show 
      // up on the next test

      if ( (4 == (nInBoxStat & 4)) && ( (nInBoxStat & 1) == 1) )
        nInBoxStat = nInBoxStat & 2;
    }
  if ( nInBoxStat > 4 ) return (nInBoxStat);
      
  // bottom edge
  // ========

  m = x2-x1;
  if ( m == 0.0 )
    {
      nInBoxStat = 16;
      return (nInBoxStat);
      //C     call errout 
      //C     $      ('ERROR in INBOX, box points have the same X value',0,0)
      //C     write (6,*) '(x1, y1) = ',x1, y1
      //C     write (6,*) '(x2, y2) = ',x2, y2
    }
  else
    {
      m = (y2-y1) / m;
      b = y1 - m*x1;
      //C     nInBoxStat = nInBoxStat + 8*fabsf(yn .lt. m*xn+b)

      // As with the previous test, if it is both below and above,
      // it must be to the left or right.  Keep the right or left
      // flag by anding with 5

      if (yn < m*xn+b) nInBoxStat = nInBoxStat + 8;
      if ( (2 == (nInBoxStat & 2)) && (8 == (nInBoxStat & 8)) ) 
        nInBoxStat = nInBoxStat & 5;
    }
      
  return (nInBoxStat);
}

int
Ccalibrate::nInterpolate(tagCorrect_parameters *ptCorrectParams,
                         Cimage *poImgDistorXINT, Cimage *poImgDistorYINT,
                         Cimage *poImageIn, Cimage *poImageOut)
{
// This now calculates the TRANSFORM (m_poImgTransform) image.
/*
C     Correct (or distort) an entire image using interpolation
C     lookup tables
C
C     For each pixel in the input image, calculate where the
C     left, right, top and bottom edges will fall in the output
C     image.  Then divide intensity up by the relative overlap
C     of the transformed input pixel on the output image.
C
C
      subroutine interpolate ( input, output, 
     $     xasize, yasize, xsize, ysize,
     $     x_int, xxasize, xyasize, xxsize, xysize,
     $     y_int, yxasize, yyasize, yxsize, yysize, 
     $     xstart, xstep, ystart, ystep, pscale, ioffst, work )

C-----------------------------------------------------------------------      
C     to convert from (x,y) to (x',y')
C
C      i = (x-xstart) / xstep + 1
C      j = (y-ystart) / ystep + 1
C      x' = x_inv (i, j) * pscale   // DANGER should be x_inv(i-1,j-1) in C++ version!
C      y' = y_inv (i, j) * pscale   // DANGER should be y_inv(i-1,j-1) in C++ version!
C
C     where i,j are not necessarily integers, and interpolation
C     is required if not.
C
C-----------------------------------------------------------------------      
*/
  int    i, j;
  int    ii, jj;
  float  fPixOut;
  int    ixs, iys, ixe, iye;
  float  valu, valc;
  float  xp, yp;
  int    ix, iy;
  int    ioffst, icount;

  int    xasize, yasize, xsize, ysize;
  int    xxasize, xyasize, xxsize, xysize;
  int    yxasize, yyasize, yxsize, yysize;
  int    xstart, xstep, ystart, ystep;
  float  pscale;
  int    error;

  float    area;
  float    pix_per_area, xtemp, ytemp;
  float    fSatVal;
  float    xn, yn;
  int    xs, ys;
  float  *pfY0;
  float  x0, fIntBadFlag;
  float fCalib;
  float fCalibMax = 0.0;
  float fPost, fPostMin, fPostMax, fPostArea;

  LONG lIntBadFlag;
  LONG lYINT[8];
  LONG lXINT[8];

  cout << "\n***   Correct an image for spatial distortions..." << endl;

  if ( (NULL == poImageIn) || (NULL == poImageOut) )
    {
      return (-1);
    }
  else if (!poImageIn->bIsAvailable() || !poImageOut->bIsAvailable())
    {
      return (-2);
    }
  else if (NULL == poImgDistorXINT)
    {
      return (-3);
    }
  else if (!poImgDistorXINT->bIsAvailable())
    {
      return (-4);
    }
  else if (NULL == poImgDistorYINT)
    {
      return (-3);
    }
  else if (!poImgDistorYINT->bIsAvailable())
    {
      return (-4);
    }

  // Get the nonunf info

  // Danger here, you need a nonunf AND a dark image to do this, even
  // if you are not going to correct for nonunf and/or dark!

  // To save memory, delete member dark image and nonunf images if available
  // since they will not be used anymore directly.
 
  if (NULL != m_poImgNonunf)
    {
      delete m_poImgNonunf;
      m_poImgNonunf = NULL;
    }

  if (NULL != m_poImgDark)
    {
      delete m_poImgDark;
      m_poImgDark = NULL;
    }

  float   fNonunf, fBadNonunf;

  if (NULL == m_poNonunf) 
    m_poNonunf = new Cnonunf(sTransSymbol(m_sNonunfName), sTransSymbol(m_sDarkName));

  if (!m_poNonunf->bIsAvailable())
    {
      delete m_poNonunf;
      cout << "ERROR getting Nonunf object!" << endl;
      return (3);
    }

  m_poNonunf->nList();
  fBadNonunf = m_poNonunf->fGetBadFlag();

  icount      = 0;
  fIntBadFlag = (float) ptCorrectParams->badflag;
  lIntBadFlag = (long)  ptCorrectParams->badflag;

  fSatVal = poImageIn->fGetSatValue();

  (void) poImageIn->nGetDimensions(&xasize, &yasize);
  xsize   = xasize;
  ysize   = yasize;

  (void) poImgDistorXINT->nGetDimensions(&xxasize, &xyasize);
  xxsize  = xxasize;
  yysize  = xyasize;
  yxasize = xxasize;
  yyasize = xyasize;

  xstart  = ptCorrectParams->xint_start;
  ystart  = ptCorrectParams->yint_start;
  xstep   = ptCorrectParams->xint_step;
  ystep   = ptCorrectParams->yint_step;
  pscale  = ptCorrectParams->pscale;

  if (3 < m_nVerbose)
    {
      cout <<   "      Image Size            : " << xsize   << ", " << ysize
           << "\n      Image Array Size      : " << xasize  << ", " << yasize
           << "\n      Int Table Size        : " << xxsize  << ", " << yysize
           << "\n      Int Table Array Size  : " << xyasize << ", " << yyasize
           << "\n      Int X start, step     : " << xstart  << ", " << xstep
           << "\n      Int Y start, step     : " << ystart  << ", " << ystep
           << "\n      Pixel scale, 1/pscale : " << pscale  << ", " << 1.0 / pscale
           << endl;
    }

  ioffst = 0;
  if (  (NULL != m_poImgNonunf) && m_poImgNonunf->bIsAvailable())
    {
      m_poImgNonunf->m_oHeader.nGetValue(D_K_NonunfOffset, &ioffst);
    }
  cout << "      Image offset to use (from NONUNF header): " << ioffst << endl;

  // See if the post transform image is available and use it!

  fPostMin          = 5.0;
  fPostMax          = 0.0;
  Cimage *poImgPost = NULL;
  if ("" != m_sPostName)
    {
      poImgPost = new Cimage(m_sPostName);
      if (!poImgPost->bIsAvailable())
        {
          cout << "ERROR reading posttransform info!";
          delete poImgPost;
          poImgPost = NULL;
        }
    }

  // Start by setting output image to the offset value
  // (This could be done faster another way, but would be
  //  specific for image data types)

  float fOffset;
  fOffset = (float)ioffst;
  for (j = 0; j < yasize; j++)
    {
      for (i = 0; i < xasize; i++)
        {
          (poImageOut->*poImageOut->prnSetPixel)(i, j, fOffset);
        }
    }

  xs = min ( xsize, xstep*xxsize + xstart);
  ys = min ( ysize, ystep*yysize + ystart);

/*
C     Define lower edges for first row.  The lower edges for the
C     rest of the rows will just be equal to the upper edge of
C     the previous row.
C
C     Initialize y starting positions by assuming that the lower height
C     of the box is the same as the upper height of the box
C     This is to prevent getting a bad value from trying
C     an illegal position
*/
  int   nOffset;
  float fPixIJ;

  pfY0    = new float [xs+1];
  for (i = 1; i <= xs; i++)
    {
      xp = float(i-xstart)/(float)xstep + 1.0;
      yp = (1.0+0.5-(float)ystart)/(float)xstep + 1.0;
      ix = (int)xp;
      iy = (int)yp;
      if (   (ix >= 1) && (ix < yxasize) 
          && (iy >= 1) && (iy < yyasize) ) 
        {
          nOffset  = ix +  iy * yxasize;
//          cout << "ix, iy, nOffset: " << ix << ", " << iy << ", " << nOffset << endl;
          lYINT[0] = poImgDistorYINT->lGetPixel(nOffset);  // (ix,   iy);
          lYINT[1] = poImgDistorYINT->lGetPixel(nOffset-1); //ix-1, iy);
          lYINT[2] = poImgDistorYINT->lGetPixel(nOffset-yxasize); //ix, iy-1);
          lYINT[3] = poImgDistorYINT->lGetPixel(nOffset-yxasize-1);//ix-1, iy-1);
/*
////////////////////////////////
  lYINT[4] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix,   iy);
  lYINT[5] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix-1, iy);
  lYINT[6] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix,   iy-1);
  lYINT[7] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix-1, iy-1);
          if (lYINT[0] != lYINT[4])
            cout << "YINT mismatch for get pix ix, iy\n";
          if (lYINT[1] != lYINT[5])
            cout << "YINT mismatch for get pix ix-1, iy\n";
          if (lYINT[2] != lYINT[6])
            cout << "YINT mismatch for get pix ix, iy-1\n";
          if (lYINT[3] != lYINT[7])
            cout << "YINT mismatch for get pix ix-1, iy-1\n";
/////////////////////////////////
*/
          if (   (lIntBadFlag == lYINT[0])
              || (lIntBadFlag == lYINT[1])
              || (lIntBadFlag == lYINT[2])
              || (lIntBadFlag == lYINT[3]) )
            {
              pfY0[i] = fIntBadFlag;
            }
          else
            {
              valu = (float)((yp-(float)iy)*
                   ((xp-(float)ix) * (float)lYINT[0]                 //y_int(ix+1,iy+1) 
                   + (1.0+(float)ix-xp) * (float)lYINT[1]) +         //y_int(ix,iy+1)) +
                   (1.0+(float)iy-yp)*
                   ((xp-(float)ix)* (float)lYINT[2]                  //y_int(ix+1,iy)
                   + (1.0+(float)ix-xp)* (float)lYINT[3]))           //y_int(ix,iy)
                    * pscale;
               yp = (float)(1-ystart)/(float)ystep + 1.0;
               iy = int(yp);
               valc = (float)((yp-(float)iy)*
                   ((xp-(float)ix)* (float)lYINT[0]                 //y_int(ix+1,iy+1)
                   + (1.0+(float)ix-xp) *(float)lYINT[1]) +            //y_int(ix,iy+1)) +
                   (1.0+(float)iy-yp)*
                   ((xp-(float)ix) * (float)lYINT[2]                //y_int(ix+1,iy)
                   + (1.0+(float)ix-xp)* (float)lYINT[3]))           // y_int(ix,iy)))
                             * pscale;
               pfY0[i] = 2.0 * valc  -  valu;
            }
        }
      else
        {
          pfY0[i] = fIntBadFlag;
        }
    }  // end i loop

  float fTemp;

  Creflnlist *poReflnlistTable;
  Crefln     *poRefln;
  poReflnlistTable = new Creflnlist();

  int nTableSize;
  int nInOffset, nOutOffset;
  int nFI_nInOffset;
  int nFI_nOutOffset;
  int nRefsWritten = 0;
  int nFile;

  // Use H for nFI_nInOffset, and K for nFI_nOutOffset

  nFI_nInOffset    =  poReflnlistTable->nGetFieldIndex(Creflnlist::ms_snH);
  nFI_nOutOffset   =  poReflnlistTable->nGetFieldIndex(Creflnlist::ms_snK);

  poRefln    = new Crefln(poReflnlistTable);

  bool bCreateTransform  = m_bCreateTransform;
  tagTransform_parameters tTransformParams;
  Cstring sScrName;
  if ( !bCreateTransform || ("" != sGetEnv("DTREK_NOTRANSFORM")))
    {
#ifndef CORRECT_ONLY
      cout << "No transform image will be calculated.\n" << flush;
#endif
      bCreateTransform = FALSE;
    }
  else if (bCreateTransform)
    {
      // Initialize building a transform image

      if (NULL != m_poImgTransform) 
        delete m_poImgTransform;

      // Allocate all memory needed for reflns in poReflnlistTable.

      nTableSize = max(xasize * 256, 1024 * 1024);
      nTableSize = max(xasize * 512, 1024 * 1024);
      nTableSize = max(xasize * 256, 1024 * 512);
      nTableSize = max(xasize * 256, 4 * 1024 * 512);
      nTableSize = max(xasize * 256, 4 * 1024 * 1024);

      cout << "\nSize of temporary transform reflnlist: " << nTableSize << '\n' << flush;
      int nStat = poReflnlistTable->nExpand(nTableSize);

      nFile = 3;
      sScrName = sTransSymbol("$(TRANSFORM_SCR)");
      int nBytes;
      nBytes = sScrName.length();

      (void) dskbow(&nFile, sScrName.string(), &nBytes, &nTableSize, &nStat);

      if (0 != nStat)
        {
          cout << "     Error opening file " << sScrName << "Error is: " << nStat << "\n";
          return (nStat);
        }
    }

  // Start of main loop

  for (j = 1; j <= ys; j++)
    {
/*
C        Define left edges for first pixel in this row.  The left
C        edges for the rest of the row will just be equal to the
C        right edge of the previous pixel.
C
C        as above, assume that the distance from the center to
C        the left edge is the same as the distance from the
C        center to the right edge
*/
         yp = float(j-ystart)/float(ystep) + 1.0;
         xp = float(1.0+0.5-(float)xstep)/float(xstep) + 1.0;
         ix = (int)xp;
         iy = (int)yp;
         if (   (ix >= 1) && (ix < xxasize) 
             && (iy >= 1) && (iy < xyasize) ) 
           {
             nOffset =         ix +  iy * xxasize;
             lXINT[0] = poImgDistorXINT->lGetPixel(nOffset);
             lXINT[1] = poImgDistorXINT->lGetPixel(nOffset-1);
             lXINT[2] = poImgDistorXINT->lGetPixel(nOffset-xxasize);
             lXINT[3] = poImgDistorXINT->lGetPixel(nOffset-xxasize-1);
/*
////////////////////////////////
  lXINT[4] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix,   iy);
  lXINT[5] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix-1, iy);
  lXINT[6] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix,   iy-1);
  lXINT[7] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix-1, iy-1);
          if (lXINT[0] != lXINT[4])
            cout << "XINTa mismatch for get pix ix, iy\n";
          if (lXINT[1] != lXINT[5])
            cout << "XINTa mismatch for get pix ix-1, iy\n";
          if (lXINT[2] != lXINT[6])
            cout << "XINTa mismatch for get pix ix, iy-1\n";
          if (lXINT[3] != lXINT[7])
            cout << "XINTa mismatch for get pix ix-1, iy-1\n";
/////////////////////////////////
*/
             if (   (lIntBadFlag == lXINT[0])
                 || (lIntBadFlag == lXINT[1])
                 || (lIntBadFlag == lXINT[2])
                 || (lIntBadFlag == lXINT[3]) )
               {
                 x0 = fIntBadFlag;
               }
             else
               {
                 
               valu = (float)((yp-(float)iy)*
                      ((xp-(float)ix)* (float)lXINT[0]             //x_int(ix+1,iy+1)
                       + (1.0+(float)ix-xp)* (float)lXINT[1]) +      // x_int(ix,iy+1)) +
                   (1.0+(float)iy-yp)*
                   ((xp-(float)ix)*        (float)lXINT[2]         //x_int(ix+1,iy)
                    + (1.0+(float)ix-xp)*     (float)lXINT[3]))      //x_int(ix,iy))) 
                 * pscale;
               xp = (float)(1-xstart)/(float)xstep + 1.0;
               ix = (int)xp;
               valc = (float)((yp-(float)iy)*
                   ((xp-(float)ix)* (float)lXINT[0]             //x_int(ix+1,iy+1)
                    + (1.0+(float)ix-xp)*  (float)lXINT[1]) +     //x_int(ix,iy+1)) +
                   (1.0+(float)iy-yp)*
                       ((xp-(float)ix)* (float)lXINT[2]         //x_int(ix+1,iy)
                        + (1.0+(float)ix-xp)* (float)lXINT[3]))    //x_int(ix,iy)))
                 * pscale;
               x0 = 2.0 * valc  -  valu;
               }
           }
         else
           {
             x0 = fIntBadFlag;
           }

         // Now go across the row

         for (i = 1; i <= xs; i++)
           {
             // Calculate right edge

             xp = float(i+0.5-xstart)/float(xstep) + 1.0;
             yp = float(j-ystart)/float(ystep) + 1.0;
             ix = (int)xp;
             iy = (int)yp;

             if (   (ix >= 1) && (ix < xxasize) 
                 && (iy >= 1) && (iy < xyasize) ) 
               {
                 nOffset =         ix +  iy * xxasize;
                 lXINT[0] = poImgDistorXINT->lGetPixel(nOffset); 
                 lXINT[1] = poImgDistorXINT->lGetPixel(nOffset-1);
                 lXINT[2] = poImgDistorXINT->lGetPixel(nOffset-xxasize);
                 lXINT[3] = poImgDistorXINT->lGetPixel(nOffset-xxasize-1);
/*
////////////////////////////////
  lXINT[4] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix,   iy);
  lXINT[5] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix-1, iy);
  lXINT[6] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix,   iy-1);
  lXINT[7] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix-1, iy-1);
          if (lXINT[0] != lXINT[4])
            cout << "XINTb mismatch for get pix ix, iy\n";
          if (lXINT[1] != lXINT[5])
            cout << "XINTb mismatch for get pix ix-1, iy\n";
          if (lXINT[2] != lXINT[6])
            cout << "XINTb mismatch for get pix ix, iy-1\n";
          if (lXINT[3] != lXINT[7])
            cout << "XINTb mismatch for get pix ix-1, iy-1\n";
/////////////////////////////////
*/
                 if (  (lIntBadFlag == lXINT[0])
                    || (lIntBadFlag == lXINT[1])
                    || (lIntBadFlag == lXINT[2])
                    || (lIntBadFlag == lXINT[3]) )
                   {
                     xn = fIntBadFlag;
                   }
                 else
                   {
                     xn = (float)((yp-(float)iy)*  
                           ((xp-(float)ix) * (float)lXINT[0]   //x_int(ix+1,iy+1)
                            + (1.0+(float)ix-xp) * (float)lXINT[1]) + //x_int(ix,iy+1)) +
                      (1.0+(float)iy-yp)*
                           ((xp-(float)ix) * (float)lXINT[2]      //x_int(ix+1,iy)
                            + (1.0+(float)ix-xp) * (float)lXINT[3])) //x_int(ix,iy)))
                       * pscale;
                   }
               }
            else
              {
                xn = fIntBadFlag;
              }
            

             // Calculate upper edge

             xp = float(i-xstart)/float(xstep) + 1.0;
             yp = ((float)j+0.5-(float)ystart)/float(ystep) + 1.0;
             ix = (int)xp;
             iy = (int)yp;

             if (   (ix >= 1) && (ix < yxasize) 
                 && (iy >= 1) && (iy < yyasize) ) 
               {
                 nOffset  = ix +  iy * yxasize;
                 lYINT[0] = poImgDistorYINT->lGetPixel(nOffset);
                 lYINT[1] = poImgDistorYINT->lGetPixel(nOffset-1);
                 lYINT[2] = poImgDistorYINT->lGetPixel(nOffset-xxasize);
                 lYINT[3] = poImgDistorYINT->lGetPixel(nOffset-xxasize-1);
/*
////////////////////////////////
  lYINT[4] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix,   iy);
  lYINT[5] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix-1, iy);
  lYINT[6] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix,   iy-1);
  lYINT[7] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix-1, iy-1);
          if (lYINT[0] != lYINT[4])
            cout << "YINTb mismatch for get pix ix, iy\n";
          if (lYINT[1] != lYINT[5])
            cout << "YINTb mismatch for get pix ix-1, iy\n";
          if (lYINT[2] != lYINT[6])
            cout << "YINTb mismatch for get pix ix, iy-1\n";
          if (lYINT[3] != lYINT[7])
            cout << "YINTb mismatch for get pix ix-1, iy-1\n";
/////////////////////////////////
*/

                 if (   (lIntBadFlag == lYINT[0])
                     || (lIntBadFlag == lYINT[1])
                     || (lIntBadFlag == lYINT[2])
                     || (lIntBadFlag == lYINT[3]) )
                   {
                     yn = fIntBadFlag;
                   }
                 else
                   {
                     yn = (float)((yp-(float)iy)*
                                  ((xp-(float)ix)*  (float)lYINT[0]
                                   + (1.0+(float)ix-xp) * (float)lYINT[1]) + 
                                  (1.0+(float)iy-yp)*
                                  ((xp-(float)ix) * (float)lYINT[2]    
                                   + (1.0+(float)ix-xp) * (float)lYINT[3]))
                       * pscale;
                   }
               }
             else
               {
                 yn = fIntBadFlag;
               }
/*
C     divide up the pixel
C     ix, iy are the index in the output image
C
C================================================================
C     July '94 - found a bug!  With radial distortion, an input
C     pixel always went to less than one pixel in the output image.
C     Without image intensifiers, the image could expand.  I had
C     assumed that the image always shrunk.  (It wasn't completely
C     stupid - making this assumption saved compute time.  It was
C     just stupid to forget.)
C
C     Ok, we have to move the input pixel into the area
C     between (x0,y0) -> (xn, yn)
C
C     Recall:
C
C     1) INTFLAG means that the position in invalid
C
C     2) pixels defined as going from i-0.5 -> i+0.5
C
C                               pixel (ixe,iye)
C     +-----------------------+
C     |       |       |(xn,yn)|
C     |  +=============+      |
C     |  |             |      |
C     |  |             |      |
C     ---|             |-------
C     |  |             |      |
C     |  |             |      |
C     |  |             |      |
C     |  |             |      |
C     ---|             |-------
C     |  +=============+      |
C     |(x0,y0)|       |       |
C     |       |       |       |
C     |   ^   |       |       |
C     +---|-------------------+
C     ^   |
C     |   --  pixel (ixs,iys)
C     |
C     |-- (ixs-0.5, iys-0.5)
C================================================================
C
C
C
*/
             if (   (x0 != fIntBadFlag) && (pfY0[i] != fIntBadFlag)
                 && (xn != fIntBadFlag) && (yn      != fIntBadFlag)
                 && (xn > 0.0) && (yn   > 0.0)
                 && (x0 > 0.0) && (pfY0[i] > 0.0) )
               {
                 area = (xn - x0) * (yn - pfY0[i]);
//                 if (0.0 >= area)
//                   cout << "area <= 0 for input pixel i,j: " << i << ", " << j << endl;
                 fNonunf = (m_poNonunf->*m_poNonunf->m_prfNonunfFactor)(i-1, j-1);
                 if ( (fNonunf == fBadNonunf) || (0.0 == fNonunf) )
                   area = 0.0;
                 if (0.0 < area)
                   {
                     // DO NOT USE nint function?!

                     ixs = nint(x0);        // Truncate down
                     iys = nint(pfY0[i]);   // Truncate down
                     ixe = nint(xn);        // Truncate down
                     iye = nint(yn);        // Truncate down
//                     ixe = int(xn+1.0);    // Round-up  (truncate down is probably OK)
//                     iye = int(yn+1.0);    // Round-up
/*
c+jwp 12-sep-1995
c if any pixel in the input is flagged as bad non-uniformity or is
c saturated (i.e. has a value of -1 (= 65535) then the output pixel is set
c to 65535, too.
c
*/
                     fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(i-1, j-1);
                     if ( (-1.0 == fPixIJ) || (fSatVal <= fPixIJ) )
                       {
                         // Input is bad, so flag pix_per_area as bad for 
                         // later on

                         //c+jwp 13-May-1996     pix_per_area = -1.0 - fOffset
                         pix_per_area = -1.0 - fOffset;
                         //c                     pix_per_area = 0.0
                         //c-jwp 13-May-1996
                         //c
                       }
                     else
                       {
                         // i4 = jzext(input(i,j))
                         // pix_per_area = (i4 - ioffst) / area

                         pix_per_area = (fPixIJ - fOffset) / area;
                         // the above may be negative!
                       }
//c-jwp
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     nInOffset = (xasize * (j-1))  + (i-1);
                                             
                     for (jj = iys; jj <= iye; jj++)
                       {
                         if ( (0 < jj) && (jj <= yasize) )
                           {
                             ytemp = min((float)jj+0.5, yn)
                                   - max((float)jj-0.5, pfY0[i]);
                             if (8 < m_nVerbose)
                               {
                                 if (0.0 >= ytemp)
                                   {
                                     cout << "ytemp <= 0 for output ii,jj: " 
                                          << ii << ", " << jj << "     "
                                          << ytemp << ", " << yn << ", " << pfY0[i] << endl;
                                   }
                               }
                             if (0.0 < ytemp)
                               {
                                 for (ii = ixs; ii <= ixe; ii++)
                                   {
                                     if ( (0 < ii) && (ii <= xasize) )
                                       {
                                         xtemp = min((float)ii+0.5, xn) 
                                               - max((float)ii-0.5, x0);
                                         if (8 < m_nVerbose)
                                           {
                                             if (0.0 >= xtemp)
                                               {
                                                 cout << "xtemp <= 0 for output ii,jj: " << ii << ", " << jj << "    "
                                                      << xtemp << ", " << xn << ", " << x0 << endl;
                                               }
                                           }
                                         if (0.0 < xtemp)
                                           {
                                             fPixOut
                    = (poImageOut->*poImageOut->prfGetPixel)(ii-1, jj-1);
                                     //c+jwp 12-sep-1995

                                             if (   (pix_per_area == (-1.0-fOffset))
                                                 || (fPixOut == -1.0) )
                                               {
//c flag output pixel as bad and keep it that way

                                                 fPixOut = -1.0;
                                               }
                                             else
                                               {
/*
c-jwp
c
c add the contribution of this input pixel to the output pixel scaled by
c the (area overlapped with output) / (area of input pixel).
c remember that output always starts with the ioffst in it, so it should
c have little chance to go below 0.  If it does, we print a warning message.
c In the meantime, add 0.5 to the floating point value, so that it rounds
c up instead of truncates.
c
*/
                                                 fTemp = xtemp * ytemp;

                                                 if (0.0001 < fTemp)
                                                   {
                                                     fPixOut = fPixOut 
                                                 + (pix_per_area*xtemp*ytemp);

/*
    cout << setw(6) << i  << setw(6) << j  << ' ' << setw(10) << area
         << setw(6) << ii << setw(6) << jj << ' ' << setw(10) << fTemp
         << '\n';
*/
                                                     if (bCreateTransform)
                                                       {
                                                         nOutOffset = (xasize * (jj-1)) + (ii-1);
                                                         fPostArea = area;
                                                         if (NULL != poImgPost)
                                                           {
                                                             // If it is available
                                                             // rescale by posttransform info 
                                                         
                                                             fPost = poImgPost->fGetPixel(nOutOffset);
                                                             if (0.0 < fPost)
                                                               {
                                                                 fPostArea = fPostArea / fPost;
                                                                 if (fPost < fPostMin)
                                                                   fPostMin = fPost;
                                                                 else if (fPost > fPostMax)
                                                                   fPostMax = fPost;
                                                               }
                                                             else 
                                                               {
                                                                 if (8 < m_nVerbose)
                                                                   {
                                                                     cout << "fPost problem: " << fPost 
                                                                          << ": " << ii-1 << ", " << jj-1
                                                                          << endl;
                                                                   }
                                                                 fTemp = 0.0;
                                                               }
                                                           }
                                                         fPostArea = fPostArea / fNonunf;

                                                         poRefln->vSetH(nInOffset);
                                                         poRefln->vSetK(nOutOffset);
                                                         poRefln->vSetIntensity(fPostArea);
                                                         poRefln->vSetSigmaI(fTemp);
                                                         fCalib = fTemp / fPostArea;
                                                         fCalibMax = max(fCalibMax, fCalib);
                                                         // TODO: if fCalib <= 0 do not insert?
                                                         if (0.0 < fCalib)
                                                           {
                                                             poReflnlistTable->nInsert(poRefln);
                                                             if (nTableSize 
                                                                 == poReflnlistTable->nGetNumReflns())
                                                               (void) nPreTransform(CAL_DUMPSOME, nFile,
                                                                                    poReflnlistTable, 
                                                                                    &nRefsWritten);
                                                           }
                                                       }

                                                     if (fPixOut >= fSatVal)
                                                       {
                                                         // it became legitimately saturated

                                                         fPixOut = -1.0;
                                                       }
                                                     else if ( fPixOut < -1.0)
                                                       {
                                                         // oops, how did it get so negative?
                                                         // better flag it as bad

                                                         if (icount < 100)
                                                           {
                                                             // allow only 100 error
                                                             // messages of this type
                                                             cout << "WARNING, negative "
                                                                  << "output at " << ii
                                                                  << ", " << jj << endl;
                                                             icount++;
                                                           }
                                                         fPixOut = -1.0;
                                                       }
                                                   }
                                               }
                                             (poImageOut->*poImageOut->prnSetPixel)(ii-1, jj-1, fPixOut);
                                             //c+jwp 12-sep-1995
                                           }
                                       } // if
                                   } //enddo
                               }
                           }
                       } //enddo
                   }
               }

             // Set the left edge for the next pixel and the upper
             // edge for the next row

             x0      = xn;
             pfY0[i] = yn;
           }// end i loop

       } // end j loop

  delete [] pfY0;

/*
c+jwp 13-May-1996
c Now subtract off the ioffset from everywhere, except the bad pixels
c Also truncate to 0
*/
  for (j = 1; j <= ys; j++)
    {
      for (i = 1; i <= xs; i++)
        {
          fPixIJ = (poImageOut->*poImageOut->prfGetPixel)(i-1, j-1);
          if ( (-1.0 != fPixIJ) && (fSatVal > fPixIJ) )
            {
              fPixIJ = fPixIJ - fOffset;
              if (0.0 > fPixIJ) fPixIJ = 0.0;
              (poImageOut->*poImageOut->prnSetPixel)(i-1, j-1, fPixIJ);
            }
        }
    }

  // Do we need to build the transform image?

  if ( (0 >= nRefsWritten) && (0 >= poReflnlistTable->nGetNumReflns()) )
    {
      if (bCreateTransform)
        cout << "SEVERE ERROR, no transformation information available!"
             << endl << flush;
      else
        cout << "INFO, No transformation information available!"
             << endl << flush;
    }
  else
    {
      // Make a transform image

      if (0 < poReflnlistTable->nGetNumReflns())
        (void) nPreTransform(CAL_DUMPALL, nFile, poReflnlistTable, &nRefsWritten);
      cout << "\nNumber of reflns in scratch file: "
           << nRefsWritten << endl;
      if ("" != sGetEnv("DTEXPOSE_SCALE"))
	{
	  float fExposeScale = 1.0;
	  fExposeScale = atof(sGetEnv("DTEXPOSE_SCALE").string());
	  cout << "Expose scale factor is " << fExposeScale << endl;
	  //fCalibMax = fCalibMax * fExposeScale;
	}
      cout << "Maximum Calib transform scale factor: " << fCalibMax << endl;

      // Save some memory

      if (NULL != poReflnlistTable)
        {
          delete poReflnlistTable;
          poReflnlistTable = NULL;
        }
      if (NULL != poImgPost)
        {
          cout << "MinPost , MaxPost:  " << fPostMin << ", " << fPostMax
               << endl;
          delete poImgPost;
          poImgPost = NULL;
        }

      // Create transform image 10% larger than inputs to allow for
      // bad output pixels that are not in the input information

      int nExpandDiv = 10;
      if ("" != sGetEnv("DTCAL_EXPAND"))
	{
	  nExpandDiv = atoi(sGetEnv("DTCAL_EXPAND").string());
	  if (0 >= nExpandDiv) nExpandDiv = 4;
	  cout << "DTCAL_EXPAND is " << nExpandDiv << endl;
	}

      int nMaxContribs = nRefsWritten + nRefsWritten / nExpandDiv;
      //      int nMaxContribs = nRefsWritten + nRefsWritten / 20;

      // cout << "nRefsWritten: " << nRefsWritten << endl;

      // nMaxContribs = nRefsWritten + nRefsWritten / 1;

      // Shouldn't the next few lines be in ::nBuildTransform(...)?

      m_nDimTrf = 2;
      cout << "TRANSFORM image size: " << m_nDimTrf << " by " 
           << nMaxContribs << '\n' << flush;
      m_poImgTransform = new Cimage(m_nDimTrf, nMaxContribs, eImage_uI2);

      // Make sure scale factors will fit in 2-bytes (0 - 65535), but
      // in reality (2-65535) since 0 and 1 are used as special flags

      fCalibMax = 65534.9 / fCalibMax;
      if (fCalibMax > 65536.0) fCalibMax = 65536.0;

      fCalibMax = 65536.0 / fCalibMax;

      cout << "CURRENT 1: " << fCalibMax << endl;
      if ("" != sGetEnv("DTEXPOSE_SCALE"))
	{
	  float fExposeScale = 1.0;
	  fExposeScale = atof(sGetEnv("DTEXPOSE_SCALE").string());
	  cout << "Expose scale factor: " << fExposeScale << endl;
	  fCalibMax = fCalibMax * fExposeScale * fExposeScale;
	  cout << "NEW Transform scale factor: " << fCalibMax << endl;
	}

      cout << "REPLACING 1: " << fCalibMax << endl;
      m_poImgTransform->m_oHeader.nReplaceValue("TRANSFORM_SCALE", 
                                                fCalibMax, 9);

      int nStat;

      // Close the scratch file

      (void) dskbcw(&nFile, &nStat);

//      cout << "About to call nBuildTransform...\n" << flush;

#ifndef CORRECT_ONLY
      nStat = nBuildTransform(sScrName, nRefsWritten, xasize, yasize);
#else
      nStat = -1;
#endif

      // Delete the scratch file

//      (void) nFileDelete(sScrName);
    }

  if (NULL != poImgPost)
    {
      cout << "MinPost , MaxPost:  " << fPostMin << ", " << fPostMax
           << endl;
      delete poImgPost;
      poImgPost = NULL;
    }

  if (NULL != poReflnlistTable)
    {
      delete poReflnlistTable;
      poReflnlistTable = NULL;
    }
  delete poRefln;
  return (0);
}

int
Ccalibrate::nReadDistor(void)
{
  // Read spatial distortion info from image files on disk

  int nStat;
  int i;

  Cimage **appoImage[4];

  appoImage[0] = &m_poImgDistorXINV;
  appoImage[1] = &m_poImgDistorYINV;
  appoImage[2] = &m_poImgDistorXINT;
  appoImage[3] = &m_poImgDistorYINT;

  nStat = 0;
  for (i = 0; (i < 4) && (0 == nStat); i++)
    {
      if (NULL != *appoImage[i])
        {
          delete *appoImage[i];
          *appoImage[i] = NULL;
        }

      *appoImage[i] = new Cimage (m_sDistorName + ms_asExt[i]);
      if (!(*appoImage[i])->bIsAvailable())
        {
          cout << "ERROR reading distor file " << m_sDistorName << ms_asExt[i]
               << endl;
          delete *appoImage[i];
          *appoImage[i] = NULL;
          nStat = i + 1;
        }
    }

  Cimage_header *poCalpar;
  poCalpar = new Cimage_header(m_sDistorName + ms_asExt[4]);
  if (!poCalpar->bIsAvailable())
    {
      cout << "ERROR reading distor file " << m_sDistorName << ms_asExt[i]
           << endl;
      nStat = nStat + 1;
    }

  if (0 == nStat)
    {
      int nTemp;
      Cstring sWarning = "WARNING in nReadDistor(), problem reading keyword: ";
      nTemp = poCalpar->nGetValue(D_K_CalibPixelSize, &m_fpixel_size);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibPixelSize << endl;
      nStat += nTemp;

      nTemp = poCalpar->nGetValue(D_K_CalibPscale,
                                  &m_tCorrectParams.pscale);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibPscale << endl;
      nStat += nTemp;

      nTemp = poCalpar->nGetValue(D_K_CalibXintStart,
                                  &m_tCorrectParams.xint_start);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibXintStart << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibYintStart,
                                  &m_tCorrectParams.yint_start);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibYintStart << endl;
      nStat += nTemp;
      //  nTemp = poCalpar->nGetValue(D_K_CalibXinvStart,
      //                 "XINV_START");
      //  nTemp = poCalpar->nGetValue(D_K_CalibYinvStart,
      //                 "YINV_START");
      nTemp = poCalpar->nGetValue(D_K_CalibXintStep,
                                  &m_tCorrectParams.xint_step);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibXintStep << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibYintStep,
                                  &m_tCorrectParams.yint_step);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibYintStep << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibXBeam,
                                  &m_tDistorParams.x_beam);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibXBeam << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibYBeam,
                                  &m_tDistorParams.y_beam);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibYBeam << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibBadFlag,
                                  &m_tCorrectParams.badflag);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibBadFlag << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibXScale,
                                  &m_tDistorParams.x_scale);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibXScale << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibYScale,
                                  &m_tDistorParams.y_scale);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibYScale << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibSpacing,
                                  &m_tDistorParams.spacing);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibSpacing << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibXCenter,
                                  &m_tDistorParams.x_center);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibXCenter << endl;
      nTemp = poCalpar->nGetValue(D_K_CalibYCenter,
                                  &m_tDistorParams.y_center);
//      nStat += nTemp;
      if (0 != nTemp)
        cout << sWarning << D_K_CalibYCenter << endl;
//      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibRadialA,
                                  &m_tDistorParams.a);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibRadialA << endl;
//      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibRadialB,
                                  &m_tDistorParams.b);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibRadialB << endl;
//      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibRadialC,
                                  &m_tDistorParams.c);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibRadialC << endl;
//      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibRadialA1,
                                  &m_tDistorParams.a1);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibRadialA1 << endl;
//      nStat += nTemp;

      nTemp = poCalpar->nGetValue(D_K_CalibRatio,
                                  &m_tDistorParams.ratio);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibRatio << endl;
//      nStat += nTemp;

      nTemp = poCalpar->nGetValue(D_K_CalibXinvStep, 
                                  &m_tCorrectParams.xinv_step);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibXinvStep << endl;
      nStat += nTemp;

      nTemp = poCalpar->nGetValue(D_K_CalibYinvStep,
                                  &m_tCorrectParams.yinv_step);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibYinvStep << endl;
      nStat += nTemp;

      nTemp = poCalpar->nGetValue(D_K_CalibXSize,
                                  &m_tDistorParams.x_size);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibXSize << endl;
      nStat += nTemp;
      nTemp = poCalpar->nGetValue(D_K_CalibYSize,
                                  &m_tDistorParams.y_size);
      if (0 != nTemp)
        cout << sWarning << D_K_CalibYSize << endl;
      nStat += nTemp;
/*
         nTemp = poCalpar->nGetValue(D_K_CalibXCenter,
         "X_CENTER");
         nTemp = poCalpar->nGetValue(D_K_CalibYCenter,
         "Y_CENTER");
         nTemp = poCalpar->nGetValue(D_K_CalibXPtCenter,
         "X_PT_CENTER");
         nTemp = poCalpar->nGetValue(D_K_CalibYPtCenter,
         "Y_PT_CENTER");

         nTemp = poCalpar->nGetValue(D_K_CalibVertSlope,
         "VER_SLOPE");
         nTemp = poCalpar->nGetValue(D_K_CalibHorzSlope,
         "HORZ_SLOPE");
         nTemp = poCalpar->nGetValue(D_K_CalibXMaskPoints,
         "X_MASK_POINTS");
         nTemp = poCalpar->nGetValue(D_K_CalibYMaskPoints,
         "Y_MASK_POINTS");
         nTemp = poCalpar->nGetValue(D_K_CalibLeftMaskPoint,
         "LEFT___MASK_POINT");
         nTemp = poCalpar->nGetValue(D_K_CalibCenterMaskPoint,
         "CENTER_MASK_POINT");
         nTemp = poCalpar->nGetValue(D_K_CalibRightMaskPoint,
         "RIGHT__MASK_POINT");
         nTemp = poCalpar->nGetValue(D_K_CalibBottomMaskPoint,
         "BOTTOM_MASK_POINT");
         nTemp = poCalpar->nGetValue(D_K_CalibTopMaskPoint,
         "TOP____MASK_POINT");
         */

      if (0 != nStat)
        {
          cout << "ERROR getting distor keywords from " << m_sDistorName 
               << ms_asExt[4]<< endl;
        }
    }
  if (0 != nStat)
    {
      cout << "ERROR reading spatial distortion information!" << endl;
    }
  delete poCalpar;
  return (nStat);
}

int
Ccalibrate::nGoInterp(void)
{
  // TODO:  This routine should call nCorrectImage!
  int nStat;
  if (NULL != m_poImgMask)
    delete m_poImgMask;

  m_poImgMask = new Cimage(m_sMaskName);

  if (!m_poImgMask->bIsAvailable())
    {
      cout << "ERROR, mask image unavailable!" << endl;
      nStat = -1;
      return (nStat);
    }
  
  Cimage *poImgMaskCor;
  poImgMaskCor = new Cimage(*m_poImgMask);
  if (0 == nReadDistor())
    {
      nStat = nInterpolate(&m_tCorrectParams,
                           m_poImgDistorXINT, m_poImgDistorYINT,
                           m_poImgMask,
                           poImgMaskCor);
      poImgMaskCor->nWrite("maskcor.img");
    }
  delete poImgMaskCor;
  return (nStat);
}

int
Ccalibrate::nWriteDistor(void)
{
  // Write out the spatial distortion information

  int nStat, nStatL;

  Cimage *apoImage[4];

  apoImage[0] = m_poImgDistorXINV;
  apoImage[1] = m_poImgDistorYINV;
  apoImage[2] = m_poImgDistorXINT;
  apoImage[3] = m_poImgDistorYINT;

  int i;

  nStat = 0;
  for (i = 0; (i < 4) && (0 == nStat); i++)
    {
      if ( (NULL != apoImage[i]) && (apoImage[i]->bIsAvailable()) )
        nStat = apoImage[i]->nWrite(m_sDistorName + ms_asExt[i]);
      else
        nStat = 2;
      if (0 != nStat)
        {
          cout << "ERROR writing " << m_sDistorName << ms_asExt[i]
               << endl;
        }
    }
  if (0 == nStat)
    {
      if (NULL == m_poHeaderSave)
        {
          m_poHeaderSave = new Cimage_header();
        }
      nUpdateHeader(m_poHeaderSave);
      if ( (NULL != m_poHeaderSave) && m_poHeaderSave->bIsAvailable() )
        {
          nStat = m_poHeaderSave->nWrite(m_sDistorName + ms_asExt[4]);
        }
      else
        {
          nStat = 3;
        }
      if (0 != nStat)
        {
          cout << "ERROR writing " << m_sDistorName << ms_asExt[4]
               << endl;
        }
    }
  return (nStat);
}

int
Ccalibrate::nGoNonunf(void)
{
  int nStat = 0;
  int i;

  cout << "\n***Calculate the non-uniformity correction image..." << endl;

  // 1. Make sure distor files are available
  // 2. Read in DARK image, find bad pixels in the DARK image
  // 3. Read in FLOOD image, find bad pixels in the FLOOD image
  // 4. Read in the bad pixel list
  // 5. Subtract DARK from FLOOD
  // 6. Calculate the reference image (or read it in)
  // 7. Distor the reference image
  // 8. Scale FLOOD to reference image to create NONUNF image
  // 9. Write out NONUNF image

  //****************** Find Bad Pixels *************************

  if (NULL == m_poImgDark)
    {
      cout << "***   Reading dark image... " << m_sDarkName << endl;
      m_poImgDark = new Cimage(m_sDarkName);
    }
  if (!m_poImgDark->bIsAvailable())
    {
      cout << "ERROR, dark image unavailable!" << endl;
      nStat = -1;
      return (nStat);
    }

  // Create an empty bad pixel array.
  // Get the image size needed from the dark image
  // TODO:  Maybe read in a bad pixel image from disk
  
  if (NULL != m_poImgBadPix)
    {
      delete m_poImgBadPix;
    }

  int nDim0, nDim1;
  m_poImgDark->nGetDimensions(&nDim0, &nDim1);  // Should we use flood?
  m_poImgBadPix = new Cimage (nDim0, nDim1, eImage_I2);

  (void) m_poImgBadPix->nSetNextPixel(0, 0);
  for (i = 0; i < nDim0 * nDim1; i++) 
    m_poImgBadPix->vSetNextPixel((short int) 0);

  int a2nTiles[2] = { 0, 0 };
  int a5nNumBads[5];

  cout << "***   Finding bad pixels in the dark image: " << m_sDarkName << endl;
  nFlagBadPixels(m_poImgDark, m_a2nhorizontal_limits,
                 m_a2nvertical_limits, m_a2nsearch_center,
                 (int)m_fsearch_radius,
                 m_a3nmin_pixel[DARKIMG],
                 m_a3nmax_pixel[DARKIMG],
                 m_a2fpixel_sd[DARKIMG], a2nTiles, a5nNumBads);

  if (NULL == m_poImgFlood)
    {
      cout << "***   Reading flood image... " << m_sFloodName << endl;
      m_poImgFlood = new Cimage(m_sFloodName);
    }
  if (!m_poImgFlood->bIsAvailable())
    {
      cout << "ERROR, flood image unavailable!" << endl;
      nStat = -1;
      return (nStat);
    }

  cout << "***   Finding bad pixels in the flood image: " << m_sFloodName << endl;
  if ( (m_a3nmin_pixel[FLOODIMG] * 3) <  m_a3nmax_pixel[FLOODIMG])
    cout << "\n***WARNING: the minimum allowed pixel value " << m_a3nmin_pixel[FLOODIMG] << " in the flood image"
         << "\n            is less than one-third the maximum allowed value " << m_a3nmax_pixel[FLOODIMG] << ".\n";
  else if ( (m_a3nmin_pixel[FLOODIMG] * 2) <  m_a3nmax_pixel[FLOODIMG])
    cout << "\n***INFO: the minimum allowed pixel value " << m_a3nmin_pixel[FLOODIMG] << " in the flood image"
         << "\n         is less than one-half the maximum allowed value " << m_a3nmax_pixel[FLOODIMG] << ".\n";

  nStat = nFlagBadPixels(m_poImgFlood, m_a2nhorizontal_limits,
                 m_a2nvertical_limits, m_a2nsearch_center,
                 (int)m_fsearch_radius,
                 m_a3nmin_pixel[FLOODIMG],
                 m_a3nmax_pixel[FLOODIMG],
                 m_a2fpixel_sd[FLOODIMG], a2nTiles, a5nNumBads);

  if (0 == nStat)
    {
      m_poImgBadPix->vSetState(eImage_available_state);
    }

  if (3 < m_nVerbose)
    {
      m_poImgBadPix->nWrite("badpix.img");
    }

  //****************** Read Bad Pixel List *************************

  (void) nReadBadPixelList(m_sBadpixName);

  if (3 < m_nVerbose)
    {
      m_poImgBadPix->nWrite("badpix2.img");
    }

  //****************** Dark Subtract Flood *************************

  cout << "\n***   Correct flood image for dark image..." << endl;
  nStat = nSubtract(m_poImgFlood, m_poImgDark, m_fdarksub_scale,
                    m_fdarksub_const, (float)m_a3nmin_pixel[2]);
  if (0 != nStat)
    {
      cout << "ERROR performing FLOOD = FLOOD - DARK!" << endl;
      return (nStat);
    }

  // TODO: Mark bad pixels in FLOOD - DARK (set bad pixels to ?)


  //****************** Calculate Reference Image *************************

  if (NULL != m_poImgNonunf)
    {
      delete m_poImgNonunf;
    }

  if ("REFER" == m_sReferName)
    {
      cout << "\nINFO:  calculating reference image: refer.img.\n";
      nStat = nCalcReference(&m_tCorrectParams,
                             &m_tDistorParams,
                             &m_tMaskParams);

      m_poImgRefer->m_oHeader.nReplaceValue(D_K_SaturatedValue, m_fReferenceMax);
      m_poImgRefer->nWrite("refer.img");
    }
  else
    {
      // Read reference image from disk

      if (NULL != m_poImgRefer)
        {
          delete m_poImgRefer;
          m_poImgRefer = NULL;
        }
      cout << "\nINFO:  reading reference image: " << m_sReferName << endl;
      m_poImgRefer = new Cimage(m_sReferName);
    }
  if (!m_poImgRefer->bIsAvailable())
    {
      cout << "ERROR: Reference image is unavailable!\n" 
           << "       Possible cause: beam center is on a BAD pixel position.\n"
           << flush;
      exit (1);
    }

  //****************** Distort Reference Image *************************

  // Need distor files!

  m_poImgNonunf = new Cimage (*m_poImgRefer);
  
  nStat = nDistortReference(m_poImgRefer, m_poImgNonunf);
  if (0 == nStat)
    {
//      cout << "      Setting nonunf image to available.\n";
      m_poImgNonunf->vSetState(eImage_available_state);
    }

  m_poImgNonunf->m_oHeader.nReplaceValue(D_K_SaturatedValue, m_fReferenceMax);
  m_poImgNonunf->nWrite("referdist.img");
  m_poImgNonunf->m_oHeader.nDelete(D_K_SaturatedValue);

  //****************** Scale Flood to Reference Image *************************

  nStat = nScaleImages();

  // ****************** Write NONUNF File *************************

  if ( (0 == nStat) || (3 < m_nVerbose) )
    {
      m_poImgNonunf->m_oHeader.nReplaceValue("NONUNF_MINUSED", m_a3nmin_pixel[FLOODIMG]);
      nStat = m_poImgNonunf->nWrite(sTransSymbol(m_sNonunfName));
    }

  cout << "...Done with calculating NONUNF file." << endl;
  return (nStat);
}

int
Ccalibrate::nFlagBadPixels(Cimage *poImageIn,
                           const int a2nHorzLim[2],
                           const int a2nVertLim[2],
                           const int a2nCenter[2],
                           const int nRadius,
                           const int nMinVal, const int nMaxVal,
                           const float fSigma,
                           const int a2nTilesIn[2], int a5nNumBads[5])

{
#define QOUT   0
#define QBRITE 1
#define QDARK  2
#define Q2SD   3
  // Flag bad pixels in the input image, by marking the equivalent
  // pixel position in m_poImgBadPix.
  // Bad pixels 
/*
      SUBROUTINE FLGBAD(IMAGE, FLGPIX,
     $     ISIZE1, ISIZE2, IUSE1, IUSE2, HORLIM, VERLIM,
     &     CENTER, RADIUS, MINVAL, MAXVAL, SIGMA, 
     &     TILE, BADS, WORK)
C
C     FLGBAD -- FLaG BAD pixels in a given image.  Bad pixels are determined
C     ======    by the following criteria:
C               1.  Larger than a given RADIUS from a given CENTER
C               2.  Below a given MINVAL
C               3.  Above a given MAXVAL
C               4.  More than SIGMA away from an average value calculated
C                   over TILE(1) by TILE(2) pixels.
C
C  IMAGE is UNSIGNED!!!
C
C     27-Nov-1990         J. W. Pflugrath         Cold Spring Harbor Laboratory
C       Created.
C
C     18-Sep-1992         Marty Stanton
C       Modified to flag out of horizontal/vertical bounds
C       Modified to flag SD differently than below MINVAL, above MAXVAL
C
*/
/*
      INTEGER   QMNPNT
      PARAMETER (QMNPNT = 50)
      INTEGER   ISIZE1, ISIZE2, IUSE1, IUSE2, CENTER(2)
      INTEGER   HORLIM(2), VERLIM(2),
     &          RADIUS, MINVAL, MAXVAL, TILE(2), BADS(5), WORK
      REAL      SIGMA
C
      INTEGER I1, I2, J1, J2, IMNVAL, IMXVAL, NUM, IRAD, IRADSQ,
     &        IR2SQ, I4J1J2, ILOOP, ITEST, JTEST
C***********************************************************************
      REAL    RSUM1, RSUM2, RSUM3, AVG, SD, AVGOLD, SDOLD, TEMP, TEST
C
*/
  int i, nStat;
  int nRadSq;
  int nIR2SQ;
  int nIRAD;
  float fPixIJ;
  float fAvg, fSD;
  float fAvgPrev, fSDPrev, fTest;
  float fSum1, fSum2, fNum;
  int   nLoop;
  short int iMark, iMarkPrev;
  float fMaxVal, fMinVal;
  float fMaxLim, fMinLim;
  int i1, i2, j1, j2;
  int nDim0, nDim1;
  int a2nTiles[2];

  for (i = 0; i < 5; i++)
    a5nNumBads[i] = 0;
  nRadSq  = nRadius * nRadius;
  fMaxVal = (float)nMaxVal;
  fMinVal = (float)nMinVal;
  fMaxLim = (float)9999999999.0;
  fMinLim = (float)-999999999.0;

  poImageIn->nGetDimensions(&nDim0, &nDim1);

  a2nTiles[0] = a2nTilesIn[0];
  a2nTiles[1] = a2nTilesIn[1];
  if (0 == a2nTiles[0])
    a2nTiles[0] = max(32, nDim0 / 16);
  if (0 == a2nTiles[1])
    a2nTiles[1] = max(32, nDim1 / 16);
  cout << "\n***   Finding bad pixels:"
       << "\n      Horz limit:     " << a2nHorzLim[0] << ", " << a2nHorzLim[1]
       << "\n      Vert limit:     " << a2nVertLim[0] << ", " << a2nVertLim[1]
       << "\n      Center, radius: " << a2nCenter[0] << ", " << a2nCenter[1] << ", "
                                 << nRadius
       << "\n      Min, max value: " << fMinVal << ", " << fMaxVal
       << "\n      Tile size:      " << a2nTiles[0] << ", " << a2nTiles[1]
       << "\n      SD limit:       " << fSigma
       << endl;

  for (i2 = 0; i2 < nDim1; i2 = i2 + a2nTiles[1])
    {
      for (i1 = 0; i1 < nDim0; i1 = i1 + a2nTiles[0])
        {
          // Calculate average background and standard deviation for a
          // small area on the image.
          // Skip this section if ASIG <= 0 which indicates that
          // only MINVAL and MAXVAL should be used as the criteria.

          if (0.0 < fSigma)
            {
              fSum1 = fSum2 = fNum = 0.0;
              for (j2 = max(i2, a2nVertLim[0]); 
                   j2 < min(min(nDim1, i2+a2nTiles[1]), a2nVertLim[1]+1); j2++)
                //  DO 500 J2 = MAX(I2,VERLIM(1)),
                //              MIN(IUSE2, I2+TILE(2)-1, VERLIM(2))
                {
                  nIR2SQ = j2 - a2nCenter[1];
                  nIR2SQ = nIR2SQ * nIR2SQ;

                  for (j1 = max(i1, a2nHorzLim[0]);
                       j1 < min(min(nDim0, i1 + a2nTiles[0]), a2nHorzLim[1]+1); j1++)
                    {
                      // Final test of limits:
                      
                      nIRAD = j1 - a2nCenter[0];
                      nIRAD = nIRAD * nIRAD +  nIR2SQ;
                      if (nIRAD <= nRadSq)
                        {
                          // Compute sum and sum of squares, use 
                          // floating point so as not to overflow
                          // integer math.

                          fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(j1, j2);
                          fSum2 += (fPixIJ * fPixIJ);
                          fSum1 += fPixIJ;
                          fNum  += 1.0;
                        }
                    }
                }

              // Now compute average and gaussian statistics sd, but require
              // a minimum number of points in the calculation.

              fSD = 0.0;
              if (50 <= fNum)
                {
                  fAvg = fSum1 / fNum;
                  fSD  = (fSum2 - (fSum1 * fSum1 / fNum ) ) / (fNum - 1.0);

                  // SD may be negative (because of round-off errors) so be safe

                  fSD = sqrtf(fabsf(fSD));
                }

              if (0.0 < fSD)
                {
                  // Now compute average and sd, but exclude those points that
                  // are more than 3 sigma from the previous calculated avg and sd.
                  // At the same time allow for iterations.

                  // Beginning of UNTIL, DO loop.  Test is made at bottom of loop.
                  nLoop = 0;
                  do
                    {
                      nLoop++;
                      fTest    = 3.0 * fSD;
                      fAvgPrev = fAvg;
                      fSDPrev  = fSD;
                      fSum1 = 0.0;
                      fSum2 = 0.0;
                      fNum  = 0.0;
                      for (j2 = max(i2, a2nVertLim[0]); 
                           j2 < min(min(nDim1, i2+a2nTiles[1]), a2nVertLim[1]+1);
                           j2++)
                        {
                          // Preliminary test of limits:
                          
                          nIR2SQ = j2 - a2nCenter[1];
                          nIR2SQ = nIR2SQ * nIR2SQ;
                          for (j1 = max(i1, a2nHorzLim[0]);
                               j1 < min(min(nDim0, i1 + a2nTiles[0]), a2nHorzLim[1]+1); j1++)
                            {

                              // Final test of limits:
                      
                              nIRAD = j1 - a2nCenter[0];
                              nIRAD = nIRAD * nIRAD +  nIR2SQ;
                              if (nIRAD <= nRadSq)
                                {
                                  // Compute sum and sum of squares, use 
                                  // floating point so as not to overflow
                                  // integer math.

                                  fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(j1, j2);
                                  if (fabsf(fPixIJ-fAvg) <= fTest)
                                    {
                                      // Include in new average and sd calculation.

                                      fSum2 += (fPixIJ * fPixIJ);
                                      fSum1 += fPixIJ;
                                      fNum  += 1.0;
                                    }
                                }
                            }
                        }

                      // Now compute average and gaussian statistics sd, 
                      // but require
                      // a minimum number of points in the calculation.

                      fSD = 0.0;
                      if (fNum >= 50)
                        {
                          fAvg = fSum1 / fNum;
                          fSD  = (fSum2 - (fSum1 * fSum1 / fNum ) ) 
                               / (fNum - 1.0);

                          // SD may be negative (because of round-off errors) so be safe

                          fSD = sqrtf(fabsf(fSD));
                        }

                      fTest = fabsf(fAvg - fAvgPrev);

                      // Convergence is attained when the average differs by
                      // less than 1 sd from the previous average or when the
                      // average has been calculate 3 times.

                      } while ( (fTest > fSD) && (nLoop < 3) && (0.0 != fSD) );

                  fMaxLim = fAvg + fSigma * fSD;
                  fMinLim = fAvg - fSigma * fSD;
                }
            }

          // Now go through the array to look for pixels outside the 
          // seach limits OR
          // above the threshold specified by IMXVAL or MAXVAL OR
          // below the threshold specified by IMNVAL or MINVAL.

// 1200        CONTINUE

          for (j2 = i2; j2 < min(nDim1, i2 + a2nTiles[1]); j2++)
            {
              nIR2SQ = j2 - a2nCenter[1];
              nIR2SQ = nIR2SQ * nIR2SQ;
              for (j1 = i1; j1 < min(nDim0, i1 + a2nTiles[0]); j1++)
                {
                  // Test of circle limits:

                  iMark = 0;
                  nIRAD = j1 - a2nCenter[0];
                  nIRAD = nIRAD * nIRAD +  nIR2SQ;
                  if (   (nIRAD > nRadSq)
                      || (j2 < a2nVertLim[0])
                      || (j2 > a2nVertLim[1])
                      || (j1 < a2nHorzLim[0])
                      || (j1 > a2nHorzLim[1]) )
                    {
                      // Out of bounds 

                      a5nNumBads[QOUT]++;
                      iMark = QOUT + 1;
                    }
                  else
                    {
                      fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(j1, j2);
/************************		      
		      if ( (j1 == 3005) && (j2 == 1178) )
			{
			  cout << "3005, 1178 val: " << fPixIJ << endl;
			  cout << "fMaxVal, fMinVal: " << fMaxVal << ", "
                               << fMinVal << endl;
			  cout << "fMaxLim, fMinLim: " << fMaxLim << ", "
                               << fMinLim << endl;
			}
**********************/
                      if (fPixIJ > fMaxVal)
                        {
                          // Pixel too bright, so flag.
                          
                          a5nNumBads[QBRITE]++;
                          iMark = QBRITE + 1;
                        }
                      else if (fPixIJ < fMinVal)
                        {
                          // Pixel is too dark, so flag.

                          a5nNumBads[QDARK]++;
                          iMark = QDARK + 1;
                        }
                      else if (   (fPixIJ > fMaxLim)
                               || (fPixIJ < fMinLim) )
                        {
                          // Pixel out of SD bounds, so flag.
                              
                          a5nNumBads[Q2SD]++;
                          iMark = Q2SD + 1;
                        }
                    }

                  // Set mark in bad pixel array, but preserve previous value

                  if (0 != iMark)
                    {
                      (void)m_poImgBadPix->nGetPixel(j1, j2, &iMarkPrev);
                      (void)m_poImgBadPix->nSetPixel(j1, j2, (short int)(iMark | iMarkPrev));
                    }
                }
            }
        }
    }

  int a5nPercent[5];

  a5nPercent[QOUT]   = 100 * a5nNumBads[QOUT]   / (nDim0*nDim1);
  a5nPercent[QBRITE] = 100 * a5nNumBads[QBRITE] / (nDim0*nDim1);
  a5nPercent[QDARK]  = 100 * a5nNumBads[QDARK]  / (nDim0*nDim1);
  a5nPercent[Q2SD]   = 100 * a5nNumBads[Q2SD]   / (nDim0*nDim1);
  cout << "\n Total number of pixels in the image:       " << setw(10) << nDim0 * nDim1;
  cout << "\n Number of bad pixels found: "
       << "\n                             Out of bounds: " << setw(10) << a5nNumBads[QOUT]
       << "  ( " << setw(3) << a5nPercent[QOUT] << "% )"
       << "\n Too intense (above max value of " << setw(8) << fMaxVal << "): " << setw(10)  << a5nNumBads[QBRITE] 
       << "  ( " << setw(3) << a5nPercent[QBRITE] << "% )"
       << "\n Too dark    (below min value of " << setw(8) << fMinVal << "): " << setw(10)  << a5nNumBads[QDARK]
       << "  ( " << setw(3) << a5nPercent[QDARK] << "% )"
       << "\n Too high a standard dev from local average: " << setw(9)  << a5nNumBads[Q2SD]
       << "  ( " << setw(3) << a5nPercent[Q2SD] << "% )\n"
       << endl;
  if (   (5 < a5nPercent[QBRITE])
      || (5 < a5nPercent[QDARK])
      || (5 < a5nPercent[QOUT])
      || (5 < a5nPercent[Q2SD]) )
    cout << "WARNING, check to see why more than 5% of the pixels were rejected!\n\n" << flush;
  return (0);
}

int
Ccalibrate::nReadBadPixelList(Cstring& rsFilename)
{
  // Read the file rsFilename for a list of bad pixel commands
  // Valid commands are ROW, COLUMN, LINE, RECT, PIXEL

  int     nStat;
  Cstring sLine;
  Cstring asTokens[10];
  int     nNumToks;
  int     i, j, ii, jj;
  short int iMark, iMarkPrev;

  nStat = 1;
  ifstream oIn( sTransSymbol(rsFilename).string());
  if ((oIn.rdbuf()->is_open()))
    {
      getline(oIn, sLine);
      while (oIn && !oIn.eof())
        {
          // Successful getline

          if (3 < m_nVerbose)
            {
              cout << "Read line: " << sLine << endl;
            }
          nNumToks = split(sLine, asTokens, 10, " ,\t");
          if (0 == nNumToks)
            asTokens[0] = "";
          asTokens[0].downcase();
          if (   (0 == nNumToks)
              || (0 == asTokens[0].index('!')) )
            {
              // Skip the commented line or empty line
            }
          else if (   (asTokens[0] == "col")
                   || (asTokens[0] == "column") )
            {
              // column i, start, end, flag
              if (3 < nNumToks)
                {
                  if (4 < nNumToks)
                    iMark = (short int) atoi(asTokens[4].string());
                  else
                    iMark = 1;
                  i = atoi(asTokens[1].string());
                  for (j =  atoi(asTokens[2].string());
                       j <= atoi(asTokens[3].string()); j++)
                    {
                      (void)m_poImgBadPix->nGetPixel(i, j, &iMarkPrev);
                      (void)m_poImgBadPix->nSetPixel(i, j, 
                                              (short int)(iMark | iMarkPrev));
                    }
                }
            }
          else if (   (asTokens[0] == "line")
                   || (asTokens[0] == "row") )
            {
              // row j, start, end | flag
              if (3 < nNumToks)
                {
                  if (4 < nNumToks)
                    iMark = (short int) atoi(asTokens[4].string());
                  else
                    iMark = 1;
                  j = atoi(asTokens[1].string());
                  for (i = atoi(asTokens[2].string()); 
                       i <= atoi(asTokens[3].string()); i++)
                    {
                      (void)m_poImgBadPix->nGetPixel(i, j, &iMarkPrev);
                      (void)m_poImgBadPix->nSetPixel(i, j, 
                                              (short int)(iMark | iMarkPrev));
                    }
                }
            }
          else if (   (asTokens[0] == "rectangle")
                   || (asTokens[0] == "rect") )
            {
              // Rect orig0, orig1, ext0, ext1 | flag

              if (4 < nNumToks)
                {
                  if (5 < nNumToks)
                    iMark = (short int) atoi(asTokens[5].string());
                  else
                    iMark = 1;
                  int ii, jj;
                  jj = atoi(asTokens[2].string());
                  for (j = 0; j < atoi(asTokens[4].string()); j++, jj++)
                    {
                      ii = atoi(asTokens[1].string());
                      for (i = 0; i < atoi(asTokens[3].string()); i++, ii++)
                        {
                          (void)m_poImgBadPix->nGetPixel(ii, jj, &iMarkPrev);
                          (void)m_poImgBadPix->nSetPixel(ii, jj, 
                                              (short int)(iMark | iMarkPrev));
                        }
                    }
                }

            }
          else if (asTokens[0] == "insidecircle")
            {
              // Circle center0, center1, radius0 | flag

              if (3 < nNumToks)
                {
                  if ( 4 < nNumToks)
                    iMark = (short int) atoi(asTokens[4].string());
                  else
                    iMark = 1;
                  int ii, jj;
		  int icen, jcen, iradsq;
		  icen   = atoi(asTokens[1].string());
		  jcen   = atoi(asTokens[2].string());
		  iradsq = atoi(asTokens[3].string());
		  cout << "center, radius: " << icen << ", " << jcen << ", " 
                       << iradsq << endl;

		  iradsq = iradsq * iradsq;
		  int nDim0, nDim1;
		  (void)m_poImgBadPix->nGetDimensions(&nDim0, &nDim1);
                  for (j = 0; j < nDim1; j++)
                    {
                      jj = (j-jcen) * (j-jcen);
                      for (i = 0; i < nDim0; i++)
                        {
			  ii = (i-icen) * (i-icen);
			  if (iradsq > (ii + jj) )
			    {
			      (void)m_poImgBadPix->nGetPixel(i, j, &iMarkPrev);
			      (void)m_poImgBadPix->nSetPixel(i, j, 
                                              (short int)(iMark | iMarkPrev));
			    }
                        }
                    }
                }
            }
          else if (asTokens[0] == "outsidecircle")
	    {
              // Circle center0, center1, radius0 | flag

              if (3 < nNumToks)
                {
                  if ( 4 < nNumToks)
                    iMark = (short int) atoi(asTokens[4].string());
                  else
                    iMark = 1;
                  int ii, jj;
		  int icen, jcen, iradsq;
		  icen   = atoi(asTokens[1].string());
		  jcen   = atoi(asTokens[2].string());
		  iradsq = atoi(asTokens[3].string());
		  cout << "center, radius: " << icen << ", " << jcen << ", " 
                       << iradsq << endl;

		  iradsq = iradsq * iradsq;

		  int nDim0, nDim1;
		  (void)m_poImgBadPix->nGetDimensions(&nDim0, &nDim1);
                  for (j = 0; j < nDim1; j++)
                    {
                      jj = (j-jcen) * (j-jcen);
                      for (i = 0; i < nDim0; i++)
                        {
			  ii = (i-icen) * (i-icen);
			  if (iradsq < (ii + jj) )
			    {
			      (void)m_poImgBadPix->nGetPixel(i, j, &iMarkPrev);
			      (void)m_poImgBadPix->nSetPixel(i, j, 
                                              (short int)(iMark | iMarkPrev));
			    }
                        }
                    }
                }
            }
          else if (asTokens[0] == "pixel")
            {
              // pixel i j | flag
              if (2 < nNumToks)
                {
                  if (3 < nNumToks)
                    iMark = (short int) atoi(asTokens[3].string());
                  else
                    iMark = 1;
                  i = atoi(asTokens[1].string());
                  j = atoi(asTokens[2].string());
                  (void)m_poImgBadPix->nGetPixel(i, j, &iMarkPrev);
                  (void)m_poImgBadPix->nSetPixel(i, j, (short int)(iMark | iMarkPrev));
                }
            }
          else
            {
              cout << "WARNING: unknown bad pixel keyword: " << asTokens[0] << endl;
            }
          getline(oIn, sLine);
        }

      if (oIn.eof()) 
        {
          cout << "INFO end-of-file on " << sTransSymbol(rsFilename) << endl;
          nStat = 0;
        }
      else if (!oIn) 
        {
          cout << "ERROR reading bad pixel file!\n" << endl;
          nStat = 1;
        }
      else
        {
          cout << "Duh? reading bad pixel file!\n" << endl;
          nStat = 2;
        }
      }
  oIn.close();                 // Make sure we always close the stream!

//+27-Aug-2001
  // Kludge to get bad pixels from another image.  Any pixel greater than
  // a certain value in the image will be a bad pixel

  sLine = sGetEnv("JIM_BADPIX");
  if ("" != sLine)
    {
      Cimage *poImage = NULL;
      nNumToks = split(sLine, asTokens, 10, " ,\t");
      poImage = new Cimage(asTokens[0]);
      if (!poImage->bIsAvailable())
	{
	  cout << "ERROR: " << asTokens[0] << " badpix image is unavailable!\n" 
               << flush;
	  delete poImage;
	  poImage = NULL;
	}
      else
	{
	  int nDim0, nDim1;
      short int iMaxPix;
	  //short int iMinPix, 
	  m_poImgBadPix->nGetDimensions(&nDim0, &nDim1);
	  poImage->nSetNextPixel(0,0);
	  m_poImgBadPix->nSetNextPixel(0,0);
	  iMaxPix = (short int) atoi(asTokens[1]);
//	  iMinPix = (short int) atoi(asTokens[2]);
	  cout << "\n     INFO: pixels with values >= " << iMaxPix 
//               << "\n              and with values <= " << iMinPix
               << "\n              in the image " << asTokens[0] 
               << "\n              will be flagged as bad.\n" << flush;
	  ii = 0;
	  for (j = 0; j < nDim1; j++)
	    {
	      for (i = 0; i < nDim0; i++)
		{
		  // This could be done faster if we know the data type of
		  // the arrays.  The JIM_BADPIX image is supposed to be a
		  // NONUNF file, 

		  iMark = poImage->iGetNextPixel();
		  if (iMark >= iMaxPix)
		    {
		      (void)m_poImgBadPix->nSetPixel(i, j, (short int)(5));
		      ii++;
		    }
		  //else if (iMark <= iMinPix)
		  //  {
		  //  }
		}
	    }
	  cout << "     There were " << ii << " pixels with values >= " 
               << iMaxPix
	       << "\n        in the image file " << asTokens[0] << endl << flush;
	}
    }

//-27-Aug-2001
  return (nStat);
}

int
Ccalibrate::nCalcReference(tagCorrect_parameters *ptCorrectParams,
                           tagDistort_parameters *ptDistorParams,
                           tagMask_parameters    *ptMaskParams)
{
  int nStat;

  if (NULL != m_poImgRefer)
    delete m_poImgRefer;

  int nDim0, nDim1;

  cout << "\n***   Calculate reference (expected flood) image from geometry..." 
       << endl;

  m_poImgFlood->nGetDimensions(&nDim0, &nDim1);
//+JWP 2007-10-12 use unsigned short for the reference image 
//                DO NOT inherit data_type from flood image
  m_poImgRefer = new Cimage (nDim0, nDim1, eImage_uI2);

  if (   (NULL == m_poImgDistorXINT) || (!m_poImgDistorXINT->bIsAvailable())
      || (NULL == m_poImgDistorYINT) || (!m_poImgDistorYINT->bIsAvailable()) )
    {
      nStat = nReadDistor();
      if (0 != nStat)
        {
          cout << "ERROR Reading spatial distortion info!" << endl;
          return (nStat);
        }
    }
      
  float input_xbeam, input_ybeam;
  float ctod, ctod2;
  float pixsiz, pixsiz2;
  int   i, j;
  float r2;

  // Calculate where beam would be on undistorted image

//  nStat = nPxToPx(ptDistorParams->x_beam, ptDistorParams->y_beam,
//  nStat = nPxToPx(ptMaskParams->beam[0], ptMaskParams->beam[1],
  nStat = nPxToPx(m_a2fbeam_position[0], m_a2fbeam_position[1],
                  ptCorrectParams,
                  &input_xbeam, &input_ybeam);

//  ctod    = ptMaskParams->source;
  ctod    = m_fxtod_distance;
  ctod2   = ctod * ctod;
//  pixsiz  = ptDistorParams->spacing / ptDistorParams->x_scale;
  pixsiz  = m_fpixel_size;
  pixsiz2 = pixsiz * pixsiz;

  if ( (3 < m_nVerbose) || (0 != nStat) )
    {
      // Print out this info if verbosity is high enough or there is an error

      cout << "\nSource to Detector Distance (mm) : " << ctod
	   << "\n(Info) Mask to Det Distance (mm) : " << m_fmasktod_distance
           << "\nPixel Size (mm) :                  " << pixsiz

           << "\nBeam Position raw:                 " << m_a2fbeam_position[0]
           << ", " <<  m_a2fbeam_position[1]
//           << "\nBeam Position raw:                 " << ptDistorParams->x_beam
//           << ", " << ptDistorParams->y_beam 
           << "\nBeam Position undistorted image:   " << input_xbeam
           << ", " << input_ybeam
           << endl;
    }

  // This makes a cosine fall-off reference image.
  // Pixel value at shortest distance to detector is normalized to 10000.
  // TODO:  What about difference pixel sizes in the two directions?

  float fI2, fJ2;
  short int iRef;
  m_poImgRefer->nSetNextPixel(0,0);
  for (j = 0; j < nDim1; j++)
    {
      fJ2 = (float)j - input_ybeam;
      fJ2 = fJ2 * fJ2;
      for (i = 0; i < nDim0; i++)
        {
          fI2 = (float)i - input_xbeam;
          fI2 = fI2 * fI2;
          r2  = (fI2 + fJ2) * pixsiz2;
          r2  = m_fReferenceMax * ( ctod2 / (r2 + ctod2) );
//          iRef = (short int) nint(r2);
          iRef = (short int) r2;
          m_poImgRefer->vSetNextPixel(iRef);
        }
    }
  if (0 == nStat)
    {
      // This nStat comes from the return status of the above ::nPxToPx(...)
      // call that determines the direct beam position in pixels on the 
      // undistorted image.  If the direct beam position cannot be determined
      // because the spatial distortion info will not allow it (i.e. beam is on
      // a bad pixel, then do NOT change make the image "available".

      m_poImgRefer->vSetState(eImage_available_state);
      
    }
  else
    {
      cout << "WARNING: Cannot convert direct beam position from\n"
           << "         uncorrected to corrected pixel positions.\n"
           << "         Reference image will be unavailable.\n" << flush;
    }
    
  return (nStat);
}

int
Ccalibrate::nPxToPx(const float xi, const float yi,
                    tagCorrect_parameters *ptCorrectParams,
                    float *pfXc, float *pfYc)
{
  int nStat;
  
  int         ixlook, iylook;
  float       xlook, ylook;
  float       dx, dy, xc, yc;

  nStat = 0;

  int nDim0, nDim1;
  m_poImgDistorXINT->nGetDimensions(&nDim0, &nDim1);

  // Calculate position in lookup table
  // We need to add 1 to each xi, yi to account for differences
  //       between C++ and Fortran array origins (?)

  xlook  = (xi - ptCorrectParams->xint_start)
           / float(ptCorrectParams->xint_step) + 1;
  ylook  = (yi - ptCorrectParams->yint_start)
           / float(ptCorrectParams->yint_step) + 1;
  ixlook = int(xlook);
  iylook = int(ylook);
  dx     = xlook - float(ixlook);
  dy     = ylook - float(iylook);

  // Check if all the interpolation points are in bounds.  If
  // not set nStat to -1 and return
  // If all interpolation points are inbounds, make sure that
  // all positions contain values

  if (   (ixlook < 1) || (ixlook > nDim0)
      || (iylook < 1) || (iylook > nDim1) )
    {
      nStat = -1;
      *pfXc = 0.0;
      *pfYc = 0.0;
    }
  else
    {
      LONG    lPix[8];  
      LONG    lIntBadFlag;

      lIntBadFlag = (long) ptCorrectParams->badflag;

      lPix[0] = lIntBadFlag;
      lPix[1] = lIntBadFlag;
      lPix[2] = lIntBadFlag;
      lPix[3] = lIntBadFlag;
      lPix[4] = lIntBadFlag;
      lPix[5] = lIntBadFlag;
      lPix[6] = lIntBadFlag;
      lPix[7] = lIntBadFlag;

      // All corners must have values in them

      m_poImgDistorXINT->nGetPixel(ixlook-1, iylook-1, &lPix[0]);
      m_poImgDistorXINT->nGetPixel(ixlook,   iylook-1, &lPix[1]);
      m_poImgDistorXINT->nGetPixel(ixlook-1, iylook,   &lPix[2]);
      m_poImgDistorXINT->nGetPixel(ixlook,   iylook,   &lPix[3]);
      m_poImgDistorYINT->nGetPixel(ixlook-1, iylook-1, &lPix[4]);
      m_poImgDistorYINT->nGetPixel(ixlook,   iylook-1, &lPix[5]);
      m_poImgDistorYINT->nGetPixel(ixlook-1, iylook,   &lPix[6]);
      m_poImgDistorYINT->nGetPixel(ixlook,   iylook,   &lPix[7]);
/*    else if ( x_int( ixlook,iylook     ) .eq. INTFLAG .or.
     $          x_int( ixlook+1,iylook   ) .eq. INTFLAG .or.
     $          x_int( ixlook,iylook+1   ) .eq. INTFLAG .or.
     $          x_int( ixlook+1,iylook+1 ) .eq. INTFLAG .or.
     $          y_int( ixlook,iylook     ) .eq. INTFLAG .or.
     $          y_int( ixlook+1,iylook   ) .eq. INTFLAG .or.
     $          y_int( ixlook,iylook+1   ) .eq. INTFLAG .or.
     $          y_int( ixlook+1,iylook+1 ) .eq. INTFLAG  )    then
         istat = -1
C         write (6,*) 'Bad Value in table'
         return
      end if
*/
      if (   (lIntBadFlag == lPix[0])
          || (lIntBadFlag == lPix[1])
          || (lIntBadFlag == lPix[2])
          || (lIntBadFlag == lPix[3])
          || (lIntBadFlag == lPix[4])
          || (lIntBadFlag == lPix[5])
          || (lIntBadFlag == lPix[6])
          || (lIntBadFlag == lPix[7]) )
        {
          nStat = -1;
          *pfXc = 0.0;
          *pfYc = 0.0;
        }
      else
        {
          // Everything seems OK - do the interpolation

          xc = lPix[0]
               + dx * ( lPix[1] - lPix[0] )
               + dy * ( lPix[2] - lPix[0] )
               + dx * dy * ( lPix[0] + lPix[3] - lPix[1] - lPix[2] );

          yc = lPix[4]
               + dx * ( lPix[5] - lPix[4] )
               + dy * ( lPix[6] - lPix[4] )
               + dx * dy * ( lPix[4] + lPix[7] - lPix[5] - lPix[6] );

          // The interpolation tables contain integers which must be
          // multiplied by 'pscale' to get position

          *pfXc = xc * ptCorrectParams->pscale;
          *pfYc = yc * ptCorrectParams->pscale;
          nStat = 0;
        }
    }
  return (nStat);
}
int
Ccalibrate::nDistortReference(Cimage *poImageIn, Cimage *poImageOut)
{
  // Distort reference image - after distorting, distorted
  // reference is stored in BOTH the REFERENCE and NONUNF images
  // (The above statement is NOT TRUE!)
  // In reality, this routine computes the ratio of the pixel area in the
  // raw image to the area in the undistorted image.  That is, it
  // normalizes the reference image based on the actual pixel area in the
  // raw (undistorted image).

  float x_diff, y_diff;
  int   cortyp;            /*1:radial 2:inter*/
  int   i;
  int   nStat1;
  int   nStat;
  float fPixIJ;
  int   ii, jj;
  float size;

  x_diff = 0.0; y_diff = 0.0;

  if (   (NULL == m_poImgDistorXINT) || (!m_poImgDistorXINT->bIsAvailable())
         || (NULL == m_poImgDistorYINT) || (!m_poImgDistorYINT->bIsAvailable()) )
    {
      nStat = nReadDistor();
      if (0 != nStat)
        {
          cout << "ERROR Reading spatial distortion info!" << endl;
          return (nStat);
        }
    }

  cout << "\n***   Modify reference image for spatial distortion effects..." 
       << endl;

//  cout << "xstep in nDistorRef: " << m_tCorrectParams.xint_step << endl;
//  cout << "ystep in nDistorRef: " << m_tCorrectParams.yint_step << endl;
  cortyp = 2;
  if ( m_nflood_radial == 1 || m_nflood_interpolate == 1 ) 
    {
      if ( m_nflood_interpolate == 1 )
        cortyp = 2;
      else if ( m_nflood_radial == 1 )
        cortyp = 1;
    }

  if (0 == cortyp)
    {
      // Do nothing, reference image does not need to be distorted
    }
  else
    {
/*
      subroutine mkflat ( 
     $     input, output, 
     $     axsize, aysize, ixsize, iysize,
     $     x_diff, y_diff, cortyp,
     $     dparams, mparams, cparams,
     $     x_int, y_int,
     $     iaxsize, iaysize, iixsize, iiysize )
      
      
      integer   axsize, aysize, ixsize, iysize
      integer*2 input(axsize, aysize), output(axsize, aysize)
      real      x_diff, y_diff
      integer   cortyp
      record    /distort_parameters/ dparams
      record    /mask_parameters/ mparams
      record    /correct_parameters/ cparams
      integer   iaxsize, iaysize, iixsize, iiysize
      integer*4 x_int(iaxsize, iaysize), y_int(iaxsize, iaysize)

      integer*4   long
      integer     i, j, ii, jj

      real    pxl,pyl,  pxr,pyr, pxt,pyt
      real    pxb(3072), pyb(3072)
      real    pxc, pyc
      real    size
*/
      int j;
      int nDim0, nDim1;
      int nNumTooBig;
      float *pfXB;
      float *pfYB;
      float pxl, pyl;
      float pxr, pyr, pxt, pyt;

      (void) poImageIn->nGetDimensions(&nDim0, &nDim1);
      pfXB = new float [max(nDim0,nDim1)];
      pfYB = new float [max(nDim0,nDim1)];

      // For each position on the distorted detector image, calculate the
      // size of the projection of that pixel onto the corrected image.

      // Set the whole output image to all 1s.

      fPixIJ = 1.0;
      for (j = 0; j < nDim1; j++)
        {
          for (i = 0; i < nDim0; i++)
            {
              (poImageOut->*poImageOut->prnSetPixel)(i, j, fPixIJ);
            }
        }

      nNumTooBig = 0;
      cortyp     = 2;
      j          = 2;
      
      for (i = 1; i <= nDim0-1; i++)  // Watch limits?
        {
          if (1 == cortyp)
            {
              nStat1 = nRadialPxToPx((float)i, (float)j-0.5, &pfXB[i], &pfYB[i]);
            }
          else if (2 == cortyp)
            {
              nStat1 = nPxToPx((float)i, (float)j-0.5, &m_tCorrectParams,
                              &pfXB[i], &pfYB[i]);
            }
        }
         
      for (j = 1; j <= nDim1-1; j++)  // Double check limits
        {
          if (1 == cortyp)
            {
              nStat = nRadialPxToPx(2.0-0.5, (float)j, &pxl, &pyl);
            }
          else if (2 == cortyp)
            {
              nStat = nPxToPx(2.0-0.5, (float)j, &m_tCorrectParams, &pxl, &pyl);
            }

          // do i = 2, ixsize-1
          for (i = 1; i <= nDim0-1; i++) // Watch limits!
            {
              if (1 == cortyp)
                {
                  nStat = nRadialPxToPx((float)i + 0.5, (float)j, &pxr, &pyr);
                  if (0 == nStat)
                    nStat = nRadialPxToPx((float)i, (float)j + 0.5, &pxt, &pyt);
                }
              else if (2 == cortyp)
                {
                  nStat = nPxToPx((float)i+0.5, (float)j, &m_tCorrectParams, 
                                  &pxr, &pyr);
                  if (0 == nStat)
                    nStat = nPxToPx((float)i, (float)j+0.5, &m_tCorrectParams, 
                                    &pxt, &pyt);
                }
              fPixIJ = 1.0;
              if (   (pxl > 0.0) && (pyl > 0.0)
                  && (pxr > 0.0) && (pyr > 0.0)
                  && (pxt > 0.0) && (pyt > 0.0)
                  && (pfXB[i] > 0.0) && (pfYB[i] > 0.0)
                  && (0 == nStat) )
                
                {
                  // TODO:  We should have a lower limit for the pixel size,
                  // that is, output pixels that are less than say 1/10th the
                  // area of the input pixel should be flagged as bad?

                  size = (pxr - pxl) * (pyt - pfYB[i]);
/*
                  if (0.0 >= size)
                    {
                      cout << "WARNING, pixel size <= 0.0 for pixel i,j: " << i << ", " 
                           << j << endl;
                    }
*/
                  ii = nint((pxr+pxl)/2.0f + x_diff);
                  jj = nint((pfYB[i]+pyt)/2.0f + y_diff);
                  if (   (size > 0.001) 
                      && (ii >=  1) && (ii <= nDim0)
                      && (jj >=  1) && (jj <= nDim1) )
                    {
                      fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(ii, jj)
                               * size;
                      // Debug, put the size in the pixel value to look at it
//                      fPixIJ = 100.0
//                               * size;
                      // Mark as bad pixels which distort to give counts
                      // more than m_fReferenceDistortMaxFactor
                      // times the maximum in the input reference image

                      if (  (m_fReferenceDistortMaxFactor * m_fReferenceMax) 
                          < fPixIJ)
                        {
                          fPixIJ = 1.0;
                          nNumTooBig++;
                        }
                    }
                }
              (poImageOut->*poImageOut->prnSetPixel)(i, j, fPixIJ);
              pfXB[i] = pxt;
              pfYB[i] = pyt;
              pxl     = pxr;
              pyl     = pyr;
            }
        }

      // Marty Stanton's original code sets a border of 5 pixels 
      // to 0
/*
      do j = 1, 5
         do i = 1, ixsize
            output(i,j) = 0
         end do
      end do
      do j = iysize-4, iysize
         do i = 1, ixsize
            output(i,j) = 0
         end do
      end do
      do j = 1, iysize
         do i = 1, 5
            output(i,j) = 0
         end do
      end do
      do j = 1, iysize
         do i = ixsize-4, ixsize
            output(i,j) = 0
         end do
      end do
*/
      delete [] pfXB;
      delete [] pfYB;

      cout << "Number of distorted pixels that are too big: " 
           <<  nNumTooBig
           << endl;
    }
  if (5 < m_nVerbose)
    cout << "... returning from nDistortReference!" << endl;
  return (0);
}

int
Ccalibrate::nDistortReference2(Cimage *poImageIn, Cimage *poImageOut)
{
  // Distort reference image - after distorting, distorted
  // reference is stored in BOTH the REFERENCE and NONUNF images
  // (The above statement is NOT TRUE!)
  // In reality, this routine computes the ratio of the pixel area in the
  // raw image to the area in the undistorted image.  That is, it
  // normalizes the reference image based on the actual pixel area in the
  // raw (undistorted image).

  float x_diff, y_diff;
  int   cortyp;            /*1:radial 2:inter*/
  int   i;
  int   nStat1;
  int   nStat;
  float fPixIJ;
  int   ii, jj;
  float size;
  float fSizePrev = 0.0;

  x_diff = 0.0; y_diff = 0.0;

  cout << "\n***   II. Modify reference image for spatial distortion effects..." 
       << endl;

//  cout << "xstep in nDistorRef: " << m_tCorrectParams.xint_step << endl;
//  cout << "ystep in nDistorRef: " << m_tCorrectParams.yint_step << endl;
  cortyp = 2;
  if ( m_nflood_radial == 1 || m_nflood_interpolate == 1 ) 
    {
      if ( m_nflood_interpolate == 1 )
        cortyp = 2;
      else if ( m_nflood_radial == 1 )
        cortyp = 1;
    }

  if (0 == cortyp)
    {
      // Do nothing, reference image does not need to be distorted
    }
  else
    {
      int j;
      int nDim0, nDim1;
      int nNumTooBig;
      float *pfXB;
      float *pfYB;
      float pxl, pyl;
      float pxr, pyr, pxt, pyt;

      (void) poImageIn->nGetDimensions(&nDim0, &nDim1);
      pfXB = new float [max(nDim0,nDim1)];
      pfYB = new float [max(nDim0,nDim1)];

      // For SOME positions on the distorted detector image, calculate the
      // size of the projection of that pixel onto the corrected image.

      fPixIJ = 1.0;

      nNumTooBig = 0;
      cortyp     = 2;

      int a2nStart[2], a2nEnd[2], a2nDir[2];
      int a2nMin[2], a2nMax[2];
      
// Modify the algorithm to start near a point (center?)
// and work to the four corners.

      int nQuadrant;
      for (nQuadrant = 1; nQuadrant < 2; nQuadrant++)  // Start with just 1 quadrant
	{
	  if (1 == nQuadrant)
	    {
	      a2nStart[0] = 3150; // Starting pixel in the output image
	      a2nStart[1] = 1050; // Starting pixel in the output image
	      a2nDir[0]   = -1;     // -1 or 1
	      a2nDir[1]   =  1;     // -1 or 1
	      a2nMin[0]   = 2088;   //
	      a2nMax[0]   = a2nStart[0];
	      a2nMin[1]   = a2nStart[1];
	      a2nMax[1]   = 2088;
	      cout << "INFO: " << a2nMin[0] << " <= i <= " << a2nMax[0] << endl;
	      cout << "INFO: " << a2nMin[1] << " <= j <= " << a2nMax[1] << endl;
	      cout << "INFO: Start i, j: " << a2nStart[0] << ", " << a2nStart[1] << "\n"
                   << "INFO: Step  i, j: " << a2nDir[0] << ", " << a2nDir[1] << endl;
	    }
	  else
	    {
	      // Other corners are not yet supported

	      return (0);
	    }
	  bool bInBounds;

	  i = a2nStart[0];
	  j = a2nStart[1];
	  bInBounds = ( (a2nMin[0] <= i) && (a2nMax[0] >= i) && (a2nMin[1] <= j) && (a2nMax[1] >= j) );
	  while (bInBounds) //for (i = nJWPiLow; i <= nDim0-1; i++)  // Watch limits?
	    {
	      //cout << "i, j: " << i << ", " << j << endl;
	      if (1 == cortyp)
		{
		  nStat1 = nRadialPxToPx((float)i, (float)j-0.5f, &pfXB[i], &pfYB[i]);
		}
	      else if (2 == cortyp)
		{
		  nStat1 = nPxToPx((float)i, (float)j-0.5f, &m_tCorrectParams,
				   &pfXB[i], &pfYB[i]);
		  if (0 != nStat)
		    {
		      //cout << "Non-zero nStat when loading edge in nDistorReference for index i = " << i << " !\n";
		      pfXB[i] = 0.0;
		      pfYB[i] = 0.0;
		    }
		}
	      i = i + a2nDir[0];
	      bInBounds = ( (a2nMin[0] <= i) && (a2nMax[0] >= i) && (a2nMin[1] <= j) && (a2nMax[1] >= j) );
	    }

	  i = a2nStart[0];
	  j = a2nStart[1];
	  //cout << "i, j: " << i << ", " << j << endl;
	  bInBounds = ( (a2nMin[0] <= i) && (a2nMax[0] >= i) && (a2nMin[1] <= j) && (a2nMax[1] >= j) );
	  while (bInBounds) //for (j = nJWPjLow; j <= nDim1-1; j++)  // Double check limits
	    {
	      // For each j

	      if (1 == cortyp)
		{
		  nStat = nRadialPxToPx((float)i-0.5f, (float)j, &pxl, &pyl);
		}
	      else if (2 == cortyp)
		{
		  nStat = nPxToPx((float)i-0.5f, (float)j, &m_tCorrectParams, &pxl, &pyl);
		}
	      i = a2nStart[0]; // Reset i
	      while (bInBounds) // for (i = nJWPiLow; i <= nDim0-1; i++) // Watch limits!
		{
		  // For each i

		  //cout << "i, j: " << i << ", " << j << endl;
		  if (1 == cortyp)
		    {
		      nStat = nRadialPxToPx((float)i + 0.5f, (float)j, &pxr, &pyr);
		      if (0 == nStat)
			nStat = nRadialPxToPx((float)i, (float)j + 0.5f, &pxt, &pyt);
		    }
		  else if (2 == cortyp)
		    {
		      nStat = nPxToPx((float)i+0.5f, (float)j, &m_tCorrectParams, 
				      &pxr, &pyr);
		      if (0 == nStat)
			nStat = nPxToPx((float)i, (float)j+0.5f, &m_tCorrectParams, 
					&pxt, &pyt);
		    }
		  fPixIJ = 1.0;
		  if (   (pxl > 0.0) && (pyl > 0.0)
		      && (pxr > 0.0) && (pyr > 0.0)
		      && (pxt > 0.0) && (pyt > 0.0)
		      && (pfXB[i] > 0.0) && (pfYB[i] > 0.0)
		      && (0 == nStat) )
		    {
		      // TODO:  We should have a lower limit for the pixel size,
		      // that is, output pixels that are less than say 1/10th the
		      // area of the input pixel should be flagged as bad?

		      size = (pxr - pxl) * (pyt - pfYB[i]);
		      size = size * (float)a2nDir[0] * (float)a2nDir[1];
		      if (0.0 >= size)
			{
			  cout << "WARNING, pixel size " << size << " <= 0.0 for pixel i,j: " << i << ", " 
			       << j << "   prevsize: " << fSizePrev << endl;
			  cout << "pxr - pxl  * pyt - pfYB[i] " 
			       << pxr << ", "
			       << pxl << ", "
			       << pyt << ", "  << pfYB[i] << endl;
			  // Better to use previous nearby size than a bogus size
			  size = -size;
			}
		      else
			fSizePrev = size;

		      ii = nint((pxr+pxl)/2.0f + x_diff);
		      jj = nint((pfYB[i]+pyt)/2.0f + y_diff);
		      if (   (size > 0.001) 
			     //&& (ii >=  1) && (ii <= nDim0)
			     //&& (jj >=  1) && (jj <= nDim1) )
			     && (ii >=  0) && (ii < nDim0)
			     && (jj >=  0) && (jj < nDim1) )

			{
			  fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(ii, jj)
			    * size;
			  // Debug, put the size in the pixel value to look at it
			  fPixIJ = 1000.0
			    * size;
			  // Mark as bad pixels which distort to give counts
			  // more than m_fReferenceDistortMaxFactor
			  // times the maximum in the input reference image

			  if (  (m_fReferenceDistortMaxFactor * m_fReferenceMax) 
				< fPixIJ)
			    {
			      fPixIJ = 1.0;
			      nNumTooBig++;
			    }
			}
		    }
		  (poImageOut->*poImageOut->prnSetPixel)(i, j, fPixIJ);
		  pfXB[i] = pxt;
		  pfYB[i] = pyt;
		  pxl     = pxr;
		  pyl     = pyr;

		  i = i + a2nDir[0];
		  //cout << "A: i, j: " << i << ", " << j << endl;
		  // Check only i
		  bInBounds = ( (a2nMin[0] <= i) && (a2nMax[0] >= i) );
		}
	      j = j + a2nDir[1];
	      //cout << "B: i, j: " << i << ", " << j << endl;
	      //cout << "INFO: " << a2nMin[0] << " <= i <= " << a2nMax[0] << endl;
	      //cout << "INFO: " << a2nMin[1] << " <= j <= " << a2nMax[1] << endl;
	      //cout << "INFO: Start i, j: " << a2nStart[0] << ", " << a2nStart[1] << "\n"
	      // << "INFO: Step  i, j: " << a2nDir[0] << ", " << a2nDir[1] << endl;
	      // Check only j
	      bInBounds = ( (a2nMin[1] <= j) && (a2nMax[1] >= j) );
	      // if (bInBounds) cout << "bInBounds is TRUE\n" << endl;
	    }

	} // End of different corner

      delete [] pfXB;
      delete [] pfYB;

      cout << "Number of distorted pixels that are too big: " 
           <<  nNumTooBig
           << endl;
    }
  if (5 < m_nVerbose)
    cout << "... returning from nDistortReference2!" << endl;
  return (0);
}

int
Ccalibrate::nRadialPxToPx(const float xi, const float yi,
                          float *pfXo, float *pfYo)
{
  float x, y, r;

  if (5 < m_nVerbose)
    cout << "nRadialPxToPx called!\n";
  y    = yi - m_tDistorParams.y_center;
  x    = xi - m_tDistorParams.x_center;
  r    = x * x  +  y * y;
  r    = m_tDistorParams.a * m_tDistorParams.a1
       + m_tDistorParams.b * sqrtf(r) 
       + m_tDistorParams.c * r;
  *pfYo = y * r / m_tDistorParams.ratio + m_tDistorParams.y_center;
  *pfXo = x * r + m_tDistorParams.x_center;
  return (0);
}

int
Ccalibrate::nScaleImages(void)
{
  // Scale the distorted reference image found in m_poImgNonunf to
  // the flood field image found in m_poImgFlood.  Bad pixels are found
  // in image m_poImgBadpix
/*
C
C
C     Scale reference to flood image to get nonunf
C
C     On input, nonunf contains the reference image
C
      subroutine  scale (nonunf, flood, flgpix,
     $     xsize, ysize, xpix, ypix, flags, work )

      implicit        none

      include         'caldefs'

      integer         xpix, ypix, xsize, ysize
      integer*2       flood (xsize, ysize)
      integer*2       nonunf (xsize, ysize)
      byte            flgpix (xsize, ysize)
      integer         flags(6)
      real            min_ratio, max_ratio, r, f, m

      integer         i, j, n
      integer         iyoff, ixoff
      integer         work
      real            sc
      integer         icount
*/
  int   nStat;
  float min_ratio, max_ratio;
  float fPixFlood;
  float fPixRefer;
  float sc, r;

  nStat     = 0;
  max_ratio = 0.0;
  min_ratio = 10000000.0;

  cout << "\n***   Scale flood image to modified reference image..." << endl;

  if ( (NULL == m_poImgNonunf) || (!m_poImgNonunf->bIsAvailable()) )
    {
      cout << "ERROR in Ccalibrate::nScaleImages: Nonunf image not available!"
           << endl;
      return (-1);
    }
  if ( (NULL == m_poImgFlood) || (!m_poImgFlood->bIsAvailable()) )
    {
      cout << "ERROR in Ccalibrate::nScaleImages: Flood image not available!"
           << endl;
      return (-1);
    }

  if ( (NULL == m_poImgBadPix) || (!m_poImgBadPix->bIsAvailable()) )
    {
      cout << "ERROR in Ccalibrate::nScaleImages: BadPix image not available!"
           << endl;
      return (-1);
    }

  int i, j;
  int nDim0, nDim1;
  short int iBadPixFlag;
  int nNewBadCount = 0;
  int nMinMaxR;

  m_poImgFlood->nGetDimensions(&nDim0, &nDim1);
  m_poImgNonunf->nSetNextPixel(0,0);
  m_poImgFlood->nSetNextPixel(0,0);
  m_poImgBadPix->nSetNextPixel(0,0);
  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
        {
          // This could be done faster if we know the data type of
          // the arrays

          fPixFlood   = (m_poImgFlood->*m_poImgFlood->prfGetPixel)(i, j);
          fPixRefer   = (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(i, j);
          iBadPixFlag = m_poImgBadPix->iGetNextPixel();
          if (   (0 == iBadPixFlag) 
              && (0.0 != fPixFlood)
              && (2.0 < fPixRefer) ) // Maybe we should flag bad pixels in refer
            {
               r = fPixRefer / fPixFlood;
	       nMinMaxR = 0;
               if ( r > max_ratio ) 
		 {
		   // Save the max_ratio found

		   nMinMaxR  = 1;
		   max_ratio = r;
		 }
               if ( r < min_ratio ) 
		 {
		   // Save the min_ratio found, but only down to 2 (see above if)

		   nMinMaxR  = -1;
		   min_ratio = r;
		 }
	       if (0 != nMinMaxR)
		 cout << "INFO max or min found at " << i << ", " << j << ": "
                      << r << endl;
	    }
          if ( (1.0 >= fPixRefer) && (0 == iBadPixFlag) )
            {
	      // Set this as a bad pixel, but why?

              m_poImgBadPix->nSetPixel(i, j, (short int) 4);
              nNewBadCount++;
              if (5 < m_nVerbose)
                {
                  if ( (100 > nNewBadCount) || (0 == (nNewBadCount % 5000) ) )
                    {
                      cout << "New bad pixel in reference image: " 
                           << i << ", " << j
                           << "  Flood, refer: " << fPixFlood << ", " 
                           << fPixRefer
                           << endl;
                    }
                }
            }
        }
    }

  sc = 32000. / max_ratio;
  cout << "\n      Max reference/flood ratio: " << max_ratio
       << "\n      Min reference/flood ratio: " << min_ratio
       << "\n      Scale factor             : " <<  sc
       << "\n      Number of new bad pixels : " << nNewBadCount
       << endl;

  // Compute a scale factor that makes the numbers in the nonunf img fall
  // in a range from 100 to 32000.  
  // Well isn't that what sc = 32000. / max_ratio above does???

  int icount = 0;
  float fPixNonunf;
  short int iBadPix;

  m_poImgBadPix->nSetNextPixel(0,0);
  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
        {
          // This could be done faster if we know the data type of
          // the arrays

          fPixFlood   = (m_poImgFlood->*m_poImgFlood->prfGetPixel)(i, j);
          iBadPixFlag = m_poImgBadPix->iGetNextPixel();

          // In the original code, if a pixel was bad (by iBadPixFlag)
          // then iBadPixFlag was placed into the output

//          fPixNonunf  = 2.0;
          fPixNonunf  = (float)iBadPixFlag;
          if ( (0 == iBadPixFlag) && (0.0 != fPixFlood) )
            {
              fPixRefer   = (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(i, j);
              r           = fPixRefer / fPixFlood;
              if (r < max_ratio)
                {
                  fPixNonunf = r * sc;
                }
              else
                {
                  if (100 > icount)
                    {
                      cout << "\n      Max ratio found at pixel : " << i << ", " << j
                           << "\n      Flood value              : " << fPixFlood
                           << "\n      Reference value          : " << fPixRefer
                           << "\n      Ratio Refer/Flood        : " <<  r
                           << "\n      Bad pixel flag           : " << iBadPixFlag
                           << endl;
                    }
                  if (r != max_ratio)
		    {
		      // Careful here.  Do we really want to set these pixels as bad?
		      // fPixNonunf = 4.0; //  nonunf(i,j) = flags(1)
		    }
                  icount++;
                }
            }
          (m_poImgNonunf->*m_poImgNonunf->prnSetPixel)(i, j, fPixNonunf + (float)0.5);
        }
    }

  // Calculate an average around the center in order to compute the
  // denominator of a correction factor:  fact = nonunf(i,j) / denom
  // Store denom in flags(5) so calling program has access to it.
  // If the center has alot of bad pixels, then shift away from center
  // by iyoff, ixoff until away from bad pixels

  int n;
  int ixoff, iyoff;

  ixoff = 0;
  iyoff = 0;
  n     = 0;
  while (  (200 > n) && (ixoff < nDim0) && (iyoff < nDim1) )
    {
      n  = 0;
      sc = 0.0;
      for (j = (nDim1 / 2) - 10 + iyoff; j <= (nDim1 / 2) + 10 + iyoff; j++)
        {
          for (i = (nDim0 / 2) - 10 + ixoff; i <= (nDim0 / 2) + 10 + ixoff; i++)
            {
              (void) m_poImgBadPix->nGetPixel(i, j, &iBadPixFlag);
              fPixNonunf = (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(i, j);
              if ( (0 == iBadPixFlag) && (0.0 < fPixNonunf) )
                {
                  sc = sc + fPixNonunf;
                  n++;
                }
            }
        }
      if (0 < n)
        {
          m_poImgNonunf->m_oHeader.nReplaceValue(D_K_NonunfFlag1, (int)1);
          m_poImgNonunf->m_oHeader.nReplaceValue(D_K_NonunfFlag2, (int)2);
          m_poImgNonunf->m_oHeader.nReplaceValue(D_K_NonunfFlag3, (int)3);
          m_poImgNonunf->m_oHeader.nReplaceValue(D_K_NonunfFlag4, (int)4);
          m_poImgNonunf->m_oHeader.nReplaceValue(D_K_NonunfOffset,
                                                 m_fdarksub_const);
          sc = sc / (float)n;
          m_poImgNonunf->m_oHeader.nReplaceValue(D_K_NonunfDenominator,
                                                   nint(sc));
//          flags(5) = nint(sc);
        }
      else
        {
//          flags(5) = -1;
        }

      // Move offsets just in case we had too many bad pixels where we just were
      // We should check that we have not moved to far

      if (3 < m_nVerbose)
        {
          cout << "\n      Denominator calculation"
               << "\n      X,Y offsets from center  : " << ixoff << ", " << iyoff 
               << "\n      Number of pixels used    : " << n
               << "\n      Average (= denominator)  : " << nint(sc)
               << endl;
        }
      ixoff = ixoff + nDim0 / 5;
      iyoff = iyoff + nDim1 / 5;
    }

//+jwp
  // THE FOLLOWING DOES NOT TREAT INTERNAL EDGES OF MULTI-MODULE DETECTORS!

  // Remove some edge effects if 0 < m_a3nmax_pixel[2]
  // 1. On j = nDim1-1 edge, set first non-bad pixel to be bad
  // 2. On i = nDim0-1 edge, set first non-bad pixel to be bad
  // 3. Set second pixel to be bad if m_a3nmaxPixel[2] is > 1.
  
  int nNumBad = 0;
  if (0 >= m_a3nmax_pixel[2])
    {
      cout << "INFO: No high edge treatment performed." << endl;
    }
  else
    {
      cout << "Performing high edge treatment type: " << m_a3nmax_pixel[2] << endl;
      for (i = 0; i < nDim0; i++)
        {
          j = nDim1 - 1;
          while (   (0 < j)
                 && (4.0 >= (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(i, j)))
            {
              j--;
            }
          (m_poImgNonunf->*m_poImgNonunf->prnSetPixel)(i, j,   (float)4.0);
          if (1 < m_a3nmax_pixel[2])
            (m_poImgNonunf->*m_poImgNonunf->prnSetPixel)(i, j-1, (float)4.0);
        }
      for (j = 0; j < nDim1; j++)
        {
          i = nDim0 - 1;
          while (   (0 < i)
                 && (4.0 >= (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(i, j)))
            {
              i--;
            }
          (m_poImgNonunf->*m_poImgNonunf->prnSetPixel)(i, j, (float)4.0);
          if (1 < m_a3nmax_pixel[2])
            (m_poImgNonunf->*m_poImgNonunf->prnSetPixel)(i-1, j, (float)4.0);
        }
    }

  // More edge badness.
  // On the other edges (i=0 and j=0), if the second pixel is good,
  // but the third pixel is bad, then make the second pixel bad.

  for (j = 0; j < nDim1; j++)
    {
      if (  (4.0 >= (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(2, j))
          && (4.0 < (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(1, j)) )
        {
          (m_poImgNonunf->*m_poImgNonunf->prnSetPixel)(1, j, (float)4.0);
          nNumBad++;
      }
    }

  for (i = 0; i < nDim0; i++)
    {
      if (  (4.0 >= (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(i, 2))
          && (4.0 < (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(i, 1)) )
        {
          (m_poImgNonunf->*m_poImgNonunf->prnSetPixel)(i, 1, (float)4.0);
          nNumBad++;
        }
    }
//-jwp

  cout << "Number of new bad edge pixels: " << nNumBad << endl;
  if (3 < m_nVerbose)
    cout << "... returning from nScaleImages!" << endl;
  return (nStat);
}

int
Ccalibrate::nCorrectImage(const int nCorrectMask, Cimage *poImageIn,
                          Cimage *poImageOut)
{
  // Correct a single image for dark, non-uniformity and/or spatial distortion.

  int nStat = 1;

  if ( (NULL == poImageIn) || (!poImageIn->bIsAvailable()) )
    {
      cout << "ERROR, INPUT image unavailable!" << endl;
      nStat = -1;
      return (nStat);
    }

  if ( (NULL == poImageOut) || (!poImageOut->bIsAvailable()) )
    {
      cout << "ERROR, OUTPUT image unavailable!" << endl;
      nStat = -1;
      return (nStat);
    }

  cout << "\n      Corrections to make: >>> ";
  if (0 != (nCorrectMask & CAL_DO_DARK) )
    cout << " dark";
  if (0 != (nCorrectMask & CAL_DO_TRANSFORM) )
    cout << " transform";
  if (0 != (nCorrectMask & CAL_DO_NONUNF) )
    cout << " nonunf";
  if (0 != (nCorrectMask & CAL_DO_SPATIAL) )
    cout << " spatial_distortion";
  cout << " <<<" << endl;

  nStat = 0;
  if (   (0 != (nCorrectMask & CAL_DO_DARK) )
      && (0.0 == m_fdarksub_scale) )
    {
      cout << "WARNING dark scale is 0, so NO dark subtraction performed!"
           << endl;
    }

  // Special case of DO_TRANSFORM.
  // This means DO_NONUNF and DO_SPATIAL are automatically excluded.

  if (   (0 != (nCorrectMask & CAL_DO_TRANSFORM) )
      && (0 == nStat) )
    {
      // 1.  Apply non-uniformity and spatial correction in one go
      //      with the transform information.  Deal with possible dark.

      Cimage *poImageDark = NULL;

      if (   (0 != (nCorrectMask & CAL_DO_DARK) )
          && (0.0 != m_fdarksub_scale) )
        {
          // Need to worry about the dark image
          if (NULL == m_poImgDark)
            {
              cout << "Reading dark image... " << m_sDarkName << endl;
              m_poImgDark = new Cimage(m_sDarkName);
            }
          if (!m_poImgDark->bIsAvailable())
            {
              cout << "ERROR, dark image unavailable!" << endl;
              nStat = -1;
              return (nStat);
            }
          poImageDark = m_poImgDark;
        }
      cout << "nCorrectImage just before ::nTransform(...)\n" << endl;
      nStat = nTransform(poImageIn, poImageOut, poImageDark);
      nStat = 0;
      if (0 != nStat)
        {
          cout << "ERROR performing transform!" << endl;
        }
      ///////////////
      return (nStat);
      ///////////////
    }

  // Get here only if NOT using transform image

#ifdef CORRECT_ONLY
  cout << "\nWARNING, newer & better algorithms exist for correcting images!\n";
#endif
  if (   (0 != (nCorrectMask & CAL_DO_DARK) )
      && (0.0 != m_fdarksub_scale) 
      && (0 == nStat) )
    {
      // 1. Subtract DARK image

      cout << "\n***   Correct image for dark image..." << endl;

      if (NULL == m_poImgDark)
        {
          cout << "Reading dark image... " << m_sDarkName << endl;
          m_poImgDark = new Cimage(m_sDarkName);
        }
      if (!m_poImgDark->bIsAvailable())
        {
          cout << "ERROR, dark image unavailable!" << endl;
          nStat = -1;
          return (nStat);
        }
      Cimage *poImageResult;
      if (nCorrectMask == CAL_DO_DARK)
        poImageResult = poImageOut;  // Ensures result is in poImageOut
      else
        poImageResult = NULL;         // Ensures result is in poImageIn
      nStat = nSubtract(poImageIn, m_poImgDark, m_fdarksub_scale,
                        m_fdarksub_const, (float)m_a3nmin_pixel[2],
                        poImageResult);
      if (0 != nStat)
        {
          cout << "ERROR performing Image - DARK!" << endl;
          return (nStat);
        }
    }

  if (   (0 != (nCorrectMask & CAL_DO_NONUNF) )
      && (0 == nStat) )
    {
      // 2. Apply non-uniformity correction

      Cimage *poImageResult;
      if (0 == (nCorrectMask & CAL_DO_SPATIAL) )
        poImageResult = poImageOut;   // Ensures result is in poImageOut
      else
        poImageResult = NULL;         // Ensures result is in poImageIn
      nStat = nNormalize(poImageIn, poImageResult);
    }

  if (   (0 != (nCorrectMask & CAL_DO_SPATIAL) )
      && (0 == nStat) )
    {
      // 3. Apply spatial distortion correction

      if (   (NULL == m_poImgDistorXINT) || (!m_poImgDistorXINT->bIsAvailable())
          || (NULL == m_poImgDistorYINT) || (!m_poImgDistorYINT->bIsAvailable()) )
        {
          nStat = nReadDistor();
          if (0 != nStat)
            {
              cout << "ERROR Reading spatial distortion info!" << endl;
              return (nStat);
            }
        }
      nStat = nInterpolate(&m_tCorrectParams, 
                           m_poImgDistorXINT, m_poImgDistorYINT,
                           poImageIn, poImageOut);
      
      (void) nUpdateHeaderCorrectedImage(poImageIn, poImageOut);
    }
  return (nStat);
}

int
Ccalibrate::nGoCorrect(const Cstring sImgNameIn, const Cstring sImgNameOut)
{
  int     nStat = 1;
  Cimage *poImageIn  = NULL;
  Cimage *poImageOut = NULL;

  cout << "\n***Correcting image... " << sImgNameIn << " to " << sImgNameOut
       << " ..." << endl;

  poImageIn = new Cimage(sImgNameIn);
  if (poImageIn->bIsAvailable())
    {
      // Check if the corrected image is the same size as the input
      // image

      if ( (0 >= m_a2nTOutDim[0]) || (0 >= m_a2nTOutDim[1]) )
        {
          poImageIn->nGetDimensions(&m_a2nTOutDim[0], &m_a2nTOutDim[1]);
        }

      if (   (0 == (m_nCorrectionMask & CAL_DO_TRANSFORM) )
          || (   (m_a2nTOutDim[0] == poImageIn->nGetDimension(0))
              && (m_a2nTOutDim[1] == poImageIn->nGetDimension(1))) )
        {
          // NOT using transform algorithm, so MUST be same size
          // OR it is the same size
          // so simply create the image from the input image

          poImageOut = new Cimage (*poImageIn);
        }
      else
        {
          // Different size, so create the output image in special way.
          // This is probably only valid with the TRANSFORM algorithm.

          cout << "INFO! Output image not same size as input image!\n"
               << flush;
          poImageOut = new Cimage(m_a2nTOutDim[0], m_a2nTOutDim[1],
                                  eImage_uI2);
          Cimage_header oNewHeader;
          oNewHeader = poImageIn->m_oHeader;
          oNewHeader.nReplaceValue(Cimage_header::ms_sSize1, m_a2nTOutDim[0]);
          oNewHeader.nReplaceValue(Cimage_header::ms_sSize2, m_a2nTOutDim[1]);
          poImageOut->nSetHeader(&oNewHeader);
          poImageOut->vSetState(eImage_available_state);
        }
      nStat = nCorrectImage(m_nCorrectionMask, poImageIn, poImageOut);
      if (0 == nStat)
        {
          nStat = poImageOut->nWrite(sImgNameOut);
          if (0 != nStat)
            {
              cout << "ERROR writing image " << sImgNameOut << " !" << endl;
            }
        }
      else
        {
          cout << "ERROR correcting image " << sImgNameIn << " !" << endl;
        }
      delete poImageOut;
      
      // These next 2 lines reset the output image file dimensions
      // Therefor if -tlimits was used, it has to be used for each -gocorrect
      // option

      m_a2nTOutDim[0] = 0;
      m_a2nTOutDim[1] = 0;
    }
  else
    {
      cout << "ERROR: Input image " << sImgNameIn << " unavailable!" << endl;
    }

  delete poImageIn;
  return (nStat);
}

int
Ccalibrate::nGoCorrectScan(const Cstring sScanNameIn,
                           const Cstring sScanNameOut)
{
  int     i;
  int     nStat      = 1;
  Cimage *poImageIn  = NULL;
  Cimage *poImageOut = NULL;
  Cscan  *poScanIn   = NULL;
  Cscan  *poScanOut   = NULL;
  Cstring sImageNameIn;
  Cstring sImageNameOut;
  int     nNumImages = 0;
  
  cout << "\n***Correcting scan template of images... " << sScanNameIn 
       << " to " << sScanNameOut
       << " ..." << endl;

  sImageNameIn = sScanNameIn;
  poScanIn     = new Cscan(sScanNameIn, m_a2nSeqNum[0], 1);

  // Find first image in scan by looking at the disk files

  poScanIn->vSetSeqStart(m_a2nSeqNum[0]);
  poScanIn->vSetImgNum(0);
  for (i = m_a2nSeqNum[0]; (i <= m_a2nSeqNum[1]) && (1 == nStat); i++)
    {
      poScanIn->nGetImageName(&sImageNameIn);
      if (bFileExists(sImageNameIn))
        {
          nStat = 0;
        }
      else
        {
          poScanIn->vNextSeqNum();
        }
    }
  if (0 != nStat)
    {
      cout << "ERROR no image match scan template: " << sScanNameIn << endl;
      delete poScanIn;
      return (nStat);
    }

  sImageNameOut = sScanNameOut;
  poScanOut     = new Cscan(sScanNameOut, m_a2nSeqNum[0], 1);

  // Loop through all images that can be found

  Cimage_header oNewHeader;
  poImageIn = new Cimage(sImageNameIn);

  while (poImageIn->bIsAvailable()  && (0 == nStat) )
    {
      // Check if the corrected image is the same size as the input image

      if ( (0 >= m_a2nTOutDim[0]) || (0 >= m_a2nTOutDim[1]) )
        {
          poImageIn->nGetDimensions(&m_a2nTOutDim[0], &m_a2nTOutDim[1]);
        }

      if (   (0 == (m_nCorrectionMask & CAL_DO_TRANSFORM) )
          || (   (m_a2nTOutDim[0] == poImageIn->nGetDimension(0))
              && (m_a2nTOutDim[1] == poImageIn->nGetDimension(1))) )
        {
          // NOT using transform algorithm, so MUST be same size
          // OR it is the same size
          // so simply create the image from the input image

          if (NULL == poImageOut)
            poImageOut = new Cimage (*poImageIn);
        }
      else
        {
          // Different size, so create the output image in special way.
          // This is probably only valid with the TRANSFORM algorithm.

          cout << "INFO! Output image not same size as input image!\n"
               << flush;
          if (NULL == poImageOut)
            poImageOut = new Cimage(m_a2nTOutDim[0], m_a2nTOutDim[1],
                                    eImage_uI2);
        }

      // Copy input image header to output image header

      oNewHeader = poImageIn->m_oHeader;
      oNewHeader.nReplaceValue(Cimage_header::ms_sSize1, m_a2nTOutDim[0]);
      oNewHeader.nReplaceValue(Cimage_header::ms_sSize2, m_a2nTOutDim[1]);
      poImageOut->nSetHeader(&oNewHeader);
      poImageOut->vSetState(eImage_available_state);

      nStat = nCorrectImage(m_nCorrectionMask, poImageIn, poImageOut);
      if (0 == nStat)
        {
          poScanOut->vSetSeqNum(poScanIn->nGetSeqNum());
          poScanOut->nGetImageName(&sImageNameOut);
          nStat  = poImageOut->nWrite(sImageNameOut);
          if (0 != nStat)
            {
              cout << "ERROR writing image " << sImageNameOut << " !" << endl;
            }
          else
            {
              cout << "SUCCESS writing image " << sImageNameOut << " !" << endl;
              nNumImages++;
            }
        }
      else
        {
          cout << "ERROR correcting image " << sImageNameIn << " !" << endl;
        }
      // Get next image

      poScanIn->vNextSeqNum();
      if (m_a2nSeqNum[1] >= poScanIn->nGetSeqNum())
        {
          nStat = poScanIn->nGetImage(poImageIn);
        }
      else
        poImageIn->vSetState(eImage_unknown_state); // Force end of while
    }
  cout << "Number of images corrected: " << nNumImages << endl;

  if (NULL != poImageIn)
    delete poImageIn;
  if (NULL != poImageOut)
    delete poImageOut;
  if (NULL != poScanIn)
    delete poScanIn;
  if (NULL != poScanOut)
    delete poScanOut;
  return (nStat);
}

int
Ccalibrate::nNormalize(Cimage *poImageIn, Cimage *poImageOut)
{
  // Apply the non-uniformity correction to the input image *poImageIn.
  // Place the result in *poImageOut if (NULL != poImageOut), otherwise
  // place the result in *poImageIn.

  int nStat = 0;
  cout << "\n***   Correct image for non-uniformity of response..." << endl;

  if ( (NULL == poImageIn) || !poImageIn->bIsAvailable())
    {
      cout << "ERROR in nNormalize, input image not available!" << endl;
      nStat = -1;
      return (nStat);
    }
  if ( (NULL != poImageOut) && !poImageOut->bIsAvailable())
    {
      cout << "ERROR in nNormalize, output image not available!" << endl;
      nStat = -1;
      return (nStat);
    }
  if (NULL == m_poImgNonunf)
    {
      cout << "Reading nonunf image... " << m_sNonunfName << endl;
      m_poImgNonunf = new Cimage(m_sNonunfName);
    }
  if (!m_poImgNonunf->bIsAvailable())
    {
      cout << "ERROR, nonunf image unavailable!" << endl;
      nStat = -1;
      return (nStat);
    }

  Cimage *poImageResult;

  if (NULL == poImageOut)
    poImageResult = poImageIn;
  else
    poImageResult = poImageOut;

  float fDenom;
  float fOffset;

  nStat = m_poImgNonunf->m_oHeader.nGetValue(D_K_NonunfDenominator, &fDenom);
  if (0 != nStat)
    {
      cout << "ERROR, getting nonunf denominator!" << endl;
      nStat = -1;
      return (nStat);
    }

  nStat = m_poImgNonunf->m_oHeader.nGetValue(D_K_NonunfOffset, &fOffset);
  if (0 != nStat)
    {
      fOffset = 0.0;
      nStat   = 0;
    }

  cout << "      Denominator used: " << fDenom
       << "\n      Offset used     : " << fOffset
       << endl;

  int   i, j;
  int   nDim0, nDim1;
  float fPixNonunf;
  float fPixIJ;
  float fSatVal;
  float afFlags[4] = { 1., 2., 3., 4. };

  fSatVal = poImageIn->fGetSatValue();

  m_poImgNonunf->nGetDimensions(&nDim0, &nDim1);
  m_poImgNonunf->m_oHeader.nGetValue(D_K_NonunfFlag1, &afFlags[0]);
  m_poImgNonunf->m_oHeader.nGetValue(D_K_NonunfFlag2, &afFlags[1]);
  m_poImgNonunf->m_oHeader.nGetValue(D_K_NonunfFlag3, &afFlags[2]);
  m_poImgNonunf->m_oHeader.nGetValue(D_K_NonunfFlag4, &afFlags[3]);

  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
        {
          fPixNonunf =         (m_poImgNonunf->*m_poImgNonunf->prfGetPixel)(i, j);
          fPixIJ     =         (poImageIn->*poImageIn->prfGetPixel)(i, j);
          if (   (afFlags[0] == fPixNonunf) 
              || (afFlags[1] == fPixNonunf) 
              || (afFlags[2] == fPixNonunf) 
              || (afFlags[3] == fPixNonunf) )
            {
              // Change pixels with bad non-uniformity to have a value of 
              // fOffset, so they do not contribute later on

              fPixIJ = fOffset;
            }
          else if (fSatVal != fPixIJ)
            {
              // Otherwise undo any applied offset before correcting 
              // for nonunf and then put the offset back in

              fPixIJ = (fPixIJ - fOffset) * fPixNonunf / fDenom;
              fPixIJ = fPixIJ + fOffset;
              if (fSatVal < fPixIJ)
                fPixIJ = fSatVal;
              else if (0.0 > fPixIJ) 
                fPixIJ = 0.0;   // Truncate to 0.0
            }
          (poImageResult->*poImageResult->prnSetPixel)(i, j, fPixIJ);
        } // end i loop
    } // end j loop
  return (nStat);
}

int
Ccalibrate::nMakeMarks(void)
{
  //    Based on the flood field, figure out which pixels are good and
  //    to which module they belong.
  //    Marty Stanton
  //    July 1994
      
  int i, j, ii, jj;
  int nDim0, nDim1;
  int nModDim0, nModDim1;
  int nNumGood = 0;
  int nNumBad  = 0;
  int nStat    = 0;

  cout << "\n***Make a mark image from flood image..." << endl;

  if (NULL == m_poImgFlood)
    {
      cout << "***   Reading flood image... " << m_sFloodName << endl;
      m_poImgFlood = new Cimage(m_sFloodName);
    }
  if (!m_poImgFlood->bIsAvailable())
    {
      cout << "ERROR, flood image unavailable!" << endl;
      nStat = -1;
      return (nStat);
    }

  if (NULL != m_poImgMark)
    delete m_poImgMark;

  m_poImgMark = new Cimage(*m_poImgFlood);

  m_poImgFlood->nGetDimensions(&nDim0, &nDim1);

  nModDim0 = nDim0 / m_a2nNumModules[0];
  nModDim1 = nDim1 / m_a2nNumModules[1];
  float fPixIJ;

  for (j = 0; j < nDim1; j++)
    {
      for (i = 0; i < nDim0; i++)
        {
          fPixIJ = (m_poImgFlood->*m_poImgFlood->prfGetPixel)(i, j);
          if (fPixIJ >= (float)m_a3nmin_pixel[FLOODIMG])
            {
              nNumGood++;
              fPixIJ =   (float)(m_a2nNumModules[1] * (j / nModDim1)
                                                    + (i / nModDim0) 
                                                    + 1) * 1000.0;
            }
          else
            {
              nNumBad++;
              fPixIJ = 0.0;
            }
          (m_poImgMark->*m_poImgMark->prnSetPixel)(i, j, fPixIJ);
        }
    }

  cout <<   "      Image size        : " << nDim0 << ", " << nDim1
       << "\n      Number of modules : " << m_a2nNumModules[0] << ", "
                                         << m_a2nNumModules[1]
       << "\n      Module size       : " << nModDim0 << ", " << nModDim1
       << "\n      Min pixel value   : " << m_a3nmin_pixel[FLOODIMG]
       << "\n      Number bad pixels : " << nNumBad
       << "\n      Number good pixels: " << nNumGood;

  // Go through mark image. 
  // If you don't have at least 3 good neighbors, mark a good pixel as bad

  nNumBad = 0;
  for (j = 0; j < nDim1; j++)
    {
      if ( (3 < m_nVerbose) && (0 == (j % 100) ) )
        cout << "          ... examining " << j << " of " << nDim1 << "..." << endl;
      for (i = 0; i < nDim0; i++)
        {
          fPixIJ   = (m_poImgMark->*m_poImgMark->prfGetPixel)(i, j);
          if (0.0 != fPixIJ)
            {
              nNumGood = 0;
              for (jj = max(0, j-1); (jj <= min(j+1, nDim1-1))
                                  && (4 > nNumGood); jj++)
                {
                  for (ii = max(0, i-1); (ii <= min(i+1, nDim0-1))
                                      && (4 > nNumGood); ii++)
                    {
                      fPixIJ   = (m_poImgMark->*m_poImgMark->prfGetPixel)(ii, jj);
                      if (0.0 < fPixIJ)
                        nNumGood++;
                    }
                }
              if (3 >= nNumGood)
                {
                  (m_poImgMark->*m_poImgMark->prnSetPixel)(i, j, 0.0);
                  nNumBad++;
                }
            }
        }
    }
  cout << "\n      Number of new bad : " << nNumBad << endl;
  m_poImgMark->vSetState(eImage_available_state);
//  cout << "About to write mark image!\n" << flush;
  nStat = m_poImgMark->nWrite(m_sMarkName);
  return (nStat);
}

int
Ccalibrate::nMergeTables(void)
{
  // Merge interpolation tables
  
  //  0. Initial module spatial distortion info must be available.
  //  1. If does not exist, create big spatial distortion tables of 
  //     size nNumMod0 by nNumMod1 modules.
  //  2. Set all values in the big tables to BAD_PIXEL.
  //  3. Read a distor file of a reference module
  //  4. Get module size in pixels (2 nums),
  //     get reference pixel position (2nums), pixel size to use in big tables
  // Start_loop:
  //  5. Get module spd basename, fiducial pixel posn, mm shift, pixel shift
  //  6. Read spd info of module
  //  7. Get shift of module from somewhere
  //  8. Shift the mm values of the module 
  //  9. Merge module with a pixel shift into the big tables
  // 10. Loop back
  // End_loop
  // 11. Fixup boundaries of modules
  // 12. Read mark image and spatially undistor it.
  // 13. Create inverse tables
  //  ?.
  //  ?. Write out merged tables

  int i, j;
  int nStat = 0;
  int nDim0, nDim1;
  int nDimMerge0, nDimMerge1;
  LONG    lIntBadFlag;
  LONG    lPixVal;
  LONG    lMin0PixVal, lMax0PixVal, lMin1PixVal, lMax1PixVal;
  LONG    lShiftMM0, lShiftMM1;

  int ipx, jpx, ioff, joff;
  int nMod;
  float fShiftMM0, fShiftMM1;
  float fShiftMM0Unscaled, fShiftMM1Unscaled;
  float fShiftPx0, fShiftPx1;
  float fShiftMMOv0, fShiftMMOv1;
  float fFidPx0, fFidPx1;
  float fFidMM0, fFidMM1;
  float fRefMM0, fRefMM1;
  float fRefPx0, fRefPx1;
  float fFidRefMM0, fFidRefMM1;
  float fMin0, fMin1, fExt0, fExt1;

  LONG lTest;
  int  nNumMod0, nNumMod1; 
  int  nModSizePX0, nModSizePX1;
  Cspatial *poSpatial;
  float fRefPixSize;
  float fOutPixSize;
  
  Cstring sDistorName;

  cout << "\n***Merge module distor tables into single set..." << endl;

  if (   (NULL == m_poImgDistorXINT) || (!m_poImgDistorXINT->bIsAvailable())
      || (NULL == m_poImgDistorYINT) || (!m_poImgDistorYINT->bIsAvailable()) )
    {
      nStat = nReadDistor();
      if (0 != nStat)
        {
          cout << "ERROR Reading spatial distortion info!" << endl;
          return (nStat);
        }
    }

  // Delete any existing merged distor tables

  if (NULL != m_poImgMergeDistorXINT)
    {
      delete m_poImgMergeDistorXINT;
    }
  if (NULL != m_poImgMergeDistorXINV)
    {
      delete m_poImgMergeDistorXINV;
    }
  if (NULL != m_poImgMergeDistorYINT)
    {
      delete m_poImgMergeDistorYINT;
    }
  if (NULL != m_poImgMergeDistorYINV)
    {
      delete m_poImgMergeDistorYINV;
    }

  cout << "\n***   Opening merge info file... " << m_sMergeInfoName << endl;

  ifstream oIn( sTransSymbol(m_sMergeInfoName).string());

  // Do we need all these extra () in the next line?

  if (!(oIn.rdbuf()->is_open())) 
    {
      cout << "ERROR opening merge info file " << m_sMergeInfoName << "!\n" << endl;
      nStat = 3;
      return (nStat);
    }

  oIn >> nNumMod0 >> nNumMod1; 
 
  // Create the merged distor tables

  m_poImgDistorXINT->nGetDimensions(&nDim0, &nDim1);
  nDimMerge0 = nDim0 * nNumMod0;
  nDimMerge1 = nDim1 * nNumMod1;
  cout << "\n***   Creating merged spatial distor tables for a " 
       << nNumMod0 << " by " << nNumMod1 << " modular detector." 
       << "\n      Size of tables: " << nDimMerge0 << " by "
       << nDimMerge1 << endl;

  m_poImgMergeDistorXINT = new Cimage(nDimMerge0, nDimMerge1, eImage_I4);
  m_poImgMergeDistorXINV = new Cimage(nDimMerge0, nDimMerge1, eImage_I4);
  m_poImgMergeDistorYINT = new Cimage(nDimMerge0, nDimMerge1, eImage_I4);
  m_poImgMergeDistorYINV = new Cimage(nDimMerge0, nDimMerge1, eImage_I4);

  // Set all values in the big tables BADFLAG
  
  m_poImgMergeDistorXINT->nSetNextPixel(0,0);
  m_poImgMergeDistorXINV->nSetNextPixel(0,0);
  m_poImgMergeDistorYINT->nSetNextPixel(0,0);
  m_poImgMergeDistorYINV->nSetNextPixel(0,0);

  lIntBadFlag = (long) m_tCorrectParams.badflag;
//  cout << "bad flag is " << lIntBadFlag << endl;
  for (i = 0; i < (nDimMerge0 * nDimMerge1); i++)
    {
      m_poImgMergeDistorXINT->vSetNextPixel(lIntBadFlag);
      m_poImgMergeDistorXINV->vSetNextPixel(lIntBadFlag);
      m_poImgMergeDistorYINT->vSetNextPixel(lIntBadFlag);
      m_poImgMergeDistorYINV->vSetNextPixel(lIntBadFlag);
    }

  // Get reference distortion information

//  cout << "\nEnter module size in pixels: ";
  oIn >> nModSizePX0 >> nModSizePX1;
//  cout << "\nEnter reference position in pixels: ";
  oIn >> fRefPx0 >> fRefPx1;
//  cout << "\nEnter pixel size to use in output tables: ";
  oIn >> fOutPixSize;
//  cout << "\nEnter overall shift of mask in mm: "
  oIn >> fShiftMMOv0 >> fShiftMMOv1;

  fRefPixSize = m_fpixel_size;
  if (0.0 >= fOutPixSize)
    fOutPixSize = fRefPixSize;
  nStat = nStantonPixeltoMM(fRefPx0, fRefPx1, &fRefMM0, &fRefMM1);
  cout << "\nReference position in pixels: " << fRefPx0 << ", " << fRefPx1
       << "\nReference position in mms:    " << fRefMM0 << ", " << fRefMM1
       << "\nReference pixel size:         " << fRefPixSize
       << "\nMerged table pixel size:      " << fOutPixSize
       << "\nMerged table overall mm shift:" << fShiftMMOv0 
                                     << ", " << fShiftMMOv1
       << endl;

  // Skip a line in the merge info file that contains column info
  // sDistorName is used as a scratch var here

  getline(oIn, sDistorName);
  if (0 >= sDistorName.length())
    getline(oIn, sDistorName);

  // Now loop over each module

  for (nMod = 0; nMod < (nNumMod0 * nNumMod1) && (0 == nStat); nMod++)
    {
      cout << "\n***********************************************"
           << "\nModule " << nMod + 1;
//      cout << "\nEnter distor file basename: ";
      oIn >> sDistorName;
      
      cout << "\nModule " << nMod +1 << " distor basename:    " 
           << sDistorName << endl;
      m_sDistorName = sDistorName;
      nStat = nReadDistor();
      if (0 != nStat)
        {
          cout << "ERROR Reading spatial distortion info!" << endl;
          return (nStat);
        }
      
//      cout << "\nEnter fiducial position in pixels: ";
      oIn >> fFidPx0 >> fFidPx1;
//      cout << "\nEnter fiducial X,Y shift in mm: ";
      oIn >> fShiftMM0 >> fShiftMM1;

//+JWP 2008-09-03
// Scale the fiducial X,Y shift by the necessary factor due to the 
// fact that the mask is not up against the detector

      fShiftMM0Unscaled = fShiftMM0;
      fShiftMM1Unscaled = fShiftMM1;

      fShiftMM0 = fShiftMM0 * m_fxtod_distance
                             / (m_fxtod_distance - m_fmasktod_distance);
      fShiftMM1 = fShiftMM1 * m_fxtod_distance 
                             / (m_fxtod_distance - m_fmasktod_distance);

//      cout << "\nEnter module pixel shift in pixels: ";
      oIn >> fShiftPx0 >> fShiftPx1;
//      cout << "\nEnter minimum origin mm, extents mm in OUTPUT table: ";
      oIn >> fMin0 >> fMin1 >> fExt0 >> fExt1;

      nStat = nStantonPixeltoMM(fFidPx0, fFidPx1, &fFidMM0, &fFidMM1);
      cout << "\nModule " << nMod + 1 << " pixel size in mm:             " 
           << m_fpixel_size
           << "\nModule " << nMod + 1 << " requested pixel shift:        " 
           << fShiftPx0 << ", " << fShiftPx1
           << "\nModule " << nMod + 1 << " fiducial  position in pixels: " 
           << fFidPx0 << ", " << fFidPx1
           << "\nModule " << nMod + 1 << " fiducial  position in mm:     " 
           << fFidMM0 << ", " << fFidMM1
           << "\nModule " << nMod + 1 << " requested unscaled mm shift:  " 
           << fShiftMM0Unscaled << ", " << fShiftMM1Unscaled
           << "\nModule " << nMod + 1 << " requested actual mm shift:    " 
           << fShiftMM0 << ", " << fShiftMM1;
      fShiftMM0 = fShiftMM0 + fRefMM0 - fFidMM0;
      fShiftMM1 = fShiftMM1 + fRefMM1 - fFidMM1;
      cout << "\nModule " << nMod + 1 << " actual    mm shift:           " 
           << fShiftMM0 << ", " << fShiftMM1
           << "\nModule " << nMod + 1 << " allowed minimum mm:           "
           << fMin0 << ", " << fMin1
           << "\nModule " << nMod + 1 << " allowed extents mm:           "
           << fExt0 << ", " << fExt1
           << endl;

      // Shift the XINT and YINT tables

      nStat = nShiftMM(fShiftMM0 + fShiftMMOv0, fShiftMM1 + fShiftMMOv1, fOutPixSize);

      nStat = nStantonPixeltoMM(fFidPx0, fFidPx1, &fFidRefMM0, &fFidRefMM1);

      cout << "Module " << nMod + 1 << " new fiducial position in mm:  " 
           << fFidRefMM0 << ", " << fFidRefMM1
           << "\nModule " << nMod + 1 << " new pixel size in mm:         " 
           << m_fpixel_size
           << endl;

      // This expects xint_start, yint_start to be 1!

      if (   (1 != m_tCorrectParams.xint_start)
          || (1 != m_tCorrectParams.yint_start) )
        {
          cout << "WARNING! xint_start and/or yint_start is not 1!" << endl;
        }
      ioff = nint(fShiftPx0 / m_tCorrectParams.xint_step);
      joff = nint(fShiftPx1 / m_tCorrectParams.yint_step);

      // Now copy the shifted distor table into the merged table
      // One might not use the entire spatial distortion table

      nDim0 = nModSizePX0 / m_tCorrectParams.xint_step;
      nDim1 = nModSizePX1 / m_tCorrectParams.yint_step;

//      fMin0 = (fMin0 + fShiftMMOv0) / (m_tCorrectParams.pscale * fOutPixSize);
//      fMin1 = (fMin1 + fShiftMMOv1) / (m_tCorrectParams.pscale * fOutPixSize);
      fMin0 = fMin0 / (m_tCorrectParams.pscale * fOutPixSize);
      fMin1 = fMin1 / (m_tCorrectParams.pscale * fOutPixSize);
      fExt0 = fExt0 / (m_tCorrectParams.pscale * fOutPixSize);
      fExt1 = fExt1 / (m_tCorrectParams.pscale * fOutPixSize);

      lMin0PixVal = (LONG) fMin0;
      lMax0PixVal = (LONG) (fMin0 + fExt0 + 0.5);
      lMin1PixVal = (LONG) fMin1;
      lMax1PixVal = (LONG) (fMin1 + fExt1 + 0.5);
      if (6 < m_nVerbose)
        {
          cout << "Min, Max 0; Min, Max 1: "
               << lMin0PixVal << ", "
               << lMax0PixVal << ", "
               << lMin1PixVal << ", "
               << lMax1PixVal << endl;
        }

      for (j = 0; j < nDim1; j++)
        {
          jpx = j + joff;
          for (i = 0; i < nDim0; i++)
            {
              ipx = i + ioff;
              nStat = m_poImgDistorXINT->nGetPixel(i, j, &lPixVal);
//
// TODO: May not want to mark boundaries as bad!
//
              if (0 != nStat)
                {
                  lPixVal = lIntBadFlag;
                }
/*
              else if (   (0 == i) || ((nDim0-1) == i) 
                       || (0 == j) || ((nDim1-1) == j) )
                {
                  // Make sure boundaries between modules are marked bad

                  lPixVal = lIntBadFlag;
                }
*/
////#define DOMARKS 1
#ifndef DOMARKS
              else if (   (lPixVal != lIntBadFlag)
                       && ( (lPixVal < lMin0PixVal) || (lPixVal > lMax0PixVal) ) )
                {
                  // If outside of mm boundaries for this module, mark bad

                  lPixVal = lIntBadFlag;
                }
#endif

              if (lPixVal != lIntBadFlag)
                m_poImgMergeDistorXINT->nSetPixel(ipx, jpx, lPixVal);

              nStat = m_poImgDistorYINT->nGetPixel(i, j, &lPixVal);
//
// TODO: May not want to mark boundaries as bad?  Or mark only one side
//       of the border as bad?
//
              if (0 != nStat)
                {
                  lPixVal = lIntBadFlag;
                }
/*
              else if (   (0 == i) || ((nDim0-1) == i) 
                       || (0 == j) || ((nDim1-1) == j) )
                {
                  // Make sure boundaries between modules are marked bad

                  lPixVal = lIntBadFlag;
                }
*/
#ifndef DOMARKS
              else if (   (lPixVal != lIntBadFlag)
                       && ( (lPixVal < lMin1PixVal) || (lPixVal > lMax1PixVal) ) )
                {
                  // If outside of mm boundaries for this module, mark bad

                  lPixVal = lIntBadFlag;
                }
#endif
              if (lPixVal != lIntBadFlag)
                m_poImgMergeDistorYINT->nSetPixel(ipx, jpx, lPixVal);
            }
        }
      if (0 != nStat)
	{
	  cout << "WARNING, nStat != 0, for the above module.  Reset to 0!\n" << endl;
	  nStat = 0;
	}
    }
  oIn.close();                 // Make sure we always close the stream!

  cout << "\n***********************************************"
       << "\n***  ... done merging tables.\n" << endl;

//+jwp 14-Nov-1998 MAY NOT NEED THE FOLLOWING SECTION ANYMORE!
//                 MAY NOT NEED MARK IMAGE ANYMORE!

#ifdef DOMARKS

  // Try to resolve mm overlap among merged spatial distortion tables
  // If there is mm overlap, then perhaps the pixel size is wrong

  int j1, j2;
  int jlo, jhi;
  long lPixVal2;
  long lPixValSav1, lPixValSav2;
  for (j = nDim1-1; j < nDim1 * (nNumMod1-1); j = j+nDim1)
    {
//      cout << "fixup j is: " << j << endl;
      jlo = max(0, j - 20);
      jhi = min(j + 20, nDim1 * nNumMod1);

      for (i = 0; i < (nDim0 * nNumMod0); i++)
        {
          j1 = j;
          j2 = j + 1;
          lPixValSav1 = lIntBadFlag;
          lPixValSav2 = lIntBadFlag;

          bool bLoop = TRUE;
          while (bLoop)
            {
              lPixVal = lIntBadFlag;
              while ( (j1 > jlo) && (lPixVal == lIntBadFlag) )
                {
                  j1--;
                  m_poImgMergeDistorYINT->nGetPixel(i, j1, &lPixVal);
                }
              lPixVal2 = lIntBadFlag;
              while ( (j2 < jhi) && (lPixVal2 == lIntBadFlag) )
                {
                  j2++;
                  m_poImgMergeDistorYINT->nGetPixel(i, j2, &lPixVal2);
                }
//              cout << "l1, l2: " << lPixVal << ", " << lPixVal2 << endl;
              bLoop = (lPixVal2 < lPixVal)&& (j1 > jlo) && (j2 < jhi);
              if (bLoop)
                {
                  lPixValSav1 = lPixVal;
                  lPixValSav2 = lPixVal2;
                  m_poImgMergeDistorYINT->nSetPixel(i, j1, lIntBadFlag);
                  m_poImgMergeDistorYINT->nSetPixel(i, j2, lIntBadFlag);
                }
              else
                {
//                  cout << "i = " << i << "j1+1, j2-1: " << j1+1 << ", "
//                       << j2-1 << " sav1,2: " << lPixValSav1 << ", " 
//                       << lPixValSav2 << endl;
//                  lPixValSav1 = lPixValSav1 + lPixValSav2;
//                  lPixValSav1 = lPixValSav1 / 2;
/*
//1 set one pixel bad
                  if (abs(lPixVal-lPixValSav1) < abs(lPixVal2-lPixVal))
                    m_poImgMergeDistorYINT->nSetPixel(i, j2, lIntBadFlag);
                  else
                    m_poImgMergeDistorYINT->nSetPixel(i, j1, lIntBadFlag);
*/
/*
//2 Do nothing
*/
/*
//3 restore one of the previous ones
                  if (abs(lPixVal-lPixValSav1) < abs(lPixVal2-lPixVal))
                    m_poImgMergeDistorYINT->nSetPixel(i, j1+1, lPixValSav1);
                  else
                    m_poImgMergeDistorYINT->nSetPixel(i, j2-1, lPixValSav2);
*/
//4 restore both of the previous ones
                  m_poImgMergeDistorYINT->nSetPixel(i, j1+1, lPixValSav1);
                  m_poImgMergeDistorYINT->nSetPixel(i, j2-1, lPixValSav2);

                }
            }
        }
    }

  int ihi, ilo;
  int i1, i2;
  for (i = nDim0-1; i < nDim0 * (nNumMod0-1); i = i+nDim0)
    {
//      cout << "fixup i is: " << i << endl;
      ilo = max(0, i - 20);
      ihi = min(i + 20, nDim0 * nNumMod0);

      for (j = 0; j < (nDim1 * nNumMod1); j++)
        {
          i1 = i;
          i2 = i + 1;
          lPixValSav1 = lIntBadFlag;
          lPixValSav2 = lIntBadFlag;

          bool bLoop = TRUE;
          while (bLoop)
            {
              lPixVal = lIntBadFlag;
              while ( (i1 > ilo) && (lPixVal == lIntBadFlag) )
                {
                  i1--;
                  m_poImgMergeDistorXINT->nGetPixel(i1, j, &lPixVal);
                }
              lPixVal2 = lIntBadFlag;
              while ( (i2 < ihi) && (lPixVal2 == lIntBadFlag) )
                {
                  i2++;
                  m_poImgMergeDistorXINT->nGetPixel(i2, j, &lPixVal2);
                }
//              cout << "l1, l2: " << lPixVal << ", " << lPixVal2 << endl;
              bLoop = (lPixVal2 < lPixVal)&& (i1 > ilo) && (i2 < ihi);
              if (bLoop)
                {
                  lPixValSav1 = lPixVal;
                  lPixValSav2 = lPixVal2;
                  m_poImgMergeDistorXINT->nSetPixel(i1, j, lIntBadFlag);
                  m_poImgMergeDistorXINT->nSetPixel(i2, j, lIntBadFlag);
                }
              else
                {
//                  cout << "j = " << j << "i1+1, i2-1: " << i1+1 << ", "
//                       << i2-1 << " sav1,2: " << lPixValSav1 << ", " 
//                       << lPixValSav2 << endl;
/*
1
                  if (abs(lPixVal-lPixValSav1) < abs(lPixVal2-lPixVal))
                    m_poImgMergeDistorXINT->nSetPixel(i2, j, lIntBadFlag);
                  else
                    m_poImgMergeDistorXINT->nSetPixel(i1, j, lIntBadFlag);
*/
/*
2 Do nothing (this looks good)
*/

/*
//3 restore previous one
                  if (abs(lPixVal-lPixValSav1) < abs(lPixVal2-lPixVal))
                    m_poImgMergeDistorXINT->nSetPixel(i1+1, j, lPixValSav1);
                  else
                    m_poImgMergeDistorXINT->nSetPixel(i2-1, j, lPixValSav2);
*/
//4 Restore both previous ones
                  m_poImgMergeDistorXINT->nSetPixel(i1+1, j, lPixValSav1);
                  m_poImgMergeDistorXINT->nSetPixel(i2-1, j, lPixValSav2);
                }
            }
        }
    }
#endif

  if (0 == nStat)
    {
      // Set image state to available

      m_poImgMergeDistorXINT->vSetSatValue((float)lIntBadFlag);
      m_poImgMergeDistorYINT->vSetSatValue((float)lIntBadFlag);
      m_poImgMergeDistorXINT->vSetState(eImage_available_state);
      m_poImgMergeDistorYINT->vSetState(eImage_available_state);
    }

#ifdef DOMARKS

  // TODO: Mark border pixels as bad.

  // Do something with the mark image.

  Cimage *poImgMark, *poImgMarkD;
  poImgMark  = new Cimage(m_sMarkName);
  poImgMarkD = new Cimage(*poImgMark);
  poImgMarkD->vSetSatValue(32767.0);
  poImgMarkD->vSetState(eImage_available_state);
  
  // nMoveMark undistorts the mark image spatially, but keeps the pixel
  // values intact.

  nStat = nMoveMark(&m_tCorrectParams, 
                    m_poImgMergeDistorXINT, m_poImgMergeDistorYINT,
                    poImgMark, poImgMarkD);
  cout << "nstat from nMoveMark: " << nStat << endl;
//  poImgMarkD->m_oHeader.nDeleteMask("SCAN*");
  poImgMarkD->nWrite(m_sMarkName + "u");
  if (0 == nStat)
    nStat = nFillHoles(poImgMarkD);
  poImgMarkD->nWrite(m_sMarkName + "f");
  delete poImgMark;

#endif
  // Start reverse interpolation table

  cout << "***      Starting REVERSE interpolation tables..." << endl;

  LONG lPix[8];
  int  xintstart, xintstep;
  int  yintstart, yintstep;
  int  xinvstart, xinvstep;
  int  yinvstart, yinvstep;

  float pscale;

  int       timx, timy;
  float     x1, x2, x3, x4, y1, y2, y3, y4;
  int       nxmin, nxmax, nymin, nymax;
  int       ii, jj;
  float     a, b, c, m;
  float     dxpos, dypos;
  float     xp, yp;
  int       cmod, dmod;
  short int iPixIJ;
  int       indexi, indexj;
  float     xm_size, ym_size;

  pscale    = m_tCorrectParams.pscale;
  xintstart = m_tCorrectParams.xint_start;
  xintstep  = m_tCorrectParams.xint_step;
  yintstart = m_tCorrectParams.yint_start;
  yintstep  = m_tCorrectParams.yint_step;
  xinvstart = m_tCorrectParams.xint_start;
  xinvstep  = m_tCorrectParams.xint_step;
  yinvstart = m_tCorrectParams.yint_start;
  yinvstep  = m_tCorrectParams.yint_step;

  xm_size   = (float)nModSizePX0;
  ym_size   = (float)nModSizePX1;

  m_poImgMergeDistorXINT->nGetDimensions(&timx, &timy);

  if (3 < m_nVerbose)
    {
      cout << "\nReverse table params: "
           << "\nX,Y int start, step: " << xintstart << ", " << yintstart
                                        << ", " << xintstep << ", " << yintstep
           << "\nX,Y inv start, step: " << xinvstart << ", " << yinvstart
                                        << ", " << xinvstep << ", " << yinvstep
           << "\nPixel scale        : " << pscale
           << "\nModule size pixels: " << xm_size << ", " << ym_size
           << "\nInverse array size:  " << timx << ", " << timy
           << endl;
    }

  for (j = 1; j <= timy-1; j++)
    {
      for (i = 1; i <= timx-1; i++)
        {
          // All corners must have values in them

          m_poImgMergeDistorXINT->nGetPixel(i-1, j-1, &lPix[0]);
          m_poImgMergeDistorXINT->nGetPixel(i,   j-1, &lPix[1]);
          m_poImgMergeDistorXINT->nGetPixel(i-1, j,   &lPix[2]);
          m_poImgMergeDistorXINT->nGetPixel(i,   j,   &lPix[3]);
          m_poImgMergeDistorYINT->nGetPixel(i-1, j-1, &lPix[4]);
          m_poImgMergeDistorYINT->nGetPixel(i,   j-1, &lPix[5]);
          m_poImgMergeDistorYINT->nGetPixel(i-1, j,   &lPix[6]);
          m_poImgMergeDistorYINT->nGetPixel(i,   j,   &lPix[7]);

          if (   (lIntBadFlag != lPix[0])
              && (lIntBadFlag != lPix[1])
              && (lIntBadFlag != lPix[2])
              && (lIntBadFlag != lPix[3])
              && (lIntBadFlag != lPix[4])
              && (lIntBadFlag != lPix[5])
              && (lIntBadFlag != lPix[6])
              && (lIntBadFlag != lPix[7]) )
            {
/*
            if ( xint(i,j)     .le. ndim1 .and.
     $           xint(i+1,j)   .le. ndim1 .and.
     $           xint(i,j+1)   .le. ndim1 .and.
     $           xint(i+1,j+1) .le. ndim1 .and.
     $           yint(i,j)     .le. ndim2 .and.
     $           yint(i+1,j)   .le. ndim2 .and.
     $           yint(i,j+1)   .le. ndim2 .and.
     $           yint(i+1,j+1) .le. ndim2 ) then
*/            

              // Determine min, max positions that might be
              // within this box

              x1 = (float)lPix[0] * pscale;
              x2 = (float)lPix[1] * pscale;
              x3 = (float)lPix[2] * pscale;
              x4 = (float)lPix[3] * pscale;
              y1 = (float)lPix[4] * pscale;
              y2 = (float)lPix[5] * pscale;
              y3 = (float)lPix[6] * pscale;
              y4 = (float)lPix[7] * pscale;
            
              nxmin = (int) (min(x1, min(x2, min(x3, x4 ))));
              nxmax = (int) (max(x1, max(x2, max(x3, x4 )))) + 1;
              nymin = (int) (min(y1, min(y2, min(y3, y4 ))));
              nymax = (int) (max(y1, max(y2, max(y3, y4 )))) + 1;
               
              for (jj = nymin; jj <= nymax; jj++)
                {
                  for (ii = nxmin; ii <= nxmax; ii++)
                    {
                      // Is it a step pixel ?

                      if (   (0 == ((ii - xinvstart) % xinvstep))
                          && (0 == ((jj - yinvstart) % yinvstep)) )
                        {
                          // Is it really within this box ?

                          if (0 == nInBox( float(ii), float(jj), 
                                           x1, y1, x2, y2, 
                                           x3, y3, x4, y4 ))
                            {
                              // Do a simple linear interpolation

                              a = ((float)ii-x3) / (x4-x3) 
                                   - ((float)ii-x1)/(x2-x1);
                              c = ((float)ii-x1) / (x2-x1);
                              m = ((float)jj-y2) / (y4-y2) 
                                    - ((float)jj-y1)/(y3-y1);
                              b = ((float)jj-y1) / (y3-y1);
                              xp = (a * b  + c) / (1.0 - a * m);
                              yp = m * xp + b;
/*
    (I think)    

 Distorted:      (xintstep*(i+xp-1)+xintstart, yintstep*(j+yp-1)+yintstart)


 Corrected       (ii,jj)

     The last question to answer is if this position in corrected
     space really came from this position in distorted, or if
     it came from another

     To test, 1) figure out which module the distorted position
                 corresponds to
              2) compare with module that the corrected position
                 must come from (mark file)
     
*/
                              dxpos = xintstep * (i+xp-1) + xintstart;
                              dypos = yintstep * (j+yp-1) + yintstart;
                              dmod = int( dxpos / xm_size) + 
                                     nNumMod0 * int( dypos / ym_size) + 1;
#ifdef DOMARKS
                              poImgMarkD->nGetPixel(ii, jj, &iPixIJ);
#else
                              iPixIJ = dmod * 1000;
#endif
                              cmod = nint((float)iPixIJ / 1000.f);
                              if ( (5 < m_nVerbose) && (cmod != dmod) && (0 != cmod) )
                                {
                                  cout << "\nDx,ypos: " << dxpos << ", " << dypos
                                       << "\nii,jj: " << ii << ", " << jj
                                       << "\ndmod,cmod: " << dmod << ", " << cmod
                                       << "\ni, j:      " << i << ", " << j
                                       << "\nxp, yp:    " << xp << ", " << yp
                                       << endl;
                                }
                              else
                                {
                                  indexi = (ii - xinvstart) / xinvstep + 1;
                                  indexj = (jj - yinvstart) / yinvstep + 1;

                                  // Remember that nSetPixel checks bounds

                                  m_poImgMergeDistorXINV->nSetPixel(indexi-1, indexj-1,
                                   (LONG)nint((xintstep*((float)i+xp-1.0)
                                              +xintstart) / pscale));
                                  m_poImgMergeDistorYINV->nSetPixel(indexi-1, indexj-1,
                                   (LONG)nint((yintstep*((float)j+yp-1.0)
                                              +yintstart) / pscale));
                                }
/*
                                 xinv  (indexi, indexj) = 
     $                                nint((xintstep*(i+xp-1)+xintstart)
     $                                / packscale)
                                 yinv  (indexi, indexj) = 
     $                                nint((yintstep*(j+yp-1)+yintstart)
     $                                / packscale)
*/
                            }
                        }
                    }
                }
            }
        }
    }

  // TODO: We probably want to fill in "holes" between modules in the
  //       XINV, YINV arrays.

  if (0 == nStat)
    {
      m_poImgMergeDistorXINV->vSetSatValue((float)lIntBadFlag);
      m_poImgMergeDistorYINV->vSetSatValue((float)lIntBadFlag);
      m_poImgMergeDistorXINV->vSetState(eImage_available_state);
      m_poImgMergeDistorYINV->vSetState(eImage_available_state);
    }

  if (0 == nStat)
    {
//      cout << "...About to write merged tables!\n" << endl;
      nStat = m_poImgMergeDistorXINT->nWrite("merge" + ms_asExt[2]);
      if (0 == nStat)
        nStat = m_poImgMergeDistorYINT->nWrite("merge" + ms_asExt[3]);
      if (0 == nStat)
        nStat = m_poImgMergeDistorXINV->nWrite("merge" + ms_asExt[0]);
      if (0 == nStat)
        nStat = m_poImgMergeDistorYINV->nWrite("merge" + ms_asExt[1]);
      if (0 == nStat)
        {
          if (NULL == m_poHeaderSave)
            {
              m_poHeaderSave = new Cimage_header();
            }
          m_a2nNumModules[0]     = nNumMod0;
          m_a2nNumModules[1]     = nNumMod1;
          m_tDistorParams.x_beam = m_a2fbeam_position[0];
          m_tDistorParams.y_beam = m_a2fbeam_position[1];
          m_tDistorParams.x_size = nNumMod0 * nModSizePX0;
          m_tDistorParams.y_size = nNumMod1 * nModSizePX1;
          m_fpixel_size          = fOutPixSize;
          m_tCorrectParams.x_scale = 1.0 / m_fpixel_size;
          m_tCorrectParams.y_scale = 1.0 / m_fpixel_size;
          nUpdateHeader(m_poHeaderSave);
          if ( (NULL != m_poHeaderSave) && m_poHeaderSave->bIsAvailable() )
            {
              nStat = m_poHeaderSave->nWrite("merge" + ms_asExt[4]);
            }
          else
            {
              nStat = 3;
            }
          if (0 != nStat)
            {
              cout << "ERROR writing " << "merge" << ms_asExt[4]
                << endl;
            }
        }
    }

#ifdef DOMARKS
  if (NULL != poImgMarkD)
    delete poImgMarkD;
#endif
  return (nStat);
}

int 
Ccalibrate::nMoveMark(tagCorrect_parameters *ptCorrectParams,
                      Cimage *poImgDistorXINT, Cimage *poImgDistorYINT,
                      Cimage *poImageIn, Cimage *poImageOut)
{
  // This is the same code as ::nInterpolate, EXCEPT instead of interpolating
  // we just copy the input pixel to the undistorted output pixel.  The
  // first input pixel to make it into the output pixel wins!

  int    i, j;
  int    ii, jj;
  float  fPixOut;
  int    ixs, iys, ixe, iye;
  float  valu, valc;
  float  xp, yp;
  int    ix, iy;
  int    ioffst;

  int    xasize, yasize, xsize, ysize;
  int    xxasize, xyasize, xxsize, xysize;
  int    yxasize, yyasize, yxsize, yysize;
  int    xstart, xstep, ystart, ystep;
  float  pscale;
  int    error;

  float    area;
  float    pix_per_area, xtemp, ytemp;
  float    fSatVal;
  float    xn, yn;
  int    xs, ys;
  float  *pfY0;
  float  x0, fIntBadFlag;

  LONG lIntBadFlag;
  LONG lYINT[8];
  LONG lXINT[8];

  cout << "\n***   Correct a mark image for spatial distortions..." << endl;

  if ( (NULL == poImageIn) || (NULL == poImageOut) )
    {
      return (-1);
    }
  else if (!poImageIn->bIsAvailable() || !poImageOut->bIsAvailable())
    {
      return (-2);
    }
  else if (NULL == poImgDistorXINT)
    {
      return (-3);
    }
  else if (!poImgDistorXINT->bIsAvailable())
    {
      return (-4);
    }
  else if (NULL == poImgDistorYINT)
    {
      return (-3);
    }
  else if (!poImgDistorYINT->bIsAvailable())
    {
      return (-4);
    }

  fIntBadFlag = (float) ptCorrectParams->badflag;
  lIntBadFlag = (long)  ptCorrectParams->badflag;

  fSatVal = poImageIn->fGetSatValue();

  (void) poImageIn->nGetDimensions(&xasize, &yasize);
  xsize   = xasize;
  ysize   = yasize;

  (void) poImgDistorXINT->nGetDimensions(&xxasize, &xyasize);
  xxsize  = xxasize;
  yysize  = xyasize;
  yxasize = xxasize;
  yyasize = xyasize;

  xstart  = ptCorrectParams->xint_start;
  ystart  = ptCorrectParams->yint_start;
  xstep   = ptCorrectParams->xint_step;
  ystep   = ptCorrectParams->yint_step;
  pscale  = ptCorrectParams->pscale;

  if (3 < m_nVerbose)
    {
      cout <<   "      Image Size            : " << xsize   << ", " << ysize
           << "\n      Image Array Size      : " << xasize  << ", " << yasize
           << "\n      Int Table Size        : " << xxsize  << ", " << yysize
           << "\n      Int Table Array Size  : " << xxasize << ", " << xyasize
           << "\n      Int X start, step     : " << xstart  << ", " << xstep
           << "\n      Int Y start, step     : " << ystart  << ", " << ystep
           << "\n      Pixel scale, 1/pscale : " << pscale  << ", " << 1.0 / pscale
           << endl;
    }

  // Start by setting output image to the offset value
  // (This could be done faster another way, but would be
  //  specific for image data types)

  float fOffset;
  fOffset = 0.0;  // (float)ioffst;
  for (j = 0; j < yasize; j++)
    {
      for (i = 0; i < xasize; i++)
        {
          (poImageOut->*poImageOut->prnSetPixel)(i, j, fOffset);
        }
    }

  xs = min ( xsize, xstep*xxsize + xstart);
  ys = min ( ysize, ystep*yysize + ystart);

/*
C     Define lower edges for first row.  The lower edges for the
C     rest of the rows will just be equal to the upper edge of
C     the previous row.
C
C     Initialize y starting positions by assuming that the lower height
C     of the box is the same as the upper height of the box
C     This is to prevent getting a bad value from trying
C     an illegal position
*/
  int   nOffset;
  float fPixIJ;
  pfY0    = new float [xs+1];
  for (i = 1; i <= xs; i++)
    {
      xp = (i-xstart)/(float)xstep + 1.0;

//      yp = (1+0.5-ystart)/(float)xstep + 1.0;
      yp = (1.0-(float)ystart)/(float)xstep + 1.0;
      ix = (int)xp;
      iy = (int)yp;
      if (   (ix >= 1) && (ix < yxasize) 
          && (iy >= 1) && (iy < yyasize) ) 
        {
          nOffset  = ix +  iy * yxasize;
//          cout << "ix, iy, nOffset: " << ix << ", " << iy << ", " << nOffset << endl;
          lYINT[0] = poImgDistorYINT->lGetPixel(nOffset);  // (ix,   iy);
          lYINT[1] = poImgDistorYINT->lGetPixel(nOffset-1); //ix-1, iy);
          lYINT[2] = poImgDistorYINT->lGetPixel(nOffset-yxasize); //ix, iy-1);
          lYINT[3] = poImgDistorYINT->lGetPixel(nOffset-yxasize-1);//ix-1, iy-1);

          if (   (lIntBadFlag == lYINT[0])
              || (lIntBadFlag == lYINT[1])
              || (lIntBadFlag == lYINT[2])
              || (lIntBadFlag == lYINT[3]) )
            {
              pfY0[i] = fIntBadFlag;
            }
          else
            {
              valu = (float)((yp-iy)*
                   ((xp-ix) * lYINT[0]             //y_int(ix+1,iy+1) 
                   + (1+ix-xp) * lYINT[1]) +         //y_int(ix,iy+1)) +
                   (1+iy-yp)*
                   ((xp-ix)* lYINT[2]               //y_int(ix+1,iy)
                   + (1+ix-xp)* lYINT[3]))         //y_int(ix,iy)
                    * pscale;
               yp = (1-ystart)/ystep + 1.0;
               iy = int(yp);
               valc = (float)((yp-iy)*
                   ((xp-ix)* lYINT[0]                 //y_int(ix+1,iy+1)
                   + (1+ix-xp)*lYINT[1]) +            //y_int(ix,iy+1)) +
                   (1+iy-yp)*
                   ((xp-ix)* lYINT[2]                //y_int(ix+1,iy)
                   + (1+ix-xp)* lYINT[3]))           // y_int(ix,iy)))
                             * pscale;
               pfY0[i] = 2.0 * valc  -  valu;
            }
        }
      else
        {
          pfY0[i] = fIntBadFlag;
        }
    }  // end i loop

  // Start of main loop

  poImageIn->nSetNextPixel(0,0);
  for (j = 1; j <= ys; j++)
    {
/*
C        Define left edges for first pixel in this row.  The left
C        edges for the rest of the row will just be equal to the
C        right edge of the previous pixel.
C
C        as above, assume that the distance from the center to
C        the left edge is the same as the distance from the
C        center to the right edge
*/
//         xp = (1+0.5-xstep)/float(xstep) + 1.0;
//         yp = (j-ystart)/float(ystep) + 1.0;
         xp = (1.0-(float)xstep)/float(xstep) + 1.0;
         yp = float(j-ystart)/float(ystep) + 1.0;
         ix = (int)xp;
         iy = (int)yp;
         if (   (ix >= 1) && (ix < xxasize) 
             && (iy >= 1) && (iy < xyasize) ) 
           {
             nOffset =         ix +  iy * xxasize;
             lXINT[0] = poImgDistorXINT->lGetPixel(nOffset);
             lXINT[1] = poImgDistorXINT->lGetPixel(nOffset-1);
             lXINT[2] = poImgDistorXINT->lGetPixel(nOffset-xxasize);
             lXINT[3] = poImgDistorXINT->lGetPixel(nOffset-xxasize-1);
/*
////////////////////////////////
  lXINT[4] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix,   iy);
  lXINT[5] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix-1, iy);
  lXINT[6] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix,   iy-1);
  lXINT[7] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix-1, iy-1);
          if (lXINT[0] != lXINT[4])
            cout << "XINTa mismatch for get pix ix, iy\n";
          if (lXINT[1] != lXINT[5])
            cout << "XINTa mismatch for get pix ix-1, iy\n";
          if (lXINT[2] != lXINT[6])
            cout << "XINTa mismatch for get pix ix, iy-1\n";
          if (lXINT[3] != lXINT[7])
            cout << "XINTa mismatch for get pix ix-1, iy-1\n";
/////////////////////////////////
*/
             if (   (lIntBadFlag == lXINT[0])
                 || (lIntBadFlag == lXINT[1])
                 || (lIntBadFlag == lXINT[2])
                 || (lIntBadFlag == lXINT[3]) )
               {
                 x0 = fIntBadFlag;
               }
             else
               {
                 
               valu = (float)((yp-iy)*
                      ((xp-ix)* lXINT[0]             //x_int(ix+1,iy+1)
                       + (1+ix-xp)* lXINT[1]) +      // x_int(ix,iy+1)) +
                   (1+iy-yp)*
                   ((xp-ix)*        lXINT[2]         //x_int(ix+1,iy)
                    + (1+ix-xp)*     lXINT[3]))      //x_int(ix,iy))) 
                 * pscale;
               xp = (1-xstart)/xstep + 1.0;
               ix = (int)xp;
               valc = (float)((yp-iy)*
                   ((xp-ix)* lXINT[0]             //x_int(ix+1,iy+1)
                    + (1+ix-xp)*  lXINT[1]) +     //x_int(ix,iy+1)) +
                   (1+iy-yp)*
                       ((xp-ix)* lXINT[2]         //x_int(ix+1,iy)
                        + (1+ix-xp)* lXINT[3]))    //x_int(ix,iy)))
                 * pscale;
               x0 = 2.0 * valc  -  valu;
               }
           }
         else
           {
             x0 = fIntBadFlag;
           }

         // Now go across the row

         for (i = 1; i <= xs; i++)
           {
             // Calculate right edge

//             xp = (i+0.5-xstart)/float(xstep) + 1.0;
//             yp = (j-ystart)/float(ystep) + 1.0;
             xp = float(i-xstart)/float(xstep) + 1.0;
             yp = float(j-ystart)/float(ystep) + 1.0;
             ix = (int)xp;
             iy = (int)yp;

             if (   (ix >= 1) && (ix < xxasize) 
                 && (iy >= 1) && (iy < xyasize) ) 
               {
                 nOffset =         ix +  iy * xxasize;
                 lXINT[0] = poImgDistorXINT->lGetPixel(nOffset); 
                 lXINT[1] = poImgDistorXINT->lGetPixel(nOffset-1);
                 lXINT[2] = poImgDistorXINT->lGetPixel(nOffset-xxasize);
                 lXINT[3] = poImgDistorXINT->lGetPixel(nOffset-xxasize-1);
/*
////////////////////////////////
  lXINT[4] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix,   iy);
  lXINT[5] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix-1, iy);
  lXINT[6] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix,   iy-1);
  lXINT[7] = (LONG) (poImgDistorXINT->*poImgDistorXINT->prfGetPixel)(ix-1, iy-1);
          if (lXINT[0] != lXINT[4])
            cout << "XINTb mismatch for get pix ix, iy\n";
          if (lXINT[1] != lXINT[5])
            cout << "XINTb mismatch for get pix ix-1, iy\n";
          if (lXINT[2] != lXINT[6])
            cout << "XINTb mismatch for get pix ix, iy-1\n";
          if (lXINT[3] != lXINT[7])
            cout << "XINTb mismatch for get pix ix-1, iy-1\n";
/////////////////////////////////
*/
                 if (  (lIntBadFlag == lXINT[0])
                    || (lIntBadFlag == lXINT[1])
                    || (lIntBadFlag == lXINT[2])
                    || (lIntBadFlag == lXINT[3]) )
                   {
                     xn = fIntBadFlag;
                   }
                 else
                   {
                     xn = (float)((yp-iy)*  
                           ((xp-ix) * lXINT[0]   //x_int(ix+1,iy+1)
                            + (1+ix-xp) * lXINT[1]) + //x_int(ix,iy+1)) +
                      (1+iy-yp)*
                           ((xp-ix) * lXINT[2]      //x_int(ix+1,iy)
                            + (1+ix-xp) * lXINT[3])) //x_int(ix,iy)))
                       * pscale;
                   }
               }
            else
              {
                xn = fIntBadFlag;
              }
            

             // Calculate upper edge

             xp = (i-xstart)/float(xstep) + 1.0;
//             yp = ((float)j+0.5-(float)ystart)/float(ystep) + 1.0;
             yp = float(j-ystart)/float(ystep) + 1.0;
             ix = (int)xp;
             iy = (int)yp;

             if (   (ix >= 1) && (ix < yxasize) 
                 && (iy >= 1) && (iy < yyasize) ) 
               {
                 nOffset  = ix +  iy * yxasize;
                 lYINT[0] = poImgDistorYINT->lGetPixel(nOffset);
                 lYINT[1] = poImgDistorYINT->lGetPixel(nOffset-1);
                 lYINT[2] = poImgDistorYINT->lGetPixel(nOffset-xxasize);
                 lYINT[3] = poImgDistorYINT->lGetPixel(nOffset-xxasize-1);
/*
////////////////////////////////
  lYINT[4] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix,   iy);
  lYINT[5] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix-1, iy);
  lYINT[6] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix,   iy-1);
  lYINT[7] = (LONG) (poImgDistorYINT->*poImgDistorYINT->prfGetPixel)(ix-1, iy-1);
          if (lYINT[0] != lYINT[4])
            cout << "YINTb mismatch for get pix ix, iy\n";
          if (lYINT[1] != lYINT[5])
            cout << "YINTb mismatch for get pix ix-1, iy\n";
          if (lYINT[2] != lYINT[6])
            cout << "YINTb mismatch for get pix ix, iy-1\n";
          if (lYINT[3] != lYINT[7])
            cout << "YINTb mismatch for get pix ix-1, iy-1\n";
/////////////////////////////////
*/

                 if (   (lIntBadFlag == lYINT[0])
                     || (lIntBadFlag == lYINT[1])
                     || (lIntBadFlag == lYINT[2])
                     || (lIntBadFlag == lYINT[3]) )
                   {
                     yn = fIntBadFlag;
                   }
                 else
                   {
                     yn = (float)((yp-iy)*
                                  ((xp-ix)*  lYINT[0]
                                   + (1+ix-xp) * lYINT[1]) + 
                                  (1+iy-yp)*
                                  ((xp-ix) * lYINT[2]    
                                   + (1+ix-xp) * lYINT[3]))
                       * pscale;
                   }
               }
             else
               {
                 yn = fIntBadFlag;
               }

             if (   (x0 != fIntBadFlag) && (pfY0[i] != fIntBadFlag)
                 && (xn != fIntBadFlag) && (yn      != fIntBadFlag) )
               {
                 area = (xn - x0) * (yn - pfY0[i]);
                 if (0.0 >= area)
                   cout << "area <= 0 for input pixel i,j: " << i << ", " << j << endl;
                 if (0.0 < area)
                   {
                     ixs = nint(x0);
                     iys = nint(pfY0[i]);
                     ixe = nint(xn);
                     iye = nint(yn);

                     fPixIJ = (poImageIn->*poImageIn->prfGetPixel)(i-1, j-1);
//                     for (jj = iys; jj <= iye; jj++)
                     for (jj = iys-1; jj <= iye+1; jj++)
                       {
                         if ( (0 < jj) && (jj <= yasize) )
                           {
//                             for (ii = ixs; ii <= ixe; ii++)
                             for (ii = ixs-1; ii <= ixe+1; ii++)
                               {
                                 if ( (0 < ii) && (ii <= xasize) )
                                   {
                                     fPixOut
                    = (poImageOut->*poImageOut->prfGetPixel)(ii-1, jj-1);
                                     if (0.0 == fPixOut)
                                       (poImageOut->*poImageOut->prnSetPixel)(ii-1, jj-1, fPixIJ);
                                   } // if
                               } //enddo ii
                           }
                       } //enddo jj
                   }
               }

             // Set the left edge for the next pixel and the upper
             // edge for the next row

             x0      = xn;
             pfY0[i] = yn;
           }// end i loop
       } // end j loop

  delete [] pfY0;
  return (0);
}

int 
Ccalibrate::nFillHoles(Cimage *poImage)
{

  int i, j, ii, jj;
  int nIter;
  int nNumNewGood;
  int nNumGood, nNumBad;
  int nDim0, nDim1;
  int nNeighbors;
  short int iPixIJ;

  cout << "\n***   Filling holes in a mark image..." << endl;

  if ( (NULL == poImage) || !poImage->bIsAvailable())
    {
      return (-1);
    }

  poImage->nGetDimensions(&nDim0, &nDim1);
  nIter       = 0;
  nNumNewGood = 10000;
  nNeighbors  = 7;
  
  while ( (3 > nIter) || ((10 < nNumNewGood) && (10 >= nIter)) ) 
    {
      nIter++;
      nNumNewGood = 0;
      nNumGood    = 0;
      nNumBad     = 0;

      for (j = 0; j < nDim1; j++)
        {
          for (i = 0; i < nDim0; i++)
            {
              poImage->nGetPixel(i, j, &iPixIJ);
              if (0 == iPixIJ)
                {
                  nNumBad++;
                  nNumGood = 0;
                  // Average the pixels around this spot

                  float fSum = 0;
                  for (jj = max(0, j-1); jj <= min(j+1, nDim1-1); jj++)
                    {
                      for (ii = max(0, i-1); ii <= min(i+1, nDim0-1); ii++)
                        {
                          poImage->nGetPixel(ii, jj, &iPixIJ);
                          if (0 < iPixIJ)
                            {
                              nNumGood++;
                              fSum += (float)iPixIJ;
                            }
                        }
                    }
                  if (nNumGood >= nNeighbors)
                    {
                      fSum = fSum / (float)nNumGood;
                      if (nint(fSum/1000.0) * 1000 == fSum)
                        {
                          nNumNewGood++;
                          iPixIJ = (short int) fSum;
                          poImage->nSetPixel(i, j, iPixIJ);
                        }
                    }
                }
            }
        }                  
      if (3 < m_nVerbose)
        {
          cout << "\n================"
               << "\nIteration     " << nIter
               << "\nNew good      " << nNumNewGood
               << "\nNumber of bad " << nNumBad
               << endl;
        }
    }
  return (0);
}

int
Ccalibrate::nShiftMM(const float fXmm, const float fYmm, const float fNewPixSize)
{
  // Shift a XINT and YINT spatial distortion tables by a 
  // set amount of mm's
  // If fPixSize is not 0.0, then apply a new pixel size, too.

  // TODO: After applying a shift, then the XINV and YINV tables
  //       should be recalculated!

  int nStat;
  int i;
  if (3 < m_nVerbose)
    {
      cout << "\n***Shift distor tables ..." << endl;
      cout << "Input: Xmm, Ymm, pixel size: " << fXmm << ", " << fYmm 
           << ", " << fNewPixSize << endl;
    }

  if (   (NULL == m_poImgDistorXINT) || (!m_poImgDistorXINT->bIsAvailable())
      || (NULL == m_poImgDistorYINT) || (!m_poImgDistorYINT->bIsAvailable()) )
    {
      nStat = nReadDistor();
      if (0 != nStat)
        {
          cout << "ERROR Reading spatial distortion info!" << endl;
          return (nStat);
        }
    }

  int nDim;
  LONG lIntBadFlag;
  LONG lPixVal;
  LONG lPixShift;

  float pscale;
  pscale = m_tCorrectParams.pscale;
  if (3 < m_nVerbose)
    {
      cout << "\n      Pixel scale, 1/pscale: " << pscale  << ", "
           << 1.0 / pscale << "\nPixel size (mm): " <<  m_fpixel_size << endl;
    }
  lIntBadFlag = (long) m_tCorrectParams.badflag;
  lPixShift = (LONG) nint(fXmm 
                              / (m_tCorrectParams.pscale * m_fpixel_size));

  pscale = 0.0;
  if (0.0 < fNewPixSize)
    pscale = m_fpixel_size / fNewPixSize;

  nDim = m_poImgDistorXINT->nGetDimension();
  for (i = 0; i < nDim; i++)
    {
      lPixVal = m_poImgDistorXINT->lGetPixel(i);
      if (lPixVal != lIntBadFlag)
        {
          lPixVal += lPixShift;
          if (0.0 < pscale)
            lPixVal = nint( (float)lPixVal  * pscale);
          m_poImgDistorXINT->vSetPixel(i, lPixVal);
        }
    }

  lPixShift = (LONG) nint(fYmm 
                              / (m_tCorrectParams.pscale * m_fpixel_size));

  nDim = m_poImgDistorYINT->nGetDimension();
  for (i = 0; i < nDim; i++)
    {
      lPixVal = m_poImgDistorYINT->lGetPixel(i);
      if (lPixVal != lIntBadFlag)
        {
          lPixVal += lPixShift;
          if (0.0 < pscale)
            lPixVal = nint( (float)lPixVal  * pscale);
          m_poImgDistorYINT->vSetPixel(i, lPixVal);
        }
    }

  if (0.0 < fNewPixSize)
    m_fpixel_size = fNewPixSize;
  
  return (0);
}

int 
Ccalibrate::nStantonPixeltoMM(const float fPx0, const float fPx1,
                              float *pfXmm, float *pfYmm)
{
  float f0tmp, f1tmp, f0look, f1look, f0delta, f1delta;
  int   n0look, n1look;
  int   n0Interp[4], n1Interp[4];
  float f0Interp[4], f1Interp[4];
  int   nStat;
  
  // *In the Cimage class, elements are numbered
  //   from 0, not 1!

  f0look = ((fPx0+1.0) - (float) m_tCorrectParams.xint_start)
                        / (float) m_tCorrectParams.xint_step;
  f1look = ((fPx1+1.0) - (float) m_tCorrectParams.yint_start) 
                        / (float) m_tCorrectParams.yint_step;

  n0look = (int) f0look;
  n1look = (int) f1look;

  f0delta = f0look - (float) n0look;
  f1delta = f1look - (float) n1look;

  // Check that the interpolation points are all in bounds and not flagged bad
  
  nStat = m_poImgDistorXINT->nGetPixel( n0look, n1look, &n0Interp[0]);
  if ( (0 != nStat) || (n0Interp[0] == m_tCorrectParams.badflag) ) return (1);

  nStat = m_poImgDistorXINT->nGetPixel( n0look+1, n1look, &n0Interp[1]);
  if ( (0 != nStat) || (n0Interp[1] == m_tCorrectParams.badflag) ) return (1);

  nStat = m_poImgDistorXINT->nGetPixel( n0look, n1look+1, &n0Interp[2]);
  if ( (0 != nStat) || (n0Interp[2] == m_tCorrectParams.badflag) ) return (1);

  nStat = m_poImgDistorXINT->nGetPixel( n0look+1, n1look+1, &n0Interp[3]);
  if ( (0 != nStat) || (n0Interp[3] == m_tCorrectParams.badflag) ) return (1);

  nStat = m_poImgDistorYINT->nGetPixel( n0look, n1look, &n1Interp[0]);
  if ( (0 != nStat) || (n1Interp[0] == m_tCorrectParams.badflag) ) return (1);

  nStat = m_poImgDistorYINT->nGetPixel( n0look+1, n1look, &n1Interp[1]);
  if ( (0 != nStat) || (n1Interp[1] == m_tCorrectParams.badflag) ) return (1);

  nStat = m_poImgDistorYINT->nGetPixel( n0look, n1look+1, &n1Interp[2]);
  if ( (0 != nStat) || (n1Interp[2] == m_tCorrectParams.badflag) ) return (1);

  nStat = m_poImgDistorYINT->nGetPixel( n0look+1, n1look+1, &n1Interp[3]);
  if ( (0 != nStat) || (n1Interp[3] == m_tCorrectParams.badflag) ) return (1);

  // Now we have the coordinates of the 4 corners to interpolate from,
  //   so convert to float and do the interpolation

  for (int i = 0; i < 4; i++)
    {
      f0Interp[i] = (float) n0Interp[i];
      f1Interp[i] = (float) n1Interp[i];
    }

  f0tmp = f0Interp[0] +
            f0delta * (f0Interp[1] - f0Interp[0])  +
            f1delta * (f0Interp[2] - f0Interp[0])  +
            f0delta * f1delta * (f0Interp[0] + f0Interp[3] - 
                                 f0Interp[1] - f0Interp[2]);

  f1tmp = f1Interp[0] +
            f0delta * (f1Interp[1] - f1Interp[0])  +
            f1delta * (f1Interp[2] - f1Interp[0])  +
            f0delta * f1delta * (f1Interp[0] + f1Interp[3] - 
                                 f1Interp[1] - f1Interp[2]);

  // Apply fDirVec here

  // Correct for m_a2x2fDirVec next
      
//  *pfXmm = m_a2x2fDirVec[0][0] * f0tmp  +  m_a2x2fDirVec[1][0] * f1tmp;
//  *pfYmm = m_a2x2fDirVec[0][1] * f0tmp  +  m_a2x2fDirVec[1][1] * f1tmp;

    *pfXmm = f0tmp;
    *pfYmm = f1tmp;

  // Then scale to proper mm with appropriate beam center offset
  
  *pfXmm = *pfXmm * m_tCorrectParams.pscale * m_fpixel_size; // - m_fXCenter;
  *pfYmm = *pfYmm * m_tCorrectParams.pscale * m_fpixel_size; // - m_fYCenter;

  return (0);
}

int
Ccalibrate::nTransform(Cimage *poImageIn, Cimage *poImageOut, 
                       Cimage *poImageDark)
{
  // Transform the input image into the output image by using
  // the m_poImgTransform information which applies the non-uniformity
  // and spatial distortion corrections simultaneously.
  //
  // The transform image m_poImgTransform has the following properties:
  // There are pairs of values:  (input_pixel_offset, 65536 * scale_factor)
  // sorted in the order of which output_pixel they contribute to.  
  // A scale_factor=0 means the output_pixel is bad and/or has no input_pixels
  // that contribute to it.  If an output_pixel has contributions from more
  // than one input_pixel, then the contributions appear in the image from
  // lowest scale_factor to highest scale_factor.  Thus, if a scale_factor is
  // lower than the previous scale_factor, that means a new output_pixel is
  // being worked on.  A scale_factor of 1 is an indication to go to the
  // next output_pixel and the associated input_pixel does not contribute
  // to the output pixel.  In this way, the m_poImgTransform does not hold
  // any explicit information about which output_pixel the input_pixel(s)
  // contribute to.
  
  // If NULL != poImageDark, then *poImageIn is modified during the transform.

  // CAVEAT:  This ONLY works for (unsigned short int) images!

  static float s_fPixSize0 = 0.0;
  static float s_fPixSize1 = 0.0;

  int nStat = 0;
  if (NULL == m_poImgTransform)
    {
      cout << "Reading transform image (a) ... " << m_sTransformName << endl;
      m_poImgTransform = new Cimage(m_sTransformName);
      s_fPixSize0 = 0.0;
      s_fPixSize1 = 0.0;
    }
  //+jwp 2006-06-02
  else
    {
      //delete m_poImgTransform;
      //m_poImgTransform = NULL;
      cout << "NOT Reading transform image (b) ... " << m_sTransformName << endl;
      //m_poImgTransform = new Cimage(m_sTransformName);
      s_fPixSize0 = 0.0;
      s_fPixSize1 = 0.0;
    }
  //=jwp 2006-06-02
  if (!m_poImgTransform->bIsAvailable())
    {
      cout << "ERROR, transform image not available!" << endl;
      nStat = 3;
    }
  else
    {
      // Do the transformation

      int nDim0, nDim1;
      poImageIn->nGetDimensions(&nDim0, &nDim1);

      float fCalibScale = 1.0;
      int nScale;
      int nScalePrev     = -1;
      int nSaturated     = 0;
      int nOutPedestal   = 20;
      int nInPedestal    = 20;
      int nBadValue      = 0;
      int nInOffsetPixel = 0;
      int nNumNegative   = 0;

      
      int nNumTruncUp    = 0;
      int nNumTruncDown  = 0;
      int nNumBad        = 0;
      int nNumAvail      = 0;
      long long int nOutput = 0;
      long long int nInput  = 0;
      int nDark          = 0;
      int nOutDim0, nOutDim1;
      int a2nTOutDim[2];
      short int iInOffsetDelta;
      bool bSaturated = FALSE;
      bool bHaveSPD   = FALSE;
      bool bDark_RAXIS_compression = FALSE;
      bool bInp_RAXIS_compression  = FALSE;
      bool bOut_RAXIS_compression  = FALSE;

      float a4fSpdInfo[4];
      Cstring sPrefix;

      float fSumPix, fPix;
      float fBadValue;
      int   nNumPix, nOffset;
      unsigned short int uiPix;
      int   nPixVal;

#ifdef ANL_TIMER
      anl_reset_timer(0, "Ccal::nTransform, Tran");
      anl_start_timer(0);
#endif          

      // Get saturated value and number of input values for the transform

      nSaturated = (int) poImageIn->fGetSatValue();
//+-jwp 2007-Oct-12      nSaturated = min(65535, nSaturated);
      int nInp_CompressInfo0, nInp_CompressInfo1;
      int nOut_CompressInfo0, nOut_CompressInfo1;
      int nDark_CompressInfo0, nDark_CompressInfo1;
      bInp_RAXIS_compression = (1 < poImageIn->m_fCompressInfo[0]);
      nInp_CompressInfo0 = (int) poImageIn->m_fCompressInfo[0];
      nInp_CompressInfo1 = (int) poImageIn->m_fCompressInfo[1];
      
      bOut_RAXIS_compression = (1 < poImageOut->m_fCompressInfo[0]);
      nOut_CompressInfo0 = (int) poImageOut->m_fCompressInfo[0];
      nOut_CompressInfo1 = (int) poImageOut->m_fCompressInfo[1];

      if (bInp_RAXIS_compression != bOut_RAXIS_compression)
	{
	  cout << "ERROR pixel value compression scheme not same for IN and OUT!\n" << flush;
	}
//-jwp
      if (0 != m_poImgTransform->m_oHeader.nGetValue("CALIB_PEDESTAL", &nOutPedestal))
        (void) m_poImgTransform->m_oHeader.nGetValue("TRANSFORM_PEDESTAL", &nOutPedestal);

      (void) m_poImgTransform->m_oHeader.nGetValue("TRANSFORM_BADVALUE", &nBadValue);
      (void) m_poImgTransform->m_oHeader.nGetValue("TRANSFORM_SCALE", &fCalibScale);
      int nNumBadOut = 0;
      (void) m_poImgTransform->m_oHeader.nGetValue("TRANSFORM_BI_OUTPUT", &nNumBadOut);

      if ( (0.0 >= s_fPixSize0) || (0.0 >= s_fPixSize1) )
        {
          // Get the transformed pixel size
          // Use static vars so this is done only once.

          (void) m_poImgTransform->m_oHeader.nGetValue(Cdetector::ms_sDetectorNames,
                                                       1, &sPrefix);
       
          s_fPixSize0 = 0.10;
          s_fPixSize1 = 0.10;
          nStat = m_poImgTransform->m_oHeader.nGetValue(sPrefix
                                 + Cspatial::ms_sSpatialDistortionInfo,
                                                        4, a4fSpdInfo);
          if (0 == nStat)
            {
              s_fPixSize0 = a4fSpdInfo[2];
              s_fPixSize1 = a4fSpdInfo[3];
            }

          // Try to read spd tables if they haven't been read already

          if (   (NULL == m_poImgDistorXINT) || (!m_poImgDistorXINT->bIsAvailable())
              || (NULL == m_poImgDistorYINT) || (!m_poImgDistorYINT->bIsAvailable()) )
            {
              if (bFileExists(m_sDistorName + ms_asExt[0]))
                  nStat = nReadDistor();
              else
                nStat = -1;
            }
          else
            {
              nStat = 0;
            }
          if (0 == nStat)
            bHaveSPD = TRUE;
          nStat = 0;
        }

      poImageOut->nGetDimensions(&nOutDim0, &nOutDim1);
      cout << "Expected output size from Transform header: "
           << m_a2nTOutDim[0] << " by " << m_a2nTOutDim[1]
           << "\n  Actual output size from output    header: "
           << nOutDim0 << " by " << nOutDim1 << endl;
      if (   (nOutDim0 != m_a2nTOutDim[0])
          || (nOutDim1 != m_a2nTOutDim[1]) )
        cout << "ERROR in nTransform, mismatch of output image dimensions!\n"
             << flush;

      cout << "Transform scale factor: " << fCalibScale << endl;

      if (NULL == poImageDark)
        {
          nInPedestal = 0; // Or read from the input image?
	  cout << "WARNING! poImageDark is NULL!\n" << endl;
        }
      else
        {
          // Subtract dark image, keeping saturated pixels saturated
          // TODO: Special problem with pedestal!  Should it be added
          //       back in during dark subtraction?  Then removed later?
	  
#ifdef ANL_TIMER
          anl_reset_timer(1, "Ccal::nTransform, Dark");
          anl_start_timer(1);
#endif          
          nNumAvail = poImageDark->nGetDimension();
	  bDark_RAXIS_compression = (1 < poImageDark->m_fCompressInfo[0]);
	  nDark_CompressInfo0 = (int) poImageDark->m_fCompressInfo[0];
	  nDark_CompressInfo1 = (int) poImageDark->m_fCompressInfo[1];
          poImageIn->nSetNextPixel(0, 0);
          poImageDark->nSetNextPixel(0, 0);
          nNumTruncUp   = 0;
          nNumTruncDown = 0;
          nNumBad       = 0;
	  cout << "INFO: nInPedestal used in ::nTransform\n"
               << "      for-dark-subtraction (added)\n"
               << "      AND later before applying transform (subtracted): " << nInPedestal << endl;
	  cout << "INFO: nOutPedestal used in ::nTransform AFTER applying transform (added): " << nOutPedestal << endl;
          while (0 < nNumAvail)
            {
              nOutput =  (int)poImageIn->uiGetNextPixelNoInc();
	      if (bInp_RAXIS_compression && (32767 < nOutput))
		nOutput = ((short int)nOutput + nInp_CompressInfo1) * nInp_CompressInfo0;
              nDark   =  (int)poImageDark->uiGetNextPixel();
	      if (bDark_RAXIS_compression && (32767 < nDark))
		nDark = ((short int) nDark + nDark_CompressInfo1) * nDark_CompressInfo0;		  
              if (nSaturated > nOutput)
                {
                  nOutput = nOutput - nDark + nInPedestal;
                  if (0 > nOutput)
                    {
                      nOutput = 0;
                      nNumTruncUp++;
                    }
                  else if (nSaturated <= nOutput)
                    {
                      nOutput = nSaturated;
                      nNumTruncDown++;
                    }
                }
	      if (bInp_RAXIS_compression && (32767 < nOutput)) // Compress back
		{
		  nOutput = (unsigned short int)((float)nOutput / (float)nInp_CompressInfo0) + nInp_CompressInfo1;
		}
              poImageIn->vSetNextPixel( (unsigned short int) nOutput);
              nNumAvail--;
            }
#ifdef ANL_TIMER
          anl_stop_timer(1);
#endif
          cout << "\nFor dark subtraction: "
               << "\n     Pedestal used for dark subtraction:            " << nInPedestal
               << "\n     Number sat. or truncated down to " << nSaturated
               << "\n                                                   : " << nNumTruncDown
               << "\n     Number truncated up to 0 after adding pedestal: " << nNumTruncUp
               << endl << endl;
        }

      // If spatial distortion info is in the header, then
      // convert beam position to corrected position

      nStat = poImageIn->m_oHeader.nGetValue(Cdetector::ms_sDetectorNames, 1,
                                            &sPrefix);
      if (0 == nStat)
        {
          if (0 != poImageIn->m_oHeader.nGetValue(sPrefix
                             + Cspatial::ms_sSpatialBeamPosn,
                                                  2, a4fSpdInfo))
            {
              nStat = poImageIn->m_oHeader.nGetValue(sPrefix 
                                 + Cspatial::ms_sSpatialDistortionInfo,
                                                     4, a4fSpdInfo);
              if (0 != nStat)
                {
                  // No beam position available anywhere, so use 1/2 image size

                  a4fSpdInfo[0] = (float) (nOutDim0 / 2);
                  a4fSpdInfo[1] = (float) (nOutDim1 / 2);
                }
            }

          // At this point we have a valid direct beam position in raw pixels
          // so convert to corrected pixels if we have spd tables,

          if (bHaveSPD)
            {
              float f0, f1;

              nStat = nPxToPx(a4fSpdInfo[0], a4fSpdInfo[1], &m_tCorrectParams, 
                              &f0, &f1);
//            cout << "result from nPxtoPx: " << nStat << ",\nin: " 
//                   << a4fSpdInfo[0] << ", " << a4fSpdInfo[1] << ", out: "
//                   << f0 << ", " << f1 << endl << flush;
              if (0 == nStat)
                {
                  a4fSpdInfo[0] = f0;
                  a4fSpdInfo[1] = f1;
                }
              else
                {
                  a4fSpdInfo[0] = m_a2fbeam_position[0];
                  a4fSpdInfo[1] = m_a2fbeam_position[1];
                }
            }
          else
            {
              a4fSpdInfo[0] = m_a2fbeam_position[0];
              a4fSpdInfo[1] = m_a2fbeam_position[1];
            }

          // Even if we cannot convert the beam center, 
          // change the spd and nonunf info.

          a4fSpdInfo[2] = s_fPixSize0;
          a4fSpdInfo[3] = s_fPixSize1;

          poImageOut->m_oHeader.nReplaceValue(sPrefix
                                        + Cspatial::ms_sSpatialBeamPosn,
                                              2, a4fSpdInfo);
          poImageOut->m_oHeader.nReplaceValue(sPrefix 
                                        + Cspatial::ms_sSpatialDistortionInfo,
                                              4, a4fSpdInfo);
                  
          poImageOut->m_oHeader.nReplaceValue(sPrefix 
                                        + Cspatial::ms_sSpatialDistortionType,
                                              Cspatial::ms_sSpatialTypeSimple);
          poImageOut->m_oHeader.nReplaceValue(sPrefix 
                                        + Cnonunf::ms_sNonunfType,
                                           Cnonunf::ms_sNonunfStateSimpleMask);
          poImageOut->m_oHeader.nReplaceValue(sPrefix 
                                        + Cnonunf::ms_sNonunfInfo,
                                           "$(NONUNF_MASK)");
        }
      int nNumWritten = 0;
      nNumAvail       = 0;
      nNumTruncUp     = 0;
      nNumTruncDown   = 0;
      nNumBad         = 0;
      int nInFile     = 0;
      nStat      = m_poImgTransform->m_oHeader.nGetValue("CALIB_USED", &nNumAvail);
      nInFile = nNumAvail;
      m_poImgTransform->nGetDimensions(&nDim0, &nDim1);
      if (nNumAvail > nDim1)
	{
	  cout << "WARNING! number available less than transform dimension, so reset!\n";
	  nNumAvail = nDim1;
	}

      int i;
      int nJ0, nJ1, nJV, nJF;
/*
      int nJ0Prev, nJ1Prev;
      int nTest;
      nJ0Prev = -1;
      nJ1Prev = -1;
      m_nDimTrf = m_poImgTransform->nGetDimension(0);
*/
      // If interpolation of bad output pixels is to be done, 
      //  then prepare for that

      if (0 < nNumBadOut)
        {
          if (0 == m_nNumBadInterpOut)
            {
              m_nNumBadInterpOut = nNumBadOut;
              m_poImgBadInterp   = new Cimage(4, m_nNumBadInterpOut, eImage_I4);
              (void) m_poImgBadInterp->nSetNextPixel(0, 0);
            }
          else 
            {
              nNumBadOut = 0;
            }
//        cout << "nNumBadOut, m_nNumBadInterpOut: "
//               << nNumBadOut << ", " << m_nNumBadInterpOut << endl << flush;
        }

      // Prime the algorithm

      fCalibScale = fCalibScale / 65536.0;

      poImageOut->nSetNextPixel(0, 0);
      m_poImgTransform->nSetNextPixel(0, 0);
      nOutput         =  0;
      nInOffsetPixel  =  0;
      nScalePrev      = -1; // Cannot be 0 first time through

      // Loop through all available input_pixel contributions

      while (0 < nNumAvail)
        {
          iInOffsetDelta  = m_poImgTransform->iGetNextPixel();
	  nScale          = (int)m_poImgTransform->uiGetNextPixel();

	  //+jwp 2006-06-02
	  if ( (0 >= iInOffsetDelta) && (65535 == nScale) )
	    {
	      //+jwp 2009-01-09, this is now very unlikely and probably can only 
              // occur around the edges of the active area of the detector or in the
              // seams of a multimodule detector.
              //-jwp 2009-01-09      
	      // This is impossible and probably indicates a change in the
	      // m_poImgTransform.m_Next_Pixel pointer by 1 byte instead
	      // of 2 bytes somewhere.  Generally, if you get here you 
	      // will have a core dump later on 
	      // when deleting the m_poImageDark

	      cout << "nNumAvail: " << nNumAvail << " at: ";
	      cout << nInFile - nNumAvail << endl;
	      cout << "AB: " << iInOffsetDelta << ", " << nScale << endl;
	      //+jwp 2007-12-05
	      // Change because the true bug was found elsewhere
	      //iInOffsetDelta = 0;
	      //nScale         = 65535;
	      //-jwp 2007-12-05
	    }
	  //-jwp 2006-06-02

/*
          if (2 < m_nDimTrf)
            nJ0 = (int) m_poImgTransform->uiGetNextPixel();
          if (3 < m_nDimTrf)
            nJ1 = (int) m_poImgTransform->uiGetNextPixel();
*/
          nInOffsetPixel += iInOffsetDelta;
          nNumAvail--;
          
          // A value of 0 for iInOffsetDelta (if nScale not 0 nor 1)
          // A value of -32768 for iInOffsetDelta (if nScale not 0 nor 1)
          // indicates next pixel pair has the actual nInOffsetPixel

//          if (   (0 == iInOffsetDelta) 
//              && (1 < nScale)
//              && (0 < nNumAvail) )
//          if (0 == iInOffsetDelta) // fewest
          if (-32768 == iInOffsetDelta) // fewest
            if (1 < nScale)        // next fewest
              if (0 < nNumAvail)   // really never happens
            {
              // This means the next "pixels" contain the actual nInOffsetPixel
              // value
	      if (0 != nScalePrev)
		{
		  // Write out error message if this is not a bad pixel
		  cout << "SPECIAL! nNumAvail: " << nNumAvail << ", " << nScalePrev << endl;
		}
              
              nInOffsetPixel  = 65536 * m_poImgTransform->iGetNextPixel();
              nInOffsetPixel += (int)m_poImgTransform->uiGetNextPixel();
/*
              if (2 < m_nDimTrf)
                nJ0 = (int) m_poImgTransform->uiGetNextPixel();
              if (3 < m_nDimTrf)
                nJ1 = (int) m_poImgTransform->uiGetNextPixel();
*/
              nNumAvail--;
            }

          if (0 == nScalePrev)
            {
              // Output pixel is flagged as bad

              nNumBad++;
              nOutput = nBadValue;  // nBadValue does not need compression in any case
              poImageOut->vSetNextPixel( (unsigned short int) nOutput);
              if (0 != nNumBadOut)
                {
                  // Do this only the 1st image through here 
                  // AND if interp required

                  // Store list of bad pixels in the output image

                  // debug line:
//                if (nNumBadOut == m_nNumBadInterpOut)
//                  cout << "FIRST bad interpout!\n" << flush;

                  nNumBadOut--;
                  m_poImgBadInterp->vSetNextPixel(
                          (LONG) (nNumWritten % nOutDim0)); // px0
                  m_poImgBadInterp->vSetNextPixel(
                          (LONG) (nNumWritten / nOutDim0)); // px1
                  m_poImgBadInterp->vSetNextPixel(
                          (LONG)0); // value
                  m_poImgBadInterp->vSetNextPixel(
                          (LONG)0); // flag?
                }
              nNumWritten++;
              nOutput    = 0;
              bSaturated = FALSE;
            }
          else if (nScale < nScalePrev)
            {
              // Write out previous output pixel, but round-up,
              // scale properly and add in the pedestal.
              
              //          nOutput = ((nOutput + 32768) >> 16) + nOutPedestal;
              //          nOutput = ((nOutput+32768) / 65536) + nOutPedestal;

              nOutput = int(((double)nOutput * fCalibScale) + 0.5) + nOutPedestal;
              if ( bSaturated || (nOutput >= nSaturated) )
                {
                  nNumTruncDown++;
                  nOutput = nSaturated;
                }
              else if (0 <= nOutput)
                {
                  // No nothing, this should be the vast majority of the pixels
                }
	      else if (-16384 > nOutput)
                {
                  // Occurs on overflows into the sign bit of nOutput

                  nNumTruncDown++;
                  nOutput = nSaturated;
                }
              else if (0 > nOutput)
                {
                  nNumTruncUp++;
                  nOutput = 0;
                }
	      if (bOut_RAXIS_compression && (32767 < nOutput)) // Compress back
		{
		  nOutput = (unsigned short int)((float)nOutput / (float)nOut_CompressInfo0) + nOut_CompressInfo1;
		}
	      poImageOut->vSetNextPixel( (unsigned short int) nOutput);
/*
//              if (5000 > nNumWritten || (0 == (nNumWritten % 1000)))
              if (   (nJ0Prev  != (nNumWritten % nOutDim0) + m_a2nTOutOrig[0])
                  || (nJ1Prev  != (nNumWritten / nOutDim0) + m_a2nTOutOrig[1]) )
                {
                  cout << "Ac " << nNumWritten % nOutDim0  + m_a2nTOutOrig[0]
                       << ", " 
                       << nNumWritten / nOutDim0 + m_a2nTOutOrig[1]
                       << "  Ex " << nJ0Prev << ", " << nJ1Prev << "  W,O: " 
                       << nNumWritten << "  " << nOutput << endl;
                }

              nJ0Prev    = nJ0;
              nJ1Prev    = nJ1;
*/
              nNumWritten++;
              nOutput    = 0;
              bSaturated = FALSE;
            }
          if (1 < nScale)
            {
              // This pixel contributes to currently being summed output pixel

//              nInput = (int) poImageIn->uiGetPixel(nInOffsetPixel) - nInPedestal;
              nInput = (int) poImageIn->uiGetPixel(nInOffsetPixel);
	      if (bInp_RAXIS_compression && (32767 < nInput))
		nInput = ((short int)nInput + nInp_CompressInfo1) * nInp_CompressInfo0;
              if (nInput >= nSaturated)
                {
                  bSaturated = TRUE;
                }
              else
                {
                  nInput -= nInPedestal;
                  if (0 > nInput)
		    {
		      nNumNegative++;
                      //cout << "B NEGATIVE where not expected: " 
		      //   << nInput << endl;
		      //nInput = 0;  +JWP 2009-01-13 leave negative
		    }
                }
              nOutput += (nScale * nInput);
            }
          nScalePrev = nScale;
        } // end while

      // If requested, apply any interpolation of bad output pixels
      // after all output pixels are available

      if (NULL != m_poImgBadInterp)
        {
	  //cout << "::nTransform, doing interpolation!\n" << flush;
	  //m_poImgBadInterp->nWrite("badinterp.img");

          nInterpolationFunction(nOutDim0, nOutDim1, nBadValue,
                                 m_nNumBadInterpOut, m_poImgBadInterp, 
                                 poImageOut);

          // Now that all interpolated values are calculated,
          // copy them into the output image
          int nPix;
          m_poImgBadInterp->nSetNextPixel(0, 0);          
          for (i = 0; i < m_nNumBadInterpOut; i++)
            {
              nJ0     = (int) m_poImgBadInterp->lGetNextPixel(); // px0
              nJ1     = (int) m_poImgBadInterp->lGetNextPixel(); // px1
//+JWP 2008-11-11
// Colin noticed that uiPix is not set, but nPix is, so ...
// a few lines below poImageOut->vSetPixel(nOffset, uiPix) is meaningless
//
              //uiPix   = (unsigned short int) m_poImgBadInterp->lGetNextPixel(); // value
              nPix    = m_poImgBadInterp->lGetNextPixel(); // value
              nJF     = (int) m_poImgBadInterp->lGetNextPixel(); // flag
              if (1 < nPix)
                {
                  nOffset = nJ0 + (nOutDim0 * nJ1);

                  if (bOut_RAXIS_compression && (32767 < nPix)) // Compress back
                     uiPix = (unsigned short int)((float)nPix / (float)nOut_CompressInfo0) + nOut_CompressInfo1;
		  else
		    uiPix = (unsigned short int) nPix;
                  (void)poImageOut->vSetPixel(nOffset, uiPix);
/*
		  //(void)poImageOut->nSetPixel(nJ0, nJ1, uiPix);
                  cout << "OKpix nOffset, px0, px1: " << nOffset << ", "
                       << nOffset % nOutDim0 << ", " << nOffset / nOutDim0
                       << "   :   " << nJ0 << ", " << nJ1 << ", " << nPix
                       << endl << flush;

                }
              else
                {
                  cout << "baduipix nOffset, px0, px1: " << nOffset << ", "
                       << nOffset % nOutDim0 << ", " << nOffset / nOutDim0
                       << "   :   " << uiPix
                       << endl << flush;
*/
                }
            }
        } // end do interpolation if

/*
      cout << "just written " << nNumWritten
           << "   " << nNumWritten % nOutDim0 << ", " << nNumWritten / nOutDim0
           << " last read " << nJ0 << ", " << nJ1 
           << "\nlast written " << nJ0Prev << ", " << nJ1Prev << endl;
*/
      if (nNumWritten < poImageOut->nGetDimension())
        {
          // One should NOT get here!

          cout << "\nWARNING in ::nTransform, output is lite: Dim, Out: "
               << poImageOut->nGetDimension() << "  " << nNumWritten << endl;

          // So fill out with zeroes.
	  
	  nOutput = 0;
	  nNumWritten = poImageOut->nGetDimension() - nNumWritten;
	  while (nNumWritten > 0)
	    {
	      // No need for compress test since nOutput == 0
	      poImageOut->vSetNextPixel( (unsigned short int) nOutput);
	      nNumWritten--;
	    }
        }
        
      cout << "\nFor spatial and nonunf correction: "
           << "\n     Pedestal added to output pixel: " << nOutPedestal     
           << "\n     Number input pixels < 0       : " << nNumNegative     
           << "\n     Number truncated down to 65535: " << nNumTruncDown
           << "\n     Number truncated up   to     0: " << nNumTruncUp
           << "\n     Number flagged as bad:          " << nNumBad 
           << "  (" << nNumBad / (nOutDim0 * nOutDim1) * 100 << "%)\n"
           << endl;
#ifdef ANL_TIMER
      anl_stop_timer(0);
      if (NULL != poImageDark) anl_print_timer(1);
      anl_print_timer(0);
#endif
    }
  return (nStat);
}

int
Ccalibrate::nPostTransform(void)
{
  // Create a posttransform scaling image to do some post-scaling 
  // on the transform information.  The image is call posttransform.img
  // and can be used a subsequent run of ::nInterpolate.

  int nStat = 0;

  // Get reference flood image which should have been previously written
  // to disk.  (In future, make want to call ::nCalcReference again?)

  if (NULL != m_poImgRefer)
    {
      delete m_poImgRefer;
      m_poImgRefer = NULL;
    }
  cout << "Reading reference image... " << m_sReferName << endl;
  m_poImgRefer = new Cimage(m_sReferName);
  if (!m_poImgRefer->bIsAvailable())
    {
      cout << "ERROR, reference image not available!" << endl;
      nStat = 3;
    }

  // Get corrected flood image.

  Cimage *poImgRefCor;
  poImgRefCor = new Cimage(m_sFloodName);
  if (!poImgRefCor->bIsAvailable())
    {
      cout << "ERROR, corrected flood image unavailable!" << endl;
      nStat = 3;
    }

  if (0 == nStat)
    {
      // Have two images available, compute scale factors and store in
      // a new image.

      float fScale;
      float fPxScale;
      float fRef      = 0.0;
      float fCor      = 0.0;
      float fMinScale = 0.8;
      int   i, j, nDim0, nDim1;
      int   n0, n1;
      int   nOff;

      poImgRefCor->nGetDimensions(&nDim0, &nDim1);
//      cout << "dimensions of poImgRefCor: " << nDim0 << ", " << nDim1 << endl;
      n0  = nDim0 / 2;
      n1  = nDim1 / 2;

      Cimage *poImgPost;
      poImgPost = new Cimage(nDim0, nDim1, eImage_realIEEE);

      float fOffsetRef = 0.0f;
      float fOffsetCor = 0.0f;
      float fCorPx;
      float fNumPx = 0.0f;

//TODO: Should not include BAD PIXELS in the scale factor calculation
//      Should use more pixels in scale factor calculation

      poImgRefCor->m_oHeader.nGetValue("CALIB_PEDESTAL", &fOffsetCor);
      cout << "Offset(pedestal) used: " << fOffsetCor << endl;
      nOff = min(nDim0, nDim1) / 5;
      for (j = n1-nOff; j < n1+nOff; j++)
        {
          for (i = n0-nOff; i < n0+nOff; i++)
            {
              fCorPx = (poImgRefCor->*poImgRefCor->prfGetPixel)(i, j)
                          - fOffsetCor;
              if (0.0f < fCorPx)
                {
                  fRef += ( (m_poImgRefer->*m_poImgRefer->prfGetPixel)(i, j)
                            - fOffsetRef);
                  fCor   += fCorPx;
                  fNumPx += 1.0f;
                }
            }
        }
      if ((float)fNumPx < (0.5f * 4.0f * (float)nOff))
        cout << "WARNING only " << fNumPx << " used for scaling.\n" << flush;
//      cout << "n0, n1, nOff, fRef = " << n0 << ", " << n1 << ", " << nOff
//           << ", " << fRef << endl;
      fScale = fCor / fRef;
      cout << "Scale factor for cor / ref: " << fScale << endl;

      int nNumBad   = 0;
      int nNumAvail = poImgRefCor->nGetDimension();
      poImgPost->nSetNextPixel(0, 0);
      for (j = 0; j < nDim1; j++)
	{
	  for (i = 0; i < nDim0; i++)
	    {
	      fRef = (m_poImgRefer->*m_poImgRefer->prfGetPixel)(i, j) - fOffsetRef;
	      fCor = (poImgRefCor->*poImgRefCor->prfGetPixel)(i, j) - fOffsetCor;
	      fPxScale = 0.0;
	      if (0.0 != fCor)
		{
		  fPxScale = fScale * fRef / fCor;
		}

	      // Scale should be very close to 1.0, if not make it 0
	      // because something is drastically wrong with the original correction

	      if ( (fMinScale > fPxScale) || ((1.0 / fMinScale) < fPxScale) )
		{
		  //            if ( (50 > nNumBad) || (0 == (nNumBad % 1000)) )
		  //              cout << "nN, PostSc: " << nNumBad << "  " << fPxScale << endl;
		  fPxScale = 0.0;
		}
	      if (0.0 == fPxScale) 
		nNumBad++;
	      poImgPost->vSetNextPixel(fPxScale);
	      nNumAvail--;
	    }
	}

      // Write it out

      cout << "Number of bad pixels: " << nNumBad << endl;
      nStat = poImgPost->nWrite(m_sPostName);
      delete poImgPost;
    }
  delete poImgRefCor;
  return (nStat);
}

int
Ccalibrate::nBuildTransform(const Cstring &rsFilename, const int nNumToRead,
                            const int nDim0, const int nDim1)
{
  // Build a transform image from the data in the binary scratch file

  int i;
  int nStat;
  int nCalib;
  int nMaxNumPixels;
  int nNumInOut = 0;
  float fCalib;
  float fCalibMax;
  Cstring sErrorNoSpace
             = "ERROR not enough space allocated for the transform image!";

#ifdef CORRECT_ONLY
  return (-1);
#else
  if (NULL == m_poImgTransform) 
    {
      cout << "ERROR, no transform image available in nBuildTransform!\n" << endl;
      nStat = -1;
      return (nStat);
    }

  float fSum        = 0.0;
  unsigned short int uiCalib;

  Creflnlist *poOutList;
  Crefln     *poOutRef;
  int        *pnOutIdx;
  float       fMaxScale      = 0.0;
  float       fMinScale      = 65535.0;
  int         nNumPartial    = 0;
  int         nNumTooBig     = 0;
  int         nOutput        = 0;
  int         nNumRead       = 0;
  int         nNumWritten    = 0;
  int         nNumAdded      = 0;
  int         nNumError      = 0;
  int         nNumBad        = 0;
  int         nNumOutRange   = 0;
  int         nInOffsetPrev  = 0;
  int         nOutOffsetPrev = 0;
  int         nInOffsetDelta = 0;
  int         nOutOffsetDelta= 0;
  unsigned short int uiScalePrev    = 0;

  // To begin with, allow up to 18 input pixels per output pixels

  int         nNumAllocated = 18; 
  int         nOutOffset;
  int         nInOffset;
  int         nFile, nBytes;
  int         nOutI, nOutJ;
  tagTransformRef tRef;

  // Copy detector info from the flood image to the transform image.

  Cimage_header *poHeader;
  poHeader = new Cimage_header(m_sFloodName);
  if (poHeader->bIsAvailable())
    {
      Cstring sDetName;
      Cstring sTemp;
      nStat = poHeader->nGetValue(D_K_DetectorNames, &sDetName);
      if (0 != nStat)
        sDetName = "CCD_";
      sTemp = sDetName + '*';

//      cout << "copy mask is: " << sTemp << endl;
      m_poImgTransform->m_oHeader.nCopyMask(*poHeader, sTemp);
      sTemp = "DETECTOR*";
      m_poImgTransform->m_oHeader.nCopyMask(*poHeader, sTemp);
      m_poImgTransform->m_oHeader.nReplaceValue(sDetName
                                                + D_K_SpatialDistortionType,
                                                D_K_SpatialTypeSimple);

      // Convert beam position in raw image pixels to the position
      // in the transformed image pixels
      
      float fBeam0, fBeam1; 

      nPxToPx(m_tDistorParams.x_beam, m_tDistorParams.y_beam,
              &m_tCorrectParams, &fBeam0, &fBeam1);
      sTemp = Cstring(fBeam0) + ' ' + Cstring(fBeam1);
      m_poImgTransform->m_oHeader.nReplaceValue(sDetName
                                                + D_K_SpatialBeamPosition,
                                                sTemp);
      sTemp += Cstring(' ') + Cstring(m_fpixel_size) + ' ' + Cstring(m_fpixel_size);
      m_poImgTransform->m_oHeader.nReplaceValue(sDetName
                                                + D_K_SpatialDistortionInfo,
                                                sTemp);
      nStat = 0;
    }
  else
    {
      cout << "WARNING, flood image header (from " << m_sFloodName
           << ") not available to update transform header." << endl;
    }
  delete poHeader;

  if ( (0 >= m_a2nTOutDim[0]) || (0 >= m_a2nTOutDim[1]) )
    {
      cout << "ERROR! Transform output dimensions are bogus!\n" << endl;
    }

  // nStat is a dummy argument in the next line

  (void) m_poImgTransform->nGetDimensions(&nStat, &nMaxNumPixels);
  m_poImgTransform->m_oHeader.nReplaceValue("OUTPUT_SIZE", 2, m_a2nTOutDim);

  cout << "Creating TRANSFORM for output size of " << m_a2nTOutDim[0] 
       << " by " << m_a2nTOutDim[1] << " pixels\n" << flush;

  // Prepare a bad pixel image of same size as output image

  if (NULL != m_poImgBadPix)
    delete m_poImgBadPix;
  m_poImgBadPix = new Cimage(m_a2nTOutDim[0], m_a2nTOutDim[1], eImage_I2);
  m_poImgBadPix->nSetNextPixel(0, 0);
  for (i = 0; i < m_poImgBadPix->nGetDimension(); i++)
    m_poImgBadPix->vSetNextPixel((short int) 0);

  // Try to read a bad pixel file

  if (0 == nReadBadPixelList(m_sTBadpixName))
    {
      if (3 < m_nVerbose)
        {
          m_poImgBadPix->nWrite("tbadpix.img");
        }
    }
  else
    cout << "No additional bad pixels marked for the output image.\n" << endl;

  m_poImgBadPix->nSetNextPixel(0,0);
  short int iBadPixFlag = 0;
  
  nStat      = 0;

  // Open the scratch file

  nFile = 2;
  nBytes = rsFilename.length();
  (void) dskbor(&nFile, rsFilename.string(), &nBytes, &nStat);
  if (0 != nStat)
    {
      cout << "     Error opening file " << rsFilename
           << " Error is: " << nStat << "\n";
      return (nStat);
    }

  m_nDimTrf = m_poImgTransform->nGetDimension(0);

  fCalibMax = 1.0;
  m_poImgTransform->m_oHeader.nGetValue("TRANSFORM_SCALE", &fCalibMax);
  cout << "Transform scale factor: " << fCalibMax << endl;

  fCalibMax = 65536.0 / fCalibMax;

  poOutList = new Creflnlist();
  poOutRef  = new Crefln(poOutList);
  pnOutIdx  = new int [nNumAllocated];

  nBytes = sizeof(tagTransformRef);

  int nMin0, nMin1, nMax0, nMax1, nODim0, nODim1;
  bool bInOutRange = FALSE;

  nMin0  = m_a2nTOutOrig[0];
  nMin1  = m_a2nTOutOrig[1];
  nMax0  = m_a2nTOutOrig[0] + m_a2nTOutDim[0] - 1;
  nMax1  = m_a2nTOutOrig[1] + m_a2nTOutDim[1] - 1;
  nODim0 = m_a2nTOutDim[0];
  nODim1 = m_a2nTOutDim[1];

  cout << "\n***Build output transform for origin, extents: " 
       << nMin0 << ", " << nMin1
       << ", " << nMax0 << ", " << nMax1 << endl;
  cout << "         Output dimensions: " << nODim0 << ", " << nODim1 << endl << flush;
  nOutput = nMin0 + (nMin1 * nDim0); 

  while ( (nOutput < (nDim0 * nDim1)) && (0 == nStat) )
    {
      // Read next pixel contribution to get 
      //    (nInOffset, nOutOffset, areas)
      // Compute scale factor from ratio of areas

      if (nNumRead < nNumToRead)
        {
          // Use this little 'do while' to skip output pixels that are
          // outside the tlimits 

          do 
            {
              (void) dskbr(&nFile, (char*)&tRef, &nBytes, &nStat);
              if (0 != nStat)
                {
                  cout << "ERROR in ::nBuildTransform in dskbr(), nStat = "
                       << nStat << endl;
                }
              nNumRead++;
              nInOffset  = tRef.nInOffset;
              nOutOffset = tRef.nOutOffset;
              fCalib     = tRef.fCalibNumer / tRef.fCalibDenom
                             * fCalibMax;
              nOutI = nOutOffset % nDim0;
              nOutJ = nOutOffset / nDim0;
            } while ((0 == nStat)
                     && (nNumRead < nNumToRead)
                     &&
                     (   (nOutI < nMin0)
                      || (nOutI > nMax0)
                      || (nOutJ < nMin1)
                      || (nOutJ > nMax1) ) );

//          cout << "I,J: " << nOutI << ", " << nOutJ << endl;
          if ( (0 != nStat) || (nNumRead >= nNumToRead) )
            {
              // Set out nOutOffset to max possible nOutput + 1, 
              // so we know we are done with inputs
              
              nOutOffset = nDim0 * nDim1;
              // Should be no problem, but just in case ....
              uiCalib = 0; fCalib = 0.0;
            }
          else
            {
              fMaxScale = max(fMaxScale, fCalib);
              fMinScale = min(fMinScale, fCalib);

              //  Test range so fCalib fits in unsigned short int
              
              nCalib = nint(fCalib);

              if ( (2 > nCalib) || (65535 < nCalib) )
                {
                  // If there are any out of range you are in trouble!
                  // However, these usually occur around the edges of the active area,
                  // so they are inconsequential in the grand scheme of things.

		  if (200 > nNumOutRange)
		    cout << "INFO: fCalib out of range for output pixel (" 
                         << nOutOffset % nDim0 << ", " << nOutOffset / nDim0
                         << "): "  << fCalib << "  int: " << nCalib << endl;
		  else if (200 == nNumOutRange)
		    cout << "WARNING: more than 200 out of range for output pixels\n" 
                         << "      further WARNINGS suppressed!\n" << endl;
                  nNumOutRange++;


                  // Reset the range to values we think do the least damage.

                  if (2 > nCalib)
                    nCalib = 2;
//+JWP 2007-12-03 (bug on next line was repaired!)
                  else if (65535 < nCalib) 
		    {
		      nCalib = 65534;
		  //nCalib = 65535;
		    }
//-JWP 2007-12-03
                }

              uiCalib = (unsigned short int) nCalib;

              if ( (7 < m_nVerbose) && (0 == (nNumRead % 300000)) )
                {
                  // Some debug stuff
                  cout << "nNR, In, Out, Calib: " << nNumRead << ", "
                       << nInOffset << ", " << nOutOffset << ", "
                       << uiCalib << endl;
                }
//              cout << "Nread, NoutO: " << nNumRead << ", " << nOutOffset
//                   << endl;
            }
        }
      else
        {
          // Set out nOutOffset to max possible nOutput + 1, 
          // so we know we are done with inputs
              
          nOutOffset = nDim0 * nDim1;
          // Should be no problem, but just in case ....
          uiCalib = 0; fCalib = 0.0;
        }
      if (nOutOffset < nOutput)
        {
          // This is an error!

          cout << "SEVERE ERROR: nOutOffset < nOutput: " 
               << nOutOffset << " < " << nOutput
               << nOutOffset % nDim0 << ", " << nOutOffset / nDim0 
               << "  :  "
               << nOutput % nDim0 << ", " << nOutput / nDim0 
               << "\n   ***Consider increasing nTableSize elsewhere." 
               << endl;
          nStat = 1;
        }
      else if (nOutOffset != nOutput)
        {
          // Need to output contributions to previous
          // output pixel

          // Also check if this output pixel is BAD?

          iBadPixFlag = m_poImgBadPix->iGetNextPixel();

          if (   (0 == poOutList->nGetNumReflns())
              || (0 != iBadPixFlag))
            {
              //  If output pixel is BAD, then come to here also
              //  in that case be sure to poOutList->vDeleteAll() if not empty

              if (0 < poOutList->nGetNumReflns())
                poOutList->vDeleteAll();

              // Place output pixel with no info in output image
              // Write out (0,0);

              uiScalePrev = 0;
              nNumBad++;
              if (nNumWritten <= nMaxNumPixels) 
                {
                  m_poImgTransform->nSetPixel(0, nNumWritten,
                                              (unsigned short int) 0);
                  m_poImgTransform->nSetPixel(1, nNumWritten,
                                              uiScalePrev);
/*
                  if (2 < m_nDimTrf)
                    m_poImgTransform->nSetPixel(2, nNumWritten,
                                              (unsigned short int) (nOutput % nDim0));
                  if (3 < m_nDimTrf)
                    m_poImgTransform->nSetPixel(3, nNumWritten,
                                              (unsigned short int) (nOutput / nDim0));
                  if (4 < m_nDimTrf)
                    m_poImgTransform->nSetPixel(4, nNumWritten,
                                              (unsigned short int) 0);
                  if (5 < m_nDimTrf)
                    m_poImgTransform->nSetPixel(5, nNumWritten,
                                              (unsigned short int) 0);
*/
                }
              else
                {
                  if (50 > nNumError)
                    cout << sErrorNoSpace << endl;
                  nNumError++;

                  //                  nStat = -2;
                }
              nNumWritten++;
            }
          else
            {
              // More than one input pixel contributes to the output pixel
              // Check if the sum of areas is close to 1.0

              fSum = 0.0;
              for (i = 0; i < poOutList->nGetNumReflns(); i++)
                {
                  fSum += poOutList->poGetRefln(i)->fGetSigmaI();
                }
              if (1.015 < fSum)
                {
                  if ( (6 < m_nVerbose) && (50 > nNumTooBig) )
                    cout << "WARNING: (fSum > 1.015) = " << fSum << endl;
                  iBadPixFlag = 1;
                  nNumTooBig++;
                }
              if (   (fSum < 0.985)
                  || (0 != iBadPixFlag) )
                {
                  // Reject this output reflection!

                  if ( (6 < m_nVerbose) && (0 == iBadPixFlag) )
                    {
                      cout << "Partial output pixel, so rejected!" 
                           << "  (" <<  nOutput % nDim0 << ", " 
                           << nOutput / nDim0 << ")\n";
                    }

                  uiScalePrev = 0;
                  nNumBad++;
                  if (0 == iBadPixFlag)  // Count only if partial
                    nNumPartial++;
                  if (nNumWritten <= nMaxNumPixels) 
                    {
                      m_poImgTransform->nSetPixel(0, nNumWritten,
                                                  (unsigned short int) 0);
//+JWP 2007-12-03
		      if (uiScalePrev >= 65534)
			{
			  cout << "WARNING, tell Jim X: uiScaleprev: " << uiScalePrev 
                               << " numwritten: " << nNumWritten << endl;

			  uiScalePrev = 65534;
			}
//-JWP 2007-12-03
                      m_poImgTransform->nSetPixel(1, nNumWritten,
                                                  uiScalePrev);
/*
                      if (2 < m_nDimTrf)
                        m_poImgTransform->nSetPixel(2, nNumWritten,
                                                  (unsigned short int)
                                                  (nOutput % nDim0));
                      if (3 < m_nDimTrf)
                        m_poImgTransform->nSetPixel(3, nNumWritten,
                                                  (unsigned short int)
                                                  (nOutput / nDim0));
                      if (4 < m_nDimTrf)
                        m_poImgTransform->nSetPixel(4, nNumWritten,
                                                  (unsigned short int) 0);
                      if (5 < m_nDimTrf)
                        m_poImgTransform->nSetPixel(5, nNumWritten,
                                                  (unsigned short int) 0);
*/
                    }
                  else
                    {
                      if (50 > nNumError)
                        cout << sErrorNoSpace << endl;
                      nNumError++;
//                      nStat = -2;
                    }
                  nNumWritten++;
                }
              else
                {
                  // Sort on the Scale 

                  if (nNumAllocated < poOutList->nGetNumReflns())
                    {
                      delete [] pnOutIdx;
                      nNumAllocated = 2 * poOutList->nGetNumReflns();
                      pnOutIdx = new int [nNumAllocated];
                    }
                  poOutList->vSort(eReflnField_int_type, 
                                   poOutList->m_nFI_nL, pnOutIdx);

                  // Write out the contributions (InDelta, fScale)

                  if (   (uiScalePrev
                          <= (unsigned short int)poOutList->poGetRefln(pnOutIdx[0])->nGetL())
                      && (0 != uiScalePrev) )
                    {
                      if (7 < m_nVerbose)
			{
			  cout << "WARNING, prev scale <= curr scale: "
                               << uiScalePrev << " < " 
                               << 
                             (unsigned short int)poOutList->poGetRefln(pnOutIdx[0])->nGetL()
                               << "  (" <<  nOutput % nDim0 << ", " 
                               << nOutput / nDim0 << ")\n";
                        }

                      // Write out a bogus pixel so above does not happen
                      nNumAdded++;
                      uiScalePrev = 1;
                      if (nNumWritten <= nMaxNumPixels) 
                        {
                          m_poImgTransform->nSetPixel(0, nNumWritten,
                                                      (unsigned short int) 0);
                          m_poImgTransform->nSetPixel(1, nNumWritten,
                                                      uiScalePrev);
/*
                          if (2 < m_nDimTrf)
                            m_poImgTransform->nSetPixel(2, nNumWritten,
                                                      (unsigned short int) 
                                                      (nOutput % nDim0));
                          if (3 < m_nDimTrf)
                            m_poImgTransform->nSetPixel(3, nNumWritten,
                                                      (unsigned short int)
                                                      (nOutput / nDim0));
                          if (4 < m_nDimTrf)
                            m_poImgTransform->nSetPixel(4, nNumWritten,
                                                      (unsigned short int) 0);
                          if (5< m_nDimTrf)
                            m_poImgTransform->nSetPixel(5, nNumWritten,
                                                      (unsigned short int) 0);
*/
                        }
                      else
                        {
                          if (50 > nNumError)
                            cout << sErrorNoSpace << endl;
                          nNumError++;
//                          nStat = -2;
                        }
                      nNumWritten++;
                    }
                  for (i = 0; i < poOutList->nGetNumReflns(); i++)
                    {
                      uiScalePrev 
                        = (unsigned short int)
                          poOutList->poGetRefln(pnOutIdx[i])->nGetL();
                      nInOffsetDelta 
                        = poOutList->poGetRefln(pnOutIdx[i])->nGetH() 
                          - nInOffsetPrev;
		      
                      bInOutRange = FALSE;

                      if (   (-32767 > nInOffsetDelta)
                          || ( 32767 < nInOffsetDelta) )
//                          || (    0 == nInOffsetDelta) )
                        {
                          // Input pixel offset delta is out of range
                          // so used a value of -32768 which is a flag,
                          // then  uiScalePrev
                          // Next output pixel holds the actual input Pixel
                          // NOT the offsetDelta!

                          if (   (100 > nNumInOut) 
                              || (0 == (nNumInOut % 5000) ) )
                            cout << "INFO In off delta is large:"
                                 << nInOffsetDelta << ", " 
                                 << poOutList->poGetRefln(pnOutIdx[i])->nGetH() 
                                 << ", " << nInOffsetPrev << "  :  ("
                                   // nInOffset not what you want here
                                 << nInOffset % nDim0 << ", " 
                                 << nInOffset / nDim0
                                 << ")  (" << nInOffsetPrev % nDim0 
                                 << ", " << nInOffsetPrev / nDim0 << ")"
                                 <<  endl;

//                          nInOffsetDelta = 0;
                          nInOffsetDelta = -32768;
                          bInOutRange = TRUE;
                          nNumInOut++;
                        }
                      if (nNumWritten <= nMaxNumPixels) 
                        {
//+JWP 2007-11-09			  
                          // Sometimes fortuitously an input pixel is entirely included in
                          // an output pixel (?), so the next if-block prevents a problem
                          // in the ::nTransform() later on.

			  if ( (-1 == nInOffsetDelta) && (65535 == uiScalePrev) )
			    {
			      uiScalePrev = 65534;
			      //cout << "INFO: Special uiScalePrev, please report to Jim!\n";
			    }
			  if (uiScalePrev == 65535)
			    {
			      cout << "ABBBuiScaleprev, offset: " << uiScalePrev << ", "
                                   << nInOffsetDelta
				   << " numwritten: " << nNumWritten << endl;

			      uiScalePrev = 65535;
			    }
//-JWP 2007-11-09			  
                          m_poImgTransform->nSetPixel(0, nNumWritten,
                                                      (unsigned short int) nInOffsetDelta);
                          m_poImgTransform->nSetPixel(1, nNumWritten,
                                                      uiScalePrev);
/*
                          if (2 < m_nDimTrf)
                            m_poImgTransform->nSetPixel(2, nNumWritten,
                                                      (unsigned short int)
                                                      (nOutput % nDim0));
                          if (3 < m_nDimTrf)
                            m_poImgTransform->nSetPixel(3, nNumWritten,
                                                      (unsigned short int)
                                                      (nOutput / nDim0));
                          if (4 < m_nDimTrf)
                            m_poImgTransform->nSetPixel(4, nNumWritten,
                                                      (unsigned short int)
                                                      (poOutList->poGetRefln(pnOutIdx[i])->nGetH() % nDim0));
                          if (5 < m_nDimTrf)
                            m_poImgTransform->nSetPixel(5, nNumWritten,
                                                      (unsigned short int) 
                                                      (poOutList->poGetRefln(pnOutIdx[i])->nGetH() / nDim0));
*/
                        }
                      else
                        {
                          if (50 > nNumError)
                            cout << sErrorNoSpace << endl;

                          nNumError++;
//                          nStat = -2;
                        }
                      nInOffsetPrev  = poOutList->poGetRefln(pnOutIdx[i])->nGetH();
                      nNumWritten++;
                      if (bInOutRange)
                        {
                          // Need to output the actual nInOffset for
                          // the just written pixel (held in nInOffsetPrev)

                          unsigned short int uiHi, uiLo;
                          if (nNumWritten <= nMaxNumPixels) 
                            {
                              uiHi = nInOffsetPrev / 65536;
                              uiLo = nInOffsetPrev % 65536;
                              m_poImgTransform->nSetPixel(0, nNumWritten, uiHi);
                              m_poImgTransform->nSetPixel(1, nNumWritten, uiLo);
/*
                              if (2 < m_nDimTrf)
                                m_poImgTransform->nSetPixel(2, nNumWritten,
                                                      (unsigned short int)
                                                      (nOutput % nDim0));
                              if (3 < m_nDimTrf)
                                m_poImgTransform->nSetPixel(3, nNumWritten,
                                                      (unsigned short int)
                                                      (nOutput / nDim0));
                              if (4 < m_nDimTrf)
                                m_poImgTransform->nSetPixel(4, nNumWritten,
                                                      (unsigned short int)
                                                      (poOutList->poGetRefln(pnOutIdx[i])->nGetH() % nDim0));
                              if (5 < m_nDimTrf)
                                m_poImgTransform->nSetPixel(5, nNumWritten,
                                                      (unsigned short int) 
                                                      (poOutList->poGetRefln(pnOutIdx[i])->nGetH() / nDim0));
*/
                            }
                          else
                            {
                              if (50 > nNumError)
                                cout << sErrorNoSpace << endl;
                              nNumError++;
//                              nStat = -2;
                            }
                          nNumWritten++;
                        }
                    }
                }
              poOutList->vDeleteAll();
            }
//          cout << "AnOut: " << nOutput % nDim0
//               << ", " << nOutput / nDim0 << "   "
//               << nOutOffset % nDim0 << ", " 
//               << nOutOffset / nDim0 << ", " << nNumRead
//               << endl;

          nOutput++;

          // At this point we have written out all 
          // info related to output pixel nOutput

          int nOI, nOJ;
          while (nOutput < nOutOffset)
            {
              // Fill in any missing output pixels between nOutput 
              // and the new nOutOffset.
              // But keep in bounds of transformed output image

              nOI = nOutput % nDim0;
              nOJ = nOutput / nDim0;
              if (   (nOI < nMin0)
                  || (nOI > nMax0)
                  || (nOJ < nMin1)
                  || (nOJ > nMax1) )
                {
                  // Out of bounds, so do nothing
                  // Since nNumWritten will be incremented below, 
                  // decrement here to keep count correct

                  nNumWritten--;
                }
              else if (nNumWritten <= nMaxNumPixels) 
                {
                  // Place output pixel with no info in output image
                  // Write out (0,0);

                  // Move uiScalePrev = 0 inside this loop, since if
                  // no filler or dummy pixels are written out, we
                  // need to maintain the true uiScalePrev.

                  iBadPixFlag = m_poImgBadPix->iGetNextPixel();
                  uiScalePrev = 0;
                  nNumBad++;
//                  cout << "BnOut: " << nOI
//                       << ", " << nOJ << endl;
                  m_poImgTransform->nSetPixel(0, nNumWritten,
                                              (unsigned short int) 0);
                  m_poImgTransform->nSetPixel(1, nNumWritten,
                                              uiScalePrev);
/*
                  if (2 < m_nDimTrf)
                    m_poImgTransform->nSetPixel(2, nNumWritten,
                                              (unsigned short int) (nOutput % nDim0));
                  if (3 < m_nDimTrf)
                    m_poImgTransform->nSetPixel(3, nNumWritten,
                                              (unsigned short int) (nOutput / nDim0));
                  if (4 < m_nDimTrf)
                    m_poImgTransform->nSetPixel(4, nNumWritten,
                                              (unsigned short int) 0);
                  if (5 < m_nDimTrf)
                    m_poImgTransform->nSetPixel(5, nNumWritten,
                                              (unsigned short int) 0);
*/
                }
              else
                {
                  if (50 > nNumError) 
                    cout << sErrorNoSpace << endl;
                  nNumError++;
//                  nStat = -2;
                }
              nNumWritten++;
              nOutput++;
            }
        }

      // At this point nOutput should = nOutOffset

      if (nOutput != nOutOffset)
        {
          cout << "WARNING! nOutput != nOutOffset!" << endl;
        }

      // Store contributions to output pixel

      poOutRef->vSetH(nInOffset);
      poOutRef->vSetK(nOutOffset);
      poOutRef->vSetL((int)uiCalib);
      poOutRef->vSetSigmaI(tRef.fCalibNumer);
      poOutList->nInsert(poOutRef);

    } // end of while

  // Fill out to end of array with zeroes

  // Write out a final "flush" pixel (may want to write 0,0 instead of 1,0?!)

  m_poImgTransform->nSetPixel(0, nNumWritten,
                              (unsigned short int) 1);
  m_poImgTransform->nSetPixel(1, nNumWritten,
                              (unsigned short int) 0);
  nNumWritten++;
  cout << "Total <pixels> used in transform file: " << nNumWritten << endl;
  m_poImgTransform->m_oHeader.nReplaceValue("CALIB_USED", nNumWritten);
  while (nNumWritten < nMaxNumPixels)
    {
      m_poImgTransform->nSetPixel(0, nNumWritten,
                                  (unsigned short int) 1);
      m_poImgTransform->nSetPixel(1, nNumWritten,
                                  (unsigned short int) 0);
/*
      if (2 < m_nDimTrf)
        m_poImgTransform->nSetPixel(2, nNumWritten,
                                  (unsigned short int) 0);
      if (3 < m_nDimTrf)
        m_poImgTransform->nSetPixel(3, nNumWritten,
                                  (unsigned short int) 0);
      if (4 < m_nDimTrf)
        m_poImgTransform->nSetPixel(4, nNumWritten,
                                  (unsigned short int)1);
      if (5 < m_nDimTrf)
        m_poImgTransform->nSetPixel(5, nNumWritten,
                                  (unsigned short int)1);
*/
      nNumWritten++;
    }
      
  cout << "MinScale, MaxScale: " << fMinScale 
       << ", " << fMaxScale
       << endl;

  cout << "Number written:          " << nNumWritten << endl;
  cout << "Number read:             " << nNumRead    << endl;
  cout << "Number < 0.985:          " << nNumPartial << endl;
  cout << "Number > 1.015:          " << nNumTooBig  << endl;
  cout << "Number dummies:          " << nNumAdded   << endl;
  cout << "Number scales outrange:  " << nNumOutRange << endl;
  cout << "Number offsets outrange: " << nNumInOut << endl;
  cout << "Number write errors:     " << nNumError << endl;
  cout << "Number output bads:      " << nNumBad << endl;
  cout << "Number of reflns in scratch file: "
       << nNumToRead << endl << flush;
  if (0 == nStat) 
    {
      m_poImgTransform->m_oHeader.nReplaceValue("Create_date", sGetDate());
      m_poImgTransform->m_oHeader.nReplaceValue("Create_time", sGetTime());
      m_poImgTransform->m_oHeader.nReplaceValue("TRANSFORM_BI_OUTPUT", nNumBad);

      // Add the TLIMITs that were used, even if they were not specified
      int a4nTemp[4];
      a4nTemp[0] = m_a2nTOutOrig[0];
      a4nTemp[1] = m_a2nTOutOrig[1];
      a4nTemp[2] = m_a2nTOutDim[0];
      a4nTemp[3] = m_a2nTOutDim[1];
      m_poImgTransform->m_oHeader.nReplaceValue("TRANSFORM_LIMITS", 4, a4nTemp);
      m_poImgTransform->vSetState(eImage_available_state);
      nStat = m_poImgTransform->nWrite(m_sTransformName);
    }

  delete poOutList;
  delete poOutRef;
  delete [] pnOutIdx;
  if (nOutput != (nDim0 * nDim1))
    cout << "WARNING!  Ending nOutput: " << nOutput << endl;

  // Close the scratch file

  (void) dskbcr(&nFile, &nBytes);  // nBytes is temp here

  delete m_poImgBadPix;
  m_poImgBadPix = NULL;

  if (3 < m_nVerbose)
    {
      cout << "INFO:L  end of ::nBuildTransform, nStat is " 
           << nStat << endl << flush;
    }
  cout << flush;
  return (nStat);
#endif
}

int
Ccalibrate::nPreTransform(const int nMode, int nFile,
                          Creflnlist *poReflnlistTable, int *pnNumWritten)
{
  // Since building the transform requires too much memory for very
  // large detector, use a scratch file to store the sorted preliminary
  // information.
  // nMode   == CAL_DUMPALL write out all reflns in poReflnlistTable
  //         != write out 1/2 of reflns in poReflnlistTable
  // nFile   the file number used in the dskb* routines
  // poReflnlistTable  hold input pixel, output pixel and calibration
  //                   numerator and denominator values
  // pnNumWritten      TOTAL number of reflections written to disk

  int i;
  int nStat;
  int nNumToWrite;
  int *pnIndex;
  int *pnDelFlag;
  int nFI_nOutOffset;
  Crefln *poRefln;

  tagTransformRef tRef;
  int nBytes;

  // Sort the reflection list and write out some or all of the reflections 
  // to a special binary scratch file.

//  cout << "About to new pnIndex in nPreTransform with "
//       << poReflnlistTable->nGetNumReflns() << " reflns.\n" << flush;
  pnIndex   = new int [poReflnlistTable->nGetNumReflns()];
  pnDelFlag = new int [poReflnlistTable->nGetNumReflns()];

  nFI_nOutOffset =  poReflnlistTable->nGetFieldIndex(Creflnlist::ms_snK);
//  cout << "About to sort in nPreTransform\n" << flush;
  poReflnlistTable->vSort(eReflnField_int_type, nFI_nOutOffset, pnIndex);
//  cout << "Done with sort in nPreTransform\n" << flush;

  nNumToWrite = poReflnlistTable->nGetNumReflns();
//  cout << "Number of reflections in poReflnListTable at beginning: "
//       << nNumToWrite << endl;

  for (i = 0; i < nNumToWrite; i++)
    pnDelFlag[i] = 0;

  if (CAL_DUMPALL != nMode)
//    nNumToWrite = 3 * nNumToWrite / 4;
    nNumToWrite = 1 * nNumToWrite / 2;
  //    nNumToWrite = 1 * nNumToWrite / 3;
    
  poRefln    = new Crefln(poReflnlistTable);

  nBytes = sizeof(tagTransformRef);

  int nNumWritten = 0;
  for (i = 0; i < nNumToWrite; i++)
    {
      poRefln          = poReflnlistTable->poGetRefln(pnIndex[i]);
      tRef.nInOffset   = poRefln->nGetH();
      tRef.nOutOffset  = poRefln->nGetK();
      tRef.fCalibNumer = poRefln->fGetSigmaI();
      tRef.fCalibDenom = poRefln->fGetIntensity();

      // Note: file nFile must be opened and closed by the calling routine!

      (void) dskbw(&nFile, (char *)&tRef, &nBytes, &nStat);
      if (0 != nStat)
        {
          cout << "ERROR writing scratch file in nPreTransform!\n" << flush;
        }

      nNumWritten++;

      // Mark this refln for deletion

      pnDelFlag[pnIndex[i]] = -1;
    }

  *pnNumWritten += nNumWritten;

//  cout << "\nTotal number written: " << *pnNumWritten << endl;

  // Delete the reflns that are marked for deletion

  nNumWritten = poReflnlistTable->nDelete(-1, pnDelFlag);
//  cout << "Number deleted: " << nNumWritten << endl << flush;

//  cout << "Done with deletion!\n" << flush;

  delete [] pnIndex;
  delete [] pnDelFlag;
  if (3 < m_nVerbose)
    {
      cout << "Number of reflections in poReflnListTable at end: "
           << poReflnlistTable->nGetNumReflns() << endl;
      cout << "Returning from nPreTransform, nStat = " << nStat << endl << flush;
    }
  else
    {
      cout << '.' << flush;
    }
  return (nStat);
}

int
Ccalibrate::nMakeT(void)
{
  // Special test routine

  float fCalibMax = 2.2627;

  cout << "Maximum Calib transform scale factor: " << fCalibMax << endl;

  fCalibMax = 65534.9 / fCalibMax;
  if (fCalibMax > 65536.0) fCalibMax = 65536.0;

  fCalibMax = 65536.0 / fCalibMax;

  // Create transform image 3% larger than inputs to allow for
  // bad output pixels that are not in the input information

  int nMaxContribs = 37005158;

//  cout << "About to new transform image in nInterpolate size, j: "
//       << nMaxContribs << '\n' << flush;

  // Shouldn't the next few lines be in ::nBuildTransform(...)?

  if (NULL != m_poImgTransform) delete m_poImgTransform;

  m_nDimTrf = 2;
  m_poImgTransform = new Cimage(m_nDimTrf, nMaxContribs, eImage_uI2);

  // Make sure scale factors will fit in 2-bytes (0 - 65535), but
  // in reality (2-65535) since 0 and 1 are used as special flags

  cout << "Transform scale factor: " << fCalibMax << endl;
  if ("" != sGetEnv("DTEXPOSE_SCALE"))
    {
      float fExposeScale = 1.0;
      fExposeScale = atof(sGetEnv("DTEXPOSE_SCALE").string());
      cout << "Expose scale factor: " << fExposeScale << endl;
      fCalibMax = fCalibMax * fExposeScale;
      cout << "NEW Transform scale factor: " << fCalibMax << endl;
    }
  cout << "REPLACING 2: " << fCalibMax << endl;
  m_poImgTransform->m_oHeader.nReplaceValue("TRANSFORM_SCALE", 
                                            fCalibMax, 9);

  int nStat;

//  cout << "About to call nBuildTransform...\n" << flush;

  Cstring sScrName = "TRANSFORM_SCR";
  int xasize, yasize;
  xasize = 4968;
  yasize = 4608;
//int nRefsWritten = 37005158; //??
  int nRefsWritten = 35243008; // 37005158
  cout << "Number of reflns in scratch file: " << nRefsWritten << endl;

  nStat = nBuildTransform(sScrName, nRefsWritten, xasize, yasize);

  // Delete the scratch file

//      (void) nFileDelete(sScrName);

  return (nStat);
}

int
Ccalibrate::nUpdateHeaderCorrectedImage(Cimage *poImageIn, Cimage *poImageOut)
{
  // If spatial distortion info is in the header, then
  // convert beam position to corrected position

  int nStat;
  Cstring sPrefix;
  float a4fSpdInfo[4];
  int nOutDim0, nOutDim1;

  if (!poImageIn->bIsAvailable())
    return (-1);
  if (!poImageOut->bIsAvailable())
    return (-2);

  nStat = poImageIn->m_oHeader.nGetValue(Cdetector::ms_sDetectorNames, 1,
                                         &sPrefix);
  
  a4fSpdInfo[0] = m_a2fbeam_position[0];
  a4fSpdInfo[1] = m_a2fbeam_position[1];
  a4fSpdInfo[2] = m_fpixel_size;
  a4fSpdInfo[3] = m_fpixel_size;

  poImageOut->m_oHeader.nReplaceValue(sPrefix
                                      + Cspatial::ms_sSpatialBeamPosn,
                                      2, a4fSpdInfo);
  poImageOut->m_oHeader.nReplaceValue(sPrefix 
                                      + Cspatial::ms_sSpatialDistortionInfo,
                                      4, a4fSpdInfo);
                  
  poImageOut->m_oHeader.nReplaceValue(sPrefix 
                                      + Cspatial::ms_sSpatialDistortionType,
                                      Cspatial::ms_sSpatialTypeSimple);
  poImageOut->m_oHeader.nReplaceValue(sPrefix 
                                      + Cnonunf::ms_sNonunfType,
                                      Cnonunf::ms_sNonunfStateSimpleMask);
  poImageOut->m_oHeader.nReplaceValue(sPrefix 
                                      + Cnonunf::ms_sNonunfInfo,
                                      "$(NONUNF_MASK)");
  return (nStat);
}
