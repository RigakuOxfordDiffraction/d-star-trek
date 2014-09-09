#ifndef DT_CCALIBRATE_H
#define DT_CCALIBRATE_H
//
// Copyright (c) 1998 Molecular Structure Corporation
//
// Ccalibrate.h        Initial author: J.W. Pflugrath           02-June-1998
//    This file is the header file for class Ccalibrate
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

//+Include files

#include "Dtrek.h"
#include "dtreksys.h"
#include "dtrekdefs.h"
#include "Cimage.h"
#include "Cimage_header.h"
#include "Creflnlist.h"
#include "Cnonunf.h"

#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
#include "CXprop.h"
#endif

#include "dtrekvec.h"

//+Definitions, constants, and initialization of static member variables

// Status codes for points

#define QBAD      0
#define QGOOD     1
#define QAWAY     2
#define QMOVE     4
#define QLOW      8
#define QREPEAT  16 
#define QINTERPOLATED  513
#define MASK_SPACING 1.0
#define DARKIMG   0
#define FLOODIMG  1
#define FLOOD_DARKIMG 2
#define CTR       300
#define MASK_SIZE 600
#define DO_DUMP   1
#define DO_UPDATE 1
#define DO_RADIAL 1
#define INTFLAG 100000000

#define CAL_DO_DARK      1
#define CAL_DO_NONUNF    2
#define CAL_DO_SPATIAL   4
#define CAL_DO_TRANSFORM 8
#define CAL_DUMPSOME     0
#define CAL_DUMPALL      1

typedef struct _tagPoints
{
  float xo [ MASK_SIZE][ MASK_SIZE];
  float yo [ MASK_SIZE][ MASK_SIZE];
  float xc [ MASK_SIZE][ MASK_SIZE];
  float yc [ MASK_SIZE][ MASK_SIZE];
  float xu [ MASK_SIZE][ MASK_SIZE];
  float yu [ MASK_SIZE][ MASK_SIZE];
  float x_diff [ MASK_SIZE][ MASK_SIZE];
  float y_diff [ MASK_SIZE][ MASK_SIZE];
  int   status [MASK_SIZE][MASK_SIZE];
  float inten  [MASK_SIZE][ MASK_SIZE];
  int   xb, xe, pts_x;
  int   yb, ye, pts_y;
} tagPoints;

/*
 *     o : observed position
 *     e : expected position
 *     u : unpinned position
 * _diff : differences between position predicted by unpinning,
 *            and the observed position.
 */

typedef struct _tagMask_parameters
{
  int center[2];
  int horizontal[2];
  int vertical[2];
  int search_center[2], resrad, resrad2, lim, dump, pkwid, work;
  int misses, back, radial;
  float dis, move, away, pkrat, mask_angle;
  float beam[2], maskd, source;
  float maskparams[2], mask_spacing, pixel_size;
} tagMask_parameters;

/*
* (cent_to_maxval)    away : max dis centroid from max val
* (background_pixels) back : number of pixels to use for background
* (beam_position)     beam : beam position
* (center_peak)       center : center position
* (peak_distance)     dis : minimum separation between spots
* (0)                 dump : dump intermediate results
* (min_peak)          lim : minimumpeak height
* (???)               maskd : mask to detector distance
* (mask_angle)        mask_angle : angle of mask
* (bad_peaks)         misses : max misses in a line
* (cent_to_pred)      move : max. calculated-average position
* (0)                 pkrat : fraction of max peak value to use as background
* (peak_size)         pkwid : peak size
* (qradial)           radial : radial/noradial flag for interpolation routine
* (search_radius)     resrad, resrad2 : radial searching limits
* (vertical_limits)   vertical : y searching limits
* (xtod_distance)     source : source to detector distance
* (horizontal_limits) horizontal : x searching limits
* (1)                 work : plot on the workstation
* (mask_angle, 1)     maskparams : user vals for mask angle, x/y ratio
*/

/*
 *    Data structure for calibration parameters in calibration
 *     routines.
 */
typedef struct _tagDistort_parameters
{
  float x_center;
  float y_center;
  float x_pt_center;
  float y_pt_center;
  float x_scale;
  float y_scale;
  float ratio;
  float ver_slope;
  float horz_slope;
  float a1;
  float a;
  float b;
  float c;
  float spacing;
  float x_beam;
  float y_beam;
  int   x_size;
  int   y_size;
} tagDistort_parameters;

typedef struct _tagCorrect_parameters
{
  float x_center;
  float y_center;
  float x_pt_center;
  float y_pt_center;
  float x_scale;
  float y_scale;
  float ratio;
  float ver_slope;
  float horz_slope;
  int   xint_step, yint_step;
  int   xinv_step, yinv_step;
  int   xint_start, yint_start;
  int   xinv_start, yinv_start;
  int   badflag;
  float pscale;
} tagCorrect_parameters;

typedef struct _tagTransform_parameters
{
  float fMaxScale;
  float fMinScale;
  float fCalibMax;
  int   nNumPartial;
  int   nOutput;
  int   nNumRead;
  int   nNumWritten;
  int   nNumAdded;
  int   nNumOutRange;
  int   nInOffsetPrev;
  int   nOutOffsetPrev;
  int   nInOffsetDelta;
  int   nOutOffsetDelta;
  unsigned short int   uiScalePrev;
} tagTransform_parameters;

typedef struct _tagTransformRef
{
  int   nInOffset;
  int   nOutOffset;
  float fCalibNumer;
  float fCalibDenom;
}  tagTransformRef;

//+Code begin

class Ccalibrate {

public:

  static Cstring ms_sDtcalibFiles;
  static Cstring ms_sDtcalibOptions;

private:

  int          m_nVerbose;
  int          m_a2nDim[2];
  int          m_a2nTOutDim[2];
  int          m_a2nTOutOrig[2];
  int          m_nNumBadInterpOut;
  int          m_nDimTrf;
  int          m_nCorrectionMask;
#if ((!defined(SSI_PC)) && (!defined(NO_X_WINDOWS)))
  CXprop      *m_poXprop;
#endif
  Cstring      m_sOutHeader;   // Output header filename;
  Cstring      m_sNonunfName;
  Cstring      m_sDistorName;
  Cstring      m_sMaskName;
  Cstring      m_sReferName;
  Cstring      m_sDarkName;
  Cstring      m_sFloodName;
  Cstring      m_sSaveName;
  Cstring      m_sBadpixName;
  Cstring      m_sTBadpixName;
  Cstring      m_sDumpName;
  Cstring      m_sMarkName;
  Cstring      m_sCalibOptions;
  Cstring      m_sCalibFiles;
  Cstring      m_sTransformName;
  Cstring      m_sPostName;
  Cstring      m_sMergeInfoName;
  Cstring      m_sMaskPeaksName;

  Cimage        *m_poImgTransform;
  Cimage        *m_poImgNonunf;
  Cimage        *m_poImgDistorXINT;
  Cimage        *m_poImgDistorXINV;
  Cimage        *m_poImgDistorYINT;
  Cimage        *m_poImgDistorYINV;
  Cimage        *m_poImgMask;
  Cimage        *m_poImgRefer;
  Cimage        *m_poImgDark;
  Cimage        *m_poImgFlood;
  Cimage        *m_poImgBadPix;
  Cimage        *m_poImgMark;
  Cimage        *m_poImgBadInterp;
  Cimage_header *m_poHeaderSave;
  Creflnlist    *m_poReflnlistMask;

  Cimage        *m_poImgMergeDistorXINT;
  Cimage        *m_poImgMergeDistorXINV;
  Cimage        *m_poImgMergeDistorYINT;
  Cimage        *m_poImgMergeDistorYINV;
  Cimage_header *m_poHeadMergeDistor;

  Cnonunf       *m_poNonunf;

  int          m_nflood_radial;
  int          m_nflood_interpolate;
  int          m_nflood_geometry;
  int          m_nflood_film;
  float        m_fdarksub_const;
  float        m_fdarksub_scale;
  float        m_a2fpixel_sd[2];
  int          m_a3nmax_pixel[3];
  int          m_a3nmin_pixel[3];
  int          m_nradial_distortion;
  int          m_nbad_peaks;
  float        m_fmask_angle;
  int          m_nmin_peak;
  float        m_fcent_to_maxval;
  float        m_fcent_to_pred;
  int          m_nbackground_pixels;
  int          m_npeak_size;
  float        m_fpeak_distance;
  int          m_a2ncenter_peak[2];
  float        m_fsearch_radius;
  int          m_a2nsearch_center[2];
  int          m_a2nvertical_limits[2];
  int          m_a2nhorizontal_limits[2];
  float        m_fxtod_distance;
  float        m_fmasktod_distance;
  float        m_a2fbeam_position[2];
  float        m_fmask_spacing;   // This means mask must be SQUARE
  float        m_fMaskPeaksDivisor;
  float        m_fpixel_size;
  int          m_ndump;
  int          m_nedges;

  int          m_a2nSeqNum[2];
  int          m_a2nNumModules[2];
  bool         m_bCreateTransform;

  static Cstring ms_asExt[5];

  tagCorrect_parameters  m_tCorrectParams;
  tagDistort_parameters  m_tDistorParams;
  tagMask_parameters     m_tMaskParams;
  tagPoints              m_tMaskPoints;

  float        m_fReferenceMax;
  float        m_fReferenceDistortMaxFactor;
  std::ofstream  *m_poOutDump;
//

public:

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

Ccalibrate (Cimage_header& oHeader);

~Ccalibrate ();

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

//

int nInitValues(const int nDim0 = 1024, const int nDim1 = 1024);
int nInitValues(Cimage_header& oHeader);
int nList(const int nFlag = 1);
int nUpdateHeader(Cimage_header *poHeader, const Cstring &sPre="");
inline void vSetVerboseLevel(const int nLev) { m_nVerbose = nLev; }
int nSetup(void);
int nDoCalibrate(const Cstring& sCalibOptions);

inline Cstring sGetOutHeader(void) { return (m_sOutHeader); }
int nNotifyDisplay(const Cstring sImageIn = "",
		   const Cstring sReflnlist = "", const Cstring sHeader = "");
int nGoDistor(void);
int nGoNonunf(void);
int nSubtract(Cimage *poImage1, Cimage *poImage2, const float fScale=1.0,
              const float fOffset = 0.0, const float fMin = 0.0,
	      Cimage *poImageOut=NULL);
 
int nGetMask(Cimage *poMask, tagMask_parameters tMaskParams);
int nCalcInterpTable(void);
int nFindMaskCenter(tagPoints *ptMaskPoints,
		    tagDistort_parameters *ptDistorParams,
		    tagMask_parameters *ptMaskParams,
		    float *pfXcent, float *pfYcent);
int nFindRadial(tagPoints *ptMaskPoints,tagDistort_parameters *ptDistorParams,
		tagMask_parameters *ptMaskParams);

int
nInterpolate(tagCorrect_parameters *ptCorrectParams,
	     Cimage *poImgDistorXINT, Cimage *poImgDistorYINT,
	     Cimage *poImageIn, Cimage *poImageOut);

int nWriteDistor(void);
int nReadDistor(void);
int nFlagBadPixels(Cimage *poImageIn,
		   const int a2nHorzLim[2],
		   const int a2nVertLim[2],
		   const int a2nCenter[2],
		   const int nRadius,
		   const int nMinVal, const int nMaxVal,
		   const float fSigma,
		   const int a2nTiles[2], int a5nNumBads[5]);


int nReadBadPixelList(Cstring& rsFilename);
int nCorrectImage(const int nCorrectMask, Cimage *poImageIn, Cimage *poImageOut);

 int nUpdateHeaderCorrectedImage(Cimage *poImageIn, Cimage *poImageOut);
int nShiftMM(const float fXmm, const float fYmm, const float fPixSize=0.0);
int nStantonPixeltoMM(const float fPx0, const float fPx1,
		      float *pfXmm, float *pfYmm);

int nTransform(Cimage *poImageIn, Cimage *poImageOut, Cimage *poImageDark);
int nPostTransform(void);
int nBuildTransform(const Cstring &rsFilename, const int nNumToRead, 
		    const int nDim0, const int nDim1);

int nMakeT(void);

 private: // Some of the above function will be made private

int nFindCenter(const int nSize,        // box size to find center in
                const float fMove,      // max distance point can move
                const int nMinLim,      // min peak value
                const float fAway,       // max dist centroid can be from peak
                Cimage *poMask,         // Image to look for peaks in
                float *pfX, float *pfY, // position to find center at
                float *pfIntensity);    // intensity of peak

int nSearchCenter(Cimage *poMask,
                  tagMask_parameters tMaskParams,
                  float *pfX, float *pfY, float *pfInten);

void vDumpString(const Cstring& rsString);
void vDumpPoint(const int i, const int j, const float fX,
		const float fY, const float fInten, const int nStatus);

int
nPolyFit(const float *pafX, const float *pafY, const float *pafW,
	 const int nNumPoints, const int nDegree, 
	 float *pafC, float *pfS2, float *pfR2);

int
nSineFit(const float *pfX, const float *pfY, float *pfW,
         const int nNumPoints, const int nDegree, 
         const int nNumIter, const float fCutoff,
         float *pfC, float *pfS2, float *pfR2);

void vCalcResiduals(const float *pfX, const float *pfY, 
		   const int nNumPoints, const int nDegree,
		   const float *pfC, float *pfW);
float fCalcAverage(const int nNumPoints, const float *pfW);
void vCalcWeights(const int nNumPoints, const float fAverage,
		  const float fCutoff, float *pfW);

int nCalcInterpolate(tagPoints *ptMaskPoints,
		     tagDistort_parameters *ptDistorParams,
		     tagMask_parameters *ptMaskParams,
		     tagCorrect_parameters *ptCorrectParams,
		     Cimage *poImgMask,
		     Cimage *poImgDistorXINT,
		     Cimage *poImgDistorXINV,
		     Cimage *poImgDistorYINT,
		     Cimage *poImgDistorYINV,
		     int nEdges);

void vCalcExpected(tagPoints             *ptMaskPoints,
		   tagDistort_parameters *ptDistorParams,
		   tagMask_parameters    *ptMaskParams,
		   tagCorrect_parameters *ptCorrectParams);

void
vCalcUnpinned(tagPoints             *ptMaskPoints,
	      tagDistort_parameters *ptDistorParams,
	      tagMask_parameters    *ptMaskParams);

void
vReverseUnpinPoint(tagDistort_parameters *ptDistorParams,
		   const float xu, const float yu, 
		   float *pfxo, float *pfyo);
void
vCalcMaskPosition(tagCorrect_parameters *ptCorrectParams,
		  const float x, const float y, float *pfXm, float *pfYm);

void
vCompleteArray(Cimage *poImageIn, Cimage *poImageOut,
	       const int xstart, const int ystart, const int xend,
	       const int yend, const float fFlag);

void vFitLine(const float *pfX, const float *pfY, const float *pfW,
	      const int nNumPoints, float *pfC);
void vFitPoly2(const float *pfX, const float *pfY, const float *pfW,
	       const int nNumPoints, float *pfC);
  
float fBSpline1 (const float G[4], const float T);

float fBSpline2 ( const float lt[4][4], const float x, const float y );

int nInBox(const float xn, const float yn,
	   const float x1, const float y1, const float x2, const float y2,
	   const float x3, const float y3, const float x4, const float y4);


int nGoInterp(void);
int nGoCorrect(const Cstring sImgNameIn, const Cstring sImgNameOut);
int nGoCorrectScan(const Cstring sImgNameIn, const Cstring sImgNameOut);
int nNormalize(Cimage *poImageIn, Cimage *poImageOut=NULL);
int nCalcReference(tagCorrect_parameters *ptCorrectParams,
		   tagDistort_parameters *ptDistorParams,
		   tagMask_parameters    *ptMaskParams);

int nPxToPx(const float xi, const float yi,
            tagCorrect_parameters *ptCorrectParams,
            float *pfXc, float *pfYc);
int nDistortReference(Cimage *poImageIn, Cimage *poImageOut);
int nDistortReference2(Cimage *poImageIn, Cimage *poImageOut);

int nRadialPxToPx(const float xi, const float yi, float *pfX0, float *pfYo);
int nScaleImages(void);
int nMakeMarks(void);
int nMergeTables(void);
int nMoveMark(tagCorrect_parameters *ptCorrectParams,
	      Cimage *poImgDistorXINT, Cimage *poImgDistorYINT,
	      Cimage *poImageIn, Cimage *poImageOut);
int nFillHoles(Cimage *poImage);
int nPreTransform(const int nMode, int nFile, Creflnlist *poReflnlist, int *pnNumWritten);

};  // end of class Ccalibrate

#endif   // DT_CCALIBRATE_H
