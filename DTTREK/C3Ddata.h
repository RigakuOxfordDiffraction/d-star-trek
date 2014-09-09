#ifndef DT_C3DDATA_H
#define DT_C3DDATA_H
//
// Copyright (c) 1995 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// C3Ddata.h        Initial author: J.W. Pflugrath           24-Mar-1995
//    This file is the header file for class C3Ddata
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
#include "Cimage.h"
#include "Crefln.h"
#include "CinterpMesh.h"

class Cprofit2DPerSlice;

//+Forward class declarations

class Cimage;   // Forward declaration of Cimage class.
class Cstat;
class C3Ddata;
class Cdetector;

//+Code begin

enum  e3Ddata_types
{
  e3Ddata_ushort,
  e3Ddata_int,     // signed int
  e3Ddata_float,
  e3Ddata_double
};

enum e3Ddata_specialtypes
{
    eSpecial_TrimEdges,
    eSpecial_DoNotIntegrateGarbage,
    eSpecial_PhiHKLLimits,
    eSpecial_UseMosaicityBoundsOnStrong
};

struct DTREK_EXPORT C3DdataInput {
  
  // Crystalographic data.

  float   m_fDetectorGain;          // INPUT  Detector gain.  Used in computation of sigma values.
  float   m_a2fShiftKa1[2];			// INPUT  Shift of Ka1
  float	  m_a2fShiftKa2[2];			// INPUT  Shift of Ka2
  float   m_a2x2fEllipsoidAIn[2][2];// INPUT Mask ellipsoid supplied by user. (2x2 A matrix).
  float   m_a2fEllipsoidbIn[2];     // INPUT Mask ellipsoid supplied by user. (b vector);
  float   m_fEllipsoidcIn;          // INPUT Mask ellipsoid supplied by user. (c constant);

  // User supplied settings.
  float   m_fUseBackgroundBorder;   // INPUT  Amount of pixels to use on edge of shoebox for background computation.  Very important for speed issues.
  float   m_fOnEdgeMaxPercent;      // INPUT  Maximum total integrated intensity that can lie on shoebox 0 or 1 edges.
  float   m_fMinPercentBackground;  // INPUT  Minimum percentage of pixels that can be in the background.
  float   m_fSigmaAboveBackground;  // INPUT  Sigma above background used during background searching.
  float   m_fSigmaBelowBackground;  // INPUT  Sigma below background used during background searching.
  int     m_nBackgroundBorder;      // INPUT  Pixels to use in boarder.
  float   m_fFattening;             // INPUT  How much to stretch a 2D peak profile.
  float   m_fFractionIntense;       // INPUT  Used to remove 'tails' from a peak for width purposes (does not affect integrated intensity or profile fitting).  Used to load m_nIntensePeakXXXXX variables.  
  float   m_fMosaicityLowestSpotIoverSig;          //INPUT   Lowest I/sig(I) allowed for a slice that contributes to mosaicity.
  int     m_nMinPeakRad;            // INPUT  Minimum peak radius that will be used.  This is also 'smeared' by the ka1/ka2 characteristics.  
  int     m_nIntegrateWeakFlag;     // INPUT  0 == Integrate off-centroid slices only if strong.  1== Integrate all predicted mosaicity slices 2 == Integrate everything (including padding)
  bool    m_bFitToPlanar;           // INPUT  Should we fit to a planar background.
  bool    m_bPrint;                 // INPUT  Does a few extra calculations.
  bool    m_bSpotChase;             // INPUT  Do we chase the starting spot centroid?
  bool    m_bUsePerSliceBackground; // INPUT  Do we use a per-slice background, or an overall background?
  bool    m_bProfileFitWeak;		// INPUT  Determines if we profile fit weak spots.
  bool    m_bIncludeSat;            // INPUT  Do we include saturated reflections?
  float   m_fMinProfitPeakIntensity;// INPUT  Minimum intensity to contribute to profile fitting.
  int     m_a3nHKL[3];              // INPUT  Only needed if we are enqueueing spots into the profile fitter.


  void  vInit();
  C3DdataInput();
  void  vCopy(const C3DdataInput&);
};

class DTREK_EXPORT C3DdataOutput {
    friend class C3Ddata;
    friend class Cprofit2DPerSlice;
    friend class Cprofit2DSliceModel;

  // Output variables for various routines.
  float   m_fAvgBackgroundAvg;                      // OUTPUT Background average
  float   m_fAvgBackgroundSigma;                    // OUTPUT and its standard deviation
  float   m_a3fAvgBackground[3];                    // OUTPUT Background plane ([0]*delta0 + [1]*delta1 + [2])
  float   m_a3fCentroid[3];                         // OUTPUT Centroid of peak in 3D array
  float   m_fRotSigma;                              // OUTPUT Sigma for rotation centroid.  In units of slices.
  float   m_fPeakIntensity;                         // OUTPUT Peak intensity
  float   m_fPeakIntensitySigma;                    // OUTPUT and its standard deviation (also includes background variance).
  float   m_fPeakBackgroundSigma;                   // OUTPUT sigma for peak intensity due to background.
  float   m_a3fPeakSize[3];                         // OUTPUT Size of peak in pixels
  float   m_fPeakMinValue;                          // OUTPUT Min pixel value in peak
  float   m_fPeakMaxValue;                          // OUTPUT Max pixel value in peak
  float   m_fPeakSharpness;                         // OUTPUT Describes how sharp a peak is.
  float   m_a3fPeakMaxValue[3];                     // OUTPUT Localtion of Max pixel value.
  int     m_nPeakStart;                             // OUTPUT Starting slice with calculated data.
  int     m_nPeakEnd;                               // OUTPUT Ending slice with calculated data.
  int     m_nIntensePeakStart;                      // OUTPUT Starting slice of most intense region of shoebox.
  int     m_nIntensePeakEnd;                        // OUTPUT Ending slice of most intense region of shoebox.
  int     m_nCentroid;                              // OUTPUT Centroid slice.  
  int     m_a2x2nEllipsoidRange[2][2];              // OUTPUT Ellipsoids on all slices are bounded here (provides a tight bound for optimized loops).
  int     m_nSatCount;                              // OUTPUT Count of saturated pixels.
  float   m_a2x2fEllipsoidA[2][2];                  // OUTPUT Mask ellipsoid (2x2 A matrix).
  float   m_a2fEllipsoidb[2];                       // OUTPUT Mask ellipsoid (b vector).
  float   m_fEllipsoidc;                            // OUTPUT Mask ellipsoid (c constant);

  float   m_a3fAbsCentroid[3];                      // MEMORY Offset centroid.  Set by user during bootstrap call.
  float   m_fRotWidth;                              // MEMORY Rot width specified by user during bootstrap call.


public:

  float fGetPerPixelBackground(void)                    { return m_fAvgBackgroundAvg; };
  float fGetPerPixelBackgroundSigma(void)       { return m_fAvgBackgroundSigma; };
  float fGetPeakBackgroundSigma(void) { return m_fPeakBackgroundSigma; };
  float fGetPeakIntensity(void) { return (m_fPeakIntensity); }
  float fGetPeakIntensitySigma(void) { return (m_fPeakIntensitySigma); }
  float fGetPeakMaxValue(void) { return (m_fPeakMaxValue); }
  float fGetPeakMinValue(void) { return (m_fPeakMinValue); }
  float fGetRotationSigma(void) { return (m_fRotSigma); };
  float fGetPeakSharpness(void) { return (m_fPeakSharpness); };
  float fGetPeakMax(int nDim)   { return m_a3fPeakMaxValue[nDim]; };
  float* pfGetEllipsoidA(void) { return &m_a2x2fEllipsoidA[0][0]; };
  float* pfGetEllipsoidb(void) { return &m_a2fEllipsoidb[0]; };
  float* pfGetAvgBackground()  { return &m_a3fAvgBackground[0]; };
  float& fGetEllipsoidc(void)  { return m_fEllipsoidc; };
  int   nGetEllipsoidRange(int n0,int n1) { return m_a2x2nEllipsoidRange[n0][n1]; };
  
  void  vGetCentroid(float *pfCentroid);
  void  vGetPeakSize(float *pfSize);
  void  vGetPeakStartEnd(int *pnStart, int *pnEnd);
  void  vGetIntensePeakStartEnd(int *pnStart,int *pnEnd);
  void  vGetSatCount(int* pnSatCount) { *pnSatCount = m_nSatCount; };


  void vInit();
  void vCopy(const C3DdataOutput& oOther); 
  C3DdataOutput();
};

const  int    m_nScratchFloatSize = 19;
const  int    m_nScratchFloatIntensity = 0;                       //[0] Per slice intensity
const  int    m_nScratchFloatSigma = 1;                             //[1] Per slice sigma value.
const  int    m_nScratchFloatIntensityProfit = 2;   //[2] Per slice intensity (Profit)
const  int    m_nScratchFloatSigmaProfit = 3;               //[3] Per slice sigma value (Profit)
const  int    m_nScratchFloatCompete = 4;                   //[4][5] Per slice competing pixel
const  int    m_nScratchFloatBackground = 6;                //[6][7][8] Per slice background
const  int    m_nScratchFloatStatus = 9;                    //[9] Per slice status. (see status values above).
const  int    m_nScratchFloatStat = 10;                         //[10] Per slice status. 
const  int    m_nScratchFloatEllipsoidA = 11;       //[11][12][13][14]  Per slice ellipsoid A
const  int    m_nScratchFloatEllipsoidb = 15;       //[15][16] Per slice ellipsoid b.
const  int    m_nScratchFloatEllipsoidc = 17;       //[17] Per slice ellipsoid c.
const  int    m_nScratchFloatSigmaBackground = 18;  //[18] Per slice background sigma.


class DTREK_EXPORT C3Ddata:public C3DdataInput,public C3DdataOutput  {

public:

  e3Ddata_types  m_eType;           // Simple data type, 2-bytes or 4-bytes
  Cimage *m_poImage;                // Pointer to image so conversion methods available
  int     m_nExt[3];                // Extent of data in 3 directions
  int     m_nAlloc;                 // Size in floats of allocated data
  int     m_nFilled[3];             // Amount filled in pfData in 3 directions;
  int     m_nOrigOffset[3];         // Offset in pxls & images from image and scan origin
  float  *m_pfData;                 // Pointer to the float data
  double *m_pflData;                // Pointer to the double data
  unsigned short int *m_puiData;    // Pointer to the unsigned short integer data
  int    *m_pnData;                 // Pointer to the signed int (4-byte) data

  static int    ms_nBorderMask;
  static int    ms_nHighMask;
  static int    ms_nLowMask;
  static int    ms_nSatMask;
  static int    ms_nBadMask;
  static int    ms_nReservedMask;
  static int    ms_nEllipsoidMask;
  static int    ms_nCentroidPixelMask;
  static int    ms_nContiguousHighMask; // Used when building the ellipsoid.
  static int    ms_nMinPeakPixels;      // Used as a temporary when calculating the contiguous region.

  static int    ms_nBoundaryPeakPixels;

  static int    ms_nSpecialOptions;     // Special processing options mask.

  static int	ms_nPerSliceInMosaicity;
  static int    ms_nPerSliceCentroid;
  static int	ms_nPerSliceIntegrateWeak;
  static int    ms_nPerSlicePostRefineNotInStatisticalMosaicity;
  static int    ms_nPerSlicePostRefineNotInRmergeMosaicity;

  static Cprofit2DPerSlice* m_poProfit;

  private:

  friend class Cprofit2DPerSlice;
  friend class Cprofit2DSliceModel;

  static int    ms_nNumObjects;
  static int   *ms_pnMask;
  static float *ms_pfMask;
  static float *ms_pfPrint;
  static int    ms_nAllocSmooth;
  static int   *ms_apnScratchInt[2];
  static int    ms_nAllocScratchInt;
         float *m_apfScratchFloat[m_nScratchFloatSize];
  

  static int    *ms_pnFreeSize;
  static float **ms_ppfFreeFloat;
  static int     ms_nFreeCount;
  static int     ms_nSoftCount;
  static int     ms_nHardCount;

//

////////////////////////////////////////////////////////////////////////
//+Constructors, destructors and assignments

public:


C3Ddata ();                          // Construct an empty 3Ddata object

C3Ddata(const int nInExt0, const int nInExt1, const int nInExt2,
        const int nOrig0, const int nOrig1, const int nOrig2,
	const e3Ddata_types eType = e3Ddata_float);

C3Ddata(const C3Ddata& oOther);

~C3Ddata ();

inline void vSetSigmaAbove(const float fValue) { m_fSigmaAboveBackground = fValue; }

inline int operator > (const C3Ddata& oIn)
{
  return (m_nAlloc > oIn.m_nAlloc);
}


C3Ddata& operator+=(C3Ddata& o3DdataIn);  // WARNING:  This operator does not ADD values... it REPLACES them.

////////////////////////////////////////////////////////////////////////
//+Member Function prototypes
 public:

inline int nGetExt(const int nWhich)    { return (m_nExt[nWhich]);}
inline int nGetOffset(const int nWhich) { return (m_nOrigOffset[nWhich]);}
inline void vSetOffsets(const int nOff0, const int nOff1, const int nOff2)
                { m_nOrigOffset[0] = nOff0; m_nOrigOffset[1] = nOff1;
		  m_nOrigOffset[2] = nOff2; }
inline int* pnGetMask(void) { return ms_pnMask; };

void vGetShiftedCentroid(float a3fCent[3],double fIntensity,int nSlice,bool bSet,int a2nAdjCent[2]);
void vGetPerSliceEllipsoidRelObsCent(int nSlice,float* a2x2fEllipsoidA,float* a2fEllipsoidb,float* pfEllipsoidc,float* pfAbsObsCentroidAlternate = NULL);

bool  bIntersect(const C3Ddata& oDataOther,int a6nIntersection[6]);

inline void vSetExtents(const int nExt0, const int nExt1, const int nExt2)
         { m_nExt[0] = nExt0; m_nExt[1] = nExt1; m_nExt[2] = nExt2; }
inline void vEmpty(void) { m_nFilled[0] = 0; m_nFilled[1] = 0; m_nFilled[2] = 0;}
inline void vSetType(const e3Ddata_types eType)
                           { m_eType = eType; }
inline e3Ddata_types eGetType(void) { return (m_eType); }

// The inline function below do not check what the m_eType is, so consequently
// there could be problems!

inline float* pfGetPerSliceIntensity() { return m_apfScratchFloat[m_nScratchFloatIntensity]; };
inline float* pfGetPerSliceSigma() { return m_apfScratchFloat[m_nScratchFloatSigma]; };
inline float* pfGetPerSliceProfitIntensity() { return m_apfScratchFloat[m_nScratchFloatIntensityProfit]; };
inline float* pfGetPerSliceProfitSigma() { return m_apfScratchFloat[m_nScratchFloatSigmaProfit]; };
inline float* pfGetPerSliceBackgroundSigma() { return m_apfScratchFloat[m_nScratchFloatSigmaBackground]; };
inline int	  nGetPerSliceProfitFlag(int nImage) { return (int) m_apfScratchFloat[m_nScratchFloatStatus][nImage]; };

inline float* pfGetPerSliceBackground(int nSlice) { return m_apfScratchFloat[m_nScratchFloatBackground] + nSlice*3; };
inline float* pfGetPrintMask() { return ms_pfPrint; };

inline float& fGetValue(const int n0, const int n1, const int n2)
  { return (m_pfData[n2 * (m_nExt[0] * m_nExt[1]) + (n1 * m_nExt[0]) + n0]); }

inline int& nGetValue(const int n0, const int n1, const int n2)
  { return (m_pnData[n2 * (m_nExt[0] * m_nExt[1]) + (n1 * m_nExt[0]) + n0]); }

inline double& lfGetValue(const int n0, const int n1, const int n2)
  { return (m_pflData[n2 * (m_nExt[0] * m_nExt[1]) + (n1 * m_nExt[0]) + n0]); }

inline unsigned short int& uiGetValue(const int n0, const int n1, const int n2)
  { return (m_puiData[n2 * (m_nExt[0] * m_nExt[1]) + (n1 * m_nExt[0]) + n0]); }

inline void vSetValue(const float fValue, const int n0, const int n1, const int n2)
  { m_pfData[n2 * (m_nExt[0] * m_nExt[1]) + (n1 * m_nExt[0]) + n0] = fValue; }

inline void vSetValue(const int nValue, const int n0, const int n1, const int n2)
  { m_pnData[n2 * (m_nExt[0] * m_nExt[1]) + (n1 * m_nExt[0]) + n0] = nValue; }

inline void vSetValue(const double fValue, const int n0, const int n1, const int n2)
  { m_pflData[n2 * (m_nExt[0] * m_nExt[1]) + (n1 * m_nExt[0]) + n0] = fValue; }

inline void vSetValue(const unsigned short int uiValue, const int n0, const int n1, const int n2)
  { m_puiData[n2 * (m_nExt[0] * m_nExt[1]) + (n1 * m_nExt[0]) + n0] = uiValue; }

inline void vSetCompetingSpot(int nSlice,float fPix0,float fPix1) { m_apfScratchFloat[m_nScratchFloatCompete + 0][nSlice] = fPix0; m_apfScratchFloat[m_nScratchFloatCompete + 1][nSlice] = fPix1; };
inline void vGetCompetingSpot(int nSlice,float* pfPix0,float* pfPix1) { *pfPix0 = m_apfScratchFloat[m_nScratchFloatCompete + 0][nSlice]; *pfPix1 = m_apfScratchFloat[m_nScratchFloatCompete + 1][nSlice]; };

inline int  nGetLayersFilled(int nLayer = 2) { return (m_nFilled[nLayer]); }
inline void vSetFilled(int n0,int n1,int n2) { m_nFilled[0] = max(m_nFilled[0],n0); m_nFilled[1] = max(n1,m_nFilled[1]); m_nFilled[2] = max(n2,m_nFilled[2]); };
inline bool bIsFilled(int nIndex = 2) { return (m_nFilled[nIndex] >= m_nExt[nIndex]); }


// Book-keeping.

int  nFill2D(Cimage *poImage, const int nWhich = -1);
void vZero(void);
int  nList(const int nLayer);
int  nWrite(const Cstring& rsFilename, const int nDir=0, Crefln *poRefln=NULL);
int  nCopyFrom(const C3Ddata& o3DdataIn);
int  nCopyExtra(const C3Ddata& o3DdataIn); 
int  nInitValues(void);

int  nConvertToFloat(Cimage *poImageIn = NULL, C3Ddata *po3DOther = NULL);

void vFreeInit(void);
void vGetMemory(const e3Ddata_types eType= e3Ddata_float);
void vFreeMemory(double **ppdLocation);
void vCheckMemory();    // Conditional call to vGetMemory()

// Calculation routines.

int  nCalcMask(const float a3fAbsCentroid[3],const float fRotWidth);
int  nBackgroundAvgSD(int nSlice,float* pfBackground,float* pfAverage,float* pfDeviation,float* pfMaxValue,float* pfPeakCentroid,float* pfTotalCentroid);
int  nFindBestCentroidSlice(const float a3fCent[3],const float* pfPerSliceIntensity,int nMosPeakStartSearch,int nMosPeakEndSearch);
int  nCalcContiguous(int nSlice,int nFlagsIn,int nFlagsOut,int a2nStart[2]);
int  nCalcEllipsoid(int nSlice,int a2nSliceCent[2],float a2x2fEllipsoidA[2][2],float a2fEllipsoidb[2],float& fEllipsoidc);
int  nCalcEllipsoidBitmask(int nSlice,int a2nSliceCent[2],int nBitsToSet,float a2x2fEllipsoidA[2][2],float a2fEllipsoidb[2],float fEllipsoidc,int a2x2nSearch[2][2],float* pfEllipsoidalBackground = NULL);
int  nCheckFuzzySpotShape(int nSlice);
int  nLoadPrint(bool bPrintHKLBounds);
int  nPrintMask(bool bPrintFull = TRUE,int nMaskToPrint =  0xffffffff);

int  nCalcGetPeakInfo(float *pfCentroid,
		      float *pfInt,  float *pfSigmaI,
		      float *pfAvgOut = NULL, float *pfSDOut = NULL,
		      float *pfSize = NULL);

int  nCalcGetMax(int *pn0Max, int *pn1Max, int *pn2Max, float *pfValue);

int nAddRogue(int nRogueImage);
void vTermRogue();

int nCountZeros(long *plNumberOfZeroPixels, const float fBadValue = 0.0);
int nResetLowValues(const float fMinValue, const float fNewValue);

void vGet2DPeakAreaPixelsCount(int& nPeakAreaPixels, int& nEllipsPixels, int& nPeakBorderPixels);

private:

void vAllocScratch(int nAlgorithm);                     // Allocate free variables for algorithms.
void vFreeScratch();                                    // Free variables.

void vFreeDelete(void);

};  // end of class C3Ddata



const int C3Ddata_Print_Div_0 = 5;    // Print Regions:  Used when printing out spot sizes.
const int C3Ddata_Print_Div_1 = 5;

// These are used to represent "Regions" of the image.
// All fields must be doubles.
struct DTREK_EXPORT C3DdataDetectorArea 
{
    double a2x2fEllipsoidA[2][2];
    double a2fEllipsoidb[2];

    double fEllipsoidc;
    
    double fIntensity;
    double fSigma;
    
    double fBackground;
    double fBackgroundSigma;
    
    double fSize[3];
    double fShoebox[2];

    void vDup();
    void vZero();
    
    int  nReadWriteHeader(bool bRead,int nIndex,Cstring& sOut,Cstring& sName);
};

class DTREK_EXPORT C3DdataDetectorProfile 
{
public:
    C3DdataDetectorProfile(int nDim0 = 0,int nDim1 = 0);

public:
    C3DdataDetectorArea m_aoDetectorAreasPrint[C3Ddata_Print_Div_0][C3Ddata_Print_Div_1];
    C3DdataDetectorArea m_aoDetectorAreasVarPrint[C3Ddata_Print_Div_0][C3Ddata_Print_Div_1];

    int  nAddPeakArea(int nPix0,int nPix1,C3Ddata* poShoebox);
    int  nGetPeakArea(int nPix0,int nPix1,C3DdataDetectorArea& oArea,C3DdataDetectorArea* poAreaVariance = NULL,int nLevelToUse = -1);
    int  nGetPeakArea(C3DdataDetectorArea& oArea);
    int  nInitPeakAreas(int nDim0,int nDim1);
    int  nPrintPeakAreas(bool bComputeOnly,double fRotIncrement);
    double fGetAverageLevel() { return m_oDetectorAreas.fGetAverageLevel(); }

    int nInitValues(Cimage_header& oHeader);
    int nUpdateHeader(Cimage_header *poHeader);
    bool bIsInitialized() { return (m_bMeshInitialized); }
    bool bHasPrintInfo() { return (m_bHasPrintAreaData); }
  
private:
    int               m_nDim0;
    int               m_nDim1;
    CinterpMesh       m_oDetectorAreas;
    bool              m_bMeshInitialized;          // Is the CinterpMesh object initialized
    bool              m_bHasPrintAreaData;         // Are the printing tables initialized.
    int               m_nContrib;                  // # contributors to the last call to nGetPeakArea().
};


#endif   // DT_C3DDATA_H
