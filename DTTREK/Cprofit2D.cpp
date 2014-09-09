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
// Cprofit2D.cc     Initial author: T.J. Niemeyer       18-Sep-2000
//    Contains all routines for 2D profile fitting.
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

#include "Cprofit2D.h"
#include "Cspatial.h"
#include "dtsvd.h"
#include "CScanBitmap.h"
#include <string.h>
#include "Cstat.h"
#include "Csource.h"


#include <time.h>

int g_nDebug = 0;
//////////////////////////////////////
//////////////////////////////////////
/////  tProfitRec Structure //////////
//////////////////////////////////////
//////////////////////////////////////

int tProfitRec::ms_a3nDim[3] = {0,0,0};
int tProfitRec::ms_a2nCOG[2] = {0,0};
itr<int> tProfitRec::ms_anInUse;

tProfitRec::tProfitRec(int nDim[2]) {
  ms_a3nDim[0] = nDim[0];
  ms_a3nDim[1] = nDim[1];
  ms_a3nDim[2] = nDim[0]*nDim[1];
  ms_a2nCOG[0] = nDim[0]/2;
  ms_a2nCOG[1] = nDim[1]/2;
  m_pfX        = new double[ms_a3nDim[2]];
  m_pfVar      = new double[ms_a3nDim[2]];
  bHasX        = TRUE;
  bHasX2       = TRUE;
  vClear();
};
tProfitRec::tProfitRec() {
  m_pfX   = new double[ms_a3nDim[2]];
  m_pfVar = new double[ms_a3nDim[2]];
  bHasX   = TRUE;
  bHasX2  = TRUE;
  vClear();
};

tProfitRec::tProfitRec(double* pfX,double* pfVar,double fContrib) {
  m_pfX      = pfX;
  m_pfVar    = pfVar;
  m_fContrib = fContrib;
  bHasX      = FALSE;
  bHasX2     = FALSE;
};

tProfitRec::~tProfitRec() {
  if (bHasX)
    {
      delete [] m_pfX;
      bHasX = FALSE;
    }
  if (bHasX2)
    {
      delete [] m_pfVar;
      bHasX2 = FALSE;
    }
};

void tProfitRec::vClear()
{
  int nx;
  for (nx=0;nx<ms_a3nDim[2];nx++)
    {
      m_pfX[nx] = 0.0;
      m_pfVar[nx] = 0.0;
    }
  m_fContrib = 0.0;
  m_fNominalVar = 0.0;
};

int tProfitRec::nAdd(tProfitRec& oOther) {
  int nx;

  for (nx=0;nx<ms_a3nDim[2];nx++)
    {
      m_pfX[nx]   += oOther.m_pfX[nx];
      m_pfVar[nx] += oOther.m_pfVar[nx];
    }
  m_fContrib += oOther.m_fContrib;
  return 0;
};

int tProfitRec::nCalcCOG(double* pfCOG) {
  int     n0,n1;
  double  fSum;
  double *plf0;

  pfCOG[0] = 0.0;
  pfCOG[1] = 0.0;
  fSum     = 0.0;
  plf0     = m_pfX;

  // Figure out the COG.
  for (n1=0;n1<ms_a3nDim[1];n1++)
    {
      for (n0=0; n0<ms_a3nDim[0]; n0++, plf0++)
        {
          pfCOG[0] += (double)n0 * (*plf0);
          pfCOG[1] += (double)n1 * (*plf0);
          fSum     += (*plf0);
        }
    }
  pfCOG[0] /= fSum;
  pfCOG[1] /= fSum;
  return 0;    
};


int tProfitRec::nLoad(float* pfDataX,
                      float* pfDataVar,
                      float fNominalVar,
                      int nDim[3],
                      double* pfCOG)
{
    int n0,n1;
    int nn0,nn1;
    int nx,ny;
    int nGrid;
    double f0;
    double fCOG[2];
    int nCOG[2];
    double fSum,fSumAbs;
    int nSliceSize;
    float *pf0,*pf1;
    int *pnInUse;
    double a4fGrid[4];
    int    a2x4nGrid[2][4] = {{0,1,0,1},{ 0,0,1,1}};

    /*  The pixel shift matrix maps the input data centroid onto the profile 
        (*this) centroid.
        This is done so that ms_a2nCOG[*] is really the spot centroid 
        without a fractional center
    */
    
    vClear();
    
    // Find the COG of the data.

    nSliceSize = nDim[0]*nDim[1];
    fCOG[0] = 0.0;
    fCOG[1] = 0.0;
    fSum    = 0.0;
    fSumAbs = 0.0;
    pf0     = pfDataX;

    // Figure out the COG.
    for (n1=0; n1 < nDim[1]; n1++) 
      {
        for (n0=0; n0 < nDim[0]; n0++,pf0++) 
          {
            if (*pf0>g_fMaskOut) 
              {
                fCOG[0] += (double)n0 * (*pf0);
                fCOG[1] += (double)n1 * (*pf0);
                fSum    += *pf0;
                fSumAbs += (*pf0);
              } 
            else if (*pf0 < g_fMaskOut)
              {
                if ("" != sGetEnv("DTREK_VERBOSE"))
                  {
                    printf("WARNING:  Pixels less than '%f' encountered in 2D profile .dat file.\n",g_fMaskOut);
                    pf1 = pfDataX;
                    int nn1, nn0;
                    for (nn1=0; nn1 < nDim[1]; nn1++) 
                      {
                        for (nn0=0; nn0 < nDim[0]; nn0++, pf1++) 
                          {
                            printf("%.1f ",*pf1);
                          }
                        printf("\n");
                      }
                    printf("\n");
                  }
              }
          }
      }

    m_fIntensity = fSum;

    if ( (m_fIntensity == 0.0) || (fSumAbs<=1.0) )
      return 1;

    fCOG[0] /= fSumAbs;
    fCOG[1] /= fSumAbs;

    // The user might have specified a COG.
    
    if (pfCOG)
      {
        fCOG[0] = pfCOG[0];
        fCOG[1] = pfCOG[1];
      }

    nCOG[0] = (int) fCOG[0];
    nCOG[1] = (int) fCOG[1];
    
    // Compute the grid shifts.
    for (nx=0; nx<4; nx++)
      {
        a4fGrid[nx] = 1.0;
        for (ny=0;ny<2;ny++)
          {
            if (0 == a2x4nGrid[ny][nx])
              a4fGrid[nx] *= (fCOG[ny] - nCOG[ny]);
            else
              a4fGrid[nx] *= (1.0 - (fCOG[ny] - nCOG[ny]));
          }
      }

    // Now add the normalized profile @ the COG.

    nx = 0;
    pf0 = pfDataX;
    pf1 = pfDataVar;
    ms_anInUse.setsize(ms_a3nDim[0]*ms_a3nDim[1]);
    pnInUse = &ms_anInUse[0];
    memset(pnInUse,0,ms_a3nDim[0]*ms_a3nDim[1]*sizeof(int));

    for (n1=0; n1 < nDim[1]; n1++)
      {
        for (n0=0; n0 < nDim[0]; n0++,pf0++,pf1++)
          {
            nn0 = ms_a2nCOG[0] + n0 - nCOG[0] -1;
            nn1 = ms_a2nCOG[1] + n1 - nCOG[1] -1;
            if (   (nn0>=0) && (nn0<ms_a3nDim[0]-1)
                && (nn1>=0) && (nn1<ms_a3nDim[1]-1) )
              {                
                if (*pf0 > g_fMaskOut)
                  {
                    f0 = *pf0;

                    // Must center the profile at the centroid, 
                    // which is probably a floating point value.

                    for (nGrid=0;nGrid<4;nGrid++)
                      {
                        m_pfX[nx = nn0 + a2x4nGrid[0][nGrid] + (nn1 + a2x4nGrid[1][nGrid])* ms_a3nDim[0]] += f0/m_fIntensity*a4fGrid[nGrid];
                        m_pfVar[nx] += (*pf1)/(m_fIntensity*m_fIntensity)*a4fGrid[nGrid];
                        pnInUse[nx] = 1;
                      }
                  }               
              }
          }
      }
    
    m_fContrib    = 1.0;
    m_fNominalVar = fNominalVar / (m_fIntensity * m_fIntensity);
    return 0;
};


// Use this to save the data to a location in an image.
int tProfitRec::nSave(Cimage& oImage,
                      int a2nPointer[2],
                      double fScale)
{
    int    n0,n1;
    double *pfX;
    double fSum;
    
    if (0.0 > fScale)
      {
        // Normalize data so that we have a constant value.
        pfX  = m_pfX;
        fSum = 0.0;
        for (n1=0; n1 < ms_a3nDim[1]; n1++)
          {
            for (n0=0; n0 < ms_a3nDim[0]; n0++,pfX++)
              {
                fSum += (*pfX) / m_fContrib;
              }
          }
        if (fSum <= 0.0)
          fScale = 0.0;
        else
          fScale = (-fScale) / fSum;
      }
                                
    pfX = m_pfX;
    for (n1=0; n1 < ms_a3nDim[1]; n1++)
      {
        for (n0=0; n0 < ms_a3nDim[0]; n0++,pfX++)
          {
            (oImage.*oImage.prnSetPixel)(a2nPointer[0]+n0,a2nPointer[1]+n1,fScale*max(0.0,(*pfX))/m_fContrib);
          }
      }
    return 0;    
}


// Try to merge the two profiles using a weighting scheme.
int tProfitRec::nScale(tProfitRec& oRaw,
                       double& fIntensity,
                       double& fVariance,
                       double& fCorr)
{
    int n0,n1;
    int nx;
    double f0,f1,f2;
    double fK;
    double* pfX,*pfVarRaw,*pfXRaw;
    double fSumDenom,fSumNumer;
    double fTotalVariance;
    bool bDebug = FALSE;
    double fNorm;
    double fMinToUseInProfit;

    // First calculate fK, the scale factor.

    pfX       = m_pfX;
    pfVarRaw  = oRaw.m_pfVar;
    pfXRaw    = oRaw.m_pfX;
    fSumDenom = 0.0;
    fSumNumer = 0.0;
    double fSumVarDenom = 0.0;

    // Find the maximum value of the profile.
    // Only use those values in the profile
    // that are somewhat above this maximum.

    fMinToUseInProfit = 0.0;
    fNorm = 0.0;
    for (nx = ms_a3nDim[0]*ms_a3nDim[1],n0 = 0; n0 < nx; n0++)
      {
        fNorm += pfX[n0];
      }
    fMinToUseInProfit *= 0.1;
    fNorm = 1.0 / fNorm;

    for (nx=0,n1=0;n1<ms_a3nDim[1];n1++)
      {
        for (n0=0;n0<ms_a3nDim[0];n0++,nx++,pfX++,pfVarRaw++,pfXRaw++)
          {
            double fRaw;

            // Alternate Weighting scheme: f0 = 1.0;
            // Alternate Weighting scheme: f0 = *pfVarRaw
            if ((*pfX) > fMinToUseInProfit)
              {

//TJN 8-Apr-2005 make next statement always TRUE
		/*

                if ((*pfXRaw) > 0.0)
                  fRaw = (*pfXRaw);
                else
                  fRaw = 0.0;
		*/
                  fRaw = (*pfXRaw);
//-TJN
                                
                f0 = 1.0;
                //f0 = *pfVarRaw;
                //if (f0 == 0.0)
                //      f0 = oRaw.m_fNominalVar;

                fSumNumer    += (*pfX) * (fRaw) / (f0); 
                fSumDenom    += (*pfX) * (*pfX) / (f0); 
                fSumVarDenom += (*pfX) * (*pfX) / ((*pfVarRaw)?(*pfVarRaw):oRaw.m_fNominalVar);
              }
          }
      }
    if (fSumDenom==0.0)
      {
        fIntensity = 0.0;
        fVariance  = 1e20;
        return 1;
      }

    // fK should be the profile fitted intensity.

    fK =  (fSumNumer / fSumDenom) / fNorm;

    // Next, calculate the total variance and correlation.

    double fCount1   = 0.0;
    double fCount2   = 0.0;
    double fCount1_2 = 0.0;
    double fCount2_2 = 0.0;
    double fCount12  = 0.0;
    int    nCount    = 0;

    fTotalVariance   = 0.0;
    pfX = m_pfX;
    pfVarRaw = oRaw.m_pfVar;
    pfXRaw = oRaw.m_pfX;
    for (nx=0,n1=0; n1 < ms_a3nDim[1]; n1++)
      {
        for (n0=0; n0 < ms_a3nDim[0]; n0++,nx++,pfX++,pfVarRaw++,pfXRaw++)
          {
            // Use all the pixels in the profile to estimate variance.
            if ( (*pfX) > fMinToUseInProfit)
              {
                double fRaw;
                double fProfit;
                fProfit = (*pfX)*fNorm*fK;
                if ( (*pfXRaw) > 0.0) 
                  fRaw = (*pfXRaw);
                else
                  fRaw = 0.0;
                f0 = fRaw - fProfit;
                fTotalVariance += f0 * f0;
                                
                if (fProfit != 0.0)
                  {
                    // Only correlate non-zero pixels.
                    fCount1   += fRaw;
                    fCount1_2 += fRaw * fRaw;
                    fCount2   += fProfit;
                    fCount2_2 += fProfit * fProfit;
                    fCount12  += fRaw * fProfit;
                    nCount++;
                  }
              }
          }
      }

    if (nCount == 0) 
      fCorr = -1.0;
    else
      {
        f0 = sqrt(max(0.0,nCount*fCount1_2 - fCount1*fCount1));
        f1 = sqrt(max(0.0,nCount*fCount2_2 - fCount2*fCount2));
        f2 = nCount * fCount12  -  fCount1 * fCount2;
        fCorr = f2 / f0 / f1;        
      }

    if (bDebug)
      {
        Cstring sTemp = "profile";
        nDiagnostic(sTemp,oRaw,*this,fK);
      }

    fVariance = fTotalVariance * oRaw.m_fIntensity * oRaw.m_fIntensity;
    //fVariance = oRaw.m_fIntensity*oRaw.m_fIntensity/fSumVarDenom/fNorm/fNorm;
    fIntensity = fK * oRaw.m_fIntensity;
    return 0;
}


int tProfitRec::nDiagnostic(Cstring& sName,
                            tProfitRec& oRaw,
                            tProfitRec& oAverage,
                            double fScale)
{
  int nDim0,nDim1;
  int n0,n1;
  int nProfile;
  int nx;
  static int nProfileCount = 0;
  Cstring sOutputName;
        
        
  nDim0 = ms_a3nDim[0];
  nDim1 = ms_a3nDim[1];
  Cimage oImage(nDim0,nDim1*2,eImage_uI2);
  for (nProfile = 0; nProfile < 2; nProfile++)
    {
      for (n1 = 0, nx = 0; n1 < nDim1; n1++)
        {
          for (n0 = 0; n0 < nDim0; n0++,nx++)
            {
              if (nProfile == 0)
                oImage.nSetPixel(n0,n1,(unsigned short int) (1000.0*oAverage.m_pfX[nx]/oAverage.m_fContrib));
              else
                oImage.nSetPixel(n0,n1+nDim1,(unsigned short int) (1000.0*oRaw.m_pfX[nx]/oRaw.m_fContrib*fScale));
            }
        }
    }
  sOutputName = sName.before(".img");
  sOutputName += nProfileCount++;
  sOutputName += ".img";
  oImage.nWrite(sOutputName);
  
  return 0;
}


int tProfitRec_nAddFunctionPerSlice(int nDataPoints,
                                    double* pfInData,
                                    double* pfOutData)
{
  tProfitRec& oIn = *((tProfitRec*) pfInData);
  
  tProfitRec oOut(pfOutData, pfOutData+ oIn.ms_a3nDim[2],
                  *(pfOutData+2*oIn.ms_a3nDim[2]));
  oOut.nAdd(oIn);       
  pfOutData[2*oIn.ms_a3nDim[2]] = oOut.m_fContrib;

  return 0;
}


//////////////////////////////////////
//////////////////////////////////////
/////////  Cprofit2DPerSlice /////////
//////////////////////////////////////
//////////////////////////////////////


Cprofit2DPerSlice::Cprofit2DPerSlice()
{
  m_nProfitNumReflections      = 10;
  m_nProfitMaxImageRange       = 10;
  m_nCurrentImage              = -1;
  m_fMinIoverSigForContributor = 4.0;
  m_a2nProfileDim[0]           = 30;
  m_a2nProfileDim[1]           = 30;
  m_poTempProfitRec            = NULL;

  m_nMaxProfilesLevel          = 8;
  m_nMaxMeshesToKeep           = 3;
  m_pcBitmaps                  = NULL;
  m_nBitmapAlloc               = 0;
  m_nBitmapUsed                = 0;

  m_nMinImage                  = -1;
  m_nMaxImage                  = -1;
};

Cprofit2DPerSlice::~Cprofit2DPerSlice() {
  int nx;
  for (nx=0;nx < m_apoProfiles.size(); nx++)
    {
      if (NULL != m_apoProfiles[nx]) 
        {
          delete m_apoProfiles[nx];
          m_apoProfiles[nx] = NULL;
        }
    }
  if (NULL != m_poTempProfitRec)
    {
      delete m_poTempProfitRec;
      m_poTempProfitRec = NULL;
    }
  if (NULL != m_pcBitmaps)
    {
      delete[] m_pcBitmaps;
      m_pcBitmaps = NULL;
    }
}


int Cprofit2DPerSlice::nLoadProfile(int nSlice,
                                    C3Ddata& oData,
                                    bool bAddToLast,
                                    bool bLoadProfitRec,
                                    double* pfCOG)
{
  int n0,n1,nx;
  int nSliceOffset;
  int nSliceSize;
  float *afSliceBackground;
  float *pfXBuf;
  float *pfVarBuf;
  double f0,f1,f2;

  nSliceSize   = oData.m_nExt[0] * oData.m_nExt[1];
  nSliceOffset = nSlice * nSliceSize;
  afSliceBackground = oData.m_apfScratchFloat[m_nScratchFloatBackground];

  // Make sure the array is ready.

  m_afFloatBuf.setsize(nSliceSize * 2);
  pfXBuf   = &m_afFloatBuf[0];
  pfVarBuf = &m_afFloatBuf[nSliceSize];

  for (n1=0,nx=0; n1 < oData.m_nExt[1]; n1++)
    {
      for (n0=0; n0 < oData.m_nExt[0]; n0++,nx++)
        {
          // if this pixel is in the mask..               
          if (oData.ms_pnMask[nSliceOffset + nx] & oData.ms_nEllipsoidMask)
            {
              // We have a mask pixel.  
              // Check for saturated or nonunf pixels.

              if (oData.ms_pnMask[nx + nSliceOffset] & oData.ms_nSatMask)
                {
                  // If we are bootstrapping, then don't flag as bad 
                  // Otherwise, this reflection should get rejected, 
                  // no matter what slice we are on.
                  return 1;
                } 
              else if (oData.ms_pnMask[nx + nSliceOffset] & oData.ms_nBadMask)
                {
                  // Always return an error for this as well.
                  return 1;
                } 
              else
                {
                  if (oData.m_bFitToPlanar)
                    {
                      f1 = (afSliceBackground + nSlice*3)[0]*(float)n0/(float)oData.m_nExt[0] 
                        + (afSliceBackground + nSlice*3)[1]*(float)n1/(float)oData.m_nExt[1] 
                        + (afSliceBackground + nSlice*3)[2];
                                                
                    }
                  else
                    f1 = oData.m_a3fAvgBackground[2];
                  f0 = oData.m_pfData[nx + nSliceOffset] - f1;
                  f2 = oData.m_fAvgBackgroundSigma*oData.m_fAvgBackgroundSigma + max(0.0,f0)*oData.m_fDetectorGain;
                  if ((bAddToLast) && (pfXBuf[nx] != g_fMaskOut))
                    {
                      pfXBuf[nx]   += (float) f0;
                      pfVarBuf[nx] += (float) f2;
                    }
                  else
                    {
                      pfXBuf[nx]   = (float) f0;
                      pfVarBuf[nx] = (float) f2;
                    }
                }
            } 
          else
            {
              if (!bAddToLast) 
                pfXBuf[nx] = g_fMaskOut;
            }
        }
    }

  if (!m_poTempProfitRec)
    m_poTempProfitRec = new tProfitRec(m_a2nProfileDim);

  // Load the data into a profile
  if (bLoadProfitRec)
    {
      if (m_poTempProfitRec->nLoad(pfXBuf, pfVarBuf,
                                   oData.m_fAvgBackgroundSigma
                                    * oData.m_fAvgBackgroundSigma,
                                   oData.m_nExt,pfCOG))
        return 1;
    }
  return 0;
}


int Cprofit2DPerSlice::nAdd(int nSliceStart,int nSliceEnd,C3Ddata& oData)
{
  int nMeshBin;
  float* afSliceIntensity;
  float* afSliceSigma;
  int nSlice;
  int nSlicesLoaded;

  // We want to extract the pixels used in the background 
  // into an array with background removed.
  // This will be loaded into a tProfitRec.
        
    afSliceIntensity = oData.m_apfScratchFloat[m_nScratchFloatIntensity];
    afSliceSigma = oData.m_apfScratchFloat[m_nScratchFloatSigma];
    
    // Only profit for slices that are strong enough.
    for (nSlice = nSliceStart; nSlice <= nSliceEnd; nSlice++)
      {
        if ( afSliceIntensity[nSlice] 
            / max(1.0,afSliceSigma[nSlice]) >= m_fMinIoverSigForContributor) 
          break;
      }
    nSliceStart = nSlice;
    if (nSlice>nSliceEnd)
      return 1;

    for (nSlice = nSliceEnd; nSlice >= nSliceEnd; nSlice--)
      {
        if (afSliceIntensity[nSlice]
            / max(1.0,afSliceSigma[nSlice])>= m_fMinIoverSigForContributor) 
          break;
      }
    nSliceEnd = nSlice;

    nSlicesLoaded = 0;
    for (nSlice = nSliceStart; nSlice <= nSliceEnd; nSlice++)
      {
        if (nLoadProfile(nSlice,oData,nSlicesLoaded>0,(nSlice==nSliceEnd)))
          return 1;
        nSlicesLoaded++;
      }

    nSlice = (nSliceStart + nSliceEnd) / 2;
    // Now load the interpolated mesh with the new result.
    nMeshBin = nDetermineMesh(0,nSlice +  m_nCurrentImage); 
    if (nMeshBin == -1)
      return 1;
    m_apoProfiles[nMeshBin]->nAdd(max(0,
                                      min(m_nDim0-1,oData.m_nOrigOffset[0] + oData.m_nExt[0]/2)),
                                  max(0,
                                      min(m_nDim1-1,oData.m_nOrigOffset[1] + oData.m_nExt[1]/2)),
                                  (double*) m_poTempProfitRec,
                                  &tProfitRec_nAddFunctionPerSlice);
    nMeshBin = nDetermineMesh(1,nSlice + m_nCurrentImage);
    if (nMeshBin != -1)
      {
        m_apoProfiles[nMeshBin]->nAdd(max(0,min(m_nDim0-1,oData.m_nOrigOffset[0] + oData.m_nExt[0]/2)),max(0,min(m_nDim1-1,oData.m_nOrigOffset[1] + oData.m_nExt[1]/2)),(double*) m_poTempProfitRec,&tProfitRec_nAddFunctionPerSlice);
      }
    return 0;
}

int Cprofit2DPerSlice::nDetermineMesh(int nRelIndex,int nCurrentImageNumber)
{
  int nIndex;
  int nDeleteIndex;

  if (nCurrentImageNumber == -1)
    nCurrentImageNumber = m_nCurrentImage;
  if (m_nMaxImage == -1)
    {
      m_nMaxImage = m_nMinImage = m_nCurrentImage;
    }
  m_nMaxImage = max(m_nCurrentImage,m_nMaxImage);
  m_nMinImage = min(m_nCurrentImage,m_nMinImage);

  if ((nCurrentImageNumber <= -1) || (nRelIndex >= 2))
    return -1;

  while (  (!m_apoProfiles.size()) 
           || 
           (nCurrentImageNumber>m_anImageEnd.last() - m_nProfitMaxImageRange/2))
    {
      // Allocate a new mesh.
      if (!m_apoProfiles.size())
        {
          m_anImageStart + nCurrentImageNumber;
          m_anImageEnd + (nCurrentImageNumber + m_nProfitMaxImageRange - 1);
        } 
      else
        {
          m_anImageStart + (m_anImageEnd.last() + 1 - m_nProfitMaxImageRange/2);
          m_anImageEnd + (m_anImageStart.last() + m_nProfitMaxImageRange - 1);
        }
      m_apoProfiles + new CinterpMesh;
      m_apoProfiles.last()->vSetInterpConditions(m_nProfitNumReflections);
      m_apoProfiles.last()->vInit(m_nMaxProfilesLevel,
                                  (1+m_a2nProfileDim[0]*m_a2nProfileDim[1]*2),
                                  0,m_nDim0,0,m_nDim1);
    }
  nIndex = (nCurrentImageNumber - m_anImageStart[0])/(m_nProfitMaxImageRange/2);
  if (nIndex < 0)
    {
      // This can happen sometimes.
      nIndex = 0;
    }
    
  // Delete unwanted meshes.
  for (nDeleteIndex = nIndex - m_nMaxMeshesToKeep; 
       (nDeleteIndex >=0) && (m_apoProfiles[nDeleteIndex]); nDeleteIndex--)
    {
      delete m_apoProfiles[nDeleteIndex];
      m_apoProfiles[nDeleteIndex] = NULL;
    }

  // Make sure that the mesh needed really exists.
  while ((nIndex < m_apoProfiles.size()) && (!m_apoProfiles[nIndex]))
    nIndex++;
  if (nIndex >= m_apoProfiles.size())
    {
      printf("WARNING:  Missmatch of images during profile fitting.\n");
      return -1;
    }

  if (nRelIndex == 1)
    {
      if (   (nIndex - 1>=0) 
          && (nCurrentImageNumber >= m_anImageStart[nIndex-1])
          && (nCurrentImageNumber <= m_anImageEnd[nIndex-1]))
        nIndex--;
      else
        return -1;
      
    } 

  if ((nIndex < m_apoProfiles.size()) && (m_apoProfiles[nIndex]))
    return nIndex;
  return -1;
}

int Cprofit2DPerSlice::nFit(int nSlice,
                            C3Ddata& oData,
                            double a2fOverallCent[2],
                            float& fIntensityOut,
                            float& fSigmaOut)
{
  const int     nTestProfileCount = 4;
  int           anTestProfiles[nTestProfileCount];
  double        afTestProfilesIntensity[nTestProfileCount]; // Intensity from fitting to  each profile.
  double        afTestProfilesVariance[nTestProfileCount];  // Variance due to each profile.
  double        afTestProfilesCorr[nTestProfileCount];      // Correlation coeffs. due to each profile.

  int nMeshBin;
  int nProf;
  int nProfilesTested;
  int nx,ny;
  double *pfDataX,*pfDataX2,*pfDataCounts;
  double fIntensity,fVariance,fCorr;

  nMeshBin = nDetermineMesh(0,nSlice + m_nCurrentImage);
  if (nMeshBin == -1)
    return 1;
  if (-1 != (nx = nDetermineMesh(1,nSlice + m_nCurrentImage)))
    nMeshBin = nx;
  CinterpMesh& oMesh = *m_apoProfiles[nMeshBin];

  if (nLoadProfile(nSlice,oData,false,true,&a2fOverallCent[0]))
    return 1;

  if (-1 == oMesh.nGet(
                max(0,min(m_nDim0-1,oData.m_nOrigOffset[0] + oData.m_nExt[0]/2)),
                max(0,min(m_nDim1-1,oData.m_nOrigOffset[1] + oData.m_nExt[1]/2)),NULL,NULL))
    return 1;
        
  // We will be testing several profiles to see which gives the best fit.

  for (ny = 0; ny < nTestProfileCount; ny++)
    {
      if (oMesh.m_anInterpOffset.size() - ny - 1>=0)
        anTestProfiles[ny] = oMesh.m_anInterpOffset[oMesh.m_anInterpOffset.size() - ny - 1];
      else
        anTestProfiles[ny] = -1;
    }

  nProfilesTested = 0;
        
  for (nProf=0;nProf<nTestProfileCount;nProf++)
    {
      if (anTestProfiles[nProf]!=-1)
        {
          pfDataX = oMesh.pfGetItem(anTestProfiles[nProf]);
          pfDataX2 = pfDataX + m_a2nProfileDim[0]*m_a2nProfileDim[1];
          pfDataCounts = pfDataX2 + m_a2nProfileDim[0]*m_a2nProfileDim[1];
          tProfitRec oProf(pfDataX,pfDataX2,*pfDataCounts);
          if (oProf.nScale(*m_poTempProfitRec,fIntensity,fVariance,fCorr))
            {
              afTestProfilesVariance[nProf] = 1e20;
              afTestProfilesCorr[nProf] = -1.0;
            } 
          else 
            {
              afTestProfilesIntensity[nProf] = fIntensity; 
              afTestProfilesVariance[nProf] = fVariance;  
              afTestProfilesCorr[nProf] = fCorr;
              nProfilesTested++;
            }
        } 
      else 
        {
          afTestProfilesVariance[nProf] = 1e20;
          afTestProfilesCorr[nProf] = -1.0;
        }
    }
  // Did we have any profiles that worked?
  if (!nProfilesTested)
    return 1;

  // Search through the profiles, and choose the one that fits with the lowest variance.                
        
  for (nProf=0,nx=0;nx<nTestProfileCount;nx++)
    {
      if (afTestProfilesCorr[nx]>afTestProfilesCorr[nProf])
        {
          nProf = nx;
        }
    }

  // Return the fit for this profile.
  fIntensityOut = afTestProfilesIntensity[nProf];
  fSigmaOut     = sqrt(afTestProfilesVariance[nProf]);
  return 0;
}


int Cprofit2DPerSlice::nGetBitmapPixel(int nOffset,
                                       int a3nPixel[3],
                                       int a3nExt[3])
{
  int nx;

  nx  = a3nPixel[2]*a3nExt[0]*a3nExt[1] + a3nPixel[1]*a3nExt[0] + a3nPixel[0];
  return ((m_pcBitmaps[nOffset + nx/8] >> (nx % 8)) & 1);
}

int Cprofit2DPerSlice::nAllocBitmaps(int nDesiredSize)
{
  unsigned char* pcNewBitmaps;
  int nOldSize;
  int nNewSize;
  if (m_nBitmapUsed+nDesiredSize>m_nBitmapAlloc) 
    {
      nOldSize = m_nBitmapUsed;
      nNewSize = max(100000,nOldSize*2);
      pcNewBitmaps = new unsigned char[nNewSize];
      memset(pcNewBitmaps,0,nNewSize);
      if (nOldSize) 
        memcpy(pcNewBitmaps, m_pcBitmaps, nOldSize);
      m_nBitmapAlloc = nNewSize;
      if (NULL != m_pcBitmaps)
        delete [] m_pcBitmaps;
      m_pcBitmaps = pcNewBitmaps;
    }
  return 0;
}

int Cprofit2DPerSlice::nWriteScanBitmap(Cstring& sName,Cimage_header* poHeader)
{
  int nImageNumber;
  int nRef;
  int nOffset;
  int nValue;
  int a3nExt[3];
  int a3nOrigOffset[3];
  int a3nOffset[3];
      
  int nx;
  CScanBitmap oScanBitmap;

  printf("Writing scan bitmap.\n");

  nFileAppendVersion(sName, TRUE);

  if ((m_nMaxImage == -1) || (!m_pcBitmaps))
    return 0;

  if (oScanBitmap.nOpenOutput(sName,m_nDim0,m_nDim1,m_nMinImage,m_nMaxImage,3,0))
    return 1;

  for (nImageNumber = m_nMinImage; nImageNumber <= m_nMaxImage; nImageNumber++)
    {
      oScanBitmap.vClearLocalBitmap();
      for (nRef = 0; nRef < m_anBitmapOffset.size(); nRef++)
        {
          for (nx=0; nx < 3; nx++)
            {
              a3nExt[nx] = m_a3nExt[nx][nRef];
              a3nOrigOffset[nx] = m_a3nOrigOffset[nx][nRef];
            }
          nOffset = m_anBitmapOffset[nRef];
           
          if (   (a3nOrigOffset[2]<=nImageNumber) 
              && (a3nOrigOffset[2]+a3nExt[2]>nImageNumber))
            {
              a3nOffset[2] = (nImageNumber - a3nOrigOffset[2]);                

              int anUsedValues[8];

              for (nx=0;nx<8;nx++)
                anUsedValues[nx] = 0;
                
              for (a3nOffset[1] = 0;(a3nOffset[1]<a3nExt[1]);a3nOffset[1]++)
                {
                  for (a3nOffset[0] = 0;(a3nOffset[0]<a3nExt[0]);a3nOffset[0]++)
                    {
                      anUsedValues[oScanBitmap.nGetLocalBitmap(a3nOffset[0] + a3nOrigOffset[0],a3nOffset[1] + a3nOrigOffset[1])]=1;
                    }
                }
                // Get an unused value.
              for (nValue=1;anUsedValues[nValue] && (nValue+1<=6);nValue++);
              for (a3nOffset[1] = 0;(a3nOffset[1]<a3nExt[1]);a3nOffset[1]++)
                {
                  for (a3nOffset[0] = 0;(a3nOffset[0]<a3nExt[0]);a3nOffset[0]++)
                    {
                      if (nGetBitmapPixel(nOffset,a3nOffset,a3nExt))
                        {
                          if (oScanBitmap.nGetLocalBitmap(a3nOffset[0] + a3nOrigOffset[0],a3nOffset[1] + a3nOrigOffset[1]))
                            oScanBitmap.vSetLocalBitmap(a3nOffset[0] + a3nOrigOffset[0],a3nOffset[1] + a3nOrigOffset[1],7);
                          else
                            oScanBitmap.vSetLocalBitmap(a3nOffset[0] + a3nOrigOffset[0],a3nOffset[1] + a3nOrigOffset[1],nValue);
                        }
                    }                        
                }
            }
        }
      // oScanBitmap.nWriteMaskedImage("test.img",NULL);

      if (oScanBitmap.nWriteNextBitmap())
        return 1;
    }
  printf("Finished writing scan bitmap.\n");
  if (poHeader) 
    oScanBitmap.nUpdateHeader(*poHeader);

  if (oScanBitmap.nCloseOutput())
    return 1;
  return 0;
}

int  Cprofit2DPerSlice::nAddBitmap(C3Ddata& oData)
{
  int a3nPixel[3];
  int a3nExt[3];
  int a3nOrigOffset[3];
  int nStart,nEnd;
  int nSlice,nSlicePoint;
  int n0,n1;
  int nBitPoint;
  int nx;
    

  oData.vGetPeakStartEnd(&nStart,&nEnd);
  m_a3nExt[0] + (a3nExt[0] = (oData.nGetEllipsoidRange(0,1) - oData.nGetEllipsoidRange(0,0) + 1));
  m_a3nExt[1] + (a3nExt[1] = (oData.nGetEllipsoidRange(1,1) - oData.nGetEllipsoidRange(1,0) + 1));
  m_a3nExt[2] + (a3nExt[2] = (nEnd - nStart + 1));
  m_a3nOrigOffset[0] + (a3nOrigOffset[0] = oData.nGetEllipsoidRange(0,0) + oData.m_nOrigOffset[0]);
  m_a3nOrigOffset[1] + (a3nOrigOffset[1] = oData.nGetEllipsoidRange(1,0) + oData.m_nOrigOffset[1]);
  m_a3nOrigOffset[2] + (a3nOrigOffset[2] = m_nCurrentImage + nStart);
  m_anBitmapOffset + m_nBitmapUsed;

  nBitPoint = a3nExt[0]*a3nExt[1]*a3nExt[2];
  nAllocBitmaps(nBitPoint/8 + ((nBitPoint % 8)!=0));

  nBitPoint = 0;

  for (a3nPixel[2] = nStart; a3nPixel[2] <= nEnd; a3nPixel[2] ++)
    {
      nSlicePoint = 0;
      nSlice = a3nPixel[2];
      for (a3nPixel[1] = 0; a3nPixel[1] < a3nExt[1]; a3nPixel[1] ++)
        {
          n0 = a3nOrigOffset[0] - oData.m_nOrigOffset[0];
          n1 = a3nPixel[1] + a3nOrigOffset[1] - oData.m_nOrigOffset[1];
          nx = n0 + n1*oData.m_nExt[0] + 
            nSlice*oData.m_nExt[0]*oData.m_nExt[1];
          for (a3nPixel[0] = 0; a3nPixel[0] < a3nExt[0]; a3nPixel[0]++,nx++,n0++,nSlicePoint++)
            {
              if (oData.pnGetMask()[nx] & C3Ddata::ms_nEllipsoidMask)
                {
                  m_pcBitmaps[nBitPoint/8+m_nBitmapUsed] |= 1 << (nBitPoint % 8);
                }
              nBitPoint++;
            }
        }
    }
    m_nBitmapUsed += (nBitPoint/8) + ((nBitPoint % 8)!=0);
    return 0;
}

int Cprofit2DPerSlice::nDiagnostic(Cstring& sName,int nImage)
{
    /*
    int nx,ny;
    int n0,n1;
    int nn0,nn1;
    int nStat;
    int nInterp;
        int nMeshIndex;
        int a2nMeshIndex[2];
        int nPix0,nPix1;
    float fStep0,fStep1;
        double* pfDataX,*pfDataX2,*pfDataCounts,*pfData;
        const int nBestPrint = 2;

    int     anBestPrintInterp[nBestPrint];
        int     anBestPrintMesh[nBestPrint];
    double  afBestPrintValues[nBestPrint];


    
    Cimage oImage(m_nDim0,m_nDim1,eImage_uI2);


        a2nMeshIndex[0] = nDetermineMesh(0,nImage);
        if (a2nMeshIndex[0] < 0) 
                return 1;
        a2nMeshIndex[1] = nDetermineMesh(1,nImage);
        if (a2nMeshIndex[1] < 0)
                a2nMeshIndex[1] = a2nMeshIndex[0];

        CinterpMesh* a2poMesh[2] = { m_apoProfiles[a2nMeshIndex[0]],m_apoProfiles[a2nMeshIndex[1]] };

    // Zero the image.
    for (n0=0;n0<m_nDim0;n0++) {
        for (n1=0;n1<m_nDim1;n1++) {
            oImage.nSetPixel(n0,n1,(unsigned short int) 0);
        };
    };
    fStep0 = m_nDim0/10.0;
    fStep1 = m_nDim1/10.0;
    for (n0=(int) (fStep0/2);n0+fStep0/2<m_nDim0;n0+=(int) fStep0) {
        for (n1=(int) (fStep1/2);n1+fStep1/2<m_nDim1;n1+=(int) fStep1) {
            

                        for (nMeshIndex = 0; nMeshIndex < 2; nMeshIndex++) {
                                nStat = a2poMesh[nMeshIndex]->nGet(n0,n1,NULL);
                                if (nStat != -1) {
                                        // Choose the best solutions at this location (have the greatest # of contributors.
                    for (nx=0;nx<nBestPrint;nx++)
                        afBestPrintValues[nx] = 0.0;

                    pfDataX = a2poMesh[nMeshIndex]->pfGetItem(a2poMesh[nMeshIndex]->m_anInterpOffset.last());
                    pfDataX2 = pfDataX + m_a2nProfileDim[0]*m_a2nProfileDim[1];
                    pfDataCounts = pfDataX2 + m_a2nProfileDim[0]*m_a2nProfileDim[1];
                    
                    int   nLowest;
                    for (nx=0,nLowest=0;nx<nBestPrint;nx++) {
                        if (afBestPrintValues[nx]<afBestPrintValues[nLowest])
                            nLowest = nx;
                    };
                    if ((*pfDataCounts>=1.0) && (*pfDataCounts > afBestPrintValues[nLowest])) {
                        afBestPrintValues[nLowest] = *pfDataCounts;
                        anBestPrintInterp[nLowest] = a2poMesh[nMeshIndex]->m_anInterpOffset.last();
                        anBestPrintMesh[nLowest] = nMeshIndex;
                    };
                                };
                        };
                        
                        for (nx=0;nx<nBestPrint;nx++) {
                                if (afBestPrintValues[nx]>0.0) {
                                        CinterpMesh& oMesh = *a2poMesh[anBestPrintMesh[nx]];
                                        nInterp = anBestPrintInterp[nx];
                                        oMesh.nGet(n0,n1,NULL);
                                        pfDataX = oMesh.pfGetItem(nInterp);
                                        pfDataX2 = pfDataX + m_a2nProfileDim[0]*m_a2nProfileDim[1];
                                        pfDataCounts = pfDataX2 + m_a2nProfileDim[0]*m_a2nProfileDim[1];
                                        pfData = pfDataX;
                                        tProfitRec oTemp(pfDataX,pfDataX2,*pfDataCounts);
                                        for (nPix1=0;nPix1<m_a2nProfileDim[0];nPix1++) {
                                                for (nPix0=0;nPix0<m_a2nProfileDim[0];nPix0++,pfData++) {
                                                        nn1 = n1+nPix1;
                                                        nn0 = n0+nPix0+nx*(1+m_a2nProfileDim[0]);
                                                        if ((nn0>=0) && (nn1>=0) && (nn0<m_nDim0) && (nn1<m_nDim1))
                                                                (oImage.*oImage.prnSetPixel)(nn0,nn1,10000*(*pfData)/(*pfDataCounts));
                                                };
                                        };
                                };
                        };

        };
    };
        Cstring sOutputName;
        sOutputName = sName.before(".img");
        sOutputName += nImage;
        sOutputName += ".img";

    oImage.nWrite(sOutputName);
    */
     return 0;
};

//////////////////////////////////////
//////////////////////////////////////
//////// Cprofit2DSlices /////////////
//////////////////////////////////////
//////////////////////////////////////


// TJN says Try playing with fSigmaProfit vs fSigma


#define DONTUSEPROFIT  ((fIntensityProfit == 0.0) \
                      || (fSigmaProfit > 4.0*fSigma))
//
//                      || (fIntensity>m_fMaxIntensityToProfit)) )

//Remove these 2 things from previous line per TJN 8-Apr-2005
//&& ((!(nPartialStat & C3Ddata::ms_nPerSliceIntegrateWeak) || (fIntensityProfit == 0.0))))

Cprofit2DSlices::Cprofit2DSlices()
{
  m_poIntegrate = NULL;
  m_poPartial = NULL;
  m_poScaleList = NULL;
  m_poAxialRefs = NULL;
  m_poSource = NULL;
  m_poCrystal = NULL;
  m_fBroadening = 1.0;

  m_nMinRefsPerProfile = 500;
  m_nResoGroups = g_nMaxResoGroups;
  m_nNumProfiles = 0;
  m_nProfilePoints = 0;

  m_fGainEst = -1.0;
  m_fMaxIntensityToProfit = 1000000;

  m_fStandardProfileStep = 0.01;
  m_fStandardProfileMax = 2.0;
  m_fIntersectMax = 0.03;

  m_nTotalIntersections = 0;
}

int Cprofit2DSlices::nLoadHeader(Cstring& sHeader)
{
  Cimage_header oHeader(sHeader);
  if (!oHeader.bIsAvailable())
    return 1;
  return nLoadHeader(oHeader);
}

int Cprofit2DSlices::nLoadHeader(Cimage_header& oHeader)
{
  if (m_poSource)
    delete m_poSource;
  m_poSource = new Csource(oHeader);
  if (!m_poSource->bIsAvailable())
    return 1;
  m_poCrystal = new Ccrystal(oHeader);
  if (!m_poCrystal->bIsAvailable())
    return 1;
  return 0;
}

Cstring Cprofit2DSlices::sCreatePartName(Cstring& sName,int nPart)
{
  Cstring sTemp;
  Cstring sTemp2;

  if (nPart == -1)
    return sName;
  else
    {
      sTemp = "Part";
      if (nPart == -2)
        sTemp += "*";
      else
        sTemp += nPart;
      sTemp += "_";
      sTemp += sFileGetBasename(sName);
      sTemp2 = sFileGetDirectory(sName);
      return sFileBuildName(sTemp2, sTemp);
    }
}

int Cprofit2DSlices::nWriteParts(bool bWriteIntegrate,bool bWritePartial,bool bWriteEllipsoid)
{
  int nPart;
  Cstring sTemp;

  if (m_anIntegrateFilePartStart.size()<=1)
    return 0;
  if (m_poIntegrate)
    delete m_poIntegrate;
  m_poIntegrate = NULL;
  if (m_poPartial)
    delete m_poPartial;
  m_poPartial = NULL;
  if (m_poScaleList)
    delete m_poScaleList;
  m_poScaleList = NULL;
  
  if (bWriteIntegrate)
    {
      Creflnlist oIntegrate(sCreatePartName(m_sIntegrate,0));
      for (nPart = 1; nPart < m_anIntegrateFilePartStart.size(); nPart++)
        {
          Creflnlist oIntegrateNext(sCreatePartName(m_sIntegrate,nPart));
          oIntegrate.nInsertListFrom(oIntegrateNext);
        }
      oIntegrate.nWrite(m_sIntegrate);
      (void) nFileDelete(sCreatePartName(m_sIntegrate,-2));
    };

  if (bWritePartial)
    {
      Creflnlist oPartial(sCreatePartName(m_sPartial,0));
      for (nPart = 1; nPart < m_anIntegrateFilePartStart.size(); nPart++)
        {
          Creflnlist oPartialNext(sCreatePartName(m_sPartial,nPart));
          oPartial.nInsertListFrom(oPartialNext);
        }
      oPartial.nWrite(m_sPartial);
      (void) nFileDelete(sCreatePartName(m_sPartial,-2));
    }


  if (bWriteEllipsoid)
    {
      Creflnlist oEllipsoids(sCreatePartName(m_sEllipsoid,0));
      for (nPart = 1; nPart < m_anIntegrateFilePartStart.size(); nPart++)
        {
          Creflnlist oEllispoidsNext(sCreatePartName(m_sEllipsoid,nPart));
          oEllipsoids.nInsertListFrom(oEllipsoids);
        }
      oEllipsoids.nWrite(m_sEllipsoid);
      (void) nFileDelete(sCreatePartName(m_sEllipsoid,-2));
    }
    return 0;
}
        
int Cprofit2DSlices::nLoad(Cstring& sIntegrateList,
                           Cstring& sPartialList,
                           Cstring& sScaleList,
                           int nFilePart)
{
    int nx;
    int nIntRef;
    double fRotMid;
   
    if (m_poIntegrate)
      delete m_poIntegrate;
    m_poIntegrate = new Creflnlist(sIntegrateList,(nFilePart==-1)?-1:(m_anIntegrateFilePartEnd[nFilePart] - m_anIntegrateFilePartStart[nFilePart] + 1),(nFilePart == -1)?-1:(m_anIntegrateFilePartStart[nFilePart]));
    m_sIntegrate = sIntegrateList;
    if (!m_poIntegrate->bIsAvailable())
      return 1;
    if (m_poPartial)
      delete m_poPartial;
    if (m_poScaleList)
      delete m_poScaleList;
    m_poScaleList = NULL;
    m_sPartial = sPartialList;
    m_poPartial = new Creflnlist(sPartialList,(nFilePart==-1)?-1:(m_anPartialFilePartEnd[nFilePart] - m_anPartialFilePartStart[nFilePart] + 1),(nFilePart == -1)?-1:(m_anPartialFilePartStart[nFilePart]));
    if (!m_poPartial->bIsAvailable())
      return 1;
    if (sScaleList.length() && (nFilePart == -1))
      {
        m_poScaleList = new Creflnlist(sScaleList);
        if (!m_poScaleList->bIsAvailable())
          return 1;
      }
    
    m_sEllipsoid = sTransSymbol(sDtrekGetPrefix());
    m_sEllipsoid += "ellipsoids.ref";

    m_nIntRefs = m_poIntegrate->nGetNumReflns();

    if  (   (m_poIntegrate->m_nFI_fObsRotMid < 0)
         || (m_poIntegrate->m_nFI_fObsRotWidth < 0) 
         || (m_poIntegrate->m_nFI_fCalcRotMid < 0)
         || (m_poPartial->m_nFI_fCalcRotMid < 0)
         || (m_poPartial->m_nFI_fObsRotMid < 0) )
      {
        printf("ERROR: The proper fields were not in the reflection list for truncaction.\n");
        return 1;
      }

    // Find refinement batches.
    int nCurrentRefineBatch = -999; //some crazy value to initialize

    -m_anRefineBatchStart;       // Start/End of refinement batches used for profile fitting.
    -m_anRefineBatchEnd;

    if (m_nIntRefs)
      {
        nCurrentRefineBatch = (*m_poIntegrate)[0].nGetField(m_poIntegrate->m_nFI_nRefineBatch);
        -m_anRefineBatchStart + 0;
        -m_anRefineBatchEnd + 0;
      }
    for (nIntRef = 0; nIntRef < m_nIntRefs; nIntRef++)
      {
        if (nCurrentRefineBatch == (*m_poIntegrate)[nIntRef].nGetField(m_poIntegrate->m_nFI_nRefineBatch))
          m_anRefineBatchEnd.last() = nIntRef;
        else
          {
            nCurrentRefineBatch = (*m_poIntegrate)[nIntRef].nGetField(m_poIntegrate->m_nFI_nRefineBatch);
            m_anRefineBatchEnd + nIntRef;
            m_anRefineBatchStart + nIntRef;
          }
      }

    double fCosTheta;
    double fReso;
    double fDstarSq;
    double fDstar;
    double fWavelength;
    double fTheta;

    // Assign various variables.

    -m_afLorentz;
    -m_afTheta;
    -m_afDStar;

    fWavelength = m_poSource->m_poWavelength->fGetWavelength();

    for (nIntRef = 0; nIntRef < m_nIntRefs; nIntRef++)
      {
        m_afLorentz + (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fLorentz);
        fReso = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fResolution);
        fDstar     = fWavelength/fReso;
        fDstarSq   = fDstar*fDstar;
        fCosTheta  = sqrt(fDstarSq - 0.25 * fDstarSq * fDstarSq);
        fCosTheta /= fDstar;
        fTheta     = acos(max(-1.0,min(1.0,fCosTheta)));
        m_afTheta + fTheta;
        m_afDStar + fDstar;
        fRotMid = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fCalcRotMid);
      }
    m_pfLorentz = & m_afLorentz[0];
    m_pfTheta   = & m_afTheta[0];
    m_pfDStar   = & m_afDStar[0];

    // Determine the resolution rings.
    itr<int> anThetaSort;
    int nRef,nRefSort;

    // If we break the data up into the given resolution bins, will we still have enough reflections per group?
    if (m_nIntRefs/m_nResoGroups < m_nMinRefsPerProfile)
      {
        // No?  Decrease the # of groups.
        m_nResoGroups = m_nIntRefs/m_nMinRefsPerProfile;
      }
    if (m_nResoGroups<1)
      m_nResoGroups = 1;

    // Sort based on the Theta value.
    g_pfCmpDoubles = m_pfTheta;
    anThetaSort.setsize(m_afTheta.size());
    for (nx = 0; nx < m_afTheta.size(); nx++)
        anThetaSort[nx] = nx;
    
    if( m_afTheta.size() > 0 ) // safety check
        qsort(&anThetaSort[0], m_afTheta.size(), sizeof(int), double_cmp_rel);
    
    // Determine cutoffs for each resolution bin.
    -m_afShellThetaGroupMin;
    -m_afShellThetaGroupMax;
    -m_afShellResoGroupMin;
    -m_afShellResoGroupMax;

    nx = 0;

    for (nRefSort=0;nRefSort<m_nIntRefs;nRefSort++) {
        nRef=anThetaSort[nRefSort];

        if ((nx == 0) && (m_afShellThetaGroupMin.size() == m_afShellThetaGroupMax.size())) {
            m_afShellThetaGroupMin + m_pfTheta[nRef];
            m_afShellResoGroupMin + (fWavelength/m_pfDStar[nRef]);
        };
        if (1) { // Put test here of accepting a reflection.
            if (nx == 0) {
                m_afShellThetaGroupMax + m_pfTheta[nRef];
                m_afShellResoGroupMax + (fWavelength/m_pfDStar[nRef]);
            }

            m_afShellThetaGroupMax.last() = m_pfTheta[nRef];
            m_afShellResoGroupMax.last() = (fWavelength/m_pfDStar[nRef]);
            nx++;
            if ((nx >= m_nIntRefs/m_nResoGroups) && (m_afShellThetaGroupMin.size()!=m_nResoGroups))
                nx = 0;
        };
    };
    // Now go back and assign resolution bins.
    m_anResoBatch.setsize(m_nIntRefs);
    m_pnResoBatch = & m_anResoBatch[0];
    for (nRef = 0; nRef < m_nIntRefs; nRef++) {
        for (nx = 0; nx < m_afShellThetaGroupMin.size();nx++) {
            if ((m_pfTheta[nRef]>= m_afShellThetaGroupMin[nx]) && (m_pfTheta[nRef] <= m_afShellThetaGroupMax[nx]))
                break;
        };
        m_pnResoBatch[nRef] = min(m_afShellThetaGroupMax.size()-1,nx);
    };
    if (nLoadProfileBatches())
        return 1;

        return 0;
};

int Cprofit2DSlices::nLoadProfileBatches()
{
    int    nBatch;
    int    nRef;
    int    nRefs = 0;
    int    nResoGroup;
    
    // Determine appropriate widths for profiling batches.
    int    nNextRefineBatchStart = 0;
    int    nProfileBatchStart,nProfileBatchEnd;
    int    nRefineBatchStart,nRefineBatchEnd;

    m_a2anProfileBatch[0].setsize(m_nIntRefs);
    m_a2anProfileBatch[1].setsize(m_nIntRefs);
    m_anBestProfileBatchi.setsize(m_nIntRefs);
    m_a2pnProfileBatch[0] = & m_a2anProfileBatch[0][0];
    m_a2pnProfileBatch[1] = & m_a2anProfileBatch[1][0];
    m_pnBestProfileBatchi = & m_anBestProfileBatchi[0];
    for (nRef = 0; nRef < m_nIntRefs; nRef++)
      {
        m_a2pnProfileBatch[0][nRef] = -1;
        m_a2pnProfileBatch[1][nRef] = -1;
        m_pnBestProfileBatchi[nRef] = -1;
      }
    
    m_nNumProfiles = 0;
    
    for (nResoGroup = 0; nResoGroup < m_nResoGroups; nResoGroup++)
      {
        nNextRefineBatchStart = 0;
        
        -m_anProfileBatchStart[nResoGroup];
        -m_anProfileBatchEnd[nResoGroup];
        -m_anProfileBatchContrib[nResoGroup];
        -m_anProfileBatchIndex[nResoGroup];
        for (nBatch = 0;1;nBatch++)
          {
            nRefineBatchStart = nNextRefineBatchStart;
            nRefineBatchEnd = nNextRefineBatchStart;        
            nProfileBatchStart = m_anRefineBatchStart[nRefineBatchStart];
            nProfileBatchEnd = m_anRefineBatchEnd[nRefineBatchEnd];
            
            while (1)
              {
                // Count how many reflections are in the profile batch.
                nRefs = 0;
                for (nRef = nProfileBatchStart; nRef <= nProfileBatchEnd;nRef++)
                  {
                    if (m_pnResoBatch[nRef] == nResoGroup)
                      nRefs++;
                  }
                if (nRefs >= m_nMinRefsPerProfile) 
                  break;
                if (nRefineBatchEnd + 1 < m_anRefineBatchStart.size())
                  {
                    nRefineBatchEnd++;
                    nProfileBatchEnd = m_anRefineBatchEnd[nRefineBatchEnd];
                  }
                else
                  break;
              }
            if (nRefs < m_nMinRefsPerProfile)
              {
                if (m_anProfileBatchEnd[nResoGroup].size())
                  {
                    m_anProfileBatchEnd[nResoGroup].last() = nProfileBatchEnd;
                  } 
                else
                  {
                    m_anProfileBatchStart[nResoGroup] + nProfileBatchStart;
                    m_anProfileBatchEnd[nResoGroup] + nProfileBatchEnd;
                    m_anProfileBatchContrib[nResoGroup] + 0;
                    m_anProfileBatchIndex[nResoGroup] + (m_nNumProfiles++);
                  }
                break;
              } 
            else
              {
                m_anProfileBatchStart[nResoGroup] + nProfileBatchStart;
                m_anProfileBatchEnd[nResoGroup] + nProfileBatchEnd;
                m_anProfileBatchContrib[nResoGroup] + 0;
                m_anProfileBatchIndex[nResoGroup] + (m_nNumProfiles++);
              }      
            nNextRefineBatchStart = max(nRefineBatchStart + 1,(nRefineBatchStart + nRefineBatchEnd)/2);
            if (nNextRefineBatchStart>=m_anRefineBatchStart.size())
              break;
          }
        
        for (nRef = 0; nRef < m_nIntRefs; nRef++)
          {
            if (m_pnResoBatch[nRef] == nResoGroup)
              {
                // Find the profile batch in which the reflection participates.
                for (nBatch = 0; nBatch < m_anProfileBatchStart[nResoGroup].size(); nBatch++)
                  {
                    if ((nRef >= m_anProfileBatchStart[nResoGroup][nBatch]) && (nRef <= m_anProfileBatchEnd[nResoGroup][nBatch]))
                      {
                        if (m_a2pnProfileBatch[0][nRef] == -1)
                          {
                            m_a2pnProfileBatch[0][nRef] = nBatch;
                            m_anProfileBatchContrib[nResoGroup][nBatch]++;
                          } 
                        else if (m_a2pnProfileBatch[1][nRef] == -1)
                          {
                            m_a2pnProfileBatch[1][nRef] = nBatch;
                            m_anProfileBatchContrib[nResoGroup][nBatch]++;
                          }
                      }
                  }
              }
          }
      }
    return 0;
}


Cprofit2DSlices::~Cprofit2DSlices()
{
  if (m_poIntegrate)
    delete m_poIntegrate;
  if (m_poPartial)
    delete m_poPartial;
  if (m_poAxialRefs)
    delete m_poAxialRefs;
  if (m_poSource)
    delete m_poSource;
  if (m_poCrystal)
    delete m_poCrystal;
}


int Cprofit2DSlices::nCalcRemoveNoise(bool bFirstPart,bool bLastPart)
{
  Cstat oStat;
  Cstat oStat2;
  double f0;
  int nPass;
  int nStat;
  int nRef=0;
  int nIntRef,nPartialRef;
  int nBin;

  int nPartialRefs;

  
  int nRefsRemovedBadPosition                 = 0;
  int nRefsRemovedBadProfitFit                = 0;
  int nRefsRemovedOverlap                     = 0;    // Previous removal.
  int nRefsForOutput                          = 0;
  int nRefsAdjustedForNoise                   = 0;
  int nRefsWithMixedProfitIntegratedIntensity = 0;
  int nRefsWithOnlyProfitIntensity            = 0;
  int nRefsWithOnlyIntegratedIntensity        = 0;
  int nRefsWithNoPartialInformation           = 0;
  int nRefsWithIntensityBeyondMosaicityBounds = 0;
  int nRefsWithIntensityOnlyAtCentroid        = 0;
  int nProfitRefsUsed                         = 0;
  int nProfitRefsExamined                     = 0;

  static int s_nRefsRemovedBadPosition;
  static int s_nRefsRemovedBadProfitFit;
  static int s_nRefsRemovedOverlap;
  static int s_nRefsForOutput;
  static int s_nRefsAdjustedForNoise;
  static int s_nRefsWithMixedProfitIntegratedIntensity;
  static int s_nRefsWithOnlyProfitIntensity;
  static int s_nRefsWithOnlyIntegratedIntensity;
  static int s_nRefsWithNoPartialInformation;
  static int s_nRefsWithIntensityBeyondMosaicityBounds;
  static int s_nRefsWithIntensityOnlyAtCentroid;
  static int s_nProfitRefsUsed;
  static int s_nProfitRefsExamined;
  

  double fObs0,fObs1,fObsRot,fCalc0,fCalc1,fCalcRot;
  double fPxDiff,fRotDiff;
  double fPxAverage  = 0.0, fPxDeviation = 0.0;
  double fRotAverage,fRotDeviation;
  double fPxDeviationSigma = 5.0;


  if ((!m_poPartial) || (!m_poPartial->bIsAvailable()))
    return 1;
  nPartialRefs = m_poPartial->nGetNumReflns();
  m_nIntRefs = m_poIntegrate->nGetNumReflns();

  nRefsRemovedOverlap = m_nTotalIntersections;

  // Find the average deviation PxObs from PxCalc and RotObs from RotCalc.
  // We will want to reject these reflections.
        
  for (nPass = 0; nPass < 2; nPass++)
    {
      oStat.vClearRaw();
      oStat2.vClearRaw();
      for (nIntRef = 0; nIntRef < m_nIntRefs; nIntRef++)
        {
          fObs0 = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fObsPx0);
          fObs1 = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fObsPx1);
          fObsRot = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fObsRotMid);
          fCalc0 = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fCalcPx0);
          fCalc1 = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fCalcPx1);
          fCalcRot = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fCalcRotMid);
          fPxDiff = sqrt((fObs0 - fCalc0)*(fObs0 - fCalc0) + (fObs1 - fCalc1)*(fObs1 - fCalc1));
          fRotDiff = fObsRot - fCalcRot;
          if (nPass == 0)
            {
              oStat.vAddRaw(fPxDiff);
              oStat2.vAddRaw(fRotDiff);
            }
          else
            {
              if ((fPxDiff - fPxAverage > fPxDeviationSigma*fPxDeviation) && ((*m_poIntegrate)[nIntRef].m_nSelect==0))
                {
                  (*m_poIntegrate)[nIntRef].m_nSelect = 1;
                  nRefsRemovedBadPosition++;
                }
            }
        }
      if (nPass == 0)
        {
          fPxAverage = oStat.fAverageRaw();
          fPxDeviation = oStat.fStandardDeviationRaw();
          fRotAverage = oStat2.fAverageRaw();
          fRotDeviation = oStat2.fStandardDeviationRaw();
        }
    }


  if (m_fGainEst < 0.0)
    {
      // Only estimate the gain if we didn't have an explicit value.

      if (m_nIntRefs > 100)
        {
          // Estimate the gain.
          itr<double> afBackground;
          itr<double> afVariance;
          double*     pfBackground;
          double*     pfVariance;
          double*     a2fFunctions[2];
          double      afArgs[2];
          double      afSigmas[2];
          double      fGain;
          double      fPedestal;
          double      fDetectorNoise;
          afBackground.setsize(m_nIntRefs);
          pfBackground = & afBackground[0];
          afVariance.setsize(m_nIntRefs);
          pfVariance = & afVariance[0];
          a2fFunctions[0] = pfBackground;
          a2fFunctions[1] = ms_pfConst;
          for (nIntRef = 0;  nIntRef < m_nIntRefs; nIntRef++)
            {
              pfBackground[nIntRef] = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fBackground);
              pfVariance[nIntRef] = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fBackgroundSigma);
              pfVariance[nIntRef] *= pfVariance[nIntRef];
            }
          if ((!nSolveLeastSquares(2,m_nIntRefs,pfVariance,NULL,&a2fFunctions[0],afArgs,afSigmas,NULL,NULL)) && (afArgs[0]>=0.0))
            {
              fGain = afArgs[0];
              fPedestal = pfBackground[(int) (0.02*m_nIntRefs)];
              fDetectorNoise = sqrt(max(0.0,fGain*fPedestal + afArgs[1]));
              
              // Sort the background variables, and find the one that is 5% above the lowest.  This is probably
              // We assume this is very close to a zero photon reading.
              qsort(pfBackground,m_nIntRefs,sizeof(double),double_cmp);
                
              m_fGainEst = fGain;
              if ("" != sGetEnv("DTREK_VERBOSE"))
                {
                  printf("Estimation of Detector Gain\n");
                  printf("----------------------------------------------\n");
                  printf("Gain                 %8.2lf\n",m_fGainEst);
                  printf("SQRT(Gain)           %8.2lf\n",sqrt((double)m_fGainEst));
                  printf("Min Background       %8.2lf\n",fPedestal);
                  printf("Detector Noise       %8.2lf\n",fDetectorNoise);
                  printf("----------------------------------------------\n");
                }
                
              nStat = 0;
            } 
          else
            nStat = 1;
        } 
      else 
        nStat = 1;
    }
  else 
    nStat = 0;
    
  if (nStat)
    {
      printf("Estimation of Detector Gain\n");
      printf("----------------------------------------------\n");
      printf("WARNING:             Could not estimate detector gain.\n");
      printf("Gain                 %8.2lf\n",m_fGainEst);
      printf("SQRT(Gain)           %8.2lf\n",sqrt((double)m_fGainEst));
      printf("----------------------------------------------\n");
      nStat = 0;
    }


  // The code in this block was an effort to reject very badly fitting profiles.
  // The bisection algorithm rejects an exact amount of data.
  
  if ((0) && (m_nIntRefs >= 200))
    {
      double fProfitSum,fIntSum;
      double fProfitIntensity,fIntIntensity;
      itr<int> anIntensitySort;
      itr<double> afIntensity;
      itr<double> afX;
      itr<double> afC;
      itr<double> afX2;
      itr<double> afAvg;
      itr<double> afDev;
      itr<double> afDiff;
      int nx;
      int nIntRefSort;
      const int nMinRefsPerBin = 200;
      
      // Sort based on intensity.
      anIntensitySort.setsize(m_nIntRefs);
      for (nx = 0; nx < m_nIntRefs; nx++)
        anIntensitySort[nx] = nx;
      for (nIntRef = 0; nIntRef < m_nIntRefs; nIntRef++)
        {
          afIntensity + (double) (*m_poIntegrate)[nIntRef].fGetIntensity();
        }
      afDiff.setsize(afIntensity.size());
      g_pfCmpDoubles = &afIntensity[0];
      qsort(&anIntensitySort[0],m_nIntRefs,sizeof(int),double_cmp_rel);

      afX.setsize(max(1,m_nIntRefs/nMinRefsPerBin));
      afX2.setsize(max(1,m_nIntRefs/nMinRefsPerBin));
      afC.setsize(max(1,m_nIntRefs/nMinRefsPerBin));
      for (nx = 0; nx < afX.size(); nx++)
        {
          afX[nx] = 0.0;
          afX2[nx] = 0.0;
          afC[nx] = 0.0;
        }

      for (nIntRefSort = 0; nIntRefSort < m_nIntRefs; nIntRefSort++)
        {
          nIntRef = anIntensitySort[nIntRefSort];
          fProfitSum = 0.0;
          fIntSum = 0.0;
          for (nPartialRef = m_anPartialStart[nIntRef]; nPartialRef <= m_anPartialEnd[nIntRef]; nPartialRef++)
            {
              fProfitIntensity = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fProfitIntensity);
              fIntIntensity = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fIntensity);
              if (fProfitIntensity != 0.0)
                {
                  fProfitSum += fProfitIntensity;
                  fIntSum += fIntIntensity;
                } 
              else
                break;
            }
          if (nPartialRef > m_anPartialEnd[nIntRef])
            {
              f0 = fIntSum - fProfitSum;
              nx = min(afX.size() - 1, nIntRefSort/nMinRefsPerBin);
              afX[nx] += f0;
              afX2[nx] += f0*f0;
              afC[nx] += 1.0;
              afDiff[nIntRef] = f0;
            } 
          else
            afDiff[nIntRef] = 0.0;
        }

      for (nx = 0; nx < afX.size(); nx++)
        {
          if (afC[nx] >= 20)
            {
              afAvg[nx] = afX[nx]/afC[nx];
              afDev[nx] = sqrt((afX2[nx]*afC[nx] - afX[nx]*afX[nx])/(afC[nx]*(afC[nx]-1)));
            } 
          else
            {
              afDev[nx] = 0.0;
              afAvg[nx] = 0.0;
            }
        }     

      int nBisectionStat = 0;
      int nBisectionState = 0;
      int nRefsExamined;
      double fSigmaLevel;
      double fFractionRejected = 0.0;
      double fTargetFractionRejected = 0.01;

      do
        {
          nBisectionStat = nBisection(nBisectionState,
                                      2.0,
                                      10.0,
                                      fSigmaLevel,
                                      fFractionRejected,
                                      fTargetFractionRejected,
                                      fTargetFractionRejected*0.2);
            
          nRefsRemovedBadProfitFit = 0;
          nRefsExamined = 0;
          for (nIntRefSort = 0; nIntRefSort < m_nIntRefs; nIntRefSort++)
            {
              nIntRef = anIntensitySort[nIntRefSort];
              nx = min(afX.size() - 1, nIntRefSort/nMinRefsPerBin);
              if ((afDev[nx]>0.0) && (afDiff[nx]!=0.0))
                {
                  if (ABS(afDiff[nIntRef] - afAvg[nx]) > fSigmaLevel*afDev[nx]) {
                    nRefsRemovedBadProfitFit++;                    
                    if (nBisectionStat)
                      (*m_poIntegrate)[nIntRef].m_nSelect = 1;
                  } 
                  else
                    nRefsExamined++;
                }
            }
          fFractionRejected = ((double) nRefsRemovedBadProfitFit)/nRefsExamined;
        } while (!nBisectionStat);
    }

  // Determine the level at which profile fitting should be attempted.

  if ((0) && (nPartialRefs > 100))
    {
      itr<double> afIntensity;
      itr<double> afIntensityCutoff;
      itr<double> afAverageVariance;
      itr<double> afAverageProfitVariance;
      itr<double> afAverageIntVariance;
      int nMinRefsPerBin = 500;

      double fSumIntProfitVar;
      int    nSumIntProfitVar;
      double fSumIntVar;
      int    nSumIntVar;
      double fIntIntensity;
      double fProfitIntensity;
      double fIntSigma;
      itr<double> afVarPlotIntensity;
      itr<double> afVarPlotDiff;


      // Sort based on intensity.
      afIntensity.setsize(nPartialRefs);
      -afIntensity;
      for (nPartialRef = 0;  nPartialRef < nPartialRefs; nPartialRef++)
        {
          afIntensity + (double) (*m_poPartial)[nPartialRef].fGetIntensity();
        };
      
      if( afIntensity.size() > 0 ) // safety check
        qsort(&afIntensity[0],nPartialRefs,sizeof(double),double_cmp);

      if (nPartialRefs/nMinRefsPerBin>25)
        nMinRefsPerBin = nPartialRefs/25;
      // Using our sorted intensity, make cutoff bins for our calculations.
      nRef = nMinRefsPerBin;
      do
        {
          afIntensityCutoff + afIntensity[min(afIntensity.size() - 1,nRef)];
          nRef += nMinRefsPerBin;
        } while (nRef < afIntensity.size());
      afIntensityCutoff.last() = afIntensity[afIntensity.size() - 1];


      for (nBin = 0; nBin < afIntensityCutoff.size(); nBin++)
        {
          if (m_anPartialStart[nIntRef] == -1)
            {
              continue;
            }
          fSumIntProfitVar = 0.0;
          nSumIntProfitVar = 0;
          fSumIntVar = 0.0;
          nSumIntVar = 0;
          for (nIntRef = 0; nIntRef < m_nIntRefs; nIntRef++)
            {
              for (nPartialRef = m_anPartialStart[nIntRef]; nPartialRef <= m_anPartialEnd[nIntRef]; nPartialRef++)
                {
                  fProfitIntensity = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fProfitIntensity);
                  fIntIntensity = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fIntensity);
                  fIntSigma = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fSigmaI);
                    
                  if (fIntIntensity<=afIntensityCutoff[nBin])
                    {
                      if ((nBin == 0) || (fIntIntensity > afIntensityCutoff[nBin-1]))
                        {
                          fSumIntProfitVar += (fProfitIntensity - fIntIntensity)*(fProfitIntensity - fIntIntensity);
                          nSumIntProfitVar ++;
                          fSumIntVar += fIntSigma*fIntSigma;
                          nSumIntVar ++;
                        }
                    }
                }
            }
          if (nSumIntVar >= nMinRefsPerBin*0.5)
            {
              afAverageIntVariance + (fSumIntVar/nSumIntVar);
              afAverageProfitVariance + ((fSumIntProfitVar - fSumIntVar)/nSumIntProfitVar);
            } 
          else
            {
              afAverageIntVariance + 0.0;
              afAverageProfitVariance + 0.0;
            }
        }

        // Determine the maximum intensity to profile fit.
        for (nBin = afAverageIntVariance.size() - 1; nBin >= 0; nBin--)
          {
            if (afAverageIntVariance[nBin] >= afAverageProfitVariance[nBin])
              break;
          }
        if (nBin>=0)
          m_fMaxIntensityToProfit = afIntensityCutoff[nBin];
        else
          m_fMaxIntensityToProfit = -10000.0;
        
        if ("" != sGetEnv("DTREK_VERBOSE"))
          {
            printf("\nMaximum Intensity to profile fit = %.2lf\n\n",m_fMaxIntensityToProfit);
            
            printf("Estimated profit/non-profit variances vs Intensity\n");
            printf("-----------------------------------------\n");
            printf(" %10s %10s %13s\n","Intensity","sigma(int)","sigma(profit)");
            printf("-----------------------------------------\n");
            for (nBin = 0; nBin < afAverageIntVariance.size(); nBin++)
              {
                printf("%10.2lf %10.2lf %13.2lf\n",(double) afIntensityCutoff[nBin],(double) sqrt(max(0.0,afAverageIntVariance[nBin])),(double) sqrt(max(0.0,afAverageProfitVariance[nBin])));
              }
            printf("-----------------------------------------\n");
          }
    }

    
    const int nIntensityBins = 7;
    const double fMaxBin = 70000.0;
    double afIntensityMax[nIntensityBins];
    double afIntensitySigma[nIntensityBins];
    double afIoverSigma[nIntensityBins];
    int    anIntensitySigma[nIntensityBins];
    int    anIntensityProfited[nIntensityBins];
    int    anIntensityNotProfited[nIntensityBins];

    for (nBin = 0; nBin < nIntensityBins; nBin++) {
        afIntensityMax[nBin] = fMaxBin*exp(-(double)(nIntensityBins - 2 - nBin));
        afIntensitySigma[nBin] = 0.0;
        afIoverSigma[nBin] = 0.0;
        anIntensitySigma[nBin] = 0;
        anIntensityProfited[nBin] = 0;
        anIntensityNotProfited[nBin] = 0;
    };

        const bool bDebug = false;
        FILE* pFDebugOut = (bDebug)?fopen("noise.txt","wt"):NULL;
        double fAverage;
        double fDeviation;


        for (nPass = 0; nPass < 2; nPass++)
          {

            oStat.vClearRaw();
            nProfitRefsUsed = 0;
            nProfitRefsExamined = 0;

            for (nIntRef = 0; nIntRef < m_nIntRefs; nIntRef++)
              {
                double fBackground;
                double fSigma;
                double fIntensity;
                double fSigmaProfit;
                double fIntensityProfit;
                double fSigmaBackground;
                double fStatistic;
                double fSumIntegrated;
                double fSumVariance;
                double fSumBackgroundVariance;
                int    nSumVariance;
                int    nSumIntegrated;
                int    nSumProfit;
                int    nSumOutsideMosaicity;
                int    nSumOutsideCentroid;
                int    nPartialStat;
                bool   bUsedProfit;
                bool   bIsTooWeak;
                double fObsRotMid;
                double fMinObsRot,fMaxObsRot,fMinObsRotWeak,fMaxObsRotWeak;
                double fObsRotWidth,fObsRotWidthWeak;
                int    nCountObsRot,nCountObsRotWeak;

                fSumIntegrated = 0.0;
                fSumVariance = 0.0;
                fSumBackgroundVariance = 0.0;
                nSumVariance = 0;
                nSumIntegrated = 0;
                nSumProfit = 0;
                nSumOutsideMosaicity = 0;
                nSumOutsideCentroid = 0;
                fMinObsRot = fMaxObsRot = fMinObsRotWeak = fMaxObsRotWeak = -9999.0;
                fObsRotWidth = fObsRotWidthWeak = 0.0;
                nCountObsRot = nCountObsRotWeak = 0;
                        
                        
                fBackground = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fBackground);

                        
                if (m_anPartialStart[nIntRef] == -1)
                  {
                    nRefsWithNoPartialInformation++;
                    continue;
                  };

                if ((*m_poIntegrate)[nIntRef].m_nSelect)
                  {
                    // Clear all the partials for this reflection since it has been rejected.
                    for (nPartialRef = m_anPartialStart[nIntRef]; nPartialRef <= m_anPartialEnd[nIntRef]; nPartialRef++)
                      {
                        (*m_poPartial)[nPartialRef].m_nSelect = 1;
                      }
                    continue;
                  }
                        
                for (nPartialRef = m_anPartialStart[nIntRef]; nPartialRef <= m_anPartialEnd[nIntRef]; nPartialRef++)
                  {
                    fIntensityProfit = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fProfitIntensity);
                    fSigmaProfit = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fProfitSigmaI);
                    fIntensity = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fIntensity);
                    fSigma = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fSigmaI);
                    fSigmaBackground = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fBackgroundSigma);
                    nPartialStat = (*m_poPartial)[nPartialRef].nGetField(m_poPartial->m_nFI_nPartialFlag);
                    (*m_poPartial)[nPartialRef].m_nSelect = 1;


                    if (DONTUSEPROFIT)
                      {
                        bUsedProfit = false;
                      } 
                    else
                      {
                        fIntensity = fIntensityProfit;
                        fSigma = fSigmaProfit;

                        nProfitRefsExamined++;
                        bUsedProfit = true;
                      }

                    fStatistic = fIntensity/max(1.0,fBackground);
                    if (nPass == 0)
                      {
                        if (fIntensity <= 0.0)
                          {
                            oStat.vAddRaw(((nRef % 2)*2 - 1)*fStatistic);
                            if (pFDebugOut) 
                              fprintf(pFDebugOut,"%d %lf\n",nPartialRef,fStatistic);
                          }
                      } 
                    else if (nPass == 1) 
                      {

                        //bIsTooWeak = (fStatistic < 2.5*fDeviation + fAverage);
                        bIsTooWeak = true;
                        if( 0 == (nPartialStat & 
                                 (C3Ddata::ms_nPerSlicePostRefineNotInStatisticalMosaicity | 
                                 C3Ddata::ms_nPerSlicePostRefineNotInRmergeMosaicity)) )
                            bIsTooWeak = false;

                        f0 = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fObsRotMid);
                        if (fMinObsRot == -9999.0)
                          {
                            fMinObsRot = fMaxObsRot = f0;
                            nCountObsRot = 1;
                          } 
                        else
                          {
                            fMinObsRot = min(f0,fMinObsRot);
                            fMaxObsRot = max(f0,fMaxObsRot);
                            nCountObsRot++;
                          }
                        if ((!bIsTooWeak) && (fMinObsRotWeak == -9999.0))
                          {
                            fMinObsRotWeak = fMaxObsRotWeak = (*m_poPartial)[nPartialRef].fGetField(m_poPartial->m_nFI_fObsRotMid);
                            nCountObsRotWeak = 1;
                          } 
                        else if (!bIsTooWeak)
                          {
                            fMinObsRotWeak = min(f0,fMinObsRotWeak);
                            fMaxObsRotWeak = max(f0,fMaxObsRotWeak);
                            nCountObsRotWeak++;
                          }


                        if (   (nPartialStat & C3Ddata::ms_nPerSliceCentroid)
                            || !bIsTooWeak) 
                          {
                            fSumIntegrated += fIntensity;
                            fSumVariance += fSigma*fSigma;
                            fSumBackgroundVariance += fSigmaBackground*fSigmaBackground;
                            nSumVariance++;
                            (*m_poPartial)[nPartialRef].m_nSelect = 0;
                        
                            for (nBin = 0; (afIntensityMax[nBin]< fIntensity)
                                            && (nBin + 1 < nIntensityBins);
                                 nBin++);

                            afIntensitySigma[nBin] += fSigma;
                            afIoverSigma[nBin] += fIntensity/max(1.0,fSigma);
                            anIntensitySigma[nBin]++;
                            anIntensityProfited[nBin] += (bUsedProfit);
                            anIntensityNotProfited[nBin] += (!bUsedProfit);
                            

                            if (bUsedProfit) 
                              nSumProfit++;
                            else 
                              nSumIntegrated++;
                            if (!(nPartialStat & C3Ddata::ms_nPerSliceInMosaicity))
                              nSumOutsideMosaicity++;
                            if (!(nPartialStat & C3Ddata::ms_nPerSliceCentroid))
                              nSumOutsideCentroid++;
                          }
                      }
                  }
                if (nPass == 1)
                  {
                    if (nSumVariance != (m_anPartialEnd[nIntRef] - m_anPartialStart[nIntRef]+1))
                      nRefsAdjustedForNoise++;
                    if (nSumIntegrated && nSumProfit)
                      nRefsWithMixedProfitIntegratedIntensity++;
                    if (nSumIntegrated && !nSumProfit)
                      nRefsWithOnlyIntegratedIntensity++;
                    if (!nSumIntegrated && nSumProfit)
                      nRefsWithOnlyProfitIntensity++;
                    if (nSumOutsideMosaicity)
                      nRefsWithIntensityBeyondMosaicityBounds++;
                    if (!nSumOutsideCentroid)
                      nRefsWithIntensityOnlyAtCentroid++;
                    nProfitRefsUsed += nSumProfit;

                    (*m_poIntegrate)[nIntRef].vSetField(0,(float) fSumIntegrated);
                    (*m_poIntegrate)[nIntRef].vSetField(1,(float) sqrt(fSumVariance));
                    (*m_poIntegrate)[nIntRef].vSetField(m_poIntegrate->m_nFI_fBackgroundSigma,(float) sqrt(fSumBackgroundVariance));


                    // Adjust the rotation width to represent the truncated intensity.  This is important because we
                    // want to correctly plot ellipsoids.
                    if (nCountObsRotWeak < 2)
                      (*m_poIntegrate)[nIntRef].vSetField(m_poIntegrate->m_nFI_fObsRotWidth,(float) 0.00);
                    else
                      {
                        fObsRotMid = (*m_poIntegrate)[nIntRef].fGetField(m_poIntegrate->m_nFI_fObsRotMid);
                        fObsRotWidth = (fMaxObsRot - fMinObsRot)/(nCountObsRot - 1);
                        fObsRotWidthWeak = (fMaxObsRotWeak - fMinObsRotWeak)/(nCountObsRotWeak - 1);
                        fMinObsRot -= fObsRotWidth/2.0;
                        fMaxObsRot += fObsRotWidth/2.0;
                        fMinObsRotWeak -= fObsRotWidth/2.0;
                        fMaxObsRotWeak += fObsRotWidth/2.0;

                        fObsRotWidth = 2.0*min(ABS(fObsRotMid - fMinObsRot),ABS(fObsRotMid - fMaxObsRot));
                        fObsRotWidthWeak = 2.0*min(ABS(fObsRotMid - fMinObsRotWeak),ABS(fObsRotMid - fMaxObsRotWeak));

                        (*m_poIntegrate)[nIntRef].vSetField(m_poIntegrate->m_nFI_fObsRotWidth,(float) fObsRotWidthWeak);
                      }
                  }
              }
            fAverage = oStat.fAverageRaw();
            fDeviation = oStat.fStandardDeviationRaw();


            if (pFDebugOut)
              {
                fclose(pFDebugOut);
                pFDebugOut = NULL;
              }
          }

        for (nIntRef = 0; nIntRef < m_nIntRefs; nIntRef++)
          {
            if ((*m_poIntegrate)[nIntRef].m_nSelect == 0)
              nRefsForOutput++;
          }
    
        if ("" != sGetEnv("DTREK_VERBOSE"))
          {
            printf("Analysis of profile fitted and non-profile fitted data\n");
            printf("-------------------------------------------------------------------------------\n");
        
            printf("Intensity");
            for (nBin = 0; nBin < nIntensityBins; nBin++)
              {
                char pcBuf[50];
                sprintf(pcBuf,"%2s%-.0lf",(nBin == nIntensityBins-1)?" >":"<=",(afIntensityMax[min(nBin,nIntensityBins-2)]));
                printf("%9s ",pcBuf);
              }
            printf("\n");
        
            printf("<Sigma>  ");
            for (nBin = 0; nBin < nIntensityBins; nBin++)
              {
                printf("  %7.0lf ",afIntensitySigma[nBin]/max(1,anIntensitySigma[nBin]));
              }
            printf("\n");
            printf("<I/Sigma>");
            for (nBin = 0; nBin < nIntensityBins; nBin++)
              {
                printf("  %7.1lf ",afIoverSigma[nBin]/max(1,anIntensitySigma[nBin]));
              }
            printf("\n");    
            printf("#Profited");
            for (nBin = 0; nBin < nIntensityBins; nBin++)
              {
                printf("  %7d ",anIntensityProfited[nBin]);
              }
            printf("\n");    
            printf("#Integ.  ");
            for (nBin = 0; nBin < nIntensityBins; nBin++)
              {
                printf("  %7d ",anIntensityNotProfited[nBin]);
              }
            printf("\n");    
            printf("-------------------------------------------------------------------------------\n");
            printf("\n\n");
          }

        if (bFirstPart)
          {
            s_nRefsRemovedBadPosition                   = 0;
            s_nRefsRemovedBadProfitFit                  = 0;
            s_nRefsRemovedOverlap                       = 0;// Previous removal.
            s_nRefsForOutput                            = 0;
            s_nRefsAdjustedForNoise                     = 0;
            s_nRefsWithMixedProfitIntegratedIntensity   = 0;
            s_nRefsWithOnlyProfitIntensity              = 0;
            s_nRefsWithOnlyIntegratedIntensity          = 0;
            s_nRefsWithNoPartialInformation             = 0;
            s_nRefsWithIntensityBeyondMosaicityBounds   = 0;
            s_nRefsWithIntensityOnlyAtCentroid          = 0;
            s_nProfitRefsUsed                           = 0;
            s_nProfitRefsExamined                       = 0;
          } 
        s_nRefsRemovedBadPosition     += nRefsRemovedBadPosition;
        s_nRefsRemovedBadProfitFit    += nRefsRemovedBadProfitFit;
        s_nRefsRemovedOverlap         += nRefsRemovedOverlap;
        s_nRefsForOutput              += nRefsForOutput;
        s_nRefsAdjustedForNoise       += nRefsAdjustedForNoise;
        s_nRefsWithMixedProfitIntegratedIntensity += nRefsWithMixedProfitIntegratedIntensity;
        s_nRefsWithOnlyProfitIntensity += nRefsWithOnlyProfitIntensity;
        s_nRefsWithOnlyIntegratedIntensity += nRefsWithOnlyIntegratedIntensity;
        s_nRefsWithNoPartialInformation       += nRefsWithNoPartialInformation;
        s_nRefsWithIntensityBeyondMosaicityBounds += nRefsWithIntensityBeyondMosaicityBounds;
        s_nRefsWithIntensityOnlyAtCentroid        += nRefsWithIntensityOnlyAtCentroid;
        s_nProfitRefsUsed                         += nProfitRefsUsed;
        s_nProfitRefsExamined                     += nProfitRefsExamined;
    

        if (bLastPart)
          {
            // Print statistics gleaned from merging the data.
            printf("Analysis of partials contributing to integrated data\n");
            printf("-----------------------------------------------------------------------------\n");
            printf("Reflections removed with highly deviant pixel positions:             %8d\n",s_nRefsRemovedBadPosition);
            printf("Reflections removed due to highly deviant profile fit:               %8d\n",s_nRefsRemovedBadProfitFit);
            printf("Reflections removed due to overlaps:                                 %8d\n",s_nRefsRemovedOverlap);
            printf("Maximum deviation allowed between calculated and observed positions: %8.1lf\n",fPxAverage + fPxDeviationSigma*fPxDeviation);
            printf("Reflections with noisy tails removed:                                %8d\n",s_nRefsAdjustedForNoise);
            printf("Reflections containing both profiled and non-profiled data:          %8d\n",s_nRefsWithMixedProfitIntegratedIntensity);
            printf("Reflections containing only profiled data:                           %8d\n",s_nRefsWithOnlyProfitIntensity);
            printf("Reflections containing only non-profiled data:                       %8d\n",s_nRefsWithOnlyIntegratedIntensity);
            printf("Reflections which lacked adequate profile information:               %8d\n",s_nRefsWithNoPartialInformation);
            printf("Reflections containing intensity beyond predicted mosaicity bounds:  %8d\n",s_nRefsWithIntensityBeyondMosaicityBounds);
            printf("Reflections containing intensity only in the centroid slice:         %8d\n",s_nRefsWithIntensityOnlyAtCentroid);
            printf("-----------------------------------------------------------------------------\n");
            printf("Total profiles examined:                                             %8d\n",s_nProfitRefsExamined);
            printf("Total profiles used:                                                 %8d\n",s_nProfitRefsUsed);
            printf("Total reflections to be output:                                      %8d\n",s_nRefsForOutput);
            printf("-----------------------------------------------------------------------------\n");    
          }

        return 0;
}

int Cprofit2DSlices::nPrintAxial(bool bPrint)
{
    // Do some things with oAxialReflnlist
    Crefln*     poRefln;
    int*        pnHKL;
    int         nLoop;
    int         nRef;
    
    char*       pcSelect = NULL;

    int         nIntRefs = m_poIntegrate->nGetNumReflns();
    
    // Reset all the deletion flags.
    for (nRef = 0; nRef < nIntRefs; nRef++)
      (*m_poIntegrate)[nRef].m_nSelect = 1;
    // Reset all the deletion flags.
    for (nRef = 0; nRef < nIntRefs; nRef++)
      (*m_poIntegrate)[nRef].m_nSelect = 0;


    pcSelect = new char [m_poIntegrate->m_nTotalFieldsPlus1];
    memset(pcSelect,0,(m_poIntegrate->m_nTotalFieldsPlus1));    
    
    // Reset all the deletion flags.
    for (nRef = 0; nRef < nIntRefs; nRef++)
      (*m_poIntegrate)[nRef].m_nSelect = 1;


    for (nLoop=0;nLoop<3;nLoop++)
      {
        for (nRef = 0; nRef < m_poIntegrate->nGetNumReflns(); nRef++)
          {
            poRefln = m_poIntegrate->poGetRefln(nRef);
            pnHKL = poRefln->pnGetHKL();
            if ((pnHKL[(nLoop+1) % 3]==0) && (pnHKL[(nLoop+2) % 3]==0))
              {
                poRefln->m_nSelect = 0;
              } 
          }
      } 
    if (!m_poAxialRefs)
      {
        m_poAxialRefs = new Creflnlist(*m_poIntegrate);
      }
    m_poAxialRefs->nInsertListFrom(*m_poIntegrate,pcSelect);
    // Reset all the deletion flags.
    for (nRef = 0; nRef < nIntRefs; nRef++)
      (*m_poIntegrate)[nRef].m_nSelect = 0;



    if (bPrint)
      {
        int nx;
        Cstring sTemp;
	int nPrinted = 0;
        
        printf("\nh00, 0k0, 00l (axial) reflections. Examine them for systematic absences."
            "\n--------------------------------------------------------------"
            "\n    h     k     l Intensity SigmaInt I/sigI  Pix0  Pix1  Angle"
            "\n--------------------------------------------------------------\n");
        
        for (nLoop=0; nLoop < 3; nLoop++)
          {
            for (nx = 0; nx < m_poAxialRefs->nGetNumReflns(); nx++)
              {
                
                poRefln = m_poAxialRefs->poGetRefln(nx);
                pnHKL = poRefln->pnGetHKL();
                if ((pnHKL[(nLoop+1) % 3]==0) && (pnHKL[(nLoop+2) % 3]==0))
                  {
                    printf("%5d %5d %5d %9.0f %8.0f %6.1f %5.0f %5.0f %6.2f\n",
                           poRefln->nGetH(),
                           poRefln->nGetK(),
                           poRefln->nGetL(),
                           poRefln->fGetIntensity(),
                           poRefln->fGetSigmaI(),
                           poRefln->fGetIntensity() / max(1.0,poRefln->fGetSigmaI()),
                           poRefln->fGetField(m_poIntegrate->m_nFI_fObsPx0),
                           poRefln->fGetField(m_poIntegrate->m_nFI_fObsPx1),
                           poRefln->fGetField(m_poIntegrate->m_nFI_fObsRotMid));
		    nPrinted++;
                  }
              }
          }
        printf("--------------------------------------------------------------\n");        
        if (0 < nPrinted)
	  {
	    sTemp = sDtrekGetPrefix()+"axial.ref";
	    m_poAxialRefs->nWrite(sTemp);
	    printf("INFO: Use dtdisplay and File>Read refln list ... to view these axial reflns.\n\n");
	  }
	else
	  {
	    printf("INFO: There were no axial reflns in the integrated results.\n\n");
	  }
      }
    delete [] pcSelect;
    return 0;
}


int Cprofit2DSlices::nWrite(int nPart,
                            bool bWriteIntegrate,
                            bool bWriteProfit,
                            bool bWriteEllipsoid)
{
  char* pcSelect;

  pcSelect = new char [m_poIntegrate->m_nTotalFieldsPlus1];
  memset(pcSelect,0,(m_poIntegrate->m_nTotalFieldsPlus1));

  if ((bWriteIntegrate) && (m_poIntegrate))
    m_poIntegrate->nWrite(sCreatePartName(m_sIntegrate,nPart),NULL,pcSelect);
  if ((bWriteProfit) && (m_poPartial))
    m_poPartial->nWrite(sCreatePartName(m_sPartial,nPart));

  delete [] pcSelect;

  if (bWriteEllipsoid)
    {
      int nx,ny;
      (void) m_poPartial->nExpandGetField(m_poPartial->ms_sfEllipseMajorAxis);
      (void) m_poPartial->nExpandGetField(m_poPartial->ms_sfEllipseMinorAxis);
      (void) m_poPartial->nExpandGetField(m_poPartial->ms_sfEllipseAxisMajorOffset);
      (void) m_poPartial->nExpandGetField(m_poPartial->ms_sfObsPx0);
      (void) m_poPartial->nExpandGetField(m_poPartial->ms_sfObsPx1);
      (void) m_poPartial->nExpandGetField(m_poPartial->ms_sfObsRotWidth);
      pcSelect = new char [m_poPartial->m_nTotalFieldsPlus1];
      memset(pcSelect,0,(m_poPartial->m_nTotalFieldsPlus1));

      // The partial's list was lacking the ObsPx0 and ObsPx1 fields.  Set these.
      for (nx=0;nx<m_poIntegrate->nGetNumReflns();nx++)
        {
          if (m_anPartialStart[nx]>=0)
            {
              for (ny = m_anPartialStart[nx];ny <= m_anPartialEnd[nx]; ny++)
                {
                  (*m_poPartial)[ny].vSetField(m_poPartial->m_nFI_fObsPx0,(*m_poIntegrate)[nx].fGetField(m_poIntegrate->m_nFI_fObsPx0));
                  (*m_poPartial)[ny].vSetField(m_poPartial->m_nFI_fObsPx1,(*m_poIntegrate)[nx].fGetField(m_poIntegrate->m_nFI_fObsPx1));
                }
            }
        }
        
      for (nx=0;nx<m_poPartial->nGetNumReflns();nx++)
        {
          (*m_poPartial)[nx].nConvertToAxialEllipse();
          (*m_poPartial)[nx].vSetField(m_poPartial->m_nFI_fObsRotWidth,(float) 0.01);
        }
      m_poPartial->nWrite(sCreatePartName(m_sEllipsoid,nPart),NULL,pcSelect);
        
      delete [] pcSelect;
    }
  return 0;
}

int Cprofit2DSlices::nSyncLists(bool bSyncWithHKL)
{
  int* pnPartialStart;
  int* pnPartialEnd;
  int nRef,nRefPartial,nRefScale;
  int nStat;

  if ((!m_poIntegrate) || (!m_poPartial))
    return 1;

  m_anPartialStart.setsize(m_poIntegrate->nGetNumReflns());
  pnPartialStart = &m_anPartialStart[0];
  m_anPartialEnd.setsize(m_poIntegrate->nGetNumReflns());
  pnPartialEnd = &m_anPartialEnd[0];

  // Set all reflections so that they are not pointing to any list of partials.
  for (nRef = 0; nRef < m_poIntegrate->nGetNumReflns(); nRef++)
    {
      pnPartialStart[nRef] = -1;
      pnPartialEnd[nRef] = -1;
    }
  int nRefs;
  int nRefsNotFound = 0;
  for (nRef = 0, nRefPartial = 0; (nRef < m_poIntegrate->nGetNumReflns()); nRef++)
    {
      nRefs = 0;
      m_anPartialStart[nRef] = nRefPartial;
      if (nRefPartial >=m_poPartial->nGetNumReflns())
        break;
      while (   (nRefPartial < m_poPartial->nGetNumReflns())
             && ((*m_poIntegrate)[nRef].nPackHKL() == ((*m_poPartial)[nRefPartial].nPackHKL())))
        {
          nRefs++;
          nRefPartial++;
        }
      m_anPartialEnd[nRef] = nRefPartial-1;
      if (!nRefs)
        {
          nRefPartial++;
          nRef--;
          nRefsNotFound++;
        }                      
    }
  while (nRef < m_poIntegrate->nGetNumReflns())
    {
      m_anPartialStart[nRef] = -1;
      m_anPartialEnd[nRef] = -1;
      nRef++;
    }
  if (nRefsNotFound)
    printf("WARNING:  %d partial reflections did not have an associated integrated reflection.\n",nRefsNotFound);
    
  if (m_poScaleList)
    {
      // Build list of scale factors already computed for each reflection.
      // To do this, we must synchronize the scaled reflections with their dtintegrate counterparts.
        
      int nFI_fScale = m_poScaleList->nGetFieldIndex((Cstring) "fScale");
        
      if (nFI_fScale >= 0)
        {
          m_afPerReflnScale.setsize(m_poIntegrate->nGetNumReflns());
          for (nRef = 0; nRef < m_poIntegrate->nGetNumReflns(); nRef++) 
            m_afPerReflnScale[nRef] = 1.0;
          nRef = 0;
          nStat = 0;
          for (nRefScale = 0; nRefScale < m_poScaleList->nGetNumReflns(); nRefScale++)
            {
              if (nRef >=m_poIntegrate->nGetNumReflns())
                {
                  printf("ERROR:  Could not sync. scale file with integrate list.\n");
                  nStat = 1;
                  break;
                }
              for (; nRef < m_poIntegrate->nGetNumReflns(); nRef++)
                {
                  if ((*m_poIntegrate)[nRef].nPackHKL() == (*m_poScaleList)[nRefScale].nPackHKL())
                    {
                      m_afPerReflnScale[nRef] = (*m_poScaleList)[nRefScale].fGetField(nFI_fScale);
                      nRef++;
                      break;
                    } 
                  else
                    {
                      m_afPerReflnScale[nRef] = 1.0;
                    }
                }
            }

          if (nRefScale < m_poScaleList->nGetNumReflns())
            nStat = 1;
          // Free up the scale list.
          delete m_poScaleList;
          m_poScaleList = NULL;
          if (nStat)
            return nStat;
        } 
      else 
        {
          printf("ERROR:  The input scale list is invalid!\n");
          return 1;
        }
    }
  else
    m_afPerReflnScale.clear();

  return 0;
}


int Cprofit2DSlices::nCalcWidthsHKL(double fFixedMosaicity)
{
  int *pnSort;
  int nLastHKL,nThisHKL=0;
  int nStartIndex,nEndIndex;
  int nRef,nRefSort;
  int nRefIndex;
  int nResoGroup;
  int nPass;

  double f0,f1;
  int nx,ny;

  itr<double> afSliceMinMosaicity;        // One entry per slice for a given HKL group.  The minimum mosaicity required for the slice to be predicted.
  itr<int>    anSliceSort;                // One entry per slice for a given HKL group.  The index into m_poPartials
  itr<int>    anRefSort;                  // One entry per slice for a given HKL group.  The index into m_poIntegrate (relative to nStartIndex.  Thus it's >=0)
  itr<double> afRefIntensity;             // One entry per refln for a given HKL group.  Used in calculation of variance.
  itr<int>    anRefIntensityIsValid;      // One entry per refln for a given HKL group.  If we already scaled are data and are profile fitting a second time, this is used.

  const double fHKLMosaicityStep = 0.05;
  const double fHKLMosaicityMax = 10.0;
  itr<double> m_afHKLMosaicity[g_nMaxResoGroups];
  itr<double> m_afHKLMosaicityContrib[g_nMaxResoGroups];
  double m_afMosaicity[g_nMaxResoGroups];
  double m_afMosaicityMults[g_nMaxResoGroups];
  double m_afMosaicityContrib[g_nMaxResoGroups];

  double fExpansionFactor;
  int nPartialStart,nPartialEnd,nPartialBins,nPartialSlice;
  double fObsRotStart,fObsRotEnd,fObsRotWidth,fCalcRotMid;

  // Initialize the mosaicity histogram.
  for (nResoGroup = 0; nResoGroup < g_nMaxResoGroups ; nResoGroup++)
    { 
      m_afHKLMosaicity[nResoGroup].setsize((int) (fHKLMosaicityMax/fHKLMosaicityStep));
      m_afHKLMosaicity[nResoGroup].dupc(0);
      m_afHKLMosaicityContrib[nResoGroup].setsize((int) (fHKLMosaicityMax/fHKLMosaicityStep));
      m_afHKLMosaicityContrib[nResoGroup].dupc(0);
      m_afMosaicity[nResoGroup] = 0.0;
      m_afMosaicityMults[nResoGroup] = 0.0;
      m_afMosaicityContrib[nResoGroup] = 0.0;
    }

  if (m_poIntegrate->nReduce(*m_poCrystal, false)) 
    return 1;
  pnSort = m_poIntegrate->pnGetSortIndex();
    
  for (nPass = 0; nPass < 2; nPass++)
    {
      nLastHKL = (*m_poIntegrate)[pnSort[0]].nGetField(m_poIntegrate->m_nFI_nPackedHKL);
      nStartIndex = 0;
        
      for (nRefSort=0;nRefSort<=m_nIntRefs;nRefSort++)
        {
          nRef=pnSort[nRefSort];
            
          if (nRefSort<m_nIntRefs)
            nThisHKL = (*m_poIntegrate)[nRef].nGetField(m_poIntegrate->m_nFI_nPackedHKL);
          if ((nRefSort!=m_nIntRefs) && (nThisHKL==nLastHKL))
            continue;
          nEndIndex=nRefSort-1;
            
          if (nEndIndex >= nStartIndex)
            {
              nResoGroup = m_pnResoBatch[pnSort[nStartIndex]];
                
              -afSliceMinMosaicity;
              -anSliceSort;
              -anRefSort;
              // Go through each slice, and find the minimum mosaicity needed to produce that slice.
              for (nRefIndex = nStartIndex; nRefIndex <= nEndIndex; nRefIndex++)
                {
                  if ((m_anPartialStart[pnSort[nRefIndex]] < m_anPartialEnd[pnSort[nRefIndex]]) && (m_anPartialStart[pnSort[nRefIndex]]!=-1))
                    { 
                      fExpansionFactor = m_pfLorentz[pnSort[nRefIndex]]*m_pfDStar[pnSort[nRefIndex]];
                      nPartialStart = m_anPartialStart[pnSort[nRefIndex]];
                      nPartialEnd = m_anPartialEnd[pnSort[nRefIndex]];
                      nPartialBins = nPartialEnd - nPartialStart + 1;
                        
                      fObsRotStart = (*m_poPartial)[nPartialStart].fGetField(m_poPartial->m_nFI_fObsRotMid);
                      fObsRotEnd = (*m_poPartial)[nPartialEnd].fGetField(m_poPartial->m_nFI_fObsRotMid);
                      fObsRotWidth = (fObsRotEnd - fObsRotStart)/(nPartialBins - 1);
                      fCalcRotMid = (*m_poPartial)[nPartialEnd].fGetField(m_poPartial->m_nFI_fCalcRotMid);
                        
                      for (nPartialSlice = nPartialStart; nPartialSlice <= nPartialEnd; nPartialSlice++)
                        {
                          f0 = ((*m_poPartial)[nPartialSlice].fGetField(m_poPartial->m_nFI_fObsRotMid) - fObsRotWidth/2.0);
                          f1 = ((*m_poPartial)[nPartialSlice].fGetField(m_poPartial->m_nFI_fObsRotMid) + fObsRotWidth/2.0);
                          if ((fCalcRotMid >= f0) && (fCalcRotMid <= f1))
                            {
                              afSliceMinMosaicity + 0.0;
                              anSliceSort + nPartialSlice;
                              anRefSort + (nRefIndex - nStartIndex);
                            } 
                          else
                            {
                              f0 = min(ABS(fCalcRotMid - f0), ABS(fCalcRotMid - f1));
                              f0 /= fExpansionFactor;
                              afSliceMinMosaicity + f0;
                              anSliceSort + nPartialSlice;
                              anRefSort + (nRefIndex - nStartIndex);
                            }
                        }
                    }
                  else
                    {
                      afSliceMinMosaicity + 0.0;
                      anSliceSort + m_anPartialStart[pnSort[nRefIndex]];
                      anRefSort + (nRefIndex - nStartIndex);
                    }
                }
                
                if (nPass == 1)
                  {
                    int nPartialStat;
                    // Compute the refln. widths based on best R-merge.
                    
                    for (nx = 0; nx < anSliceSort.size(); nx++)
                      {
                        ny = anSliceSort[nx];
                        nPartialStat = (*m_poPartial)[ny].nGetField(m_poPartial->m_nFI_nPartialFlag);
                        // See if this slice is noise or not.
                        f0 = m_afMosaicity[nResoGroup];
                        if (fFixedMosaicity)
                          f0 = fFixedMosaicity;
                        if ((afSliceMinMosaicity[nx]>f0) && (afSliceMinMosaicity[nx]>0.0))
                          nPartialStat |= C3Ddata::ms_nPerSlicePostRefineNotInRmergeMosaicity;
                        else
                          nPartialStat &= ~C3Ddata::ms_nPerSlicePostRefineNotInRmergeMosaicity;
                        (*m_poPartial)[ny].vSetField(m_poPartial->m_nFI_nPartialFlag,nPartialStat);
                      }
                  } 
                else
                  {
                    
                    if (afSliceMinMosaicity.size() && (anRefSort.size()>=4))
                      {
                        // Sort the mosaicity slices.
                        -g_apvSwapPointers + ((void*) &afSliceMinMosaicity[0]) + ((void*) & anSliceSort[0]) + ((void*) & anRefSort[0]);
                        -g_anSwapSizes + sizeof(double) + sizeof(int) + sizeof(int);
                        qsortswap(&afSliceMinMosaicity[0],afSliceMinMosaicity.size(),sizeof(double),double_cmp,qsort_swap_arrays);
                        
                        // Fill in the values for each reflection as they would be if we didn't truncate anything.
                        -afRefIntensity;
                        -anRefIntensityIsValid;
                        afRefIntensity.setsize(nEndIndex - nStartIndex + 1);
                        anRefIntensityIsValid.setsize(nEndIndex - nStartIndex + 1);
                        afRefIntensity.dupc(0);
                        anRefIntensityIsValid.dupc(0);
                        for (nx = 0; nx < anSliceSort.size(); nx++)
                          {
                            if (!m_afPerReflnScale.size())
                              {
                                anRefIntensityIsValid[anRefSort[nx]] = 1;
                                afRefIntensity[anRefSort[nx]] += (*m_poPartial)[anSliceSort[nx]].fGetIntensity();
                              } 
                            else 
                              {
                                if (m_afPerReflnScale[pnSort[anRefSort[nx]+ nStartIndex]]>0.0) 
                                  anRefIntensityIsValid[anRefSort[nx]] = 1;                        
                                afRefIntensity[anRefSort[nx]] += (*m_poPartial)[anSliceSort[nx]].fGetIntensity()*m_afPerReflnScale[pnSort[anRefSort[nx] + nStartIndex]];
                              }
                          }

                        int nCumul;
                        int nReject;
                        int nBin;
                        double fCumul;
                        double fRmerge;
                        double fIntensityBefore,fIntensityAfter;
                        nCumul = 0;
                        fCumul = 0.0;
                        for (nx = 0; nx < nEndIndex - nStartIndex + 1; nx++)
                          {
                            f0 = afRefIntensity[nx];
                            if (anRefIntensityIsValid[nx])
                              {
                                nCumul++;
                                fCumul += f0;
                              }
                          }
                        if ((nCumul >= 2) && (fCumul>0.0))
                          {
                            // Now, we must iteratively compute the R-merge by successively rejecting slices.
                            // Start with the highest mosaicity and work downward.
                            for (nReject = anSliceSort.size() - 1; (afSliceMinMosaicity[nReject]>0.0); nReject--)
                              {
                                // Compute the current R-merge.
                                fRmerge = 0.0;
                                for (nx = 0; nx < nEndIndex - nStartIndex + 1; nx++)
                                  {
                                    if (anRefIntensityIsValid[nx])
                                      {
                                        fRmerge += ABS(fCumul/nCumul - afRefIntensity[nx]);
                                      }
                                  }
                                fRmerge /= fCumul;
                                // Put R-merge into the proper mosaicity bin.
                                nBin = ((int) (afSliceMinMosaicity[nReject]/fHKLMosaicityStep));
                                if ((fRmerge > 0.2) && (nReject == anSliceSort.size() - 1))
                                  break;
                                m_afHKLMosaicity[nResoGroup][nBin] += fRmerge;
                                m_afHKLMosaicityContrib[nResoGroup][nBin] ++;

                                // Now remove this slice.
                                fIntensityBefore = afRefIntensity[anRefSort[nReject]];
                                fIntensityAfter = fIntensityBefore - (*m_poPartial)[anSliceSort[nReject]].fGetIntensity();
                                fCumul += (fIntensityAfter - fIntensityBefore);
                                afRefIntensity[anRefSort[nReject]] = fIntensityAfter;
                              }
                            m_afMosaicityMults[nResoGroup]++;
                            m_afMosaicityContrib[nResoGroup]+= nCumul;                                
                          }
                      }
                  }
            }  
            
            
          nStartIndex=nRefSort;
          if (nRefSort<m_nIntRefs)
            {
              nLastHKL =  (*m_poIntegrate)[nRef].nGetField(m_poIntegrate->m_nFI_nPackedHKL);
            }
        }
    }

  Cstring sName;
  itr<double> afTemp;
  for (nResoGroup = 0; nResoGroup < g_nMaxResoGroups; nResoGroup++)
    {
      sName = "RmergevsResoPlot";
      sName += nResoGroup;
      sName += ".txt";
      -afTemp;
      for (nx = 0; nx < m_afHKLMosaicity[nResoGroup].size(); nx++)
        afTemp + (m_afHKLMosaicity[nResoGroup][nx]/max(1.0,m_afHKLMosaicityContrib[nResoGroup][nx]));        
      nPlot(sName.string(),afTemp,0,fHKLMosaicityStep);
    }

  // Compute the mosaicity that yields the minimum R-merge for each resolution shell.
  printf("Estimated mosaicity per resolution\n");
  printf("----------------------------------------------\n");
  printf("  Resolution     Num   Avg  Rmerge Mosaicity\n");
  printf("     range     mults  mult   shell     shell\n");
  printf("----------------------------------------------\n");
  for (nResoGroup = 0; nResoGroup < g_nMaxResoGroups; nResoGroup++)
    {
      double fMinRmerge;
        
      ny = 0;
      fMinRmerge = 1.0;
      for (nx = 0; nx < m_afHKLMosaicity[nResoGroup].size(); nx++)
        {
          if ((m_afHKLMosaicityContrib[nResoGroup][nx]>1) && ((f0 = m_afHKLMosaicity[nResoGroup][nx]/max(1.0,m_afHKLMosaicityContrib[nResoGroup][nx])) < fMinRmerge))
            {
              fMinRmerge = f0;
              ny = nx;
            }
        }
      m_afMosaicity[nResoGroup] = fHKLMosaicityStep*(ny + 1);
      printf(" %5.2lf - %4.2lf %6d %5.2lf %7.3lf %9.3lf\n",
             m_afShellResoGroupMin[nResoGroup],
             m_afShellResoGroupMax[nResoGroup],
             (int) m_afMosaicityMults[nResoGroup],
             (double) (m_afMosaicityContrib[nResoGroup]/max(1.0,m_afMosaicityMults[nResoGroup])),
             (double) (fMinRmerge),
             m_afMosaicity[nResoGroup]            
             );
    }
  printf("-----------------------------------------------\n");

  return 0;
}

double Cprofit2DSlices::fEvaluate(itr<double>& afCalcValues,
                                  itr<double>& afObsValues,
                                  itr<double>& afSigmaValues)
{
  double f0;
  int nPoint;

  double fExpansionFactor;
  double fShift;
  double fCalcRotMid;
  double fObsRotStart, fObsRotEnd, fObsRotWidth;
  double fProfileRot;
  double fProfileShift;
  double fTotalIntensity,fTotalFitIntensity;
  int    nRefineBatch;
  int    nProfileBatch,nProfileBatchi;
  int    nProfileIndex;
  int    nResoGroup;
  int    nIntRef;
  int    nPartialStart,nPartialEnd,nPartialBins;
  int    nBin;
  int    nPass;
  int    nx;
  int    nPartialStat;
  itr<double> afFractionContrib;
  double*     pfFractionContrib;


  double fFractionCounts;
  double fFractionCountsNotAccounted;

  int nSumChiSq;      // The overall contributors to the chi-square computation.
  double fSumChiSq;

  double a2fChiSq[2]; // Contributors in each of the two test profile batches.
  int    a2nProfileBatch[2];
  int    nBatchesTested;

  // Initialize all profiles to 0.0

  m_afStandardProfileYNext.setsize(m_nProfilePoints*m_nNumProfiles);
  memset(&m_afStandardProfileYNext[0],0,sizeof(double)*m_nNumProfiles*m_nProfilePoints);
  m_pfStandardProfileYNext = & m_afStandardProfileYNext[0];
  m_afStandardProfileWidthY.setsize(m_nProfilePoints*m_nResoGroups);
  memset(&m_afStandardProfileWidthY[0],0,sizeof(double)*m_nProfilePoints*m_nResoGroups);
  m_pfStandardProfileWidthY = &m_afStandardProfileWidthY[0];

  // Initialize all per-refinement batch parameters.
  m_afRefineBatchChiSq.setsize(m_anRefineBatchStart.size());
  m_pfRefineBatchChiSq = & m_afRefineBatchChiSq[0];
  memset(m_pfRefineBatchChiSq,0,sizeof(double)*m_anRefineBatchStart.size());
  m_anRefineBatchContrib.setsize(m_anRefineBatchStart.size());
  m_pnRefineBatchContrib = & m_anRefineBatchContrib[0];
  memset(m_pnRefineBatchContrib,0,sizeof(int)*m_anRefineBatchStart.size());
  m_afRefineBatchObsShift.setsize(m_anRefineBatchStart.size());
  m_pfRefineBatchObsShift = & m_afRefineBatchObsShift[0];
  memset(m_pfRefineBatchObsShift,0,sizeof(double)*m_anRefineBatchStart.size());
  m_afRefineBatchContrib.setsize(m_anRefineBatchStart.size());
  m_pfRefineBatchContrib = & m_afRefineBatchContrib[0];
  memset(m_pfRefineBatchContrib,0,sizeof(double)*m_anRefineBatchStart.size());

  fSumChiSq = 0.0;
  nRefineBatch = 0;
  nSumChiSq = 0;
  -afObsValues;
  -afCalcValues;
  -afSigmaValues;

  for (nx = 0; nx < m_poPartial->nGetNumReflns(); nx++)
    {
      nPartialStat = (*m_poPartial)[nx].nGetField(m_poPartial->m_nFI_nPartialFlag);
      nPartialStat = nPartialStat & ~(C3Ddata::ms_nPerSlicePostRefineNotInStatisticalMosaicity);
      (*m_poPartial)[nx].vSetField(m_poPartial->m_nFI_nPartialFlag,nPartialStat);
    }
    
  // Loop through every reflection that is covered by this profile.
  for (nIntRef = 0; nIntRef < m_nIntRefs; nIntRef++)
    {
        // Update to the correct refinement batch.
      while (    (nRefineBatch < m_anRefineBatchEnd.size())
              && ((nIntRef > m_anRefineBatchEnd[nRefineBatch]) 
              || (nIntRef < m_anRefineBatchStart[nRefineBatch])) )
        nRefineBatch++;
      if (nRefineBatch >= m_anRefineBatchEnd.size())
        return 0.0;
      fShift = m_pfRefineParamShift[nRefineBatch];
      nResoGroup = m_pnResoBatch[nIntRef];
      nBatchesTested = 0;
      fProfileShift = 0.0;
      fTotalIntensity = 0.0;

        
      for (nProfileBatchi = 0; nProfileBatchi < 2; nProfileBatchi++)
        {
          nProfileBatch = m_a2pnProfileBatch[nProfileBatchi][nIntRef];
          if (nProfileBatch == -1)
            continue;
          nProfileIndex = m_anProfileBatchIndex[nResoGroup][nProfileBatch];
          
          // Must have at least two slices to do this.
          if (    (m_anPartialStart[nIntRef] < m_anPartialEnd[nIntRef])
               && (m_anPartialStart[nIntRef]!=-1))
            {
              nPartialStart = m_anPartialStart[nIntRef];
              nPartialEnd = m_anPartialEnd[nIntRef];
              nPartialBins = nPartialEnd - nPartialStart + 1;
              
              fObsRotStart = (*m_poPartial)[nPartialStart].fGetField(m_poPartial->m_nFI_fObsRotMid);
              fObsRotEnd = (*m_poPartial)[nPartialEnd].fGetField(m_poPartial->m_nFI_fObsRotMid);
              fObsRotWidth = (fObsRotEnd - fObsRotStart)/(nPartialBins - 1);
              fCalcRotMid = (*m_poPartial)[nPartialEnd].fGetField(m_poPartial->m_nFI_fCalcRotMid);
              fCalcRotMid +=  fShift;
              fObsRotStart -= fObsRotWidth/2.0;
              fObsRotEnd += fObsRotWidth/2.0;
                
              // Make enough room, and clear the contributions to each partial.
              afFractionContrib.setsize(nPartialBins);
              pfFractionContrib = &afFractionContrib[0];
              memset(pfFractionContrib,0,sizeof(double)*(nPartialBins));
                
                
                
              // Calculate the expansion factor as a product of lorentz factor
              // and the theta dependent expansion factor.
              
              fExpansionFactor = m_pfLorentz[nIntRef]*m_pfDStar[nIntRef];
              // Build fractional contributors.
              fFractionCounts = 0.0;
              fFractionCountsNotAccounted = 0.0;
              for (nPoint = 0; nPoint < m_nProfilePoints; nPoint++)
                {
                  fProfileRot = fCalcRotMid + (m_pfStandardProfileX[nPoint])*fExpansionFactor/m_fBroadening;
                  nBin = (int) ((fProfileRot - fObsRotStart)/fObsRotWidth);
                  if ((fProfileRot >= fObsRotStart) && (nBin < nPartialBins))
                    {
                      pfFractionContrib[nBin] += m_pfStandardProfileY[nProfileIndex*m_nProfilePoints + nPoint];
                      fFractionCounts += m_pfStandardProfileY[nProfileIndex*m_nProfilePoints + nPoint];

                    } 
                  else
                    fFractionCountsNotAccounted += m_pfStandardProfileY[nProfileIndex*m_nProfilePoints + nPoint];
                }
                
              // Normalize pfFractionContrib[];
              for (nBin = 0; nBin < nPartialBins; nBin++)
                pfFractionContrib[nBin]/=max(0.00001,fFractionCounts + fFractionCountsNotAccounted);                
              f0 = fFractionCounts + fFractionCountsNotAccounted;
              if (f0)
                {
                  fFractionCountsNotAccounted /= f0;
                  fFractionCounts /= f0;
                }
                
              // Fit the normalized profile to the data with sigma weighting.
              double fSumNumer;
              double fSumDenom;
              double fIntensity,fIntensityProfit;
              double fSigma,fSigmaProfit;
              double fK;
              double fChiSq=0.0;
              double fFitIntensity;
              double fSumIntensityToLeft;
              int    nMaxWidth;
                
              // We are minimizing sum((( fIntensity - K * fProfileFraction)/fSigma)^2)
                
              fK = 0.0;
              fTotalIntensity = 0.0;       
              fTotalFitIntensity = 0.0;
              fSumIntensityToLeft = 0.0;
              nMaxWidth = 0;
              for (nPass = 0; nPass < 2; nPass++)
                {
                  fSumNumer = 0.0;
                  fSumDenom = 0.0;
                  fChiSq = 0.0;
                  fSumIntensityToLeft = 0.0;
                    
                  for (nBin = 0; nBin < nPartialBins; nBin++)
                    {
                        
                      fIntensityProfit = (*m_poPartial)[nPartialStart + nBin].fGetField(m_poPartial->m_nFI_fProfitIntensity);
                      fSigmaProfit = (*m_poPartial)[nPartialStart + nBin].fGetField(m_poPartial->m_nFI_fProfitSigmaI);
                      fIntensity = (*m_poPartial)[nPartialStart + nBin].fGetField(m_poPartial->m_nFI_fIntensity);
                      fSigma = (*m_poPartial)[nPartialStart + nBin].fGetField(m_poPartial->m_nFI_fSigmaI);
                      nPartialStat = (*m_poPartial)[nPartialStart + nBin].nGetField(m_poPartial->m_nFI_nPartialFlag);                         
                        
                      if (!DONTUSEPROFIT)
                        {
                          fIntensity = fIntensityProfit;
                          fSigma = fSigmaProfit;
                        }
                      if (nPass == 0)
                        {
                          // Calculate parameters for doing fitting.
                          fSumDenom += pfFractionContrib[nBin]*pfFractionContrib[nBin]/fSigma/fSigma;
                          fSumNumer += pfFractionContrib[nBin]*fIntensity/fSigma/fSigma;                   
                          fTotalIntensity += fIntensity;
                        } 
                      else
                        {
                          // Calculate the (normalized) chi^2 of the fit.
                          fFitIntensity = fK*pfFractionContrib[nBin];
                          fTotalFitIntensity += fFitIntensity;
                          f0 = (fFitIntensity - fIntensity)/fSigma;
                          fChiSq += f0*f0;
                          afCalcValues + fFitIntensity;
                          afObsValues + fIntensity;
                          afSigmaValues + fSigma;

                          if (    (fTotalIntensity > 0.0) 
                               && (fIntensity > 0.0) 
                               && (fSumIntensityToLeft + fIntensity > fTotalIntensity/2.0) 
                               && (fProfileShift== 0.0))
                            {
                              fProfileShift = fObsRotStart +  nBin*fObsRotWidth + fObsRotWidth*(fTotalIntensity/2.0 - fSumIntensityToLeft)/fIntensity - fCalcRotMid;
                            }
                          fSumIntensityToLeft += fIntensity;


                          // See if this slice is noise or not.
                          if (    (min(fFitIntensity,fIntensity) < fSigma)
                               && (!(nPartialStat & C3Ddata::ms_nPerSliceCentroid))) 
                            (*m_poPartial)[nPartialStart + nBin].vSetField(m_poPartial->m_nFI_nPartialFlag,nPartialStat | C3Ddata::ms_nPerSlicePostRefineNotInStatisticalMosaicity);
                          else if (nPartialStat & C3Ddata::ms_nPerSlicePostRefineNotInStatisticalMosaicity)
                            (*m_poPartial)[nPartialStart + nBin].vSetField(m_poPartial->m_nFI_nPartialFlag,nPartialStat & ~C3Ddata::ms_nPerSlicePostRefineNotInStatisticalMosaicity);

                        }
                        
                      // Add to the next standard profile.
                        if (   (nPass == 1) 
                            && (fK > 0.0) 
                            && (fTotalIntensity > 0.0))
                          {
                            double fNormXStart;
                            double fNormXEnd;
                            int nStart;
                            int nEnd;
                            int nx;

                            fNormXStart = (nBin*fObsRotWidth - fCalcRotMid + fObsRotStart)/fExpansionFactor;
                            fNormXEnd = ((nBin + 1)*fObsRotWidth - fCalcRotMid + fObsRotStart)/fExpansionFactor;
                            nStart = (int) ((fNormXStart - m_pfStandardProfileX[0]) / (m_pfStandardProfileX[1] - m_pfStandardProfileX[0]));
                            nEnd = (int) ((fNormXEnd - m_pfStandardProfileX[0]) / (m_pfStandardProfileX[1] - m_pfStandardProfileX[0]));
                            nStart = max(0,min(m_nProfilePoints - 1,nStart));
                            nEnd = max(0,min(m_nProfilePoints - 1,nEnd));
                            
                            // Get the shape of the last profile at this point.  Scale it accordingly.
                            f0 = 0.0;
                            for (nx = nStart; nx <= nEnd; nx++) 
                              f0 += m_pfStandardProfileY[nProfileIndex*m_nProfilePoints + nx];
                            if (1 /*f0 == 0.0 */)
                              {
                                for (nx = nStart; nx <= nEnd; nx++)
                                  {
                                    m_pfStandardProfileYNext[nProfileIndex*m_nProfilePoints + nx] += max(0.0,fIntensity/fTotalIntensity); //max(0.0,fFitIntensity/fK);
                                  }
                              } 
                            else
                              {
                                for (nx = nStart; nx <= nEnd; nx++)
                                  {
                                    m_pfStandardProfileYNext[nProfileIndex*m_nProfilePoints + nx] += max(0.0,fIntensity/fTotalIntensity)*m_pfStandardProfileY[nProfileIndex*m_nProfilePoints + nx]/f0; 
                                  }
                              }
                            
                            nMaxWidth = max(nMaxWidth,nEnd - nStart);
                          }
                    }
                  if (nPass==0)
                    {
                      if (fSumDenom == 0.0)
                        fK = 0.0;
                      else
                        fK =fSumNumer/fSumDenom;
                    } 
                  else
                    fChiSq /= (nPartialBins - 1);
                }

                // Record this width in the profile.
                if (nMaxWidth>0)
                  {
                    for (nx = -nMaxWidth/2;nx <= nMaxWidth/2; nx++)
                      m_pfStandardProfileWidthY[nResoGroup*m_nProfilePoints + m_nProfilePoints/2 + nx] += 1.0;                            
                  }

                a2fChiSq[nBatchesTested] = fChiSq;
                a2nProfileBatch[nBatchesTested] = nProfileBatch;
                nBatchesTested++;
            }
        }
      if (nBatchesTested)
        {
          nx = ((nBatchesTested == 1) || (a2fChiSq[0] < a2fChiSq[1]))?0:1;
          m_pnBestProfileBatchi[nIntRef] = nx;
                           
          m_pfRefineBatchChiSq[nRefineBatch] += a2fChiSq[nx];
          m_pnRefineBatchContrib[nRefineBatch] ++;

            
          if ((fProfileShift) && (fTotalIntensity > 0.0))
            {
              m_pfRefineBatchObsShift[nRefineBatch] += fProfileShift*fTotalIntensity;
              m_pfRefineBatchContrib[nRefineBatch] += fTotalIntensity;
            }
          fSumChiSq += a2fChiSq[nx];
          nSumChiSq ++;                           
        }
    }

  fSumChiSq/= nSumChiSq;
  m_fBestChiSq = fSumChiSq;

  for (nRefineBatch = 0; nRefineBatch < m_anRefineBatchStart.size(); nRefineBatch++)
    {
      if (m_pnRefineBatchContrib[nRefineBatch] > 10)
        m_pfRefineBatchChiSq[nRefineBatch] /= m_pnRefineBatchContrib[nRefineBatch];
      else
        m_pfRefineBatchChiSq[nRefineBatch] = 0.0;
    }

  return fSumChiSq;   
}


int Cprofit2DSlices::nBuildStandardProfile(bool bInit)
{
  double f0;
  int nx,ny;
  volatile bool bDebug = true;

  // Build the standard profile X (axis) array.
  if (bInit)
    {
      -m_afStandardProfileX;
      for (f0 = -m_fStandardProfileMax; f0 <= m_fStandardProfileMax; f0 += m_fStandardProfileStep)
        m_afStandardProfileX + f0;
      m_nProfilePoints = m_afStandardProfileX.size();
    }
  m_pfStandardProfileX = &m_afStandardProfileX[0];
    
  if (bInit)
    {
      double fM;
      double fCoeff;
      double fX,fX0;
      double fFWHM;
      fM = 1.0;
      fCoeff = 4.0*(exp(log(2.0)/fM)-1.0);
      fX0 = 0.0;
      fFWHM = 1.0;
      m_afStandardProfileY.setsize(m_nNumProfiles*m_nProfilePoints);
      m_pfStandardProfileY = &m_afStandardProfileY[0];
      for (ny = 0; ny < m_nNumProfiles; ny++)
        {
          for (nx=0;nx < m_nProfilePoints;nx++)
            {
              fX = m_pfStandardProfileX[nx];
              f0 = 1.0 + fCoeff*(fX - fX0)*(fX - fX0)/fFWHM;
              f0 = exp(-fM*log(f0));
              m_pfStandardProfileY[nx + ny*m_nProfilePoints] = f0;
            }
        }
    } 
  else 
    {
      for (ny = 0; ny < m_nNumProfiles; ny++)
        {
          for (nx=0;nx < m_nProfilePoints;nx++)
            {
              m_pfStandardProfileY[ny*m_nProfilePoints + nx] = m_pfStandardProfileYNext[ny*m_nProfilePoints + nx]; // / m_pfStandardProfileYNextCounts[ny*m_nProfilePoints + nx];
            }
        }
    }
  return 0;
}

int Cprofit2DSlices::nSuggest()
{
  int nx;

  -m_afRefineParamShift;
  for (nx=0; nx< m_anRefineBatchStart.size(); nx++)
    {
      m_afRefineParamShift + 0.0;
    }
  m_pfRefineParamShift = &m_afRefineParamShift[0];
  return 0;
}

int Cprofit2DSlices::nPrintMosaicityModel(int nPart)
{
  int nx,ny;
  bool bMakePlots = FALSE;
  bool bPrintMosaicity = TRUE;

  if (bPrintMosaicity)
    {
      const int nShells = 6;
      int nShell;
      double afShellFraction[nShells] = {0.7,0.5,0.3,0.1,0.05,0.01};
      double afAreaMosaicity[nShells];
      double afMaxMosaicity[nShells];
      Cstat oStatArea;
      Cstat oStatMax;
      int nMax,nBound,nStart,nEnd;
      double fIntegrated,fTotal;
      double fMosaicity;
        
      for (nShell = 0; nShell < nShells; nShell++)
        {
          oStatArea.vClearRaw();
          oStatMax.vClearRaw();
          afAreaMosaicity[nShell] = 0.0;
          afMaxMosaicity[nShell]  = 0.0;
          for (ny = 0; ny < m_nNumProfiles; ny++)
            {
              nMax = 0;
              fTotal = 0.0;
              for (nx=0;nx < m_nProfilePoints;nx++)
                {
                  if (   m_pfStandardProfileY[nx + ny*m_nProfilePoints]
                      >  m_pfStandardProfileY[nMax + ny*m_nProfilePoints])
                    {
                      nMax = nx;
                    }
                  fTotal += m_pfStandardProfileY[nx + ny*m_nProfilePoints];
                }
              fIntegrated = m_pfStandardProfileY[nMax + ny*m_nProfilePoints];
              for (nBound = 1; 
                   ((fTotal - fIntegrated)/fTotal > afShellFraction[nShell]);
                   nBound++)
                {
                  if (nMax + nBound < m_nProfilePoints)
                    fIntegrated += m_pfStandardProfileY[nMax + nBound + ny*m_nProfilePoints];
                  if (nMax - nBound >= 0)
                    fIntegrated += m_pfStandardProfileY[nMax - nBound + ny*m_nProfilePoints];
                }
              nBound--;
              
              for (nStart = nMax; nStart >0; nStart--)
                {
                  if (m_pfStandardProfileY[nStart + ny*m_nProfilePoints] < afShellFraction[nShell]*m_pfStandardProfileY[nMax + ny*m_nProfilePoints]) 
                    break;
                }
              for (nEnd = nMax; nEnd < m_nProfilePoints - 1; nEnd++)
                {
                  if (m_pfStandardProfileY[nEnd + ny*m_nProfilePoints] < afShellFraction[nShell]*m_pfStandardProfileY[nMax + ny*m_nProfilePoints]) 
                    break;
                }
                
              fMosaicity = (m_pfStandardProfileX[1] - m_pfStandardProfileX[0])*(nBound*2 + 1);
              oStatArea.vAddRaw(fMosaicity);
              fMosaicity = (m_pfStandardProfileX[1] - m_pfStandardProfileX[0])*(nEnd - nStart + 1);
              oStatMax.vAddRaw(fMosaicity);
            }
            
          afAreaMosaicity[nShell] = oStatArea.fAverageRaw();
          afMaxMosaicity[nShell]  = oStatMax.fAverageRaw();
        }
        
      printf("\nRocking curve width analysis");
      if (nPart!=-1)
        {
          printf(" (for reflections %d through %d only)\n",m_anIntegrateFilePartStart[nPart],m_anIntegrateFilePartEnd[nPart]);
        } 
      else
        printf("\n");
      printf("------------------------------------\n");
      printf(" Full-width         (Degrees)\n");
      printf(" at %% of     Total_area   Max_height\n");
      printf("------------------------------------\n");
      for (nShell = 0; nShell < nShells; nShell++) 
        {
          printf("  %4.0lf%%         %.3lf        %.3lf\n",
                 100.0*afShellFraction[nShell], afAreaMosaicity[nShell],
                 afMaxMosaicity[nShell]);
        }
      printf("------------------------------------\n\n");
    }
    
  if (bMakePlots)
    {
      FILE* pFOut;
      Cstring sTemp;
      Cstring sFileOut;
      int nReso,nBatch;
      int nIndex;
      double f0;
      sTemp = sTransSymbol(sDtrekGetPrefix());
        
      sFileOut = sTemp;
      sFileOut += "profile_conv.txt";
      pFOut = fopen(sFileOut.string(),"wt");
      // Normalize the width profile.
      for (nx=0;nx<m_nResoGroups;nx++)
        {
          f0 = 0.0;
          for (ny = 0; ny < m_nProfilePoints;ny++)
            f0 += m_pfStandardProfileWidthY[nx*m_nProfilePoints + ny];
          for (ny = 0; ny < m_nProfilePoints;ny++)
            fprintf(pFOut,"%lf %lf\n",m_pfStandardProfileX[ny],m_pfStandardProfileWidthY[nx*m_nProfilePoints + ny]/f0);
          fprintf(pFOut,"\n\n");
        }
      fclose(pFOut);
        
      for (nReso = 0; nReso < m_nResoGroups;nReso++)
        {
          if (m_anProfileBatchIndex[nReso].size())
            {
              sFileOut = sTemp;
              sFileOut += "profile_reso";
              sFileOut += nReso;
              sFileOut += ".txt";
              printf("INFO: Writing profile file '%s'\n",sFileOut.string());
              pFOut = fopen(sFileOut.string(),"wt");
              for (nBatch = 0; nBatch < m_anProfileBatchIndex[nReso].size(); nBatch++)
                {
                  nIndex = m_anProfileBatchIndex[nReso][nBatch];
                    
                  f0 = 0.0;
                  for (nx=0;nx< m_nProfilePoints;nx++) 
                    f0 += m_pfStandardProfileY[nIndex*m_nProfilePoints + nx];
                  f0 = max(1.0,f0);
                    
                  for (nx=0;nx< m_nProfilePoints;nx++)
                    {
                      fprintf(pFOut,"%lf %lf\n",m_pfStandardProfileX[nx],m_pfStandardProfileY[nIndex*m_nProfilePoints + nx]/f0);
                    }
                  fprintf(pFOut,"\n");
                }
              fclose(pFOut);
            }
        }
    }
    
  return 0;
}

int Cprofit2DSlices::nRefineMosaicityModel()
{
  itr<double> afValueY0;
  itr<double> afValueYObs;
  itr<double> afSigmaY0;
  int nx;
  double fChiSq0;
  double fShift=0.0;
  int nRefineBatch;
  
  if (nBuildStandardProfile(true))
    return 1;
  if (nSuggest())
    return 1;
    
  printf("\nWorking ... ");
  fflush (stdout);
  for (nx = 0; nx< 4; nx++)
    {
      fChiSq0 = fEvaluate(afValueY0,afValueYObs,afSigmaY0);
      nBuildStandardProfile(false);
      fChiSq0 = fEvaluate(afValueY0,afValueYObs,afSigmaY0);
      for (nRefineBatch =0 ;nRefineBatch < m_afRefineBatchChiSq.size(); nRefineBatch++)
        {
          if (m_pfRefineBatchContrib[nRefineBatch]>0.0) 
            fShift = m_pfRefineBatchObsShift[nRefineBatch] /m_pfRefineBatchContrib[nRefineBatch];
          m_pfRefineParamShift[nRefineBatch] += fShift;
        }
    }
  printf(" Almost done.\n");
  fflush (stdout);
  return 0;
}


int Cprofit2DSlices::nFindIntersections()
{
  double      a2x2fEllipse1A[2][2];
  double      a2fEllipse1b[2];
  double      fEllipse1c;
  double      a2x2fEllipse2A[2][2];
  double      a2fEllipse2b[2];
  double      fEllipse2c;
  double      fIntensity;
  double      fFraction;
  double      fMaxIntersection,fIntersection;
  int         a3nHKL[3];          // HKL of the reflection.
  int         a3nRelHKL[3];       // HKL of neighboring reflections.
  int         nPackedHKL;
  int         nRef1,nRef2,nRef2i;
  int         nSlice1,nSlice2;
  int         nPartialStat1,nPartialStat2;
  int*        pnSort;
  Crefln*     poRefln;
  int         nx;
  double      f0; 
  double      fObsRotMid1,fObsRotMid2;
  double      fFind1Pix0,fFind1Pix1;
  double      fFind2Pix0,fFind2Pix1;
  bool        bCent1In2,bCent2In1;
  const       int g_nMaxHKLSearch = 1;

  int         nTotalIntersects;
  int         a20nAvgIntersect[20];

  m_poIntegrate->nExpandGetField(m_poPartial->ms_snPackedHKL);
  for (nRef1 = 0; nRef1 < m_poIntegrate->nGetNumReflns();nRef1++) 
    (*m_poIntegrate)[nRef1].vSetField(m_poIntegrate->m_nFI_nPackedHKL,(*m_poIntegrate)[nRef1].nPackHKL());
  m_poIntegrate->vSort(eReflnField_int_type,m_poIntegrate->m_nFI_nPackedHKL,NULL);
  pnSort = m_poIntegrate->pnGetSortIndex();    

  // Clear all intersections.
  nTotalIntersects = 0;
  for (nx = 0; nx < 20; nx++)
    {
      a20nAvgIntersect[nx] = 0;
    }

  for (nRef1 = 0; nRef1 < m_poIntegrate->nGetNumReflns(); nRef1++)
    {
      fIntensity = (*m_poIntegrate)[nRef1].fGetIntensity();
      if ((m_anPartialStart[nRef1] >= 0) && (fIntensity>0.0))
        {
          fFind1Pix0 = (*m_poIntegrate)[nRef1].fGetField(m_poIntegrate->m_nFI_fObsPx0);
          fFind1Pix1 = (*m_poIntegrate)[nRef1].fGetField(m_poIntegrate->m_nFI_fObsPx1);
          fMaxIntersection = 0.0;
          for (nSlice1 = m_anPartialStart[nRef1]; nSlice1 <= m_anPartialEnd[nRef1]; nSlice1++)
            {
              (*m_poPartial)[nSlice1].nPutGetEllipse(false,&a2x2fEllipse1A[0][0],&a2fEllipse1b[0],&fEllipse1c);
              fObsRotMid1 = (*m_poPartial)[nSlice1].fGetField(m_poPartial->m_nFI_fObsRotMid);
              nPartialStat1 = (*m_poPartial)[nSlice1].nGetField(m_poPartial->m_nFI_nPartialFlag);
              poRefln = m_poIntegrate->poGetRefln(nRef1);
              poRefln->nUnPackHKL(poRefln->nGetField(m_poIntegrate->m_nFI_nPackedHKL),a3nHKL);
              for (a3nRelHKL[0] = a3nHKL[0] - g_nMaxHKLSearch; 
                   a3nRelHKL[0] <= a3nHKL[0] + g_nMaxHKLSearch; a3nRelHKL[0] ++)
                {
                  for (a3nRelHKL[1] = a3nHKL[1] - g_nMaxHKLSearch; 
                       a3nRelHKL[1] <= a3nHKL[1] + g_nMaxHKLSearch; a3nRelHKL[1] ++)
                    {
                      for (a3nRelHKL[2] = a3nHKL[2] - g_nMaxHKLSearch; 
                           a3nRelHKL[2] <= a3nHKL[2] + g_nMaxHKLSearch; 
                           a3nRelHKL[2] ++) 
                        {
                          if (   (a3nRelHKL[0] != a3nHKL[0])
                              || (a3nRelHKL[1] != a3nHKL[1])
                              || (a3nRelHKL[2] != a3nHKL[2]))
                            {
                              nPackedHKL = poRefln->nPackHKL(a3nRelHKL);
                              nRef2i = m_poIntegrate->nFindFirst(m_poIntegrate->m_nFI_nPackedHKL,nPackedHKL);
                              if (nRef2i>=0)
                                {
                                  while ((nRef2i< m_poIntegrate->nGetNumReflns()) && ((*m_poIntegrate)[pnSort[nRef2i]].nGetField(m_poIntegrate->m_nFI_nPackedHKL) == nPackedHKL))
                                    {
                                      nRef2 = pnSort[nRef2i]; 
                                      fIntersection = 0.0;                                    
                                      if (m_anPartialStart[nRef2]>=0)
                                        {
                                          fFind2Pix0 = (*m_poIntegrate)[nRef2].fGetField(m_poIntegrate->m_nFI_fObsPx0);
                                          fFind2Pix1 = (*m_poIntegrate)[nRef2].fGetField(m_poIntegrate->m_nFI_fObsPx1);
                                            
                                          for (nSlice2 = m_anPartialStart[nRef2]; nSlice2 <= m_anPartialEnd[nRef2]; nSlice2++)
                                            {
                                              (*m_poPartial)[nSlice2].nPutGetEllipse(false,&a2x2fEllipse2A[0][0],&a2fEllipse2b[0],&fEllipse2c);
                                              fObsRotMid2 = (*m_poPartial)[nSlice2].fGetField(m_poPartial->m_nFI_fObsRotMid);
                                              nPartialStat2 = (*m_poPartial)[nSlice2].nGetField(m_poPartial->m_nFI_nPartialFlag);

                                              // Check to see if we are on the same image.
                                              if (ABS(fObsRotMid1 - fObsRotMid2) < 0.001)
                                                {
                                                  // Check the intersection of ellipsoids.
                                                  bCent1In2 = ((f0 = fEvalEllipse(fFind1Pix0 - fFind2Pix0,fFind1Pix1 - fFind2Pix1,&a2x2fEllipse2A[0][0],&a2fEllipse2b[0],&fEllipse2c))<=1.0);
                                                  bCent2In1 = ((f0 = fEvalEllipse(fFind2Pix0 - fFind1Pix0,fFind2Pix1 - fFind1Pix1,&a2x2fEllipse1A[0][0],&a2fEllipse1b[0],&fEllipse1c))<=1.0);
                                                  if ((bCent1In2) && (!(nPartialStat1 & C3Ddata::ms_nPerSlicePostRefineNotInStatisticalMosaicity)))
                                                    fIntersection = max(fIntersection,(*m_poPartial)[nSlice1].fGetIntensity());
                                                  if ((bCent2In1) && (!(nPartialStat2 & C3Ddata::ms_nPerSlicePostRefineNotInStatisticalMosaicity)))
                                                    fIntersection = max(fIntersection,(*m_poPartial)[nSlice1].fGetIntensity());
                                                }
                                            }
                                        }
                                      fMaxIntersection = max(fMaxIntersection,fIntersection);
                                      nRef2i++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
          fFraction = fMaxIntersection/fIntensity;
          if (fFraction > 0.0)
            {
              nx = (int) (log(fFraction*1000.0)/log(2.0));
              nx = min(11-1,max(0,nx));
              a20nAvgIntersect[nx]++;
            }

          if (fFraction > m_fIntersectMax)
            {
              // Found an intersection.
              (*m_poIntegrate)[nRef1].m_nSelect = 1;
              nTotalIntersects++;
            }
        }
    }
    
    
    if (nTotalIntersects)
      {
        printf("\nTotal spatial overlaps (intersections): %d\n",nTotalIntersects);    
        printf("--------------------------------------------------------------------\n");
        printf(" <=0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.256 0.512 >=1.0\n");
        printf("   ");
        for (nx=0;nx<11;nx++)
            printf("%5d ",a20nAvgIntersect[nx]);
        printf("\n");
        printf("--------------------------------------------------------------------\n\n");       
    }
    m_nTotalIntersections = nTotalIntersects;

    delete[] m_poIntegrate->m_pnSortIndex;
    m_poIntegrate->m_pnSortIndex = NULL;


    return 0;
};

int Cprofit2DSlices::nSplitFileProcessSize(Cstring& sIntegrateHeader,int nMaxSize) {
    Cimage_header oHeader(sIntegrateHeader);
    if (!oHeader.bIsAvailable())
        return -1;
    else
        return nSplitFileProcessSize(oHeader,nMaxSize);
};



int Cprofit2DSlices::nSplitFileProcessSize(Cimage_header& oHeader,int nMaxSize) {
    itr<int> anCutPoints;
    int nNumCuts;
    int nPart,nLastPart;
    int nAvgIntCutSize;

    int nLastIntegrateStart,nLastPartialStart;
    int nx;
    int nTotalInts,nAvgPartialCutSize;
    
    if (oHeader.nGetValue((Cstring) "DTPROFIT_CUTPOINTS",1,&anCutPoints[0])) 
        return -1;
    if ((anCutPoints[0] % 2) || (anCutPoints[0]/2 < 2))
        return -1;
    anCutPoints.setsize(anCutPoints[0] + 1);
    if (oHeader.nGetValue((Cstring) "DTPROFIT_CUTPOINTS",(int) (1 + anCutPoints[0]),&anCutPoints[0])) 
        return -1;
    

    -m_anIntegrateFilePartStart;
    -m_anIntegrateFilePartEnd;
    -m_anPartialFilePartStart;
    -m_anPartialFilePartEnd;

    // Build the 'start' and 'stop' values for each of the segments that we will be proccessing.
    // Discover the average size of the cutpoints.
    nTotalInts = (anCutPoints[anCutPoints.size() - 2]);
    nAvgPartialCutSize = (anCutPoints.last() - anCutPoints[2])/(anCutPoints[0]/2 - 1);
    nAvgIntCutSize = (anCutPoints[3] - anCutPoints[1]);
    if (nTotalInts<=nMaxSize)
        return -1;
    
    nNumCuts = (int)ceil((double) nTotalInts / (double)nMaxSize);

    nLastIntegrateStart = 0;
    nLastPartialStart = 0;
    nLastPart = 0;
    
    for (nx = 1 + 2; nx < anCutPoints.size(); nx+=2) {        
        nPart = (int) ((((double) anCutPoints[nx])/nTotalInts)*nNumCuts);
        if ((nPart > nLastPart) || (nx+2 >= anCutPoints.size())) {
            m_anIntegrateFilePartStart + nLastIntegrateStart;
            m_anPartialFilePartStart + nLastPartialStart;
            m_anIntegrateFilePartEnd + ((anCutPoints[nx] - 1) + ((nx+2 >= anCutPoints.size())?nAvgIntCutSize:0));
            m_anPartialFilePartEnd + ((anCutPoints[nx+1] - 1) + ((nx+2 >= anCutPoints.size())?(nAvgPartialCutSize*3):0));
            nLastIntegrateStart = (anCutPoints[nx]);
            nLastPartialStart = (anCutPoints[nx+1]);
            nLastPart = nPart;
        }
    };


    return m_anIntegrateFilePartStart.size();
};
