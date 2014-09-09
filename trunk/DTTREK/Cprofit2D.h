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
#ifndef DT_CPROFIT2D_H
#define DT_CPROFIT2D_H

#include "Dtrek.h"
#include "Cimage.h"

#include "Cnonunf.h"
#include "dtarray.h"

#include "CinterpMesh.h"

class Cprofit2DSliceModel;
extern int g_nDebug;
const float g_fMaskOut = -99999;
const int g_nMaxReflnsToProfit = 100000;

class DTREK_EXPORT tProfitRec {
    double* m_pfX;                  // Sum of X
    double* m_pfVar;                // Variances associated with each pixel
    double  m_fContrib;             // Contribution count.
    double  m_fIntensity;           // When calling nLoad(float* ...) this is loaded.
	double  m_fNominalVar;			// Set by a call to nLoad().
    bool bHasX;                     // At least one constructor allows user to specify these pointers.
    bool bHasX2;

public:

    // Do NOT use these functions to modify the data.  Only to read it for printing purposes.
    double* pfGetData()             { return m_pfX; };
    double  fGetDataCounts()        { return m_fContrib; };

    static int ms_a3nDim[3];        // All spots get the same profile dimensions.
    static int ms_a2nCOG[2];        // All spots use a COG pixel to line up profiles.
	static itr<int> ms_anInUse;		// Loaded when nLoad() is called.

    void vClear();
    int nLoad(float* pfDataX,float* pfDataVar,float fNominalVar,int nDim[3],double* pfCOG);
    int nSave(Cimage& oImage,int a2nPointer[2],double fScale);
    int nAdd(tProfitRec&);
    double fGetIntensity() { return m_fIntensity; };
    int nScale(tProfitRec& oRaw,double& fIntensity,double& fVariance,double& fCorr);
	int nDiagnostic(Cstring& sName,tProfitRec& oRaw,tProfitRec& oAverage,double fScale);
	int nCalcCOG(double* pfCOG);
   
    friend int tProfitRec_nAddFunctionPerSlice(int nDataPoints,double* pfInData,double* pfOutData);
    friend class Cprofit2DSliceModel;


    tProfitRec();
    tProfitRec(int nDim[2]);
    tProfitRec(double* pfX,double* pfX2,double fContrib);
    ~tProfitRec();
};


class DTREK_EXPORT Cprofit2DPerSlice {
	itr<CinterpMesh*> m_apoProfiles;
	itr<int>		 m_anImageStart;
	itr<int>		 m_anImageEnd;
    int              m_nMaxImage;
    int              m_nMinImage;
	int				 m_nCurrentImage;
	int				 m_nProfitNumReflections;
	int				 m_nProfitMaxImageRange;
    int              m_nMaxMeshesToKeep;


	// Temporary buffers to avoid constant reallocation.
	tProfitRec*		 m_poTempProfitRec;
	itr<float>		 m_afFloatBuf;

    // Bitmap related information and functions.
    unsigned char*   m_pcBitmaps;
    int              m_nBitmapUsed;
    int              m_nBitmapAlloc;
    itr<int>         m_a3nExt[3];
    itr<int>         m_a3nOrigOffset[3];
    itr<int>         m_anBitmapOffset;
    int              nAllocBitmaps(int nDesiredSize);
    int              nGetBitmapPixel(int nOffset,int a3nPixel[3],int a3nExt[3]);


	int nLoadProfile(int nSlice,C3Ddata& oData,bool bAddToLast,bool bLoadProfitRec,double* pfCOG = NULL);
public:
	int m_nDim0,m_nDim1;
	int	m_a2nProfileDim[2];
	int m_nMaxProfilesLevel;
    double m_fMinIoverSigForContributor;

	void vSetCurrentImage(int nImage) { m_nCurrentImage = nImage; };
	void vSetInterpConditions(int nProfitNumReflections,int nProfitMaxImageRange) { m_nProfitNumReflections = nProfitNumReflections; m_nProfitMaxImageRange = 2*max(1,((1+nProfitMaxImageRange)/2)); };
	int  nWriteScanBitmap(Cstring& sName,Cimage_header* poHeader);
	int  nAddBitmap(C3Ddata& oData);
	
	int nDiagnostic(Cstring& sName,int nImage);
	int nAdd(int nSliceStart,int nSliceEnd,C3Ddata& oData);  
	int nFit(int nSlice,C3Ddata& oData,double a2fOverallCent[2],float& fIntensity,float& fSigma); 
	int nDetermineMesh(int nIndex,int nCurrentImageNumber = -1);

	Cprofit2DPerSlice();
	~Cprofit2DPerSlice();
};

const int g_nMaxResoGroups = 8;

class DTREK_EXPORT Cprofit2DSlices {
	Creflnlist* m_poIntegrate;
	Creflnlist* m_poPartial;
    Creflnlist* m_poScaleList;
    Creflnlist* m_poAxialRefs;
    Csource*    m_poSource;
    Ccrystal*   m_poCrystal;
	Cstring		m_sIntegrate;
	Cstring		m_sPartial;
    Cstring     m_sEllipsoid;
    int         m_nIntRefs;
	
    itr<double> m_afPerReflnScale;           // For each reflection, a scale factor.
    itr<int>    m_anPartialStart;            // For each reflection in m_poIntegrate, the Start/End of partials in m_poPartial
	itr<int>	m_anPartialEnd;

        
    itr<int>    m_anRefineBatchStart;       // Start/End of refinement batches used for profile fitting.
    itr<int>    m_anRefineBatchEnd;
    itr<double> m_afRefineParamShift;
    double*     m_pfRefineParamShift;
    itr<double> m_afRefineBatchChiSq;
    double*     m_pfRefineBatchChiSq;
    itr<int>    m_anRefineBatchContrib;
    int*        m_pnRefineBatchContrib;
    itr<double> m_afRefineBatchObsShift;
    double*     m_pfRefineBatchObsShift;
    itr<double> m_afRefineBatchContrib;
    double*     m_pfRefineBatchContrib;
    
    itr<int>    m_anProfileBatchContrib[g_nMaxResoGroups];
    itr<int>    m_anProfileBatchStart[g_nMaxResoGroups];      // Start/End of profile batches.
    itr<int>    m_anProfileBatchEnd[g_nMaxResoGroups];
    itr<int>    m_anProfileBatchIndex[g_nMaxResoGroups];

    itr<double> m_afLorentz;                // Convenient variables calculated here.
    double*     m_pfLorentz;
    itr<double> m_afTheta;                  // All of these are PER REFLECTION.
    double*     m_pfTheta;
    itr<double> m_afDStar;
    double*     m_pfDStar;
    itr<int>    m_a2anProfileBatch[2];
    int*        m_a2pnProfileBatch[2];
    itr<int>    m_anBestProfileBatchi;
    int*        m_pnBestProfileBatchi;
    itr<int>    m_anResoBatch;
    int*        m_pnResoBatch;

    itr<double> m_afShellThetaGroupMin;     // Resolution shell start/end
    itr<double> m_afShellThetaGroupMax;     // Resolution shell start/end
    itr<double> m_afShellResoGroupMin;      // Resolution shell start/end
    itr<double> m_afShellResoGroupMax;      // Resolution shell start/end


    itr<double> m_afStandardProfileX;       // Contains range of points out to 3.0 times the FWHM.
    double*     m_pfStandardProfileX;

    itr<double> m_afStandardProfileY;       // Standard profile derived from the peak profile.  
    double*     m_pfStandardProfileY;
    itr<double> m_afStandardProfileYNext;   // Standard profile for next iteration.
    double*     m_pfStandardProfileYNext;
    itr<double> m_afStandardProfileWidthY;  // Average width function used for deconvolution.
    double*     m_pfStandardProfileWidthY;


	int			nRemoveReflns(itr<int> anRefNumbers);
    
    double fEvaluate(itr<double>& afCalcValues,itr<double>& afObsValues,itr<double>& afSigmaValues);
    int nSuggest();
    int nLoadProfileBatches();
    int nBuildStandardProfile(bool bInit);

    itr<int> m_anIntegrateFilePartStart;// INPUT: Created by nSplitFileProcessSize()
    itr<int> m_anIntegrateFilePartEnd;  // INPUT: Created by nSplitFileProcessSize()
    itr<int> m_anPartialFilePartStart;  // INPUT: Created by nSplitFileProcessSize()
    itr<int> m_anPartialFilePartEnd;    // INPUT: Created by nSplitFileProcessSize()
    Cstring  sCreatePartName(Cstring& sName,int nPart);


public:
    double m_fStandardProfileStep;  // INPUT: Step size for standard profile.
    double m_fStandardProfileMax;   // INPUT: Maximum +/- range for standard profile.
    double m_fBroadening;           // INPUT: Broadening factor.

    double m_fGainEst;              // OUTPUT: Estimation of gain from background intensities and variances.
    double m_fMaxIntensityToProfit; // OUTPUT: Maximum intensity to profile fit.  
    double m_fBestChiSq;            // OUTPUT: Best chi^2 from fEvaluate().

    int    m_nMinRefsPerProfile;    // INPUT:  Minimum # of reflections per profile.
    int    m_nResoGroups;           // INPUT:  # of resolution bins for postprefinement.
    int    m_nNumProfiles;          // OUTPUT: Calculated # of profiles from m_nMinRefsPerProfile and m_nResoGroups.
    int    m_nProfilePoints;        // OUTPUT: Calculated # of points in standard profile from m_fStandardProfileStep and m_fStandardProfileMax
    double m_fIntersectMax;         // INPUT:  Maximum intersection fraction permitted.
    int    m_nTotalIntersections;   // OUTPUT: Total # of intersections.


	int nSyncLists(bool bSyncWithHKL);
	int nLoad(Cstring& sIntegrateList,Cstring& sProfitList,Cstring& sScaleFile,int nFilePart);
    int nSplitFileProcessSize(Cstring& sIntegrateHeader,int nMaxSize);
    int nSplitFileProcessSize(Cimage_header& oHeader,int nMaxSize);
    
    int nLoadHeader(Cstring& seHeader);
    int nLoadHeader(Cimage_header& oHeader);
	int nCalcRemoveNoise(bool bFirstPart,bool bLastPart);
	int nWrite(int nFilePart,bool bWriteIntegrate,bool bWritePartial,bool bWriteEllipsoid);
    int nPrintAxial(bool bPrint);
    void vSetIntegrateOutputName(Cstring& sName) { m_sIntegrate = sName; };
    int nRefineMosaicityModel();    
    int nPrintMosaicityModel(int nPart);
    int nFindIntersections();
    int nCalcWidthsHKL(double fFixMosaicity);
    int nWriteParts(bool bWriteIntegrate,bool bWritePartial,bool bWriteEllipsoid);


	Cprofit2DSlices();
	~Cprofit2DSlices();
};

const int g_nSliceModelContours = 30;
const int g_nMaxSpotWidth = 60;

#endif
