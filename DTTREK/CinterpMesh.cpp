//
// Copyright (c) 2001 Molecular Structure Corporation
//
// CinterMesh.cc        Initial author: T.J. Niemeyer  Jan-2001
//  This file contains the member functions of class CinterpMesh.cc
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

#include "CinterpMesh.h"
#include <memory.h>
#include "Cimage.h"

int nFindMinimum(int nCount,int nToFind,double* pfOrder,int* pnIndices) {
    int nHighest;
    int nx,ny;

    for (nx=0;nx<nToFind;nx++)
        pnIndices[nx] = -1;

    for (nx=0;nx<nCount;nx++) {
        nHighest = 0;
        for (ny=0;ny<nToFind;ny++) {
            if (pnIndices[ny]==-1) {
                nHighest = ny;
                break;
            } else if (pfOrder[pnIndices[nHighest]]<pfOrder[pnIndices[ny]]) {
                nHighest = ny;
            };
        };
        if ((pnIndices[nHighest]==-1) || (pfOrder[pnIndices[nHighest]]>pfOrder[nx])) {
            pnIndices[nHighest] = nx;
        };
    };
    return 0;
};
        


int CinterpMesh::ms_nLevelPointers[g_nMaxInterpLevel+2];
double CinterpMesh::ms_fLevelGridStep[g_nMaxInterpLevel+2];
int CinterpMesh::ms_nLevelGridStep[g_nMaxInterpLevel+2];
bool CinterpMesh::ms_bInitialized = FALSE;
  
CinterpMesh::CinterpMesh()
{
    int nx;
    m_pfItems = NULL;
    m_pnCounts = NULL;
	m_pnFirst = NULL;
    m_nMinimumToInterp = /*5*/ 1;
    m_nWidth = 1;
    m_nMaxLevel = 0;
	m_bSaveOriginals = false;
	m_fSigmaReject = 0.0;

    
    if (!ms_bInitialized) {
        // Fill in the static level pointers.
        ms_nLevelPointers[0] = 0;
        ms_fLevelGridStep[0] = 1.0;
        ms_nLevelGridStep[0] = 1;
        for (nx=1;nx<g_nMaxInterpLevel+2;nx++) {
            ms_nLevelGridStep[nx] = ms_nLevelGridStep[nx-1] + 1;
            ms_fLevelGridStep[nx] = 1.0/(ms_nLevelGridStep[nx]);
            ms_nLevelPointers[nx] = ms_nLevelPointers[nx-1] + ms_nLevelGridStep[nx-1]*ms_nLevelGridStep[nx-1];
        };
        ms_bInitialized = TRUE;
    };    
};

CinterpMesh::~CinterpMesh() {
    if (m_pfItems)
        delete[] m_pfItems;
    m_pfItems = NULL;
    if (m_pnCounts)
        delete[] m_pnCounts;
    m_pnCounts = NULL;
	if (m_pnFirst)
		delete[] m_pnFirst;
	m_pnFirst = NULL;
};

void CinterpMesh::vClear() {
    if (m_pnCounts)
        memset(m_pnCounts,0,sizeof(int)*ms_nLevelPointers[m_nMaxLevel+1]);
    if (m_pfItems)
        memset(m_pfItems,0,sizeof(double)*ms_nLevelPointers[m_nMaxLevel+1]*m_nWidth);
	if (m_pnFirst) {
		int nx,ny;
		ny = ms_nLevelPointers[m_nMaxLevel+1];
		for (nx=0;nx<ny;nx++)
			m_pnFirst[nx] = -1;
	};

    m_nNumContrib = 0;
	m_afOrigItems.clear();

};
  
void CinterpMesh::vInit(int nMaxLevel,int nWidth,double fDim0Min,double fDim0Max,double fDim1Min,double fDim1Max) {
    m_fDim0Range[0] = fDim0Min;
    m_fDim0Range[1] = fDim0Max;
    m_fDim0Range[2] = m_fDim0Range[1] - m_fDim0Range[0];
    m_fDim1Range[0] = fDim1Min;
    m_fDim1Range[1] = fDim1Max;
    m_fDim1Range[2] = m_fDim1Range[1] - m_fDim1Range[0];
    m_nMaxLevel = min(nMaxLevel,g_nMaxInterpLevel);
    m_nWidth = nWidth;

    if (m_pfItems)
        delete[] m_pfItems;
	m_pfItems = NULL;
    if (m_pnCounts)
        delete[] m_pnCounts;
	m_pnCounts = NULL;
	if (m_pnFirst)
		delete[] m_pnFirst;
	m_pnFirst = NULL;

	if (m_bSaveOriginals) {
		m_pnFirst = new int[ms_nLevelPointers[m_nMaxLevel+1]];
		m_afTempBuf.setsize(m_nWidth*4);
		m_pfTempBuf = &m_afTempBuf[0];
	} else {
		if (m_nWidth!=0) 
			m_pfItems = new double [m_nWidth*ms_nLevelPointers[m_nMaxLevel+1]];
		else
			m_pfItems = NULL;
	};
	
    m_pnCounts = new int [ms_nLevelPointers[m_nMaxLevel+1]];
    vClear();
};

void CinterpMesh::vSetInterpConditions(int nMinimumToInterp,double fSigmaReject) {
    m_nMinimumToInterp = nMinimumToInterp;
	m_fSigmaReject = fSigmaReject;
	m_bSaveOriginals = (m_fSigmaReject>0.0);
};



int CinterpMesh::nAdd(double fDim0,double fDim1,double* pfItems,nInterpCallback pxCallback) {
    int nLevelRange0,nLevelRange1;
    int nOffset;
    int nx,ny;
	int nLevel;
	double* pfInsert;

    // The loop below might actually add a point that is out of range:
    if ((fDim0>m_fDim0Range[1]) || (fDim0<m_fDim0Range[0]) ||
        (fDim1>m_fDim1Range[1]) || (fDim1<m_fDim1Range[0]))
        return 1;
    for (nLevel=0;nLevel<=m_nMaxLevel;nLevel++) {
        nLevelRange0 = (int) ((fDim0-m_fDim0Range[0])/(ms_fLevelGridStep[nLevel]*m_fDim0Range[2]));
        nLevelRange1 = (int) ((fDim1-m_fDim1Range[0])/(ms_fLevelGridStep[nLevel]*m_fDim1Range[2]));
        
        nOffset = (ms_nLevelPointers[nLevel] +  nLevelRange0 + nLevelRange1*ms_nLevelGridStep[nLevel]);
        
        if (pxCallback) {
            pxCallback(m_nWidth,pfItems,&m_pfItems[nOffset*m_nWidth]);
            m_pnCounts[nOffset]++;
        } else if (m_bSaveOriginals) {
			if (nLevel == 0) {
				m_afOrigItems.setsize(m_afOrigItems.size() + m_nWidth);
				pfInsert = &m_afOrigItems[m_afOrigItems.size() - m_nWidth];
				for (nx = 0; nx < m_nWidth; nx++)
					pfInsert[nx] = pfItems[nx];
			};
			if (m_pnFirst[nOffset]==-1) {
				m_pnFirst[nOffset] = m_anNextItem.size();
				m_anNextItem + (-1);
				m_anOrigItem + (m_afOrigItems.size()/m_nWidth - 1);
			} else {
				m_anNextItem + m_pnFirst[nOffset];
				m_pnFirst[nOffset] = (m_anNextItem.size() - 1);
				m_anOrigItem + (m_afOrigItems.size()/m_nWidth - 1);
			};

			m_pnCounts[nOffset]++;
		} else {       
            for (ny = 0; ny<m_nWidth; ny++) {
                if (m_pnCounts[nOffset]==0) 
                    m_pfItems[nOffset*m_nWidth + ny] = pfItems[ny];
                else 
                    m_pfItems[nOffset*m_nWidth + ny] += pfItems[ny];                                
            };
            m_pnCounts[nOffset]++;
        };                
    };
    m_nNumContrib++;


    return 0;
};

int CinterpMesh::nGet(double fDim0,
                      double fDim1,
                      double* pfItems,
                      double* pfItemsSigma,
                      int nLevelToUse) 
{
    int nLevel;
    int nLevelRange0;
    int nLevelRange1;
    int nOffset;
    int nx;
    int nCounts;

    if (pfItems)
        memset(pfItems,0,m_nWidth*sizeof(double));
   
    // The loop below might actually add a point that is out of range:
    if ((fDim0>m_fDim0Range[1]) || (fDim0<m_fDim0Range[0]) ||
        (fDim1>m_fDim1Range[1]) || (fDim1<m_fDim1Range[0]))
        return -1;

    -m_anInterpOffset;
    -m_anInterpCounts;

    for (nx=0;nx<=m_nMaxLevel;nx++) {

        nLevelRange0= (int) ((fDim0-m_fDim0Range[0])/(ms_fLevelGridStep[nx]*m_fDim0Range[2]));
        nLevelRange1= (int) ((fDim1-m_fDim1Range[0])/(ms_fLevelGridStep[nx]*m_fDim1Range[2]));       

        nOffset = ms_nLevelPointers[nx] +  nLevelRange0 + nLevelRange1*ms_nLevelGridStep[nx];
        nCounts = m_pnCounts[nOffset];

        if (nCounts<m_nMinimumToInterp) 
            break;
        m_anInterpOffset + nOffset;
        m_anInterpCounts + nCounts;

        if ((nLevelToUse != -1) && (nx + 1> nLevelToUse)) {
            nx++;
            break;
        };
        
    };

    // If there were not even enough on the top level, then we are done.
    if (nx==0) 
        return -1;

    nLevel = nx - 1;
    if ((nLevelToUse != -1) && (nLevel<nLevelToUse))
        return -1;

    nOffset = m_anInterpOffset[nLevel];
    nCounts = m_anInterpCounts[nLevel];

    if (pfItems) {
		pfCalcGetItem(nOffset,pfItems,pfItemsSigma);
    };
    return nLevel;
};

double* CinterpMesh::pfCalcGetItem(int nOffset,double* pfBuffer,double* pfVarBuffer) 
{
	if (m_bSaveOriginals) 
    {
		int nPass;
		int nCount;
		int nIndex,nOrigIndex;
		double* pfOrigItems;
		int* pnNextItem;
		int* pnOrigItem;
		double* pfC,*pfX,*pfX2,*pfVar;
		double f0;
		int nx;

		pfOrigItems = &m_afOrigItems[0];
		pnNextItem = &m_anNextItem[0];
		pnOrigItem = &m_anOrigItem[0];
		pfC = m_pfTempBuf;
		pfX = m_pfTempBuf + m_nWidth*1;
		pfX2 = m_pfTempBuf + m_nWidth*2;
		pfVar = m_pfTempBuf + m_nWidth*3;


		for (nPass = 0; nPass < 2; nPass++) {
			nCount = 0;
			memset(pfC,0,sizeof(double)*m_nWidth);
			memset(pfX,0,sizeof(double)*m_nWidth);
			memset(pfX2,0,sizeof(double)*m_nWidth);
			for (nIndex = m_pnFirst[nOffset];nIndex != -1;nIndex = pnNextItem[nIndex]) {
				nOrigIndex = pnOrigItem[nIndex];
				nCount++;
				for (nx = 0;nx < m_nWidth; nx++) {
					f0 = pfOrigItems[nOrigIndex*m_nWidth + nx];
					if ((nPass == 0) || ((f0 - pfBuffer[nx])*(f0 - pfBuffer[nx]) < m_fSigmaReject*m_fSigmaReject*pfVar[nx])) {
						pfC[nx] += 1.0;
						pfX[nx] += f0;
						pfX2[nx] += f0*f0;
					};
				};
			};
			if (nCount != m_pnCounts[nOffset]) {
				printf("ERROR:  In CinterpMesh.\n");
				return pfBuffer;
			};
			if ((nPass == 0) && (nCount>=5)) {
				for (nx = 0; nx < m_nWidth; nx++) {
                    if (nCount>1)
					    pfVar[nx] = (pfX2[nx]*nCount - pfX[nx]*pfX[nx])/nCount/(nCount - 1);
                    else
                        pfVar[nx] = 0.0;
                    if (nCount > 0)
					    pfBuffer[nx] = pfX[nx]/nCount;
                    else
                        pfBuffer[nx] = 0.0;
				};
			} else  {
                for (nx = 0; nx < m_nWidth; nx++) {
					pfBuffer[nx] = pfX[nx]/max(1.0,pfC[nx]);
                    if (pfVarBuffer) {
                        if (pfC[nx]>1.0)
                            pfVarBuffer[nx] = (pfX2[nx]*nCount - pfX[nx]*pfX[nx])/pfC[nx]/(pfC[nx] - 1);
                        else
                            pfVarBuffer[nx] = 0.0;
                    };
                };
				break;
			}
		};
	} else {
		int ny;
		for (ny=0;ny<m_nWidth;ny++) {
			pfBuffer[ny] = m_pfItems[nOffset*m_nWidth + ny]/ m_pnCounts[nOffset];
		};
	};
	return pfBuffer;
};

int
CinterpMesh::nTest(Cstring& sTestImage) {
    double nDim0[2];
    double nDim1[2];
    double fData;
    int nx,ny;
    int nStat;
    int nMaxLevel;

    nDim0[0] = 0;
    nDim0[1] = 1000;
    nDim1[0] = 0;
    nDim1[1] = 1000;

    nMaxLevel = 8;
    vInit(nMaxLevel,1,0,1000,0,1000);
    vSetInterpConditions(3);

    // This draws a parabola.
    if (0) {
        for (ny=0;ny<10;ny++) {
            for (nx=0;nx<nDim0[1];nx++) {
                fData = 1000*min(ny,10-ny);
                nStat = nAdd(nx,ny+(nx-nDim0[1]/2)*(nx-nDim0[1]/2)/(nDim0[1]*nDim0[1]/4)*(nDim1[1]*0.8) + nDim1[1]*0.1,&fData);
                if (nStat)
                    return 1;
            };
        };
    };
    // This draws a Serpenski gasket.
    if (1) {
        float a2fVertex[2];
        float a3x2fVertices[3][2];
        a3x2fVertices[0][0] = nDim0[0];
        a3x2fVertices[0][1] = nDim1[0];
        a3x2fVertices[1][0] = nDim0[1];
        a3x2fVertices[1][1] = nDim1[0];
        a3x2fVertices[2][0] = 0.5*(nDim0[0] + nDim0[1]);
        a3x2fVertices[2][1] = nDim1[1];
        
        a2fVertex[0] = a3x2fVertices[0][0];
        a2fVertex[1] = a3x2fVertices[0][1];
        
        for (nx=0;nx<10000;nx++) {
            fData = 1000;
            nAdd(a2fVertex[0],a2fVertex[1],&fData);
            ny = rand() % 3;
            a2fVertex[0] = (a2fVertex[0] + a3x2fVertices[ny][0])*0.5;
            a2fVertex[1] = (a2fVertex[1] + a3x2fVertices[ny][1])*0.5;
        };
    };
    
        

    // This places the figure in an image.
    Cimage oImage((int) nDim0[1],(int) nDim1[1],eImage_uI2);
    for (nx=0;nx<nDim0[1];nx++) {
        for (ny=0;ny<nDim1[1];ny++) {
            nStat = nGet(nx,ny,&fData,NULL);
            (oImage.*oImage.prnSetPixel)(nx,ny,max(0,0));
        };
    };
    oImage.nWrite(sTestImage);
    return 0;
};

int CinterpMesh::nPrint(Cstring& sTestImage,int nDim0,int nDim1) {
    double fData;
    double f0,f1;
    int nx,ny;
    int nStat;

    if (m_nWidth!=1) {
        printf("ERROR: CinterpMesh::nPrint() called with non unity field width.\n");
        return 0;
    };
    Cimage oImage((int) nDim0,(int) nDim1,eImage_uI2);
    for (nx=0;nx<nDim0;nx++) {
        for (ny=0;ny<nDim1;ny++) {
            f0 = ((float) nx)/nDim0;
            f1 = ((float) ny)/nDim0;
            nStat = nGet(f0*m_fDim0Range[2] + m_fDim0Range[0],f1*m_fDim1Range[2] + m_fDim1Range[0],&fData,NULL);
            (oImage.*oImage.prnSetPixel)(nx,ny,max(0.0,fData));
        };
    };
    oImage.nWrite(sTestImage);
    return 0;

};


double CinterpMesh::fGetAverageLevel() {
    int nx,ny,nz;
    double fSumNumer,fSumDenom;
    double f0;

    if ((!m_pnCounts) || (!m_pfItems))
        return 0.0;
    fSumNumer = 0.0;
    fSumDenom = 0.0;
    for (nx = 0;nx < m_nMaxLevel; nx++) {
        nz = ms_nLevelGridStep[nx]*ms_nLevelGridStep[nx];
        f0 = 0.0;
        for (ny = 0; ny < nz; ny ++ )
            f0 += (m_pnCounts[ms_nLevelPointers[nx] + ny]>=m_nMinimumToInterp);
        f0 /= nz;
        fSumNumer += nx*f0;
        fSumDenom += f0;
    };
    if (!fSumDenom)
        return 0.0;
    f0 = fSumNumer/fSumDenom;
    return f0;
};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/////////////////// CSphericalMesh /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



CSphericalMesh::CSphericalMesh() {
    vInit(3,1,90.0,5.0,5.0);
};

void CSphericalMesh::vClear() {
    m_nNumContrib = 0;
    m_nNumContribAtLastRecompute = 0;
    m_afOrigItems.clear();
};

void CSphericalMesh::vInit(int nHarmonicOrder,int nWidth,double fMaxPolar,double fPolarStep,double fAziStep) {

    m_nWidth = nWidth;
    m_nDim0Azi = max(2,(int) floor(0.5 + 360.0/fAziStep));
    m_fStep0Azi = fAziStep;
    m_nDim1Polar = max(2,(int) floor(0.5 + fMaxPolar/fPolarStep));
    m_fStep1Polar = fPolarStep;

    m_oHarmonics.vSetHarmonics(nHarmonicOrder,0,2);
    m_oHarmonics.vSetCentroSymmetric(1);
    vClear();
};

int CSphericalMesh::nAdd(double fPolar,double fAzi,double fSigma,double* pfItems) {
    int nx;
    m_afAzi + fAzi;
    m_afPolar + fPolar;
    m_afOrigItemsSigma + fSigma;
    for (nx = 0; nx < m_nWidth;nx++)
        m_afOrigItems + pfItems[nx];
    m_nNumContrib++;
    return 0;
};

int CSphericalMesh::nRecompute()
{
    int nVar;
    double* pfValues;
    double* pfValuesIn;
    double* pfValuesOut;
    int nPoint0,nPoint1;
    double fAzi,fPolar;
    int nPoint;
    int nPoints;
    itr<double> afAzi;
    itr<double> afPolar;
    itr<double> afRotAngle;
    
    nPoints = m_nNumContrib;
    if (!nPoints)
        return 0;
    m_afValues.setsize(nPoints);

    m_afMeshValues.setsize(m_nDim0Azi*m_nDim1Polar*m_nWidth);

    for (nPoint1 = 0; nPoint1 < m_nDim1Polar; nPoint1++) {
        fPolar = (nPoint1 + 0.5)*m_fStep1Polar;
        for (nPoint0 = 0; nPoint0 < m_nDim0Azi; nPoint0++) {
            fAzi = (nPoint0 + 0.5)*m_fStep0Azi;
            afAzi + fAzi;
            afPolar + fPolar;
        };
    };

    for (nVar = 0; nVar < m_nWidth; nVar++) {
        pfValues = & m_afValues[0];
        pfValuesIn = & m_afOrigItems[0] + nVar;
        for (nPoint = 0; nPoint < nPoints; nPoint++) {
            *pfValues = *pfValuesIn;
            pfValues++;
            pfValuesIn += m_nWidth;
        };
        m_afValues.setsize(nPoints);
        m_oHarmonics.nScale(m_afValues,m_afOrigItemsSigma,m_afPolar,m_afAzi);

        m_oHarmonics.nSimulate(m_afValues,afPolar,afAzi,afRotAngle);

        pfValues = & m_afValues[0];
        pfValuesOut = & m_afMeshValues[0] + nVar;

        for (nPoint = 0; nPoint < m_afValues.size(); nPoint++)
        {
            *pfValuesOut = *pfValues;
            
            //RB 3/27/06 I am commenting these lines out because they don't do anything
            //int nx;
            //if (*pfValues < 0.0)
            //    nx=nx;
            
            pfValues++;
            pfValuesOut += m_nWidth;
        };
    };
    
    m_nNumContribAtLastRecompute = m_nNumContrib;
    
    return 0;
}

int CSphericalMesh::nGet(double fPolar,double fAzi,double* pfItems) {
    int nVar;
    int nMesh;
    double* pfInterp;

    if (bReadyToAdd())
        return 1;

    nMesh = ((int) floor((fAzi/m_fStep0Azi) + 0.5)) + m_nDim0Azi*((int) floor((fPolar/m_fStep1Polar) + 0.5));

    pfInterp = &m_afMeshValues[0] + nMesh*m_nWidth;

    for (nVar = 0; nVar < m_nWidth; nVar++) {
        pfItems[nVar] = pfInterp[nVar];
    };

    return 0;
};

