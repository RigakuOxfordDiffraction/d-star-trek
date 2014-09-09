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
// CBackground.cc     Initial author: Thaddeus Niemeyer        15-Nov-2001
//    Background tools
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


#include "CBackground.h"

#include "Cstat.h"
#include <string.h>
#include "raxis.h"
#include "dtarray.h"
#include "Cdetector.h"
#include "Csource.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

Cbackground::Cbackground(Cimage& oImage):m_oCi(oImage) {
    m_pfBackground = NULL;
    m_pcBackground = NULL;
    m_nDim[0] = m_oCi.nGetDimension(0);
    m_nDim[1] = m_oCi.nGetDimension(1);
    m_nMaskValue = 127;
    m_nTempValue = 255;
}

unsigned char Cbackground::ms_cZero = 0;
unsigned char&   Cbackground::rcGetBackground(int px0,int px1)     { if (!m_pcBackground) return ms_cZero; else return m_pcBackground[((px0)+(px1)*(m_nDim[0]))]; };
unsigned char&   Cbackground::rcGetBackground(int px0,int px1,unsigned char* pcUser)     { return pcUser[((px0)+(px1)*(m_nDim[0]))]; };


Cbackground::~Cbackground()
{
    if( m_pfBackground )
    {
        delete[] m_pfBackground;
        m_pfBackground = NULL;
    }
    
    if( m_pcBackground )
    {
        delete[] m_pcBackground;
        m_pcBackground = NULL;
    }
}

int Cbackground::nAllocateBackground() {
    if (m_pcBackground)
        return 1;
    m_pcBackground= new unsigned char[m_nDim[0]*m_nDim[1]];
    memset(m_pcBackground,0,m_nDim[0]*m_nDim[1]);
    return 0;
}


int Cbackground::nClearMasks(int nNewValue) {
    int a2n[2];
    nAllocateBackground();
    for (a2n[0]=0;a2n[0]<m_nDim[0];a2n[0]++) {
        for (a2n[1]=0;a2n[1]<m_nDim[1];a2n[1]++) {
            if (rcGetBackground(a2n[0],a2n[1]) == m_nMaskValue)
                rcGetBackground(a2n[0],a2n[1]) = nNewValue;
        }
    }
    return 0;
}

int Cbackground::nClearAll() {
    int a2n[2];
    nAllocateBackground();
    for (a2n[0]=0;a2n[0]<m_nDim[0];a2n[0]++) {
        for (a2n[1]=0;a2n[1]<m_nDim[1];a2n[1]++) {
            rcGetBackground(a2n[0],a2n[1]) = 0;
        }
    }
    return 0;
}


int Cbackground::nMaskRegion(int na0,int na1,int nb0,int nb1,bool bInvert) {
    int a2n[2];
    if (na0>nb0) 
        std::swap(na0,nb0);
    if (na1>nb1)
        std::swap(na1,nb1);
    na0 = min(m_nDim[0],max(0,na0));
    nb0 = min(m_nDim[0],max(0,nb0));
    na1 = min(m_nDim[1],max(0,na1));   
    nb1 = min(m_nDim[1],max(0,nb1));
    nAllocateBackground();
    if (bInvert) {
        for (a2n[0]=0;a2n[0]<m_nDim[0];a2n[0]++) {
            for (a2n[1]=0;a2n[1]<m_nDim[1];a2n[1]++) {
                if ((a2n[0]<na0) || (a2n[0]>nb0) || (a2n[1]<na1) || (a2n[1]>nb1))
                    rcGetBackground(a2n[0],a2n[1]) = m_nMaskValue;
            }
        }
    } else {
        for (a2n[0]=na0;a2n[0]<=nb0;a2n[0]++) {
            for (a2n[1]=na1;a2n[1]<=nb1;a2n[1]++) {
                rcGetBackground(a2n[0],a2n[1]) = m_nMaskValue;
            }
        }
    }
    return 0;
}

/*  This function iteratively passes over the image and adds 2 3 ... nLevels to pixels that are adjacent
    at these distances.

*/

int Cbackground::nBuildAdjacentToPeaks(int nLevels) {
    int nLevel;
    int npx[2];
    int nx,ny;

    if (!m_pcBackground)
        return 1;
    for (nLevel=2;nLevel<=nLevels;nLevel++) {
        printf("\rPass %d of %d",nLevel,nLevels);
        // At each pass, we only consider pixels that have a value of (nLevel-1)
        for (npx[0]=1;npx[0]<m_nDim[0]-1;npx[0]++) {
            for (npx[1]=1;npx[1]<m_nDim[1]-1;npx[1]++) {
                if (rcGetBackground(npx[0],npx[1])==nLevel-1) {
                    // See all pixels in the neighborhood that are 0 to the value nLevel.
                    for (nx=-1;nx<=1;nx++) {
                        for (ny=-1;ny<=1;ny++) {
                            if (rcGetBackground(npx[0]+nx,npx[1]+ny)==0)
                                rcGetBackground(npx[0]+nx,npx[1]+ny) = nLevel;
                        }
                    }
                }
            }
        }
    }
    printf("\n");
    return 0;
}

int Cbackground::nFatten(int nFattenCount) {
    int npx[2];
    nBuildAdjacentToPeaks(nFattenCount+1);
    for (npx[0]=1;npx[0]<m_nDim[0]-1;npx[0]++) {
        for (npx[1]=1;npx[1]<m_nDim[1]-1;npx[1]++) {
            if ((rcGetBackground(npx[0],npx[1])>1) && (rcGetBackground(npx[0],npx[1])<=nFattenCount+1))
                rcGetBackground(npx[0],npx[1]) = 1;
        }
    }

    return 0;
}




int Cbackground::nFindSpot(int nPix0,int nPix1,itr<int>& anPix0,itr<int>& anPix1,itr<int>& anStack0,itr<int>& anStack1,int& nCount,int& nMaxCount) {
	int nx,ny;
	int n0,n1;

	anStack0.clear();
	anStack1.clear();
	anPix0.clear();
	anPix1.clear();

	anStack0.push(nPix0);
	anStack1.push(nPix1);
	nCount = 0;
    nMaxCount = 0;
	while (!anStack0.isempty()) {
		n0 = anStack0.pop();
		n1 = anStack1.pop();
		if ((nx=rcGetBackground(n0,n1))!=m_nTempValue) {
			anPix0 + n0;
			anPix1 + n1;
            nMaxCount = max(nx,nMaxCount);
			rcGetBackground(n0,n1) = m_nTempValue;
			nCount += nx;
		}
		
		for (nx=-1;nx<=1;nx++) {
			for (ny=-1;ny<=1;ny++) {
				if ((n0 + nx>=0) && (n0 + nx < m_nDim[0]) &&
					(n1 + ny>=0) && (n1 + ny < m_nDim[1]) &&
					(rcGetBackground(n0 + nx,n1 + ny)>=1) && (rcGetBackground(n0 + nx,n1 + ny)<m_nMaskValue)) {
					anStack0.push(n0 + nx);
					anStack1.push(n1 + ny);
				}
			}
		}
	}

	return 0;
}

int Cbackground::nFindSpots(itr<int>& anPix0,itr<int>& anPix1,itr<int>& anDim0,itr<int>& anDim1,itr<int>& anTotal,itr<int>& anCount,itr<int>& anMax,itr<int>* panCumulPix0,itr<int>* panCumulPix1) {
	int n0,n1;
	int nx;
	itr<int> anStack0;
	itr<int> anStack1;
	itr<int> anLocalPix0;
	itr<int> anLocalPix1;

	int n0Avg,n1Avg,nTotal,nCount,nMaxCount;
	int nPix0,nPix1;
	int nMax0,nMin0,nMax1,nMin1;

    -anTotal;	
	-anPix0;
	-anPix1;
	-anDim0;
	-anDim1;
    -anCount;
    -anMax;
    if (panCumulPix0)
        -(*panCumulPix0);
    if (panCumulPix1)
        -(*panCumulPix1);
	for (n0 = 0; n0 < m_nDim[0]; n0++) {
        for (n1 = 0; n1 < m_nDim[1] ;n1++) {
			if ((rcGetBackground(n0,n1)>=1) && (rcGetBackground(n0,n1)<m_nMaskValue)) {
				nFindSpot(n0,n1,anLocalPix0,anLocalPix1,anStack0,anStack1,nTotal,nMaxCount);
				n0Avg = 0;
				n1Avg = 0;
				nCount = 0;
				for (nx = 0; nx < anLocalPix0.size(); nx++) {
					n0Avg += (nPix0 = anLocalPix0[nx]);
					n1Avg += (nPix1 = anLocalPix1[nx]);
					if (nCount == 0) {
						nMax0 = nMin0 = nPix0;
						nMax1 = nMin1 = nPix1;
					} else {
						nMax0 = max(nPix0,nMax0);
						nMax1 = max(nPix1,nMax1);
						nMin0 = min(nPix0,nMax0);
						nMin1 = min(nPix1,nMax1);
					}
					nCount++;
				}
				anPix0 + (n0Avg/nTotal);
				anPix1 + (n1Avg/nTotal);
				anDim0 + (nMax0 - nMin0);
				anDim1 + (nMax1 - nMin1);
                anTotal + nTotal;
                anCount + nCount;
                anMax + nMaxCount;
                if (panCumulPix0)
                    (*panCumulPix0) + anLocalPix0;
                if (panCumulPix1)
                    (*panCumulPix1) + anLocalPix1;
			}
		}
	}
    
        
    for (n0 = 0; n0 < m_nDim[0]; n0++) {
        for (n1 = 0; n1 < m_nDim[1] ;n1++) {
            if (rcGetBackground(n0,n1)==m_nTempValue)
                rcGetBackground(n0,n1) = 1;
        }
    }



	return 0;
}


int Cbackground::nPrintBackground(Cstring& sName) {
    int npx[2];

    unsigned short int uiPix;

    printf("Printing background pixel map...\n");
    uiPix=0;

    for (npx[0]=0;npx[0]<m_nDim[0];npx[0]++) {
        for (npx[1]=0;npx[1]<m_nDim[1];npx[1]++) {
            m_oCi.nSetPixel(npx[0],npx[1],(unsigned short int) rfGetBackground(npx[0],npx[1]));
        }
    }
    
    if (m_oCi.nWrite(sName)) 
        cout << "Could not write " << sName << "\n";

return 0;
}

int Cbackground::nPrintPeaks(Cstring& sName)
{
    int npx[2];
    Cimage oImageOut(m_oCi);

    for (npx[0]=0;npx[0]<m_nDim[0];npx[0]++) {
        for (npx[1]=0;npx[1]<m_nDim[1];npx[1]++) {
            if ((rcGetBackground(npx[0],npx[1])>=1) && (rcGetBackground(npx[0],npx[1])<m_nMaskValue))
                oImageOut.nSetPixel(npx[0],npx[1],(unsigned short int) rcGetBackground(npx[0],npx[1]));
            else
                oImageOut.nSetPixel(npx[0],npx[1],(unsigned short int) (m_oCi.*m_oCi.prfGetPixel)(npx[0],npx[1]));
        }
    }
    if (oImageOut.nWrite(sName)) 
        cout << "Could not write " << sName << "\n";

    return 0;
}

int Cbackground::nPrintAdjacentToPeaks(Cstring& sName) {
    int npx[2];
    int nMaxVal;

    // Find the maximum value in the mask.
    nMaxVal = 0;
    for (npx[0]=0;npx[0]<m_nDim[0];npx[0]++) {
        for (npx[1]=0;npx[1]<m_nDim[1];npx[1]++) {
            nMaxVal = max(nMaxVal,(int) rcGetBackground(npx[0],npx[1]));
        }
    }

    for (npx[0]=0;npx[0]<m_nDim[0];npx[0]++) {
        for (npx[1]=0;npx[1]<m_nDim[1];npx[1]++) {
            if (rcGetBackground(npx[0],npx[1]))
                m_oCi.nSetPixel(npx[0],npx[1],(unsigned short int) rcGetBackground(npx[0],npx[1]));
        }
    }
    if (m_oCi.nWrite(sName)) 
        cout << "Could not write " << sName << "\n";

    return 0;
}

int Cbackground::nExcludeZingers(int nMinAdjacent,int nNewValue) {
    unsigned char* pcExclude;
    int npx[2];
    int nAdjacent;
    int nx,ny;

    pcExclude = new unsigned char[m_nDim[0]*m_nDim[1]];
    memset(pcExclude,0,m_nDim[0]*m_nDim[1]);
    
    for (npx[0]=0;npx[0]<m_nDim[0];npx[0]++) {
        for (npx[1]=0;npx[1]<m_nDim[1];npx[1]++) {
            if (bInRange(npx[0],npx[1]) && (rcGetBackground(npx[0],npx[1])==1)) {
                // Count up pixels in adjacent boxes.
                nAdjacent = 0;
                for (nx=-1;nx<=1;nx++) {
                    for (ny=-1;ny<=1;ny++) {
                        if (bInRange(npx[0]+nx,npx[1]+ny) && (rcGetBackground(npx[0]+nx,npx[1]+ny)==1))
                            nAdjacent++;
                    }
                }
                if (nAdjacent<nMinAdjacent) 
                    rcGetBackground(npx[0],npx[1],pcExclude)=nNewValue;
            }
        }
    }
    for (npx[0]=0;npx[0]<m_nDim[0];npx[0]++) {
        for (npx[1]=0;npx[1]<m_nDim[1];npx[1]++) {

        }
    }

    delete[] pcExclude;
    return 0;
}

int Cbackground::nExcludeAdjacentToExcluded(int nSearchValue) {
    bool bChangesMade;
    bool bFound;
    int npx[2];
    int nx,ny;

    printf("Removing peaks located next to the boarder.\n");
    do {
        bChangesMade = FALSE;

        for (npx[0]=0;npx[0]<m_nDim[0];npx[0]++) {
            for (npx[1]=0;npx[1]<m_nDim[1];npx[1]++) {
                if (rcGetBackground(npx[0],npx[1])==nSearchValue) {
                    bFound = FALSE;
                    for (nx=-1;nx<=1;nx++) {
                        for (ny=-1;ny<=1;ny++) {
                            if (!bInRange(npx[0]+nx,npx[1]+ny)) {
                                bFound = TRUE;
                                break;
                            }
                        }
                        if (bFound)
                            break;
                    }
                    if (bFound) {
                        bChangesMade = TRUE;
                        rcGetBackground(npx[0],npx[1]) = m_nMaskValue;
                    }
                }
            }
        }
    } while (bChangesMade);

    return 0;
}


float Cbackground::fGetSmoothBackground(int nPix0,int nPix1)
{
    int a2x2nPixBounds[2][2];
    float fRatio;
    float fValue;
    int nx,ny;

    fValue = 0.0;
    a2x2nPixBounds[0][0] = nPix0/m_nBackgroundParams[0];
    a2x2nPixBounds[0][1] = min(m_nBackgroundParams[2] - 1,a2x2nPixBounds[0][0]+1);
    a2x2nPixBounds[1][0] = nPix1/m_nBackgroundParams[1];
    a2x2nPixBounds[1][1] = min(m_nBackgroundParams[3] - 1,a2x2nPixBounds[1][0]+1);
    
	float	fTemp1 = 0.0f;
	float	fTemp2 = 0.0f;
	
	for (nx=0;nx<2;nx++)
	{
        for (ny=0;ny<2;ny++)
		{
            fTemp1 = (float)(nPix0 - m_nBackgroundParams[0]*a2x2nPixBounds[0][!nx]);
			fTemp2 = (float)(nPix1 - m_nBackgroundParams[1]*a2x2nPixBounds[1][!ny]);
			
			fRatio = fabs(fTemp1)/(float)m_nBackgroundParams[0] * fabs(fTemp2)/(float)m_nBackgroundParams[1];

			fValue += fRatio * m_pfBackground[a2x2nPixBounds[1][ny]*m_nBackgroundParams[2] + a2x2nPixBounds[0][nx]];
        }
    }
    return fValue;
}

// rb: this routine seems to bin all "peak" (i.e. non-backround) pixels in terms of their two-theta values

int Cbackground::nIntegrate(int                 nBoxSize,
                            itr<double>&        af2Theta,
                            double              f2ThetaStep,
                            Cdetector*          poDetector,
                            Csource*            poSource,
                            int                 nCompression)
{
    // Set size to "Theta arrays"
    int     nBins = (int) (180.0 / f2ThetaStep) + 1;
        
    itr<int>                an2ThetaCounts;
        an2ThetaCounts.setsize(nBins);
        
    itr<double>             af2Theta2;
        af2Theta2.setsize(nBins);
        
    itr<double>             af2ThetaAvg;
        af2ThetaAvg.setsize(nBins);
        
    itr<double>             af2ThetaDev;
        af2ThetaDev.setsize(nBins);

    af2Theta.setsize(nBins);
    ///////////////////////////////////////////

    nAllocateBackground();      // rb: this routine actully allocates not the background, but the storage for background - m_pcBackground
    
    if( nBuildPeakPixels(nBoxSize,true) )   //rb:  this routine classifies pixels as either "peak" or "background"
                                            // and fills in m_pcBackground and m_pfBackground arrays accordingly
      return 1;
        
    float           a3fS0[3];
    poSource->vCalcGetS0(&a3fS0[0]);
    double          fWavelength = poSource->fGetWavelength();  
    vMulVec3DScalar(a3fS0,1.0/fWavelength,a3fS0);
        
    
    int             nBin = 0;
    double          fValue = 0.0;
    double          fReso = 0.0;
    double          f2Theta = 0.0;
        
    for(int nPass = 0; nPass < 2; nPass ++) 
      {
        if(nPass == 1)
          {
            for (nBin = 0; nBin < nBins; nBin++)
              {
                if (an2ThetaCounts[nBin]>=10)
                  {
                    af2ThetaAvg[nBin] = af2Theta[nBin]/an2ThetaCounts[nBin];
                    af2ThetaDev[nBin] = sqrt( (af2Theta2[nBin]*an2ThetaCounts[nBin] - af2Theta[nBin]*af2Theta[nBin])
                                              / an2ThetaCounts[nBin] / (an2ThetaCounts[nBin]-1) );
                  }
                else 
                  {
                    af2ThetaAvg[nBin] = 0.0;
                    af2ThetaDev[nBin] = 0.0;
                  }
              }
          } 
        else 
          {
            memset(&af2ThetaAvg[0],0,sizeof(double)*nBins);
            memset(&af2ThetaDev[0],0,sizeof(double)*nBins);
          }
                
        memset(&af2Theta[0],0,sizeof(double)*nBins);
        memset(&af2Theta2[0],0,sizeof(double)*nBins);
        memset(&an2ThetaCounts[0],0,sizeof(int)*nBins);
                                
        for(int n1 = 0; n1 < m_nDim[1]; n1++) 
          {
            for(int n0 = 0; n0 < m_nDim[0]; n0++)
              {
                fValue = m_pfBackground[n0 + n1*m_nDim[0]]; 
                // Or we can use the pixel values themselves as in ... (m_oCi.*m_oCi.prfGetPixel)(n0,n1);
                if ((fValue>0) && (fValue < 65500))
                  {
                    fReso = poDetector->fCalcGetResolution(n0*nCompression,n1*nCompression,a3fS0);
                    f2Theta = ((asin(min(1.0,fWavelength/2.0/fReso)))*2.0/Gs_dRADIANS_PER_DEGREE);
                    if ((f2Theta>=0.0) && (f2Theta<180))
                      {
                        nBin = int(f2Theta/f2ThetaStep);
                        if ((nPass==0) || (af2ThetaDev[nBin]==0.0) || (ABS(fValue - af2ThetaAvg[nBin])<3.0*af2ThetaDev[nBin]))
                          {
                            an2ThetaCounts[nBin] ++;
                            af2Theta[nBin]  += fValue;
                            af2Theta2[nBin] += fValue * fValue;
                          }
                      }
                  }
              }                         
          }
      }
        
    for (nBin = 0; nBin < nBins; nBin++)
      {
        if (an2ThetaCounts[nBin])
          af2Theta[nBin] /= an2ThetaCounts[nBin];
      }

    return 0;
}

int Cbackground::nFind2D(Creflnlist& oList,
                         const double fRot,
                         double fSignificanceLevel,
                         Cdetector* poDetector,
                         Csource* poSource) 
{
    int nx,ny;
    int nSpot;
    itr<int> anPos0;
    itr<int> anPos1;
	itr<int> anDim0;
	itr<int> anDim1;
    itr<int> anCount;
    itr<int> anTotal;
    itr<int> anMax;
    itr<int> anCumulPix0;
    itr<int> anCumulPix1;
    itr<int> anDelete;
    float fReso;
    float a3fS0[3];

    oList.vDeleteAll();
    Crefln oRefln(&oList);
    if (poSource) 
        poSource->vCalcGetS0(&a3fS0[0]);  // Required for resolution calc below

    nAllocateBackground();
    if (nBuildPeakPixels(30))
        return 1;

    //nPrintPeaks((Cstring) "test_before.img");

	if (nFindSpots(anPos0,anPos1,anDim0,anDim1,anTotal,anCount,anMax,&anCumulPix0,&anCumulPix1))
		return 1;
    
    
    // Do a further restriction by removing entries that don't have a significant # of pixels.
        // Compute the standard deviation and average for all pixels
    double fSumX;
    double fSumX2;
    double fSum;
    double fAverage;
    double fDeviation;
    int nPass;
    int nOffset;
    int nTotal,nCount,nMax;

    

    fAverage = 0.0;
    fDeviation = 0.0;
    for (nPass = 0; nPass < 2; nPass++) {
        fSumX = 0.0;
        fSumX2 = 0.0;
        fSum = 0.0;
        nOffset = 0;

        for (nx = 0; nx < anTotal.size(); nx++) {
            nTotal = anTotal[nx];
            nCount = anCount[nx];
            nMax = anMax[nx];
            fSumX += nMax;
            fSumX2 += nMax*nMax;
            fSum ++;
            if (nPass == 1) {
                if (nMax - fAverage >= fDeviation*2.0) {
                    anDelete + 0;
                } else {
                    // Remove the pixels as well.
                    for (ny = 0; ny < nCount; ny++)
                        rcGetBackground(anCumulPix0[nOffset + ny],anCumulPix1[nOffset + ny]) = 0;
                    anDelete + 1;
                }

                    
            }
            nOffset += nCount;
        }
        if (!fSum)
            return 0;
        fAverage = fSumX/fSum;
        fDeviation = sqrt((fSumX2*fSum-fSumX*fSumX)/fSum/(fSum-1));
    }

    
    //nPrintPeaks((Cstring) "test_after.img");


    for (nSpot = 0; nSpot < anPos0.size(); nSpot++) {
        
        if (!anDelete[nSpot]) {
            
            
            if (oRefln.bEllipsoidAvailable(false)) {
                double a2x2fEllipsoidA[2][2];
                double a2fEllipsoidb[2];
                double fEllipsoidc;
                double fRadius;
                fRadius = max(anDim0[nSpot],anDim1[nSpot])/2.0;
                nBuildEllipse(fRadius,&a2x2fEllipsoidA[0][0],&a2fEllipsoidb[0],&fEllipsoidc);
                oRefln.nPutGetEllipse(true,&a2x2fEllipsoidA[0][0],&a2fEllipsoidb[0],&fEllipsoidc);
            }
            
            
            
            oRefln.vSetField(oList.m_nFI_fObsPx0,
                anPos0[nSpot]+ (float) 0.5);
            oRefln.vSetField(oList.m_nFI_fObsPx1,
                anPos1[nSpot]+ (float) 0.5);
            oRefln.vSetIntensity(0.0);
            oRefln.vSetSigmaI(0.0);
            oRefln.vSetField(oList.m_nFI_fObsRotMid, (float) fRot);
            if ((poDetector) && (poSource)) {
                fReso = poDetector->fCalcGetResolution(
                    oRefln.fGetField(oList.m_nFI_fObsPx0),
                    oRefln.fGetField(oList.m_nFI_fObsPx1),
                    a3fS0);
                oRefln.vSetField(oList.m_nFI_fDetResolution,(float) fReso);
            }
            
            
            
            oList.nInsert(&oRefln);
        }
    }
    
    return 0;
}


int Cbackground::nExcludeMinimum(double fMinimum) {
    double f0;
    int n0,n1;

    for (n0 = 0; n0 < m_nDim[0]; n0++) {
        for (n1 = 0; n1 < m_nDim[1] ;n1++) {
            f0 = (m_oCi.*m_oCi.prfGetPixel)(n0,n1);
            if (f0<fMinimum)
                rcGetBackground(n0,n1) = 1;
        }
    }
    return 0;
}




int Cbackground::nBuildPeakPixels(int nBoxDim, bool bSaveBackground)
{
  int             nPixCount = m_nDim[0]*m_nDim[1];
    
  itr<int>        anSumXHoriz;                           
   anSumXHoriz.setsize(nPixCount);
  int*            pnSumXHoriz = &anSumXHoriz[0]; // NOTE : THE ORDER IS IMPORTANT: FIRST SET SIZE, THEN TAKE ADDRESS
                                                       // THE REASON: THE ARRAY MAY GET REALLOCATED, WHEN SETTING SIZE 
  itr<int>        anSumX2Horiz;
   anSumX2Horiz.setsize(nPixCount);
  int*            pnSumX2Horiz = &anSumX2Horiz[0];
        
  itr<int>        anSumCHoriz;
   anSumCHoriz.setsize(nPixCount);
  int*            pnSumCHoriz = &anSumCHoriz[0];
    
  itr<int>        anRollingX;
   anRollingX.setsize(nBoxDim);
  int*            pnRollingX = &anRollingX[0];
        
  itr<int>        anRollingX2;
   anRollingX2.setsize(nBoxDim);
  int*            pnRollingX2 = &anRollingX2[0];
        
    itr<int>        anRollingSel;
        anRollingSel.setsize(nBoxDim);
        int*            pnRollingSel = &anRollingSel[0];
   
    itr<float>      afPixel;
        afPixel.setsize(nPixCount);
    float*          pfPixel = &afPixel[0];
    //////////////////////////////////////////////      

    if( bSaveBackground )
      {
        if( m_pfBackground )
          delete [] m_pfBackground;
        
        m_pfBackground = new float[nPixCount];
        memset(m_pfBackground,0,sizeof(float)*nPixCount);
      }

    if (!m_pcBackground)
      nAllocateBackground();
    memset(m_pcBackground,0,m_nDim[0]*m_nDim[1]);

    int         nSelect = 0;
    int         nValue = 0;
    float       fValue = 0.0;
    int         nSumHorizPoint = 0;
    int         nPixelPoint = 0;
        
    int         nMinContrib = max(3, nBoxDim*nBoxDim/10);
    
    int         nRollingX = 0;
    int         nRollingX2 = 0;
    int         nRollingC = 0;
    int         nRollingPoint = 0;
    int         n1 = 0;

    for(int nPass = 0; nPass < 2; nPass++) 
      {
        // Build up horizontal summations.
        nPixelPoint = 0;         
        nSumHorizPoint = 0;
        for (n1 = 0,nPixelPoint = 0; n1 < m_nDim[1]; n1++) 
          {
            // nPixelPoint is a pixel coordinate in a 1-dimensional pixel array
            // n1, n0 are pixel coordinates in a 2-dimensinal image
            
            nRollingPoint = 0;
            nRollingX = 0;
            nRollingX2 = 0;
            nRollingC = 0;
            for (int n0 = 0; n0 < m_nDim[0]; n0++,nPixelPoint++)
              {
                /////////////////////////////////////////////////
                // rb: determine fvalue 
                if (nPass == 0)
                  {
                    fValue = (m_oCi.*m_oCi.prfGetPixel)(n0,n1);
                    pfPixel[nPixelPoint] = fValue;
                  }
                else
                  fValue = pfPixel[nPixelPoint];
                nValue = (int)fValue;
                /////////////////////////////////////////////////

                nSelect = !rcGetBackground(n0,n1) && fValue > 0 && fValue < 65500 ? 1 : 0;  //rb: is 65500 a saturated value???
                pnRollingSel[nRollingPoint] = nSelect;
                                
                if( nSelect ) 
                  {
                    pnRollingX[nRollingPoint] = nValue;
                    pnRollingX2[nRollingPoint] = nValue*nValue;
                    nRollingX += nValue;
                    nRollingX2 += nValue*nValue;
                    nRollingC++;
                  } 
                else 
                  {
                    pnRollingX[nRollingPoint] = 0;
                    pnRollingX2[nRollingPoint] = 0;
                  }

                if (n0>=nBoxDim-1)
                  {
                    pnSumXHoriz[nSumHorizPoint] = nRollingX;
                    pnSumX2Horiz[nSumHorizPoint] = nRollingX2;
                    pnSumCHoriz[nSumHorizPoint] = nRollingC;
                    nSumHorizPoint++;
                                        
                    nRollingPoint = ((nRollingPoint + 1) % nBoxDim);   // rb: this forces the "rolling point" to "oscillate" in [0, nBoxDim-1]
                    
                    if (pnRollingSel[nRollingPoint])
                      {
                        nRollingX -= pnRollingX[nRollingPoint];
                        nRollingX2 -= pnRollingX2[nRollingPoint];
                        nRollingC--;
                      }
                  } 
                else
                  nRollingPoint = ((nRollingPoint + 1) % nBoxDim);
              }
          }
                
        double      fAverage = 0.0;
        double      fDeviationSq = 0.0;
        double      f0 = 0.0;
                
        for (int n0 = 0; n0 < m_nDim[0]; n0++)
          {
            nRollingPoint = 0;
            nRollingX = 0;
            nRollingX2 = 0;
            nRollingC = 0;
            nSumHorizPoint = max(0,min(m_nDim[0] - nBoxDim,n0 - nBoxDim/2));
            nRollingPoint = 0;
            for (int nx = 0; nx < nBoxDim; nx++)
              {
                nRollingX += pnSumXHoriz[nSumHorizPoint];
                nRollingX2 += pnSumX2Horiz[nSumHorizPoint];
                nRollingC += pnSumCHoriz[nSumHorizPoint];
                                
                pnRollingX[nRollingPoint] = pnSumXHoriz[nSumHorizPoint];
                pnRollingX2[nRollingPoint] = pnSumX2Horiz[nSumHorizPoint];
                pnRollingSel[nRollingPoint] = pnSumCHoriz[nSumHorizPoint];
                                
                nSumHorizPoint += (m_nDim[0] - nBoxDim + 1);
                                
                nRollingPoint = ((nRollingPoint + 1) % nBoxDim);
              }
                        
            for (n1 = 0,nPixelPoint = n0; n1 < m_nDim[1]; n1++,nPixelPoint+=m_nDim[0])
              {
                if (nRollingC >= nMinContrib)
                  {
                    fAverage = ((double) nRollingX)/nRollingC;
                    fDeviationSq = (((double) nRollingX2)*nRollingC-((double) nRollingX)*nRollingX)/nRollingC/(nRollingC-1);
                    
                    f0 = (afPixel[nPixelPoint] - fAverage);
                                        
                    if ((f0*f0>fDeviationSq*2.0*2.0) && (f0>0.0))
                                                rcGetBackground(n0,n1) = 1;
                    if (bSaveBackground)
                      m_pfBackground[n0 + n1*m_nDim[0]] = (float) fAverage;
                  } 
                else
                  {
                    rcGetBackground(n0,n1) = 0;
                  }
                                        
                if ((n1 >= nBoxDim/2) && (n1 + nBoxDim/2 < m_nDim[1]))
                  {
                    nRollingX += pnSumXHoriz[nSumHorizPoint];
                    nRollingX2 += pnSumX2Horiz[nSumHorizPoint];
                    nRollingC += pnSumCHoriz[nSumHorizPoint];
                    nRollingX -= pnRollingX[nRollingPoint];
                    nRollingX2 -= pnRollingX2[nRollingPoint];
                    nRollingC -= pnRollingSel[nRollingPoint];
                    pnRollingX[nRollingPoint] = pnSumXHoriz[nSumHorizPoint];
                    pnRollingX2[nRollingPoint] = pnSumX2Horiz[nSumHorizPoint];
                    pnRollingSel[nRollingPoint] = pnSumCHoriz[nSumHorizPoint];
                    nRollingPoint = ((nRollingPoint + 1) % nBoxDim);
                    nSumHorizPoint += (m_nDim[0] - nBoxDim + 1);
                  }
              }
          }
      }
    return 0;
}
