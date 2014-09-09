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
// CGonioMask.cpp           Initial author: Thaddeus Niemeyer Jan 2003
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


#include "CGonioMask.h"
#include "Cgoniometer.h"
#include "Cscan.h"
#include "Cimage.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif

static bool bGetImagePixelSize(Cimage_header* poHeader, double daPixSize[2]);

CGonioMask::CGonioMask() {
    m_nStage = m_nOmegaStage;
    m_fScale = 1.0;
    m_poDetector = NULL;
    m_poCrysGonio = NULL;
    m_poScan = NULL;
    m_pcMask = NULL;
    m_dDetDistance = 100.0;
    m_dDetSwing    = 0.0;
};

CGonioMask::~CGonioMask() {
    if (NULL != m_poDetector)
      {
        delete m_poDetector;
	m_poDetector = NULL;
      }
    if (NULL != m_poCrysGonio)
      {
        delete m_poCrysGonio;
	m_poCrysGonio = NULL;
      }
    if (NULL != m_poScan)
      {
        delete m_poScan;
	m_poScan = NULL;
      }
    if (NULL != m_pcMask)
      {
        delete[] m_pcMask;
	m_pcMask = NULL;
      }
    
    for(int ii=0; ii < m_vpShadows.size(); ii++)
    {
        delete (CGonioShadow*)m_vpShadows[ii];
        m_vpShadows[ii] = NULL;
    }
    m_vpShadows.clear();
};

int CGonioMask::nAddBox(itr<double>& afPosX,itr<double>& afPosY,itr<double>& afPosZ)
{
    
    int nVariable = 0;
    int nx;
    int nPoints;
    itr<double>* afPos[3] = { &afPosX,&afPosY,&afPosZ };

    nPoints = afPosX.size();

    if ((afPosX.size()!=afPosY.size()) || (afPosY.size()!=afPosZ.size())) 
      return 1;
    if (afPosX.size() != 4)
      return 1;
    m_anStage + m_nStage;
    for (nx=0;nx<4;nx++)
      {
	m_a4anBoxX[nx] + (int) (m_fScale*afPosX[nx]);
	m_a4anBoxY[nx] + (int) (m_fScale*afPosY[nx]);
	m_a4anBoxZ[nx] + (int) (m_fScale*afPosZ[nx]);
      }
    return 0;
};

void CGonioMask::vSetOmegaStage() {
    m_nStage = m_nOmegaStage;
};
void CGonioMask::vSetChiStage() {
    m_nStage = m_nChiStage;
};
void CGonioMask::vSetPhiStage() {
    m_nStage = m_nPhiStage;
};

void CGonioMask::vSetScale(double fScale) {
    m_fScale = fScale;
};

int CGonioMask::nReadGonio(double fOmega,double fChi,double fPhi) {
    float a3fValues[3];

    a3fValues[0] = fOmega;
    a3fValues[1] = fChi;
    a3fValues[2] = fPhi;

    if (!m_poCrysGonio) {
        printf("ERROR:  Need to extract detector information from an actual image!\n");
        return 1;
    };
    m_poCrysGonio->nSetDatum(3,&a3fValues[0]);
    return 0;
};

int CGonioMask::nReadGonio(Cstring& sImageName) {
    Cimage_header oHeader(sImageName);
	return nReadGonio(oHeader);
};

int CGonioMask::nReadGonio(Cimage_header& oHeader) {
    Cstring sDetName;
    int nx;
    
    if (!oHeader.bIsAvailable())
        return 1;

    sDetName = "";    
    (void) oHeader.nGetValue(Cdetector::ms_sDetectorNames, 1, &sDetName);    

    // The next line is a real killer if the spatial distortion uses interpolation tables
    // since they will be read in each time the next statement is executed.

    m_poDetector = new Cdetector(oHeader, sDetName, TRUE, FALSE);
    if (!m_poDetector->bIsAvailable())
        return 1;
    m_poCrysGonio = new Cgoniometer(oHeader, Ccrystal::ms_sCrystalPrefix);
    if (!m_poCrysGonio->bIsAvailable())
        return 1;
    m_poScan = new Cscan(oHeader);
    if (!m_poScan->bIsAvailable())
        return 1;
    nx = m_poDetector->m_a2nDim[0]*m_poDetector->m_a2nDim[1];
    m_pcMask = new unsigned char[nx];
    return 0;
};
int CGonioMask::nWriteMask(Cstring& sNameOut,Cstring& sNameIn,int nMaskValue) {
    Cimage* poImage;

    if (sNameIn.length())
        poImage = new Cimage(sNameIn);
    else
        poImage = new Cimage(m_poDetector->m_a2nDim[0],m_poDetector->m_a2nDim[1],eImage_uI2);

    Cimage& oImage = *poImage;

    int n0,n1;
    int nValue;

    
    for (n1=0;n1<m_poDetector->m_a2nDim[1];n1++) {
        for (n0=0;n0<m_poDetector->m_a2nDim[0];n0++) {
            if (m_pcMask[n0 + n1*m_poDetector->m_a2nDim[0]] || (!sNameIn.length())) {
                if (nMaskValue < 0)
                    nValue = m_pcMask[n0 + n1*m_poDetector->m_a2nDim[0]];
                else
                    nValue = nMaskValue;
                oImage.nSetPixel(n0,n1,(unsigned short int) nValue);
            };
        };
    };
    oImage.nWrite(sNameOut);
    delete poImage;
    return 0;
};

int
CGonioMask::nApplyMask(Cimage& oImage)
{
  int n0, n1;
  for (n1=0; n1<m_poDetector->m_a2nDim[1]; n1++)
    {
      for (n0=0; n0<m_poDetector->m_a2nDim[0]; n0++)
	{
	  if (m_pcMask[n0 + n1*m_poDetector->m_a2nDim[0]])
	    oImage.nSetPixel(n0,n1,(unsigned short int) 0);
        }
    }
  return 0;
}


int CGonioMask::nWriteTestMask(Cstring& sInputImage,double fValue) {
    Cstring sOutputImage;

    sOutputImage = sInputImage.before('.');
    sOutputImage += "_masked.";
    sOutputImage += sInputImage.after('.');
    if (nReadGonio(sInputImage))
        return 1;
    if (nBuildMask()) 
        return 1;
    if (nWriteMask(sOutputImage,sInputImage,(int) fValue))
        return 1;
    return 0;
};



int CGonioMask::nBuildMask()
{
    int nBox;
    int nEdge;
    int nCount;
    int nPoints;
    int nStart,nEnd;
    int nStat;
    const double fMaxEdgeStep = 1;	// In MM.
    const int nMaxPixelStep = 30;	// Maximum pixels to step for a triangular piece.
    double a3fPoint[3];		// Point along edge.
    double a3fPointStart[3];	// Edge Start
    double a3fPointEnd[3];	// Edge End
    double a3fNormDiff[3];	// Normalized a3fPointEnd - a3fPointStart
    double fLenDiff,fDiff;
    double a3fRotPoint[3];	// Rotated point.
    double a3fPlatePoint[3];	// Point projected onto plate.
    int	a2nPixel[2];

    bool bGoodPixelFound;	// We load anPoint0 and anPoint1 with possibly off the plate pixels. 
    bool bGoodCenterFound;      // Did we find the center pixel with ease.
    itr<int> anPoint[2];
    itr<int> anOnPlate;
    itr<int> anCorner;
    int a2nCenter[2];		// Center of fan.
    int nCenterOnPlate;

    double a3x3fOmegaMat[3][3];
    double a3x3fOmegaChiMat[3][3];
    double a3x3fOmegaChiPhiMat[3][3];
    double a3x3fDDInv[3][3];
    double a3x3fDD[3][3];
    double a3fDN[3];
    double a3fTemp[3];
    double a3fVTS[3];
    float a3fRots[3];
    float a3x3fTemp[3][3];
    float fPix0,fPix1;
    int nx;
    int nRot;

    int nHorizCoord, nVertCoord;
    int a2nHoriz[2];
    int a2nVert[2];     // [0] = bottom most.  [1] == top most.

    m_poDetector->vUpdateDDDN();
    (void) m_poDetector->nCalcGetDDDN(&a3x3fDD[0][0], a3fDN);   // Get detector information and normal vector.
    fInvMat3D(&a3x3fDD[0][0], &a3x3fDDInv[0][0]);
    
    // Determine horizontal direction.
    {
      double        fMaxX = -99999.0;
      float         fVDC[3];
      double        a3fTemp1[3];
      double        a3fTemp2[3];
      double        fMaxXIndex = 0.0;
      int           nMaxXIndex = -1;
      int           nIndex = 0;
      int           nPix0 = 0;
      int           nPix1 = 0;
      
      
      for (nIndex = 0; nIndex < 4; nIndex++)
	{
	  if (nIndex == 0)
	    {
	      nPix0 = m_poDetector->m_a2nDim[0]/2;
	      nPix1 = 0;
            } 
	  else if (nIndex == 1)
	    {
	      nPix0 = m_poDetector->m_a2nDim[0]/2;
	      nPix1 = m_poDetector->m_a2nDim[1]-1;
            } 
	  else if (nIndex == 2)
	    {
	      nPix0 = m_poDetector->m_a2nDim[0]-1;
	      nPix1 = m_poDetector->m_a2nDim[1]/2;
            } 
	  else 
	    {
	      nPix0 = 0;
	      nPix1 = m_poDetector->m_a2nDim[1]/2;
            }       	
	  nStat   = m_poDetector->m_poSpatial->nPixeltoMM(nPix0, nPix1, &fVDC[0], &fVDC[1], &fVDC[2]);
	  vCopyVec3D(fVDC,a3fTemp1);
	  if (!nStat)
	    {
	      vMulMat3DVec3D(a3x3fDD, a3fTemp1,a3fTemp2);  // Calc scattered beam wavevector
	      if ((nMaxXIndex == -1) || (a3fTemp2[0] > fMaxXIndex))
		{
		  fMaxXIndex = a3fTemp2[0];
		  nMaxXIndex = nIndex;
                }
            }
        }
      if (nMaxXIndex == -1)
	return 1;
      else if (nMaxXIndex == 0)
	{
	  nHorizCoord = 0;
	  nVertCoord = 1;
	  a2nHoriz[0] = 0;
	  a2nHoriz[1] = m_poDetector->m_a2nDim[0] - 1;
	  a2nVert[0] = 0;
	  a2nVert[1] = m_poDetector->m_a2nDim[1] - 1;
        } 
      else if (nMaxXIndex == 1)
	{
	  nHorizCoord = 0;
	  nVertCoord = 1;
	  a2nHoriz[0] = 0;
	  a2nHoriz[1] = m_poDetector->m_a2nDim[0] - 1;
	  a2nVert[0] = m_poDetector->m_a2nDim[1] - 1;
	  a2nVert[1] = 0;
        } 
      else if (nMaxXIndex == 2)
	{
	  nHorizCoord = 1;
	  nVertCoord = 0;
	  a2nHoriz[0] = 0;
	  a2nHoriz[1] = m_poDetector->m_a2nDim[1] - 1;
	  a2nVert[0] = m_poDetector->m_a2nDim[0] - 1;
	  a2nVert[1] = 0;
        } 
      else
	{
	  nHorizCoord = 1;
	  nVertCoord = 0;
	  a2nHoriz[0] = 0;
	  a2nHoriz[1] = m_poDetector->m_a2nDim[1] - 1;
	  a2nVert[0] = 0;
	  a2nVert[1] = m_poDetector->m_a2nDim[0] - 1;
        }
    }

    memset(m_pcMask,0,m_poDetector->m_a2nDim[0]*m_poDetector->m_a2nDim[1]);
    nCount = 1;
    
    // Loop through start and end rotation positions.
    for (nRot = 0; nRot < 2; nRot ++)
      {
        // Determine the rotation matrices.
        m_poCrysGonio->nGetDatum(3,&a3fRots[0]);
        
        nx = m_poCrysGonio->nGetNum(m_poScan->m_poRotation->sGetName());
        if ((nx<0) || (nx>=3))
	  {
            printf("ERROR:  Could not find rotation axis '%s' for goniometer shadow masking!\n",
		   m_poScan->m_poRotation->sGetName().string());
            return 1;
	  }
        if (nRot == 0) 
	  a3fRots[nx] += m_poScan->m_poRotation->fGetRotStart();
        else
	  a3fRots[nx] += m_poScan->m_poRotation->fGetRotEnd();
        
        m_poCrysGonio->vCalcGetRotMatrix(
					 &a3x3fTemp[0][0],-1,
					 a3fRots[0],0.0,0.0);
        vCopyMat3D(&a3x3fTemp[0][0],&a3x3fOmegaMat[0][0]);
        m_poCrysGonio->vCalcGetRotMatrix(
					 &a3x3fTemp[0][0],-1,
					 a3fRots[0],a3fRots[1],0.0);
        vCopyMat3D(&a3x3fTemp[0][0],&a3x3fOmegaChiMat[0][0]);
        m_poCrysGonio->vCalcGetRotMatrix(
					 &a3x3fTemp[0][0],-1,
					 a3fRots[0],a3fRots[1],a3fRots[2]);
        vCopyMat3D(&a3x3fTemp[0][0],&a3x3fOmegaChiPhiMat[0][0]);
        
        for (nBox = 0; nBox < m_anStage.size(); nBox++)
	  {
            -anPoint[0];
            -anPoint[1];
            -anOnPlate;
            -anCorner;
            bGoodPixelFound = false;
            bGoodCenterFound = false;
            for (nEdge = 0; nEdge < 5; nEdge++)
	      {
                if (nEdge < 4)
		          {
                            a3fPointStart[0] = m_a4anBoxX[nEdge][nBox]; 
                            a3fPointStart[1] = m_a4anBoxY[nEdge][nBox]; 
                            a3fPointStart[2] = m_a4anBoxZ[nEdge][nBox]; 
                            a3fPointEnd[0] = m_a4anBoxX[(nEdge + 1) % 4][nBox]; 
                            a3fPointEnd[1] = m_a4anBoxY[(nEdge + 1) % 4][nBox]; 
                            a3fPointEnd[2] = m_a4anBoxZ[(nEdge + 1) % 4][nBox]; 
                            vSubVec3DVec3D(a3fPointEnd,a3fPointStart,a3fNormDiff);
                            fLenDiff = fNormVec3D(a3fNormDiff);
		          } 
		            else 
		              {
                                fLenDiff = 0.0;
                                vZeroMat(3,1,a3fPointStart);
                                vZeroMat(3,1,a3fPointEnd);
                                for (nx=0;nx<4;nx++)
		                  {
                                    a3fPointStart[0] += m_a4anBoxX[nx][nBox]
			                               *(((rand() % 100)-50)/100.0 + 1.0);
                                    a3fPointStart[1] += m_a4anBoxY[nx][nBox]
			                               *(((rand() % 100)-50)/100.0 + 1.0);
                                    a3fPointStart[2] += m_a4anBoxZ[nx][nBox]
			                               *(((rand() % 100)-50)/100.0 + 1.0);
		                  }
                                vMulVec3DScalar(a3fPointStart,1.0/4.0,a3fPointStart);
		              }
                
          for (fDiff = 0; fDiff < fLenDiff + fMaxEdgeStep; fDiff += fMaxEdgeStep)
		  {
              vMulVec3DScalar(a3fNormDiff,min(fLenDiff,fDiff),a3fTemp);
              vAddVec3DVec3D(a3fPointStart,a3fTemp,a3fPoint);
                    
              // Rotate the points according to the matrix that is needed.
              if (m_anStage[nBox]==m_nOmegaStage) 
                  vMulMat3DVec3D(a3x3fOmegaMat,a3fPoint,a3fRotPoint);
              else if (m_anStage[nBox]==m_nChiStage)
		          vMulMat3DVec3D(a3x3fOmegaChiMat,a3fPoint,a3fRotPoint);
              else if (m_anStage[nBox]==m_nPhiStage)
		          vMulMat3DVec3D(a3x3fOmegaChiPhiMat,a3fPoint,a3fRotPoint);
                    
            // We have the rotated coordinates.  
            // Figure out where the rotated coordinate hits the detector plate.
            fNormVec3D(a3fRotPoint);
            
            vCopyVec3D(a3fRotPoint,a3fTemp);
            
            bool        bScaleSOnPlateSuccess = (0 == m_poDetector->nScaleSOnPlate(a3fRotPoint, a3fVTS));
            
            if( bScaleSOnPlateSuccess )
            {
                vMulMat3DVec3D(a3x3fDDInv, a3fVTS, a3fPlatePoint);
            
                nStat = 0;
            }
            else
                nStat = 1;

            if( 0 == nStat && (fDot3D(a3fRotPoint,a3fVTS) <= 0.0 || fDot3D(a3fTemp,a3fRotPoint) < 0.0) )
		    {
                nStat = 1;
		    }
                    
            // Find the pixel coordinates.                    
            // Don't let the Cdetector routines clip the coordinates.
            if( 0 == nStat && (0 != m_poDetector->m_poSpatial->nMMtoPixel(a3fPlatePoint[0],
						                                                  a3fPlatePoint[1],
						                                                  a3fPlatePoint[2],
						                                                 &fPix0,
                                                                         &fPix1,
						                                                  FALSE)) )
		    {
                nStat = 1;
		    }
                    
            if( nStat )
		    {
                anOnPlate + 1; // RB 2/3/2006 This looks suspicious. Maybe Thad meant "anOnPlate + 0" ?                   
		    } 
		    else 
		    {
                if (nEdge == 4)
			    {
                          bGoodCenterFound = TRUE;
                          a2nCenter[0] = (int) fPix0;
                          a2nCenter[1] = (int) fPix1;
                          nCenterOnPlate = !
			              ((a2nCenter[0]<0) || (a2nCenter[0]>=m_poDetector->m_a2nDim[0]) ||
			              (a2nCenter[1]<0) || (a2nCenter[1]>=m_poDetector->m_a2nDim[1]));
                      
                } 
			    else 
                {
                          a2nPixel[0] = (int) fPix0;
                          a2nPixel[1] = (int) fPix1;
                      
                          // Use pixels even if they are off the plate.
                          if (   (a2nPixel[0] < 0) 
				                   || (a2nPixel[0] >= m_poDetector->m_a2nDim[0])
				                   || (a2nPixel[1] < 0) 
                                 || (a2nPixel[1] >= m_poDetector->m_a2nDim[1]))
                          {
                              anOnPlate + 0;
                          }
			              else
                          {
                              bGoodPixelFound = TRUE;
                              anOnPlate + 1;
                          }
                }
		    }
                    
               // Add this point to the array.

               anPoint[0] + a2nPixel[0];
               anPoint[1] + a2nPixel[1];
		  }
          
          if( anOnPlate.size() )
		  {
                    while (anCorner.size() < anOnPlate.size())
		                anCorner + 0;
                    
                    anCorner.last() = 1;
		  }
	      }
            if ((anPoint[0].size()>=3) && (bGoodPixelFound))
	      {
                nPoints = anPoint[0].size();
                for (nStart = 0; nStart < nPoints;nStart++)
		  {
                    nEnd = (nStart + 1) % nPoints;
                    if ((!anOnPlate[nStart]) && (!anOnPlate[nEnd]))
		      continue;

                    int nPointHorizStart,nPointHorizEnd,nPointHoriz;
                    int nPointVert;
		    int nStartVert;
                    nPointHorizStart = min(anPoint[nHorizCoord][nStart],
					   anPoint[nHorizCoord][nEnd]);
                    nPointHorizEnd = max(anPoint[nHorizCoord][nStart],
					 anPoint[nHorizCoord][nEnd]);
                    for (nPointHoriz = nPointHorizStart; 
			 nPointHoriz <= nPointHorizEnd; nPointHoriz++)
		      {
                        if ( (nPointHoriz < m_poDetector->m_a2nDim[nHorizCoord])
			    && (nPointHoriz >=0))
			  {
			    nStartVert = (int) ((anPoint[nVertCoord][nStart]*((double) nPointHorizEnd - nPointHoriz)
						 + anPoint[nVertCoord][nEnd]*((double) nPointHoriz - nPointHorizStart))/(nPointHorizEnd - nPointHorizStart ));
			    if (a2nVert[0]>a2nVert[1])
			      {
				if (nStartVert > a2nVert[0])
				  continue;
				else if (nStartVert<a2nVert[1])
				  nStartVert = a2nVert[1];
			      } 
			    else 
			      { 
				// (a2nVert[0]<a2nVert[1])
				if (nStartVert < a2nVert[0])
				  continue;
				else if (nStartVert > a2nVert[1])
				  nStartVert = a2nVert[1];
			      }

                            if (a2nVert[0]<nStartVert)
			      {
                                for (nPointVert = nStartVert; 
				     nPointVert >= a2nVert[0]; nPointVert--)
				  {
				    if (nVertCoord)
				      {
					nx = m_pcMask[m_poDetector->m_a2nDim[0]*nPointVert + nPointHoriz];
					if (nx)
					  break;
					else 
					  m_pcMask[m_poDetector->m_a2nDim[0]*nPointVert + nPointHoriz] = 1;
				      } 
				    else 
				      {
					nx = m_pcMask[m_poDetector->m_a2nDim[0]*nPointHoriz + nPointVert];
					if (nx)
					  break;
					else
					  m_pcMask[m_poDetector->m_a2nDim[0]*nPointHoriz + nPointVert] = 1;
				      }
				  }
			      } 
			    else 
			      {
                                for (nPointVert = nStartVert; 
				     nPointVert <= a2nVert[0]; nPointVert++)
				  {
				    if (nVertCoord)
				      {
					nx = m_pcMask[m_poDetector->m_a2nDim[0]*nPointVert + nPointHoriz];
					if (nx)
					  break;
					else 
					  m_pcMask[m_poDetector->m_a2nDim[0]*nPointVert + nPointHoriz] = 1;
				      } 
				    else 
				      {
					nx = m_pcMask[m_poDetector->m_a2nDim[0]*nPointHoriz + nPointVert];
					if (nx)
					  break;
					else
					  m_pcMask[m_poDetector->m_a2nDim[0]*nPointHoriz + nPointVert] = 1;
				      }
				  }
			      }
			  }
                        
		      }
		  }
	      }
	  }
      }
    return 0;
}

int 
CGonioMask::nInitShadows(Cimage *poImage)
{
    // This routine and the following one (nMaskShadows(...)) implements
    // a straightforward algorithm to mask out shadows on images.  Shadows
    // are specified by header keywords of the form CRYSTAL_GONIO_SHADOW_# where 
    // # is an integer starting at 1.  All mask keywords are checked to see if they
    // are present, so shadows do not need to be numbered consecutively nor
    // contiguously.  Each shadow is defined by a CGonioShadow object, defined
    // in the CGonioShadow.h header.

    // One can define circles and polygon shadows and store them
    // in the image header or in environment variables.  A first shadow is defined
    // here at the present date 15-Apr-2004.
    // RAPID users: A set of shadows will be defined for you.
    // NOTE:  Shadows are applied in the Cimage routine during integration
    //        and depending on the DTREK_IMAGE_APPLYMASK environment variable.
    // For example, shadows are not normally applied in dtdisplay, so you cannot
    // see them there ('setenv DTREK_IMAGE_APPLYMASK anything' to apply them).
  
    int       nStat = nReadGonio(poImage->m_oHeader);
    if (0 == nStat) 
    {
        m_dDetDistance = m_poDetector->fGetDistance();
        m_dDetSwing    = m_poDetector->fGetSwing();
    }

    Cimage_header*    poHeader = &poImage->m_oHeader;

    Cstring           sTemp("CRYSTAL_GONIO_NUM_SHADOWS");
    int               nNumShadows = 0;
    nStat = poHeader->nGetValue(sTemp, &nNumShadows);

    bool            bImagePixelSizeAvailable = false;
    double          daPixSize[2] = {0.0, 0.0};
    if( 0 < nNumShadows )
    {
        bImagePixelSizeAvailable = bGetImagePixelSize(poHeader, daPixSize);

        Cstring       strShadowData("");
        for (int nx = 0; nx < nNumShadows; nx++)
        {
            // Look for header keywords (and possible environment variables)
            sTemp = "CRYSTAL_GONIO_SHADOW_" + Cstring(nx+1);
            if( 0 != poHeader->nGetValue(sTemp, &strShadowData) )
            {
                // In case the keyword was not read properly
	            cout << "WARNING expected keyword " << sTemp << " has problems.\n";
            }
            else
            {
                // Create a shadow object and stuff it into the vector
                m_vpShadows.push_back(new CGonioShadow(strShadowData, daPixSize));
            }
        }
    }

    return 0;
}

int 
CGonioMask::nMaskShadows(Cimage* poImage)
{
    //cout << "nMaskShadows ...\n" << flush;
    
    if( !poImage ) 
      return 0;

    Cimage_header*  poHeader = &poImage->m_oHeader;
    double          a3fGonioValues[3];
    //
    int             nStat = poHeader->nGetValue(Cstring(D_K_CrystalPrefix)+D_K_GonioValues, 3, &a3fGonioValues[0]);

    if( 0 != nStat ) 
      return nStat;

    // GET PHI
    Cgoniometer     oGon(*poHeader, D_K_CrystalPrefix);

    int             nPhiValueIndex = oGon.nGetNum("phi");
    double          dImagePhiValue = a3fGonioValues[nPhiValueIndex];

    // GET CHI (KAPPA)
    int             nChiValueIndex = oGon.nGetNum("chi");
    double          dImageChiKappaValue = a3fGonioValues[nChiValueIndex];

    // GET OMEGA
    double          daRotationRange[2] = {0.0, 0.0};
    nStat = poHeader->nGetValue(Cstring(D_K_ScanPrefix)+D_K_Rotation, 2, daRotationRange);
    if(0 != nStat)
      return (nStat);
    //////////////////////////////////////////////////////////////////////////////////////////
    
    CGonioShadow*       poShadow = NULL;
    double      dRotationStepForShadowPurposes = 1.0;  // trial
    double      dCurrentRotationValue = 0.0;
    for(int nx = 0; nx < m_vpShadows.size(); nx++)
    {
        poShadow = m_vpShadows[nx];

        if( !poShadow )
          continue;

        if( !poShadow->bIsAvailable() )
          continue;

        if( !poShadow->bIsHardwareLimitsApply(dImageChiKappaValue, DT_GH_HARDWARE_LIMITS_CHI_KAPPA) ||
            !poShadow->bIsHardwareLimitsApply(dImagePhiValue, DT_GH_HARDWARE_LIMITS_PHI) ||
            !poShadow->bIsHardwareLimitsApply(m_dDetDistance, DT_GH_HARDWARE_LIMITS_DET_DISTANCE) ||
            !poShadow->bIsHardwareLimitsApply(m_dDetSwing, DT_GH_HARDWARE_LIMITS_DET_SWING) )
            continue;
        
        dCurrentRotationValue = daRotationRange[0];
        while( dCurrentRotationValue <= daRotationRange[1] )
        {
            if( poShadow->bIsHardwareLimitsApply(dCurrentRotationValue, DT_GH_HARDWARE_LIMITS_OMEGA) )
            {
                // Satifies Phi, ChiKappa and Omega criteria, so act on this shadow
                // Adjust pixel coordinates to account for 
                // this image's omega value.  
                // Do this for all polygon points
                // and for the polygon limits of the circle.
                poShadow->bRotate(dCurrentRotationValue, poImage->nGetDimension(0));
                poShadow->bErasePixels(poImage);
            }

            dCurrentRotationValue += dRotationStepForShadowPurposes;
        }
    }
    
    return (0);
}

static bool bGetImagePixelSize(Cimage_header* poHeader, double daPixSize[2])
{
    if( !poHeader )
        return false;

    Cstring             strDetectorPrefix(""); 
    if ( 0 != poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;

    float      afDims[2] = {0.0, 0.0}; 
    if( 0 != poHeader->nGetValue(strDetectorPrefix+D_K_DetectorDimensions, 2, afDims) )
        return false;

    float      afSizes[2] = {0.0, 0.0}; 
    if( 0 != poHeader->nGetValue(strDetectorPrefix+D_K_DetectorSize, 2, afSizes) )
        return false;

    if( afDims[0] == 0.0 || afDims[1] == 0.0 || afSizes[0] == 0.0 ||afDims[1] == 0.0 )
        return false; // safety

    daPixSize[0] = (double)(afSizes[0] / afDims[0]);
    daPixSize[1] = (double)(afSizes[1] / afDims[1]);
    
    return true;
}




