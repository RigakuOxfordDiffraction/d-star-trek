//
// Copyright (c) 2000 Molecular Structure Corporation
//
// dtprofit.cc     Initial author: T.J. Niemeyer         Jul-2003
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

#include "CExperimentalStrategy.h"
#include "CCommandLine.h"
#include "CCellReduce.h"

double CExperimentalStrategy::ms_fDefaultSpotSize[3] = {20.0, 15.0, 10.0}; 

char* CExperimentalStrategy::ms_pcHelp = 
       " -reso fResoMin fResoMax\n"
       "                    Specifies the resolution range for completeness\n"
       "                    calculations.  The default resolution range is to\n"
       "                    use '-reso edges'.  If -predicted is used,\n"
       "                    the resolution in the input file determines the\n"
       "                    upper and lower bounds.\n\n"
       " -reso edges\n"
       "                    Set upper resolution to be equal to the highest\n"
       "                    *minimum* resolution found on any of the four\n"
       "                    edges (sides) of the image.  If the calculated beam\n"
       "                    center is not located on the image, then also compute\n"
       "                    the lower resolution to be the lowest *maximum* \n"
       "                    resolution found on any of the four edges.\n\n"
       " -reso corners\n"
       "                    Set upper resolution to be equal to the highest\n"
       "                    *maximum* resolution found on any of the four\n"
       "                    edges (sides) of the image.  If the calculated beam\n"
       "                    center is not located on the image, then also compute\n"
       "                    the lower resolution to be the lowest *minimum* \n"
       "                    resolution found on any of the four edges.\n\n"
       " -reso upperEqualsHighest\n"
       "                    Synonym for '-reso corners'.\n\n"
       " -reso minarc fMinArcLength\n"
       "                    Set the upper resolution to be equal to the highest\n"
       "                    value for which diffracted Debye rings have at least\n"
       "                    fMinArcLength degrees of data.  If the calculated beam\n"
       "                    center is not located on the image, then also compute\n"
       "                    the lower resolution range to be the lowest value for\n"
       "                    which diffracted Debye rings have at least fMinArcLength\n"
       "                    degrees of data.\n"
       "                    WARNING:  This option will fail if requested conditions\n"
       "                              cannot be found on the detector.\n\n"
       " -reso 2theta fTargetReso f2ThetaMin f2ThetaMax\n"
       "                    This option will search for a 2Theta swing between\n"
       "                    f2ThetaMin and f2ThetaMax that obtains an upper\n"
       "                    resolution of fTargetReso.  The user can still control how\n"
       "                    the resolution is computed by specifying a separate -reso\n"
       "                    command ('-reso edges', '-reso corners', '-reso upperEqualsHighest',\n" 
       "                    or '-reso minarc') anywhere else on the command line.\n"
       "                    WARNING:  This option will fail if requested conditions\n"
       "                              cannot be found on the detector.\n\n"
       " -reso distance fTargetReso fDistanceMin fDistanceMax\n"
       "                    This option will search for a detector distance between\n"
       "                    fDistanceMin and fDistanceMax that obtains an upper\n"
       "                    resolution of fTargetReso.  The user can still control how\n"
       "                    the resolution is computed by specifying a seperate -reso\n"
       "                    command ('-reso edges', '-reso corners', '-reso upperEqualsHighest',\n" 
       "                    or '-reso minarc') anywhere else on the command line.\n"
       "                    anywhere else on the command line.\n"
       "                    WARNING:  This option will fail if requested conditions\n"
       "                              cannot be found on the detector.\n\n"
       " -reso printonly\n"
       "                    Only print resolution info, then terminate dtstrategy.\n\n"
       " -distance fDistance\n"
       "                    Set the crys-2-detector distance.\n\n"
       " -swing2theta f2Theta\n"
       "                    Set the 2Theta swing angle\n\n";

CExperimentalStrategy::CExperimentalStrategy(Cimage_header& oHeader):m_oHeader(oHeader),m_oDetector(oHeader),m_oSource(oHeader)
{
    m_sOutputGraph = "";
    m_bPrint = false;
    m_nMethod = eResoInner;
    m_fDebyeRingCoverage = 90.0;
};

bool CExperimentalStrategy::bIsAvailable() 
{
    return ((m_oDetector.bIsAvailable()) && (m_oSource.bIsAvailable()));
}

int CExperimentalStrategy::nBuildResoShells(double fResoLow,double fResoHigh,int nNumShells,double* pfShellStart) {
    double fSlopeReso;
    double fLow;
    int nx;

    if (fResoLow < fResoHigh)
        std::swap(fResoLow,fResoHigh);

    fSlopeReso = (1.0/(fResoHigh * fResoHigh * fResoHigh) - 1.0/(fResoLow * fResoLow * fResoLow))/nNumShells;
    fLow = fSlopeReso;
    for (nx = 0; nx < nNumShells; nx++) 
    {
        pfShellStart[nx] = 1.0 / max(1e-10, pow(fLow, 0.333333));
        fLow += fSlopeReso;
    };
    pfShellStart[nNumShells] = fLow;
    return 0;
};

double CExperimentalStrategy::fConvertToDegrees(double fReso) {
    // Change to degrees.
    // lambda = 2dsin(theta)
    // So, 2theta = asin(lambda/(2*d))*2.0
    double fLambda = m_oSource.m_poWavelength->fGetWavelength();
    return asin(fLambda/(2.0*fReso))*2.0/Gs_dRADIANS_PER_DEGREE;
    
};

int CExperimentalStrategy::nCalcReso(double& fResoHigh,double& fResoLow,double f2Theta,double fDistance)
{
    float a3fS0[3];
    float f0,f1;
    double a4fTemp[4];
    double fWavelength;
    int nx,ny;
    int nStat;
    double f2ThetaOrig;
    double fDistOrig;
    itr<double> afChiCoverage;
    itr<double> afResolution;
    float fResoMin,fResoMax,fResoMaxEdge;

    
    f2ThetaOrig = m_oDetector.m_poGoniometer->fGetSwing();
    fDistOrig = m_oDetector.m_poGoniometer->fGetDistance();
    if (f2Theta > -999.0)
        m_oDetector.m_poGoniometer->nSetSwing(f2Theta);
    if (fDistance > -999.0)
        m_oDetector.m_poGoniometer->nSetDistance(fDistance);
        
    
    m_oSource.vCalcGetS0(a4fTemp);            
    vCopyVec3D(a4fTemp,a3fS0);    
    fWavelength = m_oSource.fGetWavelength();  
    vMulVec3DScalar(a3fS0,1.0/fWavelength,a3fS0);
    ny = min(m_oDetector.m_a2nDim[0],m_oDetector.m_a2nDim[1]);
    nStat = 0;
    
    if (m_oDetector.nGetResolution(a3fS0,&fResoMin,&fResoMax,&fResoMaxEdge))
        nStat = 2;
    else {
        
        if (m_nMethod == eResoOuter) {
            fResoLow = fResoMin;
            fResoHigh = fResoMax;
        } else if (m_nMethod == eResoInner) {
            fResoLow = fResoMin;
            fResoHigh = fResoMaxEdge;
        } else if (m_nMethod == eResoChiCoverage) {
            int nBinsReso = 50;
            
            afChiCoverage.setsize(nBinsReso);
            afResolution.setsize(nBinsReso+1);
            
            fResoHigh = fResoLow = fResoMin;
            fResoMin = min(100.0,fResoMin);
            nBuildResoShells(fResoMax,fResoMin,nBinsReso,&afResolution[0]);
            afResolution.setsize(afChiCoverage.size());
            
            if( afChiCoverage.size() == 0 ||       // just a safety check in case the code is changed later
                m_oDetector.nGetResolution(a3fS0,
                                           0,
                                           ny/50,
                                           0.0,
                                           &f0,
                                           &f1,
                                           afChiCoverage.size(),
                                           &afResolution[0],
                                           &afChiCoverage[0]) )
            {
                printf("ERROR:  Could not determine resolution bounds!\n");
                nStat = 2;
            }
            else
            {
                // Now, starting at the lower resolution, continue until we find the place were we no longer
                // have the chi coverage we need.
                bool bResoFound = false;

                for (nx = 0; nx < afResolution.size(); nx++) {
                    if (afChiCoverage[nx] >= min(350.0,m_fDebyeRingCoverage)) {
                        fResoHigh =  afResolution[nx];
                        bResoFound = true;
                    } else if (bResoFound)
                        break;
                };
            };
            
            //nPlot("resoplot.txt",afChiCoverage,&afResolution);
        };
    };

  
    m_oDetector.m_poGoniometer->nSetSwing(f2ThetaOrig);
    m_oDetector.m_poGoniometer->nSetDistance(fDistOrig);

    return nStat;

};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CExperimentalStrategy::fCalcBest2ThetaSwingForTargetHighReso(double f2ThetaMin, 
                                                                    double f2ThetaMax, 
                                                                    double fTargetHighReso, 
                                                                    double& fBestPossibleHighReso)
{
    // Check if the problem could be solved in principle, i.e. if the high reso target can be achieved at the highest detector position
    double      fResoLow  = -1.0;
    double      fResoHigh = -1.0;
    
    fBestPossibleHighReso = -1.0; // initialize
    
    if( 0 != nCalcReso(fResoHigh, fResoLow, f2ThetaMax) ) 
    {
        return -999.0;
    }
    
    fBestPossibleHighReso = fResoHigh;
    if( fTargetHighReso < fResoHigh )
    {
        return -999.0;
    }
    
    // Else find the best 2Theta swing iteratively
    double      f2ThetaStep = max(1.0, floor((f2ThetaMax - f2ThetaMin)/20.0)); 

    bool        b2ThetaWorks = false;
    double      f2ThetaResult = 0.0;
    
    double      f2Theta = 0.0;
    
    for(f2Theta = f2ThetaMin; f2ThetaMax - f2Theta > f2ThetaStep / 10.0; f2Theta += f2ThetaStep)
    {
        if( 0 != nCalcReso(fResoHigh, fResoLow, f2Theta) ) 
            continue;
        
        if( fTargetHighReso >= fResoHigh && fTargetHighReso <= fResoLow ) 
        {
            f2ThetaResult = f2Theta;
            b2ThetaWorks = true;
            
            fBestPossibleHighReso = fResoHigh;

            break;
        }
    }
    
    if (b2ThetaWorks)
        return  f2ThetaResult;
    
    return -999.0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CExperimentalStrategy::fCalcBest2ThetaSwingForTargetLowReso(double f2ThetaMin, 
                                                                   double f2ThetaMax, 
                                                                   double fTargetLowReso,
                                                                   double& fBestPossibleLowReso)
{
    if( fTargetLowReso > DTREK_DET_RESO_MAX )
    {
        fTargetLowReso = DTREK_DET_RESO_MAX;
    }

    // Check if the problem could be solved in principle, i.e. if the low reso target can be achieved at the lowest detector position
    double      fResoLow  = -1.0;
    double      fResoHigh = -1.0;
    
    fBestPossibleLowReso = -1.0; // initialize
    
    if( 0 != nCalcReso(fResoHigh, fResoLow, f2ThetaMin) ) 
    {
        return -999.0;
    }
    
    fBestPossibleLowReso = fResoLow;
    if( fTargetLowReso > fResoLow )
    {
        return -999.0;
    }

    // Else find the best 2Theta swing iteratively
    double      f2ThetaStep = max(1.0, floor((f2ThetaMax - f2ThetaMin)/20.0)); 

    bool        b2ThetaWorks = false;
    double      f2ThetaResult = 0.0;
    
    double      f2Theta = 0.0;
    
    for(f2Theta = f2ThetaMax; f2Theta - f2ThetaMin > f2ThetaStep / 10.0; f2Theta -= f2ThetaStep)
    {
        if( 0 != nCalcReso(fResoHigh, fResoLow, f2Theta) ) 
            continue;
        
        if( fTargetLowReso >= fResoHigh && fTargetLowReso <= fResoLow ) 
        {
            f2ThetaResult = f2Theta;
            b2ThetaWorks = true;
            
            fBestPossibleLowReso = fResoLow;
            
            break;
        }
    }
    
    if (b2ThetaWorks)
        return  f2ThetaResult;
    
    return -999.0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CExperimentalStrategy::fCalcBestDistance(double fDistanceMin,double fDistanceMax,double fTargetReso)
{
    double      fResoLow = 0.0 , fResoHigh = 0.0;

    double          fDistanceStep = max(1.0,floor((fDistanceMax - fDistanceMin)/20.0));
    double          fDistanceResult = 0.0;
    bool            bDistanceWorks = false;
    
    for (double fDistance = fDistanceMax; fDistance - fDistanceMin > fDistanceStep / 10.0; fDistance -= fDistanceStep)
    {
        if( nCalcReso(fResoHigh,fResoLow,-999.0,fDistance) ) 
            continue;

        if( fTargetReso >= fResoHigh && fTargetReso <= fResoLow )
        {
            fDistanceResult = fDistance;
            bDistanceWorks = true;
            break;
        }
    }
    
    if( bDistanceWorks )
        return fDistanceResult;
    
    return -999.0;
}

double CExperimentalStrategy::fCalcBestDistanceSpotSeparation(double fPixelSep, 
                                                              int nModeAxis, 
                                                              int nModeSpot,
                                                              SPOT_SIZE_INFO* pstReturnSpotSizeInfo)
{
    int         nx = 0, ny = 0;
    
    double      a3fAxis[3];
    
    double      a3x3fDist[3][3];  // [0] = Spot [1] = Crystal
    
    char        acTemp[3][200];
    
    double      fMMSep;
    double      fLambda;
    double      f2ThetaRad;
    double      fDistance;
    
    Ccrystal    oCrystal(m_oHeader);
    if( !oCrystal.bIsAvailable() ) 
        return -999.0;
    
    /////////////////////////////////////////////////////////////////////////
    // REDUCE CELL
    double          a6fCell[6];
    oCrystal.vGetCell(a6fCell);
    printf("\nInput cell  : a=%.2f b=%.2f c=%.2f alpha=%.2f beta=%.2f gamma=%.2f\n", a6fCell[0],a6fCell[1],a6fCell[2],a6fCell[3],a6fCell[4],a6fCell[5]);

    CCellReduce     oReduce;
    oReduce.nReduceCell(a6fCell);

    oCrystal.vSetCell(a6fCell);
    printf("Reduced cell: a=%.2f b=%.2f c=%.2f alpha=%.2f beta=%.2f gamma=%.2f\n", a6fCell[0],a6fCell[1],a6fCell[2],a6fCell[3],a6fCell[4],a6fCell[5]);
    ///////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Determine spot size
    double      a3fSpotDim[3]={ms_fDefaultSpotSize[0], ms_fDefaultSpotSize[1], ms_fDefaultSpotSize[2]}; // Initialize max, average and min spot sizes
    
    SPOT_SIZE_INFO  stSpotSizeInfo={0.0};
    
    if( bGetSpotSizeFromHeader(stSpotSizeInfo) )
        stSpotSizeInfo.vToArray(a3fSpotDim);
    else
    {
        printf("\nWarning: the header does not contain spot size information.\nAssuming default spot size.\n");
    }

    if( pstReturnSpotSizeInfo )
        pstReturnSpotSizeInfo->vFromArray(a3fSpotDim);
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    if( nModeAxis == eComputeDistanceSpotSeparationDefault) 
        nModeAxis = eComputeDistanceSpotSeparationWorstCase;
    
    if (nModeSpot == eComputeDistanceSpotSeparationDefault) 
        nModeSpot = eComputeDistanceSpotSeparationAverageCase;

    a3fAxis[0] = max(oCrystal.fGetCell(0), max(oCrystal.fGetCell(1), oCrystal.fGetCell(2)) );

    a3fAxis[1] = (oCrystal.fGetCell(0) + oCrystal.fGetCell(1) + oCrystal.fGetCell(2)) / 3.0;
    
    a3fAxis[2] = min(oCrystal.fGetCell(0),min(oCrystal.fGetCell(1),oCrystal.fGetCell(2)));

    fLambda = m_oSource.m_poWavelength->fGetWavelength();
   
    // Calculate the ideal distance. 
    for (nx = 0; nx < 3; nx++) 
    {
        for (ny = 0; ny < 3; ny++) 
        {
            // Using lambda = 2d sin(theta), we get 2theta = asin(lambda/(2*d))*2.
        
            // So, the MM shift on the plate is Shift = Distance * sin(2theta);

            // Calculate the MM separation from the pixel separation.
            fMMSep = max(m_oDetector.m_poSpatial->fGetNominalPixelSize(0)*( fPixelSep + 1.2 * a3fSpotDim[nx] ),
                         m_oDetector.m_poSpatial->fGetNominalPixelSize(1)*( fPixelSep + 1.2 * a3fSpotDim[nx] )  );
            
            f2ThetaRad = asin(fLambda/( 2.0 * a3fAxis[ny])) * 2.0;   // RB: Thad is actually calculating the 2-Theta for a lowest order reflection.
                                                                     // So he is figuring out the distance of a lowest-order reflection from the center.
            
            fDistance = fMMSep / sin(f2ThetaRad);
            
            a3x3fDist[nx][ny] = fDistance;
        }
    }

    sprintf(&acTemp[0][0],"(%.0f A)",a3fAxis[0]);
    sprintf(&acTemp[1][0],"(%.0f A)",a3fAxis[1]);
    sprintf(&acTemp[2][0],"(%.0f A)",a3fAxis[2]);

    printf("\nDetector distance for %d pixel spot separation\n",(int) (fPixelSep));
    printf("-----------------------------------------------------\n");
    printf("                 Long Axis  Average Axis  Short Axis\n");
    printf(" Spot Diameter   %9s  %12s  %10s\n",acTemp[0],acTemp[1],acTemp[2]);
    printf("-----------------------------------------------------\n");
    sprintf(&acTemp[0][0],"Large (%3.0f Px)", a3fSpotDim[0]);
    sprintf(&acTemp[1][0],"Avg   (%3.0f Px)", a3fSpotDim[1]);
    sprintf(&acTemp[2][0],"Small (%3.0f Px)", a3fSpotDim[2]);

    for (nx = 0; nx < 3; nx++) 
    {
        printf(" %s", &acTemp[nx][0]);
        printf("  %9.1f  %12.1f  %10.1f\n", a3x3fDist[nx][0], a3x3fDist[nx][1], a3x3fDist[nx][2]);
    }

    printf("-----------------------------------------------------\n");
    fDistance = a3x3fDist[nModeSpot][nModeAxis];
    
    printf("Suggested distance: %.1f mm\n", fDistance);

    return fDistance;
}

int CExperimentalStrategy::nCalcBest2ThetaSwings(double fInitSwing, double fDegreesOverlaping, itr<double>& af2Thetas) 
{
    itr<double>         a2f2Theta[2];
    itr<double>         af2ThetaSwing;
    itr<double>         a2fReso[2];
    
    char pcBuf1[30];
    char pcBuf2[30];

    double f0,f1,f2;
    double fShift,fLastShift;
    int nx,ny;

    if (fInitSwing<=-999.0)
        fInitSwing = m_oDetector.m_poGoniometer->fGetSwing();

    af2ThetaSwing + fInitSwing;
    
    nCalcReso(f0, f1, af2ThetaSwing.last());
    
    a2fReso[0] + f0;
    a2fReso[1] + f1;
    
    a2f2Theta[0] + fConvertToDegrees(f0);
    a2f2Theta[1] + fConvertToDegrees(f1);

    // Now compute the swings.
    while( a2fReso[1].last() < 100.0 )
    {
        fShift = 0.0;
        do 
        {
            fLastShift = fShift;
            
            if( af2ThetaSwing.last() < 0.0)
                fShift+=1.0;
            else
                fShift-=1.0;
            
            if( nCalcReso(f0, f1, af2ThetaSwing.last()+fShift) )
                return 2;
            
            double ff0 = fConvertToDegrees(f0);
            double ff1 = fConvertToDegrees(f1);

            if ((af2ThetaSwing.last()<=0.0) && (af2ThetaSwing.last()+fShift>=0.0))
            {
                fLastShift = -af2ThetaSwing.last();
                break;
            }
            
            if ((af2ThetaSwing.last()>=0.0) && (af2ThetaSwing.last()+fShift<=0.0))
            {
                fLastShift = -af2ThetaSwing.last();
                break;
            }
            
            f2 = min(ABS(fConvertToDegrees(f0)-a2f2Theta[1].last()),ABS(fConvertToDegrees(f1) - a2f2Theta[0].last()));
        
        } while( f2 > fDegreesOverlaping);



        if (fLastShift!=0.0)
            af2ThetaSwing + (af2ThetaSwing.last() + fLastShift);
        else
            af2ThetaSwing + (af2ThetaSwing.last() + fShift);
        
        if (nCalcReso(f0,f1,af2ThetaSwing.last()))
            return 2;
           
        a2fReso[0] + f0;
        a2fReso[1] + f1;
        a2f2Theta[0] + fConvertToDegrees(f0);
        a2f2Theta[1] + fConvertToDegrees(f1);
    }

    printf("\nMultiple 2theta swings needed for full coverage\n");
    printf("------------------------------------------------------------\n");
    printf("    2Theta  Resolution Range  Resolution Range  Scan Overlap\n");
    printf(" #   (deg)        (Angstrom)             (deg)         (deg)\n");
    printf("------------------------------------------------------------\n");
    
    for (nx = 0; nx < af2ThetaSwing.size(); nx++)
    {
        sprintf(pcBuf1,"%.2lf - %.2lf",min(999.99,a2fReso[1][nx]),a2fReso[0][nx]);
        sprintf(pcBuf2,"%5.1lf - %5.1lf",a2f2Theta[1][nx],a2f2Theta[0][nx]);
        printf("%2d  %6.1lf  %16s  %16s  ",
            nx + 1,af2ThetaSwing[nx],pcBuf1,pcBuf2);

        itr<int> anOverlaps;
        itr<double> afOverlaps;
            
        // Now, check for overlaps between the data points.
        for (ny = 0; ny < af2ThetaSwing.size(); ny++)
        {
            if ((ny!=nx) && ((ABS(a2f2Theta[1][nx] - a2f2Theta[0][nx]) + ABS(a2f2Theta[1][ny] - a2f2Theta[0][ny]))>(max(max(a2f2Theta[1][ny],a2f2Theta[0][ny]),max(a2f2Theta[1][nx],a2f2Theta[0][nx])) - max(min(a2f2Theta[1][ny],a2f2Theta[0][ny]),min(a2f2Theta[1][nx],a2f2Theta[0][nx]))))) {
                f0 = min(ABS(a2f2Theta[1][nx] - a2f2Theta[0][ny]),ABS(a2f2Theta[0][nx] - a2f2Theta[1][ny]));
                afOverlaps + f0;
                anOverlaps + ny;
            }
        }
        
        if (afOverlaps.size() > 0)
        {
            -g_apvSwapPointers + ((void*) &afOverlaps[0]) + ((void*) & anOverlaps[0]);
            -g_anSwapSizes + sizeof(double) + sizeof(int);
            
            qsortswap((void*) &afOverlaps[0],afOverlaps.size(),sizeof(double),double_cmp,qsort_swap_arrays);
            
            for (ny = afOverlaps.size()-1; ny >= max(0,afOverlaps.size()-1-2); ny--)
            {            
                if (ny < afOverlaps.size()-1)
                    printf("\n%48s","");
                printf("With %d: %4.1lf",anOverlaps[ny]+1,afOverlaps[ny]);
            }
        }
        
        printf("\n");
    }

    printf("------------------------------------------------------------\n");

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CExperimentalStrategy::nParseResoCommandLine(int argc,char** argv, 
                                                 double& fResoHigh, double& fResoLow,
                                                 double& f2ThetaSwing, double& fDistance) 
{
    int         nArg = -2;
    
    double      f0 = 0.0;
    double      f1 = 0.0;

    const char* apcOptions[] = 
    {
        "-reso","n.outer",
        "-reso","n.corners",
        "-reso","n.upperEqualsHighest",
        "-reso","n.inner",
        "-reso","n.edges",
        "-reso","f.minarc",
        "-reso","f3.2theta",
        "-reso","f2.2theta",
        "-reso","f1.2theta",
        "-reso","f3.distance",
        "-reso","f2.distance",
        "-reso","f1.distance",
        "-reso","*n.distance","*n.worstcase","*n.bestcase","*n.averagecase","*f1.spotseparation",
        "-reso","f","f",
        "-reso","n.print",
        "-reso","n.printonly",
        "-reso","f1.multi2theta",
        "-reso","n.multi2theta",
        "-swing2theta","f",
        "-distance","f",
        NULL
    };

    int         nStat = 0;
    double      fTargetResoDist;
    double      a2fRangeDist[2] = {0.0};
    double      fTargetReso2Theta = 0.0;
    double      a2fRange2Theta[2] = {0.0};
    bool        bComputeReso = false;
    bool        bCompute2Theta = false;
    bool        bComputeDistance = false;
    bool        bComputeResoOnly = false;
    bool        bComputeDistanceSpotSeparation = false;
    int         nComputeDistanceSpotSeparation = 0;
    int         nComputeDistanceSpotSeparationModeAxis = 0;
    int         nComputeDistanceSpotSeparationModeSpot = 0;
    bool        bCompute2ThetaMultipleSwing = false;
    double      fCompute2ThetaMultipleSwingOverlap;


    g_pcHelp = ms_pcHelp;
    
    //nPrintItemSyntax(apcOptions, Cstring("-reso"));

    while( 0 == (nStat = nReadCommandLine(nArg, argc, argv, apcOptions, 0)) )  // the last argument passed, 0 means: do not skip any args
    {
        if (g_sOptionArg == "-distance")
        {
            m_oDetector.m_poGoniometer->nSetDistance(g_afOptionBuf[0]);
            fDistance = g_afOptionBuf[0];
        } 
        else if (g_sOptionArg == "-swing2theta")
        {
            m_oDetector.m_poGoniometer->nSetSwing(g_afOptionBuf[0]);
            f2ThetaSwing = g_afOptionBuf[0];            
        } 
        else if ((g_sOptionArg == "-reso-print") || (g_sOptionArg == "-reso-printonly"))
        {
            bComputeResoOnly = true;
        } 
        else if (g_sOptionArg == "-reso") 
        {
            fResoHigh = min(g_afOptionBuf[0],g_afOptionBuf[1]);
            fResoLow = max(g_afOptionBuf[0],g_afOptionBuf[1]);
        } 
        else if (g_sOptionArg == "-reso-outer")
        {
            bComputeReso = true;
            m_nMethod = eResoOuter;
        }
        else if (g_sOptionArg == "-reso-corners" || g_sOptionArg == "-reso-upperEqualsHighest") 
        {
            bComputeReso = true;
            m_nMethod = eResoOuter;
        } 
        else if (g_sOptionArg == "-reso-inner") 
        {
            bComputeReso = true;
            m_nMethod = eResoInner;
        } 
        else if (g_sOptionArg == "-reso-edges") 
        {
            bComputeReso = true;
            m_nMethod = eResoInner;
        } 
        else if (g_sOptionArg == "-reso-minarc")
        {
            bComputeReso = true;
            m_nMethod = eResoChiCoverage;
            m_fDebyeRingCoverage = g_afOptionBuf[0];
        }
        else if (g_sOptionArg == "-reso-2theta") 
        {
            bCompute2Theta = true;
            fTargetReso2Theta = g_afOptionBuf[0];
            
            if( g_afOptionBuf.size() >=2 )
                a2fRange2Theta[0] = g_afOptionBuf[1];
            else
                a2fRange2Theta[0] = 0.0;
            
            if (g_afOptionBuf.size() >=3 )
                a2fRange2Theta[1] = g_afOptionBuf[2];
            else
                a2fRange2Theta[1] = 90 + a2fRange2Theta[0];
        } 
        else if (g_sOptionArg == "-reso-distance")
        {
            if (g_afOptionBuf.size()>=1)
            {
                bComputeDistance = true;
                fTargetResoDist = g_afOptionBuf[0];
                if (g_afOptionBuf.size()>=2)
                    a2fRangeDist[0] = g_afOptionBuf[1];
                else
                    a2fRangeDist[0] = 10.0;
                if (g_afOptionBuf.size()>=3)
                    a2fRangeDist[1] = g_afOptionBuf[2];
                else
                    a2fRangeDist[1] = 100.0 + a2fRangeDist[0];
            }
            else 
            {
                bComputeDistanceSpotSeparation = true;
                nComputeDistanceSpotSeparationModeSpot = eComputeDistanceSpotSeparationDefault;
                nComputeDistanceSpotSeparationModeAxis = eComputeDistanceSpotSeparationDefault;
                nComputeDistanceSpotSeparation = 10;
            };
        } 
        else if (g_sOptionArg == "-reso-worstcase")
        {
            if (!bComputeDistanceSpotSeparation) 
                nStat = 2;
            if (nComputeDistanceSpotSeparationModeAxis == eComputeDistanceSpotSeparationDefault)
                nComputeDistanceSpotSeparationModeAxis = eComputeDistanceSpotSeparationWorstCase;
            else
                nComputeDistanceSpotSeparationModeSpot = eComputeDistanceSpotSeparationWorstCase;
        } 
        else if (g_sOptionArg == "-reso-bestcase") 
        {
            if (!bComputeDistanceSpotSeparation) 
                nStat = 2;
            if (nComputeDistanceSpotSeparationModeAxis == eComputeDistanceSpotSeparationDefault)
                nComputeDistanceSpotSeparationModeAxis = eComputeDistanceSpotSeparationBestCase;
            else
                nComputeDistanceSpotSeparationModeSpot = eComputeDistanceSpotSeparationBestCase;
            
        } 
        else if (g_sOptionArg == "-reso-averagecase") 
        {
            if (!bComputeDistanceSpotSeparation) 
                nStat = 2;
            if (nComputeDistanceSpotSeparationModeAxis == eComputeDistanceSpotSeparationDefault)
                nComputeDistanceSpotSeparationModeAxis = eComputeDistanceSpotSeparationAverageCase;
            else
                nComputeDistanceSpotSeparationModeSpot = eComputeDistanceSpotSeparationAverageCase;

        } 
        else if (g_sOptionArg == "-reso-spotseparation") 
        {
            if (!bComputeDistanceSpotSeparation) 
                nStat = 2;
            nComputeDistanceSpotSeparation = (int)g_afOptionBuf[0];
        } 
        else if (g_sOptionArg == "-reso-multi2theta") 
        {
            bCompute2ThetaMultipleSwing = true;
            if (g_afOptionBuf.size()>=1)
                fCompute2ThetaMultipleSwingOverlap = g_afOptionBuf[0];
            else
                fCompute2ThetaMultipleSwingOverlap = 20.0;
        }

    }
    
    if (nStat != 2)
        nStat = 0;
    else
        printf("ERROR:  Could not parse resolution options?\n");   

    // Consistency check.
    if( bCompute2Theta && (bComputeDistance || bComputeDistanceSpotSeparation) ) 
    {
        printf("ERROR:  2theta determination and distance determination are exclusive\n");
        printf("ERROR:  Please choose one or the other.\n");
        nStat = 1;
    }


    if( !nStat ) 
    {
        if( bCompute2Theta ) 
        {
            double      fTemp = 0.0;
            f0 = fCalcBest2ThetaSwingForTargetHighReso(a2fRange2Theta[0],a2fRange2Theta[1],fTargetReso2Theta, fTemp);        
            if (f0 > -999.0) 
            {
                f2ThetaSwing = f0;        
                if (nCalcReso(fResoHigh,fResoLow,f2ThetaSwing))
                    nStat = 2;
            } 
            else 
            {
                nCalcReso(fResoHigh,fResoLow);
                nStat = 2;
            }
        } 
        else if( bComputeDistance || bComputeDistanceSpotSeparation ) 
        {
            if( bComputeDistance )
                f0 = fCalcBestDistance(a2fRangeDist[0],a2fRangeDist[1],fTargetResoDist);
            else
                f0 = fCalcBestDistanceSpotSeparation(nComputeDistanceSpotSeparation,
                                                     nComputeDistanceSpotSeparationModeAxis,
                                                     nComputeDistanceSpotSeparationModeSpot);
            
            if( f0 > -999.0 )
            {
                fDistance = f0;        
                
                if( nCalcReso(fResoHigh,fResoLow,-999.0,fDistance) ) 
                    nStat = 2;                
            } 
            else 
            {
                nCalcReso(fResoHigh,fResoLow);
                
                nStat = 2;
            }
        }
        else if (bComputeReso)
        {
            if (nCalcReso(fResoHigh,fResoLow)) 
                nStat = 2;
        } 

        if (bCompute2ThetaMultipleSwing)
        {
            itr<double>         afSwings;
            // If we have multiple 2theta swings, then compute all swing values.
            // No further processing will be done with these swings.
            if( !nCalcBest2ThetaSwings(f2ThetaSwing, fCompute2ThetaMultipleSwingOverlap, afSwings) )
            {
                if (nCalcReso(fResoHigh,fResoLow,f2ThetaSwing))
                    nStat = 2;
            }
        }
        
        printf("---------------------------------------------------------------\n");
        printf("INFO:  Resolution range from %.2lf to %.2lf Angstrom selected.\n",fResoLow,fResoHigh);
        if ( (f2ThetaSwing > 0.0) && bCompute2Theta)
        {
            printf("INFO:  The input crystal to detector distance is %.2f mm.\n",
           m_oDetector.fGetDistance());
            printf("INFO:  The calculated detector swing angle is %.2lf deg.\n",f2ThetaSwing);
        }
        else if (bCompute2Theta) 
            printf("ERROR: No 2theta swing from %.2lf to %.2lf would attain a resolution of %.2lf\n",a2fRange2Theta[0],a2fRange2Theta[1],fTargetReso2Theta);
        if ( (fDistance > 0.0) && bComputeDistance)
        {
            printf("INFO:  The input detector swing angle is %.2f deg.\n",
           m_oDetector.fGetSwing());
            printf("INFO:  The calculated crystal to detector distance is %.2lf mm.\n",fDistance);
        }
        else if (bComputeDistance)
            printf("ERROR: No crystal-detector from %.2lf to %.2lf would attain a resolution of %.2lf\n",a2fRangeDist[0],a2fRangeDist[1],fTargetResoDist);
        printf("---------------------------------------------------------------\n");
    }

    if( nStat ) 
    {
        printf("ERROR:  Problem during resolution calculation?\n");
    }
    
    
    if( bComputeResoOnly  &&  !nStat )
        nStat = 1;

    return nStat;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Discover the max, min, and average spot dimensions.
bool CExperimentalStrategy::bGetSpotSizeFromHeader(SPOT_SIZE_INFO& stSpotSizeInfo)  
{
    C3DdataDetectorProfile      oDetectorProfiles;

    if( 0 != oDetectorProfiles.nInitValues(m_oHeader) )
        return false;

    stSpotSizeInfo.fMaxSpotSize = 0.0;
    stSpotSizeInfo.fAverageSpotSize = 0.0;
    stSpotSizeInfo.fMinSpotSize = oDetectorProfiles.m_aoDetectorAreasPrint[0][0].fSize[0];

    for(int nx = 0; nx < C3Ddata_Print_Div_0; nx ++)
    {
        for(int ny = 0; ny < C3Ddata_Print_Div_1; ny++)
        {
            stSpotSizeInfo.fMaxSpotSize     = max(stSpotSizeInfo.fMaxSpotSize, oDetectorProfiles.m_aoDetectorAreasPrint[nx][ny].fSize[0]);
        
            stSpotSizeInfo.fAverageSpotSize += oDetectorProfiles.m_aoDetectorAreasPrint[nx][ny].fSize[0];
        
            stSpotSizeInfo.fMinSpotSize     = min(stSpotSizeInfo.fMinSpotSize, oDetectorProfiles.m_aoDetectorAreasPrint[nx][ny].fSize[0]);
        }
    }

    stSpotSizeInfo.fAverageSpotSize /= (C3Ddata_Print_Div_0 * C3Ddata_Print_Div_1);
    
    return true;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


