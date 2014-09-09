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
// CSphericalHarmonic.cc        Initial author: Thaddeus Niemeyer 12-24-02
//  This file contains new spherical harmonic routines.
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

#include "CSphericalHarmonic.h"
#include "dtsvd.h"
#include "dtrekvec.h"
#include "Cstat.h"

// Azimuthal angle (Azi) is the longitudinal coordiante [0 to 2*pi].
// Polar angle (Polar) is the coladitudinal coordinate [0 to pi].

CSphericalHarmonic::CSphericalHarmonic() {
    vSetHarmonics(1,1,2);
    vSetScaleInfo(3.0,50.0,5,0.005);
    m_nScans = 1;                      // Don't know this till we call nScale().
	vSetPrint(1);
    vSetCentroSymmetric(1);
};

int CSphericalHarmonic::nCalcIncidentHarmonics(double fRot) {
    int nx;

    m_nIncidentCoeffs = m_nIncidentHarmonicOrder*2 + 1;
    m_afIncidentHarmonicCoeffs.setsize(m_nIncidentCoeffs);
    m_pfIncidentHarmonicCoeffs = &m_afIncidentHarmonicCoeffs[0];
    for (nx = 0; nx < m_nIncidentCoeffs; nx++) {
        if (nx <= m_nIncidentCoeffs/2)
            m_pfIncidentHarmonicCoeffs[nx] = cos(Gs_dRADIANS_PER_DEGREE*fRot*nx);
        else
            m_pfIncidentHarmonicCoeffs[nx] = sin(Gs_dRADIANS_PER_DEGREE*fRot*(nx - m_nIncidentCoeffs/2));
    };
    return 0;
};

int CSphericalHarmonic::nCalcHarmonics(double fPolar,double fAzi) {
    double afSinPolar[g_nMaxHarmonicOrder+1];
    double afCosPolar[g_nMaxHarmonicOrder+1];
    double afSinAzi[g_nMaxHarmonicOrder+1];
    double afCosAzi[g_nMaxHarmonicOrder+1];
    double afPolarHarmonics[g_nMaxHarmonicOrder+1][g_nMaxHarmonicOrder+1];

    int nx,ny;
    int nCoeff;

    afCosPolar[1] = cos(Gs_dRADIANS_PER_DEGREE*fPolar);
    afSinPolar[1] = sin(Gs_dRADIANS_PER_DEGREE*fPolar);
    for (nx=2;nx <= m_nHarmonicOrder; nx++) {
        afCosPolar[nx] = afCosPolar[1]*afCosPolar[nx-1];
        afSinPolar[nx] = afSinPolar[1]*afSinPolar[nx-1];
    };
    for (nx = 1; nx <= m_nHarmonicOrder; nx++) {
        afCosAzi[nx] = cos(Gs_dRADIANS_PER_DEGREE*nx*fAzi);
        afSinAzi[nx] = sin(Gs_dRADIANS_PER_DEGREE*nx*fAzi);
    };

    m_nCoeffs = 0;

    while (1) {
        afPolarHarmonics[0][0] = 1.0;
        m_nCoeffs += 1;
        if (m_nHarmonicOrder < 1)
            break;
        afPolarHarmonics[1][0] = afCosPolar[1];
        afPolarHarmonics[1][1] = afSinPolar[1];
        m_nCoeffs += 3;
        if (m_nHarmonicOrder < 2)
            break;
        afPolarHarmonics[2][0] = 3*afCosPolar[2] - 1.0;
        afPolarHarmonics[2][1] = afSinPolar[1]*afCosPolar[1];
        afPolarHarmonics[2][2] = afSinPolar[2];
        m_nCoeffs += 5;
        if (m_nHarmonicOrder < 3)
            break;
        afPolarHarmonics[3][0] = 5.0*afCosPolar[3] - 3.0*afCosPolar[1];
        afPolarHarmonics[3][1] = (5.0*afCosPolar[2] - 1.0)*afSinPolar[1];
        afPolarHarmonics[3][2] = afSinPolar[2]*afCosPolar[1];
        afPolarHarmonics[3][3] = afSinPolar[3];
        m_nCoeffs += 7;
        if (m_nHarmonicOrder < 4)
            break;
        afPolarHarmonics[4][0] = 35.0*afCosPolar[4] - 30.0*afCosPolar[2] + 3.0;
        afPolarHarmonics[4][1] = (7.0*afCosPolar[3] - 3.0*afCosPolar[1])*afSinPolar[1];
        afPolarHarmonics[4][2] = (7.0*afCosPolar[2] - 1.0)*afSinPolar[2];
        afPolarHarmonics[4][3] = afSinPolar[3]*afCosPolar[1];
        afPolarHarmonics[4][4] = afSinPolar[4];
        m_nCoeffs += 9;
        if (m_nHarmonicOrder < 5)
            break;
        afPolarHarmonics[5][0] = 63.0*afCosPolar[5] - 70.0*afCosPolar[3] + 15.0*afCosPolar[1];
        afPolarHarmonics[5][1] = (21.0*afCosPolar[4] - 14.0*afCosPolar[2] + 1.0)*afSinPolar[1];
        afPolarHarmonics[5][2] = (3.0*afCosPolar[3] - afCosPolar[1])*afSinPolar[2];
        afPolarHarmonics[5][3] = (9.0*afCosPolar[2] - 1.0)*afSinPolar[3];
        afPolarHarmonics[5][4] = afSinPolar[4]*afCosPolar[1];
        afPolarHarmonics[5][5] = afSinPolar[5];
        m_nCoeffs += 11;
        if (m_nHarmonicOrder < 6)
            break;
    };
    m_afHarmonicCoeffs.setsize(m_nCoeffs);
    m_pfHarmonicCoeffs = &m_afHarmonicCoeffs[0];

    nCoeff = 0;
    for (nx= 0; nx <= m_nHarmonicOrder;nx++) {
        for (ny = -nx; ny <= nx; ny++) {
            if (ny < 0) 
                m_pfHarmonicCoeffs[nCoeff++] = afPolarHarmonics[nx][-ny]*afSinAzi[-ny];
            else if (ny > 0)
                m_pfHarmonicCoeffs[nCoeff++] = afPolarHarmonics[nx][ny]*afCosAzi[ny];
            else
                m_pfHarmonicCoeffs[nCoeff++] = afPolarHarmonics[nx][ny];
        };
    };
    return 0;
};
int CSphericalHarmonic::nFit(itr<double>& afValues,itr<double>& afValuesSigma,itr<double>& afPolar,itr<double>& afAzi,itr<double>& afScanRot,itr<int>& anUse) {
    double* pfPolar;
    double* pfAzi;
    double* pfValues;
    double* pfSigma;

    itr<double> afFuncs;
    itr<double*> apfFuncs;
    double** ppfFuncs;
    
    double fPolar,fAzi;
	double f0;
    int nPoint;
    int nPointsUsed;
    int nIncident;
    int nHarmonic;
    int nx;
    int nSizeMult;
    int nSym;
	int nStat;

    itr<double> afValuesUsed;
    itr<double> afSigmaUsed;

    if (!afValues.size())
        return 1;


    pfValues = &afValues[0];
    pfSigma = (afValuesSigma.size())?(&afValuesSigma[0]):NULL;
    pfPolar = &afPolar[0];
    pfAzi = &afAzi[0];

    // Discover the # of harmonic terms.
    if (nCalcHarmonics(0.0,0.0))
        return 1;
    if (nCalcIncidentHarmonics(0.0))
        return 1;

    nSizeMult = 1 + max(0,min(1,(m_nCentroSymmetric)));
    m_nTotalCoeffs = m_nIncidentCoeffs*m_nCoeffs;

    // Possibly a very large array!
    afFuncs.setsize(m_nTotalCoeffs*afValues.size()*nSizeMult);

    for (nx = 0; nx < m_nTotalCoeffs; nx++) {
        apfFuncs[nx] = &afFuncs[nx*afValues.size()*nSizeMult];
    };
    ppfFuncs = &apfFuncs[0];

    // Fill in the points.
    afValuesUsed.setsize(afValues.size()*nSizeMult);
    afSigmaUsed.setsize(afValues.size()*nSizeMult);
    -afValuesUsed;
    -afSigmaUsed;

    for (nPoint = 0,nPointsUsed = 0; nPoint < afValues.size(); nPoint++) {
        if ((!anUse.size()) || (anUse[nPoint])) {
            for (nSym = 0; nSym < nSizeMult; nSym++) {
                afValuesUsed + pfValues[nPoint];
                afSigmaUsed + ((pfSigma)?(pfSigma[nPoint]):1.0);
                fPolar = pfPolar[nPoint];
                fAzi = pfAzi[nPoint];
                if (nSym) {
                    fPolar = fPolar + 180.0;
                    fAzi = 180.0 - fAzi;
                };                
                if (nCalcHarmonics(fPolar,fAzi))
                    return 1;
                if (nCalcIncidentHarmonics((afScanRot.size())?(afScanRot[nPoint]):0.0))
                    return 1;
                for (nIncident = 0; nIncident < m_nIncidentCoeffs; nIncident++) {
                    for (nHarmonic = 0; nHarmonic < m_nCoeffs; nHarmonic++) {
                        ppfFuncs[nHarmonic*m_nIncidentCoeffs + nIncident][nPointsUsed] = m_pfHarmonicCoeffs[nHarmonic]*m_pfIncidentHarmonicCoeffs[nIncident];
                    };
                };
                nPointsUsed++;
            };
        };
    };
    if (nPointsUsed < m_nTotalCoeffs)
        return 1;
    m_afFit.setsize(m_nTotalCoeffs);
    m_afFitSigma.setsize(m_nTotalCoeffs);
    pfValues = &afValuesUsed[0];
    pfSigma = &afSigmaUsed[0];

    // Run the linear least squares.
    -m_anOutlierReject;
	f0 = ms_fLeastSquaresOutlierRejectIoverSig;
	ms_fLeastSquaresOutlierRejectIoverSig = m_fRejectSigma;
    nStat = nSolveLeastSquares(m_nTotalCoeffs,nPointsUsed,pfValues,pfSigma,ppfFuncs,&m_afFit[0],&m_afFitSigma[0],NULL,NULL,&m_anOutlierReject);
	ms_fLeastSquaresOutlierRejectIoverSig = f0;
    if ((nSizeMult > 1) && (m_anOutlierReject.size() % nSizeMult == 0)) {
        for (nPoint = 0; nPoint < m_anOutlierReject.size()/nSizeMult; nPoint++)
            m_anOutlierReject[nPoint] = m_anOutlierReject[nPoint*nSizeMult];
    };
    return 0;
};

void CSphericalHarmonic::vVec3DToSpherical(double* pfVec,double& fPolar,double& fAzi,double* pfAxisMatIn) {
    double a3x3fDefaultAxis[3][3] = {{ 1.0,0.0,0.0},{ 0.0,1.0,0.0},{ 0.0,0.0,1.0}};
    double* pfAxisMat;
    double a3fProjVec[3];
    int nx;
    double f0,f1;
    if (pfAxisMatIn)
        pfAxisMat = pfAxisMatIn;
    else
        pfAxisMat = & a3x3fDefaultAxis[0][0];

    f0 = fDot3D(pfVec,pfAxisMat + 6);
    fPolar = acos(max(-1.0,min(1.0,f0)))/Gs_dRADIANS_PER_DEGREE;
    for (nx=0;nx<3;nx++) 
        a3fProjVec[nx] = pfVec[nx] - f0*pfAxisMat[6 + nx];
    fNormVec3D(a3fProjVec);
    f0 = fDot3D(a3fProjVec,pfAxisMat);
    f1 = fDot3D(a3fProjVec,pfAxisMat + 3);
    fAzi = acos(max(-1.0,min(1.0,f0)))/Gs_dRADIANS_PER_DEGREE;
    if (f1 < 0.0)
        fAzi = 360.0 - fAzi;

};
void CSphericalHarmonic::vSphericalToVec3D(double fPolar,double fAzi,double* pfVec,double* pfAxisMatIn) {
    double a3x3fDefaultAxis[3][3] = {{ 1.0,0.0,0.0},{ 0.0,1.0,0.0},{ 0.0,0.0,1.0}};
    double* pfAxisMat;
    double f0;
    int nx;
    if (pfAxisMatIn)
        pfAxisMat = pfAxisMatIn;
    else
        pfAxisMat = & a3x3fDefaultAxis[0][0];

    fPolar *= Gs_dRADIANS_PER_DEGREE;
    fAzi *= Gs_dRADIANS_PER_DEGREE;
    
    f0 = cos(fAzi)*sin(fPolar);
    for (nx=0;nx<3;nx++)
        pfVec[nx] = pfAxisMat[nx]*f0;
    f0 = sin(fAzi)*sin(fPolar);
    for (nx=0;nx<3;nx++)
        pfVec[nx] += pfAxisMat[3 + nx]*f0;
    f0 = cos(fPolar);
    for (nx=0;nx<3;nx++)
        pfVec[nx] += pfAxisMat[6 + nx]*f0;
};

void CSphericalHarmonic::vSetGroupInfo(int& nCode,int nGroup,int nScan,int nReject) {

	nCode = (int) (((unsigned int) (nGroup + nReject*g_nGroupMax)) +  ((unsigned int) nScan)*((unsigned int) g_nGroupMax)*g_nSphericalFlagMax); 
};

void CSphericalHarmonic::vGetGroupInfo(int nCode,int& nGroup,int& nScan,int& nReject) {
	unsigned int nCodeu;

	nCodeu = (unsigned int) nCode;

	nScan = (nCodeu / (g_nGroupMax*g_nSphericalFlagMax));
	nCodeu -= nScan*(g_nGroupMax*g_nSphericalFlagMax);
	nReject = (nCodeu/g_nGroupMax);
	nCodeu -= nReject*g_nGroupMax;
	nGroup = (int) nCodeu;
};

int CSphericalHarmonic::nGetSetHarmonicCoeffs(bool bGet,itr<double>& afHarmonicCoeffs) {
    if (bGet) {
        -afHarmonicCoeffs + m_afHarmonicCoeffs;
    } else {
        -m_afHarmonicCoeffs + afHarmonicCoeffs;
    };
    return 0;
};

int CSphericalHarmonic::nSimulate(itr<double>& afValues,
                                  itr<double>& afPolar,
                                  itr<double>& afAzi,
                                  itr<double>& afRotAngle) 
{
    double* pfValues;
    double* pfPolar;
    double* pfAzi;
    double* pfFit;
    int nx;
    int nIncident,nHarmonic;
    double f0;
    afValues.setsize(afPolar.size());

    pfValues = &afValues[0];
    pfPolar = &afPolar[0];
    pfAzi = &afAzi[0];
    pfFit = & m_afFit[0];
    for (nx=0; nx < afPolar.size();nx++) {
        if (nCalcHarmonics(pfPolar[nx],pfAzi[nx]))
            return 1;
        if (nCalcIncidentHarmonics((afRotAngle.size())?(afRotAngle[nx]):0.0))
            return 1;

        f0 = 0.0;
        for (nIncident = 0; nIncident < m_nIncidentCoeffs; nIncident++) {
            for (nHarmonic = 0; nHarmonic < m_nCoeffs; nHarmonic++) {
                f0 += m_pfHarmonicCoeffs[nHarmonic]*m_pfIncidentHarmonicCoeffs[nIncident]*pfFit[nHarmonic*m_nIncidentCoeffs + nIncident];
            };
        };
        pfValues[nx] = f0;
    };

    return 0;
};

int CSphericalHarmonic::nList() {
    int nx;
    
    for (nx = 0; nx < m_nCoeffs*m_nIncidentCoeffs;nx++) {
        printf("%3d = %5.4lf\n",
        nx,m_afFit[nx]);
    };
    printf("\n");


    return 0;
};

int CSphericalHarmonic::nTest() {
    itr<double> afAzi;
    itr<double> afPolar;
    itr<double> afScanRot;
    itr<double> afAziSim;
    itr<double> afPolarSim;
    itr<double> afValues;
    itr<double> afValuesOut;
    itr<double> afValuesSigma;
    itr<int>    anInUse;
    const int nTestPoints = 500;
    int nPoint;
    int nx,ny;
    
    double a3x3fSimulateAxes[3][3];

    while (1) {
    // Make a set of simulation axes for our function.
    for (nx=0;nx<3;nx++)
        a3x3fSimulateAxes[0][nx] = (rand() % 10000) - 5000;
    fNormVec3D(a3x3fSimulateAxes[0]);
    vBuildBasis3D(a3x3fSimulateAxes[0],a3x3fSimulateAxes);
    vNormMat3D(a3x3fSimulateAxes);

    nCalcHarmonics(0.0,0.0);
    nCalcIncidentHarmonics(0.0);
    -m_afFit;
    // Generate a test set of harmonic coefficients.
    for (nx = 0; nx < m_nCoeffs; nx++)  {
        for (ny = 0; ny < m_nIncidentCoeffs;ny++) {
            m_afFit + (((rand() % 10)) + 1);
        };
    };

    -afScanRot;
    -afAzi;
    -afPolar;
    -afAziSim;
    -afPolarSim;
    // Generate random points on a sphere.
    for (nPoint = 0; nPoint < nTestPoints; nPoint++) {
        double a3fVec[3];
        double fPolar,fAzi;
        for (nx=0; nx< 3;nx++)
            a3fVec[nx] = (rand() % 10000) - 5000;
        fNormVec3D(a3fVec);
        vVec3DToSpherical(a3fVec,fPolar,fAzi);
        afAzi + fAzi;
        afPolar + fPolar;
        afScanRot + (double) (rand() % 360);
        vVec3DToSpherical(a3fVec,fPolar,fAzi,&a3x3fSimulateAxes[0][0]);
        afAziSim + fAzi;
        afPolarSim + fPolar;

    };
    // Simulate the data.
    nSimulate(afValues,afPolarSim,afAziSim,afScanRot);
    // Make fake sigma estimates.
    -afValuesSigma;
    for (nx=0;nx< afValues.size(); nx++)
        afValuesSigma + 1.0;

    nList();
    
    nFit(afValues,afValuesSigma,afPolar,afAzi,afScanRot,anInUse);

    nSimulate(afValuesOut,afPolar,afAzi,afScanRot);

    for (nx=0;nx<20;nx++) 
        printf("%.5lf %.5lf\n",afValues[nx],afValuesOut[nx]);

    nList();

    };

    return 0;
};

int CSphericalHarmonic::nScale(itr<double>& afValues,itr<double>& afValuesSigma,itr<double>& afPolar,itr<double>& afAzi) {
    itr<int> anInUse;
    int nPass;
    const int nMaxPasses = 2;
    itr<double> afScanRot;      // This remains of zero length.
    itr<double> afSimulateValues;
    int nPoints;
    nPoints = afValues.size();
    
    for (nPass = 0; nPass < nMaxPasses; nPass++) {
        if (nFit(afValues,afValuesSigma,afPolar,afAzi,afScanRot,anInUse))
            return 1;
        if (nSimulate(afSimulateValues,afPolar,afAzi,afScanRot))
            return 1;
    };
    return 0;

};


int CSphericalHarmonic::nScale(itr<double>& afValues,itr<double>& afValuesSigma,itr<double>& afPolar,itr<double>& afAzi,itr<int>& anGroup,itr<double>& afScanRot,itr<double>& afScaleFactorsOut) {
    
    // Per scan arrays.
    itr<double> afDiffValues;
    itr<double> afDiffValuesSigma;
    itr<double> afDiffValuesPolar;
    itr<double> afDiffValuesAzi;
    itr<double> afDiffValuesScanRot;

    itr<double> afLastCorrectFactor;
    itr<double> afThisCorrectFactor;
    itr<int>    anDiffInUse;

    int nLoop;
    int nScan,nScans;
    int nRef;
    int nRefGroup,nRefScan,nRefReject;

    // Discover the # of scans.
    nScans = 1;
    for (nRef = 0; nRef < anGroup.size(); nRef++) {
        vGetGroupInfo(anGroup[nRef],nRefGroup,nRefScan,nRefReject);
        nScans = max(1 + nRefScan,nScans);
    };

    // Set the last correction factor array to be 1.0 everywhere.
    afThisCorrectFactor.setsize(afValues.size());
    for (nRef = 0; nRef < anGroup.size(); nRef++) {
        afThisCorrectFactor[nRef] = 1.0;
        afLastCorrectFactor[nRef] = 1.0;
    };
    nLoop = 0;


    do {

        for (nScan = 0; nScan < nScans; nScan++) {
            -afDiffValues;
            -afDiffValuesSigma;
            -afDiffValuesPolar;
            -afDiffValuesAzi;
            -afDiffValuesScanRot;


            double fGroupSum;
            double fSumInvSigma;
            double fLastCorrectSum;
            int nGroupContrib;
            int nGroupContribInScan;
            int nLastRefGroup;
            int nLastRef;
            int nRef2;

            nLastRefGroup = -1;
            nGroupContrib = 0;
            nGroupContribInScan = 0;

            for (nRef = 0; nRef <= anGroup.size(); nRef++) {
                if (nRef < anGroup.size())
                    vGetGroupInfo(anGroup[nRef],nRefGroup,nRefScan,nRefReject);
                if ((nRefGroup != nLastRefGroup) || (nRef == anGroup.size())) {
                    
                    // Process the last group of reflections.
                    if ((nGroupContrib >= 2) && (nGroupContribInScan >= 1)) {
                        fLastCorrectSum /= nGroupContrib;
                        fGroupSum/= (fSumInvSigma*fLastCorrectSum);
                        for (nRef2 = nLastRef; nRef2 < nRef; nRef2++) {
                            vGetGroupInfo(anGroup[nRef2],nRefGroup,nRefScan,nRefReject);
                            if ((nRefScan == nScan) && (!(nRefReject & g_nSphericalFlagUserReject)) && (!(nRefReject & g_nSphericalFlagExclude))) {
                                afDiffValues + (afValues[nRef2]/fGroupSum);
                                afDiffValuesSigma + (afValuesSigma[nRef2]/fGroupSum);
                                afDiffValuesPolar + afPolar[nRef2];
                                afDiffValuesAzi + afAzi[nRef2];
                                afDiffValuesScanRot + afScanRot[nRef2];
                            };
                        };
                    } else {
						if (nGroupContrib >= 1) {
							for (nRef2 = nLastRef; nRef2 < nRef; nRef2++) {
								vGetGroupInfo(anGroup[nRef2],nRefGroup,nRefScan,nRefReject);
								if (nScan == nRefScan) {
									// Set the exclusion status of this reflection.
									nRefReject |= g_nSphericalFlagTooLowRedundancy;
									vSetGroupInfo(anGroup[nRef2],nRefGroup,nRefScan,nRefReject);
								};
							};
						};
					};

                    fGroupSum = 0.0;
                    nGroupContrib = 0;
                    nGroupContribInScan = 0;
                    fLastCorrectSum = 0.0;
                    fSumInvSigma = 0.0;
                    nLastRefGroup = nRefGroup;
                    nLastRef = nRef;

                    if (nRef == anGroup.size())
                        break;

                    // Get the info again.
                    vGetGroupInfo(anGroup[nRef],nRefGroup,nRefScan,nRefReject);
                } 

                if (!(nRefReject & g_nSphericalFlagUserReject)) {
                    fGroupSum += afValues[nRef]/afLastCorrectFactor[nRef]/(afValuesSigma[nRef]*afValuesSigma[nRef]);
                    fLastCorrectSum += 1.0/afLastCorrectFactor[nRef];
                    fSumInvSigma += 1.0/(afValuesSigma[nRef]*afValuesSigma[nRef]);
                    nGroupContrib ++;

					if (afValues[nRef]/max(1.0,afValuesSigma[nRef]) < m_fSigma) {
						nRefReject |= g_nSphericalFlagExclude;
						// Set the exclude status.
						vSetGroupInfo(anGroup[nRef],nRefGroup,nRefScan,nRefReject);
					} 
                    
                    if ((nScan == nRefScan) && (!(nRefReject & g_nSphericalFlagExclude)))
                        nGroupContribInScan++;
                };
            };      
               

            if (nFit(afDiffValues,afDiffValuesSigma,afDiffValuesPolar,afDiffValuesAzi,afDiffValuesScanRot,anDiffInUse))
                return 1;

            -afDiffValuesPolar;
            -afDiffValuesAzi;
            -afDiffValuesScanRot;
            // Construct the list of absorbtion factors for each reflection in the scan.
            for (nRef = 0; nRef < anGroup.size(); nRef++) {
                vGetGroupInfo(anGroup[nRef],nRefGroup,nRefScan,nRefReject);
                if (nRefScan == nScan) {
                    afDiffValuesPolar + afPolar[nRef];
                    afDiffValuesAzi + afAzi[nRef];
                    afDiffValuesScanRot + afScanRot[nRef];
                };
            };
            if (nSimulate(afDiffValues,afDiffValuesPolar,afDiffValuesAzi,afDiffValuesScanRot))
                return 1;
            double fMaxAbsorb,fMinAbsorb;
            double fPercentChange;
            double f0,f1;

            fPercentChange = 0.0;
            for (nRef = 0,nRef2= 0; nRef < anGroup.size(); nRef++) {
                vGetGroupInfo(anGroup[nRef],nRefGroup,nRefScan,nRefReject);
                if (nRefScan == nScan) {
                    f0 = afDiffValues[nRef2];
                    f1 = 1.0;
                    afThisCorrectFactor[nRef] = f0;
                    if ((nRef2 == 0) || (f0*f1>fMaxAbsorb))
                        fMaxAbsorb = f0*f1;
                    if ((nRef2 == 0) || (f0*f1<fMinAbsorb))
                        fMinAbsorb = f0*f1;
					if (m_anOutlierReject[nRef2]) 
						nRefReject |= g_nSphericalFlagOutlierReject;
					else
						nRefReject &= ~g_nSphericalFlagOutlierReject;
					vSetGroupInfo(anGroup[nRef],nRefGroup,nRefScan,nRefReject);
                    nRef2++;
                    fPercentChange += ABS(afThisCorrectFactor[nRef] - afLastCorrectFactor[nRef]);
					
                };
            };
            fPercentChange /= nRef2;

    
			if (m_nPrintFlag)
				printf("Scan %d Maximum Absorb = %.3lf Minimum Absorb = %.3lf Change = %.2lf%%\n",nScan+1,fMaxAbsorb,fMinAbsorb,fPercentChange*100.0);
           
        };

        // Copy the last set of correction factors over.
        afLastCorrectFactor.setsize(afValues.size());
        for (nRef = 0; nRef < anGroup.size(); nRef++) {
            afLastCorrectFactor[nRef] = afThisCorrectFactor[nRef];
        };
        


        nLoop++;

    } while (nLoop < m_nMaxIterations);


    afScaleFactorsOut.setsize(afThisCorrectFactor.size());
    for (nRef = 0; nRef < anGroup.size(); nRef++) {
        afScaleFactorsOut[nRef] = afLastCorrectFactor[nRef];
    };
    return 0;
};
