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
// Ccrystal.cc            Initial author: J.W. Pflugrath           01-May-1995
//    This file contains the member functions of class Ccrystal.
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

//+ToDo
//
//   Error messages need implementing
//
//+Include files

#include "Dtrek.h"
#include "Ccrystal.h"            // Class definition and prototypes
                                 //  Ccrystal.h includes others
#include "Crefln.h"
#include "dtrekdefs.h"

#if !defined(VC6)
using std::cout;
using std::endl;
using std::flush;
#endif


//+Definitions, constants, and initialization of static member variables

Cstring Ccrystal::ms_sCrystalPrefix        = D_K_CrystalPrefix;
Cstring Ccrystal::ms_sCrystalKey           = D_K_CrystalKey;
Cstring Ccrystal::ms_sCrystalDescription   = D_K_CrystalDescription;

// These strings are not determined until the twin componant is known.  
// They are loaded in every call to nInitValues(Cimage_header&)
// The initial values loaded here are for the first crystal componant.

Cstring Ccrystal::ms_sCrystalXUnitCell     = D_K_CrystalPrefix D_K_CrystalXUnitCell;
Cstring Ccrystal::ms_sCrystalXMosaicity    = D_K_CrystalPrefix D_K_CrystalXMosaicity;     
Cstring Ccrystal::ms_sCrystalXOrientAngles = D_K_CrystalPrefix D_K_CrystalXOrientAngles;  
Cstring Ccrystal::ms_sCrystalXOrientVectors= D_K_CrystalPrefix D_K_CrystalXOrientVectors; 
Cstring Ccrystal::ms_sCrystalXSpacegroup   = D_K_CrystalPrefix D_K_CrystalXSpacegroup;
Cstring Ccrystal::ms_sCrystalXNumTwinLaws  = D_K_CrystalPrefix D_K_CrystalXNumTwinLaws;
Cstring Ccrystal::ms_sCrystalXTwinLaws	   = D_K_CrystalPrefix D_K_CrystalXTwinLaws;
Cstring Ccrystal::ms_sCrystalXTwinFraction = D_K_CrystalPrefix D_K_CrystalXTwinFraction;
Cstring Ccrystal::ms_sCrystalXRecipShift   = D_K_CrystalPrefix D_K_CrystalXRecipShift;

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Ccrystal::Ccrystal()
{
  (void) nInitValues();
}

Ccrystal::Ccrystal(Cimage_header& oHeader,int nTwin)
{
  (void) nInitValues(oHeader,nTwin);
}

Ccrystal::Ccrystal(Ccrystal& oOther) {
    Cimage_header oHeader;
    int nOrigTwin;

    nOrigTwin = oOther.m_nTwin;
    oOther.m_nTwin = 1;
    
    (void) oOther.nUpdateHeader(&oHeader);
    nInitValues(oHeader,1);
    oOther.m_nTwin = nOrigTwin;
    m_nTwin = nOrigTwin;
    m_nTwins = oOther.m_nTwins;
	m_nMaxTwinLaws = oOther.m_nMaxTwinLaws;
};

Ccrystal::~Ccrystal()
{
  if (NULL != m_poSpacegroup)
    {
      delete m_poSpacegroup;
      m_poSpacegroup = NULL;
    }
  m_eThe_State = eCrystal_unknown_state;
}

int Ccrystal::nInitValues()
{
  int i, j;       // Loop counters

  m_eThe_State = eCrystal_unknown_state;
  for (i = 0; i < 6; i++)
    {
      m_fCell[i]      = 0.0;
      m_fRecipCell[i] = 0.0;
    }

  for (j = 0; j < 3; j++)
    {
      for (i = 0; i < 3; i++)
        {
          m_fOrientVectors[j][i] = 0.0;    // Define axes for fOrientAngles as
          m_fOrientMatrix[j][i]  = 0.0;    // Define axes for fOrientAngles as
          m_fRotMatrix[j][i]     = 0.0;    // Define axes for fOrientAngles as
          m_fBMatrix[j][i]       = 0.0;
        }
      m_fOrientAngles[j]     = 0.0;
      m_fRotMatrix[j][j]     = 1.0;
      m_fOrientMatrix[j][j]  = 1.0;
      m_fOrientVectors[j][j] = 1.0;
      m_fBMatrix[j][j]       = 1.0;
	  for (i = 0; i < g_nMaxTwinLaws + 1; i++)
		  m_aa3fRecipShift[i][j] = 0.0;
	  for (i = 0; i < g_nMaxTwinLaws; i++)
		  m_aa3fTwinLaw[i][j] = 0.0;
    }

  m_fVolume        = 0.0f;
  m_fRecipVol      = 0.0f;
  m_nTwinLaws      = 0;
  for (i = 0; i < g_nMaxTwinLaws; i++)
	  m_afTwinFraction[i] = 0.0;

  m_fExposureTime  = 0.0f;
  m_sDescription   = "unknown";
  m_sCrystal_key   = "";

  m_fMosaicity = 0.3f;
  m_fEffictiveSpectralDispersion = 0.0;
  m_fMosaicityOffset = 0.0;
  
  m_nTwin = 1;
  m_nTwins = 1;
  m_nMaxTwinLaws = 0;

  m_poSpacegroup   = new Cspacegroup ();

  return (0);
}

int Ccrystal::nInitValues(Cimage_header& oHeader,int nTwin)
{
    int    nStat, nTemp;
    Cstring sTemp;
    double fTemp;
    float a3fMosHeader[4];
    
    
    (void) nInitValues();
    m_nTwin = nTwin;
    
    
    if (!oHeader.bIsAvailable())
    {
        cout << "Ccrystal ERROR: image header is not valid."
            << "  Cannot construct crystal!\n";
        nStat = 1;
    }
    else if (0 == oHeader.nGetValue(ms_sCrystalKey, &sTemp))
    {
        cout << "Csource ERROR: cannot use database key yet!\n";
        nStat = 2;
    }
    else
    {
        // Discover how many twins are in the crystal.  We will count the CRYSTAL_X_UNIT_CELL keyword's
        
        nFindNumTwins(oHeader);        
        if (nTwin>m_nTwins) {
            m_eThe_State = eCrystal_unknown_state;
            return (1);
        };

        
        nLoadKeywords(m_nTwin);
        
        // Try to get required information from the header
        
        nStat = 0;
        
        nTemp = oHeader.nGetValue(ms_sCrystalXUnitCell, 6, m_fCell);
        if (0 != nTemp)
        {
            cout << "Ccrystal WARNING: reading " << ms_sCrystalXUnitCell
                << " keyword!\n";
            nStat++;
        }
        else
        {
            m_fVolume = fCalcVolume();
        }
        
        nTemp = oHeader.nGetValue(ms_sCrystalXOrientAngles, 3, m_fOrientAngles);
        if (0 != nTemp)
        {
            cout << "Ccrystal WARNING: reading " << ms_sCrystalXOrientAngles
                << " keyword!\n";
            nStat++;
        }
        
        // Try to read mosaicity.  A single value is treated differently.
        
        fTemp = 0.0;
        nTemp = oHeader.nGetValue(ms_sCrystalXMosaicity,3,&a3fMosHeader[0]);
        if (!nTemp) {
            vSetMosaicity(a3fMosHeader[0],a3fMosHeader[1],a3fMosHeader[2]);
        } else {
            nTemp = oHeader.nGetValue(ms_sCrystalXMosaicity, &fTemp);
            if (0.0 >= fTemp) 
                fTemp = 0.3f;
            
            if (!nTemp) 
                vSetMosaicity(fTemp);
            else {
                cout << "Ccrystal WARNING: reading " << ms_sCrystalXMosaicity
                    << " keyword!\n";
                nStat++;
            }
        };
        // It is not required to have the following 2 keywords, so
        // do not flag if not found.
        
        nTemp = oHeader.nGetValue(ms_sCrystalXOrientVectors, 9,
            &m_fOrientVectors[0][0]);
        nTemp = oHeader.nGetValue(ms_sCrystalDescription, &m_sDescription);
        
        int nSpacegroup;
        nTemp = oHeader.nGetValue(ms_sCrystalXSpacegroup, &nSpacegroup);
        if (0 != nTemp)
        {
            cout << "Ccrystal WARNING: reading " << ms_sCrystalXSpacegroup
                << " keyword!\n";
            nStat++;
        }
        else
        {
            // Delete old spacegroup object before constructing new one
            
            if (NULL != m_poSpacegroup) delete m_poSpacegroup;
            m_poSpacegroup = new Cspacegroup(nSpacegroup);
        }
        
        // Try to read twin law stuff.
        nTemp = oHeader.nGetValue(ms_sCrystalXNumTwinLaws,1,&m_nTwinLaws);
        if (nTemp) {
            m_nTwinLaws = 0;
        } else if (m_nTwinLaws == 1) {
            nTemp = oHeader.nGetValue(ms_sCrystalXTwinLaws,m_nTwinLaws*3,&m_aa3fTwinLaw[0][0]);
            if (!nTemp)
                nTemp = oHeader.nGetValue(ms_sCrystalXTwinFraction,m_nTwinLaws,&m_afTwinFraction[0]);
            if (nTemp) {
                cout << "Ccrystal WARNING: reading " << ms_sCrystalXMosaicity
                    << " keyword!\n";
            };

		} else {
			cout << "Ccrystal WARNING:  There should be only 0 or 1 twin laws!\n";						
        }
        oHeader.nGetValue(ms_sCrystalXRecipShift,(m_nTwinLaws+1)*3,&m_aa3fRecipShift[0][0]);
		
    }
    
    if (0 != nStat)
    {
        // If there are errors, no sense have this around
        m_eThe_State = eCrystal_unknown_state;
        return (1);
    }
    else
    {
        m_eThe_State = eCrystal_notrefined_state;
        return (0);
    }
}


int Ccrystal::nFindNumTwins(Cimage_header& oHeader) {
    int nTemp;
    int nStat;
	int nx;
    Cstring sTemp;
	Cstring sPrefix;
    char pcBuf[10];
    double fCell[6];
    
    nTemp = 0; 
	m_nMaxTwinLaws = 0;
    do {
        nTemp++;
        sTemp = ms_sCrystalPrefix;
        if (nTemp!=1) {
            sprintf(pcBuf,"%d",nTemp);
            sTemp += pcBuf;
            sTemp += "_";
        };
		sPrefix = sTemp;
        sTemp += D_K_CrystalXUnitCell;
        nStat = oHeader.nGetValue(sTemp, 6, fCell);
		sTemp = sPrefix;
		sTemp += D_K_CrystalXNumTwinLaws;
		if (!oHeader.nGetValue(sTemp,1,&nx)) {
			m_nMaxTwinLaws = max(m_nMaxTwinLaws,nx);
		};


    } while (!nStat);
	m_nTwins = nTemp - 1;
	return 0;
};

int Ccrystal::nList(const int nFlag)
{
  printf("\nCrystal listing:\n\n");
  printf(" Unit cell lengths: %9.4f %9.4f %9.4f\n",
         m_fCell[0], m_fCell[1], m_fCell[2]);
  printf(" Unit cell  angles: %9.4f %9.4f %9.4f\n",
         m_fCell[3], m_fCell[4], m_fCell[5]);
  printf(" Unit cell  volume: %.3f\n", m_fVolume);
  if (2 > nFlag)
  {
      printf("Orientation angles: %9.4f %9.4f %9.4f\n",
             m_fOrientAngles[0], m_fOrientAngles[1], m_fOrientAngles[2]);
      printf("         Mosaicity:  %8.3f\n", fGetMosaicity());
      printf("       Description: %s\n", (const char *)m_sDescription);
  }
  
  fflush(stdout);
  
  m_poSpacegroup->nList(nFlag);
  
  printf("\n");
  
  return 0;
}

double
Ccrystal::fGetMosaicity() {
   return m_fMosaicity;
};

void
Ccrystal::vGetMosaicity(double* pfMosConst) {
    pfMosConst[0] = m_fMosaicity;
    pfMosConst[1] = m_fEffictiveSpectralDispersion;
    pfMosConst[2] = m_fMosaicityOffset;
};

void
Ccrystal::vGetMosaicity(float* pfMosConst) {
    pfMosConst[0] = m_fMosaicity;
    pfMosConst[1] = m_fEffictiveSpectralDispersion;
    pfMosConst[2] = m_fMosaicityOffset;
};


void
Ccrystal::vSetMosaicity(const double fMosaicity)
{
  m_fMosaicity = fMosaicity;
  m_fEffictiveSpectralDispersion = 0.0;
  m_fMosaicityOffset = 0.0;
  return;
}


void
Ccrystal::vSetMosaicity(const double fMosaicity,double fEffSpectralDispersion,double fMosaicityOffset)
{
  m_fMosaicity = fMosaicity;
  m_fEffictiveSpectralDispersion = fEffSpectralDispersion;
  m_fMosaicityOffset = fMosaicityOffset;
  return;
}


int
Ccrystal::nCalcGetRealMatrix(double* pfMatrix) {
    float a3x3fMatrix[3][3];
    int nStat;

    nStat = nCalcGetRealMatrix(&a3x3fMatrix[0][0]);
    vZeroMat3D(pfMatrix);
    if (nStat)
        return nStat;
    vCopyMat3D(&a3x3fMatrix[0][0],pfMatrix);
    return 0;
};

int
Ccrystal::nCalcGetRealMatrix(float* pfMatrix)
{
  // Calculate the real (direct cell) matrix from the real cell parameters
  //
  // D =  | a           0                                          0       |
  //      | b cos(gam)  b sin(gam)                                 0       |
  //      | c cos(bet)  c(cos(alp) - cos(bet)cos(gam))/sin(gam)    c *     |
  //                           sqrt(1 - [cos2(alp) + cos2(bet) + cos2(gam) |
  //                                     - 2(cos(alp)cos(bet)cos(gam))])   |
  //                                              / sin(gam)               |
  // Put real axis a along lab X.
  // See Prince, p. 18 and Int. Tables. Vol A, section 5.
  // WARNING WATCH OUT FOR ROW/COLUMN consistency!  Some methods may want
  // the transpose of this matrix.

  double fCosAlp, fCosBet, fCosGam, fSinGam, fTemp;

  fCosAlp     =  cos(m_fCell[3] * Gs_dRADIANS_PER_DEGREE);
  fCosBet     =  cos(m_fCell[4] * Gs_dRADIANS_PER_DEGREE);
  fCosGam     =  cos(m_fCell[5] * Gs_dRADIANS_PER_DEGREE);
  fSinGam     =  sin(m_fCell[5] * Gs_dRADIANS_PER_DEGREE);
  if (0.0 == fSinGam)
    return (-2);

  pfMatrix[0] = m_fCell[0];
  pfMatrix[1] = m_fCell[1] * fCosGam;
  pfMatrix[2] = m_fCell[2] * fCosBet;
  pfMatrix[3] = 0.0;
  pfMatrix[4] = m_fCell[1] * fSinGam;
  pfMatrix[5] = m_fCell[2] * (fCosAlp - fCosBet * fCosGam) / fSinGam;
  pfMatrix[6] = 0.0;
  pfMatrix[7] = 0.0;
  fTemp       = 1.0f - (  fCosAlp * fCosAlp
                        + fCosBet * fCosBet
                        + fCosGam * fCosGam
                        - 2.0f * fCosAlp * fCosBet * fCosGam);

  if (0.0f < fTemp)
    {
      pfMatrix[8] = m_fCell[2] * (double)sqrt(fTemp) / fSinGam;
      return (0);
    }
  else
    {
      // This is a bogus real cell

//      cout << "fTemp > 0 in nCalcGetRealCell: " << fTemp << endl << flush;
      pfMatrix[8] = 0.0f;
      return (-1);
    }
}

int Ccrystal::nCalcBMatrix(const double fWavelength)
{
  // Calculate the Busing and Levy (1967) B matrix:
  //
  // B =  | a*    b* cos(gam*)   c* cos(bet*)          |
  //      | 0     b* sin(gam*)  -c* sin(bet*) cos(alp) |
  //      | 0         0          c* sin(bet*) sin(alp) |
  //
  // -or- (more useful when taking derivatives of B wrt reciprocal cell params)
  // B =  | a*    b* cos(gam*)   c* cos(bet*)                                  |
  //      | 0     b* sin(gam*)  -c* [cos(bet*)cos(gam*) - cos(alp*)]/sin(gam*) |
  //      | 0         0          c* sqrt[1 - [cos(alp*)cos(alp*) +             |
  //                                          cos(bet*)cos(bet*) +
  //                                          cos(gam*)cos(gam*) -
  //                                    2 cos(alp*)cos(bet*)cos(gam*))]]
  //                                / sin(gam*)
  //
  // normalized by fWavelength

  int nStat;
  nStat = nCalcRecipCell();
  
  if (0 != nStat) return (nStat);  // Unit cell not available!

  m_fBMatrix[0][0] =  m_fRecipCell[0] * fWavelength;
  m_fBMatrix[0][1] =  0.0;
  m_fBMatrix[0][2] =  0.0;
  m_fBMatrix[1][0] =  m_fRecipCell[1] * cos(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE)
                                  * fWavelength;
  m_fBMatrix[1][1] =  m_fRecipCell[1] * sin(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE)
                                  * fWavelength;
  m_fBMatrix[1][2] =  0.0;
  m_fBMatrix[2][0] =  m_fRecipCell[2] * cos(m_fRecipCell[4] * Gs_dRADIANS_PER_DEGREE)
                                  * fWavelength;
  
  
  m_fBMatrix[2][1] = -m_fRecipCell[2] * sin(m_fRecipCell[4] * Gs_dRADIANS_PER_DEGREE)
                                  * cos(m_fCell[3] * Gs_dRADIANS_PER_DEGREE)
                                  * fWavelength;
  m_fBMatrix[2][2] =  m_fRecipCell[2] * sin(m_fRecipCell[4] * Gs_dRADIANS_PER_DEGREE)
                                  * sin(m_fCell[3] * Gs_dRADIANS_PER_DEGREE)
                                  * fWavelength;
  

  m_fBMatrix[2][1] = -m_fRecipCell[2] * fWavelength * (cos(m_fRecipCell[4] * Gs_dRADIANS_PER_DEGREE) * cos(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE) - cos(m_fRecipCell[3] * Gs_dRADIANS_PER_DEGREE)) / sin(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE);

  m_fBMatrix[2][2] =  m_fRecipCell[2] * fWavelength * sqrt(max(0.0,1.0 - (cos(m_fRecipCell[3] * Gs_dRADIANS_PER_DEGREE) * cos(m_fRecipCell[3] * Gs_dRADIANS_PER_DEGREE) +
                                                                     cos(m_fRecipCell[4] * Gs_dRADIANS_PER_DEGREE) * cos(m_fRecipCell[4] * Gs_dRADIANS_PER_DEGREE) +
                                                                     cos(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE) * cos(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE) - 
                                                            2 * cos(m_fRecipCell[3] * Gs_dRADIANS_PER_DEGREE) * cos(m_fRecipCell[4] * Gs_dRADIANS_PER_DEGREE) * cos(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE)))) / sin(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE); 
  
  return (0);
}

int Ccrystal::nCalcRecipCell(void)
{
  // Calculate reciprocal cell from real unit cell!
  // These are NOT normalized to any wavelength.

  double dTemp;
  double dCosAlp, dCosBet, dCosGam;
  double dSinAlp, dSinBet, dSinGam;

  dCosAlp = cos(m_fCell[3] * Gs_dRADIANS_PER_DEGREE);
  dCosBet = cos(m_fCell[4] * Gs_dRADIANS_PER_DEGREE);
  dCosGam = cos(m_fCell[5] * Gs_dRADIANS_PER_DEGREE);

  dSinAlp = sin(m_fCell[3] * Gs_dRADIANS_PER_DEGREE);  // maybe change to 1-cos**2?
  dSinBet = sin(m_fCell[4] * Gs_dRADIANS_PER_DEGREE);
  dSinGam = sin(m_fCell[5] * Gs_dRADIANS_PER_DEGREE);

  m_fVolume = fCalcVolume();
  if (   (0.0 >= m_fVolume)
      || (0.0 == dSinAlp)
      || (0.0 == dSinBet)
      || (0.0 == dSinGam) )
    return (2);

  m_fRecipCell[0] = m_fCell[1] * m_fCell[2] * dSinAlp / m_fVolume;
  m_fRecipCell[1] = m_fCell[0] * m_fCell[2] * dSinBet / m_fVolume;
  m_fRecipCell[2] = m_fCell[0] * m_fCell[1] * dSinGam / m_fVolume;

  dTemp = ((dCosBet * dCosGam) - dCosAlp) / (dSinBet * dSinGam);
  m_fRecipCell[3] = acos(dTemp) / Gs_dRADIANS_PER_DEGREE;

  dTemp = ((dCosAlp * dCosGam) - dCosBet) / (dSinAlp * dSinGam);
  m_fRecipCell[4] = acos(dTemp) / Gs_dRADIANS_PER_DEGREE;

  dTemp = ((dCosAlp * dCosBet) - dCosGam) / (dSinAlp * dSinBet);
  m_fRecipCell[5] = acos(dTemp) / Gs_dRADIANS_PER_DEGREE;

  m_fRecipVol = 1.f / m_fVolume;

  return (0);
}

int Ccrystal::nCalcOrientMatrix(const double fWavelength)
{
  int nStat;
  nStat = nCalcBMatrix(fWavelength);
  if (0 != nStat) return (nStat);
  nStat = nCalcRotMatrix();
  if (0 != nStat) return (nStat);

  vMulMat3DMat3D(m_fRotMatrix, m_fBMatrix, m_fOrientMatrix);

  return (0);
}

int Ccrystal::nCalcRotMatrix(void)
{
  // Convert the orient angles of the crystal (no goniostat involved!)
  // to a rotation matrix.

 vConv3Rot3Vec3DMat3D(m_fOrientAngles[0], m_fOrientAngles[1], m_fOrientAngles[2],
                       &m_fOrientVectors[0][0], &m_fOrientVectors[1][0],
                       &m_fOrientVectors[2][0], &m_fRotMatrix[0][0]);
  return (0);
}


void Ccrystal::vGetRecipCell(float *pfRecip)
{
  for (int i = 0; i < 6; i++) pfRecip[i] = m_fRecipCell[i];
}

void Ccrystal::vGetRecipCell(double *pfRecip)
{
  for (int i = 0; i < 6; i++) pfRecip[i] = m_fRecipCell[i];
}

void Ccrystal::vGetBMatrix(float *pfMat)
{
  double *pfTemp;
  pfTemp = &m_fBMatrix[0][0];
  for (int i = 0; i < 9; i++) pfMat[i] = *pfTemp++;
}

void Ccrystal::vGetBMatrix(double *pfMat)
{
  double *pfTemp;
  pfTemp = &m_fBMatrix[0][0];
  for (int i = 0; i < 9; i++) pfMat[i] = *pfTemp++;
}

void Ccrystal::vGetOrientMatrix(float *pfMat)
{
  double *pfTemp;
  pfTemp = &m_fOrientMatrix[0][0];
  for (int i = 0; i < 9; i++) pfMat[i] = *pfTemp++;
}

void Ccrystal::vGetOrientMatrix(double *pfMat)
{
  double *pfTemp;
  pfTemp = &m_fOrientMatrix[0][0];
  for (int i = 0; i < 9; i++) pfMat[i] = *pfTemp++;
}


int Ccrystal::nSetOrientMatrix(float *pfMat, const double fWavelength) {
    double a3x3fTempBuf[3][3];
    vCopyMat3D(pfMat,&a3x3fTempBuf[0][0]);
    return nSetOrientMatrix(&a3x3fTempBuf[0][0],fWavelength);
};

int Ccrystal::nSetOrientMatrix(double *pfMat, const double fWavelength)
{
  int nStat;
  double fAngle;
  double a3fTemp[3];

  double *pfTemp;
  pfTemp = &m_fOrientMatrix[0][0];
  for (int i = 0; i < 9; i++) *pfTemp++ = pfMat[i]/fWavelength;

  // Now build unit cell and rotation

  double a3x3fRealMat[3][3];
  m_fRecipVol = fInvMat3D(&m_fOrientMatrix[0][0],
                        &a3x3fRealMat[0][0]);

  vTranMat3D(a3x3fRealMat);
  nStat = 1;
  if (0.0f != m_fRecipVol)
    {
      // Get real cell from the Matrix
                
      m_fCell[0] = fLenVec3D(a3x3fRealMat[0]);
      m_fCell[1] = fLenVec3D(a3x3fRealMat[1]);
      m_fCell[2] = fLenVec3D(a3x3fRealMat[2]);

      fAngle  = fDot3D(a3x3fRealMat[1], a3x3fRealMat[2])
                / (m_fCell[1] * m_fCell[2]);
      fAngle = min(fAngle, 1.0f);
      fAngle = max(fAngle, -1.0f);
      m_fCell[3] = acos(fAngle) / Gs_dRADIANS_PER_DEGREE;

      fAngle  = fDot3D(a3x3fRealMat[0], a3x3fRealMat[2])
                / (m_fCell[0] * m_fCell[2]);
      fAngle = min(fAngle, 1.0f);
      fAngle = max(fAngle, -1.0f);
      m_fCell[4] = acos(fAngle) / Gs_dRADIANS_PER_DEGREE;
      fAngle  = fDot3D(a3x3fRealMat[0], a3x3fRealMat[1])
                / (m_fCell[0] * m_fCell[1]);
      fAngle = min(fAngle, 1.0f);
      fAngle = max(fAngle, -1.0f);
      m_fCell[5] = acos(fAngle) / Gs_dRADIANS_PER_DEGREE;

/*
      printf("\n: "
             " %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f\n",
             m_fCell[0], m_fCell[1], m_fCell[2],
             m_fCell[3], m_fCell[4], m_fCell[5]);
*/
      // Now determine orient angles

      double a3x3fInvBMat[3][3];
      (void) nCalcBMatrix();
      m_fRecipVol = fInvMat3D(&m_fBMatrix[0][0], &a3x3fInvBMat[0][0]);
      vMulMat3DMat3D(m_fOrientMatrix, a3x3fInvBMat, m_fRotMatrix);
                                
      // Compute determinant, should be 1.0 since it is a rotation matrix
      
//      vTranMat3D(m_fRotMatrix);
      fAngle = fInvMat3D( (double *)m_fRotMatrix, (double *)a3x3fInvBMat);

    
      if (fabs(1.0f - fAngle) > 0.001f)
        {
          cout << "ERROR: Determinant of new rotation matrix is: " << fAngle << '\n';
          return (1);
        }
    
      // vDeconvMat3D3XYZ absolutely NEEDS some starting angles

      a3fTemp[0] = m_fOrientAngles[0];
      a3fTemp[1] = m_fOrientAngles[1];
      a3fTemp[2] = m_fOrientAngles[2];
      
      vDeconvMat3D3XYZ(m_fRotMatrix, 
		       &a3fTemp[0],&a3fTemp[1],&a3fTemp[2]);
      m_fOrientAngles[0] = a3fTemp[0];
      m_fOrientAngles[1] = a3fTemp[1];
      m_fOrientAngles[2] = a3fTemp[2];
      nStat = 0;
/*
      cout << "Orient angles are: " 
           << m_fOrientAngles[0] << ", " << m_fOrientAngles[1] << ", "
           << m_fOrientAngles[2] << endl;
*/
    }

  return (nStat);
}

void Ccrystal::vGetRotMatrix(float *pfMat)
{
  double *pfTemp;
  pfTemp = &m_fRotMatrix[0][0];
  for (int i = 0; i < 9; i++) pfMat[i] = *pfTemp++;
}
void Ccrystal::vGetRotMatrix(double *pfMat)
{
  double *pfTemp;
  pfTemp = &m_fRotMatrix[0][0];
  for (int i = 0; i < 9; i++) pfMat[i] = *pfTemp++;
}


double Ccrystal::fCalcVolume()
{
  double fTemp;
  double  fCosAlp, fCosBet, fCosGam;

  fCosAlp = cos(m_fCell[3] * Gs_dRADIANS_PER_DEGREE);
  fCosBet = cos(m_fCell[4] * Gs_dRADIANS_PER_DEGREE);
  fCosGam = cos(m_fCell[5] * Gs_dRADIANS_PER_DEGREE);

  fTemp         = 1.0f -  fCosAlp * fCosAlp
                      -  fCosBet * fCosBet
                      -  fCosGam * fCosGam
                      +  2.0f *fCosAlp * fCosBet * fCosGam;

  if (fTemp <= 0.0) return (0.0);

  m_fVolume       = m_fCell[0] * m_fCell[1] * m_fCell[2] * sqrt(fTemp);

  return (m_fVolume);
}

// Function figures out how many cell parameters should be refined.  It returns
// the (normalized) refinement vectors in the fNormDirVecs array.


int Ccrystal::nGetDerivVecs(int& nNumberRefinedCellParameters,float fNormDirVecs[6][6]) {
    int nx;

    // Now rearrange matrices according to the crystal system if required
  // and report the number of refinable cell parameters

   vZeroMat(6,6,&fNormDirVecs[0][0]);

    if (m_poSpacegroup->eGetClass() == eSpacegroup_triclinic_class)
    {
      nNumberRefinedCellParameters = 6;  // a*, b*, c*, alp*, bet*, gam*
      for (nx=0; nx<6;nx++)
          fNormDirVecs[nx][nx]=1.0;

    }
  else if (m_poSpacegroup->eGetClass() == eSpacegroup_monoclinic_class)
    {
      nNumberRefinedCellParameters = 4;    // a*, b*, c*, bet*
      fNormDirVecs[0][0]=1.0;
      fNormDirVecs[1][1]=1.0;
      fNormDirVecs[2][2]=1.0;
      fNormDirVecs[3][4]=1.0; 
    }
  else if (m_poSpacegroup->eGetClass() == eSpacegroup_orthorhombic_class)
    {
      nNumberRefinedCellParameters = 3;   // a*, b*, c*
      fNormDirVecs[0][0]=1.0;
      fNormDirVecs[1][1]=1.0;
      fNormDirVecs[2][2]=1.0;
    }
  else if (   (m_poSpacegroup->eGetClass() == eSpacegroup_tetragonal_class)
      || (m_poSpacegroup->eGetClass() == eSpacegroup_trigonal_class)
	   || (m_poSpacegroup->eGetClass() == eSpacegroup_hexagonal_class) )
    {
      // If tetragonal, trigonal or hexagonal
      nNumberRefinedCellParameters = 2;  // (a*,c*), b*

      fNormDirVecs[0][0]=1.0;
      fNormDirVecs[0][1]=1.0;
      fNormDirVecs[1][2]=1.0;
    }
  else if (m_poSpacegroup->eGetClass() == eSpacegroup_cubic_class)
    {
      // If cubic
      nNumberRefinedCellParameters = 1;  // (a*,b*,c*)

      fNormDirVecs[0][0]=1.0;
      fNormDirVecs[0][1]=1.0;
      fNormDirVecs[0][2]=1.0;
    }
  else if (m_poSpacegroup->eGetClass() == eSpacegroup_rhombohedral_class)
    {
      // If rhombohedral
      nNumberRefinedCellParameters = 2;  // (a*,b*,c*) (alp*,bet*,gam*)
      fNormDirVecs[0][0]=1.0;
      fNormDirVecs[0][1]=1.0;
      fNormDirVecs[0][2]=1.0;
      fNormDirVecs[1][3]=1.0;
      fNormDirVecs[1][4]=1.0;
      fNormDirVecs[1][5]=1.0;
    }
  return (0);
};

int Ccrystal::nCalcGetBDeriv(float fMat[7][3][3],
			     int *pnNumberRefinedCellParameters,
                 const float fWavelength) {
    double fMat_[7][3][3];
    int nStat;
    nStat = nCalcGetBDeriv(fMat_,pnNumberRefinedCellParameters,fWavelength);
    vCopyVecND(7*3*3,&fMat_[0][0][0],&fMat[0][0][0]);
    return nStat;
};


int Ccrystal::nCalcGetBDeriv(double fMat[7][3][3],
			     int *pnNumberRefinedCellParameters,
			     const double fWavelength)

{
  // Calculate and return in the argument list the following:
  // fMat[0][][]  -  The B matrix (see nCalcBMatrix)
  // fMat[1][][]  -  Derivative of B with respect to a*
  // fMat[2][][]  -  Derivative of B with respect to b*
  // fMat[3][][]  -  Derivative of B with respect to c*
  // fMat[4][][]  -  Derivative of B with respect to alpha*
  // fMat[5][][]  -  Derivative of B with respect to beta*
  // fMat[6][][]  -  Derivative of B with respect to gamma*

  int nStat;         // General status flag
  int i, j;          // Loop counters
  double  fTemp;
  double  fCosAlpStar, fCosBetStar, fCosGamStar;
  double  fSinAlpStar, fSinBetStar, fSinGamStar;

  // Calculate the B Matrix and reciprocal cell values

  nStat = nCalcBMatrix(fWavelength);
  if (0 != nStat) return (nStat);

  // Get the B matrix into fMat[0][][]

  vGetBMatrix(&fMat[0][0][0]);

  // Compute some sines and cosines for later use

  fCosAlpStar = cos(m_fRecipCell[3] * Gs_dRADIANS_PER_DEGREE);
  fCosBetStar = cos(m_fRecipCell[4] * Gs_dRADIANS_PER_DEGREE);
  fCosGamStar = cos(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE);

  fSinAlpStar = sin(m_fRecipCell[3] * Gs_dRADIANS_PER_DEGREE);  // Change 1-cos**2??
  fSinBetStar = sin(m_fRecipCell[4] * Gs_dRADIANS_PER_DEGREE);
  fSinGamStar = sin(m_fRecipCell[5] * Gs_dRADIANS_PER_DEGREE);

  // Compute derivative of B with respect to a*

  vZeroMat3D(&fMat[1][0][0]);
  fMat[1][0][0] = fMat[0][0][0] / (m_fRecipCell[0] );

  // Compute derivative of B with respect to b*

  vZeroMat3D(&fMat[2][0][0]);
  fMat[2][1][0] = fMat[0][1][0] / (m_fRecipCell[1] );
  fMat[2][1][1] = fMat[0][1][1] / (m_fRecipCell[1] );

  // Compute derivative of B with respect to c*

  vZeroMat3D(&fMat[3][0][0]);
  fMat[3][2][0] = fMat[0][2][0] / (m_fRecipCell[2] );
  fMat[3][2][1] = fMat[0][2][1] / (m_fRecipCell[2] );
  fMat[3][2][2] = fMat[0][2][2] / (m_fRecipCell[2] );

  // Computer derivative of B with respect to alpha*

// Does fTemp = m_fRecipVolume?
  fTemp = (double)sqrt(1.0f - (fCosAlpStar * fCosAlpStar +
                      fCosBetStar * fCosBetStar +
                      fCosGamStar * fCosGamStar -
                      2.0f * fCosAlpStar * fCosBetStar * fCosGamStar));

  vZeroMat3D(&fMat[4][0][0]);
  fMat[4][2][1] = -(m_fRecipCell[2] * fWavelength) / fSinGamStar * fSinAlpStar;
  fMat[4][2][2] =  /* - (This negative was removed by tjn on 08/22/00 */ 
                  fMat[0][2][1] / fTemp * fSinAlpStar;
  //            = (m_fRecipCell[2] * fWavelength) / fSinGamStar  /
  //                     fTemp * fSinAlpStar *
  //                     (fCosAlpStar - (fCosBetStar * fCosGamStar));

  // Computer derivative of B with respect to beta*

  vZeroMat3D(&fMat[5][0][0]);
  fMat[5][2][0] = -(m_fRecipCell[2] * fWavelength ) * fSinBetStar;
  fMat[5][2][1] =  (m_fRecipCell[2] * fWavelength ) * fSinBetStar
                                                 * fCosGamStar / fSinGamStar;

  fMat[5][2][2] =  (m_fRecipCell[2] * fWavelength ) / fSinGamStar / fTemp
                                                 * fSinBetStar *
                           (fCosBetStar - (fCosAlpStar * fCosGamStar));

  // Computer derivative of B with respect to gamma*

  vZeroMat3D(&fMat[6][0][0]);
  fMat[6][1][0] = -fMat[0][1][1];  // d(b*cos(gam*))/dgam* = -b*sin(gam*)
  fMat[6][1][1] =  fMat[0][1][0];  // d(b*sin(gam*))/dgam* = b*cos(gam*)
  fMat[6][2][1] =   m_fRecipCell[2] * fWavelength  * (fCosBetStar /
                    (fSinGamStar * fSinGamStar) - fCosAlpStar *
                     fCosGamStar / (fSinGamStar * fSinGamStar));

  fMat[6][2][2] = 
      m_fRecipCell[2] * fWavelength *(
      ( fCosGamStar - fCosAlpStar*fCosBetStar)/fTemp - fTemp*fCosGamStar/fSinGamStar/fSinGamStar);

  int nx;
  for (nx=4;nx<=6;nx++) {
      vMulVecNDScalar(9,&fMat[nx][0][0],Gs_dRADIANS_PER_DEGREE,&fMat[nx][0][0]);
  };
            
      
/*  This was the previous calculation for this derivative.  It is almost correct I think except for a sign error.
    Changed on 08/22/00
      -m_fRecipCell[2]  / fSinGamStar *
                  (fTemp * fCosGamStar / fSinGamStar
                   + (fSinGamStar * fCosGamStar - fCosAlpStar * fCosBetStar *
                   fSinGamStar) / fTemp);

  */

  // Remember:
  //  (more useful when taking derivatives of B wrt reciprocal cell params)
  // B =  | a*    b* cos(gam*)   c* cos(bet*)                                  |
  //      | 0     b* sin(gam*)  -c* [cos(bet*)cos(gam*) - cos(alp*)]/sin(gam*) |
  //      | 0         0          c* sqrt[1 - [cos(alp*)cos(alp*) +             |
  //                                          cos(bet*)cos(bet*) +
  //                                          cos(gam*)cos(gam*) -
  //                                    2 cos(alp*)cos(bet*)cos(gam*))]
  //                                / sin(gam*)
  //

  // Now rearrange matrices according to the crystal system if required
  // and report the number of refinable cell parameters

  if (m_poSpacegroup->eGetClass() == eSpacegroup_triclinic_class)
    {
      *pnNumberRefinedCellParameters = 6;  // a*, b*, c*, alp*, bet*, gam*
    }
  else if (m_poSpacegroup->eGetClass() == eSpacegroup_monoclinic_class)
    {
      // If Monoclinic, check sine of beta* and gamma* to determine setting
      int nOtherNinety = 5;                     // Assume beta* is not 90.
      if (fSinGamStar < 1.0) nOtherNinety = 6;  // Oops, it was gam* not 90.

      for (i = 0; i < 3; i++)
	{
	  for (j = 0; j < 3; j++)
	    {
	      fMat[4][i][j] = fMat[4][i][j] + fMat[nOtherNinety][i][j];
	    }
	}
      *pnNumberRefinedCellParameters = 4;  // a*, b*, c*, bet*
    }
  else if (m_poSpacegroup->eGetClass() == eSpacegroup_orthorhombic_class)
    {
      *pnNumberRefinedCellParameters = 3;  // a*, b*, c*
    }
  else if (   (m_poSpacegroup->eGetClass() == eSpacegroup_tetragonal_class)
           || (m_poSpacegroup->eGetClass() == eSpacegroup_trigonal_class)
	   || (m_poSpacegroup->eGetClass() == eSpacegroup_hexagonal_class) )
    {
      // If tetragonal, trigonal or hexagonal

      for (i = 0; i < 3; i++)
	{
	  for (j = 0; j < 3; j++)
	    {
	      fMat[1][i][j] = fMat[1][i][j] + fMat[2][i][j];
	      fMat[2][i][j] = fMat[3][i][j];
	    }
	}
      *pnNumberRefinedCellParameters = 2;  // a*, c*
    }
  else if (m_poSpacegroup->eGetClass() == eSpacegroup_cubic_class)
    {
      // If cubic, there is only 1 refinable parameter, that is a*
      for (i = 0; i < 3; i++)
	{
	  for (j = 0; j < 3; j++)
	    {
	      fMat[1][i][j] = fMat[1][i][j] + fMat[2][i][j] + fMat[3][i][j];
	    }
	}
      *pnNumberRefinedCellParameters = 1;  // a*
    }
  else if (m_poSpacegroup->eGetClass() == eSpacegroup_rhombohedral_class)
    {
      // If tetragonal, trigonal or hexagonal

      for (i = 0; i < 3; i++)
	{
	  for (j = 0; j < 3; j++)
	    {
	      fMat[1][i][j] = fMat[1][i][j] + fMat[2][i][j] + fMat[3][i][j];
	      fMat[2][i][j] = fMat[4][i][j] + fMat[5][i][j] + fMat[6][i][j];
	    }
	}
      *pnNumberRefinedCellParameters = 2;  // a*, alp*
    }
  return (0);
}


int Ccrystal::nCalcGetRotDeriv(double fMat[3][3][3]) {
    double a3x3fRotMatTemp[3][3];

     vConv3Rot3Vec3DMat3D(
               m_fOrientAngles[0], m_fOrientAngles[1], m_fOrientAngles[2],
		       &m_fOrientVectors[0][0], &m_fOrientVectors[1][0],
		       &m_fOrientVectors[2][0], &a3x3fRotMatTemp[0][0],
               &fMat[0][0][0],
               &fMat[1][0][0],
               &fMat[2][0][0]               
               );

     return 0;
};

int Ccrystal::nCalcGetRealCell(const float *pfRecipCell,
                               float *pfRealCell,
                               const double fWavelength,
                               const float *pfRecipCellSig,
                               float *pfRealCellSig)

{
  // Calculate the real space unit cell from reciprocal space constants
  // If sigmas on the recip cell are given, then calculate sigmas for the
  // real cell.
  // Input reciprocal cell angles are in degrees.  The lengths are
  // unitless and must be divided by the wavelength

  int i;                 // Loop counter

  double  fCosAlpStar, fCosBetStar, fCosGamStar;
  double  fCos[3];
  double  fSin[3];
  double  fSinAlpStar, fSinBetStar, fSinGamStar;

  // Compute some sines and cosines for later use

  fCos[0] = fCosAlpStar = cos((double)pfRecipCell[3] * Gs_dRADIANS_PER_DEGREE);
  fCos[1] = fCosBetStar = cos((double)pfRecipCell[4] * Gs_dRADIANS_PER_DEGREE);
  fCos[2] = fCosGamStar = cos((double)pfRecipCell[5] * Gs_dRADIANS_PER_DEGREE);

  fSin[0] = fSinAlpStar = sqrt(1.0f - fCosAlpStar * fCosAlpStar);
  fSin[1] = fSinBetStar = sqrt(1.0f - fCosBetStar * fCosBetStar);
  fSin[2] = fSinGamStar = sqrt(1.0f - fCosGamStar * fCosGamStar);

  double fTemp =
    1.0f - fCosAlpStar * fCosAlpStar - fCosBetStar * fCosBetStar
        - fCosGamStar * fCosGamStar
        + 2.0f * fCosAlpStar * fCosBetStar * fCosGamStar;

  if (0.0 >= fTemp)
    {
      cout << "Error in Ccrystal::nCalcGetRealCell, check recip cell angles!\n";
      return (1);
    }

  fTemp        = sqrt(fTemp);

  // Calculate real cell lengths

  pfRealCell[0] = fWavelength * fSinAlpStar / (pfRecipCell[0] * fTemp);
  pfRealCell[1] = fWavelength * fSinBetStar / (pfRecipCell[1] * fTemp);
  pfRealCell[2] = fWavelength * fSinGamStar / (pfRecipCell[2] * fTemp);

  // Calculation of direct cell angles in degrees

  pfRealCell[3] = (float)acos((double) (fCosBetStar * fCosGamStar - fCosAlpStar) /
                          (fSinBetStar * fSinGamStar)) / Gs_dRADIANS_PER_DEGREE;
  pfRealCell[4] = (float)acos((double) (fCosAlpStar * fCosGamStar - fCosBetStar) /
                          (fSinAlpStar * fSinGamStar)) / Gs_dRADIANS_PER_DEGREE;
  pfRealCell[5] = (float)acos((double) (fCosAlpStar * fCosBetStar - fCosGamStar) /
                          (fSinAlpStar * fSinBetStar)) / Gs_dRADIANS_PER_DEGREE;

  // Calculate standard deviations if sigmas are not NULL

  if ( (pfRecipCellSig != NULL) && (pfRealCellSig != NULL) ) {
    double fTemp3 = fTemp * fTemp * fTemp;
    int j, k;

    // Standard deviations of direct length constants

    double fT1, fT2, fT3, fT4;
    for (i = 0; i < 3; i++) {
      j = (i+1) % 3;             // 1, 2, 0
      k = (i+2) % 3;             // 2, 0, 1

      fT1 = 1.0f / (pfRecipCell[i] * pfRecipCell[i] * fTemp)
             * fSin[i] * pfRecipCellSig[i];

      fT2 = 1.0f / pfRecipCell[i] * pfRecipCellSig[i]
             * (fCos[i] / fTemp  +  fSin[i] / fTemp
             *(-fCos[i] * fSin[i]  +  fSin[i] * fCos[j] * fCos[k]));

      fT3 = 1.0f / (pfRecipCell[i] * fTemp) * fSin[i] * pfRecipCellSig[j]
            *(-fCos[j] * fSin[j]  +  fCos[i] * fSin[j] * fCos[k]);

      fT4 = 1.0f / (pfRecipCell[i] * fTemp) * fSin[i] * pfRecipCellSig[k]
            *(-fCos[k] * fSin[k]  +  fCos[i] * fCos[j] * fSin[k]);


      pfRealCellSig[i] =  sqrt(fT1 * fT1  +  fT2 * fT2  +
                                fT3 * fT3  +  fT4 * fT4);
     }

  //  Standard deviations of direct cell angles

    for (i = 0; i < 3; i++) {
      j = (i+1) % 3;             // 1, 2, 0
      k = (i+2) % 3;             // 2, 0, 1

      fTemp3 = (fCos[j] * fCos[k]  -  fCos[i]) / (fSin[j] * fSin[k]);

      fTemp3 = -1.0f / sqrt(1.0f - (fTemp3 * fTemp3));

      fT1 = fTemp3 * (fSin[i] / (fSin[j] * fSin[k]))
             * pfRecipCellSig[i+3] * Gs_dRADIANS_PER_DEGREE;

      fT2 = fTemp3 * (-fCos[k] / (fSin[j] * fSin[j] * fSin[k])
             + fCos[i] * fCos[j] / (fSin[j]* fSin[j] * fSin[k]))
             * pfRecipCellSig[j+3] * Gs_dRADIANS_PER_DEGREE;

      fT3 = fTemp3 * (-fCos[j] / (fSin[j] * fSin[k] * fSin[k])
             + fCos[i] * fCos[k] / (fSin[j] * fSin[k] * fSin[k]))
             * pfRecipCellSig[k+3] * Gs_dRADIANS_PER_DEGREE;

      pfRealCellSig[i+3] = sqrt(fT1 * fT1  +  fT2 * fT2  +  fT3 * fT3)
                                 / Gs_dRADIANS_PER_DEGREE;
    }

  }

  return (0);
}

int Ccrystal::nUpdateHeader(Cimage_header* poHeader)
{

  // Update the header with the current values in this crystal

  int nStat;
  int nx;
  double a3fMosCoeffs[3];

  a3fMosCoeffs[0] = m_fMosaicity;
  a3fMosCoeffs[1] = m_fEffictiveSpectralDispersion;
  a3fMosCoeffs[2] = m_fMosaicityOffset;

  if (m_nTwin<=0) {
	  nFindNumTwins(*poHeader);
	  m_nTwin = m_nTwins + 1;
  };
  nLoadKeywords(m_nTwin);
  nStat =         poHeader->nReplaceValue(ms_sCrystalXUnitCell, 6, m_fCell, 4);
  nStat = nStat + poHeader->nReplaceValue(ms_sCrystalXMosaicity, 3,&a3fMosCoeffs[0], 4);
  nStat = nStat + poHeader->nReplaceValue(ms_sCrystalXOrientAngles, 3, m_fOrientAngles, 4);
  nStat = nStat + poHeader->nReplaceValue(ms_sCrystalXSpacegroup, m_poSpacegroup->nGet());
  nStat = nStat + poHeader->nReplaceValue(ms_sCrystalDescription, m_sDescription);
  if (m_nTwinLaws) {
	  nStat = nStat + poHeader->nReplaceValue(ms_sCrystalXNumTwinLaws,m_nTwinLaws);
	  nStat = nStat + poHeader->nReplaceValue(ms_sCrystalXTwinLaws,m_nTwinLaws*3,&m_aa3fTwinLaw[0][0]);
	  nStat = nStat + poHeader->nReplaceValue(ms_sCrystalXTwinFraction,m_nTwinLaws,&m_afTwinFraction[0]);  
  } else {
	  poHeader->nDelete(ms_sCrystalXNumTwinLaws);
  };
  for (nx = 0; nx < (m_nTwinLaws+1)*3;nx++) {
	  if (*(&(m_aa3fRecipShift[0][0]) + nx) != 0.0) {
		  // Only write these guys out if you really need to.  No need to clutter up the header file otherwise
  		nStat = nStat + poHeader->nReplaceValue(ms_sCrystalXRecipShift,(m_nTwinLaws+1)*3,&(m_aa3fRecipShift[0][0]));
		// And break out since we just wrote them.
		break;
	  };
  };


  return (nStat);
}

int Ccrystal::nUpdateHeaderSigmaValues(Cimage_header* poHeader,
									   float* pfCellSigmas,
									   float* pfRotSigmas,
									   float* pfMosaicitySigmas,
									   float* pfRecipShiftSigmas,
									   float* pfTwinLaw1Sigmas,
									   float* pfTwinFraction1Sigmas) {
    // Place sigmas in header.  This should be encapsulated in the classes
    // and not here!!!
    
	int nx;
    Cstring sSigmas = "_SIGMA";
	nLoadKeywords(m_nTwin);
    
    poHeader->nReplaceValue(ms_sCrystalXUnitCell + sSigmas, 6,
        pfCellSigmas, 5);
    
    poHeader->nReplaceValue(ms_sCrystalXOrientAngles + sSigmas, 3,
        pfRotSigmas, 4);
	
	poHeader->nReplaceValue(ms_sCrystalXMosaicity + sSigmas,3,
        pfMosaicitySigmas, 4);

	
	if (m_nTwinLaws) {
		poHeader->nReplaceValue(ms_sCrystalXTwinLaws + sSigmas,3*m_nTwinLaws,
			pfTwinLaw1Sigmas, 4);
		poHeader->nReplaceValue(ms_sCrystalXTwinFraction + sSigmas,m_nTwinLaws,
			pfTwinFraction1Sigmas, 4);
	};


	for (nx = 0; nx < (m_nTwinLaws+1)*3;nx++) {
		if (*(&(m_aa3fRecipShift[0][0]) + nx) != 0.0) {
			// Only write these guys out if you really need to.  No need to clutter up the header file otherwise
			poHeader->nReplaceValue(ms_sCrystalXRecipShift + sSigmas,(m_nTwinLaws+1)*3,pfRecipShiftSigmas);
			// And break out since we just wrote them.
			break;
		};
	};

	
	return 0;
};

	    



int Ccrystal::nLoadKeywords(int nTwin) {
    Cstring sTemp;
    char pcBuf[10];
    
    sTemp = ms_sCrystalPrefix;
    if (nTwin!=1) {
        sprintf(pcBuf,"%d",nTwin);
        sTemp += pcBuf;
        sTemp += "_";
    };
    
    ms_sCrystalXUnitCell = sTemp;
    ms_sCrystalXMosaicity = sTemp;
    ms_sCrystalXOrientAngles = sTemp;
    ms_sCrystalXOrientVectors = sTemp;
    ms_sCrystalXSpacegroup = sTemp;
	ms_sCrystalXNumTwinLaws  = sTemp;
	ms_sCrystalXTwinLaws = sTemp;
	ms_sCrystalXTwinFraction = sTemp;
	ms_sCrystalXRecipShift = sTemp;
    
    ms_sCrystalXUnitCell += D_K_CrystalXUnitCell;      
    ms_sCrystalXMosaicity += D_K_CrystalXMosaicity;
    ms_sCrystalXOrientAngles += D_K_CrystalXOrientAngles;
    ms_sCrystalXOrientVectors += D_K_CrystalXOrientVectors;
    ms_sCrystalXSpacegroup  += D_K_CrystalXSpacegroup;   
	ms_sCrystalXNumTwinLaws  += D_K_CrystalXNumTwinLaws;
	ms_sCrystalXTwinLaws += D_K_CrystalXTwinLaws;
	ms_sCrystalXTwinFraction += D_K_CrystalXTwinFraction;
	ms_sCrystalXRecipShift   += D_K_CrystalXRecipShift;

    return 0;
};

void
Ccrystal::vGetCell(float *pfCell)
{
  int i;
  for (i = 0; i < 6; i++)
    pfCell[i] = m_fCell[i];
}
void
Ccrystal::vGetCell(double *pfCell)
{
  int i;
  for (i = 0; i < 6; i++)
    pfCell[i] = m_fCell[i];
}


void
Ccrystal::vGetCell(float *pfA, float *pfB, float *pfC,
		   float *pfAlpha, float *pfBeta, float *pfGamma)
{
  *pfA     = m_fCell[0];
  *pfB     = m_fCell[1];
  *pfC     = m_fCell[2];
  *pfAlpha = m_fCell[3];
  *pfBeta  = m_fCell[4];
  *pfGamma = m_fCell[5];
}

void
Ccrystal::vSetCell(const float *pfCell)
{
  int i;
  for (i = 0; i < 6; i++)
    m_fCell[i] = pfCell[i];
  (void) nCalcRecipCell();  // Do this to assure consistency
  // Should we calculate Bmatrix, Rotation matrix, orientation matrix?
}

void
Ccrystal::vSetCell(const double *pfCell)
{
  int i;
  for (i = 0; i < 6; i++)
    m_fCell[i] = pfCell[i];
  (void) nCalcRecipCell();  // Do this to assure consistency
  // Should we calculate Bmatrix, Rotation matrix, orientation matrix?
}

void
Ccrystal::vSetRecipCell(const float* pfRecipCell)
{
  int i;
  float f0;
  for (i = 0; i < 6; i++)
    m_fCell[i] = pfRecipCell[i];
  (void) nCalcRecipCell();  // Since the transformations changing reciprocal cell to real cell are exaclty the
                            // same as to change from real cell to reciprocal cell, we will use this function,
                            // and then swap the results.

  f0 = m_fRecipVol;
  m_fRecipVol = m_fVolume;
  m_fVolume = f0;
  for (i=0;i<6;i++) {
      f0 = m_fCell[i];
      m_fCell[i] = m_fRecipCell[i];
      m_fRecipCell[i] = f0;
  };
  // Should we calculate Bmatrix, Rotation matrix, orientation matrix?
}

void
Ccrystal::vSetRecipCell(const double* pfRecipCell)
{
  int i;
  double f0;
  for (i = 0; i < 6; i++)
    m_fCell[i] = pfRecipCell[i];
  (void) nCalcRecipCell();  // Since the transformations changing reciprocal cell to real cell are exaclty the
                            // same as to change from real cell to reciprocal cell, we will use this function,
                            // and then swap the results.

  f0 = m_fRecipVol;
  m_fRecipVol = m_fVolume;
  m_fVolume = f0;
  for (i=0;i<6;i++) {
      f0 = m_fCell[i];
      m_fCell[i] = m_fRecipCell[i];
      m_fRecipCell[i] = f0;
  };
  // Should we calculate Bmatrix, Rotation matrix, orientation matrix?
}

void
Ccrystal::vSetCell(const double fA, const double fB, const double fC,
		   const double fAlpha, const double fBeta, const double fGamma)
{
  m_fCell[0] = fA;
  m_fCell[1] = fB;
  m_fCell[2] = fC;
  m_fCell[3] = fAlpha;
  m_fCell[4] = fBeta;
  m_fCell[5] = fGamma;
  (void) nCalcRecipCell();  // Do this to assure consistency Wavelength == 1!
  // Should we calculate Bmatrix, Rotation matrix, orientation matrix?
  // Probably should give a warning if the spacegroup does not match!
}

void
Ccrystal::vGetOrientAngles(float *pfOrientAngles)
{
  int i;
  for (i = 0; i < 3; i++)
    pfOrientAngles[i] = m_fOrientAngles[i];
}

void
Ccrystal::vGetOrientAngles(double *pfOrientAngles)
{
  int i;
  for (i = 0; i < 3; i++)
    pfOrientAngles[i] = m_fOrientAngles[i];
}

void
Ccrystal::vGetOrientAngles(float *pfOr1, float *pfOr2, float *pfOr3)
{
    *pfOr1 = m_fOrientAngles[0];
    *pfOr2 = m_fOrientAngles[1];
    *pfOr3 = m_fOrientAngles[2];
}
void
Ccrystal::vGetOrientAngles(double *pfOr1, double *pfOr2, double *pfOr3)
{
    *pfOr1 = m_fOrientAngles[0];
    *pfOr2 = m_fOrientAngles[1];
    *pfOr3 = m_fOrientAngles[2];
}

void
Ccrystal::vSetOrientAngles(const float *pfOrientAngles)
{
  int i;
  for (i = 0; i < 3; i++)
    m_fOrientAngles[i] = pfOrientAngles[i];
  (void) nCalcRotMatrix();
}

void
Ccrystal::vSetOrientAngles(const double *pfOrientAngles)
{
  int i;
  for (i = 0; i < 3; i++)
    m_fOrientAngles[i] = pfOrientAngles[i];
  (void) nCalcRotMatrix();
}

void
Ccrystal::vSetOrientAngles(const double fOr1, const double fOr2, const double fOr3)
{
  m_fOrientAngles[0] = fOr1;
  m_fOrientAngles[1] = fOr2;
  m_fOrientAngles[2] = fOr3;

  (void) nCalcRotMatrix();
// Probably should recalculate other matrixes at this point too!
}





double Ccrystal::fCalcTransForOrientMatrix(Ccrystal& oCrystalOther,double* pfTransMat) {
    double a3x3fCellOther[3][3];
    oCrystalOther.vGetOrientMatrix(&a3x3fCellOther[0][0]);
    return fCalcTransForOrientMatrix(a3x3fCellOther,pfTransMat);
};

double Ccrystal::fCalcTransForOrientMatrix(double a3x3fOrientMatrixOther[3][3],double* pfTransMat) {
    double a3x3fOrientMatrix[3][3];
    double a3x3fInvOrientMatrix[3][3];
    double a3x3fTransMatrix[3][3];
    double a3x3fRoundedTransMatrix[3][3];
        
    double afValidRoundoff[] = {3.0,-3.0,2.0,-2.0,1.0,-1.0,-0.5,0.5,-0.33333,0.33333,2,-2,3,-3,0.0,10.0};

    int nx,ny,nz;
    double f0,f1;


    vGetOrientMatrix(&a3x3fOrientMatrix[0][0]);
    fInvMat3D(&a3x3fOrientMatrix[0][0],&a3x3fInvOrientMatrix[0][0]);
    vMulMat3DMat3D(a3x3fInvOrientMatrix,a3x3fOrientMatrixOther,a3x3fTransMatrix);

    for (nx=0;nx<3;nx++) {
        for (ny=0;ny<3;ny++) {
            f0 = 0.1;
            a3x3fRoundedTransMatrix[nx][ny] = 10.0;
            for (nz = 0;afValidRoundoff[nz]!=10.0;nz++) {
                f1 = ABS(a3x3fTransMatrix[nx][ny] - afValidRoundoff[nz]);
                if (f1<f0) {
                    f0 = f1;
                    a3x3fRoundedTransMatrix[nx][ny] = afValidRoundoff[nz];
                };
            };
            if (a3x3fRoundedTransMatrix[nx][ny]==10.0)
                return 10.0;
        };
    };
    f0 = ABS(fDetMat3D(&a3x3fRoundedTransMatrix[0][0]));
    if (pfTransMat) {
        vCopyMat3D(&a3x3fRoundedTransMatrix[0][0],pfTransMat);
    };
    return (ABS(f0-1.0));
};



int Ccrystal::nCountUnique(int nNumReso,                                
                                double* pfResoLower,
                                double* pfResoUpper,
                                int*    pnResoUnique,        
                                const int nAnomFlag)
{
    
    int nStat;
    int nx;
    double f0;
    int h, k, l;

    unsigned char* pcCovered;
    Crefln  *poRefln;
    float    a3x3fBMat[3][3];      // Local var for crystal B matrix
    float    a3fX[3];              // Local var for recip lattice vector Bh = x
    float    a3fH[3];              // hkl of refln
    float    fResoMax, fResoMin;   // Min and max resolution in Angstrom units
    float    fDstarSq;
    int      nHmin, nHmax, nKmin, nKmax, nLmin, nLmax;
       
    // Get the crystal B matrix into local variable
    
    nCalcBMatrix();
    vGetBMatrix(&a3x3fBMat[0][0]);
        
    
    
    
    
    // Get min and max resolution no matter which order the calling program
    // passed them in.

    for (nx=0;nx<nNumReso;nx++) {    
        if (nx==0) {
            fResoMax = min(pfResoLower[nx],pfResoUpper[nx]);
            fResoMin = max(pfResoLower[nx],pfResoUpper[nx]);
        } else {
            fResoMax = min(fResoMax,min(pfResoLower[nx],pfResoUpper[nx]));
            fResoMin = max(fResoMin,min(pfResoLower[nx],pfResoUpper[nx]));
        };
        if (pfResoUpper[nx]>pfResoLower[nx])
            std::swap(pfResoUpper[nx],pfResoLower[nx]);
        pnResoUnique[nx] = 0;
    };
    
    // Get possible min, max hkl for the reflection range
    
    nHmax = (int) ( fGetCell(0) / fResoMax) + 1;
    nKmax = (int) ( fGetCell(1) / fResoMax) + 1;
    nLmax = (int) ( fGetCell(2) / fResoMax) + 1;

   
    nHmin = -nHmax;
    nKmin = -nKmax;
    nLmin = -nLmax;
    
    
       
    nHmin = 0;

    // Also, the hkl limits could be better calculated depending on the
    // spacegroup.  Below is a crude attempt at this.
    
    if (  (eSpacegroup_orthorhombic_class == m_poSpacegroup->eGetClass())
        ||(eSpacegroup_tetragonal_class == m_poSpacegroup->eGetClass())
        ||(eSpacegroup_cubic_class == m_poSpacegroup->eGetClass()) )
    {

        // tjn:  Sorry, this causes problems.
        //        nKmin = 0;
        //      nLmin = 0;
    }
    if (eSpacegroup_monoclinic_class == m_poSpacegroup->eGetClass())
    {
        // tjn:  This probably won't cause problems, but I'm removing it just to be sure.
        // nKmin = 0;
    }
    
    
    // Allocate space for up to a hemisphere of hkls within the reflection range
    
    nStat =   ((nHmax - nHmin) + 1)        // Use nStat as temp variable here
        * ((nKmax - nKmin) + 1)
        * ((nLmax - nLmin) + 1);

    

    pcCovered = new unsigned char[nStat];
    poRefln = new Crefln(NULL);
    memset(pcCovered,0,nStat);
    
    

    int nFlag;
    int nCentPhase;
    int nFplusminus;
    int a3nReducedHKL[3];

    
    for (h = nHmin; h <= nHmax; h++)
    {
        a3fH[0] = float(h);
        a3fH[1] = 0.0;
        a3fH[2] = 0.0;
        poRefln->vSetH(h);
        for (k = nKmin; k <= nKmax; k++)
        {
            a3fH[1] = float(k);
            poRefln->vSetK(k);
            for (l = nLmin; l <= nLmax; l++)
            {
                a3fH[2] = float(l);
                poRefln->vSetL(l);        

                nCentPhase =  m_poSpacegroup->nReduceHKL(poRefln, a3nReducedHKL,
                    &nFplusminus);

                if (0 <= nCentPhase) {

                    if (nAnomFlag == 0)
                        nFlag = 3;
                    else if (nCentPhase > 0)
                        nFlag = 3;
                    else if (nFplusminus==1)
                        nFlag = 1;
                    else
                        nFlag = 2;
                    
                    
                    // Due to our simplifications, we might need to flip some values in a3nReducedHKL[]
                    if ((nHmin == 0) && (a3nReducedHKL[0]<0)) {
                        // Flip to the other one.
                        a3nReducedHKL[0] = -a3nReducedHKL[0];
                        a3nReducedHKL[1] = -a3nReducedHKL[1];
                        a3nReducedHKL[2] = -a3nReducedHKL[2];
                    };
                    if ((a3nReducedHKL[0]<nHmin) ||
                        (a3nReducedHKL[1]<nKmin) ||
                        (a3nReducedHKL[2]<nLmin)) {
                        printf("ERROR:  While determining %% completeness!\n");
                    };

                    // This sometimes happens in hexagonal cells.  It's not an error.
                    if ((a3nReducedHKL[0]>nHmax) ||
                        (a3nReducedHKL[1]>nKmax) ||
                        (a3nReducedHKL[2]>nLmax)) 
                        continue;

                    
                    nStat = (nHmax - nHmin + 1)*(nKmax - nKmin +1)*(a3nReducedHKL[2] - nLmin) +
                        (nHmax - nHmin + 1)*(a3nReducedHKL[1] - nKmin) + 
                        (a3nReducedHKL[0] - nHmin);                
                    
                    if (!(pcCovered[nStat] & nFlag)) {
                        // This unique HKL was not set yet, so we count it as unique.
                        // Find it's resolution bin.
                        
                        
                        // If it already had a value (which is possible if I+ != I-),
                        // then we agree to use the bin that was already computed.
                        if (pcCovered[nStat]) {
                            nx = pcCovered[nStat] >> 2;
                            pnResoUnique[nx]++;
                        } else {
                            // Otherwise, we need to compute the resolution bin.
                            vMulMat3DVec3D(a3x3fBMat, a3fH, a3fX);
                            fDstarSq = fDot3D(a3fX, a3fX);
                            f0 = 1.0/sqrt(max(0.0000001,fDstarSq));                        
                            for (nx=0;nx<nNumReso;nx++) {
                                if ((f0>pfResoUpper[nx]) && (f0<pfResoLower[nx])) {
                                    pnResoUnique[nx]++;
                                    break;
                                };                            
                            };
                        };
                        if (nx != nNumReso) {
                            pcCovered[nStat] |= nFlag;
                            // Encode the resolution bin used in pcCovered[nStat].
                            pcCovered[nStat] = (pcCovered[nStat] & 3) + (nx << 2);
                        };
                    };                
                };
            }
        }
    }


    // Traverse through each HKL again.  Add counts for each non-zero resolution bin
    // which does not have a value of 3 (indicating that I+ and I- were accounted for).  
    // Note:  We are assuming nothing about the intersection of the sets {I+} and {I-} with
    //        the assymetric unit or the evaluated set of reflections.


    for (h = nHmin; h <= nHmax; h++)
    {
        for (k = nKmin; k <= nKmax; k++)
        {
            for (l = nLmin; l <= nLmax; l++)
            {                
                
                nStat = (nHmax - nHmin + 1)*(nKmax - nKmin +1)*(l - nLmin) +
                    (nHmax - nHmin + 1)*(k - nKmin) + 
                    (h - nHmin);
                nx = pcCovered[nStat] & 3;
                if ((nx) && (nx!=3)) {
                    nx = pcCovered[nStat] >> 2;
                    pnResoUnique[nx]++;
                };
            };
        };
    };                   


    
    delete[] pcCovered;
    delete poRefln;

    return 0;       
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function scales the cell lengths by a linear scale.  The return value is the resultant volume scale factor.
// This function is needed by strategy. 
bool Ccrystal::bLinearScaleToVolumeScale(double dLinearScale, double& dVolumeScale)
{
    if( dLinearScale <= 0.0 )
    {
        dVolumeScale = 1.0;
        
        return false;    // Invalid input. Do nothing. 
    }

    ////////////////////////////////////////////////////
    if( dLinearScale == 1.0 )
    {
        dVolumeScale = 1.0;
        return true; // nothing to do
    }
    
    double      dOriginalVolume = fCalcVolume();
    
    // Adjust crystal cell lengths by scale factor
    double      a6fCell[6];
    vGetCell(a6fCell);

    for (int ii = 0; ii < 3; ii++)
    {
        a6fCell[ii] *= dLinearScale;
    }

    vSetCell(a6fCell);

    dVolumeScale = dOriginalVolume / fCalcVolume();
    
    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////


void Ccrystal::vSetTwinLaw(int nTwinLawIndex,double* pfTwinLaw,bool bIsMatrix) {
    
    if ((nTwinLawIndex <= g_nMaxTwinLaws) && (nTwinLawIndex>0)) {
		m_nTwinLaws = max(nTwinLawIndex,m_nTwinLaws);
        if (bIsMatrix) {
			double a3x3fMatrix[3][3];
			vCopyMat3D(pfTwinLaw,&a3x3fMatrix[0][0]);
            vDeconvMat3D3XYZ(a3x3fMatrix, 
                &(m_aa3fTwinLaw[nTwinLawIndex-1][0]),&(m_aa3fTwinLaw[nTwinLawIndex-1][1]),&(m_aa3fTwinLaw[nTwinLawIndex-1][2]));

        } else {
            vCopyVec3D(pfTwinLaw,&(m_aa3fTwinLaw[nTwinLawIndex-1][0]));
		};
    };
};

int Ccrystal::nGetTwinLaw(int nTwinLawIndex,double* pfTwinLaw,bool bIsMatrix) {
    if ((nTwinLawIndex <= g_nMaxTwinLaws) && (nTwinLawIndex>0)) {
        if (bIsMatrix) {
            vConv3Rot3Vec3DMat3D(
                m_aa3fTwinLaw[nTwinLawIndex-1][0],m_aa3fTwinLaw[nTwinLawIndex-1][1],m_aa3fTwinLaw[nTwinLawIndex-1][2],
                &m_fOrientVectors[0][0], &m_fOrientVectors[1][0],
                &m_fOrientVectors[2][0], pfTwinLaw);
        } else {
              vCopyVec3D(&m_aa3fTwinLaw[nTwinLawIndex-1][0],pfTwinLaw);
		};
    } else {
		double a3x3fTempMat[3][3];
		if (bIsMatrix) {
			vIdentMat3D(a3x3fTempMat);
			vCopyMat3D(&a3x3fTempMat[0][0],pfTwinLaw);
		} else {
			vZeroMat(3,1,pfTwinLaw);
		};
	};
	return 0;
};

int Ccrystal::nClearTwinLaw(int nTwinLawIndex) {
    int nx; 
    if ((nTwinLawIndex<=m_nTwinLaws) && (nTwinLawIndex>0)) {
        for (nx = nTwinLawIndex; nx+1 <= g_nMaxTwinLaws; nx++) {
            vCopyVec3D(&m_aa3fTwinLaw[nx][0],&m_aa3fTwinLaw[nx-1][0]);
            m_afTwinFraction[nx-1] = m_afTwinFraction[nx];
            vCopyVec3D(&m_aa3fRecipShift[nx+1][0],&m_aa3fRecipShift[nx][0]);
        };
        m_nTwinLaws--;
    };
	return 0;
};

int Ccrystal::nGetRecipShift(int nTwinLawIndex,double* pfRecipShift)
{
    if( nTwinLawIndex <= m_nTwinLaws )
    {
        vCopyVec3D(&m_aa3fRecipShift[nTwinLawIndex][0], pfRecipShift);
    }
	
    return 0;
}

void Ccrystal::vSetRecipShift(int nTwinLawIndex,double* pfRecipShift)
{
    if( nTwinLawIndex <= m_nTwinLaws )
    {
        vCopyVec3D(pfRecipShift, &m_aa3fRecipShift[nTwinLawIndex][0]);
    }
}



int Ccrystal::nCalcTwinLawDeriv(int nTwinLawIndex,double fMat[3][3][3]) {
	double a3x3fRotMatTemp[3][3];

     vConv3Rot3Vec3DMat3D(
               m_aa3fTwinLaw[nTwinLawIndex-1][0], m_aa3fTwinLaw[nTwinLawIndex-1][1], m_aa3fTwinLaw[nTwinLawIndex-1][2],
		       &m_fOrientVectors[0][0], &m_fOrientVectors[1][0],
		       &m_fOrientVectors[2][0], &a3x3fRotMatTemp[0][0],
               &fMat[0][0][0],
               &fMat[1][0][0],
               &fMat[2][0][0]               
               );

     return 0;
};
int Ccrystal::nGetTwinLawRotAndVector(int nTwinLawIndex,int nMaxHKLToSearch,double* pfRot,double* pfRotErr,double* pfVector,int* pnHKL) {
    int a3nHKL[3];
    int a3nMinHKL[3];
    double a3x3fOrientMat[3][3];
    double a3x3fTwinLaw[3][3];
    double a3fTwinLawRot[3];
    double a3fRecipRot[3];
    double a3fDiff[3];
    double fTwinLawRot;
    double fDiff;
    double fMinDiff;
    int nx;
    int nNegative,nPositive;

    vGetOrientMatrix(&a3x3fOrientMat[0][0]);
    if (nGetTwinLaw(nTwinLawIndex,&a3x3fTwinLaw[0][0],true)) 
        return 1;
    nDeconvMat3DVec3DRot1(a3x3fTwinLaw,a3fTwinLawRot,&fTwinLawRot);

    fMinDiff = 1.0;
    a3nMinHKL[0] = 0;
    for (a3nHKL[0] = -nMaxHKLToSearch; a3nHKL[0] <= nMaxHKLToSearch; a3nHKL[0]++) {
        for (a3nHKL[1] = -nMaxHKLToSearch; a3nHKL[1] <= nMaxHKLToSearch; a3nHKL[1]++) {
            for (a3nHKL[2] = -nMaxHKLToSearch; a3nHKL[2] <= nMaxHKLToSearch; a3nHKL[2]++) {
                if (!bDivides(3,&a3nHKL[0])) {
                    // Build the rotation vector.
                    vZeroMat(3,1,&a3fRecipRot[0]);
                    for (nx = 0; nx < 3; nx++) {
                        double a3fTempVec[3];
                        vMulVec3DScalar(a3x3fOrientMat[nx],a3nHKL[nx],a3fTempVec);
                        vAddVec3DVec3D(a3fRecipRot,a3fTempVec,a3fRecipRot);
                    };
                    fNormVec3D(a3fRecipRot);
                    vSubVec3DVec3D(a3fRecipRot,a3fTwinLawRot,a3fDiff);
                    fDiff = fLenVec3D(a3fDiff);
                    if (fDiff < fMinDiff) {
                        fMinDiff = fDiff;
                        a3nMinHKL[0] = a3nHKL[0];
                        a3nMinHKL[1] = a3nHKL[1];
                        a3nMinHKL[2] = a3nHKL[2];
                    };
                    
                };
            };
        };
    };
    *pfRot = fTwinLawRot;
    *pfRotErr = fMinDiff;
    nNegative = 0;
    nPositive = 0;
    if (a3nMinHKL[0] < 0)
        nNegative ++;
    if (a3nMinHKL[1] < 0)
        nNegative ++;
    if (a3nMinHKL[2] < 0)
        nNegative ++;
    if (a3nMinHKL[0] > 0)
        nPositive ++;
    if (a3nMinHKL[1] > 0)
        nPositive ++;
    if (a3nMinHKL[2] > 0)
        nPositive ++;
    if (nNegative > nPositive) {
        a3nMinHKL[0] *= -1;
        a3nMinHKL[1] *= -1;
        a3nMinHKL[2] *= -1;
        *pfRot = 360 - *pfRot;
    };


    pnHKL[0] = a3nMinHKL[0];
    pnHKL[1] = a3nMinHKL[1];
    pnHKL[2] = a3nMinHKL[2];
    return 0;
};

// This function assumes that the header has been loaded with residual information, (perhaps after a call to nCalcRecipLatticePoints()
int Ccrystal::nGetSetResidValues(Cimage_header& oHeader,int nTwinID,int nTwinLaw,double* pfValues,int nAddDelete) {
    const int nBufSize = 200;
    double afResidBuf[nBufSize + 1];
    int nNumValues;
    int nPointer;
    int nx;

    if (nAddDelete == -1) {
        oHeader.nDelete((Cstring) D_K_CrystalResidInfo);
        return 0;
    };

    for (nPointer = 0; nPointer < nBufSize;nPointer++)
        afResidBuf[nPointer] = -1;

    oHeader.nGetValue((Cstring) D_K_CrystalResidInfo,200,&afResidBuf[0]);
    for (nNumValues = 0; (afResidBuf[nNumValues] != -1); nNumValues++);

    if (nAddDelete == 1) {
        afResidBuf[nNumValues++] = nTwinID;
        afResidBuf[nNumValues++] = nTwinLaw;
        afResidBuf[nNumValues++] = g_nResidInfoRMSEntries;
        for (nx = 0; nx <g_nResidInfoRMSEntries; nx++)
            afResidBuf[nNumValues++] = pfValues[nx];
        oHeader.nReplaceValue((Cstring) D_K_CrystalResidInfo,nNumValues,&afResidBuf[0]);
        return 0;
    };
    
    for (nPointer = 0; nPointer < g_nResidInfoRMSEntries; nPointer++)
        pfValues[nPointer] = -1;

    for (nPointer = 0; nPointer < nNumValues;) {
        if (afResidBuf[nPointer + 2] <= 0)
            break;
        if ((afResidBuf[nPointer] = nTwinID) && (afResidBuf[nPointer + 1] == nTwinLaw)) {
            for (nx = 0; nx < afResidBuf[nPointer + 2]; nx++)
                pfValues[nx] = afResidBuf[nPointer + 2 + 1 + nx];
            return 0;
        };
        nPointer += 2;
        nPointer += (int) afResidBuf[nPointer] + 1;
    };
    return 1;
};

double 
Ccrystal::dCalcGetResolution(Crefln *poRefln)
{
  // Calculate the resolution in Angstrom of a refln with a valid hkl
  // For this to work, the B matrix MUST already be available

  double dDstar;
  double a3dHKL[3], a3dDstar[3];
  a3dHKL[0] = (double) poRefln->nGetH();
  a3dHKL[1] = (double) poRefln->nGetK();
  a3dHKL[2] = (double) poRefln->nGetL();

  vMulMat3DVec3D(m_fBMatrix, &a3dHKL[0], &a3dDstar[0]);
  dDstar = fLenVec3D(a3dDstar);
  if (0.0 < dDstar)
    return ( 1.0 / dDstar);
  else
    {
      cout << "Ccrystal::dCalcGetResolution() programmer error!\n" << endl;
      return (0.0);
    }
}
