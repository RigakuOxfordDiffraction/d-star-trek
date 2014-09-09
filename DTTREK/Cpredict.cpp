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
// Cpredict.cc            Initial author: J.W. Pflugrath           01-May-1995
//  This file contains the member functions of class Cpredict which implements
//    the reflection predicting encapsulation of d*TREK.
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

#include "Dtrek.h"
#include "Cpredict.h"         // Class definition and prototypes

using namespace std;

//+Definitions, constants, and initialization of static member variables

//+Code begin

//+Public functions

// Constructors, destructors and assignments

Cpredict::Cpredict()
{
    nInitValues();
}

Cpredict::Cpredict(Cimage_header& oHeader, 
                   Creflnlist*    poReflnlistIn,
                   int            nScanNumber, 
                   const Cstring& sStrategyPrefix)
{
    nInitValues();

    m_sStrategyPrefix = sStrategyPrefix;
    m_nScanNumber = nScanNumber;

    nInitValues(oHeader);

    m_poReflnlist   = poReflnlistIn;

    m_bNewReflnlist = FALSE;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Cpredict::Cpredict(Csource*      poSourceIn, 
                   Cdetector*    poDetectorIn,
                   Ccrystal*     poCrystalIn, 
                   Cgoniometer*  poCrysGonioIn,
                   Crotation*    poRotationIn,
                   Creflnlist*   poReflnlistIn)
{
  nInitValues();

  // There should be some error checking below to make sure that the
  // pointers actually point to real objects.

  m_poSource      = poSourceIn;
  m_poCrysGonio   = poCrysGonioIn;
  m_nNumDetectors = 1;
  m_nNumCrystals  = 1;
  m_nMaxTwinLaws  = 0;

  m_ppoCrystal    = new Ccrystal*[1];
  *m_ppoCrystal   = poCrystalIn;
  m_pbActive      = new bool[1];
  *m_pbActive     = TRUE;
 
  m_ppoDetector   = new Cdetector* [m_nNumDetectors];
  *m_ppoDetector  = poDetectorIn;

  // The following may be NULL pointers.

  m_poRotation    = poRotationIn;
  m_poReflnlist   = poReflnlistIn;
}

Cpredict::~Cpredict()
{
  if (m_bNewDetNames && (NULL != m_psDetectorNames) )
    {
      delete [] m_psDetectorNames;
      m_psDetectorNames = NULL;
      m_bNewDetNames = FALSE;
    }

  if (m_bNewReflnlist && (NULL != m_poReflnlist) )
    {
      delete m_poReflnlist;
      m_poReflnlist = NULL;
      m_bNewReflnlist = FALSE;
    }

  if (m_bNewDetector && (NULL != m_ppoDetector) )
    {
      for (int i = 0; i < m_nNumDetectors; i++)
        {
          delete *(m_ppoDetector+i);
          *(m_ppoDetector+i) = NULL;
        }
      delete [] m_ppoDetector;
      m_ppoDetector = NULL;
      m_bNewDetector = FALSE;
    }
  else if (NULL != m_ppoDetector)
    {
      delete [] m_ppoDetector;
      m_ppoDetector = NULL;
    }

  if (m_bNewCrystal && (NULL != m_ppoCrystal) )
    {
      for (int i = 0; i < m_nNumCrystals; i++)
      {
          delete *(m_ppoCrystal+i);
          *(m_ppoCrystal+i) = NULL;
      }
      delete[] m_ppoCrystal;
      m_ppoCrystal = NULL;
      m_bNewCrystal = FALSE;
    }
  else if (NULL != m_ppoCrystal)
  {
      delete[] m_ppoCrystal;
      m_ppoCrystal = NULL;
  }
  if (m_pbActive) {
      delete[] m_pbActive;
      m_pbActive = NULL;
  };

  
  if (m_bNewSource && (NULL != m_poSource) )
    {
      delete m_poSource;
      m_poSource = NULL;
      m_bNewSource = FALSE;
    }


  if (m_bNewCrysGonio && (NULL != m_poCrysGonio) )
    {
      delete m_poCrysGonio;
      m_poCrysGonio = NULL;
      m_bNewCrysGonio = FALSE;
    }

  if (m_bNewRotation && (NULL != m_poRotation) )
    {
      delete m_poRotation;
      m_poRotation = NULL;
      m_bNewRotation = FALSE;
    }
  m_nNumDetectors  = 0;
}

int Cpredict::nInitValues(void)
{
  m_nStat            = 0;
  m_nDisplay         = 0;
  m_nNumDetectors    = 0;
  m_nNumCrystals     = 0;
  m_nWhichDetector   = 0;
  m_nWhichCrystal    = 0;
  m_nWhichCrystalTwinLaw = 0;

  m_psDetectorNames  = NULL;
  m_pbActive         = NULL;
  m_poReflnlist      = NULL;
  m_ppoDetector      = NULL;
  m_poSource         = NULL;
  m_ppoCrystal       = NULL;
  m_poCrysGonio      = NULL;
  m_poRotation       = NULL;

  m_bNewDetNames     = FALSE;
  m_bNewDetector     = FALSE;
  m_bNewSource       = FALSE;
  m_bNewCrystal      = FALSE;
  m_bNewCrysGonio    = FALSE;
  m_bNewRotation     = FALSE;
  m_bNewReflnlist    = FALSE;
  
  m_bUseEdgeReso     = FALSE;
  m_fReflnWidthMax   = 8.0;
  m_fReflnWidthMaxRad = m_fReflnWidthMax * fRADIANS_PER_DEGREE;
  m_fResolutionMin   = 0.0;  // Default resolution is the ENTIRE detector!
  m_fResolutionMax   = 0.0;  // Default resolution is the ENTIRE detector!
  m_fMaxLorentz      = 50.0;
  if ("" != sGetEnv("DTREK_MAXLORENTZ"))
    m_fMaxLorentz = atof(sGetEnv("DTREK_MAXLORENTZ").string());

  m_nNonunfFlag      = 0;
  m_fScaleLength     = 1.0;
  m_bCalcImageDrift  = FALSE;
  m_bCalcMosaicityLinearCoeffs = FALSE;
  m_bCalcMinMosaicityForRefln = FALSE;
  m_bDetectorIsCylindrical     = FALSE;

  m_bDoFastPrediction = FALSE;
  m_fPrePredictRotWidth = 0.0;

  m_wVerbose = 0U;

  m_nScanNumber = 0;

  m_sStrategyPrefix = "";
  
  for(int ii= 0; ii < 3; ii++)
      m_a3fLocalDN[ii] = 0.0;

  return 0;
}

int Cpredict::nInitValues(Cimage_header& oHeader)
{
  int     i = 0;       // Loop counter
  int     nTemp = 0;
  Cstring sTemp("");
  Ccrystal*         poCrystal = NULL;

  if (!oHeader.bIsAvailable())
    {
      cout << "Cpredict ERROR: image header is not valid."
           << "  Cannot construct predict!\n";
      m_nStat = 1;
    }
  else
    {
      // Try to get required information from the header

      m_nStat = 0;
      
      poCrystal     = new Ccrystal(oHeader);
      m_nNumCrystals  = poCrystal->nNumTwins();
          m_nMaxTwinLaws  = poCrystal->nMaxTwinLaws();
      if (m_nNumCrystals >= 1) {
          m_ppoCrystal  = new Ccrystal* [m_nNumCrystals];
          m_pbActive    = new bool[m_nNumCrystals];
          m_ppoCrystal[0] = poCrystal;
          m_pbActive[0]   = TRUE;
      } else
          m_nStat = 1;
      for (i = 1; i< m_nNumCrystals; i++) {
          poCrystal = new Ccrystal (oHeader,i+1);
          m_ppoCrystal[i] = poCrystal;
          m_pbActive[i] = TRUE;
      };

      m_bNewCrystal   = TRUE;

      m_poCrysGonio = new Cgoniometer(oHeader, sBuildStrategyScanPrefix(m_sStrategyPrefix, Ccrystal::ms_sCrystalPrefix, m_nScanNumber));
      
      m_bNewCrysGonio = TRUE;

      // Get the number of detectors

      nTemp = oHeader.nGetValue(Cdetector::ms_sDetectorNumber, &m_nNumDetectors);

      if (0 == nTemp)
        {
          // Get the detector names

          m_psDetectorNames = new Cstring [m_nNumDetectors];
          m_bNewDetNames = TRUE;
          nTemp = oHeader.nGetValue(Cdetector::ms_sDetectorNames, m_nNumDetectors,
                                    m_psDetectorNames);
          if (0 == nTemp)
            {
              m_ppoDetector = new Cdetector* [m_nNumDetectors];
              for (i = 0; i < m_nNumDetectors; i++)
                {
                  m_ppoDetector[i] = new Cdetector (oHeader, m_psDetectorNames[i]);
                }
              m_bNewDetector  = TRUE;
            }
          else
            m_nStat++;
        }
      else
        m_nStat++;

      m_poSource      = new Csource(oHeader);
      m_bNewSource    = TRUE;
      
      m_poRotation    = new Crotation(oHeader, sBuildStrategyScanPrefix(m_sStrategyPrefix,
                                                                        D_K_ScanPrefix, 
                                                                        m_nScanNumber));
      m_bNewRotation  = TRUE;
    }
  return (m_nStat);
}

int Cpredict::nList()
{
  float fTempW,fREdge,fRMin,fRMax;
  float fTempS0[3];
  int nx;

  cout << "Predict listing:\n";
  for (nx=0;nx<m_nNumCrystals;nx++) {
      if (m_pbActive[nx]) {
          (void) m_ppoCrystal[nx]->nList();
          double a3x3fOMat[3][3];
          m_ppoCrystal[nx]->nCalcOrientMatrix();
          m_ppoCrystal[nx]->vGetOrientMatrix(&a3x3fOMat[0][0]);
          if ("" != sGetEnv("DTREK_PREDICT_MATRICES"))
            {
              cout << "Crystal Orientation ";
              vListMat3D(&a3x3fOMat[0][0]);
              cout << "Crystal B Matrix ";
              vListMat3D(&m_ppoCrystal[nx]->m_fBMatrix[0][0]);
              cout << "Crystal Rot Matrix ";
              vListMat3D(&m_ppoCrystal[nx]->m_fRotMatrix[0][0]);
            }
      };
  };

  (void) m_poCrysGonio->nList();
  // List the crystal goniometer matrix separately
  double a3x3dRotMatrix[3][3];
  m_poCrysGonio->vCalcGetRotMatrix(&a3x3dRotMatrix[0][0], -1);
  if ("" != sGetEnv("DTREK_PREDICT_MATRICES"))
    {
      cout << "Crystal Goniometer ";
      vListMat3D(&a3x3dRotMatrix[0][0]);
    }

  (void) m_poSource->nList();
  fTempW = m_poSource->fGetWavelength();
  m_poSource->vCalcGetS0(&fTempS0[0]);
  for (int i = 0; i < m_nNumDetectors; i++)
    {
      (void) m_ppoDetector[i]->nList();
      (void) m_ppoDetector[i]->nGetResolution(fTempS0, &fRMin, &fRMax, &fREdge);
      cout << "DetResolution min:  " << fRMin  * fTempW << endl;
      cout << "DetResolution max:  " << fRMax  * fTempW << endl;
      cout << "DetResolution edge: " << fREdge * fTempW << endl;
    }
  // Get rotation axis direction
  int nAxis;
  float a3fTemp[3];
  nAxis = m_poCrysGonio->nGetNum(m_poRotation->sGetName());
  if (0 > nAxis)
    {
      nAxis = 0;
      cout << "WARNING Cpredict::nSetupRotation, axis name '"
           << m_poRotation->sGetName() << "' not a valid axis\n"
           << " for this crystal goniometer, using '" 
           << m_poCrysGonio->sGetName(nAxis) << "' instead.\n";
    }
  if (m_poCrysGonio->nGetRotVector(nAxis,&a3fTemp[0])) 
    {
      cout << "ERROR Cpredict::nSetupRotation calculating rotation axis vector.\n" << flush;
    }
  //+jwp 9-Dec-2003
  if (NULL != m_poRotation)
    {
      // Copy the actual rotation vector into the Crotation object
      m_poRotation->vSetVector(&a3fTemp[0]);
    }
  //-jwp 9-Dec-2003
  (void) m_poRotation->nList();
  cout << "Predict resol min:          " << m_fResolutionMin << endl;
  cout << "Predict resol max:          " << m_fResolutionMax << endl;
  cout << "Predict max Lorentz factor: " << m_fMaxLorentz << endl << flush;
  return 0;
}


int Cpredict::nPredict(Creflnlist& oListIn,double fPad) {
        double afGonio[20];
        double afGonioRef[20];
        int        nRotAxis;
        int    anGonioRotIndices[20];
        const int nMaxRotIndex = 5;
        int    nRots,nRotIndex;
        double fRotStart,fRotEnd;
        double fRotStartRef,fRotEndRef;
        double fRotInc = m_poRotation->fGetIncrement();
        bool   bNewRotRange;
        int        nStat;
        int        nRef;
        

        anGonioRotIndices[0] = oListIn.m_nFI_fGonio1;
        anGonioRotIndices[1] = oListIn.m_nFI_fGonio2;
        anGonioRotIndices[2] = oListIn.m_nFI_fGonio3;
        anGonioRotIndices[3] = oListIn.m_nFI_fGonio4;
        anGonioRotIndices[4] = oListIn.m_nFI_fGonio5;

    for (nRotIndex = 0; nRotIndex < nMaxRotIndex; nRotIndex++) {
                afGonio[nRotIndex] = -999.0;
        afGonioRef[nRotIndex] = -999.0;
    };
        fRotStart = -999.0;
        fRotEnd = -999.0;

    // Add width information for backwards compatiblity.
    if (oListIn.m_nFI_fObsRotWidth<0)
        (void) oListIn.nExpandGetField(oListIn.ms_sfObsRotWidth);
    for (nRef = 0; nRef < oListIn.nGetNumReflns(); nRef++)
      {
        if (oListIn[nRef].fGetField(oListIn.m_nFI_fObsRotWidth) < 0.0)
          oListIn[nRef].vSetField(oListIn.m_nFI_fObsRotWidth,(float) fRotInc);
      }
    
    for (nRef = 0; nRef < oListIn.nGetNumReflns(); nRef++)
      {
        for (nRots = 0,nRotIndex = 0; nRotIndex < nMaxRotIndex; nRotIndex++)
          {
            if (anGonioRotIndices[nRotIndex]>=0)
              {
                afGonioRef[nRots] = oListIn[nRef].fGetField(anGonioRotIndices[nRotIndex]);
                nRots++;
              } 
            else if (nRots < m_poCrysGonio->m_nNumRotValues)
              {
                afGonioRef[nRots] = m_poCrysGonio->m_pfDatumValue[nRots];
                nRots++;
              } 
          }
        if ((oListIn.m_nFI_fObsRotMid>=0) && (oListIn.m_nFI_fObsRotWidth>=0))
          {
            fRotStartRef = oListIn[nRef].fGetField(oListIn.m_nFI_fObsRotMid)
                           - 0.5*oListIn[nRef].fGetField(oListIn.m_nFI_fObsRotWidth);
            fRotEndRef = oListIn[nRef].fGetField(oListIn.m_nFI_fObsRotMid)
                           + 0.5*oListIn[nRef].fGetField(oListIn.m_nFI_fObsRotWidth);
          }
        else
          return 1;  

        if (oListIn.m_nFI_nGonioRotAxis >= 0)
          nRotAxis = oListIn[nRef].nGetField(oListIn.m_nFI_nGonioRotAxis) - 1;
        else
          nRotAxis = m_poCrysGonio->nGetNum(m_poRotation->sGetName());


        // See if we have had a change in variables.
        for (nRotIndex = 0; nRotIndex < nMaxRotIndex; nRotIndex++)
          {
            if (ABS(afGonio[nRotIndex] - afGonioRef[nRotIndex])>0.01)
              break;
          }

        if (   (    (fRotStartRef > fRotStart - 0.01) 
                 && (fRotStartRef < fRotEnd + 0.01) )
            || ( (fRotEndRef > fRotStart - 0.01) 
                 && (fRotEndRef < fRotEnd + 0.01) )) 
          {
            bNewRotRange = false;
            fRotStart = min(fRotStartRef,fRotStart);
            fRotEnd = max(fRotEndRef,fRotEnd);
          } 
        else
          bNewRotRange = true;

        if ((   (nRotIndex != nMaxRotIndex)
             || (nRef + 1 == oListIn.nGetNumReflns())
             || (bNewRotRange))
             && (nRef != 0)) 
          {
            // We should run the prediction algorithm on this set of images.
            // First we must update the goniometer to the correct settings.
            for (nRotIndex = 0; nRotIndex < m_poCrysGonio->m_nNumRotValues; nRotIndex++)
              m_poCrysGonio->m_pfDatumValue[nRotIndex] = afGonio[nRotIndex];
            m_poRotation->vSetName(m_poCrysGonio->m_psName[nRotAxis]);
            m_poRotation->nSetRotRange(fRotStart-fPad,fRotEnd+fPad);

            // Run the prediction algorithm.
            nStat  = nPredict();
            if (nStat)
              return nStat;

            bNewRotRange = true;
          }

        // Copy goniometer settings over.
        for (nRotIndex = 0; nRotIndex < nMaxRotIndex; nRotIndex++) 
          afGonio[nRotIndex] = afGonioRef[nRotIndex];
        // Update the rotation start and end.
        if (bNewRotRange)
          {
            fRotStart = fRotStartRef;
            fRotEnd = fRotEndRef;
          } 
        else 
          {
            fRotStart = min(fRotStartRef,fRotStart);
            fRotEnd = max(fRotEndRef,fRotEnd);
          }
      }
    return 0;
}

int Cpredict::nPredict(Crotation *poRotationIn, Creflnlist *poReflnlistIn)
{
    int     nStat = 0;
    int     nDet = 0;
    int     nThisBatchesFirstRefln = 0;
    int     nFirstBatchFirstRefln = 0;
    int     i = 0, j = 0;               // Loop counter
    int     nRef = 0;                   // Temporary variable to hold number of reflections

    if( poRotationIn != NULL )
    {
        // A valid pointer to a Crotation object was passed.  Use this new
        // rotation object, but first get rid of old one if necessary.
        if( (m_bNewRotation) && (m_poRotation != NULL) )
        {
            delete m_poRotation;
            m_bNewRotation = FALSE;
        }
        
        m_poRotation = poRotationIn;
    }

    if( poReflnlistIn != NULL )
    {
        // A valid pointer to a Creflnlist object was passed.  Use this new
        // reflnlist object, but first get rid of old one if necessary.
        if( (m_bNewReflnlist) && (m_poReflnlist != NULL) )
        {
            delete m_poReflnlist;
            m_bNewReflnlist = FALSE;
        }
        
        m_poReflnlist = poReflnlistIn;
    }

    nStat = nExpandReflnlist();

    if( m_bDoFastPrediction && m_nNumCrystals == 1 && m_ppoCrystal[0]->m_nTwinLaws == 0 )
    {

        nSetupSourceCrystal(0,0);
        
        nSetupRotation(m_poRotation->fGetRotStart(),m_poRotation->fGetRotEnd());
        
        vCross3D(m_afRotVec, m_afS0, m_afE3crossS0);    // Member variable for later use
        
        nSetupDetector(0);
        
        if (!m_afPrePredictRotMid.size())
        {
            if (nPredictWholeDataSet(m_fWholeRotStart,m_fWholeRotEnd))
                return 1;
        }
        
        if (nPredictReflnsFast(m_poRotation->fGetRotStart(),m_poRotation->fGetRotEnd()))
            return 1;
        
        m_poReflnlist->vSetState(eReflnlist_predicted_state);
        
        return 0;
    }

    for(int nCrystal = 0; nCrystal < m_nNumCrystals; nCrystal++)
    {
        for(int nTwinLaw = 0; nTwinLaw <= m_ppoCrystal[nCrystal]->m_nTwinLaws; nTwinLaw++)
        {
            // Are we predicting for this crystal component?

            nFirstBatchFirstRefln = m_poReflnlist->nGetNumReflns();

            if( !m_pbActive[nCrystal] )
                continue;

            nStat = nSetupSourceCrystal(nCrystal,nTwinLaw);   // Setup source and crystal
            
            if (0 != nStat)
            {
                cout << "ERROR Cpredict:nPredict, setup of source and crystal failed: " << nStat << endl;
                
                return nStat;
            }

            // Loop over required rotation range in at most 4 degree steps

            m_fRotOverallEnd   = m_poRotation->fGetRotEnd();
            m_fRotOverallStart = m_poRotation->fGetRotStart();

            double      fRotS = 0.0;

            // fRotS will be set inside the loop from fRotE, so set fRotE here

            double      fRotE = min(m_fRotOverallStart, m_fRotOverallEnd) - m_fReflnWidthMax;
            
            m_fRotOverallEnd = max(m_fRotOverallStart, m_fRotOverallEnd);

            do 
            {
                nThisBatchesFirstRefln = m_poReflnlist->nGetNumReflns();

                fRotS = fRotE;                 // Set start to the last end
                
                fRotE = fRotS + 4.0;           // Set end 4 degrees further on or the end

                // Do not cut off the begining.
                if (fRotE > m_fRotOverallEnd + m_fReflnWidthMax) 
                    fRotE = m_fRotOverallEnd + m_fReflnWidthMax;

                nStat = nSetupRotation(fRotS, fRotE);  // Setup the final rotation matrices

                vCross3D(m_afRotVec, m_afS0, m_afE3crossS0);    // Member variable for later use
                nStat = nSortAxes();          // Sort cell axes for efficiency in prediction

                // Loop over each detector

                for(nDet = 0; (nDet < m_nNumDetectors) && (!nStat); nDet++)
                {
                    nStat = nSetupDetector(nDet);  // Setup the detector
                    
                    // Need fDstarSqMax, fDstarSqMin in next routine JWP
                    if (!nStat)                  
                        nStat = nPredictReflns(fRotS, fRotE); // predict! prdwtx
                }

                if( 0 != nStat )
                    return 1;

                // Clean-up list because of the bit of sloppiness in nPredictReflns

                // Delete reflections which are larger or smaller than the total reflection width.
                nRef = m_poReflnlist->nGetNumReflns();
                
                for(i = nRef-1; i >= nThisBatchesFirstRefln; i--)
                {
                    if(   (m_poReflnlist->poGetRefln(i)->fGetField(m_poReflnlist->m_nFI_fCalcRotEnd)   <= m_fRotOverallStart)
                       || (m_poReflnlist->poGetRefln(i)->fGetField(m_poReflnlist->m_nFI_fCalcRotStart) >  m_fRotOverallEnd) )
                    {
                        m_poReflnlist->nDelete(i);
                    }
                }              
            } while (fRotE < m_fRotOverallEnd + m_fReflnWidthMax);

            // remove all duplicates.
            int*        pnSort    = new int[m_poReflnlist->nGetNumReflns() - nFirstBatchFirstRefln];
            int*        pnSortI   = new int[m_poReflnlist->nGetNumReflns() - nFirstBatchFirstRefln];
            
            for(i = nFirstBatchFirstRefln,j = 0; i < m_poReflnlist->nGetNumReflns(); i++,j++)
            {
                pnSort[j] = (*m_poReflnlist)[i].nPackHKL();
                pnSortI[j] = j;
            }

            /////////////////////////////////////////////////////////////////
            int*        pnDelete  = new int[m_poReflnlist->nGetNumReflns()];
            for (i = 0; i < m_poReflnlist->nGetNumReflns(); i++)
                pnDelete[i] = 0;
            /////////////////////////////////////////////////////////////////

            g_pnCmpInts = pnSort;
            qsort(&pnSortI[0], m_poReflnlist->nGetNumReflns() - nFirstBatchFirstRefln, sizeof(int), int_cmp_rel_2);

            for(i = 0; i < m_poReflnlist->nGetNumReflns() - nFirstBatchFirstRefln; i++)
                pnSortI[i] += nFirstBatchFirstRefln;

            double      fCalcRotMid_i        = 0.0;
            double      fCalcRotMid_j        = 0.0;
            
            double      fCalcRotWidth_i      = 0.0;
            double      fCalcRotWidth_j      = 0.0;

            int         nPackedHKL_i = 0;
            int         nPackedHKL_j = 0;

            for (i = nFirstBatchFirstRefln; i < m_poReflnlist->nGetNumReflns()-1; i++)
            {
                fCalcRotMid_i   = (*m_poReflnlist)[pnSortI[i-nFirstBatchFirstRefln]].fGetField(m_poReflnlist->m_nFI_fCalcRotMid);
                fCalcRotWidth_i = (*m_poReflnlist)[pnSortI[i-nFirstBatchFirstRefln]].fGetField(m_poReflnlist->m_nFI_fCalcRotWidth);

                nPackedHKL_i    = (*m_poReflnlist)[pnSortI[i-nFirstBatchFirstRefln]].nPackHKL();

                for(j=i+1; j < m_poReflnlist->nGetNumReflns(); j++)
                {
                    nPackedHKL_j = (*m_poReflnlist)[pnSortI[j-nFirstBatchFirstRefln]].nPackHKL();
                    
                    if( nPackedHKL_i == nPackedHKL_j )
                    {
                        fCalcRotMid_j   = (*m_poReflnlist)[pnSortI[j-nFirstBatchFirstRefln]].fGetField(m_poReflnlist->m_nFI_fCalcRotMid);
                        fCalcRotWidth_j = (*m_poReflnlist)[pnSortI[j-nFirstBatchFirstRefln]].fGetField(m_poReflnlist->m_nFI_fCalcRotWidth);

                        if( ABS(fCalcRotMid_i - fCalcRotMid_j) < (fCalcRotWidth_i + fCalcRotWidth_j) * 0.5 )
                        {
                            pnDelete[pnSortI[j-nFirstBatchFirstRefln]] = -1;
                        } 
                        // TJN: Don't break here.  The sorting on packed HKL might mix up the rotation midpoints.
                    }
                    else
                        break;
                }
            }

            m_poReflnlist->nDelete(-1, pnDelete);
            
            delete[] pnSort;
            delete[] pnSortI;
            delete[] pnDelete;
        }
    }
    
    m_poReflnlist->vSetState(eReflnlist_predicted_state);

    return 0;
}

int Cpredict::nSolveQuadratic(const double fA, const double fB, const double fC,
                              double *pfV)
{
  // This solves the special quadratic equation: AX^2 + 2BX + C = 0
  // used by the nPredictReflns routine.  The (up to 2) solutions are returned
  // in pfV[0] and pfV[1].  The number of solutions is returned by the routine.

  double fD;

  if (fA != 0.0) {
    fD = fB * fB - fA * fC;
    if (fD == 0.0) {
      *pfV = -fB / fA;
      return (1);
    }
    else if (fD > 0.0) {
      fD = sqrt(fD);
      pfV[0] = (-fB + fD)  / fA;
      pfV[1] = (-fB - fD)  / fA;
      return (2);
    }
  }
  else if (fB != 0.0) {

    //  Here if A = 0 and B not = 0

    *pfV = -0.5f * fC / fB;
    return (1);
  }
  return (0);
}


void Cpredict::vSwap2Values(int* pnInt1, int* pnInt2)
{
  register int nTemp;
  nTemp   = *pnInt1;
  *pnInt1 = *pnInt2;
  *pnInt2 = nTemp;
}
void Cpredict::vSwap2Values(double* pfFloat1, double* pfFloat2)
{
  register double fTemp;
  fTemp     = *pfFloat1;
  *pfFloat1 = *pfFloat2;
  *pfFloat2 = fTemp;
}

int Cpredict::nOutsideLimits(const int nMin, const int nMax,
                            int* pnBegin,  int* pnEnd)
{
  //  Checks whether loop indices are outside limits.
  //  Used in REEKE algorithm - called by nPredictReflns

  //  Returns 0,  limits are ok , things may have been adjusted!
  //          1,  limits are bad
  //

  if ( (*pnBegin <= nMax) && (*pnEnd >= nMin) ) {
//    *pnBegin = *pnBegin >? nMin;
//    *pnEnd   = *pnEnd <? nMax;
    *pnBegin = max(*pnBegin , nMin);
    *pnEnd   = min(*pnEnd , nMax);
    return (0);
  }
  return (1);
}

void Cpredict::vAdjustLoopLimits(const int nLoop,
                                 const double fQ1, const double fQ2,
                                 const double fR1, const double fR2,
                                 int* pnBegin, int* pnEnd)
{
  // Used in REEKE algorithm - called by nPredictReflns
  // Adjust beginning and ending reflection index loop limits

  pnBegin[0] = nIntMinusOne(fQ1);
  pnEnd[0]   = nIntPlusOne(fQ2);
  if (nLoop > 1) {
    pnBegin[1] = nIntMinusOne(fR1);
    pnEnd[1]   = nIntPlusOne(fR2);
  }
}

double Cpredict::fCalcGetPolarz(void)
{
  // Calculates and returns polarization correction factor
  //
  // Uses member variables:
  //  m_afXR[]  Dimensionless reciprocal lattice coordinate of
  //            reflection in the laboratory frame when reflection
  //            lies on Ewald sphere (at rotation angle ROT)
  //  m_afS0[]
  //  m_afPol[]
  //
  // Calculates:
  //           Polarization correction factor
  //  True intensity is proportional to LOPOL*(measured intensity)
  //  Polarization correction "POL" for twice-reflected beam from
  //  Azaroff, ACTA CRYST.(1955) 8, 701-704 as formulated by PHILLIPS,
  //  et al., ACTA CRYST.(1977) A33, 445-455 (appendix).
  //   Kabsch, (1988) J. Appl.Cryst p. 916 (THIS IS THE BEST ONE TO READ!)

  // General specification of polarization.  If degree of polarization = 0.5,
  // then beam is unpolarized.  Note:  The following code yields
  //
  //              (1 + costwotheta**2)
  //        POL = ------------------, if beam is unpolarized! (so don't add it!)
  //                    2

  //  The intensity of scattering is proportional to sinsquaredC, where
  //  C is the angle between the electric field vector of the incident
  //  beam and the direction of the diffracted beam wavevector.

  //  Let A be the angle between the diffracted beam wavevector S and
  //  the polarization plane.  Let B be the angle beteen the diffracted
  //  beam wavevector S and the normal to the polarization plane POLNRM.
  //  The polarization plane is perpendicular to both the primary beam
  //  wavevector S0 and the normal to the polarization plane, so its 'direction'
  //  can be defined by SN = S0 x POLNRM, which has been calculated in
  //  nSetupSourceCrystal.
  //  If the beam is polarized by fraction POLARZ(1) in the polarization
  //  plane and (1 - POLARZ(1)) in the normal to the polarization plane, we
  //   have:
  //   POL =
  //   sinsquaredC = POLARZ(1) * sinsquaredA + (1 - POLARZ(1)) * sinsquaredB
  //    (convert sines to cosines by using identity: sinX**2 = 1 - cosX**2)
  //              = POLARZ(1) * (1 - cossquaredA) + (1 - POLARZ(1)) *
  //                                                (1 - cossquaredB)
  //
  //               = 1 - [POLARZ(1) * cossquaredA + (1 - POLARZ(1)) *
  //                                                 cossquaredB]
  //  but cosA = (S . SN) / (|S| |SN|), and cosB = (S . POLNRM) / (|S| |POLNRM|)
  //  We have though:  |S0| = |S| = |SN| = |POLNRM| = 1, that is, all
  //  these vectors are unit length (in this routine anyways!).
  //  So POL = 1 - [POLARZ(1) * (S.SN)**2 + (1-POLARZ(1)) * (S.POLNRM)**2]
  //  S = XR - S0 in our geometry, so now make the calculations:


  double fS[3], fS1, fS2;
  vSubVec3DVec3D(m_afXR, m_afS0, fS);
  fS1 = fDot3D(fS, m_afSN);
  fS2 = fDot3D(fS, &m_afPol[1]);
  return (1.0f - ( m_afPol[0] * fS1*fS1  + (1.0f - m_afPol[0]) * fS2*fS2));
}

int Cpredict::nCalcRecipCoords(
        const double fX[3],        // Input dimensionless recip. lattice
                                  //   coordinate in lab frame
        const double fRotStartRad, // Starting rotation range within which we are
                                  //   are predicting (in RADIANS)
        const double fRotEndRad)   // Ending rotation range within which we are
                                  //   are predicting (in RADIANS)
{
  // Calculate reciprocal lattice coordinates and rotation values
  // of a reflection.  Things are returned in member variables.

  // Note: in comments below:
  //     E3 is the rotation axis vector
  //     EPS is crystal effective mosaic spread in radians (REPS is in degrees)

  // Also use: fReflnWidthMaxRad, fCrysMosaicityRad, fSourceDispersion
  //           fDstarSqMax, fDstarSqMin, fRotVec, fS0

  // Local variables

  double fE1S0, fE2S0, fE3S0, fCEA, fCEB, fCEC, fXCE1, fXCE3, fARG1, fROTC,fROTCA,fROTCB,f0;
  double fROTC1[3][3], fDSCOTH;
  double fE1[3], fE2[3], fT1, fT2;
  double fRS[3];
  double fDE2;
  double fCEACEB;
  double fCrysMosaicityRad;
  double a3fCrysMosaicity[3];

  // Calculate resolution as Dstar squared and check if within limits
  //    Remember this is dimensionless (i.e. wavelength = 1.0!!!!)

  m_fReflnDstarSq = fDot3D(fX, fX);

  // Check resolution limits

  if ( (m_fReflnDstarSq > m_fDstarSqMax) || (m_fReflnDstarSq < m_fDstarSqMin) )
    return (1);

  // Rotate lattice point (by fROTC) to bring reflection onto the
  // Ewald sphere.
  // For this purpose calculate Rot-value at the
  // reflecting position. Apply this rotation to the reciprocal
  // space coordinates of the reflection (at the start of the rotation)
  // to get the coordinates of the reflection on Ewald-sphere (fXR).

  // We first set up 3 mutually orthogonal unit vectors.
  // fE3 already defined along rotation axis as fRotVec
  // fE2 perp. to fE3 and fX
  // fE1 perp. to fE3 and fE2 (to give a right handed set)
  //
  //    E2 = [E3 x XC] / ||E3 x XC||
  //    --    --   --      --   --

  vCross3D(m_afRotVec, fX, fE2);
  fDE2 = fNormVec3D(fE2);
  if (fDE2 <= 0.0) return (2);         // Error if length of fE2 is not positive

  //  E1 = [E2 x fRotVec]

  vCross3D(fE2, m_afRotVec, fE1);

  // Determination of fROTC
  //
  // Coefficients for ROT equation CEA, CEB, CEC
  //
  // CEA*cos(ROTc) + CEB*sin(ROTc) = CEC
  //
  // CEA = (X.E1) * (E1.S0)
  //        - --     -- --
  // CEB = (X.E1) * (E2.S0)
  //        - --     -- ---
  // CEC = 0.5 * (X.E1)**2 + 0.5 * (X.E3)**2  - (X.E3) * (E3.S0)
  //              - --              - --         - --     -- --

  /*   
     For an alternate derivation/solution of this problem, see Cimagepredict.cc.
     I attempted implementing the code in Cimagepredict.cc here, but I was
     unable to improve the accuracy of these numbers.  This was much to my suprise,
     since this derivation uses an acos() below, where fARG1 is often near 1.0.  
     I guess the acos() is more accurate than one might expect! 
     However, I made one small change:
     Often NEITHER solution is correct, and BOTH should be discarded.  The code
     was assuming that the second solution was good, if the first solution was 
     bad.

    tjn-07-2000 

  */


  

  fXCE1 = fDot3D(fX, fE1);
  fXCE3 = fDot3D(fX, m_afRotVec);
  fE1S0 = fDot3D(fE1, m_afS0);
  fE2S0 = fDot3D(fE2, m_afS0);
  fE3S0 = fDot3D(m_afRotVec, m_afS0);

  fCEA = fXCE1 * fE1S0;
  fCEB = fXCE1 * fE2S0;
  fCEC = 0.5f * fXCE1* fXCE1  +  0.5f * fXCE3 * fXCE3  -  fXCE3 * fE3S0;

  fCEACEB = fCEA * fCEA  +  fCEB * fCEB;

// Skip if reciprocal lattice point (X) lies on the rotation axis

  if (fCEACEB <= 0.0) return (4);
  fARG1 = fCEC / sqrt(fCEACEB);

  if (fARG1 > 1.0)
    fARG1 = 1.0;
  else if (fARG1 < -1.0)
    fARG1 = -1.0;



// Equation for fReflnRot0 has two solutions.  The one of
// interest would lie between fRotStartRad and fRotEndRad.

  fT1 = acos(fARG1);
  fT2 = atan2(fCEB, fCEA);

  fROTCA = fT1  +  fT2;
  fROTCB = -fT1 + fT2;
  
  if ((min(ABS(fROTCA + fRotStartRad - fRotStartRad),ABS(fROTCA + fRotStartRad - fRotEndRad))) <
      (min(ABS(fROTCB + fRotStartRad - fRotStartRad),ABS(fROTCB + fRotStartRad - fRotEndRad)))) 
      fROTC = fROTCA;
  else
      fROTC = fROTCB;

  // fROTC = fT1  +  fT2;
  // m_fReflnRot0  = fRotStartRad + fROTC;
  //if ( (m_fReflnRot0 < (fRotStartRad - m_fReflnWidthMaxRad)) ||
  //     (m_fReflnRot0 >= (fRotEndRad + m_fReflnWidthMaxRad)) )
  //  fROTC = -fT1  +  fT2;
  

  m_fReflnRot0  = fRotStartRad + fROTC;

  
  // Don't exclude based on large length.
  /*
  if ( (m_fReflnRot0 < (fRotStartRad - m_fReflnWidthMaxRad)) ||
       (m_fReflnRot0 >= (fRotEndRad + m_fReflnWidthMaxRad)) )
       return (7);
   */
  

  // Calculate coordinates XR when spot lies on sphere

  vConvRotVec3DMat3D(fROTC/fRADIANS_PER_DEGREE, m_afRotVec, &fROTC1[0][0]);
  vMulMat3DVec3D(fROTC1, &fX[0], m_afXR);    // m_afXR is member variable


    // Calculate the rotated shift vector.
  vMulMat3DVec3D(m_afGonMatrix,m_afRecipShift,fRS);
  vConvRotVec3DMat3D((fROTC + fRotStartRad)/fRADIANS_PER_DEGREE, m_afRotVec, &fROTC1[0][0]);
  vMulMat3DVec3D(fROTC1,fRS,m_afRecipShiftRotated);


  // Calculate central ROT-value, ROT range, and ROT start
  //
  // Determining reflection range from Eq. VII.17 of Greenhough & Helliwell
  //  modified for inclined geometry.

  // ROTrange = Lorentz factor * ( effective mosaic spread * Dstar *
  //                    cos theta + spectral dispersion * Dstar * sin theta)
  //
  // Formula for Lorentz-factor taken from J.R. Milch and T.C. Minor,
  //  J. Appl. Cryst. (1974) 7, 502 - 505
  //                         -
  //  Lorentz factor = | 1. / ( XR . [E3 X S0] ) |
  //  fLorentz is the Lorentz factor (and not its inverse)
  //  thus true intensity Itrue = Iobs/Lorentz

  fT1 = fDot3D(m_afXR, m_afE3crossS0);
  if (fT1 == 0.0) return (5);
  m_fLorentz = fabs(1.0f / fT1);

  // Now for the reflecting range
  // ROTrange = [eps*dstar*cos(theta) +
  //       dispersion*tan(theta)*dstar*cos(theta)]*Lorentz_factor
  //
  // EPS is the full rocking range (in radians)
  // and dstar*cos(theta) = SQRT(dstar**2 - 0.25*dstar**4)
  // lambda = 2*d*cos(2*theta) = 

  fDSCOTH =  sqrt(m_fReflnDstarSq  -  0.25f * m_fReflnDstarSq * m_fReflnDstarSq);

  // Note that REPS in degrees, EPS (= m_fCrysMosaicityRad) in radians
  //
  // We can write:
  //   fReflnWidth = ( (EPS + (0.5 * m_fSourceDispersion * m_fReflnDstarSq/DSCOTH) )
  //                   * m_fLorentz * DSCOTH
  //   fReflnWidth = ( (EPS + ( m_fSourceDispersion * tan(theta) )
  //                   * m_fLorentz * DSCOTH
  //
  // where on the equatorial plane DSCOTH * m_fLorentz = 1.0 and
  // as m_fSourceDispersion*tan(theta) is small then fReflnWidth ~ EPS.
  // At 2 Angstrom using Copper Kalpha theta=22.67 degrees
  // and m_fSourceDispersion=0.0025 so that
  //     m_fSourceDispersion *tan(theta)=0.06 degrees
  //

  //fDSCOTH = fabs(2.0 - m_fReflnDstarSq ) * 0.5;

  // Must call this here because we need the pixel coordinates if have detector dependent mosaicity information available on the plate.
  // if (nCalcDetCoords())
  //    return 1;

  
  // Use this line for detector dependent mosaicity:  
  // fCrysMosaicityRad = m_ppoCrystal[m_nWhichCrystal]->fGetMosaicity(m_fPxNorm0,m_fPxNorm1)  * fRADIANS_PER_DEGREE;
  // Otherwise, use this line:

  m_ppoCrystal[m_nWhichCrystal]->vGetMosaicity(a3fCrysMosaicity);
  fCrysMosaicityRad = a3fCrysMosaicity[0] * fRADIANS_PER_DEGREE;
  m_fReflnWidth = m_fLorentz*(fCrysMosaicityRad*fDSCOTH
      + a3fCrysMosaicity[1]*sqrt(m_fReflnDstarSq-fDSCOTH*fDSCOTH))  + a3fCrysMosaicity[2]*fRADIANS_PER_DEGREE;
  
  m_fReflnStart = m_fReflnRot0  - (m_fReflnWidth * 0.5f);
  m_fReflnEnd   = m_fReflnStart + m_fReflnWidth;

  if (m_bCalcMosaicityLinearCoeffs) {
      m_fMosCoeffA = m_fLorentz*(fDSCOTH)*Gs_dRADIANS_PER_DEGREE;
      m_fMosCoeffB = m_fLorentz*(m_fSourceDispersion*sqrt(m_fReflnDstarSq-fDSCOTH*fDSCOTH));
  } else if (m_bCalcMinMosaicityForRefln) {
      

      // Perhaps the centroid is already in the region? If so, this reflection will always
      // get predicted regardless of mosaicity.
      if ( (m_fReflnRot0 >= (fRotStartRad)) &&
          (m_fReflnRot0 <= (fRotEndRad)))
          m_fMosCoeffA = 0.0;
      else {
          f0 = 2.0*min(ABS(m_fReflnRot0 - fRotStartRad),ABS(m_fReflnRot0 - fRotEndRad));
          // We need to have a reflection who's width is at least f0.
          f0 -= m_fLorentz*(a3fCrysMosaicity[1]*sqrt(m_fReflnDstarSq-fDSCOTH*fDSCOTH))  + a3fCrysMosaicity[2]*fRADIANS_PER_DEGREE;
          f0 /= m_fLorentz*fDSCOTH;
          if (f0>0.0)
              m_fMosCoeffA = f0/Gs_dRADIANS_PER_DEGREE;
          else
              m_fMosCoeffA = 10.0;
      };
  };
  
  // Return error if reflection with is bogusly large

  if (m_fReflnWidth > 90.0 * fRADIANS_PER_DEGREE) 
      return (6);

  return (0);
}

/*  This function is now disabled.  It's original purpose was to compute the change in
    spot centroid as a function of image.  This contribution is usually quite small.
    It was removed because it was causing problems when integrating wide sliced images.

*/

int Cpredict::nCalcDetCoords(void)
{
  // Calculate the detector coordinates of a reflection, given
  // dimensionless reciprocal lattice coordinates at diffracting
  // condition, detector info and fS0;

  //  Calculates spot coordinates on the detector and tests
  //  whether spot is within detector limits.
  //  Return 0 on success, non-zero on failure

  //  Calculates image drift parameters if the user has requested them.

  //  This subroutine created using code written by Albrecht
  //  Messerschmidt.    November 1987

  //
  // Uses member variables:
  // m_afXR[3]           Dimensionless reciprocal lattice coordinate of this
  //                      reflection in the laboratory frame.
  // fDN[3]           Vector defining detector normal in lab frame.
  // ffInvDD[3][3]    Matrix used in calculation of detector coords
  //
  // Calculates member variables:
  //
  //  fPx0, fPx1  Reflection detector pixel coordinates in FAST, SLOW directions
  //  fXmm, fYmm  Reflection mm coordinates in lab X, Y directions
  //  nNonunfFlag 0 if pixel has OK non-uniformity, undefined if pixel is not on
  //                the detector
  //
  //  Local variables

  double fXRdotDN;
  double fTS[3], fVDCC[3];

  // Convert the coordinates of the vector XR in reflection
  // position in diffractometer coordinate system into
  // mm-detector-coordinates.

  // fVDCC - calculated detector coordinates on inclined
  //         detector.
  //
  //  Now we need  TS = (XR - S0) * DIST**2 / (XR.DN - S0.DN)
  //               --    --   --               -- --   -- --
  // See Proceedings Phase III, pp. 57-64

  fXRdotDN = fDot3D(m_afXR, m_afDN);

  // Added to make sure that 'ghost reflections' don't get predicted.  This was a problem
  // with large 2theta angle experiments.
  if ((m_ppoDetector[m_nWhichDetector]->m_poSpatial->eSpatialType()==eSpatial_flat_geometry) && (fXRdotDN < m_fS0dotDN))
    return  1;

  vSubVec3DVec3D(m_afXR, m_afS0, m_afXRminusS0);   // Also used by oblique incidence calc
  
  if( 0 != m_ppoDetector[m_nWhichDetector]->nScaleSOnPlateItr(m_afXRminusS0, 
                                                        fTS, 
                                                        m_afRecipShiftRotated) )
    return 1;

    // The calculated detector coordinates VDCC are derived
    // by    VDCC = m_afInvDD * TS
    //       ----   =========   --
    vMulMat3DVec3D(m_afInvDD, fTS, fVDCC);

    m_fXmm = fVDCC[0];
    m_fYmm = fVDCC[1];
    m_fZRelativeCoordinate = fVDCC[2];

    // Convert mm to pixels to see if reflection falls on the detector:
    float     f0 = 0.0f;
    float     f1 = 0.0f;

    if( 0 != m_ppoDetector[m_nWhichDetector]->m_poSpatial->nMMtoPixel(m_fXmm, m_fYmm, m_fZRelativeCoordinate, &f0,&f1) )
        return 1;

    m_fPx0 = f0;
    m_fPx1 = f1;

    m_fPxNorm0 = f0/m_ppoDetector[m_nWhichDetector]->m_poSpatial->nGetMaxPixel(0);
    m_fPxNorm1 = f1/m_ppoDetector[m_nWhichDetector]->m_poSpatial->nGetMaxPixel(1);

    if( m_bDetectorIsCylindrical )
        m_ppoDetector[m_nWhichDetector]->vGetLocalNormal(m_a3fLocalDN);
  
    return 0;
}

int Cpredict::nCalcDriftDetCoords(void)
{
// Calculate the image drift if the user has requested it.
  // This is used by the integration to calculate shifted centroids without calling predict twice.
    
    
    // Try positive and negative change.
    double f0;
    double a3x3fTemp[3][3];
    double a3fTemp[3];
    int nStat;
    float fPx0_1,fPx1_1;
    float fPx0_2,fPx1_2;
    double a3fXR[3];

    m_fPx0Drift = -999;
    m_fPx1Drift = -999;
    
    f0 = m_fReflnEnd - m_fReflnRot0;
    vCopyVec3D(m_afXR,a3fXR);
    vConvRotVec3DMat3D(f0/Gs_dRADIANS_PER_DEGREE, m_afRotVec, &a3x3fTemp[0][0]);
    vMulMat3DVec3D(a3x3fTemp,m_afXR,a3fTemp);
    vCopyVec3D(a3fTemp,m_afXR);
    nStat = nCalcDetCoords();
    vCopyVec3D(a3fXR,m_afXR);
    if (nStat)
        return 1;
    fPx0_1 = m_fPx0;
    fPx1_1 = m_fPx1;

    f0 = m_fReflnStart - m_fReflnRot0;
    vConvRotVec3DMat3D(f0/Gs_dRADIANS_PER_DEGREE, m_afRotVec, &a3x3fTemp[0][0]);
    vMulMat3DVec3D(a3x3fTemp,m_afXR,a3fTemp);
    vCopyVec3D(a3fTemp,m_afXR);
    nStat = nCalcDetCoords();
    vCopyVec3D(a3fXR,m_afXR);
    if (nStat)
        return 1;
    fPx0_2 = m_fPx0;
    fPx1_2 = m_fPx1;

    
    m_fPx0Drift = (fPx0_1 - fPx0_2)/(m_fReflnEnd - m_fReflnStart)*Gs_dRADIANS_PER_DEGREE;
    m_fPx1Drift = (fPx1_1 - fPx1_2)/(m_fReflnEnd - m_fReflnStart)*Gs_dRADIANS_PER_DEGREE;
    return (0);
};


int Cpredict::nPredictWholeDataSet(double fRotStart,double fRotEnd)
{
    int nHmax,nKmax,nLmax;
    int nHmin,nKmin,nLmin;
    int nH,nK,nL;

    double f0;
    float a3fFloatBuf[3];
    float a3x3fFloatBuf[3][3];
    int nx;
    double fRotStartRad,fRotEndRad,fRotStartOffsetRad;
    int nMaxSize;
    double fRotTol;
    double a3fMosCoeffs[3];


    m_fPrePredictRotWidth = fDetermineRockingWidth(0.95);
    
    nHmax = (int) ( m_ppoCrystal[0]->fGetCell(0) / m_fResolutionMax) + 1;
    nKmax = (int) ( m_ppoCrystal[0]->fGetCell(1) / m_fResolutionMax) + 1;
    nLmax = (int) ( m_ppoCrystal[0]->fGetCell(2) / m_fResolutionMax) + 1;
    nHmin = -nHmax;
    nKmin = -nKmax;
    nLmin = -nLmax;

    m_ppoCrystal[0]->vGetMosaicity(a3fMosCoeffs);
    fRotTol = m_fPrePredictRotWidth*a3fMosCoeffs[0];


    nMaxSize =   (int) (((nHmax - nHmin) + 1)        // Use nStat as temp variable here
        * ((nKmax - nKmin) + 1)
        * ((nLmax - nLmin) + 1) * 2.0*(fRotEnd - fRotStart + 2.0*fRotTol)/360.0);
    -m_afPrePredictRotMid;
    -m_anPrePredictPackedHKL;

    m_afPrePredictRotMid.setsize(nMaxSize);
    m_anPrePredictPackedHKL.setsize(nMaxSize);
    float* pfPrePredictRotMid = & m_afPrePredictRotMid[0];
    int* pnPrePredictPackedHKL = & m_anPrePredictPackedHKL[0];

    fRotStartRad = (fRotStart - fRotTol)*Gs_dRADIANS_PER_DEGREE;
    fRotEndRad = (fRotEnd + fRotTol)*Gs_dRADIANS_PER_DEGREE;
    fRotStartOffsetRad = 0.0;

    int nEntry;
    int nRotAxis;
    double fWavelength;
    double fDstarSq;
    double a3fHKL[3];
    double a3fXRC[3];
    double a3x3fGUBMat[3][3];
    double a3x3fUBMat[3][3];
    double a3x3fGMat[3][3];
    double a3fE1[3];
    double a3fE2[3];
    double a3fE3[3];
    double a3fS0[3];
    double fPercentCompleted,fNewPercentCompleted;

    // Compute Goniometer * U * B matrix.
    nRotAxis = m_poCrysGonio->nGetNum(m_poRotation->sGetName());
    m_poCrysGonio->vCalcGetRotMatrix(& a3x3fFloatBuf[0][0],nRotAxis); 
    vCopyMat3D(&a3x3fFloatBuf[0][0],&a3x3fGMat[0][0]);
    m_ppoCrystal[0]->nCalcOrientMatrix(fWavelength = m_poSource->fGetWavelength());
    m_ppoCrystal[0]->vGetOrientMatrix(&a3x3fUBMat[0][0]);
    vMulMat3DMat3D(a3x3fGMat,a3x3fUBMat,a3x3fGUBMat);
    // Compute rotational vector.
    if (m_poCrysGonio->nGetRotVector(nRotAxis,&a3fFloatBuf[0])) 
        return 1;
    vCopyVec3D(&a3fFloatBuf[0],&a3fE3[0]);
    m_poSource->vCalcGetS0(&a3fS0[0]);

    printf("Computing all reflections for entire data set from %.2lf to %.2lf degrees.\n",
        fRotStart,fRotEnd);
    fPercentCompleted = -1.0;


    nEntry = 0;
    for (nH = nHmin; nH <= nHmax; nH++)
    {
        fNewPercentCompleted = 5.0*floor((100.0*(nH - nHmin))/(nHmax - nHmin)/5.0);
        if (fNewPercentCompleted > fPercentCompleted) {
            printf("\r%.2lf percent done.",fNewPercentCompleted);
            fPercentCompleted = fNewPercentCompleted;
        };

        a3fHKL[0] = float(nH);
        a3fHKL[1] = 0.0;
        a3fHKL[2] = 0.0;
        for (nK = nKmin; nK <= nKmax; nK++)
        {
            a3fHKL[1] = float(nK);
            for (nL = nLmin; nL <= nLmax; nL++)
            {
                a3fHKL[2] = float(nL);
                
                vMulMat3DVec3D(a3x3fGUBMat, a3fHKL, a3fXRC);
                
                fDstarSq = fDot3D(a3fXRC, a3fXRC);
                f0 = fWavelength/sqrt(max(0.0000001,fDstarSq));
                
                
                if ((f0 >= m_fResolutionMax) && (f0 < m_fResolutionMin)) {
                    
                    int nPackedHKL;

                    nPackedHKL = ((nH + 512) << 21)
                        | ((nK + 512) << 11)
                        | ((nL + 512) <<  1);

                    double fXCE1;
                    double fXCE3;
                    double fE1S0;
                    double fE2S0;
                    double fE3S0;
                    double fCEA;
                    double fCEB;
                    double fCEC;
                    double fCEACEB;
                    double fT1;
                    double fT2;
                    double fDE2;
                    double fROTC;
                    double fARG1;

                    vCross3D(a3fE3, a3fXRC, a3fE2);
                    fDE2 = fNormVec3D(a3fE2);
                    if (fDE2 <= 0.0) 
                        continue;
                    //  E1 = [E2 x fRotVec]
                    vCross3D(a3fE2, a3fE3, a3fE1);
                    fXCE1 = fDot3D(a3fXRC, a3fE1);
                    fXCE3 = fDot3D(a3fXRC, m_afRotVec);
                    fE1S0 = fDot3D(a3fE1, a3fS0);
                    fE2S0 = fDot3D(a3fE2, a3fS0);
                    fE3S0 = fDot3D(m_afRotVec, a3fS0);
                    
                    fCEA = fXCE1 * fE1S0;
                    fCEB = fXCE1 * fE2S0;
                    fCEC = 0.5f * fXCE1* fXCE1  +  0.5f * fXCE3 * fXCE3  -  fXCE3 * fE3S0;
                    
                    fCEACEB = fCEA * fCEA  +  fCEB * fCEB;
                    
                    // Skip if reciprocal lattice point (X) lies on the rotation axis
                    
                    if (fCEACEB <= 0.0) return (4);
                    fARG1 = fCEC / sqrt(fCEACEB);
                    
                    if (fARG1 > 1.0)
                        fARG1 = 1.0;
                    else if (fARG1 < -1.0)
                        fARG1 = -1.0;
                    
                    fT1 = acos(fARG1);
                    fT2 = atan2(fCEB, fCEA);
                    
                    if (nEntry + 4 >= m_afPrePredictRotMid.size()) {
                        m_afPrePredictRotMid.setsize((int) (m_afPrePredictRotMid.size()*1.5));
                        m_anPrePredictPackedHKL.setsize((int) (m_anPrePredictPackedHKL.size()*1.5));
                        pfPrePredictRotMid = & m_afPrePredictRotMid[0];
                        pnPrePredictPackedHKL = & m_anPrePredictPackedHKL[0];
                    };
                    
                    // Equation for fReflnRot0 has up to four solutions.  The one of
                    // interest would lie between fRotStartRad and fRotEndRad.
                    for (nx = 0; nx < 2; nx++) {
                        if (nx == 0) 
                            fROTC = fRotStartOffsetRad + fT1  +  fT2;
                        else
                            fROTC = fRotStartOffsetRad +  -fT1 + fT2;
                        if ((fROTC + Gs_dPI*2.0 <= fRotEndRad) && (fROTC + Gs_dPI*2.0 >= fRotStartRad)) {
                            pfPrePredictRotMid[nEntry] = (fROTC + Gs_dPI*2.0)/Gs_dRADIANS_PER_DEGREE;
                            pnPrePredictPackedHKL[nEntry] = nPackedHKL;
                            nEntry++;
                        };
                        
                        if ((fROTC - Gs_dPI*2.0 <= fRotEndRad) && (fROTC - Gs_dPI*2.0 >= fRotStartRad)) {
                            pfPrePredictRotMid[nEntry] = (fROTC - Gs_dPI*2.0)/Gs_dRADIANS_PER_DEGREE;
                            pnPrePredictPackedHKL[nEntry] = nPackedHKL;
                            nEntry++;
                        };
                        if ((fROTC >= fRotStartRad) && (fROTC <= fRotEndRad)) {
                            pfPrePredictRotMid[nEntry] = (fROTC)/Gs_dRADIANS_PER_DEGREE;
                            pnPrePredictPackedHKL[nEntry] = nPackedHKL;
                            nEntry++;
                        };
                        

                        /*
                        if (nEntry >= m_afPrePredictRotMid.size()) {
                            m_afPrePredictRotMid.setsize(min(nMaxSize,(nEntry + 1)*2));
                            m_anPrePredictPackedHKL.setsize(min(nMaxSize,(nEntry + 1)*2));
                        };
                        
                        m_afPrePredictRotMid[nEntry] = fROTC;
                        m_anPrePredictPackedHKL[nEntry] = nPackedHKL;
                        nEntry++;
                        */

                    };
                };
            };
        };
    };

    // A special purpose implementation of quick sort is used placed here.  We don't want to allocate another array because we
    // might have an extremely large # of data points at this stage.

    -g_anSwapSizes + sizeof(float) + sizeof(int);
    -g_apvSwapPointers + (void*) &m_afPrePredictRotMid[0] + (void*) &m_anPrePredictPackedHKL[0];
    printf("\nSorting predicted reflections.");
    qsortswap((void*) &m_afPrePredictRotMid[0],nEntry,sizeof(float),float_cmp,qsort_swap_arrays);
    
    printf("\nDone.\n");
    
    cout << flush;

    m_afPrePredictRotMid.setsize(nEntry);
    m_anPrePredictPackedHKL.setsize(nEntry);


    return 0;
};


double Cpredict::fDetermineRockingWidth(double fFractionUsed) {
        // Determine the 'mosaicity spread' over the surface of the plate.
    // We sample the plate at even positions, and project an S vector back to diffracting condition.
    // The effective width at a 1.0 degree mosaicity is calculated.
    // We want to find the width which is >= 95 percent of our sampled spots.
    // This then gives us an estimate of how far 'out' to carry the prediction calculations.

    double f0;
    int nPix0,nPix1;
    int nIndex0,nIndex1;
    double fStep0,fStep1;
    double a3fVDC[3];
    float  a3fVDC1[3];
    double a3fS[3];
    double a3fX[3];
    double fT1;
    double fCrysMosaicityRad;
    double fDSCOTH;
    double a3fCrysMosaicity[3];
    itr<double> afWidthArray;

    fStep0 = m_ppoDetector[0]->m_a2nDim[0]/10;
    fStep1 = m_ppoDetector[0]->m_a2nDim[0]/10;
    fCrysMosaicityRad = 1.0*Gs_dRADIANS_PER_DEGREE;
    m_ppoCrystal[m_nWhichCrystal]->vGetMosaicity(a3fCrysMosaicity);
        
    for (nIndex0 = 0; 1; nIndex0++) {
        nPix0 = (int) (nIndex0*fStep0);
        if (nPix0 >= m_ppoDetector[0]->m_a2nDim[0])
            break;
        for (nIndex1 = 0; 1; nIndex1++) {
            
            
            nPix1 = (int) (nIndex1*fStep1);
            if (nPix1 >= m_ppoDetector[0]->m_a2nDim[1])
                break;
            if (m_ppoDetector[0]->m_poSpatial->nPixeltoMM(
                nPix0,nPix1,&a3fVDC1[0],&a3fVDC1[1],&a3fVDC1[2]))
                continue;
            vCopyVec3D(a3fVDC1,a3fVDC);
            vMulMat3DVec3D(m_afDD, a3fVDC, a3fS); // Scattered beam wavevector
            fNormVec3D(a3fS);
            vAddVec3DVec3D(m_afS0, a3fS,a3fX);            
            m_fReflnDstarSq = fDot3D(a3fX, a3fX);
            f0 = m_fWavelength/sqrt(max(0.0000001,m_fReflnDstarSq));                        
            
            if ((f0 >= m_fResolutionMax) && (f0 < m_fResolutionMin)) {
                
                
                fDSCOTH =  sqrt(m_fReflnDstarSq  -  0.25f * m_fReflnDstarSq * m_fReflnDstarSq);
                
                fT1 = fDot3D(a3fX, m_afE3crossS0);
                if (fT1 == 0.0) 
                    continue;
                m_fLorentz = fabs(1.0f / fT1);
                if (m_fLorentz < 10.0) {
                    
                    m_fReflnWidth = m_fLorentz*(fCrysMosaicityRad*fDSCOTH
                        + a3fCrysMosaicity[1]*sqrt(m_fReflnDstarSq-fDSCOTH*fDSCOTH))  + a3fCrysMosaicity[2]*fRADIANS_PER_DEGREE;
                    afWidthArray + (m_fReflnWidth/Gs_dRADIANS_PER_DEGREE);
                };
            };
        };
    };
    if (!afWidthArray.size())
        return 1;
    qsort(&afWidthArray[0],afWidthArray.size(),sizeof(double),double_cmp);
    return (afWidthArray[(int) (afWidthArray.size()*fFractionUsed)]);
};




int Cpredict::nPredictReflns(const double fRotStart, const double fRotEnd)
{
  // This is the inner loop of predicting reflections.  All matrices
  // for the Source, Crystal, Crystal Goniometer and Detector have been setup
  // This is from the PRDWTX routine in MADNES with many people to thank
  // for: Peter Brick, Andrew Leslie, Jean Claude Thierry, Alan Wonacott,
  // Albrecht Messerschmidt, Jim Pflugrath, George Reeke and Bob Sweet.

  int i, j, k, l, ip, iq, ir, iloop;                 // Loop counters
  double fP0[2][3][3], fP[2][3][3];
  double fPP[2][4][3], fT[4][4];

  Crefln oRefln(m_poReflnlist);   // Create a reflection object

  int   nIPLB,  nIPLE;
  double fQminE, fQmaxE;
  int   nIQMIN, nIQMAX, nIQLE, nIQLB;
  int   nIRMIN, nIRMAX;
  int   nIER;
  int   nNRT;
  int   nJ1, nJ2;

  double fS1P2R2;
  double fRT[2], fRS1[2], fRS0[2], fRU2[2], fRU1[2], fRU0[2], fQU[2];
  double fP1H[2], fP2H[2], fP3H[2], fRLB[2], fRLE[2];
  double fP1H0[3], fPX0[3], fX12[3];
  double fSca1;
  double fStartRad, fEndRad;

  int   nNRTA[2];
  int   nHKL[3];
  int   nIRLE[2], nIRLB[2];
  int   nLOOP;
  double fFP, fFR, fFQ;
  double fFP2, fS5P, fS6P, fPX, fPY, fPZ;
  double fT1, fT2;                            // Temporary variables
  double fS1, fS2, fS3, fS4, fS5, fS6;
  double fPA, fPB, fPC, fPD;

  double fRotS, fRotE;

  // Allow inputs in any order, these should follow order in ::nSetupRotation

  fRotS = min(fRotStart, fRotEnd);
  fRotE = max(fRotStart, fRotEnd);

  // Set detector number field of reflection template that will be inserted
  // into the resultant poReflnlist.

  oRefln.vSetField(m_poReflnlist->m_nFI_nDetNum, m_nWhichDetector);

  // Compute start and end of rotation range in radians:

  fStartRad = fRotS * fRADIANS_PER_DEGREE;

  fEndRad   = fRotE * fRADIANS_PER_DEGREE;

  // Re-order rotation matrices for efficient looping

  for (k = 0; k < 3; k++) {
    for (j = 0; j < 3; j++) {
       fP[0][k][j] = m_afRotStartMatrix[m_anSortedAxes[k]][j];
      fP0[0][k][j] = m_afRotStartMatrix[m_anSortedAxes[k]][j];
       fP[1][k][j] = m_afRotPastEndMatrix[m_anSortedAxes[k]][j];
      fP0[1][k][j] = m_afRotEndMatrix[m_anSortedAxes[k]][j];
    }
  }

  // Set fPP (Pprime) as a 3x4(x2) matrix. Translation to move origin
  // to center of the Ewald sphere.  Remember that m_afS0 is a unit vector
  // antiparallel to the X-ray beam.

  for (i = 0; i < 2; i++) {
    for (k = 0; k < 3; k++) {
      fPP[i][3][k] = -m_afS0[k];     // JWP: Be sure to check sign of m_afS0!
      for (j = 0; j < 3; j++) {
        fPP[i][j][k] = fP[i][j][k];
      }
    }
  }

  //  Calculate the reciprocal metric tensor fT = fPP(trans). fPP
  //  The 3*3 portion of matrix is same at beginning and end of rotation
  //  Note that RT, RS1, RU2 and QU don't depend on L and thus don't need
  //  to be calculated twice.
  //  fPP is used here to calculate fT but is not used again

  for (l = 0; l < 2; l++) {
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        fT1 = 0.0;
        for (k = 0; k < 3; k++) {
          fT1 = fT1 + fPP[l][i][k] * fPP[l][j][k];
        }  // endl of k loop
        fT[j][i] = fT1;
      }  // end of j loop
    }  // end of i loop

    fRT[l]  =         fT[2][1] * fT[2][1]  -  fT[2][2] * fT[1][1];
    fRS1[l] =         fT[2][0] * fT[2][1]  -  fT[2][2] * fT[1][0];
    fRS0[l] =         fT[2][1] * fT[3][2]  -  fT[3][1] * fT[2][2];
    fRU2[l] =         fT[2][0] * fT[2][0]  -  fT[0][0] * fT[2][2];
    fRU1[l] = 2.0f * (fT[2][0] * fT[3][2]  -  fT[3][0] * fT[2][2]);
    fRU0[l] =         fT[3][2] * fT[3][2];
    fQU[l]  =         fT[2][2];

  }  // end of l loop

  fS1 = fT[0][0];
  fS2 = fT[1][1];
  fS3 = fT[2][2];
  fS4 = fT[1][2] + fT[2][1];
  fS5 = fT[2][0] + fT[0][2];
  fS6 = fT[1][0] + fT[0][1];

  fPA = fS4 * fS4  -  4.0f * fS2 * fS3;
  fPB = fS4 * fS5  -  2.0f * fS3 * fS6;
  fPC = fS5 * fS5  -  4.0f * fS1 * fS3;
  fPD = 4.0f * fS3 * m_fDstarSqMax;

  double fPmin  = (double)  1.E20;
  double fPmax  = (double) -1.E20;
  double fPminE = (double)  1.E20;
  double fPmaxE = (double) -1.E20;
  double fPMN[2], fPMX[2];
  double fSinTh = 0.5 * sqrt(m_fDstarSqMax);
  double fSinSq = 0.25 * m_fDstarSqMax;
  double fSin2T = sin(2.0 * atan2(fSinTh, sqrt(1.0-fSinSq)));
  double fVec[3], fVecX;
  double fAvec, fAvecX;
  double fA1, fPX1, fPX2, fP1, fP2;
  double fIsign;


  // Ewald sphere limit and limiting sphere both computed here.
  // Note that in order to allow for a general X-ray beam direction,
  // the algorithm for calculating the limiting sphere index limits
  // has been changed from that proposed by Reeke, and is now
  // essentially the same as the algorithm used for the Ewald sphere.

  for (i = 0; i < 2; i++) {

    // Determine if the reciprocal lattice vector closest to the X-ray
    // beam is parallel or antiparallel.  Remember S0 is antiparallel to
    // X-ray beam!

    fSca1 = fDot3D(&fP[i][0][0], m_afS0);
    fIsign = (double) (-(int) (1.1 * fSca1 / fabs(fSca1))); // Transfer opp.
                                                //  sign of fSca1 for use later
    for (j = 0; j < 3; j++) {
      nJ1 =  (j+1) % 3;      // ftn: 2, 3, 1    C: 1, 2, 0
      nJ2 =  (j+2) % 3;      // ftn: 3, 1, 2 C: 2, 0, 1
      fVec[j] = fP[i][1][nJ1] * fP[i][2][nJ2]  -  fP[i][2][nJ1] * fP[i][1][nJ2];
    }

    // fVec is the vector normal to the reciprocal lattice planes with
    // ip loop index zero.

    fT1 = 0.0;
    fT2 = 0.0;
    for (j = 0; j < 3; j++) {
      fT1 = fT1  +  fVec[j] * fVec[j];
      fT2 = fT2  +  fP[i][0][j] * fVec[j];
    }

    // fT2 is the volume of the reciprocal space unit cell. Because of
    // the permutation of axes (by IAX), this may be negative.

    fAvec = sqrt(fT1);
    fT2 = 1.0f / fabs(fT2);

    // Limiting sphere limits

    fVecX  = fDot3D(m_afS0, fVec);
    fAvecX =  fabs(fVecX);
    fA1 = fT1  -  fAvecX * fAvecX;
    if (fA1 > 0.0)
      fA1 = sqrt(fA1);
    else
      fA1 = 0.0;

    fPX1 = -fIsign * fT2 * (2.f * fSinSq * fAvecX  +  fSin2T * fA1);
    fPX2 = -fIsign * fT2 * (2.f * fSinSq * fAvecX  -  fSin2T * fA1);
//    fP1  = fPX1 <? fPX2;  // Min of fPX1, fPX2
//    fP2  = fPX1 >? fPX2;  // Max of fPX1, fPX2
    fP1  = min(fPX1 , fPX2);  // Min of fPX1, fPX2
    fP2  = max(fPX1 , fPX2);  // Max of fPX1, fPX2

    // Ewald sphere limits

    fT1    = -fIsign * fT2 * (fAvec + fAvecX);
    fT2    =  fIsign * fT2 * (fAvec - fAvecX);
//    fPminE =  fT1 <? fT2;
//    fPmaxE =  fT1 >? fT2;
    fPminE =  min(fT1 , fT2);
    fPmaxE =  max(fT1 , fT2);

    // Select limits on basis of fIsign

    if (fIsign > 0.0) {
//      fPMN[i] = fP1 >? fPminE;
//      fPMX[i] = fP2 >? fPmaxE;
      fPMN[i] = max(fP1 , fPminE);
      fPMX[i] = max(fP2 , fPmaxE);
    }
    else {
//      fPMN[i] = fP1 <? fPminE;
//      fPMX[i] = fP2 <? fPmaxE;
      fPMN[i] = min(fP1 , fPminE);
      fPMX[i] = min(fP2 , fPmaxE);
    }
  } // end of i loop

  // Select overall limits

//  fPmin = fPMN[0] <? fPMN[1];
//  fPmax = fPMX[0] >? fPMX[1];
  fPmin = min(fPMN[0] , fPMN[1]);
  fPmax = max(fPMX[0] , fPMX[1]);

  nIPLB = nIntMinusOne(fPmin);
  nIPLE = nIntPlusOne(fPmax);

  for (ip = nIPLB; ip <= nIPLE; ip++) {

    fFP  = (double) ip;
    fFP2 = fFP * fFP;
    fS5P = fS5 * fFP;
    fS6P = fS6 * fFP;
    fS1P2R2 = fFP2 * fS1  -  m_fDstarSqMax;

    // Calculate looping limits for iq limiting sphere

    nNRT = nSolveQuadratic(fPA, fFP * fPB, (fFP2 * fPC  +  fPD), fVec);

    if (nNRT == 1) {
      nIQMIN = nIntMinusOne(fVec[0]);
      nIQMAX = nIntPlusOne(fVec[0]);
    }
    else if (nNRT > 1) {
//      nIQMIN = nIntMinusOne(fVec[0] <? fVec[1]);
//      nIQMAX = nIntPlusOne(fVec[0] >? fVec[1]);
      nIQMIN = nIntMinusOne(min(fVec[0] , fVec[1]));
      nIQMAX = nIntPlusOne(max(fVec[0] , fVec[1]));
    }
    if (nNRT >= 1) {

      // Ewald sphere

      fQminE = (double)  1.E20;
      fQmaxE = (double) -1.E20;

      for (i = 0; i < 3; i++) {
        fP1H0[i] = fP0[0][0][i] * fFP;
      }
      for (i = 0; i < 2; i++) {
        fP1H[i] = fP[i][0][0] * fFP  -  m_afS0[0];
        fP2H[i] = fP[i][0][1] * fFP  -  m_afS0[1];
        fP3H[i] = fP[i][0][2] * fFP  -  m_afS0[2];
        fT1     = fRU0[i]  +  fFP * (fRU1[i]  +  fFP * fRU2[i]);

        nNRTA[i] = nSolveQuadratic(fRT[i], (fFP*fRS1[i] + fRS0[i]), fT1, fVec);

        if (nNRTA[i] == 1) {
//          fQminE = fQminE <? fVec[0];
//          fQmaxE = fQmaxE >? fVec[0];
          fQminE = min(fQminE , fVec[0]);
          fQmaxE = max(fQmaxE , fVec[0]);
        }
        else if (nNRTA[i] > 1) {
//          fQminE = fQminE <? fVec[0] <? fVec[1];
//          fQmaxE = fQmaxE >? fVec[0] >? fVec[1];
          fQminE = min(fQminE , min(fVec[0] , fVec[1]));
          fQmaxE = max(fQmaxE , max(fVec[0] , fVec[1]));
        }
      } // end of i loop

      if ( (nNRTA[0] != 0) || (nNRTA[1] != 0) ) {

        // Choose most restrictive limits

        nIQLB = nIntMinusOne(fQminE);
        nIQLE = nIntPlusOne(fQmaxE);

        nIER = nOutsideLimits(nIQMIN, nIQMAX, &nIQLB, &nIQLE);
        if (nIER == 0) {

          // Loop over IQ

          for (iq = nIQLB; iq <= nIQLE; iq++) {

            fFQ = (double) iq;

            // Get looping limits for ir
            // d* limit

            fT1 = fS1P2R2  +  fFQ * (fFQ * fS2  +  fS6P);

            nNRT = nSolveQuadratic(fS3, (0.5f * (fS4 * fFQ  + fS5P)), fT1, fVec);

            if (nNRT == 1) {
              nIRMIN = nIntMinusOne(fVec[0]);
              nIRMAX = nIntPlusOne(fVec[1]);
            }
            else if (nNRT > 1) {
//              nIRMIN = nIntMinusOne(fVec[0] <? fVec[1]);
//              nIRMAX = nIntPlusOne(fVec[0] >? fVec[1]);
              nIRMIN = nIntMinusOne(min(fVec[0] , fVec[1]));
              nIRMAX = nIntPlusOne(max(fVec[0] , fVec[1]));
            }
            if (nNRT >= 1) {

              // Ewald sphere limits

              for (i = 0; i < 2; i++) {
                fPX = fP1H[i]  +  fP[i][1][0] * fFQ;
                fPY = fP2H[i]  +  fP[i][1][1] * fFQ;
                fPZ = fP3H[i]  +  fP[i][1][2] * fFQ;
                fT1 = fPX * fPX  +  fPY * fPY  +  fPZ * fPZ  -  1.0f;
                fT2 = fP[i][2][0] * fPX + fP[i][2][1] * fPY + fP[i][2][2] * fPZ;

                nNRTA[i] = nSolveQuadratic(fQU[i], fT2, fT1, fVec);

                if (nNRTA[i] == 1) {
                  fRLB[i] = fVec[0];
                  fRLE[i] = fVec[0];
                }
                else if (nNRTA[i] > 1) {
//                  fRLB[i] = fVec[0] <? fVec[1];
//                  fRLE[i] = fVec[0] >? fVec[1];
                  fRLB[i] = min(fVec[0] , fVec[1]);
                  fRLE[i] = max(fVec[0] , fVec[1]);
                }
              } // end i loop

              // Use quadratic roots to set up loop limits depending on sphere
              // intersections and adjust loops to other limits.

              // JWP: I gave up trying to decipher the original Fortran GOTOs,
              //       so they appear below.  I apologize!

              nLOOP = 1;
              if (nNRTA[0] > 0) goto L270;
              if (nNRTA[1] <= 0) goto L450;   // Could be replace with continue;
              nNRTA[0] = nNRTA[1];
              fRLB[0] = fRLB[1];
              fRLE[0] = fRLE[1];
              goto L300;

            L270: ;
              if (nNRTA[1] <= 0)
                {
                  goto L300;
                }
              if (fRLB[0] <= fRLB[1]) goto L280;

              vSwap2Values(&nNRTA[0], &nNRTA[1]);
              vSwap2Values(&fRLB[0],  &fRLB[1]);
              vSwap2Values(&fRLE[0],  &fRLE[1]);
            L280: ;
              if (nNRTA[0] == 1) goto L290;
              if (nNRTA[1] == 2) goto L310;
            L290: ;
              nLOOP = 2;
            L300: ;
              // This adjust loops has problems with UMC sometimes
              vAdjustLoopLimits(nLOOP, fRLB[0], fRLE[0], fRLB[1],
                                fRLE[1], &nIRLB[0], &nIRLE[0]);
              goto L330;

            L310: ;
              if (fRLE[1] > fRLE[0]) goto L320;
              nLOOP = 2;
              vAdjustLoopLimits(nLOOP, fRLB[0], fRLB[1], fRLE[1],
                                fRLE[0], &nIRLB[0], &nIRLE[0]);
              goto L330;

            L320: ;
              if (fRLB[1] >= fRLE[0]) goto L290;
              nLOOP = 2;
              vAdjustLoopLimits(nLOOP, fRLB[0], fRLB[1], fRLE[0],
                                fRLE[1], &nIRLB[0], &nIRLE[0]);

              //  Check ir limits

            L330: ;
              nIER = nOutsideLimits(nIRMIN, nIRMAX, &nIRLB[0], &nIRLE[0]);

              if (nIER != 0) goto L350;
              if (nLOOP == 1) goto L370;
            L340: ;
              vSwap2Values(&nIRLB[0], &nIRLB[1]);
              vSwap2Values(&nIRLE[0], &nIRLE[1]);
              nIER = nOutsideLimits(nIRMIN, nIRMAX, &nIRLB[0], &nIRLE[0]);
              if (nIER == 0) goto L370;
              goto L360;

            L350: ;
              if (nLOOP == 1) goto L450;   // Could be replaced with continue;
              nLOOP = 1;
              goto L340;


            L360: ;
              if (nLOOP ==  1) goto L450;   // Could be replaced with continue;
              vSwap2Values(&nIRLB[0], &nIRLB[1]);
              vSwap2Values(&nIRLE[0], &nIRLE[1]);
              nLOOP = 1;

              //  Loop over IR

              L370: ;

              for (i = 0; i < 3; i++) {
                fPX0[i] = fP1H0[i]  +  fP0[0][1][i] * fFQ;
              }

//+jwp 26-Sep-1995
//  Some problems cropped up double-predicting reflections
//  So try to overcome them here, always have nIRLB[0] <= nIRLE[0]
//                                        and nIRLB[1] <= nIRLE[1]

              if (nLOOP == 2)
                {
                  if ( (nIRLE[0] >= nIRLB[1]) && (nIRLE[0] <= nIRLE[1]) )
                    {
                      // Loops have overlapping range, so adjust into one range

                      nLOOP = 1;
                      nIRLE[0] = nIRLE[1];
                    }
                  else if ( (nIRLE[1] >= nIRLB[0]) && (nIRLE[1] <= nIRLE[0]) )
                    {
                      nLOOP = 1;
                      nIRLB[0] = nIRLB[1];
                    }
                }
//-jwp 26-Sep-1995

              for (iloop = 1; iloop <= nLOOP; iloop++) {
                for (ir = nIRLB[iloop-1]; ir <= nIRLE[iloop-1]; ir++) {

                    nHKL[m_anSortedAxes[0]] = ip;
                    nHKL[m_anSortedAxes[1]] = iq;
                    nHKL[m_anSortedAxes[2]] = ir;

                    // fX12 will contain the reciprocal lattice coords of the
                    // reflection at the start of the rotation range.
                    fFR = (double) ir;                    
                    for (i = 0; i < 3; i++) {
                        fX12[i] = fPX0[i]  +  fP0[0][2][i] * fFR;
                    }
                    if (!nPredictInsertRefln(nHKL,fStartRad,fEndRad,fX12,oRefln)) {
                        m_poReflnlist->nInsert(&oRefln);
                    };
                }  // end ir loop
              }  // end iloop loop
            }  // nNRT endif
          L450: ;
          }  // iq loop
        }  // nNRT endif
      } // nNRTA[0] && nNRTA[1] endif
    } // nNRT endif
  } // end ip loop 460
  return (0);

}

int Cpredict::nPredictReflnsFast(const double fRotStart, const double fRotEnd) {
    float fRotS,fRotE,fRotM;
    double fRotTol;
    Crefln oRefln(m_poReflnlist);   // Create a reflection object
    float* pf0;
    double fStartRad,fEndRad;
    double a3fMosCoeffs[3];
    double a3x3fGUBMat[3][3];
    double a3x3fGMat[3][3];
    double a3x3fUBMat[3][3];
    double a3fHKL[3];
    int    nHKL[3];
    double f0;
    float* pfPrePredictRotMid;
    int*   pnPrePredictPackedHKL;
    int nRotAxis;
    int nMidIndex,nIndex;
    int nPass;
    int nStep;
    float a3x3fFloatBuf[3][3];
    double a3fXRC[3];
    double fRotStartRef,fRotEndRef;
    double fAverageShift;
    int    nAverageShift;
  
    fRotS = min(fRotStart, fRotEnd);
    fRotE = max(fRotStart, fRotEnd);
    fRotM = (fRotS + fRotE)*0.5;
    fStartRad = fRotS * fRADIANS_PER_DEGREE;    
    fEndRad   = fRotE * fRADIANS_PER_DEGREE;

    nRotAxis = m_poCrysGonio->nGetNum(m_poRotation->sGetName());
    m_poCrysGonio->vCalcGetRotMatrix(& a3x3fFloatBuf[0][0],nRotAxis); 
    vCopyMat3D(&a3x3fFloatBuf[0][0],&a3x3fGMat[0][0]);
    m_ppoCrystal[0]->nCalcOrientMatrix(m_poSource->fGetWavelength());
    m_ppoCrystal[0]->vGetOrientMatrix(&a3x3fUBMat[0][0]);
    vMulMat3DMat3D(a3x3fGMat,a3x3fUBMat,a3x3fGUBMat);

    
    // Set detector number field of reflection template that will be inserted
    // into the resultant poReflnlist.
    
    oRefln.vSetField(m_poReflnlist->m_nFI_nDetNum, m_nWhichDetector);

    m_ppoCrystal[0]->vGetMosaicity(a3fMosCoeffs);
    fRotTol = m_fPrePredictRotWidth*a3fMosCoeffs[0];
    g_fCmpFloatsTol = fRotTol;
    
    if( 0 == m_afPrePredictRotMid.size() )
        return 1;

    pfPrePredictRotMid = &m_afPrePredictRotMid[0];
    
    pnPrePredictPackedHKL = & m_anPrePredictPackedHKL[0];
    
    pf0 = (float*) bsearch(&fRotM,pfPrePredictRotMid,m_afPrePredictRotMid.size(),sizeof(float),float_cmp_tol);
    
    if (!pf0)
        return 1;
    nMidIndex = pf0 - pfPrePredictRotMid;
    fAverageShift = 0.0;
    nAverageShift = 0;

    for (nPass = 0; nPass < 2; nPass++) {
        nStep = (nPass==0)?-1:1;
        for (nIndex = nMidIndex + nPass;(nIndex >= 0) && (nIndex < m_afPrePredictRotMid.size()); nIndex += nStep) {
            // We might have zeroed out this reflection to speed things up later.
            if (pnPrePredictPackedHKL[nIndex]) {
                // Check to see that the reflection is in the rotation rocking width.
                // But make sure to account for the average shift encountered due to refinements made during integration.
                // The average shift is cumulative and is always computable by using fAverageShift/max(1,nAverageShift)
                f0 = pfPrePredictRotMid[nIndex] + fAverageShift/max(1,nAverageShift);
                if ((f0 >= fRotStart - fRotTol) && (f0 <= fRotEnd + fRotTol)) {
                    a3fHKL[0] = nHKL[0] = ((pnPrePredictPackedHKL[nIndex] & 0x7fe00000) >> 21) -  512;  //(1023L << 21) = 0x7fe00000
                    a3fHKL[1] = nHKL[1] = ((pnPrePredictPackedHKL[nIndex] &   0x1ff800) >> 11) -  512;  //(1023L << 11) =   0x1ff800
                    a3fHKL[2] = nHKL[2] = ((pnPrePredictPackedHKL[nIndex] &      0x7fe) >>  1) -  512;  //(1023L <<  1) =    0x7fe
                                        
                    vMulMat3DVec3D(m_afRotStartMatrix, a3fHKL, a3fXRC);
                    if (!nPredictInsertRefln(nHKL,fStartRad,fEndRad,a3fXRC,oRefln)) {
                        
                        fRotStartRef = oRefln.fGetField(m_poReflnlist->m_nFI_fCalcRotStart);
                        fRotEndRef = oRefln.fGetField(m_poReflnlist->m_nFI_fCalcRotEnd);
                        
                        if (max(fRotEnd,fRotEndRef) - min(fRotStart,fRotStartRef) <
                            fRotEnd -fRotStart + fRotEndRef - fRotStartRef) {   
                            // For each reflection that is added to the list, we keep a record of the average shift encountered.
                            // This shift must be used to offset which reflections are in the list.
                            fAverageShift += (fRotStartRef + fRotEndRef)*0.5 - pfPrePredictRotMid[nIndex] ;
                            nAverageShift ++;
                            m_poReflnlist->nInsert(&oRefln);
                        };
                    } else
                        pnPrePredictPackedHKL[nIndex] = 0;
                } else
                    break;
            };

        };
    };
    
    return 0;
        
}

int Cpredict::nPredictInsertRefln(int nHKL[3], 
                              double fStartRad, 
                              double fEndRad, 
                              double fX12[3], 
                              Crefln& oRefln)
{
    int     nIER = 0;
    bool    bNonunfReject = false;
    double  fSca1 = 0.0;

    float    fObliqueIncidenceCorrection = 0.0f;
    
    // Test for any systematic absences
    if( !m_ppoCrystal[m_nWhichCrystal]->m_poSpacegroup->bIsExtinct(nHKL[0], nHKL[1], nHKL[2]) )
    {
        // Test for valid reflection within d* limit
        // calculate refln starting rot angle value and
        // reflecting range value.
        // Returns reflection coords at diffraction condition
        // on the Ewald sphere in m_afXR and the
        // reflection's rotation start, end, width, and
        // if refln OK the reciprocal Lorentz factor in
        // member variables:

        nIER = nCalcRecipCoords(fX12, fStartRad, fEndRad);                  
        
        if( nIER == 0 )
        {
            
            // Test whether reflection is within detector limits
            // calculates member variables
            
            // The function nCalcDriftDetCoords() must be called BEFORE calling
            // nCalcDetCoords().  nCAlcDriftDetCoords() modifies the m_XXXX variables
            // by making a sub-call to nCalcDetCoords().
            
            nIER = nCalcDetCoords();                    

            if (nIER == 0)
            {
                // Acceptable reflection!
                oRefln.vSetH(nHKL[0]);
                oRefln.vSetK(nHKL[1]);
                oRefln.vSetL(nHKL[2]);
                
                
                // See if the predicted position falls in a
                // bad nonunformity region.
                // If so, either keep it in the list but flag it,
                // or  reject it.
                oRefln.vSetField(m_poReflnlist->m_nFI_nNonunfFlag, (int) 0);
                
                bNonunfReject = m_ppoDetector[m_nWhichDetector]->m_poNonunf->bPixelIsBad((int) m_fPx0, (int) m_fPx1);
                
                if( bNonunfReject )
                {
                    // In bad non-uniformity region decide what to do
                    if (0 == m_nNonunfFlag)
                    {
                        // Keep it, but flagged as bad
                        oRefln.vSetField(m_poReflnlist->m_nFI_nNonunfFlag, (int) 1);
                        
                        bNonunfReject = FALSE;
                    }
                }
                
                if( !bNonunfReject )
                {
                    // If NOT a reject proceed to calculate rest of reflections
                    // parameters and insert into the list
                    oRefln.vSetField(m_poReflnlist->m_nFI_fLorentz, (float)m_fLorentz);
                    
                    if( m_fLorentz > m_fMaxLorentz )
                    {
                        // Keep it, but flagged as bad
                        oRefln.vSetField(m_poReflnlist->m_nFI_nNonunfFlag, (int) 1);
                    }
                    
                    // Calculate polarization correction
                    oRefln.vSetField(m_poReflnlist->m_nFI_fPolarz, (float)fCalcGetPolarz());
                    
                    
                    /////////////////////////////////////////////////////////////////////////////
                    // Oblique incidence correction = 1/cos(angle) between detector_normal and
                    // scattered_beam_wavevector
                    // 1/cos(angle) = | DN | / ( DN.(XR - S0) )
                    //                  --       --  --   --
                    // Where XR-S0 = S, and |S| = 1.0
                    // For a flat detector, the detector_normal is a constant vector.
                    // For cylindrical detector, each point will have its own normal.

                    vSubVec3DVec3D(m_afXR, m_afS0, m_afXRminusS0);
		    
                    // The reason we treat cylindrical and flat detectors separately below
                    // is because we don't want to call fLenVec3D() for the flat detector over and over
                    // again as that value is a constant.
                    if( m_bDetectorIsCylindrical )
                      fObliqueIncidenceCorrection = fLenVec3D(m_a3fLocalDN) / fDot3D(m_a3fLocalDN, m_afXRminusS0);
		            else
                      fObliqueIncidenceCorrection = m_fLenDN / fDot3D(m_afDN, m_afXRminusS0);
                    
                    oRefln.vSetField(m_poReflnlist->m_nFI_fOblique, fObliqueIncidenceCorrection); 
                    ////////////////////////////////////////////////////////////////////////////////////////

                    // Resolution calculation
                    m_fReflnResolution = m_fWavelength / sqrt(max(m_fReflnDstarSq,(double)1.E-8));
                    
                    oRefln.vSetField(m_poReflnlist->m_nFI_fResolution, (float) m_fReflnResolution);
                    
                    // Predicted rot start, end, mid, width
                    m_fReflnRot0  = m_fReflnRot0  / fRADIANS_PER_DEGREE;
                    m_fReflnStart = m_fReflnStart / fRADIANS_PER_DEGREE;
                    m_fReflnEnd   = m_fReflnEnd   / fRADIANS_PER_DEGREE;
                    m_fReflnWidth = m_fReflnWidth / fRADIANS_PER_DEGREE;
                    
                    oRefln.vSetField(m_poReflnlist->m_nFI_fCalcRotMid, (float) m_fReflnRot0);
                    oRefln.vSetField(m_poReflnlist->m_nFI_fCalcRotStart, (float) m_fReflnStart);
                    oRefln.vSetField(m_poReflnlist->m_nFI_fCalcRotEnd, (float) m_fReflnEnd);
                    oRefln.vSetField(m_poReflnlist->m_nFI_fCalcRotWidth, (float) m_fReflnWidth);
                    
                    // Mosaicity parameters.
                    if ((m_bCalcMosaicityLinearCoeffs) || (m_bCalcMinMosaicityForRefln)) 
                        oRefln.vSetField(m_poReflnlist->m_nFI_fCalcMosCoeffA,(float) m_fMosCoeffA);
                    if (m_bCalcMosaicityLinearCoeffs)
                        oRefln.vSetField(m_poReflnlist->m_nFI_fCalcMosCoeffB,(float) m_fMosCoeffB);
                    
                    // Detector coordinates in pixels
                    oRefln.vSetField(m_poReflnlist->m_nFI_fCalcPx0, (float) m_fPx0);
                    oRefln.vSetField(m_poReflnlist->m_nFI_fCalcPx1, (float) m_fPx1);
                    oRefln.vSetField(m_poReflnlist->m_nFI_fCalcXmm, (float) m_fXmm);
                    oRefln.vSetField(m_poReflnlist->m_nFI_fCalcYmm, (float) m_fYmm);
                    
                    if( m_poReflnlist->m_nFI_fCalcZmm >= 0 ) 
                        oRefln.vSetField(m_poReflnlist->m_nFI_fCalcZmm, (float) m_fZRelativeCoordinate);
                    
                    //  Keep reciprocal lattice coordinates when reflection
                    //   is on sphere (at Laue condition).
                    //  NOTE: These are dimensionless!
                    
                    oRefln.vSetField(m_poReflnlist->m_nFI_fRecipCoord0, (float) m_afXR[0]);
                    oRefln.vSetField(m_poReflnlist->m_nFI_fRecipCoord1, (float) m_afXR[1]);
                    oRefln.vSetField(m_poReflnlist->m_nFI_fRecipCoord2, (float) m_afXR[2]);
                    
                    if (m_poReflnlist->m_nFI_nTwinID != -1)
                        oRefln.vSetField(m_poReflnlist->m_nFI_nTwinID,m_nWhichCrystal+1 + m_nWhichCrystalTwinLaw*1000);
                    
                    
                    //  Set rejection flag if reflection is too wide in Phi
                    //  (But keep reflection in list!):
                    //
                    // IREFL(QSEG,NR) = 0
                    // IF (PHIRNG .GT. fReflnWidthMax) IREFL(QSEG, NR) = QWMAX
                    
                    // Finally insert this reflection into the reflection list
                    return 0;
                    
                } // bNonunfReject endif
            }  // nIER nCalcDetCoords endif
       }  // nIER nCalcRecipCoords endif
                  
    }  // !bIsExtinct endif
                  
    return 1;
}


int Cpredict::nSortAxes(void)
{
  double fVec1[3], fVec2[3];
  int   i;                // Loop counter
  double fDC[3];           // Direction cosines
  double fTemp1, fTemp2;   // Temp variables
  int   nI1, nI2, nI3;    // More temp variables

  // Construct 3 orthonormal reference vectors: S0 (the source), fVec1, fVec2

  vCross3D(m_afS0, m_afRotVec, fVec1); // fVec1 is perpend. to S0 and RotVec
  vCross3D(fVec1, m_afS0, fVec2);   // fVec2 is perpend. to S0 and fTemp1

  for (i = 0; i < 3; i++) {
    fTemp1 = fDot3D(&m_afRotStartMatrix[i][0], m_afS0);
    fTemp2 = fDot3D(&m_afRotStartMatrix[i][0], &m_afRotStartMatrix[i][0]);
    fDC[i] = fabs(fTemp1) / sqrt(fTemp2);    // Normalize here
  }

  // Get index of longest direction cosine

  nI3 = nGetIndexLargest3D(fDC);

  for (i = 0; i < 3; i++) {
    fTemp1 = fDot3D(&m_afRotStartMatrix[i][0], fVec2);
    fTemp2 = fDot3D(&m_afRotStartMatrix[i][0], &m_afRotStartMatrix[i][0]);
    fDC[i] = fabs(fTemp1) / sqrt(fTemp2);    // Normalize here
  }

  nI1 = nGetIndexLargest3D(fDC);

  if (nI1 == nI3) nI1 = (nI1+1) % 3;
  nI2 = 3 - nI1 - nI3;                // Set nI2 to the one leftover

  m_anSortedAxes[0] = nI3;
  m_anSortedAxes[1] = nI2;
  m_anSortedAxes[2] = nI1;

  return (0);
}

int Cpredict::nGetIndexLargest3D(const double fDC[3])
{
  int i, j;
  register double fTemp = fDC[0];
  j = 0;
  for (i = 0; i < 3; i++) {
    if (fTemp < fDC[i]) {
      j = i;
      fTemp = fDC[i];
    }
  }
  return (j);
}


int Cpredict::nSetupRotation(const double fRotStart, const double fRotEnd)
{
  // Setup rotation matrices for a particular rotation.  This is done
  // after the crystal and goniometer matrices have been setup.
  // This might be a function of the rotation class.

  double fRotS, fRotE;
  float a3fTemp[3];

  // Allow inputs in any order

  fRotS = min(fRotStart, fRotEnd);
  fRotE = max(fRotStart, fRotEnd);

  // Get rotation axis direction
  int nAxis;
  static int s_nNumWarn = 0;
  nAxis = m_poCrysGonio->nGetNum(m_poRotation->sGetName());
  if (0 > nAxis)
    {
      nAxis = 0;
      s_nNumWarn++;
      if (2 > s_nNumWarn)
        cout << "WARNING Cpredict::nSetupRotation, axis name '"
             << m_poRotation->sGetName() << "' not a valid axis\n"
             << " for this crystal goniometer, using '" 
             << m_poCrysGonio->sGetName(nAxis) << "' instead.\n";
    }
  if (m_poCrysGonio->nGetRotVector(nAxis,&a3fTemp[0])) 
    {
      cout << "ERROR Cpredict::nSetupRotation calculating rotation axis vector.\n" << flush;
      return 1;
    }
  //+jwp 9-Dec-2003
  if (NULL != m_poRotation)
    {
      // Copy the actual rotation vector into the Crotation object
      m_poRotation->vSetVector(&a3fTemp[0]);
      //      cout << "Predict rotation axis ";
      //      vListVec3D(a3fTemp);
    }
  //-jwp 9-Dec-2003
  vCopyVec3D(a3fTemp,m_afRotVec);

  // Setup rotation matrices at start and end of input range
  // Multiply these by the goniometer and crystal matrix

  double fTemp[3][3];
  vConvRotVec3DMat3D(fRotS, m_afRotVec, &fTemp[0][0]);
  vMulMat3DMat3D(fTemp, m_afGonCrysMatrix, m_afRotStartMatrix);

  vConvRotVec3DMat3D(fRotE,  m_afRotVec, &fTemp[0][0]);
  vMulMat3DMat3D(fTemp, m_afGonCrysMatrix, m_afRotEndMatrix);

  vConvRotVec3DMat3D(fRotEnd+(m_fReflnWidthMax*0.5f),  m_afRotVec, &fTemp[0][0]);
  vMulMat3DMat3D(fTemp, m_afGonCrysMatrix, m_afRotPastEndMatrix);

  return (0);
}

int Cpredict::nSetupDetector(const int nWhich)
{
  // Setup Detector things for a prediction run

  int nStat;
  float fResoMaxEdge;

  if (nWhich > m_nNumDetectors) return (1);
  m_nWhichDetector = nWhich;

  // Get detector DD matrix and DNvec

  (void) m_ppoDetector[m_nWhichDetector]->nCalcGetDDDN(&m_afDD[0][0], m_afDN);
  m_fDet  = fInvMat3D(&m_afDD[0][0], &m_afInvDD[0][0]);    // Invert DD matrix
  if (m_fDet == 0.0) return (2);

  // Get Min and Max resolutions on detector.  Since m_afS0 is
  //  already normalized

  nStat = m_ppoDetector[m_nWhichDetector]->nGetResolution(pfCast3D(m_afS0),
                                      &m_fDetResolMin, &m_fDetResolMax,
                                                          &fResoMaxEdge);

  // If the user set resolution range is less than 0, use the
  // detector max resolution to the edge (not to corner)
  
  //  cout << "INFO: Det reso: min, max, edge: " 
  //       << m_fDetResolMin << ", " << m_fDetResolMax << ", " 
  //       << fResoMaxEdge << endl;
  if (0.0 >= m_fResolutionMin)
    {
      m_fResolutionMin = m_fDetResolMin * m_fWavelength;
    }
  
  if (0.0 >= m_fResolutionMax)
    {
      if ( (m_bUseEdgeReso) || (0.0 > m_fResolutionMax) )
        {
          m_fResolutionMax = fResoMaxEdge * m_fWavelength;
        }
      else
        {
          m_fResolutionMax = m_fDetResolMax * m_fWavelength;
        }
    }
  else if (m_fResolutionMax < (m_fDetResolMax * m_fWavelength))
    {
      m_fResolutionMax = m_fDetResolMax * m_fWavelength;
      //      cout << "INFO: Max resolution restricted to max detector resolution of "
      //           << m_fResolutionMax << endl;
    }
  //  cout << "INFO: Predict min reso is: " << m_fResolutionMin << endl;
  //  cout << "INFO: Predict max reso is: " << m_fResolutionMax << endl;

  if( 0 != nStat ) 
      return 2;  // Problem with detector object

  // fDetResolMin,Max already calculated with normalized m_afS0
  // Compute (d*)^2 min and max (assuming wavelength = 1.0)

  m_fDstarSqMin = max( (1.0f/(m_fDetResolMin * m_fDetResolMin)),
                       (m_fWavelength * m_fWavelength)
                       / (m_fResolutionMin * m_fResolutionMin) );
  m_fDstarSqMax = min( (1.0f/(m_fDetResolMax * m_fDetResolMax)),
                       (m_fWavelength * m_fWavelength)
                       / (m_fResolutionMax * m_fResolutionMax) );
  m_fMax2Theta = acos(min(1.0,max(-1.0,1.0 - 0.5/(m_fDetResolMax*m_fDetResolMax))));

  // Compute some things that don't change for each reflection:

  m_fS0dotDN  = fDot3D(m_afS0, m_afDN);                    // MUST have m_afS0 first!
  m_fDNdotDN  = fDot3D(m_afDN, m_afDN);
  m_fLenDN    = fLenVec3D(m_afDN);

  m_fDetDistance = m_ppoDetector[m_nWhichDetector]->fGetDistance();
  m_bDetectorIsCylindrical
    = (eSpatial_cylindrical_geometry == m_ppoDetector[0]->m_poSpatial->eSpatialType());

  return (0);
}

int Cpredict::nSetupSourceCrystal(const int nWhichTwin,const int nWhichTwinLaw)
{
  // Setup Source and Crystal things for a prediction run

  int nStat;
  float a4fTemp[4];
  double a3fMosCoeffs[3];

  // Source stuff:

  m_poSource->vCalcGetS0(a4fTemp);               // Get normalized source dir. vector
  vCopyVec3D(a4fTemp,m_afS0);
  m_fWavelength = m_poSource->fGetWavelength(); // Get source wavelength;
  m_poSource->vGetPolarz(a4fTemp);              // Get normalized source polarization
  vCopyVecND(4, a4fTemp, m_afPol);
  //  vCopyVec3D(a4fTemp,m_afPol);
  vCross3D(&m_afS0[0], &m_afPol[1], m_afSN);        // SN = S0 X POL (perpend to both)
  m_fDet = fNormVec3D(m_afSN);                  //    and normalize it
  if (m_fDet <= 0.0) return (3);

  // Get source spectral dispersion

  m_fSourceDispersion = m_poSource->fGetSpectralDispersion();

  // Crystal stuff:

  m_nWhichCrystal = nWhichTwin;
  m_nWhichCrystalTwinLaw = nWhichTwinLaw;

  if ( (1.0 != m_fScaleLength)  && (0.0 < m_fScaleLength) )
    {
      // If m_fScaleLength is not 1.0, scale crystal cell lengths by
      // m_fScaleLength, then reset m_fScaleLength to 1.0 to prevent any
      // future redundant rescaling.  Scaling the cell lengths is used in
      // some strategy calculations to increase the speed of prediction

      int i;
      double a6fCell[6];
      m_ppoCrystal[m_nWhichCrystal]->vGetCell(a6fCell);
      for (i = 0; i < 3; i++)
        {
          a6fCell[i] = m_fScaleLength * a6fCell[i];
        }
      m_fScaleLength = 1.0;
      m_ppoCrystal[m_nWhichCrystal]->vSetCell(a6fCell);
    }
  // Calculate crystal orientation matrix

  nStat = m_ppoCrystal[m_nWhichCrystal]->nCalcOrientMatrix(m_fWavelength);
  if (0 != nStat) return (4);

  // Get crystal orientation matrix with crystal missetting or
  //   orientation angles applied

  m_ppoCrystal[m_nWhichCrystal]->vGetOrientMatrix(&m_afCrysMatrix[0][0]);

  // Multiply by the twin law if applicable.
  if (m_nWhichCrystalTwinLaw) {
          double a3x3fTwinLaw[3][3];
          double a3x3fTemp[3][3];
          m_ppoCrystal[m_nWhichCrystal]->nGetTwinLaw(m_nWhichCrystalTwinLaw,&a3x3fTwinLaw[0][0],true);
          vMulMat3DMat3D(a3x3fTwinLaw,m_afCrysMatrix,a3x3fTemp);
          vCopyMat3D(&a3x3fTemp[0][0],&m_afCrysMatrix[0][0]);
  };
  m_ppoCrystal[m_nWhichCrystal]->nGetRecipShift(m_nWhichCrystalTwinLaw,&m_afRecipShift[0]);

  if (!m_poCrysGonio->bIsAvailable())
    {
      // This means the crystal goniometer axes are not defined

      cout << "ERROR Cpredict::nSetupSourceCrystal: bogus crystal goniometer!\n"
           << flush;
      return (5);
      // Possible workaround is to make m_afGonMatrix == Identity matrix

    }
  else
    {
      // Get crystal goniometer rotation matrix at datum

      double a3x3fTemp[3][3];      
      int   nAxis;

      vZeroMat(3,3,&a3x3fTemp[0][0]);
      nAxis = m_poCrysGonio->nGetNum(m_poRotation->sGetName());
      m_poCrysGonio->vCalcGetRotMatrix(&a3x3fTemp[0][0], nAxis);
      vCopyMat3D(&a3x3fTemp[0][0], &m_afGonMatrix[0][0]);
/*
      cout << "DEBUG in nSetupSourceCrystal\n";
      vListMat3D((double*)&m_afGonMatrix[0][0]);
      vListMat3D((double*)&a3x3fTemp[0][0]);
*/
    }

  // Multiply goniometer rotation matrix by crystal orientation matrix

  vMulMat3DMat3D(m_afGonMatrix, m_afCrysMatrix, m_afGonCrysMatrix);

  // Get crystal mosaicity in RADIANS:

  m_ppoCrystal[m_nWhichCrystal]->vGetMosaicity(a3fMosCoeffs);
    
  
  if( m_wVerbose & DTREK_PREDICT_VERBOSE_MOSAICITY )
  {
        // RB: Need to restore the previous precision.  
        int       nSavedIOSPrecision = cout.precision();
  
        cout << "\nComponent " << (m_nWhichCrystal+1) << " crystal mosaicity: " 
            << setprecision(3) << (a3fMosCoeffs[0] + a3fMosCoeffs[2]) << endl << flush;

        cout << setprecision(nSavedIOSPrecision);
  }
  
  return 0;
}

int Cpredict::nExpandReflnlist(Creflnlist *poReflnlistIn)
{

// Add required fields to the reflection list if they are not already there.
// The fields are named by the static member variables of the Creflnlist class
//

  Creflnlist *poReflnlistL;   // Local reflection list

  if (NULL == poReflnlistIn)
    {
      poReflnlistL = m_poReflnlist;     // Use on member variable
    }
   else
    {
      poReflnlistL = poReflnlistIn;   // Use passed reflection list
    }

  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_snDetNum);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_snNonunfFlag);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcPx0);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcPx1);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcRotMid);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcRotStart);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcRotEnd);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcRotWidth);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcXmm);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcYmm);
  if ((m_ppoDetector) && (m_ppoDetector[0]->m_poSpatial->eSpatialType() == eSpatial_cylindrical_geometry)) {
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcZmm);
  };

  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfPolarz);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfLorentz);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfOblique);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfPartiality);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfResolution);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfRecipCoord0);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfRecipCoord1);
  (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfRecipCoord2);

  if ((m_nNumCrystals > 1) || (m_nMaxTwinLaws > 0))
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_snTwinID);

  if ((m_bCalcMosaicityLinearCoeffs) || (m_bCalcMinMosaicityForRefln))
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcMosCoeffA);
  if (m_bCalcMosaicityLinearCoeffs) 
      (void) poReflnlistL->nExpandGetField(poReflnlistL->ms_sfCalcMosCoeffB);

  return (0);

}

double 
Cpredict::dGetResolution(const int nMinMax)
{
  if (0 == nMinMax)
    return (m_fResolutionMin);
  else
    return (m_fResolutionMax);
}
void 
Cpredict::vSetResolution(const double fResoMin, const double fResoMax)
{
  // Set the min and max resolution in Angstroms.  Min is closest to
  // beam, max is as far away from beam, thus Min has a higher value
  // in Angstroms than Max

  // Note: the actual resolution used will likely be altered
  // by the detector position

  m_fResolutionMin = max(fResoMin, fResoMax);
  m_fResolutionMax = min(fResoMin, fResoMax);
}

Cnonunf* Cpredict::poGetNonunf(void)
{
  Cnonunf *poNonunf = NULL;
  if (NULL != m_ppoDetector)
    {
      if (NULL != m_ppoDetector[0])
        {
          if (m_ppoDetector[0]->bIsAvailable())
            {
              poNonunf = m_ppoDetector[0]->m_poNonunf;
              if (poNonunf->bIsAvailable())
                return (poNonunf);
              poNonunf = NULL;
            }
        }
    }
  return (poNonunf);
}

int Cpredict::nCalcNearestCentroids(Creflnlist& oList,
                                    int nReflnStart,
                                    int nReflnEnd,
                                    double fObsRotStart,
                                    double fObsRotEnd,
                                    int nMinSearchRadius,
                                    float* pfCentroid0,
                                    float* pfCentroid1)
{
    int nDivision = 20;
    int nx;
    int nDivisionSq;
    itr<int> anOccupied;
    itr<int> anLink;
    itr<int> anDiv0;
    itr<int> anDiv1;
    itr<int> anCompare;
    int nRef,nRef2;
    double fPix0,fPix1;
        double fRotStart,fRotEnd;
    int nDiv0,nDiv1;
    int nDim0,nDim1;
    int a2nDelta[2];

    

    if (nReflnEnd < 0)
        nReflnEnd = oList.nGetNumReflns() - 1;

    nDim0 = m_ppoDetector[0]->m_a2nDim[0];
    nDim1 = m_ppoDetector[0]->m_a2nDim[1];

    nDivision = min(nDim1/nMinSearchRadius,nDim0/nMinSearchRadius);
    nDivisionSq = nDivision*nDivision;

    for (nx = 0; nx< nDivisionSq; nx++)
        anOccupied[nx] = -1;
    for (nRef = nReflnStart; nRef <= nReflnEnd; nRef++) {
        fPix0 = oList[nRef].fGetField(oList.m_nFI_fCalcPx0);
        fPix1 = oList[nRef].fGetField(oList.m_nFI_fCalcPx1);
        nDiv0 = max(0,min(nDivision-1,(int) (nDivision*(fPix0)/nDim0)));
        nDiv1 = max(0,min(nDivision-1,(int) (nDivision*(fPix1)/nDim1)));
        anDiv0[nRef - nReflnStart] = nDiv0;
        anDiv1[nRef - nReflnStart] = nDiv1;
        if (anOccupied[nDiv0 + nDivision*nDiv1]==-1) {
            anOccupied[nDiv0 + nDivision*nDiv1] = nRef;
            anLink[nRef-nReflnStart] = -1;
        } else {
            anLink[nRef-nReflnStart] = anOccupied[nDiv0 + nDivision*nDiv1];
            anOccupied[nDiv0 + nDivision*nDiv1] = nRef;
        };
    };

    memset(pfCentroid0,0,sizeof(float)*oList.nGetNumReflns());
    memset(pfCentroid1,0,sizeof(float)*oList.nGetNumReflns());

    // Next, do comparisions for all reflections that are near the given reflection.
    for (nRef = nReflnStart; nRef <= nReflnEnd; nRef++) {
        nDiv0 = anDiv0[nRef - nReflnStart];
        nDiv1 = anDiv1[nRef - nReflnStart];
        -anCompare;
        for (a2nDelta[0] = -1; a2nDelta[0] <= 1; a2nDelta[0]++) {
            for (a2nDelta[1] = -1; a2nDelta[1] <= 1; a2nDelta[1]++) {
                if ((a2nDelta[0] + nDiv0 >=0) && (a2nDelta[0] + nDiv0 < nDivision) &&
                    (a2nDelta[1] + nDiv1 >=0) && (a2nDelta[1] + nDiv1 < nDivision)) {
                    for (nRef2 = anOccupied[(a2nDelta[0]+nDiv0) + (a2nDelta[1]+nDiv1)*nDivision]; nRef2 != -1; nRef2 = anLink[nRef2 - nReflnStart]) {
                        if (nRef != nRef2) 
                            anCompare + nRef2;
                    };
                };
            };
        };
        fPix0 = oList[nRef].fGetField(oList.m_nFI_fCalcPx0);
        fPix1 = oList[nRef].fGetField(oList.m_nFI_fCalcPx1);
                fRotStart = oList[nRef].fGetField(oList.m_nFI_fCalcRotStart);
                fRotEnd = oList[nRef].fGetField(oList.m_nFI_fCalcRotEnd);


        // We have a list of reflections that are 'near' nRef.
        double fMinDist,fDist;
        double fPix0Comp,fPix1Comp;
                double fRotStartComp,fRotEndComp;
                int        nMinRef;
                
        int nCompare;

                fMinDist = 1e9;
                nMinRef = -1;
        for (nCompare = 0; nCompare < anCompare.size(); nCompare++) {
            nRef2 = anCompare[nCompare];


                        fRotStartComp = oList[nRef2].fGetField(oList.m_nFI_fCalcRotStart);
                        fRotEndComp = oList[nRef2].fGetField(oList.m_nFI_fCalcRotEnd);
                        // See if this reflection has an overlap with the range specified.
                        if (max(fObsRotEnd,fRotEndComp) - min(fObsRotStart,fRotStartComp) <
                                fObsRotEnd -fObsRotStart + fRotEndComp - fRotStartComp) {                               
                                fPix0Comp = oList[nRef2].fGetField(oList.m_nFI_fCalcPx0);
                                fPix1Comp = oList[nRef2].fGetField(oList.m_nFI_fCalcPx1);
                                fDist = (fPix0Comp - fPix0)*(fPix0Comp - fPix0) +
                                        (fPix1Comp - fPix1)*(fPix1Comp - fPix1);
                                if (fDist < fMinDist) {
                                        fMinDist = fDist;
                                        nMinRef = nRef2;
                                };
                        };
        };
                if (nMinRef != -1) {
                        fPix0Comp = oList[nMinRef].fGetField(oList.m_nFI_fCalcPx0);
                        fPix1Comp = oList[nMinRef].fGetField(oList.m_nFI_fCalcPx1);
                        pfCentroid0[nRef] = fPix0Comp;
                        pfCentroid1[nRef] = fPix1Comp;                  
                } else {
                        pfCentroid0[nRef] = 0.0;
                        pfCentroid1[nRef] = 0.0;
                };
    };

        return 0;
          
  
};
