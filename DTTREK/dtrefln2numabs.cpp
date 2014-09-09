//
// Copyright (c) 2006 Rigaku
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtrefln2numabs.cc     Initial author: J.W. Pflugrath           21-Feb-2000
//    This reads a d*TREK reflnlist file and converts it to a file suitable
//    for input to T. Higashi's numabs program.
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
//    dtrefln2numabs reads from the command line the reflnlist filename
//    if present (default: *dtprofit.ref) and the name of the numabs output
//    file to write if present (default: numabs.hkl).
//
//    Syntax:   % dtrefln2numabs [input_reflnlist_file [output_numabs_file]]
//    Example:  % dtrefln2numabs numabs.hkl
//
//+ToDo
//
//   Better error messages need implementing.
//
//+Include files

#include "Dtrek.h"
#include "dtrekdefs.h"
#include "Cstring.h"
#include "Crefln.h"
#include "Creflnlist.h"
#include "dskio.h"
#include "Cimage_header.h"
#include "dtrekvec.h"

#define DT_MAX_ORIENT  100
#define NUMABS_HKL_LIMIT 511

typedef struct _tagNumabsRec
{
  int     nReflnNum;           // Refln counter for given orientation and batch
  int     nH;                  // Miller index h
  int     nK;                  // Miller index k
  int     nL;                  // Miller index l
  float   fIntensity;          // Observed estimated intensity
  float   fSigmaI;             // Observed estimated standard dev. of intensity
  float   a3fSource[3];        // Source wave vector
  float   a3fScatt[3];         // Scattered beam wave vector
  float   fTransFactor;        // Place holder for transmission factor
  int  nWrite(FILE* f);
} tagNumabsRec;

int tagNumabsRec::nWrite(FILE* f) {
  int nStat;
	(void) fprintf(f,"%d %d %d %f %f %f %f %f %f %f %f %f\n",
	nH, nK, nL, fIntensity, fSigmaI, 
        a3fSource[0], a3fSource[1], a3fSource[2],
	a3fScatt[0],  a3fScatt[1],  a3fScatt[2], fTransFactor);
	nStat = ferror(f);  // bug, no clearerr() call
	return (nStat);
};

#ifdef SSI_PC
#include "CrclHelper.h"
#include "CrclIncludes.h"
#define main dtrefln2numabs_main
#endif

int main(int   argc,    // Command line argument count
         char *argv[])  // Pointers to command line args

{
  int      nStat;
  Cstring  sIn            = "dtprofit.ref";
  Cstring  sOut           = "numabs.hkl";
  Creflnlist *poReflnlist = NULL;
  Crefln     *poRefln;
  FILE       *pFTextOutFile = NULL;
  tagNumabsRec tNumabsRec;

#if !defined(VC6)
    using std::cout;
    using std::flush;
    using std::endl;
#endif
  cout << "\ndtrefln2numabs: Copyright (c) 2006 Rigaku\nModified for TEXT output\n";
  cout << D_K_DTREKVersion << endl << flush;

#ifdef SSI_PC
  Cstring sCCAppVersion = (const char*) argv[argc-1];
  if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1) {
	cout << sCCAppVersion;
	argc--; 
	cout << CCrclHelper::GetInstance()->GetTimestamp() << endl;
	argv[argc] = NULL;
  }
#endif

  // Copy command line to output log
  cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << endl << flush;

  // Parse command line arguments, if any

  argc--; argv++; // Skip program name

  sIn = sDtrekGetPrefix() + sIn;
  if (argc >= 1)
    {
      sIn = (const char*) argv[0];
      argc--;
      argv++;
      if (argc >= 1)
	{
	  sOut = (const char*) argv[0];
	argc--;
	argv++;
	}
    }

  pFTextOutFile = fopen(sOut.string(), "w+t");

  cout << "dtrefln2numabs: convert "
       << sIn  << " to "
       << sOut << '\n';

  // Read the reflnlist

  cout << "Reading " << sIn << " ..." << endl << flush;
  poReflnlist = new Creflnlist (sIn);

  if (!poReflnlist->bIsAvailable())
    {
      cout << "dtrefln2numabs: reflnlist " << sIn << " unavailable!\n" << flush;
      delete poReflnlist;
	  #ifdef SSI_PC
		  return 1;
	  #else
	  	  exit (1);
	  #endif
    }

  if (1 > poReflnlist->nGetNumReflns())
    {
      cout << "dtrefln2numabs: reflnlist " << sIn << " has no reflections!\n" << flush;
      delete poReflnlist;
	  #ifdef SSI_PC
		  return 1;
	  #else
	  	  exit (1);
	  #endif
    }

  int nFI_fSvec[3];
  int nFI_fS0vec[3];
  int nFI_nErrorFlag;
  nFI_fSvec[0]  = poReflnlist->nGetFieldIndex(D_K_IntegratefSvec0);
  nFI_fSvec[1]  = poReflnlist->nGetFieldIndex(D_K_IntegratefSvec1);
  nFI_fSvec[2]  = poReflnlist->nGetFieldIndex(D_K_IntegratefSvec2);

  nFI_fS0vec[0] = poReflnlist->nGetFieldIndex(D_K_IntegratefS0vec0);
  nFI_fS0vec[1] = poReflnlist->nGetFieldIndex(D_K_IntegratefS0vec1);
  nFI_fS0vec[2] = poReflnlist->nGetFieldIndex(D_K_IntegratefS0vec2);
  nFI_nErrorFlag= poReflnlist->nGetFieldIndex(D_K_IntegratenErrorFlag);
  if (0 > nFI_nErrorFlag)
    {
      nFI_nErrorFlag= poReflnlist->nGetFieldIndex(D_K_ProfitnBadFlag);
    }

  // Check that the reflnlist has the fields we need ...

  nStat = 0;
  if (0 > poReflnlist->m_nFI_nH)
    nStat++;
  if (0 > poReflnlist->m_nFI_nK)
    nStat++;
  if (0 > poReflnlist->m_nFI_nL)
    nStat++;
  if (0 > poReflnlist->m_nFI_fIntensity)
    nStat++;
  if (0 > poReflnlist->m_nFI_fSigmaI)
    nStat++;
  //  if (0 > poReflnlist->m_nFI_sBatch)
  //    nStat++;
  if (0 > nFI_nErrorFlag)
    nStat++;
  if (0 > nFI_fSvec[0])
    nStat++;
  if (0 > nFI_fSvec[1])
    nStat++;
  if (0 > nFI_fSvec[2])
    nStat++;
  if (0 > nFI_fS0vec[0])
    nStat++;
  if (0 > nFI_fS0vec[1])
    nStat++;
  if (0 > nFI_fS0vec[2])
    nStat++;

  if (0 != nStat)
    {
	  cout << "dtrefln2numabs: reflnlist does not have all the required fields!\nRequired fields are: nH, nK, nL, fIntensity, fSigmaI, nErrorFlag, fSvec vector, and fS0vec vector.\nPlease try using valid reflection file as input.\n" << flush;
      delete poReflnlist;
	  #ifdef SSI_PC
		  return 2;
	  #else
	  	  exit (2);
	  #endif
    }

  int i;
  int nWritten  = 0;
  int nRejected = 0;
  int nExcluded = 0;

  // Open output file

  int nLength   = sOut.length();
  int nFile     = 1;

  if (0 != nStat)
    {
      cout << "ERROR opening output file " << sOut << "!\n" << flush;
      return (nStat);
    }

  nStat                 = 0;
  int     nReflnCount   = 0;
  int     nErrorHKL     = 0;
  Cstring sBatch;
  
  tNumabsRec.fTransFactor = 1.0;

    // Get crystal orientation rotation matrix.
    float angles[3], rotMatrix[3][3], rotMatrix2[3][3], s0[3], s[3], s0p[3], sp[3];
    float x[3], xp[3];
    Cimage_header header( (Cstring)"input.head" );
    //header.nGetValue( D_K_CrystalOrientAngles, 3, angles );
    Ccrystal crystal( header );
    //crystal.vSetOrientAngles( angles );
    crystal.nCalcRotMatrix();
    crystal.vGetRotMatrix( &(rotMatrix[0][0]) );
    crystal.vGetRotMatrix( &(rotMatrix2[0][0]) );
    vTranMat3D( rotMatrix2 );

    //printf("\nOriginal Orientation Angles: %.4f %.4f %.4f\n", angles[0], angles[1], angles[2]);
    printf("U:  %3.4f %3.4f %3.4f\n    %3.4f %3.4f %3.4f\n    %3.4f %3.4f %3.4f\n\n", rotMatrix[0][0],
        rotMatrix[1][0], rotMatrix[2][0], rotMatrix[0][1], rotMatrix[1][1], rotMatrix[2][1], rotMatrix[0][2], 
        rotMatrix[1][2], rotMatrix[2][2]);

  for (i = 0; (0 == nStat) && (i < poReflnlist->nGetNumReflns()); i++)
    {
        poRefln               = poReflnlist->poGetRefln(i);

        // Convert s0 and s to Higashi's definition.
        for( int m=0; m<3; m++ ) {
            s0[m] = poRefln->fGetField(nFI_fS0vec[m]);
            s[m] = poRefln->fGetField(nFI_fSvec[m]);
        }
        vMulMat3DVec3D( rotMatrix, s0, s0p );
        vMulMat3DVec3D( rotMatrix, s, sp );

        //vAddVec3DVec3D( s0, s, x );
        //vAddVec3DVec3D( s0p, sp, xp );


      tNumabsRec.nH          = poRefln->nGetH();
      tNumabsRec.nK          = poRefln->nGetK();
      tNumabsRec.nL          = poRefln->nGetL();
      tNumabsRec.fIntensity  = poRefln->fGetIntensity();
      tNumabsRec.fSigmaI     = poRefln->fGetSigmaI();
      tNumabsRec.a3fSource[0]= s0p[0];//poRefln->fGetField(nFI_fS0vec[2]);
      tNumabsRec.a3fSource[1]= s0p[1];//poRefln->fGetField(nFI_fS0vec[1]);
      tNumabsRec.a3fSource[2]= s0p[2];//poRefln->fGetField(nFI_fS0vec[0]);
      tNumabsRec.a3fScatt[0] = sp[0];//-poRefln->fGetField(nFI_fSvec[2]);
      tNumabsRec.a3fScatt[1] = sp[1];//-poRefln->fGetField(nFI_fSvec[1]);
      tNumabsRec.a3fScatt[2] = sp[2];//-poRefln->fGetField(nFI_fSvec[0]);

        if( i<10 ) {
            printf("S0: %.4f %.4f %.4f  --->  %.4f %.4f %.4f\n", s0[0], s0[1], s0[2], s0p[0], s0p[1], s0p[2]);
            printf("S:  %.4f %.4f %.4f  --->  %.4f %.4f %.4f\n", s[0], s[1], s[2], sp[0], sp[1], sp[2]);
            /*printf("S0: %.4f %.4f %.4f  --->  %.4f %.4f %.4f\n", tNumabsRec.a3fSource[0], tNumabsRec.a3fSource[1], 
                tNumabsRec.a3fSource[2], tNumabsRec.a3fSource[0], tNumabsRec.a3fSource[1], tNumabsRec.a3fSource[2]);
            printf("S:  %.4f %.4f %.4f  --->  %.4f %.4f %.4f\n", tNumabsRec.a3fScatt[0], tNumabsRec.a3fScatt[1], 
            tNumabsRec.a3fScatt[2], tNumabsRec.a3fScatt[0], tNumabsRec.a3fScatt[1], tNumabsRec.a3fScatt[2]);*/
            //printf("|x| = |xp|    --->  %.4f = %.4f\n\n", fLenVec3D(x), fLenVec3D(xp));
        }

      if (   (-999. < tNumabsRec.a3fSource[0])
	  && (-999. < tNumabsRec.a3fSource[1])
	  && (-999. < tNumabsRec.a3fSource[2])
	  && ( 999. > tNumabsRec.a3fScatt[0])
	  && ( 999. > tNumabsRec.a3fScatt[1])
	  && ( 999. > tNumabsRec.a3fScatt[2])
	  && ( -NUMABS_HKL_LIMIT < tNumabsRec.nH)
	  && (  NUMABS_HKL_LIMIT > tNumabsRec.nH)
	  && ( -NUMABS_HKL_LIMIT < tNumabsRec.nK)
	  && (  NUMABS_HKL_LIMIT > tNumabsRec.nK)
	  && ( -NUMABS_HKL_LIMIT < tNumabsRec.nL)
	  && (  NUMABS_HKL_LIMIT > tNumabsRec.nL)
	  )
	{
	  // Legitimate reflection at first glance

	  tNumabsRec.nReflnNum    = nReflnCount++;

	  // Should this refln be excluded (by making reflnnum <0)?

/*
	  if (   (0.0 >= tNumabsRec.fIntensity)
	      || (0.0 >= tNumabsRec.fSigmaI)
	      || (0   != poRefln->nGetField(nFI_nErrorFlag))
	      )
	    {
	      tNumabsRec.nReflnNum = -tNumabsRec.nReflnNum;
	      nExcluded++;
	    }
    else */
	  if (pFTextOutFile)
	    {
	      nStat = tNumabsRec.nWrite(pFTextOutFile);
	      if (0 != nStat)
		{
		  cout << "ERROR writing to output file!: " << nStat << flush;
		}
	      else
		{
		  nWritten++;
		}
	    }
	}
      else
	{
	  if (   ( -NUMABS_HKL_LIMIT > tNumabsRec.nH)
	      || (  NUMABS_HKL_LIMIT < tNumabsRec.nH)
	      || ( -NUMABS_HKL_LIMIT > tNumabsRec.nK)
	      || (  NUMABS_HKL_LIMIT < tNumabsRec.nK)
	      || ( -NUMABS_HKL_LIMIT > tNumabsRec.nL)
	      || (  NUMABS_HKL_LIMIT < tNumabsRec.nL) )
	    {
	      nErrorHKL++;
	    }
	  nRejected++;
	}
    } // end i loop

  // Close the output file

  if (pFTextOutFile) fclose(pFTextOutFile);

  cout << "\nReflections in input file: " << poReflnlist->nGetNumReflns()
       << "\nReflections excluded:      " << nExcluded
       << "\nReflections rejected:      " << nRejected
       << "\nReflections written:       " << nWritten
       << endl << flush;

  delete poReflnlist;

  if (0 < nErrorHKL)
    {
      cout << "WARNING! There were " << nErrorHKL << " reflections with "
           << "HKL values outside the range of -" << NUMABS_HKL_LIMIT 
           << " < h,k, or l < " << NUMABS_HKL_LIMIT 
           << "\nallowed by numabs.\n" << flush;
    }

  return (nStat);
}
