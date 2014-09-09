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
// dtnumabs2refln.cc     Initial author: J.W. Pflugrath           11-Aug-1997
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
//    dtnumabs2refln reads from the command line the reflnlist filename
//    if present (default: *dtprofit.ref), the name of the numabs output
//    file to read if present (default: numabs.hkl_new) and writes a new
//    d*TREK reflnlist file (default: *dtnumabs.ref).
//
//    Syntax:   % dtnumabs2refln [input_reflnlist_file [output_file]]
//    Example:  % dtnumabs2refln dtprofit.ref numabs.hkl_new dtnumabs.ref
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
#include "CrclHelper.h"

#define DT_MAX_ORIENT   100

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
  int  nRead(FILE* f);
} tagNumabsRec;

int tagNumabsRec::nRead(FILE* f) {
  int nStat;
	(void) fscanf(f,"%d %d %d %f %f %f %f %f %f %f %f %f\n",
        &nH, &nK, &nL, &fIntensity, &fSigmaI, 
        &(a3fSource[0]), &(a3fSource[1]), &(a3fSource[2]),
        &(a3fScatt[0]),  &(a3fScatt[1]),  &(a3fScatt[2]), &fTransFactor);
	nStat = ferror(f);  // bug, no clearerr() call
	return (nStat);
};

#ifdef SSI_PC
#include "CrclIncludes.h"
#define main dtnumabs2refln_main
#endif

int main(int   argc,    // Command line argument count
         char *argv[])  // Pointers to command line args

{
    int      nStat;
    Cstring  sTemplateFile            = "dtprofit.ref";
    Cstring  sInNumabs      = "numabs.hkl_new";
    Cstring  sOut           = "dtnumabs.ref";
    Creflnlist *poReflnlist = NULL;
    Creflnlist *poReflnlistOut = NULL;
    Crefln     *poRefln;

#if !defined(VC6)
    using std::cout;
    using std::flush;
    using std::endl;
#endif
    cout << "\ndtnumabs2refln: Copyright (c) 2006 Rigaku\n" << flush;
    cout << D_K_DTREKVersion << endl;

#ifdef SSI_PC
    Cstring sCCAppVersion = (const char*) argv[argc-1];
	if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1) {
		cout << sCCAppVersion << endl;
		argc--;
		cout << CCrclHelper::GetInstance()->GetTimestamp() << endl;
		argv[argc] = NULL;
	}
#endif

    // Copy command line to output log
    cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << endl << flush;

    // Parse command line arguments, if any

    argc--; argv++; // Skip program name

    sTemplateFile  = sDtrekGetPrefix() + sTemplateFile;
    sOut = sDtrekGetPrefix() + sOut;

    if( argc < 3 ) {
        cout << "Usage: dtnumabs2refln sInputNumabsFile sOutputDtrekFile sOutputTemplateFile\n";
		#ifdef SSI_PC
			return 1;
		#else
			exit (1);
		#endif
    }
    // First argument is numabs hkl input file.
    sInNumabs = (const char*) argv[0];
    // Second argument is d*Trek output file.
    sOut = (const char*) argv[1];
    // Third argument is template d*Trek reflection file.
    sTemplateFile = (const char*) argv[2];

    cout << "dtnumabs2refln: convert "
        << sInNumabs  << " to "
        << sOut << '\n';

    // Read the reflnlist

    cout << "Reading " << sTemplateFile << " ..." << endl << flush;
    poReflnlist = new Creflnlist (sTemplateFile);

    if (!poReflnlist->bIsAvailable()) {
        cout << "dtnumabs2refln: template reflnlist " << sTemplateFile << " unavailable!\n" << flush;
        delete poReflnlist;
		#ifdef SSI_PC
			return 1;
		#else
			exit (1);
		#endif
    }

    if (1 > poReflnlist->nGetNumReflns()) {
        cout << "dtnumabs2refln: reflnlist " << sTemplateFile << " has no reflections!\n" << flush;
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
    if (0 > nFI_nErrorFlag) {
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
    if (0 > poReflnlist->m_nFI_sBatch)
        nStat++;
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
        cout << "dtnumabs2refln: template reflnlist does not have all the required fields!\n" << flush;
        delete poReflnlist;
		#ifdef SSI_PC
			return 2;
		#else
			exit (2);
		#endif
    }

    int i;
    int nSize;
    int nRead     = 0;
    int nRejected = 0;
    int nExcluded = 0;
    float fScale;

    // Open input numabs file
    int nLength   = sInNumabs.length();
    int nFile     = 1;
    FILE *pInNumabs = fopen( sInNumabs, "r" );
    //(void) dskbor(&nFile, sInNumabs.string(), &nLength, &nStat);

    if ( !pInNumabs ) {
        cout << "ERROR opening input file " << sInNumabs << "!\n" << flush;
        delete poReflnlist;
        return (nStat);
    }

    // Create the output reflnlist and add 2 new columns
    poReflnlistOut      = new Creflnlist(*poReflnlist);
    int nFI_fTgeos      = poReflnlistOut->nExpandGetField("fTgeos");
    int nFI_fScale      = poReflnlistOut->nExpandGetField("fScale");

    // Copy the old list to the new (empty) list and then use the new list
    // from now on.
    poReflnlistOut->nInsertListFrom(*poReflnlist);
    delete poReflnlist;
    poReflnlist = poReflnlistOut;
    poReflnlistOut = NULL;

    int nFI_fOtherInt   = poReflnlist->nGetFieldIndex(D_K_ProfitfOtherInt,
        eReflnField_float_type);
    int nFI_fOtherSig   = poReflnlist->nGetFieldIndex(D_K_ProfitfOtherSig,
        eReflnField_float_type);

    char buffer[512];
    int tempInt[3];
    float tempFloat[9];
    for (i = 0; (0 == nStat) && (i < poReflnlist->nGetNumReflns()); i++)
    {
        poRefln                  = poReflnlist->poGetRefln(i);

#define NUMABS_HKL_LIMIT  512

        if (   (-999. < poRefln->fGetField(nFI_fS0vec[2]))
            && (-999. < poRefln->fGetField(nFI_fS0vec[1]))
            && (-999. < poRefln->fGetField(nFI_fS0vec[0]))
            && ( 999. > -poRefln->fGetField(nFI_fSvec[2]))
            && ( 999. > -poRefln->fGetField(nFI_fSvec[1]))
            && ( 999. > -poRefln->fGetField(nFI_fSvec[0]))
            && ( -NUMABS_HKL_LIMIT < poRefln->nGetH())
            && (  NUMABS_HKL_LIMIT > poRefln->nGetH())
            && ( -NUMABS_HKL_LIMIT < poRefln->nGetK())
            && (  NUMABS_HKL_LIMIT > poRefln->nGetK())
            && ( -NUMABS_HKL_LIMIT < poRefln->nGetL())
            && (  NUMABS_HKL_LIMIT > poRefln->nGetL()) )
        {
            // Legitimate reflection at first glance
            
            if( !fgets( buffer, 512, pInNumabs ) ) {
                cout << "ERROR reading numabs hkl file!\n" << flush;
            }
            else {
                nRead++;
                sscanf( buffer, "%d %d %d %f %f %f %f %f %f %f %f %f ", &tempInt[0], &tempInt[1], &tempInt[2], 
                    &tempFloat[0], &tempFloat[1], &tempFloat[2], &tempFloat[3], &tempFloat[4], &tempFloat[5], 
                    &tempFloat[6], &tempFloat[7], &tempFloat[8] );
                poRefln->vSetIntensity( poRefln->fGetIntensity() / tempFloat[8] );
                poRefln->vSetSigmaI( poRefln->fGetSigmaI() / tempFloat[8] );
                poRefln->vSetField( nFI_fTgeos, tempFloat[8] );
                poRefln->vSetField( nFI_fScale, 1.0f );
            }
        }
        else {
            nRejected++;
        }
    } // end i loop
    
    // Close the output file

    //(void) dskbcr(&nFile, &nStat);
    fclose( pInNumabs );

    // Write out the reflnlist

    cout << "Writing " << sOut << " ..." << endl;
    nStat = poReflnlist->nWrite(sOut);

    cout << "\nReflections in input file: " << poReflnlist->nGetNumReflns()
        << "\nReflections excluded:      " << nExcluded
        << "\nReflections rejected:      " << nRejected
        << "\nReflections written:       " << nRead
        << endl << flush;

    delete poReflnlist;
    return (nStat);
}
