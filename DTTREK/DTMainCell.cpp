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
// DTMainCell.cpp       Initial author: RB          10-Sep-2004
// This file contains the member functions of class CDTMainCell

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

#include "DTMainCell.h"

#include "Cstring.h"
#include "Dtrek.h"
#include "Cimage.h"
#include "dtrekvec.h"
#include "dtsvd.h"
#include "Crefln.h"
#include "Creflnlist.h"
#include "Cindex.h"

#include "dtcell.h"

#include "Cspacegroup_check.h"

#include <string.h>
#include <stdio.h>
#include <ctype.h>


#ifdef SSI_PC
#include "CrclHelper.h"
#include "CrclIncludes.h"
#endif

#if !defined(VC6)
using std::cin;
using std::cout;
using std::endl;
using std::flush;
#endif

// RB: commenting this out, because of a compile warning under Windows
#ifndef WIN32
    #define ERROR {  cout << "Bad Logic Line " << __LINE__ << " of file '" << __FILE__ "'\n"; GETCH; exit(0); };
#endif

//template <class T> void swap(T& x,T& y) { T z; z=x; x=y; y=z; };
#define STR(x) #x

#define strlwr(x) { int nx; for (nx=0; nx < (int)strlen(x);nx++) x[nx]=tolower(x[nx]);};
#define strupr(x) { int nx; for (nx=0; nx < (int)strlen(x);nx++) x[nx]=toupper(x[nx]);};


char* cpHelp =
"dtcell Usage:                                                        \n"
"dtcell [input.head] [input.ref] [options]                            \n"
"Command line options: Description:                                   \n"
"                                                                     \n"
"    -spacegroup sIdent            Set the spacegroup to sIdent.      \n"
"                                  sIdent is a spacegroup number or a \n"
"                                  spacegroup identification string.  \n"
"                                  Still runs systematic absence check\n"
"                                  but uses result to resolve         \n"
"                                  ambiguities. Laue check disabled.  \n"
"                                                                     \n"
"    -reduce                       Calls cell reduction code to reduce\n"
"                                  input cell.  Disables laue check   \n"
"                                  and absence check.                 \n"
"                                                                     \n"
"    -lauegroup sIdent             Set the Laue group to sIdent.     \n"
"                                  (See -stringhelp).                 \n"
"                                  Disables Laue check.               \n"
"                                                                     \n"
"    -maxlauegroup sIdent          Provide the highest laue group to \n"
"                                  check.  Laue check will not use    \n"
"                                  groups with higher symmetry.       \n"
"                                                                     \n"
"    -maxreflections nLines        Max # of reflections to use.       \n"
"                                                                     \n"
"    -prompt                       Prompt at ambiguities (Default).   \n"
"                                                                     \n"
"    -noprompt                     Do not prompt for ambiguities.     \n"
"                                  Terminate at uncertainties.        \n"
"                                                                     \n"
"    -centric -acentric            Set centric or acentric.  Disables \n"
"                                  centricity check.                  \n"
"                                                                     \n"
"    -transform sIdent             Transform the input reflection list\n"
"                                  using the given transform. sIdent  \n"
"                                  can also be a target reflection    \n"
"                                  file.                              \n"
"                                  (See -stringhelp).                 \n"
"                                                                     \n"
"    -lauescale [nMinRefsPerBatch] Batch scale reflection data during \n"
"                                  Laue check. Use no fewer than      \n"
"                                  nMinRefsPerBatch reflns. per batch.\n"
"                                  Default: 50 per batch.             \n"
"                                                                     \n"
"    -sigma fSigma                 Minimum I/sigma for a reflection to\n"
"                                  be used in Laue check.             \n"
"                                  Default: 3                         \n"
"                                                                     \n"
"    -asigma fASigma               Maximum I/sigma for a reflection to\n"
"                                  be considered absent.              \n"
"                                  Default: 4                         \n"
"                                                                     \n"
"    -maxrmerge fPoor [fMaxRatio] [fIdeal]                            \n"
"                                  Specifies Laue group R-merge test. \n"
"                                  Any Laue group with R-merge larger \n"
"                                  or equal to fPoor is rejected.     \n"
"                                  Any Laue group with R-merge smaller\n"
"                                  or equal to fIdeal is accepted.    \n"
"                                  All other Laue groups are treated  \n"
"                                  based on on the ratio of their     \n"
"                                  R-merge to the R-merge of the Laue \n"
"                                  group -1. If that ratio is larger  \n"
"                                  than fMaxRatio, the group is rejected,\n"
"                                  else the group is accepted.        \n"
"                                  Defaults: 0.20 1.5 0.08 resp.      \n"
"                                                                     \n"
"    -anom                         Laue group determination assumes   \n"
"                                  I+ != I-                           \n"
"                                                                     \n"
"    -chiral                       Use only spacegroups that are      \n"
"                                  compatible with chiral structures. \n"
"                                                                     \n"
"    -chiralunknown                Use all spacegroups (opposite of   \n"
"                                  -chiral option).                   \n" 
"                                                                     \n"
"    -cell a b c alp bet gam       Sets the cell to a given cell.     \n"
"                                  (Used only in the Laue check)      \n"
"                                  Default: (10,10,10,90,90,90)       \n"
"                                                                     \n"
"    -maxresid fResidMax           Largest permissible least squares  \n"
"                                  residual for a reduced cell display\n"
"                                  Default: 5.0.                      \n"
"                                                                     \n"
"    -rejects [nMaxRejects]        Writes out dtcell_rejects.ref      \n"
"                                  containing information about       \n"
"                                  rejected reflections.              \n"
"                                  If nMaxRejects is specified, the   \n"
"                                  maximum number of reflections for  \n"
"                                  any unique zone/systematic absence \n"
"                                  category will be nMaxRejects.      \n"
"                                  nMaxRejects defaults to 20.        \n"
"                                                                     \n"
"    -orient theta1 theta2 theta3  Sets the orientation of the cell.  \n"
"                                                                     \n"
"    -sigmas s1 s2 s3 s4 s5 s6     Sets cell parameter sigma values.  \n"
"                                                                     \n"
"    -lattice STRING               Sets the Bravais lattice type.     \n"
"                                  This will be used in cell reduction\n"
"                                  to convert to a primitive cell if  \n"
"                                  the input cell is centered.        \n"
"                                                                     \n"
"    -ref OUTPUT.REF               Explicitly supply an output .ref   \n"
"                                  file for file transformations and  \n"
"                                  reindexing operations.             \n"
"                                  Enables re-indexing.               \n"
"                                  Default: *dtcell.ref               \n"
"                                                                     \n"
"    -head OUTPUT.HEAD             Explicitly supply an output .head  \n"
"                                  file for output.                   \n"
"                                  Default: *dtcell.head              \n"
"                                                                     \n"
"    -twin nTwinID                 Restrict calculations to a         \n"
"                                  particular twin.                   \n"
"                                                                     \n"
"    -twinlaw fRot [sRotAxes]      Create/add twin law.  The sRotAxes \n"
"                                  string should contain one or two   \n"
"                                  of 'a*', 'b*', or 'c*'			  \n"
"    -stringhelp                   Prints help on Laue,Transform.     \n"
"                                  and Lattice strings used in        \n"
"                                  -lauegroup -transform and          \n"
"                                  -lattice.                          \n"
;


char* cpTransHelp =
"    Transformation Strings.  (Ex:  Flip a and b with \"-transform AB\")    \n"
"                                                                           \n"
"    String Choice   Action                                                 \n"
"    -------------   ------                                                 \n"
"    I     | ABC     Identity transformation.                               \n"
"    P1    | CAB     Axis Permutation (abc) -> (cab).                       \n"
"    P2    | BCA     Axis Permutation (abc) -> (bca).                       \n"
"    AB    | BAC     Flip a and b axis (Right handed).                      \n"
"    BC    | ACB     Flip b and c axis (Right handed).                      \n"
"    AC    | CBA     Flip a and c axis (Right handed).                      \n"
"    C2-C1           Monoclinc Cell type 2 to type 1.                       \n"
"    C3-C1           Monoclinc Cell type 3 to type 1.                       \n"
"    C1-C2           Monoclinc Cell type 1 to type 2.                       \n"
"    C3-C2           Monoclinc Cell type 3 to type 2.                       \n"
"    C1-C3           Monoclinc Cell type 1 to type 3.                       \n"
"    C2-C3           Monoclinc Cell type 2 to type 3.                       \n"
"    mI-mC           I centered to C centered monoclinic.                   \n"
"    mC-mI           C centered to I centered monoclinic.                   \n"
"    hR-hP           Rhombohedral to Hexagonal.                             \n"
"    hP-hR           Hexagonal to Rhombohedral.                             \n"
"    hP-oC1          Hexagonal to Orthorhombic centered. (Cell choice 1)    \n"
"    hP-oC2          Hexagonal to Orthorhombic centered. (Cell choice 2)    \n"
"    hP-oC3          Hexagonal to Orthorhombic centered. (Cell choice 3)    \n"
"    oC1-hP          Orthorhombic centered to Hexagonal. (Cell choice 1)    \n"
"    oC2-hP          Orthorhombic centered to Hexagonal. (Cell choice 2)    \n"
"    oC3-hP          Orthorhombic centered to Hexagonal. (Cell choice 3)    \n"
"                                                                           \n"
"                        --- OR ---                                         \n"
"                                                                           \n"
"    Another reflection file can act as a transform:                        \n"
"    Example:  dtcell file1.ref -transform file2.ref                        \n"
"                                                                           \n"
"    This finds a transform and scale factor for file1.ref that minimizes   \n"
"    the R-merge between the data sets.                                     \n"
"                                                                           \n"
"                        --- OR ---                                         \n"
"                                                                           \n"
"    h00 h01 h02 h10 h11 h12 h20 h21 h22 [fScale]                           \n"
"    Where the number h00 .. h22 are the coefficeints of an hkl             \n"
"    transformation matrix:                                                 \n"
"                                                                           \n"
"   [h'] [  h00  h01  h02 ] [h]                                             \n"
"   [k']=[  h10  h11  h12 ]*[k]                                             \n"
"   [l'] [  h20  h21  h22 ] [l]                                             \n"
"                                                                           \n"
"   Each of number h00 .. h22 can be in the form of a floating point or     \n"
"   a fraction.  (Examples:  1.0 0.5 -2/3 1.0/3.0 )                         \n"
"                                                                           \n"
"   If the number fScale is provided, the scale factor is applied to each   \n"
"   intensity and sigma value in the reindexed reflection file.             \n"
"                                                                           \n"
;

char* cpLaueHelp =
"    Laue Strings (Ex: -lauegroup 4/mmm)                                    \n"
"                                                                           \n"
"    Laue Strings                                                           \n"
"    ------------                                                           \n"
"    -1     2/m   mmm     4/m   4/mmm   -3                                  \n"
"    -3m1  -31m   6/m   6/mmm   m-3   m-3m                                  \n"
"                                                                           \n"
;

char* cpLatticeHelp =
"    Lattice Strings                                                        \n"
"    ---------------                                                        \n"
"    aP mP mC oP oC oF oI hR hP tP tI cP cI cF                              \n" 
;

CDTMainCell::CDTMainCell()
{
    m_bCheckCentric     = true;           
    m_bCheckLaue        = true;         
    m_bSpacegroupCheck  = true;
    m_bReduce           = false;
}
///////////////////////////////////////////////////////////////////////////
int CDTMainCell::nExecute(unsigned int argc, char* argv[])
{
    Creflnlist  *poReflnlistIn=NULL;
	Cimage_header *poHeader=NULL;

	Cstring sTemp;
    Claue oLaue;
    int nx,ny;
    float fTemp;
    double f0;
    char cpBuf[100];
    bool bTemp;

    vDtrekSetModuleName("dtcell");
    vPrintCopyrightInfo();
    //    cout << "\ndtcell:  (c) 2006-1995 Rigaku\n";
    //    cout << D_K_DTREKVersion << endl << flush;

#ifdef SSI_PC
	Cstring sCCAppVersion = (const char*) argv[argc-1];
	if( sCCAppVersion.find("CrystalClear") != -1 || sCCAppVersion.find("SCXmini") != -1) {
		cout << sCCAppVersion;
		argc -= 2;
		cout << CCrclHelper::GetInstance()->GetTimestamp() << endl;
		argv[argc] = NULL;
	}
#endif

    cout << "Command line:\n " << sGetCommandLine(argc, argv, 71) << endl << endl << flush;

    sTemp = (const char*) argv[1];
    if ("-help" == sTemp || "-h" == sTemp || argc<2) {
        printf(cpHelp);
        return 0;
    };

	argc--; argv++; // Skip program name
	int nArg = 0;
	int nStat= 0;

	int		  nMaxReflectionsToRead = -1;	// Maximum # of reflections to read in.
    int       nTransformIndex;              // What is the index of the transform desired?  -1==User, -2==w.r.t. reflection list.
    Creflnlist* poRefTransform = NULL;      // If we are transforming w.r.t. a reference reflection list, what is the reflection list?
    double     fTransformScale = 1.0;       // Scale factor for reindexing transformation.
    double     a3x3fUserTrans[3][3];        // User specified transformation matrix.
    bool      bDoTransform =  FALSE;        // If we are to do a transform, all other options are ignored.
    bool      bPrompt       = TRUE;         // Prompt the user during ambiguity step.    
    bool      bReindexingDone = FALSE;      // Returned by nChooseSpaceGroup() and eLaueCheck()
	bool      bAddTwinLaw      = FALSE;		// Are we adding a twin law
    Cstring   sOutputHeader = sDtrekGetPrefix() + "dtcell.head"; // Output header name.
    Cstring   sHeaderName;                  // Input header name.
    Cstring   sOutputRef= sDtrekGetPrefix() + "dtcell.ref";     // Output reflection file.
    Cstring   sRefName;                     // Input reflection file.
	Cstring   sTwinLaw;						// Twin law for making such.
	double    fTwinLawRot;					// Rotation for twin laws.
    
    int       nTwinID = 1;                  // Twin ID.

#ifdef SSI_PC
    bool      bGetInputCell = FALSE;          // Call to get input cell to reduce
	bool	  bAutoMode = FALSE;			  // Flag whether to indicate automode operation
#endif

    // Prescan of the -twin option.

    for (nArg = 0; nArg < (int)argc; nArg++) {
        if (!strcmp(argv[nArg],"-twin")) {
            argv[nArg][0] = 0;
            if ((nArg+1 >= (int)argc) || (1!=sscanf(argv[nArg+1],"%d",&nTwinID))) {
                printf("Incorrect syntax with option '-twin'\n");
                goto main_termination;
            };
            argv[nArg+1][0] = 0;
            oLaue.nTwinID = nTwinID;
            nArg++;
        };
    };
    nArg = 0;
    
    // See if a header has been supplied.  If so, we read the cell information from this header.

    if (nArg < (int)argc) {
        if ((strlen(argv[nArg])>5) && (!strcmp(&argv[nArg][strlen(argv[nArg])-5],".head"))) {
            sHeaderName=argv[nArg];
            printf("Header \"%s\" opened.\n",sHeaderName.string());
            nArg++;
            if (poHeader)
                delete poHeader;
            poHeader = new Cimage_header(sHeaderName);   // Create image header
            if (!poHeader->bIsAvailable()) {
                printf("Could not access header \"%s\"\n",sHeaderName.string());
                goto main_termination;
            };
            
            // Read in the cell parameters from the image header.
            if (oLaue.nSetCrystal(*poHeader))
                goto main_termination;
            oLaue.nSetChirality(TRUE);
        };
    };


	// Get the reflection list

    sRefName="";
	if (nArg < (int)argc) {
        // If the first argument looks like an option, then skip this.
        if (argv[nArg][0]!='-') {
            sRefName = argv[nArg];
            nArg++;
        };
    } else {  printf(cpHelp);
            goto main_termination;};

    // Read in the reflection list.

    if (sRefName.string()[0]!=0) {
        cout << "Reading Reflection file " << sRefName << " ..." << endl << flush;
	    poReflnlistIn=new Creflnlist(sRefName,nMaxReflectionsToRead);
	    if (!poReflnlistIn->bIsAvailable()  || (1 > poReflnlistIn->nGetNumReflns()) ) {
            printf("Problem with input reflection list: %s",sRefName.string());
	        goto main_termination;
        }
        cout << flush;

        // Process the crystal object if it is available.
        // But only if there is not a previous header available.
        if (!poHeader) {
            Cimage_header* poRefHeader;
            Ccrystal* poCrystal;
            poCrystal = poReflnlistIn->poGetCrystal();
            if (poCrystal) {
                poRefHeader = poReflnlistIn->poGetHeader();
                if (oLaue.nSetCrystal(*poRefHeader)) {
                    printf("WARNING:  Problem with embedded crystal object in reflection list '%s'.\n"
                        "Ignoring the embedded crystal!\n",sRefName.string());
                    delete poRefHeader;
                } else {
                    oLaue.nSetChirality(TRUE);
                    if (poHeader)
                        delete poHeader;
                    poHeader = poRefHeader;
                };
                delete poCrystal;
            };
        };
    };

    if ((poReflnlistIn) && (poReflnlistIn->m_nFI_nTwinID>=0)) {
        for (nx=0;nx<poReflnlistIn->nGetNumReflns();nx++) {
            if ((*poReflnlistIn)[nx].nGetField(poReflnlistIn->m_nFI_nTwinID)!=nTwinID) {
                poReflnlistIn->nDelete(nx);
                nx--;
            };
        };
    };            


    ////////////////////////////////////
    ////// Begin of reading in options
    ////////////////////////////////////
       

    while (nArg< (int)argc) {
        printf("Command line string: >>%s<<\n",argv[nArg]);


        if (argv[nArg][0]==0) {
        } else if (!strcmp(argv[nArg],"-spacegroup")) {
            nArg++;

            if (nArg==argc) {
                printf("Could not parse -spacegroup\n");
                goto main_termination;
            };
            if (oLaue.m_eLaue != LAUE_UNKNOWN) {
                printf("Laue type already set!\n");
                goto main_termination;
            };
            // Check for string name first.
            strlwr(argv[nArg]);
            for (nx=0; nx<230; nx++) {
                strcpy(cpBuf,g_cpSpaceGroupNames[nx]);
                strlwr(cpBuf);
                if (!strcmp(cpBuf,argv[nArg]))
                    break;
            };
            if (nx == 230) {
                // Need to parse the space group number instead.
                if (1!= sscanf(argv[nArg],"%d",&ny)) {
                    printf("Could not parse -spacegroup number\n");
                    goto main_termination;
                };
                if (ny>230) ny=230;
                if (ny<1) ny=1;
                oLaue.m_nSpaceGroup=ny;
            } else 
                oLaue.m_nSpaceGroup=nx+1;

            m_bSpacegroupCheck = TRUE;
            // We don't need to check the Laue group.
            m_bCheckLaue=FALSE;
            
            // But we do have to have a valid Laue group if spacegroups
            // were not disabled with -spacegroup -1
            
            oLaue.m_eLaue=oLaue.eLauefromSpacegroup(oLaue.m_nSpaceGroup);
            // Set the lattice number.
            Cspacegroup oSpace;
            oSpace.vSet(oLaue.m_nSpaceGroup);
            oLaue.a2cLattice[0]=oSpace.cGetClass();
            oLaue.a2cLattice[1]=g_cpSpaceGroupNames[oLaue.m_nSpaceGroup-1][0];
            printf("Space group set to: %s (#%d)\nLaue check disabled.\nLaue group set to: %s\nLattice set to %s\n",
                g_cpSpaceGroupNames[oLaue.m_nSpaceGroup-1],oLaue.m_nSpaceGroup,g_cpLaueNameSmall[oLaue.m_eLaue],oLaue.a2cLattice);
        } else if (!strcmp(argv[nArg],"-reduce")) {
            m_bReduce = TRUE;

        } else if (!strcmp(argv[nArg],"-lauegroup")) {

            if (oLaue.m_eLaue != LAUE_UNKNOWN) {
                printf("Laue group already set!\n");
                goto main_termination;
            };
            nArg++;
            if (nArg==argc) {
                printf("Could not parse -lauegroup\n");
                goto main_termination;
            };
            strlwr(argv[nArg]);
            for (nx=0; nx<12; nx++) 
                if (!strcmp(g_cpLaueNameSmall[nx],argv[nArg]))
                    break;
                if (nx == 12) {
                    printf("Laue group name %s not found.\n\n",argv[nArg]);
                    printf(cpLaueHelp);
                    goto main_termination;
                } else {
                    oLaue.m_eLaue=(eLaueType) nx;
                    printf("Laue group set to: %s\nLaue check disabled.\n",g_cpLaueNameSmall[nx]);
                    m_bCheckLaue=FALSE;
                };
        } else if (!strcmp(argv[nArg],"-maxlauegroup")) {
            nArg++;
            if (nArg==argc) {
                printf("Could not parse -maxlauegroup\n");
                goto main_termination;
            };
            strlwr(argv[nArg]);
            for (nx=0; nx<12; nx++) 
                if (!strcmp(g_cpLaueNameSmall[nx],argv[nArg]))
                    break;
                if (nx == 12) {
                    printf("Laue group name %s not found.\n\n",argv[nArg]);
                    printf(cpLaueHelp);
                    goto main_termination;
                } else {
                    oLaue.eMaxLaueCheck=(eLaueType) nx;
                    printf("Maximum Laue group to check set to: %s\n",g_cpLaueNameSmall[nx]);
                };

        } else if (!strcmp(argv[nArg],"-prompt")) {
            bPrompt=TRUE;
            printf("Prompting ON.\n");
        } else if (!strcmp(argv[nArg],"-noprompt")) {
            printf("Prompting OFF.\n");
            bPrompt=FALSE;
		} else if (!strcmp(argv[nArg],"-twinlaw")) {
			if (!poHeader) {
				printf("Cannot use -twinlaw command without a header file.\n");
				goto main_termination;
			} else if ((nArg + 1 < (int)argc) && (1 == sscanf(argv[nArg+1],"%lf",&f0))) {
				nArg++;
				fTwinLawRot = f0;
				if ((nArg + 1 < (int)argc) && (argv[nArg+1][0]!='-')) {
					sTwinLaw = argv[nArg + 1];
					nArg++;
				} else
					sTwinLaw = "";
			} else {
				printf("Cannot parse -twinlaw\n");
				goto main_termination;
			};
			bAddTwinLaw = TRUE;
			if (sTwinLaw.length())
				printf("Twin Law '%s' to be applied.\n",sTwinLaw.string());
			else
				printf("Twin Laws to be removed.\n");
			
        } else if (!strcmp(argv[nArg],"-centric")) {
            oLaue.m_bIsCentric=TRUE;
            m_bCheckCentric=FALSE;
            printf("Space group is Centric.\nCentricity check disabled.\n");
        } else if (!strcmp(argv[nArg],"-acentric")) {
            oLaue.m_bIsCentric=FALSE;
            m_bCheckCentric=FALSE;
            printf("Space group is Acentric.\nCentricity check disabled.\n");                
        } else if (!strcmp(argv[nArg],"-chiral")) {
            oLaue.m_bIsChiral = true;
            printf("Space groups displayed will be compatible with chiral structures.\n");

        } else if (!strcmp(argv[nArg],"-chiralunknown")) {
            oLaue.m_bIsChiral = false;
        } else if (!strcmp(argv[nArg],"-transform")) {
            nArg++;
            fTransformScale = 1.0;
            if (nArg==argc) {
                printf("Could not parse -transform\n");
                goto main_termination;
            };
            sTemp = argv[nArg];
            sTemp.upcase();
            for (nx=0,ny=0;nx<NUM_REINDEX_TRANSFORMS; nx++) {
                while (g_cpTransformStrings[ny])
                    if (sTemp==g_cpTransformStrings[ny])
                        break;
                    else
                        ny++;
                if (g_cpTransformStrings[ny])
                    break;
                ny++;
            };
            if (nx==NUM_REINDEX_TRANSFORMS) {


                if (bFileExists(sTemp=argv[nArg])) {
                    // The user wants to reindex w.r.t. a file.
                    printf("Reading %s ...\n",sTemp.string());
                    poRefTransform = new Creflnlist(sTemp);
                    if (!poRefTransform->bIsAvailable()) {
                        printf("Could not open transformation '%s'\n",sTemp.string());
                        goto main_termination;
                    };
                    nTransformIndex = -2;
                    printf("Reindexing w.r.t. reflection list '%s'\n",sTemp.string());
                } else {
                    // The user must have provided the set of numbers to reindex with.
                    if (nArg+9> (int)argc) {
                        printf("Could not parse -transform option.\n");
                        goto main_termination;
                    };
                    for (nx=0;nx<9;nx++) {
                        sTemp = argv[nArg+nx];

                        fTemp = (float)(*(nx+&a3x3fUserTrans[0][0]));
                        if (!sTemp.bParseFloat(fTemp)) {
                            printf("Could not parse -transform (near %s)\n",sTemp.string());
                            goto main_termination;
                        };
                        *(nx+&a3x3fUserTrans[0][0]) = fTemp;
                    };
                    vTranMat3D(a3x3fUserTrans);
                    nArg+=9-1;
                    if ((nArg+1< (int)argc) && (1==sscanf(argv[nArg+1],"%lf",&f0))) {
                        fTransformScale = f0;
                        nArg++;
                    };
                    nTransformIndex=-1;
                    printf("Reindexing with user specified transformation.\n");
                };
            } else {
                nTransformIndex=nx;
                printf("Reindexing transformation set to: %s\n",g_cpTransforms[nTransformIndex]);
            };

            bDoTransform=TRUE;
            
            
        } else if (!strcmp(argv[nArg],"-asigma")) {
            nArg++;
            if ((nArg==argc) || (1!=sscanf(argv[nArg],"%lf",&f0))) {
                printf("Could not parse -asigma\n");
                goto main_termination;
            };
            oLaue.fAbsentSigma=f0;
            printf("Absence I/Sigma(I) set to %f\n",oLaue.fAbsentSigma);

        } else if (!strcmp(argv[nArg],"-sigma")) {
            nArg++;
            if ((nArg==argc) || (1!=sscanf(argv[nArg],"%lf",&f0))) {
                printf("Could not parse -sigma\n");
                goto main_termination;
            };
            oLaue.fLaueSigma=f0;
            printf("Laue I/Sigma(I) set to %f\n",oLaue.fLaueSigma);
        } else if (!strcmp(argv[nArg],"-lauescale")) {
            oLaue.bScaleForLaue = TRUE;
			if ((nArg+1< (int)argc) && (1==sscanf(argv[nArg+1],"%d",&nx))) {
				oLaue.nScaleForLaueMinRefs = nx;
				nArg++;
			};
		} else if (!strcmp(argv[nArg],"-maxreflections")) {
			if ((nArg+1>= (int)argc) || (1!=sscanf(argv[nArg+1],"%d",&nMaxReflectionsToRead))) {
				printf("Could not parse -maxreflections.\n");
				goto main_termination;
			};
        } else if (!strcmp(argv[nArg],"-maxrmerge")) {
            nArg++;
            if ((nArg>= (int)argc) || (1!=sscanf(argv[nArg],"%lf",&f0))) {
                printf("Could not parse -maxrmerge\n");
                goto main_termination;
            };
            oLaue.fRMergeAwful=f0;
            printf("Maximum R-merge is %f\n",oLaue.fRMergeAwful);
            if ((nArg+1< (int)argc) && (1==sscanf(argv[nArg+1],"%lf",&f0))) {
                oLaue.fMaxRMergeRatio = f0;
                nArg++;
                printf("Maximum R-merge ratio is %f\n",oLaue.fMaxRMergeRatio);
            };            
            if ((nArg+1< (int)argc) && (1==sscanf(argv[nArg+1],"%lf",&f0))) {
                oLaue.fRMergeGreat = f0;
                nArg++;
                printf("Minimum R-merge is %f\n",oLaue.fRMergeGreat);
            };            
            
        } else if (!strcmp(argv[nArg],"-anom")) {
            oLaue.bAnom = TRUE;
            printf("Laue check assumes I+ != I-\n");
        } else if (!strcmp(argv[nArg],"-cell")) {
            nArg++;
            for (nx=0;nx<6; nx++) {
                if ((nArg+nx>= (int)argc) || (1!=sscanf(argv[nArg+nx],"%lf",&f0))) {
                    printf("Could not parse -cell\n");
                    goto main_termination;
                } else
                    oLaue.fCellParams[nx] = f0;
            };

            printf("Cell:        [ %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f]\n"
                   "Orientation: [ %6.2f %6.2f %6.2f ]\n",                    
                   oLaue.fCellParams[0],oLaue.fCellParams[1],oLaue.fCellParams[2],oLaue.fCellParams[3],oLaue.fCellParams[4],oLaue.fCellParams[5],
                   oLaue.fOrientAngles[0],oLaue.fOrientAngles[1],oLaue.fOrientAngles[2]                   
                   );
            nArg+=5;

        } else if (!strcmp(argv[nArg],"-maxresid")) {
            nArg++;
            if ((nArg>= (int)argc) || (1!=sscanf(argv[nArg],"%lf",&f0))) {
                printf("Could not parse -maxresid");
                goto main_termination;
            };
            oLaue.fResidMax = f0;
        } else if (!strcmp(argv[nArg],"-rejects")) {
            if ((nArg+1==argc) || (1!=sscanf(argv[nArg+1],"%d",&nx)))
                nx = -1;
            else
                nArg++;
            if (nx == -1)
                nx = 20;
            oLaue.bCollectRejects = TRUE;
            oLaue.nMaxCollectRejects = nx;

        } else if (!strcmp(argv[nArg],"-orient")) {
            nArg++;
            for (nx=0;nx<3; nx++) {
                if ((nArg+nx>= (int)argc) || (1!=sscanf(argv[nArg+nx],"%lf",&f0))) {
                    printf("Could not parse -orient\n");
                    goto main_termination;
                } else
                    oLaue.fOrientAngles[nx] = f0;
            }
            printf("Cell orientation set to [ %6.2f %6.2f %6.2f ]\n",
                oLaue.fOrientAngles[0],oLaue.fOrientAngles[1],oLaue.fOrientAngles[2]);
            nArg+=2;

        } else if (!strcmp(argv[nArg],"-sigmas")) {
            nArg++;
            for (nx=0;nx<6;nx++) {
                if ((nArg+nx>= (int)argc) || (1!=sscanf(argv[nArg+nx],"%lf",&f0))) {
                    printf("Could not parse -sigmas\n");
                    goto main_termination;
                } else
                    oLaue.fCellSigmas[nx] = f0;
            };
            printf("Cell sigmas set to [ %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f ]\n",
                oLaue.fCellSigmas[0],oLaue.fCellSigmas[1],oLaue.fCellSigmas[2],
                oLaue.fCellSigmas[3],oLaue.fCellSigmas[4],oLaue.fCellSigmas[5]
                );
            nArg+=5;

        } else if (!strcmp(argv[nArg],"-lattice")) {
            nArg++;
            for (nx=0;nx<14;nx++) 
                if (!strcmp(argv[nArg],g_pcValidLattices[nx]))
                    break;
            if (nx == 14) {
                printf("Could not identify lattice '%s'\n",argv[nArg]);
                printf(cpLatticeHelp);
                goto main_termination;
            };
            strcpy(oLaue.a2cLattice,argv[nArg]);            

        } else if (!strcmp(argv[nArg],"-ref")) {
            nArg++;
            if (nArg==argc) {
                printf("Could not parse -outputref\n");
                goto main_termination;
            };
            sOutputRef=argv[nArg];
            printf("Output File Name set to \"%s\".\nRe-indexing enabled.\n",sOutputRef.string());
        } else if (!strcmp(argv[nArg],"-head")) {
            nArg++;
            if (nArg==argc) {
                printf("Could not parse -head\n");
                goto main_termination;
            };
            sOutputHeader=argv[nArg];
            printf("Output Header Name set to \"%s\".\n",sOutputHeader.string());
        
        } else if (!strcmp(argv[nArg],"-stringhelp")) {
            printf(cpTransHelp);
            printf(cpLaueHelp);
            printf(cpLatticeHelp);
            goto main_termination;
        } else if (!strcmp(argv[nArg],"-centricityonly")) {
            // UNDOCUMENTED OPTION.  For use only with Crystal-Clear
          m_bCheckCentric = TRUE;
          m_bCheckLaue    = FALSE;
          m_bSpacegroupCheck = FALSE;
        } else if (!strcmp(argv[nArg],"-laueonly")) {
            // UNDOCUMENTED OPTION.  For use only with Crystal-Clear
          m_bCheckCentric = FALSE;
          m_bCheckLaue = TRUE;
          m_bSpacegroupCheck = FALSE;
        } 
#ifdef SSI_PC		  
        else if (!strcmp(argv[nArg],"-automode"))
		{
            bAutoMode = TRUE;
		}
#endif
        else {
            printf("Unrecognized option >>%s<<\n",argv[nArg]);
            goto main_termination;
        };


        nArg++;

    }; // while

    

    ////////////////////////////////////
    ////// End of reading in options
    ////////////////////////////////////


    if (m_bReduce) {
        // Disable other checks.
        m_bCheckCentric = FALSE;
        m_bCheckLaue = FALSE;
        m_bSpacegroupCheck = FALSE;
    };

    printf("//////////////////////////////\n");
    printf("Finished reading command line.\n\n");

    printf("Input Cell:        [ %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f ]\n",
        oLaue.fCellParams[0],oLaue.fCellParams[1],oLaue.fCellParams[2],
        oLaue.fCellParams[3],oLaue.fCellParams[4],oLaue.fCellParams[5]
        );
    printf("Input Sigmas:      [ %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ]\n",
        oLaue.fCellSigmas[0],oLaue.fCellSigmas[1],oLaue.fCellSigmas[2],
        oLaue.fCellSigmas[3],oLaue.fCellSigmas[4],oLaue.fCellSigmas[5]
        );

    printf("Input Orientation: [ %6.2f %6.2f %6.2f ]\n\n",
        oLaue.fOrientAngles[0],oLaue.fOrientAngles[1],oLaue.fOrientAngles[2]
        ); 

restart_reduce:

    // If we are doing a cell reduction, then do it.
    if (m_bReduce) {
        CCellReduce  oReduce;


        int          iCell,iOrient;
        int          nRelist;
        double         a3x3fRealOrientMatIn[3][3];
        double         a3x3fRealOrientMatOut[3][3];   // This is used for the fFindOrigVecs, and then passed to nEquivRotations.  After that it is not used.
        double         a3x3fOrientMatIn[3][3];
        double         a3x3fOrientCenteredMatIn[3][3];// Centered version... used in calculating reindex matrices.      
        double         a3x3fRealMatOut[3][3];         // This guy does not have orientation information.  Loaded from cell selections.
        double         pfOrientAngles[3];             // Used whenever orientation angles are desired.
        
        Ccrystal     oCrystalIn;
        Ccrystal     oCrystalOut;
        
        oCrystalIn.vSetOrientAngles(oLaue.fOrientAngles);
        oCrystalIn.vSetCell(oLaue.fCellParams);
		oCrystalIn.nCalcOrientMatrix();
		oCrystalIn.vGetOrientMatrix(&a3x3fOrientCenteredMatIn[0][0]);

		oReduce.vConvertToPrimitive(oCrystalIn,oLaue.a2cLattice);	
		
		// Get the real orientation matrix of the original cell.
        // Cannot use nCalcGetRealMatrix() .. it does not encode orientation information.

        oCrystalIn.nCalcOrientMatrix();
        oCrystalIn.vGetOrientMatrix(&a3x3fOrientMatIn[0][0]);
        fInvMat3D(&a3x3fOrientMatIn[0][0],&a3x3fRealOrientMatIn[0][0]);
        vTranMat3D(a3x3fRealOrientMatIn);


        oReduce.m_fResidualMax = oLaue.fResidMax;   

        Crefln oRefln(&oReduce.oGetSolutions());
        
        oReduce.nAddCrystal(oCrystalIn);
        oReduce.nCalcGetBestLattices();

        // Allow user to choose best solution.
        do {
            double a6fUserCell[6];

            oReduce.nListResults(&iCell);
               
            printf("Enter one of the following:\n\t"
		   "a)          Solution number (1-14)\n\t"
		   "b)          New max residual (>=15)\n\t"
                   "c)          Other cell (6 cell parameters)\n\t"
		   "Enter)      Choose first cell in list.\n\t"
                   "Q,X)        Exit.\n\t"
		   "L)          List 44 lattice characters. (or 'L #X' to choose #X)"
                   "\n\tR #)        Set residual to #"
                   "\n\nChoice>");
                   
#ifdef SSI_PC
			CCrclHelper::GetInstance()->vBuildCmdCellChooseSoln( oReduce, oLaue);
#endif
//            if (bPrompt)
//                getline(cin,sTemp);
//            else 
//                sTemp = "";
            if (bPrompt)
			{
#ifdef SSI_PC
				if( bAutoMode )
				{
					sTemp = ""; // Choose first cell in list
				}
				else
				{
					getline(cin,sTemp);
				}
#else
				getline(cin,sTemp);
#endif
			}
            else
                sTemp = "";

#ifdef SSI_PC
			if (0 < sTemp.length()) {
				cout << sTemp << "\n";
			}
			else {
				cout << "<Enter>" << "\n";
			}
#else
            printf("\n");
#endif
	        sTemp.upcase();
            if (sTemp.string()[0]==0) {
                // Use the iCell value selected.
            } else if ((sTemp=="Q") || (sTemp=="X")) {
                iCell=-1;
                break;
            } else if (!strcmp(sTemp.string(),"L")) {

                oReduce.nListResultsL();
                printf("Hit any key to continue...\n");
                getline(cin,sTemp);
                iCell=100;
            } else if (sTemp.string()[0]=='L') {
                // Must load the given lattice character into the cell.
                // Just use the first cell as a dummy.
                iCell=100;
                // Load the lattice number.
                if ((1==sscanf(sTemp.string()+1,"%d",&nx)) && (nx>=1) && (nx<=44)) {
                    iCell = oReduce.nAddSolutionL(nx);
                    if (iCell == -1) {
                        printf("Invalid Lattice!\n");
                    };
                };
            } else if (sTemp.string()[0]=='R') {
                if (1==sscanf(sTemp.string()+1,"%lf",&f0)) {
                    oReduce.m_fResidualMax = f0;
                };
              iCell = 100;
              continue;

            } else if (6==sscanf(sTemp.string(),"%lf %lf %lf %lf %lf %lf",
                    &a6fUserCell[0],&a6fUserCell[1],&a6fUserCell[2],
                    &a6fUserCell[3],&a6fUserCell[4],&a6fUserCell[5])) {

                    iCell=100;

               
                    double  a3x3fOrientMatOut[3][3];           // a3x3fOrientEquivMatOut used elsewhere.
                    double a6fUserCellFound[6];
                    double  a3x3fTrans[3][3];

                    // Try to find the cell, and print out the residual if found.
                    oCrystalOut.vSetCell(a6fUserCell);
                    oCrystalOut.nCalcGetRealMatrix(&a3x3fRealMatOut[0][0]);
                    vTranMat3D(a3x3fRealMatOut);
                    f0=oReduce.fFindOrigVecs(a3x3fRealMatOut, a3x3fRealOrientMatIn, a3x3fRealOrientMatOut);
                    vTranMat3D(a3x3fRealOrientMatOut);            
                    // The orientation matrix is in reciprocal space.
                    fInvMat3D(&a3x3fRealOrientMatOut[0][0],&a3x3fOrientMatOut[0][0]);

                    oCrystalOut.nSetOrientMatrix(&a3x3fOrientMatOut[0][0]);
                    oCrystalOut.vGetOrientAngles(pfOrientAngles);
                    oCrystalOut.vGetCell(a6fUserCellFound);

                    printf("Spec. Cell:  [%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f]\n"
                           "Found Cell:  [%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f]\n"
                           "Orientation: [%4.2f %4.2f %4.2f ]\n"
                           "%% Error:      %4.2f.\n",
                           a6fUserCell[0],a6fUserCell[1],a6fUserCell[2],
                           a6fUserCell[3],a6fUserCell[4],a6fUserCell[5],
                           a6fUserCellFound[0],a6fUserCellFound[1],a6fUserCellFound[2],
                           a6fUserCellFound[3],a6fUserCellFound[4],a6fUserCellFound[5],
                           pfOrientAngles[0],pfOrientAngles[1],pfOrientAngles[2],
                           100.0f*f0
                           );
                    // Print the reindexing matrix.
                    if (oLaue.nCalcReindex(a3x3fOrientCenteredMatIn,a3x3fOrientMatOut,a3x3fTrans)) {
                            printf("Reindexing matrix:  \n"
                                "[h']   [%5.2f %5.2f %5.2f]   [h]\n"
                                "[k'] = [%5.2f %5.2f %5.2f] * [k]\n"
                                "[l']   [%5.2f %5.2f %5.2f]   [l]\n\n",
                                a3x3fTrans[0][0], a3x3fTrans[1][0], a3x3fTrans[2][0],
                                a3x3fTrans[0][1], a3x3fTrans[1][1], a3x3fTrans[2][1],
                                a3x3fTrans[0][2], a3x3fTrans[1][2], a3x3fTrans[2][2]);
                    } else {
                        printf("No reindexing matrix available.\n");
                    };

                    printf("\nUse this cell? Y/[N]: ");
#ifdef SSI_PC
					CCrclHelper::GetInstance()->vBuildCmdCellUse(a6fUserCell, a6fUserCellFound, pfOrientAngles, f0,a3x3fTrans);
#endif
                    getline(cin,sTemp);
                    strupr(sTemp.string());
#ifdef SSI_PC
					cout << sTemp << "\n";
#endif
                    if (!strcmp(sTemp.string(),"Y")) {
                        vCopyVecND(6,a6fUserCell,oLaue.fCellParams);
                        oCrystalOut.vGetOrientAngles(oLaue.fOrientAngles);
                        oLaue.m_nSpaceGroup = 1;
                        
                        if (oLaue.nCalcReindex(a3x3fOrientCenteredMatIn,a3x3fOrientMatOut,a3x3fTrans)) {
                            if (poReflnlistIn) {                               
                                if (oLaue.nReindexError(a3x3fTrans,*poReflnlistIn,bPrompt)) {
                                    poReflnlistIn->nReindex(&a3x3fTrans[0][0]);               
                                    bReindexingDone=TRUE;
                                } else {
                                    printf("Terminating.\n");
                                    goto main_termination;
                                };
                            };
                        } else {
                            printf("No reindex matrix available.\nNo re-indexing Done.\n");
                        };
                        goto end_of_reduce;
                    };                    
                    continue;

            } else if (oReduce.nParseSelection(sTemp,iCell))
              continue;
           
            // Allow re-listing if input > oReduce.oGetSolutions().nGetNumReflns().
            if (oReduce.oGetSolutions().nGetNumReflns() < iCell) {
                nRelist = iCell;
                oLaue.fResidMax = (double) nRelist;
                oReduce.m_fResidualMax = oLaue.fResidMax;
            };
        } while (iCell>oReduce.oGetSolutions().nGetNumReflns());

        if (iCell<=0) {
            printf("Terminating.\n");
            goto main_termination;
        };

        // Load the new cell.
        oReduce.nLoadSolution(iCell-1, &oCrystalOut);

        // Get the real orientation matrix of the new cell.
        oCrystalOut.nCalcGetRealMatrix(&a3x3fRealMatOut[0][0]);
        vTranMat3D(a3x3fRealMatOut);
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Find the real orientation matrix.
        oReduce.fFindOrigVecsA(a3x3fRealMatOut, a3x3fRealOrientMatIn, a3x3fRealOrientMatOut);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Allow choices of orientation matrix and reindexing matrix.
        double  a3x3fEquivRealOrientMatOut[24][3][3];           // These are returned by the nEquivRotations routine.
        double  a3x3fEquivOrientMatOut[24][3][3];               // We convert them to a reciprocal space orientation matrix.
        double  a3x3fEquivTrans[24][3][3];                      // Then, we find the reindexing transformation 
        int     nNumEquivOrientations = 0;
        int     nEquivCt = 0;
        int     nBestSelection = -1;                            // The selection that has angles closest to the input cell.
        double  fBestSelection = 1.0e20;                        // Sum of the absolute values of angle differences.
         
        nNumEquivOrientations = oReduce.nEquivRotations(a3x3fRealOrientMatOut,
                                                        a3x3fEquivRealOrientMatOut,
                                                        oReduce.m_a44tLattChar[oReduce.m_nLattNum].sLatticeType);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        for(nEquivCt=0; nEquivCt < nNumEquivOrientations; nEquivCt++)
        {
            vTranMat3D(a3x3fEquivRealOrientMatOut[nEquivCt]);
            
            // The orientation matrix is in reciprocal space.
            fInvMat3D(&a3x3fEquivRealOrientMatOut[nEquivCt][0][0], &a3x3fEquivOrientMatOut[nEquivCt][0][0]);
                       
            oLaue.nCalcReindex(a3x3fOrientCenteredMatIn, a3x3fEquivOrientMatOut[nEquivCt], a3x3fEquivTrans[nEquivCt]);
        }
        
        do
        {
            Ccrystal oCrystal;
#ifdef SSI_PC
			CCrclHelper::GetInstance()->vBuildCmdCellSetOrients( nNumEquivOrientations, a3x3fEquivTrans, a3x3fEquivOrientMatOut);
#endif
            for (nEquivCt=0; nEquivCt <nNumEquivOrientations; nEquivCt ++)
            {
                for (nx=0;nx<3;nx++)
                {
                    char* cphkl="hkl";
                    char cpNumField[5];
                    
                    if (nx==0)
                        sprintf(cpNumField,"#%-2d",nEquivCt+1);
                    else
                        strcpy(cpNumField,"   ");
                    
                    printf("%s [%c']%c[ ",cpNumField,cphkl[nx],(nx==1)?('='):(' '));
                    for (ny=0;ny<3;ny++) 
                        printf("%5.2f ",a3x3fEquivTrans[nEquivCt][ny][nx]);
                    printf("]%c[%c]",(nx==1)?('*'):(' '), cphkl[nx]);                        
                    if (nx==1) {
                        // Load a crystal, and extract the orientation matrix.
                        oCrystal.nSetOrientMatrix(&a3x3fEquivOrientMatOut[nEquivCt][0][0]);
                        oCrystal.vGetOrientAngles(pfOrientAngles);
                        
                        printf("  Orientation = [ %4.2f %4.2f %4.2f ]\n",
                               pfOrientAngles[0],pfOrientAngles[1],pfOrientAngles[2]);
                    }
                    else 
                        printf("\n");
                }
                printf("\n");

                f0=fabs(pfOrientAngles[0]-oLaue.fOrientAngles[0])+
                    fabs(pfOrientAngles[1]-oLaue.fOrientAngles[1])+
                    fabs(pfOrientAngles[2]-oLaue.fOrientAngles[2]);
                if (f0<fBestSelection) {
                    fBestSelection=f0;
                    nBestSelection=nEquivCt;
                }
            }


            printf("Enter rotation selection (0=Abort Default=Choose #%d): ",nBestSelection+1);
            if (bPrompt)
			{
#ifdef SSI_PC
				if( bAutoMode )
				{
					sTemp.Format("%d", nBestSelection+1);
				}
				else
				{
					getline(cin,sTemp);
				}
#else
				getline(cin,sTemp);
#endif
			}
            else
                sTemp = "";
#ifdef SSI_PC
		    if (0 < sTemp.length())
			    cout << sTemp << "\n";
		    else
			    cout << "<Enter>" << "\n";
#else
            printf("\n");
#endif
            if (sTemp.string()[0]==0)
            {
                iOrient=nBestSelection+1;
            }
            else
            {
                sscanf(sTemp.string(),"%d",&iOrient);
            }
            
            iOrient--;
            
            if (iOrient==-1)
            {
                goto restart_reduce;
            }
        } while ((iOrient<0) || (iOrient>=nNumEquivOrientations));

        // Clear the lattice so that it dosen't cause problems later.
        oLaue.a2cLattice[0] = 0;

        if (poReflnlistIn) {
            if (oLaue.nReindexError(a3x3fEquivTrans[iOrient],*poReflnlistIn,bPrompt)) {
                poReflnlistIn->nReindex(&a3x3fEquivTrans[iOrient][0][0]);
                bReindexingDone=TRUE;
            } else {
                printf("Terminating.\n");
                goto main_termination;
            };
        };
 
        // We now have the oriented matrix in a3x3fEquivOrientMatOut.  Load the cell with this UB matrix.
        oCrystalOut.nSetOrientMatrix(&a3x3fEquivOrientMatOut[iOrient][0][0]);
               
        // Get the orientation angles of this cell.
        oCrystalOut.vGetOrientAngles(&oLaue.fOrientAngles[0],&oLaue.fOrientAngles[1],&oLaue.fOrientAngles[2]);
        oLaue.m_nSpaceGroup = oCrystalOut.m_poSpacegroup->nGet();

        // Reload the solution.  We want to print out the ideal cell, with the non-ideal orientation angles.
        oReduce.nLoadSolution(iCell-1, &oCrystalOut);
        oCrystalOut.vGetCell(oLaue.fCellParams);
        printf("Crystal Cell:        [%4.2f %4.2f %4.2f  %4.2f %4.2f %4.2f]\n",
            oLaue.fCellParams[0],oLaue.fCellParams[1],oLaue.fCellParams[2],
            oLaue.fCellParams[3],oLaue.fCellParams[4],oLaue.fCellParams[5]
            );
        printf("Crystal Orientation: [%4.2f %4.2f %4.2f]\n",oLaue.fOrientAngles[0],oLaue.fOrientAngles[1],oLaue.fOrientAngles[2]);
        
    };
    end_of_reduce:

	if (bAddTwinLaw) {
        CCellReduce  oReduce;
		Ccrystal oCrystalIn(*poHeader);
		if (oReduce.nAddCreateTwinLaw(oCrystalIn,sTwinLaw,fTwinLawRot)) {
			goto main_termination;
		};
		oCrystalIn.vSetTwin(nTwinID);
		oCrystalIn.nUpdateHeader(poHeader);
		poHeader->nWrite(sOutputHeader);
		goto main_termination;

	};
    
    // If we are transforming, then that is all we are doing.  We can just transform and exit.
    if (bDoTransform) {
        double fReindexMat[3][3];            

        if (nTransformIndex==-1) {
            printf("Reindexing using user specified transform.\n");
            vCopyMat3D(&a3x3fUserTrans[0][0],&fReindexMat[0][0]);
        } else if (nTransformIndex==-2) {
            if (oLaue.m_nSpaceGroup<=0) {
                printf("Spacegroup must be known to use -transform sRefFile\n"
                    "Please provide header or use -spacegroup option.\n"
                    );
                goto main_termination;
            };                
            printf("Reindexing using reference file.\n");
            if (oLaue.nReindexRelative(*poReflnlistIn,*poRefTransform,fReindexMat,fTransformScale)) {
                printf("Terminating.\n");
                goto main_termination;
            };
        } else {
            printf("Reindexing using transform %s\n",g_cpTransforms[nTransformIndex]);
            vCopyMat3D(g_fReindexTransMats[nTransformIndex],&fReindexMat[0][0]);
        }

        if (!poReflnlistIn) {
            printf("Reflection list must be specified to use -transform sRefFile\n");
            goto main_termination;
        };
        

        printf("Reindexing matrix:  \n"
            "[h']   [%5.2f %5.2f %5.2f]   [h]\n"
            "[k'] = [%5.2f %5.2f %5.2f] * [k]\n"
            "[l']   [%5.2f %5.2f %5.2f]   [l]\n\n",
            fReindexMat[0][0], fReindexMat[1][0], fReindexMat[2][0],
            fReindexMat[0][1], fReindexMat[1][1], fReindexMat[2][1],
            fReindexMat[0][2], fReindexMat[1][2], fReindexMat[2][2]);      
        if (!oLaue.nReindexError(fReindexMat,*poReflnlistIn,bPrompt)) {
            printf("Terminating.\n");
            goto main_termination;
        };
        if (fTransformScale != 1.0) {
            printf("Applying scale factor ...\n");
            for (ny=0;ny<poReflnlistIn->nGetNumReflns();ny++) {
                (*poReflnlistIn)[ny].vSetIntensity((*poReflnlistIn)[ny].fGetIntensity()*(float)fTransformScale);
                (*poReflnlistIn)[ny].vSetSigmaI((*poReflnlistIn)[ny].fGetSigmaI()*(float)fTransformScale);
            };
        };


        bReindexingDone=TRUE;
        poReflnlistIn->nReindex(&fReindexMat[0][0]);

        oLaue.nCalcCell(fReindexMat);
        if ((bPrompt) && (bReindexingDone)) {
            printf("Write \"%s\" to \"%s\" ? [Y]/N  ",
                sRefName.string(),sOutputRef.string());
#ifdef SSI_PC
			CCrclHelper::GetInstance()->vBuildCmdCellWriteRefFile();
#endif
            getline(cin,sTemp);
#ifdef SSI_PC
			if (0 < sTemp.length()) {
				cout << sTemp << "\n";
			}
			else {
				cout << "<Enter>" << "\n";
			}
#endif
            if (!((sTemp.length()==0) || (sTemp=="Y") || (sTemp=="y")))
            {
                printf("Reflections NOT updated.\n");
                if (!poHeader)
                    poHeader = new Cimage_header;
                oLaue.nWriteHeader(&sOutputHeader,*poHeader);
#ifdef SSI_PC                
                vNotifyClientAboutFileWritten(sOutputHeader, DTREK_CCHI_BCCOF_FILE_HDR);
#endif
                goto main_termination;
            }
        }


        if( bReindexingDone )
        {
            printf("Writing out reflection file \"%s\" ...\n",sOutputRef.string());
            if (!poHeader)
                poHeader = new Cimage_header;
            oLaue.nWriteHeader(NULL,*poHeader);
            poReflnlistIn->nPutCrystal(*poHeader);
            poReflnlistIn->nWrite(sOutputRef);
#ifdef SSI_PC
			vNotifyClientAboutFileWritten(sOutputRef, DTREK_CCHI_BCCOF_FILE_REF);
#endif
        }
        
        if( !poHeader )
            poHeader = new Cimage_header;
        
        oLaue.nWriteHeader(&sOutputHeader,*poHeader);

#ifdef SSI_PC
	    vNotifyClientAboutFileWritten(sOutputHeader, DTREK_CCHI_BCCOF_FILE_HDR);
#endif

        goto main_termination;
    }

    if( poReflnlistIn )
    {
        ////////////////////////////////////////////////////////////////////////////////
        //  DO LAUE CHECK 
        /////////////////////////////////////////////////////////////////////////////////
        if( m_bCheckLaue )
        {
            oLaue.nLaueCheck(*poReflnlistIn, bTemp, bPrompt);
            bReindexingDone = bTemp || bReindexingDone;
			
            if( oLaue.m_eLaue == LAUE_UNKNOWN )
            {
				printf("Could not determine Laue group.\n");
				goto main_termination;
			}
		}
        
        ////////////////////////////////////////////////////////////////////////////////
        //  DO CENTRICITY CHECK 
        /////////////////////////////////////////////////////////////////////////////////
        if( m_bCheckCentric )
        {
            oLaue.nCentricTest(*poReflnlistIn);
        }
        
        ////////////////////////////////////////////////////////////////////////////////
        //  DO SPACEGROUP CHECK 
        /////////////////////////////////////////////////////////////////////////////////
        if( m_bSpacegroupCheck && oLaue.m_eLaue != LAUE_UNKNOWN ) 
        {
            if( poHeader )
                poHeader->nClearDictionary("DTCELL_INCONCLUSIVE_HKL_CONDITIONS", "*");

            printf("Checking systematic absences...\n");
           
            nx = oLaue.nChooseSpaceGroup(*poReflnlistIn, bTemp, bPrompt);

            if( nx != 0 )
            {
                printf("Terminating.\n");
                goto main_termination;
            }
            
            bReindexingDone = bTemp || bReindexingDone;
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( bReindexingDone )
        {
            bool    bConfirmedWritingReindexedReflnList = true;
        
            if( bPrompt )
            {
			    // This print statement required by StructureStudio.
			    printf("\nReindexing was done.  The reindexed information should be written to %s.\n", sOutputRef.string());
                printf("\nWrite reindexed information to \"%s\"  [Y]/N ", sOutputRef.string());
#ifdef SSI_PC
		        CCrclHelper::GetInstance()->vBuildCmdCellWriteRefFile();
#else
		        fflush(stdout);
#endif
                getline(cin,sTemp);
#ifdef SSI_PC
                printf("\n%s\n", sTemp.string());
#endif
                bConfirmedWritingReindexedReflnList = sTemp.length() == 0 || sTemp=="Y" || sTemp=="y";
            } 
            
            if( bConfirmedWritingReindexedReflnList )
            {
                printf("\nWriting out reflection file \"%s\" ...\n\n",sOutputRef.string());
                
                if (!poHeader)
                    poHeader = new Cimage_header;
                
                oLaue.nWriteHeader(NULL,*poHeader);
                poReflnlistIn->nPutCrystal(*poHeader);
                poReflnlistIn->nWrite(sOutputRef);
#ifdef SSI_PC
			    vNotifyClientAboutFileWritten(sOutputRef, DTREK_CCHI_BCCOF_FILE_REF);
#endif
            }
            else
            {
#ifdef SSI_PC
                printf("\nN\n");
#endif
                printf("\nReflection file not updated.\n");
            }
        }
        else
        {
            printf("\nNo reindexing done.\nReflection file not updated.\n\n");
        }
    }
    
    // Always write out a header.
    if( !poHeader )
        poHeader = new Cimage_header();
    
    oLaue.nWriteHeader(&sOutputHeader,*poHeader);

#ifdef SSI_PC
	// Let CC know we have written the header file
    vNotifyClientAboutFileWritten(sOutputHeader, DTREK_CCHI_BCCOF_FILE_HDR);
#endif
    
main_termination:
	if( poReflnlistIn )
        delete poReflnlistIn;
	
    if( poHeader )
        delete poHeader;
    
    if( poRefTransform )
        delete poRefTransform;
    
    GETCH;

#ifdef WIN32
    #ifndef SSI_PC
        #ifdef _DEBUG
            _CrtDumpMemoryLeaks();
        #endif
    #endif
#endif
    
    printf ("\ndtcell: Done.\n");
    return 0;
}

#ifdef SSI_PC
// Let CrystalClear know we have written an output file. 
// sText - file name; wCtrl - file type.
// Some CrystalClear dialogs do not need the message, so do not send it.
// Some CrystalClear dialogs do need the message, plus they need d*TREK to wait for a confirmation,
// so CrystalClear could do some maintenance while d*TREK is waiting.
void CDTMainCell::vNotifyClientAboutFileWritten(Cstring& sText, DTREK_WORD wCtrl)
{
    bool    bConfirm = false;
    if(    (m_bReduce && !m_bCheckLaue && !m_bSpacegroupCheck && !m_bCheckCentric)    // Reduce Cell
        || (m_bCheckLaue && !m_bSpacegroupCheck && !m_bCheckCentric)                  // Laue Check
        || (!m_bCheckLaue && m_bSpacegroupCheck && !m_bCheckCentric)                  // Space Group Check
        || (!m_bCheckLaue && !m_bSpacegroupCheck && m_bCheckCentric)                  // Centricity Check
      )
      bConfirm = true; // In these cases CrystalClear will confirm it knows the header file has been written    
    
    bool    bRet = CCrclHelper::GetInstance()->bBuildCmdCellOutputFile(sText.string(), wCtrl);
    if( !bRet )
        return;
    
    if( !bConfirm )
    {
        CCrclHelper::GetInstance()->SendTclCmd();
        return;
    }

    // Send message and wait for a confirmation from client
    Cstring     sTemp("");
    while( true )
    {
        getline(cin,sTemp);
        
        if( DTREK_CCHI_BCCOF_FILE_REF == wCtrl )
        {
            if( sTemp == "REF OK" )
                break;
        }
        else if( DTREK_CCHI_BCCOF_FILE_HDR == wCtrl )
        {
            if( sTemp == "HDR OK" )
                break;
        }
    }
}
#endif

void CDTMainCell::vError(const int nErrorNum, const Cstring& sMessage)
{
}


