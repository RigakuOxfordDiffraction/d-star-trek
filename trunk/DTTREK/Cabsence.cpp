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
// Cindex.cc            Initial author: Thaddeus J Niemeyer     11-Nov-1999
//  This file contains the member functions of class Cindex which implements
//    the autoindexing encapsulation of d*TREK.
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

#include "Cstring.h"
#include "Dtrek.h"
#include "Cimage.h"
#include "dtrekvec.h"
#include "dtsvd.h"
#include "Crefln.h"
#include "Creflnlist.h"

#define CCELL_DEFINE_DATA
#include "dtcell.h"

#include "Cspacegroup_check.h"

#include <string.h>
#include <stdio.h>


#ifdef SSI_PC
#include "CrclHelper.h"
#endif

#if !defined(VC6)
using std::cin;
using std::cout;
using std::endl;
using std::flush;
#endif

#if 0
#define ERROR {  printf("Bad Logic Line %d of file '%s'\n", __LINE__ , __FILE__ ); GETCH; exit(0); };
//template <class T> void swap(T& x,T& y) { T z; z=x; x=y; y=z; };
#define STR(x) #x
int g_anCentricSpacegroups[] = { 2,2,  10,15,  47,74,  83,88,  123,142,  123,142,  147,148,  164,167,  162,163,  175,176,  191,194,  200,206,  221,230,  -1 };
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Claue::nChooseSpaceGroup(Creflnlist& oRef,
                             bool& bReIndexingDone,
                             bool bPrompt)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set up control structure for space group check
    SPACEGROUP_CHECK_CTRL_INFO       stInfo;
   
    stInfo.eLaue = m_eLaue;
    
    stInfo.nReflnCountTol = 5;
    stInfo.nReflnCountTol_IOverSigma = 3;
    stInfo.fIOverSigmaTol = fAbsentSigma;

    stInfo.nRejects_MaxSaveCount = nMaxCollectRejects;
    stInfo.bCollectRejects = bCollectRejects;

    stInfo.eCentricity = m_bIsCentric ? 
                         SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric : 
                         SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Acentric;

    stInfo.eChirality  = m_bIsChiral  ? 
                         SPACEGROUP_CHECK_CTRL_INFO::enChirality_Chiral  : 
                         SPACEGROUP_CHECK_CTRL_INFO::enChirality_Unknown;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Instantiate the spacegroup check class
    Cspacegroup_check   oCheck(stInfo);
    
    // Run the space group check
    oCheck.bTestReflnList(oRef);
    oCheck.vTestSpaceGroupTableEntries();
//mrp    oCheck.vSelectCandidateSpaceGroups();
	oCheck.vSelectCandidateSpaceGroups(m_nSpaceGroup); //mrp

    m_strInconclusiveConditions = oCheck.sGetListInconclusiveConditionsForCandidateSpaceGroups();

    // Fill up Texsan output (parities) classes
    for (int ii = 0; ii < TEXSANCLASSES; ii++) 
    {
        m_nTexsanTotal[ii]    = oCheck.m_stTexsanStats.m_nTexsanTotal[ii];
        m_nTexsanObsvd[ii]    = oCheck.m_stTexsanStats.m_nTexsanObsvd[ii];
        m_fTexsanIoverSig[ii] = oCheck.m_stTexsanStats.m_fTexsanIoverSig[ii];
    }    
    ///////////////////////////////////////////////////////////////////////////////////

    // Printing flags
#ifdef SSI_PC
    // John wants all tables printed out in one pass.
    bool        bPrintAllAbsences = true;
    bool        bPrintSelectedAbsences = true;
    bool        bPrintTexsan = true;
    bool        bPrintSelectedSpacegroups = true;
#else
    bool        bPrintAllAbsences = false;
    bool        bPrintSelectedAbsences = true;
    bool        bPrintTexsan = false;
    bool        bPrintSelectedSpacegroups = true;
#endif

    bReIndexingDone = false;

    // Figure out the spacegroup.  We might need to prompt.
    if( oCheck.nGetSelectedCandidateSpaceGroupIndex() < 0 )
    {
        printf("No spacegroups available for Laue group %s\n"
               "(given centricity, chirality and absences).\n"
               "Consider increasing max allowed Rmerge with -maxrmerge.\n", g_cpLaueName[m_eLaue]);

        bPrintAllAbsences = true;
    }
    
    Cstring     sOut("");                                   //  Data sent to CrystalClear.
    Cstring     sTemp("");
    
    int         nUserSpaceGroup = -1;                       // Space group number (International Tables) selected by the user
    int         iSpaceGroupSel = -1;                        // Index of the selected space group among candidates
    int         nSelected_SpaceGroupIntTablesNumber = -1;   // Space group number (International Tables) of the selected space group (among candidates)   

    SPACEGROUP_CHECK_CTRL_INFO::eCentricityType     eSelectedSpaceGroupCentricity = SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Unknown;
    bool        bIsUserCentricity = false;
    
    bool        bUserSpaceGroupEntered = false;

    do
    {
        bUserSpaceGroupEntered = false;

        if( bPrintTexsan )                
            nPrintTexsan();

        if( bPrintAllAbsences )
        {
            sOut = "";
            oCheck.vPrintAllConditionsRestrictions(sOut);
#ifdef SSI_PC
            if( sOut != "" )
            {
                CCrclHelper::GetInstance()->vBuildCmdCellAllAbsenceTable(sOut);
                printf("Sending tcl command\n");
                CCrclHelper::GetInstance()->SendTclCmdRespond();
            }
#endif                         
        }
                      
        // Print out the restrictions, and absences for the selected line of the table.
        if( bPrintSelectedAbsences ) 
        {
            sOut = "";
            oCheck.vPrintSystematicAbsencesForSelectedSpaceGroup(sOut);
#ifdef SSI_PC
            if( sOut != "" )
            {
                CCrclHelper::GetInstance()->vBuildCmdCellAbsenceTable(sOut);
                printf("Sending tcl command\n");
                CCrclHelper::GetInstance()->SendTclCmdRespond();
            }
#endif              
        }

        if( bPrintSelectedSpacegroups )
        {
            sOut = "";
            oCheck.vPrintSelectedSpaceGroups(sOut);

#ifdef SSI_PC
            if( sOut != "" )
            {
                CCrclHelper::GetInstance()->vBuildCmdCellSpacegroupTable( sOut );
                CCrclHelper::GetInstance()->SendTclCmdRespond();
            }
#endif
        }
        
        nSelected_SpaceGroupIntTablesNumber = oCheck.nGetSelectedCandidateSpaceGroupIntTablesNumber();    
        eSelectedSpaceGroupCentricity = oCheck.eGetSelectedCandidateSpaceGroupCentricity();
        bIsUserCentricity = (  m_bIsCentric && eSelectedSpaceGroupCentricity == SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric) ||
                            ( !m_bIsCentric && eSelectedSpaceGroupCentricity == SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Acentric);
#ifdef SSI_PC
        bool bWarn = (nSelected_SpaceGroupIntTablesNumber == -1);
#endif
        if( -1 == nSelected_SpaceGroupIntTablesNumber )
        {
            printf("WARNING:  Did not find any space groups that were %s\n", m_bIsCentric ? "centric" : "acentric");
        }
        else if( !bIsUserCentricity && bPrintSelectedSpacegroups ) 
        {
            printf("WARNING:  Could not find %s spacegroups.\n", !m_bIsCentric ? "acentric" : "centric");
        }
        
        if( m_strInconclusiveConditions != "" )
        {
            printf("\nWARNING: Spacegroup selection is inconclusive, because the following\n"
                     "         reflection conditions have too few reflections in the input file:\n\n"
                     "         %s\n\n",  m_strInconclusiveConditions.string());
        }

        //////////////////////////////////////////////////////////////////////////////////////////
        printf("\nOptions:\n"
               "T)      Print reflection parities (Texsan output)\n"
               "S)      Print found spacegroups\n"
               "A)      Print Table of absences applicable to found spacegroups\n"
               "B)      Print Table of all absences applicable to Laue class\n"
               "Q)      Quit\n"
               "#)      Type spacegroup number\n");

        if( -1 != nSelected_SpaceGroupIntTablesNumber )
            printf("Enter)  Select spacegroup %d\n\n", nSelected_SpaceGroupIntTablesNumber);
        
        printf("Choice: ");
        ////////////////////////////////////////////////////////////////////////////////////////////

#ifdef SSI_PC
        CCrclHelper::GetInstance()->vBuildCmdCellSpacegroupSelected(nSelected_SpaceGroupIntTablesNumber,
                                                                    m_strInconclusiveConditions,
                                                                    bWarn);
#else
		fflush(stdout);
#endif

        if( !bPrompt )
            sTemp = "";
        else
            getline(cin, sTemp);
        
#ifdef SSI_PC
        sTemp.upcase();
#endif

        printf("%s\n", sTemp.string());

        bPrintTexsan = false;
        bPrintAllAbsences = false;
        bPrintSelectedAbsences = false;
        bPrintSelectedSpacegroups = false;
        
        if (sTemp == "T") 
        {
            bPrintTexsan = true;
            continue;
        }
        else if (sTemp=="S")
        {
            bPrintSelectedSpacegroups=true;
            continue;
        }
        else if (sTemp=="A")
        {
            bPrintSelectedAbsences=true;
        }
        else if (sTemp=="B")
        {
            bPrintAllAbsences = true;
        } 
        else if (sTemp=="Q")
        {
            m_nSpaceGroup = -1;
            m_strInconclusiveConditions = "";
            return 1;
        } 
        else if (sTemp.string()[0]==0)
        {
            nUserSpaceGroup = nSelected_SpaceGroupIntTablesNumber;
            
            bUserSpaceGroupEntered = true;
        }
        else
        {
            if( 1!=sscanf(sTemp.string(),"%d", &nUserSpaceGroup) ) 
                continue;
            
            bUserSpaceGroupEntered = true;
        }
        
        if( nUserSpaceGroup >= 1 && nUserSpaceGroup<= 230 )
        {
            if( -1 == (iSpaceGroupSel = oCheck.nGetCandidateSpaceGroupIndexBySGNumber(nUserSpaceGroup)) )
            {
                m_nSpaceGroup = nUserSpaceGroup; 
                
                printf("Spacegroup #%d not found in list of candidates... no reindexing information\n"
                       "available.\n", nUserSpaceGroup);
                
                return 0;
            }
            else // if user selection is found among the space groups candidates
            {
                m_nSpaceGroup = nUserSpaceGroup;

                if( SPACEGROUP_CHECK_CTRL_INFO::enCentricity_Centric == oCheck.eGetCandidateSpaceGroupCentricity(iSpaceGroupSel) )
                    m_bIsCentric = true;
                else
                    m_bIsCentric = false;
            }
        }
        else if( bUserSpaceGroupEntered )
        {
            printf("\n\nERROR: spacegroup number must be between 1 and 230.\n\n");
        }
    } while( iSpaceGroupSel == -1 );
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    
    oCheck.vSetSelectedCandidateSpaceGroupIndexBySGNumber(m_nSpaceGroup);

    // We need to determine whether a transformation should be carried out.
    // This will only be done if we are in orthorhombic or monoclinic spacegroups
    int         nSelSGNumberIntTables = oCheck.nGetSelectedCandidateSpaceGroupIntTablesNumber();
    Cstring     sSelSGPresentation = oCheck.sGetSelectedCandidateSpaceGroupPresentation();

    if( !(sSelSGPresentation == g_cpSpaceGroupNames[nSelSGNumberIntTables-1]) )
    {
        eIndexTransform     eITrans = oCheck.eFindTranslationPermForSelectedSpaceGroup();

#ifdef SSI_PC
		CCrclHelper::GetInstance()->vBuildCmdCellSpacegroupTransform(g_cpTransforms[eITrans]);
#endif

        printf("\nTransformation %s is required to move to standard presentation.\n"
               "Do you want to apply this transformation now? [Y]/N ", g_cpTransforms[eITrans]);
        
        if( bPrompt )
        {
            getline(cin, sTemp);
#ifdef SSI_PC
            printf("\n%s\n", sTemp.string());
#endif
        }

        if( sTemp.length() == 0 || sTemp == "Y" || sTemp=="y" || !bPrompt )
        {
            printf("\n\nReindexing reflections...\n\n");
            
            oRef.nReindex(g_fReindexTransMats[eITrans]);
            
            // Must also calculate the new cell induced.
            double      fReindexMat[3][3];
            
            vCopyMat3D(g_fReindexTransMats[eITrans], &fReindexMat[0][0]);
            
            nCalcCell(fReindexMat);
            
            bReIndexingDone = true;
        }
    }
    else
    {
#ifdef SSI_PC
		// Let CrystalClear know that no reindexing is required
        CCrclHelper::GetInstance()->vBuildCmdCellSpacegroupTransform("I");
        CCrclHelper::GetInstance()->SendTclCmd();
#endif
    }
    
    return 0;
}

