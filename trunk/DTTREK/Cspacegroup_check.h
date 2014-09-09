//
// Copyright (c) 2007 Rigaku Americas Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// Cspacegroup_check.h   Initial author: RB           01-Aug-2007

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
#ifndef DT_CSPACEGROUPCHECK_H
#define DT_CSPACEGROUPCHECK_H

#include "Cstring.h"
#include "Creflnlist.h"
#include "dtcell.h"

struct SPACEGROUP_CHECK_CTRL_INFO
{
    eLaueType eLaue;
    
    int     nReflnCountTol;
    int     nReflnCountTol_IOverSigma;
    double  fIOverSigmaTol;
    
    bool    bCollectRejects;
    int     nRejects_MaxSaveCount;

    enum eCentricityType
    {
        enCentricity_Unknown=-1,
        enCentricity_Centric,
        enCentricity_Acentric
    }     eCentricity;

    enum eChiralityType
    {
        enChirality_Unknown=-1,
        enChirality_Chiral,
    }     eChirality;

    void vInit()
    {
        eLaue = LAUE_UNKNOWN;

        nReflnCountTol = 0;
        nReflnCountTol_IOverSigma = 0;
        fIOverSigmaTol = 0.0;
        
        bCollectRejects = false;
        nRejects_MaxSaveCount = 0;

        eCentricity = enCentricity_Unknown;
        eChirality = enChirality_Unknown;
    }
    
    SPACEGROUP_CHECK_CTRL_INFO& operator=(const SPACEGROUP_CHECK_CTRL_INFO& another)
    {
        eLaue                        = another.eLaue;

        nReflnCountTol               = another.nReflnCountTol;                      
        nReflnCountTol_IOverSigma    = another.nReflnCountTol_IOverSigma;
        fIOverSigmaTol               = another.fIOverSigmaTol;
                                          
        bCollectRejects              = another.bCollectRejects;
        nRejects_MaxSaveCount        = another.nRejects_MaxSaveCount;

        eCentricity                  = another.eCentricity;
        eChirality                   = another.eChirality;
        
        return *this;
   }
};

#define MAX_NUMBER_HKL_CONDITIONS_IN_TABLE      7

typedef enum enCond
{ 
    enCond_HKL=0, 
    enCond_0KL, 
    enCond_H0L,
    enCond_HK0, 
    enCond_H00,
    enCond_0K0, 
    enCond_00L, 
    enCond_HHL, 
    enCond_HH0, 
    enCond_HKIL, 
    enCond_H_H0L, 
    enCond_HH_2_HL, 
    enCond_000L,
    enCond_Unknown
};

typedef enum enRestr
{ 
    enRestr_H,
    enRestr_K,
    enRestr_L,
    enRestr_KL,
    enRestr_HL,
    enRestr_HK,
    enRestr_4H,
    enRestr_4K,
    enRestr_4L,
    enRestr_4KL,
    enRestr_4HL,
    enRestr_4HK,
    enRestr_HKL,
    enRestr_W,
    enRestr_X,
    enRestr_Y,
    enRestr_Z,
    enRestr_3H,
    enRestr_3L,
    enRestr_6L,
    enRestr_3HL,
    enRestr_Unknown
};

#define  DTREK_SGC_REFLN_RESTR_UNKNOWN    0x00000000

#define  DTREK_SGC_REFLN_RESTR_H    0x00000001
#define  DTREK_SGC_REFLN_RESTR_K    0x00000002
#define  DTREK_SGC_REFLN_RESTR_L    0x00000004

#define  DTREK_SGC_REFLN_RESTR_HK   (DTREK_SGC_REFLN_RESTR_H | DTREK_SGC_REFLN_RESTR_K)
#define  DTREK_SGC_REFLN_RESTR_KL   (DTREK_SGC_REFLN_RESTR_K | DTREK_SGC_REFLN_RESTR_L)
#define  DTREK_SGC_REFLN_RESTR_HL   (DTREK_SGC_REFLN_RESTR_H | DTREK_SGC_REFLN_RESTR_L)

#define  DTREK_SGC_REFLN_RESTR_4    0x00000040

#define  DTREK_SGC_REFLN_RESTR_4H  (DTREK_SGC_REFLN_RESTR_4 | DTREK_SGC_REFLN_RESTR_H)
#define  DTREK_SGC_REFLN_RESTR_4K  (DTREK_SGC_REFLN_RESTR_4 | DTREK_SGC_REFLN_RESTR_K)
#define  DTREK_SGC_REFLN_RESTR_4L  (DTREK_SGC_REFLN_RESTR_4 | DTREK_SGC_REFLN_RESTR_L)

#define  DTREK_SGC_REFLN_RESTR_4HK (DTREK_SGC_REFLN_RESTR_4 | DTREK_SGC_REFLN_RESTR_HK)
#define  DTREK_SGC_REFLN_RESTR_4KL (DTREK_SGC_REFLN_RESTR_4 | DTREK_SGC_REFLN_RESTR_KL)
#define  DTREK_SGC_REFLN_RESTR_4HL (DTREK_SGC_REFLN_RESTR_4 | DTREK_SGC_REFLN_RESTR_HL)

#define  DTREK_SGC_REFLN_RESTR_HKL  0x00001000

// W:    -nH      + nL = 3n
#define  DTREK_SGC_REFLN_RESTR_W    0x00002000

// X: 2 * nH      + nL = 4n
#define  DTREK_SGC_REFLN_RESTR_X    0x00004000

// Y:    -nH + nK + nL = 3n
#define  DTREK_SGC_REFLN_RESTR_Y    0x00008000

// Z:     nH - nK + nL = 3n
#define  DTREK_SGC_REFLN_RESTR_Z    0x00010000

#define  DTREK_SGC_REFLN_RESTR_3H   0x00020000
#define  DTREK_SGC_REFLN_RESTR_3L   0x00040000

#define  DTREK_SGC_REFLN_RESTR_6    0x00080000
#define  DTREK_SGC_REFLN_RESTR_6L  (DTREK_SGC_REFLN_RESTR_6 | DTREK_SGC_REFLN_RESTR_3L)

#define  DTREK_SGC_REFLN_RESTR_3HL  0x00100000

///////////////////////////////////////////////////////////////////////////////////////
struct TEXSAN_STATS
{
    int m_nTexsanTotal[TEXSANCLASSES];
    int m_nTexsanObsvd[TEXSANCLASSES];
    double m_fTexsanIoverSig[TEXSANCLASSES];
    
    TEXSAN_STATS(){vInit();}

    void vInit()
    {
        for(int ii=0; ii < TEXSANCLASSES; ii++)
        {
            m_nTexsanTotal[ii] = 0;
            m_nTexsanObsvd[ii] = 0;
            m_fTexsanIoverSig[ii] = 0.0;
        }
    }

    void vAdd(int nIndex, bool bAbsent, double fIoverSig)
    {
        m_nTexsanTotal[nIndex]++; 
        
        m_nTexsanObsvd[nIndex] += bAbsent ? 0  : 1;  
        
        m_fTexsanIoverSig[nIndex] += fIoverSig;
    }
    
    void vTally()
    {
        for(int ii=0; ii < TEXSANCLASSES; ii++) 
        {
            if( m_nTexsanTotal[ii] > 0 )
                m_fTexsanIoverSig[ii] /= m_nTexsanTotal[ii];
        }
    }
};
///////////////////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CSpaceGroupTableInfo
{
public:
    CSpaceGroupTableInfo();

public:
    DTREK_DWORD dwGetRestrictionCode(Cstring& sName);
    enRestr  eGetRestrictionEnum(DTREK_DWORD dwCode);
    int  nGetRestrictionCount(){return (int)m_mapRestrictionEnum.size();}
    int nGetRestrictionParity(DTREK_DWORD dwCode);
    Cstring& rsGetTableByLaueType(eLaueType eLaue);
    enCond eGetConditionValue(Cstring& sName);
    int nGetConditionIndex(enCond eID);
    std::map<Cstring, DTREK_DWORD>& rmapGetRestrictionNamesIDs(){return *(&m_mapRestrictions);}
    eLaueType eGetLaueClassEnumFromName(Cstring& sName);
    Cstring sGetLaueName(eLaueType eType);
    int nGetTotalNumberOfFullListRestrictions(){return (int)m_vecRestrictions.size();}
    bool bGetRestrictionNameByIndex(int iIndex, Cstring& sRestrName);
    bool bGetRestrictionLongNameByName(Cstring& sRestrName, Cstring& sRestrLongName);
    bool bGetRestrictionInvertedLongNameByName(Cstring& sRestrName, Cstring& sRestrInvertedLongName);

public:
    static int m_nSpaceGroupFreq[]; // Frequencies of spacegroups in Cambridge database

private:
    Cstring     m_sTriclinic_table;
    Cstring     m_sMonoclinic_table;
    Cstring     m_sOrthorhombic_table;
    Cstring     m_sTetragonal_table;
    Cstring     m_sTrigonal_table;
    Cstring     m_sHexagonal_table;
    Cstring     m_sCubic_table;

    std::vector<Cstring>                    m_vecLaueClasses;
    Cstring                                 m_sSpaceGroupFreq;
    std::map<Cstring, enCond>               m_mapConditions;
    std::map<Cstring, DTREK_DWORD>    m_mapRestrictions;
    std::map<Cstring, Cstring>              m_mapRestrictionInvertedLongNames;
    std::map<Cstring, Cstring>              m_mapRestrictionLongNames;
    std::map<DTREK_DWORD, enRestr>    m_mapRestrictionEnum;
    std::map<DTREK_DWORD, int>        m_mapRestrictionParity;
    std::vector<Cstring>                    m_vecRestrictions;
};
///////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
    enUnknown,
    enInconclusive,
    enConclusive,
    enFailed,
    enSuperseded
}enumTestStatus;
class CReflectionRestriction;
class CReflectionCondition;
class DTREK_EXPORT CReflectionRestrictionTableEntry
{
public:
    CReflectionRestrictionTableEntry();
    CReflectionRestrictionTableEntry(Cstring& sName);
    CReflectionRestrictionTableEntry(const CReflectionRestrictionTableEntry& oAnother);
    
    ~CReflectionRestrictionTableEntry(){}
    
    CReflectionRestrictionTableEntry& operator=(const CReflectionRestrictionTableEntry& oAnother);

public:
    void vInit();

    enumTestStatus  eGetTestStatus()const{return m_enTestStatus;}
    void vSetSuperseded(){m_enTestStatus=enSuperseded;}
    bool bIsConclusive()   {return (m_enTestStatus==enConclusive);}
    bool bIsInconclusive() {return (m_enTestStatus==enInconclusive);}
    bool bIsFailed()       {return (m_enTestStatus==enFailed);}
    DTREK_DWORD   dwGetCode(){return m_dwCode;}
    int nGetParity(){return m_nParity;}
    void vGetRestrictionNames(std::vector<Cstring>& asNames);
    bool bSetStatusFromRestrictionsStats(std::map<Cstring, CReflectionRestriction>&   m_mapRestrictions);
    void vSelectRestrictionsInCondition(CReflectionCondition* poCondition);
    void vFillAbsensesArrayOldStyle(CReflectionCondition*, int nRestrictions[4]);

protected:
    DTREK_DWORD   dwGetRestrictionCode(Cstring& sName);
    int    nGetRestrictionParity(DTREK_DWORD dwCode);

protected:
    Cstring m_sName;
    std::vector<Cstring>    m_asNames;
    DTREK_DWORD       m_dwCode;
    int                     m_nParity;
    enumTestStatus          m_enTestStatus;
};
///////////////////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CReflectionRestriction : public CReflectionRestrictionTableEntry
{
public:
    CReflectionRestriction();
    CReflectionRestriction(Cstring& sName, SPACEGROUP_CHECK_CTRL_INFO& m_stCtrlInfo);
    CReflectionRestriction(const CReflectionRestriction& oAnother);
    ~CReflectionRestriction(){}

    CReflectionRestriction& operator=(const CReflectionRestriction& oAnother);

public:
    void vInit();

    bool bTestReflection(int iIndex, int nH, int nK, int nL, double fIntensity, double fSigma);
    void vTallyUpStatistics();

    void vOutputRejects(int iConditionIndex, Cstring& sConditionName, Creflnlist& oOriginalReflnList, Creflnlist& oRejectsList);
    
    int nGetTotalNumberOfCompliantReflns()   {return (m_nHKLCompliantIntegratedReflns    + m_nHKLCompliantNonIntegratedReflns);   }
    int nGetTotalNumberOfNonCompliantReflns(){return (m_nHKLNonCompliantIntegratedReflns + m_nHKLNonCompliantNonIntegratedReflns);} 

    int nGetNumberOfCompliantIntegratedReflns()   {return m_nHKLCompliantIntegratedReflns;}
    int nGetNumberOfNonCompliantIntegratedReflns(){return m_nHKLNonCompliantIntegratedReflns;} 

    double fGetIoverSigmaOfCompliantIntegratedReflns(){return m_fIOverSigmaHKLCompliantIntegratedReflns;} 
    double fGetIoverSigmaOfNonCompliantIntegratedReflns(){return m_fIOverSigmaHKLNonCompliantIntegratedReflns;}

    int nGetNumberOfNonCompliantIntegratedSignificantReflns(){return m_nHKLNonCompliantIntegratedSignificantReflns;}
    int nGetNumberOfCompliantIntegratedSignificantReflns(){return m_nHKLCompliantIntegratedSignificantReflns;}

    void vSelect(bool bSelect){m_bSelected = bSelect;}
    bool bIsSelected(){return m_bSelected;}
    
    void   vGetInvertedLongName(Cstring& sInvName);

    void vFillRestrictionIndexOldStyle(int& nRestrictionIndex);

private:
    bool   bIsHKLCompliant(int nH, int nK, int nL);

    bool   bTestHKLCompliantReflection(double fIntensity, double fSigma);
    bool   bTestHKLNonCompliantReflection(int iIndex, double fIntensity, double fSigma);

    int    nGetTableRestrictionIndex();
    int    nGetTableRestrictionCount();
private:
    SPACEGROUP_CHECK_CTRL_INFO m_stCtrlInfo;
    
    // Number of reflections compliant with a restriction. 
    // For example, if restriction is h=2n, it is the number of reflections with h=2n 
    int     m_nHKLCompliantIntegratedReflns;       // Reflections with valid Intensity and Sigma
    int     m_nHKLCompliantIntegratedSignificantReflns;
    int     m_nHKLCompliantNonIntegratedReflns;    // Reflections w/o valid intensity or sigma

    // Number of reflections non-compliant with a restriction. 
    // For example, if restriction is h=2n, it is the number of reflections with h!=2n 
    int     m_nHKLNonCompliantIntegratedReflns;       // Reflections with valid Intensity and Sigma
    int     m_nHKLNonCompliantIntegratedSignificantReflns;       
    int     m_nHKLNonCompliantNonIntegratedReflns;    // Reflections w/o valid intensity or sigma

    double  m_fIOverSigmaHKLCompliantIntegratedReflns;
    double  m_fIOverSigmaHKLNonCompliantIntegratedReflns;  

    std::vector<int>    m_anRejects;

    bool   m_bSelected;
};

class DTREK_EXPORT CReflectionCondition
{
public:
    CReflectionCondition(Cstring& strID, SPACEGROUP_CHECK_CTRL_INFO& m_stCtrlInfo);
    ~CReflectionCondition(){}

public:
    bool bIsReflectionCondition(int& nH, int& nK, int& nL);

    bool bTestReflection(int iIndex, Crefln& oRefln);

    void vTallyUpStatistics();
    
    enCond  eGetID(){return m_eID;}
    Cstring sGetID(){return m_strID;}
    
    bool bAddSpaceGroupRestriction(Cstring& strRestrictionName);
    bool bAddRestriction(Cstring& strRestrictionName);
    void vMakeFullRestrictionList();

    void vOutputRejects(Creflnlist& oOriginalReflnList, Creflnlist& oRejectsList);

    enumTestStatus eGetRestrictionTableEntryStatus(Cstring& strRestrictionName);

    int nGetNumberOfReflns(){return m_nNumberOfReflns;}

    CReflectionRestriction*            poGetRestrictionPtrByName(Cstring& sName);
    CReflectionRestrictionTableEntry*  poGetRestrictionTableEntryPtrByName(Cstring& sName);

    void vSelectRestriction(Cstring& sName);
    void vUnselectAllRestrictions();
    int nGetNumberOfSelectedRestrictions();

    CReflectionRestriction*            poGetFirstSelectedRestrictionPtr();
    
    bool bIsInconclusiveForCandidateSpaceGroup(){return m_bInconclusiveForCandidateSpaceGroup;}
    void vSetInconclusiveForCandidateSpaceGroup(){m_bInconclusiveForCandidateSpaceGroup = true;}

public:
    Cstring                     m_strID;

private:
    void vSupersedRestrictionTableEntries();

private:
    SPACEGROUP_CHECK_CTRL_INFO      m_stCtrlInfo;

    enCond  m_eID;
    
    int     m_nNumberOfReflns;

    std::map<Cstring, CReflectionRestriction>            m_mapRestrictions;
    std::map<Cstring, CReflectionRestrictionTableEntry>  m_mapRestrictionTableEntries;

    bool    m_bInconclusiveForCandidateSpaceGroup;

};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CSpaceGroupCheckSolution
{
public:
    CSpaceGroupCheckSolution(int iIndex, 
                             int nNumber, 
                             Cstring& sName, 
                             eLaueType eLaue,
                             bool bCentric,
                             int nFreq,
                             bool bNonStandard);
    
    CSpaceGroupCheckSolution(const CSpaceGroupCheckSolution& oAnother);
    ~CSpaceGroupCheckSolution(){}

public:
    int nGetCSDBFrequency(){return m_nFrequency;}
    void vSetProbability(double fProb){m_fProbability = fProb;}
    double fGetProbability(){return m_fProbability;}
    int nGetIndexInLaueClassTable(){return m_nIndexInLaueClassTable;}
    bool bIsCentric(){return m_bCentric;}
    bool bIsChiral(){return m_bChiral;}
    bool bIsStandard(){return !m_bNonStandard;}
    int nGetIntTablesNumber(){return m_nIntTablesNumber;}
    eLaueType eGetLaueType(){return m_eLaue;}
    
    Cstring sGetName(){return m_sName;}
    Cstring sGetPresentation(){return m_sPresentation;}

    Cstring sGetIntTablesName();

private:
    bool bIsChiralCompatibleSpacegroup(int nSpaceGroup_IntTables);

private:
    int m_nIndexInLaueClassTable; 
    int m_nIntTablesNumber;
    
    Cstring m_sName;
    Cstring m_sPresentation;  // RB: The difference between m_sName and m_sPresentation is
                              // that m_sName comes directly from Thad's spacegroup check table, so m_sName 
                              // for some monoclinic cells has cell choice information in parentheses, e.g. (C2), (C3). 
                              // m_sPresentation is just a spacegroup name without cell choice information

    eLaueType m_eLaue;
    bool    m_bCentric;
    bool    m_bChiral;
    
    int     m_nFrequency;
    double  m_fProbability;
    bool    m_bNonStandard;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
class DTREK_EXPORT CSpaceGroupTableEntry
{
public:
    CSpaceGroupTableEntry(int nIndex){m_eTestStatus = enUnknown;
                                      m_nConclusivePassedTests = 0;
                                      m_bSelected = false;
                                      m_nIndexInLaueClassTable = nIndex;
                                      m_strAbsenceTest = "";}

    ~CSpaceGroupTableEntry(){}

public:
    void vAddConditionRestriction(std::pair<Cstring, Cstring>& pair){m_mapConditionsRestrictions.insert(pair);}
    
    bool bIsConclusive()  {return m_eTestStatus == enConclusive;}
    bool bIsInconclusive(){return m_eTestStatus == enInconclusive;}

    int  nGetNumberOfConclusivePassedTests(){return m_nConclusivePassedTests;}
    void vSetSelected(bool bSel=true){m_bSelected=bSel;}
    bool bIsSelected(){return m_bSelected;}

    bool bTest(std::vector<CReflectionCondition*>& vecConditionPtrs);

    int nGetIndexInLaueClassTable(){return m_nIndexInLaueClassTable;}

    void vSelectRestrictionsInConditions(std::vector<CReflectionCondition*>& vecConditionPtrs);
    enumTestStatus eGetConditionStatus(CReflectionCondition* poCond); 
    Cstring sGetConditionRestriction(CReflectionCondition* poCond);

    void vFillAbsensesArrayOldStyle(std::vector<CReflectionCondition*>& vecConditionPtrs, 
                                    int nAbsences[MAX_NUMBER_HKL_CONDITIONS_IN_TABLE][4]);

    bool bGenerateSpaceGroupSolutions(Cstring& sLaueClassAndCentricityInfo, Cstring& sSpaceGroupsInfo);
    int  nAddSpaceGroupSolutionsToVector(std::vector<CSpaceGroupCheckSolution>& vecSpaceGroupSolutions,
                                         eLaueType eLaue,
                                         bool bCentric,
                                         bool bSelectChiralOnly);

    bool bIsSpaceGroupInStandardSetting(int nSG);

    Cstring sGetAbsenceTestSummary(){return m_strAbsenceTest;}
    void vMarkInconclusiveCondtions(std::vector<CReflectionCondition*>& vecConditionPtrs);

private:
    eLaueType eGetLaueClassEnumFromName(Cstring& sName);

private:
    std::map<Cstring, Cstring>      m_mapConditionsRestrictions;
    Cstring                         m_strAbsenceTest;
    std::vector<CSpaceGroupCheckSolution>   m_vecSpaceGroupSolutions;
    enumTestStatus                  m_eTestStatus;
    int                             m_nConclusivePassedTests;
    bool                            m_bSelected;
    int                             m_nIndexInLaueClassTable; 
};
//////////////////////////////////////////////////////////////////////
class DTREK_EXPORT Cspacegroup_check
{
public:
    Cspacegroup_check(SPACEGROUP_CHECK_CTRL_INFO& stInfo);
    ~Cspacegroup_check();

public:
    bool bTestReflnList(Creflnlist& oReflnList);
    void vTestSpaceGroupTableEntries();
//mrp    void vSelectCandidateSpaceGroups();
    void vSelectCandidateSpaceGroups(const int nFavoredSpaceGroupNumber=0); //mrp
    
    int nGetTotalNumberOfConditions(){return (int)m_vecReflnConditionPtrs.size();}
    void vGetConditionNameByIndex(int iIndex, Cstring& sName){sName = m_vecReflnConditionPtrs[iIndex]->sGetID();}
    int nGetNumberOfReflnInConditionByIndex(int iIndex){return m_vecReflnConditionPtrs[iIndex]->nGetNumberOfReflns();}
    CReflectionCondition* poGetConditionPtrByIndex(int iIndex){return m_vecReflnConditionPtrs[iIndex];}
    
    int nGetTotalNumberOfFullListRestrictions();
    void vGetFullListRestrictionNameByIndex(int iIndex, Cstring& sRestrName);
    void vGetFullListRestrictionLongNameByName(Cstring& sRestrName, Cstring& sRestrLongName);
    void vGetFullListRestrictionInvertedLongNameByName(Cstring& sRestrName, Cstring& sRestrInvertedLongName);

    int nGetSelectedCandidateSpaceGroupIndex(){return m_nSelectedCandidateSpaceGroupIndex;}
    void vSetSelectedCandidateSpaceGroupIndex(int nIndex){m_nSelectedCandidateSpaceGroupIndex = nIndex;}

    void vPrintSystematicAbsencesForSelectedSpaceGroup(Cstring& sOut);
    void vPrintAllConditionsRestrictions(Cstring& sOut);
    void vPrintSelectedSpaceGroups(Cstring& sOut);
    
    int  nGetSelectedCandidateSpaceGroupIntTablesNumber();

    SPACEGROUP_CHECK_CTRL_INFO::eCentricityType eGetSelectedCandidateSpaceGroupCentricity();
    int  nGetFirstCandidateSpaceGroupIntTablesNumber();
    
    int  nGetCandidateSpaceGroupIndexBySGNumber(int nSG);
    
    SPACEGROUP_CHECK_CTRL_INFO::eCentricityType eGetCandidateSpaceGroupCentricity(int iIndex);
    
    void vSetSelectedCandidateSpaceGroupIndexBySGNumber(int nSG);
    
    Cstring sGetSelectedCandidateSpaceGroupName();
    Cstring sGetSelectedCandidateSpaceGroupPresentation();
    
    int nGetIndexOfSpaceGroupTableEntryBySelectedCandidateIndex();
    
    int nGetIndexOfSpaceGroupTableEntryContainingSpaceGroupInStandardSetting(int nSG);

    void vFillAbsensesArrayOldStyleFromSpaceGroupTableEntry(int iEntry, int nAbsences[MAX_NUMBER_HKL_CONDITIONS_IN_TABLE][4]);
    
    Cstring sGetListInconclusiveConditionsForCandidateSpaceGroups();
    eIndexTransform eFindTranslationPermForSelectedSpaceGroup();

public:    
    TEXSAN_STATS    m_stTexsanStats;

private:
    SPACEGROUP_CHECK_CTRL_INFO                     m_stCtrlInfo;

    std::vector<CReflectionCondition*>             m_vecReflnConditionPtrs;

    std::vector<CSpaceGroupTableEntry*>            m_vecSpaceGroupTableEntryPtrs; 
    std::vector<CSpaceGroupCheckSolution>          m_vecCandidateSpaceGroups; 
    
    int                                            m_nSelectedCandidateSpaceGroupIndex;

private:
    void vLoadTableByLaueType();
    void vLoadTable(Cstring& sInput);
    
    void vAddCondition(Cstring& strCondition);

    CReflectionCondition*   poGetConditionByID(enCond eID);
    
    bool bAddRestrictionToCondition(int iConditionIndex, Cstring& strRestrictionName);

    void vMakeFullRestrictionListForConditionHKL();

    void vOutputRejects(Creflnlist& oOriginalReflnList);

    void vDoTexsanStats(Crefln& oRefln);

    int nGetTotalNumberOfSelectedRestrictionsInAllConditions();

    int nGetSpaceGroupSolutionsFromSpaceGroupTableEntry(int iIndex, 
                                                        SPACEGROUP_CHECK_CTRL_INFO::eCentricityType eType);
    void vFindProbabilitiesForCandidateSpaceGroups();
    void vMarkInconclusiveConditionsForCandidateSpaceGroupTableEntries(std::vector<CSpaceGroupTableEntry*>& vecCandidates);
    
    int nFindBestCandidateSpaceGroup(SPACEGROUP_CHECK_CTRL_INFO::eCentricityType eType);
    
    void vSelectRestrictionsFromSelectedCandidateSpaceGroup();

    Cstring sGetAbsenceTestSummary(int iCandidateSpaceGroupIndex);

    void vDeselectCandidateSpaceGroupTableEntriesByComparison(std::vector<CSpaceGroupTableEntry*>& vecCandidates);
    
    Cstring sGetLaueName();
};
#endif   // DT_CSPACEGROUPCHECK_H

