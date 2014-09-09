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

// CrclHelper.h: interface for the CCrclHelper class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CRCLHELPER_H__2AD4EDA1_E0BF_11D3_80F5_005004D38A88__INCLUDED_)
#define AFX_CRCLHELPER_H__2AD4EDA1_E0BF_11D3_80F5_005004D38A88__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "Cstring.h"

class Cindex;
class CCellReduce;
class Crefln;
class Claue;
class CCoreTclParam;
class CCoreSocket;
struct tagLaueTable;

#ifdef EXPORT_DTREK_DLL
#define DTREK_EXPORT _declspec(dllexport)
#elif defined IMPORT_DTREK_DLL
#define DTREK_EXPORT _declspec(dllimport)
#else
#define DTREK_EXPORT
#endif

#ifndef VC9
#else
        // don't do this in a header file! see: http://www.parashift.com/c++-faq-lite/coding-standards.html#faq-27.5
        //using namespace std;
#endif


/*
 * Class CCrclHelper
 * This is a singleton class used to provide interface for d*Trek to CrystalClear
 * communication.
 */
class DTREK_EXPORT CCrclHelper  
{
public:
    static CCrclHelper* GetInstance();
    static void DeleteInstance();

public:
    //virtual void SetInteractive( bool b ) { m_bInteractive = b; }
    virtual void RedirectStdout( const char* sFilename );
    virtual void RestoreStdout();
    virtual FILE* GetLogFileStream();
    virtual void LogPrint( const char* str );
	virtual Cstring GetTimestamp();
#ifndef VC9
    virtual Cstring SendTclCmd( istream_withassign stream );
#else
    virtual Cstring SendTclCmd( std::istream& stream );
#endif
	//virtual void SendTclCmd( istream_withassign stream, Cstring &response );
    virtual void SendTclCmd();
    virtual void SendTclCmdRespond();

    virtual void vUtilCmdInit( const char* pcCommand = "" );
    virtual void vUtilCmdSetCommand( const char* pcCommand );
    virtual void vUtilCmdSetParam( const char* pcName, const char* pcString );
    virtual void vUtilCmdSend();

    // Tcl command builders
    virtual void vBuildCmdCellReindex( int nTotRefs, int nOffRefs, float fAvgRef, float fAvgISig );
    virtual void vBuildCmdCellWriteRefFile();
	virtual void vBuildCmdCellUse( double pfUserCell[6], double pfUserCellFound[6],
		double pfOrientAngles[3], double fErr, double a3x3fMat[3][3] );
    virtual void vBuildCmdCellChooseSoln( CCellReduce& oReduce, Claue& oLaue );
    virtual void vBuildCmdReduceCellInput( Claue& oLaue );
    virtual void vBuildCmdCellSetOrients( int numOrients, double a3x3fMat[24][3][3], double a3x3fOrientMat[24][3][3] );
    virtual void vBuildCmdCellNewCellFromLaue(double oldCell[6], double newCell[6], double oldRot[3], double newRot[3], int newLaueGroup);
    virtual void vResetLaueTable();
    virtual void vAddLaueTableRow( tagLaueTable &row );
    virtual void vBuildCmdLaueTable( const char *sSelectedLaue, const char *sSelectedLattice, float fRmerge );
    virtual void vBuildCmdCellReindexForAxis( char axis );
	virtual void vBuildCmdCellCentricityTable( double af10MeasuredPercent[], double af10MeasuredDeviation[], double fAvgDeviation, bool centric );
    virtual void vBuildCmdCellAbsenceTable( Cstring &str );
    virtual void vBuildCmdCellAllAbsenceTable( Cstring &str );
    virtual void vBuildCmdCellSpacegroupTable( Cstring &str );
    virtual void vBuildCmdCellSpacegroupSelected( int sg, Cstring& strInconclusive, bool warning );
    virtual void vBuildCmdCellSpacegroupTransform( const char *transform );
    virtual void vBuildCmdCellNumReindexRelativeChecks( int numChecks );
    virtual void vBuildCmdCellReindexRelativeStep();
    virtual void vAddReindexRelativeMat( double a3x3fSelectedMat[3][3], float fSelectedScale, float fRmerge );
    virtual void vBuildCmdCellReindexRelativeMats();
    virtual void vBuildCmdCellLeftHandedTransform(bool bYes);
    
#define DTREK_CCHI_BCCOF_FILE_REF       0x0001
#define DTREK_CCHI_BCCOF_FILE_HDR       0x0002
     virtual bool bBuildCmdCellOutputFile(const char *File, unsigned short wFileType);

    virtual void vBuildCmdYesNoQuestion( const char *msg );
    virtual void vSendMessage( const char *msg );
    virtual void vBuildCmdCellValidLaueSels( bool *bValidSel, int num );

    // Scale average tables
    // Rmerge vs Resln
    virtual void vResetRmergeReslnRow();
    virtual void vAddRmergeReslnRow( double resoLow, double resoHigh, Cstring &avgCounts,
        int rej, int mults, Cstring &ISig, Cstring &ISig2, Cstring &chiSq, Cstring &Eadd,Cstring &rmergeShell, Cstring &rmergeCumul );
    virtual void vCompleteRmergeReslnRow();

    // Completeness vs Resln
    virtual void vResetCompReslnRow();
    virtual void vAddCompReslnRow( double resoLow, double resoHigh, int calcUnique, int obs,
        int rej, int mults, int single, int unique, Cstring &avgMult, Cstring &shell, Cstring &cumul );
    virtual void vCompleteCompReslnRow();

    // Rmerge vs Batch ID
    virtual void vResetRmergeBatchRow();
    virtual void vAddRmergeBatchRow( Cstring &name, int avgCounts, int obs, int rej, int ovlps, int single,
        Cstring &ISig, Cstring &ChiSq, Cstring &rmergeBatch, Cstring &rmergeCumul );
    virtual void vCompleteRmergeBatchRow();

    // Rmerge vs I/Sig
    virtual void vResetRmergeISigRow();
    virtual void vAddRmergeISigRow( double rangeMin, double rangeMax, Cstring counts, int obs, int rej, //int ovlps,
        int mults, Cstring &ISig, Cstring &ISig2, Cstring &ChiSqNorm, Cstring &rmergeShell, Cstring &rmergeCumul );
    virtual void vCompleteRmergeISigRow();

    virtual Cstring sSendIndexUseFoundCell( double a6fCell[], double a6fFoundCell[], float orient[], double error );

    // Display udpate stuff
    virtual void vSendUpdateDisplay( const char* sCmd, const char* sImage, const char* sReflnlist, bool bObs, bool bCalc, bool bDiff );
    virtual void vSendFindUpdateDisplay( const char* sImage, const char* sReflnlist );
    virtual void vSendRefineUpdateDisplay( const char* sImage, const char* sReflnlist );
    virtual void vSendRefineUpdateDisplay( const char* sReflnlist );
    virtual void vSendPredictUpdateDisplay( const char* sReflnlist );
    virtual void vSendPredictUpdateDisplay( const char* sImage, const char* sReflnlist );
    virtual void vSendIntegrateUpdateDisplay( const char* sImage, const char* sReflnlist, bool bObs, bool bCalc, bool bDiff );

    virtual void vSendCollectionStrategy( double dReslnMin, double dReslnMax, double dRotStart, double dRotEnd );

    virtual Cstring sSendIndexTable( CCellReduce& oReduce );
    virtual void vSendIndexOrientTable( CCellReduce& oReduce, int nChoose );

    virtual Cstring sBuildIndexTableRow( Crefln* poSoln, int nIndex, Cstring sSolnNum, int nLattice, int nLattSymm, 
        int nA, int nB, int nC, int nVolume, int nAlpha, int nBeta, int nGamma );

	virtual Cstring sBuildReduceCellInputRow(double a, double b, double c, 
					double alpha, double beta, double gamma,
					double Rot1, double Rot2, double Rot3, double vol );

    virtual Cstring GetInstallDir();
	virtual Cstring sPromptIncompletenessStrategy(double dCompleteness, double dRedundancy);

    /*
     * SetSocket
     * CrclHelper does NOT take ownership of the socket.  Some other object must be
     * responsible for deletion of the socket.
     */
    //virtual void SetSocket( CCoreSocket* pSocket ) { m_pSocket = pSocket; }

    //protected methods
protected:
    virtual char* GetAllAbsenceRowData( Cstring &cmd, char *pc, const char* delim, int numCols );


protected:
    FILE* m_fpDLL;
	int m_copy;
    /*Cstring m_csLastCmd;
    bool m_bInteractive;
    tagLaueTable* m_vLaueTable;
    int m_iNumRowsInLaueTable;
    int m_iNumReindexRelMats;
    int m_iNumRmergeReslnRows;
    int m_iNumCompReslnRows;
    int m_iNumRmergeBatchRows;
    int m_iNumRmergeISigRows;*/

    // A utility command for when d*Trek modules need to do extensive building
    // of a tcl command using CCoreTclParam.  Since we want CCoreTclParam
    // out of the d*Trek code, we provide this along with helper methods.
    //CCoreTclParam* m_pUtilCmd;

    // The socket to send commands across with.
    //CCoreSocket* m_pSocket;

// Protected constructor/destructor
protected:
	CCrclHelper();
	virtual ~CCrclHelper();

// The instance
protected:
    static CCrclHelper* m_pInstance;
};

#endif // !defined(AFX_CRCLHELPER_H__2AD4EDA1_E0BF_11D3_80F5_005004D38A88__INCLUDED_)
