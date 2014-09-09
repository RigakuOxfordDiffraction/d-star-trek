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

// CrclHelper.cpp: implementation of the CCrclHelper class.
//
//////////////////////////////////////////////////////////////////////

#include "CrclHelper.h"
#include "dtreksys.h"
#include <time.h>
#include <io.h>


#if !defined(VC6)
using std::ios;
#endif

//#include "Cindex.h"
//#include "CoreTclParam.h"
//#include "CoreSocket.h"
//#include "dtcell.h" // for tagLaueTable
//#include <limits.h> // for INT_MAX
//#include <assert.h>

// Static initializer
CCrclHelper* CCrclHelper::m_pInstance = NULL;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


CCrclHelper* CCrclHelper::GetInstance()
{
    if(!m_pInstance) {
        m_pInstance = new CCrclHelper();
    }
    return m_pInstance;
}

void CCrclHelper::DeleteInstance()
{
    if(m_pInstance) {
        delete m_pInstance;
        m_pInstance = NULL;
    }
}

CCrclHelper::CCrclHelper()
{
    m_fpDLL = NULL;
}

CCrclHelper::~CCrclHelper()
{
    if( m_fpDLL ) {
        fclose( m_fpDLL );
    }
}

void CCrclHelper::RedirectStdout( const char* sFilename )
{
    if( m_fpDLL ) {
        fclose( m_fpDLL );
    }

    ios::sync_with_stdio();

	// Version the log file.
    if(bFileExists(sFilename) && lFileGetSize(sFilename) > 0 ) {
        nFileAppendVersion(sFilename, true);
    }

	// Redirect stdout to a log file.
	m_copy = _dup(1);
	m_fpDLL = freopen(sFilename, "w", stdout);

    if( m_fpDLL ) {
        // Make sure that printf calls are not buffered.
        setvbuf(m_fpDLL, NULL, _IONBF, 0);
    }
    else {
        fprintf( stderr, "Error redirecting stdout.\n" );
    }
}

void CCrclHelper::RestoreStdout()
{
    if( m_fpDLL ) {
        fclose( m_fpDLL );
        m_fpDLL = NULL;

		_dup2(m_copy,1);
    }
}

FILE* CCrclHelper::GetLogFileStream()
{
    return m_fpDLL;
}

void CCrclHelper::LogPrint( const char* str )
{
    printf( str );
}

Cstring CCrclHelper::GetTimestamp()
{
	struct tm when;
	time_t now;
	time( &now );
	when = *localtime( &now );
   	Cstring sTimestamp;
	sTimestamp.Format("This log file created %s", asctime( &when ) );
	return sTimestamp;
}

#ifndef VC9
Cstring CCrclHelper::SendTclCmd( istream_withassign stream ) {return "";}
#else
Cstring CCrclHelper::SendTclCmd( std::istream& stream ) {return "";}
#endif

//void CCrclHelper::SendTclCmd( istream_withassign stream, Cstring &response ) {}
void CCrclHelper::SendTclCmd() {}
void CCrclHelper::SendTclCmdRespond() {}
void CCrclHelper::vUtilCmdInit( const char* pcCommand ) {}
void CCrclHelper::vUtilCmdSetCommand( const char* pcCommand ) {}
void CCrclHelper::vUtilCmdSetParam( const char* pcName, const char* pcString ) {}
void CCrclHelper::vUtilCmdSend() {}
void CCrclHelper::vBuildCmdCellReindex( int nTotRefs, int nOffRefs, float fAvgRef, float fAvgISig ) {}
void CCrclHelper::vBuildCmdCellWriteRefFile() {}
void CCrclHelper::vBuildCmdCellUse( double pfUserCell[6], double pfUserCellFound[6],
		double pfOrientAngles[3], double fErr, double a3x3fMat[3][3] ) {}
void CCrclHelper::vBuildCmdCellChooseSoln( CCellReduce& oReduce, Claue& oLaue ) {}
void CCrclHelper::vBuildCmdReduceCellInput( Claue& oLaue ) {}
void CCrclHelper::vBuildCmdCellSetOrients( int numOrients, double a3x3fMat[24][3][3], double a3x3fOrientMat[24][3][3] ) {}
void CCrclHelper::vBuildCmdCellNewCellFromLaue(double oldCell[6], double newCell[6], double oldRot[3], double newRot[3], int newLaueGroup) {}
void CCrclHelper::vResetLaueTable() {}
void CCrclHelper::vAddLaueTableRow( tagLaueTable &row ) {}
void CCrclHelper::vBuildCmdLaueTable( const char *sSelectedLaue, const char *sSelectedLattice, float fRmerge ) {}
void CCrclHelper::vBuildCmdCellReindexForAxis( char axis ) {}
void CCrclHelper::vBuildCmdCellCentricityTable( double af10MeasuredPercent[], double af10MeasuredDeviation[], double fAvgDeviation, bool centric ) {}
void CCrclHelper::vBuildCmdCellAbsenceTable( Cstring &str ) {}
void CCrclHelper::vBuildCmdCellAllAbsenceTable( Cstring &str ) {}
char* CCrclHelper::GetAllAbsenceRowData( Cstring &cmd, char *pc, const char* delim, int numCols ) { return NULL; }
void CCrclHelper::vBuildCmdCellSpacegroupTable( Cstring &str ) {}
void CCrclHelper::vBuildCmdCellSpacegroupSelected( int sg, Cstring& strInconclusive, bool warning ) {}
void CCrclHelper::vBuildCmdCellSpacegroupTransform( const char *transform ) {}
void CCrclHelper::vBuildCmdCellNumReindexRelativeChecks( int numChecks ) {}
void CCrclHelper::vBuildCmdCellReindexRelativeStep() {}
void CCrclHelper::vAddReindexRelativeMat( double a3x3fSelectedMat[3][3], float fSelectedScale, float fRmerge ) {}
void CCrclHelper::vBuildCmdCellReindexRelativeMats() {}

bool CCrclHelper::bBuildCmdCellOutputFile(const char *File, unsigned short wFileType){return false;}
void CCrclHelper::vBuildCmdCellLeftHandedTransform(bool bYes){}

void CCrclHelper::vBuildCmdYesNoQuestion( const char *msg ) {}
void CCrclHelper::vSendMessage( const char *msg ) {}
void CCrclHelper::vBuildCmdCellValidLaueSels( bool *bValidSel, int num ) {}
void CCrclHelper::vResetRmergeReslnRow() {}
void CCrclHelper::vAddRmergeReslnRow( double resoLow, double resoHigh, Cstring &avgCounts,
    int rej, int mults, Cstring &ISig, Cstring &ISig2, Cstring &Eadd, Cstring &chiSq, Cstring &rmergeShell, Cstring &rmergeCumul ) {}
void CCrclHelper::vCompleteRmergeReslnRow() {}
void CCrclHelper::vResetCompReslnRow() {}
void CCrclHelper::vAddCompReslnRow( double resoLow, double resoHigh, int calcUnique, int obs,
        int rej, int mults, int single, int unique, Cstring &avgMult, Cstring &shell, Cstring &cumul ) {}
void CCrclHelper::vCompleteCompReslnRow() {}
void CCrclHelper::vResetRmergeBatchRow() {}
void CCrclHelper::vAddRmergeBatchRow( Cstring &name, int avgCounts, int obs, int rej, int ovlps, int single,
    Cstring &ISig, Cstring &ChiSq, Cstring &rmergeBatch, Cstring &rmergeCumul ) {}
void CCrclHelper::vCompleteRmergeBatchRow() {}
void CCrclHelper::vResetRmergeISigRow() {}
void CCrclHelper::vAddRmergeISigRow( double rangeMin, double rangeMax, Cstring counts, int obs, int rej, //int ovlps,
    int mults, Cstring &ISig, Cstring &ISig2, Cstring &ChiSqNorm, Cstring &rmergeShell, Cstring &rmergeCumul ) {}
void CCrclHelper::vCompleteRmergeISigRow() {}
Cstring CCrclHelper::sSendIndexUseFoundCell( double a6fCell[], double a6fFoundCell[], float orient[], double error ) { return ""; }
void CCrclHelper::vSendUpdateDisplay( const char* sCmd, const char* sImage, const char* sReflnlist, bool bObs, bool bCalc, bool bDiff ) {}
void CCrclHelper::vSendFindUpdateDisplay( const char* sImage, const char* sReflnlist ) {}
void CCrclHelper::vSendRefineUpdateDisplay( const char* sImage, const char* sReflnlist ) {}
void CCrclHelper::vSendRefineUpdateDisplay( const char* sReflnlist ) {}
void CCrclHelper::vSendPredictUpdateDisplay( const char* sReflnlist ) {}
void CCrclHelper::vSendPredictUpdateDisplay( const char* sImage, const char* sReflnlist ) {}
void CCrclHelper::vSendIntegrateUpdateDisplay( const char* sImage, const char* sReflnlist, bool bObs, bool bCalc, bool bDiff ) {}
void CCrclHelper::vSendCollectionStrategy( double dReslnMin, double dReslnMax, double dRotStart, double dRotEnd ) {}
Cstring CCrclHelper::sSendIndexTable( CCellReduce& oReduce ) { return ""; }
void CCrclHelper::vSendIndexOrientTable( CCellReduce& oReduce, int nChoose ) {}
Cstring CCrclHelper::sBuildIndexTableRow( Crefln* poSoln, int nIndex, Cstring sSolnNum, 
    int nLattice, int nLattSymm, int nA, int nB, int nC, int nVolume, int nAlpha, 
    int nBeta, int nGamma ) { return ""; }
Cstring CCrclHelper::sBuildReduceCellInputRow(double a, double b, double c, 
					double alpha, double beta, double gamma,
					double Rot1, double Rot2, double Rot3, double vol ) { return ""; }
Cstring CCrclHelper::GetInstallDir() { return ""; }
Cstring CCrclHelper::sPromptIncompletenessStrategy(double dCompleteness, double dRedundancy) { return ""; }
