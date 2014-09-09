//
// Copyright (c) 2007 Rigaku
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMainScreener.cpp  Initial author: RB           16-Mar-2007
//
//  This file contains the member functions of class CDTMainScreener
//

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

#include "DTMainScreener.h"
#include "ResoBins.h"
#include "DTMainStrategy.h"
#include "DTMainPredict.h"

int CDTMainScreener::nExecute(unsigned int argc, char* argv[])
{
    vSetCommandArguments(argc, argv);

    if( bIsHelpRequest() )
    {
        DTREK_ERROR(0, "");  // display general help
    }
    
    CDTMainPredict      oMainPredict;
    double  fRotMin = -90.;
    double  fRotMax = 72.;
    
    Cstring     sOption(fRotMin);
    sOption += ", ";
    Cstring     sRotMax(fRotMax);
    sOption += sRotMax;


    bUpdateCommandOption("-rot", sOption);
    
    //bUpdateCommandOption("-rot", Cstring("-180, 180"));


    oMainPredict.nExecute(m_nDTREKCommandArgumentCount, m_ppcDTREKCommandArguments);

//#if 0
    //if( argc < 2 )
    //{
        //DTREK_ERROR(1, "No header filename specified\n");
    //}
    //
    //Cstring     sError("");
    //char        acTemp[DTREK_MAX_LINE_LENGTH];
    //
    //// Create a header object
    //Cstring     sHeader(argv[1]);  // assuming that the first argument is a header
    //if( !bCreateHeader(sHeader) )
    //{
        //sError = "Invalid header file ";
        //sError += Cstring(argv[1]);
        //
        //DTREK_ERROR(1, sError);
    //}
//
    //Cpredict*       poPredict = new Cpredict(*m_poHeader, NULL, 0, "");
//
    //if( !poPredict->bIsAvailable() )
    //{
        //sprintf(acTemp, "Could not create predict object.");
        //sError = acTemp;
//
        //delete  poPredict;
        //poPredict = NULL;
//
        //DTREK_ERROR(1, sError);
//
        //return false;
    //}
//
    //if( 0 != poPredict->nSetupSourceCrystal() )
    //{
        //sprintf(acTemp, "Could not set-up the source and crystal info in the predict object.");
        //sError = acTemp;
//
        //delete  poPredict;
        //poPredict = NULL;
        //
        //DTREK_ERROR(1, sError);
    //}
    //
    //Creflnlist*     poReflnList = new Creflnlist();
    //if( 0 != poPredict->nPredict(NULL, poReflnList) ) 
    //{
        //sprintf(acTemp, "Failure during reflection prediction");
        //sError = acTemp;
//
        //delete poPredict;
        //poPredict = NULL;
//
        //delete poReflnList;
        //poReflnList = NULL;
//
        //DTREK_ERROR(1, sError);
    //}
    
    Creflnlist*     poReflnList = new Creflnlist("dtpredict.ref");

    int     nNumberOfPredictedReflns = poReflnList->nGetNumReflns();
    if( 0 == nNumberOfPredictedReflns )
    {
        printf("\nWARNING: Zero reflections predicted for scan.\n");
    }
    
    //poReflnList->nWrite("screener.ref");

    //CDTMainStrategy     oMainStrategy;
    //double  fTemp1=0.0, fTemp2=0.0;
    //oMainStrategy.vDoExternalReflnListCalcSetDStarCubed(m_poPreviousReflnList, fTemp1, fTemp2);
    //
    //oMainStrategy.vDoExternalReflnListSortReduce(m_poPreviousReflnList);
   
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate bins for overall resolution range
    CResoStrategyBins*      poResoBins = new CResoStrategyBins(900.0, 0.8, 4); 

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Now do the actual statistics

    poReflnList->vToStratResoBins(poResoBins);
    
    double      fCoeff = 18.0/(fRotMax - fRotMin);

    CResoStrategyBin*       pBin = NULL;
    printf("\n\n");
    printf(" #     Low     High     Num   Num\n");
    printf("      reso     reso  reflns  rejs\n");
    printf("\n");
    for(int ii=0; ii < poResoBins->nGetNumberOfBins(); ii++)
    {
        pBin = poResoBins->pGetBin(ii);

        printf("%2d %8.2f %8.2f  %5d %5d\n",
                              ii+1,
                              pBin->dGetLowReso(), 
                              pBin->dGetHighReso(), 
                              (int)((double)(pBin->nGetNumRefls()) * fCoeff), 
                              (int)((double)pBin->nGetNumRejs()  * fCoeff));
    }



    //delete poPredict;
    delete poReflnList;
    delete poResoBins;
//#endif

    return 0;
}

void CDTMainScreener::vError(const int nErrorNum, const Cstring& sMessage)
{
  if (0 != nErrorNum)
  {
    std::cout << sMessage << std::endl;
  }
  
  std::cout << "\ndtextractheader - Usage:\n"
       << "      dtextractheader imagefile newheaderfile\n\n";
}
