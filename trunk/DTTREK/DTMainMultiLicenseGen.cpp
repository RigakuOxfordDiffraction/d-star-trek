//
// Copyright (c) 2007 Rigaku
//
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// DTMainMultiLicenseGen.cpp   Initial author: RB           07-Mar-2007
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
 
#include "DTMainMultiLicenseGen.h"
#include "DTMainLicenseGen.h"

using namespace std;

const char*     c_pcMainMultiLicenseGenTerm                   = "-term";  
const char*     c_pcMainMultiLicenseGenInfo                   = "-info";  
const char*     c_pcMainMultiLicenseGenLinux                  = "-linux";  

int CDTMainMultiLicenseGenerator::nExecute(unsigned int argc,    // Command line argument count
                                           char      *argv[])    // Pointers to command line args
{
    vDtrekSetModuleName("dtlicensegenWIN");
    cout << "\ndtlicensegenWIN:  Copyright (c) 2007 Rigaku\n";
    cout << D_K_DTREKVersion << endl;

    // Copy command line to output log
    cout << "\nCommand line:\n" << sGetCommandLine(argc, argv, 71) << endl << endl << flush;
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    vSetCommandArguments(argc, argv);
    
    Cstring     sError("");
    
    // Read the license term
    Cstring     sTerm("");
    if( !bGetCommandOption(c_pcMainMultiLicenseGenTerm, sTerm) )
    {
        DTREK_ERROR(1, "ERROR - no license term specified\n");
    }
    
    // Read the hardware information file path
    Cstring     sInfoFile("");
    if( !bGetCommandOption(c_pcMainMultiLicenseGenInfo, sInfoFile) )
    {
        DTREK_ERROR(1, "No hardware information file specified\n");
    }
    
    // See if this is a linux (Jim's) info file
    bool    bLinux = false;
    bLinux = bGetCommandOption(c_pcMainMultiLicenseGenLinux);


    FILE*       fp = fopen(sInfoFile.string(), "rt");
    if( !fp )
    {
        sError = "Cannot open hardware information file ";
        sError += Cstring(sInfoFile);

        DTREK_ERROR(1, sError);
    }
    const  int     cnSize = 1024;
    char        pcBuf[cnSize];
    Cstring     sLineRead("");
    
    CDTMainLicenseGenerator     oMainLicenseGenerator;
    int         nRet = -1;
    Cstring     sCommand("");
    Cstring     sLicense("");
    
    std::vector<Cstring>    vecLicenses;
    if( !bLinux )
    {
		// Assuming a format like:
		// CPU ID: 138264140 1659169773		
        while( fgets(pcBuf, sizeof(pcBuf), fp) )
        {
            sLineRead = pcBuf;
            
            if( !sLineRead.contains("CPU ID:") )
            {
                cout << pcBuf << endl << flush;
                continue;
            }

            std::vector<Cstring>       vecTokens;
            sLineRead.nListToVector(vecTokens, " ");
            
            sCommand = "dtlicensegen.exe ";

            sCommand += sTerm;
            sCommand += " ";
            sCommand += vecTokens[2];
            sCommand += " ";
            sCommand += vecTokens[3];
            
            vCommandLineToCommandArgumentArray(sCommand.string());
            
            oMainLicenseGenerator.vInit();
            
            cout << endl << flush;

            nRet =  oMainLicenseGenerator.nExecute(m_nDTREKCommandArgumentCount, m_ppcDTREKCommandArguments);
            
            oMainLicenseGenerator.vGetLicenseString(sLicense);

            vecLicenses.push_back(sLicense);
        }
    }
    else
    {
        //Assuming a format like:
		//Linux idesc 2.6.9-42.0.10.ELsmp #1 SMP Fri Feb 16 17:17:21 EST 2007 i686 i686 i386 GNU/Linux
		//d*TREK version:   d*TREK version 9.7LDz -- Jun 14 2007  This Host name: idesc   This CPU ID,    Host ID:
		//110843483 1330121899
		bool        bReadyToRead = false;
        while( fgets(pcBuf, sizeof(pcBuf), fp) )
        {
            sLineRead = pcBuf;
            if( bReadyToRead )
            {
                // remove all "weird" characters that might have come from Linux or e-mail program
                for(int ii=0; ii < cnSize; ii++)
                {
                    if( !( (sLineRead.GetAt(ii) >= '0' && sLineRead.GetAt(ii) <= 'Z') ||
                            sLineRead.GetAt(ii) == '\n' ) ) 
                        sLineRead.SetAt(ii, ' ');
                }

                std::vector<Cstring>       vecTokens;
                char* c_psep = " ";
                sLineRead.nListToVector(vecTokens, c_psep);
                
                int     nn = (int)vecTokens.size();

                sCommand = "dtlicensegen.exe ";

                sCommand += sTerm;
                sCommand += " ";
                sCommand += vecTokens[0]; // because the first token is a white space
                sCommand += " ";
                sCommand += vecTokens[1];
                
                vCommandLineToCommandArgumentArray(sCommand.string());
                
                oMainLicenseGenerator.vInit();
                
                cout << endl << flush;

                nRet =  oMainLicenseGenerator.nExecute(m_nDTREKCommandArgumentCount, m_ppcDTREKCommandArguments);
                
                oMainLicenseGenerator.vGetLicenseString(sLicense);

                vecLicenses.push_back(sLicense);
                
                bReadyToRead = false;
            }
            else if( sLineRead.contains("CPU ID,") )
            {
                cout << pcBuf << endl << flush;
                bReadyToRead = true;
                continue;
            }
        }
    }

    fclose(fp);

    fp = fopen(sInfoFile.string(), "at");
    if( !fp )
    {
        sError = "Cannot open hardware information file for update ";
        sError += Cstring(sInfoFile);

        DTREK_ERROR(1, sError);
    }
    
    int     ii = 0;
    for(ii=0; ii < (int)vecLicenses.size(); ii++)
    {
        fprintf(fp, "%s\n", vecLicenses[ii].string());
    }
    
    fclose(fp);

    fp = fopen("license_key.txt", "wt");
    if( !fp )
    {
        sError = "Cannot open license info file.";
        sError += Cstring(sInfoFile);

        DTREK_ERROR(1, sError);
    }

    for(ii=0; ii < (int)vecLicenses.size(); ii++)
    {
        fprintf(fp, "%s\n", vecLicenses[ii].string());
    }
    
    fclose(fp);
    
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////
void CDTMainMultiLicenseGenerator::vError(const int nErrorNum, const Cstring& sMessage)
{
    if( 0 != nErrorNum ) 
    {
        Cstring     sErrorMsg = "\nERROR: ";
        sErrorMsg += sMessage;
        sErrorMsg += "\n";

        printf("%s\n",sErrorMsg.string());
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////








