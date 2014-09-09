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
// dtimageheaderedit.cc        Initial author: Thaddeus Niemeyer   Jan-2003
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
#include "CCommandline.h"
#include "Cstring.h"
#include "dtarray.h"
#include "Cimage.h"
#include <io.h>
#include "Cdetector.h"



int nReplaceVariables(Cstring& sKeyword,itr<Cstring>& asVarNames,itr<Cstring>& asVarValues) {
    int nVar;
    Cstring sTemp;

    // Replace all variables in the string.
    do {
        for (nVar = 0; nVar < asVarNames.size(); nVar++) {
            if (sKeyword.contains(asVarNames[nVar])) {
                sTemp = sKeyword.before(asVarNames[nVar]);
                sTemp += asVarValues[nVar];
                sTemp += sKeyword.after(asVarNames[nVar]);
                sKeyword = sTemp;
                break;
            };
        };
    } while (nVar < asVarNames.size());
    return 0;
};

int main(int argc,char** argv) {
	


	const char* apcOptions[] = {
        "-START","s[sFileSpec]","r1",
			"-distance", "f[fDistance]",
			"-2theta","f[f2Theta]",
			"-beamcenter","f[fBeam0]","f[fBeam1]",
            "-set","s[sVariableName]","s[sHeaderItemName]","s[sStringPattern]","s[sStringPatternOut]",
            "-setadd","s[sVariableName]","s[sExp1]","s[sExp2]",
            "-setmul","s[sVariableName]","s[sExp1]","s[sExp2]",
			"-add","s[sHeaderItemName]","s[sStringOut]",
			"-replace","s[sHeaderItemName]","s[sOut]",
			"-replace","s[sHeaderItemName]","*s[sStringPatternIn]","s[sStringPatternOut]",
			"-remove","s[sHeaderItemName]",
			"-list","s[sHeaderItemName]","r1",
			"-prefix","s[sPrefixForOutputFilename]",
			NULL
	};
	const char cAddOnlyFlag = '\x1';
	const char cRemoveOnlyFlag = '\x2';
	const char cListOnlyFlag = '\x3';
    const Cstring sAddString = "_add_";
    const Cstring sMultString = "_mult_";
	int nArg;
	int nStat;
	int nx,ny;
	Cstring sTemp;
	
	itr<Cstring> asFilesToModify;

	itr<Cstring> asOptionNameToFind;
	itr<Cstring> asOptionDataIn;
	itr<Cstring> asOptionDataOut;

    itr<Cstring> asVariable;
    itr<Cstring> asVariableValue;
    itr<Cstring> asVariableExtractKeyword;
    itr<Cstring> asVariableExtractPatternIn;        // Or equal to First expression if asVariableExtractKeyword[*]==""
    itr<Cstring> asVariableExtractPatternOut;       // Or equal to Second expression if asVariableExtractKeyword[*]==""
	Cstring sPrefix;
	
	sPrefix = "";
	-asFilesToModify;
	nArg = -1;

    // Define variables that are always available.
    asVariable + (Cstring) "%DET%";
    asVariableExtractKeyword + Cdetector::ms_sDetectorNames;
    asVariableExtractPatternIn + (Cstring) "*";
    asVariableExtractPatternOut + (Cstring) "%1";
    asVariableValue + (Cstring) "";


	while (0 == (nStat = nReadCommandLine(nArg,argc,argv,apcOptions))) {
		if (g_sOptionArg == "-START") {
			for (nx=0;nx<g_asOptionBuf.size();nx++) {
				_finddata_t oFind;
				int nHandle;
				Cstring sDirectory;
				sDirectory = sFileGetDirectory(g_asOptionBuf[nx]);
				nHandle = _findfirst(g_asOptionBuf[nx],&oFind);
				ny = 0;
				if (-1 != nHandle) {
					do {
						asFilesToModify + sFileBuildName(sDirectory,(Cstring) oFind.name);
						ny++;
					} while (-1 != _findnext(nHandle,&oFind));
					_findclose(nHandle);
					
				} else {
					printf("ERROR:  Could not find any matches to file spec. '%s'\n",g_asOptionBuf[nx].string());
					nStat = 2;
					break;
				};
				if (!nStat) {
					printf("Found %d matches to file spec. '%s'\n",ny,g_asOptionBuf[nx].string());
				};
			};		
		} else if (g_sOptionArg == "-replace") {
			asOptionNameToFind + g_asOptionBuf[0];
			if (g_asOptionBuf.size()==3) {
				asOptionDataIn + g_asOptionBuf[1];
				asOptionDataOut + g_asOptionBuf[2];
			} else {
				asOptionDataIn + (Cstring) "*";
				asOptionDataOut + g_asOptionBuf[1];
			};
		} else if (g_sOptionArg == "-add") {
			asOptionNameToFind + g_asOptionBuf[0];
			asOptionDataIn + cAddOnlyFlag;
			asOptionDataOut + g_asOptionBuf[1];
		} else if (g_sOptionArg == "-list") {
			for (nx = 0; nx < g_asOptionBuf.size(); nx++) {
				asOptionNameToFind + g_asOptionBuf[nx];
				asOptionDataIn + cListOnlyFlag;
				asOptionDataOut + (Cstring) "";
			};
        } else if (g_sOptionArg == "-set") {
            if (g_asOptionBuf[0].contains('%'))
                asVariable + g_asOptionBuf[0];
            else {
                sTemp = "%";
                sTemp += g_asOptionBuf[0];
                sTemp += "%";
                asVariable + sTemp;
            };
            asVariableExtractKeyword + g_asOptionBuf[1];
            asVariableExtractPatternIn + g_asOptionBuf[2];
            asVariableExtractPatternOut + g_asOptionBuf[3];
            asVariableValue + (Cstring) "";
        } else if (g_sOptionArg == "-setadd") {
            if (g_asOptionBuf[0].contains('%'))
                asVariable + g_asOptionBuf[0];
            else {
                sTemp = "%";
                sTemp += g_asOptionBuf[0];
                sTemp += "%";
                asVariable + sTemp;
            };
            asVariableExtractKeyword + sAddString;
            asVariableExtractPatternIn + g_asOptionBuf[1];
            asVariableExtractPatternOut + g_asOptionBuf[2];
            asVariableValue + (Cstring) "";
        } else if (g_sOptionArg == "-setmul") {
            if (g_asOptionBuf[0].contains('%'))
                asVariable + g_asOptionBuf[0];
            else {
                sTemp = "%";
                sTemp += g_asOptionBuf[0];
                sTemp += "%";
                asVariable + sTemp;
            };
            asVariableExtractKeyword + sMultString;
            asVariableExtractPatternIn + g_asOptionBuf[1];
            asVariableExtractPatternOut + g_asOptionBuf[2];
            asVariableValue + (Cstring) "";

		} else if (g_sOptionArg == "-remove") {
			asOptionNameToFind + g_asOptionBuf[0];
			asOptionDataIn + cRemoveOnlyFlag;
			asOptionDataOut + (Cstring) "";
		} else if (g_sOptionArg == "-distance") {
			asOptionNameToFind + (Cstring) "%DET%GONIO_VALUES";
			asOptionDataIn + (Cstring) "*.*.*.*.*.*";
			sTemp = (Cstring) "%1 %2 %3 %4 %5 ";
			sTemp += g_afOptionBuf[0];
			asOptionDataOut + sTemp;
		} else if (g_sOptionArg == "-2theta") {
			asOptionNameToFind + (Cstring) "%DET%GONIO_VALUES";
			asOptionDataIn + (Cstring) "*.*.*.*.*.*";
			sTemp = (Cstring) "%1 ";
			sTemp += g_afOptionBuf[0];
			sTemp += " %3 %4 %5 %6";
			asOptionDataOut + sTemp;
		} else if (g_sOptionArg == "-beamcenter") {
			asOptionNameToFind + (Cstring) "%DET%SPATIAL_BEAM_POSITION";
			sTemp = g_afOptionBuf[0];
			sTemp += " ";
			sTemp += g_afOptionBuf[1];
			asOptionDataIn + (Cstring) "*";
			asOptionDataOut + sTemp;

			asOptionNameToFind + (Cstring) "%DET%SPATIAL_DISTORTION_INFO";
			sTemp += " %3 %4";
			asOptionDataIn + (Cstring) "*.*.*.*";
			asOptionDataOut + sTemp;
		} else if (g_sOptionArg == "-prefix") {
			sPrefix = g_asOptionBuf[0];
		} else
			nStat = 2;
	};
	if (nStat==2)
		return 1;
	nStat = 0;

	if (asOptionNameToFind.size()) {
        int nVar;
		int nFile;
		int nOption;
		int nKeyword;
        double f0,f1;
		bool bFound;
		bool bEditOption;
		Cstring sIn;
		Cstring sOut;
		Cstring sOutName;
		Stringreg oReg;
		Cstring sKeyword;
		Cstring sKeywordPat;
		Cstring sTemp;
		itr<Cstring> asKeywords;
		Cimage_header* poHeader;
		Cimage* poImage;

		for (nFile = 0; (!nStat) && (nFile < asFilesToModify.size());nFile++) {
			poHeader = NULL;
			poImage = NULL;

			printf("Processing file '%s'\n",asFilesToModify[nFile].string());
			if (asFilesToModify[nFile].contains(".head")) {
				poHeader = new Cimage_header(asFilesToModify[nFile]);
				cout << flush;
				fflush(stdout);
				if (!poHeader->bIsAvailable()) {
					printf("ERROR:  Could not open header file '%s'\n",asFilesToModify[nFile].string());
					nStat = 1;
				};
			} else {
				poImage = new Cimage(asFilesToModify[nFile]);
				cout << flush;
				fflush(stdout);
				if (!poImage->bIsAvailable()) {
					printf("ERROR: Could not open file '%s'\n",asFilesToModify[nFile].string());
					nStat = 1;
				} 
			};
		
			if (!nStat) {
				Cimage_header& oHeader = (poHeader)?(*poHeader):poImage->m_oHeader;
				bEditOption = false;


                // Find variable values.
                for (nVar = 0; nVar < asVariable.size(); nVar++) {
                    if (asVariableExtractKeyword[nVar] == sAddString) {
                        sTemp = asVariableExtractPatternIn[nVar];
                        nReplaceVariables(sTemp,asVariable,asVariableValue);
                        if (1==sscanf(sTemp.string(),"%lf",&f0)) {
                            sTemp = asVariableExtractPatternOut[nVar];
                            nReplaceVariables(sTemp,asVariable,asVariableValue);
                            if (1==sscanf(sTemp.string(),"%lf",&f1)) {
                                f0 += f1;
                                sTemp = f0;
                                asVariableValue[nVar] = sTemp;
                            };
                        };                        
                    } else if (asVariableExtractKeyword[nVar] == sMultString) {
                        sTemp = asVariableExtractPatternIn[nVar];
                        nReplaceVariables(sTemp,asVariable,asVariableValue);
                        if (1==sscanf(sTemp.string(),"%lf",&f0)) {
                            sTemp = asVariableExtractPatternOut[nVar];
                            nReplaceVariables(sTemp,asVariable,asVariableValue);
                            if (1==sscanf(sTemp.string(),"%lf",&f1)) {
                                f0 *= f1;
                                sTemp = f0;
                                asVariableValue[nVar] = sTemp;
                            };
                        };                        
                    } else if (oHeader.nGetValue(asVariableExtractKeyword[nVar],&sTemp)) {
                        printf("WARNING:  Could not load variable %s\n",asVariable[nVar].string());
                    } else if (oReg.nParse(sTemp,asVariableExtractPatternIn[nVar])) {
					    printf("WARNING:  Data for variable '%s' does not match pattern '%s'\n",asVariable[nVar].string(),asVariableExtractPatternIn[nVar].string());
                    } else if (oReg.nPrint(asVariableExtractPatternOut[nVar],sTemp)) {
                        printf("WARNING:  Could not output variable '%s' using format '%s'\n",asVariable[nVar].string(),asVariableExtractPatternOut[nVar].string());
                    } else {
                        asVariableValue[nVar] = sTemp;
                    };
                };


				for (nOption = 0; nOption < asOptionNameToFind.size(); nOption++) {
					sKeywordPat = asOptionNameToFind[nOption];
                    
                    nReplaceVariables(sKeywordPat,asVariable,asVariableValue);				
                    -asKeywords;
					oHeader.nFindKeywordMask(sKeywordPat,asKeywords);
					if (!asKeywords.size())
						asKeywords + asOptionNameToFind[nOption];
					for (nKeyword = 0; nKeyword < asKeywords.size(); nKeyword++) {
						sKeyword = asKeywords[nKeyword];
						if (oHeader.nGetValue(sKeyword,&sIn)) {
							sIn = "";
							bFound = false;
						} else
							bFound = true;
                        for (nx = 0; (nx < sIn.length()) && (NULL != strchr("\t\n\r ",sIn.GetAt(nx))); nx++);
                        sIn = sIn.after(nx-1);
                        for (nx = sIn.length()-1; (nx >= 0) && (NULL != strchr("\t\n\r ",sIn.GetAt(nx))); nx++);
                        sIn = sIn.before(nx + 1);

						if (asOptionDataIn[nOption] == cListOnlyFlag) {
							printf("%30s %s\n",sKeyword.string(),(bFound)?sIn.string():"N/A");
							nStat = 1;
						} else 	if (asOptionDataIn[nOption] == cRemoveOnlyFlag) {
							if (!bFound) {
								printf("WARNING:  Option '%s' does not exist.  Cannot remove from header.\n",sKeyword.string());
							} else {
								bEditOption = true;
								oHeader.nDelete(sKeyword);
								printf("INFO: Deleting option '%s' from header.\n",sKeyword.string());
							};
							nStat = 1;
						} else if (asOptionDataIn[nOption] == cAddOnlyFlag) {
							if (bFound) {
								printf("WARNING:  Option '%s' already exists. Please use -replace if you want to modify.\n",sKeyword.string());
								nStat = 1;
							};
							sOut = asOptionDataOut[nOption];
						} else if (oReg.nParse(sIn,asOptionDataIn[nOption])) {
							printf("WARNING:  Data for option '%s' does not match pattern '%s'\n",sKeyword.string(),asOptionDataIn[nOption].string());
							nStat = 1;
							break;
						} else {
                            sTemp = asOptionDataOut[nOption];
                            nReplaceVariables(sTemp,asVariable,asVariableValue);				
							if (oReg.nPrint(sTemp,sOut)) {
								printf("WARNING:  Could not construct output for option '%s'\n",sKeyword.string());
								nStat = 1;
								break;
							};
						};
						if (nStat) {
							nStat = 0;
							continue;
						} else if (oHeader.nReplaceValue(sKeyword,sOut)) {
							printf("WARNING:  Could not update option '%s'\n",sKeyword.string());
							nStat = 1;
						} else {
							bEditOption = true;
							printf("INFO: Updating option '%s'.\n",sKeyword.string());
							printf("      Was: %s\n",(bFound)?sIn.string():"N/A");
							printf("      Now: %s\n",sOut.string());
						};
					};
				};
			};
			if ((!nStat) && (bEditOption)) {
				if (sPrefix)
					sOutName = sPrefix;
				else
					sOutName = "";
				sOutName += sFileGetBasename(asFilesToModify[nFile]);
				sOutName = sFileBuildName(sFileGetDirectory(asFilesToModify[nFile]),sOutName);
                if (poImage) {
                    poImage->m_bWriteRAXIS = false;
                    poImage->nWrite(sOutName);
                };
				if (poHeader)
					poHeader->nWrite(sOutName);
			};
            if (poImage) 
				delete poImage;

			if (poHeader)
				delete poHeader;
		};
	} else {
		printf("INFO:  Nothing to do!\n");
	};

	return nStat;
    
};