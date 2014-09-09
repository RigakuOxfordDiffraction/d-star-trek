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
// CCommandLine.cc        Initial author: Thaddeus Niemeyer   Sep-2002
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

#include "CCommandLine.h"
#include "dtreksys.h"

char*                   g_pcHelp = NULL;
bool                    m_bCommandLinePassive = false;
itr<int>                g_anOptionBuf;
itr<double>             g_afOptionBuf;
std::vector<Cstring>    g_asOptionBuf;
Cstring                 g_sOptionArg;
char**                  g_ppcOptionArgs = NULL;
int                     g_nArgc = 0;
int                     g_nAllocArgc = 0;

int nFreeCommandLine()
{
    g_pcHelp = NULL;
    m_bCommandLinePassive = 0;
    g_anOptionBuf.clear();
    g_afOptionBuf.clear();
    g_asOptionBuf.clear();
    g_sOptionArg = "";
    
    if( g_ppcOptionArgs )
    {
        for (int nx=0; nx < g_nAllocArgc; nx++)
            delete[] g_ppcOptionArgs[nx];
        
        delete[] g_ppcOptionArgs;
        
        g_ppcOptionArgs = NULL;
    }
    
    g_nArgc = 0;
    g_nAllocArgc = 0;

    return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int nReadCommandLine(int& nArg, unsigned int argc, char* argv[], const char* apcOptions[], int nArgsSkip)
{
	Cstring         sTemp("");
	Cstring         sTemp2("");
    
    int             nArgStart = 0;
    int             nArgEnd = 0;
    int             nArgCount = 0;    // For each argument, the list of following arguments are indexed by these.
    
    int             nx = 0, ny = 0;
    double          f0 = 0.0, f1 = 0.0;
	
    int             nOptionLang = 0;
	
	g_anOptionBuf.clear();
	g_afOptionBuf.clear();
	g_asOptionBuf.clear();
	
	//argc--; 
	//argv++; // Skip program name
	
    argc -= nArgsSkip; 
	argv += nArgsSkip;

    if( nArg < 0 )
    {
        nFreeCommandLine();
        
        m_bCommandLinePassive = (nArg == -2);
            
		nArg = 0;

		// If there are no arguments at all, print the help string.
		if( argc == 0 )
        {
            Cstring         sTemp("");
			
            nPrintItemSyntax(apcOptions, sTemp);
            
            if( !m_bCommandLinePassive && g_pcHelp )
			    printf(g_pcHelp);
			
            return 1;
		}
        ////////////////////////////////////////////////////////////////////
        
        if( !m_bCommandLinePassive && argv[0][0] == '@' )
        {
            // Read arguments in from a file since this option was requested.
            
            Cstring sFile = &argv[0][1];
            FILE* pFArgs;
            char* pcTok;
            char* pcBuf;
            int nPass;
            int nArgCount;
            int nComment;
            for (nPass = 0; nPass < 2; nPass++) {
                nArgCount = 0;
                pFArgs = fopen(sFile.string(),"rt");
                if( !pFArgs) {
                    printf("ERROR:  Could not open argument file '%s'\n",sFile.string());
                    return 2;
                }
                nx = lFileGetSize(sFile);
                pcBuf = new char[nx+1];
                ny = fread(pcBuf,1,nx,pFArgs);
                pcBuf[ny] = 0;
                for (nx=0,nComment = 0;nx<ny;nx++) {
                    if( pcBuf[nx] == '#') 
                        nComment = 1;
                    else if( (pcBuf[nx] == '\n') || (pcBuf[nx] == '\r'))
                        nComment = 0;
                    if( nComment)
                        pcBuf[nx] = ' ';
                }
                if( nPass) {
                    g_ppcOptionArgs[nArgCount] = new char[30];
                    strcpy(g_ppcOptionArgs[nArgCount],"-START");
                    nArgCount++;
                }

                for (pcTok = strtok(pcBuf,"\t\n\r ");pcTok;pcTok = strtok(NULL,"\t\n\r ")) {
                    if( nPass) {
                        g_ppcOptionArgs[nArgCount] = new char[strlen(pcTok)+1];
                        strcpy(g_ppcOptionArgs[nArgCount],pcTok);
                    }
                    nArgCount++;
                }
                delete[] pcBuf;


                if( !nPass) {
                    g_nArgc = nArgCount + 1; 
                    g_nAllocArgc = nArgCount+10;
                    g_ppcOptionArgs = new char*[g_nAllocArgc];
                }
                fclose(pFArgs);                
            }
            for (;nArgCount<g_nAllocArgc;nArgCount++) {
                g_ppcOptionArgs[nArgCount] = new char[10];
                strcpy(g_ppcOptionArgs[nArgCount],"");
            }

        }
        else
        {
            g_nArgc = argc + 1;
            g_nAllocArgc = g_nArgc + 10;
            
            g_ppcOptionArgs = new char*[g_nAllocArgc];
            
            g_ppcOptionArgs[0] = new char[30];
            
            strcpy(g_ppcOptionArgs[0], "-START");
            
            for (nx=0; nx < argc; nx++)
            {
                g_ppcOptionArgs[nx+1] = new char[strlen(argv[nx])+1];
                
                strcpy(g_ppcOptionArgs[nx+1],argv[nx]);
            }
            
            for(nx++; nx < g_nAllocArgc; nx++)
            {
                g_ppcOptionArgs[nx] = new char[10];
                
                strcpy(g_ppcOptionArgs[nx], "");
            }
        }
        
        if( !m_bCommandLinePassive )
        {
            printf("Command line:\n ");
        
            for (nx=0,ny=0; nx < g_nArgc; nx++)
            {
                if( (strcmp(g_ppcOptionArgs[nx],"-START")) && (strcmp(g_ppcOptionArgs[nx],"--")))
                {
                    if( ny + strlen(g_ppcOptionArgs[nx])>=72)
                    {
                        printf(" \\\n   ");
                        ny = 0;
                    }
                    
                    printf("%s ",g_ppcOptionArgs[nx]);
                    
                    ny += strlen(g_ppcOptionArgs[nx]) + 1;
                }
            }		
            
            printf("\n");
        }
	}
	
    if( nArg == g_nArgc )
		return 1;
	
	if( g_ppcOptionArgs[nArg][0] != '-' && g_ppcOptionArgs[nArg][0] != '+' )
    {
        if( !m_bCommandLinePassive )
		    printf("ERROR: '-' or '+' must preceed option %s\n", g_ppcOptionArgs[nArg]);
		
        return 2;
	}
	
    while( nArg < g_nArgc && !strcmp(g_ppcOptionArgs[nArg],"--") )
        nArg++;

    if( nArg == g_nArgc )
		return 1;

	g_sOptionArg = g_ppcOptionArgs[nArg];
	
    bool        bBackTrack = false;
    bool        bIsMinus = false;
	
	nArg++;
    nArgStart = nArg;
    nArgEnd = nArgStart - 1;
    
    for(; nArg < g_nArgc; nArg++)
    {
		bIsMinus = false;
		
		//if( ((g_ppcOptionArgs[nArg][0]=='-') || (g_ppcOptionArgs[nArg][0]=='+')) && 
        //    (strspn(g_ppcOptionArgs[nArg],"-0123456789.+:") != strlen(g_ppcOptionArgs[nArg])))
        
		sTemp = g_ppcOptionArgs[nArg];

        if( bIsDTCommandLineOptionString(sTemp.string()) )
        {
            // Make sure it's a bona-fide option in the loop below.
			
            bIsMinus = false;
            
            for(nOptionLang = 0; NULL != apcOptions[nOptionLang]; nOptionLang++)
            {
                if( apcOptions[nOptionLang][0] == '-' && sTemp == apcOptions[nOptionLang] ) 
                {
                    bIsMinus = true;
                    break;
                }
            }
            
            if( sTemp.GetAt(0) == '-' && m_bCommandLinePassive )
                bIsMinus = true;
		} 
        else if( 0 == strcmp(g_ppcOptionArgs[nArg], "--") )
			bIsMinus = true;
		
		if( bIsMinus)                  
			break;
        else if( g_sOptionArg.find('-') == 0 ||  g_sOptionArg.find('+') == 0 )
        {
			nArgEnd = nArg;
		} 
        else
        {
            if( !m_bCommandLinePassive)
                printf("ERROR:  Expecting option but got '%s'\n", g_sOptionArg.string());
			
            return 2;
		}
	}
	
    nArgCount = nArgEnd - nArgStart + 1;
    
        
	Cstring sPrefix("");
	Cstring sPostfix("");
	bool	bOptional = false;
	int     nRepeats = 0;
	int		nArgRead = 0;
	int		nArgReadStart = 0;                  // Used mostly for printing.
    bool    bNameFound = 0;                     // Was the name found?
    bool    bParseError = 0;                    // A flag used to ignore the error in this option in favor of the next one.  
                                                // This allows for overloads of option names, provided they are placed consecutively on the command line.
    bool    bNonPostfixEncountered = false;     // Was a postfixed sub-option encountered.
    bool    bPrintedProcOptionString = false;   // Did we print information on this option yet?
	
    bNameFound = false;
    bPrintedProcOptionString = false;

    for(nOptionLang = 0; NULL != apcOptions[nOptionLang]; nOptionLang++)
    {
        g_anOptionBuf.clear();
        g_afOptionBuf.clear();
        g_asOptionBuf.clear();

		// Step to next option set which should start with a '-'
		while( NULL != apcOptions[nOptionLang] && 
			   ((Cstring)apcOptions[nOptionLang]).find('-') != 0 &&
               ((Cstring)apcOptions[nOptionLang]).find('+') != 0 )
        {	
               nOptionLang++;
        }
        
        if( NULL == apcOptions[nOptionLang] )
            break;

        // Print information if we found this string.
        if( g_sOptionArg == apcOptions[nOptionLang] && !bPrintedProcOptionString )
        {
            if( (nArgCount == 0) || (g_ppcOptionArgs[nArgStart][0]))
            {        
                printf("Processing Option: %s ",(nArgStart==1)?"":g_sOptionArg.string());
                
                for(nx=nArgStart; nx<=nArgEnd; nx++)
                {
                    printf("%s ",g_ppcOptionArgs[nx]);
                }
                printf("\n");
                
                bPrintedProcOptionString = true;
            }
        }

		if( g_sOptionArg == apcOptions[nOptionLang])
        {
			nOptionLang++;
			nArgRead = 0;
            bNameFound = 1;
			
            // For multiple option commands, we might have already
			// gone through and set some string entries to "".
			
			while( !g_ppcOptionArgs[nArgStart + nArgRead][0] && nArgRead < nArgCount )
            {
				nArgRead++;
            }
            
            nArgReadStart = nArgRead;
            
            bParseError = false;
            
            bNonPostfixEncountered = false;

            for(nx = nOptionLang + 1; apcOptions[nx]             && 
                                      apcOptions[nOptionLang][0] &&
				                      apcOptions[nx][0] != '-'   &&
                                      apcOptions[nx][0] != '+';                        nx++) 
            {
                sTemp = apcOptions[nx];
                
                if( !sTemp.after('.').length())
                    bNonPostfixEncountered = true;
            }
            ///////////////////////////////////////////////////////////////////////////////////

			do
            {
				sTemp = apcOptions[nOptionLang];
				sTemp = sTemp.before('[') + sTemp.after(']');
				sPrefix = sTemp.before(".");
				sPostfix = sTemp.after(".");
                sPostfix.downcase();
			
                if( sPrefix.contains('*'))
                {
					bOptional = true;
					sPrefix = sPrefix.after('*');
				} 
                else
					bOptional = false;


				if( (sPrefix.length() == 2) && (sPrefix.GetAt(1)>='0') && (sPrefix.GetAt(1)<='9'))
					nRepeats = sPrefix.GetAt(1) - '0';
				else if( sPrefix == "n")
                    nRepeats = 0;
                else
					nRepeats = 1;

                sTemp = (Cstring) g_ppcOptionArgs[nArgStart + nArgRead];
                sTemp.downcase();
				
                if( (sPostfix.length()) && 
					(nArgRead + nRepeats < nArgCount) &&
					(sPostfix == sTemp))
                {
					nArgRead++;
				} 
                else if( (sPostfix.length()) && (bOptional))
					continue;
                else if( (sPostfix.length()) && (!bOptional))
                {
                    bParseError = true;
                    break;
                }


				if( sPrefix == "n")
                {
				} 
                else if( (sPrefix.GetAt(0) == 's') && 
					(nArgRead + nRepeats <= nArgCount))
                {
					for (nx=0;nx<nRepeats;nx++)
                    {
                        g_asOptionBuf.push_back((Cstring) g_ppcOptionArgs[nArgStart+nArgRead]);
                        nArgRead++;
					}
				} 
                else if( (sPrefix.GetAt(0) == 'd') && 
					     (nArgRead + nRepeats <= nArgCount))
                {
					for (nx = 0; nx< nRepeats; nx++)
                    {
						if( (strlen(g_ppcOptionArgs[nArgStart+nArgRead])==strspn(g_ppcOptionArgs[nArgStart+nArgRead],"+-0123456789")) &&
							(1 == sscanf(g_ppcOptionArgs[nArgStart+nArgRead],"%d",&ny)))
                        {
							g_anOptionBuf + ny;
							nArgRead++;
						} 
                        else
							break;
					}
					
                    if( (nx < nRepeats) && (!bOptional)) 
                    {
						bParseError = true;
						break;
					}
				} 
                else if( (sPrefix.GetAt(0) == 'f') && (nArgRead + nRepeats <= nArgCount) )
                {
					for (nx = 0; nx < nRepeats; nx++)
                    {
						if( 1 == sscanf(g_ppcOptionArgs[nArgStart+nArgRead],"%lf",&f0))
                        {
							g_afOptionBuf + f0;
							nArgRead++;
						} 
                        else
							break;
					}
					
                    if( (nx < nRepeats) && (!bOptional))
                    {
						bParseError = true;
						break;
					}
                } 
                else if( (sPrefix.length()==2) && (sPrefix.GetAt(0) == 'r'))
                {
                        if( nArgRead < nArgCount) 
                        {
                            nOptionLang -= (1 + (sPrefix.GetAt(1)-'0'));
                            continue;                    
                        }
				} 
                else if( !bOptional) 
                {
                    bParseError = true;
                    break;
				} 
                else
					continue;
				
				if( sPostfix.length())
                {
					printf("Processing sub option: ");
				
                    for (nx=nArgReadStart;nx<nArgRead;nx++)
						printf("%s ",g_ppcOptionArgs[nArgStart + nx]);
					
                    printf("\n");
					
					g_sOptionArg += "-";
					g_sOptionArg += sPostfix;
					
                    for (nx=0;nx<nArgRead;nx++)
						strcpy(g_ppcOptionArgs[nArgStart + nx],"");
                    
                    nArg = nArgStart - 1;

                    if( (nArgRead == nArgCount) && (!bNonPostfixEncountered))
                    {
						nArg = nArgEnd + 1;
                    } 
                    else
                        nArg = nArgStart - 1;

					break;					
				} else
					nArg = nArgEnd + 1;
				
			} while (nOptionLang++,
				((apcOptions[nOptionLang]) &&
				(apcOptions[nOptionLang][0]) &&
				(apcOptions[nOptionLang][0] != '-') &&
                (apcOptions[nOptionLang][0] != '+')));

            
            if( bParseError)
            {
                nOptionLang--;
                continue;
            }
			
            if( (nArgRead != nArgCount) && (nArg != nArgStart - 1))
            {
                nOptionLang--;
                continue;
			} 
            else
				return 0;
		}
	}

    // Should not get here unless there was an error.
    if( (m_bCommandLinePassive) && (!bPrintedProcOptionString) && (!bNameFound)) 
    {
        g_sOptionArg = "";
        
        return 0;
    } 
    else 
    {
        if( bNameFound)
        {
            printf("ERROR:  Could not parse remainder of option '%s':  ",
                g_sOptionArg.string());
        
            nArgRead = 0;
            
            while (nArgRead < nArgCount)
            {
                if( g_ppcOptionArgs[nArgStart + nArgRead][0])
                    printf("%s ",g_ppcOptionArgs[nArgStart + nArgRead]);
                
                nArgRead++;
            }
            printf("\n");
            
            if( nPrintErrorHelp(g_sOptionArg))
            {
                nPrintItemSyntax(apcOptions,g_sOptionArg);
            }
        } 
        else
            printf("ERROR: Unrecognized option '%s'\n",g_sOptionArg.string());
    }
	
	return 2;
		
}
///////////////////////////////////////////////////////////////////////////////
int nPrintErrorHelp(Cstring& sOption)
{
	Cstring sHelp;
	Cstring sFind;
	Cstring sTemp;
    bool bNextOption;
	int nx;
    int nStat;

	if( (!g_pcHelp) || (!g_pcHelp[0]))
		return 0;

	sHelp = g_pcHelp;
	sFind = "\n";
	sFind += sOption;
    nx = 0;
    nStat = 1;
    while (nx != -1) {
        nx= sHelp.find(sFind);
            
        if( nx != -1) {
            printf("\n");
            sHelp = sHelp.after(nx);
            do {
                sTemp = sHelp.before('\n');
                printf("%s\n",sTemp.string());
                sHelp = sHelp.after('\n');
                bNextOption = ((sHelp.length() == 0) || (sHelp.GetAt(0) == '-'));                   
            } while (!bNextOption);
            nStat = 0;
        }
    }


	return nStat;
}

int nPrintItemSyntax(const char* apcOptions[],Cstring& sItemName) {
	Cstring sName;
	Cstring sOption;
	Cstring sPrefix;
	Cstring sPostfix;
	int nx;
	int nCount;
	int nFormatIndent;
	int nFormatPrinted;
	int nRepeatCount;
	int nCountf,nCountd,nCounts;
	int nRepeats;
	bool bOptional;
	bool bStartOption;

    if( sItemName.length() == 0)
	    printf("\nSummary of command line options.\n");
	nFormatIndent = 7;
	for (nCount = 0; apcOptions[nCount];nCount++) {
		if( apcOptions[nCount][0] == '-') 
			nFormatIndent = max(nFormatIndent,(int) strlen(apcOptions[nCount]));
	}

	for (nCount = 0; apcOptions[nCount];) {
		if( !strcmp(apcOptions[nCount],"-START")) {
			
			printf("Format: ");
			nFormatPrinted = 8;
			bStartOption = true;
		} else if( (sItemName.length()==0) || (sItemName == apcOptions[nCount])) {
			printf("%s",apcOptions[nCount]);
			nFormatPrinted = strlen(apcOptions[nCount]);
			bStartOption = false;
        } else {
           nCount++;
           while ((apcOptions[nCount]) && (apcOptions[nCount][0]!='-'))
               nCount++;
           continue;
        }

            
		
		nCount++;
		nCountf = nCountd = nCounts = 0;
		nRepeatCount = 0;
		while (apcOptions[nCount]) {
			sOption = apcOptions[nCount];
			if( sOption.GetAt(0) == '-')
				break;
			sName = sOption.after('[');
			sName = sName.before(']');
			sOption = sOption.before('[') + sOption.after(']');
			sPrefix = sOption.before(".");
			sPostfix = sOption.after(".");
			if( sPrefix.contains('*')) {
				bOptional = true;
				nCountf = nCountd = nCounts = 0;
				sPrefix = sPrefix.after('*');
			} else
				bOptional = false;

			// Formated printing.
			if( nFormatPrinted + sPostfix.length() + sName.length() + 2>65) {
				nFormatPrinted = 0;
				printf("\n");
			}
			while (nFormatPrinted < nFormatIndent) {
				printf(" ");
				nFormatPrinted++;
			}
			printf(" ");

			nFormatPrinted += sPostfix.length() + sName.length() + 2;

			if( (sPrefix.length() == 2) && (sPrefix.GetAt(1)>='0') && (sPrefix.GetAt(1)<='9'))
				nRepeats = sPrefix.GetAt(1) - '0';
			else
				nRepeats = 1;


			Cstring     sTemp("");
            
            if( !sName.length())
            {
				sName = "";
				
                for(nx=0; nx < nRepeats; nx++) 
                {
					if( nx >= 1 )
						sName += ' ';
					
                    if( sPrefix.contains('d') )
                    {
                        sTemp = "nArg";
                        sTemp += (++nCountd);
                        
                        sName += sTemp;
                    }
                    else if( sPrefix.contains('f') )
                    {
                        sTemp = "fArg";
                        sTemp += (++nCountf);
                        
                        sName += sTemp;
                    }
                    else if( sPrefix.contains('s'))
                    {
                        sTemp = "sArg";
                        sTemp += (++nCounts);
                        
                        sName += sTemp;
                    }
                    else
						break;
				}
			}
			
            if( sPrefix.contains('r')) {
				if( nRepeatCount == 1) {
					printf("... }");
					nCount++;
					break;
				}
				printf("{");
				nCount -= (sPrefix.GetAt(1) - '0');
				nRepeatCount++;
				continue;
			}

			if( bOptional)
				printf("[");
			if( (sPostfix.length()) && (sName.length())) {
				printf("%s %s",sPostfix.string(),sName.string());
				
			} else {
				printf("%s",sPostfix.length()?sPostfix.string():sName.string());
			}
			if( bOptional)
				printf("]");
			nCount++;
		}		
		if( bStartOption) 
			printf(" [options ...]\n\nOption:\n");
		else
			printf("\n");
	}
	printf("\n");
	return 0;
}
///////////////////////////////////////////////////////////////////////////////////
bool bIsDTCommandLineOptionString(char* pcString)
{
    int     nLength = strlen(pcString);
    
    if( nLength < 2 )
        return false;              // there must be a dash and at least one letter

    if( pcString[0] != '-' && pcString[0] != '+')
        return false;              // must start with a dash or plus sign 
    
    if( strspn(pcString, "+-0123456789.:,; ") == nLength )
        return false;             // starts with a dash or plus sign, but the rest are not alphabetical characters

    return true;
}
