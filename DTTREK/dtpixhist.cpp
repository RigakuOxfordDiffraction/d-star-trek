//
// dtpixhist.cpp     Initial author: RB           27-JULY-2004

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

#ifdef SSI_PC
#define main   dtete_main
#endif

#include "Cstring.h"
#include "Cimage_header.h"
#include "Cscan.h"

#include "DTMain.h"

using std::cout;
using std::endl;
using std::flush;
using std::cerr;

#define MAX_NUMBER_OF_BITS      20


typedef struct _tagIMAGE_SEQUENCE
{
    int     m_nSeq1;
    int     m_nSeq2;
    
    _tagIMAGE_SEQUENCE(int nSeq1, int nSeq2){m_nSeq1=nSeq1; m_nSeq2=nSeq2;}
}IMAGE_SEQUENCE;

static void vError(const int nErrorNum, const Cstring& sMessage);

static bool bSetCustomLinearBins(int nMinPixValue, 
                                 int nMaxPixValue, 
                                 int nBinWidth,
                                 std::vector<_tagPIXEL_INTENSITY_BIN>& vecBins);

static bool bSetCustomLogarithmicBins(int nMinPixValue, 
                                      int nMaxPixValue, 
                                      int nPowerBaseForBinning,
                                      std::vector<_tagPIXEL_INTENSITY_BIN>& vecBins);

#define   DTREK_PIXHIST_OCPIH_CCD            0x0001
#define   DTREK_PIXHIST_OCPIH_RAXIS_II       0x0002
#define   DTREK_PIXHIST_OCPIH_NO_SKIP_BINS   0x0004
static bool bOutputCustomPixelIntensityHistogram(std::vector<_tagPIXEL_INTENSITY_BIN>& vecBins, 
                                                 Cstring& strOutputFilePath,
                                                 DTREK_WORD wCtrl);

static bool bOutputPixelIntensityBitTestResults(std::vector<int>& vecBitBins);

static bool bIsBinMustBeEmptyOnIP(int nBegin, int nEnd, bool bRAXIS_II);

static bool bIsRAXIS_II_image(Cimage_header* poHeader);
static bool bIsCCD_image(Cimage_header* poHeader);

///////////////////////////////////////////////////////////////////////////////////////////////
int main(int   argc,     // Command line argument count
         char *argv[])   // Pointers to command line args
{
    vDtrekSetModuleName("dtpixhist");
    vPrintCopyrightInfo();
    //cout << "\ndtpixhist:  Copyright (c) 2006 Rigaku\n";
    //cout << D_K_DTREKVersion << endl;

    // Copy command line to output log
    cout << "\nCommand line:\n" << sGetCommandLine(argc, argv, 71) << endl << endl << flush;

    // check for help request
    argc--; argv++;

    if( 0 == argc )
    {
        DTREK_ERROR(1, "ERROR: dtpixhist - input image/header file is not specified.");
    }

    if( 1 <= argc && (0==strcmp(argv[0],"-h") || 0==strcmp(argv[0],"-help")) )
    {
        DTREK_ERROR(0, "");  // display help
    }
    
    argc++; argv--;  // restore original values
    /////////////////////////////////////////////////////////////////////////////////////////    
    
    Cstring                                 strHeaderPath("");
    std::vector<_tagIMAGE_SEQUENCE>         vecImageSequences;
    Cstring                                 strOutputFilePath("");
    

    Cstring                                 strInput("");

    int             nSeq1 = 0;
    int             nSeq2 = 0;
    
    int             nMinPixelValue = 0;
    int             nMaxPixelValue = (int)pow(2.0, MAX_NUMBER_OF_BITS);
    int             nBinCount = -1;
    int             nBinWidth = -1;
    bool            bPowerBins = false;
    int             nPowerBaseForBinning = -1;

    bool            bCustomHistogram = true;
    bool            bNoBinSkip = false;


    int             ii = 1;   // because ii=0 is the executable path
    while( ii < argc )
    {
        strInput = argv[ii];
        
        if( ii == 1 && ii < argc )
        {
            strHeaderPath = strInput;
            ii++;
            continue;
        }
        else if( "-seq" == strInput && ii < argc - 1 )
        {
            // get the first image number
            nSeq1 = atoi(argv[ii+1]);

            if( ii < argc - 2 ) // determine if there is a second argument to the -seq
            {
                strInput = argv[ii+2];
                if( strInput.bIsNumeric() )
                {
                    nSeq2 = atoi(strInput);
                    ii += 3;
                }
                else
                {
                    nSeq2 = nSeq1;
                    ii += 2;
                }
            }
            else
            {
                nSeq2 = nSeq1;
                ii += 2;
            }
            
            vecImageSequences.push_back(_tagIMAGE_SEQUENCE(nSeq1, nSeq2));
            continue;
        }
        else if( "-min" == strInput && ii < argc - 1  )
        {
            nMinPixelValue =  atoi(argv[ii+1]);
            ii += 2;
            continue;
        }
        else if( "-max" == strInput && ii < argc - 1  )
        {
            nMaxPixelValue =  atoi(argv[ii+1]);
            ii += 2;
            continue;
        }
        else if( "-bins" == strInput && ii < argc - 1  )
        {
            nBinCount =  atoi(argv[ii+1]);
            ii += 2;
            continue;
        }
        else if( "-binwidth" == strInput && ii < argc - 1  )
        {
            strInput = argv[ii+1];
            
            if( strInput == "pow2" )
            {
                bPowerBins = true;
                nPowerBaseForBinning = 2;
            }
            else if( strInput == "pow10" )
            {
                bPowerBins = true;
                nPowerBaseForBinning = 10;
            }
            else
                nBinWidth =  atoi(argv[ii+1]);
            
            ii += 2;
            continue;
        }
        else if( "-out" == strInput && ii < argc - 1  )
        {
            strOutputFilePath = argv[ii+1];
            ii += 2;
            continue;
        }
        else if( "-noskipbins" == strInput  )
        {
            bNoBinSkip = true;
            ii += 1;
            continue;
        }
        
        ii++;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////
    // Error checking
    ////////////////////////
    if( vecImageSequences.size() == 0 )
    {
        printf("\nImage sequence not specified or invalid.\n");
        return 1;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Set bins
    int      nBinMaxCurrent = 0;       
    int      nBinMinCurrent = 0;       
    
    std::vector<_tagPIXEL_INTENSITY_BIN>   vecBins;            
    
    if( bPowerBins && nPowerBaseForBinning > 1 )
    {
        if( nMinPixelValue <= 0 )                           // Safety check
            nMinPixelValue = 1;
        
        if( nMaxPixelValue <= 1 )
            nMaxPixelValue = (int)pow(2.0, MAX_NUMBER_OF_BITS);  // restore default, if customer entered a wrong value

        bSetCustomLogarithmicBins(nMinPixelValue, nMaxPixelValue, nPowerBaseForBinning, vecBins);
    }
    else if( (nBinWidth > 0 || nBinCount > 0) && nMaxPixelValue - nMinPixelValue > 0 )    
    {    
        if( nBinWidth < 1 && nBinCount > 0 ) // if the bin width is NOT set, but bin number IS set. 
        {
            nBinWidth = (int)((double)(nMaxPixelValue - nMinPixelValue) / (double)nBinCount);
            printf("\nCalculated integer bin width: %d\n\n", nBinWidth);
            
            if( nBinWidth < 1 )
            {
                nBinWidth = 1;
                printf("\nSetting bin width to 1\n\n");
            }
        }
        
        bSetCustomLinearBins(nMinPixelValue, nMaxPixelValue, nBinWidth, vecBins);
    }
    else if( nMaxPixelValue - nMinPixelValue > 0 )  // just force it to behave as if a pow2 option is given!
    {
        if( nMinPixelValue <= 0 )                           // Safety check
            nMinPixelValue = 1;

        if( nMaxPixelValue <= 1 )
            nMaxPixelValue = (int)pow(2.0, MAX_NUMBER_OF_BITS);  // restore default, if user entered a wrong value

        nPowerBaseForBinning = 2;
        bPowerBins = true;
        bSetCustomLogarithmicBins(nMinPixelValue, nMaxPixelValue, nPowerBaseForBinning, vecBins);
    }
    else
    {
        printf("\nCustom histogram will not be created\n");

        if( nMaxPixelValue - nMinPixelValue <= 0 )
            printf("\nMax value %d cannot be equal or smaller than Min value %d\n\n", nMaxPixelValue, nMinPixelValue);

        bCustomHistogram = false;
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Set up the MAX_NUMBER_OF_BITS bin array to check for the presence of either of the MAX_NUMBER_OF_BITS bits
    std::vector<int>   vecBitBins;            
    int  nTemp = 0;
    
    vecBitBins.push_back(nTemp); // to count values with no bits set, i.e. zero values
    for(ii=0; ii < MAX_NUMBER_OF_BITS; ii++)
        vecBitBins.push_back(nTemp);
    vecBitBins.push_back(nTemp); // to count the rest of the values
    //////////////////////////////////////////////////////////////////////

    
    // Read images and build the histogram(s)

    Cimage_header*          poHeader = new Cimage_header(strHeaderPath);
    Cscan*                  poScan = new Cscan(*poHeader);
    Cimage*                 poImage = new Cimage();

    bool                    bRAXIS_II = bIsRAXIS_II_image(poHeader);
    bool                    bCCD      = bIsCCD_image(poHeader);

    DTREK_WORD              wCtrl = DTREK_IMAGE_BIH_BITTEST;
    if( bCustomHistogram )
        wCtrl |= DTREK_IMAGE_BIH_CUSTOM;

    if( bRAXIS_II )
        wCtrl |= DTREK_IMAGE_BIH_RAXIS_II;
    else if( bCCD )
        wCtrl |= DTREK_IMAGE_BIH_CCD;

    /////////////////////////////////////////////////////////////////////
    int         nNumberOfPixValuesInOverlapRange_PMT1 = 0;
    int         nNumberOfPixValuesInOverlapRange_PMT2 = 0;

    int         nTemp_PMT1 = 0;
    int         nTemp_PMT2 = 0;

    for(ii=0; ii < vecImageSequences.size(); ii++)
    {
        for(int jj=vecImageSequences[ii].m_nSeq1; jj<=vecImageSequences[ii].m_nSeq2; jj++)
        {
            poScan->vSetSeqNum(jj);
            poScan->nGetImage(poImage);
            
            
            nTemp_PMT1 = 0;
            nTemp_PMT2 = 0;
            
            poImage->vBuildIntensityHistogram(vecBins, 
                                              vecBitBins, 
                                              nTemp_PMT1, 
                                              nTemp_PMT2, 
                                              wCtrl);
            
            nNumberOfPixValuesInOverlapRange_PMT1 += nTemp_PMT1;
            nNumberOfPixValuesInOverlapRange_PMT2 += nTemp_PMT2;
        }
    }
    
    if( bCustomHistogram )
    {
        DTREK_WORD wCtrl = 0U;
        
        if( bRAXIS_II )
            wCtrl |= DTREK_PIXHIST_OCPIH_RAXIS_II;
        else if( bCCD )
            wCtrl |= DTREK_PIXHIST_OCPIH_CCD;

        if( bNoBinSkip )
            wCtrl |= DTREK_PIXHIST_OCPIH_NO_SKIP_BINS;
        
        bOutputCustomPixelIntensityHistogram(vecBins, strOutputFilePath, wCtrl);
    }

    bOutputPixelIntensityBitTestResults(vecBitBins);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    if( !bCCD )
    {
        int     nBitShiftCoeff = bRAXIS_II ? 8 : 32; 
        
        printf("\nNumber of pixel intensities in the range %d - %d that ARE NOT multiples of %d : %d",
                PMT_OVERLAP_BEGIN,
                PMT_OVERLAP_END,
                nBitShiftCoeff,
                nNumberOfPixValuesInOverlapRange_PMT1);
             
        printf("\nNumber of pixel intensities in the range %d - %d that   ARE   multiples of %d : %d\n",
                PMT_OVERLAP_BEGIN,
                PMT_OVERLAP_END,
                nBitShiftCoeff,
                nNumberOfPixValuesInOverlapRange_PMT2);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////

    // clean-up
    delete poHeader;
    delete poScan;
    delete poImage;

    return 0;
}

bool bOutputCustomPixelIntensityHistogram(std::vector<_tagPIXEL_INTENSITY_BIN>& vecBins, 
                                          Cstring& strOutputFilePath,
                                          DTREK_WORD wCtrl)
{
    FILE*       pfOut = NULL;
    if( "" != strOutputFilePath )
    {
        pfOut = fopen(strOutputFilePath.string(),"wt");
        if( !pfOut )
        {
            printf("\n\nError: cannot open file: %s", strOutputFilePath.string());
        }
        
        return false;
    }
    //////////////////////////////////////////////////////////////////////////////
    
    bool    bNoBinSkip = 0U != (wCtrl & DTREK_PIXHIST_OCPIH_NO_SKIP_BINS);
    bool    bRAXIS_II  = 0U != (wCtrl & DTREK_PIXHIST_OCPIH_RAXIS_II);  
    bool    bCCD       = 0U != (wCtrl & DTREK_PIXHIST_OCPIH_CCD);  

    int         nTotalCount = 0;
    int         ii=0;
    for(ii=0; ii < vecBins.size(); ii++)
         nTotalCount += vecBins[ii].m_nCount;

    double      dFraction = 0.0;

    printf("\nNumber of pixels with intensity values within specified intensity ranges\n");
    printf("-------------------------------------------------------------------\n");
    printf("     >=Min\t     <=Max\t          Num pixels\t %% of total\n");
    printf("-------------------------------------------------------------------\n");
    
    printf("          \t%10d\t%20d\t%11.2f\n", vecBins[0].m_nEnd - 1 , vecBins[0].m_nCount, 
                                             (double)vecBins[0].m_nCount / max(1.0, (double)nTotalCount)*100.0 );
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool        bBinMustBeEmptyBecauseOfBitShift = false;
    bool        bBinIsInPMTOverlapRange = false;
    for(ii=1; ii < vecBins.size()-1; ii++)
    {
        if( !bCCD )
        {
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            bBinMustBeEmptyBecauseOfBitShift = bIsBinMustBeEmptyOnIP(vecBins[ii].m_nBegin, vecBins[ii].m_nEnd-1, bRAXIS_II);
            if( vecBins[ii].m_nCount != 0 && bBinMustBeEmptyBecauseOfBitShift ) 
            {
                printf("\n\nError: pixel intensities found between %d and %d inclusive\n\n", vecBins[ii].m_nBegin, vecBins[ii].m_nEnd-1);
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            bBinIsInPMTOverlapRange          = vecBins[ii].m_nBegin >= PMT_OVERLAP_BEGIN && vecBins[ii].m_nEnd-1 <= PMT_OVERLAP_END;

            if( !bNoBinSkip && (bBinMustBeEmptyBecauseOfBitShift || bBinIsInPMTOverlapRange) )
                continue;   // skip that bin
        }
        
        printf("%10d\t%10d\t%20d\t%11.2f\n", vecBins[ii].m_nBegin, vecBins[ii].m_nEnd-1, vecBins[ii].m_nCount,
                                             (double)vecBins[ii].m_nCount/max(1.0, (double)nTotalCount)*100.0);
        if( NULL != pfOut )
            fprintf(pfOut, "%10d\t%10d\t%20d\t%11.2f\n", vecBins[ii].m_nBegin, vecBins[ii].m_nEnd-1, vecBins[ii].m_nCount,
                                                         (double)vecBins[ii].m_nCount/max(1.0, (double)nTotalCount)*100.0);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if( NULL != pfOut )
        fclose(pfOut);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    printf("%10d\t          \t%20d\t%11.2f\n", vecBins[vecBins.size()-1].m_nBegin, vecBins[vecBins.size()-1].m_nCount,
                                              (double)vecBins[vecBins.size()-1].m_nCount/max(1.0, (double)nTotalCount)*100.0);
    
    printf("-------------------------------------------------------------------\n");
    printf("Total     \t          \t%20d\n", nTotalCount);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool bOutputPixelIntensityBitTestResults(std::vector<int>& vecBitBins)
{
    printf("\n\nNumber of pixels having intensity values with a specified bit set\n");
    printf("------------------------------------\n");
    printf("        Bit\t          Num pixels\n");
    printf("------------------------------------\n");
    
    printf("No bits set\t%20d\n", vecBitBins[0]);

    for(int ii=1; ii < vecBitBins.size()-1; ii++)
    {
        printf("%11d\t%20d\n", ii, vecBitBins[ii]);
    }
    printf("     Others\t%20d\n",vecBitBins[vecBitBins.size()-1]);

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////
void vError(const int nErrorNum, const Cstring& sMessage)
{
    if( 0 != nErrorNum )
    {
        cout << sMessage << endl;
    }
    
    cout << "\ndtpixhist creates two tables representing histograms of pixel intensities.\n"
       "One table is created according to user-specified parameters: maximum and\n"
       "minimum pixel intensities, number of bins, bin width, etc. Another table\n"
       "represents a test whether a particular bit is present in any of the pixel\n"
       "values.\n"
       "\n"
       "Notes on using dtpixhist with IP images:\n\n"
       "The primary purpose of dtpixhist is to check the two read-out PMTs. PMT1\n"
       "is used to provide data for weak pixels. PMT2 is used for strong pixels.\n"
       "The circuitry that decides which PMT is used for a given pixel has hysterisis,\n"
       "so there is a band between 8191 and 16383 where the output from either PMT\n"
       "might be used AND the PMT used depends on whether it is on an up slope or down\n"
       "slope. Below 8191 only PMT1 is read and above 16383 only PMT2 is read.\n"
       "In addition, above 16363 only every 32nd value (8th value for R-AXIS II) can\n"
       "be observed, because of the 5 bit shift (3 bit shift for R-AXIS-II), used to\n"
       "store-restore the numbers.\n"
       "\n"
       "If a particular bit is not set, that suggests an ADC output line maybe bad. If\n"
       "pixels above 16383 are never observed, it suggests the HV cable for PMT2 is\n"
       "broken. If there are no values below 8192 it suggests the HV cable for PMT1\n"
       "is broken. If the histogram shows an anomalous number of pixels near 8191 or\n"
       "16383 a PMT alignment, adjusting HV1 or HV2, or replacing a PMT might be in\n"
       "order.\n"
       "\n"
       "By default, dtpixhist will not list any pixel intensity bins between 8192 and\n"
       "16383. Neither will it list bins, where no pixel values can be observed\n"
       "because of the bit shift. The user can override this default behavior and have\n"
       "all bins listed.\n";
       
    cout << "\ndtpixhist - Usage:\n"
       "dtpixhist sImageFile [sHeaderFile] -seq n1 [n2] [-min nMinIntensity]\n"
       "[-max nMaxIntensity] [-bins nBinsNumber] [-binwidth nBinWidth | pow2 | pow10]\n"
       "[-out sOutFile] [-noskipbins]\n"
       "Command line arguments:\n"
       "\n"
       "sImageFile         Image file path. If this is not specified, a header file\n"
       "                   must be specified.\n"
       "\n"
       "sHeaderFile        Header file path. If this is not specified, an image file\n"
       "                   must be specified.\n"  
       "\n"
       "n1 n2              Image sequence begin and end. If n2 is not specified,\n"
       "                   it is assumed n1=n2. There could be more than one -seq\n"
       "                   options specified.\n"
       "\n"
       "nMinIntensity      Minimum pixel intensity for the histogram.\n"
       "\n"
       "nMaxIntensity      Maximum pixel intensity for the histogram.\n"
       "\n"
       "nBinsNumber        Number of equal width bins between nMinIntensity and \n"
       "                   nMaxIntensity.\n"
       "\n"
       "nBinWidth          A number, specifying all bins' equal width.\n"
       "                   NOTE: If option -binwidth is specified, do not specify\n"
       "                   option -bins.\n"
       "\n"
       "pow2               A keyword, specifying base-2 logarithmic intensity scale\n"
       "                   NOTE: If option -binwidth is specified, do not specify\n"
       "                   option -bins.\n"
       "\n"
       "pow10              A keyword, specifying base-10 logarithmic intensity scale\n"
       "                   NOTE: If option -binwidth is specified, do not specify\n"
       "                   option -bins.\n"
       "\n"
       "sOutFile           An output file name for the user-specified pixel intensity\n"
       "                   histogram.\n"
       "\n"
       "-noskipbins        Do not skip the bins, corresponding to the PMT overlap, or\n"
       "                   where no pixel values should be observed, because of the\n"
       "                   bit shift.\n"
       "\n"
       "Examples:\n"
       "\n"
       "   dtpixhist lysodata001.osc -seq 1 -binwidth pow2\n"
       "   dtpixhist lysodata001.img -seq 1 -seq 12 14 -min 0 -max 40000 -bins 400\n"
       "   dtpixhist lysodata001.osc -seq 1 100 -min 0 -max 40000 -binwidth 400\n"
       "\n"
       << flush;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
bool bSetCustomLinearBins(int nMinPixelValue, 
                          int nMaxPixelValue, 
                          int nBinWidth,
                          std::vector<_tagPIXEL_INTENSITY_BIN>& vecBins)
{
    if( nMaxPixelValue <= nMinPixelValue )
        return false;
    if( nBinWidth < 1 )
        return false;

    int     nBinMaxCurrent = nMinPixelValue;                    // Just initializing to start the while loop
    int     nBinMinCurrent = nBinMaxCurrent;                    // Just initializing to start the while loop

    // Add one bin to count values less than min
    vecBins.push_back(_tagPIXEL_INTENSITY_BIN( -1, nBinMaxCurrent));

    while(true)
    {
        nBinMinCurrent = nBinMaxCurrent;
        nBinMaxCurrent = nBinMinCurrent + nBinWidth >= nMaxPixelValue ? nMaxPixelValue : nBinMinCurrent + nBinWidth;
    
        vecBins.push_back(_tagPIXEL_INTENSITY_BIN(nBinMinCurrent, nBinMaxCurrent));
    
        if( nBinMaxCurrent >= nMaxPixelValue )
            break;
    }
    // Add one more bin to count values equal or greater than max
    vecBins.push_back(_tagPIXEL_INTENSITY_BIN(nBinMaxCurrent, -1));

    return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
bool bSetCustomLogarithmicBins(int nMinPixelValue, 
                               int nMaxPixelValue, 
                               int nPowerBaseForBinning,
                               std::vector<_tagPIXEL_INTENSITY_BIN>& vecBins)
{
    if( nMaxPixelValue <= nMinPixelValue )
        return false;
    if( nPowerBaseForBinning < 2 )
        return false;
    
    int     nBinMaxCurrent = nMinPixelValue;                    // Just initializing to start the while loop
    int     nBinMinCurrent = nBinMaxCurrent;                    // Just initializing to start the while loop

    // Add one bin to count values less than min
    vecBins.push_back(_tagPIXEL_INTENSITY_BIN( -1, nMinPixelValue));
    
    while(true)
    {
        nBinMinCurrent = nBinMaxCurrent;
        nBinMaxCurrent = nBinMinCurrent * nPowerBaseForBinning > nMaxPixelValue ? nMaxPixelValue : nBinMinCurrent * nPowerBaseForBinning;
        
        vecBins.push_back(_tagPIXEL_INTENSITY_BIN(nBinMinCurrent, nBinMaxCurrent));
        
        if( nBinMaxCurrent >= nMaxPixelValue )
            break;
    }
    
    // Add one more bin to count values equal or more than max
    vecBins.push_back(_tagPIXEL_INTENSITY_BIN(nBinMaxCurrent, -1));

    return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// For IPs: above 16363 we need only every 32nd value (8th value for R-AXIS II), all
// other values are zero because of the 5 bit shift (3 bit shift for the
// II) used to compress the pixel values.
bool bIsBinMustBeEmptyOnIP(int nBegin, int nEnd, bool bRAXIS_II)
{
    if( nEnd <= PMT_OVERLAP_END )
        return false;  // no bit shift below PMT_OVERLAP_END

    //////////////////////////////////////////////////////////////////////////////////////
    // Otherwise see if all values in the bin cannot materialize because of the bin shift
    int     nEveryBinPresent = bRAXIS_II ? 8 : 32;

    int     nQuotient_begin = nBegin / nEveryBinPresent;
    int     nRemainder_begin = nBegin - (nQuotient_begin * nEveryBinPresent);
    
    int     nQuotient_end = nEnd / nEveryBinPresent;
    int     nRemainder_end = nEnd - (nQuotient_end * nEveryBinPresent);

    if( 0 == nRemainder_begin || 0 == nRemainder_end )
        return false;  // at least one of the bin's ends IS divisable by nEveryBinPresent
    
    if( nQuotient_begin != nQuotient_end )
        return false;  // a divisable number is somewhere between the begin and end 

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool bIsRAXIS_II_image(Cimage_header* poHeader)
{
    if( !poHeader )
        return false;

    Cstring             strDetectorPrefix(""); 
    if ( 0 != poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;

    if( strDetectorPrefix.contains("CCD") )
        return false; // Should not be a CCD
    ///////////////////////////////////////////////////////////////////////////////////

    float      afDims[2] = {0.0, 0.0}; 
    if( 0 != poHeader->nGetValue(strDetectorPrefix+D_K_DetectorDimensions, 2, afDims) )
        return false;

    float      afSizes[2] = {0.0, 0.0}; 
    if( 0 != poHeader->nGetValue(strDetectorPrefix+D_K_DetectorSize, 2, afSizes) )
        return false;

    if( afDims[0] == 0.0 || afDims[1] == 0.0 || afSizes[0] == 0.0 || afDims[1] == 0.0 )
        return false; // safety

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
    const int           c_nRAXIS_II_pixel_number = 1900;
    const float         c_fRAXIS_II_size_0       = 193.23f;
    const float         c_fRAXIS_II_size_1       = 199.5f;

    
    if( afDims[0] == afDims[1] && afDims[0] == c_nRAXIS_II_pixel_number && 
        afSizes[0] == c_fRAXIS_II_size_0 && afSizes[1] == c_fRAXIS_II_size_1 )
        return true;

    if( afDims[0] == afDims[1] && afDims[0] == c_nRAXIS_II_pixel_number/2 && 
        afSizes[0] == c_fRAXIS_II_size_0*2.0f && afSizes[1] == c_fRAXIS_II_size_1*2.0f )
        return true;

    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool bIsCCD_image(Cimage_header* poHeader)
{
    if( !poHeader )
        return false;

    Cstring             strDetectorPrefix(""); 
    if ( 0 != poHeader->nGetValue(Cstring(D_K_DetectorNames), &strDetectorPrefix) )
        return false;

    if( strDetectorPrefix.contains("CCD") )
        return true;

    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
