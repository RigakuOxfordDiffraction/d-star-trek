#ifndef DIRECTBEAMINFO_H
#define DIRECTBEAMINFO_H

//# 
//# Copyright © 2000-2001 Molecular Structure Corporation
//# Copyright © 2001-2003 Rigaku/MSC, Inc.
//#                       9009 New Trails Drive
//#                       The Woodlands, TX, USA  77381
//#  
//# The contents are unpublished proprietary source
//# code of Rigaku/MSC, Inc.
//# 
//# All rights reserved   
//# 
//# DirectBeamInfo.h      Initial author: T.L.Hendrixson            Sep 2000
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
 
/****************************************************************************
 *                              Include Files                               *
 ****************************************************************************/

#include <iostream>

#include <sstream>
#include <string>

/****************************************************************************
 *                               Definitions                                *
 ****************************************************************************/

/****************************************************************************
 *                          Structure definitions                           *
 ****************************************************************************/

/****************************************************************************
 *                               Enumerations                               *
 ****************************************************************************/

/****************************************************************************
 *                                 Typedefs                                 *
 ****************************************************************************/

/****************************************************************************
 *                                Constants                                 *
 ****************************************************************************/

/****************************************************************************
 *                            External variables                            *
 ****************************************************************************/

/****************************************************************************
 *                           Forward declarations                           *
 ****************************************************************************/

class Cimage;
class C3Ddata;
class Creflnlist;
class Cfind;

/****************************************************************************
 *                           DirectBeamInfo class                           *
 ****************************************************************************/

// Description of class goes here

class DirectBeamInfo
{

//# Enumerations, structures and other type definitions

public:

   struct PeakInfo {

      PeakInfo() : m_IsUsed(false),
                   m_PixelIntensity(0),
                   m_MostIntensePixelX(-1),
                   m_MostIntensePixelY(-1),
                   m_MostIntensePixelI(0),
                   m_NetIntensity(0),
                   m_X(-1.0),
                   m_Y(-1.0),
                   m_Intensity(0.0),
                   m_SigmaI(0.0),
                   m_FWHM_X(0.0),
                   m_FWHM_Y(0.0),
                   m_SigmaWidthX(0.0),
                   m_SigmaWidthY(0.0),
                   m_Rank(0.0),
                   m_RankI(0.0),
                   m_RankR(0.0),
                   m_nBoxSize(0)
      {}

      PeakInfo &operator=(const PeakInfo &rhs)
      {
         if(this == &rhs)
            return *this;

         m_IsUsed = rhs.m_IsUsed;
         m_PixelIntensity = rhs.m_PixelIntensity;
         m_MostIntensePixelX = rhs.m_MostIntensePixelX;
         m_MostIntensePixelY = rhs.m_MostIntensePixelY;
         m_MostIntensePixelI = rhs.m_MostIntensePixelI;
         m_NetIntensity = rhs.m_NetIntensity;

         m_X = rhs.m_X;
         m_Y = rhs.m_Y;
         m_Intensity = rhs.m_Intensity;
         m_SigmaI = rhs.m_SigmaI;
         m_FWHM_X = rhs.m_FWHM_X;
         m_FWHM_Y = rhs.m_FWHM_Y;
         m_SigmaWidthX = rhs.m_SigmaWidthX;
         m_SigmaWidthY = rhs.m_SigmaWidthY;

         m_Rank = rhs.m_Rank;
         m_RankI = rhs.m_RankI;
         m_RankR = rhs.m_RankR;
         m_nBoxSize = rhs.m_nBoxSize;

         return *this;
      }

      bool m_IsUsed;

      long m_PixelIntensity;
      long m_MostIntensePixelX;
      long m_MostIntensePixelY;
      long m_MostIntensePixelI;
      long m_NetIntensity;

      double m_X;
      double m_Y;
      double m_Intensity;
      double m_SigmaI;
      double m_FWHM_X;
      double m_FWHM_Y;
      double m_SigmaWidthX;
      double m_SigmaWidthY;
      double m_Rank;
      double m_RankI;
      double m_RankR;
      int   m_nBoxSize;
   };

protected:

private:

//# Routines

public:

//# Constructor

      // Create an instance of the class.
   DirectBeamInfo(void);  

//# Destructor  

      // The destructor for this class.  This will free any memory that  
      // was allocated by the class. 
	   // <hr>
   ~DirectBeamInfo();  

//# Operators

   DirectBeamInfo &operator<<(std::ostream &out);

//# Methods

      // <P>
      //    This routine evaluates the direct beam image <VAR>pImage</VAR>.
      //    The parameters <VAR>StartLine</VAR> and <VAR>NumberOfLines</VAR>
      //    can be used to restrict the analysis of the direct beam to a
      //    smaller portion of the image.  If these parameters are not
      //    specified, then the entire image will be examined.
      // </P>
      // <P>
      //    During the evaluation, information about the direct beam, such
      //    as position, size and intensity, is calculated and an output
      //    string containing an ASCII representation of the direct beam and
      //    the information determined is generated (but not written to any
      //    output device).
      // </P>
   void Evaluate(Cimage *pImage, const int StartLine = 0,
                 const int NumberOfLines = 0);

      // This routine places the integrated intensity, I, and the error
      // in the intensity, &sigma;(I), into the parameters <VAR>Intensity</VAR>
      // and <VAR>SigmaI</VAR>, respectively.
   void GetIntensity(double &Intensity, double &SigmaI);

      // This routine places the X (fastest index) and Y (slowest index)
      // positions of the most intense pixel in the direct beam in the 
      // parameters <VAR>XPosition</VAR> and <VAR>YPosition</VAR>,
      // respectively.  The values of the positions are in pixels.
   void GetMostIntensePixelPosition(long &XPosition, long &YPosition);

      // This routine returns the intensity of the most intense pixel in
      // the direct beam.
   long GetMostIntensePixelIntensity(void);

      // This routine returns the net intensity of the direct beam.  The
      // net intensity is calulated by doing a box sum over the direct
      // beam and subtracting off the total background.  The background 
      // for the box is calculated by using the outer ring of pixels in
      // the box.
   long GetNetIntensity(void);

      // This routine returns the output string.
   std::string GetOutputString(void);

      // This routine returns the intensity of the pixel closest to the
      // position of the direct beam.
   long GetPixelIntensity(void);

      // This routine returns a PeakInfo structure containing all the peak       
      // information.  Note that most all of this information can be 
      // obtained by using other <CODE>GetXXX()</CODE> routines.
   PeakInfo GetPeakInfo(void);

      // This routine places the X (fastest index) and Y (slowest index)
      // positions of the direct beam in the parameters <VAR>XPosition</VAR>
      // and <VAR>YPosition</VAR>, respectively.  The values of the
      // positions are in pixels.
   void GetPosition(double &XPosition, double &YPosition);

      // This routine returns a string containing a summary of the ranking
      // results for the candidate beam positions.
   std::string GetRankingSummary(void);

      // This routine places the width of the direct beam into the 
      // parameters <VAR>XWidth</VAR> and <VAR>YWidth</VAR> for the X
      // (fastest index) and Y (slowest index), respectively.
      // The width is calculated using the full width at half maximum pixel
      // value for the peak.
   void GetWidth_FWHM(double &XWidth, double &YWidth);

      // This routine places the width of the direct beam into the 
      // parameters <VAR>XWidth</VAR> and <VAR>YWidth</VAR> for the X
      // (fastest index) and Y (slowest index), respectively.
      // The width is calculated using the point at which the intensity
      // in the peak is 3&sigma; above the background.
   void GetWidth_3Sigma(double &XWidth, double &YWidth);

      // This routine sets the identifier used in the output string (in
      // the first line).
   void SetIdentifier(const std::string &Identifier);

      // This routine sets the filename of the mask that will be used when
      // determining the direct beam position.  <VAR>Filename</VAR> is the
      // filename of the image file that will be used as the mask.  If
      // <VAR>Filename</VAR> is an empty string, then a mask will not be
      // used in the determination.
   void SetMaskFilename(const std::string &Filename);

      // This routine sets the weights for the ranking scheme used when
      // determining the direct beam position and returns a boolean
      // indicating if the weights could be set or not.  <VAR>Intensity</VAR> 
      // is the weight for the intensity term and <VAR>Resolution</VAR> is
      // the weight for the resolution term.  If the value returned is
      // <I>true</I>, then the weights could be set.  If the value returned
      // is <I>false</I>, then the weights could not be set.
   bool SetRankingWeights(const double &Intensity, const double &Resolution);

protected:

//# These variables are defined here because the derived classes need to  
//# be able to access them. 

private:

//# There is no copy constructor  
      // A non-existant copy constructor.  This is defined as a private
      // member function, but not implemented, so that user code cannot 
      // contain the following:
      // <PRE><CODE>
      // DirectBeamInfo instance1;
      // DirectBeamInfo instance2(instance1);
      // </CODE></PRE>
      // in an attempt to create a copy of an instance of this class.
   DirectBeamInfo(const DirectBeamInfo &Copy);  

//# There is no assignment operator  
      // A non-existant assignment operator.  This is defined as a private
      // member function, but not implemented, so that user code cannot 
      // contain the following:
      // <PRE><CODE>
      // DirectBeamInfo instance1,instance2;
      // instance2 = instance1;
      // </CODE></PRE>
      // in an attempt to create a copy of an instance of this class.
   DirectBeamInfo &operator=(const DirectBeamInfo &AssignEqual);  


      // This routine analyzes the image to locate the position of the
      // direct beam and to generate information about the shape and
      // intensity of the direct beam.  The parameters <VAR>StartLine</VAR>
      // and <VAR>NumberOfLines</VAR> can be used to restrict the
      // analysis of the direct beam to a smaller portion of the image.
   void Analyze(const int StartLine, const int NumberOfLines);

      // This routine calculates the size of the direct beam in the 
      // X and Y positions, using a 3sigma criterion and full width at
      // half maximum determination.
   void CalculateDirectBeamSize(void);

      // This routine calculates the full width at half-maximum (FWHM)
      // for the direct beam in the X and Y positions.
   void CalculateFWHM(C3Ddata &BeamData, const float BkgdAvg);

      // This routine calculates the net intensity of the direct beam
      // using a box sum algorithm.
   void CalculateNetIntensity(C3Ddata &BeamData);

      // This routine calculates the size for the box that will be used
      // in the profile output.  On return, <VAR>Xmin</VAR> and 
      // <VAR>Ymin</VAR> will contain the minimum positions for the profile 
      // box, in X and Y respectively, and <VAR>Xsize</VAR> and 
      // <VAR>Ysize</VAR> will contain the size of the profile box, in X 
      // and Y respectively..
   void CalculateProfileBox(int &Xmin, int &Ymin, int &Xsize, int &Ysize);

      // This routine finds, using <VAR>Find</VAR>, the most intense spot 
      // in the reflexion list, <VAR>pRefnlist</VAR>, and stores the spot 
      // information as the current direct beam information.
   bool FindMostIntenseSpot(Creflnlist &Reflnlist, Cfind &Find);

      // This routine returns a peak containing bogus information.  The
      // peak will be roughly centered on the image.
   PeakInfo GenerateBogusPeak(void);

      // This routine generates the output string containing the direct
      // beam profile and information about the direct beam.  Note that
      // if the image containing the direct beam is smaller than 7 pixels
      // high, the profile will probably look a little strange.
   void GenerateOutputString(void);

      // This routine returns the detector prefix contained in the header
      // of the image being analyzed.
   std::string GetDetectorPrefix(void);

      // This routine returns the base intensity ranking for the peak 
      // described by <VAR>Peak</VAR>, based on the maximum intensity
      // <VAR>MaxIntensity</VAR>.
   double GetIntensityRank(const PeakInfo &Peak, const double MaxIntensity);

      // This routine determines the centroids of all the peaks found in
      // the reflexion list <VAR>Reflnlist</VAR>, using the peak finding
      // object <VAR>Find</VAR>.  Once the centroids have been determined,
      // the peak information for each peak in the reflexion list is placed
      // in an element of <VAR>pPeaks</VAR>.  It is up to the calling routine
      // to ensure that <VAR>pPeaks</VAR> contains enough storage space to
      // hold information on all of the peaks.  Upon return, the largest peak
      // intensity is placed in <VAR>MaxIntensity</VAR>.
   void GetPeakCentroids(Creflnlist &Reflnlist, Cfind &Find, PeakInfo *pPeak, 
                         double &MaxIntensity);

      // This routine returns the base resolution (typically in terms of 
      // distance from the center of the image) ranking for the peak
      // described by <VAR>Peak</VAR>. 
   double GetResolutionRank(const PeakInfo &Peak);

      // This routine locates the position of the direct beam in the image, 
      // recording its position, intensity and sigma(I).  If the image
      // being evaluated does not exist, then the values returned will be
      // pseudo-randomly generated fake values.  If no candidate for the 
      // direct beam can be found, then the routine will return 
      // <VAR>false</VAR>, otherwise it will return <VAR>true</VAR>.
   bool LocateDirectBeam(void);

      // This routine returns a pair containing the lines that make up the
      // edges of the image.  Lines that only contain 0 are assumed to be 
      // padding and not part of the actual image.  The first element of the 
      // pair returned is the top line in the actual image, and the second
      // element of the pair returned is the bottom line in the actual image.
   std::pair<int,int> LocateImageEdges(Cimage &Image);

      // This routine pads a portion of the edges of the image, if they 
      // contain 0 values, with an average background value.  This helps
      // locate the direct beam when it is near the edge of the image.
      //
      // Note that this routine is defined, but is not used since I 
      // haven't decided if this class should be modifying the image
      // since it doesn't create its own copy (which could be very
      // expensive in memory terms).
   void PadImageEdges(void);

      // This routine determines the position of the half-height of one
      // side of a peak profile of length <VAR>Length</VAR>.  It is
      // assumed that the maximum of the peak is the first array element
      // of <VAR>pProfile</VAR>.  The routine will return the postion of
      // the half height relative to the first array element, and will
      // interpolate the position if required.
   double PositionOfHalfMaximum(const double *pProfile, const int Length);

      // This routine returns a string containing a text representation of
      // the peak profile for the line <VAR>Line</VAR> of the profile.
      // <VAR>Data</VAR> is the data object containing the peak to be
      // used.  <VAR>Xsize</VAR>,<VAR>Xcenter</VAR> and <VAR>Ycenter</VAR>
      // are the size of the profile in the X direction, and the center of
      // the peak in the X and Y directions, respectively.  
      
   std::string ProfileString(C3Ddata &Data, const int Line, const int Xsize,
                             const int Xorigin, const int Yorigin);

      // This routine returns a string containing the ranking summary of
      // the peaks pointed to by the elements of <VAR>pPeaks</VAR>.  At
      // most, information on <VAR>N</VAR> peaks will be contained in the
      // summary, with the summary stopping when a peak with a rank of
      // 0.0 is encountered.
   std::string RankingSummary(PeakInfo **pPeaks, const int N);

      // This routine sorts the pointers contains in <VAR>pPeaks</VAR> 
      // based on the intensity of the peak objects that they point to.
      // The sort is in descending order, and only the inclusive range
      // [<VAR>LowerBound</VAR>,<VAR>UpperBound</VAR>] are sorted.
   void SortByRank(DirectBeamInfo::PeakInfo **pPeaks, const int LowerBound, 
                   const int UpperBound);

      // This routine contains a string containing text information used
      // in the peak profile output.  <VAR>Line</VAR> indicates the current
      // line in the profile output, relative to the line in the profile 
      // output that contains the center of the peak.
   std::string TextInfo(const int Line);


//# Member variables

public:

protected:

private:

   static int m_BoxSize;
   static int m_BoxWidth;

   double m_WtI;
   double m_WtR;

   PeakInfo m_PeakInfo;

   Cimage *m_pImage;

   std::string m_Identifier;
   std::string m_MaskFilename;
   std::string m_OutputString;
   std::string m_RankingSummary;
};

inline std::ostream &operator<<(std::ostream &o, DirectBeamInfo &a)
{ a.operator<<(o); return o; }

#endif /* DIRECTBEAMINFO_H */

