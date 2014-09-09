//
// Copyright (c) 1997 Molecular Structure Corporation
//
// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No.
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// dtrekdefs.h     Initial author: N. Chen           Summer 1997
//    This header file defines most strings for d*TREK.
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

#ifndef DT_DTREKDEFS_H
#define DT_DTREKDEFS_H

#ifdef OSF1
#define DTREK_PLATFORM "O"
#elif defined(BNL_LICENSE)
#define DTREK_PLATFORM "BNL"
#elif defined(LINUX)
#define DTREK_PLATFORM "L"
#elif defined(IRIX)
#define DTREK_PLATFORM "I"
#elif defined (WIN32)
    #ifdef _DEBUG
        #ifdef VC9
            #define DTREK_PLATFORM "W9D"
        #else
            #define DTREK_PLATFORM "WD"
        #endif
    #else
        #ifdef VC9
            #define DTREK_PLATFORM "W9R"
        #else
            #define DTREK_PLATFORM "WR"
        #endif
    #endif 
#elif defined(OSX)
#define DTREK_PLATFORM "D"
#elif defined(SUNOS)
#define DTREK_PLATFORM "S"
#else
#define DTREK_PLATFORM "U"
#endif

#ifdef SSI_PC
    #ifdef VC9
        #ifdef _DEBUG
            #define D_K_DTREKVersion     ("d*TREK version 9.9.9.10 W9DSSI -- " __DATE__)
        #else
            #define D_K_DTREKVersion     ("d*TREK version 9.9.9.10 W9RSSI -- " __DATE__)
        #endif
    #else
        #ifdef _DEBUG
            #define D_K_DTREKVersion     ("d*TREK version 9.9.9.10 WDSSI -- " __DATE__)
        #else
            #define D_K_DTREKVersion     ("d*TREK version 9.9.9.10 WRSSI -- " __DATE__)
        #endif
    #endif
#else
    #ifdef DEMO_LICENSE
        #define D_K_DTREKVersion     ("d*TREK version 9.9.9.10" DTREK_PLATFORM "Dz" " -- " __DATE__)
    #else
        #define D_K_DTREKVersion     ("d*TREK version 9.9.9.10" DTREK_PLATFORM " -- " __DATE__)
    #endif
#endif


// Installation
#define D_K_DTREK_ROOT                      "DTREK_ROOT"

// Cimage_header Strings
#define D_K_Comment                        "COMMENT"
#define D_K_Comment2                       "COMMENT2"
#define D_K_DtdisplayOrientation           "DTDISPLAY_ORIENTATION"
#define D_K_DtdisplayTile                  "DTDISPLAY_TILE"
#define D_K_HeaderBytes                    "HEADER_BYTES"
#define D_K_ByteOrder                      "BYTE_ORDER"
#define D_K_ByteOrderOriginal              "BYTE_ORDER_ORIGINAL"
#define D_K_Compression                    "COMPRESSION"
#define D_K_Dim                            "DIM"
#define D_K_DataType                       "Data_type"
#define D_K_Filename                       "FILENAME"
#define D_K_NumReflns                      "NUM_REFLNS"
#define D_K_Raxis                          "RAXIS"
#define D_K_RxPrefix                       "RX_"
#define D_K_RaxisCompressionRatio          "RAXIS_COMPRESSION_RATIO"
#define D_K_RaxisReadLines                 "RAXIS_READ_LINES"
#define D_K_RaxisReadStart                 "RAXIS_READ_START"
#define D_K_RaxisRecordSize                "RAXIS_RECORD_SIZE"
#define D_K_SaturatedValue                 "SATURATED_VALUE"
#define D_K_MinRawPixOKValue               "DTREK_NONUNF_OKVALUE"
#define D_K_Size1                          "SIZE1"
#define D_K_Size2                          "SIZE2"
#define D_K_Type                           "TYPE"
#define D_K_LittleEndian                   "little_endian"
#define D_K_BigEndian                      "big_endian"
#define D_K_OtherEndian                    "Other_endian"
#define D_K_DtrekCBFPixDataType            "DTREK_CBF_PIXDATATYPE"

#define D_K_CoordNames                     "COORD_NAMES"
#define D_K_CoordUnits                     "COORD_UNITS"

#define D_K_OriginalImageFormat            "ORIGINAL_IMAGE_FORMAT"

// Cimage Strings
#define D_K_SignedChar                     "signed char"
#define D_K_UnsignedChar                   "unsigned char"
#define D_K_ShortInt                       "short int"
#define D_K_LongInt                        "long int"
#define D_K_LongIntMarty                   "long_integer"
#define D_K_SignedLongMarty                "signed_long"
#define D_K_MadMarty                       "mad"
#define D_K_UnsignedShortInt               "unsigned short int"
#define D_K_UnsignedShortIntMarty          "unsigned_short"
#define D_K_UnsignedLongInt                "unsigned long int"
#define D_K_FloatIEEE                      "float IEEE"
#define D_K_Compressed                     "Compressed"
#define D_K_OtherType                      "Other_type"
#define D_K_WinBMP                         "WinBMP"
#define D_K_BitmapSize                     "BitmapSize"
#define D_K_BitmapType                     "BitmapType"
#define D_K_BitmapRLE                      "BitmapRLE"


// ?
#define D_K_Cog1                           "COG_1"
#define D_K_Cog2                           "COG_2"
#define D_K_Fint                           "FINT"

// C3Ddata Strings
#define D_K_ShoeboxSize                    "SHOEBOX_SIZE"
#define D_K_ShoeboxOrig                    "SHOEBOX_ORIG"
#define D_K_ShoeboxDir                     "SHOEBOX_DIR"

// Ccrystal Strings
#define D_K_CrystalDescription             "CRYSTAL_DESCRIPTION"
#define D_K_CrystalResidInfo               "CRYSTAL_RESID_INFO"
#define D_K_CrystalKey                     "CRYSTAL_KEY"

#define D_K_CrystalPrefix                  "CRYSTAL_"

#define D_K_CrystalMolecularFormula        "CRYSTAL_MOLECULAR_FORMULA"
#define D_K_CrystalMorphologyColor         "CRYSTAL_MORPHOLOGY_COLOR"
#define D_K_CrystalMorphologyMount         "CRYSTAL_MORPHOLOGY_MOUNT"
#define D_K_CrystalMorphologyShape         "CRYSTAL_MORPHOLOGY_SHAPE"
#define D_K_CrystalMorphologySize          "CRYSTAL_MORPHOLOGY_SIZE"
#define D_K_CrystalOtherInfo               "CRYSTAL_OTHER_INFO"
#define D_K_CrystalOffsetSubG              "CRYSTAL_OFFSETS"

// Options added to .REF header info
#define D_K_HeaderCrystalInfo              D_K_CrystalMask
#define D_K_HeaderDtfindInfo               "DTFIND*"
#define D_K_HeaderDontAdd                  "SIZE"

#define D_K_CrystalMask                    "CRYSTAL*MOSAICITY;CRYSTAL*ORIENT*;CRYSTAL*SPACEGROUP;CRYSTAL*UNIT_CELL;CRYSTAL*NAME;CRYSTAL*SIGMA;CRYSTAL*TWINLAWS;CRYSTAL*TWINFRACTIONS;CRYSTAL*RECIPSHIFT"

#define D_K_CrystalXMosaicity              "MOSAICITY"
#define D_K_CrystalXOrientAngles           "ORIENT_ANGLES"
#define D_K_CrystalXOrientVectors          "ORIENT_VECTORS"
#define D_K_CrystalXSpacegroup             "SPACEGROUP"
#define D_K_CrystalXUnitCell               "UNIT_CELL"
#define D_K_CrystalXName                   "NAME"
#define D_K_CrystalXTwinLaws               "TWINLAWS"
#define D_K_CrystalXTwinFraction           "TWINFRACTION"
#define D_K_CrystalXRecipShift             "RECIP_SHIFT"
#define D_K_CrystalXNumTwinLaws            "NUM_TWIN_LAWS"

// Cdetector Strings
#define D_K_DetectorPrefix                 "DETECTOR_"
#define D_K_DetectorADCOffset              "DETECTOR_ADC_OFFSET"
#define D_K_DetectorDescription            "DETECTOR_DESCRIPTION"
#define D_K_DetectorDimensions             "DETECTOR_DIMENSIONS"
#define D_K_DetectorKey                    "DETECTOR_KEY"
#define D_K_DetectorIdentification         "DETECTOR_IDENTIFICATION"
#define D_K_DetectorNames                  "DETECTOR_NAMES"
#define D_K_DetectorNumber                 "DETECTOR_NUMBER"
#define D_K_DetectorRefineFlags            "DETECTOR_REFINE_FLAGS"
#define D_K_DetectorOptions                "DETECTOR_OPTIONS"
#define D_K_DetectorSize                   "DETECTOR_SIZE"
#define D_K_DetectorTemperature            "DETECTOR_TEMPERATURE"
#define D_K_DetectorTranslation            "DETECTOR_TRANSLATION"
#define D_K_DetectorType                   "DETECTOR_TYPE"
#define D_K_DetectorVectors                "DETECTOR_VECTORS"

#define D_K_DetectorBiasVoltage            "DETECTOR_BIAS_VOLTAGE"

#define D_K_CylindricalDetectorRadius      "DETECTOR_CYLINDER_RADIUS"

#define D_K_UnbinnedDimensions             "UNBINNED_DIMENSIONS"
#define D_K_UnbinnedBeamPosition           "UNBINNED_BEAM_POSITION"

// Crotation Strings
#define D_K_Rotation                       "ROTATION"
#define D_K_RotAxisName                    "ROTATION_AXIS_NAME"
#define D_K_RotVector                      "ROTATION_VECTOR"
#define D_K_RotLimits                      "ROTATION_LIMITS"

// Cgoniometer Strings
#define D_K_GonioDescription               "GONIO_DESCRIPTION"
#define D_K_GonioLabCoordSystem            "GONIO_LABORATORY_COORDINATE_SYSTEM"
#define D_K_GonioKey                       "GONIO_KEY"
#define D_K_GonioNames                     "GONIO_NAMES"
#define D_K_GonioNumValues                 "GONIO_NUM_VALUES"
#define D_K_GonioUnits                     "GONIO_UNITS"
#define D_K_GonioValues                    "GONIO_VALUES"
#ifdef SSI_PC
#define D_K_GonioValuesHardware            "GONIO_VALUES_HARDWARE"
#endif
#define D_K_GonioValuesMin                 "GONIO_VALUES_MIN"
#define D_K_GonioValuesMax                 "GONIO_VALUES_MAX"
#define D_K_GonioValuesSigma               "GONIO_VALUES_SIGMA"
#define D_K_GonioVectors                   "GONIO_VECTORS"
#define D_K_GonioOffsetTableGonio          "GONIO_OFFSET_TABLE_GONIO"
#define D_K_GonioOffsetTableOffset         "GONIO_OFFSET_TABLE_OFFSET"
#define D_K_GonioOffsetTableEntries        "GONIO_OFFSET_TABLE_ENTRIES"
#define D_K_GonioOffsetTableRotAxis        "GONIO_OFFSET_TABLE_ROT_AXIS"
#define D_K_GonioCollisionOffsets          "GONIO_COLLISION_OFFSET"
#define D_K_GonioScanAxisNames             "GONIO_SCAN_AXES"

// Cnonunf Strings
#define D_K_NonunfDenominator              "NONUNF_DENOMINATOR"
#define D_K_NonunfOffset                   "NONUNF_OFFSET"
#define D_K_NonunfFlag1                    "NONUNF_FLAG1"
#define D_K_NonunfFlag2                    "NONUNF_FLAG2"
#define D_K_NonunfFlag3                    "NONUNF_FLAG3"
#define D_K_NonunfFlag4                    "NONUNF_FLAG4"
#define D_K_NonunfInfo                     "NONUNF_INFO"
#define D_K_NonunfKey                      "NONUNF_KEY"
#define D_K_NonunfNumerator                "NONUNF_NUMERATOR"
#define D_K_NonunfType                     "NONUNF_TYPE"

#define D_K_NonunfStateDarkDCoffsetNonunf  "Dark_dcoffset_nonunf"
#define D_K_NonunfStateDarkNonunf          "Dark_nonunf"
#define D_K_NonunfStateNone                "None"
#define D_K_NonunfStateSimpleScale         "Simple_scale"
#define D_K_NonunfStateSimpleMask          "Simple_mask"
#define D_K_NonunfStateUnknown             "Unknown"

// Cspatial Strings
#define D_K_SpatialBeamPosn                "SPATIAL_BEAM_POSITION"
#define D_K_SpatialBeamPosition            "SPATIAL_BEAM_POSITION"
#define D_K_SpatialDistortionInfo          "SPATIAL_DISTORTION_INFO"
#define D_K_SpatialDistortionType          "SPATIAL_DISTORTION_TYPE"
#define D_K_SpatialDistortionVectors       "SPATIAL_DISTORTION_VECTORS"
#define D_K_SpatialKey                     "SPATIAL_DISTORTION_KEY"
#define D_K_SpatialXYSigma                 "SPATIAL_XY_SIGMA"

#define D_K_SpatialTypeComplex             "Complex_spatial"
#define D_K_SpatialTypeInterp              "Interp_spatial"
#define D_K_SpatialTypeNone                "None"
#define D_K_SpatialTypeSimple              "Simple_spatial"
#define D_K_SpatialTypeUnknown             "Unknown"

// Calibrate Strings
#define D_K_CalibPixelSize                 "PIXEL_SIZE"
#define D_K_CalibPscale                    "PSCALE"
#define D_K_CalibXintStart                 "XINT_START"
#define D_K_CalibYintStart                 "YINT_START"
#define D_K_CalibXinvStart                 "XINV_START"
#define D_K_CalibYinvStart                 "YINV_START"
#define D_K_CalibXintStep                  "XINT_STEP"
#define D_K_CalibYintStep                  "YINT_STEP"
#define D_K_CalibXinvStep                  "XINV_STEP"
#define D_K_CalibYinvStep                  "YINV_STEP"
#define D_K_CalibXSize                     "X_SIZE"
#define D_K_CalibYSize                     "Y_SIZE"
#define D_K_CalibXBeam                     "X_BEAM"
#define D_K_CalibYBeam                     "Y_BEAM"
#define D_K_CalibBadFlag                   "BAD_FLAG"
#define D_K_CalibXCenter                   "X_CENTER"
#define D_K_CalibYCenter                   "Y_CENTER"
#define D_K_CalibXPtCenter                 "X_PT_CENTER"
#define D_K_CalibYPtCenter                 "Y_PT_CENTER"
#define D_K_CalibXScale                    "X_SCALE"
#define D_K_CalibYScale                    "Y_SCALE"
#define D_K_CalibRatio                     "RATIO"
#define D_K_CalibVertSlope                 "VER_SLOPE"
#define D_K_CalibHorzSlope                 "HORZ_SLOPE"
#define D_K_CalibRadialA1                  "RADIAL_A1"
#define D_K_CalibRadialA                   "RADIAL_A"
#define D_K_CalibRadialB                   "RADIAL_B"
#define D_K_CalibRadialC                   "RADIAL_C"
#define D_K_CalibSpacing                   "SPACING"
#define D_K_CalibXMaskPoints               "X_MASK_POINTS"
#define D_K_CalibYMaskPoints               "Y_MASK_POINTS"
#define D_K_CalibLeftMaskPoint             "LEFT___MASK_POINT"
#define D_K_CalibCenterMaskPoint           "CENTER_MASK_POINT"
#define D_K_CalibRightMaskPoint            "RIGHT__MASK_POINT"
#define D_K_CalibBottomMaskPoint           "BOTTOM_MASK_POINT"
#define D_K_CalibTopMaskPoint              "TOP____MASK_POINT"
#define D_K_CalibNonunfFile                "NONUNF"
#define D_K_CalibDistorBasename            "DISTOR"
#define D_K_CalibDarkFile                  "DARK"
#define D_K_CalibMaskFile                  "MASK"
#define D_K_CalibFloodFile                 "FLOOD"
#define D_K_CalibReferenceFile             "REFER"
#define D_K_CalibBadpixelFile              "BADPIX"
#define D_K_CalibDumpFile                  "DUMP"
#define D_K_CalibPrefix                    "CALIB_"
#define D_K_Calib_number_modules           "CALIB_number_modules"
#define D_K_Calib_flood_radial             "CALIB_flood_radial"
#define D_K_Calib_flood_interpolate        "CALIB_flood_interpolate"
#define D_K_Calib_flood_geometry           "CALIB_flood_geometry"
#define D_K_Calib_flood_film               "CALIB_flood_film"
#define D_K_Calib_darksub_const            "CALIB_darksub_const"
#define D_K_Calib_darksub_scale            "CALIB_darksub_scale"
#define D_K_Calib_pixel_sd                 "CALIB_pixel_sd"
#define D_K_Calib_max_pixel                "CALIB_max_pixel"
#define D_K_Calib_min_pixel                "CALIB_min_pixel"
#define D_K_Calib_radial_distortion        "CALIB_radial_distortion"
#define D_K_Calib_bad_peaks                "CALIB_bad_peaks"
#define D_K_Calib_mask_angle               "CALIB_mask_angle"
#define D_K_Calib_min_peak                 "CALIB_min_peak"
#define D_K_Calib_cent_to_maxval           "CALIB_cent_to_maxval"
#define D_K_Calib_cent_to_pred             "CALIB_cent_to_pred"
#define D_K_Calib_background_pixels        "CALIB_background_pixels"
#define D_K_Calib_peak_size                "CALIB_peak_size"
#define D_K_Calib_peak_distance            "CALIB_peak_distance"
#define D_K_Calib_center_peak              "CALIB_center_peak"
#define D_K_Calib_search_radius            "CALIB_search_radius"
#define D_K_Calib_search_center            "CALIB_search_center"
#define D_K_Calib_vertical_limits          "CALIB_vertical_limits"
#define D_K_Calib_horizontal_limits        "CALIB_horizontal_limits"
#define D_K_Calib_xtod_distance            "CALIB_xtod_distance"
#define D_K_Calib_masktod_distance         "CALIB_masktod_distance"
#define D_K_Calib_beam_position            "CALIB_beam_position"
#define D_K_Calib_mask_spacing             "CALIB_mask_spacing"
#define D_K_Calib_pixel_size               "CALIB_pixel_size"
#define D_K_Calib_edges                    "CALIB_edges"
#define D_K_Calib_dump_mode                "CALIB_dump_mode"
#define D_K_DtcalibFiles                   "DTCALIB_FILES"
#define D_K_DtcalibOptions                 "DTCALIB_OPTIONS"

// Cscan Strings
#define D_K_ScanCrysDatum                  "SCAN_CRYS_RELZERO"
#define D_K_ScanDetDatum                   "SCAN_DET_RELZERO"
#define D_K_ScanDezinger                   "SCAN_DEZINGER_IMAGE"
#define D_K_ScanShutterMode                "SCAN_SHUTTER_MODE"
#define D_K_ScanDetBinMode                 "SCAN_DET_BINMODE"
#define D_K_ScanDetectorOptions            "SCAN_DETECTOR_OPTIONS"
#define D_K_ScanKey                        "SCAN_KEY"
#define D_K_ScanMode                       "SCAN_MODE"

#define D_K_ScanPrefix                     "SCAN_"

#define D_K_ScanSeqInfo                    "SCAN_SEQ_INFO"
#define D_K_ScanTemplate                   "SCAN_TEMPLATE"
#define D_K_ScanWavelength                 "SCAN_WAVELENGTH"
#define D_K_ScanWavelengthOpts             "SCAN_WAVELENGTH_OPTIONS"
#define D_K_ScanSeqSelected                "SCAN_SEQ_SELECTED"
#define D_K_ScanComment                    "Scan_Comment"

#define D_K_ScanModeStillC                 "Still_Closed"
#define D_K_ScanModeStillO                 "Still_Open"
#define D_K_ScanModeScanC                  "Scan_Closed"
#define D_K_ScanModeScanO                  "Scan_Open"
#define D_K_ScanModeUnknown                "Scan_Unknown"

// Csource Strings
#define D_K_SourceCrossfire                "SOURCE_CROSSFIRE"
#define D_K_SourceFocus                    "SOURCE_FOCUS"
#define D_K_SourceIntensity                "SOURCE_INTENSITY"
#define D_K_SourcePolarz                   "SOURCE_POLARZ"
#define D_K_SourcePrefix                   "SOURCE_"
#define D_K_SourceRefineFlags              "SOURCE_REFINE_FLAGS"
#define D_K_SourceSize                     "SOURCE_SIZE"
#define D_K_SourceKey                      "SOURCE_KEY"
#define D_K_SourceSpectralDispersion       "SOURCE_SPECTRAL_DISPERSION"
#define D_K_SourceValues                   "SOURCE_VALUES"
#define D_K_SourceValuesSigma              "SOURCE_VALUES_SIGMA"
#define D_K_SourceVectors                  "SOURCE_VECTORS"
#define D_K_SourceWavelength               "SOURCE_WAVELENGTH"
#define D_K_SourceVoltage                  "SOURCE_VOLTAGE"
#define D_K_SourceAmperage                 "SOURCE_AMPERAGE"

// Cwavelength Strings
#define D_K_Wavelength                     "WAVELENGTH"
#define D_K_WavelengthKey                  "WAVELENGTH_KEY"

// Strategy Strings
#define D_K_StrategyTestPrefix             "STRATEGY_INPUT_"
#define D_K_StrategyInputPrefix            "STRATEGY_INPUT_"
#define D_K_StrategyResultPrefix           "STRATEGY_OUTPUT_"
#define D_K_DtstrategyPercentComplete      "DTSTRATEGY_PERCENT_COMPLETE"
#define D_K_DtstrategyRedundancy           "DTSTRATEGY_REDUNDANCY"
#define D_K_DtstrategyRedundancyDev        "DTSTRATEGY_REDUNDANCY_DEVIATION"

// Crefine Strings
#define D_K_DtrefineFiles                  "DTREFINE_FILES"
#define D_K_DtrefineOptions                "DTREFINE_OPTIONS"
#define D_K_DtrefineRmsMm                  "DTREFINE_RMS_MM"
#define D_K_DtrefineRmsDeg                 "DTREFINE_RMS_DEG"
#define D_K_DtrefineLambda                 "DTREFINE_LAMBDA"
#define D_K_DtrefineLambdaProfile          "DTREFINE_LAMBDA_PROFILE"
#define D_K_DtrefineResidProfile           "DTREFINE_RESID_PROFILE"
#define D_K_DtrefineProfileCount           "DTREFINE_PROFILE_COUNT"
#define D_K_DtrefineResoUsed               "DTREFINE_RESO_USED"
#define D_K_DtrefineReflectionNumbers      "DTREFINE_REFLN_NUMBERS"
#define D_K_DtrefineTwinsRefined           "DTREFINE_TWINS_REFINED"
#define D_K_DtrefineTwinRejected             "DTREFINE_TWIN_REJECTED"
#define D_K_DtrefineReflectionNumbersPhotons "DTREFINE_REFLN_NUMBERS_PHOTONS"

// Cfind Strings
#define D_K_DtfindOptions                  "DTFIND_OPTIONS"
#define D_K_DtfindSpotSize                 "DTFIND_SPOTSIZE"
#define D_K_DtfindSeqOptions               "DTFIND_SEQ_OPTIONS"
#define D_K_DtfindAvgSpotsPerImage         "DTFIND_AVG_SPOTS_PER_IMAGE"
#define D_K_DtfindTotalPhotons             "DTFIND_TOTAL_PHOTONS"
#define D_K_DtfindPhotonValues             "DTFIND_PHOTON_VALUES"
#define D_K_DtfindPhotonIoverSig           "DTFIND_PHOTON_IOVERSIG"
#define D_K_DtfindFractionSaturated        "DTFIND_FRACTION_SATURATED_SPOTS"

// Cindex Strings
#define D_K_DtindexOptions                 "DTINDEX_OPTIONS"

// Cintegrate Strings
#define D_K_IntegratenDataFlag             "nDataFlag"
#define D_K_IntegratenErrorFlag            "nErrorFlag"
#define D_K_IntegratefDataAddr             "fDataAddr"
#define D_K_IntegratefWidth0               "fWidth0"
#define D_K_IntegratefWidth1               "fWidth1"
#define D_K_IntegratefSvec0                "fSvec0"
#define D_K_IntegratefSvec1                "fSvec1"
#define D_K_IntegratefSvec2                "fSvec2"
#define D_K_IntegratefS0vec0               "fS0vec0"
#define D_K_IntegratefS0vec1               "fS0vec1"
#define D_K_IntegratefS0vec2               "fS0vec2"
#define D_K_IntegrateOblique               "DTINTEGRATE_OBLIQUE"

// Cprofit Strings
#define D_K_ProfitnBadFlag                 "nBadFlag"
#define D_K_ProfitfOtherInt                "fOtherInt"
#define D_K_ProfitfOtherSig                "fOtherSig"

// Project info
#define D_K_ProjectTitle                   "PROJECT_TITLE"

// Optics info
#define D_K_OpticsType                     "OPTICS_TYPE"
#define D_K_OpticsAngle                    "OPTICS_ANGLE"
#define D_K_OpticsCollimator               "OPTICS_COLLIMATOR"
#define D_K_OpticsFilter                   "OPTICS_FILTER"

// Data collection strings (dtlaunch)
#define D_K_DataCollectionDetBinMode       "DATA_COLLECTION_DET_BIN_MODE"
#define D_K_DataCollectionDezinger         "DATA_COLLECTION_DEZINGER"
#define D_K_DataCollectionDistance         "DATA_COLLECTION_DISTANCE"
#define D_K_DataCollectionImageWidth       "DATA_COLLECTION_IMAGE_WIDTH"
#define D_K_DataCollectionTemperature      "DATA_COLLECTION_TEMPERATURE"
#define D_K_DataCollectionTime             "DATA_COLLECTION_TIME"
#define D_K_DataCollection2Theta           "DATA_COLLECTION_TWO_THETA"
#define D_K_DataCollectionWavelength       "DATA_COLLECTION_WAVELENGTH"

// Indexing wedge strings (dtlaunch)
#define D_K_IndexWedges                    "INDEX_WEDGES"
#define D_K_IndexImagesPerWedge            "INDEX_IMAGES_PER_WEDGE"
#define D_K_IndexTimePerImage              "INDEX_TIME_PER_IMAGE"
#define D_K_IndexImageWidth                "INDEX_IMAGE_WIDTH"

// backward compatability for earlier dtlaunch and teXsan for Windows95
#define D_K_DataCollectionFrameWidth       "DATA_COLLECTION_FRAME_WIDTH"
#define D_K_IndexFramesPerWedge            "INDEX_FRAMES_PER_WEDGE"
#define D_K_IndexTimePerFrame              "INDEX_TIME_PER_FRAME"
#define D_K_IndexFrameWidth                "INDEX_FRAME_WIDTH"

// Miscellaneous (detector) options
#define D_K_BinModeOption                "-bin_mode"
#define D_K_DezingerOption               "-dezinger"
#define D_K_DezingerMultiplierOption     "-dezinger_multiplier"
#define D_K_ImageKindOption              "-image_kind"
#define D_K_ReadoutSpeedOption           "-readout_speed"
#define D_K_PixelSizeOption              "-pixel_size"
#define D_K_SubtractClosedOption         "-subtract_closed"
#define D_K_SaveRawOption                "-save_raw"
#define D_K_RawFilenameOption            "-raw_filename"
#define D_K_TriggeredExposureOption      "-trigger"

//Collision 
#define D_K_CollisionInfo                "COLLISION_INFO"
#define D_K_CollisionInfoFilePath        "DTREK_COLLISION_INFO_FILE_PATH"

// Stress keywords
#define D_K_StressNumOfZPositions        "STRESS_NUM_Z_POSITIONS"
#define D_K_StressZPositions             "STRESS_Z_POSITIONS"

#define D_K_StressSlitWidth              "STRESS_SLIT_WIDTH"
#define D_K_StressInclinationType        "STRESS_INCLINATION_TYPE"

#define D_K_StressPsiAngles              "STRESS_PSI_ANGLES"

#define D_K_StressPeakAngle              "STRESS_PEAK_ANGLE"
#define D_K_StressYoungModulus           "STRESS_YOUNG_MODULUS"
#define D_K_StressPoissonRatio           "STRESS_POISSON_RATIO"
#define D_K_StressStressConstant         "STRESS_STRESS_CONSTANT"

#define D_K_ExposureTimeEvaluationInfo   "EXPOSURE_TIME_EVALUATION_INFO"
#define D_K_ProposedExposureTime         "PROPOSED_EXPOSURE_TIME"

typedef unsigned short DTREK_WORD;
typedef unsigned long  DTREK_DWORD;

#ifdef EXPORT_DTREK_DLL
#define DTREK_EXPORT _declspec(dllexport)
#elif defined IMPORT_DTREK_DLL
#define DTREK_EXPORT _declspec(dllimport)
#else
#define DTREK_EXPORT
#endif

#ifdef DTREK_WIN_EXE
    #define DTREK_WIN_DLL_DATA_EXPORT _declspec (dllimport)
#else
    #define DTREK_WIN_DLL_DATA_EXPORT
#endif // DTREK_WIN_EXE

#endif  //!DT_DTREKDEFS_H

