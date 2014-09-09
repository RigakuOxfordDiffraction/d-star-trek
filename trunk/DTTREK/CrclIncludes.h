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

#ifdef EXPORT_DTREK_DLL
#define DTREK_EXPORT _declspec(dllexport)
#elif defined IMPORT_DTREK_DLL
#define DTREK_EXPORT _declspec(dllimport)
#else
#define DTREK_EXPORT
#endif

DTREK_EXPORT int dtfind_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtindex_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtrefine_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtpredict_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtorient_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtstrategy_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtmultistrategy_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtintegrate_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtprofit_main(int argc, char *argv[]);
DTREK_EXPORT int dtcell_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtaxial_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtrefln2numabs_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtnumabs2refln_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtreflnmerge_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtscaleaverage_main(unsigned int argc, char *argv[]);
DTREK_EXPORT int dtranker_main(unsigned int argc, char *argv[]);

