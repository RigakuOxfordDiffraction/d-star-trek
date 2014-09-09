#ifndef DT_ANLTIMER_H
#define DT_ANLTIMER_H

// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// anlTimer.h           Initial author: M.W.Westbrook et al.
//    This file contains the functions for timer routines.
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

//+Include files

#include<sys/times.h>
#include<sys/time.h>

//+External prototypes

//+Definitions and constants

#define TICKS 100
#define MAXIMUM_TIMERS  255

#define CPU_LABEL "CPU="
#define SYS_LABEL "SYS="
#define USR_LABEL "USR="
#define ELAPSED_LABEL "ELA="

struct time_info {
       int timer_no;
       int counter;
       char description[50];
       float start_usr;
       float start_sys;
       float start_cpu;
       float start_elapsed;
       float stop_usr;
       float stop_sys;
       float stop_cpu;
       float stop_elapsed;
       float accu_usr;
       float accu_sys;
       float accu_cpu;
       float accu_elapsed;
       float min_usr;
       float min_sys;
       float min_cpu;
       float min_elapsed;
       float max_usr;
       float max_sys;
       float max_cpu;
       float max_elapsed;
};

//+Code begin

void anl_set_processor (int index);

void anl_reset_timer(int no, char * description);
void anl_start_timer(int timer_no);
void anl_stop_timer(int timer_no); 

void anl_print_timer(int timer_no); 
void anl_print_timer(int no);
void anl_print_last_interval(int timer_no); 

void anl_print_accu_timer(int timer_no); 
void anl_print_timer_info();


void anl_timer_add(int no,int a,int b,char *description);

float anl_last_cpu_time (int no);
float anl_last_sys_time (int no);
float anl_last__usr_time (int no);
float anl_last_elapsed_time (int no);

float anl_get_accu_cpu_time(int no);
float anl_get_accu_sys_time(int no);
float anl_get_accu_usr_time(int no);

#endif   // DT_ANLTIMER_H
