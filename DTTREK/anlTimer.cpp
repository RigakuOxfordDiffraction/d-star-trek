// RESTRICTED RIGHTS NOTICE SHORT FORM (JUNE 1987)
//
// Use, reproduction, or disclosure is subject to restrictions set
// forth in Contract No. W-31-109-ENG-38 and Contract No. 
// 943072401 with the University of Chicago, Operator of
// Argonne National Laboratory.
//
// anlTimer.cc           Initial author: M.W.Westbrook et al.
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

#if !defined(VC8) && !defined(__APPLE__)
    #include <iostream.h>
    #include <iomanip.h>
#else
    #include <iostream>
    #include <iomanip>
#endif

#include "anlTimer.h"

int anl_processor = 0;

/* to reduce overhead of local variables */

struct tms stop;
struct tms start;
struct timeval wallclock_start;
struct timeval wallclock_stop;

struct time_info timer[MAXIMUM_TIMERS];

void anl_set_processor (int index) {
  anl_processor = index;
}
 
/* ------------------------------------------------------------ */

void anl_add_timer (int no,int a,int b,char * description) {

 timer[no].timer_no=no;
 strcpy(timer[no].description,description);

 timer[no].start_usr     = timer[a].start_usr   + timer[b].start_usr;
 timer[no].start_cpu     = timer[a].start_cpu   + timer[b].start_cpu;
 timer[no].start_sys     = timer[a].start_sys   + timer[b].start_sys;
 timer[no].start_elapsed = timer[a].start_elapsed + timer[b].start_elapsed;
 timer[no].stop_usr      = timer[a].stop_usr    + timer[b].stop_usr;
 timer[no].stop_cpu      = timer[a].stop_cpu    + timer[b].stop_cpu;
 timer[no].stop_sys      = timer[a].stop_sys    + timer[b].stop_sys;
 timer[no].stop_elapsed  = timer[a].stop_elapsed+ timer[b].stop_elapsed;
 timer[no].accu_usr      = timer[a].accu_usr    + timer[b].accu_usr;
 timer[no].accu_cpu      = timer[a].accu_cpu    + timer[b].accu_cpu;
 timer[no].accu_sys      = timer[a].accu_sys    + timer[b].accu_sys;
 timer[no].accu_elapsed  = timer[a].accu_elapsed+ timer[b].accu_elapsed;
 timer[no].min_usr       = timer[a].min_usr     + timer[b].min_usr;
 timer[no].min_cpu       = timer[a].min_cpu     + timer[b].min_cpu;
 timer[no].min_sys       = timer[a].min_sys     + timer[b].min_sys;
 timer[no].min_elapsed   = timer[a].min_elapsed + timer[b].min_elapsed;
 timer[no].max_usr       = timer[a].max_usr     + timer[b].max_usr;
 timer[no].max_cpu       = timer[a].max_cpu     + timer[b].max_cpu;
 timer[no].max_sys       = timer[a].max_sys     + timer[b].max_sys;
 timer[no].max_elapsed   = timer[a].max_elapsed + timer[b].max_elapsed;
}

/* ------------------------------------------------------------ */
void anl_reset_timer (int no,char *description) {

  timer[no].timer_no = no;
  timer[no].counter = 0;
  
  strcpy(timer[no].description,description);
  
  timer[no].start_usr = 0.0;
  timer[no].start_sys = 0.0;
  timer[no].start_cpu = 0.0;
  timer[no].start_elapsed= 0.0;
  
  timer[no].stop_usr = 0.0;
  timer[no].stop_sys = 0.0;
  timer[no].stop_cpu = 0.0;
  timer[no].stop_elapsed= 0.0;
  
  timer[no].accu_usr = 0.0;
  timer[no].accu_sys = 0.0;
  timer[no].accu_cpu = 0.0;
  timer[no].accu_elapsed= 0.0;

  timer[no].min_usr = 0.0;
  timer[no].min_sys = 0.0;
  timer[no].min_cpu = 0.0;
  timer[no].min_elapsed= 0.0;
  
  timer[no].max_usr = 0.0;
  timer[no].max_sys = 0.0;
  timer[no].max_cpu = 0.0;
  timer[no].max_elapsed= 0.0;
}
/* ------------------------------------------------------------ */

/* ------------------------------------------------------------ */
float anl_get_accu_cpu_time (int no) {  return timer[no].accu_cpu; }
/* ------------------------------------------------------------ */
float anl_get_accu_sys_timer (int no) {  return timer[no].accu_sys; }
/* ------------------------------------------------------------ */
float anl_get_accu_usr_timer (int no) {  return timer[no].accu_usr; }
/* ------------------------------------------------------------ */

/* ------------------------------------------------------------ */

float anl_last_cpu_time (int no)  { 
  float cpu_tm;
  cpu_tm = timer[no].stop_cpu-timer[no].start_cpu;
  return cpu_tm;

}

/* ------------------------------------------------------------ */

float anl_last_sys_time(int no) {
 float sys_tm;
 sys_tm = timer[no].stop_sys-timer[no].start_sys;
 return sys_tm;
}

/* ------------------------------------------------------------ */

float anl_last_usr_time (int no) {
 float usr_tm;
 usr_tm = timer[no].stop_usr-timer[no].start_usr;
 return usr_tm;
}

/* ------------------------------------------------------------ */

float anl_last_elapsed_time (int no) {
 float elapsed_tm;
 elapsed_tm = timer[no].stop_elapsed-timer[no].start_elapsed;
 return elapsed_tm;
}

/* ------------------------------------------------------------ */

void anl_start_timer (int no) { 

  gettimeofday(&wallclock_start);
  timer[no].start_elapsed = (float)wallclock_start.tv_sec;

  times(&start);
  timer[no].start_usr = (float)start.tms_utime/TICKS;
  timer[no].start_sys = (float)start.tms_stime/TICKS;
  timer[no].start_cpu = 
    timer[no].start_usr + timer[no].start_sys; 

}

/* ------------------------------------------------------------ */

void anl_stop_timer (int no) {

 float  Elapsed,Usr,Cpu,Sys;


 times(&stop);
 gettimeofday(&wallclock_stop);

 timer[no].counter++;
/* 
cout << "(stop_timer)TIMER[" << no << "].counter= " << timer[no].counter << "\n";
*/
 timer[no].stop_usr = (float)stop.tms_utime/TICKS;
 timer[no].stop_sys = (float)stop.tms_stime/TICKS;
 timer[no].stop_cpu = timer[no].stop_usr + timer[no].stop_sys; 
 timer[no].stop_elapsed = (float)wallclock_stop.tv_sec;


 Usr     = anl_last_usr_time(no);
 Cpu     = anl_last_cpu_time(no);
 Sys     = anl_last_sys_time(no);
 Elapsed = anl_last_elapsed_time(no);

 timer[no].accu_usr = timer[no].accu_usr + Usr;
 timer[no].accu_sys = timer[no].accu_sys + Sys;
 timer[no].accu_cpu = timer[no].accu_cpu + Cpu;
 timer[no].accu_elapsed= timer[no].accu_elapsed + Elapsed;

 if(timer[no].counter==1){

   timer[no].max_usr = timer[no].min_usr = Usr;
   timer[no].max_cpu = timer[no].min_cpu = Cpu;
   timer[no].max_sys = timer[no].min_sys = Sys;

 } else {

   /* set minimum */

   if (timer[no].min_usr > Usr) {
     timer[no].min_usr = Usr;
   }
   if (timer[no].min_cpu > Cpu) {
     timer[no].min_usr=Cpu;
   }
   if (timer[no].min_sys > Sys) {
    timer[no].min_sys=Sys;
   }
   if (timer[no].min_elapsed > Elapsed) {
    timer[no].min_elapsed = Elapsed;
   }

   /* set maximum */

   if (timer[no].max_usr < Usr) {
     timer[no].max_usr = Usr;
   }
   if (timer[no].max_cpu < Cpu) {
     timer[no].max_usr = Cpu;
   }
   if (timer[no].max_sys < Sys) {
    timer[no].max_sys = Sys;
   }
   if (timer[no].max_elapsed < Elapsed) {
    timer[no].max_elapsed = Elapsed;
   }
   
 } /* end for else */

 }


/* ------------------------------------------------------------ */

void anl_print_accu_timer(int no) {

  cout << "[" << anl_processor << "] " 
       << "TIMER[" << no << "] ACC: "
       << timer[no].description << ":"
       << " " << CPU_LABEL << timer[no].accu_cpu
       << " " << USR_LABEL << timer[no].accu_usr
       << " " << SYS_LABEL << timer[no].accu_sys
       << " " << ELAPSED_LABEL << timer[no].accu_elapsed
       << "\n";

}

/* ------------------------------------------------------------ */

void anl_print_last_interval(int no) {

  cout << "[" << anl_processor << "] " 
       << "TIMER[" << no << "] LST: "
       << timer[no].description << ":"
       << " " << CPU_LABEL<< timer[no].stop_cpu-timer[no].start_cpu
       << " " << USR_LABEL<< timer[no].stop_usr-timer[no].start_usr
       << " " << SYS_LABEL<< timer[no].stop_sys-timer[no].start_sys
       << " " << ELAPSED_LABEL
       << timer[no].stop_elapsed-timer[no].start_elapsed
       << "\n";
}

/* ------------------------------------------------------------ */

void anl_timer_info() {
  
  cout 
    << "/*****************************************************/\n"
    <<  " The following are the descriptions about Timer Library\n"

    <<  " [USR] User time is the time spent "
    <<  "by the program running on user mode.\n"

    <<  " [SYS] System time is the time spent by"
    <<  "by the program runnng on system mode.\n"

    <<  " [CPU] CPU time is the sum of user time and system time.\n"
    <<  "/****************************************************/\n";

}

/* ------------------------------------------------------------ */

void anl_print_timer (int no) {

  cout << "\n";

  cout << "[" << anl_processor << "] " 
       << "TIMER[" << no << "] ACC: "
       << timer[no].description << ":"
       << " " << CPU_LABEL << timer[no].accu_cpu
       << " " << USR_LABEL << timer[no].accu_usr
       << " " << SYS_LABEL << timer[no].accu_sys
       << " " << ELAPSED_LABEL << timer[no].accu_elapsed
       << "\n";

  cout << "[" << anl_processor << "] " 
       << "TIMER[" << no << "] MIN: "
       << timer[no].description << ":"
       << " " << CPU_LABEL << timer[no].min_cpu
       << " " << USR_LABEL << timer[no].min_usr
       << " " << SYS_LABEL << timer[no].min_sys
       << " " << ELAPSED_LABEL << timer[no].min_elapsed
       << "\n";

  cout << "[" << anl_processor << "] " 
       << "TIMER[" << no << "] MAX: "
       << timer[no].description << ":"
       << " " << CPU_LABEL << timer[no].max_cpu
       << " " << USR_LABEL << timer[no].max_usr
       << " " << SYS_LABEL << timer[no].max_sys
       << " " << ELAPSED_LABEL << timer[no].max_elapsed
       << "\n";

  if (0 != timer[no].counter)
    {
      cout << "[" << anl_processor << "] " 
           << "TIMER[" << no << "] AVE: "
           << timer[no].description << ":"
	   << " " << CPU_LABEL << timer[no].accu_cpu/timer[no].counter
           << " " << USR_LABEL << timer[no].accu_usr/timer[no].counter
           << " " << SYS_LABEL << timer[no].accu_sys/timer[no].counter
           << " " << ELAPSED_LABEL << timer[no].accu_elapsed/timer[no].counter
           << "\n";
    }

/*
  cout << "(print_timer)TIMER[" << no << "].counter= " << timer[no].counter << "\n";
*/
}

/* ------------------------------------------------------------ */
