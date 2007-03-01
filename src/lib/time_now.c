/*cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis  
cdis 
cdis*/
/*

	ROUTINE TIME_NOW IS CAST FOR FORTRAN COMPATIBILITY

	It can be called from Fortran using:

	integer*4 c_time
	call time_now(c_time)

	The routine returns the current system's c_time.


c Author:  Dan Birkenheuer
c Date:     April 1, 1994

*/
#include <config.h>
#include <stdio.h>
#include <stddef.h>

/*
LW commented out ifdefs for sys/time.h and time.h....just include time.h
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifndef TM_IN_SYS_TIME
#include <time.h>
#endif
*/

#include <time.h>

#ifdef FORTRANUNDERSCORE
#define time_now time_now_
#endif

#ifdef FORTRANDOUBLEUNDERSCORE
#define time_now time_now__
#endif

#ifdef __STDC__
void time_now (long *c_time)
#else
void time_now (c_time)
long *c_time;
#endif

{
  unsigned long now;

  now = time(NULL);
  *c_time = (long) now;

	return;
     
}
