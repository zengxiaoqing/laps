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
#include <config.h>

/*
LW Comment out ifdefs for sys/time.h and time.h below....just include time.h
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifndef TM_IN_SYS_TIME
#include <time.h>
#endif
*/
#include <time.h>

#include <stdio.h>

#ifdef FORTRANUNDERSCORE
#define i4time_now_gg i4time_now_gg_
#endif
#ifdef FORTRANDOUBLEUNDERSCORE
#define i4time_now_gg i4time_now_gg__
#endif
#ifdef FORTRANCAPS
#define i4time_now_gg I4TIME_NOW_GG
#endif

#if defined(__alpha) || defined(__IP21)
	int i4time_now_gg ()
#else
	long i4time_now_gg ()
#endif
{
#if defined(__alpha) || defined(__IP21)
	int i_i4time;
#else
	long i_i4time;
#endif
	time_t utime;
	
	utime = time((time_t)0);

#if defined(__alpha) || defined(__IP21)
	i_i4time = (int)utime;
#else
	i_i4time = (long)utime;
#endif
	i_i4time += 315619200; /* Convert i4time to Jan 1 1960 reference */
	return i_i4time;
}
