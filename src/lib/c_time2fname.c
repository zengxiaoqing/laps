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
	Subroutine c_time2fname (Ictime, c_time)

This routine is designed to be callable from FORTRAN as illustrated above. 
The input value is an integer*4 type variable (Ictime) which is the
conventional c-time string.  The routine computes a CHARACTER*9 return string
(c_time) that contains the filename in the form YYJJJHHMM.

Author:   Dan Birkenheuer
Last modified: 3/24/94

*/

#include <stdio.h>
#include <time.h>

#ifdef FORTRANUNDERSCORE
#define c_time2fname  c_time2fname_
#endif
#ifdef FORTRANDOUBLEUNDERSCORE
#define c_time2fname  c_time2fname__
#endif
#ifdef FORTRANCAPS
#define c_time2fname  C_TIME2FNAME
#endif

#ifdef __STDC__
void c_time2fname (long *ictime, char *c_time)
#else
void c_time2fname (ictime, c_time)
long *ictime;
char *c_time;
#endif


{   
    struct tm *timePtr;
    time_t valTime;

    /* doesn't do any input checks */
    valTime = (time_t) *ictime;

    timePtr = gmtime(&valTime);

    sprintf(c_time,"%02d%03d%02d%02d", timePtr->tm_year,
	timePtr->tm_yday+1,
        timePtr->tm_hour, timePtr->tm_min);

	return;

}
