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
/*      waiting_c (iseconds)



This routine will pause execution for the requested number of seconds.  It
is designed to be FORTRAN callable.

Author:  Dan Birkenheuer
Date:    5/17/95

From FORTRAN use:



integer isecs

isecs = 10 ! or some number of seconds.
call waiting_c (isecs)



This routine returns no parameters.  If the user sends a negative value for
isecs, the wait will be for the absolute value of isecs.


*/


#include <config.h>

#include <unistd.h>


#ifdef FORTRANUNDERSCORE
#define waiting_c waiting_c_
#endif

#ifdef FORTRANDOUBLEUNDERSCORE
#define waiting_c waiting_c__
#endif

#ifdef __STDC__
void waiting_c ( int *isec)
#else
void waiting_c (isec)
int *isec;
#endif
{
int istatus;
unsigned int  secs_to_wait;
secs_to_wait = abs(*isec);
istatus = sleep (secs_to_wait);
return;
}


