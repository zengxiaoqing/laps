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
#include <stdio.h>
#include <stdlib.h>

#ifdef FORTRANUNDERSCORE
#define fstrncpy fstrncpy_
#define nstrncpy nstrncpy_
#endif
#ifdef FORTRANDOUBLEUNDERSCORE
#define fstrncpy fstrncpy__
#define nstrncpy nstrncpy__
#endif
#ifdef FORTRANCAPS
#define fstrncpy FSTRNCPY
#define nstrncpy NSTRNCPY
#endif

/*************************************************************************
*       fstrncpy  -  copy function used to copy strings with embedded blanks
*************************************************************************/
#ifdef __STDC__
void fstrncpy (char *target,char *source,long maxlen)
#else
void fstrncpy (target,source,maxlen)
char *target;           /* space to be copied into */
char *source;           /* string to be copied     */
long maxlen;            /* max length of *source   */
#endif
{
        while (maxlen-- && *source != '\0')
          *target++ = *source++;
        *target = '\0';
}
 
/*************************************************************************
*       nstrncpy  -  copy function used to copy strings terminated with blanks
*************************************************************************/
#ifdef __STDC__
void nstrncpy (char *target,char *source,long maxlen)
#else
void nstrncpy (target,source,maxlen)
char *target;           /* space to be copied into */
char *source;           /* string to be copied     */
long maxlen;            /* max length of *source   */
#endif
{
        while (maxlen-- && *source != ' ')
          *target++ = *source++;
        *target = '\0';
}
