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
#ifndef _fort2c_str_
#define _fort2c_str_
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


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
/***********************************************************************/
#ifdef __STDC__
void downcase_c(char *instr,char *outstr)
#else
void downcase_c(instr,outstr)
char *instr;
char *outstr;
#endif
{
        char *inptr, *tempout;
        short i,len_in;
        char outchar[2];

        outchar[1] = '\0';
	len_in = strlen(instr);
        tempout = (char *) malloc(len_in + 1); 
        tempout[0] = '\0';
        tempout[1] = '\0';
        inptr = instr;
 
        for (i = 0; i < len_in; i++) {
          outchar[0] = tolower(*inptr);
          if (i == 0)
            strncpy(tempout,outchar,1);
          else
            strncat(tempout,outchar,1);
          inptr++;
        }
        
        strcpy(outstr,tempout);
        free(tempout);
}
/***********************************************************************/
#ifdef __STDC__
void upcase_c(char *instr,char *outstr)
#else
void upcase_c(instr,outstr)
char *instr;
char *outstr;
#endif
{
        char *inptr, *tempout;
        short i,len_in;
        char outchar[2];

        outchar[1] = '\0';
	len_in = strlen(instr);
        tempout = (char *) malloc(len_in + 1); 
        tempout[0] = '\0';
        tempout[1] = '\0';
        inptr = instr;
 
        for (i = 0; i < len_in; i++) {
          outchar[0] = toupper(*inptr);
          if (i == 0)
            strncpy(tempout,outchar,1);
          else
            strncat(tempout,outchar,1);
          inptr++;
        }
        
        strcpy(outstr,tempout);
        free(tempout);
}
/***********************************************************************/
#endif
