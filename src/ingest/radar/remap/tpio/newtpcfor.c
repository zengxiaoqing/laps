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
 *	Newtpcfor...
 *
 *	Fortran to C interface for TPIO package.
 */

#ifndef	LINT
static char newtpcfor_id[] = "@(#)newtpcfor.c	2.1	2/15/95";
#endif	/* LINT */

#include <stdio.h>
#include <malloc.h>
#include "newtpcfor.h"

tpopen_(iunit,dev,mode,iflag)
int *iunit,*mode,*iflag;
char *dev;
{
	int i,n;
	char *ptr, *ptr2, *dev2;

#ifdef	COMPILE_IN_THIS_CODE
#ifdef	SUN
/*
 *	Fix a Sun bug.
 *
 *	If a fortran program sends to stdout, and is piped to a non-existent
 *	program (i.e., a typo in the piped program), then the user may get a
 *	segmentation violation or an open pipe complaint.  If trace(1) is
 *	run, then one of the complaints will occur.  As a workaround, we disable
 *	signal 13 (SIGPIPE).
 */
	signal(13,1);
#endif
#endif

/*
 *	We need to clear out the spaces from the device name, as it seems
 *	to be possible to get them.
 */

	n = strlen(dev) + 1;
	ptr = dev;
	ptr2 = malloc(n);
	if (ptr2 == NULL)
	{
		perror("tpopen:  malloc");
		*iflag = -1;
	}
	bzero(ptr2,n);

	dev2 = ptr2;
	for(i=0; i<n; i++)
	{
		if (*ptr != ' ')
			*ptr2 = *ptr;
		else
			*ptr2 = NULL;
		*ptr++;
		*ptr2++;
	}
	*ptr2 = NULL;

	*iflag = c_tpopen(*iunit,dev2,*mode);
	return;
}

tpread_(iunit,ibuf,ibuflen,lenr)
int *iunit,*ibuflen,*lenr;
char *ibuf;
{
	*lenr = c_tpread(*iunit,ibuf,*ibuflen);
/*
 *	Now, make things compatible with the current version.
 */
	if (*lenr == 0)
		return(1);
	if (*lenr == -1)
	{
		*lenr = 0;
		return(2);
	}
	return(0);
}

tpwrite_(iunit,ibuf,ibuflen,lenr)
int *iunit,*ibuflen,*lenr;
char *ibuf;
{
	*lenr = c_tpwrite(*iunit,ibuf,*ibuflen);
/*
 *	Now, make things compatible with the current version.
 */
	if (*lenr == 0)
		return(1);
	if (*lenr == -1)
	{
		*lenr = 0;
		return(2);
	}
	return(0);
}

tpclose_(iunit)
int *iunit;
{
	c_tpclose(*iunit);
}

tprew_(iunit)
int *iunit;
{
	c_tprew(*iunit);
}

tpbsf_(iunit,n)
int *iunit,*n;
{
	c_tpbsf(*iunit,*n);
}

tpbsr_(iunit,n)
int *iunit,*n;
{
	c_tpbsr(*iunit,*n);
}

tpfsf_(iunit,n)
int *iunit,*n;
{
	c_tpfsf(*iunit,*n);
}

tpfsr_(iunit,n)
int *iunit,*n;
{
	c_tpfsr(*iunit,*n);
}

tpweof_(iunit)
int *iunit;
{
	c_tpweof(*iunit);
}

tprewoff_(iunit)
int *iunit;
{
	c_tprewoff(*iunit);
}

tpstatus_(iunit,rc)
int *iunit,*rc;
{
	*rc = c_tpstatus(*iunit);
}

tpseek_(iunit,k)
int *iunit,*k;
{
	c_tpseek(*iunit,*k);
}
