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
 *	Merge routines for SSAP.  UNTESTED!
 */

#ifndef	LINT
static char merge[] = "@(#)merge.c	2.2	2/14/95";
#endif	/* LINT */

#include <config.h>
#include <stdio.h>
#include "defs.h"
#include "merge.h"

dataload8(rt_data, c_data, min, res, ng)
float rt_data[];
unsigned char c_data[];
int min, res, ng;
{
	unsigned char *ucptr;
	int i, max, value;

	max = min + res * 255;

	ucptr = c_data;

#ifdef	XXX
	for(i=0; i<NUM_DATA; i++)
#endif
	for(i=0; i<ng; i++)
	{
		if (i < ng)
		{
			value = (int) (rt_data[i] * 100);
			if (value < min || value > max)
				*ucptr++ = MISSING_VALUE;
			else
				*ucptr++ = (value - min) / res;
		}
		else
			*ucptr++ = MISSING_VALUE;
	}
	return;
}

dataunload8(rt_data, c_data, min, res, ng)
float rt_data[];
unsigned char c_data[];
int min, res, ng;
{
	int i;

	for(i=0; i<ng; i++)
	{
		if ( c_data[i] == 0 )
			rt_data[i] = ALG_MISS;
		else
			rt_data[i] = (float) 
				( (short) c_data[i] * res + min ) * 0.01;
	}

	return;
}

alloc_dbz_radial(ng, fg, gs, az, data)
int ng, fg, gs, az;
unsigned char data[];
{
	unsigned char *ptr, *dptr;
	unsigned int size;
	int i;

	if (rt_num_radials >= rt_current_alloc)
	{
		size = (rt_current_alloc + INCREMENT) * 4;
		if (rt_current_alloc == 0)
			mptr = (unsigned char **) malloc(size);
		else
			mptr = (unsigned char **) realloc(mptr,size);
		if (mptr == NULL)
		{
			perror("alloc_dbz_radial:  malloc/realloc");
			return;
		}
		rt_current_alloc += INCREMENT;
/*
		printf("rt_current_alloc %d\n",rt_current_alloc);
		fflush(stdout);
 */
	}
	ptr = (unsigned char *) malloc(dz_size + ng);

	if (ptr == NULL)
	{
		perror("alloc_dbz_radial:  malloc");
		return;
	}
/*
 *	Save the data.
 */

	dz_radial = (struct dz_radial *) ptr;
	dz_radial->first_gate = fg;
	dz_radial->gate_spacing = gs;
	dz_radial->number_of_gates = ng;
	dz_radial->azimuth = az;

	dptr = &ptr[dz_size];
	for(i=0; i<ng; i++)
		*dptr++ = data[i];

	mptr[rt_num_radials] = ptr;
	rt_num_radials++;
	return;
}

free_dbz_radials()
{
	int i, n ;

	if (rt_using_dbz)
	{
/*
		printf("freeing up radials\n");
 */
		rt_using_dbz = 0;
	}
	else
		return;

	if (rt_num_radials == 0)
		return;

	for(i=0; i<rt_num_radials; i++)
	{
		n = cfree(mptr[i]);
		if (n == 0)
			perror("free_dbz_radials:  free");
	}
	rt_num_radials = 0;
	return;
}

restore_dbz_radial(ng, fg, gs, az, data)
int ng, fg, gs, az;
unsigned char data[];
{
	int i, k, kstart;
	int del, delold, kdel;
	int myrange, hisrange;
	int max;
	unsigned char *ptr;
	unsigned char *dptr;
	unsigned char *eptr;

	rt_using_dbz = 1;

/*
 *	Initialize everything to the missing value.
 */

	ptr = data;
	for(i=0; i<ng; i++)
		*ptr++ = MISSING_VALUE;

	kdel = -1;
	delold = 999;

/*
 *	Find the closest radial.
 */

	for(i=0; i<rt_num_radials; i++)
	{
		dz_radial = (struct dz_radial *) mptr[i];

		del = az - dz_radial->azimuth;
		if (del < 0)
			del = -del;
		if (del > 30000)
			del = 36000 - del;

		if (del < delold)
		{
			delold = del;
			kdel = i;
		}
	}

/*
 *	Is the nearest radial close enough?
 */

	if (delold > AZI_ACCEPT)
	{
#ifdef	DEBUG_0
		printf("restore_dbz_radial: failed to find a radial for %d\n",az);
#endif	/* DEBUG_0 */
		return;
	}

/*
 *	It is.
 */

	dz_radial = (struct dz_radial *) mptr[kdel];

/*
 *	Check to see if everything looks ok!
 */
	if ( (dz_radial->gate_spacing == 4 * gs) && 
		(dz_radial->first_gate -3 * gs == fg) )
	{
		ptr = mptr[kdel];
		dptr = &ptr[dz_size];
		for(i=0; i<dz_radial->number_of_gates; i++)
			data[i] = dptr[i];
		max = dz_radial->number_of_gates * 4;
		if (max > ng)
			max = ng;
		expand(data, ng);
	}
	else
	{
		printf("restore_dbz_radial: test fails!\n");
		if  (dz_radial->gate_spacing == 4 * gs) 
			printf("restore_radial: test 1 ok\n");
		else
			printf("restore_radial: test 1 fail\n");
		if (dz_radial->first_gate - 3 * gs == fg)
			printf("restore_radial: test 2 ok\n");
		else
			printf("restore_radial: test 2 fail\n");
		printf("restore_radial:  method fails!\n");
	}

	return;
}
