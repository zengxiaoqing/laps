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
 *	Read_88d_a2_ ...
 *
 *	Read a realtime A2 data record and pass it back to the user to run
 *	in his/her algorithm.  All the current calling programs are in
 *	Fortran.
 */

#ifndef	LINT
static char read_88d_a2_id[] = "@(#)read88da2.c	2.2	2/14/95";
#endif	/* LINT */

#include <config.h>
#include <stdio.h>
#include <strings.h>
#include "merge.h"

#define	MAX_DATA	920
#define	NSIZE		MAX_DATA * 2
#define	ALG_MISS	999.
#define	REFL_MIN	-25.
#define	DBZ		0
#define	VEL		1

#define	DBZ_PRESENT	get_status(DBZ)
#define	VEL_PRESENT	get_status(VEL)

static int old_scan = -1;
static int new_scan;
static int endflag = 0;

read_88d_a2_(azm, elev, swp_nbr, ndate, ntime, ng_dbz, ng_vel, unyq, refl,
	vel, uendv, idata_count, radial_nbr, sw, vcp, ibegt, iendt,
	first_rad, rtortp, rth_vel, scan_number, gs_dbz, gs_vel )
float	*azm;				/* azimuth angle */
float	*elev;				/* elevation angle */
long	*swp_nbr;			/* tilt number */
long	*ndate;				/* date info */
long	*ntime;				/* time info */
long	*ng_dbz;			/* number of reflectivity gates */
long	*ng_vel;			/* number of velocity/sw gates */
float	*unyq;				/* nyquist velocity */
short	refl[MAX_DATA/2];		/* reflectivity  */	
float	vel[MAX_DATA];			/* velocity data */
long	*uendv;				/* end of volume scan flag */
long	*idata_count;			/* ignored */
long	*radial_nbr;			/* radial number, ignored */
short	sw[MAX_DATA];			/* sw data (short?) */
long	*vcp;				/* VCP number */
long	*ibegt;				/* beginning time, ignored */
long	*iendt;				/* ending time, ignored */
long	*first_rad;			/* first radial, ignored */
long	*rtortp;			/* 0 for playback, 1 for RT */
long	*rth_vel;			/* reflectivity threshold */
long	*scan_number;			/* volume scan number */
float   *gs_dbz;			/* dbz gate spacing in km */
float   *gs_vel;			/* vel gate spacing in km */
{
	float tmp_dbz[NSIZE], tmp_vel[NSIZE], tmp_snr[NSIZE], tmp_sw[NSIZE];
	unsigned char c_dbz[NSIZE];
	int i;
	int k;
	int merging;

	if ( endflag )
	{
		*uendv = 1;
		endflag = 0;
		return;
	}

	*azm = (float) get_azi() / 100.0;
	*elev = (float) get_elev() / 100.0;
	*swp_nbr = get_tilt();
	*ndate = 10000 * get_month() + 100 * get_day() + get_year();
	*ntime = 10000 * get_hour() + 100 * get_min() + get_sec();
	*ng_dbz = get_number_of_gates(DBZ);
	*ng_vel = get_number_of_gates(VEL);
	*unyq = (float) get_nyquist() / 100.0;
	*vcp = get_vcp();
/*
 *	DON'T USE THIS. SSAMAIN SETS "RTORTP" AND PASSES THE INFO TO US!
	*rtortp = get_rt_mode();
 */
	*scan_number = get_scan();
	*gs_dbz = ( float ) get_gate_spacing( DBZ ) / 1000.0;
	*gs_vel = ( float ) get_gate_spacing( VEL ) / 1000.0;
/*
 *	The next three variables are not important.  They are set to zero.
 */
	*idata_count = 0;
	*radial_nbr = 0;
	*first_rad = 0;

	get_data(tmp_dbz, tmp_vel, tmp_snr, tmp_sw, NSIZE);

/*
 *	Save reflectivity data or retrieve reflectivity data, as needed.
 *
 *	Yes, I know this is *VERY* ugly because this is interfacing to routines
 *	that were designed for LL code.
 */

	merging = 0;

	if ( !VEL_PRESENT )
	{
		free_dbz_radials();
/*
 *	Scale floating point to character.
 */
		dataload8( tmp_dbz, c_dbz, DZ_MIN, DZ_RES, *ng_dbz);
		alloc_dbz_radial( *ng_dbz , get_first_gate( DBZ ),
			get_gate_spacing( DBZ ), get_azi(), c_dbz );
	}
	else
	{
		if ( DBZ_PRESENT )
			free_dbz_radials();
		else
		{
			restore_dbz_radial( *ng_vel, get_first_gate( VEL ),
				get_gate_spacing( VEL ), get_azi(), c_dbz );
/*
 *	Scale character to floating point.  (See, I told you it was ugly).
 */
			dataunload8( tmp_dbz, c_dbz, DZ_MIN, DZ_RES, *ng_vel );
			merging = 1;
			*ng_dbz = *ng_vel;
		}
	}

/*
 *	Handle reflectivity.  Algorithms want raw data, so that means we
 *	have to do the following:
 *
 *		(1)  Threshold the raw data.
 *		(2)  Pass data as it arrives from tilts 1 and 3.
 *		(3)  Delete data for tilts 2 and 4.
 *		(4)  Convert data from float to short integers
 */

/*
 *	Threshold the data fields.
 *
 *	Note:  this code assumes that DATA_MISSING and RANGE_FOLDED_DATA are
 *	       less than -700.0.
 */

	if (*ng_dbz > 0)
	{
		for(i=0; i<*ng_dbz; i++)
		{
			if (tmp_dbz[i] < -700.0)
				tmp_dbz[i] = ALG_MISS;
			if (tmp_dbz[i] >= 777.0)
				tmp_dbz[i] = REFL_MIN;
		}
	}

	if (*ng_vel > 0)
	{
		for(i=0; i<*ng_vel; i++)
		{
			if (tmp_vel[i] < -700.0)
				tmp_vel[i] =  ALG_MISS;
			if (tmp_sw[i] < -700.0)
				tmp_sw[i] = ALG_MISS;
/*
 *	Threshold using reflectivity.  Keep in mind that the gate spacing
 *	for reflectivity is 1km, while the gate spacing for velocity is
 *	0.25km.
 *
 *	New note:  the above applies only when "merging" is not set.  If
 *	"merging" is set, then the data is 0.25km spacing.
 */
			if ( merging == 0 )
				k = i / 4;
			else
				k = i;
			if (*ng_dbz > 0 && tmp_dbz[k] < *rth_vel)
			{
				tmp_vel[i] = ALG_MISS;
				tmp_sw[i] = ALG_MISS;
			}
		}
	}

/*
 *	Handle reflectivity.
 */

	if ( *ng_dbz > 0 && !merging )
	{
		if (*ng_dbz > MAX_DATA/2)
			*ng_dbz = MAX_DATA/2;
		for(i=0; i<*ng_dbz; i++)
			refl[i] = (short) tmp_dbz[i];
	}

/*
 *	Handle velocity/sw data.  Just copy the data.
 */

	if (*ng_vel > 0)
	{
		if (*ng_vel > MAX_DATA)
			*ng_vel = MAX_DATA;
		for(i=0; i<*ng_vel; i++)
		{
			vel[i] = tmp_vel[i];
			sw[i] = (short) tmp_sw[i];
		}
	}

	read_radial();

	new_scan = get_scan();
	if (old_scan == -1)
		old_scan = new_scan;
	if (old_scan != new_scan)
		endflag = 1;
	*uendv = 0;

	old_scan = new_scan;

/*
 *	Restore dbz status.
 */

	if ( merging )
		*ng_dbz = 0;

	return;
}

rt_alg_init_( s )
char *s;
{
	char *ptr;
/*
 *	Check for white space and stick in a NULL.
 */
	ptr = index( s, ' ');
	if (ptr != NULL)
		*ptr = NULL;
	radar_init( s );
	read_radial();
}

rt_set_radar_name_( s )
char *s;
{
	set_radar_name( s );
}
