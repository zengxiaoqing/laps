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
 *	Precip.c...
 *
 *	This package of routines will be interfaced to the precip algorithm.
 *	It is designed to make tape playback and realtime data look like a
 *	freshly loaded A2 tape.
 */

#ifndef	LINT
static char precip[] = "@(#)precip.c	2.2	2/14/95";
#endif	/* LINT */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "a2io.h"

/*
 *	mode values:
 *		0 ... never called, send ARCHIVE2
 *		1 ... send eof
 *		2 ... send header
 *		3 ... send data
 */

int mode = 0;

#ifdef	RTS
char vs_info[24];
#define	DATA_NUM	2432
#endif	/* RTS */

char buf[MAX_RECORD_SIZE];

#define	WARN_FLAG	500
int warn = 0;

int eof = 0;

int scan_number = -99;

char *ptr;
char *dptr;
struct a2header *a2header;
struct a2data *a2data;
struct a2scan *a2scan;

precip_open()
{
	char *ptr;
#ifndef	RTS

	ptr = getenv("PRECIP_TAPE");
	if ( ptr == NULL )
	{
		printf("You must set the PRECIP_TAPE environment parameter.\n");
		exit(1);
	}

	radar_init( ptr );
#else

	ptr = getenv("RT_MODE");

	if ( ptr == NULL )
		radar_init( "Raw");
	else
		radar_init( ptr );
#endif
	return;
}

precip_read( s )
char *s;
{
	int n;
	int rc;

	if ( mode == 0 )
	{
#ifdef	RTS
/*
 *	Search for a valid volume scan number record.
 */
		for(;;)
		{
			rc = read_record( buf );

			if ( is_vol_scan() )
				break;
			warn++;
			if ( warn >= WARN_FLAG )
			{
				printf("Looking for volume scan record\n");
				warn = 0;
			}	
		}
#endif	/* RTS */
		bcopy( "ARCHIVE2", s, 8 );
		mode++;
		return( 8 );
	}
	else if ( mode == 1 )
	{
		mode++;
		return( 0 );
	}
	else if ( mode == 2 )
	{
#ifndef	RTS
		if ( eof > 1 )
			return( 0 );
		while(( rc = read_record( buf )) != 24 )
		{
			eof_check( rc );
			if ( eof > 1 )
				return( 0 );
		}
		bcopy( buf, s, 24 );
		buf[12] = NULL;
		scan_number = atoi( &buf[9] );
#else
		bcopy( vs_info, s, 24 );
#endif	/* RTS */
		mode++;
		return( 24 );
	}
	else if ( mode == 3 )
	{
#ifndef	RTS
		if ( eof > 1 )
			return( 0 );
reread:
		rc = read_record( buf );
/*
 *	Make sure we didn't get an i/o error.  If we did get an i/o error,
 *	the bcopy() call below would get a segmentation violation.
 *
 *	We'll don't worry about too many i/o errors, as tpread() [called by
 *	read_record()] does that for us.
 */

		if ( rc == -1 )
			goto reread;

		eof_check( rc );

		if ( eof >  0 )
		{
			mode = 2;
			return( 0 );
		}
		else
		{
			bcopy( buf, s, rc );
			return( rc );
		}
#else
		rc = read_record( buf );
		if ( is_vol_scan() )
		{
			mode = 2;
			return( 0 );
		}
		else
		{
			bzero( s, DATA_NUM );
			bcopy( buf, &s[12], rc );
/*
			if (rc + 16 != DATA_NUM)
			{
				printf("rc %d DATA_NUM %d\n",rc,DATA_NUM);
			}
 */
			return( DATA_NUM );
		}
#endif	/* RTS */
	}
	else
	{
		printf("precip_read:  help:  reached impossible code\n");
		return( -1 );
	}
/*NOTREACHED*/
}

#ifndef	RTS
eof_check( n )
int n;
{
	if ( n == 0 )
	{
		eof++;
		if ( eof > 1 )
			printf("precip_read:  double eof detected\n");
	}
	else
		eof = 0;
	return;
}
#endif

#ifdef	RTS

#define	SECONDS_PER_DAY		(24 * 60 * 60)

struct fake_header {
	char s[12];
	long date_info;
	long time_info;
	long dummy;
} fake_header;

is_vol_scan()
{
	long now;

	ptr = buf;
	a2header = ( struct a2header * ) ptr;
	dptr = &buf[ sizeof( struct a2header ) ];

/*
 *	There are differences between realtime data and tape playback data.
 *	The tape data seems to be the correct version.  Make the realtime
 *	stream look like playback.
 */
	if ( a2header->message_type == A2_DATA_TYPE )
	{
		a2data = ( struct a2data *) dptr;
		if ( a2data->radial_status == 3 && 
			a2data->elevation_number != 1 )
				a2data->radial_status = 0;
		a2data->dbz_ptr = a2data->arc_dbz_ptr;
		a2data->vel_ptr = a2data->arc_vel_ptr;
		a2data->spw_ptr = a2data->arc_spw_ptr;
		return( 0 );
	}
	if ( a2header->message_type != A2_VOLSCAN_TYPE )
		return( 0 );

	a2scan = ( struct a2scan * ) dptr;

	scan_number = a2scan->vsn;
	sprintf( fake_header.s,"ARCHIVE2.%.3d", a2scan->vsn );

	now = time( 0 );

	fake_header.date_info = now / SECONDS_PER_DAY;
/*
 *	Change UNIX date to NEXRAD date.
 */
	fake_header.date_info += 1;
	fake_header.time_info = now % SECONDS_PER_DAY;
/*
 *	NEXRAD wants milliseconds.
 */
	fake_header.time_info *= 1000;
	fake_header.dummy     = 0;

	bcopy( (char *) &fake_header, vs_info, 24 );
	return( 1 );
}
#endif	/* RTS */

get_precip_scan_( num )
long int *num;
{
	*num = scan_number;
	return;
}
