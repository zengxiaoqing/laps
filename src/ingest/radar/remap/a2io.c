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
 *	a2io...
 *
 *	Decode A2 data.
 */
 
#ifndef	LINT
static char a2io[] = "%W%	%G%";
#endif	LINT

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include "c_buff.h"
#include "a2io.h"
#include "debug.h"
#include "defs.h"
#include "vcp.h"

static char buf[145920];

#define	A_SCALE		180. / 32768.

/*
 *	I/O mode dependencies.
 */

#ifndef	RTS

#define	REC_SIZE	2432
#define	MAX_FAIL	10
int luin = 0;
static int failure = 0;

#else

#define	REC_SIZE	2416

#define	OLD_SCAN_METHOD

#ifdef	OLD_SCAN_METHOD
static unsigned short old_elevation_number = 999;
#endif	OLD_SCAN_METHOD

struct c_buff *wideband_init_reader();
struct c_buff *my_cbuff;

int pointer;

#endif

#define	SECONDS_PER_DAY		24 * 60 * 60

#define	BEAMWIDTH	0.95

int radar_name_set = 0;
char radar_name[4];
char site_name[8];
int lat, lon;

/*
 *	We need to do our own scan number.
 */

static int scan_number = 0;

static int new = 0;
static int left = 0;

static struct tm *gmt;
char gmt_data[sizeof(struct tm)];

unsigned short int *vcp;

struct a2header *a2header;
struct a2data *a2data;

int a2size;

char *ptr;
char *dptr;

/*
 *	Warning message frequency.
 */

#define	WARN_FREQ	500

int warn_flag = 0;
int scan_flag = 0;

a2_init( s )
char *s;
{

/*
 *	Retrieve radar info.
 */

	get_radar_info();

	a2size = sizeof( struct a2header );

#ifdef	RTS
	
	while( ( my_cbuff = wideband_init_reader( 0x88d2, &pointer )) == NULL )
	{
		printf("Waiting for access to the circular buffer...\n");
		sleep(4);
	}
	printf("\nConnection to the circular buffer has been established.\n\n");
#ifndef	OLD_SCAN_METHOD
	scan_number = get_vol_cnt();
#endif	OLD_SCAN_METHOD
#else	RTS
	if (c_tpopen(luin,s,0) < 0)
	{
		perror( s );
		exit(1);
	}

#endif	RTS
	return;

}

read_radial()
{
	int rc;
	static int old_seq = -1;
	int diff;

reread:	;

/*
 *	Do we need to read more data?
 */

	if (left == 0)
	{
		rc = read_record();

		if ( rc == -1 )				/* i/o error */
		{
			left = 0;
			goto reread;
		}
		else if ( rc == 0 )			/* eof */
		{
			left = 0;
			goto reread;
		}
#ifndef	RTS
		else if ( rc == 8 )			/* tape header */
		{
			adjust_a2_pointer( "tape header", 8 );
			goto reread;
		}
		else if ( rc == 24 )			/* file header */
		{
			buf[12] = NULL;
			scan_number = atoi( &buf[9] );
			adjust_a2_pointer( "file header", 24 );
			goto reread;
		}
		else if ( rc < REC_SIZE )
#else	RTS
		else if ( rc < a2size )
#endif	RTS
		{
			printf("read_radial:  Funny record length of %d\n",rc);
			adjust_a2_pointer( "funny record", rc );
			goto reread;
		}

		left = rc;
		new = 0;
	}

	ptr = &buf[new];
	a2header = (struct a2header *) ptr;

#ifdef	MISPLACED
	dptr = &buf[new + sizeof(struct a2header)];
	a2data = (struct a2data *) dptr;

/*
 *	Make sure the number of gates info is right.  It *might* be wrong.
 *	Do NOT use a2data->dbz/vel/sw_ptr!
 */

	if (a2data->new_dbz_ptr == 0)
		a2data->num_gates_dbz = 0;
	if (a2data->new_vel_ptr == 0)
		a2data->num_gates_vel = 0;
#endif	MISPLACED

#ifdef	DEBUG_1
	printf("type %d ",a2header->message_type);
	printf("seq #%x\n",a2header->seq);
#endif	DEBUG_1

	if (old_seq >= 0)
	{
/*
 *	Repair an incorrect data value if necessary.
 */
		if (a2header->message_type != 1 && a2header->seq == 0)
		{
#ifdef	DEBUG_1
			printf("read_radial:  fixing seq number...\n");
#endif	DEBUG_1
			a2header->seq = old_seq + 1;
		}
		diff = a2header->seq - old_seq;
		if (old_seq > a2header->seq)
			diff += 0x7fff;
/*
 *	Don't complain too much, as we may not see all non-data messages.
 */
		if (diff > 5)
			printf("read_radial:  %d radials missing\n",diff);
	}
	old_seq = a2header->seq;
	if (old_seq == 0x7fff)
		old_seq = 0;

/*
 *	Decode data types.  Those types above 200 are local additions.
 */

	if ( a2header->message_type == A2_DATA_TYPE )
	{
		dptr = &buf[new + a2size];
		a2data = (struct a2data *) dptr;

/*
 *	Make sure the number of gates info is right.  It *might* be wrong.
 *	Do NOT use a2data->dbz/vel/sw_ptr!
 */

		if (a2data->new_dbz_ptr == 0)
			a2data->num_gates_dbz = 0;
		if (a2data->new_vel_ptr == 0)
			a2data->num_gates_vel = 0;
	}

	else if ( a2header->message_type == A2_VOLSCAN_TYPE )
	{
		scan_number = (short) &buf[new + a2size];
#ifdef	DEBUG_1
		printf("Decode scan_number is %d\n",scan_number);
#endif	DEBUG_1
#ifdef	OLD_SCAN_METHOD
		printf("OLD_SCAN_METHOD defined, still got scan_number of %d\n",
			scan_number);
#endif	OLD_SCAN_METHOD
		adjust_a2_pointer( "decode scan number", a2size+2 );
		goto reread;
	}

	else
	{
#ifdef	DEBUG_1
		printf("read_radial:  wrong type...trying again...\n");
		printf("read_radial:  not processing message_type of %d\n",
			a2header->message_type;
#endif	DEBUG_1
		adjust_a2_pointer( "skipping uneeded data", REC_SIZE );
		goto reread;
	}

#ifdef	OLD_SCAN_METHOD
/*
 *	Is this a new volume scan?  (RTS only)
 */

#ifdef	RTS
	if (old_elevation_number == 999)
	{
		scan_number = 1;
		old_elevation_number = a2data->elevation_number;
	}
	if (old_elevation_number > a2data->elevation_number)
		scan_number++;
	old_elevation_number = a2data->elevation_number;
#endif	RTS
#endif	OLD_SCAN_METHOD

/*
 *	I'm not sure why I have to make these first gate corrections...
 */

	if (get_status(DBZ))
		a2data->range_first_gate_dbz += (0.5 * a2data->gate_size_dbz);
	if (get_status(VEL))
		a2data->range_first_gate_vel += (0.5 * a2data->gate_size_vel);

/*
 *	Find VCP info.
 */

	switch (a2data->vcp)
	{
		case 11:
			vcp = vcpat11;
			break;
		case 21:
			vcp = vcpat21;
			break;
		case 31:
			vcp = vcpat31;
			break;
		case 32:
			vcp = vcpat32;
			break;
		case 300:
			vcp = vcpat300;
			break;
		default:
			if (warn_flag == 0)
				printf("read_radial:  %d:  unknown vcp number\n",
					a2data->vcp);
			vcp = NULL;
			warn_flag++;
			if (warn_flag > WARN_FREQ)
				warn_flag = 0;
			break;	
	}

/*
 *	Convert milli-seconds to seconds.
 */
 
	a2data->zulu_time /= 1000;

/*
 *	Add the number of days.
 */

	a2data->zulu_time += a2data->mod_julian_date * SECONDS_PER_DAY;

/*
 *	For UNIX, time=0 is December 31, 1969, while it is January 1, 1970
 *	for NEXRAD systems, so we need to subtract off one day.
 */

	a2data->zulu_time -= SECONDS_PER_DAY;
	gmt = gmtime(&a2data->zulu_time);
	bcopy(gmt,gmt_data,sizeof(struct tm));
	gmt = (struct tm *) gmt_data;

/*
 *	Check to see if we have a valid scan number.  If not, we try again.
 */

	if ( scan_number == 0 )
	{
		if (scan_flag == 0)
			printf("read_radial:  looking for valid scan number\n");
		scan_flag++;
		if (scan_flag > WARN_FREQ)
			scan_flag = 0;
		adjust_a2_pointer( "done decoding data", REC_SIZE );
		goto reread;
	}

	adjust_a2_pointer( "done decoding data", REC_SIZE );

	return;
}

/*
 *	Read a data record.
 */

read_record()
{
	int i, rc;
	static int old_rc = -1;
	static int data_flag = 1;

#ifndef	RTS
	rc = c_tpread(luin, buf, 145920);
	i = c_tpstatus(luin);
	if (i != 0)
		failure++;
	if (failure >= MAX_FAIL)
	{
		printf("read_record:  Too many read failures.  Exit.\n");
		exit(1);
	}
	if (old_rc == 0 && rc == 0)
	{
		printf("read_record:  double eof\n");
		exit(0);
	}
	if (buf[0] == 'U' && buf[1] == 'F')
	{
		printf("read_record:  Can't process Universal Format data.  Exit.\n");
		exit(1);
	}

	old_rc = rc;
#else	RTS
reread:	;
	rc = wideband_read( my_cbuff, buf, 145920, &pointer );
	if (rc < 0)
	{
		perror( "read_record:  wideband_read" );
		wideband_done_reader();
		exit(1);
	}
	if (rc == 0)					/* no data */
	{
		if ( data_flag == 0 )
			printf("Waiting for data from the circular buffer\n");
		data_flag++;
		if ( data_flag > WARN_FREQ )
			data_flag = 0;
		usleep( 100000 );
		goto reread;
	}
	else
		data_flag = 1;
#endif	RTS

#ifdef	DEBUG_0
	if (rc == 0)
		printf("read_record:  eof\n");
	else if (rc == -1)
		printf("read_record:  i/o error?\n");
	else
	{
		printf("read_record:  %d byte record",rc);
		if ( rc < a2size )
			printf(" --- bad record length");
		printf("\n");
	}
#endif	DEBUG_0

	left = rc;
	new = 0;
	return( rc );
}

/*
 *	Function data_wait() only gets called when we are broadcasting data.
 *	It really should be in "make_frame.c" or "ll_radial.c", but was put
 *	here as this file is still being edited quite a bit.
 *
 *	We use this to simulate pseudo-realtime for playback.
 */

/*
 *	SS10 needs 10000.
 *	IPC needs 25000.
 *	Use 70000 to simulate realtime.
 */

data_wait()
{
#ifndef	RTS
/*
	usleep(10000);
	usleep(25000);
	usleep(70000);
 */
	usleep(50000);
#endif	RTS
}

/*
 *      The altitude info (station elevation) may now be in the data file.
 */

get_altitude()
{
        return( alt );
}

get_azi()
{
	return( (int) ( (float) (a2data->azimuth) * A_SCALE  * 100) );
}

get_beamwidth()
{
	return( (int) ( BEAMWIDTH * 100 ) );
}

/*
 *	Need to initialize to missing, and check to see the number of fields
 *	returned?
 */

get_data(dbz,vel,snr,sw,n)
float dbz[], vel[], snr[], sw[];
int n;
{
	int i;
	int k;
	float *zptr, *vptr, *sptr;
	unsigned char value;

/*
 *	Initialize.
 */

	for(i=0; i<n; i++)
	{
		dbz[i] = MISSING_DATA;
		vel[i] = MISSING_DATA;
		snr[i] = DEFAULT_SNR;
		sw[i]  = MISSING_DATA;
	}

/*
 *	Load reflectivity.
 */

	if (get_status(DBZ) != 0)
	{
		k = n;
		if (a2data->num_gates_dbz < n)
			k = a2data->num_gates_dbz;
		zptr = dbz;
		for(i=0; i<k; i++)
		{
			value = (unsigned) dptr[a2data->dbz_ptr+i];
			if ( value == 0 )
				*zptr = MISSING_DATA;
			else if ( value == 1 )
				*zptr = RANGE_FOLDED_DATA;
			else
				*zptr = 0.5 * (float) value - 33.0;
			zptr++;
		}
	}

/*
 *	Load velocity and spectral width.
 */

	if (get_status(VEL) != 0)
	{
		k = n;
		if (a2data->num_gates_vel < n)
			k = a2data->num_gates_vel;
		vptr = vel;
		sptr = sw;
		for(i=0; i<k; i++)
		{
			value = (unsigned) dptr[a2data->vel_ptr+i];

			if ( value == 0 )
			{
				*vptr = MISSING_DATA;
				*sptr = MISSING_DATA;
			}
			else if ( value == 1 )
			{
				*vptr = RANGE_FOLDED_DATA;
				*sptr = RANGE_FOLDED_DATA;
			}

			else
			{
				if (a2data->vel_res == 2)
					*vptr = 0.5 * (float) value - 64.5;
				else if (a2data->vel_res == 4)
					*vptr = (float) value - 129.0;
				else
				{
					*vptr = MISSING_DATA;
					printf("get_data:  %d:  bad vel_res\n",
						a2data->vel_res);
				}
				value = (unsigned) dptr[a2data->sw_ptr+i];
				*sptr = 0.5 * (float) value - 64.5;
			}
			vptr++;
			*sptr++;
		}
	}
	return;
}

get_day()
{
	return( gmt->tm_mday );
}

get_elev()
{
	return( (int) ( (float) (a2data->elevation) * A_SCALE * 100) );
}

get_first_gate(n)
int n;
{
	if (get_status(n) == 0)
		return( MISSING );
	if (n == DBZ)
		return( a2data->range_first_gate_dbz );
	else if (n == VEL || n == SW)
		return( a2data->range_first_gate_vel );
	else
		return( 0 );
}

get_fixed_angle()
{
	int magic;
	int afix;

/*
 *	If we don't have a valid VCP number, then just send back the current
 *	elevation angle.
 */

	if ( vcp == NULL )
		return( get_elev() );

	magic = (a2data->elevation_number - 1) * 17 + 6;
	afix = (float) vcp[magic] * A_SCALE * 100;

	return ( (int) afix );
}

get_frequency()
{
	return( MISSING );
}

get_gate_spacing(n)
{
	if (get_status(n) == 0)
		return( MISSING );
	if (n == DBZ)
		return( a2data->gate_size_dbz );
	else if (n == VEL || n == SW)
		return( a2data->gate_size_vel );
	else
		return( 0 );
}

get_hour()
{
	return( gmt->tm_hour );
}

get_latitude()
{
	return( lat );
}

get_longitude()
{
	return( lon );
}

get_min()
{
	return( gmt->tm_min );
}

get_month()
{
	return( gmt->tm_mon + 1 );
}

get_nyquist()
{
	if (get_status(VEL) == 0)
		return( MISSING );
	return( (int) a2data->nyquist );
}

get_number_of_gates(n)
int n;
{
	if (n == DBZ)
		return( a2data->num_gates_dbz );
	else if (n == VEL || n == SW)
		return( a2data->num_gates_vel );
	else
		return( 0 );
}

/*
 *	Check this out later...
 */

get_polarization()
{
	return( 1 );
}

/*
 *	I'm not sure what this should be...
 */

get_power()
{
	return( MISSING );
}

get_prf()
{
	return( 1000 );
/*
	return ( MISSING );
 */
}

/*
 *	Fill it later...
 */

get_prt()
{
	if (get_status(VEL) == 0)
		return( MISSING );
/*
	return( (int) (uf_doparm[VEL][10] * 1000) );
 */
	return( MISSING );
}
	
get_proj_name(s)
char *s;
{
	bzero(s,16);
	strncpy(s,"REALTIME NEXRAD",15);
	return;
}

/*
 *	I'm not sure what this is...
 */

get_pulse_width()
{
	return( MISSING );
}

get_radar_name(s)
char *s;
{
	bzero(s,16);
	strncpy(s,radar_name,4);
	return;
}

get_scan()
{
	return( scan_number );
}

/*
 *	NEEDS CHECKED OUT!
 */

get_scan_dir()
{
	int magic;
	int afix;

/*
 *	If we have an invalid VCP number, just punt.
 */
	if ( vcp == NULL )
		return( 0 );

	magic = (a2data->elevation_number - 1) * 17 + 11;
	afix = (float) vcp[magic] * A_SCALE;

	if (afix > 0.0)
		return( -1 );
	else if (afix < 0.0)
		return( 1 );
	else
		return( 0 );
}

get_sec()
{
	return( gmt->tm_sec );
}

get_site_name(s)
char *s;
{
	bzero(s,16);
	strncpy(s,site_name,8);
	return;
}

get_status(n)
int n;
{
	int val;

	if (n == DBZ)
		val = a2data->num_gates_dbz;
	else if (n == VEL || n == SW)
		val = a2data->num_gates_vel;
	else
		val = 0;
	
	if (val > 0)
		return( 1 );
	else
		return( 0 );
}

get_tilt()
{
	return( a2data->elevation_number );
}

get_tilt_type()
{
	return( 1 );
}

get_vcp()
{
	return( a2data->vcp );
}

get_year()
{
	return( gmt->tm_year );
}

adjust_a2_pointer( msg, offset )
char *msg;
int offset;
{
	left -= offset;
	new += offset;

	if ( left < 0 )
	{
		printf("read_radial: logic fault: %s\n",msg);
		left = 0;
	}
}

set_radar_name(s)
char *s;
{
	strncpy(radar_name,s,4);
	radar_name_set++;
}

#define	SIZE	80

FILE *fp;

get_radar_info()
{
	char line[SIZE];
	char *default_name = "radarinfo.dat";
	char name[4], site[8];
	char *ptr;
	int lat1,lat2,lat3,lon1,lon2,lon3;

/*
 *	Initialize info to missing.  Do NOT initialize "radar_name" here as
 *	the user may have already called set_radar_name().
 */

	strncpy( site_name, "Unknown ", 8 );
	lat = 0;
	lon = 0;

/*
 *	Retrieve the radar name.
 *
 *	First, look and see if "set_radar_name()" has been called.  If it
 *	hasn't, we look for an environmental variable.
 */

	if ( radar_name_set == 0 )
	{
		ptr = getenv("RADARNAME");
		if (ptr == NULL)
		{
			printf("WARNING:  no call to set_radar_name, and no environmental variable RADARNAME\n");
			strncpy(radar_name,"UNKN",4);
			goto finish;
		}
		else
			strncpy(radar_name,ptr,4);
	}

/*
 *	Open the radar info file.
 *
 *	Look for an environment variable.
 *
 *	If the environmental variable isn't set,  look for the radar file in
 *	the current directory.  If it is set, and the file doesn't exist,
 *	complain and still look for the file in the current directory.
 */

	ptr = getenv("RADARFILE");
	if ( ptr != NULL )
	{
		fp = fopen(ptr,"r");
		if (fp == NULL)
			perror(ptr);
	}

	if (fp == NULL)
	{
		fp = fopen(default_name,"r");
		if (fp == NULL)
		{
			perror(default_name);
			printf("WARNING:  can't find a valid '%s' file\n",
				default_name);
			goto finish;
		}
	}

	fgets(line,SIZE,fp);			/* eat header line */

	while( fgets(line,SIZE,fp) != NULL )
	{
		sscanf(line,"%4s %8c %d %d %d %d %d %d",
			name,site,&lat1,&lat2,&lat3,&lon1,&lon2,&lon3);
		if (strncmp( radar_name, name, 4 ) == 0)
		{
			strncpy(site_name, site, 8);
			break;
		}
	}

	fclose(fp);

	lat = dms_df( lat1, lat2, lat3 );
	lon = dms_df( lon1, lon2, lon3 );

finish:	;

	printf("Radar:	%-4.4s\nSite:	%-8.8s\nLat.:	%f\nLong.:	%f\n\n",
		radar_name, site_name, lat/100000., lon/100000.);

	return;
}

dms_df(a,b,c)
int a,b,c;
{
	float tmp;
	int ans;

	tmp = ( (float) a + (float) b / 60.0 + (float) c / 3600.0 );

	ans = (int) ( tmp * 100000 );

	return( ans );
}
