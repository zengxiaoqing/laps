/*
 *	a2io...
 *
 *	Decode A2 data.
 */
 
#ifndef	LINT
static char a2io[] = "%W%	%G%";
#endif	/* LINT */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include <errno.h>
#ifdef	RTS
#include "c_buff.h"
#endif	/* RTS */
#include "a2io.h"
#include "defs.h"
#include "vcp.h"

static char buf[MAX_RECORD_SIZE];

#define	A_SCALE		a_scale
float a_scale = 18000. / 32768.;

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

struct c_buff *wideband_init_reader();
struct c_buff *my_cbuff;

int pointer;

#endif

#define	SECONDS_PER_DAY		24 * 60 * 60

/*
 *	Complain if it appears that more than MAXDIFF radials have been
 *	dropped.
 */

#define	MAXDIFF		5

static int dbl_eof = 0;

int radar_name_set = 0;
char radar_name[4];
char site_name[8];
int lat, lon;
int lat1,lat2,lat3,lon1,lon2,lon3;
int alt;

/*
 *	We need to do our own scan number.
 */

static int scan_num = 0;

static int old_scan_num = -1;
static int old_tilt_num = -1;

static int scan_status, tilt_status;

static int new = 0;
static int left = 0;

static struct tm *gmt;
char gmt_data[sizeof(struct tm)];

unsigned short int *vcp;

struct a2header *a2header;
struct a2data *a2data;
struct a2scan *a2scan;
struct a2site *a2site;

int a2size;

static int vs_timestamp = 0;

char *ptr;
char *dptr;

/*
 *	Warning message frequency.
 */

#define	WARN_FREQ	500

int warn_flag = 0;
int scan_flag = 0;

#define	MAX_ELEV	2100		/* max legal elevation angle x 100 */

unsigned short elev_save = 0;		/* same type as a2data->elevation */
unsigned short elev_thresh;		/* ditto */

/*
 *	Allow the user to control the data processed.
 */

#ifndef	RTS
static int a2_start_scan = 1;
static int a2_end_scan = 9999;
static int a2_num_scan = 9999;
#endif	/* RTS */

/*
#ifndef	RTS
#define	FIX_BAD_SCAN
#else
#undef	FIX_BAD_SCAN
#endif
 */

#ifdef	FIX_BAD_SCAN
#define	SCAN_OFFSET	1000
static int scan_fix = 0;
static int use_scan_fix = 0;
#endif

/*
 *	a2_init is really a poor name selection since generic names are being
 *	used for the rest of the code.
 */

a2_init( s )
char *s;
{
	radar_init( s );
}

radar_init( s )
char *s;
{
#ifdef	RTS
	int circ_buf_name;
#endif	/* RTS */

/*
 *	Retrieve radar info.
 */

	get_radar_info();

	a2size = sizeof( struct a2header );

/*
 *	Open the input device, and prepare to do i/o.
 */

#ifdef	RTS

/*
 *	Set line buffering for output so that we can log the information.
 */

	fflush( stdout );
	setlinebuf( stdout );

	if ( s[0] == 'R' )
		circ_buf_name = CIRC_RAW;
	else if ( s[0] == 'E' )
		circ_buf_name = CIRC_CREMS_DEALIAS;
	else
	{
		circ_buf_name = CIRC_RAW;
		printf("WARNING:  undefined circular buffer '%c', using RAW.\n",
			s[0]);
	}
	while( ( my_cbuff = wideband_init_reader( circ_buf_name, &pointer )) 
		== NULL )
	{
		printf("Waiting for access to the circular buffer...\n");
		sleep(4);
	}
	printf("\nConnection to ");
	if ( circ_buf_name == CIRC_RAW )
		printf("the raw ");
	else if ( circ_buf_name == CIRC_CREMS_DEALIAS )
		printf("the crems/dealias ");
	else
		printf("an unknown ");
	printf("circular buffer has been established.\n\n");
	scan_num = 0; /* getvolid(); */

#ifdef	DEBUG_A2IO
	printf("radar_init:  getvolid() returns %d\n",scan_num);
#endif	/* DEBUG_A2IO */

/*
 *	I have been told that a "scan_num" of "-1" can occur if nothing is
 *	running.  I'll change this to "0" so my code can go into search mode.
 */
	if ( scan_num < 0 )
		scan_num = 0;
#else	/* RTS */

/*
 *	c_tpopen() tells the user what went wrong, so we just exit without
 *	saying anything.
 */
	if (c_tpopen(luin,s,0) < 0)
		exit(1);

/*
 *	Access environmental info.
 */

	parse_user_settings();

#endif	/* RTS */

/*
 *	Make pointers point to something so we don't core dump on program
 *	abuse (access data without a success read_radial() call).  THE INFO
 *	RETURNED WILL BE GARBAGE!
 */

	ptr = &buf[0];
	a2header = ( struct a2header *) ptr;
	a2data = ( struct a2data *) ptr;

/*
 *	Make gmt point to something useful (actually useless) so the program
 *	doesn't core dump if someone tries to access date info without going
 *	through a success read_radial() call.
 */
	bzero( gmt_data, sizeof(struct tm));
	gmt = (struct tm *) gmt_data;

/*
 *	Compute maximum elevation threshold in the units the data comes in.
 */

	elev_thresh = (unsigned short) ( (float) MAX_ELEV / A_SCALE );

	return;

}

#ifndef	RTS

/*
 *	Parse user settings
 */

parse_user_settings()
{
	char *myptr;
	int flag = 0;
	int error;

	error = 0;
	myptr = getenv("A2_START_SCAN");
	if ( myptr != NULL )
	{
		a2_start_scan = atoi( myptr );
		printf("A2_START_SCAN:	%d\n",a2_start_scan);
		if ( a2_start_scan < 1 )
			error++;
		flag++;
	}
	myptr = getenv("A2_END_SCAN");
	if ( myptr != NULL )
	{
		a2_end_scan = atoi( myptr );
		printf("A2_END_SCAN:	%d\n",a2_end_scan);
		if ( a2_end_scan < 1 )
			error++;
		flag++;
	}
	myptr = getenv("A2_NUM_SCAN");
	if ( myptr != NULL )
	{
		a2_num_scan = atoi( myptr );
		printf("A2_NUM_SCAN:	%d\n",a2_num_scan);
		if ( a2_num_scan < 1 )
			error++;
		flag++;
	}

	if ( error )
	{
		printf("FATAL ERROR:  ");
		printf("One or more of the above values are less than or equal to zero.\n");
		exit( 1 );
	}

	if ( a2_end_scan < a2_start_scan )
	{
		printf("FATAL ERROR:  ");
		printf("ending scan (%d) is less than starting scan (%d).\n",
			a2_end_scan, a2_start_scan );
		exit( 1 );
	}

	if ( flag )
		printf("\n");

	return;
}
#endif

/*
 *	Make the next radial of data to the user.  The i/o itself is done in
 *	read_record().
 *
 *	The user must call other routines to retrieve the info.
 */

read_radial()
{
	int rc;
	static int old_seq = -1;
	int diff;
#ifndef	RTS
	int eot;
#endif

reread:	;

/*
 *	Do we need to read more data?  We probably need to for tape data,
 *	as it is usually blocked together.
 */

	if (left == 0)
	{
		rc = read_record( buf );

		if ( rc == -1 )				/* i/o error */
		{
			left = 0;
			goto reread;
		}
		else if ( rc == 0 )			/* eof */
		{
			left = 0;
			if ( dbl_eof )
				return( 1 );
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
			scan_num = atoi( &buf[9] );
			adjust_a2_pointer( "file header", 24 );
#ifdef	DEBUG_A2IO
	printf("scan_num %d a2_start_scan %d a2_end_scan %d a2_num_scan %d\n",
		scan_num,a2_start_scan,a2_end_scan,a2_num_scan);
#endif
#ifdef	FIX_BAD_SCAN
			use_scan_fix = 0;
			scan_fix = scan_num;
#endif
/*
 *	If we haven't reached the starting location, reset our search code.
 *	This may produce extra "looking for valid scan number" messages,
 *	though this might be fixable.
 */
			if ( scan_num < a2_start_scan )
			{
				printf("Found volume scan %d\n", scan_num );
				scan_num = 0;
			}
/*
 *	See if we need to make any decisions.
 */
			if ( scan_num != 0 )
			{
				--a2_num_scan;
				eot = 0;
				if ( a2_num_scan < 0 )
					eot++;
				if ( a2_end_scan < scan_num )
					eot++;
				if ( eot )
				{
					dbl_eof = 1;
					return( 1 );
				}
			}
			goto reread;
		}
/*
 *	Complain if Exabyte is freaking out.  For unknown reasons, some Exabyte
 *	drives get into a state where all reads return two less bytes then they
 *	should.
 */
		else if ( rc == 6 || rc == 22 )
		{
printf("Data records appear to be two bytes short.  Take the following action:\n");
printf("1.  Dismount your tape.\n");
printf("2.  Verify that it is an A2 tape.\n");
printf("3a. If it isn't, get the correct tape.\n");
printf("3b. If it is, turn off your drive for 10 seconds, then turn it back on.\n");
			exit( 1 );
		}
		else if ( rc < REC_SIZE )
#else	/* RTS */
		else if ( rc < a2size )
#endif	/* RTS */
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

#ifdef	DEBUG_A2IO
	if (a2header->message_type != A2_DATA_TYPE )
	{
		printf("read_radial:  type %d ",a2header->message_type);
		printf("read_radial:  seq #%x\n",a2header->seq);
	}
#endif	/* DEBUG_A2IO */

	if (old_seq >= 0)
	{
/*
 *	Repair an incorrect data value if necessary.
 *
 *	Code change:  set "a2header->seq" to "old_seq" instead of "old_seq+1".
 *	This is necessary as 202 messages follow 201 messages, and cause the
 *	a complaint.
 */
		if (a2header->message_type != 1 && a2header->seq == 0)
		{
#ifdef	DEBUG_A2IO
			printf("read_radial:  fixing seq number...\n");
#endif	/* DEBUG_A2IO */
			a2header->seq = old_seq;
		}
		diff = a2header->seq - old_seq;
		if (old_seq > a2header->seq)
			diff += 0x7fff;
/*
 *	Don't complain too much, as we may not see all non-data messages.
 */
		if (diff > MAXDIFF)
			printf("read_radial:  %d radials missing\n",diff);
	}
/*
 *	Check for a "bug" case.  Tape file headers have a zero for the sequence
 *	number.  This will cause "old_seq" to get set to 0.  When we read the
 *	first valid radial, we will get a bogus "radials missing" message.
 */
#ifndef	RTS
	if ( old_seq < 0 && a2header->seq == 0 )
		a2header->seq = old_seq;
#endif	/* RTS */
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
 *	Make sure the number of gates info is right.  It *might* be wrong
 *	on some of the older Archive II tapes.
 *
 *	Do NOT use a2data->dbz/vel/spw_ptr!
 */

		if (a2data->arc_dbz_ptr == 0)
			a2data->num_gates_dbz = 0;
		if (a2data->arc_vel_ptr == 0)
			a2data->num_gates_vel = 0;
	}

#ifdef	RTS
	else if ( a2header->message_type == A2_VOLSCAN_TYPE )
	{
		dptr = &buf[new + a2size];
		a2scan = (struct a2scan *) dptr;
		scan_num = a2scan->vsn;
#ifdef	DEBUG_A2IO
		printf("read_radial:  scan_num is %d\n",scan_num);
#endif	/* DEBUG_A2IO */
		adjust_a2_pointer( "decode scan number", rc );
		goto reread;
	}
	else if ( a2header->message_type == A2_SITE_TYPE )
	{
		dptr = &buf[new + a2size];
		a2site = (struct a2site *) dptr;
		vs_timestamp = a2site->reset_tm;
		adjust_a2_pointer( "decode A2_SITE_TYPE", rc );
		goto reread;
	}
#endif	/* RTS */

	else
	{
#ifdef	DEBUG_A2IO
		printf("read_radial:  not processing message_type of %d\n",
			a2header->message_type);
#endif	/* DEBUG_A2IO */

#ifndef	RTS
		adjust_a2_pointer( "skipping unused data", REC_SIZE );
#else
		adjust_a2_pointer( "skipping unused data", rc );
#endif	/* RTS */
		goto reread;
	}

/*
 *	Check to see if we have a valid scan number.  If not, we try again.
 *
 *	If we are in search mode ( a2_start_scan ), then don't
 *	complain.
 */

	if ( scan_num == 0 )
	{
#ifndef	RTS
		if ( a2_start_scan > 1 )
			scan_flag = 1;
#endif	/* RTS */
		if (scan_flag == 0)
			printf("read_radial:  looking for valid scan number\n");
		scan_flag++;
		if (scan_flag > WARN_FREQ)
			scan_flag = 0;
		adjust_a2_pointer( "done decoding data", REC_SIZE );
		goto reread;
	}

/*
 *	Set the "status_scan" and "status_tilt" flags.
 */

	if ( old_scan_num != scan_num )
	{
		scan_status = 1;
		tilt_status = 1;
	}
	else if ( old_tilt_num != a2data->elevation_number )
	{
		scan_status = 0;
		tilt_status = 1;
	}
	else
	{
		scan_status = 0;
		tilt_status = 0;
	}

/*
 *	Check for an error state, which seems to occur on Build 8 tapes.
 *	The error is that new volume scans can start without eof and without
 *	the 24-byte header record.
 */

	if ( tilt_status && a2data->elevation_number == 1 )
	{
		if ( !scan_status )
		{
			printf("read_radial:  new volume scan without new header.\n");
#ifdef	FIX_BAD_SCAN
			use_scan_fix++;
			scan_fix += SCAN_OFFSET;
#endif
		}
	}

	old_scan_num = scan_num;
	old_tilt_num = a2data->elevation_number;

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
	bcopy( (char *) gmt, gmt_data, sizeof(struct tm));
	gmt = (struct tm *) gmt_data;

/*
 *	Check for a bogus elevation angle.  This is an indication of a radar
 *	hardware problem.
 */

	if ( a2data->elevation > elev_thresh )
	{
		printf("read_radial:  bogus elevation angle of %.2f degrees.\n",
			(float)get_elev()/100.);
		a2data->elevation = elev_save;
		printf("read_radial:  setting elevation angle to %.2f degrees.\n",
			(float)get_elev()/100.);
	}

	elev_save = a2data->elevation;
	
	adjust_a2_pointer( "done decoding data", REC_SIZE );

	return( 0 );
}

/*
 *	Read a data record.  No data unblocking (tape mode only) is done
 *	here.
 *
 *	The return value is the number of bytes read, a zero indicating
 *	an eof, or a -1, indicating an error.  If this routine is called
 *	after a double eof has occurred, then it is consided a fatal error
 *	and the program exits.
 */

read_record( mybuf )
char *mybuf;
{
	int i, rc;
	static int old_rc = -1;
#ifdef	RTS
	static int data_flag = 1;
#endif	/* RTS */

#ifndef	RTS
/*
 *	Fail if we're called and have already hit a double eof.  If this occurs,
 *	this means the user program fails to check for this case.
 */
	if ( dbl_eof )
	{
		printf("read_record:  user program ignored end of data set.  Exit.\n");
		exit( 1 );
	}
/*
 *	Read data, and keep track of i/o errors.  User program frequently do not
 *	keep track of i/o errors and fail if there is a data problem.
 */
	rc = c_tpread(luin, mybuf, MAX_RECORD_SIZE);
	i = c_tpstatus(luin);
	if (i != 0)
		failure++;
	if (failure >= MAX_FAIL)
	{
		printf("read_record:  Too many read failures.  Exit.\n");
		exit(1);
	}
/*
 *	Keep track of the number of consecutive eofs.
 */
	if (old_rc == 0 && rc == 0)
	{
		dbl_eof++;
		left = 0;	/* shouldn't be necessary, but be safe */
		return( rc );
	}
/*
 *	Some of us have been know to grab the universal version of A2 tapes.
 */
	if (mybuf[0] == 'U' && mybuf[1] == 'F' && rc >= 2 )
	{
		printf("read_record:  Can't process Universal Format data.\n");
		printf("read_record:  Try linking your program with -luftp.\n");
		exit(1);
	}

	old_rc = rc;
#else	/* RTS */
reread:	;
	rc = wideband_read( my_cbuff, mybuf, MAX_RECORD_SIZE, &pointer );
	if (rc < 0)
	{
		if ( errno == EBADMSG )
		{
			printf("read_record:  data lost\n");
			pointer = sync_c_buff( my_cbuff );
			printf("read_record:  realtime data resynced\n");
			goto reread;
		}
		else
		{
			perror( "read_record:  wideband_read" );
			wideband_done_reader( my_cbuff );
			exit(1);
		}
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
#endif	/* RTS */

	left = rc;
	new = 0;
	return( rc );
}

/*
 *	The altitude info (station elevation) may now be in the data file.
 */

get_altitude()
{
	return( alt );
}

get_azi()
{
	return( (int) ( (float) (a2data->azimuth) * A_SCALE ) );
}

/*
 *	This is just a constant for 88Ds.
 */

get_beamwidth()
{
	return( 95 );
}

/*
 *	Get_data just returns the four basic fields, and is left in for
 *	compability with the previous release.  Users are encouraged to use
 *	get_data_field().
 *	one data field data a time.
 */

get_data(dbz,vel,snr,spw,n)
float dbz[], vel[], snr[], spw[];
int n;
{
	get_data_field( DBZ, dbz, n );
	get_data_field( VEL, vel, n );
	get_data_field( SNR, snr, n );
	get_data_field( SPW, spw, n );
}

/*
 *	Retrieve a decode data field.
 *
 *	Return values:
 *		0		requested field returned
 *		MISSING		requested field not present this radial,
 *				or not present at all
 *
 *	In the last two cases, all data values are set to MISSING_DATA.
 *
 *	Arguments are:
 *		num		data field number from defs.h
 *		data		real array for the data
 *		n		the number of data points
 */

get_data_field( num, data, n )
int num;
float data[];
int n;
{
	int i, k;
	int num_gates;
	int start_pointer;
	float *fptr;
	float f1, f2;
	char *vptr;
	unsigned char value;

/*
 *	Initialize data.
 */

/*
 *	Special case:  SNR.  We don't even do a get_status() which will tell
 *	us it isn't here.
 */

	f1 = MISSING_DATA;

	if ( num == SNR )
		f1 = DEFAULT_SNR;

	fptr = data;
	for( i=0; i<n; i++)
		*fptr++ = f1;

	if ( f1 == DEFAULT_SNR )
		return( 0 );

	if ( get_status( num ) == 0 )
		return( MISSING );

	if ( num == DBZ )
	{
		num_gates = a2data->num_gates_dbz;
		start_pointer = a2data->arc_dbz_ptr;
		f1 = 0.5;
		f2 = 33.0;
	}
	else if ( num == VEL || num == RVEL )
	{
		num_gates = a2data->num_gates_vel;
		start_pointer = a2data->arc_vel_ptr;
		if ( a2data->vel_res == 2 )
		{
			f1 = 0.5;
			f2 = 64.5;
		}
		else if ( a2data->vel_res == 4 )
		{
			f1 = 1.0;
			f2 = 129.0;
		}
		else
		{
printf("get_data_field:  %d:  can't decode velocity this radial, bad vel_res.\n",
				a2data->vel_res);
			return( MISSING );
		}
	}
	else if ( num == SPW )
	{
		num_gates = a2data->num_gates_vel;
		start_pointer = a2data->arc_spw_ptr;
		f1 = 0.5;
		f2 = 64.5;
	}
	else
	{
		printf("get_data_field:  %d:  'num' is corrupt.  help!\n",num);
		return( MISSING );
	}

	k = n;
	if ( num_gates < n )
		k = num_gates;
	fptr = data;
	vptr = &dptr[start_pointer];

	for(i=0; i<k; i++)
	{
		value = (unsigned) *vptr;
		if ( value == 0 )
			*fptr = MISSING_DATA;
		else if ( value == 1 )
			*fptr = RANGE_FOLDED_DATA;
		else
			*fptr = f1 * (float) value - f2;
		fptr++;
		vptr++;
	}
	return( 0 );
}

/*
 *	Retrieve a raw data field.
 *
 *	Return values:
 *		0		requested field returned
 *		MISSING		requested field not present this radial,
 *				or not present at all
 *
 *	In the last two cases, all data values are set to MISSING_DATA.
 *
 *	Arguments are:
 *		num		data field number from defs.h
 *		data		real array for the data
 *		n		the number of data points
 *		flag		possible scale info
 *
 *	THIS CALL IS DATA SPECIFIC.  IT MAY BE DIFFERENT IN OTHER RADAR
 *	FORMAT DECODERS!
 */

get_data_field_raw( num, data, n, flag )
int num;
unsigned short int data[];
int n;
int *flag;
{
	int i, k;
	int num_gates;
	int start_pointer;
	char *vptr;
	unsigned short int *iptr;
	unsigned char value;

/*
 *	Initialize data to missing.
 */

	iptr = data;
	for( i=0; i<n; i++)
		*iptr++ = 0;

	if ( get_status( num ) == 0 )
		return( MISSING );

	if ( num == DBZ )
	{
		num_gates = a2data->num_gates_dbz;
		start_pointer = a2data->arc_dbz_ptr;
		*flag = 0;
	}
	else if ( num == VEL || num == RVEL )
	{
		num_gates = a2data->num_gates_vel;
		start_pointer = a2data->arc_vel_ptr;
		*flag = a2data->vel_res;
	}
	else if ( num == SPW )
	{
		num_gates = a2data->num_gates_vel;
		start_pointer = a2data->arc_spw_ptr;
		*flag = a2data->vel_res;
	}
	else
	{
		printf("get_data_field:  %d:  'num' is corrupt.  help!\n",num);
		return( MISSING );
	}

	k = n;
	if ( num_gates < n )
		k = num_gates;
	iptr = data;
	vptr = &dptr[start_pointer];

	for(i=0; i<k; i++)
	{
		value = (unsigned) *vptr;
/*
		*iptr++ = (unsigned short) value;
 */
		*iptr++ = (unsigned short) value;
		vptr++;
	}
	return( 0 );
}

get_day()
{
	return( gmt->tm_mday );
}

/*
 *	Current elevation angle.  Use get_fixed_angle() to get the fixed
 *	elevation angle.
 *
 */

get_elev()
{
	return( (int) ( (float) (a2data->elevation) * A_SCALE ) );
}

/*
 *	Provide users with some friendly conversions.
 */

get_field_name( n, s  )
int n;
char *s;
{
	if ( n == DBZ )
		strcpy(s,"DBZ");
	else if ( n == VEL )
		strcpy(s,"VEL");
	else if ( n == SPW )
		strcpy(s,"SPW");
	else if ( n == SNR )
		strcpy(s,"SNR");
	else if ( n == RVEL )
		strcpy(s, "RVEL" );
	else
		strcpy(s,"???");
	return;
}

get_field_num( s )
char *s;
{
	if ( !strcmp(s, "DBZ") )
		return( DBZ );
	else if ( !strcmp(s, "VEL" ) )
		return( VEL );
	else if ( !strcmp(s, "SPW" ) )
		return( SPW );
	else if ( !strcmp(s, "SNR" ) )
		return( SNR );
	else if ( !strcmp(s, "RVEL" ) )
		return( RVEL );
	else
		return(-1);
}

get_first_gate(n)
int n;
{
	if (get_status(n) == 0)
		return( MISSING );
	if (n == DBZ)
		return( a2data->range_first_gate_dbz );
	else if ( n == VEL || n == SPW || n == RVEL )
		return( a2data->range_first_gate_vel );
	else
		return( 0 );
}

/*
 *	Fixed elevation angle.  Used get_elev() to get the elevation angle
 *	of the current radial.
 */

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
	afix = (float) vcp[magic] * A_SCALE ;

	return ( (int) afix );
}

/*
 *	I don't know how to retrieve this.
 */

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
	else if ( n == VEL || n == SPW || n == RVEL )
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

get_latitude_dms( d, m, s )
int *d, *m, *s;
{
	*d = lat1;
	*m = lat2;
	*s = lat3;
}

get_longitude()
{
	return( lon );
}

get_longitude_dms( d, m, s )
int *d, *m, *s;
{
	*d = lon1;
	*m = lon2;
	*s = lon3;
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
	else if (n == VEL || n == SPW || n == RVEL )
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

/*
 *	Ditto.  I seem to recall there is some equation in one of the
 *	packages that wants a number, so we'll give it one.
 */

get_prf()
{
	return( 1000 );
}

/*
 *	Fill it later...
 */

get_prt()
{
/*
	if (get_status(VEL) == 0)
		return( MISSING );
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

get_radial_status()
{
	return( a2data->radial_status );
}

/*
 *	Give the user a mechanism to determine if we are running in realtime
 *	mode or tape playback mode.
 */

get_rt_mode()
{
#ifdef	RTS
	return( 1 );
#else
	return( 0 );
#endif
}

get_scan()
{
#ifdef	FIX_BAD_SCAN
	if ( use_scan_fix )
		return( scan_fix );
#endif
	return( scan_num );
}

get_scan_status()
{
	return( scan_status );
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
	else if (n == VEL || n == SPW || n == RVEL)
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

get_tilt_status()
{
	return( tilt_status );
}

get_tilt_type()
{
	return( 1 );
}

get_timestamp()
{
	return( a2data->zulu_time );
}

/*
 *	This routine will go away.
 *	Don't believe it.  I decided to keep it.
 */

get_vcp()
{
	return( a2data->vcp );
}

get_vs_timestamp()
{
	return( vs_timestamp );
}

get_year()
{
	return( gmt->tm_year );
}

/*
 *	Keep track of where we are in the data, so that the next time that
 *	read_radial() is called, the right things happen.
 */

adjust_a2_pointer( msg, offset )
char *msg;
int offset;
{
	left -= offset;
	new += offset;

/*
 *	If "offset" is bad, make it fault.
 */
	if ( offset <= 0 )
		left = -1;

	if ( left < 0 )
	{
		printf("read_radial:  logic fault: %s\n",msg);
		left = 0;
	}
}

set_radar_name(s)
char *s;
{
	strncpy(radar_name,s,4);
	radar_name_set++;
}

/*
 *	Radar name, site name, latitude, and longitude are all important
 *	radar parameters.  Unfortunately, none of them are in the actual
 *	data themselves.
 *
 *	The user may set the radar name by:
 *	o  set_radar_name( some_name );	[in code]
 *	o  setenv RADARNAME some_name	[ignored if set_radar_name() used]
 *
 *	The program retrieves the info from the file:
 *	o  setenv RADARFILE some_file	[data file]
 *	o  from "./radarinfo.dat"	[if RADARFILE not used]
 */

#define	SIZE	120

FILE *fp;

get_radar_info()
{
	int rc;
	char line[SIZE];
	char *default_name = "radarinfo.dat";
	char name[4], site[8];
	char *ptr;

/*
 *	Initialize info to missing.  Do NOT initialize "radar_name" here as
 *	the user may have already called set_radar_name().
 */

	strncpy( site_name, "Unknown ", 8 );
	lat = 0;
	lon = 0;
	alt = MISSING;

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
			printf("WARNING:  radar name set to 'UNKN'\n");
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
			printf("WARNING:  can't find a valid '%s' file\n",
				default_name);
			goto finish;
		}
	}

	fgets(line,SIZE,fp);			/* eat header line */

	while( fgets(line,SIZE,fp) != NULL )
	{
		rc = sscanf(line,"%4s %8c %d %d %d %d %d %d %d",
			name,site,&lat1,&lat2,&lat3,&lon1,&lon2,&lon3,&alt);
		if ( rc < 8 || rc > 9 )
		{
			printf("get_radar_info:  %d arguments returned.  ");
			printf("Radarinfo file is corrupt.  Exit.\n");
			exit( 1 );
		}
		if ( rc == 8 )
			alt = MISSING;
#ifdef	DEBUG_RADARINFO
		printf("name %4s site %8.8s lat1 %d lat2 %d lat3 %d ",
			name, site, lat1, lat2, lat3 );
		printf("lon1 %d lon2 %d lon3 %d alt %d\n",
			lon1, lon2, lon3, alt );
#endif
		if (strncmp( radar_name, name, 4 ) == 0)
		{
			strncpy(site_name, site, 8);
			break;
		}
	}

	fclose(fp);

	lat = dms_df( lat1, lat2, lat3 );
	lon = dms_df( lon1, lon2, lon3 );

	if ( lat == 0 )
		printf("WARNING:  radar '%s' not found in data file\n",
			radar_name);

finish:	;

	if ( lat == 0 )
		printf("WARNING:  latitude and longitude set to 0\n");
#ifdef	DEBUG_A2IO
	printf("\nRadar:	%-4.4s\nSite:	%-8.8s\nLat.:	%f\nLong.:	%f\n\n",
		radar_name, site_name, lat/100000., lon/100000.);
#endif	/* DEBUG_A2IO */

	return;
}

/*
 *	Convert degrees, minutes, seconds to an integer (floating point
 *	value times 100000).
 */

dms_df(a,b,c)
int a,b,c;
{
	float tmp;
	int ans;

	tmp = ( (float) a + (float) b / 60.0 + (float) c / 3600.0 );

	ans = (int) ( tmp * 100000 );

	return( ans );
}
