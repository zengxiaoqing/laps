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

	Filename:	wideband_scan.c

	Description:

	This file contains the wideband_scan program source.  The
	purpose of the program is to access NEXRAD wideband messages on 
	the wideband circular buffer, and to report scan (VCP) information
	to standard out.

	Routines:	main
			interrupt_handler

	Date		Author		Action
	----		------		------
	02-27-94	B.H. Stevens	Original version.

 **********************************************************************/

#include <stdio.h>
#include <fcntl.h>
#include <signal.h>

#include "mess.h"

#define TRUE	(1==1)
#define FALSE	(1==0)

#define DEGREES_PER_BIT	180.0/32768.0

extern int 	errno;

static int	terminate_flag;

/**********************************************************************


	Routine:	main

	Description:

	(See file description above.)

	Inputs:		argc		(not used)
			argv		(not used)

	Outputs:	

	Returns:	0		all terminations
								      */

main( argc, argv )
	int	argc;
	char	*argv;

{	int 	nbytes;
	char 	buffer[2500];	
	void 	interrupt_handler();
	int 	status;

struct c_buff *server;
struct c_buff *wideband_init_reader();

int lowwater;

	struct DRD_mess {
		struct wsr88d_mess head;
		struct wsr88d_drd  body;
		};

	struct DRD_mess *mess;
	short 	type, sequence;
	int	last_tilt = 0, tilt, last_scan = 0, az_count = 0;

	terminate_flag = FALSE;

	signal( SIGTERM, interrupt_handler ) ;
	signal( SIGINT, interrupt_handler ) ;

/*	Attach to the circular buffer shared memory */

	while( terminate_flag == FALSE && (server = wideband_init_reader(0x88d2,&lowwater)) == '\0' ) {
		fprintf( stderr, "Waiting for access to wideband data source\n" );
		sleep(3);
		}

	if ( terminate_flag == FALSE ) {
		fprintf( stderr,"Wideband data source accessed\n" );
		}

	nbytes = 0;
	
	mess = (struct DRD_mess *)buffer;

/*	Read the circular buffer forever */

	while( terminate_flag == FALSE && nbytes >= 0 ) {
		errno = 0;
		nbytes = wideband_read(server, buffer, sizeof(buffer),&lowwater );

/*		Analyze message contents */

		if ( nbytes > 0 && nbytes >= (mess->head.hd_size)*2 ) {
			type = (short)mess->head.hd_type;
			sequence = mess->head.hd_sequence;

			if ( type == MESS_DRD ) {

				tilt = mess->body.drd_el_number;

				if ( tilt != last_tilt ) {

					if ( az_count > 0 ) {
						++az_count;
						fprintf( stdout, "Azcnt=%3d\n", az_count );
						fflush( stdout );
						az_count = 0;
						}

					az_count = 1;
					last_tilt = tilt;

					if ( tilt == 1 || last_scan == 0 ) ++last_scan;

					fprintf( stdout, "Volcnt=%3d VCP=%3d Elnum=%3d ",
					         last_scan,
						 mess->body.drd_vl_pattern,
						 last_tilt );
					fflush( stdout );
					}
				else
					++az_count;
				}
			}
		else if ( nbytes < 0 ) 
			fprintf( stderr,"wideband_read status=%d\n", nbytes );
		else
			usleep( 100000 );
		}
	
	if ( az_count > 0 ) fprintf( stdout, "\n" );

/*	Release (detach from) the shared memory */

	status = wideband_done_reader(server);

	if ( status >= 0 ) 
		fprintf( stderr, "Wideband data source released\n" );		
/*	else
		fprintf( stderr,"wideband_done_reader status = %d\n", status );
*/

	exit(0);
}

/**********************************************************************


	Routine:	interrupt handler

	Description:

	This routine handles INT and TERM interrupts.  It sets a
	"terminate" flag so the normal program polling loop will clean
	up properly.

	Inputs:		sig
			code
			scp
			addr

	Outputs:	(none)

	Sets:		terminate_flag

	Returns:	(none)

								      */
		
void interrupt_handler( sig, code, scp, addr )
	int	sig, code ;
	struct  sigcontext *scp ;
	char 	*addr ;
{
	/* fprintf(stderr,"In interrupt handler\n"); */

	terminate_flag = TRUE;

	return;
	}
