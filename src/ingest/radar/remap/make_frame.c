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
 *	make_frame...
 *
 *	Build and broadcast frames of data.
 */

#ifndef	LINT
static char make_frame_id[] = "@(#)make_frame.c	1.4	3/9/94";
#endif	LINT

#include <stdio.h>
#include "rt_packet.h"
#include "rt_config.h"

char packet[1520];
int port;

int rec_size 	= sizeof(struct record_id);
int frame_size 	= sizeof(struct frame_header);
int rad_size 	= sizeof(struct radial_header);
int datah_size 	= sizeof(struct data_header);
int data_size 	= sizeof(struct data_radial);
int crc_size 	= sizeof(struct crc);

int pcount;

/*
 *	Create all the frames for the current radial, and broadcast them.
 */

load_frames()
{
	int i, nframes;

	nframes = compute_frames();

	for(i=0; i<nframes; i++)
	{
		make_frame_header(i+1,nframes);
		make_frame(i);

		if (i == 0)
			make_first_frame();
		else
			make_second_frame();
		broadcast();
	}

/*
 *	Pause between frames.
 */
	data_wait();

	return;
}

/*
 *	Figure out how many frames we'll need to broadcast the data.
 */

compute_frames()
{
	int n, nframes;

	n = radial_header.number_of_gates;;
	nframes = 0;

	while(n > 0)
	{
		nframes++;
		if (nframes == 1)
			n -= FRAME_TYPE1_NUM_DATA;
		else
			n -= FRAME_TYPE2_NUM_DATA;
	}
	return(nframes);
}

/*
 *	Prepare the frame header info.
 */

make_frame_header(n1,n2)
int n1,n2;
{
	if (n1 == 1)
		frame_header.frame_type = FRAME_TYPE1;
	else
		frame_header.frame_type = FRAME_TYPE2;

	frame_header.frame_seq_number++;
	frame_header.num_frames_per_radial = n2;
	frame_header.frame_number = n1;
}

/*
 *	Prepare the frame data.
 */

make_frame(n)
int n;
{
	int i;
	int begin, kount;

	if (n == 0)
	{
		begin = 1;
		kount = FRAME_TYPE1_NUM_DATA;
	}
	else
	{
		begin = FRAME_TYPE1_NUM_DATA + 
			(n - 1) * FRAME_TYPE2_NUM_DATA + 1;
		kount = FRAME_TYPE2_NUM_DATA;
	}

	if (begin + kount - 1 > radial_header.number_of_gates)
		kount = radial_header.number_of_gates - begin + 1;
	data_header.starting_gate = begin;
	data_header.number_of_gates = kount;

	if (begin == 1)
		return;

	for(i=begin-1; i<begin+kount; i++)
	{
		data_radial[i-begin].prod_set[0] = data_radial[i].prod_set[0];
		data_radial[i-begin].prod_set[1] = data_radial[i].prod_set[1];
		data_radial[i-begin].prod_set[2] = data_radial[i].prod_set[2];
		data_radial[i-begin].prod_set[3] = data_radial[i].prod_set[3];
	}
}

/*
 *	Load the first frame.
 */

make_first_frame()
{
	pcount = rec_size;

	load_frame_header();
	load_radial_header();
	load_data_header();

}

/*
 *	Load the second (and later) frame.
 */

make_second_frame()
{
	pcount = rec_size;

	load_frame_header();
	load_data_header();

}

/*
 *	Load the record id info.
 */

load_record_id()
{
	bcopy((char *)&record_id,&packet[pcount],rec_size);
	pcount += rec_size;
}

/*
 *	Load the frame header.
 */

load_frame_header()
{
	bcopy((char *)&frame_header,&packet[pcount],frame_size);
	pcount += frame_size;
}

/*
 *	Load the radial header.
 */

load_radial_header()
{
	bcopy((char *)&radial_header,&packet[pcount],rad_size);
	pcount += rad_size;
}

/*
 *	Load the data.
 */

load_data_header()
{
	int add_size;

	add_size = data_size * data_header.number_of_gates;

	bcopy((char *)&data_header,&packet[pcount],datah_size);
	pcount += datah_size;

	bcopy((char *)data_radial,&packet[pcount],add_size);
	pcount += add_size;

	bcopy((char *)&crc,&packet[pcount],crc_size);
	pcount += crc_size;
}

/*
 *	Initialize the CFT output data stream.  This requires the following
 *	line in the file pointed by the environmental variable ALG_SERVICES:
 *
 *	cftdata		everyone	NIT	4000	le1
 *
 *	Change "le1" to the appropriate ethernet interface.  The following
 *	line must be in /etc/ethers (or in NIS master file):
 *
 *	ff:ff:ff:ff:ff:ff	everyone
 *
 */

start_cft()
{
	port = sc_server_init("cftdata");
}

/*
 *	Load the record id info, and broadcast the packet.
 */

broadcast()
{
	int struct_size;

	struct_size = pcount;
	pcount = 0;

	make_record_id(struct_size);
	load_record_id();

	sc_server_send_data(port,packet,struct_size);
}
