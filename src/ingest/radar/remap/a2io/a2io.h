/*
 *	A2io.h...
 *
 *	A2 parameters...
 */

#ifndef	LINT
static char a2io_id[] = "@(#)a2io.h	2.3	12/8/95";
#endif	/* LINT */

/*
 *	Data structures
 */

struct a2header {
#ifndef	RTS
	char ctm[12];					/* not used */
#endif	/* RTS */
	unsigned short message_size;			/* not used */
	unsigned char rda_channel;			/* not used */
	unsigned char message_type;			/* was short */
	short seq;
	unsigned short julian_date;			/* not used */
	unsigned long milsec;				/* not used */
	unsigned short num_mes_seg;			/* not used */
	unsigned short cur_mes_seg;			/* not used */
} ;

struct a2data {
	/*unsigned*/ long zulu_time;
	/*unsigned*/ short mod_julian_date;
	unsigned short unamb_range;
	unsigned short azimuth;
	unsigned short azimuth_number;
	unsigned short radial_status;			/* might be wrong! */
	unsigned short elevation;
	unsigned short elevation_number;
	short range_first_gate_dbz;
	short range_first_gate_vel;
	unsigned short gate_size_dbz;
	unsigned short gate_size_vel;
	unsigned short num_gates_dbz;
	unsigned short num_gates_vel;
	unsigned short cut_sector_number;
	long calibration_constant;
	unsigned short dbz_ptr;				/* do not use! */
	unsigned short vel_ptr;				/* do not use! */
	unsigned short spw_ptr;				/* do not use! */
	unsigned short vel_res;
	unsigned short vcp;
	unsigned short dummy[4];
	unsigned short arc_dbz_ptr;			/* these are correct */
	unsigned short arc_vel_ptr;			/* these are correct */
	unsigned short arc_spw_ptr;			/* these are correct */
	unsigned short nyquist;
	short atmos;
	unsigned short tover;
} ;

struct a2scan {
	short vsn;
	short pad;
};

struct a2site {
	unsigned char st_name[64];
	unsigned char st_site[32];
	unsigned long reset_tm;
};

/*
 *	Max buffer sizes.  Put here so other programs can include this info.
 */

#ifndef	RTS
#define	MAX_RECORD_SIZE	145920
#else
#define	MAX_RECORD_SIZE	8192
#endif	/* RTS */

/*
 *	Data types.
 */

#define	A2_DATA_TYPE		1
#define	A2_VOLSCAN_TYPE		201
#define	A2_SITE_TYPE		202

/*
 *	Circular buffer types.
 */

#define	CIRC_RAW		0x88d2
#define	CIRC_CREMS_DEALIAS	0x88e0

/*
 *	Default SNR
 */

#define	DEFAULT_SNR		40.0
