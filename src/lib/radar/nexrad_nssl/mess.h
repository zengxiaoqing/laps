/**************************************************************************************************
 *						messages					  *
 *												  *
 *  Version  Description								  Date	  *
 *  -------  -------------------------------------------------------------------------  --------  *
 *  1.00.00  Preliminary program development						10.13.92  *
 *  1.01.00  Additions for ConCurrent interface						01.12.93  *
 *  1.02.00  Modifications for IPC Handler						02.02.93  *
 *  1.02.01  Use only one byte for message type (pending NEXRAD release)		02.25.94  *
 *    1.03   Additions for NSSL message types						05.03.94  *
 *												  *
 **************************************************************************************************/

/**************************************************************************************************
						Definitions
 **************************************************************************************************/

#define	MESS_DRD		0x0001			/* Digital Radar Data			  */
#define	MESS_RDA		0x0002			/* RDA Status Data			  */
#define MESS_CONSOLE_RDA_RPG	0x0004			/* Console Message, RDA to RPG		  */
#define MESS_CONSOLE_RPG_RDA	0x000a			/* Console Message, RPG to RDA		  */
#define MESS_LOOP_RPG_RDA	0x000b			/* Loop Back Test RDA to RPG		  */
#define MESS_LOOP_RDA_RPG	0x000c			/* Loop Back Test RPG to RDA		  */

#define MESS_CTU		0x0014			/* Console to user message		  */
#define MESS_UTC		0x0015			/* User to console message		  */
#define MESS_LBR		0x0016			/* Loop Back Messages from RDA		  */
#define MESS_LBU		0x0017			/* Loop Back Messages from User		  */

#define NSSL_VOL_ID		0x00c9			/* NSSL Message, Volume ID		  */

#define WSR88D_DRD_FINE		0.5			/* Fine doppler Vel resolution		  */
#define WSR88D_DRD_COARSE	1.0			/* Coarse doppler Vel resolution	  */

#define NO_ALARM		0x0000			/* No alarms				  */


/**************************************************************************
				structures
 **************************************************************************/

struct wsr88d_mess {			/* struct for message header	  */
	unsigned short	hd_size;	/* message segment size 	  */
	unsigned short	hd_type;	/* message type			  */
	unsigned short	hd_sequence;	/* ID sequence (Mod 0x7fff)	  */
	unsigned short  hd_date;	/* Julian Date from Jan 1, 1970	  */
	unsigned short	hd_time[2];	/* milliseconds from midnight GMT */
	unsigned short	hd_count;	/* number of message segments	  */
	unsigned short	hd_segment;	/* segment number of this message */
	};

struct wsr88d_x25mess {			/* struct for X.25 message hdr	  */
	unsigned short	undocumented;	/* undocumented halfword	  */
	struct wsr88d_mess header;	/* message header		  */
	};

struct wsr88d_drd {			/* struct for digital radar data  */
	unsigned long	drd_time;	/* data collection time millsec   */
	unsigned short	drd_date;	/* Julian Date from Jan 1,1970    */
	unsigned short	drd_range;	/* Unambiguous range interval size*/
	unsigned short	drd_az_angle;	/* azimuth angle (horz rast scan) */
	unsigned short  drd_az_number;	/* azimuth number within cut	  */
	unsigned short	drd_rd_status;	/* radial status		  */
	unsigned short	drd_el_angle;	/* elevation angle		  */
	unsigned short	drd_el_number;	/* elevation number		  */
	unsigned short	drd_sv_range;	/* surveillance range		  */
	unsigned short	drd_dp_range;	/* doppler range		  */
	unsigned short	drd_sv_interval;/* surveillance gate size	  */
	unsigned short	drd_dp_interval;/* doppler gate size		  */
	unsigned short	drd_sv_count;	/* number of surveillance gates	  */
	unsigned short	drd_dp_count;	/* number of doppler gates	  */
	unsigned short	drd_scan_prof;	/* Pointer within RDA Control	  */
	unsigned long	drd_cal_const;	/* calibration constant for Z	  */
	unsigned short	drd_sv_pointer; /* first location of survey data  */
	unsigned short	drd_ve_pointer; /* first location of velocity	  */
	unsigned short	drd_wt_pointer; /* first location of width	  */
	unsigned short	drd_ve_resol;	/* doppler velocity resolution	  */
	unsigned short  drd_vl_pattern;	/* volume coverage pattern	  */
	unsigned char   drd_spare_a[8];	/* reserved for V+V simulation	  */
	unsigned short	drd_a2_ref;	/* Archive II pointer (not used)  */
	unsigned short	drd_a2_vel;	/* Archive II pointer (not used)  */
	unsigned short	drd_a2_width;	/* Archive II pointer (not used)  */
	unsigned short	drd_nyquist;	/* nyquist velocity		  */
	unsigned short	drd_atmos;	/* Atmospheric Attenuation Fact	  */
	unsigned short	dtd_tover;	/* threshold parameter		  */
	unsigned char	drd_spare_b[33];/* reserved for future use	  */
	unsigned short  drd_rd_data[2300]; /* data segment for radial data  */
	};

struct wsr88d_rda {			/* structure for RDA status data  */
	unsigned short	rda_status;	/* RDA status code		  */
	unsigned short	rda_operation;	/* RDA operational status code	  */
	unsigned short	rda_control;	/* RDA control status code	  */
	unsigned short	rda_genstatus;	/* RDA aux power generator status */
	unsigned short	rda_xmitpow;	/* RDA average transmitter power  */
	unsigned short	rda_refcorrect;	/* RDA reflectivity calib corr	  */
	unsigned short	rda_dtenables;	/* RDA data transmission enable   */
	unsigned short	rda_vol_pattern;/* RDA Volume coverage pattern    */
	unsigned short	rda_cntrl_auth;	/* RDA Control authorization	  */
	unsigned short	rda_infer_det;	/* RDA Interference Detection Rate*/
	unsigned short	rda_oper_mode;	/* RDA Operational Mode		  */
	unsigned short	rda_infer_sup;	/* RDA Interference Suppresion	  */
	unsigned short	rda_archive_st;	/* RDA archive II status	  */
	unsigned short	rda_archive_lt; /* RDA archive II remaining CAP	  */
	unsigned short	rda_alarm_sum;	/* RDA alarm summary		  */
	unsigned short	rda_comm_ack;	/* RDA command acknowledgement	  */
	unsigned short	rda_spare_a[10];/* RDA spare data space		  */
	unsigned short	rda_alarms[14];	/* RDA alarms			  */
	};

struct wsr88d_console {			/* Console messages 		  */
	unsigned short	con_size;	/* number of bytes in message	  */
	unsigned char	con_mess[402];	/* message buffer		  */
	};

struct wsr88d_loopback {		/* Loopback message		  */
	unsigned short	loop_size;	/* number of bytes in message	  */
	unsigned short	loop_mess[1199];/* message buffer		  */
	};

struct nssl_vol_id {			/* structure for NSSL Volume id	  */
	unsigned short	vol_number;	/* volume number		  */
	};
