/******************************************************************************************************************
 *						riscprg.h							  *
 *														  *
 *  Version  Description										  Date    *
 *  -------  -----------------------------------------------------------------------------------------  --------  *
 *    1.00   Preliminary program development								01.12.93  *
 *    1.01   Mods for new circular buffers								02.08.94  *
 *    2.00   Modification for release 3 
 *														  *
 ******************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <sys/param.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include "c_buff.h"

/**** #include "x.25.h"  ****/
#include "mess.h"

#ifndef _RIDDS_H_

#define _RIDDS_H_

#define MESS_ELEMENT		     5
#define MESS_WIDTH		    80

#define STARTUP_RESTART		0x0001					/* restart system from workstation fail	  */
#define STARTUP_WIDEBAND	0x0002					/* start up wideband services		  */
#define STARTUP_STATUS		0x0004					/* start up status services		  */
#define STARTUP_CLIENT		0x0800					/* start up client instead of servers	  */

#define RADAR_MESS		0x0001					/* radar messages include DRD and RDA 	  */
#define LOCAL_MESS		0x0002					/* local messages sent to the que	  */
#define ADMIN_MESS		0x0004					/* local control messages sent to que	  */

struct configuration {							/* strucure used for config information	  */
	char	host[64];						/* name of ingest machine		  */
	char	current_host[64];					/* name of current machine		  */
	char	default_dir[256];

	char	radar_name[64];						/* radar name (generally WFO office site  */
	char	radar_call_sign[16];					/* radar call sign			  */
	double	radar_lat;						/* radar position in lat		  */
	double	radar_lon;						/* radar position in lon		  */
	double	radar_elev;						/* radar elevation in meters		  */

	int	shm_key_comm;						/* shared memory key value for x25 comms  */
	int	shm_key_widebd;						/* shared memory key for raw data	  */
	int	shm_key_edited;						/* shared memory key for edited data	  */
	int	shm_key_tape_status;					/* shared memory key for tape status	  */
	int	shm_key_global_status;					/* shared memory key for global status	  */
	int	shm_key_local_status;					/* shared memory key for local status	  */
	int	shm_key_rpg_comm;					/* shared memory key for rpg control	  */
	int	shm_key_tape_comm;					/* shared memory key for tape control	  */

	int	shm_size_comm;						/* size of shared memory for x25 comms	  */
	int	shm_size_widebd;					/* size of wideband data distr buffer	  */
	int	shm_size_edited;					/* size of edited data distr buffer	  */

	int	supervisor_mode;					/* run as the root or superuser		  */

	char	hsi_init_command[256];					/* command to init HSI device		  */
	char	hsi_iflayer_command[256];				/* command to attach HSI device to layer  */

	char 	config_fl[256];						/* path name to configuration file	  */
	char	error_log_fl[256];					/* path name to error log file		  */
	char	watchdog_device[256];					/* path name for watchdog device	  */

	char	widebd_broadcast_addr[32];				/* broadcast address for the widebd data  */
	int	widebd_broadcast_port;					/* port address for the widebd data	  */
	char	edited_broadcast_addr[32];				/* broadcast address for edited data	  */
	int	edited_broadcast_port;					/* port address for edited data		  */
	char	status_broadcast_addr[32];				/* broadcast address for status data	  */
	int	status_broadcast_port;					/* port address for status data		  */

	int	volume_reset_flag;					/* flag value to allow automatic volume # */
	time_t	volume_reset_time;					/* time of day to reset volume number	  */

	int	tape_sync_flag;						/* sync the tape change with volume reset */
	int	tape_stragety;						/* tape archive stragety		  */
	int	tape_device_count;					/* number of tape devices		  */
	char	tape_devices[4][256];					/* tape devices				  */
	int	tape_size[4];						/* tape size in megabytes		  */
	char	tape_log_fl[256];					/* path name to tape index log file	  */

	char	restart_vol_id[256];					/* filename for volume restart		  */
	char	restart_vol_alt[256];					/* filename for alternate volume restart  */
	char	reset_prog[256];					/* filename to run when restart occurs	  */
	char	startup_prog[256];					/* filename to run at statup		  */
	char	restart_prog[256];					/* filename to run when volume is reset	  */

	time_t	disconnect_wait;					/* seconds to wait before reconnect	  */
	time_t	disconnect_limit;					/* seconds to allow the retries in	  */
	int	disconnect_retries;					/* retries within disconnect limit	  */

	time_t	heartbeat;						/* heartbeat time in seconds		  */

	int	short_vol_count;
	};

struct system_status {
	time_t	system_time;						/* Current System Time (UTC)		  */
	time_t	radar_time;						/* Current Radar Time, from DRD mess	  */
	time_t	volume_time;						/* Start time of the current volume	  */
	int	volume;							/* current volume number		  */
	int	radar_vcp;						/* current vcp setting			  */
	int	azimuth;						/* current antenna azimuth		  */

	int	radial;

	int	elevation;						/* current antenna elevation		  */
	int	tilt;							/* current tilt number			  */
	int	connection_state;					/* state of connection ot wsr88d host	  */
	int	watts;							/* power output				  */

	short	rda_status;						/* RDA current status (from RDA Packet)	  */
	short	ops_status;						/* Operability status			  */
	short	cntrl_status;						/* control status			  */
	short	generator_state;					/* power generator state		  */
	short	rda_data_streams;					/* RDA data output enabled		  */
	short	isu_state;						/* Interference Suppresion State	  */
	short	id_rate;						/* Interference Detection Rate		  */
	short	aii_state;						/* Archive II state			  */
	short	aii_space;						/* remaining archive II space in scans	  */
	short	rda_alarm_sum;						/* RDA alarm summary			  */

	int	x25_sent;						/* number of x25 packets sent		  */
	int	x25_seen;						/* number of x 25 packets seen		  */
	int	messages_sent;						/* number of messages sent		  */
	int	messages_seen;						/* number of messages seen		  */

	int	la_state;						/* local archive state			  */
	int	la_space;						/* local archive space in volume scans	  */
	
	int	system_status;						/* status summary of system		  */
	
	int	console_mess_cnt;					/* id number for current console message  */
	int	error_mess_cnt;						/* id number for current error number	  */
	int	tape_mess_cnt;
	char	console_mess[MESS_ELEMENT*MESS_WIDTH];			/* console_message			  */
	char	error_mess[MESS_ELEMENT*MESS_WIDTH];			/* error_message			  */
	char	tape_mess[MESS_ELEMENT*MESS_WIDTH];
	};

struct tape_status {
	time_t	tape_time;						/* tape archive time			  */
	time_t	volume_time;
	int	la_space;						/* tape archive	space available		  */
	time_t	la_time;
	};

struct mess_stat {							/* message status used by the message handler */
	int	mess_type;						/* message type				  */
	int	mess_flag;						/* message status or special use flag	  */
	time_t	mess_time;						/* message time				  */
	};

#define SYSTEM_MONITOR_SHMKEY	0x88d1				/* key for system monitor shared memory segment	  */
#define LOCAL_MONITOR_SHMKEY	0x88d3				/* key for local machine shared memory segment	  */
#define DATA_BRDCST_SHMKEY	0x88d2				/* Key for wideband message shared memory segment */

#define RAW_RADAR_BUFF_SIZE	512000				/* c_buffer size for raw radar data stream	  */
#define EDIT_RADAR_BUFF_SIZE	512000				/* c_buffer size for editted radar data stream	  */

#define BUFFER_SIZE		2560

#define EL_SYSTEM		0x0001				/* system error					  */
#define EL_X25_LOW		0x0002				/* x25 low level error (packet transport layer)	  */
#define EL_X25_HIGH		0x0004				/* X25 hi level error (packet processing layer)	  */
#define EL_X25_STATE		0x0008				/* x25 state and status changes			  */
#define EL_X25_PROGRESS		0x0010				/* x25 high level state changes			  */
#define EL_RADAR_STATUS		0x0020				/* Radar Status					  */
#define EL_RADAR_ALARMS		0x0040				/* Radar Alarms					  */
#define EL_VOLUME_STATUS	0x0080				/* Volume change and state status		  */
#define EL_RDA_MESS		0x0100				/* RDA/RPG Console Messages			  */
#define EL_SUN_CONSOLE		0x0200				/* Program gernerated console messages		  */
#define EL_DIST_STATE		0x0400				/* Data Disturbution Status and State Messages	  */
#define EL_TAPE_STATE		0x0800				/* Tape state and status messages		  */
#define EL_SYS_ADMIN_MESSAGES	0x1000				/* Misc system adminstration messages		  */
#define ALL_ERRORS		0x1fff				/* all error levels				  */

#define ER_FILE			0x0001				/* Report Error to error log file		  */
#define ER_TERMINAL		0x0002				/* Report Error to program terminal		  */
#define ER_CONSOLE		0x0004				/* Report Error to system console		  */
#define ER_DIST			0x0008				/* Report Error to data distritbuter		  */
#define ALL_STREAMS		0x000f				/* Report Errors to all data streams		  */

#define RUNNING_STATE		0				/* NexRAD programs runnning ok			  */

#define TAPE_PLAYBACK		0x010				/* Playback is the current data source		  */
char *selfield(char*,int,int);
#endif
