/******************************************************************************************************************
 *						getconfig							  *
 *														  *
 *  Description												  Date	  *
 *  --------------------------------------------------------------------------------------------------  --------  *
 *  Preliminary program development, routine to process riscrpg configuration file to structure		07.01.94  *
 *														  *
 ******************************************************************************************************************/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "ridds.h"

static char SccsId[]="@(#)getconfig.c	3.3  3/28/95  NSSL, RIDDS Configuration processing utility";

struct configuration configuration;

/*
main()
{
getconfig(&configuration,"nexrad.cfg");
}
*/

/******************************************************************************************************************
						getconfig
 ******************************************************************************************************************/

int getconfig(table,filenm)
struct configuration *table;
char *filenm;
{
char **config_table;
char **initconfig();
char **openconfig();

char *getstring();

/*
 * prcess configuration table
 */

initbuff((char *)table,sizeof(struct configuration),0);

config_table=initconfig();
config_table=openconfig(filenm,config_table);

strcpy(table->config_fl,filenm);

/*
 * process string information from configuration table
 *
 * this includes 
 *	Ingest Host Machine Name
 *	Radar Name
 *	HSI initialization strings
 *	Configuration file name, error logging file and watchdog alerting device
 *	IP addresses (simply because they are already strings) for services
 */

getstring(table->host,"INGEST_HOST",config_table);		
getstring(table->radar_name,"RADAR_NAME",config_table);	
getstring(table->radar_call_sign,"RADAR_CALL_SIGN",config_table);
getstring(table->hsi_init_command,"HSI_INIT_COMMAND",config_table);
getstring(table->hsi_iflayer_command,"HSI_IFLAYER_COMMAND",config_table);
getstring(table->error_log_fl,"ERROR_LOG_FILE",config_table);
getstring(table->watchdog_device,"WATCHDOG_DEVICE",config_table);
getstring(table->widebd_broadcast_addr,"LEVEL2_BRDCST_ADD",config_table);
getstring(table->edited_broadcast_addr,"EDDATA_BRDCST_ADD",config_table);
getstring(table->status_broadcast_addr,"STATUS_BRDCST_ADD",config_table);
getstring(table->restart_vol_id,"RESTART_VOL_FILE",config_table);
getstring(table->restart_vol_alt,"RESTART_VOL_ALT",config_table);
getstring(table->reset_prog,"RESET_PROGRAM",config_table);
getstring(table->startup_prog,"STARTUP_PROGRAM",config_table);
getstring(table->restart_prog,"RESTART_PROGRAM",config_table);

/*
 * process integer information from configuration table
 *
 * this includes
 *	shared memory keys for a variety of services
 *	sizes for shared memory segments 
 *	port addresses for IP services
 *	various flags
 *	some tape parameters
 */

getint(&(table->shm_key_comm),"HSI_SHM_KEY",config_table);
getint(&(table->shm_key_widebd),"WIDEBAND_DATA_KEY",config_table);
getint(&(table->shm_key_edited),"EDITED_DATA_KEY",config_table);
getint(&(table->shm_key_tape_status),"TAPE_STATUS_KEY",config_table);
getint(&(table->shm_key_global_status),"SYSTEM_STATUS_KEY",config_table);
getint(&(table->shm_key_local_status),"LOCAL_STATUS_KEY",config_table);
getint(&(table->shm_key_rpg_comm),"RPG_COMM_KEY",config_table);
getint(&(table->shm_key_tape_comm),"TAPE_COMM_KEY",config_table);

getint(&(table->shm_size_comm),"HSI_SHM_SIZE",config_table);
getint(&(table->shm_size_widebd),"WIDEBAND_DATA_SIZE",config_table);
getint(&(table->shm_size_edited),"EDITED_DATA_SIZE",config_table);

getint(&(table->widebd_broadcast_port),"LEVEL2_BRDCST_PORT",config_table);
getint(&(table->edited_broadcast_port),"EDDATA_BRDCST_PORT",config_table);
getint(&(table->status_broadcast_port),"STATUS_BRDCST_PORT",config_table);

getint(&(table->disconnect_retries),"DISCONNECT_RETRY",config_table);

getint(&(table->short_vol_count),"SHORT_VOL_CNT",config_table);

/* process time information from configuration table
 *
 * this includes
 *	disconnect wait
 *	disconnect limit
 *	heartbeat
 */

#if 0    /* commented out for wideband_client... rcl 4/2/96 */
gettimevalue(&(table->disconnect_wait),"DISCONNECT_WAIT",config_table);
gettimevalue(&(table->disconnect_limit),"DISCONNECT_LIMIT",config_table);
gettimevalue(&(table->heartbeat),"HEARTBEAT",config_table);
#endif

/*
 * process double information from configuration table
 *
 * this includes
 *	radar altitude
 */

getdouble(&(table->radar_elev),"RADAR_ELEVATION",config_table);

/*
 * special process information from configuration table
 *
 * this include
 *	radar latitude and longitude
 *	configuration file name
 *	volume reset flag
 *	volume reset time
 *	tape sync flag
 *	tape strategy
 *	tape device count;
 *	tape devices
 *	tape log file basename
 */

/*
 * Process Volume reset flag parameters
 */
#if 0    /* commented out for wideband_client... rcl 4/2/96 */
get_volume_values(table,config_table);
#endif

/*
 * Get tape archive parameters
 */

get_tape_values(table,config_table);

/*
 * Get Radar Positional Information
 */

get_radar_position(table,config_table);

/*
 * Get Current Host Name
 */

gethostname(table->current_host,64);
return(0);
 }

/**********************************************************************************************************
					get_radar_position
 **********************************************************************************************************/

get_radar_position(table,config)
struct configuration *table;
char **config;
{
char buffer[256];
char *tmp;
char *selfield(char*,int,int);
char *getvalue();

if ((tmp=getvalue("RADAR_POSITION",config)) == '\0') {
	table->radar_lat=0.0;
	table->radar_lon=0.0;
	}

else {
	strcpy(buffer,tmp);
	trimboth(buffer);
	table->radar_lat=atof(selfield((char *) buffer,1,(int) ','));
	table->radar_lon=atof(selfield((char *) buffer,2,(int) ','));
	}
}

/**********************************************************************************************************
					get_tape_values
 **********************************************************************************************************/

get_tape_values(table,config)
struct configuration *table;
char **config;
{
static char *changetab[]={"Independent","independent","INDEPENDENT",	/* tape is independent of volume id */
	"Dependent","dependent","DEPENDENT",				/* tape is dependent on voleme id   */
	"Disable","disable","DISABLE",					/* tape is turned off		    */
	'\0'};

static int changetokens[]={2,2,2,1,1,1,0,0,0};

static char *schemetab[]={"Single","single","SINGLE",			/* archive on a single device	  */
	"Multiple","multiple","MULTIPLE",				/* archive to multiple devices	  */
	"Jukebox","jukebox","JUKEBOX",					/* archive to a jukebox		  */
	'\0'};

static int schemetokens[]={1,1,1,2,2,2,3,3,3};

char buffer[256];
char working[256];
char *tmp;
char *getvalue();
int token;
int loop;
int fields;

/*
 * determine if the tape archive is to be used and if so, how to store volume numbers
 */

if ((tmp=getvalue("TAPE_CHANGE_FLAG",config)) == '\0') {
	table->tape_sync_flag=0;
	return;
	}

strcpy(buffer,tmp);

if ((token=strtoken(buffer,changetab,changetokens)) < 1) {
	table->tape_sync_flag=0;
	return;
	}

/*
 * determine if archive scheme, single device, multiple device or jukebox
 */

if ((tmp=getvalue("TAPE_SCHEME",config)) == '\0') {
	table->tape_stragety=1;
	}
else {
	strcpy(buffer,tmp);
	if ((token=strtoken(buffer,schemetab,schemetokens)) < 1)
		table->tape_stragety=1;
	}

/*
 * determine the devices and the number of devices to use for archiving
 */

if ((tmp=getvalue("TAPE_DEVICES",config)) == '\0') {
	table->tape_sync_flag=0;
	return;
	}

strcpy(buffer,tmp);
trimboth(buffer);
table->tape_device_count=textfields(buffer,',');
if (table->tape_device_count > 4)
	table->tape_device_count=4;

for (loop=0; loop < table->tape_device_count; loop++) {
	strcpy(table->tape_devices[loop],selfield(buffer,loop+1,','));
	trimboth(table->tape_devices[loop]);
	}


if ((tmp=getvalue("TAPE_SIZES",config)) == '\0') {		
	for (loop=0;loop < table->tape_device_count; loop++) 
		table->tape_size[loop]=2048;
	}

else {
	strcpy(buffer,tmp);
	trimboth(buffer);
	fields=textfields(buffer,',');

	for (loop=0; loop < fields && loop < 4; loop++) {
		strcpy(working,selfield(buffer,loop+1,','));
		if ((table->tape_size[loop]=atoi(working)) < 1024)
			table->tape_size[loop]=2048;
		}

	for (; loop < table->tape_device_count; loop++)
		table->tape_size[loop]=2048;
	}

/*
 * determine the log file to write tape index information to
 */

if ((tmp=getvalue("TAPE_LOG_FILE",config)) == '\0') 
	strcpy(table->tape_log_fl,"archive.log");
else
	strcpy(table->tape_log_fl,tmp);
return;
}

#if 0    /* commented out for wideband_client... rcl 4/2/96 */
/**********************************************************************************************************
					get_volume_values
 **********************************************************************************************************/

get_volume_values(table,config)
struct configuration *table;
char **config;
{
static char *stringtab[]={"ENABLE","enable","Enable","ON","on","On","DISABLE","disable","Disable","OFF","off","Off",'\0'};
static int tokens[]={1,1,1,1,1,1,0,0,0,0,0,0};

char buffer[64];
char *tmp;
char *getvalue();

int token;
time_t strtotime();

if ((tmp=getvalue("VOLUME_RESET_FLAG",config)) == '\0') {
	table->volume_reset_flag=0;
	table->volume_reset_time=(time_t)0;
	return;
	}

strcpy(buffer,tmp);

if ((token=strtoken(buffer,stringtab,tokens)) != 1) {
	table->volume_reset_flag=0;
	table->volume_reset_time=(time_t)0;
	return;
	}


if ((tmp=getvalue("VOLUME_RESET_TIME",config)) == '\0') {
	table->volume_reset_flag=0;
	table->volume_reset_time=(time_t)0;
	return;
	}
	
strcpy(buffer,tmp);
table->volume_reset_time=strtotime(buffer,6);
table->volume_reset_flag=1;
return;
}

#endif


/**********************************************************************************************************
						getdouble
 **********************************************************************************************************/

getdouble(value,para,table)
double *value;
char *para;
char **table;
{
char *getvalue();
char *tmp;
char buffer[256];
double atof();

if ((tmp=getvalue(para,table)) == '\0') {
	*value=0.0;
	return(-1);
	}
else {
	strcpy(buffer,tmp);
	*value=atof(buffer);
	return(0);
	}
}

/**********************************************************************************************************
						getint
 **********************************************************************************************************/

getint(value,para,table)
time_t *value;
char *para;
char **table;
{
char *getvalue();
char *tmp;
char buffer[256];
long strtol();

if ((tmp=getvalue(para,table)) == '\0') {
	*value=0;
	return(-1);
	}
else {
	strcpy(buffer,tmp);
	*value=(int)strtol(buffer,(char **)NULL,'\0');
	return(0);
	}
}

/**********************************************************************************************************
						gettimevalue
 **********************************************************************************************************/

gettimevalue(value,para,table)
time_t *value;
char *para;
char **table;
{
char *getvalue();
char *tmp;
char buffer[256];
time_t strtotime();

if ((tmp=getvalue(para,table)) == '\0') {
	*value=0;
	return(-1);
	}
else {
	strcpy(buffer,tmp);
	*value=strtotime(buffer,6);
	return(0);
	}
}

/**********************************************************************************************************
						getstring
 **********************************************************************************************************/

char *getstring(dest,para,table)
char *dest;
char *para;
char **table;
{
char *getvalue();
char *tmp;

if ((tmp=getvalue(para,table)) == '\0')
	strcpy(dest,"");
else
	strcpy(dest,tmp);

return(tmp);
}
