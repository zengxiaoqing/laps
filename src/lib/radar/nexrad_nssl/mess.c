/**********************************************************************************************************
 *						mess.c							  *
 *													  *
 *  Description											  Date	  *
 *  ------------------------------------------------------------------------------------------  --------  *
 *  Nexrad wideband message handling routines for the user port communications process		06.15.94  *
 *													  *
 **********************************************************************************************************/

#include <config.h>
#include <stdio.h>
#include <sys/types.h>

#include "ridds.h"

static char SccsId[]="%W%  %G% NSSL RIDDS, Message handling routines";

/*
 * Global variables
 */

extern int errno;

/**********************************************************************************************************
						handle_mess
 **********************************************************************************************************/

struct mess_stat *handle_mess(mess,status)
char *mess;
struct system_status *status;
{
static struct mess_stat m_stat;
int temp;
struct wsr88d_mess *message;
message=(struct wsr88d_mess *)mess;

switch(message->hd_type) { 
		case MESS_DRD:	m_stat.mess_type=MESS_DRD;
				handle_drd(mess,&m_stat,status);
				break;

		case MESS_RDA: 	m_stat.mess_type=MESS_RDA;
				handle_rda(mess,&m_stat,status);
				break;

		case MESS_CTU:  m_stat.mess_type=MESS_CTU; 
				handle_rda(mess,&m_stat,status);
				handle_ctu(mess,&m_stat,status);
				break;

		case MESS_UTC:	m_stat.mess_type=MESS_UTC;
				handle_utc(mess,&m_stat,status);
				break;

		case MESS_LBR: 	m_stat.mess_type=MESS_LBR;
				handle_lbr(mess,&m_stat,status);
				break;

		case MESS_LBU: 	m_stat.mess_type=MESS_LBU;
				handle_lbu(mess,&m_stat,status);
				break;

		case 201: 	m_stat.mess_type= 201;
				handle_vol(mess,&m_stat,status);
				break;

		default:	undefined_mess(mess,&m_stat,status);
				show_log(71,"Undefined Radar Message Received");
				break;
		}

return(&m_stat);
}

/******************************************************************************************************************
						handle_drd
*******************************************************************************************************************/

int handle_drd(mess,m_stat,l_stat)
struct wsr88d_mess *mess;
struct mess_stat *m_stat;
struct system_status *l_stat;
{
extern struct c_buff *server_buff;

struct wsr88d_drd  *radar_data;

time_t wsrtimetounix();
time_t temp;
time_t seconds;

double wsrangle();
double angles;
int value;

radar_data=(struct wsr88d_drd *)(((char *)mess)+sizeof(struct wsr88d_mess));

memcpy((char *)&seconds,(char *)(&(radar_data->drd_time)),sizeof(unsigned long));
m_stat->mess_time=wsrtimetounix((unsigned short)radar_data->drd_date,seconds);
m_stat->mess_flag=0;

/*
 * Get various pieces of information from the DRD message for the status message
 *
 * Currently this includes:
 *	radar azimuth
 *	radar elevation
 *	radial
 *	tilt
 *	vcp
 */

l_stat->azimuth=wsrangle(radar_data->drd_az_angle);				/* get radial azimuth		  */
l_stat->elevation=wsrangle(radar_data->drd_el_angle);				/* get radial elevation		  */

l_stat->radial=(int)radar_data->drd_az_number;					/* get radial number		  */
l_stat->tilt=(int)radar_data->drd_el_number;					/* get tilt number		  */
l_stat->radar_vcp=(int)radar_data->drd_vl_pattern;				/* get vcp			  */

l_stat->radar_time=m_stat->mess_time;

if (radar_data->drd_rd_status == 3 && l_stat->tilt == 1) {			/* new volume ?			  */
	m_stat->mess_flag=1;
	l_stat->volume_time=m_stat->mess_time;
	}

/*
 * The real-time data has a couple of problems that makes it inconsistant with the documentation and 
 * with programs down stream that espect things to match the documentation.  here we try to 'fix'
 * things (seems like the convenient thing to do)
 */

if (radar_data->drd_sv_pointer == radar_data->drd_ve_pointer) {
	radar_data->drd_sv_count=0;
	radar_data->drd_a2_ref=0;
	}
else {
	radar_data->drd_a2_ref=radar_data->drd_sv_pointer;
	}

radar_data->drd_a2_vel=radar_data->drd_ve_pointer;
radar_data->drd_a2_width=radar_data->drd_wt_pointer;

return(0);
}

/******************************************************************************************************************
						get_vcp_type
 ******************************************************************************************************************/

get_vcp_type(mess)
struct wsr88d_mess *mess;
{
struct wsr88d_drd  *radar_data;
int vcp;

radar_data=(struct wsr88d_drd *)(((char *)mess)+sizeof(struct wsr88d_mess));

vcp=(int)radar_data->drd_vl_pattern;						/* get vcp			  */

return(vcp);
}

/******************************************************************************************************************
						handle_vol
 ******************************************************************************************************************/

handle_vol(mess,m_stat,l_stat)						/* handler for NSSL Set Volume Number		  */
struct wsr88d_mess *mess;
struct mess_stat *m_stat;
struct system_status *l_stat;
{
struct nssl_vol_id *id;
extern int secondary_stream;
extern struct system_status local_stat;

id=(struct nssl_vol_id *)(((char *)mess)+sizeof(struct wsr88d_mess));
m_stat->mess_flag=id->vol_number;

return(0);
}

/******************************************************************************************************************
						handle_rda
*******************************************************************************************************************/

int handle_rda(mess,m_stat,l_stat)
struct wsr88d_mess *mess;
struct mess_stat *m_stat;
struct system_status *l_stat;
{
struct wsr88d_rda *rda_status;
extern struct system_status local_stat;

m_stat->mess_flag=0;

rda_status=(struct wsr88d_rda *)(((char *)mess)+sizeof(struct wsr88d_mess));	/* get rda status		  */

l_stat->rda_status=(int)rda_status->rda_status;					/* get rda status		  */
l_stat->ops_status=(int)rda_status->rda_operation;				/* get rda operability		  */
l_stat->cntrl_status=(int)rda_status->rda_control;				/* get rda control point	  */
l_stat->generator_state=(int)rda_status->rda_genstatus;				/* get generator state		  */
l_stat->rda_data_streams=(int)rda_status->rda_dtenables;			/* get data streams enabled	  */
l_stat->isu_state=(int)rda_status->rda_infer_sup;				/* get interfernece sup state	  */
l_stat->id_rate=(int)rda_status->rda_infer_det;					/* get interference detect rate	  */
l_stat->aii_state=(int)rda_status->rda_archive_st;				/* get AII state		  */
l_stat->aii_space=(int)rda_status->rda_archive_lt;				/* get tape space remaining	  */
l_stat->rda_alarm_sum=(int)rda_status->rda_alarm_sum;				/* get alarm summary		  */
l_stat->watts=(int)rda_status->rda_xmitpow;					/* get average transmitted power  */

if ((l_stat->generator_state & 0x0008) != 0) 
	show_log(40,"Power Change Over");

return(0);
}

/******************************************************************************************************************
						handle_ctu
 ******************************************************************************************************************/

int handle_ctu(mess,m_stat,l_stat)
struct wsr88d_mess *mess;
struct mess_stat *m_stat;
struct system_status *l_stat;
{
struct wsr88d_console *console;

char *timetostr();
char tempmess[80];
char message[2560];

console=(struct wsr88d_console *)(((char *)mess)+sizeof(struct wsr88d_mess));
*(console->con_mess+console->con_size)='\0';

strcpy(message,"UCP Message: ");
strcat(message,console->con_mess);
show_log(40,message);

strcpy(tempmess,timetostr(&(m_stat->mess_time),1));
strcat(tempmess,"  ");
strncat(tempmess,console->con_mess,60);
tempmess[61]='\0';

add_mess(MESS_ELEMENT,MESS_WIDTH,l_stat->console_mess,tempmess);
l_stat->console_mess_cnt++;

return(0);
}

/******************************************************************************************************************
						handle_utc
 ******************************************************************************************************************/

int handle_utc(mess,m_stat,l_stat)
struct wsr88d_mess *mess;
struct mess_stat *m_stat;
struct system_status *l_stat;
{
struct wsr88d_console *console;
m_stat->mess_flag=0;
return(0);
}

/******************************************************************************************************************
						handle_lbr
 ******************************************************************************************************************/

int handle_lbr(mess,m_stat,l_stat)
struct wsr88d_mess *mess;
struct mess_stat *m_stat;
struct system_status *l_stat;
{
m_stat->mess_flag=0;
return(0);
}

/******************************************************************************************************************
						handle_lbu
 ******************************************************************************************************************/

int handle_lbu(mess,m_stat,l_stat)
struct wsr88d_mess *mess;
struct mess_stat *m_stat;
struct system_status *l_stat;
{
m_stat->mess_flag=0;
return(0);
}


/******************************************************************************************************************
						undefined_mess
 ******************************************************************************************************************/

int undefined_mess(mess, m_stat, l_stat)
struct wsr88d_mess *mess;
struct mess_stat *m_stat;
struct system_status *l_stat;
{
m_stat->mess_flag=0;
return(0);
}

/******************************************************************************************************************
						write_new_vol
 ******************************************************************************************************************/

write_new_vol(buf,vol)
struct c_buff *buf;
int vol;
{
extern struct system_status local_stat;
char buffer[256];

unsigned long time_tmp;
struct wsr88d_mess *mess;
struct nssl_vol_id *volume;

mess=(struct wsr88d_mess *)buffer;
volume=(struct nssl_vol_id *)(buffer+sizeof(struct wsr88d_mess));

mess->hd_size=sizeof(struct wsr88d_mess)+sizeof(unsigned short);		/* set message size		  */
if ((mess->hd_size%2) == 1)
	mess->hd_size++;

mess->hd_type=201;								/* set message type		  */
mess->hd_sequence=0;								/* set sequence number		  */
mess->hd_date=(local_stat.volume_time/86400)+1;					/* set modified julian date	  */
time_tmp=(local_stat.volume_time%86400)*1000;					/* set time			  */
memcpy((char *)mess->hd_time,(char *)&time_tmp,sizeof(unsigned long));
mess->hd_count=1;								/* set the total segment count	  */
mess->hd_segment=1;								/* set the segment number	  */

volume->vol_number=vol;								/* set the volume number	  */

wideband_write(buf,(char *)mess,mess->hd_size*2);			/* write the message		  */
return(0);
}

/* end of file */
