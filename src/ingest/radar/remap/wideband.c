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
/**********************************************************************************************************************************
 *								wideband.c							  *
 *																  *
 *  Version  Description												  Date	  *
 *  -------  ---------------------------------------------------------------------------------------------------------  --------  *
 *    1.00   Program Development, program facilitates moving NEXRAD wideband data between processes, (c_buff)		02.13.94  *
 *    1.02   Enhancement to allow access to more that one wideband message handler					04.20.94  *
 *																  *
 **********************************************************************************************************************************/

#include "c_buff.h"
#include "riscrpg.h"

/*********************************************************************************************************************************
							wideband_init_reader
 *********************************************************************************************************************************/

struct c_buff *wideband_init_reader(key,pointer)
int key;
int *pointer;
{
struct c_buff *pnt;

if ((pnt=open_c_buff(key,0,0)) == '\0')
	return('\0');

*pointer=sync_c_buff(pnt);
return(pnt);
}

/*********************************************************************************************************************************
							wideband_init_writer
 *********************************************************************************************************************************/

struct c_buff *wideband_init_writer(key,size)
int key;
int size;
{
struct c_buff *pnt;

if ((pnt=open_c_buff(key,size,0644|IPC_CREAT)) == '\0')
	return(0);
return(pnt);
}

/*********************************************************************************************************************************
							wideband_read
 *********************************************************************************************************************************/

int wideband_read(circ_buff, msg_buffer, msg_buflen, read_pointer )
struct c_buff *circ_buff;
char *msg_buffer;
int msg_buflen;
int *read_pointer;
{
return( read_c_buff( circ_buff, msg_buffer, msg_buflen, read_pointer ) );
}

/*********************************************************************************************************************************
							wideband_write
 *********************************************************************************************************************************/

int wideband_write(circ_buff, msg_buffer, msg_len)
struct c_buff  *circ_buff;
char *msg_buffer;
int	msg_len;
{
return( write_c_buff( circ_buff, msg_buffer, msg_len ) );
}

/*********************************************************************************************************************************
							wideband_done_reader
 *********************************************************************************************************************************/

int wideband_done_reader(circ_buff)
struct c_buff *circ_buff;
{
return(close_c_buff(circ_buff));
}

/*********************************************************************************************************************************
							wideband_done_reader
 *********************************************************************************************************************************/

int	wideband_done_writer(circ_buff)
struct c_buff *circ_buff;
{
return(close_c_buff(circ_buff));
}
