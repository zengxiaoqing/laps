/**********************************************************************************************************************************
 *								wideband.c							  *
 *																  *
 *  Version  Description												  Date	  *
 *  -------  ---------------------------------------------------------------------------------------------------------  --------  *
 *    1.00   Program Development, program facilitates moving NEXRAD wideband data between processes, (c_buff)		02.13.94  *
 *    1.02   Enhancement to allow access to more that one wideband message handler					04.20.94  *
 *																  *
 **********************************************************************************************************************************/

#include <config.h>
#include "c_buff.h"
#include "ridds.h"
/***#include "riscrpg.h" ***/

static char SccsId[]="@(#)wideband.c	3.1 8/25/94 Wideband Circular buffer resources";


int wideband_timeout=5000;

/*********************************************************************************************************************************
							wideband_init_reader
 *********************************************************************************************************************************/

struct c_buff *wideband_init_reader(key,pointer)
int key;
int *pointer;
{
struct c_buff *pnt;

if ((pnt=open_c_buff(key,0,SHM_RDONLY)) == '\0')
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

if ((pnt=open_c_buff(key,size,0666|IPC_CREAT)) == '\0')
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
int size;
extern int wideband_timeout;

if ((size=read_c_buff(circ_buff,msg_buffer,msg_buflen,read_pointer)) == 0) {
	usleep(wideband_timeout);
	return(size);
	}

else if (size < 0) {
	if (errno == EBADMSG) {
		*read_pointer=sync_c_buff(circ_buff);
		return(0);
		}
	else {
		perror("wideband read");
		return(-1);
		}
	}

return(size);
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
