/******************************************************************************************************************
 *							c_buff.h						  *
 *														  *
 *  Version  Description										  Date	  *
 *  -------  -----------------------------------------------------------------------------------------  --------  *
 *    1.00   Preliminary program development, structure used to support c_buff program			04.21.94  *
 *														  *
 ******************************************************************************************************************/

#ifndef _c_buff_h
#define _c_buff_h

#define C_BUFF_MAGIC	0xff00f0a5

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

struct c_buff {						/* circular buffer structure definitions	 	  */
	int	c_shmid;				/* shared memory id associated with circular buffer	  */
	int	c_shmkey;				/* shared memory key associated with circular buffer	  */
	int	c_count;				/* number of times the writer has cycled through 	  */
	int	c_rec_count;				/* number of writes to the buffer			  */
	int	c_highwater;				/* highwater mark, current safe place to write to	  */
	int	c_buffsize;				/* sizeof the buffer in bytes				  */
	};

struct c_buffmess {					/* record header for each message written to the c_buff	  */
	unsigned char	magic[4];			/* value used to identify true beginning of recoed	  */
	unsigned short	mess_id;			/* message id or record number				  */
	unsigned short	mess_size;			/* size of the message in bytes				  */
	unsigned short	mess_checksum;			/* checksum value of the message			  */
	};

#endif

struct c_buff *open_c_buff();
int close_c_buff();

int write_c_buff();
int read_c_buff();

int sync_c_buff();
