/******************************************************************************************************************
 *							c_buff.c						  *
 *														  *
 *  Version  Description										  Date	  *
 *  -------  -----------------------------------------------------------------------------------------  --------  *
 *    1.00   Preliminary program development								02-05-93  *
 *    2.00   Enhancements added by Ben Stevens to match new c_buff package				02-08-94  *
 *    2.01   Further enhancemnts by Ben Stevens								02-29-94  *
 *    3.00   Rewrite to simplfy interface, no reason for the number of calls				04-20-94  * 
 *														  *
 ******************************************************************************************************************/

#include <config.h>
#include <errno.h>
#include "c_buff.h"

#ifdef SHM_W
#else
#define SHM_W 0200
#endif

#ifdef SHM_R
#else
#define SHM_R 0400
#endif

static char SccsId[]="@(#)c_buff.c	1.9	5/31/94";

/******************************************************************************************************************
						open_c_buff
 ******************************************************************************************************************/

struct c_buff *open_c_buff(key,size,mode)			/* open or attach to a circular buffer segment	  */
int key;							/* key for attachment to shared memory seg	  */
int size;							/* size of the segment in bytes			  */
int mode;							/* mode of the shared memory segment		  */
{
int id;
int uid;
int gid;
int mode_fl;
int granted;

char *pnt;					
struct c_buff *buff;
struct shmid_ds shrmem;
 /**** char *shmat(); ****/

uid=getuid();
gid=getgid();

granted=0;

if ((id=shmget((key_t)key,size+sizeof(struct c_buff),mode)) < 0)	/* get shared memory id			  */
	return('\0');							/* error, could not create or locate 	  */

if (shmctl(id,IPC_STAT,&shrmem) != 0)
	return('\0');

if (mode == 0)
	granted=SHM_RDONLY;

else if (uid == shrmem.shm_perm.uid) {
	mode_fl=shrmem.shm_perm.mode & 0700;
	if ((mode_fl&mode&SHM_W) == SHM_W)
		granted|=SHM_W;
	if ((mode_fl&mode&SHM_R) == SHM_R)
		granted|=SHM_W;
	}

else if (gid == shrmem.shm_perm.gid) {
	mode_fl=(shrmem.shm_perm.mode * 070) << 3;
	if ((mode_fl&mode&SHM_W) == SHM_W)
		granted|=SHM_W;
	if ((mode_fl&mode&SHM_R) == SHM_R)
		granted|=SHM_W;
	granted>>=3;
	}

else {
	mode_fl=(shrmem.shm_perm.mode * 07) << 6;
	if ((mode_fl&mode&SHM_W) == SHM_W)
		granted|=SHM_W;
	if ((mode_fl&mode&SHM_R) == SHM_R)
		granted|=SHM_W;
	granted>>=6;
	}


if ((pnt=shmat(id,0,granted))==(char *)(-1)) 					/* attach to shared memory segment	  */
	return('\0');							/* cannot attach to shared memory 	  */

if ((mode&IPC_CREAT)==IPC_CREAT) {					/* create a new circular buffer segment	  */
	buff=(struct c_buff *)pnt;					/* map pointer to c_buff structure	  */
	buff->c_shmid=id;						/* retain the shared memory id		  */
	buff->c_shmkey=key;						/* set the key, just for argument	  */
	buff->c_count=0;						/* reset the cycle counter		  */
	buff->c_rec_count=0;
	buff->c_highwater=0;						/* reset the highwater marker		  */
	buff->c_buffsize=size;						/* set the buffer size in bytes		  */
	initbuff(pnt+sizeof(struct c_buff),size,0);			/* init the circular buffer		  */
	}

return((struct c_buff *)pnt);						/* return pointer to c_buff		  */
}

/******************************************************************************************************************
						write_c_buff
 ******************************************************************************************************************/
		
int write_c_buff(circ_buff,buffer,length)
struct c_buff *circ_buff;
char *buffer;
int length;
{
struct c_buff local_buff;
char *buff;

int highwater;

struct c_buffmess header;
unsigned long magic;
unsigned long checksum();

if (length == 0)
	return(0);

if (length > 32768) {
	errno=ENOMEM;
	return(-1);
	}
	

memcpy((char *)&local_buff,(char *)circ_buff,sizeof(struct c_buff));
buff=((char *)circ_buff)+sizeof(struct c_buff);

magic=C_BUFF_MAGIC;
memcpy((char *)header.magic,(char *)&magic,sizeof(unsigned long));
header.mess_id=circ_buff->c_rec_count;
header.mess_size=length;
header.mess_checksum=(unsigned short)(checksum(buffer,length)&0x0ffff);

highwater=circ_buff->c_highwater;	

if (_write_c_seg(&local_buff,buff,(char *)&header,sizeof(struct c_buffmess)) != sizeof(struct c_buffmess)) {
	local_buff.c_highwater=highwater;
	return(-1);
	}

if (_write_c_seg(&local_buff,buff,buffer,length) < 0) {
	local_buff.c_highwater=highwater;
	return(-1);
	}

circ_buff->c_count=local_buff.c_count;
circ_buff->c_rec_count=local_buff.c_rec_count;
circ_buff->c_highwater=local_buff.c_highwater;
return(length);
}

/******************************************************************************************************************
						_write_c_seg
 ******************************************************************************************************************/

_write_c_seg(c_buff,pnt,buffer,size)
struct c_buff *c_buff;
char *pnt;
char *buffer;
int size;
{
int spacetoend;							/* distance from highwater to the eob		  */

if (size > (c_buff->c_buffsize) || size < 0) {			/* check the size of the message		  */
	errno=EINVAL;
	return(-1);
	}

spacetoend=c_buff->c_buffsize-c_buff->c_highwater;		/* determine the distance to the eob		  */

if (size > spacetoend) {					/* is the message greater than the free space	  */
	memcpy(pnt+c_buff->c_highwater,buffer,spacetoend);	/* copy first part of the message to c buffer	  */
	memcpy(pnt,buffer+spacetoend,size-spacetoend);		/* copy the second part of the message		  */
	}
else {								/* the message is less than the free space	  */
	memcpy(pnt+c_buff->c_highwater,buffer,size);		/* copy the entire message to the c buffer	  */
	}

if (c_buff->c_highwater > (c_buff->c_highwater+size)%c_buff->c_buffsize) 	/* did highwater rollover	  */
	c_buff->c_count++;							/* highwater rolled over	  */
	
c_buff->c_highwater=(c_buff->c_highwater+size)%c_buff->c_buffsize;		/* adjust highwater mark	  */
c_buff->c_rec_count++;								/* adjust record counter	  */
return(size);
}
/******************************************************************************************************************
						read_c_buff
 ******************************************************************************************************************/
		
int read_c_buff(circ_buff,buffer,length,lowwater)
struct c_buff *circ_buff;
char *buffer;
int length;
int *lowwater;
{
int tmp_low;
int tmp_length;
int lowwater_mark;

unsigned long magic;
struct c_buffmess header;
unsigned long checksum();
unsigned short tmp_check;

char tmp_buffer[32768];						/* temperal buffer for local use		  */

if (length==0)
	return(0);

lowwater_mark=(*lowwater);

if (lowwater_mark == circ_buff->c_highwater)
	return(0);

if (length > 32768) {
	errno=ENOMEM;
	return(-1);
	}

if (_read_c_seg(circ_buff,(char *)&header,sizeof(struct c_buffmess),&lowwater_mark) != sizeof(struct c_buffmess)) {
	errno=EBADMSG;
	return(-1);
	}

memcpy((char *)&magic,(char *)header.magic,sizeof(unsigned long));
if (magic != C_BUFF_MAGIC) {
	errno=EBADMSG;
	return(-1);
	}

if (_read_c_seg(circ_buff,tmp_buffer,header.mess_size,&lowwater_mark) != header.mess_size) {
	errno=EBADMSG;
	return(-1);
	}

if (length >= header.mess_size) {
	tmp_length=header.mess_size;
	errno=0;
	}
else {
	tmp_length=length;
	errno=EMSGSIZE;
	}

tmp_check=(unsigned short)(checksum(tmp_buffer,header.mess_size)&0xffff);

if (tmp_check == header.mess_checksum) {
	memcpy(buffer,tmp_buffer,tmp_length);
	*lowwater=lowwater_mark;
	return(tmp_length);
	}
else {
	errno=EBADMSG;
	return(-1);
	}
}

/******************************************************************************************************************
						_read_c_seg
 ******************************************************************************************************************/

_read_c_seg(c_buff,buffer,length,lowwater)
struct c_buff *c_buff;
char *buffer;
int length;
int *lowwater;
{
int spacetoend;							/* distance from highwater to the eob		  */
char *pnt;						/* pointer to the acual buffer			  */

if (length > (c_buff->c_buffsize) || length < 0) {		/* check the size of the message		  */
	errno=EINVAL;
	return(-1);
	}

if (*lowwater > c_buff->c_buffsize) {				/* lowwater marker out of range			  */
	errno=EINVAL;						
	return(-1);
	}

pnt=((char *)c_buff)+sizeof(struct c_buff);			/* determine where the buffer really is		  */
spacetoend=c_buff->c_buffsize-(*lowwater);			/* determine the distance to the eob		  */

if (length > spacetoend) {					/* is the message greater than the free space	  */
	memcpy(buffer,pnt+(*lowwater),spacetoend);		/* copy first part of the message to buffer	  */
	memcpy(buffer+spacetoend,pnt,length-spacetoend);		/* copy the second part of the message		  */
	}

else {								/* the message is less than the free space	  */
	memcpy(buffer,pnt+(*lowwater),length);			/* copy the entire message to the c buffer	  */
	}

	
*lowwater=((*lowwater)+length)%c_buff->c_buffsize;		/* adjust highwater mark	 		 */
return(length);
}


/******************************************************************************************************************
						close_c_buff
 ******************************************************************************************************************/
		
close_c_buff(circ_buff)
struct c_buff *circ_buff;
{
struct shmid_ds buff;
int id;

id=circ_buff->c_shmid;

if (shmdt((char *)circ_buff) != 0)
	return(-1);

if (shmctl(id,IPC_STAT,&buff) != 0)
	return(-1);

if (buff.shm_nattch == 0) {					/* nothing is attached here, delete segment	  */
	if (shmctl(id,IPC_RMID,&buff) != 0)
		return(1);
	}
return(0);
}


/******************************************************************************************************************
						sync_c_buff
 ******************************************************************************************************************/
		
sync_c_buff(circ_buff)
struct c_buff *circ_buff;
{
return(circ_buff->c_highwater);
}
