/******************************************************************************************************************
 *						ip_utilities							  *
 ******************************************************************************************************************/
#include <config.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>

#include <signal.h>

#include "mudp.h"

static char SccsId[]="%W% %G% IP utilities, NSSL";

#define IP_MAX_PACK_LENGTH 1450
#define IP_HEADER_LENGTH 16
#define IP_BUFFER_LENGTH IP_MAX_PACK_LENGTH-IP_HEADER_LENGTH
#define IP_MIN_PACK_LENGTH 128
#define IP_MIN_BUFFER_LENGTH IP_MIN_PACK_LENGTH-IP_HEADER_LENGTH

/**********************************************************************************************************
						ip_send
 **********************************************************************************************************/


ip_send(sock,buffer,length)
int sock;
char *buffer;
int length;
{
static int count=0;
int nlength;
int nbytes;
int size;

struct mudp_header *head;
char *buff;
char ip_buffer[2000];

if (length <= 0)
	return(0);

head=(struct mudp_header *)ip_buffer;
buff=ip_buffer+sizeof(struct mudp_header);

head->mess_length=length;
head->no_packets=length/(IP_MAX_PACK_LENGTH-sizeof(struct mudp_header));
if ((nbytes%(IP_MAX_PACK_LENGTH-sizeof(struct mudp_header))) != 0)
	head->no_packets++;

head->packet_no=1;		
head->packet_offset=0;
head->packet_length=0;
head->flags=0;
head->serial_no=count;

/*
 * loop throuth the number of packets required to send the message
 *
 * if the packet is to small to send over the ethernet interface, pad the packet
 * so that it will go.
 */

while((nlength=length-head->packet_offset) > 0) {
	memmove(buff,buffer+head->packet_offset,(IP_BUFFER_LENGTH < nlength) ? IP_BUFFER_LENGTH : nlength);
	head->packet_length=(IP_BUFFER_LENGTH < nlength) ? IP_BUFFER_LENGTH : nlength;

	if ((size=send(sock,ip_buffer,head->packet_length+sizeof(struct mudp_header),0)) < 0) {
		perror("send");
		return(-1);
		}


/*
 * set up the header for subsequent packets to be sent
 */

	count++;

	head->serial_no=count;
	head->packet_no++;
	head->packet_offset+=head->packet_length;
	}

return(length);
}

/**********************************************************************************************************
						ip_sendto
 **********************************************************************************************************/

ip_sendto(sock,addr,buffer,length)
int sock;
struct sockaddr *addr;
char *buffer;
int length;
{
static int count=0;
int nlength;
int nbytes;
int size;

struct mudp_header *head;
char *buff;
char ip_buffer[2000];

if (length <= 0)
	return(0);

head=(struct mudp_header *)ip_buffer;
buff=ip_buffer+sizeof(struct mudp_header);

head->mess_length=length;
head->no_packets=length/(IP_MAX_PACK_LENGTH-sizeof(struct mudp_header));
if ((nbytes%(IP_MAX_PACK_LENGTH-sizeof(struct mudp_header))) != 0)
	head->no_packets++;

head->packet_no=1;		
head->packet_offset=0;
head->packet_length=0;
head->flags=0;
head->serial_no=count;

/*
 * loop throuth the number of packets required to send the message
 *
 * if the packet is to small to send over the ethernet interface, pad the packet
 * so that it will go.
 */

while((nlength=length-head->packet_offset) > 0) {
	memmove(buff,buffer+head->packet_offset,(IP_BUFFER_LENGTH < nlength) ? IP_BUFFER_LENGTH : nlength);
	head->packet_length=(IP_BUFFER_LENGTH < nlength) ? IP_BUFFER_LENGTH : nlength;

	if ((size=sendto(sock,ip_buffer,head->packet_length+sizeof(struct mudp_header),0,addr,sizeof(*addr))) < 0) {
		perror("send");
		return(-1);
		}
/*
 * set up the header for subsequent packets to be sent
 */

	count=(count+1)%65536;

	head->serial_no=count;
	head->packet_no++;
	head->packet_offset+=head->packet_length;
	}

return(length);
}

/**********************************************************************************************************
						ip_recv
 **********************************************************************************************************/

ip_recv(sock,buffer,length)
int sock;
char *buffer;
int length;
{
static int packets=(-1);

int count;
int nlength;
int nbytes;
int rec_bytes;
int dropped;

struct mudp_header *head;
char *buff;
char ip_buffer[2000];

if (length <= 0)
	return(0);

head=(struct mudp_header *)ip_buffer;
buff=ip_buffer+sizeof(struct mudp_header);

count=0;
rec_bytes=0;
initbuff(buffer,length,0);

/*
 * loop through the packets until a message is seen
 */

while((nlength=recv(sock,ip_buffer,5120,0)) >= 0) {
	if (packets == (-1))
		packets=head->serial_no;
	else  {
		packets=(packets+1)%65536;
		if (packets != head->serial_no && packets < head->serial_no) {
			fprintf(stderr,"ip_recv: PACKETS DROPPED, Dropped: % 5d UDP packets\n",head->serial_no-packets);
			packets=head->serial_no;
			}
		if (packets != head->serial_no && packets > head->serial_no) {
			dropped=(65535-packets)+head->serial_no;
			fprintf(stderr,"ip_recv: PACKETS DROPPED, Dropped: % 5d UDP packets\n",dropped);
			packets=head->serial_no;
			}
		}

	if (nlength != head->packet_length+sizeof(struct mudp_header))
		return(0);

	if (head->packet_offset+head->packet_length > length)
		return(0);

	memmove(buffer+head->packet_offset,buff,head->packet_length);
	rec_bytes+=head->packet_length;
	count++;	
	
	if (head->packet_no == head->no_packets) {
		if (count == head->no_packets)
			return(rec_bytes);
		else
			return(0);
		}
	}
return(-1);
}

/**********************************************************************************************************
						ip_read
 **********************************************************************************************************/

ip_read(sock,buffer,length)
int sock;
char *buffer;
int length;
{
int count;
int nlength;
int nbytes;
int rec_bytes;

struct mudp_header *head;
char *buff;
char ip_buffer[5120];

if (length <= 0)
	return(0);

head=(struct mudp_header *)ip_buffer;
buff=ip_buffer+sizeof(struct mudp_header);

count=0;
rec_bytes=0;
initbuff(buffer,length,0);

/*
 * loop through the packets until a message is seen
 */

while((nlength=read(sock,ip_buffer,5120)) > 0) {
	if (nlength != head->packet_length+sizeof(struct mudp_header))
		return(0);

	if (head->packet_offset+head->packet_length > length)
		return(0);

	memmove(buffer+head->packet_offset,buff,head->packet_length);
	rec_bytes+=head->packet_length;
	count++;	
	
	if (head->packet_no == head->no_packets) {
		if (count == head->no_packets)
			return(rec_bytes);
		else
			return(0);
		}
	}
return(-1);
}

/**********************************************************************************************************
						ip_recvfrom
 **********************************************************************************************************/

ip_recvfrom(sock,addr,buffer,length)
int sock;
struct sockaddr *addr;
char *buffer;
int length;
{
int count;
int nlength;
int nbytes;

struct mudp_header *head;
char *buff;
char ip_buffer[2000];

if (length <= 0)
	return(0);

head=(struct mudp_header *)ip_buffer;
buff=ip_buffer+sizeof(struct mudp_header);

count=0;
initbuff(buffer,length,0);

/*
 * loop through the packets until a message is seen
 */

while((nlength=recvfrom(sock,ip_buffer,2000,0,addr,(unsigned long *) sizeof(*addr))) > 0) {
	memmove(buffer+head->packet_offset,buff,head->packet_length);
	count++;	
	
	if (head->packet_no == head->no_packets) {
		if (count == head->no_packets)
			return(length);
		else
			return(0);
		}
	}
return(-1);
}

