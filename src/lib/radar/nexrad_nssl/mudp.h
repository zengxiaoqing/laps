/******************************************************************************************************************
 *							mudp.h							  *
 *														  *
 *  Version  Description										  Date	  *
 *  -------  ----------------------------------------------------------------------------------------  ---------  *
 *    1.00   Preliminary Program development								11.17.93  *
 *														  *
 ******************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <netdb.h>

#include <netinet/in.h>

#define SERVER_PORT 3274				/* arbitary port number for server			  */

#define DUP_MESSAGE	0x0001				/* this message is critcal sent more than once		  */
#define TIME_MESSAGE	0x0002;				/* this message is a time stamp only			  */
#define LOCAL_MESSAGE	0x8000;				/* this is a message sent in the local stream		  */


struct mudp_header {
	unsigned long	serial_no;			/* incremental counter for all packets			  */
	unsigned short	mess_length;			/* length of entire message, even if broken up		  */
	unsigned short	no_packets;			/* number of packets to send message			  */
	unsigned short  packet_no;			/* packet number for this packet ( n of x)		  */
	unsigned short	packet_offset;			/* number of bytes previosly sent in packet		  */
	unsigned short	packet_length;			/* length of the packet minus this header		  */
	unsigned short  flags;				/* specific flags used for various indicators		  */
	};


