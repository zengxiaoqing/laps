/**************************************************************************************************
 *						messages					  *
 **************************************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

/**************************************************************************************************
						init_mess
 **************************************************************************************************/

init_mess(nel,width,buff)
int nel;
int width;
char *buff;
{
return(initbuff(buff,nel*width,0));
}

/**************************************************************************************************
						add_mess
 **************************************************************************************************/

add_mess(nel,width,buff,mess)
int nel;
int width;
char *buff;
char *mess;
{
memcpy(buff,buff+width,(nel-1)*width);
initbuff(buff+((nel-1)*width),width,0);
strncpy(buff+((nel-1)*width),mess,width);
}

/**************************************************************************************************
						del_mess
 **************************************************************************************************/

del_mess(nel,width,buff,element)
int nel;
int width;
char *buff;
int element;
{
int loop;

if (element < 0 || element > nel-1)
	return(-1);

for (loop=element; loop > 0; loop--)
	memcpy(buff+(loop*width),buff+((loop-1)*width),width);

initbuff(buff,width,0);
return(0);
}

/**************************************************************************************************
						extract_mess
 **************************************************************************************************/

extract_mess(element,mess,nel,width,buff)
int element;
char *mess;
int nel;
int width;
char *buff;
{
if (element < 0 || element > nel-1)
	return(-1);

strncpy(mess,buff+(element*width),width);
*(mess+width)='\0';
return(0);
}

