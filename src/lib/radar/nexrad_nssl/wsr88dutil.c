/**************************************************************************************************************************
 *								wsr88dutil.c						  *
 *															  *
 *  Description													  Date	  *
 *  ----------------------------------------------------------------------------------------------------------  --------  *
 *  A collection of routines usefull when processing wsr88d data						05.01.94  *
 *															  *
 **************************************************************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include <sys/types.h>
#include <time.h>

static char SccsId[]="@(#)wsr88dutil.c	3.1  2/10/95  NSSL, RIDDS Nexrad Data Processing Routines";

/**************************************************************************************************************************
						wsrtimetounix
 **************************************************************************************************************************/

time_t wsrtimetounix(wsr_date,wsr_time)
unsigned short wsr_date;
unsigned long wsr_time;
{
unsigned long seconds;
time_t temp;
char tmp[256];

temp=((time_t)wsr_date)-1;
temp*=(time_t)86400;

memcpy((char *)&seconds,(char *)&wsr_time,sizeof(unsigned long));

seconds/=1000;
temp+=seconds;

return(temp);
}

/**************************************************************************************************************************
						wsrangle
 **************************************************************************************************************************/

double wsrangle(angle)
unsigned short angle;
{
int value;
double angles;

value=(int)angle;
value>>=3;
angles=0.043945*(double)value;
angles*=100;

return(angles);
}
