/******************************************************************************************************************
 *							strtime							  *
 *														  *
 *  Version  Description										  Date	  *
 *  -------  -----------------------------------------------------------------------------------------  --------  *
 *  1.00.00  Preliminary program development								03.17.93  *
 *														  *
 ******************************************************************************************************************/

#include <config.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

static char SccsId[]="@(#)strtime.c	1.1  11/7/94, string to time and time to string functions";


/******************************************************************************************************************
						strtotime
 ******************************************************************************************************************/

time_t strtotime(buff,format)
char *buff;
int format;
{
struct tm sortanow;
char *strptime();
time_t now;
char *tmp;
char tmpbuf[256];

switch(format) {
	case 1:	tmp=strptime(buff,"%m/%d/%y %H:%M:%S",&sortanow);
		sortanow.tm_wday=0;
		sortanow.tm_yday=0;
		sortanow.tm_isdst=0;
		now=timegm(&sortanow);
		return(now);
		break;

	case 2:	tmp=strptime(buff,"%m/%d %H:%M",&sortanow);
		sortanow.tm_sec=0;
		sortanow.tm_year=70;
		sortanow.tm_wday=0;
		sortanow.tm_yday=0;
		sortanow.tm_isdst=0;
		break;

	case 3:	tmp=strptime(buff,"%m/%d/%y",&sortanow);
		sortanow.tm_sec=0;
		sortanow.tm_min=0;
		sortanow.tm_hour=0;
		sortanow.tm_wday=0;
		sortanow.tm_yday=0;
		sortanow.tm_isdst=0;
		break;

	case 4: tmp=strptime(buff,"%m/%d",&sortanow);
		sortanow.tm_sec=0;
		sortanow.tm_min=0;
		sortanow.tm_hour=0;
		sortanow.tm_year=70;
		sortanow.tm_wday=0;
		sortanow.tm_yday=0;
		sortanow.tm_isdst=0;
		break;

	case 5: tmp=strptime(buff,"%m/%y",&sortanow);
		sortanow.tm_sec=0;
		sortanow.tm_min=0;
		sortanow.tm_hour=0;
		sortanow.tm_mday=1;
		sortanow.tm_wday=0;
		sortanow.tm_yday=0;
		sortanow.tm_isdst=0;
		break;

	case 6:	tmp=strptime(buff,"%H:%M:%S",&sortanow);
		sortanow.tm_mday=1;
		sortanow.tm_mon=0;
		sortanow.tm_year=70;
		sortanow.tm_wday=0;
		sortanow.tm_yday=0;
		sortanow.tm_isdst=0;
		break;

	case 7: tmp=strptime(buff,"%H:%M",&sortanow);
		sortanow.tm_sec=0;
		sortanow.tm_mday=1;
		sortanow.tm_mon=0;
		sortanow.tm_year=70;
		sortanow.tm_wday=0;
		sortanow.tm_yday=0;
		sortanow.tm_isdst=0;
		break;

	case 8:	tmp=strptime(buff,"%y/%m/%d %H:%M:%S",&sortanow);
		sortanow.tm_wday=0;
		sortanow.tm_yday=0;
		sortanow.tm_isdst=0;
		break;
	}	
			
now=timegm(&sortanow);
return(now);
}

/******************************************************************************************************************
						datetotime
 ******************************************************************************************************************/

time_t datetotime(buff)
char *buff;
{
char buffer[32];
time_t strtotime();
time_t now;

strcpy(buffer,buff);
strcat(buffer," 00:00:00");

now=strtotime(buffer);
return(now);
}

/******************************************************************************************************************
						hourtotime
 ******************************************************************************************************************/

time_t hourtotime(buff)
char *buff;
{
char buffer[32];
time_t strtotime();
time_t now;

strcpy(buffer,"70-01-01 ");
strcat(buffer,buff);

now=strtotime(buffer);
return(now);
}

/******************************************************************************************************************
						timetostr()
 ******************************************************************************************************************/

char *timetostr(now,format)					/* convert a time value to various strings	  */
time_t *now;							/* time to be converted (GMT value)		  */
int format;							/* format code, used to determine output string	  */
{
struct tm *tm;
struct tm *gmtime();

char buffer[256];

tm=gmtime(now);

switch(format) {	
	case 1:	sprintf(buffer,"%02d/%02d/%02d ",tm->tm_mon+1,tm->tm_mday,tm->tm_year);	/* MM/DD/YY HH:MM:SS	  */
		sprintf(buffer+strlen(buffer),"%02d:%02d:%02d",tm->tm_hour,tm->tm_min,tm->tm_sec);
		return(buffer);
		break;

	case 2: sprintf(buffer,"%02d/%02d ",tm->tm_mon+1,tm->tm_mday);			/* MM/DD HH:MM		  */
		sprintf(buffer+strlen(buffer),"%02d:%02d",tm->tm_hour,tm->tm_sec);
		return(buffer);
		break;

	case 3:	sprintf(buffer,"%02d/%02d/%02d",tm->tm_mon+1,tm->tm_mday,tm->tm_year);	/* MM/DD/YY		  */
		return(buffer);
		break;

	case 4: sprintf(buffer,"%02d/%02d",tm->tm_mon+1,tm->tm_mday);			/* MM/DD		  */
		return(buffer);
		break;

	case 5: sprintf(buffer,"%02d/%02d",tm->tm_mon+1,tm->tm_year);			/* MM/YY		*/
		return(buffer);
		break;

	case 6: sprintf(buffer,"%02d:%02d:%02d",tm->tm_hour,tm->tm_min,tm->tm_sec);	/* HH:MM:SS		*/
		return(buffer);
		break;

	case 7:	sprintf(buffer,"%02d:%02d",tm->tm_hour,tm->tm_min);			/* HH:MM		*/
		return(buffer);
		break;

	case 8: sprintf(buffer,"%02d/%02d/%02d ",tm->tm_year,tm->tm_mon+1,tm->tm_mday); /* YY/MM/DD HH:MM:SS	   */
		sprintf(buffer+strlen(buffer),"%02d:%02d:%02d",tm->tm_hour,tm->tm_min,tm->tm_sec);
		return(buffer);
		break;
	}
}
