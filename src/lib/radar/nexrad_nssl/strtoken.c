/**********************************************************************************************************************************
 *								strtoken.c							  *
 *																  *
 *  Version  Description												  Date	  *
 *  -------  ---------------------------------------------------------------------------------------------------------  --------  *
 *    1.00   Preliminary program development, convert strings to token values						05.09.94  *
 *																  *
 **********************************************************************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>

static char SccsId[]="@(#)strtoken.c	1.1  8/11/94  String Utility, Douglas Rhue";

/**********************************************************************************************************************************
						strtoken
 **********************************************************************************************************************************/

strtoken(token,list,values)
char *token;
char **list;
int *values;
{
int offset;
int loop;
int count;

count=0;
for (loop=0; *(list+loop) != '\0' && (list+loop) != '\0' && strlen(*(list+loop)) > 0; loop++) {
	if (strncmp(list[loop],token,strlen(token)) == 0) {
		count++;	
		offset=loop;
		}
	}

if (count != 1)
	return(-1);

return(values[offset]);
}

