/******************************************************************************
 *				checksum.c				      *
 *									      *
 *  Version  Description					      Date    *
 *  -------  -----------------------------------------------------  --------  *
 *    1.00   Preliminary program development,                                 *
 *           various checksum calculating routines		    01.26.94  *
 *									      *
 *****************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>

/******************************************************************************
				Checksum
******************************************************************************/

unsigned long checksum(buffer,length)
unsigned char *buffer;
int length;
    {
    int loop;
    unsigned long result;

    result=0;

    for (loop=0; loop<length; loop++) 
	result+=(unsigned long)buffer[loop];

    return(result);
    }

/******************************************************************************
				modsum
******************************************************************************/

unsigned long modsum(buffer,length,mod)
    unsigned char *buffer;
    int length;
    int mod;
    {
    int loop;
    unsigned long result;

    if (mod == 0) 
	return(0);

    result=0;
    for (loop=0; loop < length; loop+=mod)
	result+=(unsigned long)buffer[loop];

    return(result);
    }

/******************************************************************************
				fdchecksum
******************************************************************************/

unsigned long fdchecksum(fd)
int fd;
    {
    char buffer[5120];

    unsigned long checksum();
    unsigned long result;

    int size;

    result=0;
    while((size=read(fd,buffer,5120)) > 0) 
	result+=checksum(buffer,size);
	
    return(result);
    }

/******************************************************************************
				fdmodsum
******************************************************************************/

unsigned long fdmodsum(fd,mod)
    int fd;
    int mod;
    {
    char buffer[5120];
    int size;

    unsigned long modsum();
    unsigned long result;

    if (mod > 512)
	return(0);

    result=0;
    while((size=read(fd,buffer,mod*10)) > 0)
	result+=modsum(buffer,size,mod);

    return(result);
    }







