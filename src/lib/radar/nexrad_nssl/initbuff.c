/********************************************************************************************************************************
 *							initbuff()								*
 *																*
 * Version  Description												  Date	  Int	*
 * -------  --------------------------------------------------------------------------------------------------	--------  ---	*
 * 1.00.00  Preliminary program development									11.14.89  DTR	*
 *																*
 * This routine will initialize a buffer of a given length with a given character. This routine is generrally used to remove	*
 * garbage from the buffer before it is to be used.										*
 *																*
 * This routine returns the number of bytes initialized.									*
 *																*
 ********************************************************************************************************************************/

#include <config.h>

int initbuff(buffer,size,space)
char *buffer;							/* address of buffer to initialize				*/
int size;							/* size of the buffer in bytes					*/
char space;							/* character to initilize the buffer with			*/
{
int loop;							/* counter for buffer size					*/

for (loop=0; loop < size; loop++) 
	*(buffer+loop) = space;					/* set point in buffer to equal to the init character		*/

return(loop+1);
}

