#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include "fill_bigfile.h"

/* accepts three filename formats:
     YY is year 
     JJJ is julian day 
     HH is hour model is run    
     MM is minutes past top of hour that model is run
     hh is hour beyond HH that data is valid (ie. 00, 01, 02 ...)
     mm is minutes past top of hour that data is valid
     OO is month
     DD is day 

     1. YYJJJHHhh  
     2. YYJJJHHMMhhmm
     3. YYYYOODD_HHhh (wfo filename format as of 7/98)
*/
int main(int argc, char *argv[])
{
       time_t refTime, validTime;
       char filename[100];

       if (argc == 2) {   /* transfer argv[1] to filename */

         strcpy(filename,argv[1]);
         filenameToUnixtime(&refTime, &validTime, filename);
         printf("%12d %12d\n",refTime, validTime); 
       }
       else {
	 fprintf(stderr,"Program fname2unixtime could not convert filetime %s to unixtime\n",
	         filename);
       }
       return;
}
