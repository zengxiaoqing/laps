#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <netcdf.h>
#include "fill_bigfile.h"

main()
{
       int xfr_status, numElements, xstatus, fn_len;
       PARAMETER_LIST_T paramInfo[MAX_PARAMETER_LIST_LEN];
       FILE *fp;
       time_t unix_filetime;
       char filename[100], fcstRunTime[20], ext[4];
       char paramString[100], inputDirname[128];
       char timestr[21], namestr[21], *nptr;
       long i4time;


/* open $TRANS_TBL in $XFR_HOME to read first extension */

       strcpy(filename,getenv("TRANS_TBL"));

       if ((fp = fopen(filename, "r")) == (FILE *) 0) {
         /* fprintf(stderr, "s", strerror(errno));
         fprintf(stderr, "\n"); */
         fprintf(stdout, "Unable to open %s file.\n", filename);
       }
       else {

/* read in parameter file for passing into processLAPS */
         xstatus = get_transfer_parameter_list(paramInfo, &numElements);

/* open laps/nest7grid/sched/systime.dat */

         strcpy(filename,getenv("LAPSTIME"));

         if ((fp = fopen(filename, "r")) == (FILE *) 0) {
        /*   fprintf(stderr, "s", strerror(errno));
           fprintf(stderr, "\n"); */
           fprintf(stdout, "Unable to open LAPS %s file.\n", filename);
         }
         else {

/* get i4time and filename of LAPS from file */

           if (!(feof(fp))) {
             fgets(timestr,20,fp);
             fgets(namestr,20,fp);
             fclose(fp);

             sscanf(timestr,"%ld", &i4time);

             nptr = namestr;
/* strip leading blanks and trailing newline from namestr */
             while (strncmp(nptr," ",1) == 0) nptr++;
             strncpy(fcstRunTime,nptr,LAPS_FN_LEN);
             fn_len = strlen(fcstRunTime);
             if (fcstRunTime[fn_len - 1] == '\n')
	       fcstRunTime[fn_len - 1] = '\0'; 
             else
               fcstRunTime[fn_len] = '\0';

/* convert i4time with base of 1/1/60 to unixtime with base of 1/1/70 */
             unix_filetime = i4time - 315619200;

             xfr_status = processLAPS(fcstRunTime, unix_filetime, 
                                      unix_filetime, numElements, paramInfo);
           }
         }
       }
       exit (0);
}
