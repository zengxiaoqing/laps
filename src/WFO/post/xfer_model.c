#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <netcdf.h>
#include "fill_bigfile.h"

/* Usage: This program expects that the following environment variables are set:
	    XFR_HOME  - location of new_xfer_rams executable.
            TRANS_TBL - full path and file name of table showing what to transfer.
	    FCSTPRD   - location of /ram and /rsf directories where input files reside
                        (or whatever directories are shown in $TRANS_TBL).
            WFO_FCSTPRD - location of bigfile output directory, which must have 
                          a template file in it.

          If you want to process exactly one file, pass in one argument:
            filename to process (with no path, and no extension) that is located 
              in $FCSTPRD

          If this program is run with no arguments, it will look in the
            $FCSTPRD directory, find the most recent model run, and process
            the analysis and all forecasts available for that model run.
*/

main(int argc, char *argv[])
{
       int xfr_status, numElements;
       PARAMETER_LIST_T paramInfo[MAX_PARAMETER_LIST_LEN];
       FILE *fp;
       time_t unix_filetime, unix_fcsttime, curr_time;
       DIR *inputDir;
       long l_closestEnt, l_startFile, l_currtime, last_hour;
       struct dirent dir_ent, *p_dirent;
       char filename[100], fcstRunTime[20], ext[4], *fn_dot;
       char currtime[20], startFilename[20], closestEnt[20];
       char paramString[100], inputDirname[128];
       char c_max_fcst_hr[4], base_filename[132];
       int itemp, found;
       long next_hour, max_fcst_hr, fcst_period;
       char c_next_hour[3], char1, char2;


/* open $TRANS_TBL in $XFR_HOME to read first extension */

       strcpy(filename,getenv("TRANS_TBL"));

       if ((fp = fopen(filename, "r")) == (FILE *) 0) {
         /* fprintf(stderr, "s", strerror(errno));
         fprintf(stderr, "\n"); */
         fprintf(stdout, "Unable to open %s file.\n", filename);
       }
       else {

/* read in parameter file for passing into processFill */
         get_transfer_parameter_list(paramInfo, &numElements);

/* get extension for first file of xfer from parameter list */
         strcpy(ext,paramInfo[0].LAPS_dir_name);

/* modified 5-11-98 to accept parameters of filename(no path), model run unixtime,
     and data unixtime. */ 
/* modified again on 7-23-98 to accept one parameter only, the filename(no path) */

         if (argc == 2) {  /* process exactly one file defined by fcstRunTime */
	   strcpy(fcstRunTime,argv[1]);
           filenameToUnixtime((time_t *) &unix_filetime, &unix_fcsttime, fcstRunTime);

           xfr_status = processFill(fcstRunTime, unix_filetime, 
                                    unix_fcsttime, numElements, paramInfo);
         }
         else {  /* do it the old way, which will loop through all files in FCSTPRD  */
/* make a directory that is $FCSTPRD/ext, read entries and determine most recent filename */
           strcpy(inputDirname,getenv("FCSTPRD"));
           strcat(inputDirname,"/");
           strcat(inputDirname,ext);

           inputDir = opendir(inputDirname);
           if (inputDir == NULL) {
             fprintf(stdout, "Unable to open directory %s.\n", inputDirname);
             fprintf(stdout, "Cannot access LAPS/ForecastModel files");
             return ERROR;
           }
           else {

/* get current UTC time to compare to, set closestEnt to 0 */
             curr_time = time(NULL);
             l_currtime = (long)curr_time;
             closestEnt[0] = '\0';
             strcat(closestEnt,"0");
             l_closestEnt = 0;

/* find the file in the directory that is closest to currtime_hr" */

             p_dirent = &dir_ent;
             p_dirent = readdir(inputDir);
             startFilename[19] = '\0';
             while (p_dirent != NULL) {
               strcpy(startFilename,p_dirent->d_name);
               if (strncmp(startFilename,".",1) == 0) {}
               else {
                 fn_dot = strchr(startFilename,46); /* 46 is ascii for "." */
                 if (fn_dot != NULL) {
                   *fn_dot = '\0';
                 }
                 filenameToUnixtime((time_t *) &l_startFile,&unix_fcsttime,startFilename);
                 if ((l_startFile < l_currtime) && (l_startFile > l_closestEnt)) {
                   strcpy(closestEnt,startFilename);
                   l_closestEnt = l_startFile;
                 }
               }
               p_dirent = readdir(inputDir);
             }

/* when get to here, should have unixtime of most recent model run in l_closestEnt */
/* loop through directory entries again and process any files that have
   a unix_filetime that matches l_closestEnt */

             rewinddir(inputDir);
             p_dirent = readdir(inputDir);
             fcstRunTime[19] = '\0';
             while (p_dirent != NULL) {
               strcpy(fcstRunTime,p_dirent->d_name);
               if (strncmp(fcstRunTime,".",1) == 0) {}
               else {
                 fn_dot = strchr(fcstRunTime,46); /* 46 is ascii for "." */
                 if (fn_dot != NULL) {
                   *fn_dot = '\0'; /* get rid of file extension */
                 }
                 filenameToUnixtime((time_t *) &unix_filetime,&unix_fcsttime,fcstRunTime);
                 if (unix_filetime == l_closestEnt) {
                   xfr_status = processFill(fcstRunTime, unix_filetime, 
                                            unix_fcsttime, numElements, paramInfo);
                 } 
               }
               p_dirent = readdir(inputDir);
             }
           }
         }
       }
       exit (0);
}
