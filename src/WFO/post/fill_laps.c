#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <netcdf.h>
#include "fill_bigfile.h"

/*************************************************************************/
int processLAPS(char *r_filename, time_t reftime, time_t valtime, 
                int numElements, PARAMETER_LIST_T *paramInfo)
{
    int i; 
    long x_out, y_out;
    int x,y,numLevels;
    int cdfId_in, cdfId_out, nc_status;
    float *data;
    short *LAPSinv; 
    float *level;
    long index;
    char timeString[MAX_TIMESTRING_LEN];
    char prevFilename[128];
    struct tm *filtime = gmtime (&reftime);
        
    fprintf(stdout,"filetime = %d\n", reftime);
    strftime(timeString, MAX_TIMESTRING_LEN, "%y%j%H00", filtime);
    fprintf(stdout,"WFO filetime %s\n", timeString);
      
    cdfId_out = openOutputWFOfile(&index, reftime, &x_out, &y_out);
    if (cdfId_out == (-1)) {
      return ERROR;
    }
    else {

      cdfId_in = -2;
      prevFilename[0] = '\0';
    
      for(i = 0; i < numElements; i++) {
        if(extractLAPSdata(&cdfId_in, prevFilename, paramInfo[i], r_filename, 
                           reftime, valtime, &x, &y, &numLevels, &level,
                           &data, &LAPSinv) != SUCCESS) {
          fprintf(stdout,"   Error from extractLAPSdata on LAPS variable %s\n",
            		paramInfo[i].LAPS_cdl_varname);
        }
        else {
          if (((long) x <= x_out) && ((long) y <= y_out)) {
          }
          else {
            fprintf(stdout,"Error in grid sizes - input = %dx%d - output = %dx%d.\n",
                    x,y,x_out,y_out);
            free(level);
            free(LAPSinv);
            free(data);
            nc_status = nc_close(cdfId_in);
            nc_status = nc_close(cdfId_out);
            return ERROR;
          }
         
          if((transferLAPStoWFO(cdfId_out,index,paramInfo[i],reftime,valtime,
              x,y,numLevels,&level,&data,&LAPSinv)) != SUCCESS) {
            fprintf(stdout,"    Error from transferLAPStoWFO on WFO variable %s\n",
                    paramInfo[i].WFO_cdl_varname);
          }
          else { 
            fprintf(stdout,"LAPS %s %s stored as %s in WFO file\n",
             	    paramInfo[i].LAPS_dir_name,
             	    paramInfo[i].LAPS_cdl_varname,
            	    paramInfo[i].WFO_cdl_varname);
          } 
        }
      }

      free(level);
      free(LAPSinv);
      free(data);
      nc_status = nc_close(cdfId_in);
      nc_status = nc_close(cdfId_out);
      return SUCCESS;    
    }   
}
        
/*************************************************************************/
