#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <netcdf.h>
#include "fill_bigfile.h"

/******************************************************************/
int processFill(char *r_filename, time_t reftime, time_t valtime, 
                int numElements, PARAMETER_LIST_T *paramInfo)
{
    int i, x, y, numLevels, found, fill_status;
    int cdfId_in, cdfId_out, dimId, varId, nc_status;
    float *level, *data;
    short *LAPSinv;
    long index;
#ifdef alpha
    int valtimeMINUSreftime, *vMr_array;
#else
    long valtimeMINUSreftime, *vMr_array;
#endif
    size_t start[1], count[1];
    long n_valtimes, x_out, y_out;
    char timeString[MAX_TIMESTRING_LEN];
    char prevFilename[128];
    struct tm *filtime = gmtime (&reftime);
        
    fprintf(stdout,"unix_filetime = %d\n", reftime);
    fprintf(stdout,"unix_fcsttime = %d\n", valtime);
    fprintf(stdout,"Forecast hour is %d\n",((valtime - reftime)/3600));
    strftime(timeString, MAX_TIMESTRING_LEN, "%y%j%H00", filtime);
    fprintf(stdout,"WFO filetime %s\n", timeString);
      
    cdfId_out = openOutputWFOfile(&index, reftime, &x_out, &y_out); 

    if (cdfId_out == (-1)) {
      return ERROR;
    }
    else {

/* determine which index to write into using valtimeMINUSreftime data from
   output file and calculation of (valtime - reftime) passed in */

#ifdef alpha
      valtimeMINUSreftime = (int)(valtime - reftime);
#else
      valtimeMINUSreftime = (long)(valtime - reftime);
#endif

      nc_status = nc_inq_dimid(cdfId_out, "n_valtimes", &dimId);
      if (nc_status != NC_NOERR) {
        fprintf(stdout,"No dimension 'n_valtimes' in output file.\n ");
        return ERROR;
      }
      else {
        nc_status = nc_inq_dimlen(cdfId_out, dimId, (size_t *) &n_valtimes);
        if (nc_status != NC_NOERR) {
          fprintf(stdout,"Cannot determine value of 'n_valtimes' from output file.\n ");
          return ERROR;
        }
        else {
#ifdef alpha
          vMr_array = (int *)malloc(n_valtimes * sizeof(int));
#else
          vMr_array = (long *)malloc(n_valtimes * sizeof(long));
#endif
          nc_status = nc_inq_varid(cdfId_out, "valtimeMINUSreftime", &varId);
          if (nc_status != NC_NOERR) {
            fprintf(stdout,"No variable 'valtimeMINUSreftime' in output file.\n ");
            free(vMr_array);
            return ERROR;
          }
          else {
            start[0] = 0;
            count[0] = n_valtimes;
#ifdef alpha
            nc_status = nc_get_vara_int(cdfId_out,varId,start,count,vMr_array);
#else
            nc_status = nc_get_vara_long(cdfId_out,varId,start,count,vMr_array);
#endif
            if (nc_status != NC_NOERR) {
              fprintf(stdout, "Cannot read 'valtimeMINUSreftime from output file.\n ");
              free(vMr_array);
              return ERROR;
            }
            else {
              found = FALSE;
              i = 0;
              index = n_valtimes;
              while ((found == FALSE) && (i < n_valtimes)) {
                if (valtimeMINUSreftime == vMr_array[i]) {
                  index = i;
                  found = TRUE;
                }          
                else {
                  i++;
                }
              } 
              if (index >= n_valtimes) {  /* time not found */
                fprintf(stdout,"Requested forecast period not available in output file.\n ");
                fprintf(stdout,"Check entries in valtimeMINUSreftime in template file.\n ");
                free(vMr_array);
                return ERROR;
              }
            }
          }
          free(vMr_array);
        }
      }

      cdfId_in = -2;
      prevFilename[0] = '\0';
    
      for(i = 0; i < numElements; i++) {
        if(extractLAPSdata(&cdfId_in, prevFilename, paramInfo[i], r_filename, 
                           reftime, valtime, &x, &y, &numLevels, &level,
                           &data, &LAPSinv) != SUCCESS) {
          fprintf(stdout,"Error from extractLAPSdata on LAPS variable %s\n",
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

          fill_status = transferLAPStoWFO(cdfId_out,index,paramInfo[i],reftime,
                                          valtime, x,y,numLevels,&level,&data,
                                          &LAPSinv);
          if (fill_status != SUCCESS) {
            if (fill_status == NOINV) {
              fprintf(stdout,"No inventory on one or more levels of LAPS variable %s.\n",
                      paramInfo[i].LAPS_cdl_varname);
            }
            else {
              fprintf(stdout,"Error from transferLAPStoWFO on WFO variable %s\n",
                      paramInfo[i].WFO_cdl_varname);
            }
          }
          else { 
            fprintf(stdout,"     LAPS %s %s stored as %s in WFO file\n",
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
        

/******************************************************************/
