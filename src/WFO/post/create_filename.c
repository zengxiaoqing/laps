#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "/usr/local/netcdf/include/netcdf.h"
#include "fill_laps.h"

#define TRUE    1
#define FALSE   0
#define SUCCESS TRUE
#define ERROR   FALSE

/******************************************************************/
main(int argc, char *argv[])
{
    char unixtime_in[MAX_TIMESTRING_LEN], WFO_filename[MAX_TIMESTRING_LEN];
    time_t unixtime_t;
    struct tm *unixtime_tm;


    strcpy(unixtime_in, argv[1]);
    unixtime_t = (time_t)atol(unixtime_in);
    unixtime_tm = gmtime (&unixtime_t);
    
    strftime(WFO_filename, MAX_TIMESTRING_LEN, "%Y%m%d_%H00", unixtime_tm);
    printf("%s\n",WFO_filename);

}
