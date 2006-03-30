#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "geodefs.h"
#include "geoLib.h"
#include <math.h>

/*  List of Error return values possible:
     -8: From pwrap_ll_to_xy.c  program not called with correct number of args
*/
 
int main(int argc, char *argv[])
{
       float lat, lon;
       float x, y;
       int type;
       int istatus;
       int success;
       char setsup_fname[256];
    /* char fname[20]="/tmp/setsup.dat"; */
       
       struct allgeoinfo geoinfo;
       struct allgeoinfo *geoinfoptr;
       struct proj_parm my_proj;

       success=0;

       if (argc == 5) {   /* Got all of the arguments.*/

         strcpy(setsup_fname,argv[4]);
         type = atoi(argv[1]);

	 /* Fill geo structure, from text file */
         get_geo_file(setsup_fname, &success);

	 /* Fill proj structure, from geoinfo */
         get_proj_parm(&my_proj);

	 if (type == 1) {
	   /* ll to xy */
           lat = atof(argv[2]);
           lon = atof(argv[3]);
           ll_to_xy(&lat,&lon, &x,&y, &my_proj);
           printf("%f,%f", x,y);

	 } else {
	   /* xy to ll */
           x = atof(argv[2]);
           y = atof(argv[3]);
           xy_to_ll(&x,&y, &lat,&lon, &my_proj);
           printf("%f,%f", lat,lon);
         }
	 return;

       } else {
        /* Program pwrap_ll_to_xy was not called with enough arguments. */
         istatus = -8;
         printf("%d:\n", istatus);
         return istatus;
       }
}
