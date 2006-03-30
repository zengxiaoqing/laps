#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "geodefs.h"
#include "geoLib.h"
#include <math.h>

/*  List of Error return values possible:
     -5 : From setsup_gui.c	proj_idx not in range 1..10
     -6 : From setsup_gui.c	cornertype not 1,2,3
     -7 : From setsup_gui.c	stdlat and stdlat2 must have same sign
     -8 : From gen_map_bkgnd.c  program not called with correct number of args
     -9 : From create_bcd_bkgnd_xy.c  empty bcd filename  passed in
    -10 : From create_bcd_bkgnd_xy.c  Could not open vector_instructions.tk for writing
    -11 : From create_bcd_bkgnd_xy.c  Bad point count in bcd filename
    -12 : From create_bcd_bkgnd_xy.c  Could not read points in bcd filename
*/

int main(int argc, char *argv[])
{
       int proj_idx;
       float stdlat;
       float stdlon;
       float stdlat2;
       float ul_lat; 
       float ul_lon;
       float lr_lat;
       float lr_lon;
       char bcd_fname[256];
       char vector_fname[256];
       char setsup_fname[256];

       struct allgeoinfo OUTgeoinfo;
       struct allgeoinfo *geoinfoptr;
       float flt1, flt2;
       int cornertype;
       int istatus;

       
       if (argc == 12) {   /* got all of the arguments */

         cornertype = 2;
         
         proj_idx = atoi(argv[1]);
         stdlat = atof(argv[2]);
         stdlon = atof(argv[3]);
         stdlat2 = atof(argv[4]);
         ul_lat = atof(argv[5]);
         ul_lon = atof(argv[6]);
         lr_lat = atof(argv[7]);
         lr_lon = atof(argv[8]);
         strcpy(bcd_fname,argv[9]);
         strcpy(vector_fname,argv[10]);
         strcpy(setsup_fname,argv[11]);

/*
         printf("proj_idx = %d\n", proj_idx);
         printf("stdlat = %f\n", stdlat );
         printf("stdlon = %f\n", stdlon );
         printf("stdlat2 = %f\n", stdlat2 );
         printf("ul_lat = %f\n", ul_lat );
         printf("ul_lon = %f\n", ul_lon );
         printf("lr_lat = %f\n", lr_lat );
         printf("lr_lon = %f\n", lr_lon );
         printf("bcd_fname = %s\n", bcd_fname );
         printf("vector_fname = %s\n", vector_fname );
         printf("setsup_fname = %s\n", setsup_fname );
         printf("before create_bcd_bkgnd_xy\n");
         printf("after create_bcd_bkgnd_xy \n");
*/
         istatus = create_bcd_bkgnd_xy(&proj_idx,&stdlat,&stdlon,&stdlat2,&ul_lat,
                             &ul_lon,&lr_lat,&lr_lon,bcd_fname,
			     &cornertype,&vector_fname,&setsup_fname);



         if (istatus <= -5) {
           printf("%d:\n", istatus);
         } else {
           get_all_parms(&OUTgeoinfo);

           flt1 = (float) OUTgeoinfo.proj_idx;
           flt2 = (float) OUTgeoinfo.altproj;
           printf("xmin %f:xmax %f:ymin %f:ymax %f\n",
            OUTgeoinfo.xmin,
            OUTgeoinfo.xmax,
            OUTgeoinfo.ymin,
            OUTgeoinfo.ymax);
           printf("%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f:%f\n",
            flt1,
            OUTgeoinfo.cenlat,
            OUTgeoinfo.cenlon,
            OUTgeoinfo.rot,
            OUTgeoinfo.cone,
            OUTgeoinfo.polesgn,
            flt2,
            OUTgeoinfo.aaa,
            OUTgeoinfo.bbb,
            OUTgeoinfo.ccc,
            OUTgeoinfo.ddd,
            OUTgeoinfo.eee,
            OUTgeoinfo.fff,
            OUTgeoinfo.ggg,
            OUTgeoinfo.hhh,
            OUTgeoinfo.iii);
         }
         return istatus;

/*
*/
       }
       else {
/*       Program gen_map_bkgnd was not called with enough arguments.  */
         istatus = -8;
         printf("%d:\n", istatus);
         return istatus;
       }
}
