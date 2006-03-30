/* Location_lines.c generates latitude and longitude lines to 
 * give reference coordinates to political and geographical 
 * maps (esp, those generated with create_bcd_background_xy.c).
 *
 * Author: Phil McDonald & Paula McCaslin   18 April 2002
 * 
 */

#include "geodefs.h"
#include "geoLib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define max(a,b) (a>b)?a:b
#define min(a,b) (a<b)?a:b


int location_lines(FILE * vfp, 
        float *lat_ur,
        float *lon_ur,
        float *lat_ll,
        float *lon_ll,
        float *xMin,
        float *xMax,
        float *yMin,
        float *yMax,
        float *scaleby_x,
        float *scaleby_y,
        float *x_frac,
        float *y_frac,
        float *line_color, 
        GeoInfo * my_proj,
        int lambrt)
{

float lon_west, lon_east, lat_north, lat_south; 
float lon_diff, lat_diff, lon_mid, lat_mid;
float ll_array[20000][4];
float text_array[100][3];
float i, j, x, y, x2, y2;
float offset_x, offset_y, offset_x2, offset_y2;
float x_max, x_min, y_max, y_min;
float lon_step, half_lon, half_lat;
int   k, indx, t_indx, iMax;
char  fancy_font[20] = "helvetica 11 bold";
float lat_max, lat_min, lon_max, lon_min;
float lab;


if (*lon_ll > *lon_ur) {
    printf ( "Possible ERROR in location_lines.c: left longitude is > right \
              longitude: %f > %f\n", *lon_ll, *lon_ur);
   if (*lon_ll == 180.) { *lon_ll = -180.; }
}

if (*lat_ll > *lat_ur) {
    printf ( "Possible ERROR in location_lines.c : south latitude is > north \
              latitude: %f > %f \n", *lat_ll, *lat_ur);
}


/* Break latitude lines into small segments resulting 
 * in smooth and curved lines of latitude. 
 */
lon_step=5; 

/* Display limited number of lat,lon lines (e.g. lat_diff) 
 * based on size of domain. 
 */
lon_diff=abs(*lon_ur - *lon_ll); 
lat_diff=abs(*lat_ur - *lat_ll); 
lon_mid=(*lon_ur + *lon_ll)/2; 
lat_mid=(*lat_ur + *lat_ll)/2; 

if (1) {
if (lon_diff > 51) {
    lon_diff=20;
} else if (lon_diff > 41) {
    lon_diff=15;
} else if (lon_diff > 31) {
    lon_diff=10;
} else if (lon_diff > 11) {
    lon_diff=5;
} else if (lon_diff > 3) {
    lon_diff=2.5;
} else {
    lon_diff=1;
}

if (lat_diff > 50) {
    lat_diff=20;
} else if (lat_diff > 21) {
    lat_diff=10;
} else if (lat_diff > 11) {
    lat_diff=5;
} else if (lat_diff > 3) {
    lat_diff=2.5;
} else {
    lat_diff=1;
}
} else {

lon_diff=10;
lat_diff=10;

if (lon_diff > 21) {
    lon_diff=10;
} else if (lon_diff > 11) {
    lon_diff=5;
} else if (lon_diff > 3) {
    lon_diff=2.5;
} else {
    lon_diff=1;
}

if (lat_diff > 21) {
    lat_diff=10;
} else if (lat_diff > 11) {
    lat_diff=5;
} else if (lat_diff > 3) {
    lat_diff=2.5;
} else {
    lat_diff=1;
}


}



/* With non-Mercator projections data in the lower right corner
 * is omitted, get_inclusive_extents corrects this. */
/* if (my_proj->proj_idx != 19) { */

   get_inclusive_extents_xy(vfp, my_proj,
                      xMax,   xMin,   yMax,   yMin,
                      lat_ur, lon_ur, lat_ll, lon_ll); 

x_max=*xMax;
x_min=*xMin;
y_max=*yMax;
y_min=*yMin;

lat_max=*lat_ur;
lon_max=*lon_ur;
lat_min=*lat_ll;
lon_min=*lon_ll;

/* Create intervals for grid lines and grid labels that 
 * are astetic by finding min and max lat,lon 
 * e.g. values ending in 0, 5 and store them in the vars
 * lon_west, lon_ease, lat_north, lat_south.  */
lon_west=-280;
lon_east= 280;
lat_north= 100;
lat_south=-100;

k=280;
for (i=-k; i <= k; i+=lon_diff) {
    if(lon_min <= i){
       lon_west=i;
       /*if(lon_west <= -180) {lon_west=-180; break;}*/
       for (i=lon_west; i <= k; i+=lon_diff) {
           if(lon_max <= i){
              lon_east=i-=lon_diff;
              break;
           }
       }
       break;
    }
}

for (i= 100; i >= -100; i-=lat_diff) {
    if(lat_max >= i){
       lat_north=i;
       for (i=lat_north; i >= -100; i-=lat_diff) {
           if(lat_min >= i){
              lat_south=i+=lat_diff;
              break;
           }
       }
       break;
    }
}

/* Find the center of the lat,lon lines in order to affix labels. */

x=(x_min+x_max+.005)/2;
y=(y_min+y_max+.005)/2;
xy_to_ll(&x_min,&y_min, &half_lat,&half_lon, my_proj);

/* If domain crosses the Equator then put labels on the 0 lat line. */

if (*lat_ll < 0 && *lat_ur > 0) { 
	half_lat=0; 

} else if (my_proj->proj_idx == 1) { 
/* Deal with polar labels that cross polar caps. */
   if ( abs(half_lat - *lat_ll) > abs(*lat_ur - *lat_ll) ) { 
       half_lat=((90 - *lat_ll)/2) + *lat_ll; 
   }
} else {
   half_lat=my_proj->cenlat;
}
half_lon=my_proj->cenlon;

indx=0; 
t_indx=0;

/* Perspective boundary around rectangular boundary. Not necessary. */
if (1) {
        /*
         * West edge.
         */
        ll_array[indx][0]=lat_max;
        ll_array[indx][1]=lon_min;
        ll_array[indx][2]=lat_min;
        ll_array[indx][3]=lon_min;
        indx++; 
        
        /*
         * East edge.
         */
        ll_array[indx][0]=lat_max;
        ll_array[indx][1]=lon_max;
        ll_array[indx][2]=lat_min;
        ll_array[indx][3]=lon_max;
        indx++; 
        
        /*
         * South edge.
         */
        for (j=lon_min; j <= lon_max; j+=lon_step) {
           ll_array[indx][0]=lat_min;
           ll_array[indx][1]=j;
           ll_array[indx][2]=lat_min;
           ll_array[indx][3]=j+lon_step;
           indx++; 
           }
           
        indx--; 
        ll_array[indx][3]=lon_max;
        indx++; 

        /*
         * North edge.
         */
        for (j=lon_min; j <= lon_max; j+=lon_step) {

           ll_array[indx][0]=lat_max;
           ll_array[indx][1]=j;
           ll_array[indx][2]=lat_max;
           ll_array[indx][3]=j+lon_step;
           indx++; 
        }
        indx--; 
        ll_array[indx][3]=lon_max;
        indx++; 

}        
iMax=indx;


/*
 * Create longitudinal reference lines (dashed) and labels.
 */

for (i=lon_west; i <= lon_east; i+=lon_diff) {

     ll_array[indx][0]=lat_max;
     ll_array[indx][1]=i;
     ll_array[indx][2]=lat_min;
     ll_array[indx][3]=i;
     indx++; 

     text_array[t_indx][0]=half_lat;
     text_array[t_indx][1]=i;
     if (i < -180) { 
        lab=i+360;
     } else if (i > 180) { 
        lab=i-360;
     } else {
        lab=i;
     }
     text_array[t_indx][2]=lab;
     t_indx++;

}

/*
 * Create latitudinal reference lines (dashed) and labels.
 */

for (i=lat_north; i >= lat_south; i-=lat_diff) {

     text_array[t_indx][0]=i;
     text_array[t_indx][1]=half_lon;
     text_array[t_indx][2]=i;
     t_indx++;

     for (j=lon_min; j <= lon_max; j+=lon_step) {

          ll_array[indx][0]=i;
          ll_array[indx][1]=j;
          ll_array[indx][2]=i;
          ll_array[indx][3]=j+lon_step;
          indx++; 
     }
     indx--; 
     ll_array[indx][3]=lon_max;
     indx++; 
}

/*
 * Write output location line vector instructions for tk.
 */


fprintf(vfp,"sub ll_instructions {\n");
fprintf(vfp,"my($can,$cx,$cy) = @_;\n");
fprintf(vfp,"");
 

/*
 * Write output lat,lon lines.
 */
/* Find Polar caps.*/
if (my_proj->proj_idx == 1 && lon_diff > 41) { };

if (lon_diff < 11 && lat_diff < 11) {
for (k=0; k < iMax; k++) {
    ll_to_xy(&(ll_array[k][0]),&(ll_array[k][1]), &x, &y,  my_proj);
    ll_to_xy(&(ll_array[k][2]),&(ll_array[k][3]), &x2,&y2, my_proj);

    offset_x =( (x  - x_min) * *scaleby_x );
    offset_x2=( (x2 - x_min) * *scaleby_x );
    offset_y =( (y_max - y ) * *scaleby_y );
    offset_y2=( (y_max - y2) * *scaleby_y );

    if (lambrt) {
       offset_x =*x_frac - offset_x;
       offset_x2=*x_frac - offset_x2;
       offset_y =*y_frac - offset_y;
       offset_y2=*y_frac - offset_y2;
    }

    fprintf (vfp, "$can->createLine($cx * %f, $cy * %f, $cx * %f, $cy * %f, \
                   -fill => 'white', -dash => '.  ', \
                   -width => 1, -tags => 'location_lines');\n", \
                   offset_x, offset_y, offset_x2, offset_y2);

}

for (k=iMax; k < indx; k++) {
    ll_to_xy(&(ll_array[k][0]),&(ll_array[k][1]), &x, &y,  my_proj);
    ll_to_xy(&(ll_array[k][2]),&(ll_array[k][3]), &x2,&y2, my_proj);

    offset_x =( (x  - x_min) * *scaleby_x );
    offset_x2=( (x2 - x_min) * *scaleby_x );
    offset_y =( (y_max - y ) * *scaleby_y );
    offset_y2=( (y_max - y2) * *scaleby_y );

    if (lambrt) {
       offset_x =*x_frac - offset_x;
       offset_x2=*x_frac - offset_x2;
       offset_y =*y_frac - offset_y;
       offset_y2=*y_frac - offset_y2;
    }

    fprintf (vfp, "$can->createLine($cx * %f, $cy * %f, $cx * %f, $cy * %f, \
                   -fill => 'white', -dash => '.  ', \
                   -width => 1, -tags => 'location_lines');\n", \
                   offset_x, offset_y, offset_x2, offset_y2);
}
}
/*
 * Write labels.
 */

for (k=0; k < t_indx; k++) {
    ll_to_xy(&(text_array[k][0]),&(text_array[k][1]), &x,&y, my_proj);

    offset_x =( (x - x_min)  * *scaleby_x );
    offset_y =( (y_max - y)  * *scaleby_y );
    if (lambrt) {
       offset_x=*x_frac - offset_x;
       offset_y=*y_frac - offset_y;
    }

    fprintf (vfp, "$can->createText($cx * %f, $cy * %f, \
                   -fill => 'white', -text => %f, \
                   -font => '%s', -tags => 'location_lines');\n", \
                   offset_x, offset_y, text_array[k][2],fancy_font);
}

}


/*
 * ----- The following subroutine allows for the -----
 * ----- inclusion of data in the lower right --------
 * ----- corner, data that otherwise would be --------
 * ----- omitted with non-Mercator projections. ------
 */

int get_inclusive_extents_ll(
        FILE * vfp, 
        GeoInfo * my_proj1,
        float *lat_ur,
        float *lon_ur,
        float *lat_ll,
        float *lon_ll,
        float *x_max, 
        float *x_min, 
        float *y_max, 
        float *y_min)
{
int k;
float xy_pair[8][2];
float x, y;
float x_mid, y_mid;
float x_ur, y_ur, x_ll, y_ll;

/*
 * With the min and max lat,lon values 
 *
 *  ---ur
 * |    |
 * |    |
 * ll--- 
 *
 * Calculate x,y values from the two lat,lon 
 * pairs above.
 */

ll_to_xy(lat_ur,lon_ur, &x_ur,&y_ur, my_proj1);
ll_to_xy(lat_ll,lon_ll, &x_ll,&y_ll, my_proj1);

xy_pair[0][0]= x_ll;  /* lower left */
xy_pair[0][1]= y_ll;
xy_pair[1][0]= x_ur;  /* upper right */ 
xy_pair[1][1]= y_ur;

/*
x_mid=(x_ur+x_ll)/2;
y_mid=(y_ur+y_ll)/2;
xy_pair[2][0]= x_mid;  midpoints 
xy_pair[2][1]= y_mid;
*/

/*
 * Find the min and max x,y area.
 */

*x_min=  1.;
*x_max= -1.;
*y_min=  1.;
*y_max= -2.;

for (k=0; k < 2; k++) {

     x=xy_pair[k][0];
     y=xy_pair[k][1];
     *y_max = max(y,*y_max);
     *y_min = min(y,*y_min);
     *x_max = max(x,*x_max);
     *x_min = min(x,*x_min);
}

return;
}

int get_inclusive_extents_xy(
        FILE * vfp, 
        GeoInfo * my_proj1,
        float *x_max, 
        float *x_min, 
        float *y_max, 
        float *y_min,
        float *lat_ur,
        float *lon_ur,
        float *lat_ll,
        float *lon_ll)
{
int k;
float lat, lon;
float xy_pair[8][2];
float x_mid, y_mid;
float nlat_ur, nlon_ur, nlat_ll, nlon_ll;
int flag;

x_mid=(*x_max+*x_min)/ 2;
y_mid=(*y_max+*y_min)/ 2;

/*
 * With the min and max x, y values from above create eight 
 * points on a rectangle (four corners and four mid points) 
 * by filling x,y values as follows:
 *
 * ul---upper---ur
 * |             |
 * left      right
 * |             |
 * ll---lower---lr
 */

xy_pair[0][0]=*x_min; /* ul corner */
xy_pair[0][1]=*y_max;
xy_pair[1][0]=*x_min; /* left midpoint */
xy_pair[1][1]= y_mid;
xy_pair[2][0]=*x_min; /* ll corner */
xy_pair[2][1]=*y_min;
xy_pair[3][0]= x_mid; /* upper midpoint */
xy_pair[3][1]=*y_max;
xy_pair[4][0]= x_mid; /* lower midpoint */
xy_pair[4][1]=*y_min;
xy_pair[5][0]=*x_max; /* ur corner */
xy_pair[5][1]=*y_max;
xy_pair[6][0]=*x_max; /* right midpoint */
xy_pair[6][1]= y_mid;
xy_pair[7][0]=*x_max; /* lr corner */
xy_pair[7][1]=*y_min;

/*
 * Calculate lat,lon values from the eight x,y 
 * values above, then find the min and max lat,lon 
 * values [xy_to_ll].
 */

nlat_ur= -90.;
nlon_ur=-180.;
nlat_ll=  90.;
nlon_ll= 180.;

for (k=0; k < 8; k++) {
     /* Convert xy to lat,lon pairs */
     xy_to_ll(&(xy_pair[k][0]),&(xy_pair[k][1]), &lat,&lon, my_proj1);
          
     if (*lon_ll > lon) { 
         flag=1;
     } else {
         flag=0;
     }

     if (lon > *lon_ll+180 && !flag) {
        lon=lon-360;
     } else if (lon < *lon_ll-180 && flag) {
        lon=lon+360;
     }

     /* Find min and max lat,lon values */
     nlat_ur = max(lat,nlat_ur);
     nlat_ll = min(lat,nlat_ll);
     if (k != 3) {
       nlon_ur = max(lon,nlon_ur);
       nlon_ll = min(lon,nlon_ll);
     }
}

*lat_ur=nlat_ur;
*lat_ll=nlat_ll;
*lon_ur=nlon_ur;
*lon_ll=nlon_ll;
return;

}
