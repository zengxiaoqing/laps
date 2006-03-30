#include "geodefs.h"
#include "geoLib.h"
#include <stdio.h>
#include <math.h>

static struct allgeoinfo geoinfo;
static float view_hgt = 35839.0;

void setsup(int * proj_idx,
            float * cenlat,
            float * cenlon,
            float * rot,
            float * x1lat1,
            float * x2lon1,
            float * y1lat2,
            float * y2lon2,
            int * cornertype,
            int * success)
/*
C SetSup is used to define map projections.
C The origin and orientation of the projection are selected by the user.
C Points on the Earth defined by latitude and longitude are transformed to
C points in projection (x,y), the plane of projection.   

C   Usage
C           Call SetSup (proj_idx,cenlat,cenlon,rot,
C                    x1lat1,x2lon1,y1lat2,y2lon2,
C                    cornertype,success)

C   On input    proj_idx
C   for setsup    proj_idx defines the projection type
C             according to the following code:
C               1  Stereographic
C               2  Orthographic
C               3  Lambert Conformal Conic with two standard
C                  parallels
C               4  Lambert Equal Area
C               5  Gnomonic
C               6  Azimuthal Equidistant
C               7  Straight Great Circles
C               8  Cylindrical Equidistant
C               9  Mercator
C                  10  Mollweide type
C                  16  Radar Az/Ran

C           cenlat,PoLon,rot
C             If (proj_idx.ne.3)
C             . cenlat and cenlon define in degrees the latitude and
C               longitude of the point on the globe which is to
C               transform to the origin of the x,y plane.
C               -90 .le. cenlat .le. 90
C                  -180 .le. cenlon .le. 180
C               degrees of latitude north of the Equator and
C               degrees of longitude east of the Greenwich Meridian
C               are positive.  If the origin is at the North Pole,
C               "north" is considered to be in the direction of
C               (cenlon+180.). If the origin is at the South Pole,
C               "north" is in the direction of cenlon.
C             . rot is the angle between the y axis and north at
C               the origin.  It is measured in degrees and is taken
C               to be positive if the angular movement from north
C               to the y axis is counter-clockwise.  For the
C               cylindrical projections (8,9,10), the axis of the
C               projection is parallel to the y axis.
C             If (proj_idx.eq.3) (Lambert Conformal Conic)
C             . cenlon = central meridian of projection in degrees.
C             . cenlat,rot are the two standard parallels in deg.

C           cornertype, x1lat1,x2lon1,y1lat2,y2lon2
C             cornertype can take the values 1 through 5 and 
C             specifies one of five options on the way in which 
C             the limits of the rectangular map are defined by the
C             parameters x1lat1, x2lon1, y1lat2, and y2lon2.

C             cornertype = 1
C               The maximum useful area produced by the projection
C               is plotted.  x1lat1, x2lon1, y1lat2, and y2lon2
C               are not used and may be set to zero.

C             cornertype = 2
C               In this case (x1lat1,x2lon1) and (y1lat2,y2lon2) are the
C               latitudes and longitudes in degrees of two points
C               which are to be at opposite corners of the map
C               (upper right and lower left, respectively).
C               Care must be taken when using cylindrical
C               projections and this option.

C             cornertype = 3
C               The minimum and maximum values of u and v are
C               specified by x1lat1 through y2lon2.  x1lat1 = xMin,
C               x2lon1 = xMax, y1lat2 = yMin, y2lon2 = yMax.
C               Knowledge of the transformation equations is
C               necessary for this option to be used.
*/
{
   float coslat,coslon,cosrot,sinlat,sinlon,sinrot;
   float chi1,chi2;
   float x1,y1,x2,y2;
   GeoInfo * geoinfoptr;
  
   *success = 0;
   if (*proj_idx<1 || *proj_idx>17 ||
       *proj_idx>10 && *proj_idx<16) return;
   if (*cornertype<1 || *cornertype>3) return;

   geoinfo.proj_idx = *proj_idx;
   geoinfo.cenlat = *cenlat;
   geoinfo.cenlon = *cenlon;
   geoinfo.rot = *rot;
   geoinfo.altproj = *proj_idx;

   if (*proj_idx==LAMBERT_CONFORMAL) {
       geoinfo.aaa = 0.0;
       if (geoinfo.rot+geoinfo.cenlat>0) {
           geoinfo.polesgn = 1.0;
       } else {
           geoinfo.polesgn = 1.0;
       }/*endif*/

       if (geoinfo.rot==geoinfo.cenlat) {
           chi1 = (90.0 - geoinfo.polesgn*geoinfo.cenlat)*DTR;
           geoinfo.cone = cos(chi1);
       } else {
           chi1 = (90.0 - geoinfo.polesgn*geoinfo.cenlat)*DTR;
           chi2 = (90.0 - geoinfo.polesgn*geoinfo.rot)*DTR;
           geoinfo.cone = log(sin(chi1)/sin(chi2))/
                          log(tan(chi1/2)/tan(chi2/2));
       }/*endif*/

       printf("cone = %f\n", geoinfo.cone);
       if (geoinfo.polesgn<0)
            geoinfo.polesgn = 1.0;


   } else {
       coslat = cos(*cenlat*DTR);
       coslon = cos(*cenlon*DTR);
       cosrot = cos(*rot*DTR);
       sinlat = sin(*cenlat*DTR);
       sinlon = sin(*cenlon*DTR);
       sinrot = sin(*rot*DTR);
       geoinfo.aaa = cosrot*coslon*sinlat-sinrot*sinlon;
       geoinfo.bbb = cosrot*sinlon*sinlat+sinrot*coslon;
       geoinfo.ccc = -cosrot*coslat;
       geoinfo.ddd = -sinrot*coslon*sinlat-cosrot*sinlon;
       geoinfo.eee = -sinrot*sinlon*sinlat+cosrot*coslon;
       geoinfo.fff = sinrot*coslat;
       geoinfo.ggg = coslon*coslat;
       geoinfo.hhh = sinlon*coslat;
       geoinfo.iii = sinlat;
       if (geoinfo.proj_idx==SATELLITE_VIEW)
           geoinfo.cone = (view_hgt+R_EARTH)/R_EARTH;
   }/*endif*/

   switch (geoinfo.proj_idx) {
      case STEREOGRAPHIC:
      case AZIMUTHAL_EQUAL_AREA:
      case GNOMONIC:
        geoinfo.xmin = geoinfo.ymin = -2.0;
        geoinfo.xmax = geoinfo.ymax = 2.0;
        break;
      case ORTHOGRAPHIC:
      case LAMBERT_CONFORMAL:
        geoinfo.xmin = geoinfo.ymin = -1.0;
        geoinfo.xmax = geoinfo.ymax = 1.0;
        break;
      case AZIMUTHAL_EQUIDISTANT:
      case MERCATOR:
        geoinfo.xmin = geoinfo.ymin = -PI;
        geoinfo.xmax = geoinfo.ymax = PI;
        break;
      case CYLINDRICAL_EQUIDISTANT:
        geoinfo.xmin = -180.0;
        geoinfo.ymin = -90.0;
        geoinfo.xmax = 180.0;
        geoinfo.ymax = 90.0;
        break;
      case MOLLWEIDE:
        geoinfo.xmin = -2.0;
        geoinfo.ymin = -1.0;
        geoinfo.xmax = 2.0;
        geoinfo.ymax = 1.0;
        break;
      case AZIMUTH_RANGE:
      case RADAR_AZ_RAN:
        geoinfo.xmin = geoinfo.ymin = 0;
        geoinfo.xmax = 600.0;
        geoinfo.ymax = 360.0;
        break;
      case SATELLITE_VIEW:
        geoinfo.xmax = asin(1/geoinfo.cone);
        if (geoinfo.xmax>PI/4) geoinfo.xmax = PI/4;
        geoinfo.ymax = geoinfo.xmax;
        geoinfo.ymin = geoinfo.xmin = -geoinfo.xmax;
        break;
   }/*end switch*/

   if (*cornertype==LAT_LON_CORNERS_AREA) {
       geoinfoptr = (GeoInfo *) &geoinfo;
       ll_to_xy(x1lat1,x2lon1,&x1,&y1,geoinfoptr);

       /*Sample all corners to get max size*/
       ll_to_xy(y1lat2,y2lon2,&x2,&y2,geoinfoptr);
       geoinfo.xmin = (x1<x2 ? x1 : x2);
       geoinfo.ymin = (y1<y2 ? y1 : y2);
       geoinfo.xmax = (x1>x2 ? x1 : x2);
       geoinfo.ymax = (y1>y2 ? y1 : y2);

       ll_to_xy(x1lat1,y2lon2,&x1,&y2,geoinfoptr);
       geoinfo.xmin = (x1<x2 ? x1 : x2);
       geoinfo.ymin = (y1<y2 ? y1 : y2);
       geoinfo.xmax = (x1>x2 ? x1 : x2);
       geoinfo.ymax = (y1>y2 ? y1 : y2);

       ll_to_xy(y1lat2,x2lon1,&x2,&y1,geoinfoptr);
       geoinfo.xmin = (x1<x2 ? x1 : x2);
       geoinfo.ymin = (y1<y2 ? y1 : y2);
       geoinfo.xmax = (x1>x2 ? x1 : x2);
       geoinfo.ymax = (y1>y2 ? y1 : y2);
   } else if (*cornertype==CARTESIAN_CORNERS_AREA) {
       geoinfo.xmin = *x1lat1;
       geoinfo.ymin = *y1lat2;
       geoinfo.xmax = *x2lon1;
       geoinfo.ymax = *y2lon2;
   }/*endif*/

/*
       printf("Inside setsup\n");
       printf("xmin = %f\n", geoinfo.xmin);
       printf("xmax = %f\n", geoinfo.xmax);
       printf("ymin = %f\n", geoinfo.ymin);
       printf("ymax = %f\n", geoinfo.ymax);

       printf("proj_idx = %d\n", geoinfo.proj_idx);
       printf("cenlat = %f\n", geoinfo.cenlat);
       printf("cenlon = %f\n", geoinfo.cenlon);
       printf("rot = %f\n", geoinfo.rot);
       printf("cone = %f\n", geoinfo.cone);
       printf("polesgn = %f\n", geoinfo.polesgn);
       printf("altproj = %d\n", geoinfo.altproj);
       printf("aaa = %f\n", geoinfo.aaa);
       printf("bbb = %f\n", geoinfo.bbb);
       printf("ccc = %f\n", geoinfo.ccc);
       printf("ddd = %f\n", geoinfo.ddd);
       printf("eee = %f\n", geoinfo.eee);
       printf("fff = %f\n", geoinfo.fff);
       printf("ggg = %f\n", geoinfo.ggg);
       printf("hhh = %f\n", geoinfo.hhh);
       printf("iii = %f\n\n\n", geoinfo.iii);
*/

/* if (geoinfo.xmin<geoinfo.xmax ||
       geoinfo.ymin<geoinfo.ymax)  */ *success = 1;

}/* end setsup */

void get_proj_parm(GeoInfo * projparm)
{
   struct allgeoinfo * allgeoinfoptr;

   allgeoinfoptr = (struct allgeoinfo *) projparm;
   *allgeoinfoptr = geoinfo;
}/* end get_proj_parm */

void get_all_parms(struct allgeoinfo *allparms)
{
   struct allgeoinfo *allgeoinfoptr;

   allgeoinfoptr = (struct allgeoinfo *) allparms;
   *allgeoinfoptr = geoinfo;
}/* end get_all_parms */

void save_geo_file(char * filename, int * success)
{
   FILE *fp;

   *success = 0;
   fp = fopen(filename, "w");
   if (fp==NULL) return;

   *success = fprintf(fp, "%f\n", geoinfo.cenlon);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.cone);
   if (*success) *success = fprintf(fp, "901\n");
   if (*success) *success = fprintf(fp, "902\n");
   if (*success) *success = fprintf(fp, "903\n");
   if (*success) *success = fprintf(fp, "904\n");
   if (*success) *success = fprintf(fp, "%f\n", geoinfo.polesgn);
   if (*success) *success = fprintf(fp, "905\n");
   if (*success) *success = fprintf(fp, "906\n");
   if (*success) *success = fprintf(fp, "907\n");
   if (*success) *success = fprintf(fp, "908\n");
   if (*success) *success = fprintf(fp, "909\n");
   if (*success) *success = fprintf(fp, "%d\n", geoinfo.proj_idx);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.xmin);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.xmax);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.ymin);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.ymax);
   if (*success) *success = fprintf(fp, "%d\n", geoinfo.altproj);
   if (*success) *success = fprintf(fp, "%f\n", geoinfo.rot);
   if (*success) *success = fprintf(fp, "910\n");
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.cenlat);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.aaa);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.bbb);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.ccc);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.ddd);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.eee);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.fff);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.ggg);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.hhh);
   if (*success) *success = fprintf(fp, "%.8f\n", geoinfo.iii);
   fclose(fp);
} /* end save_geo_file */

void set_view_hgt(float * hgt)
{
   view_hgt = *hgt;
} /* end set_view_hgt */

void save_geo_file2(char * filename)
/* The output order in save_geo_file2 differs from save_geo_file. */
{
   FILE *fp;
   int success;

   success = 0;
   fp = fopen(filename, "w");
   if (fp==NULL) return;

   success = 1;
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.xmin);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.xmax);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.ymin);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.ymax);
   if (success) success = fprintf(fp, "%d\n", geoinfo.proj_idx);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.cenlat);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.cenlon);
   if (success) success = fprintf(fp, "%f\n", geoinfo.rot);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.cone);
   if (success) success = fprintf(fp, "%f\n", geoinfo.polesgn);
   if (success) success = fprintf(fp, "%d\n", geoinfo.altproj);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.aaa);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.bbb);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.ccc);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.ddd);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.eee);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.fff);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.ggg);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.hhh);
   if (success) success = fprintf(fp, "%.8f\n", geoinfo.iii);

   fclose(fp);

} /* end save_geo_file2 */


void get_geo_file(char * filename, int * success)
{
   FILE *fp;

   *success = 0;
   fp = fopen(filename, "r");
   if (fp==NULL) return;

   fscanf(fp, "%f", &geoinfo.xmin);
   fscanf(fp, "%f", &geoinfo.xmax);
   fscanf(fp, "%f", &geoinfo.ymin);
   fscanf(fp, "%f", &geoinfo.ymax);
   fscanf(fp, "%d", &geoinfo.proj_idx);
   fscanf(fp, "%f", &geoinfo.cenlat);
   fscanf(fp, "%f", &geoinfo.cenlon);
   fscanf(fp, "%f", &geoinfo.rot);
   fscanf(fp, "%f", &geoinfo.cone);
   fscanf(fp, "%f", &geoinfo.polesgn);
   fscanf(fp, "%d", &geoinfo.altproj);
   fscanf(fp, "%f", &geoinfo.aaa);
   fscanf(fp, "%f", &geoinfo.bbb);
   fscanf(fp, "%f", &geoinfo.ccc);
   fscanf(fp, "%f", &geoinfo.ddd);
   fscanf(fp, "%f", &geoinfo.eee);
   fscanf(fp, "%f", &geoinfo.fff);
   fscanf(fp, "%f", &geoinfo.ggg);
   fscanf(fp, "%f", &geoinfo.hhh);
   fscanf(fp, "%f", &geoinfo.iii);
   fclose(fp);

} /* end get_geo_file */
