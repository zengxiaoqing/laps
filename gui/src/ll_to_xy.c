#include "geoLib.h"
#include "geodefs.h"
#include <math.h>

void ll_to_xy(float * lat, float * lon,
              float * x, float * y,
              GeoInfo * ingeoinfo)
{
   float xx,yy,zz;
   float xl,yl,zl;
   float latrad,lonrad,coslat;
   float t,r;

   struct allgeoinfo * geoinfo;

   geoinfo = (struct allgeoinfo *) ingeoinfo;

   if (geoinfo->proj_idx==CYLINDRICAL_EQUIDISTANT &&
       geoinfo->cenlat==0 && geoinfo->cenlon==0 && geoinfo->rot==0) {
       *x = *lon;
       *y = *lat;
       return;
   }/*endif*/

   *x = 1e12;
   *y = 1e12;
   if (*lat>1e11) return;

   if (geoinfo->proj_idx!=LAMBERT_CONFORMAL) {

/* the 3d vector (xl,yl,zl) is a unit vector on the earth's surface such that
   (0,0,1) is at the north pole, (1,0,0) is at the intersection of the equator
   and the Greenwich meridian and (0,1,0) is at 0N, 90E */

       latrad = *lat*DTR;
       lonrad = *lon*DTR;
       coslat = cos(latrad);
       xl = cos(lonrad)*coslat;
       yl = sin(lonrad)*coslat;
       zl = sin(latrad);

/* the 3d vector (xx,yy,zz) is a unit vector on the earth's surface such that
   (0,0,1) is at the projection center, (0,1,0) is in the positive x direction
   and (-1,0,0) is in the positive y direction */

       xx = geoinfo->aaa*xl+geoinfo->bbb*yl+geoinfo->ccc*zl;
       yy = geoinfo->ddd*xl+geoinfo->eee*yl+geoinfo->fff*zl;
       zz = geoinfo->ggg*xl+geoinfo->hhh*yl+geoinfo->iii*zl;

   }/*endif*/

   switch (geoinfo->proj_idx) {
      case STEREOGRAPHIC:
        if (zz<=-1.0) return;
        t = 1/(1+zz);
        *x = t*yy;
        *y = -t*xx;
        return;
      case ORTHOGRAPHIC:
        if (zz<0.0) return;
        *x = yy;
        *y = -xx;
        return;
      case LAMBERT_CONFORMAL:
        t = (90.0-geoinfo->polesgn*(*lat))*DTR;
        if (t>=PI) return;
        r = pow(tan(t/2),geoinfo->cone);
        t = (*lon)-geoinfo->cenlon;
        if (t<-180.0) 
            t += 360.0;
        else if (t>180.0) 
            t -= 360.0;
        t = DTR*t*geoinfo->cone;
        *x = sin(t)*r;
        *y = -cos(t)*r;
        return;
      case AZIMUTHAL_EQUAL_AREA:
        t = 1+zz;
        if (t<=0.0) return;
        t = sqrt(2/t);
        *x = t*yy;
        *y = -t*xx;
        return;
      case GNOMONIC:
        if (zz<=0.0) return;
        *x = yy/zz;
        *y = -xx/zz;
        return;
      case AZIMUTHAL_EQUIDISTANT:
        t = sqrt(xx*xx+yy*yy);
        if (t==0.0) {
            if (zz<0) return;
            *x = 0.0;
            *y = 0.0;
            return;
        }/*endif*/
        t = atan2(t,zz)/t;
        *x = t*yy;
        *y = -t*xx;
        return;
      case SATELLITE_VIEW:
        if (zz<1/geoinfo->cone) return;
        t = geoinfo->cone-zz;
        *x = atan2(yy,t);
        *y = -atan(xx/sqrt(yy*yy+t*t));
        return;
      case CYLINDRICAL_EQUIDISTANT:
        if (yy==0.0 && zz==0.0) return;
        *x = RTD*atan2(yy,zz);
        *y = -RTD*asin(xx);
        return;
      case MERCATOR:
        if (yy==0.0 && zz==0.0) return;
        if (xx<=-1.0 || xx>=1.0) return;
        *x = atan2(yy,zz);
        *y = -log(tan(asin(xx)/2+PI/4));
        return;
      case MOLLWEIDE:
        t = sqrt(zz*zz+yy*yy);
        if (t<=0.0) return;
        *x = atan2(yy,zz)*t*2/PI;
        *y = -xx;
        return;
      case AZIMUTH_RANGE:
        if (xx==0.0 && yy==0.0)
            *y = 0.0;
        else {
            *y = atan2(yy,-xx)*RTD;
            if (*y<0.0) *y += 360.0;
        }/*endif*/
        *x = atan2(sqrt(xx*xx+yy*yy),zz)*R_EARTH;
        return;
      case RADAR_AZ_RAN:
        if (xx==0.0 && yy==0.0)
            *y = 0.0;
        else {
            *y = atan2(yy,-xx)*RTD;
            if (*y<0.0) *y += 360.0;
        }/*endif*/
        *x = sqrt(xx*xx+yy*yy)*R_EARTH;
        return;
      default:
        return;
   }/*end switch*/

}
