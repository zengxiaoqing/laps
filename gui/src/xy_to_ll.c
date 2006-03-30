#include "geoLib.h"
#include "geodefs.h"
#include <math.h>

void xy_to_ll(float * x, float * y,
              float * lat, float * lon,
              GeoInfo * ingeoinfo)
{
   float xx,yy,zz;
   float xl,yl,zl;
   float t,r,coneinv;

   struct allgeoinfo * geoinfo;
   
   geoinfo = (struct allgeoinfo *) ingeoinfo;
   *lat = 1e12;
   *lon = 1e12;
   if (*x>1e11) return;

/* the 3d vector (xx,yy,zz) is a unit vector on the earth's surface such that
   (0,0,1) is at the projection center, (0,1,0) is in the positive x direction
   and (-1,0,0) is in the positive y direction */

   switch (geoinfo->proj_idx) {
      case STEREOGRAPHIC:
        t = 2/(1+ (*x)*(*x) + (*y)*(*y) );
        xx = -t*(*y);
        yy = t*(*x);
        zz = t-1;
        break;
      case ORTHOGRAPHIC:
        xx = -(*y);
        yy = (*x);
        t = 1-xx*xx-yy*yy;
        if (t<=0.0) return;
        zz = sqrt(t);
        break;
      case LAMBERT_CONFORMAL:
        if ((*x)==0 && (*y)==0) {
            *lat = 90.0*geoinfo->polesgn;
            *lon = geoinfo->cenlon;
            return;
        }/*endif*/
        r = sqrt((*x)*(*x)+(*y)*(*y));
        coneinv = 1/geoinfo->cone;
        t = coneinv*atan2((*x),-geoinfo->polesgn*(*y))*RTD;
        if (t<-180.0 || t>180.0) return;
        *lon = geoinfo->cenlon+t;
        if (*lon<-180)
            *lon += 360;
        else if (*lon>180)
            *lon -= 360;
        t = 2*atan(pow(r,coneinv));
        *lat = geoinfo->polesgn*(90-RTD*t);
        return;
      case AZIMUTHAL_EQUAL_AREA:
        r = (*x)*(*x)+(*y)*(*y);
        t = 1-r/4;
        if (t<=0) return;
        t = sqrt(t);
        xx = -t*(*y);
        yy = t*(*x);
        zz = 1-r/2;
        break;
      case GNOMONIC:
        zz = 1/sqrt(1+(*x)*(*x)+(*y)*(*y));
        xx = -zz*(*y);
        yy = zz*(*x);
        break;
      case AZIMUTHAL_EQUIDISTANT:
        r = sqrt((*x)*(*x)+(*y)*(*y));
        if (r==0) {
            xx = 0;
            yy = 0;
            zz = 1;
            break;
        } else if (r>PI)
            return;
        t = sin(r)/r;
        xx = -t*(*y);
        yy = t*(*x);
        zz = cos(r);
       break;
      case SATELLITE_VIEW:
        t = cos(*y);
        xx = -sin(*y);
        yy = t*sin(*x);
        zz = t*cos(*x);
        t = (geoinfo->cone*geoinfo->cone)*(zz*zz-1)+1;
        if (t<0) return;
        t = geoinfo->cone*zz-sqrt(t);
        xx *= t;
        yy *= t;
        zz = geoinfo->cone-zz*t;
       break;
      case CYLINDRICAL_EQUIDISTANT:
        if (geoinfo->cenlat==0 && geoinfo->cenlon==0 && geoinfo->rot==0) {
            *lat = *y;
            *lon = *x;
            return;
        }/*endif*/
        if (*x<-180.0 || *x>180.0 || *y<-90.0 || *y>90.0) return;
        r = (*y)*DTR;
        xx = sin(-r);
        t = cos(r);
        r = (*x)*DTR;
        yy = sin(r)*t;
        zz = cos(r)*t;
        break;
      case MERCATOR:
        if (*x<-PI || *x>PI) return;
        xx = sin(2*atan(exp(-(*y)))-PI/2);
        t = 1-xx*xx;
        if (t<0) return;
        t = sqrt(t);
        yy = sin(*x)*t;
        zz = cos(*x)*t;
        break;
      case MOLLWEIDE:
        t = 1-(*y)*(*y);
        if (t<=0) return;
        t = (*x)/sqrt(t);
        if (t<-2.0 || t>2.0) return;
        t = t*PI/2;
        yy = sin(t);
        zz = cos(t);
        xx = -(*y);
        t = 1-xx*xx;
        if (t<0) return;
        t = sqrt(t);
        yy *= t;
        zz *= t;
        break;
      case AZIMUTH_RANGE:
        r = (*x)/R_EARTH;
        if (r>=1.0) return;
        t = (*y)*DTR;
        zz = cos(r);
        r = sin(r);
        xx = -r*cos(t);
        yy = r*sin(t);
        break;
      case RADAR_AZ_RAN:
        r = (*x)/R_EARTH;
        if (r>=1.0) return;
        t = (*y)*DTR;
        zz = sqrt(1-r*r);
        xx = -r*cos(t);
        yy = r*sin(t);
        break;
   }/*end switch*/


/* the 3d vector (xx,yy,zz) is a unit vector on the earth's surface such that
   (0,0,1) is at the north pole, (1,0,0) is at the intersection of the equator
   and the Greenwich meridian and (0,1,0) is at 0N, 90E */

   xl = geoinfo->aaa*xx+geoinfo->ddd*yy+geoinfo->ggg*zz;
   yl = geoinfo->bbb*xx+geoinfo->eee*yy+geoinfo->hhh*zz;
   zl = geoinfo->ccc*xx+geoinfo->fff*yy+geoinfo->iii*zz;

   if (xl==0.0 && yl==0.0)
       *lon = geoinfo->cenlon;
   else
       *lon = atan2(yl,xl)*RTD;
   *lat = atan2(zl,sqrt(xl*xl+yl*yl))*RTD;

}
