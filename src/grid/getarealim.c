#include "geoLib.h"
#include "geodefs.h"
#include <stdio.h>

void get_area_limits(float * x_min, float * x_max,
                    float * y_min, float * y_max,
                    GeoInfo * projparm)
{
   *x_min = projparm->xmin;
   *x_max = projparm->xmax;
   *y_min = projparm->ymin;
   *y_max = projparm->ymax;
}/* end get_area_limits */
