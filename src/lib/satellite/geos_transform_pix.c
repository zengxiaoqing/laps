#include <math.h>


/*
  20121002 graemem@ssec.wisc.edu

  Initial transcod from GOESTransform.java, provided by Tom Rink, 20120705.
  
  Transforms GOES-R Fixed Grid Format (FGF) x and y coordinates to Earth coordinates 
  (i.e. lon and lat), and vice versa. Transformations are done by first calculating 
  intermediate satellite viewing angle coordinates, called lamda and theta.

       FGF  -->  satellite view angle coordinates  -->  Earth
       Earth  -->  satellite view angle coordinates  -->  FGF

  Entry routines are fgf_to_earth_ and earth_to_fgf_. These routines are designed to
  be callable from Fortran, independent of any framework.

  Using the GRS80 constants.
  
  
  20140509 - wstraka@ssec.wisc.edu
    
    Added in the transforms for EUMETSAT/JMA code. This means that the WGS84 constants had to be added
    into the code as well. In addition, there is an assumption
    of how the data is organized in the x direction for ABI which isn't consistant with
    the CGMS orientation. Locations are marked with WCS3.
  

 */

enum SCAN_GEOMETRIES {
  GOES,
  GEOS
};

#define PI (3.141592653589793238462643)
#define DEG_TO_RAD (PI/180.0)
#define RAD_TO_DEG (180.0/PI)

#define R_POL (6356.7523)     /* semi-minor axis (polar radius km) */
#define R_EQ (6378.1370)       /* semi-major axis (equatorial radius km) */
#define F (1.0/298.257222101)  /* flattening */
#define FP (1.0/((1.0-F)*(1.0-F)))

#define H_MSG (42164.0)
#define H_GOESR (42164.16)
/*#define H (H_GOESR) /* Geostationary Orbit Radius (spacecraft to barycenter distance) (km) */
/*#define D (H*H - R_EQ*R_EQ)*/

/* FUTURE: accept as an input if needed */
/*#define ORIGINAL_SCAN_GEOM (GOES)*/

/* not easy to generate a NaN in c, so using the Geocat double missing value */
#define MISSING_VALUE_DOUBLE (-999.0)

/**
 * Transform fractional FGF coordinates to (lamda, theta) radians for ABI.
 *
 * @param[in] fgf_x fractional FGF coordinate, zero-based
 * @param[in] fgf_y fractional FGF coordinate, zero-based
 * @param[in] scale_x scaleFactor from the x coordinate variable
 * @param[in] offset_x addOffset from the x coordinate variable
 * @param[in] scale_y scaleFactor from the y coordinate variable
 * @param[in] offset_y addOffset from the y coordinate variable
 * @param[out] lamda units: radians
 * @param[out] theta units: radians
 */
void fgf_to_sat_abi(const double fgf_x, const double fgf_y, const double scale_x, 
                const double offset_x, const double scale_y, const double offset_y, 
                double *lamda, double *theta) {
      *lamda = fgf_x*scale_x + offset_x;
      *theta = fgf_y*scale_y + offset_y;
}


/**
 * Transform fractional SEVIRI coordinates to (lamda, theta) radians.
 *
 * @param[in] SEVIRI fractional FGF coordinate, zero-based
 * @param[in] SEVIRI fractional FGF coordinate, zero-based
 * @param[in] offset_x is the COFF variable from HRIT Header
 * @param[in] scale_x is the CFAC variable from HRIT Header
 * @param[in] offset_y is the LOFF variable from HRIT Header
 * @param[in] scale_y is the LFAC variable from HRIT Header
 * @param[out] lamda units: radians
 * @param[out] theta units: radians
 *
 * REMINDER - EUMETSAT ASSUMES LOWER RIGHT as 0,0)
 */
void seviriidx_to_sat(const double fgf_x, const double fgf_y, const double scale_x, 
                const double offset_x, const double scale_y, const double offset_y, 
                double *lamda, double *theta) {
                
       double frnt_factor =  pow( 2, 16 );
                
      *lamda = (frnt_factor * ( fgf_x - offset_x )) / scale_x;     
      *theta = (frnt_factor * ( fgf_y - offset_y )) / scale_y;
}


/**
 * Transform fractional JMA coordinates to (lamda, theta) radians.
 *
 * @param[in] SEVIRI fractional FGF coordinate, zero-based
 * @param[in] SEVIRI fractional FGF coordinate, zero-based
 * @param[in] offset_x is the COFF variable from HRIT Header
 * @param[in] scale_x is the CFAC variable from HRIT Header
 * @param[in] offset_y is the LOFF variable from HRIT Header
 * @param[in] scale_y is the LFAC variable from HRIT Header
 * @param[out] lamda units: radians
 * @param[out] theta units: radians
 *
 * REMINDER - JMA ASSUMES UPPER LEFT as (0,0)
 * recall that lamda and theta are in radians, but JMA uses deg constants
 */
void jmaidx_to_sat(const double fgf_x, const double fgf_y, const double scale_x, 
                const double offset_x, const double scale_y, const double offset_y, 
                double *lamda, double *theta) {
                
       double frnt_factor =  pow( 2, 16 );
                
      *lamda = ((frnt_factor * ( fgf_x - offset_x )) / scale_x) * DEG_TO_RAD; 
      *theta = ((frnt_factor * ( fgf_y - offset_y )) / scale_y) * DEG_TO_RAD;
}


/**
 * Transform view angle coordinates in the GOES scan geometry frame to view angle coordinates
 * in the GEOS scan geometry frame.
 *
 *  @param[in, out] lamda units: radians
 *  @param[in, out] theta units: radians
 */
void goes_to_geos(double *lamda, double *theta) {
     double theta_geos = asin( sin(*theta)*cos(*lamda) );
     double lamda_geos = atan( tan(*lamda)/cos(*theta) );  
     *lamda = lamda_geos;
     *theta = theta_geos;
}

/**
 * Transform satellite view angle coordinates, known as the "intermediate" coordinates in the
 * CGMS Normalized Geostationary Projection, to geographic Earth coordinates.
 *
 * @param[in] lamda lamda (East-West) angle, units: radians
 * @param[in] theta theta (North-South) angle, units: radians
 * @param[in] sub_lon_degrees satellite subpoint longitude, units: degrees
 * @param[out] lon_degrees units: degrees
 * @param[out] lat_degrees units: degrees
 */
void sat_to_earth(const enum SCAN_GEOMETRIES scan_geom, const double lamda, const double theta, const double sub_lon_degrees,
                  double *lon_degrees, double *lat_degrees) {

  double sub_lon_radians = sub_lon_degrees * DEG_TO_RAD;

  double x = lamda;
  double y = theta;
  enum SCAN_GEOMETRIES scan = scan_geom;
  
  
  double h = H_MSG;

  if (scan_geom == GOES) { /* convert from GOES to GEOS for transfrom below */
    goes_to_geos(&x, &y);
    h = H_GOESR;   
  }
  

  double d = (h*h - R_EQ*R_EQ);
  double c1 = (h * cos(x) * cos(y)) * (h * cos(x) * cos(y));
  double c2 = (cos(y) * cos(y) + FP * sin(y) * sin(y)) * d;

  if (c1<c2) {
    *lon_degrees = MISSING_VALUE_DOUBLE;
    *lat_degrees = MISSING_VALUE_DOUBLE;
    return;
  }
     
  double s_d = sqrt(c1 - c2);
     
  double s_n = (h * cos(x) * cos(y) - s_d) / (cos(y) * cos(y) + FP * sin(y) * sin(y));
    
  double s_1 = h - s_n * cos(x) * cos(y);
  double s_2 = s_n * sin(x) * cos(y);
  double s_3 = -s_n * sin(y);

    
  double s_xy = sqrt(s_1*s_1 + s_2*s_2);
  double geographic_lon = atan(s_2/s_1) + sub_lon_radians;

  //WCS3
  double geographic_lat = atan(FP*(s_3/s_xy));

  if (scan_geom == GOES) { /* GOES has this flip of FP, while SEVIRI/MTSAT doesn't */
    geographic_lat = atan(-FP*(s_3/s_xy));
  }
  
 

  *lon_degrees = (RAD_TO_DEG*geographic_lon);
  *lat_degrees = RAD_TO_DEG*geographic_lat;

          
  // force output longitude to -180 to 180 range
  if (*lon_degrees < -180.0) *lon_degrees += 360.0;
  if (*lon_degrees > 180.0) *lon_degrees -= 360.0;
}

/**
 *  Transform fractional line/element coordinates to (longitude, latitude).
 *
 ***********************************
 * @param[in] sat - fixed grid identifier
 * @param[in] fgf_x fractional FGF coordinate, zero-based
 * @param[in] fgf_y fractional FGF coordinate, zero-based
 * For ABI
 * @param[in] scale_x scaleFactor from the x coordinate variable
 * @param[in] offset_x addOffset from the x coordinate variable
 * @param[in] scale_y scaleFactor from the y coordinate variable
 * @param[in] offset_y addOffset from the y coordinate variable
 * For EUMETSAT/JMA
 * @param[in] scale_x is the CFAC variable from HRIT Header
 * @param[in] offset_x is the COFF variable from HRIT Header
 * @param[in] scale_y is the LFAC variable from HRIT Header
 * @param[in] offset_y is the LOFF variable from HRIT Header 
 * @param[in] sub_lon_degrees satellite subpoint longitude, units: degrees
 ***********************************
 * @param[out] lon_degrees units: degrees
 * @param[out] lat_degrees units: degrees
 */
void fgf_to_earth_(const int *sat, const double *fgf_x, const double *fgf_y, const double *scale_x, 
                    const double *offset_x, const double *scale_y, const double *offset_y, 
                    const double *sub_lon_degrees, double *lon_degrees, double *lat_degrees) {

                   
  double lamda, theta;
  enum SCAN_GEOMETRIES scan_geom;

  if (*sat == 1) {
    fgf_to_sat_abi(*fgf_x, *fgf_y, *scale_x, *offset_x, *scale_y, *offset_y, &lamda, &theta);
    scan_geom = GOES;
  }
  else if (*sat == 2) {
    seviriidx_to_sat(*fgf_x, *fgf_y, *scale_x, *offset_x, *scale_y, *offset_y, &lamda, &theta);
    scan_geom = GEOS;
  }
  else if (*sat == 3) {
    jmaidx_to_sat(*fgf_x, *fgf_y, *scale_x, *offset_x, *scale_y, *offset_y, &lamda, &theta);
    scan_geom = GEOS;
  }

  sat_to_earth(scan_geom, lamda, theta, *sub_lon_degrees, lon_degrees, lat_degrees);


}


/**
 * Transform geographic Earth coordinates to satellite view angle coordinate system
 * also known as the "intermediate" coordinate system in CGMS Normalized Geostationary Projection.
 *
 * @param[in] lon_degrees longitude, units: degrees
 * @param[in] lat_degrees latitude, units: degrees
 * @param[in] sub_lon_degrees satellite subpoint longitude, units: degrees
 * @param[out] lamda the x or East-West view angle, units: radians
 * @param[out] theta the y or North_South view angle, units: radians
 */

void earth_to_sat(const enum SCAN_GEOMETRIES scan_geom, const double lon_degrees, const double lat_degrees, const double sub_lon_degrees, double *lamda, double *theta) {


  double h = H_MSG;
  enum SCAN_GEOMETRIES scan = scan_geom;

  if (scan_geom == GOES) { /* USE GOES-R Height */
    double h = H_GOESR;    
  }
    
  double geographic_lat = lat_degrees*DEG_TO_RAD;
  double geographic_lon = lon_degrees*DEG_TO_RAD;
  double sub_lon = sub_lon_degrees * DEG_TO_RAD;

  double geocentric_lat = atan(((R_POL*R_POL)/(R_EQ*R_EQ))*tan(geographic_lat));

  double r_earth = R_POL/sqrt(1.0 -((R_EQ*R_EQ - R_POL*R_POL)/(R_EQ*R_EQ))*cos(geocentric_lat)*cos(geocentric_lat));

  double r_1 = h - r_earth*cos(geocentric_lat)*cos(geographic_lon - sub_lon);
  double r_2 = -r_earth*cos(geocentric_lat)*sin(geographic_lon - sub_lon);
  double r_3 = r_earth*sin(geocentric_lat);

  *lamda = MISSING_VALUE_DOUBLE;
  *theta = MISSING_VALUE_DOUBLE;
  


  if (r_1 > h) { // often two geoid intersect points, use the closer one.
    return;
  }

  if (scan_geom == GEOS) { // GEOS (eg. SEVIRI, MSG)  CGMS 03, 4.4.3.2, Normalized Geostationary Projection
  //WCS3 - theta was -r_3 in the original EUMETSAT/JMA code

    *lamda = atan(-r_2/r_1);
    *theta = asin(-r_3/sqrt(r_1*r_1 + r_2*r_2 + r_3*r_3));
  }
  else if (scan_geom == GOES) { // GOES (eg. GOES-R ABI) 
    *lamda = asin(-r_2/sqrt(r_1*r_1 + r_2*r_2 + r_3*r_3));
    *theta = atan(r_3/r_1);
  }
 
/*  printf("%f\n",*lamda);
  printf("%f\n",*theta);*/

}

/**
 *  Transform (lamda, theta) in radians to fractional FGF coordinates for ABI.
 * @param[in] lamda the x or East-West view angle, units: radians
 * @param[in] theta the y or North_South view angle, units: radians
 * @param[in] scale_x scaleFactor from the x coordinate variable
 * @param[in] offset_x addOffset from the x coordinate variable
 * @param[in] scale_y scaleFactor from the y coordinate variable
 * @param[in] offset_y addOffset from the y coordinate variable
 * @param[out] fgf_x fractional fgf x coordinate
 * @param[out] fgf_y fractional fgf y coordinate
 */
void sat_to_fgf_abi(const double lamda, const double theta, const double scale_x, const double offset_x, const double scale_y, 
                const double offset_y, double *fgf_x, double *fgf_y) {
     *fgf_x = (lamda - offset_x)/scale_x;
     *fgf_y = (theta - offset_y)/scale_y;
}


/**
 *  Transform (lamda, theta) in radians to fractional EUMETSAT coordinates.
 * @param[in] lamda the x or East-West view angle, units: radians
 * @param[in] theta the y or North_South view angle, units: radians
 * @param[in] offset_x is the COFF variable from HRIT Header
 * @param[in] scale_x is the CFAC variable from HRIT Header
 * @param[in] offset_y is the LOFF variable from HRIT Header
 * @param[in] scale_y is the LFAC variable from HRIT Header
 * @param[out] fgf_x fractional EUMETSAT coordinate
 * @param[out] fgf_y fractional EUMETSAT coordinate
 */
void sat_to_seviriidx(const double lamda, const double theta, const double scale_x, const double offset_x, const double scale_y, 
                const double offset_y, double *fgf_x, double *fgf_y) {
     double frnt_factor =  pow( 2, -16 );
    
     *fgf_x = offset_x + (lamda *  frnt_factor * scale_x);
     *fgf_y = offset_y + (theta *  frnt_factor * scale_y);

/*     printf("fgf_x=%f  fgf_y=%f \n", *fgf_x, *fgf_y);*/


}

/**
 *  Transform (lamda, theta) in radians to fractional JMA coordinates.
 * @param[in] lamda the x or East-West view angle, units: radians
 * @param[in] theta the y or North_South view angle, units: radians
 * @param[in] offset_x is the COFF variable from HRIT Header
 * @param[in] scale_x is the CFAC variable from HRIT Header
 * @param[in] offset_y is the LOFF variable from HRIT Header
 * @param[in] scale_y is the LFAC variable from HRIT Header
 * @param[out] fgf_x fractional EUMETSAT coordinate
 * @param[out] fgf_y fractional EUMETSAT coordinate
 *
 * recall that lamda and theta are in radians, but JMA uses deg constants
 */
void sat_to_jmaidx(const double lamda, const double theta, const double scale_x, const double offset_x, const double scale_y, 
                const double offset_y, double *fgf_x, double *fgf_y) {
     double frnt_factor =  pow( 2, -16 );
    
     *fgf_x = offset_x + (lamda *  frnt_factor * scale_x * RAD_TO_DEG);
     *fgf_y = offset_y + (theta *  frnt_factor * scale_y * RAD_TO_DEG);

/*     printf("fgf_x=%f  fgf_y=%f \n", *fgf_x, *fgf_y);*/


}


/**
 * Transform Earth coordinates (lon,lat) to fractional line/element coordinates.
 *
 ***********************************
 * @param[in] sat - fixed grid identifier
 * @param[in] lon_degrees Longitude, units: degrees
 * @param[in] lat_degrees Latitude, units: degrees
 * For ABI
 * @param[in] scale_x scaleFactor from the x coordinate variable
 * @param[in] offset_x addOffset from the x coordinate variable
 * @param[in] scale_y scaleFactor from the y coordinate variable
 * @param[in] offset_y addOffset from the y coordinate variable
 * For EUMETSAT/JMA
 * @param[in] scale_x is the CFAC variable from HRIT Header
 * @param[in] offset_x is the COFF variable from HRIT Header
 * @param[in] scale_y is the LFAC variable from HRIT Header
 * @param[in] offset_y is the LOFF variable from HRIT Header 
 ***********************************
 * @param[in] sub_lon_degrees satellite subpoint longitude, units: degrees
 * @param[out] fgf_x fractional fgf x coordinate
 * @param[out] fgf_y fractional fgf y coordinate
 
 
 
 */
void earth_to_fgf_(const int  *sat, const double *lon_degrees, const double *lat_degrees, const double *scale_x, const double *offset_x, 
                   const double *scale_y, const double *offset_y, const double *sub_lon_degrees, double *fgf_x, double *fgf_y) {


  double lamda, theta;
  enum SCAN_GEOMETRIES scan_geom;
  
  if (*sat == 1) {
    scan_geom = GOES;
  }
  else if (*sat == 2) {
    scan_geom = GEOS;
  }
  else if (*sat == 3) {
    scan_geom = GEOS;
  }

  earth_to_sat(scan_geom, *lon_degrees, *lat_degrees, *sub_lon_degrees, &lamda, &theta);

  if (*sat == 1) {
    sat_to_fgf_abi(lamda, theta, *scale_x, *offset_x, *scale_y, *offset_y, fgf_x, fgf_y);
  }
  else if (*sat == 2) {
    sat_to_seviriidx(lamda, theta, *scale_x, *offset_x, *scale_y, *offset_y, fgf_x, fgf_y);
  }
  else if (*sat == 3) {
    sat_to_jmaidx(lamda, theta, *scale_x, *offset_x, *scale_y, *offset_y, fgf_x, fgf_y);
  }
  
}
