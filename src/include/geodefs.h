
#define PI 3.1415926535897931
#define TWOPI 6.2831853071795862
#define RTD 57.295779513082323
#define DTR 0.017453292519943295

#define STEREOGRAPHIC 1
#define ORTHOGRAPHIC 2
#define LAMBERT_CONFORMAL 3
#define AZIMUTHAL_EQUAL_AREA 4
#define GNOMONIC 5
#define AZIMUTHAL_EQUIDISTANT 6
#define SATELLITE_VIEW 7
#define CYLINDRICAL_EQUIDISTANT 8
#define MERCATOR 9
#define MOLLWEIDE 10
#define AZIMUTH_RANGE 16
#define RADAR_AZ_RAN 17

#define WHOLE_PROJECTION_AREA 1
#define LAT_LON_CORNERS_AREA 2
#define CARTESIAN_CORNERS_AREA 3

#define R_EARTH 6371.2

struct allgeoinfo {
   float xmin;
   float xmax;
   float ymin;
   float ymax;
   int   proj_idx;
   float cenlat;
   float cenlon;
   float rot;
   float cone;
   float polesgn;
   int   altproj;
   float aaa;
   float bbb;
   float ccc;
   float ddd;
   float eee;
   float fff;
   float ggg;
   float hhh;
   float iii;
};

