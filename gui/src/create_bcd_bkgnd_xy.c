/* create_bcd_bkgnd_xy.c generates political and geographical maps
 * using NWS AWIPS binary cartographic data files and original
 * NCAR Graphics supmap based routines.
 *
 * Copyright (C) 1999  James P. Edwards
 * Author: Phil McDonald & Paula McCaslin 18 April 2002 
 *         Adapted for GUI
 *
 */
#include "geodefs.h"
#include "geoLib.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int pipe=0;
#define min(a,b) (a<b)?a:b
#define max(a,b) (a>b)?a:b

FILE * dopen(char * fileName)
{
    int i;
    char scratch[256];
    FILE * fp;
    fp = fopen(fileName, "r");
    if (!fp) {
      printf("File %s does not exist or cannot be opened.\n", fileName);
      return 0;
    }
    
    i = strlen(fileName);
    if (i>3 && strcmp(".gz",fileName+i-3)==0) {
        pipe = 1;
        fclose(fp);
        strcpy(scratch, "gzip -d -c ");
        strcat(scratch, fileName);
        fp = popen(scratch, "r");
    } else if (i>2 && strcmp(".Z",fileName+i-2)==0) {
        pipe = 1;
        fclose(fp);
        strcpy(scratch, "uncompress -c ");
        strcat(scratch, fileName);
        fp = popen(scratch, "r");
    }/*endif*/
    return fp;
}

void dclose(FILE * fp)
{
    if (pipe)
        pclose(fp);
    else
        fclose(fp);
}

void eliminate_short_lines(struct bcdRecord *record, float cornertocorner)
{
  int i, npts;
  int cnt=0;

  if (record->npts<4) return;

  npts=record->npts-2;

  for(i=0;i<npts;i++){
    while(((record->pts[i]).lon-(record->pts[i+1]).lon)*
	  ((record->pts[i]).lon-(record->pts[i+1]).lon)+
	  ((record->pts[i]).lat-(record->pts[i+1]).lat)*
	  ((record->pts[i]).lat-(record->pts[i+1]).lat) <
	  0.0000001*cornertocorner && record->npts>3) {
      (record->pts[i+1]).lon=(record->pts[i+2]).lon;
      (record->pts[i+1]).lat=(record->pts[i+2]).lat;
      record->npts--;
      cnt++;
    }
  }
  /*
  if(cnt>2) eliminate_short_lines(record,cornertocorner);
  */
  
}


void dump_xy_to_tk_cmmds(FILE * vfp, 
                        float *x, float *y, int npts, 
			float *scaleby_x, float *scaleby_y,
		        float xMin, float yMax,
		        float x_frac, float y_frac,
			int lambrt,
			char *color[], 
			int color_index, 
			int dashed_line, 
			char *bcd_tag[])
{
    int i;
    float offset_x,offset_y;
    float dist_1,dist_2;


    if (npts<4) return;

    fprintf(vfp,"$can->createLine(");

    for(i=0;i<npts-1;i++){

      /* add code to eliminate "short vectors" */
      /* eliminate_short_lines(record,cornertocorner); */
      dist_1=( *(x+i) - *(x+i+1) )* ( *(x+i) - *(x+i+1) )+
	     ( *(y+i) - *(y+i+1) )* ( *(y+i) - *(y+i+1) );
      dist_2=0.0000000010; 
      /* dist_2=25*1E9; */ 
       

      if(dist_1 > dist_2)
      {
          offset_x=( (*(x+i) - xMin) * *scaleby_x );
          offset_y=( (yMax - *(y+i)) * *scaleby_y );
	  /* Invert lambrt projection, because LAMBERT statments in setsup.c 
          */
	  if (lambrt) {
             offset_x=x_frac - offset_x;
             offset_y=y_frac - offset_y;
	  }
          fprintf(vfp,"$cx * %f, $cy * %f, \n",offset_x,offset_y);
      }

    }
    
    offset_x=( (*(x+i) - xMin) * *scaleby_x );
    offset_y=( (yMax - *(y+i)) * *scaleby_y );
    /* Invert lambrt projection, because LAMBERT statments in setsup.c 
    */
    if (lambrt) {
       offset_x=x_frac - offset_x;
       offset_y=y_frac - offset_y;
    }


    fprintf(vfp,"$cx * %f, $cy * %f, \n",offset_x,offset_y);
    if (dashed_line == 2) {
      fprintf(vfp,"-fill => '%s', -dash => '.  ', -width => 1, -tags => '%s');\n",
              color[color_index],bcd_tag[color_index]);
    } else {
      fprintf(vfp,"-fill => '%s', -width => 1, -tags => '%s');\n",
              color[color_index],bcd_tag[color_index]);
    }
    
} /* end dump_xy_to_tk_cmmds */


int create_bcd_bkgnd_xy(int *proj_idx,
            float *cenlat,
            float *cenlon,
            float *rot,
            float *lat1,
            float *lon1,
            float *lat2,
            float *lon2,
            char  *bcd_fname,
            int   *cornertype,
            char  *vector_fname,
            char  *setsup_fname)

/*  List of Error return values possible:
     -5 : From setsup_gui.c     proj_idx not in range 1..10
     -6 : From setsup_gui.c     cornertype not 1,2,3
     -7 : From setsup_gui.c     stdlat and stdlat2 must have same sign
     -9 : From create_bcd_bkgnd_xy.c  empty bcd filename  passed in
    -10 : From create_bcd_bkgnd_xy.c  Could not open vector_instructions.tk for writing
    -11 : From create_bcd_bkgnd_xy.c  Bad point count in bcd filename
    -12 : From create_bcd_bkgnd_xy.c  Could not read points in bcd filename
*/

{
  struct proj_parm my_proj, all_my_proj;  
  FILE *ifp, *vfp;
  float cornertocorner;  
  float xMin,xMax,yMin,yMax;
  float xMin_orig,xMax_orig,yMin_orig,yMax_orig;  
  float xxx1,xxx2,yyy1,yyy2;
  int success;
  int lambrt;
  float dx, dy, dist_x, dist_y, scaleby_x, scaleby_y;
  float x_frac, y_frac, pct;
  float lat__, lon__;
  char *bcd_tag[3];
  char *color[3];
  /*
  char setsup_fname[20]="/tmp/setsup.dat"; 
  char vector_fname[30]="/tmp/vector_instructions.tk";
  */
  char *path_ptr;
  char bcd_path[256];
  char init_bcd_fname[80];
  int path_len;
  int icolor=1;
  float lon_diff;

  color[0] = malloc(30);
  color[1] = malloc(30);
  color[2] = malloc(30);
  strcpy(color[0],"white");
  strcpy(color[1],"coral");
  strcpy(color[2],"pink2");
  /*
  strcpy(color[0],"pink3");
  strcpy(color[2],"khaki1");
  */

  bcd_tag[0] = malloc(30);
  bcd_tag[1] = malloc(30);
  bcd_tag[2] = malloc(30);
  strcpy(bcd_tag[0],"latlon_lines");
  strcpy(bcd_tag[1],"map_lines");
  strcpy(bcd_tag[2],"bcd_lines");

  /* Separate bcd fname and its path, 
   * i.e. find the last slash (/) then separate the two file name. */
  strcpy(bcd_path,bcd_fname);
  path_ptr=strrchr(bcd_path,'/');
  path_ptr++;
  strcpy(init_bcd_fname,path_ptr);
  *path_ptr = '\0';


  vfp = fopen(vector_fname,"w");
  if (vfp==NULL) {
/*  Could not open vector_instructions.tk for writing. */
    success = -10;
    return success;
  }
  fprintf(vfp,"sub vector_instructions {\n");
  fprintf(vfp,"my($can,$cx,$cy) = @_;\n");
  fprintf(vfp,"");

  setsup(proj_idx,cenlat,cenlon,rot,lat1,lon1,lat2,lon2,
            cornertype,&success);
  if (success <= -5)  return success;
  get_proj_parm(&my_proj);

  if (*proj_idx == 3 && *cenlat < 0 ) {
     /* Lambert in Southern Hemisphere.  */
     /*
     get_inclusive_extents_ll(vfp, &my_proj, lat1, lon1, lat2, lon2, 
	   &xMax_orig, &xMin_orig, &yMax_orig, &yMin_orig);
     */
     lambrt=1;
     ll_to_xy(lat2,lon2, &xMax_orig,&yMax_orig, &my_proj);
     ll_to_xy(lat1,lon1, &xMin_orig,&yMin_orig, &my_proj);

  } else {
     /* All others. */
     lambrt=0;
     ll_to_xy(lat1,lon1, &xMax_orig,&yMax_orig, &my_proj);
     ll_to_xy(lat2,lon2, &xMin_orig,&yMin_orig, &my_proj);
  }


  /* Create buffer zone around the domain that is 15% bigger (1.15).*/
  pct=1.15;
  dx=(yMax_orig-yMin_orig)*(pct-1)/2;
  xMax=xMax_orig+dx;
  xMin=xMin_orig-dx;

  dy=(yMax_orig-yMin_orig)*(1-pct)/2;
  yMax=yMax_orig-dy;
  yMin=yMin_orig+dy;


  /* Get distances after new calculation then
   * set variables to scale map in canvas window 
   */
  dist_x=(xMax-xMin);
  dist_y=(yMax-yMin);

  if (dist_x >= dist_y) {
     x_frac=1.;
     y_frac=dist_y/dist_x;

  } else {
     x_frac=dist_x/dist_y;
     y_frac=1.;

  }
  
  scaleby_x= (1/ dist_x) * x_frac;
  scaleby_y= (1/ dist_y) * y_frac;

/*
 * Domain bounding box.
 * 
 *  The creating a rectangle, called bounding box, and calling
 *  subroutine &create_tags (found in "srt_bounding_box.pl") 
 *  allow the rectangle to become active and movable by the
 *  user.
 */
  xxx1=( (xMax_orig-xMin) * scaleby_x );
  xxx2=( (xMin_orig-xMin) * scaleby_x );
  yyy1=( (yMax-yMax_orig) * scaleby_y );
  yyy2=( (yMax-yMin_orig) * scaleby_y );

  fprintf(vfp,"$can->createRectangle(\n");
  fprintf(vfp,"$cx * %f, $cy * %f, $cx * %f, $cy * %f,\n",xxx1,yyy1,xxx2,yyy2);
  fprintf(vfp,"-outline => %s, -tags => 'bbox');\n",color[0]);


/* Create political and geographical maps */

  cornertocorner=(*lat1-*lat2)*(*lat1-*lat2)+(*lon1-*lon2)*(*lon1-*lon2);
  /*xy_to_ll(&(oneRecord.pts[ii].lat),&(oneRecord.pts[ii].lon), x+ii,y+ii, &my_proj); WRONG order...*/
  /* printf("cornertocorner = %f\n",cornertocorner); */


 if ( strncmp(init_bcd_fname, "US_", 3) == 0) { 

  /* Create political boundary reference lines. */
  strcpy(bcd_fname,bcd_path);
  strcat(bcd_fname,init_bcd_fname);

  bcd_vectors(bcd_fname, cornertype, vfp, 
            &xMin, &xMax, &yMin, &yMax,
            &scaleby_x, &scaleby_y, 
	    &x_frac, &y_frac, &my_proj, lambrt, 
	    color, 2, 1, bcd_tag);

  if ( strcmp(init_bcd_fname, "US_States") != 0) { 
     /* Create continent lines. */
     strcpy(bcd_fname,bcd_path);
     strcat(bcd_fname,"US_States.bcd");

     bcd_vectors(bcd_fname, cornertype, vfp, 
               &xMin, &xMax, &yMin, &yMax,
               &scaleby_x, &scaleby_y, 
	       &x_frac, &y_frac, &my_proj, lambrt, 
	       color, 1, 1, bcd_tag);
  }

 } else {

  /* Create continent lines. */
  strcpy(bcd_fname,bcd_path);
  strcat(bcd_fname,"US_States.bcd");

  bcd_vectors(bcd_fname, cornertype, vfp, 
            &xMin, &xMax, &yMin, &yMax,
            &scaleby_x, &scaleby_y, 
	    &x_frac, &y_frac, &my_proj, lambrt, 
	    color, 1, 1, bcd_tag);

  /* Create political boundary reference lines. */
  strcpy(bcd_fname,bcd_path);
  strcat(bcd_fname,init_bcd_fname);

  bcd_vectors(bcd_fname, cornertype, vfp, 
            &xMin, &xMax, &yMin, &yMax,
            &scaleby_x, &scaleby_y, 
	    &x_frac, &y_frac, &my_proj, lambrt, 
	    color, 2, 1, bcd_tag);

 };

  /* Create ocean lat/lon reference lines. */
  strcpy(bcd_fname,bcd_path);
  lon_diff=abs(*lon2 - *lon1);
  if (lon_diff > 51) {
     strcat(bcd_fname,"LatLon_20deg.bcd");
  } else {
     strcat(bcd_fname,"LatLon_10deg.bcd");
  }

  /*
  strcat(bcd_fname,"LatLon_10deg_over_Oceans.bcd");
  */

  bcd_vectors(bcd_fname, cornertype, vfp, 
            &xMin, &xMax, &yMin, &yMax,
            &scaleby_x, &scaleby_y, 
	    &x_frac, &y_frac, &my_proj, lambrt, 
	    color, 0, 2, bcd_tag);



/* Draw a Rectange around map data. */
  xxx1=( (xMax-xMin) * scaleby_x );
  xxx2=( 0           * scaleby_x );
  yyy1=( 0           * scaleby_y );
  yyy2=( (yMax-yMin) * scaleby_y );

  fprintf(vfp,"$can->createRectangle(\n");
  fprintf(vfp,"$cx * %f, $cy * %f, \
	       $cx * %f, $cy * %f,\n",xxx1, yyy1, xxx2, yyy2);
  fprintf(vfp,"-outline => %s, -width => 3, \
	       -tags => 'box_line');\n",color[1]);
 
  /* Draw a Rectange to mask external gridlines. */
  /* left side */
  xxx1=-10.;
  xxx2=( (xMin-xMin) * scaleby_x );
  yyy1=-10;
  yyy2=( (yMax-yMin) * scaleby_y );
  fprintf(vfp,"$can->createRectangle(\n");
  fprintf(vfp,"$cx * %f, $cy * %f, \
	       $cx * %f, $cy * %f,\n",xxx1, yyy1, xxx2, yyy2);
  fprintf(vfp,"-outline => $bg_canvas, \
	       -fill => $bg_canvas, -tags => 'masked_area');\n");

  /* right side */
  xxx1=( (xMax-xMin + 10) * scaleby_x );
  xxx2=( (xMax-xMin) * scaleby_x );
  yyy1=( (yMax-yMin + 10) * scaleby_y );
  yyy2=-2;
  fprintf(vfp,"$can->createRectangle(\n");
  fprintf(vfp,"$cx * %f, $cy * %f, \
	       $cx * %f, $cy * %f,\n",xxx1, yyy1, xxx2, yyy2);
  fprintf(vfp,"-outline => $bg_canvas, \
	       -fill => $bg_canvas, -tags => 'masked_area');\n");

  /* top */
  xxx1=0.;
  xxx2=( (xMax-xMin) * scaleby_x );
  yyy1=-2;
  yyy2=( (yMax-yMax) * scaleby_y );
  fprintf(vfp,"$can->createRectangle(\n");
  fprintf(vfp,"$cx * %f, $cy * %f, \
	       $cx * %f, $cy * %f,\n",xxx1, yyy1, xxx2, yyy2);
  fprintf(vfp,"-outline => $bg_canvas, \
	       -fill => $bg_canvas, -tags => 'masked_area');\n");

  /* bottom */
  xxx1=( (xMax-xMin) * scaleby_x );
  xxx2=-10.;
  yyy1=( (yMax-yMin) * scaleby_y );
  yyy2=( (yMax-yMin + 10) * scaleby_y );
  fprintf(vfp,"$can->createRectangle(\n");
  fprintf(vfp,"$cx * %f, $cy * %f, \
	       $cx * %f, $cy * %f,\n",xxx1, yyy1, xxx2, yyy2);
  fprintf(vfp,"-outline => $bg_canvas, \
	       -fill => $bg_canvas, -tags => 'masked_area');\n");

  /* fake-out */
   /*
   */
  xxx1=( (xMax-xMin) * scaleby_x );
  xxx2=( 0           * scaleby_x );
  yyy1=( 0           * scaleby_y );
  yyy2=( (yMax-yMin) * scaleby_y );
  fprintf(vfp,"$can->createRectangle(\n");
  fprintf(vfp,"$cx * %f, $cy * %f, \
	       $cx * %f, $cy * %f,\n",xxx1, yyy1, xxx2, yyy2);
  fprintf(vfp,"-outline => $bg_canvas, \
	       -fill => $bg_canvas, -tags => 'fake_area');\n");


  fprintf(vfp,"}\n");

/* Create latitude and longitude lines to give reference 
 * coordinates to political and geographical maps
 */
  location_lines(vfp, lat1, lon1, lat2, lon2,
                 &xMin, &xMax, &yMin, &yMax, 
		 &scaleby_x, &scaleby_y,
		 &x_frac, &y_frac,
                 icolor, &my_proj,
		 lambrt);

  fprintf(vfp,"$can->raise('location_lines');\n");
  fprintf(vfp,"$can->raise('box_line');\n");
  fprintf(vfp,"$can->raise('bbox');\n"); 
  fprintf(vfp,"$can->raise('rtags');\n"); 
  fprintf(vfp,"$can->raise('masked_area');\n");

  fprintf(vfp,"}\n");


  fprintf(vfp,"\n# Variables for map_utils.pl called via srt_localize.pl\n");
  fprintf(vfp,"$latsw=%f;\t",*lat2);
  fprintf(vfp,"$lonsw=%f;\n",*lon2);

  fprintf(vfp,"$x_frac=%f;\t",x_frac);
  fprintf(vfp,"$y_frac=%f;\n",y_frac);

  xxx1=xMin;
  yyy1=yMin;
  xy_to_ll(&xxx1,&yyy1, &lat__,&lon__, &my_proj);
  fprintf(vfp,"$latsw=%f;\t",lat__);
  fprintf(vfp,"$lonsw=%f;\n",lon__);

  xxx1=xMin;
  yyy1=yMax;
  xy_to_ll(&xxx1,&yyy1, &lat__,&lon__, &my_proj);
  fprintf(vfp,"$latnw=%f;\t",lat__);
  fprintf(vfp,"$lonnw=%f;\n",lon__);

  xxx1=xMax;
  yyy1=yMax;
  xy_to_ll(&xxx1,&yyy1, &lat__,&lon__, &my_proj);
  fprintf(vfp,"$latne=%f;\t",lat__);
  fprintf(vfp,"$lonne=%f;\n",lon__);

  xxx1=xMax;
  yyy1=yMin;
  xy_to_ll(&xxx1,&yyy1, &lat__,&lon__, &my_proj);
  fprintf(vfp,"$latse=%f;\t",lat__);
  fprintf(vfp,"$lonse=%f;\n\n",lon__);

  fprintf(vfp,"$scaleby_x=%f;\t", scaleby_x);
  fprintf(vfp,"$scaleby_y=%f;\n", scaleby_y);
  fprintf(vfp,"$scale_xMin=%f;\t", xMin);
  fprintf(vfp,"$scale_yMax=%f;\n", yMax);

  fclose (vfp);

/*  
save_geo_file2(setsup_fname, &success);
*/
  save_geo_file2(setsup_fname);

  free(*color);
  free(*bcd_tag);

  return 1;
  }

int bcd_vectors(
            char *bcd_fname,
            int *cornertype,
            FILE * vfp, 
            float *xMin,
	    float *xMax,
	    float *yMin,
	    float *yMax,
            float *scaleby_x,
            float *scaleby_y,
            float *x_frac,
            float *y_frac,
            GeoInfo * my_proj2,
	    int lambrt,
            char *color[], 
            int color_index, 
            int line_thickness, 
            char *bcd_tag[]) 
{
  int success;
  FILE *ifp;
  int n, allin, allout,inarea[MAXPTS],ii,i1,i2;  
  float minlat,minlon,maxlat,maxlon,minreclat,minreclon,maxreclat,maxreclon;
  struct bcdRecord oneRecord;
  struct bcdRecord clipRec;
  float ptX,ptY,dist;
  float x[MAXPTS],y[MAXPTS];
  float xx[MAXPTS],yy[MAXPTS],lat,lon;
  float w,wy;

  ifp = dopen(bcd_fname);

  if(bcd_fname==NULL){
/*  Need to enter a bcd filename. */
    success = -9;
    return success;
  }


  minlat = 90.0;
  minlon = 180.0;
  maxlat = 0.0;
  maxlon = -180.0;
  while (1) {

    if (my_fread(&n,4,1,ifp)!=1) break;
    if (n<2 || n>MAXPTS) {
/*    Bad point count in bcd filename. */
      success = -11;
      return success;
    }/*endif*/
    oneRecord.npts = n;
    n = n*2+4;
    if (my_fread(&(oneRecord.lat1),4,n,ifp)!=n) {
/*    Could not read points in bcd filename. */
      success = -12;
      return success;
    }/*endif*/

    allin = allout = 1;

    minreclat = 90.0;
    minreclon = 180.0;
    maxreclat = 0.0;
    maxreclon = -180.0;
    for (ii=0; ii<oneRecord.npts; ii++) {
      ll_to_xy(&(oneRecord.pts[ii].lat),&(oneRecord.pts[ii].lon),
	       x+ii,y+ii, my_proj2);
      if (x[ii]>=*xMin && x[ii]<=*xMax && y[ii]>=*yMin && y[ii]<=*yMax) {
	allout = 0;
	inarea[ii] = 1;
	minreclat = min(oneRecord.pts[ii].lat,minreclat);
	minreclon = min(oneRecord.pts[ii].lon,minreclon);
	maxreclat = max(oneRecord.pts[ii].lat,maxreclat);
	maxreclon = max(oneRecord.pts[ii].lon,maxreclon);
      } else {
	allin = inarea[ii] = 0;
      }
    }/*end for*/
    if (allin) {
      minlat = min(minlat,minreclat);
      minlon = min(minlon,minreclon);
      maxlat = max(maxlat,maxreclat);
      maxlon = max(maxlon,maxreclon);

      dump_xy_to_tk_cmmds(vfp,x,y,oneRecord.npts,
                      scaleby_x, scaleby_y, 
		      *xMin,*yMax,*x_frac,*y_frac, lambrt, 
		      color, color_index, line_thickness, bcd_tag);
      continue;
    } else if (allout)
      continue;
    
    clipRec.npts = 0;
    for (ii=0; ii<oneRecord.npts; ii++) {
      if (ii==0 || inarea[ii]==inarea[ii-1]) {
	if (!inarea[ii]) continue;
        xx[clipRec.npts] = x[ii];
        yy[clipRec.npts] = y[ii];
	clipRec.pts[clipRec.npts].lat = oneRecord.pts[ii].lat;
	clipRec.pts[clipRec.npts++].lon = oneRecord.pts[ii].lon;
	continue;
      }/*endif*/
      if (inarea[ii]) {
	i1 = ii;
	i2 = ii-1;
	minlat = min(oneRecord.pts[ii].lat,minlat);
	minlon = min(oneRecord.pts[ii].lon,minlon);
	maxlat = max(oneRecord.pts[ii].lat,maxlat);
	maxlon = max(oneRecord.pts[ii].lon,maxlon);
      } else {
	i1 = ii-1;
	i2 = ii;
      }/*endif*/
      if (x[i2]>*xMax)
	w = (*xMax-x[i1])/(x[i2]-x[i1]);
      else if (x[i2]<*xMin)
	w = (*xMin-x[i1])/(x[i2]-x[i1]);
      else
	w = 1.0;
      if (y[i2]>*yMax)
	wy = (*yMax-y[i1])/(y[i2]-y[i1]);
      else if (y[i2]<*yMin)
	wy = (*yMin-y[i1])/(y[i2]-y[i1]);
      else
	wy = 1.0;
      if (wy<w) w = wy;
      xx[clipRec.npts] = w*x[i2]+(1-w)*x[i1];
      yy[clipRec.npts] = w*y[i2]+(1-w)*y[i1];
      xy_to_ll(&xx[clipRec.npts],&yy[clipRec.npts],&lat,&lon, my_proj2); 
      clipRec.pts[clipRec.npts].lat = lat;
      clipRec.pts[clipRec.npts++].lon = lon;
      if (inarea[ii]) {
        xx[clipRec.npts] = x[ii];
        yy[clipRec.npts] = y[ii];
	clipRec.pts[clipRec.npts].lat = oneRecord.pts[ii].lat;
	clipRec.pts[clipRec.npts++].lon = oneRecord.pts[ii].lon;
      } else {


      dump_xy_to_tk_cmmds(vfp,x,y,oneRecord.npts,
                      scaleby_x, scaleby_y, 
		      *xMin,*yMax,*x_frac,*y_frac, lambrt, 
		      color, color_index, line_thickness, bcd_tag);
	clipRec.npts = 0;
      }/*endif*/
    }/*end for ii*/
    if (clipRec.npts>1) {
      dump_xy_to_tk_cmmds(vfp,x,y,oneRecord.npts,
                      scaleby_x, scaleby_y, 
		      *xMin,*yMax,*x_frac,*y_frac, lambrt, 
		      color, color_index, line_thickness, bcd_tag);
    }


    
  }/*end while*/
  dclose (ifp);

  if(maxlat-minlat>maxlon-minlon){
    minlon -= 0.5*((maxlat-minlat)-(maxlon-minlon));
    maxlon += 0.5*((maxlat-minlat)-(maxlon-minlon)) ;
  }else{
    minlat -= 0.5*((maxlon-minlon)-(maxlat-minlat));
    maxlat += 0.5*((maxlon-minlon)-(maxlat-minlat));
  }


  return 1;
}


/* 16Nov99  Phil McDonald */
int	my_fread (void *pdat, long size, long n, FILE *file)
{

    int		nread;



    nread = fread (pdat, size, n, file);


#ifdef BYTE_SWAP || SWAPBYTE

    if (size > 1)
    {
        unsigned char	byte, *ptmp, *pend, *pi, *pj;
	long		s1;

        s1   = size - 1;
        ptmp = pdat;
        pend = (unsigned char *) pdat + (n * size);

        while (ptmp < pend)
        {
            pi = ptmp;
            pj = ptmp + s1;
            while (pi < pj)
            {
                byte = *pi, *pi = *pj, *pj = byte;
                pi++;
                pj--;
            }

            ptmp += size;
        }
    }
#endif


    return nread;
}
