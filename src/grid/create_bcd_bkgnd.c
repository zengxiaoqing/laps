#include "geodefs.h"
#include "geoLib.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int pipe=0;
#define min(a,b) (a<b)?a:b
#define max(a,b) (a>b)?a:b
#ifdef FORTRANUNDERSCORE
#define create_bcd_bkgnd create_bcd_bkgnd_
#endif

FILE * dopen(char * fileName)
{
    int i;
    char scratch[256];
    FILE * fp;
    fp = fopen(fileName, "r");
    if (!fp) return 0;
    
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

  if (record->npts<2) return;

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



void dump_record_to_gnuplot(FILE * fp, struct bcdRecord * record, float cornertocorner)
{
    int i;

    /*
    eliminate_short_lines(record,cornertocorner);
    */
    if (record->npts<2) return;

    /* add code to eliminate "short vectors" */
    for(i=0;i<record->npts-1;i++){
      
      /*      if(((record->pts[i]).lon-(record->pts[i+1]).lon)*
         ((record->pts[i]).lon-(record->pts[i+1]).lon)+
	 ((record->pts[i]).lat-(record->pts[i+1]).lat)*
         ((record->pts[i]).lat-(record->pts[i+1]).lat) >
	 0.00001*cornertocorner) 
      */{
      
	  fprintf(fp,"set arrow from %f,%f to %f,%f nohead ls 2\n"
		  ,(record->pts[i]).lon,(record->pts[i]).lat
		  ,(record->pts[i+1]).lon,(record->pts[i+1]).lat);
	}
    }

} /* end dump_record_to_gnuplot */


void create_bcd_bkgnd(int * proj_idx,
            float * cenlat,
            float * cenlon,
            float * rot,
            float * x1lat1,
            float * x2lon1,
            float * y1lat2,
            float * y2lon2,
            int * cornertype,
            int * success)

{
  struct proj_parm my_proj;
  struct bcdRecord oneRecord;
  struct bcdRecord clipRec;
  FILE *ifp, *ofp;
  float cornertocorner;
  float x1,x2,y1,y2;
  float x[MAXPTS],y[MAXPTS];
  float xx,yy,lat,lon;
  float w,wy;
  int n, allin, allout,inarea[MAXPTS],ii,i1,i2;  
  float minlat,minlon,maxlat,maxlon,minreclat,minreclon,maxreclat,maxreclon;
  char *infname;

  infname = getenv("BCD_FNAME");

  if(infname==NULL){
    printf("specify the bcd filename in env variable BCD_FNAME");
    exit;
  }

  ifp = dopen(infname);
  

  ofp = fopen("mapdata","w");
  if (ofp==NULL) {
    fprintf(stderr,"Could not open mapdata for writing\n");
    exit(1);
  }/*endif*/
  setsup(proj_idx,cenlat,cenlon,rot,x1lat1,x2lon1,y1lat2,y2lon2,
            cornertype,success);

  
  get_proj_parm(&my_proj);
  get_area_limits(&x1, &x2, &y1, &y2, &my_proj);


  cornertocorner=(*x1lat1-*y1lat2)*(*x1lat1-*y1lat2)+(*x2lon1-*y2lon2)*(*x2lon1-*y2lon2);
  printf("cornertocorner = %f\n",cornertocorner);

  minlat = 90.0;
  minlon = 180.0;
  maxlat = 0.0;
  maxlon = -180.0;
  while (1) {

    if (fread(&n,4,1,ifp)!=1) break;
    if (n<2 || n>MAXPTS) {
      fprintf(stderr,"Bad point count %d in %s\n",n,infname);
      exit(1);
    }/*endif*/
    oneRecord.npts = n;
    n = n*2+4;
    if (fread(&(oneRecord.lat1),4,n,ifp)!=n) {
      fprintf(stderr,"Could not read %d points in %s\n",n,infname);
      exit(1);
    }/*endif*/

    allin = allout = 1;

    minreclat = 90.0;
    minreclon = 180.0;
    maxreclat = 0.0;
    maxreclon = -180.0;
    for (ii=0; ii<oneRecord.npts; ii++) {
      ll_to_xy(&(oneRecord.pts[ii].lat),&(oneRecord.pts[ii].lon),
	       x+ii,y+ii, &my_proj);
      if (x[ii]>=x1 && x[ii]<=x2 && y[ii]>=y1 && y[ii]<=y2) {
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

      dump_record_to_gnuplot(ofp,&oneRecord,cornertocorner); 
      continue;
    } else if (allout)
      continue;
    
    clipRec.npts = 0;
    for (ii=0; ii<oneRecord.npts; ii++) {
      if (ii==0 || inarea[ii]==inarea[ii-1]) {
	if (!inarea[ii]) continue;
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
      if (x[i2]>x2)
	w = (x2-x[i1])/(x[i2]-x[i1]);
      else if (x[i2]<x1)
	w = (x1-x[i1])/(x[i2]-x[i1]);
      else
	w = 1.0;
      if (y[i2]>y2)
	wy = (y2-y[i1])/(y[i2]-y[i1]);
      else if (y[i2]<y1)
	wy = (y1-y[i1])/(y[i2]-y[i1]);
      else
	wy = 1.0;
      if (wy<w) w = wy;
      xx = w*x[i2]+(1-w)*x[i1];
      yy = w*y[i2]+(1-w)*y[i1];
      xy_to_ll(&xx,&yy,&lat,&lon, &my_proj);
      clipRec.pts[clipRec.npts].lat = lat;
      clipRec.pts[clipRec.npts++].lon = lon;
      if (inarea[ii]) {
	clipRec.pts[clipRec.npts].lat = oneRecord.pts[ii].lat;
	clipRec.pts[clipRec.npts++].lon = oneRecord.pts[ii].lon;
      } else {
	dump_record_to_gnuplot(ofp,&clipRec,cornertocorner); 
	clipRec.npts = 0;
      }/*endif*/
    }/*end for ii*/
    if (clipRec.npts>1) dump_record_to_gnuplot(ofp,&clipRec,cornertocorner);
    
  }/*end while*/
  if(maxlat-minlat>maxlon-minlon){
    minlon -= 0.5*((maxlat-minlat)-(maxlon-minlon));
    maxlon += 0.5*((maxlat-minlat)-(maxlon-minlon));
  }else{
    minlat -= 0.5*((maxlon-minlon)-(maxlat-minlat));
    maxlat += 0.5*((maxlon-minlon)-(maxlat-minlat));
  }
  /*  fprintf(ofp,"set xrange [%f:%f]\n",minlon,maxlon);
      fprintf(ofp,"set yrange [%f:%f]\n",minlat,maxlat); */
  fclose (ofp);
  dclose (ifp);
}
