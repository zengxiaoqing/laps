#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <netcdf.h>

/* program in_domain.c */
/* reads in static.nest7grid LAT and LON and determines if lat_in and lon_in are in domain */
int main(int argc, char *argv[])
{

	float lat_in, lon_in, *lat, *lon, ri, rj,maxLat,maxLon,minLat,minLon;
        char  dataroot[256],fname[256];
        int   cdfId,varId,nc_status,istatus,slen,nx,ny,status,i,j;
        size_t start[4], count[4];

        if (argc == 3) {
	  lat_in = atof(argv[1]);
	  lon_in = atof(argv[2]);

/* get LAPS_DATA_ROOT environment variable */
          strcpy(dataroot,getenv("LAPS_DATA_ROOT"));
/*        fprintf(stdout,"dataroot = %s\n",dataroot);                */
          slen = strlen(dataroot); 
          if (slen > 0) {

/* read netCDF file static.nest7grid LAT and LON, x and y dimensions */
/*          strncpy(fname,dataroot,slen);                            */
            strcpy(fname,dataroot);
/*          fprintf(stdout,"fname = %s\n",fname);                    */
            strcat(fname,"/static/static.nest7grid");
/*          fprintf(stdout,"full fname = %s\n",fname);               */

            if( access(fname, F_OK) != 0 ) {
              fprintf(stdout,"The LAPS static file %s does not exist\n",fname);
              status = 0;
              fprintf(stdout,"%d\n",status);
              fprintf(stdout,"%d\n",slen);
              return;
            }

            nc_status = nc_open(fname, NC_NOWRITE , &cdfId);
            if (nc_status != NC_NOERR) {
              fprintf(stdout,"Error opening LAPS static file %s\n",fname);
              status = 0;
              fprintf(stdout,"%d\n",status);
              return;
            }

            nc_status = nc_inq_dimid(cdfId,"x", &varId);
            if (nc_status != NC_NOERR) {
              fprintf(stdout, "Unable to find dimension 'x' in file.\n");
              status = 0;
              fprintf(stdout,"%d\n",status);
              return;
            } else {
              nc_status = nc_inq_dimlen(cdfId, varId, (size_t *)&nx);
              if (nc_status != NC_NOERR) {
                fprintf(stdout, "Unable to retrieve dimension 'x' from file.\n");
                status = 0;
                fprintf(stdout,"%d\n",status);
                return;
              }
            }

            nc_status = nc_inq_dimid(cdfId,"y", &varId);
            if (nc_status != NC_NOERR) {
              fprintf(stdout, "Unable to find dimension 'y' in file.\n");
              status = 0;
              fprintf(stdout,"%d\n",status);
              return;
            } else {
              nc_status = nc_inq_dimlen(cdfId, varId, (size_t *)&ny);
              if (nc_status != NC_NOERR) {
                fprintf(stdout, "Unable to retrieve dimension 'y' from file.\n");
                status = 0;
                fprintf(stdout,"%d\n",status);
                return;
              }
            }

/* malloc lat and lon */

            lat = malloc(nx*ny*sizeof(float));
            lon = malloc(nx*ny*sizeof(float));

/* setup start and count arrays */
            start[0] = 0;
            start[1] = 0;
            start[2] = 0;
            start[3] = 0;

            count[0] = 1;
            count[1] = 1;
            count[2] = ny;
            count[3] = nx;

            nc_status = nc_inq_varid(cdfId,"lat", &varId);
            if (nc_status != NC_NOERR) {
              fprintf(stdout, "Unable to find variable 'lat' in file.\n");
              status = 0;
              free(lat);
              free(lon);
              fprintf(stdout,"%d\n",status);
              return;
            } else {
              nc_status = nc_get_vara_float(cdfId, varId, start, count, lat);
              if (nc_status != NC_NOERR) {
                fprintf(stdout, "Unable to retrieve variable 'lat' from file.\n");
                status = 0;
                free(lat);
                free(lon);
                fprintf(stdout,"%d\n",status);
                return;
              }
            }

            nc_status = nc_inq_varid(cdfId,"lon", &varId);
            if (nc_status != NC_NOERR) {
              fprintf(stdout, "Unable to find variable 'lon' in file.\n");
              status = 0;
              free(lat);
              free(lon);
              fprintf(stdout,"%d\n",status);
              return;
            } else {
              nc_status = nc_get_vara_float(cdfId, varId, start, count, lon);
              if (nc_status != NC_NOERR) {
                fprintf(stdout, "Unable to retrieve variable 'lon' from file.\n");
                status = 0;
                free(lat);
                free(lon);
                fprintf(stdout,"%d\n",status);
                return;
              }
            }

            nc_status = nc_close(cdfId);

/* call fortran subroutine latlon_to_rlapsgrid...if return value = 1, is in domain 
            latlon_to_rlapsgrid(&lat_in,&lon_in,lat,lon,&nx,&ny,&ri,&rj,&istatus);
*/

/* determine max/min lat/lon for crude bounding box */
            maxLat = -90.0;
            minLat = 90.0;
            maxLon = -1000.0;
            minLon = 1000.0;

/*            printf("nx %d ny %d\n",nx,ny);
*/
            for (j = 0; j < nx; j++) {
              for (i = 0; i < ny; i++) {
/*                printf("Lat %f Lon %f\n",*(lat+(nx*j+i)),*(lon+(nx*j+i)));
*/
                if (*(lat+nx*j+i) > maxLat) maxLat = *(lat+nx*j+i);
                if (*(lat+nx*j+i) < minLat) minLat = *(lat+nx*j+i);
                if (*(lon+nx*j+i) > maxLon) maxLon = *(lon+nx*j+i);
                if (*(lon+nx*j+i) < minLon) minLon = *(lon+nx*j+i);
              }
            }

/*            printf("Lat %f %f Lon %f %f\n",minLat,maxLat,minLon,maxLon);
*/

            free(lat);
            free(lon);

            if ((lat_in >= minLat) && (lat_in <= maxLat) &&
                (lon_in >= minLon) && (lon_in <= maxLon)) {
              status = 1;
            } else {
              status = 0;
            }

          } else {
            printf("Environment variable LAPS_DATA_ROOT not set.\n");
            status = 0;
          }
        } else {
          printf("Not enough command line arguments.\n");
          status = 0;
        }

        printf("%d\n",status);
        return;

}
