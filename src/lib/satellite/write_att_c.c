#ifdef hpux
#define _INCLUDE_POSIX_SOURCE
#endif
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <netcdf.h>
#define SYSCMD "ncgen -o %s %s"

#ifdef FORTRANDOUBLEUNDERSCORE
#define write_att_c write_att_c__
#define nstrncpy nstrncpy__
#endif

#ifdef FORTRANUNDERSCORE
#define write_att_c write_att_c_
#define nstrncpy nstrncpy_
#endif

#ifdef __STDC__
void write_att_c(char *f_filename, long *fn_len, char *f_cdlfile, long *cdl_len, 
                 char *f_sat_name, long *sn_len, long *num_atts, 
                 double *orb_att, long *status)
#else
void write_att_c(f_filename, fn_len, f_cdlfile, cdl_len, f_sat_name, sn_len, 
                 num_atts, orb_att, status)
char *f_filename; 
long *fn_len; 
char *f_cdlfile;
long *cdl_len;
char *f_sat_name; 
long *sn_len;
long *num_atts; 
double *orb_att; 
long *status;
#endif
{
        int cdfid, varid, dimid, istatus;
        static char *syscmd, *cdlfile, *filename, *sat_name;
        size_t start[1], count[1];
        size_t dim_len;

/* convert fortran f_filename into C string filename, fortran f_cdlfile
     into C string cdlfile and f_sat_name into C sat_name  */
        filename = malloc(*fn_len + 1);
        nstrncpy(filename,f_filename,*fn_len);
        cdlfile = malloc(*cdl_len + 1);
        nstrncpy(cdlfile,f_cdlfile,*cdl_len);
        sat_name = malloc(*sn_len + 1);
        nstrncpy(sat_name,f_sat_name,*sn_len);

/* check to see if cdl file is there */
        if( access(cdlfile, F_OK) != 0 ) {
          printf("The cdl file %s does not exist\n",cdlfile);
          *status = -2; /* error in file creation */
          return;
        }

/* SYSCMD contains "/usr/local/netcdf/bin/ncgen -o %s %s\0" which
           is 33 char, cdlfile, and filename  + 10 extra  */
        syscmd = malloc((strlen(SYSCMD)+*cdl_len+*fn_len+10) * sizeof(char));
        sprintf(syscmd,SYSCMD, filename, cdlfile);
        free(cdlfile);

/*  create file, then open it */
        system(syscmd);
        istatus = nc_open(filename, NC_WRITE, &cdfid);
        if (istatus != NC_NOERR) {
          *status = -1; /* error opening file  */
          free(filename);
          free(sat_name);
          free(syscmd);
          return;
        }
        else {
          free(syscmd);
          free(filename);
        }

/* check length of sat name less than dimension in file - 1 */
        istatus = nc_inq_dimid(cdfid, "namelen", &dimid);
        if (istatus != NC_NOERR) {
          printf("Cannot find dimension 'namelen' in file.\n");
          free(sat_name);
          *status = -3;
          return;
        } 

        istatus = nc_inq_dimlen(cdfid, dimid, &dim_len);
        if (istatus != NC_NOERR) {
          printf("Cannot read dimension 'namelen' in file.\n");
          free(sat_name);
          *status = -3;
          return;
        }

        if (*sn_len > (dim_len - 1)) {
          sat_name[dim_len - 1] = '\0';
          printf("Truncating  satellite name to %s \n", *sat_name); 
        }

/* write sat name out to file */
        istatus = nc_inq_varid(cdfid,"sat_name",&varid);
        if (istatus != NC_NOERR) {
          printf("Cannot find variable 'sat_name' in file.\n");
          free(sat_name);
          *status = -3;
          return;
        } 

        start[0] = 0;
        count[0] = dim_len - 1;
        istatus = nc_put_vara_text(cdfid, varid, start, count, sat_name);
        if (istatus != NC_NOERR) {
          printf("Cannot write variable 'sat_name' to file.\n");
          free(sat_name);
          *status = -3;
          return;
        } 

        free(sat_name);

/* check length of orb_att array is <= dimension in file  */
        istatus = nc_inq_dimid(cdfid, "num_att", &dimid);
        if (istatus != NC_NOERR) {
          printf("Cannot find dimension 'num_att' in file.\n");
          *status = -3;
          return;
        } 

        istatus = nc_inq_dimlen(cdfid, dimid, &dim_len);
        if (istatus != NC_NOERR) {
          printf("Cannot read dimension 'num_att' in file.\n");
          *status = -3;
          return;
        }

        if (*num_atts > dim_len) {
          printf("Number of elements in orb_att array exceeds available space in file.\n");
          *status = 2;
          return;
        }
        
/* write orbit attitude data out to file */
        istatus = nc_inq_varid(cdfid,"orb_att",&varid);
        if (istatus != NC_NOERR) {
          printf("Cannot find variable 'orb_att' in file.\n");
          *status = -3;
          return;
        } 

        start[0] = 0;
        count[0] = *num_atts;
        istatus = nc_put_vara_double(cdfid, varid, start, count, orb_att);
        if (istatus != NC_NOERR) {
          printf("Cannot write variable 'orb_att' to file.\n");
          *status = -3;
          return;
        }

/* close file and return */
        istatus = nc_close(cdfid);
        if (istatus != NC_NOERR) {
          printf("Error closing netCDF file.\n");
          *status = -3;
          return;
        }
        else {
          *status = 1;
          return;
        }
}
