#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <netcdf.h>
#define SYSCMD "ncgen -o %s %s"

#ifndef DEBUG
#define DEBUG 0
#endif

#ifdef FORTRANUNDERSCORE
#define write_qc write_qc_
#define nstrncpy nstrncpy_
#endif
#ifdef FORTRANCAPS
#define write_qc WRITE_QC
#endif
#ifdef FORTRANDOUBLEUNDERSCORE
#define nstrncpy nstrncpy__
#define write_qc write_qc__
#endif

#define SUCCESS 1
#define ERROR 0

void write_qc(char *, char *, int *, int *, long *, long *, char *, 
              char *, char *, int *, int *, int *, int *, int *,
              float *, float *, float *, int *, float *, float *, 
              float *, float *, float *, float *, int *, float *, float *,
              float *, float *, float *, float *, int *, float *, float *, 
              float *, float *, float *, float *, float *, float *, float *,
              float *, float *, float *, int *, float *, float *, float *, 
              float *, float *, float *, int *,  float *,  float *, 
              float *, float *, float *, float *, int *);
int write_val_ref(int, long *, long *);
int write_latlon(int, int *, float *, float *, float *);
int write_temp(int, int *, int *, float *, float *, float *, float *, float *, float *);
int write_dewpoint(int, int *, int *, float *, float *, float *, float *, float *, float *);
int write_wind(int, int *, int *, float *, float *, float *, float *, float *, float *, 
               float *, float *, float *, float *, float *, float *);
int write_pmsl(int, int *, int *, float *, float *, float *, float *, float *, float *);
int write_alt(int, int *, int *, float *, float *, float *, float *, float *, float *);
int write_char_var(int, int *, char *, char *, char *, int *, int *, int *);
void free_strings(char *, char *, char *); 
int get_dims(int, int *, int *, int *);

/*******************************************************************************/
void write_qc(char *f_filename, char *f_cdl_path, int *numSta, int *maxSta,
              long *i_reftime, long *i_valtime, char *f_stations, 
              char *f_provider, char *f_reptype, int *sta_len, int *pro_len, 
              int *type_len, int *fn_len, int *cdl_path_len,
              float *lat, float *lon, float *elev, int *qcstat, 
              float *t, float *tb, float *ta, float *tc, float *te,
              float *tf, int *qcstatd, float *td, float *tdb,
              float *tda, float *tdc, float *tde, float *tdf,
              int *qcstauv, float *u, float *ub, float *ua, float *uc,
              float *ue, float *uf, float *v, float *vb, float *va,
              float *vc, float *ve, float *vf, int *qcstapm, 
              float *pmsl, float *pmslb, float *pmsla, float *pmslc, 
              float *pmsle, float *pmslf, int *qcstal,  float *alt,  
              float *altb, float *alta, float *altc, float *alte,
              float *altf, int *status)
{

	char *filename, *cdlfile, *syscmd;
	int cdl_len, istatus;
	int cdfid, varid, dimid; 
	size_t start_1[1], count_1[1], mindex[1], dim_len;

/* convert fortran f_filename into C string filename  */
        filename = malloc(*fn_len + 1);
        nstrncpy(filename,f_filename,*fn_len);

/* allocate space for syscmd and cdlfile and fill up */
        /* cdl filename is "sfc_qc.cdl" */
        cdlfile = malloc((*cdl_path_len + 15) * sizeof(char));
        nstrncpy(cdlfile,f_cdl_path,*cdl_path_len);
	strcat(cdlfile,"/sfc_qc.cdl");
	cdl_len = strlen(cdlfile);

/* check to see if cdl file is there */
        if( access(cdlfile, F_OK) != 0 ) {
          printf("The cdl file %s does not exist\n",cdlfile);
          *status = -2; /* error in file creation */
          free(cdlfile);
	  free(filename);
          return;
        }

/* SYSCMD contains "/usr/local/netcdf/bin/ncgen -o %s %s\0" which
           is 33 char, cdlfile, and filename  + 10 extra  */
        syscmd = malloc((strlen(SYSCMD)+cdl_len+*fn_len+10) * sizeof(char));
        sprintf(syscmd,SYSCMD, filename, cdlfile);
        free(cdlfile);

/*    see if file is already there  */
        if( access(filename, F_OK) != 0 ) { /* file does not exist */

/*  create file, then open it */
          system(syscmd);
          istatus = nc_open(filename,NC_WRITE, &cdfid);
          if (istatus != NC_NOERR) {
            *status = -2; /* error in file creation */
            free(syscmd);
            free(filename);
            return;
          }
        }
        else { /* file is there...write over existing file */
          system(syscmd); 
          istatus = nc_open(filename,NC_WRITE, &cdfid);
          if (istatus != NC_NOERR) {  /* error opening file */
            printf("File %s exists, but cannot be opened.\n",filename);
            *status = -2; /* error in file creation*/
            free(syscmd);
            free(filename);
            return;
          }
        }
        free(syscmd);
        free(filename);

/* file is now open and ready for writing */

/* write valtime and reftime to file */
	istatus = write_val_ref(cdfid, i_valtime, i_reftime); 
	if (istatus != SUCCESS) {
	  *status = istatus;    
	  nc_close(cdfid);
          return;
        }

/* deal with the 3 character variables: f_stations, f_provider, f_reptype */

/* transfer the character data one field at a time  */
	istatus = write_char_var(cdfid, numSta, f_stations, f_provider, 
                                 f_reptype, sta_len, pro_len, type_len);
	if (istatus != SUCCESS) {
	  *status = istatus;    
	  nc_close(cdfid);
          return;
        }

/* transfer rest of data by writing slabs of data */

/* write latitude, longitude and elevation to file */
	istatus = write_latlon(cdfid, numSta,lat, lon, elev);

/* write temperature variables */
	istatus = write_temp(cdfid, numSta, qcstat, t, tb, 
                             ta, tc, tf, te);

/* write dewpoint temperature variables */
	istatus = write_dewpoint(cdfid, numSta, qcstatd, td, tdb, 
                                 tda, tdc, tdf, tde);

/* write wind variables */
	istatus = write_wind(cdfid, numSta, qcstauv, u, ub, ua, 
                             uc, uf, ue, v, vb, va, vc, vf, ve);

/* write MSL pressure variables */
	istatus = write_pmsl(cdfid, numSta, qcstapm, pmsl, pmslb, 
                             pmsla, pmslc, pmslf, pmsle);

/* write altimeter variables */
	istatus = write_alt(cdfid, numSta, qcstal, alt, altb, 
                            alta, altc, altf, alte);


	nc_close(cdfid);

	*status = SUCCESS;
	return;
}

/*******************************************************************************/
int write_val_ref(int cdfid, long *i_valtime, long *i_reftime) 
{
	int varid, istatus;
	size_t mindex[1];
	double d_val;

	mindex[0] = 1;
	
        istatus = nc_inq_varid (cdfid, "valtime", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'valtime'\n");
	  return 2;
	}

	d_val = (double)*i_valtime;
        istatus = nc_put_var1_double(cdfid, varid, mindex, &d_val);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'valtime'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "reftime", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'reftime'\n");
	  return 2;
	}

	d_val = (double)*i_reftime;
        istatus = nc_put_var1_double(cdfid, varid, mindex, &d_val);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'reftime'\n");
	  return 2;
	}
	return (SUCCESS);
}
/*******************************************************************************/
int write_latlon(int cdfid, int *numSta, float *lat, float *lon, float *elev)
{
	int varid, istatus;
	size_t start[1], count[1];

        start[0] = 0;
	count[0] = *numSta;
	
        istatus = nc_inq_varid (cdfid, "latitude", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'latitude'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, lat);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'latitude'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "longitude", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'longitude'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start,count, lon);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'longitude'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "elevation", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'elevation'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, elev);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'elevation'\n");
	  return 2;
	}

	return (SUCCESS);
}
/*******************************************************************************/
int write_temp(int cdfid, int *numSta, int *qcstat, float *t, 
               float *tb, float *ta, float *tc, float *tf, float *te)
{
	int varid, istatus;
	size_t start[1], count[1];

	start[0] = 0;
	count[0] = *numSta;
	
        istatus = nc_inq_varid (cdfid, "qcstat", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'qcstat'\n");
	  return 2;
	}

        istatus = nc_put_vara_int(cdfid, varid, start, count, qcstat);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'qcstat'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "t", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 't'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, t);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 't'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "tb", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'tb'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, tb);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'tb'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "ta", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'ta'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, ta);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'ta'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "tc", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'tc'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, tc);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'tc'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "tf", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'tf'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, tf);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'tf'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "te", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'te'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, te);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'te'\n");
	  return 2;
	}

	return (SUCCESS);
}
/*******************************************************************************/
int write_dewpoint(int cdfid, int *numSta, int *qcstatd, float *td, 
                   float *tdb, float *tda, float *tdc, float *tdf, float *tde)
{
	int varid, istatus;
	size_t start[1], count[1];

	start[0] = 0;
	count[0] = *numSta;
	
        istatus = nc_inq_varid (cdfid, "qcstatd", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'qcstatd'\n");
	  return 2;
	}

        istatus = nc_put_vara_int(cdfid, varid, start, count, qcstatd);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'qcstatd'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "td", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'td'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, td);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'td'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "tdb", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'tdb'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, tdb);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'tdb'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "tda", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'tda'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, tda);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'tda'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "tdc", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'tdc'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, tdc);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'tdc'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "tdf", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'tdf'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, tdf);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'tdf'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "tde", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'tde'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, tde);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'tde'\n");
	  return 2;
	}

	return (SUCCESS);

}
/*******************************************************************************/
int write_wind(int cdfid, int *numSta, int *qcstauv, float *u, 
               float *ub, float *ua, float *uc, float *uf, float *ue, float *v, 
               float *vb, float *va, float *vc, float *vf, float *ve)
{
	int varid, istatus;
	size_t start[1], count[1];

	start[0] = 0;
	count[0] = *numSta;

        istatus = nc_inq_varid (cdfid, "qcstauv", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'qcstauv'\n");
	  return 2;
	}

        istatus = nc_put_vara_int(cdfid, varid, start, count, qcstauv);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'qcstauv'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "u", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'u'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, u);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'u'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "ub", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'ub'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, ub);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'ub'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "ua", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'ua'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, ua);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'ua'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "uc", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'uc'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, uc);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'uc'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "uf", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'uf'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, uf);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'uf'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "ue", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'ue'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, ue);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'ue'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "v", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'v'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, v);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'v'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "vb", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'vb'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, vb);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'vb'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "va", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'va'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, va);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'va'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "vc", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'vc'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, vc);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'vc'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "vf", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'vf'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, vf);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'vf'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "ve", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 've'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, ve);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 've'\n");
	  return 2;
	}

	return (SUCCESS);

}
/*******************************************************************************/
int write_pmsl(int cdfid, int *numSta, int *qcstapm, float *pmsl, 
               float *pmslb, float *pmsla, float *pmslc, float *pmslf, float *pmsle)
{
	int varid, istatus;
	size_t start[1], count[1];

	start[0] = 0;
	count[0] = *numSta;
	
        istatus = nc_inq_varid (cdfid, "qcstapm", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'qcstapm'\n");
	  return 2;
	}

        istatus = nc_put_vara_int(cdfid, varid, start, count, qcstapm);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'qcstapm'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "pmsl", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'pmsl'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, pmsl);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'pmsl'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "pmslb", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'pmslb'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, pmslb);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'pmslb'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "pmsla", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'pmsla'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, pmsla);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'pmsla'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "pmslc", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'pmslc'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, pmslc);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'pmslc'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "pmslf", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'pmslf'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, pmslf);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'pmslf'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "pmsle", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'pmsle'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, pmsle);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'pmsle'\n");
	  return 2;
	}

	return (SUCCESS);

}
/*******************************************************************************/
int write_alt(int cdfid, int *numSta, int *qcstal, float *alt, 
              float *altb, float *alta, float *altc, float *altf, float *alte)
{
	int varid, istatus;
	size_t start[1], count[1];

	start[0] = 0;
	count[0] = *numSta;
	
        istatus = nc_inq_varid (cdfid, "qcstal", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'qcstal'\n");
	  return 2;
	}

        istatus = nc_put_vara_int(cdfid, varid, start, count, qcstal);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'qcstal'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "alt", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'alt'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, alt);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'alt'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "altb", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'altb'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, altb);
        if (istatus != NC_NOERR) {
 	  printf("Can't write variable 'altb'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "alta", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'alta'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, alta);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'alta'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "altc", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'altc'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, altc);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'altc'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "altf", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'altf'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, altf);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'altf'\n");
	  return 2;
	}

        istatus = nc_inq_varid (cdfid, "alte", &varid);
        if (istatus != NC_NOERR) {
 	  printf("Can't find variable 'alte'\n");
	  return 2;
	}

        istatus = nc_put_vara_float(cdfid, varid, start, count, alte);
        if (istatus != NC_NOERR){
 	  printf("Can't write variable 'alte'\n");
	  return 2;
	}

	return (SUCCESS);

}
/*******************************************************************************/
int write_char_var(int cdfid, int *numSta, char *f_stations, 
                   char *f_provider, char *f_reptype, int *sta_len,
                   int *pro_len, int *type_len)
{
	int file_prov, file_staid, file_reptype, recNum;
	int var_pro, var_sta, var_typ, istatus, slen;
        int wrtlen_pro, wrtlen_sta, wrtlen_typ;
	size_t start_pro[2], count_pro[2], start_sta[2], count_sta[2]; 
	size_t start_typ[2], count_typ[2];
	char *station, *repType, *provider, *pro_ptr, *sta_ptr, *typ_ptr;

/* determine dimension values of maxProviderLen, maxStaTypeLen, and maxStaIdLen
     from the output file */ 

	istatus = get_dims(cdfid, &file_prov, &file_staid, &file_reptype);
	if (istatus == -1) {
          return -1; 
        }

/* setup start_ and count_ arrays, and write out data */
	start_pro[1] = 0;
	count_pro[0] = 1;

	if (*pro_len+1 < file_prov) {
          wrtlen_pro = *pro_len+1; 
	}
	else if (*pro_len == file_prov) {
	  wrtlen_pro = *pro_len;
        }
	else {
	  wrtlen_pro = file_prov;
        }

	start_sta[1] = 0;
	count_sta[0] = 1;

	if (*sta_len+1 < file_staid) {
          wrtlen_sta = *sta_len+1; 
	}
	else if (*sta_len == file_staid) {
	  wrtlen_sta = *sta_len;
        }
	else {
	  wrtlen_sta = file_staid;
        }

	start_typ[1] = 0;
	count_typ[0] = 1;

	if (*type_len+1 < file_reptype) {
          wrtlen_typ = *type_len+1; 
	}
	else if (*type_len == file_reptype) {
	  wrtlen_typ = *type_len;
        }
	else {
	  wrtlen_typ = file_reptype;
        }

/* get var_pro, var_sta, var_typ, varids for each variable */
        istatus = nc_inq_varid (cdfid, "dataProvider", &var_pro);
        if (istatus != NC_NOERR) {
	  return -1;
	}

        istatus = nc_inq_varid (cdfid, "stationId", &var_sta);
        if (istatus != NC_NOERR) {
	  return -1;
	}

        istatus = nc_inq_varid (cdfid, "reportType", &var_typ);
        if (istatus != NC_NOERR) {
	  return -1;
	}

/* setup character strings to hold results of nstrncpy */
	station = malloc((wrtlen_sta) * sizeof(char));
	repType = malloc((wrtlen_typ) * sizeof(char));
	provider = malloc((wrtlen_pro) * sizeof(char));

	pro_ptr = f_provider;
	sta_ptr = f_stations;
	typ_ptr = f_reptype;

/* loop over numSta */
	for (recNum = 0; recNum < *numSta; recNum++) {

/* write provider */
	  start_pro[0] = recNum;
          nstrncpy(provider, pro_ptr, wrtlen_pro - 1);

	  slen = strlen(provider);
	  if (slen < wrtlen_pro) {
            count_pro[1] = slen + 1;
          } 
          else {
            count_pro[1] = wrtlen_pro;
          }
          istatus = nc_put_vara_text(cdfid,var_pro,start_pro,count_pro,provider);
          if (istatus != NC_NOERR) {
	    free_strings(station, repType, provider);
	    return 2;
	  }
	  pro_ptr += *pro_len;

/* write reportType */
	  start_typ[0] = recNum;
          nstrncpy(repType, typ_ptr, wrtlen_typ - 1);

	  slen = strlen(repType);
	  if (slen < wrtlen_typ) {
            count_typ[1] = slen + 1;
          } 
          else {
            count_typ[1] = wrtlen_typ;
          }
          istatus = nc_put_vara_text(cdfid,var_typ,start_typ,count_typ,repType);
          if (istatus != NC_NOERR) {
	    free_strings(station, repType, provider);
	    return 2;
	  }
	  typ_ptr += *type_len;

/* write stationId */
	  start_sta[0] = recNum;
          nstrncpy(station, sta_ptr, wrtlen_sta - 1);

	  slen = strlen(station);
	  if (slen < wrtlen_sta) {
            count_sta[1] = slen + 1;
          } 
          else {
            count_sta[1] = wrtlen_sta;
          }
          istatus = nc_put_vara_text(cdfid,var_sta,start_sta,count_sta, station);
          if (istatus != NC_NOERR) {
	    free_strings(station, repType, provider);
	    return 2;
	  }
	  sta_ptr += *sta_len;

        }

	free_strings(station, repType, provider);
	return (SUCCESS);
}
/*******************************************************************************/
void free_strings(char *station, char *repType, char *provider) 
{
        free(station);
        free(repType);
        free(provider);
	return;
}
/*******************************************************************************/
int get_dims(int cdfid, int *file_prov, int *file_staid, int *file_reptype)
{
	int dimid, istatus;
	size_t dim_len;

        istatus = nc_inq_dimid(cdfid,"maxProviderLen",&dimid);
        if (istatus != NC_NOERR) {
          istatus = nc_close(cdfid);
          return (-1);
        }
        istatus = nc_inq_dimlen(cdfid,dimid,&dim_len);
        if (istatus != NC_NOERR){
	  printf("Cannot read dimension 'maxProviderLen'\n");
          istatus = nc_close(cdfid);
          return (-1);
        }
	*file_prov = (int)dim_len;

        istatus = nc_inq_dimid(cdfid,"maxRepTypeLen",&dimid);
        if (istatus != NC_NOERR) {
          istatus = nc_close(cdfid);
          return (-1);
        }
        istatus = nc_inq_dimlen(cdfid,dimid,&dim_len);
        if (istatus != NC_NOERR){
	  printf("Cannot read dimension 'maxRepTypeLen'\n");
          istatus = nc_close(cdfid);
          return (-1);
        }
	*file_reptype = (int)dim_len;

        istatus = nc_inq_dimid(cdfid,"maxStaIdLen",&dimid);
        if (istatus != NC_NOERR) {
          istatus = nc_close(cdfid);
          return (-1);
        }
        istatus = nc_inq_dimlen(cdfid,dimid,&dim_len);
        if (istatus != NC_NOERR){
	  printf("Cannot read dimension 'maxStaIdLen'\n");
          istatus = nc_close(cdfid);
          return (-1);
        }
	*file_staid = (int)dim_len;

	return (SUCCESS);

}

	




