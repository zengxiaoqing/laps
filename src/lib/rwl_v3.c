#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include "/usr/local/netcdf/include/netcdf.h"
#define SYSCMD "/usr/local/netcdf/bin/ncgen -o %s %s"

#define LC3_LEVELS 42
#define LM1_LEVELS 3

#ifndef DEBUG
#define DEBUG 0
#endif

#ifdef FORTRANUNDERSCORE
#define dim_size_v3 dim_size_v3_
#define check_laps_inv check_laps_inv_
#define get_lvl_coord get_lvl_coord_v3_
#define get_index_v3 get_index_v3_
#define free_read_var free_read_var_
#define retrieve_hdr_v3 retrieve_hdr_v3_
#define retrieve_grid_v3 retrieve_grid_v3_
#define cstr_to_fstr cstr_to_fstr_
#define read_cdf_v3 read_cdf_v3_
#define write_val_ref_asctime write_val_ref_asctime_
#define check_grid_dimensions check_grid_dimensions_
#define get_static_info get_static_info_
#define free_file_var free_file_var_
#define free_write_var free_write_var_
#define free_static_var free_static_var_
#define fill_c_var fill_c_var_
#define write_hdr_v3 write_hdr_v3_
#define update_inv_v3 update_inv_v3_
#define get_level_index get_level_index_
#define update_laps_v3 update_laps_v3_
#define write_cdf_v3 write_cdf_v3_
#endif
#ifdef FORTRANCAPS
#define dim_size_v3 DIM_SIZE_V3
#define check_laps_inv CHECK_LAPS_INV
#define get_lvl_coord GET_LVL_COORD_V3
#define get_index_v3 GET_INDEX_V3
#define free_read_var FREE_READ_VAR
#define retrieve_hdr_v3 RETRIEVE_HDR_V3
#define retrieve_grid_v3 RETRIEVE_GRID_V3
#define cstr_to_fstr CSTR_TO_FSTR
#define read_cdf_v3 READ_CDF_V3
#define write_val_ref_asctime WRITE_VAL_REF_ASCTIME
#define check_grid_dimensions CHECK_GRID_DIMENSIONS
#define get_static_info GET_STATIC_INFO
#define free_file_var FREE_FILE_VAR
#define free_write_var FREE_WRITE_VAR
#define free_static_var FREE_STATIC_VAR
#define fill_c_var FILL_C_VAR
#define write_hdr_v3 WRITE_HDR_V3
#define update_inv_v3 UPDATE_INV_V3
#define get_level_index GET_LEVEL_INDEX
#define update_laps_v3 UPDATE_LAPS_V3
#define write_cdf_v3 WRITE_CDF_V3
#endif
#ifdef FORTRANDOUBLEUNDERSCORE
#define dim_size_v3 dim_size_v3__
#define check_laps_inv check_laps_inv__
#define get_lvl_coord get_lvl_coord_v3__
#define get_index_v3 get_index_v3__
#define free_read_var free_read_var__
#define retrieve_hdr_v3 retrieve_hdr_v3__
#define retrieve_grid_v3 retrieve_grid_v3__
#define cstr_to_fstr cstr_to_fstr__
#define read_cdf_v3 read_cdf_v3__
#define write_val_ref_asctime write_val_ref_asctime__
#define check_grid_dimensions check_grid_dimensions__
#define get_static_info get_static_info__
#define free_file_var free_file_var__
#define free_write_var free_write_var__
#define free_static_var free_static_var__
#define fill_c_var fill_c_var__
#define write_hdr_v3 write_hdr_v3__
#define update_inv_v3 update_inv_v3__
#define get_level_index get_level_index__
#define update_laps_v3 update_laps_v3__
#define write_cdf_v3 write_cdf_v3__
#endif

/************************************************************/
#ifdef __STDC__
int dim_size_v3 (int i_cdfid, char *d_name)
#else
int dim_size_v3 (i_cdfid, d_name)
int i_cdfid;
char *d_name;
#endif
{
        int i, i_dimid, i_status, i_dsize;
        long dsize;
 
/* read in the dimension id and size */
        if ((i_dimid = ncdimid (i_cdfid, (const char *)d_name)) == (-1))
          return -1;
 
        if ((i_status = ncdiminq(i_cdfid,i_dimid,(char *)0,&dsize)) == (-1))
          return -1;
 
        i_dsize = (int) dsize;
        return i_dsize;
}
 
/************************************************************/
#ifdef __STDC__
short check_laps_inv (int i_cdfid, int i_lindx, int i_fcindx,
                      char *g_name)
#else
short check_laps_inv (i_cdfid, i_lindx, i_fcindx, g_name)
int i_cdfid, i_lindx, i_fcindx;
char *g_name;
#endif
{
        int     *i4_ptr, i, i_indx, i_status, i_varid;
        long    start[2];
        short   i_flag;
        char    dim_name[15], var_name[15];
        void    cdf_i4times();
 
/* read the fctimes inventory array associated with this grid */
        sprintf (var_name, "%s%s", g_name, "_fcinv");
 
        i_varid = ncvarid (i_cdfid, (const char *)var_name);
 
        start[0] = i_fcindx;
        start[1] = i_lindx;
 
        i_status = ncvarget1 (i_cdfid, i_varid, (const long *)start, (void *)&i_flag);
        if (i_status == -1)
          return -1;
 
        if (i_flag == 1)
          return 1;
        else
          return -1;
}

/************************************************************/
#ifdef __STDC__
float *get_lvl_coord_v3 (int i_cdfid, char *v_name, int i_size)
#else
float *get_lvl_coord_v3 (i_cdfid, v_name, i_size)
int i_cdfid;
char *v_name;
int i_size;
#endif
 
{
        int i_status, i_varid;
        long start[1], count[1];
        static float *f_ptr;
 
/* allocate memory to hold the values associated with the dimension */
        f_ptr = (float *) malloc (i_size * sizeof(float));
 
/* get the variable id for this dimension */
        i_varid = ncvarid (i_cdfid, v_name);
 
/* read the contents of the variable into memory */
        start[0] = 0;
        count[0] = (long) i_size;
        i_status = ncvarget (i_cdfid, i_varid, (const long *)start, 
                             (const long *)count, (void *) f_ptr);
        if (i_status == (-1)){
          if (DEBUG == 1)
            printf("error reading values for array index %s\n", v_name);
          return (float *)i_status;
        }
 
        return f_ptr;
}

/************************************************************/
#ifdef __STDC__
int get_index_v3 (int i_cdfid, int i_value, char *d_name)
#else
int get_index_v3 (i_cdfid, i_value, d_name)
int i_cdfid, i_value;
char *d_name;
#endif
{
        int i, i_dimid, i_status, i_varid, i_start;
        long i_dsize;
        float *f_ptr, f_value;
 
        f_value = (float) i_value;

/* read in the dimension id and size */
        i_dimid = ncdimid (i_cdfid, (const char *)d_name);
        if (i_dimid == (-1)) return -1;
 
        if ((i_status = ncdiminq(i_cdfid,i_dimid,(char *)0,&i_dsize)) == (-1))
          return -1;

/* read the contents of the coordinate variable associated with this dimension
   from the net_cdf file */
        if ((f_ptr = get_lvl_coord_v3(i_cdfid,"level",i_dsize))==(float *)(-1))
          return -1;
 
/* locate the value in the array pointed to by i_ptr */
        for (i=0; ((i<i_dsize) && (f_value!=*(f_ptr+i))); i++);
 
/* deallocate the memory */
        free (f_ptr);
 
/* test to see if the value was in the array */
        if (i<i_dsize)
          return i;
        else
          return -1;
}

/************************************************************/
#ifdef __STDC__
void free_read_var(char *var, char *comment, char *ext, 
                   char *lvl_coord, char *units, 
                   char *comm_var, char *inv_var)
#else
void free_read_var(var, comment, ext, lvl_coord, units, 
              comm_var, inv_var)
char *var;
char *comment;
char *ext;
char *lvl_coord; 
char *units; 
char *comm_var;
char *inv_var;
#endif
{
	free(var);
	free(comment);
	free(ext);
	free(lvl_coord);
	free(units);
	free(comm_var);
	free(inv_var);
}

/************************************************************/
#ifdef __STDC__
int retrieve_hdr_v3(int cdfid, long *imax, long *jmax, long *kmax)
#else
int retrieve_hdr_v3(cdfid,imax,jmax,kmax)
int cdfid;
long *imax;
long *jmax;
long *kmax;
#endif
{
        int i_status, i_varid;
        int temp, str_len, i;
        long mindex[1];
        char *t_ptr;
 
/* turn off the error handling done by the netCDF routines */
        ncopts = NC_VERBOSE;

        mindex[0] = 0;
 
/* get the data variable id */
        if ((i_varid = ncvarid (cdfid, "imax")) == (-1))
          return -1;
 
/* read the var from the netcdf file */
        i_status = ncvarget1 (cdfid, i_varid, (const long *)mindex, (void *) imax);
        if (i_status == (-1))
          return -1;
 
/* get the data variable id */
        if ((i_varid = ncvarid (cdfid, "jmax")) == (-1))
          return -1;
 
/* read the var from the netcdf file */
        i_status = ncvarget1 (cdfid, i_varid, (const long *)mindex, (void *) jmax);
        if (i_status == (-1))
          return -1;
 
/* get the data variable id */
        if ((i_varid = ncvarid (cdfid, "kmax")) == (-1))
          return -1;
 
/* read the var from the netcdf file */
        i_status = ncvarget1 (cdfid, i_varid, (const long *)mindex, (void *) kmax);
        if (i_status == (-1))
          return -1;
 
/* normal return */
 
        return 0;
}

/************************************************************/
#ifdef __STDC__
int retrieve_grid_v3(int i_cdfid,int i_level,int record_indx,
                     char *var, float *dptr, char *cptr,
                     long *comm_len, char *lvlptr, char *uptr)
#else
int retrieve_grid_v3(i_cdfid, i_level, record_indx, var,
                           dptr, cptr, comm_len, lvlptr, uptr)
int i_cdfid;
int i_level;
int record_indx;
char *var;
float *dptr;
char *cptr;
long *comm_len;
char *lvlptr;
char *uptr;
#endif
{
        int i_status, i_invflag, i_varid, lvl_indx, x_dim, y_dim;
        long start[4],count[4],start_c[3],count_c[3];
        char var_name[13];
 
/* turn off the error handling done by the netCDF routines */
        ncopts = NC_VERBOSE;
 
/* get the level index of the data array */
        lvl_indx = get_index_v3 (i_cdfid, i_level, "z");
        if (lvl_indx == (-1)){
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: no level %d\n", i_level);
          return -1;
        }
 
/* check the inventory variable to see if this grid is available */
        i_invflag = check_laps_inv (i_cdfid, lvl_indx,
                                    record_indx, var);
        if (i_invflag == (-1)){
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: grid not available\n");
          return -1;
        }
 
/* get the x and y dimension sizes */
        x_dim = dim_size_v3 (i_cdfid, "x");
        y_dim = dim_size_v3 (i_cdfid, "y");
 
/* get the data variable id */
        if ((i_varid = ncvarid (i_cdfid, (const char *)var)) == (-1)) {
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: no grid available.\n");
          return -1;
        }
 
/* construct the arrays needed to read the grid */
        start[0] = record_indx;
        start[1] = lvl_indx;
        start[2] = 0;
        start[3] = 0;
 
        count[0] = 1;
        count[1] = 1;
        count[2] = y_dim;
        count[3] = x_dim;
 
/* read the grid from the netcdf file */
        i_status = ncvarget (i_cdfid, i_varid, (const long *)start, 
                             (const long *)count, (void *) dptr);
        if (i_status == (-1)) {
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: error retrieving data %s grid.\n",
                         *var);
          return -1;
        }
 
/* get attributes lvl_coord and LAPS_units */
        i_status = ncattget (i_cdfid, i_varid, "lvl_coord", 
                             (void *) lvlptr);
        if (i_status == (-1)) {
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: error retrieving lvl_coord.\n");
          return -1;
        }
 
        i_status = ncattget (i_cdfid, i_varid, "LAPS_units", 
                             (void *) uptr);
        if (i_status == (-1)) {
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: error retrieving LAPS_units.\n");
          return -1;
        }
 
/* setup to read the comment from the netcdf file */
 
        sprintf(var_name, "%s%s", var, "_comment");
        if ((i_varid = ncvarid (i_cdfid, (const char *)var_name)) == (-1)) {
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: no comment field available.\n");
          return -1;
        }
 
/* construct the arrays needed to read the comment */
        start_c[0] = record_indx;
        start_c[1] = lvl_indx;
        start_c[2] = 0;
 
        count_c[0] = 1;
        count_c[1] = 1;
        count_c[2] = *comm_len + 1;
 
/* read the comment from the netcdf file */
        i_status = ncvarget (i_cdfid, i_varid, (const long *)start_c, 
                             (const long *)count_c, (void *) cptr);
        if (i_status == (-1)) {
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: error retrieving comment.\n");
          return -1;
        }
 
/* normal return */
 
        return 0;
}

/************************************************************/
#ifdef __STDC__
void cstr_to_fstr(char *f_comment, char *comment, long *comm_len,
                  char *f_lvl_coord, char *lvl_coord, 
                  long *lvl_coord_len, char *f_units, char *units,
                  long *units_len, long *kdim)
#else
void cstr_to_fstr(f_comment, comment, comm_len, f_lvl_coord, 
                  lvl_coord, lvl_coord_len, f_units, units,
                  units_len, kdim)
char *f_comment; 
char *comment; 
long *comm_len;
char *f_lvl_coord; 
char *lvl_coord;
long *lvl_coord_len; 
char *f_units; 
char *units;
long *units_len; 
long *kdim;
#endif
{
	char *uptr, *cptr, *lptr;
        char *hld_unit, *hld_comm, *hld_lvl;
        int i, j, slen;

        *f_units = '\0';
        *f_comment = '\0';
        *f_lvl_coord = '\0';
         
        hld_unit = malloc(*units_len + 1);
        hld_comm = malloc(*comm_len + 1);
        hld_lvl = malloc(*lvl_coord_len + 1);
        uptr = units;
        cptr = comment;
        lptr = lvl_coord;

        for (i = 0; i < *kdim; i++) {
          slen = strlen(uptr);
          if (slen < (*units_len)) {
            strcpy(hld_unit,uptr);
            for (j = slen + 1; j <= *units_len; j++)
              strcat(hld_unit," ");
          }
          else
            strncpy(hld_unit,uptr,*units_len);
          strcat(f_units, hld_unit);
          uptr += (*units_len) + 1;

          slen = strlen(cptr);
          if (slen < (*comm_len)) {
            strcpy(hld_comm,cptr);
            for (j = slen + 1; j <= *comm_len; j++)
              strcat(hld_comm," ");
          }
          else
            strncpy(hld_comm,cptr,*comm_len);
          strcat(f_comment, hld_comm);
          cptr += (*comm_len) + 1;

          slen = strlen(lptr);
          if (slen < (*lvl_coord_len)) {
            strcpy(hld_lvl,lptr);
            for (j = slen + 1; j <= *lvl_coord_len; j++)
              strcat(hld_lvl, " ");
          }
          else
            strncpy(hld_lvl,lptr,*lvl_coord_len);
          strcat(f_lvl_coord, hld_lvl);
          lptr += (*lvl_coord_len) + 1;
        }

        free(hld_lvl);
        free(hld_comm);
        free(hld_unit);
}
/************************************************************/
#ifdef __STDC__
void read_cdf_v3 (char *f_filename, char *f_ext, char *f_var, 
                  char *f_comment, char *f_lvl_coord, char *f_units, 
                  long *var_len, long *comm_len, long *fn_length, 
                  long *ext_len, long *lvl_coord_len, long *units_len, 
                  long *i_reftime, long *i_valtime, long *iimax, 
                  long *jjmax, long *kkmax, long *kdim, long lvl[], 
                  float *data, long *called_from, long *status)
#else
void read_cdf_v3 (f_filename, f_ext, f_var, f_comment, f_lvl_coord, 
                  f_units, var_len, comm_len, fn_length, ext_len, 
                  lvl_coord_len, units_len, i_reftime, i_valtime,
                  iimax, jjmax, kkmax, kdim, lvl, data, called_from, 
                  status)

char *f_filename;
char *f_ext;
char *f_var;
char *f_comment;
char *f_lvl_coord;
char *f_units;
long *var_len;
long *comm_len;
long *fn_length;
long *ext_len;
long *lvl_coord_len;
long *units_len;
long *i_reftime;
long *i_valtime;
long *iimax;
long *jjmax;
long *kkmax;
long *kdim;
long lvl[];
float *data;
long *called_from;
long *status;
#endif
{
        int cdfid, istatus;
        int i,j, t_level, i_record, unconv_var;
        int val_id,dim_id, t_record, int_1, int_2, found,t_var_id;
        long num_record, imax, jmax, kmax, mindex[1];
        float *dptr;
        double reftime, valtime, d_valtime, timeoff, diff; 
        char *filename, *cpt, char1;
	char *ext, *var, *comment, *lvl_coord, *units; 
	char *comm_var, *inv_var;
        char *fvptr, *vptr, *lptr, *uptr, *cptr, *t_var;

/* turn off the error handling done by the netCDF routines */
        ncopts = NC_VERBOSE;
        reftime = (double) *i_reftime;

/* convert fortran file_name into C fname  */
        filename = malloc(*fn_length + 1);
        nstrncpy(filename,f_filename,*fn_length);

/* open netCDF file for reading, if file exists */
        if( access(filename, F_OK) != 0 ) {
          printf("The LAPS file %s does not exist\n",filename);
          *status = -1; /* error opening file */
          return;
        }

        cdfid = ncopen((const char*)filename,NC_NOWRITE);
        free(filename);
        if (cdfid == -1) {
          *status = -1; /* error opening file */
          return;
        }

/* set i_record to value returned by ncdiminq */

        if ((dim_id = ncdimid(cdfid,"record")) == (-1)) {
          ncclose(cdfid);
          *status = -2; /* error in reading header info */
          return;
        }
        
        istatus = ncdiminq(cdfid, dim_id, (char *) 0, &num_record);
        if ((istatus == -1) || (num_record == 0)) {  /* no data in file */
          ncclose(cdfid);
          *status = -2; /* error in reading header info */
          return;
        }
        else {
          if (num_record >= 1)
            i_record = (int) (num_record - 1);
        }

/* allocate space for c strings */
	var = malloc((*var_len + 1) * sizeof(char) * (*kdim));
        comment = malloc((*comm_len + 1) * sizeof(char) * (*kdim));
	ext = malloc((*ext_len + 1) * sizeof(char));
        lvl_coord = malloc((*lvl_coord_len + 1) * sizeof(char) * (*kdim));
        units = malloc((*units_len + 1) * sizeof(char) * (*kdim));
        
/* allocate space for comment variable: ""var"_comment\0" = var_len + 9 */
        comm_var = malloc((*var_len + 9) * sizeof(char));

/* allocate space for inventory variable: ""var"_fcinv\0" = var_len + 7  */
        inv_var = malloc((*var_len + 7) * sizeof(char));
 
/* null out arrays for lvl_coord,units and comment before using; 
   tfr values from f_var to var   */
 
        fvptr = f_var;
        vptr = var;
        lptr = lvl_coord;
        uptr = units;
        cptr = comment;

        for (i = 0; i < *kdim; i++) {
          nstrncpy(vptr,fvptr,*var_len);
          downcase_c(vptr,vptr);
          vptr += (*var_len);
          *vptr = '\0';
          vptr++;
          fvptr += (*var_len);
          *lptr = '\0';
          lptr += (*lvl_coord_len) + 1;
          *uptr = '\0';
          uptr += (*units_len) + 1;
          *cptr = '\0';
          cptr += (*comm_len) + 1;
        }
        nstrncpy(ext,f_ext,*ext_len);
        upcase_c(ext, ext);

        istatus = retrieve_hdr_v3(cdfid,&imax,&jmax,&kmax);
        if (istatus == -1) {
          ncclose(cdfid);
          free_read_var(var, comment, ext, lvl_coord, units,
                        comm_var, inv_var);
          *status = -2; /* error in reading header info */
          return;
        }

        if (imax > *iimax || jmax > *jjmax || *kkmax > *kdim) {
          ncclose(cdfid);
          free_read_var(var, comment, ext, lvl_coord, units,
                        comm_var, inv_var);
          *status = -3; /* error in grid_size */
          return;
        }
 
        if ((strcmp(ext,"LMR") == 0) || (strcmp(ext,"LF1") == 0)) {

          if ((val_id = ncvarid(cdfid,"valtime")) == (-1)) {
            ncclose(cdfid);
            free_read_var(var, comment, ext, lvl_coord, units,
                          comm_var, inv_var);
            *status = -2; /* error in reading header info */
            return;
          }
        }

        unconv_var = 0;
        vptr = var;
        dptr = data;
        if (*called_from == 0) {
          cptr = comment;
          lptr = lvl_coord;
          uptr = units;
        }
        else {
          cptr = f_comment;
          lptr = f_lvl_coord;
          uptr = f_units;
        }
        for (i = 0; i < *kkmax; i++) {
          t_level = lvl[i];

/* handle special case of i_record for lmr and lf1 files */
          if ((strcmp(ext,"LMR") == 0) || (strcmp(ext,"LF1") == 0)) {
            fvptr = vptr;  /* fvptr temp point to current var */
            fvptr++;
            char1 = *fvptr;
            int_1 = atoi((const char *)&char1) - atoi((const char *)'0');
            fvptr++;
            char1 = *fvptr;
            int_2 = atoi((const char *)&char1) - atoi((const char *)'0');
            fvptr = vptr;
            fvptr++;
            *fvptr = '\0';
            t_var = vptr;

            timeoff = (double) (((int_1*10) + int_2) * 10 * 60);
            diff = (double) (*i_valtime - *i_reftime);
            valtime = reftime + timeoff + diff;
            i = 0;
            found = 0;
            while ((i < num_record) && (found == 0)) {
              mindex[0] = i;
              istatus = ncvarget1(cdfid,val_id,(const long *)mindex,(void *)&d_valtime); 
              if (istatus == (-1)) {
                ncclose(cdfid);
                free_read_var(var, comment, ext, lvl_coord, units,
                              comm_var, inv_var);
                *status = -2; /* error in reading header info */
                return;
              }

              if (d_valtime == valtime)	{
                found = 1;
                t_record = i;
              }
              else
                i++;
            }
          }
          else {
            t_record = i_record;
            t_var = vptr;
          }


          istatus = retrieve_grid_v3(cdfid, t_level, t_record, 
                                     t_var, dptr, cptr, 
                                     comm_len, lptr, uptr); 

          if (istatus == -1)  {
            unconv_var += 1;
          }
          dptr += (*iimax)*(*jjmax);
          vptr += (*var_len) + 1;
          cptr += (*comm_len) + 1;
          uptr += (*units_len) + 1;
          lptr += (*lvl_coord_len) + 1;
 
        }
 
        ncclose(cdfid);

        if (*called_from == 0) {
          cstr_to_fstr(f_comment, comment, comm_len, f_lvl_coord, 
                       lvl_coord, lvl_coord_len, f_units, units,
                       units_len, kdim);
        }
 
        free_read_var(var, comment, ext, lvl_coord, units,
                      comm_var, inv_var);
        *status = (long) unconv_var;
        return;

}
/************************************************************/
#ifdef __STDC__
int write_val_ref_asctime(int cdfid, int i_record, double *valtime,
                          double *reftime, char *asctime,
                          long *asc_len)
#else
int write_val_ref_asctime(cdfid, i_record, valtime, reftime, 
                          asctime, asc_len)
int cdfid; 
int i_record;
double *valtime; 
double *reftime;
char *asctime;
long *asc_len;
#endif
{
        int varid, istatus;
        long mindex[1], start_2[2], count_2[2];

        mindex[0] = i_record;
        if ((varid = ncvarid(cdfid,"valtime")) == (-1)) {
          return -1;
        }
        else {
        if (DEBUG==1) printf("found varid for valtime\n");
        }

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)valtime);
        if (istatus == -1) {
          return -1;
        }
        else {
        if (DEBUG==1) printf("correctly wrote valtime\n");
        }
            
        if ((varid = ncvarid(cdfid,"reftime")) == (-1)) {
          return -1;
        }
        else {
        if (DEBUG==1) printf("found varid for reftime\n");
        }

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)reftime);
        if (istatus == -1) {
          return -1;
        }
        else {
        if (DEBUG==1) printf("correctly wrote reftime\n");
        }
            
/* write out asctime */
        if ((varid = ncvarid(cdfid,"asctime")) == (-1)) {
          return (-1);
        }
        else {
        if (DEBUG==1) printf("found varid for asctime\n");
        }

	start_2[0] = i_record;
        start_2[1] = 0;
        count_2[0] = 1;
        count_2[1] = *asc_len;
        if (DEBUG == 1) { 
          printf("start = [%d][%d] | count = [%d][%d] \n", start_2[0],start_2[1],count_2[0],count_2[1]);
          printf("asctime = [%s]\n",asctime);
        }

        istatus = ncvarput(cdfid, varid, (const long *)start_2,
                           (const long *)count_2, (void *)asctime);
        if (istatus == (-1)) {
          return (-1);
        }
        else {
        if (DEBUG==1) printf("correctly wrote asctime\n");
        }

/* normal return */

        return(0);

}
/************************************************************/
#ifdef __STDC__
int check_grid_dimensions(int cdfid, char *ext, long *imax, 
                          long *jmax, long *n_levels, long lc3_levels,
                          long lm1_levels)
#else
int check_grid_dimensions(cdfid, ext, imax, jmax, n_levels, lc3_levels,
                          lm1_levels)
int cdfid;
char *ext;
long *imax; 
long *jmax; 
long *n_levels;
long lc3_levels;
long lm1_levels;
#endif
{
        int dimid, istatus;
        long dim_val;

/* determine value in dimension x */
        if ((dimid = ncdimid(cdfid,"x")) == (-1)) {
          printf("No x dimension found in output file.\n");
          return -1;
        }
        if ((istatus = ncdiminq(cdfid,dimid,(char *)0,(long *)&dim_val)) == (-1)) {
          printf("Unable to access x dimension in output file.\n");
          return -1;
        }

        if (*imax == dim_val) {
        }
        else {
          printf("x dimension in output file does not match IMAX passed in.\n");
          return -1;
        }

/* determine value in dimension y */
        if ((dimid = ncdimid(cdfid,"y")) == (-1)) {
          printf("No y dimension found in output file.\n");
          return -1;
        }
        if ((istatus = ncdiminq(cdfid,dimid,(char *)0,(long *)&dim_val)) == (-1)) {
          printf("Unable to access y dimension in output file.\n");
          return -1;
        }

        if (*jmax == dim_val) {
        }
        else {
          printf("y dimension in output file does not match JMAX passed in.\n");
          return -1;
        }

/* determine value in dimension z */
        if ((dimid = ncdimid(cdfid,"z")) == (-1)) {
          printf("No z dimension found in output file.\n");
          return -1;
        }
        if ((istatus = ncdiminq(cdfid,dimid,(char *)0,(long *)&dim_val)) == (-1)) {
          printf("Unable to access z dimension in output file.\n");
          return -1;
        }

        if (dim_val == *n_levels) {
          /* file is LAPS standard 3D file */
        }
        else if (dim_val == 1) {
          *n_levels = 1;
        }
        else if ((strncmp(ext,"lc3",3) == 0) && (dim_val == (long) LC3_LEVELS)) {
          *n_levels = (long) LC3_LEVELS;
        }
        else if ((strncmp(ext,"lm1",3) == 0) && (dim_val == (long) LM1_LEVELS)) {
          *n_levels = (long) LM1_LEVELS;
        }
        else {
          printf("z dimension in output file does not match n_levels passed in.\n");
          return -1;
        }
    
/* normal return */

        return(0);  /* dimensions check out OK */
}
/************************************************************/
#ifdef __STDC__
int get_static_info(char *static_grid, float *Dx, float *Dy, 
                    float *La1, float *Lo1, float *LoV, 
                    float *Latin1, float *Latin2, char *map_proj, 
                    char *origin)
#else
int get_static_info(static_grid, Dx, Dy, La1, Lo1, LoV, Latin1,
                    Latin2, map_proj, origin)
char *static_grid; 
float *Dx; 
float *Dy; 
float *La1; 
float *Lo1; 
float *LoV; 
float *Latin1; 
float *Latin2; 
char *map_proj; 
char *origin;
#endif
{
        int cdfid_stat, dimid, varid, istatus, ret_status;
        long namelen, mindex[1], start[1], count[1];
        long start_g[2], count_g[2];

         
	cdfid_stat = ncopen(static_grid,NC_NOWRITE);
        if (cdfid_stat == -1) {
          printf("Unable to open static file %s.\n",static_grid);
          return -1;
        }

        ret_status = 0;

/* determine value in dimension namelen */
        if ((dimid = ncdimid(cdfid_stat,"namelen")) == (-1)) {
          namelen = 132;
          printf("Variables map_proj and origin may be incorrect.\n");
          printf("Unable to read dimension namelen from static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncdiminq(cdfid_stat,dimid,(char *)0,(long *)&namelen);
          if (istatus == -1) {
            namelen = 132;
            printf("Variables map_proj and origin may be incorrect.\n");
            printf("Unable to read dimension namelen from static.nest7grid.\n");
            ret_status = -1;
          }
        }

        mindex[0] = 0;

/* read Dx */
        if ((varid = ncvarid(cdfid_stat,"Dx")) == (-1)) {
          printf("Error reading Dx from  static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncvarget1 (cdfid_stat, varid, (const long *)mindex, (void *)Dx);
          if (istatus == -1) {
            printf("Error reading Dx from  static.nest7grid.\n");
            ret_status = -1;
          }
        }
        
/* read Dy */
        if ((varid = ncvarid(cdfid_stat,"Dy")) == (-1)) {
          printf("Error reading Dy from  static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncvarget1 (cdfid_stat, varid, (const long *)mindex, (void *)Dy);
          if (istatus == -1) {
            printf("Error reading Dy from  static.nest7grid.\n");
            ret_status = -1;
          }
        }
        
/* read La1 */
        if ((varid = ncvarid(cdfid_stat,"La1")) == (-1)) {
          printf("Error reading La1 from  static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncvarget1 (cdfid_stat, varid, (const long *)mindex, (void *)La1);
          if (istatus == -1) {
            printf("Error reading La1 from  static.nest7grid.\n");
            ret_status = -1;
          }
        }
        
/* read Lo1 */
        if ((varid = ncvarid(cdfid_stat,"Lo1")) == (-1)) {
          printf("Error reading Lo1 from  static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncvarget1 (cdfid_stat, varid, (const long *)mindex, (void *)Lo1);
          if (istatus == -1) {
            printf("Error reading Lo1 from  static.nest7grid.\n");
            ret_status = -1;
          }
        }
        
/* read LoV */
        if ((varid = ncvarid(cdfid_stat,"LoV")) == (-1)) {
          printf("Error reading LoV from  static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncvarget1 (cdfid_stat, varid, (const long *)mindex, (void *)LoV);
          if (istatus == -1) {
            printf("Error reading LoV from  static.nest7grid.\n");
            ret_status = -1;
          }
        }
        
/* read Latin1 */
        if ((varid = ncvarid(cdfid_stat,"Latin1")) == (-1)) {
          printf("Error reading Latin1 from  static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncvarget1 (cdfid_stat, varid, (const long *)mindex, (void *)Latin1);
          if (istatus == -1) {
            printf("Error reading Latin1 from  static.nest7grid.\n");
            ret_status = -1;
          }
        }
        
/* read Latin2 */
        if ((varid = ncvarid(cdfid_stat,"Latin2")) == (-1)) {
          printf("Error reading Latin2 from  static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncvarget1 (cdfid_stat, varid, (const long *)mindex, (void *)Latin2);
          if (istatus == -1) {
            printf("Error reading Latin2 from  static.nest7grid.\n");
            ret_status = -1;
          }
        }
        
        start_g[0] = 0;
        start_g[1] = 0;
        count_g[0] = 1;
        count_g[1] = namelen;

/* read grid_type */
        if ((varid = ncvarid(cdfid_stat,"grid_type")) == (-1)) {
          printf("Error reading grid_type from static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncvarget (cdfid_stat, varid, (const long *)start_g, 
                              (const long *)count_g, (void *)map_proj);
          if (istatus == -1) {
            printf("Error reading grid_type from  static.nest7grid.\n");
            ret_status = -1;
          }
        }

        start[0] = 0;
        count[0] = namelen;
        
/* read origin_name */
        if ((varid = ncvarid(cdfid_stat,"origin_name")) == (-1)) {
          printf("Error reading origin_name from  static.nest7grid.\n");
          ret_status = -1;
        }
        else {
          istatus = ncvarget (cdfid_stat, varid, (const long *)start, 
                              (const long *)count, (void *)origin);
          if (istatus == -1) {
            printf("Error reading origin_name from  static.nest7grid.\n");
            ret_status = -1;
          }
        }
        
        ncclose(cdfid_stat);
        if(ret_status == -1) {
          return(ret_status);
        }
        else {
          return(0);
        }
}
/************************************************************/
#ifdef __STDC__
void free_file_var(char *syscmd, char *filename)
#else
void free_file_var(syscmd, filename)

char *syscmd;
char *filename; 
#endif
{
          free(filename);
          free(syscmd);
          return;
}
/************************************************************/
#ifdef __STDC__
void free_write_var(char *var, char *comment, char *asctime, 
                    char *comm_var, char *inv_var)
#else
void free_write_var(var, comment, asctime, comm_var, inv_var)
char *var;
char *comment;
char *asctime;
char *comm_var;
char *inv_var;
#endif
{
	free(var);
	free(comment);
        free(asctime);
	free(comm_var);
	free(inv_var);
}

/************************************************************/
#ifdef __STDC__
void free_static_var(char *map_proj, char *origin)
#else
void free_static_var(map_proj, origin)
char *map_proj;
char *origin;
#endif
{
	free(map_proj);
	free(origin);
}

/************************************************************/
#ifdef __STDC__
void fill_c_var(long *kdim, char *f_var, char *var, long *var_len,
               char *f_comment, char *comment, long *comm_len,
               char *f_asctime, char *asctime, long *asc_len)
#else
void fill_c_var(kdim,f_var,var,var_len,f_comment,comment,comm_len,
                f_asctime, asctime,asc_len)
long *kdim;
char *f_var;
char *var;
long *var_len;
char *f_comment;
char *comment;
long *comm_len;
char *f_asctime;
char *asctime;
long *asc_len;
#endif
{
	int i;
        char *fvptr, *vptr, *fcptr, *cptr;

        fvptr = f_var;
        vptr = var;
        fcptr = f_comment;
        cptr = comment;

        for (i = 0; i < *kdim; i++) {
          nstrncpy(vptr,fvptr,*var_len);
          downcase_c(vptr,vptr);
          vptr += (*var_len) + 1;
          fvptr += (*var_len);
          fstrncpy(cptr,fcptr,*comm_len);
          cptr += (*comm_len) + 1;
          fcptr += (*comm_len);
        }

        fstrncpy(asctime,f_asctime,*asc_len);
}
/************************************************************/
#ifdef __STDC__
int  write_hdr_v3(int cdfid, char *ext, int i_record, long *imax, 
                  long *jmax, long *kmax, long *kdim, float *base, 
                  float *interval, long *n_levels, float *Dx, 
                  float *Dy, float *La1, float *Lo1, float *LoV, 
                  float *Latin1, float *Latin2, char *map_proj, 
                  char *origin) 
#else
int  write_hdr_v3(cdfid, ext, i_record, imax, jmax, kmax, kdim, base, 
                  interval,n_levels, Dx, Dy, La1, Lo1, LoV, Latin1, 
                  Latin2, map_proj, origin) 
int cdfid; 
char *ext;
int i_record; 
long *imax; 
long *jmax; 
long *kmax; 
long *kdim; 
float *base; 
float *interval; 
long *n_levels; 
float *Dx; 
float *Dy; 
float *La1; 
float *Lo1; 
float *LoV; 
float *Latin1; 
float *Latin2;
char *map_proj; 
char *origin;
#endif
{
        short Nx, Ny;
	int i, j, istatus, varid;
        long start[1], count[1], mindex[1], start_2[2], count_2[2];
        float level_val, *levels, *lp;
        
        istatus = 0;
        mindex[0] = 0;

/* write out imax */
        if ((varid = ncvarid(cdfid,"imax")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for imax\n");
       

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)imax);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote imax\n");

/* write out jmax */
        if ((varid = ncvarid(cdfid,"jmax")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for jmax\n");

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)jmax);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote jmax\n");

/* write out kmax */
        if ((varid = ncvarid(cdfid,"kmax")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for kmax\n");

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)kmax);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote kmax\n");

/* write out kdim */
        if ((varid = ncvarid(cdfid,"kdim")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for kdim\n");
        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)kdim);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote kdim\n");

/* write out levels */
/* levels written are based on value in n_levels, which has been adjusted to
     one of four options: standard LAPS 3D (nk_laps in nest7grid.parms),
                          surface = 1, or two special files LC3 and LM1,
                          which are currently set with #define at beginning of
                          this file */

        levels = (float *) malloc((*n_levels) * sizeof(float));
        lp = levels;

        if (strncmp(ext,"lc3",3) == 0) {
          for( i = 1; i <= LC3_LEVELS; i++) {
            *lp = i;
            lp++;
          }
        }
        else if (strncmp(ext,"lm1",3) == 0) {
          for( i = 1; i <= LM1_LEVELS; i++) {
            *lp = i * (-1);
            lp++;
          }
        }
        else if (*n_levels == 1) {
          *levels = 0;
        }
        else {  /* standard LAPS 3D */

          j = *n_levels - 1;
          for( i = 0; i < *n_levels; i++) {
            lp = levels + j;
            j = j - 1;
            *lp = *base;
            *base = *base - *interval;
          }
      
/*        for( i = 0; i < *n_levels; i++) {
            *lp = *base;
            lp++;
            *base = *base - *interval;
          }
*/
        }

        start[0] = 0;
        count[0] = *n_levels;

        if ((varid = ncvarid(cdfid,"level")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for level\n");

        istatus = ncvarput(cdfid,varid,(const long *)start,
                           (const long *) count,(void *)levels);
          if (istatus == (-1))
            return (-1);
        if (DEBUG==1) printf("correctly wrote level\n");

        free (levels);

        mindex[0] = 0;

/* write out Nx */
        if ((varid = ncvarid(cdfid,"Nx")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for Nx\n");

        Nx = (short) *imax;
        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)&Nx);
        if (istatus == (-1)) 
          return (-1);
        if (DEBUG==1) printf("correctly wrote Nx\n");

/* write out Ny */
        if ((varid = ncvarid(cdfid,"Ny")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for Ny\n");

        Ny = (short) *jmax;
        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)&Ny);
        if (istatus == (-1)) 
          return (-1);
        if (DEBUG==1) printf("correctly wrote Ny\n");

/* write out Dx */
        if ((varid = ncvarid(cdfid,"Dx")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for Dx\n");

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)Dx);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote Dx\n");

/* write out Dy */
        if ((varid = ncvarid(cdfid,"Dy")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for Dy\n");

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)Dy);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote Dy\n");

/* write out La1 */
        if ((varid = ncvarid(cdfid,"La1")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for La1\n");

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)La1);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote La1\n");

/* write out Lo1 */
        if ((varid = ncvarid(cdfid,"Lo1")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for Lo1\n");

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)Lo1);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote Lo1\n");

/* write out LoV */
        if ((varid = ncvarid(cdfid,"LoV")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for LoV\n");

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)LoV);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote LoV\n");

/* write out Latin1 */
        if ((varid = ncvarid(cdfid,"Latin1")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for Latin1\n");

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)Latin1);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote Latin1\n");

/* write out Latin2 */
        if ((varid = ncvarid(cdfid,"Latin2")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for Latin2\n");

        istatus = ncvarput1(cdfid,varid,(const long *)mindex,(void *)Latin2);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote Latin2\n");

/* write out map_proj */
        if ((varid = ncvarid(cdfid,"grid_type")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for grid_type\n");

	start_2[0] = i_record;
        start_2[1] = 0;
        count_2[0] = 1;
        count_2[1] = strlen(map_proj);

        istatus = ncvarput(cdfid,varid,(const long *)start_2,
                           (const long *)count_2,(void *)map_proj);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote map_proj\n");

/* write out origin */
        if ((varid = ncvarid(cdfid,"origin_name")) == (-1))
          return (-1);
        if (DEBUG==1) printf("found varid for origin_name\n");

        start[0] = 0;
        count[0] = strlen(origin);
        istatus = ncvarput(cdfid,varid,(const long *)start,
                           (const long *)count,(void *)origin);
        if (istatus == (-1))
          return (-1);
        if (DEBUG==1) printf("correctly wrote origin\n");


        if (istatus == (-1))
          return (-1);
        else
          return 0;
}
/************************************************************/
#ifdef __STDC__
int update_inv_v3(int cdfid,int i_record,int z_dim,int inv_id)
#else
int update_inv_v3(cdfid, i_record, z_dim, inv_id)
int cdfid;
int i_record;
int z_dim;
int inv_id;
#endif
{
	int istatus;
        long start[2];
        short i_flag;

	start[0] = i_record;
        start[1] = z_dim;

        i_flag = 1;
        istatus = ncvarput1(cdfid,inv_id,(const long *)start,
                            (void *)&i_flag);
        if (istatus == (-1)) 
          return istatus;
        else
          return 0;
}
/************************************************************/
#ifdef __STDC__
int get_level_index(int cdfid, char *var, char *dim,
                    float match_level)
#else
int get_level_index(cdfid,var,dim,match_level)
int cdfid; 
char *var; 
char *dim;
float match_level;
#endif
{
        int dim_id, var_id, i, istatus, found;
        long dim_val,start[1],count[1];
        float *file_lvls, *flvls;

/* determine value in dimension z */
        if ((dim_id = ncdimid(cdfid,dim)) == (-1)) return -1;
        if ((istatus = ncdiminq(cdfid,dim_id,(char *)0,(long *)&dim_val)) == (-1))
          return -1;

/* allocate space for dim_val floats, and pull levels */
        file_lvls = (float *) malloc(dim_val * sizeof(float));
        if ((var_id = ncvarid(cdfid,"level")) == (-1))
          return -1;

        start[0] = 0;
        count[0] = dim_val;
        istatus = ncvarget(cdfid,var_id,(const long *)start,
                           (const long *)count,(void *)file_lvls);
        if (istatus == (-1)) return -1;

        found = 0;
        i = 0;
        flvls = file_lvls;
        while ((found == 0) && (i < dim_val)) {
          if (*flvls == match_level) 
            found = 1;
          else {
            flvls++;
            i++;
          }
        }
        free(file_lvls);
 
        if (found == 1) 
          return i;
        else
          return -1;
        
}
/************************************************************/
#ifdef __STDC__
int update_laps_v3(int cdfid,long i_level,int i_record,
                   long *imax, long *jmax, int var_id, 
                   int inv_id,float *dptr,
                   int comm_id, char *cptr)
#else
int update_laps_v3(cdfid,i_level,i_record,imax,jmax, 
                   var_id,inv_id,dptr,comm_id,cptr)
int cdfid;
long i_level;
int i_record;
long *imax;
long *jmax;
int var_id;
int inv_id;
float *dptr;
int comm_id;
char *cptr;
#endif
{
	int x_dim, y_dim, z_dim, istatus, dim_id;
        long start[4], count[4], c_start[3], c_count[3];

        x_dim = *imax;
        y_dim = *jmax;
        
/* try to reduce overhead here */
        z_dim = get_level_index(cdfid,"level","z",(float)i_level);
        if (z_dim == (-1))
          return (-1);

	start[0] = i_record;
        start[1] = z_dim;
        start[2] = 0;
        start[3] = 0;


        count[0] = 1;
        count[1] = 1;
        count[2] = y_dim;
	count[3] = x_dim;

        istatus = ncvarput(cdfid,var_id,(const long *)start,
                           (const long *)count,(void *)dptr);
        if (istatus == (-1)) /* error writing out grid */
          return (-1);

        else {

	  istatus = update_inv_v3(cdfid,i_record,z_dim,inv_id);
          if (istatus == (-1)) /* error writing out grid */
            return (-1);

          c_start[0] = i_record;
          c_start[1] = z_dim;
          c_start[2] = 0;
 
          c_count[0] = 1;
          c_count[1] = 1;
          c_count[2] = strlen(cptr);

          istatus = ncvarput(cdfid,comm_id,(const long *)c_start,
                             (const long *)c_count,(void *)cptr);

          if (istatus == (-1)) /* error writing out grid */
            return (-1);
        }

/* normal return */

        return 0;

}
/************************************************************/
#ifdef __STDC__
void write_cdf_v3 (char *f_filename, char *f_ext, char *f_var, 
                   char *f_comment, char *f_asctime, char *f_cdl_path,
                   char *f_static_path, long *fn_length, long *ext_len, 
                   long *var_len, long *comm_len, long *asc_len, 
                   long *cdl_path_len, long *stat_len, long *i_reftime, 
                   long *i_valtime, long *imax, long *jmax, long *kmax, 
                   long *kdim, long lvl[], float *data, float *base, 
                   float *interval, long *n_levels, long *called_from, 
                   long *append, long *status)
#else
void write_cdf_v3 (f_filename, f_ext, f_var, f_comment, f_asctime, 
                   f_cdl_path, f_static_path, fn_length, ext_len, 
                   var_len, comm_len,  asc_len, cdl_path_len, stat_len, 
                   i_reftime, i_valtime, imax, jmax, kmax, kdim, lvl, 
                   data, base, interval, n_levels, called_from, 
                   append, status)
char *f_filename; 
char *f_ext; 
char *f_var; 
char *f_comment;
char *f_asctime;
char *f_cdl_path;
char *f_static_path;
long *fn_length; 
long *ext_len; 
long *var_len; 
long *comm_len; 
long *asc_len;
long *cdl_path_len;
long *stat_len;
long *i_reftime; 
long *i_valtime; 
long *imax;
long *jmax; 
long *kmax; 
long *kdim; 
long lvl[]; 
float *data; 
float *base; 
float *interval; 
short *n_levels;
long *called_from;
long *append;
long *status; 
#endif
{
	int i, j, istatus, cdl_len, cdfid, i_record, var_id;
        int dim_id, comm_id, old_record, int1, int2;
        int missing_grids, inv_id, lc3_levels, lm1_levels;
        long i_level, num_record, mindex[1],i4time;
        float *dptr, Dx, Dy, La1, Lo1, LoV, Latin1, Latin2;
        double timeoff, reftime, valtime, diff;
        static char *fvptr, *vptr, *cptr;
        static char *var, *comment, *filename, *ext; 
        char *map_proj, *origin;
        static char *asctime, *syscmd, *cdlfile, *comm_var, *inv_var;
        static char *static_grid;
        char char1;

/* turn off the error handling done by the netCDF routines */
        ncopts = NC_VERBOSE;

        lc3_levels = (long) LC3_LEVELS;
        lm1_levels = (long) LM1_LEVELS;

/* convert fortran f_file_name into C string filename and f_ext into C ext  */
        filename = malloc(*fn_length + 1);
        nstrncpy(filename,f_filename,*fn_length);
	ext = malloc((*ext_len + 1) * sizeof(char));
        nstrncpy(ext,f_ext,*ext_len);
        downcase_c(ext,ext);

/* allocate space for syscmd and cdlfile and fill up */
        /* cdl file contains "extension" + ".cdl\0" */
	cdlfile = malloc((*cdl_path_len + (*ext_len) + 5) * sizeof(char));
        nstrncpy(cdlfile,f_cdl_path,*cdl_path_len);
        strcat(cdlfile,ext);
        strcat(cdlfile,".cdl");
        cdl_len = strlen(cdlfile);

        /* SYSCMD contains "/usr/local/netcdf/bin/ncgen -o %s %s\0" which
           is 33 char, cdlfile, and filename  + 10 extra  */
        syscmd = malloc((strlen(SYSCMD)+cdl_len+*fn_length+10) * sizeof(char));
        sprintf(syscmd,SYSCMD, filename, cdlfile);
        free(cdlfile);
        
/*  For simplicity of error msgs, since we are not allowing append of files,
    the following logic is commented out  LW 9-25-97 
    see if file is already there by trying to open 
        printf("Checking to see if file %s exists:\n",filename);
        cdfid = ncopen(filename,NC_WRITE);
        if (cdfid == -1) {   file not there, create one  
          printf("File %s does not exist....creating new one.\n",filename);
*/
/*  create file */
          system(syscmd);
          cdfid = ncopen(filename,NC_WRITE);
          if (cdfid == -1) {
	    *status = -2; /* error in file creation */ 
            free_file_var(syscmd, filename);
            free(ext);
            return;
          }
/*      }
*/

/* set i_record to value returned by ncdiminq */

        if ((dim_id = ncdimid(cdfid,"record")) == (-1)) {
	  *status = -2; /* error in file creation */ 
          ncclose(cdfid);
          free_file_var(syscmd, filename);
          free(ext);
          return;
        }

        if ((istatus = ncdiminq(cdfid,dim_id,(char *)0,(long *)&num_record)) == (-1)) {
	  *status = -2; /* error in file creation */ 
          ncclose(cdfid);
          free_file_var(syscmd, filename);
          free(ext);
          return;
        }

        if (num_record == 0)  {   /* no data in file */
          i_record = (int)num_record;
        }
        else {
          if (*append == 0) {  /* only one analysis allowed in a file */
            *status = -6;
            ncclose(cdfid);
            free(ext);
            free_file_var(syscmd, filename);
            return; 
          }
        }

/* file is now open and ready for writing.    */
/* verify x and y in output file match imax and jmax */
/* verify that z is correct: 
     n_levels is the number of 3D levels listed in nest7grid.parms.
     options  for z are: n_levels (standard LAPS 3D file), z = 1 for surface,
                         and 2 special cases, LC3 = 42, LM1 = 3.
     set n_levels to match one of the four cases above, or abort   */

        istatus = check_grid_dimensions(cdfid, ext, imax, jmax, n_levels,
                                        lc3_levels, lm1_levels);
        if (istatus == -1) {
          *status = -3;
          free_file_var(syscmd, filename);
          free(ext);
          return;
        }

/* get info from static.nest7grid for writing to output files */
	static_grid = malloc((*stat_len + 20) * sizeof(char));
        nstrncpy(static_grid,f_static_path,*stat_len);
        strcat(static_grid,"static.nest7grid");

        map_proj = malloc(256 * sizeof(char));
        origin = malloc(256 * sizeof(char));

        istatus = get_static_info(static_grid, &Dx, &Dy, &La1, &Lo1, &LoV, 
                                  &Latin1, &Latin2, map_proj, origin);
        free(static_grid);

        if (istatus == -1) {
          printf("Error reading info from static.nest7grid.\n");
          printf("  Some navigation data will be missing from file.\n");
        }

/* convert i_reftime to reftime */
	reftime = (double)(*i_reftime);
        /* default case for analysis: valtime == reftime */
        valtime = (double)(*i_valtime);

/* transfer fortran strings to c strings */
	var = malloc((*var_len + 1) * sizeof(char) * (*kdim));
        comment = malloc((*comm_len + 1) * sizeof(char) * (*kdim));
        asctime = malloc((*asc_len + 1) * sizeof(char));

/* allocate space for comment variable: ""var"_comment\0" = var_len + 9 */
        comm_var = malloc((*var_len + 9) * sizeof(char));

/* allocate space for inventory variable: ""var"_fcinv\0" = var_len + 7 */
        inv_var = malloc((*var_len + 7) * sizeof(char));

        fill_c_var(kdim,f_var,var,var_len,f_comment,comment,comm_len,
                   f_asctime,asctime,asc_len);

/* write out header information into file */

	istatus = write_hdr_v3(cdfid, ext, i_record, imax, jmax, kmax, kdim, 
                               base, interval, n_levels, &Dx, &Dy, &La1, &Lo1, 
                               &LoV, &Latin1, &Latin2, map_proj, origin);

        free_static_var(map_proj, origin);
        if (istatus == -1) {
	  *status = -5; /* error writing out header */ 
          ncclose(cdfid);
          free_file_var(syscmd, filename);
          free(ext);
          free_write_var(var, comment, asctime, comm_var, inv_var);
          return;
        }
       
        old_record = -1;
	missing_grids = 0;
        vptr = var;

	for (i = 0; i < *kdim; i++) {

/* handle special case of i_record for lmr and lf1 files */
          if ((strcmp(ext,"lmr") == 0) || (strcmp(ext,"lf1") == 0)) {
            fvptr = vptr;  /* fvptr temp point to current var */
            fvptr++;
            char1 = *fvptr;
            int1 = atoi((const char *)&char1) - atoi((const char *)'0');
            fvptr++;
            char1 = *fvptr;
            int2 = atoi((const char *)&char1) - atoi((const char *)'0');
            fvptr = vptr;
            fvptr++;
            *fvptr = '\0';
            
            if (old_record == (-1)) diff = valtime - reftime;

            timeoff = (double) (((int1*10) + int2) * 10 * 60);
            valtime = reftime + timeoff + diff;

            i4time = valtime + 315619200;
            cv_i4tim_asc_lp(&i4time,f_asctime,&istatus);
            fstrncpy(asctime, f_asctime,*asc_len);
          }

/* write out valtime and reftime if old_record != i_record */
          if (i_record != old_record) {
            istatus = write_val_ref_asctime(cdfid, i_record, &valtime, &reftime, 
                          asctime, asc_len);
            if (istatus == -1) {
	      *status = -5; /* error writing out header */ 
              ncclose(cdfid);
              free_file_var(syscmd, filename);
              free(ext);
              free_write_var(var, comment, asctime, comm_var, inv_var);
              return;
            }
            old_record = i_record;
          }

          sprintf(comm_var,"%s_comment",vptr);
          sprintf(inv_var,"%s_fcinv",vptr);
	  var_id = ncvarid(cdfid,vptr);
          vptr += (*var_len + 1);
          comm_id = ncvarid(cdfid,comm_var);
          inv_id = ncvarid(cdfid,inv_var);
          if ((var_id == (-1)) || (comm_id == (-1)) || (inv_id == (-1))) {
          }
          else {
            i_level = lvl[i];
            dptr = (data + (i*(*imax)*(*jmax)));
            cptr = (comment + (i*(*comm_len + 1)));
            istatus = update_laps_v3(cdfid,i_level,i_record,imax,
                                     jmax,var_id,inv_id,dptr,
                                     comm_id,cptr);
            if (istatus == (-1)) missing_grids++;

            if ((strcmp(ext,"lmr") == 0) || (strcmp(ext,"lf1") == 0))
              i_record++;
          }
        }

        ncclose(cdfid);
        free(ext);
        free_write_var(var, comment, asctime, comm_var, inv_var);

        if (missing_grids < *kdim)
	  *status = missing_grids; 
        else { 
	  *status = -4; 
        }
        return;
}
/************************************************************/
