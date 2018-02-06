#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <netcdf.h>
#include <string.h>
#define SYSCMD "ncgen -o %s %s"

#define LC3_LEVELS 42
#define LM1_LEVELS 3

#ifndef DEBUG
#define DEBUG 0
#endif

#if(SIZEOF_SHORT==4)
#define fint4 short
#define ncputvar1 nc_put_var1_short
#define ncgetvar1 nc_get_var1_short
#elif(SIZEOF_INT==4)
#define fint4 int
#define ncputvar1 nc_put_var1_int
#define ncgetvar1 nc_get_var1_int
#elif(SIZEOF_LONG==4)
#define fint4 long
#define ncputvar1 nc_put_var1_long
#define ncgetvar1 nc_get_var1_long
#endif

#ifdef FORTRANUNDERSCORE
#define dim_size_v3 dim_size_v3_
#define check_laps_inv check_laps_inv_
#define get_lvl_coord_v3 get_lvl_coord_v3_
#define get_index_v3 get_index_v3_
#define free_read_var free_read_var_
#define retrieve_hdr_v3 retrieve_hdr_v3_
#define retrieve_grid_v3 retrieve_grid_v3_
#define cstr_to_fstr cstr_to_fstr_
#define read_cdf_v3 read_cdf_v3_
#define write_val_ref_asctime write_val_ref_asctime_
#define check_grid_dims_v3 check_grid_dims_v3_
#define get_static_info_v3 get_static_info_v3_
#define free_file_var free_file_var_
#define free_write_var free_write_var_
#define free_static_var free_static_var_
#define fill_c_var fill_c_var_
#define write_hdr_v3 write_hdr_v3_
#define update_inv_v3 update_inv_v3_
#define get_level_index get_level_index_
#define update_laps_v3 update_laps_v3_
#define write_cdf_v3 write_cdf_v3_
#define open_cdf open_cdf_
#define nstrncpy nstrncpy_
#define fstrncpy fstrncpy_
#define cv_i4tim_asc_lp cv_i4tim_asc_lp_
#endif
#ifdef FORTRANCAPS
#define dim_size_v3 DIM_SIZE_V3
#define check_laps_inv CHECK_LAPS_INV
#define get_lvl_coord_v3 GET_LVL_COORD_V3
#define get_index_v3 GET_INDEX_V3
#define free_read_var FREE_READ_VAR
#define retrieve_hdr_v3 RETRIEVE_HDR_V3
#define retrieve_grid_v3 RETRIEVE_GRID_V3
#define cstr_to_fstr CSTR_TO_FSTR
#define read_cdf_v3 READ_CDF_V3
#define write_val_ref_asctime WRITE_VAL_REF_ASCTIME
#define check_grid_dims_v3 CHECK_GRID_DIMS_V3
#define get_static_info_v3 GET_STATIC_INFO_V3
#define free_file_var FREE_FILE_VAR
#define free_write_var FREE_WRITE_VAR
#define free_static_var FREE_STATIC_VAR
#define fill_c_var FILL_C_VAR
#define write_hdr_v3 WRITE_HDR_V3
#define update_inv_v3 UPDATE_INV_V3
#define get_level_index GET_LEVEL_INDEX
#define update_laps_v3 UPDATE_LAPS_V3
#define write_cdf_v3 WRITE_CDF_V3
#define open_cdf OPEN_CDF
#define nstrncpy NSTRNCPY
#define fstrncpy FSTRNCPY
#define cv_i4_tim_asc_lp CV_I4_TIM_ASC_LP
#endif
#ifdef FORTRANDOUBLEUNDERSCORE
#define dim_size_v3 dim_size_v3__
#define check_laps_inv check_laps_inv__
#define get_lvl_coord_v3 get_lvl_coord_v3__
#define get_index_v3 get_index_v3__
#define free_read_var free_read_var__
#define retrieve_hdr_v3 retrieve_hdr_v3__
#define retrieve_grid_v3 retrieve_grid_v3__
#define cstr_to_fstr cstr_to_fstr__
#define read_cdf_v3 read_cdf_v3__
#define write_val_ref_asctime write_val_ref_asctime__
#define check_grid_dims_v3 check_grid_dims_v3__
#define get_static_info_v3 get_static_info_v3__
#define free_file_var free_file_var__
#define free_write_var free_write_var__
#define free_static_var free_static_var__
#define fill_c_var fill_c_var__
#define write_hdr_v3 write_hdr_v3__
#define update_inv_v3 update_inv_v3__
#define get_level_index get_level_index__
#define update_laps_v3 update_laps_v3__
#define write_cdf_v3 write_cdf_v3__
#define open_cdf open_cdf__
#define nstrncpy nstrncpy__
#define fstrncpy fstrncpy__
#define cv_i4_tim_asc_lp cv_i4tm_asc_lp__
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
        size_t dim_len;
 
/* read in the dimension id and size */
        i_status = nc_inq_dimid (i_cdfid, d_name, &i_dimid);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          return -1;
	} 
        i_status = nc_inq_dimlen(i_cdfid,i_dimid,&dim_len);
        if (i_status != NC_NOERR){
          printf("%s\n",nc_strerror(i_status));
          return -1;
	} 

        i_dsize = (int) dim_len;

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
        size_t  start[2];
        short   i_flag;
        char    dim_name[15], var_name[15];
        void    cdf_i4times();
 
/* read the fctimes inventory array associated with this grid */
        sprintf (var_name, "%s%s", g_name, "_fcinv");
 
        i_status = nc_inq_varid (i_cdfid, var_name, &i_varid);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          return -1;
	}
 
        start[0] = i_fcindx;
        start[1] = i_lindx;
 
        i_status = nc_get_var1_short (i_cdfid,i_varid,start,&i_flag);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          return -1;
	}
 
        if (i_flag == 1)
          return 1;
        else
          printf("warning in check_laps_inv i_flag != 1\n");
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
        size_t start[1], count[1];
        static float *f_ptr;
 
/* allocate memory to hold the values associated with the dimension */
        f_ptr = (float *) malloc (i_size * sizeof(float));
 
/* get the variable id for this dimension */
        i_status = nc_inq_varid (i_cdfid,v_name,&i_varid);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          if (DEBUG == 1)
            printf("error reading values for array index %s\n", v_name);
          i_status = -1;
          return (float *)i_status;
        }
 
/* read the contents of the variable into memory */
        start[0] = 0;
        count[0] = (size_t) i_size;
        i_status = nc_get_vara_float (i_cdfid,i_varid,start,count,f_ptr);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          if (DEBUG == 1)
            printf("error reading values for array index %s\n", v_name);
          i_status = -1;
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
        int i, i_dimid, i_status, i_varid, i_dsize;
        size_t dim_len;
        float *f_ptr, f_value;
 
        f_value = (float) i_value;

/* read in the dimension id and size */
        i_status = nc_inq_dimid (i_cdfid, d_name, &i_dimid);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          return -1;
	} 
        i_status = nc_inq_dimlen(i_cdfid,i_dimid,&dim_len);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          return -1;
	} 



        i_dsize = (int)dim_len;

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
int retrieve_hdr_v3(int cdfid, fint4 *imax, fint4 *jmax, fint4 *kmax)
#else
int retrieve_hdr_v3(cdfid,imax,jmax,kmax)
int cdfid;
fint4 *imax;
fint4 *jmax;
fint4 *kmax;
#endif
{
        int i_status, i_varid;
        int temp, str_len, i;
        size_t mindex[1];
        char *t_ptr;
 
        mindex[0] = 0;
 
/* get the data variable id */
        i_status = nc_inq_varid (cdfid, "imax",&i_varid);
        if (i_status != NC_NOERR) {
          printf("imax id %s\n",nc_strerror(i_status));
          return -1;
	}
/* read the var from the netcdf file */
        i_status = ncgetvar1(cdfid,i_varid,mindex,imax);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          return -1;
	}
 
/* get the data variable id */
        i_status = nc_inq_varid (cdfid, "jmax",&i_varid);
        if (i_status != NC_NOERR) {
          printf("jmax id %s\n",nc_strerror(i_status));
          return -1;
	}
 
/* read the var from the netcdf file */
        i_status = ncgetvar1(cdfid,i_varid,mindex,jmax);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          return -1;
	}
 
/* get the data variable id */
        i_status = nc_inq_varid (cdfid, "kmax",&i_varid);
        if (i_status != NC_NOERR) {
          printf("kmax id %s\n",nc_strerror(i_status));
          return -1;
	}
 
/* read the var from the netcdf file */
        i_status = ncgetvar1(cdfid,i_varid,mindex,kmax);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          return -1;
	}
 
/* normal return */
 
        return 0;
}

/************************************************************/
#ifdef __STDC__
int retrieve_grid_v3(int i_cdfid,int i_level,int record_indx,
                     char *var, float *dptr, char *cptr,
                     int name_len, char *lvlptr, char *uptr)
#else
int retrieve_grid_v3(i_cdfid, i_level, record_indx, var,
                           dptr, cptr, name_len, lvlptr, uptr)
int i_cdfid;
int i_level;
int record_indx;
char *var;
float *dptr;
char *cptr;
int name_len;
char *lvlptr;
char *uptr;
#endif
{
        int i_status, i_invflag, i_varid, lvl_indx, x_dim, y_dim;
        size_t start[4],count[4],start_c[3],count_c[3];
        char var_name[13];
 
/* get the data variable id */
        i_status = nc_inq_varid (i_cdfid,var,&i_varid);
        if (i_status != NC_NOERR) {
          printf("Variable %s not found in input file.\n",var);
          if (DEBUG == 1)
            printf("%s\n",nc_strerror(i_status));
            printf("cdf_retrieve_laps: no grid available.\n");
          return -1;
        }
 
/* get the level index of the data array */
        lvl_indx = get_index_v3 (i_cdfid, i_level, "z");
        if (lvl_indx == (-1)){
          printf("Level %d for variable %s not found in input file.\n",i_level, var);
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: no level %d\n", i_level);
          return -1;
        }
 
/* check the inventory variable to see if this grid is available */
        i_invflag = check_laps_inv (i_cdfid, lvl_indx,
                                    record_indx, var);
        if (i_invflag == (-1)){
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: grid %s not available.\n",var);
          return -1;
        }
 
/* get the x and y dimension sizes */
        x_dim = dim_size_v3 (i_cdfid, "x");
        y_dim = dim_size_v3 (i_cdfid, "y");
 
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
        i_status = nc_get_vara_float (i_cdfid,i_varid,start,count,dptr);
        if (i_status != NC_NOERR) {
          printf("Error reading grid %s from input file.\n",var);
          if (DEBUG == 1)
            printf("%s\n",nc_strerror(i_status));
            printf("cdf_retrieve_laps: error retrieving data %s grid.\n",
                         var);
          return -1;
        }
 
/* get attributes lvl_coord and LAPS_units */
        i_status = nc_get_att_text (i_cdfid,i_varid,"lvl_coord",lvlptr);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: error retrieving lvl_coord.\n");
          return -1;
        }
 
        i_status = nc_get_att_text (i_cdfid,i_varid,"LAPS_units",uptr);
        if (i_status != NC_NOERR) {
          printf("%s\n",nc_strerror(i_status));
          if (DEBUG == 1)
            printf("cdf_retrieve_laps: error retrieving LAPS_units.\n");
          return -1;
        }
 
/* setup to read the comment from the netcdf file */
 
        sprintf(var_name, "%s%s", var, "_comment");
        i_status = nc_inq_varid (i_cdfid,var_name,&i_varid);
        if (i_status != NC_NOERR) {
          printf("Error finding variable %s.\n",var_name);
          if (DEBUG == 1)
            printf("%s\n",nc_strerror(i_status));
            printf("cdf_retrieve_laps: no comment field available.\n");
          return -1;
        }
 
/* construct the arrays needed to read the comment */
        start_c[0] = record_indx;
        start_c[1] = lvl_indx;
        start_c[2] = 0;
 
        count_c[0] = 1;
        count_c[1] = 1;
        count_c[2] = name_len;
 
/* read the comment from the netcdf file */
        i_status = nc_get_vara_text (i_cdfid,i_varid,start_c,count_c,cptr);
        if (i_status != NC_NOERR) {
          printf("Error reading %s.\n",var_name);
          if (DEBUG == 1)
            printf("%s\n",nc_strerror(i_status));
            printf("cdf_retrieve_laps: error retrieving comment.\n");
          return -1;
        }
 
/* normal return */
 
        return 0;
}

/************************************************************/
#ifdef __STDC__
void cstr_to_fstr(char *f_comment, char *comment, fint4 *comm_len,
                  int name_len,
                  char *f_lvl_coord, char *lvl_coord, 
                  fint4 *lvl_coord_len, char *f_units, char *units,
                  fint4 *units_len, fint4 *kdim)
#else
void cstr_to_fstr(f_comment, comment, comm_len, name_len, f_lvl_coord, 
                  lvl_coord, lvl_coord_len, f_units, units,
                  units_len, kdim)
char *f_comment; 
char *comment; 
fint4 *comm_len;
int  name_len;
char *f_lvl_coord; 
char *lvl_coord;
fint4 *lvl_coord_len; 
char *f_units; 
char *units;
fint4 *units_len; 
fint4 *kdim;
#endif
{
	char *uptr, *cptr, *lptr, *pf, *pc, pc_char[1];
        char *hld_unit, *hld_comm, *hld_lvl;
        int i, j, slen;

        *f_units = '\0';
        *f_comment = '\0';
        *f_lvl_coord = '\0';
         
        hld_unit = malloc(*units_len + 1);
        pc = hld_unit + (*units_len);
        *pc = '\0';
        hld_comm = malloc(name_len);
        pc = hld_comm + (name_len - 1);
        *pc = '\0';
        hld_lvl = malloc(*lvl_coord_len + 1);
        pc = hld_lvl + (*lvl_coord_len);
        *pc = '\0';
        uptr = units;
        cptr = comment;
        lptr = lvl_coord;

        for (i = 0; i < *kdim ; i++) {
          slen = strlen(uptr);
          if (slen < (*units_len)) {
            strcpy(hld_unit,uptr);
            for (j = slen + 1; j <= *units_len; j++)
              strcat(hld_unit," ");
          }
          else {
            strncpy(hld_unit,uptr,*units_len);
          }
          if (i == (*kdim - 1)) {
            pf = f_units + ((*kdim - 1) * (*units_len));
            pc = hld_unit;
            for (j = 0; j < *units_len; j++) {
              strncpy(pc_char,pc,1);
              *pf = pc_char[0];
              pf++;
              pc++;
            }
          }
          else {
            strcat(f_units, hld_unit);
          }
          uptr += (*units_len) + 1;

          slen = strlen(cptr);
          if (slen < (*comm_len)) {
            strcpy(hld_comm,cptr);
            for (j = slen + 1; j <= *comm_len; j++)
              strcat(hld_comm," ");
          }
          else {
            strncpy(hld_comm,cptr,*comm_len);
          }
          if (i == (*kdim - 1)) {
            pf = f_comment + ((*kdim - 1) * (*comm_len));
            pc = hld_comm;
            for (j = 0; j < *comm_len; j++) {
              strncpy(pc_char,pc,1);
              *pf = pc_char[0];
              pf++;
              pc++;
            }
          }
          else {
            strcat(f_comment, hld_comm);
          }
          cptr += (*comm_len) + 1;

          slen = strlen(lptr);
          if (slen < (*lvl_coord_len)) {
            strcpy(hld_lvl,lptr);
            for (j = slen + 1; j <= *lvl_coord_len; j++)
              strcat(hld_lvl, " ");
          }
          else {
            strncpy(hld_lvl,lptr,*lvl_coord_len);
          }
          if (i == (*kdim - 1)) {
            pf = f_lvl_coord + ((*kdim - 1) * (*lvl_coord_len));
            pc = hld_lvl;
            for (j = 0; j < *lvl_coord_len; j++) {
              strncpy(pc_char,pc,1);
              *pf = pc_char[0];
              pf++;
              pc++;
            }
          }
          else {
            strcat(f_lvl_coord, hld_lvl);
          }
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
                  fint4 *var_len, fint4 *comm_len, fint4 *fn_length, 
                  fint4 *ext_len, fint4 *lvl_coord_len, fint4 *units_len, 
                  fint4 *i_reftime, fint4 *i_valtime, fint4 *iimax, 
                  fint4 *jjmax, fint4 *kkmax, fint4 *kdim, fint4 lvl[], 
                  float *data, fint4 *called_from, fint4 *status)
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
fint4 *var_len;
fint4 *comm_len;
fint4 *fn_length;
fint4 *ext_len;
fint4 *lvl_coord_len;
fint4 *units_len;
fint4 *i_reftime;
fint4 *i_valtime;
fint4 *iimax;
fint4 *jjmax;
fint4 *kkmax;
fint4 *kdim;
fint4 lvl[];
float *data;
fint4 *called_from;
fint4 *status;
#endif
{
        int cdfid, istatus, name_len, name_len_null;
        int i,j, t_level, i_record, unconv_var;
        int val_id,dim_id, t_record, int_1, int_2, found,t_var_id;
        long num_record;
        fint4 imax, jmax, kmax; 
        size_t mindex[1], dim_len;
        float *dptr;
        double reftime, valtime, d_valtime, timeoff, diff; 
        char *filename, *cpt, char1[3];
	char *ext, *var, *comment, *lvl_coord, *units; 
	char *comm_var, *inv_var;
        char *fvptr, *vptr, *lptr, *uptr, *cptr, *t_var;

        reftime = (double) *i_reftime;

/* convert fortran file_name into C fname  */
        filename = malloc(*fn_length + 1);
        nstrncpy(filename,f_filename,*fn_length);

/* open netCDF file for reading, if file exists */
        if( access(filename, F_OK) != 0 ) {
          printf("The LAPS file %s does not exist\n",filename);
          *status = -5; /* file not there */
          return;
        }

        istatus = nc_open((const char*)filename,NC_NOWRITE, &cdfid);
        free(filename);
        if (istatus != NC_NOERR) {
          printf("nc_open %s\n",nc_strerror(istatus));
          *status = -1; /* error opening file */
          return;
        }

/* set i_record to value returned by nc_inq_dimlen */

        istatus = nc_inq_dimid(cdfid,"record",&dim_id);
        if (istatus != NC_NOERR) {
          printf("%s\n",nc_strerror(istatus));
          istatus = nc_close(cdfid);
          *status = -2; /* error in reading header info */
          printf("Error reading record in header info\n");
          return;
        }
        
        istatus = nc_inq_dimlen(cdfid, dim_id, &dim_len);
        num_record = (long) dim_len;
        if ((istatus != NC_NOERR) || (num_record == 0)) { /* no data in file */
          printf("dimlen dim_id %s\n",nc_strerror(istatus));
          istatus = nc_close(cdfid);
          *status = -2; /* error in reading header info */
          return;
        }
        else {
          if (num_record >= 1)
            i_record = (int) (num_record - 1);
        }

/* determine value of namelen dimension in file */

        istatus = nc_inq_dimid(cdfid,"namelen",&dim_id);
        if (istatus != NC_NOERR) {
          printf("%s\n",nc_strerror(istatus));
          istatus = nc_close(cdfid);
          *status = -2; /* error in reading header info */
          printf("Error reading namelen in header info\n");
          return;
        }
        
        istatus = nc_inq_dimlen(cdfid, dim_id, &dim_len);
        name_len = (int) dim_len;
        name_len_null = name_len;
        if (istatus != NC_NOERR) { /* no "namelen" dimension */
          printf("%s\n",nc_strerror(istatus));
          istatus = nc_close(cdfid);
          *status = -2; /* error in reading header info */
          printf("Error reading namelen dimension in header info\n");
          return;
        }

/* allocate space for c strings */
	var = malloc((*var_len + 1) * sizeof(char) * (*kdim));
	ext = malloc((*ext_len + 1) * sizeof(char));
        lvl_coord = malloc((*lvl_coord_len + 1) * sizeof(char) * (*kdim));
        units = malloc((*units_len + 1) * sizeof(char) * (*kdim));
        
/* set length of comment to shortest of namelen dimension or comm_len */ 
        if (*comm_len  < name_len ) {
          name_len = *comm_len + 1; /* number of characters to pull */
          name_len_null = *comm_len;  /* location to write null in comment */
          comment = malloc((*comm_len + 1) * sizeof(char) * (*kdim));
        }
        else {  /* *comm_len >= name_len */
          name_len_null = name_len - 1;
          comment = malloc(name_len * sizeof(char) * (*kdim));
        }

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
          for ( j = 0; j < (*lvl_coord_len + 1); j++) {
            *lptr = '\0';
            lptr++;
          }
          for ( j = 0; j < (*units_len + 1); j++) {
            *uptr = '\0';
            uptr++;
          }
          for ( j = 0; j < name_len; j++) {
            *cptr = '\0';
            cptr++;
          }
        }
        nstrncpy(ext,f_ext,*ext_len);
        upcase_c(ext, ext);

        istatus = retrieve_hdr_v3(cdfid,&imax,&jmax,&kmax);
        if (istatus == -1) {
          istatus = nc_close(cdfid);
          free_read_var(var, comment, ext, lvl_coord, units,
                        comm_var, inv_var);
          *status = -2; /* error in reading header info */
          printf("Error reading grid dimensions in header info\n");
          return;
        }

        if (imax > *iimax || jmax > *jjmax || *kkmax > *kdim) {
          istatus = nc_close(cdfid);
          free_read_var(var, comment, ext, lvl_coord, units,
                        comm_var, inv_var);
          *status = -3; /* error in grid_size */
          return;
        }
 
        if ((strcmp(ext,"LMR") == 0) || (strcmp(ext,"LF1") == 0)) {

          istatus = nc_inq_varid(cdfid,"valtime",&val_id);
          if (istatus != NC_NOERR) {
            istatus = nc_close(cdfid);
            free_read_var(var, comment, ext, lvl_coord, units,
                          comm_var, inv_var);
            *status = -2; /* error in reading header info */
            printf("Error reading valtime in header info\n");
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
            char1[0] = *(vptr+1);
            char1[1] = *(vptr+2);
            char1[2] = '\0';
 
            int_1 = atoi(char1);
            fvptr = vptr;
            fvptr++;
            *fvptr = '\0';
            t_var = vptr;

            timeoff = (double) int_1*600;

            diff = (double) (*i_valtime - *i_reftime);
            valtime = reftime + timeoff + diff;
            i = 0;
            found = 0;
            while ((i < num_record) && (found == 0)) {
              mindex[0] = i;
              istatus = nc_get_var1_double(cdfid,val_id,mindex,&d_valtime); 
              if (istatus != NC_NOERR) {
                istatus = nc_close(cdfid);
                free_read_var(var, comment, ext, lvl_coord, units,
                              comm_var, inv_var);
                *status = -2; /* error in reading header info */
                printf("Error in nc_get_var1_double\n");
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
                                     name_len, lptr, uptr); 

          if (istatus == -1)  {
            unconv_var += 1;
          }
          dptr += (*iimax)*(*jjmax);
          vptr += (*var_len) + 1;
          cptr += name_len;
          uptr += (*units_len) + 1;
          lptr += (*lvl_coord_len) + 1;
 
        }
 
        istatus = nc_close(cdfid);

        if (*called_from == 0) {
          cstr_to_fstr(f_comment, comment, comm_len, name_len,
                       f_lvl_coord, lvl_coord, lvl_coord_len, f_units, 
                       units, units_len, kdim);
        }
 
        free_read_var(var, comment, ext, lvl_coord, units,
                      comm_var, inv_var);
        *status = (fint4) unconv_var;
        return;

}
/************************************************************/
#ifdef __STDC__
int write_val_ref_asctime(int cdfid, int i_record, double *valtime,
                          double *reftime, char *asctime,
                          fint4 *asc_len, int name_len)
#else
int write_val_ref_asctime(cdfid, i_record, valtime, reftime, 
                          asctime, asc_len, name_len)
int cdfid; 
int i_record;
double *valtime; 
double *reftime;
char *asctime;
fint4 *asc_len;
int name_len;
#endif
{
        int varid, istatus;
        size_t mindex[1], start_2[2], count_2[2];

        mindex[0] = i_record;
        istatus = nc_inq_varid(cdfid,"valtime",&varid);
        if (istatus != NC_NOERR) {
          return -1;
        }
        else {
          if (DEBUG==1) printf("found varid for valtime\n");
        }

        istatus = nc_put_var1_double(cdfid,varid,mindex, valtime);
        if (istatus != NC_NOERR) {
          return -1;
        }
        else {
        if (DEBUG==1) printf("correctly wrote valtime\n");
        }
            
        istatus = nc_inq_varid(cdfid,"reftime",&varid);
        if (istatus != NC_NOERR) {
          return -1;
        }
        else {
        if (DEBUG==1) printf("found varid for reftime\n");
        }

        istatus = nc_put_var1_double(cdfid, varid, mindex, reftime);
        if (istatus != NC_NOERR) {
          return -1;
        }
        else {
        if (DEBUG==1) printf("correctly wrote reftime\n");
        }
            
/* write out asctime */
        istatus = nc_inq_varid(cdfid,"asctime",&varid);
        if (istatus != NC_NOERR) {
          return (-1);
        }
        else {
        if (DEBUG==1) printf("found varid for asctime\n");
        }

	start_2[0] = i_record;
        start_2[1] = 0;
        count_2[0] = 1;
        if (name_len < *asc_len) {
          count_2[1] = name_len;
          asctime[name_len-1] = '\0';
        }
        else {
          count_2[1] = *asc_len;
        }
        if (DEBUG == 1) { 
          printf("start = [%d][%d] | count = [%d][%d] \n", start_2[0],start_2[1],count_2[0],count_2[1]);
          printf("asctime = [%s]\n",asctime);
        }

        istatus = nc_put_vara_text(cdfid,varid,start_2,count_2,asctime);
        if (istatus != NC_NOERR) {
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
int check_grid_dims_v3(int cdfid, char *ext, fint4 *imax, 
                          fint4 *jmax, fint4 *n_levels, fint4 lc3_levels,
                          fint4 lm1_levels)
#else
int check_grid_dims_v3(cdfid, ext, imax, jmax, n_levels, lc3_levels,
                          lm1_levels)
int cdfid;
char *ext;
fint4 *imax; 
fint4 *jmax; 
fint4 *n_levels;
fint4 lc3_levels;
fint4 lm1_levels;
#endif
{
        int dimid, istatus;
        size_t dim_len;
        long dim_val;

/*      printf("check_grid_dims_v3 location 1 \n");   */

/* determine value in dimension x */
        istatus = nc_inq_dimid(cdfid,"x",&dimid);
        if (istatus != NC_NOERR) {
          printf("No x dimension found in output file.\n");
          return -1;
        }
        istatus = nc_inq_dimlen(cdfid,dimid,&dim_len);
        if (istatus != NC_NOERR) {
          printf("Unable to access x dimension in output file.\n");
          return -1;
        }
        dim_val = (long)dim_len;

        if (*imax != dim_val) {
          printf("x dimension in output file does not match IMAX passed in rwl_v3.c: %d %ld\n",*imax,dim_val);
          return -1;
        }

/*      printf("check_grid_dims_v3 location 2 \n");   */

/* determine value in dimension y */
        istatus = nc_inq_dimid(cdfid,"y",&dimid);
        if (istatus != NC_NOERR) {
          printf("No y dimension found in output file.\n");
          return -1;
        }
        istatus = nc_inq_dimlen(cdfid,dimid,&dim_len);
        if (istatus != NC_NOERR) {
          printf("Unable to access y dimension in output file.\n");
          return -1;
        }
        dim_val = (long)dim_len;

        if (*jmax != dim_val) {
          printf("y dimension in output file does not match JMAX passed in rwl_v3.c: %d %ld\n",*jmax,dim_val);
          return -1;
        }

/*      printf("check_grid_dims_v3 location 3 \n");   */

/* determine value in dimension z */
        istatus = nc_inq_dimid(cdfid,"z",&dimid);
        if (istatus != NC_NOERR) {
          printf("No z dimension found in output file.\n");
          return -1;
        }
        istatus = nc_inq_dimlen(cdfid,dimid,&dim_len);
        if (istatus != NC_NOERR) {
          printf("Unable to access z dimension in output file.\n");
          return -1;
        }
        dim_val = (long)dim_len;

        printf("check_grid_dims_v3 location 4 \n");       

	if (dim_val == 1) {
          printf("check_grid_dims_v3 location 4a %d \n",*n_levels);
          *n_levels = 1;   
        }
        else if ((strncmp(ext,"lc3",3) == 0) && (dim_val == (long) LC3_LEVELS)) {
          printf("check_grid_dims_v3 location 4b \n");
          *n_levels = (fint4) LC3_LEVELS;
        }
        else if ((strncmp(ext,"lm1",3) == 0) && (dim_val == (long) LM1_LEVELS)) {
          printf("check_grid_dims_v3 location 4c \n");
          *n_levels = (fint4) LM1_LEVELS;
        }

        if(dim_val != *n_levels) {
          printf("Z dimension in output file does not match N_LEVELS passed in rwl_v3.c: %d %ld\n",*n_levels,dim_val);
          printf("Check CDL, levels namelist, and localization.\n");
          return -1;
        }

        printf("check_grid_dims_v3 location 5 \n");
    
/* normal return */

        return(0);  /* dimensions check out OK */
}
/************************************************************/
#ifdef __STDC__
int get_static_info_v3(char *static_grid, float *Dx, float *Dy, 
                      float *La1, float *Lo1, float *LoV, 
                      float *Latin1, float *Latin2, char *map_proj, 
                      char *origin, int name_len, char* ldf)
#else
int get_static_info_v3(static_grid, Dx, Dy, La1, Lo1, LoV, Latin1,
                       Latin2, map_proj, origin, name_len, ldf)
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
int name_len;
char *ldf;
#endif
{
        int cdfid_stat, dimid, varid, istatus, ret_status;
        long namelen;
        size_t mindex[1], start[1], count[1], dim_len;
        size_t start_g[2], count_g[2];

         
	istatus = nc_open(static_grid, NC_NOWRITE, &cdfid_stat);
        if (istatus != NC_NOERR) {
          printf("Unable to open static file %s.\n",static_grid);
          return -1;
        }

        ret_status = 0;

/* determine value in dimension namelen */
        istatus = nc_inq_dimid(cdfid_stat,"namelen",&dimid);
        if (istatus != NC_NOERR) {
          namelen = 132;
          printf("Variables map_proj and origin may be incorrect.\n");
          printf("Unable to read dimension namelen from static.%s\n", *ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_inq_dimlen(cdfid_stat,dimid,&dim_len);
          if (istatus != NC_NOERR) {
            namelen = 132;
            printf("Variables map_proj and origin may be incorrect.\n");
            printf("Unable to read dimension namelen from static.%s\n", *ldf);
            ret_status = -1;
          }
        }
/* namelen is length in static.nest7grid; name_len is length in LAPS output files */
        namelen = (long)dim_len;

        mindex[0] = 0;

/* read Dx */
        istatus = nc_inq_varid(cdfid_stat,"Dx",&varid);
        if (istatus != NC_NOERR) {
          printf("Error reading Dx from  static.%s\n",*ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_get_var1_float (cdfid_stat,varid,mindex,Dx);
          if (istatus != NC_NOERR) {
            printf("Error reading Dx from  static.%s\n", *ldf);
            ret_status = -1;
          }
        }
        
/* read Dy */
        istatus = nc_inq_varid(cdfid_stat,"Dy",&varid);
        if (istatus != NC_NOERR) {
          printf("Error reading Dy from  static.%s\n", *ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_get_var1_float (cdfid_stat,varid,mindex,Dy);
          if (istatus != NC_NOERR) {
            printf("Error reading Dy from  static.%s\n", *ldf);
            ret_status = -1;
          }
        }
        
/* read La1 */
        istatus = nc_inq_varid(cdfid_stat,"La1",&varid);
        if (istatus != NC_NOERR) {
          printf("Error reading La1 from  static.%s\n", *ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_get_var1_float (cdfid_stat,varid,mindex,La1);
          if (istatus != NC_NOERR) {
            printf("Error reading La1 from  static.%s\n", *ldf);
            ret_status = -1;
          }
        }
        
/* read Lo1 */
        istatus = nc_inq_varid(cdfid_stat,"Lo1",&varid);
        if (istatus != NC_NOERR) {
          printf("Error reading Lo1 from  static.%s\n", *ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_get_var1_float (cdfid_stat,varid,mindex,Lo1);
          if (istatus != NC_NOERR) {
            printf("Error reading Lo1 from  static.%s\n", *ldf);
            ret_status = -1;
          }
        }
        
/* read LoV */
        istatus = nc_inq_varid(cdfid_stat,"LoV",&varid);
        if (istatus != NC_NOERR) {
          printf("Error reading LoV from  static.%s\n", *ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_get_var1_float (cdfid_stat,varid,mindex,LoV);
          if (istatus != NC_NOERR) {
            printf("Error reading LoV from  static.%s\n", *ldf);
            ret_status = -1;
          }
        }
        
/* read Latin1 */
        istatus = nc_inq_varid(cdfid_stat,"Latin1",&varid);
        if (istatus != NC_NOERR) {
          printf("Error reading Latin1 from  static.%s\n", *ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_get_var1_float (cdfid_stat,varid,mindex,Latin1);
          if (istatus != NC_NOERR) {
            printf("Error reading Latin1 from  static.%s\n", *ldf);
            ret_status = -1;
          }
        }
        
/* read Latin2 */
        istatus = nc_inq_varid(cdfid_stat,"Latin2",&varid);
        if (istatus != NC_NOERR) {
          printf("Error reading Latin2 from  static.%s\n", *ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_get_var1_float (cdfid_stat,varid,mindex,Latin2);
          if (istatus != NC_NOERR) {
            printf("Error reading Latin2 from  static.%s\n", *ldf);
            ret_status = -1;
          }
        }
        
        start_g[0] = 0;
        start_g[1] = 0;
        count_g[0] = 1;
        if (name_len > namelen) {
          count_g[1] = namelen;
        }
        else {
          count_g[1] = name_len;
        }

/* read grid_type */
        istatus = nc_inq_varid(cdfid_stat,"grid_type",&varid);
        if (istatus != NC_NOERR) {
          printf("Error reading grid_type from static.%s\n", *ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_get_vara_text (cdfid_stat, varid,start_g,count_g,map_proj);
          if (istatus != NC_NOERR) {
            printf("Error reading grid_type from  static.%s\n", *ldf);
            ret_status = -1;
          }
          if (name_len < namelen) map_proj[name_len-1] = '\0';
        }

        start[0] = 0;
        if (name_len > namelen) {
          count[0] = namelen;
        }
        else {
          count[0] = name_len;
        }
        
/* read origin_name */
        istatus = nc_inq_varid(cdfid_stat,"origin_name",&varid);
        if (istatus != NC_NOERR) {
          printf("Error reading origin_name from  static.%s\n", *ldf);
          ret_status = -1;
        }
        else {
          istatus = nc_get_vara_text (cdfid_stat,varid,start,count,origin);
          if (istatus != NC_NOERR) {
            printf("Error reading origin_name from  static.%s\n", *ldf);
            ret_status = -1;
          }
          if (name_len < namelen) origin[name_len-1] = '\0';
        }
        
        istatus = nc_close(cdfid_stat);
        if (ret_status == -1) {
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
        return;
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
        return;
}

/************************************************************/
#ifdef __STDC__
void fill_c_var(fint4 *kdim, char *f_var, char *var, fint4 *var_len,
               char *f_comment, char *comment, fint4 *comm_len,
               char *f_asctime, char *asctime, fint4 *asc_len)
#else
void fill_c_var(kdim,f_var,var,var_len,f_comment,comment,comm_len,
                f_asctime, asctime,asc_len)
fint4 *kdim;
char *f_var;
char *var;
fint4 *var_len;
char *f_comment;
char *comment;
fint4 *comm_len;
char *f_asctime;
char *asctime;
fint4 *asc_len;
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
int  write_hdr_v3(int cdfid, char *ext, int i_record, fint4 *imax, 
                  fint4 *jmax, fint4 *kmax, fint4 *kdim, float *pr, 
                  fint4 *n_levels, float *Dx, 
                  float *Dy, float *La1, float *Lo1, float *LoV, 
                  float *Latin1, float *Latin2, char *map_proj, 
                  char *origin, int name_len, float *cdl_levels)
#else
int  write_hdr_v3(cdfid, ext, i_record, imax, jmax, kmax, kdim, pr, 
                  n_levels, Dx, Dy, La1, Lo1, LoV, Latin1, 
                  Latin2, map_proj, origin, name_len, cdl_levels)
int cdfid; 
char *ext;
int i_record; 
fint4 *imax; 
fint4 *jmax; 
fint4 *kmax; 
fint4 *kdim; 
float *pr; 
fint4 *n_levels; 
float *Dx; 
float *Dy; 
float *La1; 
float *Lo1; 
float *LoV; 
float *Latin1; 
float *Latin2;
char *map_proj; 
char *origin;
int  name_len;
float *cdl_levels;
#endif
{
        short Nx, Ny;
	int i, j, istatus, varid;
        size_t start[1], count[1], mindex[1], start_2[2], count_2[2];
        float level_val, *levels, *lp;
        
        istatus = 0;
        mindex[0] = 0;

/* write out imax */
        istatus = nc_inq_varid(cdfid,"imax",&varid);
        if (istatus != NC_NOERR){
          printf("%s\n",nc_strerror(istatus));
          return (-1);
	}
        if (DEBUG==1) printf("found varid for imax\n");
       
        istatus = ncputvar1(cdfid,varid,mindex,(const fint4 *)imax);
        if (istatus != NC_NOERR){
          printf("%s\n",nc_strerror(istatus));
          return (-1);
	}
        if (DEBUG==1) printf("correctly wrote imax\n");

/* write out jmax */
        istatus = nc_inq_varid(cdfid,"jmax",&varid);
        if (istatus != NC_NOERR){
          printf("%s\n",nc_strerror(istatus));
          return (-1);
	}
        if (DEBUG==1) printf("found varid for jmax\n");

        istatus = ncputvar1(cdfid,varid,mindex,(const fint4 *)jmax);
        if (istatus != NC_NOERR){
          printf("%s\n",nc_strerror(istatus));
          return (-1);
	}
        if (DEBUG==1) printf("correctly wrote jmax\n");

/* write out kmax */
        istatus = nc_inq_varid(cdfid,"kmax",&varid);
        if (istatus != NC_NOERR){
          printf("%s\n",nc_strerror(istatus));
          return (-1);
	}
        if (DEBUG==1) printf("found varid for kmax\n");

        istatus = ncputvar1(cdfid,varid,mindex,(const fint4 *)kmax);
        if (istatus != NC_NOERR){
          printf("%s\n",nc_strerror(istatus));
          return (-1);
	}
        if (DEBUG==1) printf("correctly wrote kmax\n");

/* write out kdim */
        istatus = nc_inq_varid(cdfid,"kdim",&varid);
        if (istatus != NC_NOERR){
          printf("%s\n",nc_strerror(istatus));
          return (-1);
	}
        if (DEBUG==1) printf("found varid for kdim\n");

        istatus = ncputvar1(cdfid,varid,mindex,(const fint4 *)kdim);
        if (istatus != NC_NOERR){
          printf("%s\n",nc_strerror(istatus));
          return (-1);
	}
        if (DEBUG==1) printf("correctly wrote kdim\n");

/* if ext is fua or fsf or pbl, cdl should have Nx, Ny, Dx, Dy, La1, Lo1, Latin1,
     Latin2, grid_type and origin_name set as data in the cdl.  The variable
     "level" will be filled using cdl_levels rather than from base and interval */

        if ((strncmp(ext,"fua",3) == 0) || (strncmp(ext,"fsf",3) == 0) ||
            (strncmp(ext,"pbl",3) == 0)) {
          start[0] = 0;
          count[0] = *n_levels;

          istatus = nc_inq_varid(cdfid,"level",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for level\n");

          istatus = nc_put_vara_float(cdfid,varid,start,count,cdl_levels);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote level\n");
        }
        else {
/* write out levels */
/* levels written are based on value in n_levels, which has been adjusted to
     one of four options: standard LAPS 3D (nk_laps in nest7grid.parms),
                            This then uses the values in pressure.nl for the levels
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
/* the ordering in *pr is from the bottom up, we need top down ordering in the
   file, so invert what is in *pr and put it into levels */

            j = *n_levels -1;
            for( i = 0; i < *n_levels; i++) {
              levels[i] = pr[j];
              j--;
            }
      
          }

          start[0] = 0;
          count[0] = *n_levels;

          istatus = nc_inq_varid(cdfid,"level",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for level\n");

          istatus = nc_put_vara_float(cdfid,varid,start,count,levels);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote level\n");

          free (levels);

          mindex[0] = 0;

/* write out Nx */
          istatus = nc_inq_varid(cdfid,"Nx",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for Nx\n");

          Nx = (short) *imax;
          istatus = nc_put_var1_short(cdfid,varid,mindex,&Nx);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote Nx\n");

/* write out Ny */
          istatus = nc_inq_varid(cdfid,"Ny",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for Ny\n");

          Ny = (short) *jmax;
          istatus = nc_put_var1_short(cdfid,varid,mindex,&Ny);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote Ny\n");

/* write out Dx */
          istatus = nc_inq_varid(cdfid,"Dx",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for Dx\n");

          istatus = nc_put_var1_float(cdfid,varid,mindex,Dx);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote Dx\n");

/* write out Dy */
          istatus = nc_inq_varid(cdfid,"Dy",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for Dy\n");

          istatus = nc_put_var1_float(cdfid,varid,mindex,Dy);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote Dy\n");

/* write out La1 */
          istatus = nc_inq_varid(cdfid,"La1",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for La1\n");

          istatus = nc_put_var1_float(cdfid,varid,mindex,La1);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote La1\n");

/* write out Lo1 */
          istatus = nc_inq_varid(cdfid,"Lo1",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for Lo1\n");

          istatus = nc_put_var1_float(cdfid,varid,mindex,Lo1);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote Lo1\n");

/* write out LoV */
          istatus = nc_inq_varid(cdfid,"LoV",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for LoV\n");

          istatus = nc_put_var1_float(cdfid,varid,mindex,LoV);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote LoV\n");

/* write out Latin1 */
          istatus = nc_inq_varid(cdfid,"Latin1",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for Latin1\n");

          istatus = nc_put_var1_float(cdfid,varid,mindex,Latin1);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote Latin1\n");

/* write out Latin2 */
          istatus = nc_inq_varid(cdfid,"Latin2",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for Latin2\n");

          istatus = nc_put_var1_float(cdfid,varid,mindex,Latin2);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote Latin2\n");

/* write out map_proj */
          istatus = nc_inq_varid(cdfid,"grid_type",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for grid_type\n");

	  start_2[0] = i_record;
          start_2[1] = 0;
          count_2[0] = 1;
          if (strlen(map_proj) > name_len) {
            count_2[1] = name_len;
            map_proj[name_len - 1] = '\0';
          }
          else {
            count_2[1] = strlen(map_proj);
          }

          istatus = nc_put_vara_text(cdfid,varid,start_2,count_2,map_proj);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote map_proj\n");

/* write out origin */
          istatus = nc_inq_varid(cdfid,"origin_name",&varid);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("found varid for origin_name\n");

          start[0] = 0;
          if (strlen(origin) > name_len) {
            count[0] = name_len;
            origin[name_len - 1] = '\0';
          }
          else {
            count[0] = strlen(origin);
          }

          istatus = nc_put_vara_text(cdfid,varid,start,count,origin);
          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }
          if (DEBUG==1) printf("correctly wrote origin\n");


          if (istatus != NC_NOERR){
            printf("%s\n",nc_strerror(istatus));
            return (-1);
	  }

        }
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
        size_t start[2];
        short i_flag;

	start[0] = i_record;
        start[1] = z_dim;

        i_flag = 1;
        istatus = nc_put_var1_short(cdfid,inv_id,start,&i_flag);
        if (istatus != NC_NOERR){
          printf("%s\n",nc_strerror(istatus));
          return istatus;
	}
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
        long dim_val;
        size_t start[1],count[1], dim_len;
        float *file_lvls, *flvls;

/* determine value in dimension z */
        istatus = nc_inq_dimid(cdfid,dim,&dim_id);
        if (istatus != NC_NOERR) {
          printf("%s\n",nc_strerror(istatus));
          return -1;
	}

        istatus = nc_inq_dimlen(cdfid,dim_id,&dim_len);
        if (istatus != NC_NOERR) {
          printf("%s\n",nc_strerror(istatus));
          return -1;
	}
        dim_val = (long)dim_len;

/* allocate space for dim_val floats, and pull levels */
        file_lvls = (float *) malloc(dim_val * sizeof(float));
        istatus = nc_inq_varid(cdfid,"level",&var_id);
        if (istatus != NC_NOERR) {
          printf("%s\n",nc_strerror(istatus));
          return -1;
	}

        start[0] = 0;
        count[0] = dim_val;
        istatus = nc_get_vara_float(cdfid,var_id,start,count,file_lvls);
        if (istatus != NC_NOERR) {
          printf("%s\n",nc_strerror(istatus));
          return -1;
	}

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
        else {
          printf("Level %8.2f not found.\n", match_level);
          return -1;
        }
        
}
/************************************************************/
#ifdef __STDC__
int update_laps_v3(int cdfid,fint4 i_level,int i_record,
                   fint4 *imax, fint4 *jmax, int var_id, 
                   int inv_id,float *dptr,
                   int comm_id, char *cptr,int name_len)
#else
int update_laps_v3(cdfid,i_level,i_record,imax,jmax, 
                   var_id,inv_id,dptr,comm_id,cptr, name_len)
int cdfid;
fint4 i_level;
int i_record;
fint4 *imax;
fint4 *jmax;
int var_id;
int inv_id;
float *dptr;
int comm_id;
char *cptr;
int name_len;
#endif
{
	int x_dim, y_dim, z_dim, istatus, dim_id;
        size_t start[4], count[4], c_start[3], c_count[3];

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

        istatus = nc_put_vara_float(cdfid,var_id,start,count,dptr);
        if (istatus != NC_NOERR){ /* error writing out grid */
          printf("%s\n",nc_strerror(istatus));
          return (-1);
        }

	istatus = update_inv_v3(cdfid,i_record,z_dim,inv_id);
	if (istatus == (-1)) /* error writing out grid */{
          printf("%s\n",nc_strerror(istatus));
          return (-1);
        }

	c_start[0] = i_record;
	c_start[1] = z_dim;
	c_start[2] = 0;
 
	c_count[0] = 1;
	c_count[1] = 1;
        if (strlen(cptr) > name_len) {
          c_count[2] = name_len;
          cptr[name_len - 1] = '\0';
        }
        else {
          c_count[2] = strlen(cptr);
        }

	istatus = nc_put_vara_text(cdfid,comm_id,c_start,
				   c_count,cptr);

	if (istatus != NC_NOERR){ /* error writing out grid */
          printf("%s\n",nc_strerror(istatus));
          return (-1);
        }

       

/* normal return */

        return 0;

}
/************************************************************/
#ifdef __STDC__
void write_cdf_v3 (char *f_filename, char *f_ext, char *f_var, 
                   char *f_comment, char *f_asctime, char *f_cdl_path,
                   char *f_static_path, char *f_ldf, fint4 *ldf_len, 
                   fint4 *fn_length, fint4 *ext_len, 
                   fint4 *var_len, fint4 *comm_len, fint4 *asc_len, 
                   fint4 *cdl_path_len, fint4 *stat_len, fint4 *i_reftime, 
                   fint4 *i_valtime, fint4 *imax, fint4 *jmax, fint4 *kmax, 
                   fint4 *kdim,fint4 lvl[], float *data, float *pr, 
                   fint4 *n_levels, float *cdl_levels,
                   fint4 *called_from, fint4 *append, fint4 *status)
#else
void write_cdf_v3 (f_filename, f_ext, f_var, f_comment, f_asctime, 
                   f_cdl_path, f_static_path, f_ldf,ldf_len, fn_length, 
                   ext_len, var_len, comm_len,  asc_len, cdl_path_len, 
                   stat_len, i_reftime, i_valtime, imax, jmax, kmax, 
                   kdim, lvl, data, pr, n_levels, cdl_levels,
                   called_from, append, status)
char *f_filename; 
char *f_ext; 
char *f_var; 
char *f_comment;
char *f_asctime;
char *f_cdl_path;
char *f_static_path;
char *f_ldf;
fint4 *ldf_len;
fint4 *fn_length; 
fint4 *ext_len; 
fint4 *var_len; 
fint4 *comm_len; 
fint4 *asc_len;
fint4 *cdl_path_len;
fint4 *stat_len;
fint4 *i_reftime; 
fint4 *i_valtime; 
fint4 *imax;
fint4 *jmax; 
fint4 *kmax; 
fint4 *kdim; 
fint4 lvl[]; 
float *data; 
float *pr; 
fint4 *n_levels;
float *cdl_levels;
fint4 *called_from;
fint4 *append;
fint4 *status; 
#endif
{
	int i, j, istatus, cdl_len, cdfid, i_record, var_id;
        int dim_id, comm_id, old_record, int1, int2, name_len;
        int missing_grids, inv_id, lc3_levels, lm1_levels,xyz;
        fint4 i_level, i4time;
        size_t mindex[1], dim_len;
        long num_record;
        float *dptr, Dx, Dy, La1, Lo1, LoV, Latin1, Latin2;
        double timeoff, reftime, valtime, diff;
        static char *fvptr, *vptr, *cptr;
        static char *var, *comment, *filename, *ext; 
        char *map_proj, *origin, *ldf;
        static char *asctime, *syscmd, *cdlfile, *comm_var, *inv_var;
        static char *static_grid;
        char char1[3];

        lc3_levels = (fint4) LC3_LEVELS;
        lm1_levels = (fint4) LM1_LEVELS;

/*      printf("Start of write_cdf_v3 f_filename %s      \n",f_filename);
        printf("Start of write_cdf_v3 f_static_path %s   \n",f_static_path);
        printf("Start of write_cdf_v3 f_ldf %s           \n",f_ldf);
        printf("Start of write_cdf_v3 ldf_len  %d        \n",*ldf_len);
        printf("Start of write_cdf_v3 fn_length  %d      \n",*fn_length);
        printf("Start of write_cdf_v3 cdl_path_len  %d   \n",*cdl_path_len);
        printf("Start of write_cdf_v3 stat_len  %d       \n",*stat_len);
        printf("Start of write_cdf_v3 imax/jmax  %d %d   \n",*imax,*jmax);  */

/* convert fortran f_file_name into C string filename and f_ext into C ext  */
        filename = malloc(*fn_length + 1);
        ldf = malloc(*ldf_len + 1);
        nstrncpy(filename,f_filename,*fn_length);
        nstrncpy(ldf, f_ldf, *ldf_len);
	ext = malloc((*ext_len + 1) * sizeof(char));
        nstrncpy(ext,f_ext,*ext_len);
        downcase_c(ext,ext);

/* allocate space for syscmd and cdlfile and fill up */
        /* cdl file contains "extension" + ".cdl\0" */
	cdlfile = malloc((*cdl_path_len + (*ext_len) + 5) * sizeof(char));
        nstrncpy(cdlfile,f_cdl_path,*cdl_path_len);

/* check to see if ext is v01...v19 */
        if ((strncmp(ext,"v0",2) == 0) || (strncmp(ext,"v1",2) == 0)) {
          strcat(cdlfile,"v00");
        }
        else {
          strcat(cdlfile,ext);
        }

        strcat(cdlfile,".cdl");
        cdl_len = strlen(cdlfile);

/* check to see if cdl file is there */
        if( access(cdlfile, F_OK) != 0 ) {
          printf("The cdl file %s does not exist\n",cdlfile);
          *status = -2; /* error in file creation */
          return;
        }
       
        /* SYSCMD contains "/usr/local/netcdf/bin/ncgen -o %s %s\0" which
           is 33 char, cdlfile, and filename  + 10 extra  */
        syscmd = malloc((strlen(SYSCMD)+cdl_len+*fn_length+10) * sizeof(char));
        sprintf(syscmd,SYSCMD, filename, cdlfile);
        free(cdlfile);
        
/*    see if file is already there  */
        if( access(filename, F_OK) != 0 ) { /* file does not exist */

/*  create file, then open it */
#ifdef SMS
          istatus = pcl_system(syscmd);
#else
          system(syscmd);
#endif
          istatus = nc_open(filename,NC_WRITE, &cdfid);
          if (istatus != NC_NOERR) {
	    *status = -2; /* error in file creation */ 
            printf("The file %s could not be opened\n",filename);
            printf("system command: %s\n",syscmd);
            free_file_var(syscmd, filename);
            free(ext);
            return;
          }
        }
        else { /* file is there...*/
          if ((strncmp(ext,"lga",3) == 0) || (strncmp(ext,"lgb",3) == 0) ||
              (strncmp(ext,"fua",3) == 0) || (strncmp(ext,"fsf",3) == 0)) { /* allow multiple writes to file */
            *called_from = 2;
          }
          else { /* 10/14/97 for write_laps_data with no append, write over it anyway */
            if (((*called_from == 0) || (*called_from == 1)) && (*append == 0)){
#ifdef SMS
              istatus = pcl_system(syscmd);
#else
              system(syscmd); /* added 10/14/97 */
#endif
            }
          }
          istatus = nc_open(filename,NC_WRITE, &cdfid);
          if (istatus != NC_NOERR) {  /* error opening file */
            printf("File %s exists, but cannot be opened.\n",filename);
	    *status = -2; /* error in file creation */ 
            free_file_var(syscmd, filename);
            free(ext);
            return;
          }
        }

        free_file_var(syscmd, filename);

/* set i_record to value returned by nc_inq_dimlen...if value > 0, set to (value - 1)  */

        istatus = nc_inq_dimid(cdfid,"record", &dim_id);
        if (istatus != NC_NOERR) {
	  *status = -2; /* error in file creation */ 
          istatus = nc_close(cdfid);
          free(ext);
          return;
        }

        istatus = nc_inq_dimlen(cdfid,dim_id,&dim_len);
        if (istatus != NC_NOERR) {
	  *status = -2; /* error in file creation */ 
          istatus = nc_close(cdfid);
          free(ext);
          return;
        }
        num_record = (long)dim_len;

        if (num_record == 0)  {   /* no data in file */
          i_record = (int)num_record;
        }
        else {
          if (*append == 0) {
            i_record = 0;
          }
          else {
            i_record = (int)num_record - 1;
          }
          if (*called_from == 2) { /* ok to write into existing record */
            if (*append == 0) { /* if only one analysis allowed in a file */
              if (i_record == 0) {  /* things are OK */
              }
              else {
                *status = -6;
                istatus = nc_close(cdfid);
                free(ext);
                return; 
              }
            }
            else {
              printf("Append option not currently implemented.\n");
              *status = -6;
              istatus = nc_close(cdfid);
              free(ext);
              return; 
            }
          }
          else {
            if (*append == 0) {  /* only one analysis allowed in a file */
              *status = -6;
              istatus = nc_close(cdfid);
              free(ext);
              return; 
            }
            else {
              printf("Append option not currently implemented.\n");
              printf("location 3\n");
              *status = -6;
              istatus = nc_close(cdfid);
              free(ext);
              return; 
            }
          }
        }

/* determine value of namelen dimension in output file  */

        printf("write_cdf_v3 location 3  \n");

        istatus = nc_inq_dimid(cdfid,"namelen", &dim_id);
        if (istatus != NC_NOERR) {
	  *status = -2; /* error in file creation */ 
          istatus = nc_close(cdfid);
          free(ext);
          return;
        }

        printf("write_cdf_v3 location 3a  \n");

        istatus = nc_inq_dimlen(cdfid,dim_id,&dim_len);
        if (istatus != NC_NOERR) {
	  *status = -2; /* error in file creation */ 
          istatus = nc_close(cdfid);
          free(ext);
          return;
        }
        name_len = (int)dim_len;

        printf("write_cdf_v3 location 3b \n");

/* file is now open and ready for writing.    */
/* verify x and y in output file match imax and jmax */
/* verify that z is correct: 
     n_levels is the number of 3D levels listed in nest7grid.parms.
     options  for z are: n_levels (standard LAPS 3D file), z = 1 for surface,
                         and 2 special cases, LC3 = 42, LM1 = 3.
     set n_levels to match one of the four cases above, or abort   */

        istatus = check_grid_dims_v3(cdfid, ext, imax, jmax, n_levels,
                                        lc3_levels, lm1_levels);

        printf("write_cdf_v3 location 3c \n");

        if (istatus == -1) {
          *status = -3;
          free(ext);
          return;
        }

/* write header info only if num_record == 0 */

        printf("write_cdf_v3 location 4  \n");

        if (num_record == 0) {

          map_proj = malloc(name_len * sizeof(char));
          origin = malloc(name_len * sizeof(char));

          if ((strncmp(ext,"fua",3) == 0) || (strncmp(ext,"fsf",3) == 0)  ||
              (strncmp(ext,"pbl",3) == 0)) {
/* don't read static.nest7grid...data should be in cdl file */
          }
          else {
/* get info from static.nest7grid for writing to output files */
	    static_grid = malloc((*stat_len + 20) * sizeof(char));
            nstrncpy(static_grid,f_static_path,*stat_len);
            strcat(static_grid,"static.");
            strcat(static_grid,ldf);

            istatus = get_static_info_v3(static_grid, &Dx, &Dy, &La1, &Lo1, 
                                       &LoV, &Latin1, &Latin2, map_proj, 
                                       origin, name_len, ldf);
            free(static_grid);

            if (istatus == -1) {
              printf("Error reading info from static.%s\n", ldf);
              printf("  Some navigation data will be missing from file.\n");
            }
          }

/* write out header information into file */

	  istatus = write_hdr_v3(cdfid, ext, i_record, imax, jmax, kmax, kdim, 
                                 pr, n_levels, &Dx, &Dy, &La1, &Lo1, 
                                 &LoV, &Latin1, &Latin2, map_proj, origin,name_len,
                                 cdl_levels);

          free_static_var(map_proj, origin);
          if (istatus == -1) {
	    *status = -5; /* error writing out header */ 
            istatus = nc_close(cdfid);
            free(ext);
            free_write_var(var, comment, asctime, comm_var, inv_var);
            return;
          }
        } /* end if num_record == 0  */

        printf("write_cdf_v3 location 5  \n");
       
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

        if (*called_from == 2) {
          if (num_record == 0) {
            old_record = -1;
          }
          else {
            old_record = i_record;
          }
        }
        else { 
          old_record = -1;
        }
	missing_grids = 0;
        vptr = var;

/*      printf("write_cdf_v3 location 6  \n"); */

	for (i = 0; i < *kdim; i++) {

/* handle special case of i_record for lmr and lf1 files */
          if ((strcmp(ext,"lmr") == 0) || (strcmp(ext,"lf1") == 0)) {
            fvptr = vptr;  /* fvptr temp point to current var */
            char1[0] = *(vptr+1);
            char1[1] = *(vptr+2);
            char1[2] = '\0';

            int1 = atoi(char1);
            fvptr = vptr;
            fvptr++;
            *fvptr = '\0';
            
            timeoff= (double) int1*600;
            if (old_record == (-1)) diff = valtime - reftime;
            
            valtime = reftime + timeoff + diff;

            i4time = valtime + 315619200;
/* commented out cv_i4tim_asc_lp call 4-5-00...will replicate f_asctime passed
   in to call to write_cdf_v3 in 1hr and 2hr forecasts */
/*            cv_i4tim_asc_lp(&i4time,f_asctime,&istatus);  */
            fstrncpy(asctime, f_asctime,*asc_len);
          }

/* write out valtime and reftime if old_record != i_record */
          if (i_record != old_record) {
            istatus = write_val_ref_asctime(cdfid, i_record, &valtime, &reftime, 
                          asctime, asc_len, name_len);
            if (istatus == -1) {
	      *status = -5; /* error writing out header */ 
              istatus = nc_close(cdfid);
              free(ext);
              free_write_var(var, comment, asctime, comm_var, inv_var);
              return;
            }
            old_record = i_record;
          }

          sprintf(comm_var,"%s_comment",vptr);
          sprintf(inv_var,"%s_fcinv",vptr);
	  istatus = nc_inq_varid(cdfid,vptr,&var_id);
          if (istatus != NC_NOERR) {
            printf("Variable %s not found in output file.\n",vptr);
            var_id = -1;
          }

          istatus = nc_inq_varid(cdfid,comm_var,&comm_id);
          if (istatus != NC_NOERR) {
            printf("Variable %s not found in output file.\n",comm_var);
            comm_id = -1;
          }
            
          istatus = nc_inq_varid(cdfid,inv_var,&inv_id);
          if (istatus != NC_NOERR) {
            printf("Variable %s not found in output file.\n",inv_var);
            inv_id = -1;
          }
            
          if ((var_id == (-1)) || (comm_id == (-1)) || (inv_id == (-1))) {
            printf("Grid %s not written to file.\n",vptr);
            missing_grids++;
          }
          else {
            i_level = lvl[i];
            dptr = (data + (i*(*imax)*(*jmax)));
            cptr = (comment + (i*(*comm_len + 1)));
            istatus = update_laps_v3(cdfid,i_level,i_record,imax,
                                     jmax,var_id,inv_id,dptr,
                                     comm_id,cptr,name_len);
            if (istatus == (-1)) {
              missing_grids++;
              printf("Variable with error is %s at level %d.\n",vptr,i_level);
            }

            if ((strcmp(ext,"lmr") == 0) || (strcmp(ext,"lf1") == 0))
              i_record++;
          }
          vptr += (*var_len + 1);
        }

/*      printf("write_cdf_v3 location 7  \n");  */

        istatus = nc_close(cdfid);
        if(DEBUG==1) printf("normal close in cdf_write %d\n",istatus);

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
#ifdef __STDC__
int open_cdf (int mode, char *fname,fint4 *no_laps_diag)
#else
int open_cdf (mode, fname, no_laps_diag)
int mode;
char *fname;
fint4 *no_laps_diag;
#endif
{
        int      cdfid, istatus, i;

/* open the netcdf file and get the file id and the variable id */
        istatus = nc_open (fname, mode, &cdfid);

        if (istatus != NC_NOERR) {
          if (*no_laps_diag == 0) {
            printf("open_cdf: cannot open file as netCDF %s.\n", fname); 
          }
          return -1;
        }
        else {
          return cdfid;
        }
}
/************************************************************/

