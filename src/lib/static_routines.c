#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "netcdf.h"
#include "grid_info.h"
#include "laps_grid_def.h"

#define SYSCMD "ncgen -o %s %s"

#ifdef FORTRANCAPS
#define read_cdf_static READ_CDF_STATIC
#define write_cdf_static WRITE_CDF_STATIC
#endif

#ifdef FORTRANUNDERSCORE
#define cre_static cre_static_
#define cre_loss cre_loss_
#define cdf_update_stat cdf_update_stat_
#define write_cdf_static write_cdf_static_
#define cdf_retr_grid_stat cdf_retr_grid_stat_
#define cdf_retr_hdr_stat cdf_retr_hdr_stat_
#define read_cdf_static read_cdf_static_
#define open_cdf open_cdf_
#define itoa itoa_
#define cdf_dim_size cdf_dim_size_
#define log_diag log_diag_
#define nstrncpy nstrncpy_
#define fstrncpy fstrncpy_
#endif

/*****************************************************************************
* CDF_UPDATE_STAT
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Update the LAPS netCDF static file 
*	Purpose			To coordinate the writing of static info into
*				  the appropriate netCDF file: find the 
*				  appropriate location to write the grid, 
*				  and call the routine that does the writing.
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 9/5/90
*			adapted for new WRITELAPSDATA 1/93 Linda Wharton
*			adapted from CDF_UPDATE_LAPS for new 
*			  WRT_LAPS_STATIC 4/93 Linda Wharton
*
*	Input :
*		s_field		A character string identifying the field:
*				    lat		latitude of grid point
*				    lon		longitude of grid point
*				    avg		average elev MSL of grid box
*				    std		std dev of elev's in box
*				    env  	avg + std
*				    zin		avg converted to mb pressure,
*						  then mapped 1100=0 to 100=20
*		gptr		A pointer to the grid (including product and data
*						headers).
*	Output :
*		None
*	Globals :		
*		NC_WRITE	A flag used to open the netCDF file for writing
*	Returns :
*		 0 if successful
*		-1 if an error occurs
***************************************************************************/
#ifdef __STDC__
int cdf_update_stat (int i_cdfid,char *s_field,char *gptr,char *commnt,
		     char *comm_ptr)
#else
int cdf_update_stat (i_cdfid, s_field, gptr, commnt, comm_ptr)
int i_cdfid;
char *s_field;
char *gptr;
char *commnt;
char *comm_ptr;
#endif
{
	int i_status, i_comid, i_varid;
#ifdef __alpha
	long start[4], count[4], start_c[3], count_c[3];
#else
	int start[4], count[4], start_c[3], count_c[3];
#endif

	log_diag (2, "cdf_update_stat:grid name = %s\n", s_field);

/* turn off the error handling done by the netCDF routines */
        ncopts = NC_VERBOSE;

/* get the x and y dimension sizes */
        count[0] = 1;
        count[1] = 1;
	count[2] = cdf_dim_size (i_cdfid, "y");
	count[3] = cdf_dim_size (i_cdfid, "x");

/* construct the arrays needed by the netcdf write routine */
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	start[3] = 0;

	start_c[0] = 0;
	start_c[1] = 0;
	start_c[2] = 0;
	count_c[0] = 1;
	count_c[1] = 1;
	count_c[2] = strlen(comm_ptr);

/* get the variable ids */
	log_diag (2, "Data variable name = %s\n", s_field);
	if ((i_varid = ncvarid (i_cdfid, s_field)) == (-1))
		return -1;
	if ((i_comid = ncvarid (i_cdfid, commnt)) == (-1))
		return -1;

/* write the grid to the netcdf file */
	i_status = ncvarput (i_cdfid, i_varid, start, count, gptr);
	
	if (i_status == (-1)) {
		log_diag (1, "cdf_update_stat: error during cdf write\n");
		return -1;
	}
	else {
	  i_status = ncvarput (i_cdfid, i_comid, start_c, count_c, comm_ptr);
	  
	  if (i_status == (-1))
		log_diag (1, "cdf_update_stat: error during cdf write\n");
	  else
		log_diag (2, "cdf_update_stat: cdf write ok\n");
	}

	return i_status;
}
/*************************************************************************
*	WRITE_CDF_STATIC
*	Category	Product Management
*	Group		General Purpose Database
*	Module		Write Grid data
*	Purpose		Write grid data into netCDF static file.
*
*	Designer/Programmer : Linda Wharton
*	Modifications : original 4/93
*
*	Input:
*		filname		NetCDF filename to open
*		asctime		Ascii time of file
*		var		Array of LAPS variables to write out
*		comment		Comments for each level to write
*		imax		X dimension of data in file
*		jmax		Y dimension of data in file
*		n_grids		Number of grids available in file & 
*		 		  dimension  of the VAR and COMMENT
*				  in the calling program
*		data		3D grid containing data to write
*		zin		Grid for AVS applications
*		model		Model generating static data
*		grid_spacing	size of grid box in M
*	Output:
*		status		Returns status to calling subroutine
*	Globals:
*		none
*	Returns:
*		none
*****************************************************************************/
#ifdef __STDC__
#ifdef __alpha
void write_cdf_static(char *filname,short *s_length,char *f_asctime,
                      char *f_dir, int *dir_len,
		      char *f_var,char *f_comment,char *f_ldf,
                      int *ldf_len,int *imax,int *jmax, int *n_grids,
                      float *data,float *zin,char *f_model,
		      float *grid_spacing,int *status)
#else
void write_cdf_static(char *filname,short *s_length,char *f_asctime,
                      char *f_dir, long *dir_len,
		      char *f_var,char *f_comment,char *f_ldf,
                      long *ldf_len,long *imax,long *jmax, long *n_grids,
                      float *data,float *zin,char *f_model,
		      float *grid_spacing,long *status)
#endif
#else
#ifdef __alpha
void write_cdf_static(filname,s_length,f_asctime,f_dir,dir_len,
                      f_var,f_comment,f_ldf,ldf_len,imax,jmax,
                      n_grids,data,zin,f_model, grid_spacing,
                      status)
char *filname;
short *s_length;
char *f_asctime;
char *f_dir;
int *dir_len;
char *f_var;
char *f_comment;
char *f_ldf;
int *ldf_len;
int *imax;
int *jmax;
int *n_grids;
float *data;
float *zin;
char *f_model;
float *grid_spacing;
int *status;
#else
void write_cdf_static(filname,s_length,f_asctime,f_dir,dir_len,
                      f_var,f_comment,f_ldf,ldf_len,imax,jmax,
                      n_grids,data,zin,f_model, grid_spacing,
                      status)
char *filname;
short *s_length;
char *f_asctime;
char *f_dir;
long *dir_len;
char *f_var;
char *f_comment;
char *f_ldf;
long *ldf_len;
long *imax;
long *jmax;
long *n_grids;
float *data;
float *zin;
char *f_model;
float *grid_spacing;
long *status;
#endif
#endif
{
	char prefix[5];
	char comm_var[13];                
	int out_file, istat, i_varid, i, process_vr, count;
        char model[132],asctime[18],var[8][4],comment[8][126];
	char fname[92],*ldf;
	char zin_comment[126];
	long start[1], edges[1];
	int grid_done[8]; /* 0=LAT,1=LON,2=AVG,3=ZIN,4=STD,5=ENV,6=LDF,7=USE */
        static char *syscmd, *cdlfile;
        int cdl_len;

/* turn off the error handling done by the netCDF routines */
        ncopts = NC_VERBOSE;

/****  Header info set in cre_static subroutine
	   1)  Record length (IMAX), 
	   2)  Number of lines (JMAX),
	   3)  Version number of write static routine.
     
       Header info set when header is written
	   1)  Number of grids written (n_grids),
	   2)  Ascii time,
	   3)  Model generating topo data.
****/

	for (i = 0; i < 8; i++)
	  grid_done[i] = 0;
	for (i = 0; i < 126; i++)
	   zin_comment[i] = '\0';
	  
/* convert fortran string f_var to c string var */
 
        for (i = 0; i < *n_grids; i++) {
          nstrncpy(var[i],(f_var+i*3),3);
          fstrncpy(comment[i],(f_comment+i*125),125);
        }
        fstrncpy(model,f_model,131);
        fstrncpy(asctime,f_asctime,17);

/* convert fortran file_name into C fname, and fortran f_ext into C ext  */
        nstrncpy(fname,filname,s_length);

/* make sure grids to write will fit in file */
        if (*imax <= NX && *jmax <= NY && *n_grids <= 7)
          {}
        else {
          *status = -3;	/* returns dimension error */
          return;
        }

/* allocate space for syscmd and cdlfile and fill up */
        /* cdl file contains domain name + ".cdl\0" */
        cdlfile = malloc(((*dir_len)+(*ldf_len)+5) * sizeof(char));
        nstrncpy(cdlfile,f_dir,*dir_len);
        ldf = malloc((*ldf_len + 1) * sizeof(char));
        nstrncpy(ldf,f_ldf,*ldf_len);
        strcat(cdlfile,ldf,*ldf_len);
        free(ldf);
        strcat(cdlfile,".cdl");
        cdl_len = strlen(cdlfile);
 
        /* SYSCMD contains "/usr/local/netcdf/bin/ncgen -o %s %s\0" which
           is 33 char, cdlfile, and fname  + 10 extra  */
        syscmd = malloc((strlen(SYSCMD)+cdl_len+*s_length+10) * sizeof(char));
        sprintf(syscmd,SYSCMD, fname, cdlfile);
        free(cdlfile);
 
/* create output file */
        system(syscmd);
        free(syscmd);
        out_file = ncopen(fname,NC_WRITE);
        if (out_file == -1) {
          *status = -2; /* error in file creation */
          return;
        }

 /*  update header */
    start[0] = 0;
  /* store n_grids */
    if ((i_varid = ncvarid (out_file, "n_grids")) == (-1)) {
      *status = -5;	/* returns "error writing header" */
      return;
    }
#ifdef __alpha
    ncvarput1(out_file, i_varid, (long *)0, (void *)n_grids);
#else
    ncvarput1(out_file, i_varid, (int *)0, (void *)n_grids);
#endif

  /* store grid_spacing */
    if ((i_varid = ncvarid (out_file, "grid_spacing")) == (-1)) {
      *status = -5;	/* returns "error writing header" */
      return;
    }
#ifdef __alpha
    ncvarput1(out_file, i_varid, (long *)0, (void *)grid_spacing);
#else
    ncvarput1(out_file, i_varid, (int *)0, (void *)grid_spacing);
#endif
      
  /* store asctime */
    edges[0] = 18;
    if ((i_varid = ncvarid (out_file, "asctime")) == (-1)) {
      *status = -5;	/* returns "error writing header" */
      return;
    }
    ncvarput(out_file, i_varid, start, edges, (void *)asctime);
 
  /* store model */
    edges[0] = 132;
    if ((i_varid = ncvarid (out_file, "process_name")) == (-1)) {
      *status = -5;	/* returns "error writing header" */
      return;
    }
    ncvarput(out_file, i_varid, start, edges, (void *)model);
 
/* write data to grid and write out comment */
	count = 0;
        for (i = 0; i < *n_grids; i++) {
          process_vr = 1;
          if (strncmp(var[i], "LAT",3) ==  0) {
            strcpy(prefix,"lat");
            strcpy(comm_var,"lat_comment");
            grid_done[0] = 1;
          }
          else if (strncmp(var[i], "LON",3) ==  0) {
            strcpy(prefix,"lon");
            strcpy(comm_var,"lon_comment");
            grid_done[1] = 1;
          }
          else if (strncmp(var[i], "AVG",3) ==  0) {
            strcpy(prefix,"avg");
            strcpy(comm_var,"avg_comment");
            grid_done[2] = 1;
          }
          else if (strncmp(var[i], "STD",3) ==  0) {
            strcpy(prefix,"std");
            strcpy(comm_var,"std_comment");
            grid_done[4] = 1;
          }
          else if (strncmp(var[i], "ENV",3) ==  0) {
            strcpy(prefix,"env");
            strcpy(comm_var,"env_comment");
            grid_done[5] = 1;
          }
          else if (strncmp(var[i], "LDF",3) ==  0) {
            strcpy(prefix,"ldf");
            strcpy(comm_var,"ldf_comment");
            grid_done[6] = 1;
          }
          else if (strncmp(var[i], "USE",3) ==  0) {
            strcpy(prefix,"use");
            strcpy(comm_var,"use_comment");
            grid_done[7] = 1;
          }
          else  {
            process_vr = 0;
            count += 1;
          }
          
          if (process_vr == 1) {
            istat = cdf_update_stat(out_file,prefix,
               	 		    (data + (i*(*imax)*(*jmax))),
                	 	    comm_var,comment[i]);
            if (istat == -1) {
              *status = -4;
              ncclose(out_file);
              return;
            }
          }
        }
   
/* write out zin if LAT, LON and AVG were written */
        if (grid_done[0]==1 && grid_done[1]==1 && grid_done[2]== 1) {
            strcpy(zin_comment,
"Grid contains AVG data converted from meters to pressure in mb, then\
 mapped from 1100mb=0 to 100mb=20");
            strcpy(prefix,"zin");
            strcpy(comm_var,"zin_comment");
            istat = cdf_update_stat(out_file,prefix,zin,comm_var,
            			    zin_comment);
          grid_done[3] = 1;
        }  
        
        if (grid_done[0]==1 && grid_done[1]==1 && grid_done[2]== 1
            && grid_done[3]==1) {
          *status = *n_grids - count;  /* returns number of grids written */
        }
        else
          *status = -6;  /* returns missing LAT,LON or AVG grids error  */
        
        ncclose(out_file);          	
        return;
}
/*****************************************************************************
*	CDF_RETR_GRID_STAT
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Retrieve Grid Data
*	Purpose			To retrieve a specified grid from a netcdf file
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 9/12/90
*			modified for new READLAPSDATA 1/93 Linda Wharton
*
*	Input :
*		i_cdfid		NetCDF file id for file to be read
*		s_field		A character string identifying the field
*		gptr		Pointer to location in data array to write data.
*		cptr		Pointer to location in data array to write comments
*		uptr		Pointer to location in data array to write units.
*	Output :
*		none
*	Globals:
*		NC_NOWRITE	A NetCDF global that opens the NetCDF file for reading
*						only.
*	Returns:
*		0 if normal return
*		-1 if an error occurs
*****************************************************************************/

#ifdef __STDC__
int cdf_retr_grid_stat(int i_cdfid, char *s_field, char version[], float *gptr,
		       char *cptr, char *uptr)
#else
int cdf_retr_grid_stat(i_cdfid, s_field, version, gptr, cptr, uptr)

int i_cdfid;
char *s_field;
char version[];
float *gptr;
char *cptr;
char *uptr;
#endif
{
	int i_status, j, i_varid, ver;
#ifdef __alpha
	long start[4],count[4],start_c[3],count_c[3];
	long v2_start[2],v2_count[2],v2_start_c[1],v2_count_c[1];
#else
	int start[4],count[4],start_c[3],count_c[3];
	int v2_start[2],v2_count[2],v2_start_c[1],v2_count_c[1];
#endif
	char var_name[13];

/*	printf("cdf_ret_grd: level = %d   fctime = %d   field = %s\n",
		i_level, i_fctime, s_field);  */

/* turn off the error handling done by the netCDF routines */
	ncopts = NC_VERBOSE;

/* get the x and y dimension sizes */
        i_status = strncmp(version,"V2",2);
        if (i_status == 0) {
          ver = 2;
	  v2_count[0] = cdf_dim_size (i_cdfid, "lat");
	  v2_count[1] = cdf_dim_size (i_cdfid, "lon");
	  v2_start[0] = 0;
	  v2_start[1] = 0;
        }
        else {
          ver = 3;
          count[0] = 1;
          count[1] = 1;
          count[2] = cdf_dim_size (i_cdfid, "y");
          count[3] = cdf_dim_size (i_cdfid, "x");
          start[0] = 0;
          start[1] = 0;
          start[2] = 0;
          start[3] = 0;
        }

/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", s_field);
	if ((i_varid = ncvarid (i_cdfid, s_field)) == (-1)) {
		printf("cdf_retrieve_laps: no grid available.\n");
		return -1;
	}
			   

/* read the grid from the netcdf file */
        if (ver == 2)
	  i_status = ncvarget (i_cdfid, i_varid, v2_start, v2_count, gptr);
        else
	  i_status = ncvarget (i_cdfid, i_varid, start, count, gptr);

	if (i_status == (-1)) {
	   printf("cdf_retrieve_laps: error retrieving data %s grid.\n", 
		 	 *s_field);
	   return -1;
	}

/* retrieve units */
	for (j = 0; j < 11; j++) 
	  *(uptr + j) = '\0';
	i_status = ncattget (i_cdfid, i_varid, "LAPS_units", uptr);
	if (i_status == (-1)) {
	   printf("cdf_retrieve_laps: error retrieving LAPS_units.\n");
	   return -1;
	}

/* setup to read the comment from the netcdf file */

	sprintf(var_name, "%s%s", s_field, "_comment");	
	if ((i_varid = ncvarid (i_cdfid, var_name)) == (-1)) {
	   printf("cdf_retrieve_laps: no comment field available.\n");
	   return -1;
	}
	for (j = 0; j < 132; j++) 
	  *(cptr + j) = '\0';

/* construct the arrays needed to read the grid */
        if (ver == 2) {
	  v2_start_c[0] = 0;
	  v2_count_c[0] = 126;
	  i_status = ncvarget (i_cdfid, i_varid, v2_start_c, v2_count_c, cptr);
        }
        else {
	  start_c[0] = 0;
	  start_c[1] = 0;
	  start_c[2] = 0;
	  count_c[0] = 1;
	  count_c[1] = 1;
	  count_c[2] = 126;
	  i_status = ncvarget (i_cdfid, i_varid, start_c, count_c, cptr);
        }

	if (i_status == (-1)) {
	   printf("cdf_retrieve_laps: error retrieving comment.\n");
	   return -1;
	}
	   
/* normal return */

	return 0;
}
/*****************************************************************************
*	CDF_RETR_HDR_STAT
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Retrieve Static Grid
*	Purpose			To retrieve header info from a netcdf file
*
*	Designer/Programmer : Linda Wharton
*	Modifications : original 4/93
*
*	Input :
*		i_cdfid		NetCDF file id for file to be read
*	Output :
*		imaxn		number of x gridpoints
*		jmaxn		number of y gridpoints
*		kmaxn		number of grids in the file
*		laps_dom_file	name of laps domain for this file
*		asctime		ascii time netCDF file was written
*		version		version of WRITELAPSDATA used to make grid
*		model		meteorological model in file
*		origin		location where file was created
*		num_variables   number of variables in this file
*	Globals:
*		NC_NOWRITE	A NetCDF global that opens the NetCDF file for
*				    reading only.
*	Returns:
*		-1 if an error occurs
*****************************************************************************/

#ifdef __STDC__
#ifdef __alpha
int cdf_retr_hdr_stat(int i_cdfid,int *imaxn, int *jmaxn, int *n_grids_n,
		      float *grid_spacing_n,char *asctime, char *version,
		      char *model, char *origin, long *num_variables)
#else
int cdf_retr_hdr_stat(int i_cdfid,long *imaxn, long *jmaxn, long *n_grids_n,
		      float *grid_spacing_n,char *asctime, char *version,
		      char *model, char *origin, long *num_variables)
#endif
#else
#ifdef __alpha
int cdf_retr_hdr_stat(i_cdfid,imaxn,jmaxn,n_grids_n,grid_spacing_n,asctime, 
		      version,model,origin,num_variables)
int i_cdfid;                               
int *imaxn;
int *jmaxn;
int *n_grids_n;
float *grid_spacing_n;
char *asctime;
char *version;
char *model;
char *origin;
long *num_variables;
#else
int cdf_retr_hdr_stat(i_cdfid,imaxn,jmaxn,n_grids_n,grid_spacing_n,asctime, 
		      version,model,origin,num_variables)
int i_cdfid;                               
long *imaxn;
long *jmaxn;
long *n_grids_n;
float *grid_spacing_n;
char *asctime;
char *version;
char *model;
char *origin;
long *num_variables;
#endif
#endif
{
	int i_status, i_varid, i, i_version, temp, str_len;
#ifdef __alpha
        long mindex[1], start[1], count_asc[1], count_long[1];
#else
        int mindex[1], start[1], count_asc[1], count_long[1];
#endif
	static char c_ver[5];
	char *t_ptr;

/* turn off the error handling done by the netCDF routines */
	ncopts = NC_VERBOSE;
	mindex[0] = 0;
	start[0] = 0;
	count_asc[0] = 18;
	count_long[0] = 132;
	
/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", "imax");
	if ((i_varid = ncvarid (i_cdfid, "imax")) == (-1))
		return -1;

/* read the var from the netcdf file */
	i_status = ncvarget1 (i_cdfid, i_varid, mindex, imaxn);
	if (i_status == (-1))
	   return -1;

/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", "jmax");
	if ((i_varid = ncvarid (i_cdfid, "jmax")) == (-1))
		return -1;

/* read the var from the netcdf file */
	i_status = ncvarget1 (i_cdfid, i_varid, mindex, jmaxn);
	if (i_status == (-1))
	   return -1;

/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", "n_grids");
	if ((i_varid = ncvarid (i_cdfid, "n_grids")) == (-1))
		return -1;

/* read the var from the netcdf file */
	i_status = ncvarget1 (i_cdfid, i_varid, mindex, n_grids_n);
	if (i_status == (-1))
	   return -1;

/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n",
		  "grid_spacing");
	if ((i_varid = ncvarid (i_cdfid, "grid_spacing")) == (-1))
		return -1;

/* read the var from the netcdf file */
	i_status = ncvarget1 (i_cdfid, i_varid, mindex, grid_spacing_n);
	if (i_status == (-1))
	   return -1;

/* read the version from the netcdf file */
        if (ncattget (i_cdfid, NC_GLOBAL, "version", (void *)&i_version) == -1) {
          i_varid = ncvarid(i_cdfid,"version");
          if (i_varid == -1) {
            printf("error reading version");
            return -1;
          }
          else {
            if (ncvarget1(i_cdfid, i_varid, mindex, (void *)&i_version) == -1) {
              return -1;
            }
          }
        }

/* convert integer value of version to character string */
	if (i_version < 10) {
	   strcpy(version,"V");
	   itoa(i_version,c_ver,2);
	   strcat(version,c_ver);
	}
	else if (i_version < 100) {
	   strcpy(version,"V");
	   itoa(i_version,c_ver,3);
	   strcat(version,c_ver);
	}
	else {
	   strcpy(version,"V");
	   itoa(i_version,c_ver,4);
	   strcat(version,c_ver);
	}
	
/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", "asctime");
	if ((i_varid = ncvarid (i_cdfid, "asctime")) == (-1))
		return -1;
	   
/* read the asctime from the netcdf file */
	i_status = ncvarget (i_cdfid, i_varid, start, count_asc, asctime);
	if (i_status == (-1))
	   return -1;

/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", "model");
	if ((i_varid = ncvarid (i_cdfid, "process_name")) == (-1)) {
	  i_varid = ncvarid (i_cdfid, "model");
          if (i_varid == -1) return -1;
        }
	   
/* read the model from the netcdf file */
	i_status = ncvarget (i_cdfid, i_varid, start, count_long, model);
	if (i_status == (-1)) 
	   return -1;
	else {
	   str_len = strlen(model);
	   if (str_len < 131) {
	      for (i = str_len+1; i < 132; i++)
	         strcat(model," ");
	   }
        }

/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", "origin");
	if ((i_varid = ncvarid (i_cdfid, "origin_name")) == (-1)) {
          i_varid = ncvarid (i_cdfid, "origin");
          if (i_varid == -1) return -1;
        }

/* read the asctime from the netcdf file */
	i_status = ncvarget (i_cdfid, i_varid, start, count_long, origin);
	if (i_status == (-1))
	   return -1;
	else {
	   str_len = strlen(origin);
	   if (str_len < 131) {
	      for (i = str_len+1; i < 132; i++)
	         strcat(origin," ");
	   }
        }

/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", 
		  "num_variables");
	if ((i_varid = ncvarid (i_cdfid, "num_variables")) == (-1)) {
          i_varid = ncvarid (i_cdfid, "n_laps_var");
          if (i_varid == -1) return -1;
        }

/* read the var from the netcdf file */
	  i_status = ncvarget1 (i_cdfid, i_varid, mindex, num_variables);
	if (i_status == (-1))
	   return -1;
   	   
/* normal return */

	return 0;
}

/*************************************************************************
*	READ_CDF_STATIC
*	Category	Product Management
*	Group		General Purpose Database
*	Module		write_cdf-Read_cdf_file
*	Purpose		Read LAPS data from netCDF format static file.
*
*	Designer/Programmer : Linda Wharton
*	Modifications : original 4/93
*
*	Input:
*		iimax      	Expected X dimension of data
*		jjmax      	Expected Y dimension of data
*		kkmax      	Expected number of data grids to retrieve
*		kdim 		Dimension  of the var_req, lvl_req, comment,
*				  lvl_coord, units in the calling program
*		var_req		Array of LAPS variables to retrieve from the file
*		lvl_req		Levels of LAPS variables to retrieve
*		fname		NetCDF filename to open
*		ext		Extension identifying data to be read
*	Output:
*		imax		Actual X dimension of data in file
*		jmax		Actual Y dimension of data in file
*		kmax		Actual number of grids available in file
*		lvl_coord	Level coordinates of each level read
*		units		Units of each level read
*		comment		Comments for each level read
*		asctime		Ascii time of file
*		version		Version of WRITELAPSDATA that wrote the file
*		data		3D grid containing data requested
*		process_var	Array containing 1 if variable was processed,
*				  otherwise contains 0
*		status		Returns status to calling subroutine
*	Globals:
*		none
*	Returns:
*		count of requested variables that were not retrieved
*************************************************************************/
#ifdef __STDC__
#ifdef __alpha
void read_cdf_static(char *filname, short *s_length,char *f_var,
		     char *f_comment,char *f_units,int *imax,int *jmax,
		     int *n_grids,float *data,float *grid_spacing,
		     int *no_laps_diag,int *status)
#else
void read_cdf_static(char *filname, short *s_length,char *f_var,
		     char *f_comment,char *f_units,long *imax,long *jmax,
		     long *n_grids,float *data,float *grid_spacing,
		     long *no_laps_diag,long *status)
#endif
#else
#ifdef __alpha
void read_cdf_static(filname,s_length,f_var,f_comment,f_units,imax,jmax,
		     n_grids,data,grid_spacing,no_laps_diag,status)
char *filname;
short *s_length;
char *f_var;
char *f_comment;
char *f_units;
int *imax;
int *jmax;
int *n_grids;
float *data;
float *grid_spacing;
int *no_laps_diag;
int *status;
#else
void read_cdf_static(filname,s_length,f_var,f_comment,f_units,imax,jmax,
		     n_grids,data,grid_spacing,no_laps_diag,status)
char *filname;
short *s_length;
char *f_var;
char *f_comment;
char *f_units;
long *imax;
long *jmax;
long *n_grids;
float *data;
float *grid_spacing;
long *no_laps_diag;
long *status;
#endif
#endif
{		   
	int istat, i, unconv_var, process_vr, cdfid;
        char var[8][4],comment[8][126],units[8][11];
	char prefix[5];
#ifdef __alpha
	int imaxn, jmaxn, n_grids_n;
#else
	long imaxn, jmaxn, n_grids_n;
#endif
	char model[132], asctime[18], origin[132], version[5];
	long num_variables;
	char fname[92];
	
/* turn off the error handling done by the netCDF routines */
        ncopts = NC_VERBOSE;

/* null out arrays for units,comment before using   */
 
        for (i = 0; i < *n_grids; i++) {
          strncpy(units[i],"          ",10);
          strncpy(comment[i],"                                                                                                                             ",125);
        }

/* convert fortran string f_var to c string var */
 
        for (i = 0; i < *n_grids; i++)
          nstrncpy(var[i],(f_var+i*4),3);
 
/* convert fortran filname into C fname */
        nstrncpy(fname,filname,s_length);
 
/* open file for reading */	
	cdfid = open_cdf(NC_NOWRITE,fname,no_laps_diag);
	if (cdfid == -1) {
		*status = -1;	/* error opening file */
		return;
	}

/* get header info */
	istat = cdf_retr_hdr_stat(cdfid,&imaxn,&jmaxn,&n_grids_n,
				  grid_spacing,asctime,version,
				  model,origin,&num_variables);
	if (istat == -1) {
	   *status = -5;
	   ncclose(cdfid);
	   return;
	}
	if (imaxn > *imax || jmaxn > *jmax || *n_grids > n_grids_n) {
	   *status = -3;
	   ncclose(cdfid);
	   return;
	}
	
	unconv_var = 0;
	for (i = 0; i < *n_grids; i++) {
          process_vr = 1;
          if (strncmp(var[i], "LAT",3) ==  0) {
            strcpy(prefix,"lat");
          }
          else if (strncmp(var[i], "LON",3) ==  0) {
            strcpy(prefix,"lon");
          }
          else if (strncmp(var[i], "AVG",3) ==  0) {
            strcpy(prefix,"avg");
          }
          else if (strncmp(var[i], "STD",3) ==  0) {
            strcpy(prefix,"std");
          }
          else if (strncmp(var[i], "ENV",3) ==  0) {
            strcpy(prefix,"env");
          }
          else if (strncmp(var[i], "ZIN",3) ==  0) {
            strcpy(prefix,"zin");
          }
          else if (strncmp(var[i], "LDF",3) ==  0) {
            strcpy(prefix,"ldf");
          }
          else if (strncmp(var[i], "USE",3) ==  0) {
            strcpy(prefix,"use");
          }
          else  {
            process_vr = 0;
	    unconv_var += 1;
          }
          
          if (process_vr == 1) {
   	    istat = cdf_retr_grid_stat(cdfid,prefix,version,
	   		               (data + i*(*imax)*(*jmax)),
				       comment[i],units[i]);
            if (istat == -1) {
              *status = -4;
              ncclose(cdfid);
              return;
            }
          }
        }
   	
	ncclose(cdfid);
	
        for (i = 0; i < *n_grids; i++) {
          fstrncpy((f_comment+(i*126)),comment[i],125);
          fstrncpy((f_units+(i*11)),units[i],10);
        }

	*status = unconv_var;
	return;
}
                                                        

