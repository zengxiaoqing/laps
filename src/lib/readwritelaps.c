/*cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis  
cdis 
cdis*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "netcdf.h"
#include "grid_info.h"
#include "laps_grid_def.h"


#ifdef FORTRANUNDERSCORE
#define cre_lw3 cre_lw3_
#define cre_lh1 cre_lh1_
#define cre_lh2 cre_lh2_
#define cre_lh3 cre_lh3_
#define cre_lh4 cre_lh4_
#define cre_lq3 cre_lq3_
#define cre_lsx cre_lsx_
#define cre_lwm cre_lwm_
#define cre_lt1 cre_lt1_
#define cre_lhe cre_lhe_
#define cre_liw cre_liw_
#define cre_lmt cre_lmt_
#define cre_lmr cre_lmr_
#define cre_lf1 cre_lf1_
#define cre_l1s cre_l1s_
#define cre_lps cre_lps_
#define cre_lrp cre_lrp_
#define cre_lba cre_lba_
#define cre_lc3 cre_lc3_
#define cre_lwc cre_lwc_
#define cre_lil cre_lil_
#define cre_lcb cre_lcb_
#define cre_lct cre_lct_
#define cre_lcv cre_lcv_
#define cre_lmd cre_lmd_
#define cre_lco cre_lco_
#define cre_lty cre_lty_
#define cre_lcp cre_lcp_
#define cre_lvd cre_lvd_
#define cre_lve cre_lve_
#define cre_lma cre_lma_
#define cre_lmf cre_lmf_
#define cre_z02 cre_z02_
#define cre_ram cre_ram_
#define cre_rsf cre_rsf_
#define cre_rsm cre_rsm_
#define cre_vrc cre_vrc_
#define cre_lm1 cre_lm1_
#define cre_lm2 cre_lm2_
#define cre_vrd cre_vrd_
#define cre_v_radar cre_v_radar_
#define cre_lga cre_lga_
#define cre_lgf cre_lgf_
#define cre_ln3 cre_ln3_
#define log_diag log_diag_
#define itoa itoa_
#define fill_empty_grids fill_empty_grids_
#define cdf_get_coord cdf_get_coord_
#define cdf_get_levels cdf_get_levels_
#define cdf_get_index cdf_get_index_
#define cdf_dim_size cdf_dim_size_
#define cdf_check_laps_inv cdf_check_laps_inv_
#define cdf_write_grid cdf_write_grid_
#define cdf_update_laps_inv cdf_update_laps_inv_
#define make_c_fname make_c_fname_
#define get_cdf_var get_cdf_var_
#define open_cdf open_cdf_
#define cdf_retrieve_hdr cdf_retrieve_hdr_
#define write_hdr_cdf write_hdr_cdf_
#define cdf_update_laps cdf_update_laps_
#define cdf_retrieve_laps_grid cdf_retrieve_laps_grid_
#define write_cdf_file write_cdf_file_
#define read_cdf_file read_cdf_file_
#define cdf_inquire_var cdf_inquire_var_
#define read_cdf_header read_cdf_header_
#define nstrncpy nstrncpy_
#define fstrncpy fstrncpy_
#endif  
#ifdef FORTRANCAPS
#define cre_lw3 CRE_LW3
#define cre_lh1 CRE_LH1
#define cre_lh2 CRE_LH2
#define cre_lh3 CRE_LH3
#define cre_lh4 CRE_LH4
#define cre_lq3 CRE_LQ3
#define cre_lsx CRE_LSX
#define cre_lwm CRE_LWM
#define cre_lt1 CRE_LT1
#define cre_lhe CRE_LHE
#define cre_liw CRE_LIW
#define cre_lmt CRE_LMT
#define cre_lmr CRE_LMR
#define cre_lf1 CRE_LF1
#define cre_l1s CRE_L1S
#define cre_lps CRE_LPS
#define cre_lrp CRE_LRP
#define cre_lba CRE_LBA
#define cre_lc3 CRE_LC3
#define cre_lwc CRE_LWC
#define cre_lil CRE_LI1
#define cre_lcb CRE_LCB
#define cre_lct CRE_LCT
#define cre_lcv CRE_LCV
#define cre_lmd CRE_LMD
#define cre_lco CRE_LCO
#define cre_lty CRE_LTY
#define cre_lcp CRE_LCP
#define cre_lvd CRE_LVD
#define cre_lve CRE_LVE
#define cre_lma CRE_LMA
#define cre_lmf CRE_LMF
#define cre_z02 CRE_Z02
#define cre_ram CRE_RAM
#define cre_rsf CRE_RSF
#define cre_rsm CRE_RSM
#define cre_vrc CRE_VRC
#define cre_lm1 CRE_LM1
#define cre_lm2 CRE_LM2
#define cre_vrd CRE_VRD
#define cre_v_radar CRE_V_RADAR
#define cre_lga CRE_LGA
#define cre_lgf CRE_LGF
#define cre_ln3 CRE_LN3
#define log_diag LOG_DIAG
#define itoa ITOA
#define fill_empty_grids FILL_EMPTY_GRIDS
#define cdf_get_coord CDF_GET_COORD
#define cdf_get_levels CDF_GET_LEVELS
#define cdf_get_index CDF_GET_INDEX
#define cdf_dim_size CDF_DIM_SIZE
#define cdf_check_laps_inv CDF_CHECK_LAPS_INV
#define cdf_write_grid CDF_WRITE_GRID
#define cdf_update_laps_inv CDF_UPDATE_LAPS_INV
#define make_c_fname MAKE_C_FNAME
#define get_cdf_var GET_CDF_VAR
#define open_cdf OPEN_CDF
#define cdf_retrieve_hdr CDF_RETRIEVE_HDR
#define write_hdr_cdf WRITE_HEADER_CDF
#define cdf_update_laps CDF_UPDATE_LAPS
#define cdf_retrieve_laps_grid CDF_RETRIEVE_LAPS_GRID
#define write_cdf_file WRITE_CDF_FILE
#define read_cdf_file READ_CDF_FILE
#define cdf_inquire_var CDF_INQUIRE_VAR
#define read_cdf_header READ_CDF_HEADER
#define nstrncpy NSTRNCPY
#define fstrncpy FSTRNCPY
#endif  
#ifdef FORTRANDOUBLEUNDERSCORE
#define cre_lw3 cre_lw3_
#define cre_lh1 cre_lh1__
#define cre_lh2 cre_lh2__
#define cre_lh3 cre_lh3__
#define cre_lh4 cre_lh4__
#define cre_lq3 cre_lq3__
#define cre_lsx cre_lsx__
#define cre_lwm cre_lwm__
#define cre_lt1 cre_lt1__
#define cre_lhe cre_lhe__
#define cre_liw cre_liw__
#define cre_lmt cre_lmt__
#define cre_lmr cre_lmr__
#define cre_lf1 cre_lf1__
#define cre_l1s cre_l1s__
#define cre_lps cre_lps__
#define cre_lrp cre_lrp__
#define cre_lba cre_lba__
#define cre_lc3 cre_lc3__
#define cre_lwc cre_lwc__
#define cre_lil cre_lil__
#define cre_lcb cre_lcb__
#define cre_lct cre_lct__
#define cre_lcv cre_lcv__
#define cre_lmd cre_lmd__
#define cre_lco cre_lco__
#define cre_lty cre_lty__
#define cre_lcp cre_lcp__
#define cre_lvd cre_lvd__
#define cre_lve cre_lve__
#define cre_lma cre_lma__
#define cre_lmf cre_lmf__
#define cre_z02 cre_z02__
#define cre_ram cre_ram__
#define cre_rsf cre_rsf__
#define cre_rsm cre_rsm__
#define cre_vrc cre_vrc__
#define cre_lm1 cre_lm1__
#define cre_lm2 cre_lm2__
#define cre_vrd cre_vrd__
#define cre_v_radar cre_v_radar__
#define cre_lga cre_lga__
#define cre_lgf cre_lgf__
#define cre_ln3 cre_ln3__
#define log_diag log_diag__
#define itoa itoa__
#define fill_empty_grids fill_empty_grids__
#define cdf_get_coord cdf_get_coord__
#define cdf_get_levels cdf_get_levels__
#define cdf_get_index cdf_get_index__
#define cdf_dim_size cdf_dim_size__
#define cdf_check_laps_inv cdf_check_laps_inv__
#define cdf_write_grid cdf_write_grid__
#define cdf_update_laps_inv cdf_update_laps_inv__
#define make_c_fname make_c_fname__
#define get_cdf_var get_cdf_var__
#define open_cdf open_cdf__
#define cdf_retrieve_hdr cdf_retrieve_hdr__
#define write_hdr_cdf write_hdr_cdf__
#define cdf_update_laps cdf_update_laps__
#define cdf_retrieve_laps_grid cdf_retrieve_laps_grid__
#define write_cdf_file write_cdf_file__
#define read_cdf_file read_cdf_file__
#define cdf_inquire_var cdf_inquire_var__
#define read_cdf_header read_cdf_header__
#define nstrncpy nstrncpy__
#define fstrncpy fstrncpy__
#endif  
/*************************************************************************
*	List of functions in this file:
*		cre_lw3,cre_lh1,cre_lh2,cre_lh3,cre_lh4,cre_lq3,
*               cre_lsx, cre_lwm, cre_lt1, cre_lhe, cre_liw, cre_lmt,
*               cre_lmr, cre_lf1, cre_l1s, cre_lps, cre_lrp, cre_lba,
*               cre_lc3, cre_lwc, cre_lil, cre_lcb, cre_lct, cre_lcv, 
*               cre_lmd, cre_lco, cre_lty, cre_lcp, cre_z02
*		log_diag  (does nothing, set to printf statment for debug)
*		itoa	  (converts integer to ascii)
*		fill_empty_grids
*		cdf_get_coord
*		cdf_get_index
*		cdf_dim_size
*		cdf_chk_laps_inv
*		cdf_write_grid
*		cdf_update_laps_inv
*		make_c_fname
*		get_cdf_var
*		open_cdf
*		cdf_retrieve_header
*		write_hdr_cdf
*		cdf_update_laps
*		cdf_retrieve_laps_grid
*		write_cdf_file
*		read_cdf_file
*		cdf_inquire_var 
* 		read_cdf_header 
*
*************************************************************************/
/* file to create netCDF format file with extension LW3 */
#ifdef __STDC__
int cre_lw3(char *fname)
#else
int cre_lw3(fname) 		/* create fname */
char *fname;
#endif
{   
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  u_id, v_id, w_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        u_comment_id, v_comment_id, w_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, u_fcinv_id, v_fcinv_id, 
        w_fcinv_id, origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  u_valid_range[2];
   float  v_valid_range[2];
   float  w_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 65L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   u_id = ncvardef (cdfid, "u", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   v_id = ncvardef (cdfid, "v", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   w_id = ncvardef (cdfid, "w", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   u_comment_id = ncvardef (cdfid, "u_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   v_comment_id = ncvardef (cdfid, "v_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   w_comment_id = ncvardef (cdfid, "w_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   u_fcinv_id = ncvardef (cdfid, "u_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   v_fcinv_id = ncvardef (cdfid, "v_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   w_fcinv_id = ncvardef (cdfid, "w_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, u_id, "long_name", NC_CHAR, 13, (void *)"eastward wind");
   ncattput (cdfid, u_id, "units", NC_CHAR, 13, (void *)"meters/second");
   u_valid_range[0] = -200;
   u_valid_range[1] = 200;
   ncattput (cdfid, u_id, "valid_range", NC_FLOAT, 2, (void *) u_valid_range);
/*   float_val = 1.0e+37;
   ncattput (cdfid, u_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val); */
   ncattput (cdfid, u_id, "LAPS_var", NC_CHAR, 2, (void *)"U3");
   ncattput (cdfid, u_id, "lvl_coord_3D", NC_CHAR, 4, (void *)"HPA ");
   ncattput (cdfid, u_id, "lvl_coord_sfc", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, u_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, v_id, "long_name", NC_CHAR, 14, (void *)"northward wind");
   ncattput (cdfid, v_id, "units", NC_CHAR, 13, (void *)"meters/second");
   v_valid_range[0] = -200;
   v_valid_range[1] = 200;
   ncattput (cdfid, v_id, "valid_range", NC_FLOAT, 2, (void *) v_valid_range);
/*   float_val = 1.0e+37;
   ncattput (cdfid, v_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val); */
   ncattput (cdfid, v_id, "LAPS_var", NC_CHAR, 2, (void *)"V3");
   ncattput (cdfid, v_id, "lvl_coord_3D", NC_CHAR, 4, (void *)"HPA ");
   ncattput (cdfid, v_id, "lvl_coord_sfc", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, v_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, w_id, "long_name", NC_CHAR, 10, (void *)"wind omega");
   ncattput (cdfid, w_id, "units", NC_CHAR, 14, (void *)"pascals/second");
   w_valid_range[0] = -20000;
   w_valid_range[1] = 20000;
   ncattput (cdfid, w_id, "valid_range", NC_FLOAT, 2, (void *) w_valid_range);
/*   float_val = 1.0e+37;
   ncattput (cdfid, w_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val); */
   ncattput (cdfid, w_id, "LAPS_var", NC_CHAR, 2, (void *)"OM");
   ncattput (cdfid, w_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (cdfid, w_id, "LAPS_units", NC_CHAR, 4, (void *)"PA/S");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, u_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, v_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, w_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 
        550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }
   
   {			/* store num_variables */
    static long num_variables = {21};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {63};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {63};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}

/*************************************************************************/
/* file to create netCDF format file with extension LH1 */
#ifdef __STDC__
int cre_lh1(char *fname) 		/* create fname */
#else
int cre_lh1(fname)
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  pw_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, pw_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, pw_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  pw_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 1L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pw_id = ncvardef (cdfid, "pw", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pw_comment_id = ncvardef (cdfid, "pw_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pw_fcinv_id = ncvardef (cdfid, "pw_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, pw_id, "long_name", NC_CHAR, 18, (void *)"precipitable water");
   ncattput (cdfid, pw_id, "units", NC_CHAR, 6, (void *)"meters");
   pw_valid_range[0] = 0;
   pw_valid_range[1] = 0.1;
   ncattput (cdfid, pw_id, "valid_range", NC_FLOAT, 2, (void *) pw_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, pw_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, pw_id, "LAPS_var", NC_CHAR, 2, (void *)"PW");
   ncattput (cdfid, pw_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'M';
   ncattput (cdfid, pw_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, pw_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {1};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {1};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}

/*************************************************************************/
/* file to create netCDF format file with extension LH2 */
#ifdef __STDC__
int cre_lh2(char *fname) 		/* create fname */
#else
int cre_lh2(fname)
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  lpw_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, lpw_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, lpw_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  lpw_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 3L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 3L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lpw_id = ncvardef (cdfid, "lpw", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lpw_comment_id = ncvardef (cdfid, "lpw_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lpw_fcinv_id = ncvardef (cdfid, "lpw_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, lpw_id, "long_name", NC_CHAR, 24, (void *)"layer precipitable water");
   ncattput (cdfid, lpw_id, "units", NC_CHAR, 6, (void *)"meters");
   lpw_valid_range[0] = 0;
   lpw_valid_range[1] = 0.1;
   ncattput (cdfid, lpw_id, "valid_range", NC_FLOAT, 2, (void *) lpw_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, lpw_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, lpw_id, "LAPS_var", NC_CHAR, 2, (void *)"PW");
   ncattput (cdfid, lpw_id, "lvl_coord", NC_CHAR, 3, (void *)"N/A");
   char_val = 'M';
   ncattput (cdfid, lpw_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, lpw_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {3};
    static short level[] = {1, 2, 3};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {3};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {3};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}

/*************************************************************************/
/* file to create netCDF format file with extension LH3 */
#ifdef __STDC__
int cre_lh3(char *fname) 		/* create fname */
#else
int cre_lh3(fname)
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  rh_id, rhl_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        rh_comment_id, rhl_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, rh_fcinv_id, rhl_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  rh_valid_range[2];
   float  rhl_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 42L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rh_id = ncvardef (cdfid, "rh", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rhl_id = ncvardef (cdfid, "rhl", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rh_comment_id = ncvardef (cdfid, "rh_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rhl_comment_id = ncvardef (cdfid, "rhl_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rh_fcinv_id = ncvardef (cdfid, "rh_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rhl_fcinv_id = ncvardef (cdfid, "rhl_fcinv", NC_SHORT, 2, dims);
   
   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, rh_id, "long_name", NC_CHAR, 17, (void *)"relative humidity");
   ncattput (cdfid, rh_id, "units", NC_CHAR, 7, (void *)"percent");
   rh_valid_range[0] = 0;
   rh_valid_range[1] = 100;
   ncattput (cdfid, rh_id, "valid_range", NC_FLOAT, 2, (void *) rh_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, rh_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, rh_id, "LAPS_var", NC_CHAR, 3, (void *)"RH3");
   ncattput (cdfid, rh_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, rh_id, "LAPS_units", NC_CHAR, 7, (void *)"PERCENT");
   ncattput (cdfid, rhl_id, "long_name", NC_CHAR, 17, (void *)"relative humidity from liquid");
   ncattput (cdfid, rhl_id, "units", NC_CHAR, 7, (void *)"percent");
   rh_valid_range[0] = 0;
   rh_valid_range[1] = 100;
   ncattput (cdfid, rhl_id, "valid_range", NC_FLOAT, 2, (void *) rh_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, rhl_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, rhl_id, "LAPS_var", NC_CHAR, 3, (void *)"RHL");
   ncattput (cdfid, rhl_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, rhl_id, "LAPS_units", NC_CHAR, 7, (void *)"PERCENT");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, rh_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   ncattput (cdfid, rhl_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {18};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {42};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {42};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}

/*************************************************************************/
/* file to create netCDF format file with extension LH4 */
#ifdef __STDC__
int cre_lh4(char * fname) 		/* create fname */
#else
int cre_lh4(fname)
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  tpw_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, tpw_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, tpw_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;


   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  tpw_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 1L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   tpw_id = ncvardef (cdfid, "tpw", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   tpw_comment_id = ncvardef (cdfid, "tpw_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   tpw_fcinv_id = ncvardef (cdfid, "tpw_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, tpw_id, "long_name", NC_CHAR, 35, (void *)"integrated total precipitable water");
   ncattput (cdfid, tpw_id, "units", NC_CHAR, 6, (void *)"meters");
   tpw_valid_range[0] = 0;
   tpw_valid_range[1] = 0.100;
   ncattput (cdfid, tpw_id, "valid_range", NC_FLOAT, 2, (void *) tpw_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, tpw_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, tpw_id, "LAPS_var", NC_CHAR, 3, (void *)"TPW");
   ncattput (cdfid, tpw_id, "lvl_coord", NC_CHAR, 4, (void *)"    ");
   ncattput (cdfid, tpw_id, "LAPS_units", NC_CHAR, 1, (void *)"M");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, tpw_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {1};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {1};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}

/*************************************************************************/
/* file to create netCDF format file with extension LQ3 */
#ifdef __STDC__
int cre_lq3(char *fname) 		/* create fname */
#else
int cre_lq3(fname)
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  sh_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, sh_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, sh_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;


   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  sh_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 21L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sh_id = ncvardef (cdfid, "sh", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sh_comment_id = ncvardef (cdfid, "sh_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sh_fcinv_id = ncvardef (cdfid, "sh_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, sh_id, "long_name", NC_CHAR, 17, (void *)"specific humidity");
   ncattput (cdfid, sh_id, "units", NC_CHAR, 5, (void *)"kg/kg");
   sh_valid_range[0] = 0;
   sh_valid_range[1] = 0.1;
   ncattput (cdfid, sh_id, "valid_range", NC_FLOAT, 2, (void *) sh_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, sh_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, sh_id, "LAPS_var", NC_CHAR, 2, (void *)"SH");
   ncattput (cdfid, sh_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, sh_id, "LAPS_units", NC_CHAR, 10, (void *)"          ");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, sh_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {21};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {21};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}

/*************************************************************************/
/* file to create netCDF format file with extension LSX */
#ifdef __STDC__
int cre_lsx(char *fname) 		/* create fname */
#else
int cre_lsx(fname)
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, 
        domain_len_dim, asc_len_dim;

   /* variable ids */
   int  su_id, sv_id, fp_id, st_id, std_id, vv_id, srh_id, ccg_id, mp_id, 
        ta_id, pot_id, ept_id, sp_id, vor_id, mr_id, mc_id, d_id, pta_id, ma_id, 
        li_id, spd_id, cssi_id, pbe_id, nbe_id, vis_id, fwx_id, hi_id, lvl_id, 
        imax_id, jmax_id, kmax_id, kdim_id, su_comment_id, sv_comment_id, 
        fp_comment_id, st_comment_id, std_comment_id, vv_comment_id, srh_comment_id, 
        ccg_comment_id, mp_comment_id, ta_comment_id, pot_comment_id, 
        ept_comment_id, sp_comment_id, vor_comment_id, mr_comment_id, 
        mc_comment_id, d_comment_id, pta_comment_id, ma_comment_id, 
        li_comment_id, spd_comment_id, cssi_comment_id, pbe_comment_id, 
        nbe_comment_id, vis_comment_id, fwx_comment_id, hi_comment_id,
        laps_domain_file_id, asctime_id, fctimes_id, level_id, su_fcinv_id, 
        sv_fcinv_id, fp_fcinv_id, st_fcinv_id, std_fcinv_id, vv_fcinv_id, 
        srh_fcinv_id, ccg_fcinv_id, mp_fcinv_id, ta_fcinv_id, pot_fcinv_id, 
        ept_fcinv_id, sp_fcinv_id, vor_fcinv_id, mr_fcinv_id, mc_fcinv_id, 
        d_fcinv_id, pta_fcinv_id, ma_fcinv_id, li_fcinv_id, spd_fcinv_id, 
        cssi_fcinv_id, pbe_fcinv_id, nbe_fcinv_id, vis_fcinv_id, fwx_fcinv_id, 
        hi_fcinv_id,origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short   short_val;
   float   float_val;

   /* attribute vectors */
   float  su_valid_range[2];
   float  sv_valid_range[2];
   float  fp_valid_range[2];
   float  st_valid_range[2];
   float  std_valid_range[2];
   float  vv_valid_range[2];
   float  srh_valid_range[2];
   float  ccg_valid_range[2];
   float  mp_valid_range[2];
   float  ta_valid_range[2];
   float  pot_valid_range[2];
   float  ept_valid_range[2];
   float  sp_valid_range[2];
   float  vor_valid_range[2];
   float  mr_valid_range[2];
   float  mc_valid_range[2];
   float  d_valid_range[2];
   float  pta_valid_range[2];
   float  ma_valid_range[2];
   float  li_valid_range[2];
   float  spd_valid_range[2];
   float  cssi_valid_range[2];
   float  pbe_valid_range[2];
   float  nbe_valid_range[2];
   float  vis_valid_range[2];
   float  fwx_valid_range[2];
   float  hi_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 27L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   su_id = ncvardef (cdfid, "su", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sv_id = ncvardef (cdfid, "sv", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   fp_id = ncvardef (cdfid, "fp", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   st_id = ncvardef (cdfid, "st", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   std_id = ncvardef (cdfid, "std", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vv_id = ncvardef (cdfid, "vv", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   srh_id = ncvardef (cdfid, "srh", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ccg_id = ncvardef (cdfid, "ccg", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mp_id = ncvardef (cdfid, "mp", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ta_id = ncvardef (cdfid, "ta", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pot_id = ncvardef (cdfid, "pot", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ept_id = ncvardef (cdfid, "ept", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sp_id = ncvardef (cdfid, "sp", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vor_id = ncvardef (cdfid, "vor", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mr_id = ncvardef (cdfid, "mr", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mc_id = ncvardef (cdfid, "mc", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d_id = ncvardef (cdfid, "d", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pta_id = ncvardef (cdfid, "pta", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ma_id = ncvardef (cdfid, "ma", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   li_id = ncvardef (cdfid, "li", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   spd_id = ncvardef (cdfid, "spd", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   cssi_id = ncvardef (cdfid, "cssi", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pbe_id = ncvardef (cdfid, "pbe", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   nbe_id = ncvardef (cdfid, "nbe", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vis_id = ncvardef (cdfid, "vis", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   fwx_id = ncvardef (cdfid, "fwx", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   hi_id = ncvardef (cdfid, "hi", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   su_comment_id = ncvardef (cdfid, "su_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sv_comment_id = ncvardef (cdfid, "sv_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   fp_comment_id = ncvardef (cdfid, "fp_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   st_comment_id = ncvardef (cdfid, "st_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   std_comment_id = ncvardef (cdfid, "std_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vv_comment_id = ncvardef (cdfid, "vv_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   srh_comment_id = ncvardef (cdfid, "srh_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ccg_comment_id = ncvardef (cdfid, "ccg_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mp_comment_id = ncvardef (cdfid, "mp_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ta_comment_id = ncvardef (cdfid, "ta_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pot_comment_id = ncvardef (cdfid, "pot_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ept_comment_id = ncvardef (cdfid, "ept_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sp_comment_id = ncvardef (cdfid, "sp_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vor_comment_id = ncvardef (cdfid, "vor_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mr_comment_id = ncvardef (cdfid, "mr_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mc_comment_id = ncvardef (cdfid, "mc_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d_comment_id = ncvardef (cdfid, "d_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pta_comment_id = ncvardef (cdfid, "pta_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ma_comment_id = ncvardef (cdfid, "ma_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   li_comment_id = ncvardef (cdfid, "li_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   spd_comment_id = ncvardef (cdfid, "spd_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   cssi_comment_id = ncvardef (cdfid, "cssi_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pbe_comment_id = ncvardef (cdfid, "pbe_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   nbe_comment_id = ncvardef (cdfid, "nbe_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vis_comment_id = ncvardef (cdfid, "vis_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   fwx_comment_id = ncvardef (cdfid, "fwx_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   hi_comment_id = ncvardef (cdfid, "hi_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   su_fcinv_id = ncvardef (cdfid, "su_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sv_fcinv_id = ncvardef (cdfid, "sv_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   fp_fcinv_id = ncvardef (cdfid, "fp_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   st_fcinv_id = ncvardef (cdfid, "st_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   std_fcinv_id = ncvardef (cdfid, "std_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vv_fcinv_id = ncvardef (cdfid, "vv_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   srh_fcinv_id = ncvardef (cdfid, "srh_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ccg_fcinv_id = ncvardef (cdfid, "ccg_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mp_fcinv_id = ncvardef (cdfid, "mp_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ta_fcinv_id = ncvardef (cdfid, "ta_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pot_fcinv_id = ncvardef (cdfid, "pot_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ept_fcinv_id = ncvardef (cdfid, "ept_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sp_fcinv_id = ncvardef (cdfid, "sp_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vor_fcinv_id = ncvardef (cdfid, "vor_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mr_fcinv_id = ncvardef (cdfid, "mr_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mc_fcinv_id = ncvardef (cdfid, "mc_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d_fcinv_id = ncvardef (cdfid, "d_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pta_fcinv_id = ncvardef (cdfid, "pta_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ma_fcinv_id = ncvardef (cdfid, "ma_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   li_fcinv_id = ncvardef (cdfid, "li_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   spd_fcinv_id = ncvardef (cdfid, "spd_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   cssi_fcinv_id = ncvardef (cdfid, "cssi_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pbe_fcinv_id = ncvardef (cdfid, "pbe_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   nbe_fcinv_id = ncvardef (cdfid, "nbe_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vis_fcinv_id = ncvardef (cdfid, "vis_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   fwx_fcinv_id = ncvardef (cdfid, "fwx_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   hi_fcinv_id = ncvardef (cdfid, "hi_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, su_id, "long_name", NC_CHAR, 34, (void *)"grid relative u component sfc wind");
   ncattput (cdfid, su_id, "units", NC_CHAR, 13, (void *)"meters/second");
   su_valid_range[0] = -200;
   su_valid_range[1] = 200;
   ncattput (cdfid, su_id, "valid_range", NC_FLOAT, 2, (void *) su_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, su_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, su_id, "LAPS_var", NC_CHAR, 1, (void *)"U");
   ncattput (cdfid, su_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, su_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, sv_id, "long_name", NC_CHAR, 34, (void *)"grid relative v component sfc wind");
   ncattput (cdfid, sv_id, "units", NC_CHAR, 13, (void *)"meters/second");
   sv_valid_range[0] = -200;
   sv_valid_range[1] = 200;
   ncattput (cdfid, sv_id, "valid_range", NC_FLOAT, 2, (void *) sv_valid_range);
   ncattput (cdfid, sv_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, sv_id, "LAPS_var", NC_CHAR, 1, (void *)"V");
   ncattput (cdfid, sv_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, sv_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, fp_id, "long_name", NC_CHAR, 23, (void *)"1500 m reduced pressure");
   ncattput (cdfid, fp_id, "units", NC_CHAR, 7, (void *)"pascals");
   fp_valid_range[0] = 0;
   fp_valid_range[1] = 1000000;
   ncattput (cdfid, fp_id, "valid_range", NC_FLOAT, 2, (void *) fp_valid_range);
   ncattput (cdfid, fp_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, fp_id, "LAPS_var", NC_CHAR, 1, (void *)"P");
   ncattput (cdfid, fp_id, "lvl_coord", NC_CHAR, 4, (void *)"AGL ");
   ncattput (cdfid, fp_id, "LAPS_units", NC_CHAR, 2, (void *)"PA");
   ncattput (cdfid, st_id, "long_name", NC_CHAR, 19, (void *)"surface temperature");
   ncattput (cdfid, st_id, "units", NC_CHAR, 14, (void *)"degrees kelvin");
   st_valid_range[0] = 210;
   st_valid_range[1] = 366;
   ncattput (cdfid, st_id, "valid_range", NC_FLOAT, 2, (void *) st_valid_range);
   ncattput (cdfid, st_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, st_id, "LAPS_var", NC_CHAR, 1, (void *)"T");
   ncattput (cdfid, st_id, "lvl_coord", NC_CHAR, 4, (void *)"AGL ");
   char_val = 'K';
   ncattput (cdfid, st_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, std_id, "long_name", NC_CHAR, 28, (void *)"surface dewpoint temperature");
   ncattput (cdfid, std_id, "units", NC_CHAR, 14, (void *)"degrees kelvin");
   std_valid_range[0] = 210;
   std_valid_range[1] = 366;
   ncattput (cdfid, std_id, "valid_range", NC_FLOAT, 2, (void *) std_valid_range);
   ncattput (cdfid, std_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, std_id, "LAPS_var", NC_CHAR, 2, (void *)"TD");
   ncattput (cdfid, std_id, "lvl_coord", NC_CHAR, 4, (void *)"AGL ");
   char_val = 'K';
   ncattput (cdfid, std_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, vv_id, "long_name", NC_CHAR, 17, (void *)"vertical velocity");
   ncattput (cdfid, vv_id, "units", NC_CHAR, 13, (void *)"meters/second");
   vv_valid_range[0] = 0;
   vv_valid_range[1] = 20000;
   ncattput (cdfid, vv_id, "valid_range", NC_FLOAT, 2, (void *) vv_valid_range);
   ncattput (cdfid, vv_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, vv_id, "LAPS_var", NC_CHAR, 2, (void *)"VV");
   ncattput (cdfid, vv_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, vv_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, srh_id, "long_name", NC_CHAR, 17, (void *)"relative humidity");
   ncattput (cdfid, srh_id, "units", NC_CHAR, 7, (void *)"percent");
   srh_valid_range[0] = 0;
   srh_valid_range[1] = 100;
   ncattput (cdfid, srh_id, "valid_range", NC_FLOAT, 2, (void *) srh_valid_range);
   ncattput (cdfid, srh_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, srh_id, "LAPS_var", NC_CHAR, 2, (void *)"RH");
   ncattput (cdfid, srh_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, srh_id, "LAPS_units", NC_CHAR, 7,(void *) "PERCENT");
   ncattput (cdfid, ccg_id, "long_name", NC_CHAR, 13, (void *)"cloud ceiling");
   ncattput (cdfid, ccg_id, "units", NC_CHAR, 6, (void *)"meters");
   ccg_valid_range[0] = -20000;
   ccg_valid_range[1] = 20000;
   ncattput (cdfid, ccg_id, "valid_range", NC_FLOAT, 2, (void *) ccg_valid_range);
   ncattput (cdfid, ccg_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, ccg_id, "LAPS_var", NC_CHAR, 3, (void *)"CCE");
   ncattput (cdfid, ccg_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'M';
   ncattput (cdfid, ccg_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, mp_id, "long_name", NC_CHAR, 12, (void *)"MSL pressure");
   ncattput (cdfid, mp_id, "units", NC_CHAR, 7, (void *)"pascals");
   mp_valid_range[0] = 0;
   mp_valid_range[1] = 100000;
   ncattput (cdfid, mp_id, "valid_range", NC_FLOAT, 2, (void *) mp_valid_range);
   ncattput (cdfid, mp_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, mp_id, "LAPS_var", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, mp_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, mp_id, "LAPS_units", NC_CHAR, 2, (void *)"PA");
   ncattput (cdfid, ta_id, "long_name", NC_CHAR, 21, (void *)"temperature advection");
   ncattput (cdfid, ta_id, "units", NC_CHAR, 21, (void *)"degrees Kelvin/second");
   ta_valid_range[0] = -20000;
   ta_valid_range[1] = 20000;
   ncattput (cdfid, ta_id, "valid_range", NC_FLOAT, 2, (void *) ta_valid_range);
   ncattput (cdfid, ta_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, ta_id, "LAPS_var", NC_CHAR, 3, (void *)"TAD");
   ncattput (cdfid, ta_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, ta_id, "LAPS_units", NC_CHAR, 3, (void *)"K/S");
   ncattput (cdfid, pot_id, "long_name", NC_CHAR, 21, (void *)"potential temperature");
   ncattput (cdfid, pot_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   pot_valid_range[0] = 210;
   pot_valid_range[1] = 366;
   ncattput (cdfid, pot_id, "valid_range", NC_FLOAT, 2, (void *) pot_valid_range);
   ncattput (cdfid, pot_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, pot_id, "LAPS_var", NC_CHAR, 2, (void *)"TH");
   ncattput (cdfid, pot_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (cdfid, pot_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, ept_id, "long_name", NC_CHAR, 32, (void *)"equivalent potential temperature");
   ncattput (cdfid, ept_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   ept_valid_range[0] = 210;
   ept_valid_range[1] = 366;
   ncattput (cdfid, ept_id, "valid_range", NC_FLOAT, 2, (void *) ept_valid_range);
   ncattput (cdfid, ept_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, ept_id, "LAPS_var", NC_CHAR, 3, (void *)"THE");
   ncattput (cdfid, ept_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (cdfid, ept_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, sp_id, "long_name", NC_CHAR, 16, (void *)"surface pressure");
   ncattput (cdfid, sp_id, "units", NC_CHAR, 7, (void *)"pascals");
   sp_valid_range[0] = 0;
   sp_valid_range[1] = 100000;
   ncattput (cdfid, sp_id, "valid_range", NC_FLOAT, 2, (void *) sp_valid_range);
   ncattput (cdfid, sp_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, sp_id, "LAPS_var", NC_CHAR, 2, (void *)"PS");
   ncattput (cdfid, sp_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, sp_id, "LAPS_units", NC_CHAR, 2, (void *)"PA");
   ncattput (cdfid, vor_id, "long_name", NC_CHAR, 9, (void *)"vorticity");
   ncattput (cdfid, vor_id, "units", NC_CHAR, 7, (void *)"/second");
   vor_valid_range[0] = -20000;
   vor_valid_range[1] = 20000;
   ncattput (cdfid, vor_id, "valid_range", NC_FLOAT, 2, (void *) vor_valid_range);
   ncattput (cdfid, vor_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, vor_id, "LAPS_var", NC_CHAR, 3, (void *)"VOR");
   ncattput (cdfid, vor_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, vor_id, "LAPS_units", NC_CHAR, 2, (void *)"/S");
   ncattput (cdfid, mr_id, "long_name", NC_CHAR, 12, (void *)"mixing ratio");
   ncattput (cdfid, mr_id, "units", NC_CHAR, 14, (void *)"grams/kikogram");
   mr_valid_range[0] = -20000;
   mr_valid_range[1] = 20000;
   ncattput (cdfid, mr_id, "valid_range", NC_FLOAT, 2, (void *) mr_valid_range);
   ncattput (cdfid, mr_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, mr_id, "LAPS_var", NC_CHAR, 2, (void *)"MR");
   ncattput (cdfid, mr_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, mr_id, "LAPS_units", NC_CHAR, 4, (void *)"G/KG");
   ncattput (cdfid, mc_id, "long_name", NC_CHAR, 20, (void *)"moisture convergence");
   ncattput (cdfid, mc_id, "units", NC_CHAR, 22, (void *)"grams/kilogram/seconds");
   mc_valid_range[0] = -20000;
   mc_valid_range[1] = 20000;
   ncattput (cdfid, mc_id, "valid_range", NC_FLOAT, 2, (void *) mc_valid_range);
   ncattput (cdfid, mc_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, mc_id, "LAPS_var", NC_CHAR, 3, (void *)"MRC");
   ncattput (cdfid, mc_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, mc_id, "LAPS_units", NC_CHAR, 6, (void *)"G/KG/S");
   ncattput (cdfid, d_id, "long_name", NC_CHAR, 10, (void *)"divergence");
   ncattput (cdfid, d_id, "units", NC_CHAR, 7, (void *)"/second");
   d_valid_range[0] = -20000;
   d_valid_range[1] = 20000;
   ncattput (cdfid, d_id, "valid_range", NC_FLOAT, 2, (void *) d_valid_range);
   ncattput (cdfid, d_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d_id, "LAPS_var", NC_CHAR, 3, (void *)"DIV");
   ncattput (cdfid, d_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, d_id, "LAPS_units", NC_CHAR, 2, (void *)"/S");
   ncattput (cdfid, pta_id, "long_name", NC_CHAR, 31, (void *)"potential temperature advection");
   ncattput (cdfid, pta_id, "units", NC_CHAR, 16, (void *)"kilograms/second");
   pta_valid_range[0] = -20000;
   pta_valid_range[1] = 20000;
   ncattput (cdfid, pta_id, "valid_range", NC_FLOAT, 2, (void *) pta_valid_range);
   ncattput (cdfid, pta_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, pta_id, "LAPS_var", NC_CHAR, 3, (void *)"THA");
   ncattput (cdfid, pta_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, pta_id, "LAPS_units", NC_CHAR, 3, (void *)"K/S");
   ncattput (cdfid, ma_id, "long_name", NC_CHAR, 18, (void *)"moisture advection");
   ncattput (cdfid, ma_id, "units", NC_CHAR, 21, (void *)"grams/kilogram/second");
   ma_valid_range[0] = -20000;
   ma_valid_range[1] = 20000;
   ncattput (cdfid, ma_id, "valid_range", NC_FLOAT, 2, (void *) ma_valid_range);
   ncattput (cdfid, ma_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, ma_id, "LAPS_var", NC_CHAR, 3, (void *)"MRA");
   ncattput (cdfid, ma_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, ma_id, "LAPS_units", NC_CHAR, 6, (void *)"G/KG/S");
   ncattput (cdfid, li_id, "long_name", NC_CHAR, 12, (void *)"lifted index");
   ncattput (cdfid, li_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   li_valid_range[0] = 210;
   li_valid_range[1] = 366;
   ncattput (cdfid, li_id, "valid_range", NC_FLOAT, 2, (void *) li_valid_range);
   ncattput (cdfid, li_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, li_id, "LAPS_var", NC_CHAR, 2, (void *)"LI");
   ncattput (cdfid, li_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (cdfid, li_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, spd_id, "long_name", NC_CHAR, 18, (void *)"surface wind speed");
   ncattput (cdfid, spd_id, "units", NC_CHAR, 13, (void *)"meters/second");
   spd_valid_range[0] = -20000;
   spd_valid_range[1] = 20000;
   ncattput (cdfid, spd_id, "valid_range", NC_FLOAT, 2, (void *) spd_valid_range);
   ncattput (cdfid, spd_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, spd_id, "LAPS_var", NC_CHAR, 3, (void *)"SPD");
   ncattput (cdfid, spd_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, spd_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, cssi_id, "long_name", NC_CHAR, 27, (void *)"colorado severe storm index");
   ncattput (cdfid, cssi_id, "units", NC_CHAR, 10, (void *)"          ");
   cssi_valid_range[0] = -20000;
   cssi_valid_range[1] = 20000;
   ncattput (cdfid, cssi_id, "valid_range", NC_FLOAT, 2, (void *) cssi_valid_range);
   ncattput (cdfid, cssi_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, cssi_id, "LAPS_var", NC_CHAR, 3, (void *)"CSS");
   ncattput (cdfid, cssi_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, cssi_id, "LAPS_units", NC_CHAR, 10, (void *)"          ");
   ncattput (cdfid, pbe_id, "long_name", NC_CHAR, 23, (void *)"positive buoyant energy");
   ncattput (cdfid, pbe_id, "units", NC_CHAR, 15, (void *)"joules/kilogram");
   pbe_valid_range[0] = -20000;
   pbe_valid_range[1] = 20000;
   ncattput (cdfid, pbe_id, "valid_range", NC_FLOAT, 2, (void *) pbe_valid_range);
   ncattput (cdfid, pbe_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, pbe_id, "LAPS_var", NC_CHAR, 3, (void *)"PBE");
   ncattput (cdfid, pbe_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, pbe_id, "LAPS_units", NC_CHAR, 4, (void *)"J/KG");
   ncattput (cdfid, nbe_id, "long_name", NC_CHAR, 23, (void *)"negative buoyant energy");
   ncattput (cdfid, nbe_id, "units", NC_CHAR, 15, (void *)"joules/kilogram");
   nbe_valid_range[0] = -20000;
   nbe_valid_range[1] = 20000;
   ncattput (cdfid, nbe_id, "valid_range", NC_FLOAT, 2, (void *) nbe_valid_range);
   ncattput (cdfid, nbe_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, nbe_id, "LAPS_var", NC_CHAR, 3, (void *)"NBE");
   ncattput (cdfid, nbe_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, nbe_id, "LAPS_units", NC_CHAR, 4, (void *)"J/KG");
   ncattput (cdfid, vis_id, "long_name", NC_CHAR, 10, (void *)"visibility");
   ncattput (cdfid, vis_id, "units", NC_CHAR, 6, (void *)"meters");
   vis_valid_range[0] = 0;
   vis_valid_range[1] = 50000;
   ncattput (cdfid, vis_id, "valid_range", NC_FLOAT, 2, (void *) vis_valid_range);
   ncattput (cdfid, vis_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, vis_id, "LAPS_var", NC_CHAR, 3, (void *)"VIS");
   ncattput (cdfid, vis_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'M';
   ncattput (cdfid, vis_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, fwx_id, "long_name", NC_CHAR, 11, (void *)"Fire danger");
   ncattput (cdfid, fwx_id, "units", NC_CHAR, 4, (void *)"none");
   vis_valid_range[0] = 0;
   vis_valid_range[1] = 20;
   ncattput (cdfid, fwx_id, "valid_range", NC_FLOAT, 2, (void *) vis_valid_range);
   ncattput (cdfid, fwx_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, fwx_id, "LAPS_var", NC_CHAR, 3, (void *)"FWX");
   ncattput (cdfid, fwx_id, "lvl_coord", NC_CHAR, 4, (void *)"none");
   ncattput (cdfid, fwx_id, "LAPS_units", NC_CHAR, 4,(void *) "none");
   ncattput (cdfid, hi_id, "long_name", NC_CHAR, 10, (void *)"Heat index");
   ncattput (cdfid, hi_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   hi_valid_range[0] = 210;
   hi_valid_range[1] = 366;
   ncattput (cdfid, hi_id, "valid_range", NC_FLOAT, 2, (void *) hi_valid_range);
   ncattput (cdfid, hi_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, hi_id, "LAPS_var", NC_CHAR, 2, (void *)"HI");
   ncattput (cdfid, hi_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, hi_id, "LAPS_units", NC_CHAR, 1, (void *)"K");

   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, su_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, sv_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, fp_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, st_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, std_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, vv_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, srh_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, ccg_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, mp_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, ta_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, pot_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, ept_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, sp_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, vor_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, mr_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, mc_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, d_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, pta_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, ma_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, li_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, spd_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, cssi_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, pbe_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, nbe_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, vis_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, fwx_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, hi_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {93};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {27};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {27};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}

/*************************************************************************/
/* file to create netCDF format file with extension LHE */
#ifdef __STDC__
int cre_lhe(char *fname)                /* create fname */
#else
int cre_lhe(fname)
char *fname;
#endif
{

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  hel_id, mu_id, mv_id, lvl_id, imax_id, jmax_id, kmax_id, 
        kdim_id, hel_comment_id, mu_comment_id, mv_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, 
        hel_fcinv_id, mu_fcinv_id, mv_fcinv_id, origin_id, model_id, 
        version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  hel_valid_range[2];
   float  mu_valid_range[2];
   float  mv_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 3L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   hel_id = ncvardef (ncid, "hel", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mu_id = ncvardef (ncid, "mu", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mv_id = ncvardef (ncid, "mv", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   hel_comment_id = ncvardef (ncid, "hel_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mu_comment_id = ncvardef (ncid, "mu_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mv_comment_id = ncvardef (ncid, "mv_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   hel_fcinv_id = ncvardef (ncid, "hel_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mu_fcinv_id = ncvardef (ncid, "mu_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mv_fcinv_id = ncvardef (ncid, "mv_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, hel_id, "long_name", NC_CHAR, 8, (void *)"helicity");
   ncattput (ncid, hel_id, "units", NC_CHAR, 16, (void *)"meters/second**2");
   hel_valid_range[0] = 0;
   hel_valid_range[1] = 0.1;
   ncattput (ncid, hel_id, "valid_range", NC_FLOAT, 2, (void *) hel_valid_range);
   float_val = 1e+37;
   ncattput (ncid, hel_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, hel_id, "LAPS_var", NC_CHAR, 3, (void *)"LHE");
   ncattput (ncid, hel_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, hel_id, "LAPS_units", NC_CHAR, 6, (void *)"M/S**2");
   ncattput (ncid, mu_id, "long_name", NC_CHAR, 31, (void *)"SFC-300mb mean wind-u component");
   ncattput (ncid, mu_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mu_valid_range[0] = 0;
   mu_valid_range[1] = 0.1;
   ncattput (ncid, mu_id, "valid_range", NC_FLOAT, 2, (void *) mu_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, mu_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mu_id, "LAPS_var", NC_CHAR, 2, (void *)"MU");
   ncattput (ncid, mu_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, mu_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, mv_id, "long_name", NC_CHAR, 31, (void *)"SFC-300mb mean wind-v component");
   ncattput (ncid, mv_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mv_valid_range[0] = 0;
   mv_valid_range[1] = 0.1;
   ncattput (ncid, mv_id, "valid_range", NC_FLOAT, 2, (void *) mv_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, mv_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mv_id, "LAPS_var", NC_CHAR, 2, (void *)"MV");
   ncattput (ncid, mv_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, mv_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, hel_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mu_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mv_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {21};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {132};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {132};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {3};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {3};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return (ncid);
}
/*************************************************************************/
/* file to create netCDF format file with extension LT1 */
#ifdef __STDC__
int cre_lt1(char *fname) 		/* create fname */
#else
int cre_lt1(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim,
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, 
        domain_len_dim, asc_len_dim;

   /* variable ids */
   int  z_id, t_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, t_comment_id, 
        z_comment_id, laps_domain_file_id, asctime_id, fctimes_id, level_id, 
        t_fcinv_id, z_fcinv_id, origin_id, model_id, version_id, 
        num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  z_valid_range[2];
   float  t_valid_range[2];

   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1);
   namelen_dim = ncdimdef(cdfid, "namelen", 132);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 42);
   var_len_dim = ncdimdef(cdfid, "var_len", 4);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18);

   /* define variables */
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   z_id = ncvardef (cdfid, "z", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   t_id = ncvardef (cdfid, "t", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   t_comment_id = ncvardef (cdfid, "t_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   z_comment_id = ncvardef (cdfid, "z_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   t_fcinv_id = ncvardef (cdfid, "t_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   z_fcinv_id = ncvardef (cdfid, "z_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);


   /* assign attributes */
   ncattput (cdfid, z_id, "long_name", NC_CHAR, 6, (void *)"height");
   ncattput (cdfid, z_id, "units", NC_CHAR, 6, (void *)"meters");
   z_valid_range[0] = 0;
   z_valid_range[1] = 100000;
   ncattput (cdfid, z_id, "valid_range", NC_FLOAT, 2, (void *) z_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, z_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, z_id, "LAPS_var", NC_CHAR, 2, (void *)"HT");
   ncattput (cdfid, z_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, z_id, "LAPS_units", NC_CHAR, 6, (void *)"METERS");
   ncattput (cdfid, t_id, "long_name", NC_CHAR, 11, (void *)"temperature");
   ncattput (cdfid, t_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   t_valid_range[0] = 0;
   t_valid_range[1] = 100;
   ncattput (cdfid, t_id, "valid_range", NC_FLOAT, 2, (void *) t_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, t_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, t_id, "LAPS_var", NC_CHAR, 2,(void *) "T3");
   ncattput (cdfid, t_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, t_id, "LAPS_units", NC_CHAR, 1,(void *)"K");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, t_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, z_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
   static int fctimes_start[] = {0};
   static int fctimes_edges[] = {1};
   static short fctimes[] = {0};
   ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {18};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store level */
   static int level_start[] = {0};
   static int level_edges[] = {21};
   static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
   ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
   static int origin_start[] = {0};
   static int origin_edges[] = {49};
   static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
   ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
   static int model_start[] = {0};
   static int model_edges[] = {43};
   static char model[] = {"LAPS - Local Analysis and Prediction System"};
   ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {42};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {42};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}

/*************************************************************************/
/* file to create netCDF format file with extension LWM */
#ifdef __STDC__
int cre_lwm(char *fname)                /* create fname */
#else
int cre_lwm(fname)
char *fname;
#endif
{

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  u_id, v_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        u_comment_id, v_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, u_fcinv_id, v_fcinv_id, origin_id, 
        model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  u_valid_range[2];
   float  v_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 2L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   u_id = ncvardef (ncid, "u", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   v_id = ncvardef (ncid, "v", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   u_comment_id = ncvardef (ncid, "u_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   v_comment_id = ncvardef (ncid, "v_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   u_fcinv_id = ncvardef (ncid, "u_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   v_fcinv_id = ncvardef (ncid, "v_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, u_id, "long_name", NC_CHAR, 21, (void *)"surface eastward wind");
   ncattput (ncid, u_id, "units", NC_CHAR, 13, (void *)"meters/second");
   u_valid_range[0] = -200;
   u_valid_range[1] = 200;
   ncattput (ncid, u_id, "valid_range", NC_FLOAT, 2, (void *) u_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, u_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, u_id, "LAPS_var", NC_CHAR, 2, (void *)"SU");
   ncattput (ncid, u_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, u_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, v_id, "long_name", NC_CHAR, 22, (void *)"surface northward wind");
   ncattput (ncid, v_id, "units", NC_CHAR, 13, (void *)"meters/second");
   v_valid_range[0] = -200;
   v_valid_range[1] = 200;
   ncattput (ncid, v_id, "valid_range", NC_FLOAT, 2, (void *) v_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, v_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, v_id, "LAPS_var", NC_CHAR, 2, (void *)"SV");
   ncattput (ncid, v_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, v_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, u_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, v_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {18};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {132};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {132};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {2};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {2};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LIW */
#ifdef __STDC__
int cre_liw(char *fname) 		/* create fname */
#else
int cre_liw(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  liw_id, w_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        liw_comment_id, w_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, liw_fcinv_id, w_fcinv_id, origin_id, 
        model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  liw_valid_range[2];
   float  w_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 2L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   liw_id = ncvardef (cdfid, "liw", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   w_id = ncvardef (cdfid, "w", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   liw_comment_id = ncvardef (cdfid, "liw_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   w_comment_id = ncvardef (cdfid, "w_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   liw_fcinv_id = ncvardef (cdfid, "liw_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   w_fcinv_id = ncvardef (cdfid, "w_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, liw_id, "long_name", NC_CHAR, 21, (void *)"LAPS Li * 600mb Omega");
   ncattput (cdfid, liw_id, "units", NC_CHAR, 18, (void *)"kiloPascals/second");
   liw_valid_range[0] = 0;
   liw_valid_range[1] = 0.1;
   ncattput (cdfid, liw_id, "valid_range", NC_FLOAT, 2, (void *) liw_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, liw_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, liw_id, "LAPS_var", NC_CHAR, 3, (void *)"LIW");
   ncattput (cdfid, liw_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, liw_id, "LAPS_units", NC_CHAR, 6, (void *)"K-PA/S");
   ncattput (cdfid, w_id, "long_name", NC_CHAR, 16, (void *)"LAPS 600mb Omega");
   ncattput (cdfid, w_id, "units", NC_CHAR, 14, (void *)"Pascals/second");
   w_valid_range[0] = 0;
   w_valid_range[1] = 0.1;
   ncattput (cdfid, w_id, "valid_range", NC_FLOAT, 2, (void *) w_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, w_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, w_id, "LAPS_var", NC_CHAR, 1, (void *)"W");
   ncattput (cdfid, w_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, w_id, "LAPS_units", NC_CHAR, 4, (void *)"PA/S");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, liw_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, w_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {18};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {2};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {2};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LMT */
#ifdef __STDC__
int cre_lmt(char *fname) 		/* create fname */
#else
int cre_lmt(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  etop_id, llr_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        etop_comment_id, llr_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, etop_fcinv_id, llr_fcinv_id, origin_id, 
        model_id, version_id, num_variables_id;


   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  etop_valid_range[2];
   float  llr_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 2L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   etop_id = ncvardef (cdfid, "etop", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   llr_id = ncvardef (cdfid, "llr", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   etop_comment_id = ncvardef (cdfid, "etop_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   llr_comment_id = ncvardef (cdfid, "llr_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   etop_fcinv_id = ncvardef (cdfid, "etop_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   llr_fcinv_id = ncvardef (cdfid, "llr_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, etop_id, "long_name", NC_CHAR, 23, (void *)"maximum radar echo tops");
   ncattput (cdfid, etop_id, "units", NC_CHAR, 6, (void *)"meters");
   etop_valid_range[0] = 0;
   etop_valid_range[1] = 0.1;
   ncattput (cdfid, etop_id, "valid_range", NC_FLOAT, 2, (void *) etop_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, etop_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, etop_id, "LAPS_var", NC_CHAR, 3, (void *)"LMT");
   ncattput (cdfid, etop_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (cdfid, etop_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, llr_id, "long_name", NC_CHAR, 22, (void *)"low level reflectivity");
   ncattput (cdfid, llr_id, "units", NC_CHAR, 3, (void *)"dBZ");
   llr_valid_range[0] = 0;
   llr_valid_range[1] = 0.1;
   ncattput (cdfid, llr_id, "valid_range", NC_FLOAT, 2, (void *) llr_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, llr_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, llr_id, "LAPS_var", NC_CHAR, 3, (void *)"LLR");
   ncattput (cdfid, llr_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, llr_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, etop_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, llr_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {18};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {2};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {2};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LMR */
#ifdef __STDC__
int cre_lmr(char *fname) 		/* create fname */
#else
int cre_lmr(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  mxrf_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        mxrf_comment_id, laps_domain_file_id, asctime_id, fctimes_id, 
        level_id, mxrf_fcinv_id, origin_id, model_id, version_id,
        num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  mxrf_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 3L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 3L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mxrf_id = ncvardef (cdfid, "mxrf", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mxrf_comment_id = ncvardef (cdfid, "mxrf_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mxrf_fcinv_id = ncvardef (cdfid, "mxrf_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, mxrf_id, "long_name", NC_CHAR, 29, (void *)"column max radar reflectivity");
   ncattput (cdfid, mxrf_id, "units", NC_CHAR, 3, (void *)"dBZ");
   mxrf_valid_range[0] = 0;
   mxrf_valid_range[1] = 0.1;
   ncattput (cdfid, mxrf_id, "valid_range", NC_FLOAT, 2, (void *) mxrf_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, mxrf_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, mxrf_id, "LAPS_var", NC_CHAR, 1, (void *)"R");
   ncattput (cdfid, mxrf_id, "lvl_coord", NC_CHAR, 4, (void *)"    ");
   ncattput (cdfid, mxrf_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, fctimes_id, "units", NC_CHAR, 5, (void *)"hours");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, mxrf_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {3};
    static short fctimes[] = {0, 1, 2};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {3};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {3};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LF1 */
#ifdef __STDC__
int cre_lf1(char *fname) 		/* create fname */
#else
int cre_lf1(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  mxrh_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, mxrh_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, mxrh_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  mxrh_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 3L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 3L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mxrh_id = ncvardef (cdfid, "mxrh", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mxrh_comment_id = ncvardef (cdfid, "mxrh_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mxrh_fcinv_id = ncvardef (cdfid, "mxrh_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, mxrh_id, "long_name", NC_CHAR, 30, (void *)"max radar reflectivity history");
   ncattput (cdfid, mxrh_id, "units", NC_CHAR, 3, (void *)"dBZ");
   mxrh_valid_range[0] = 0;
   mxrh_valid_range[1] = 0.1;
   ncattput (cdfid, mxrh_id, "valid_range", NC_FLOAT, 2, (void *) mxrh_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, mxrh_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, mxrh_id, "LAPS_var", NC_CHAR, 1, (void *)"H");
   ncattput (cdfid, mxrh_id, "lvl_coord", NC_CHAR, 4, (void *)"    ");
   ncattput (cdfid, mxrh_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, fctimes_id, "units", NC_CHAR, 5, (void *)"hours");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, mxrh_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {3};
    static short fctimes[] = {0, 1, 2};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {3};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {3};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension L1S */
#ifdef __STDC__
int cre_l1s(char *fname) 		/* create fname */
#else
int cre_l1s(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  s1hr_id, stot_id, pc_id, pt_id, lvl_id, imax_id, jmax_id, kmax_id, 
        kdim_id, s1hr_comment_id, stot_comment_id, pc_comment_id, 
        pt_comment_id, laps_domain_file_id, asctime_id, fctimes_id, 
        level_id, s1hr_fcinv_id, stot_fcinv_id, pc_fcinv_id, pt_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  s1hr_valid_range[2];
   float  stot_valid_range[2];
   float  pc_valid_range[2];
   float  pt_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 4L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s1hr_id = ncvardef (cdfid, "s1hr", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   stot_id = ncvardef (cdfid, "stot", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pc_id = ncvardef (cdfid, "pc", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pt_id = ncvardef (cdfid, "pt", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s1hr_comment_id = ncvardef (cdfid, "s1hr_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   stot_comment_id = ncvardef (cdfid, "stot_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pc_comment_id = ncvardef (cdfid, "pc_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pt_comment_id = ncvardef (cdfid, "pt_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s1hr_fcinv_id = ncvardef (cdfid, "s1hr_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   stot_fcinv_id = ncvardef (cdfid, "stot_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pc_fcinv_id = ncvardef (cdfid, "pc_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pt_fcinv_id = ncvardef (cdfid, "pt_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, s1hr_id, "long_name", NC_CHAR, 26, (void *)"LAPS 60 minute snow accum.");
   ncattput (cdfid, s1hr_id, "units", NC_CHAR, 6, (void *)"meters");
   s1hr_valid_range[0] = -200;
   s1hr_valid_range[1] = 200;
   ncattput (cdfid, s1hr_id, "valid_range", NC_FLOAT, 2, (void *) s1hr_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, s1hr_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s1hr_id, "LAPS_var", NC_CHAR, 3, (void *)"S01");
   ncattput (cdfid, s1hr_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (cdfid, s1hr_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, stot_id, "long_name", NC_CHAR, 29, (void *)"storm total snow accumulation");
   ncattput (cdfid, stot_id, "units", NC_CHAR, 6, (void *)"meters");
   stot_valid_range[0] = -200;
   stot_valid_range[1] = 200;
   ncattput (cdfid, stot_id, "valid_range", NC_FLOAT, 2, (void *) stot_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, stot_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, stot_id, "LAPS_var", NC_CHAR, 3, (void *)"STO");
   ncattput (cdfid, stot_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (cdfid, stot_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, pc_id, "long_name", NC_CHAR, 29, (void *)"LAPS 60 minute precip. accum.");
   ncattput (cdfid, pc_id, "units", NC_CHAR, 6, (void *)"meters");
   pc_valid_range[0] = -200;
   pc_valid_range[1] = 200;
   ncattput (cdfid, pc_id, "valid_range", NC_FLOAT, 2, (void *) pc_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, pc_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, pc_id, "LAPS_var", NC_CHAR, 3, (void *)"R01");
   ncattput (cdfid, pc_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (cdfid, pc_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, pt_id, "long_name", NC_CHAR, 26, (void *)"storm total precip. accum.");
   ncattput (cdfid, pt_id, "units", NC_CHAR, 6, (void *)"meters");
   pt_valid_range[0] = -200;
   pt_valid_range[1] = 200;
   ncattput (cdfid, pt_id, "valid_range", NC_FLOAT, 2, (void *) pt_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, pt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, pt_id, "LAPS_var", NC_CHAR, 3, (void *)"RTO");
   ncattput (cdfid, pt_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (cdfid, pt_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, s1hr_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, stot_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, pc_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, pt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {24};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {4};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {4};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LPS */
#ifdef __STDC__
int cre_lps(char *fname) 		/* create fname */
#else
int cre_lps(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  ref_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, ref_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, ref_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  ref_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 21L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ref_id = ncvardef (cdfid, "ref", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ref_comment_id = ncvardef (cdfid, "ref_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ref_fcinv_id = ncvardef (cdfid, "ref_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, ref_id, "long_name", NC_CHAR, 23, (void *)"LAPS radar reflectivity");
   ncattput (cdfid, ref_id, "units", NC_CHAR, 3, (void *)"dBZ");
   ref_valid_range[0] = 0;
   ref_valid_range[1] = 100;
   ncattput (cdfid, ref_id, "valid_range", NC_FLOAT, 2, (void *) ref_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, ref_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, ref_id, "LAPS_var", NC_CHAR, 3, (void *)"REF");
   ncattput (cdfid, ref_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, ref_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, ref_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {21};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {21};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LRP */
#ifdef __STDC__
int cre_lrp(char *fname) 		/* create fname */
#else
int cre_lrp(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  icg_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, icg_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, icg_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  icg_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 21L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   icg_id = ncvardef (cdfid, "icg", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   icg_comment_id = ncvardef (cdfid, "icg_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   icg_fcinv_id = ncvardef (cdfid, "icg_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, icg_id, "long_name", NC_CHAR, 20, (void *)"icing severity index");
   ncattput (cdfid, icg_id, "units", NC_CHAR, 4, (void *)"none");
   icg_valid_range[0] = 0;
   icg_valid_range[1] = 100;
   ncattput (cdfid, icg_id, "valid_range", NC_FLOAT, 2, (void *) icg_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, icg_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, icg_id, "LAPS_var", NC_CHAR, 3, (void *)"LRP");
   ncattput (cdfid, icg_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, icg_id, "LAPS_units", NC_CHAR, 4, (void *)"NONE");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, icg_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {21};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {21};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LBA */
#ifdef __STDC__
int cre_lba(char *fname) 		/* create fname */
#else
int cre_lba(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  bz_id, bu_id, bv_id, bw_id, lvl_id, imax_id, jmax_id, kmax_id, 
        kdim_id, bz_comment_id, bu_comment_id, bv_comment_id, bw_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, bz_fcinv_id, 
        bu_fcinv_id, bv_fcinv_id, bw_fcinv_id, origin_id, model_id, 
        version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  bz_valid_range[2];
   float  bu_valid_range[2];
   float  bv_valid_range[2];
   float  bw_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 84L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   bz_id = ncvardef (cdfid, "bz", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   bu_id = ncvardef (cdfid, "bu", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   bv_id = ncvardef (cdfid, "bv", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   bw_id = ncvardef (cdfid, "bw", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   bz_comment_id = ncvardef (cdfid, "bz_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   bu_comment_id = ncvardef (cdfid, "bu_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   bv_comment_id = ncvardef (cdfid, "bv_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   bw_comment_id = ncvardef (cdfid, "bw_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   bz_fcinv_id = ncvardef (cdfid, "bz_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   bu_fcinv_id = ncvardef (cdfid, "bu_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   bv_fcinv_id = ncvardef (cdfid, "bv_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   bw_fcinv_id = ncvardef (cdfid, "bw_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, bz_id, "long_name", NC_CHAR, 26, (void *)"non-linear balanced height");
   ncattput (cdfid, bz_id, "units", NC_CHAR, 6, (void *)"meters");
   bz_valid_range[0] = -200;
   bz_valid_range[1] = 200;
   ncattput (cdfid, bz_id, "valid_range", NC_FLOAT, 2, (void *) bz_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, bz_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, bz_id, "LAPS_var", NC_CHAR, 2, (void *)"HT");
   ncattput (cdfid, bz_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (cdfid, bz_id, "LAPS_units", NC_CHAR, 6, (void *)"METERS");
   ncattput (cdfid, bu_id, "long_name", NC_CHAR, 33, (void *)"non-linear balanced eastward wind");
   ncattput (cdfid, bu_id, "units", NC_CHAR, 13, (void *)"meters/second");
   bu_valid_range[0] = -200;
   bu_valid_range[1] = 200;
   ncattput (cdfid, bu_id, "valid_range", NC_FLOAT, 2, (void *) bu_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, bu_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   char_val = 'U';
   ncattput (cdfid, bu_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, bu_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (cdfid, bu_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, bv_id, "long_name", NC_CHAR, 34, (void *)"non-linear balanced northward wind");
   ncattput (cdfid, bv_id, "units", NC_CHAR, 13, (void *)"meters/second");
   bv_valid_range[0] = -200;
   bv_valid_range[1] = 200;
   ncattput (cdfid, bv_id, "valid_range", NC_FLOAT, 2, (void *) bv_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, bv_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   char_val = 'V';
   ncattput (cdfid, bv_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, bv_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (cdfid, bv_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, bw_id, "long_name", NC_CHAR, 25, (void *)"non-linear balanced omega");
   ncattput (cdfid, bw_id, "units", NC_CHAR, 14, (void *)"pascals/second");
   bw_valid_range[0] = -20000;
   bw_valid_range[1] = 20000;
   ncattput (cdfid, bw_id, "valid_range", NC_FLOAT, 2, (void *) bw_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, bw_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, bw_id, "LAPS_var", NC_CHAR, 2, (void *)"OM");
   ncattput (cdfid, bw_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (cdfid, bw_id, "LAPS_units", NC_CHAR, 4, (void *)"PA/S");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, bz_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, bu_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, bv_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, bw_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {24};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {84};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {84};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LC3 */
#ifdef __STDC__
int cre_lc3(char *fname) 		/* create fname */
#else
int cre_lc3(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  camt_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        camt_comment_id, laps_domain_file_id, asctime_id, fctimes_id, 
        level_id, camt_fcinv_id, origin_id, model_id, version_id,
        num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  camt_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 42L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 42L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   camt_id = ncvardef (cdfid, "camt", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   camt_comment_id = ncvardef (cdfid, "camt_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   camt_fcinv_id = ncvardef (cdfid, "camt_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, camt_id, "long_name", NC_CHAR, 22, (void *)"fractional cloud cover");
   ncattput (cdfid, camt_id, "units", NC_CHAR, 10, (void *)"fractional");
   camt_valid_range[0] = 0;
   camt_valid_range[1] = 1;
   ncattput (cdfid, camt_id, "valid_range", NC_FLOAT, 2, (void *) camt_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, camt_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, camt_id, "LAPS_var", NC_CHAR, 3, (void *)"LC3");
   ncattput (cdfid, camt_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, camt_id, "LAPS_units", NC_CHAR, 10, (void *)"FRACTIONAL");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, camt_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {42};
    static short level[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {42};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {42};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LWC */
#ifdef __STDC__
int cre_lwc(char *fname) 		/* create fname */
#else
int cre_lwc(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  lwc_id, ice_id, pcn_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        lwc_comment_id, ice_comment_id, pcn_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, lwc_fcinv_id, ice_fcinv_id, 
        pcn_fcinv_id, origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  lwc_valid_range[2];
   float  ice_valid_range[2];
   float  pcn_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 63L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lwc_id = ncvardef (cdfid, "lwc", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ice_id = ncvardef (cdfid, "ice", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pcn_id = ncvardef (cdfid, "pcn", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lwc_comment_id = ncvardef (cdfid, "lwc_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ice_comment_id = ncvardef (cdfid, "ice_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pcn_comment_id = ncvardef (cdfid, "pcn_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lwc_fcinv_id = ncvardef (cdfid, "lwc_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ice_fcinv_id = ncvardef (cdfid, "ice_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pcn_fcinv_id = ncvardef (cdfid, "pcn_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, lwc_id, "long_name", NC_CHAR, 18, (void *)"cloud liquid water");
   ncattput (cdfid, lwc_id, "units", NC_CHAR, 14, (void *)"grams/meter**3");
   lwc_valid_range[0] = 0;
   lwc_valid_range[1] = 100;
   ncattput (cdfid, lwc_id, "valid_range", NC_FLOAT, 2, (void *) lwc_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, lwc_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, lwc_id, "LAPS_var", NC_CHAR, 3, (void *)"LWC");
   ncattput (cdfid, lwc_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, lwc_id, "LAPS_units", NC_CHAR, 6, (void *)"G/M**3");
   ncattput (cdfid, ice_id, "long_name", NC_CHAR, 9, (void *)"cloud ice");
   ncattput (cdfid, ice_id, "units", NC_CHAR, 14, (void *)"grams/meter**3");
   ice_valid_range[0] = 0;
   ice_valid_range[1] = 100;
   ncattput (cdfid, ice_id, "valid_range", NC_FLOAT, 2, (void *) ice_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, ice_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, ice_id, "LAPS_var", NC_CHAR, 3, (void *)"ICE");
   ncattput (cdfid, ice_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, ice_id, "LAPS_units", NC_CHAR, 6, (void *)"G/M**3");
   ncattput (cdfid, pcn_id, "long_name", NC_CHAR, 25, (void *)"hydrometeor concentration");
   ncattput (cdfid, pcn_id, "units", NC_CHAR, 14, (void *)"grams/meter**3");
   pcn_valid_range[0] = 0;
   pcn_valid_range[1] = 100;
   ncattput (cdfid, pcn_id, "valid_range", NC_FLOAT, 2, (void *) pcn_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, pcn_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, pcn_id, "LAPS_var", NC_CHAR, 3, (void *)"PCN");
   ncattput (cdfid, pcn_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, pcn_id, "LAPS_units", NC_CHAR, 6, (void *)"G/M**3");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, lwc_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, ice_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, pcn_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store num_variables */
    static long num_variables = {21};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {63};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {63};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LIL */
#ifdef __STDC__
int cre_lil(char *fname) 		/* create fname */
#else
int cre_lil(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  ilw_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, ilw_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, ilw_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  ilw_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 1L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ilw_id = ncvardef (cdfid, "ilw", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ilw_comment_id = ncvardef (cdfid, "ilw_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ilw_fcinv_id = ncvardef (cdfid, "ilw_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, ilw_id, "long_name", NC_CHAR, 23, (void *)"integrated liquid water");
   ncattput (cdfid, ilw_id, "units", NC_CHAR, 14, (void *)"grams/meter**2");
   ilw_valid_range[0] = 0;
   ilw_valid_range[1] = 0.1;
   ncattput (cdfid, ilw_id, "valid_range", NC_FLOAT, 2, (void *) ilw_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, ilw_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, ilw_id, "LAPS_var", NC_CHAR, 3, (void *)"LIL");
   ncattput (cdfid, ilw_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, ilw_id, "LAPS_units", NC_CHAR, 6, (void *)"G/M**2");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, ilw_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {1};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {1};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LCB */
#ifdef __STDC__
int cre_lcb(char *fname) 		/* create fname */
#else
int cre_lcb(fname) 		
char *fname;
#endif
{
   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  cbas_id, ctop_id, cce_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        cbas_comment_id, ctop_comment_id, cce_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, cbas_fcinv_id, ctop_fcinv_id, cce_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  cbas_valid_range[2];
   float  ctop_valid_range[2];
   float  cce_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 3L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   cbas_id = ncvardef (ncid, "cbas", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ctop_id = ncvardef (ncid, "ctop", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   cce_id = ncvardef (ncid, "cce", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   cbas_comment_id = ncvardef (ncid, "cbas_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ctop_comment_id = ncvardef (ncid, "ctop_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   cce_comment_id = ncvardef (ncid, "cce_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   cbas_fcinv_id = ncvardef (ncid, "cbas_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ctop_fcinv_id = ncvardef (ncid, "ctop_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   cce_fcinv_id = ncvardef (ncid, "cce_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, cbas_id, "long_name", NC_CHAR, 15, (void *)"LAPS cloud base");
   ncattput (ncid, cbas_id, "units", NC_CHAR, 6, (void *)"meters");
   cbas_valid_range[0] = 0;
   cbas_valid_range[1] = 40000;
   ncattput (ncid, cbas_id, "valid_range", NC_FLOAT, 2, (void *) cbas_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, cbas_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, cbas_id, "LAPS_var", NC_CHAR, 3, (void *)"LCB");
   ncattput (ncid, cbas_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, cbas_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, ctop_id, "long_name", NC_CHAR, 14, (void *)"LAPS cloud top");
   ncattput (ncid, ctop_id, "units", NC_CHAR, 6, (void *)"meters");
   ctop_valid_range[0] = 0;
   ctop_valid_range[1] = 40000;
   ncattput (ncid, ctop_id, "valid_range", NC_FLOAT, 2, (void *) ctop_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, ctop_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, ctop_id, "LAPS_var", NC_CHAR, 3, (void *)"LCT");
   ncattput (ncid, ctop_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, ctop_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, cce_id, "long_name", NC_CHAR, 18, (void *)"LAPS cloud ceiling");
   ncattput (ncid, cce_id, "units", NC_CHAR, 6, (void *)"meters");
   cce_valid_range[0] = 0;
   cce_valid_range[1] = 40000;
   ncattput (ncid, cce_id, "valid_range", NC_FLOAT, 2, (void *) cce_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, cce_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, cce_id, "LAPS_var", NC_CHAR, 3, (void *)"CCE");
   ncattput (ncid, cce_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'M';
   ncattput (ncid, cce_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, cbas_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, ctop_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, cce_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store num_variables */
    static long num_variables = {21};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {3};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {3};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LCT */
#ifdef __STDC__
int cre_lct(char *fname)
#else
int cre_lct(fname)              /* create fname */
char *fname;
#endif
{
   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  spt_id, ptt_id, sct_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        spt_comment_id, ptt_comment_id, sct_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, spt_fcinv_id, ptt_fcinv_id, sct_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  spt_valid_range[2];
   float  ptt_valid_range[2];
   float  sct_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 3L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   spt_id = ncvardef (ncid, "spt", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ptt_id = ncvardef (ncid, "ptt", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sct_id = ncvardef (ncid, "sct", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   spt_comment_id = ncvardef (ncid, "spt_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ptt_comment_id = ncvardef (ncid, "ptt_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sct_comment_id = ncvardef (ncid, "sct_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   spt_fcinv_id = ncvardef (ncid, "spt_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ptt_fcinv_id = ncvardef (ncid, "ptt_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sct_fcinv_id = ncvardef (ncid, "sct_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, spt_id, "long_name", NC_CHAR, 16, (void *)"LAPS precip type");
   ncattput (ncid, spt_id, "units", NC_CHAR, 4, (void *)"none");
   spt_valid_range[0] = 0;
   spt_valid_range[1] = 0.1;
   ncattput (ncid, spt_id, "valid_range", NC_FLOAT, 2, (void *) spt_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, spt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, spt_id, "LAPS_var", NC_CHAR, 3, (void *)"PTY");
   ncattput (ncid, spt_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, spt_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, ptt_id, "long_name", NC_CHAR, 42, (void *)"LAPS surface precip type-LL Refl Threshold");
   ncattput (ncid, ptt_id, "units", NC_CHAR, 4, (void *)"none");
   ptt_valid_range[0] = 0;
   ptt_valid_range[1] = 0.1;
   ncattput (ncid, ptt_id, "valid_range", NC_FLOAT, 2, (void *) ptt_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, ptt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, ptt_id, "LAPS_var", NC_CHAR, 3, (void *)"PTT");
   ncattput (ncid, ptt_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, ptt_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, sct_id, "long_name", NC_CHAR, 23, (void *)"LAPS surface cloud type");
   ncattput (ncid, sct_id, "units", NC_CHAR, 4, (void *)"none");
   sct_valid_range[0] = 0;
   sct_valid_range[1] = 16;
   ncattput (ncid, sct_id, "valid_range", NC_FLOAT, 2, (void *) sct_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, sct_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, sct_id, "LAPS_var", NC_CHAR, 3, (void *)"SCT");
   ncattput (ncid, sct_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, sct_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, spt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, ptt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, sct_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {21};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {3};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {3};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return (ncid);
}
/*************************************************************************/
/* file to create netCDF format file with extension LCV */
#ifdef __STDC__
int cre_lcv(char *fname) 		/* create fname */
#else
int cre_lcv(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim,
        asc_len_dim;
        

   /* variable ids */
   int  ccov_id, csc_id, cwt_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        ccov_comment_id, csc_comment_id, cwt_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, ccov_fcinv_id, csc_fcinv_id, 
        cwt_fcinv_id, origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short short_val;
   float float_val;

   /* attribute vectors */
   float  ccov_valid_range[2];
   float  csc_valid_range[2];
   float  cwt_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1);
   namelen_dim = ncdimdef(cdfid, "namelen", 132);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 3);
   var_len_dim = ncdimdef(cdfid, "var_len", 4);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18);

   /* define variables */
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ccov_id = ncvardef (cdfid, "ccov", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   csc_id = ncvardef (cdfid, "csc", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   cwt_id = ncvardef (cdfid, "cwt", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ccov_comment_id = ncvardef (cdfid, "ccov_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   csc_comment_id = ncvardef (cdfid, "csc_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   cwt_comment_id = ncvardef (cdfid, "cwt_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ccov_fcinv_id = ncvardef (cdfid, "ccov_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   csc_fcinv_id = ncvardef (cdfid, "csc_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   cwt_fcinv_id = ncvardef (cdfid, "cwt_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);


   /* assign attributes */
   ncattput (cdfid, ccov_id, "long_name", NC_CHAR, 16, (void *)"LAPS cloud cover");
   ncattput (cdfid, ccov_id, "units", NC_CHAR, 4, (void *)"none");
   ccov_valid_range[0] = 0;
   ccov_valid_range[1] = 0.1;
   ncattput (cdfid, ccov_id, "valid_range", NC_FLOAT, 2, (void *) ccov_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, ccov_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, ccov_id, "LAPS_var", NC_CHAR, 3, (void *)"LCV");
   ncattput (cdfid, ccov_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, ccov_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (cdfid, csc_id, "long_name", NC_CHAR, 33, (void *)"Cloud analysis implied snow cover");
   ncattput (cdfid, csc_id, "units", NC_CHAR, 4, (void *)"none");
   csc_valid_range[0] = 0.0;
   csc_valid_range[1] = 1.0;
   ncattput (cdfid, csc_id, "valid_range", NC_FLOAT, 2, (void *) csc_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, csc_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, csc_id, "LAPS_var", NC_CHAR, 3, (void *)"CSC");
   ncattput (cdfid, csc_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, csc_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (cdfid, cwt_id, "long_name", NC_CHAR, 27, (void *)"Clear Sky Water Temperature");
   ncattput (cdfid, cwt_id, "units", NC_CHAR, 4, (void *)"none");
   cwt_valid_range[0] = 100;
   cwt_valid_range[1] = 366;
   ncattput (cdfid, cwt_id, "valid_range", NC_FLOAT, 2, (void *) cwt_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, cwt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, cwt_id, "LAPS_var", NC_CHAR, 3, (void *)"CWT");
   ncattput (cdfid, cwt_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, cwt_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, ccov_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, csc_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (cdfid, cwt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
   static int fctimes_start[] = {0};
   static int fctimes_edges[] = {1};
   static short fctimes[] = {0};
   ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }
   
   {			/* store num_variables */
    static long num_variables = {21};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store origin */
   static int origin_start[] = {0};
   static int origin_edges[] = {49};
   static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
   ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
   static int model_start[] = {0};
   static int model_edges[] = {43};
   static char model[] = {"LAPS - Local Analysis and Prediction System"};
   ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {3};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {3};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LMD */
#ifdef __STDC__
int cre_lmd(char *fname) 		/* create fname */
#else
int cre_lmd(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  mcd_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, mcd_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, mcd_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  mcd_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 21L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mcd_id = ncvardef (cdfid, "mcd", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mcd_comment_id = ncvardef (cdfid, "mcd_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mcd_fcinv_id = ncvardef (cdfid, "mcd_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, mcd_id, "long_name", NC_CHAR, 38, (void *)"mean volume diameter of cloud droplets");
   ncattput (cdfid, mcd_id, "units", NC_CHAR, 6, (void *)"meters");
   mcd_valid_range[0] = 0;
   mcd_valid_range[1] = 100;
   ncattput (cdfid, mcd_id, "valid_range", NC_FLOAT, 2, (void *) mcd_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, mcd_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, mcd_id, "LAPS_var", NC_CHAR, 3, (void *)"LMD");
   ncattput (cdfid, mcd_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   char_val = 'M';
   ncattput (cdfid, mcd_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, mcd_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {21};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {21};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LCO */
#ifdef __STDC__
int cre_lco(char *fname) 		/* create fname */
#else
int cre_lco(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  cw_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, cw_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, cw_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  cw_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 21L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   cw_id = ncvardef (cdfid, "cw", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   cw_comment_id = ncvardef (cdfid, "cw_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   cw_fcinv_id = ncvardef (cdfid, "cw_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, cw_id, "long_name", NC_CHAR, 19, (void *)"cloud derived omega");
   ncattput (cdfid, cw_id, "units", NC_CHAR, 14, (void *)"Pascals/second");
   cw_valid_range[0] = 0;
   cw_valid_range[1] = 100;
   ncattput (cdfid, cw_id, "valid_range", NC_FLOAT, 2, (void *) cw_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, cw_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, cw_id, "LAPS_var", NC_CHAR, 3, (void *)"COM");
   ncattput (cdfid, cw_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, cw_id, "LAPS_units", NC_CHAR, 4, (void *)"PA/S");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, cw_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {21};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {21};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LTY */
#ifdef __STDC__
int cre_lty(char *fname) 		/* create fname */
#else
int cre_lty(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  ptyp_id, ctyp_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, ptyp_comment_id, ctyp_comment_id, laps_domain_file_id, asctime_id, fctimes_id, level_id, ptyp_fcinv_id, ctyp_fcinv_id, origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  ptyp_valid_range[2];
   float  ctyp_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 42L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ptyp_id = ncvardef (cdfid, "ptyp", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ctyp_id = ncvardef (cdfid, "ctyp", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ptyp_comment_id = ncvardef (cdfid, "ptyp_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ctyp_comment_id = ncvardef (cdfid, "ctyp_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ptyp_fcinv_id = ncvardef (cdfid, "ptyp_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ctyp_fcinv_id = ncvardef (cdfid, "ctyp_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, ptyp_id, "long_name", NC_CHAR, 18, (void *)"precipitation type");
   ncattput (cdfid, ptyp_id, "units", NC_CHAR, 4, (void *)"none");
   ptyp_valid_range[0] = -200;
   ptyp_valid_range[1] = 200;
   ncattput (cdfid, ptyp_id, "valid_range", NC_FLOAT, 2, (void *) ptyp_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, ptyp_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, ptyp_id, "LAPS_var", NC_CHAR, 3, (void *)"PTY");
   ncattput (cdfid, ptyp_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, ptyp_id, "LAPS_units", NC_CHAR, 4, (void *)"NONE");
   ncattput (cdfid, ctyp_id, "long_name", NC_CHAR, 10, (void *)"cloud type");
   ncattput (cdfid, ctyp_id, "units", NC_CHAR, 4, (void *)"none");
   ctyp_valid_range[0] = -200;
   ctyp_valid_range[1] = 200;
   ncattput (cdfid, ctyp_id, "valid_range", NC_FLOAT, 2, (void *) ctyp_valid_range);
   float_val = 1e+37;
   ncattput (cdfid, ctyp_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, ctyp_id, "LAPS_var", NC_CHAR, 3, (void *)"CTY");
   ncattput (cdfid, ctyp_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (cdfid, ctyp_id, "LAPS_units", NC_CHAR, 4, (void *)"NONE");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, ptyp_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, ctyp_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store num_variables */
    static long num_variables = {18};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {42};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {42};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LCP */
#ifdef __STDC__
int cre_lcp(char *fname) 		/* create fname */
#else
int cre_lcp(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  ccpc_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        ccpc_comment_id, laps_domain_file_id, asctime_id, fctimes_id, 
        level_id, ccpc_fcinv_id, origin_id, model_id, 
        version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  ccpc_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 21L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ccpc_id = ncvardef (cdfid, "ccpc", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ccpc_comment_id = ncvardef (cdfid, "ccpc_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ccpc_fcinv_id = ncvardef (cdfid, "ccpc_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, ccpc_id, "long_name", NC_CHAR, 37, (void *)"fractional cloud cover pressure coord");
   ncattput (cdfid, ccpc_id, "units", NC_CHAR, 10, (void *)"fractional");
   ccpc_valid_range[0] = 0;
   ccpc_valid_range[1] = 100;
   ncattput (cdfid, ccpc_id, "valid_range", NC_FLOAT, 2, (void *) ccpc_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, ccpc_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, ccpc_id, "LAPS_var", NC_CHAR, 3, (void *)"LCP");
   ncattput (cdfid, ccpc_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, ccpc_id, "LAPS_units", NC_CHAR, 10, (void *)"FRACTIONAL");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, ccpc_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {21};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {21};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LVD */
#ifdef __STDC__
int cre_lvd(char *fname)        	/* create fname */
#else
int cre_lvd(fname) 		
char *fname;
#endif
{   
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim,
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim,
        asc_len_dim;

   /* variable ids */
   int  s8w_id, s8c_id, svs_id, svn_id, alb_id, s3a_id, s3c_id, s4a_id, s4c_id,
        s5a_id, s5c_id, s8a_id, sca_id, scc_id, lvl_id, imax_id, jmax_id, 
        kmax_id, kdim_id, s8w_comment_id, s8c_comment_id, svs_comment_id, 
        svn_comment_id, alb_comment_id, s3a_comment_id, s3c_comment_id, 
        s4a_comment_id, s4c_comment_id, s5a_comment_id, s5c_comment_id, 
        s8a_comment_id, sca_comment_id, scc_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, s8w_fcinv_id, s8c_fcinv_id, 
        svs_fcinv_id, svn_fcinv_id, alb_fcinv_id, s3a_fcinv_id, s3c_fcinv_id, 
        s4a_fcinv_id, s4c_fcinv_id, s5a_fcinv_id, s5c_fcinv_id, s8a_fcinv_id, 
        sca_fcinv_id, scc_fcinv_id, origin_id, model_id, version_id, 
        num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  s8w_valid_range[2];
   float  s8c_valid_range[2];
   long  svs_valid_range[2];
   long  svn_valid_range[2];
   float  alb_valid_range[2];
   float  s3a_valid_range[2];
   float  s3c_valid_range[2];
   float  s4a_valid_range[2];
   float  s4c_valid_range[2];
   float  s5a_valid_range[2];
   float  s5c_valid_range[2];
   float  s8a_valid_range[2];
   float  sca_valid_range[2];
   float  scc_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1);
   namelen_dim = ncdimdef(cdfid, "namelen", 132);
   lat_dim = ncdimdef(cdfid, "lat", NY);
   lon_dim = ncdimdef(cdfid, "lon", NX);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 14);
   var_len_dim = ncdimdef(cdfid, "var_len", 4);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18);

   /* define variables */
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s8w_id = ncvardef (cdfid, "s8w", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s8c_id = ncvardef (cdfid, "s8c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   svs_id = ncvardef (cdfid, "svs", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   svn_id = ncvardef (cdfid, "svn", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   alb_id = ncvardef (cdfid, "alb", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s3a_id = ncvardef (cdfid, "s3a", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s3c_id = ncvardef (cdfid, "s3c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s4a_id = ncvardef (cdfid, "s4a", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s4c_id = ncvardef (cdfid, "s4c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s5a_id = ncvardef (cdfid, "s5a", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s5c_id = ncvardef (cdfid, "s5c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s8a_id = ncvardef (cdfid, "s8a", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sca_id = ncvardef (cdfid, "sca", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   scc_id = ncvardef (cdfid, "scc", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s8w_comment_id = ncvardef (cdfid, "s8w_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s8c_comment_id = ncvardef (cdfid, "s8c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   svs_comment_id = ncvardef (cdfid, "svs_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   svn_comment_id = ncvardef (cdfid, "svn_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   alb_comment_id = ncvardef (cdfid, "alb_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s3a_comment_id = ncvardef (cdfid, "s3a_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s3c_comment_id = ncvardef (cdfid, "s3c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s4a_comment_id = ncvardef (cdfid, "s4a_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s4c_comment_id = ncvardef (cdfid, "s4c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s5a_comment_id = ncvardef (cdfid, "s5a_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s5c_comment_id = ncvardef (cdfid, "s5c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s8a_comment_id = ncvardef (cdfid, "s8a_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sca_comment_id = ncvardef (cdfid, "sca_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   scc_comment_id = ncvardef (cdfid, "scc_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s8w_fcinv_id = ncvardef (cdfid, "s8w_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s8c_fcinv_id = ncvardef (cdfid, "s8c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   svs_fcinv_id = ncvardef (cdfid, "svs_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   svn_fcinv_id = ncvardef (cdfid, "svn_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   alb_fcinv_id = ncvardef (cdfid, "alb_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s3a_fcinv_id = ncvardef (cdfid, "s3a_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s3c_fcinv_id = ncvardef (cdfid, "s3c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s4a_fcinv_id = ncvardef (cdfid, "s4a_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s4c_fcinv_id = ncvardef (cdfid, "s4c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s5a_fcinv_id = ncvardef (cdfid, "s5a_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s5c_fcinv_id = ncvardef (cdfid, "s5c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s8a_fcinv_id = ncvardef (cdfid, "s8a_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sca_fcinv_id = ncvardef (cdfid, "sca_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   scc_fcinv_id = ncvardef (cdfid, "scc_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, dims);

      
   /* assign attributes */
   ncattput (cdfid, s8w_id, "long_name", NC_CHAR, 40, (void *)"goes IR band-8 bright temp warmest pixel");
   ncattput (cdfid, s8w_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s8w_valid_range[0] = 0;
   s8w_valid_range[1] = 400;
   ncattput (cdfid, s8w_id, "valid_range", NC_FLOAT, 2, (void *) s8w_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, s8w_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s8w_id, "LAPS_var", NC_CHAR, 3, (void *)"S8W");
   ncattput (cdfid, s8w_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, s8w_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, s8c_id, "long_name", NC_CHAR, 40, (void *)"goes IR band-8 bright temp coldest pixel¨");
   ncattput (cdfid, s8c_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s8c_valid_range[0] = 0;
   s8c_valid_range[1] = 400;
   ncattput (cdfid, s8c_id, "valid_range", NC_FLOAT, 2, (void *) s8c_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, s8c_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s8c_id, "LAPS_var", NC_CHAR, 3, (void *)"S8C");
   ncattput (cdfid, s8c_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, s8c_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, svs_id, "long_name", NC_CHAR, 28, (void *)"goes visible satellite - raw");
   ncattput (cdfid, svs_id, "units", NC_CHAR, 6, (void *)"counts");
   svs_valid_range[0] = 0;
   svs_valid_range[1] = 100;
   ncattput (cdfid, svs_id, "valid_range", NC_LONG, 2, (void *) svs_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, svs_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, svs_id, "LAPS_var", NC_CHAR, 3, (void *)"SVS");
   ncattput (cdfid, svs_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   ncattput (cdfid, svs_id, "LAPS_units", NC_CHAR, 6, (void *)"COUNTS");
   ncattput (cdfid, svn_id, "long_name", NC_CHAR, 35, (void *)"goes visible satellite - normalized");
   ncattput (cdfid, svn_id, "units", NC_CHAR, 6, (void *)"counts");
   svn_valid_range[0] = 0;
   svn_valid_range[1] = 100;
   ncattput (cdfid, svn_id, "valid_range", NC_LONG, 2, (void *) svn_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, svn_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, svn_id, "LAPS_var", NC_CHAR, 3, (void *)"SVN");
   ncattput (cdfid, svn_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   ncattput (cdfid, svn_id, "LAPS_units", NC_CHAR, 6, (void *)"COUNTS");
   ncattput (cdfid, alb_id, "long_name", NC_CHAR, 6, (void *)"albedo");
   ncattput (cdfid, alb_id, "units", NC_CHAR, 4, (void *)"none");
   alb_valid_range[0] = -20000;
   alb_valid_range[1] = 20000;
   ncattput (cdfid, alb_id, "valid_range", NC_FLOAT, 2, (void *) alb_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, alb_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, alb_id, "LAPS_var", NC_CHAR, 3, (void *)"ALB");
   ncattput (cdfid, alb_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   ncattput (cdfid, alb_id, "LAPS_units", NC_CHAR, 10, (void *)"          ");
   ncattput (cdfid, s3a_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-3 bright temp averaged");
   ncattput (cdfid, s3a_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s3a_valid_range[0] = 0;
   s3a_valid_range[1] = 400;
   ncattput (cdfid, s3a_id, "valid_range", NC_FLOAT, 2, (void *) s3a_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, s3a_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s3a_id, "LAPS_var", NC_CHAR, 3, (void *)"S3A");
   ncattput (cdfid, s3a_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, s3a_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, s3c_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-3 bright temp filtered");
   ncattput (cdfid, s3c_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s3c_valid_range[0] = 0;
   s3c_valid_range[1] = 400;
   ncattput (cdfid, s3c_id, "valid_range", NC_FLOAT, 2, (void *) s3c_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, s3c_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s3c_id, "LAPS_var", NC_CHAR, 3, (void *)"S3C");
   ncattput (cdfid, s3c_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, s3c_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, s4a_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-4 bright temp averaged");
   ncattput (cdfid, s4a_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s4a_valid_range[0] = 0;
   s4a_valid_range[1] = 400;
   ncattput (cdfid, s4a_id, "valid_range", NC_FLOAT, 2, (void *) s4a_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, s4a_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s4a_id, "LAPS_var", NC_CHAR, 3, (void *)"S4A");
   ncattput (cdfid, s4a_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, s4a_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, s4c_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-4 bright temp filtered");
   ncattput (cdfid, s4c_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s4c_valid_range[0] = 0;
   s4c_valid_range[1] = 400;
   ncattput (cdfid, s4c_id, "valid_range", NC_FLOAT, 2, (void *) s4c_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, s4c_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s4c_id, "LAPS_var", NC_CHAR, 3, (void *)"S4C");
   ncattput (cdfid, s4c_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, s4c_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, s5a_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-5 bright temp averaged");
   ncattput (cdfid, s5a_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s5a_valid_range[0] = 0;
   s5a_valid_range[1] = 400;
   ncattput (cdfid, s5a_id, "valid_range", NC_FLOAT, 2, (void *) s5a_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, s5a_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s5a_id, "LAPS_var", NC_CHAR, 3, (void *)"S5A");
   ncattput (cdfid, s5a_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, s5a_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, s5c_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-5 bright temp filtered");
   ncattput (cdfid, s5c_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s5c_valid_range[0] = 0;
   s5c_valid_range[1] = 400;
   ncattput (cdfid, s5c_id, "valid_range", NC_FLOAT, 2, (void *) s5c_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, s5c_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s5c_id, "LAPS_var", NC_CHAR, 3, (void *)"S5C");
   ncattput (cdfid, s5c_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, s5c_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, s8a_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-8 bright temp averaged");
   ncattput (cdfid, s8a_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s8a_valid_range[0] = 0;
   s8a_valid_range[1] = 400;
   ncattput (cdfid, s8a_id, "valid_range", NC_FLOAT, 2, (void *) s8a_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, s8a_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, s8a_id, "LAPS_var", NC_CHAR, 3, (void *)"S8A");
   ncattput (cdfid, s8a_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, s8a_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, sca_id, "long_name", NC_CHAR, 36, (void *)"goes IR band-12 bright temp averaged");
   ncattput (cdfid, sca_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   sca_valid_range[0] = 0;
   sca_valid_range[1] = 400;
   ncattput (cdfid, sca_id, "valid_range", NC_FLOAT, 2, (void *) sca_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, sca_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, sca_id, "LAPS_var", NC_CHAR, 3, (void *)"SCA");
   ncattput (cdfid, sca_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, sca_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, scc_id, "long_name", NC_CHAR, 36, (void *)"goes IR band-12 bright temp filtered");
   ncattput (cdfid, scc_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   scc_valid_range[0] = 0;
   scc_valid_range[1] = 400;
   ncattput (cdfid, scc_id, "valid_range", NC_FLOAT, 2, (void *) scc_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, scc_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, scc_id, "LAPS_var", NC_CHAR, 3, (void *)"SCC");
   ncattput (cdfid, scc_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113'; 
   ncattput (cdfid, scc_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0; 
   ncattput (cdfid, s8w_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, s8c_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, svs_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, svn_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, alb_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, s3a_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, s3c_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, s4a_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, s4c_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, s5a_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, s5c_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, s8a_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, sca_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, scc_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
   static int fctimes_start[] = {0};
   static int fctimes_edges[] = {1};
   static short fctimes[] = {0};
   ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {54};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store level */
   static int level_start[] = {0};
   static int level_edges[] = {1};
   static short level[] = {0};
   ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
   static int origin_start[] = {0};
   static int origin_edges[] = {49};
   static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
   ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
   static int model_start[] = {0};
   static int model_edges[] = {43};
   static char model[] = {"LAPS - Local Analysis and Prediction System"};
   ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }


   {			/* store kmax */
    static long kmax = {14};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {14};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LVE */
#ifdef __STDC__
int cre_lve(char *fname)        	/* create fname */
#else
int cre_lve(fname) 		
char *fname;
#endif
{   
   int  cdfid;			/* netCDF id */
          
   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim,
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim,
        asc_len_dim;
        
   /* variable ids */                             
   int  d8w_id, d8c_id, svs_id, svn_id, alb_id, d3a_id, d3c_id, d4a_id, d4c_id,
        d5a_id, d5c_id, d8a_id, dca_id, dcc_id, lvl_id, imax_id, jmax_id, 
        kmax_id, kdim_id, d8w_comment_id, d8c_comment_id, svs_comment_id, 
        svn_comment_id, alb_comment_id, d3a_comment_id, d3c_comment_id, 
        d4a_comment_id, d4c_comment_id, d5a_comment_id, d5c_comment_id, 
        d8a_comment_id, dca_comment_id, dcc_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, d8w_fcinv_id, d8c_fcinv_id, 
        svs_fcinv_id, svn_fcinv_id, alb_fcinv_id, d3a_fcinv_id, d3c_fcinv_id, 
        d4a_fcinv_id, d4c_fcinv_id, d5a_fcinv_id, d5c_fcinv_id, d8a_fcinv_id, 
        dca_fcinv_id, dcc_fcinv_id, origin_id, model_id, version_id, 
        num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  d8w_valid_range[2];
   float  d8c_valid_range[2];
   long  svs_valid_range[2];
   long  svn_valid_range[2];
   float  alb_valid_range[2];
   float  d3a_valid_range[2];
   float  d3c_valid_range[2];
   float  d4a_valid_range[2];
   float  d4c_valid_range[2];
   float  d5a_valid_range[2];
   float  d5c_valid_range[2];
   float  d8a_valid_range[2];
   float  dca_valid_range[2];
   float  dcc_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1);
   namelen_dim = ncdimdef(cdfid, "namelen", 132);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 14);
   var_len_dim = ncdimdef(cdfid, "var_len", 4);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18);

   /* define variables */
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d8w_id = ncvardef (cdfid, "d8w", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d8c_id = ncvardef (cdfid, "d8c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   svs_id = ncvardef (cdfid, "svs", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   svn_id = ncvardef (cdfid, "svn", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   alb_id = ncvardef (cdfid, "alb", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d3a_id = ncvardef (cdfid, "d3a", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d3c_id = ncvardef (cdfid, "d3c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d4a_id = ncvardef (cdfid, "d4a", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d4c_id = ncvardef (cdfid, "d4c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d5a_id = ncvardef (cdfid, "d5a", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d5c_id = ncvardef (cdfid, "d5c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   d8a_id = ncvardef (cdfid, "d8a", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   dca_id = ncvardef (cdfid, "dca", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   dcc_id = ncvardef (cdfid, "dcc", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d8w_comment_id = ncvardef (cdfid, "d8w_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d8c_comment_id = ncvardef (cdfid, "d8c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   svs_comment_id = ncvardef (cdfid, "svs_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   svn_comment_id = ncvardef (cdfid, "svn_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   alb_comment_id = ncvardef (cdfid, "alb_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d3a_comment_id = ncvardef (cdfid, "d3a_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d3c_comment_id = ncvardef (cdfid, "d3c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d4a_comment_id = ncvardef (cdfid, "d4a_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d4c_comment_id = ncvardef (cdfid, "d4c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d5a_comment_id = ncvardef (cdfid, "d5a_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d5c_comment_id = ncvardef (cdfid, "d5c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   d8a_comment_id = ncvardef (cdfid, "d8a_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   dca_comment_id = ncvardef (cdfid, "dca_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   dcc_comment_id = ncvardef (cdfid, "dcc_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d8w_fcinv_id = ncvardef (cdfid, "d8w_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d8c_fcinv_id = ncvardef (cdfid, "d8c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   svs_fcinv_id = ncvardef (cdfid, "svs_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   svn_fcinv_id = ncvardef (cdfid, "svn_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   alb_fcinv_id = ncvardef (cdfid, "alb_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d3a_fcinv_id = ncvardef (cdfid, "d3a_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d3c_fcinv_id = ncvardef (cdfid, "d3c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d4a_fcinv_id = ncvardef (cdfid, "d4a_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d4c_fcinv_id = ncvardef (cdfid, "d4c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d5a_fcinv_id = ncvardef (cdfid, "d5a_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d5c_fcinv_id = ncvardef (cdfid, "d5c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   d8a_fcinv_id = ncvardef (cdfid, "d8a_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dca_fcinv_id = ncvardef (cdfid, "dca_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dcc_fcinv_id = ncvardef (cdfid, "dcc_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, dims);

      
   /* assign attributes */
   ncattput (cdfid, d8w_id, "long_name", NC_CHAR, 40, (void *)"goes IR band-8 bright temp warmest pixel");
   ncattput (cdfid, d8w_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   d8w_valid_range[0] = 0;
   d8w_valid_range[1] = 400;
   ncattput (cdfid, d8w_id, "valid_range", NC_FLOAT, 2, (void *) d8w_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, d8w_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d8w_id, "LAPS_var", NC_CHAR, 3, (void *)"D8W");
   ncattput (cdfid, d8w_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, d8w_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, d8c_id, "long_name", NC_CHAR, 40, (void *)"goes IR band-8 bright temp coldest pixel¨");
   ncattput (cdfid, d8c_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   d8c_valid_range[0] = 0;
   d8c_valid_range[1] = 400;
   ncattput (cdfid, d8c_id, "valid_range", NC_FLOAT, 2, (void *) d8c_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, d8c_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d8c_id, "LAPS_var", NC_CHAR, 3, (void *)"D8C");
   ncattput (cdfid, d8c_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, d8c_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, svs_id, "long_name", NC_CHAR, 28, (void *)"goes visible satellite - raw");
   ncattput (cdfid, svs_id, "units", NC_CHAR, 6, (void *)"counts");
   svs_valid_range[0] = 0;
   svs_valid_range[1] = 100;
   ncattput (cdfid, svs_id, "valid_range", NC_LONG, 2, (void *) svs_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, svs_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, svs_id, "LAPS_var", NC_CHAR, 3, (void *)"SVS");
   ncattput (cdfid, svs_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   ncattput (cdfid, svs_id, "LAPS_units", NC_CHAR, 6, (void *)"COUNTS");
   ncattput (cdfid, svn_id, "long_name", NC_CHAR, 35, (void *)"goes visible satellite - normalized");
   ncattput (cdfid, svn_id, "units", NC_CHAR, 6, (void *)"counts");
   svn_valid_range[0] = 0;
   svn_valid_range[1] = 100;
   ncattput (cdfid, svn_id, "valid_range", NC_LONG, 2, (void *) svn_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, svn_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, svn_id, "LAPS_var", NC_CHAR, 3, (void *)"SVN");
   ncattput (cdfid, svn_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   ncattput (cdfid, svn_id, "LAPS_units", NC_CHAR, 6, (void *)"COUNTS");
   ncattput (cdfid, alb_id, "long_name", NC_CHAR, 6, (void *)"albedo");
   ncattput (cdfid, alb_id, "units", NC_CHAR, 4, (void *)"none");
   alb_valid_range[0] = -20000;
   alb_valid_range[1] = 20000;
   ncattput (cdfid, alb_id, "valid_range", NC_FLOAT, 2, (void *) alb_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, alb_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, alb_id, "LAPS_var", NC_CHAR, 3, (void *)"ALB");
   ncattput (cdfid, alb_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   ncattput (cdfid, alb_id, "LAPS_units", NC_CHAR, 10, (void *)"          ");
   ncattput (cdfid, d3a_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-3 bright temp averaged");
   ncattput (cdfid, d3a_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   d3a_valid_range[0] = 0;
   d3a_valid_range[1] = 400;
   ncattput (cdfid, d3a_id, "valid_range", NC_FLOAT, 2, (void *) d3a_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, d3a_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d3a_id, "LAPS_var", NC_CHAR, 3, (void *)"D3A");
   ncattput (cdfid, d3a_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, d3a_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, d3c_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-3 bright temp filtered");
   ncattput (cdfid, d3c_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   d3c_valid_range[0] = 0;
   d3c_valid_range[1] = 400;
   ncattput (cdfid, d3c_id, "valid_range", NC_FLOAT, 2, (void *) d3c_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, d3c_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d3c_id, "LAPS_var", NC_CHAR, 3, (void *)"D3C");
   ncattput (cdfid, d3c_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, d3c_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, d4a_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-4 bright temp averaged");
   ncattput (cdfid, d4a_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   d4a_valid_range[0] = 0;
   d4a_valid_range[1] = 400;
   ncattput (cdfid, d4a_id, "valid_range", NC_FLOAT, 2, (void *) d4a_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, d4a_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d4a_id, "LAPS_var", NC_CHAR, 3, (void *)"D4A");
   ncattput (cdfid, d4a_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, d4a_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, d4c_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-4 bright temp filtered");
   ncattput (cdfid, d4c_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   d4c_valid_range[0] = 0;
   d4c_valid_range[1] = 400;
   ncattput (cdfid, d4c_id, "valid_range", NC_FLOAT, 2, (void *) d4c_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, d4c_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d4c_id, "LAPS_var", NC_CHAR, 3, (void *)"D4C");
   ncattput (cdfid, d4c_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, d4c_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, d5a_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-5 bright temp averaged");
   ncattput (cdfid, d5a_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   d5a_valid_range[0] = 0;
   d5a_valid_range[1] = 400;
   ncattput (cdfid, d5a_id, "valid_range", NC_FLOAT, 2, (void *) d5a_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, d5a_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d5a_id, "LAPS_var", NC_CHAR, 3, (void *)"D5A");
   ncattput (cdfid, d5a_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, d5a_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, d5c_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-5 bright temp filtered");
   ncattput (cdfid, d5c_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   d5c_valid_range[0] = 0;
   d5c_valid_range[1] = 400;
   ncattput (cdfid, d5c_id, "valid_range", NC_FLOAT, 2, (void *) d5c_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, d5c_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d5c_id, "LAPS_var", NC_CHAR, 3, (void *)"D5C");
   ncattput (cdfid, d5c_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, d5c_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, d8a_id, "long_name", NC_CHAR, 35, (void *)"goes IR band-8 bright temp averaged");
   ncattput (cdfid, d8a_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   d8a_valid_range[0] = 0;
   d8a_valid_range[1] = 400;
   ncattput (cdfid, d8a_id, "valid_range", NC_FLOAT, 2, (void *) d8a_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, d8a_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, d8a_id, "LAPS_var", NC_CHAR, 3, (void *)"D8A");
   ncattput (cdfid, d8a_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, d8a_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, dca_id, "long_name", NC_CHAR, 36, (void *)"goes IR band-12 bright temp averaged");
   ncattput (cdfid, dca_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   dca_valid_range[0] = 0;
   dca_valid_range[1] = 400;
   ncattput (cdfid, dca_id, "valid_range", NC_FLOAT, 2, (void *) dca_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, dca_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, dca_id, "LAPS_var", NC_CHAR, 3, (void *)"DCA");
   ncattput (cdfid, dca_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113';
   ncattput (cdfid, dca_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, dcc_id, "long_name", NC_CHAR, 36, (void *)"goes IR band-12 bright temp filtered");
   ncattput (cdfid, dcc_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   dcc_valid_range[0] = 0;
   dcc_valid_range[1] = 400;
   ncattput (cdfid, dcc_id, "valid_range", NC_FLOAT, 2, (void *) dcc_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, dcc_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, dcc_id, "LAPS_var", NC_CHAR, 3, (void *)"DCC");
   ncattput (cdfid, dcc_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   char_val = '\113'; 
   ncattput (cdfid, dcc_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0; 
   ncattput (cdfid, d8w_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, d8c_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, svs_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, svn_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, alb_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, d3a_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, d3c_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, d4a_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, d4c_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, d5a_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, d5c_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, d8a_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, dca_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, dcc_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
   static int fctimes_start[] = {0};
   static int fctimes_edges[] = {1};
   static short fctimes[] = {0};
   ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {54};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store level */
   static int level_start[] = {0};
   static int level_edges[] = {1};
   static short level[] = {0};
   ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
   static int origin_start[] = {0};
   static int origin_edges[] = {49};
   static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
   ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
   static int model_start[] = {0};
   static int model_edges[] = {43};
   static char model[] = {"LAPS - Local Analysis and Prediction System"};
   ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }


   {			/* store kmax */
    static long kmax = {14};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {14};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LMA */
#ifdef __STDC__
int cre_lma(char *fname) 		/* create fname */
#else
int cre_lma(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  maz_id, mat_id, marh_id, mau_id, mav_id, lvl_id, imax_id, 
        jmax_id, kmax_id, kdim_id, maz_comment_id, mat_comment_id, 
        marh_comment_id, mau_comment_id, mav_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, 
        maz_fcinv_id, mat_fcinv_id, marh_fcinv_id, mau_fcinv_id, 
        mav_fcinv_id, origin_id, model_id, version_id, 
        num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  maz_valid_range[2];
   float  mat_valid_range[2];
   float  marh_valid_range[2];
   float  mau_valid_range[2];
   float  mav_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", 11L);
   lon_dim = ncdimdef(cdfid, "lon", 11L);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 105L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   maz_id = ncvardef (cdfid, "maz", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mat_id = ncvardef (cdfid, "mat", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   marh_id = ncvardef (cdfid, "marh", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mau_id = ncvardef (cdfid, "mau", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mav_id = ncvardef (cdfid, "mav", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   maz_comment_id = ncvardef (cdfid, "maz_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mat_comment_id = ncvardef (cdfid, "mat_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   marh_comment_id = ncvardef (cdfid, "marh_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mau_comment_id = ncvardef (cdfid, "mau_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mav_comment_id = ncvardef (cdfid, "mav_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   maz_fcinv_id = ncvardef (cdfid, "maz_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mat_fcinv_id = ncvardef (cdfid, "mat_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   marh_fcinv_id = ncvardef (cdfid, "marh_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mau_fcinv_id = ncvardef (cdfid, "mau_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mav_fcinv_id = ncvardef (cdfid, "mav_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, maz_id, "long_name", NC_CHAR, 12, (void *)"MAPS heights");
   ncattput (cdfid, maz_id, "units", NC_CHAR, 6, (void *)"meters");
   maz_valid_range[0] = -200;
   maz_valid_range[1] = 200;
   ncattput (cdfid, maz_id, "valid_range", NC_FLOAT, 2, (void *) maz_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, maz_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, maz_id, "LAPS_var", NC_CHAR, 2, (void *)"HT");
   ncattput (cdfid, maz_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, maz_id, "LAPS_units", NC_CHAR, 6, (void *)"METERS");
   ncattput (cdfid, mat_id, "long_name", NC_CHAR, 16, (void *)"MAPS temperature");
   ncattput (cdfid, mat_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   mat_valid_range[0] = -200;
   mat_valid_range[1] = 200;
   ncattput (cdfid, mat_id, "valid_range", NC_FLOAT, 2, (void *) mat_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, mat_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   char_val = 'T';
   ncattput (cdfid, mat_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, mat_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, mat_id, "LAPS_units", NC_CHAR, 6, (void *)"KELVIN");
   ncattput (cdfid, marh_id, "long_name", NC_CHAR, 22, (void *)"MAPS relative humidity");
   ncattput (cdfid, marh_id, "units", NC_CHAR, 4, (void *)"none");
   marh_valid_range[0] = -200;
   marh_valid_range[1] = 200;
   ncattput (cdfid, marh_id, "valid_range", NC_FLOAT, 2, (void *) marh_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, marh_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, marh_id, "LAPS_var", NC_CHAR, 2, (void *)"RH");
   ncattput (cdfid, marh_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, marh_id, "LAPS_units", NC_CHAR, 4, (void *)"NONE");
   ncattput (cdfid, mau_id, "long_name", NC_CHAR, 21, (void *)"MAPS u wind component");
   ncattput (cdfid, mau_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mau_valid_range[0] = -200;
   mau_valid_range[1] = 200;
   ncattput (cdfid, mau_id, "valid_range", NC_FLOAT, 2, (void *) mau_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, mau_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   char_val = 'U';
   ncattput (cdfid, mau_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, mau_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, mau_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, mav_id, "long_name", NC_CHAR, 21, (void *)"MAPS v wind component");
   ncattput (cdfid, mav_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mav_valid_range[0] = -20000;
   mav_valid_range[1] = 20000;
   ncattput (cdfid, mav_id, "valid_range", NC_FLOAT, 2, (void *) mav_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, mav_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   char_val = 'V';
   ncattput (cdfid, mav_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, mav_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, mav_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, maz_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, mat_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, marh_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, mau_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, mav_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {27};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {11};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {11};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {105};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {105};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LMF */
#ifdef __STDC__
int cre_lmf(char *fname) 		/* create fname */
#else
int cre_lmf(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  maz_id, mat_id, marh_id, mau_id, mav_id, lvl_id, imax_id, 
        jmax_id, kmax_id, kdim_id, maz_comment_id, mat_comment_id, 
        marh_comment_id, mau_comment_id, mav_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, 
        maz_fcinv_id, mat_fcinv_id, marh_fcinv_id, mau_fcinv_id, 
        mav_fcinv_id, origin_id, model_id, version_id, 
        num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  maz_valid_range[2];
   float  mat_valid_range[2];
   float  marh_valid_range[2];
   float  mau_valid_range[2];
   float  mav_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 21L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", 11L);
   lon_dim = ncdimdef(cdfid, "lon", 11L);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 105L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   maz_id = ncvardef (cdfid, "maz", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mat_id = ncvardef (cdfid, "mat", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   marh_id = ncvardef (cdfid, "marh", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mau_id = ncvardef (cdfid, "mau", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mav_id = ncvardef (cdfid, "mav", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   maz_comment_id = ncvardef (cdfid, "maz_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mat_comment_id = ncvardef (cdfid, "mat_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   marh_comment_id = ncvardef (cdfid, "marh_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mau_comment_id = ncvardef (cdfid, "mau_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mav_comment_id = ncvardef (cdfid, "mav_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   maz_fcinv_id = ncvardef (cdfid, "maz_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mat_fcinv_id = ncvardef (cdfid, "mat_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   marh_fcinv_id = ncvardef (cdfid, "marh_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mau_fcinv_id = ncvardef (cdfid, "mau_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mav_fcinv_id = ncvardef (cdfid, "mav_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, maz_id, "long_name", NC_CHAR, 12, (void *)"MAPS heights");
   ncattput (cdfid, maz_id, "units", NC_CHAR, 6, (void *)"meters");
   maz_valid_range[0] = -200;
   maz_valid_range[1] = 200;
   ncattput (cdfid, maz_id, "valid_range", NC_FLOAT, 2, (void *) maz_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, maz_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, maz_id, "LAPS_var", NC_CHAR, 2, (void *)"HT");
   ncattput (cdfid, maz_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, maz_id, "LAPS_units", NC_CHAR, 6, (void *)"METERS");
   ncattput (cdfid, mat_id, "long_name", NC_CHAR, 16, (void *)"MAPS temperature");
   ncattput (cdfid, mat_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   mat_valid_range[0] = -200;
   mat_valid_range[1] = 200;
   ncattput (cdfid, mat_id, "valid_range", NC_FLOAT, 2, (void *) mat_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, mat_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   char_val = 'T';
   ncattput (cdfid, mat_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, mat_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, mat_id, "LAPS_units", NC_CHAR, 6, (void *)"KELVIN");
   ncattput (cdfid, marh_id, "long_name", NC_CHAR, 22, (void *)"MAPS relative humidity");
   ncattput (cdfid, marh_id, "units", NC_CHAR, 4, (void *)"none");
   marh_valid_range[0] = -200;
   marh_valid_range[1] = 200;
   ncattput (cdfid, marh_id, "valid_range", NC_FLOAT, 2, (void *) marh_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, marh_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, marh_id, "LAPS_var", NC_CHAR, 2, (void *)"RH");
   ncattput (cdfid, marh_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, marh_id, "LAPS_units", NC_CHAR, 4, (void *)"NONE");
   ncattput (cdfid, mau_id, "long_name", NC_CHAR, 21, (void *)"MAPS u wind component");
   ncattput (cdfid, mau_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mau_valid_range[0] = -200;
   mau_valid_range[1] = 200;
   ncattput (cdfid, mau_id, "valid_range", NC_FLOAT, 2, (void *) mau_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, mau_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   char_val = 'U';
   ncattput (cdfid, mau_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, mau_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, mau_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, mav_id, "long_name", NC_CHAR, 21, (void *)"MAPS v wind component");
   ncattput (cdfid, mav_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mav_valid_range[0] = -20000;
   mav_valid_range[1] = 20000;
   ncattput (cdfid, mav_id, "valid_range", NC_FLOAT, 2, (void *) mav_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, mav_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   char_val = 'V';
   ncattput (cdfid, mav_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, mav_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (cdfid, mav_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, maz_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, mat_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, marh_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, mau_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, mav_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {27};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {11};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {11};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {105};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {105};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension Z02 */
#ifdef __STDC__
int cre_z02(char *fname) 		/* create fname */
#else
int cre_z02(fname) 		
char *fname;
#endif
{
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  s8w_id, s8c_id, svs_id, svn_id, alb_id, lvl_id, imax_id, 
        jmax_id, kmax_id, kdim_id, s8w_comment_id, s8c_comment_id, 
        svs_comment_id, svn_comment_id, alb_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, s8w_fcinv_id, s8c_fcinv_id, 
        svs_fcinv_id, svn_fcinv_id, alb_fcinv_id, origin_id, model_id, 
        version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   long  long_val;
   double  double_val;

   /* attribute vectors */
   float  s8w_valid_range[2];
   float  s8c_valid_range[2];
   float  svs_valid_range[2];
   float  svn_valid_range[2];
   float  alb_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1L);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1L);
   namelen_dim = ncdimdef(cdfid, "namelen", 132L);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 5L);
   var_len_dim = ncdimdef(cdfid, "var_len", 4L);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5L);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11L);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126L);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12L);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s8w_id = ncvardef (cdfid, "s8w", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s8c_id = ncvardef (cdfid, "s8c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   svs_id = ncvardef (cdfid, "svs", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   svn_id = ncvardef (cdfid, "svn", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   alb_id = ncvardef (cdfid, "alb", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s8w_comment_id = ncvardef (cdfid, "s8w_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s8c_comment_id = ncvardef (cdfid, "s8c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   svs_comment_id = ncvardef (cdfid, "svs_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   svn_comment_id = ncvardef (cdfid, "svn_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   alb_comment_id = ncvardef (cdfid, "alb_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s8w_fcinv_id = ncvardef (cdfid, "s8w_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s8c_fcinv_id = ncvardef (cdfid, "s8c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   svs_fcinv_id = ncvardef (cdfid, "svs_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   svn_fcinv_id = ncvardef (cdfid, "svn_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   alb_fcinv_id = ncvardef (cdfid, "alb_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, s8w_id, "long_name", NC_CHAR, 40, (void *)"goes IR band-8 bright temp warmest pixel");
   ncattput (cdfid, s8w_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s8w_valid_range[0] = -200;
   s8w_valid_range[1] = 200;
   ncattput (cdfid, s8w_id, "valid_range", NC_FLOAT, 2, (void *) s8w_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, s8w_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, s8w_id, "LAPS_var", NC_CHAR, 3, (void *)"S8W");
   ncattput (cdfid, s8w_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (cdfid, s8w_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, s8c_id, "long_name", NC_CHAR, 40, (void *)"goes IR band-8 bright temp coldest pixel");
   ncattput (cdfid, s8c_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   s8c_valid_range[0] = -200;
   s8c_valid_range[1] = 200;
   ncattput (cdfid, s8c_id, "valid_range", NC_FLOAT, 2, (void *) s8c_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, s8c_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, s8c_id, "LAPS_var", NC_CHAR, 3, (void *)"S8C");
   ncattput (cdfid, s8c_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (cdfid, s8c_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, svs_id, "long_name", NC_CHAR, 28, (void *)"goes visible satellite - raw");
   ncattput (cdfid, svs_id, "units", NC_CHAR, 6, (void *)"counts");
   svs_valid_range[0] = -200;
   svs_valid_range[1] = 200;
   ncattput (cdfid, svs_id, "valid_range", NC_FLOAT, 2, (void *) svs_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, svs_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, svs_id, "LAPS_var", NC_CHAR, 3, (void *)"SVS");
   ncattput (cdfid, svs_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, svs_id, "LAPS_units", NC_CHAR, 6, (void *)"COUNTS");
   ncattput (cdfid, svn_id, "long_name", NC_CHAR, 35, (void *)"goes visible satellite - normalized");
   ncattput (cdfid, svn_id, "units", NC_CHAR, 6, (void *)"counts");
   svn_valid_range[0] = -200;
   svn_valid_range[1] = 200;
   ncattput (cdfid, svn_id, "valid_range", NC_FLOAT, 2, (void *) svn_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, svn_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, svn_id, "LAPS_var", NC_CHAR, 3, (void *)"SVN");
   ncattput (cdfid, svn_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, svn_id, "LAPS_units", NC_CHAR, 6, (void *)"COUNTS");
   ncattput (cdfid, alb_id, "long_name", NC_CHAR, 6, (void *)"albedo");
   ncattput (cdfid, alb_id, "units", NC_CHAR, 4, (void *)"none");
   alb_valid_range[0] = -20000;
   alb_valid_range[1] = 20000;
   ncattput (cdfid, alb_id, "valid_range", NC_FLOAT, 2, (void *) alb_valid_range);
   double_val = 1e+37;
   ncattput (cdfid, alb_id, "_FillValue", NC_DOUBLE, 1,(void *) &double_val);
   ncattput (cdfid, alb_id, "LAPS_var", NC_CHAR, 3, (void *)"ALB");
   ncattput (cdfid, alb_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (cdfid, alb_id, "LAPS_units", NC_CHAR, 10, (void *)"          ");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   long_val = 0;
   ncattput (cdfid, s8w_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, s8c_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, svs_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, svn_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);
   long_val = 0;
   ncattput (cdfid, alb_fcinv_id, "_FillValue", NC_LONG, 1,(void *) &long_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {27};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {5};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {5};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension RAM */
#ifdef __STDC__
int cre_ram(char *fname)                /* create fname */
#else
int cre_ram(fname)
char *fname;
#endif
{

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  rz_id, ru_id, rv_id, rw_id, rt_id, rmr_id, rh3_id, lwc_id, ice_id, rai_id, 
        sno_id, pic_id, ref_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        rz_comment_id, ru_comment_id, rv_comment_id, rw_comment_id, rt_comment_id, 
        rmr_comment_id, rh3_comment_id, lwc_comment_id, ice_comment_id, 
        rai_comment_id, sno_comment_id, pic_comment_id, ref_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, rz_fcinv_id, 
        ru_fcinv_id, rv_fcinv_id, rw_fcinv_id, rt_fcinv_id, rmr_fcinv_id, 
        rh3_fcinv_id, lwc_fcinv_id, ice_fcinv_id, rai_fcinv_id, sno_fcinv_id, 
        pic_fcinv_id, ref_fcinv_id, origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  rz_valid_range[2];
   float  ru_valid_range[2];
   float  rv_valid_range[2];
   float  rw_valid_range[2];
   float  rt_valid_range[2];
   float  rmr_valid_range[2];
   float  rh3_valid_range[2];
   float  lwc_valid_range[2];
   float  ice_valid_range[2];
   float  rai_valid_range[2];
   float  sno_valid_range[2];
   float  pic_valid_range[2];
   float  ref_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 21L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 273L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rz_id = ncvardef (ncid, "rz", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ru_id = ncvardef (ncid, "ru", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rv_id = ncvardef (ncid, "rv", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rw_id = ncvardef (ncid, "rw", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rt_id = ncvardef (ncid, "rt", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rmr_id = ncvardef (ncid, "rmr", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rh3_id = ncvardef (ncid, "rh3", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lwc_id = ncvardef (ncid, "lwc", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ice_id = ncvardef (ncid, "ice", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rai_id = ncvardef (ncid, "rai", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sno_id = ncvardef (ncid, "sno", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pic_id = ncvardef (ncid, "pic", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ref_id = ncvardef (ncid, "ref", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rz_comment_id = ncvardef (ncid, "rz_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ru_comment_id = ncvardef (ncid, "ru_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rv_comment_id = ncvardef (ncid, "rv_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rw_comment_id = ncvardef (ncid, "rw_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rt_comment_id = ncvardef (ncid, "rt_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rmr_comment_id = ncvardef (ncid, "rmr_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rh3_comment_id = ncvardef (ncid, "rh3_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lwc_comment_id = ncvardef (ncid, "lwc_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ice_comment_id = ncvardef (ncid, "ice_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rai_comment_id = ncvardef (ncid, "rai_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sno_comment_id = ncvardef (ncid, "sno_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pic_comment_id = ncvardef (ncid, "pic_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ref_comment_id = ncvardef (ncid, "ref_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rz_fcinv_id = ncvardef (ncid, "rz_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ru_fcinv_id = ncvardef (ncid, "ru_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rv_fcinv_id = ncvardef (ncid, "rv_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rw_fcinv_id = ncvardef (ncid, "rw_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rt_fcinv_id = ncvardef (ncid, "rt_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rmr_fcinv_id = ncvardef (ncid, "rmr_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rh3_fcinv_id = ncvardef (ncid, "rh3_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lwc_fcinv_id = ncvardef (ncid, "lwc_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ice_fcinv_id = ncvardef (ncid, "ice_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rai_fcinv_id = ncvardef (ncid, "rai_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sno_fcinv_id = ncvardef (ncid, "sno_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pic_fcinv_id = ncvardef (ncid, "pic_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ref_fcinv_id = ncvardef (ncid, "ref_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, rz_id, "long_name", NC_CHAR, 22, (void *)"LAPS Fcst model height");
   ncattput (ncid, rz_id, "units", NC_CHAR, 6, (void *)"meters");
   rz_valid_range[0] = 0;
   rz_valid_range[1] = 100000;
   ncattput (ncid, rz_id, "valid_range", NC_FLOAT, 2, (void *) rz_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rz_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rz_id, "LAPS_var", NC_CHAR, 2, (void *)"HT");
   ncattput (ncid, rz_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rz_id, "LAPS_units", NC_CHAR, 6, (void *)"METERS");
   ncattput (ncid, ru_id, "long_name", NC_CHAR, 23, (void *)"LAPS Fcst eastward wind");
   ncattput (ncid, ru_id, "units", NC_CHAR, 13, (void *)"meters/second");
   ru_valid_range[0] = -200;
   ru_valid_range[1] = 200;
   ncattput (ncid, ru_id, "valid_range", NC_FLOAT, 2, (void *) ru_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, ru_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, ru_id, "LAPS_var", NC_CHAR, 2, (void *)"U3");
   ncattput (ncid, ru_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, ru_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, rv_id, "long_name", NC_CHAR, 24, (void *)"LAPS Fcst northward wind");
   ncattput (ncid, rv_id, "units", NC_CHAR, 13, (void *)"meters/second");
   rv_valid_range[0] = -200;
   rv_valid_range[1] = 200;
   ncattput (ncid, rv_id, "valid_range", NC_FLOAT, 2, (void *) rv_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rv_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rv_id, "LAPS_var", NC_CHAR, 2, (void *)"V3");
   ncattput (ncid, rv_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rv_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, rw_id, "long_name", NC_CHAR, 15, (void *)"LAPS Fcst omega");
   ncattput (ncid, rw_id, "units", NC_CHAR, 14, (void *)"pascals/second");
   rw_valid_range[0] = -20000;
   rw_valid_range[1] = 20000;
   ncattput (ncid, rw_id, "valid_range", NC_FLOAT, 2, (void *) rw_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rw_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rw_id, "LAPS_var", NC_CHAR, 2, (void *)"OM");
   ncattput (ncid, rw_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rw_id, "LAPS_units", NC_CHAR, 4, (void *)"PA/S");
   ncattput (ncid, rt_id, "long_name", NC_CHAR, 21, (void *)"LAPS Fcst temperature");
   ncattput (ncid, rt_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   rt_valid_range[0] = 0;
   rt_valid_range[1] = 500;
   ncattput (ncid, rt_id, "valid_range", NC_FLOAT, 2, (void *) rt_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rt_id, "LAPS_var", NC_CHAR, 2, (void *)"T3");
   ncattput (ncid, rt_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   char_val = 'K';
   ncattput (ncid, rt_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rmr_id, "long_name", NC_CHAR, 22, (void *)"LAPS Fcst mixing ratio");
   ncattput (ncid, rmr_id, "units", NC_CHAR, 5, (void *)"kg/kg");
   rmr_valid_range[0] = 0;
   rmr_valid_range[1] = 0.1;
   ncattput (ncid, rmr_id, "valid_range", NC_FLOAT, 2, (void *) rmr_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rmr_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rmr_id, "LAPS_var", NC_CHAR, 3, (void *)"RMR");
   ncattput (ncid, rmr_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, rmr_id, "LAPS_units", NC_CHAR, 10, (void *)"          ");
   ncattput (ncid, rh3_id, "long_name", NC_CHAR, 27, (void *)"LAPS Fcst relative humidity");
   ncattput (ncid, rh3_id, "units", NC_CHAR, 7, (void *)"percent");
   rh3_valid_range[0] = 0;
   rh3_valid_range[1] = 100;
   ncattput (ncid, rh3_id, "valid_range", NC_FLOAT, 2, (void *) rh3_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rh3_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rh3_id, "LAPS_var", NC_CHAR, 3, (void *)"RH3");
   ncattput (ncid, rh3_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, rh3_id, "LAPS_units", NC_CHAR, 7, (void *)"PERCENT");
   ncattput (ncid, lwc_id, "long_name", NC_CHAR, 28, (void *)"LAPS Fcst cloud liquid water");
   ncattput (ncid, lwc_id, "units", NC_CHAR, 18, (void *)"kilograms/meter**3");
   lwc_valid_range[0] = 0;
   lwc_valid_range[1] = 0.1;
   ncattput (ncid, lwc_id, "valid_range", NC_FLOAT, 2, (void *) lwc_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lwc_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lwc_id, "LAPS_var", NC_CHAR, 3, (void *)"LWC");
   ncattput (ncid, lwc_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, lwc_id, "LAPS_units", NC_CHAR, 7, (void *)"KG/M**3");
   ncattput (ncid, ice_id, "long_name", NC_CHAR, 19, (void *)"LAPS Fcst cloud ice");
   ncattput (ncid, ice_id, "units", NC_CHAR, 18, (void *)"kilograms/meter**3");
   ice_valid_range[0] = 0;
   ice_valid_range[1] = 0.1;
   ncattput (ncid, ice_id, "valid_range", NC_FLOAT, 2, (void *) ice_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, ice_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, ice_id, "LAPS_var", NC_CHAR, 3, (void *)"ICE");
   ncattput (ncid, ice_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, ice_id, "LAPS_units", NC_CHAR, 7, (void *)"KG/M**3");
   ncattput (ncid, rai_id, "long_name", NC_CHAR, 14, (void *)"LAPS Fcst rain");
   ncattput (ncid, rai_id, "units", NC_CHAR, 11, (void *)"kg/meter**3");
   rai_valid_range[0] = 0;
   rai_valid_range[1] = 0.1;
   ncattput (ncid, rai_id, "valid_range", NC_FLOAT, 2, (void *) rai_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rai_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rai_id, "LAPS_var", NC_CHAR, 3, (void *)"RAI");
   ncattput (ncid, rai_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, rai_id, "LAPS_units", NC_CHAR, 7, (void *)"KG/M**3");
   ncattput (ncid, sno_id, "long_name", NC_CHAR, 14, (void *)"LAPS Fcst snow");
   ncattput (ncid, sno_id, "units", NC_CHAR, 11, (void *)"kg/meter**3");
   sno_valid_range[0] = 0;
   sno_valid_range[1] = 0.1;
   ncattput (ncid, sno_id, "valid_range", NC_FLOAT, 2, (void *) sno_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, sno_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, sno_id, "LAPS_var", NC_CHAR, 3, (void *)"SNO");
   ncattput (ncid, sno_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, sno_id, "LAPS_units", NC_CHAR, 7, (void *)"KG/M**3");
   ncattput (ncid, pic_id, "long_name", NC_CHAR, 20, (void *)"LAPS Fcst precip ice");
   ncattput (ncid, pic_id, "units", NC_CHAR, 11, (void *)"kg/meter**3");
   pic_valid_range[0] = 0;
   pic_valid_range[1] = 0.1;
   ncattput (ncid, pic_id, "valid_range", NC_FLOAT, 2, (void *) pic_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, pic_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, pic_id, "LAPS_var", NC_CHAR, 3, (void *)"PIC");
   ncattput (ncid, pic_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, pic_id, "LAPS_units", NC_CHAR, 7, (void *)"KG/M**3");
   ncattput (ncid, ref_id, "long_name", NC_CHAR, 20, (void *)"LAPS Fcst radar refl");
   ncattput (ncid, ref_id, "units", NC_CHAR, 3, (void *)"dBZ");
   ref_valid_range[0] = 0;
   ref_valid_range[1] = 100;
   ncattput (ncid, ref_id, "valid_range", NC_FLOAT, 2, (void *) ref_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, ref_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, ref_id, "LAPS_var", NC_CHAR, 3, (void *)"REF");
   ncattput (ncid, ref_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, ref_id, "LAPS_units", NC_CHAR, 3, (void *)"dBZ");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, rz_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, ru_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rv_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rw_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rmr_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rh3_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, lwc_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, ice_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rai_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, sno_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, pic_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, ref_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store num_variables */
    static long num_variables = {51};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {132};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {132};
    static char model[] = {"LAPS Forecast Model - see comment fields for model and version"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {273};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {273};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/*************************************************************************/
/* file to create netCDF format file with extension RSF */
#ifdef __STDC__
int cre_rsf(char *fname)                /* create fname */
#else
int cre_rsf(fname)
char *fname;
#endif
{                       /* create rsf.cdf */

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  rus_id, rvs_id, rts_id, rps_id, rtd_id, rh_id, lcb_id, lct_id, msl_id, 
        lil_id, tpw_id, r01_id, rto_id, s01_id, sto_id, th_id, the_id, pbe_id, 
        nbe_id, ps_id, cce_id, vis_id, lcv_id, lmt_id, spt_id, lhe_id, li_id, 
        hi_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        rus_comment_id, rvs_comment_id, rts_comment_id, rps_comment_id, rtd_comment_id, 
        rh_comment_id, lcb_comment_id, lct_comment_id, msl_comment_id, lil_comment_id, 
        tpw_comment_id, r01_comment_id, rto_comment_id, s01_comment_id, sto_comment_id, 
        th_comment_id, the_comment_id, pbe_comment_id, nbe_comment_id, ps_comment_id,
        cce_comment_id,vis_comment_id,lcv_comment_id,lmt_comment_id,spt_comment_id,
        lhe_comment_id, li_comment_id, hi_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, rus_fcinv_id, rvs_fcinv_id, rts_fcinv_id, rps_fcinv_id, 
        rtd_fcinv_id, rh_fcinv_id, lcb_fcinv_id, lct_fcinv_id, msl_fcinv_id, lil_fcinv_id, 
        tpw_fcinv_id, r01_fcinv_id, rto_fcinv_id, s01_fcinv_id, sto_fcinv_id, th_fcinv_id, 
        the_fcinv_id, pbe_fcinv_id, nbe_fcinv_id, ps_fcinv_id, cce_fcinv_id, vis_fcinv_id,
        lcv_fcinv_id, lmt_fcinv_id, spt_fcinv_id, lhe_fcinv_id, li_fcinv_id, hi_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  rus_valid_range[2];
   float  rvs_valid_range[2];
   float  rts_valid_range[2];
   float  rps_valid_range[2];
   float  rtd_valid_range[2];
   float  rh_valid_range[2];
   float  lcb_valid_range[2];
   float  lct_valid_range[2];
   float  msl_valid_range[2];
   float  lil_valid_range[2];
   float  tpw_valid_range[2];
   float  r01_valid_range[2];
   float  rto_valid_range[2];
   float  s01_valid_range[2];
   float  sto_valid_range[2];
   float  th_valid_range[2];
   float  the_valid_range[2];
   float  pbe_valid_range[2];
   float  nbe_valid_range[2];
   float  ps_valid_range[2];
   float  cce_valid_range[2];
   float  vis_valid_range[2];
   float  lcv_valid_range[2];
   float  lmt_valid_range[2];
   float  spt_valid_range[2];
   float  lhe_valid_range[2];
   float  li_valid_range[2];
   float  hi_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 28L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rus_id = ncvardef (ncid, "rus", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rvs_id = ncvardef (ncid, "rvs", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rts_id = ncvardef (ncid, "rts", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rps_id = ncvardef (ncid, "rps", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rtd_id = ncvardef (ncid, "rtd", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rh_id = ncvardef (ncid, "rh", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lcb_id = ncvardef (ncid, "lcb", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lct_id = ncvardef (ncid, "lct", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   msl_id = ncvardef (ncid, "msl", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lil_id = ncvardef (ncid, "lil", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   tpw_id = ncvardef (ncid, "tpw", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   r01_id = ncvardef (ncid, "r01", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rto_id = ncvardef (ncid, "rto", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s01_id = ncvardef (ncid, "s01", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sto_id = ncvardef (ncid, "sto", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   th_id = ncvardef (ncid, "th", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   the_id = ncvardef (ncid, "the", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pbe_id = ncvardef (ncid, "pbe", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   nbe_id = ncvardef (ncid, "nbe", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ps_id = ncvardef (ncid, "ps", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   cce_id = ncvardef (ncid, "cce", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vis_id = ncvardef (ncid, "vis", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lcv_id = ncvardef (ncid, "lcv", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lmt_id = ncvardef (ncid, "lmt", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   spt_id = ncvardef (ncid, "spt", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lhe_id = ncvardef (ncid, "lhe", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   li_id = ncvardef (ncid, "li", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   hi_id = ncvardef (ncid, "hi", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rus_comment_id = ncvardef (ncid, "rus_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rvs_comment_id = ncvardef (ncid, "rvs_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rts_comment_id = ncvardef (ncid, "rts_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rps_comment_id = ncvardef (ncid, "rps_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rtd_comment_id = ncvardef (ncid, "rtd_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rh_comment_id = ncvardef (ncid, "rh_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lcb_comment_id = ncvardef (ncid, "lcb_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lct_comment_id = ncvardef (ncid, "lct_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   msl_comment_id = ncvardef (ncid, "msl_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lil_comment_id = ncvardef (ncid, "lil_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   tpw_comment_id = ncvardef (ncid, "tpw_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   r01_comment_id = ncvardef (ncid, "r01_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rto_comment_id = ncvardef (ncid, "rto_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s01_comment_id = ncvardef (ncid, "s01_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sto_comment_id = ncvardef (ncid, "sto_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   th_comment_id = ncvardef (ncid, "th_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   the_comment_id = ncvardef (ncid, "the_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pbe_comment_id = ncvardef (ncid, "pbe_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   nbe_comment_id = ncvardef (ncid, "nbe_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ps_comment_id = ncvardef (ncid, "ps_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   cce_comment_id = ncvardef (ncid, "cce_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vis_comment_id = ncvardef (ncid, "vis_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lcv_comment_id = ncvardef (ncid, "lcv_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lmt_comment_id = ncvardef (ncid, "lmt_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   spt_comment_id = ncvardef (ncid, "spt_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lhe_comment_id = ncvardef (ncid, "lhe_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   li_comment_id = ncvardef (ncid, "li_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   hi_comment_id = ncvardef (ncid, "hi_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rus_fcinv_id = ncvardef (ncid, "rus_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rvs_fcinv_id = ncvardef (ncid, "rvs_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rts_fcinv_id = ncvardef (ncid, "rts_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rps_fcinv_id = ncvardef (ncid, "rps_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rtd_fcinv_id = ncvardef (ncid, "rtd_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rh_fcinv_id = ncvardef (ncid, "rh_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lcb_fcinv_id = ncvardef (ncid, "lcb_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lct_fcinv_id = ncvardef (ncid, "lct_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   msl_fcinv_id = ncvardef (ncid, "msl_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lil_fcinv_id = ncvardef (ncid, "lil_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   tpw_fcinv_id = ncvardef (ncid, "tpw_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   r01_fcinv_id = ncvardef (ncid, "r01_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rto_fcinv_id = ncvardef (ncid, "rto_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s01_fcinv_id = ncvardef (ncid, "s01_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sto_fcinv_id = ncvardef (ncid, "sto_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   th_fcinv_id = ncvardef (ncid, "th_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   the_fcinv_id = ncvardef (ncid, "the_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pbe_fcinv_id = ncvardef (ncid, "pbe_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   nbe_fcinv_id = ncvardef (ncid, "nbe_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ps_fcinv_id = ncvardef (ncid, "ps_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   cce_fcinv_id = ncvardef (ncid, "cce_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vis_fcinv_id = ncvardef (ncid, "vis_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lcv_fcinv_id = ncvardef (ncid, "lcv_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lmt_fcinv_id = ncvardef (ncid, "lmt_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   spt_fcinv_id = ncvardef (ncid, "spt_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lhe_fcinv_id = ncvardef (ncid, "lhe_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   li_fcinv_id = ncvardef (ncid, "li_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   hi_fcinv_id = ncvardef (ncid, "hi_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, rus_id, "long_name", NC_CHAR, 27, (void *)"LAPS Fcst sfc eastward wind");
   ncattput (ncid, rus_id, "units", NC_CHAR, 13, (void *)"meters/second");
   rus_valid_range[0] = -200;
   rus_valid_range[1] = 200;
   ncattput (ncid, rus_id, "valid_range", NC_FLOAT, 2, (void *) rus_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rus_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   char_val = 'U';
   ncattput (ncid, rus_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rus_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rus_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, rvs_id, "long_name", NC_CHAR, 28, (void *)"LAPS Fcst sfc northward wind");
   ncattput (ncid, rvs_id, "units", NC_CHAR, 13, (void *)"meters/second");
   rvs_valid_range[0] = -200;
   rvs_valid_range[1] = 200;
   ncattput (ncid, rvs_id, "valid_range", NC_FLOAT, 2, (void *) rvs_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rvs_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   char_val = 'V';
   ncattput (ncid, rvs_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rvs_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rvs_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, rts_id, "long_name", NC_CHAR, 25, (void *)"LAPS Fcst sfc temperature");
   ncattput (ncid, rts_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   rts_valid_range[0] = 0;
   rts_valid_range[1] = 500;
   ncattput (ncid, rts_id, "valid_range", NC_FLOAT, 2, (void *) rts_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rts_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   char_val = 'T';
   ncattput (ncid, rts_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rts_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   char_val = 'K';
   ncattput (ncid, rts_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rps_id, "long_name", NC_CHAR, 33, (void *)"LAPS Fcst 1500 M Reduced Pressure");
   ncattput (ncid, rps_id, "units", NC_CHAR, 7, (void *)"pascals");
   rps_valid_range[0] = 0;
   rps_valid_range[1] = 200000;
   ncattput (ncid, rps_id, "valid_range", NC_FLOAT, 2, (void *) rps_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rps_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   char_val = 'P';
   ncattput (ncid, rps_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rps_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, rps_id, "LAPS_units", NC_CHAR, 2, (void *)"PA");
   ncattput (ncid, rtd_id, "long_name", NC_CHAR, 38, (void *)"LAPS Fcst surface dewpoint temperature");
   ncattput (ncid, rtd_id, "units", NC_CHAR, 14, (void *)"degrees kelvin");
   rtd_valid_range[0] = 0;
   rtd_valid_range[1] = 500;
   ncattput (ncid, rtd_id, "valid_range", NC_FLOAT, 2, (void *) rtd_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rtd_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rtd_id, "LAPS_var", NC_CHAR, 2, (void *)"TD");
   ncattput (ncid, rtd_id, "lvl_coord", NC_CHAR, 4, (void *)"AGL ");
   char_val = 'K';
   ncattput (ncid, rtd_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rh_id, "long_name", NC_CHAR, 27, (void *)"LAPS Fcst relative humidity");
   ncattput (ncid, rh_id, "units", NC_CHAR, 7, (void *)"percent");
   rh_valid_range[0] = 0;
   rh_valid_range[1] = 100;
   ncattput (ncid, rh_id, "valid_range", NC_FLOAT, 2, (void *) rh_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rh_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rh_id, "LAPS_var", NC_CHAR, 2, (void *)"RH");
   ncattput (ncid, rh_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, rh_id, "LAPS_units", NC_CHAR, 7, (void *)"PERCENT");
   ncattput (ncid, lcb_id, "long_name", NC_CHAR, 20, (void *)"LAPS Fcst cloud base");
   ncattput (ncid, lcb_id, "units", NC_CHAR, 6, (void *)"meters");
   lcb_valid_range[0] = 0;
   lcb_valid_range[1] = 100000;
   ncattput (ncid, lcb_id, "valid_range", NC_FLOAT, 2, (void *) lcb_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lcb_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lcb_id, "LAPS_var", NC_CHAR, 3, (void *)"LCB");
   ncattput (ncid, lcb_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, lcb_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, lct_id, "long_name", NC_CHAR, 19, (void *)"LAPS Fcst cloud top");
   ncattput (ncid, lct_id, "units", NC_CHAR, 6, (void *)"meters");
   lct_valid_range[0] = 0;
   lct_valid_range[1] = 100000;
   ncattput (ncid, lct_id, "valid_range", NC_FLOAT, 2, (void *) lct_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lct_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lct_id, "LAPS_var", NC_CHAR, 3, (void *)"LCT");
   ncattput (ncid, lct_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, lct_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, msl_id, "long_name", NC_CHAR, 22, (void *)"LAPS Fcst MSL pressure");
   ncattput (ncid, msl_id, "units", NC_CHAR, 7, (void *)"pascals");
   msl_valid_range[0] = 0;
   msl_valid_range[1] = 200000;
   ncattput (ncid, msl_id, "valid_range", NC_FLOAT, 2, (void *) msl_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, msl_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, msl_id, "LAPS_var", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, msl_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, msl_id, "LAPS_units", NC_CHAR, 2, (void *)"PA");
   ncattput (ncid, lil_id, "long_name", NC_CHAR, 33, (void *)"LAPS Fcst integrated liquid water");
   ncattput (ncid, lil_id, "units", NC_CHAR, 14, (void *)"kilograms/meter**2");
   lil_valid_range[0] = 0;
   lil_valid_range[1] = 10;
   ncattput (ncid, lil_id, "valid_range", NC_FLOAT, 2, (void *) lil_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lil_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lil_id, "LAPS_var", NC_CHAR, 3, (void *)"LIL");
   ncattput (ncid, lil_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, lil_id, "LAPS_units", NC_CHAR, 6, (void *)"KG/M**2");
   ncattput (ncid, tpw_id, "long_name", NC_CHAR, 45, (void *)"LAPS Fcst integrated total precipitable water");
   ncattput (ncid, tpw_id, "units", NC_CHAR, 6, (void *)"meters");
   tpw_valid_range[0] = 0;
   tpw_valid_range[1] = 0.1;
   ncattput (ncid, tpw_id, "valid_range", NC_FLOAT, 2, (void *) tpw_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, tpw_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, tpw_id, "LAPS_var", NC_CHAR, 3, (void *)"TPW");
   ncattput (ncid, tpw_id, "lvl_coord", NC_CHAR, 4, (void *)"none");
   char_val = 'M';
   ncattput (ncid, tpw_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, r01_id, "long_name", NC_CHAR, 29, (void *)"LAPS Fcst cycle precip.accum.");
   ncattput (ncid, r01_id, "units", NC_CHAR, 6, (void *)"meters");
   r01_valid_range[0] = 0;
   r01_valid_range[1] = 200;
   ncattput (ncid, r01_id, "valid_range", NC_FLOAT, 2, (void *) r01_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, r01_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, r01_id, "LAPS_var", NC_CHAR, 3, (void *)"R01");
   ncattput (ncid, r01_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, r01_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rto_id, "long_name", NC_CHAR, 36, (void *)"LAPS Fcst storm total precip. accum.");
   ncattput (ncid, rto_id, "units", NC_CHAR, 6, (void *)"meters");
   rto_valid_range[0] = 0;
   rto_valid_range[1] = 200;
   ncattput (ncid, rto_id, "valid_range", NC_FLOAT, 2, (void *) rto_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rto_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rto_id, "LAPS_var", NC_CHAR, 3, (void *)"RTO");
   ncattput (ncid, rto_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, rto_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, s01_id, "long_name", NC_CHAR, 27, (void *)"LAPS Fcst cycle snow accum.");
   ncattput (ncid, s01_id, "units", NC_CHAR, 6, (void *)"meters");
   s01_valid_range[0] = 0;
   s01_valid_range[1] = 200;
   ncattput (ncid, s01_id, "valid_range", NC_FLOAT, 2, (void *) s01_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s01_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s01_id, "LAPS_var", NC_CHAR, 3, (void *)"S01");
   ncattput (ncid, s01_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, s01_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, sto_id, "long_name", NC_CHAR, 39, (void *)"LAPS Fcst storm total snow accumulation");
   ncattput (ncid, sto_id, "units", NC_CHAR, 6, (void *)"meters");
   sto_valid_range[0] = 0;
   sto_valid_range[1] = 200;
   ncattput (ncid, sto_id, "valid_range", NC_FLOAT, 2, (void *) sto_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, sto_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, sto_id, "LAPS_var", NC_CHAR, 3, (void *)"STO");
   ncattput (ncid, sto_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, sto_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, th_id, "long_name", NC_CHAR, 31, (void *)"LAPS Fcst potential temperature");
   ncattput (ncid, th_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   th_valid_range[0] = 0;
   th_valid_range[1] = 500;
   ncattput (ncid, th_id, "valid_range", NC_FLOAT, 2, (void *) th_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, th_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, th_id, "LAPS_var", NC_CHAR, 2, (void *)"TH");
   ncattput (ncid, th_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (ncid, th_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, the_id, "long_name", NC_CHAR, 42, (void *)"LAPS Fcst equivalent potential temperature");
   ncattput (ncid, the_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   the_valid_range[0] = 0;
   the_valid_range[1] = 500;
   ncattput (ncid, the_id, "valid_range", NC_FLOAT, 2, (void *) the_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, the_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, the_id, "LAPS_var", NC_CHAR, 3, (void *)"THE");
   ncattput (ncid, the_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (ncid, the_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, pbe_id, "long_name", NC_CHAR, 33, (void *)"LAPS Fcst positive buoyant energy");
   ncattput (ncid, pbe_id, "units", NC_CHAR, 15, (void *)"joules/kilogram");
   pbe_valid_range[0] = -20000;
   pbe_valid_range[1] = 20000;
   ncattput (ncid, pbe_id, "valid_range", NC_FLOAT, 2, (void *) pbe_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, pbe_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, pbe_id, "LAPS_var", NC_CHAR, 3, (void *)"PBE");
   ncattput (ncid, pbe_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, pbe_id, "LAPS_units", NC_CHAR, 4, (void *)"J/KG");
   ncattput (ncid, nbe_id, "long_name", NC_CHAR, 33, (void *)"LAPS Fcst negative buoyant energy");
   ncattput (ncid, nbe_id, "units", NC_CHAR, 15, (void *)"joules/kilogram");
   nbe_valid_range[0] = -20000;
   nbe_valid_range[1] = 20000;
   ncattput (ncid, nbe_id, "valid_range", NC_FLOAT, 2, (void *) nbe_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, nbe_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, nbe_id, "LAPS_var", NC_CHAR, 3, (void *)"NBE");
   ncattput (ncid, nbe_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, nbe_id, "LAPS_units", NC_CHAR, 4, (void *)"J/KG");
   ncattput (ncid, ps_id, "long_name", NC_CHAR, 26, (void *)"LAPS Fcst surface pressure");
   ncattput (ncid, ps_id, "units", NC_CHAR, 7, (void *)"pascals");
   ps_valid_range[0] = 0;
   ps_valid_range[1] = 200000;
   ncattput (ncid, ps_id, "valid_range", NC_FLOAT, 2, (void *) ps_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, ps_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, ps_id, "LAPS_var", NC_CHAR, 2, (void *)"PS");
   ncattput (ncid, ps_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, ps_id, "LAPS_units", NC_CHAR, 2, (void *)"PA");
   ncattput (ncid, cce_id, "long_name", NC_CHAR, 23, (void *)"LAPS Fcst cloud ceiling");
   ncattput (ncid, cce_id, "units", NC_CHAR, 6, (void *)"meters");
   cce_valid_range[0] = 0;
   cce_valid_range[1] = 100000;
   ncattput (ncid, cce_id, "valid_range", NC_FLOAT, 2, (void *) cce_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, cce_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, cce_id, "LAPS_var", NC_CHAR, 3, (void *)"CCE");
   ncattput (ncid, cce_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'M';
   ncattput (ncid, cce_id, "LAPS_units", NC_CHAR, 1, (void *)&char_val);
   ncattput (ncid, vis_id, "long_name", NC_CHAR, 20, (void *)"LAPS Fcst visibility");
   ncattput (ncid, vis_id, "units", NC_CHAR, 6, (void *)"meters");
   vis_valid_range[0] = 0;
   vis_valid_range[1] = 100000;
   ncattput (ncid, vis_id, "valid_range", NC_FLOAT, 2, (void *) vis_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, vis_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, vis_id, "LAPS_var", NC_CHAR, 3, (void *)"VIS");
   ncattput (ncid, vis_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'M';
   ncattput (ncid, vis_id, "LAPS_units", NC_CHAR, 1, (void *)&char_val);

   ncattput (ncid, lcv_id, "long_name", NC_CHAR, 21, (void *)"LAPS Fcst cloud cover");
   ncattput (ncid, lcv_id, "units", NC_CHAR, 4, (void *)"none");
   lcv_valid_range[0] = 0.0;
   lcv_valid_range[1] = 1.00;
   ncattput (ncid, lcv_id, "valid_range", NC_FLOAT, 2, (void *) lcv_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lcv_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lcv_id, "LAPS_var", NC_CHAR, 3, (void *)"LCB");
   ncattput (ncid, lcv_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, lcv_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, lmt_id, "long_name", NC_CHAR, 23, (void *)"LAPS Fcst max echo tops");
   ncattput (ncid, lmt_id, "units", NC_CHAR, 6, (void *)"meters");
   lmt_valid_range[0] = 0;
   lmt_valid_range[1] = 100000;
   ncattput (ncid, lmt_id, "valid_range", NC_FLOAT, 2, (void *) lmt_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lmt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lmt_id, "LAPS_var", NC_CHAR, 3, (void *)"LMT");
   ncattput (ncid, lmt_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, lmt_id, "LAPS_units", NC_CHAR, 1, (void *)&char_val);

   ncattput (ncid, spt_id, "long_name", NC_CHAR, 29, (void *)"LAPS Fcst surface precip type");
   ncattput (ncid, spt_id, "units", NC_CHAR, 4, (void *)"none");
   spt_valid_range[0] = 0.0;
   spt_valid_range[1] = 0.100;
   ncattput (ncid, spt_id, "valid_range", NC_FLOAT, 2, (void *) spt_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, spt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, spt_id, "LAPS_var", NC_CHAR, 3, (void *)"SPT");
   ncattput (ncid, spt_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, spt_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, lhe_id, "long_name", NC_CHAR, 18, (void *)"LAPS Fcst helicity");
   ncattput (ncid, lhe_id, "units", NC_CHAR, 16, (void *)"meters/second**2");
   lhe_valid_range[0] = -20000;
   lhe_valid_range[1] = 20000;
   ncattput (ncid, lhe_id, "valid_range", NC_FLOAT, 2, (void *) lhe_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lhe_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lhe_id, "LAPS_var", NC_CHAR, 3, (void *)"LHE");
   ncattput (ncid, lhe_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, lhe_id, "LAPS_units", NC_CHAR, 6, (void *)"M/S**2");
   ncattput (ncid, li_id, "long_name", NC_CHAR, 22, (void *)"LAPS Fcst lifted index");
   ncattput (ncid, li_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   li_valid_range[0] = -100;
   li_valid_range[1] = 100;
   ncattput (ncid, li_id, "valid_range", NC_FLOAT, 2, (void *) li_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, li_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, li_id, "LAPS_var", NC_CHAR, 2, (void *)"LI");
   ncattput (ncid, li_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (ncid, li_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, hi_id, "long_name", NC_CHAR, 20, (void *)"LAPS Fcst Heat index");
   ncattput (ncid, hi_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   hi_valid_range[0] = 210;
   hi_valid_range[1] = 366;
   ncattput (ncid, hi_id, "valid_range", NC_FLOAT, 2, (void *) hi_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, hi_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, hi_id, "LAPS_var", NC_CHAR, 2, (void *)"HI");
   ncattput (ncid, hi_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, hi_id, "LAPS_units", NC_CHAR, 1, (void *)"K");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, rus_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rvs_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rts_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rps_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rtd_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rh_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lcb_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lct_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, msl_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lil_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, tpw_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, r01_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rto_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, s01_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, sto_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, th_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, the_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, pbe_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, nbe_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, ps_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, cce_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, vis_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lcv_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lmt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, spt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lhe_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, li_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, hi_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store num_variables */
    static long num_variables = {96};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {29};
    static char model[] = {"LAPS Forecast Modeling System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {28};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {
    static long kdim = {28};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/*************************************************************************/
/* file to create netCDF format file with extension RSM */
#ifdef __STDC__
int cre_rsm(char *fname)                /* create fname */
#else
int cre_rsm(fname)
char *fname;
#endif
{                       /* create rsm.cdf */

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  lsm_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, lsm_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, lsm_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  lsm_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 11L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 11L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lsm_id = ncvardef (ncid, "lsm", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lsm_comment_id = ncvardef (ncid, "lsm_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lsm_fcinv_id = ncvardef (ncid, "lsm_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, lsm_id, "long_name", NC_CHAR, 23, (void *)"LAPS Fcst soil moisture");
   ncattput (ncid, lsm_id, "units", NC_CHAR, 11, (void *)"meter/meter");
   lsm_valid_range[0] = 0;
   lsm_valid_range[1] = 0.1;
   ncattput (ncid, lsm_id, "valid_range", NC_FLOAT, 2, (void *) lsm_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lsm_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lsm_id, "LAPS_var", NC_CHAR, 3, (void *)"LSM");
   ncattput (ncid, lsm_id, "lvl_coord", NC_CHAR, 2, (void *)"CM");
   ncattput (ncid, lsm_id, "LAPS_units", NC_CHAR, 3, (void *)"M/M");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, lsm_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {11};
    static short level[] = {-50, -40, -30, -20, -15, -10, -7, -4, -2, -1, 0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {29};
    static char model[] = {"LAPS Forecast Modeling System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {11};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {11};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/*************************************************************************/
/* file to create netCDF format file with extension VRC */
#ifdef __STDC__
int cre_vrc(char *fname)
#else
int cre_vrc(fname) 		/* create fname */
char *fname;
#endif
{   
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim,
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim,
        asc_len_dim;

   /* variable ids */
   int  now_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, now_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, now_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  now_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1);
   namelen_dim = ncdimdef(cdfid, "namelen", 132);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 1);
   var_len_dim = ncdimdef(cdfid, "var_len", 4);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18);

   /* define variables */
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   now_id = ncvardef (cdfid, "now", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   now_comment_id = ncvardef (cdfid, "now_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   now_fcinv_id = ncvardef (cdfid, "now_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, now_id, "long_name", NC_CHAR, 15, (void *)"NOWRAD 2D radar");
   ncattput (cdfid, now_id, "units", NC_CHAR, 3, (void *)"DBZ");
   now_valid_range[0] = -50;
   now_valid_range[1] = 100;
   ncattput (cdfid, now_id, "valid_range", NC_FLOAT, 2, (void *) now_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, now_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, now_id, "LAPS_var", NC_CHAR, 3, (void *)"REF");
   ncattput (cdfid, now_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (cdfid, now_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, now_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
   static int fctimes_start[] = {0};
   static int fctimes_edges[] = {1};
   static short fctimes[] = {0};
   ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
   static int level_start[] = {0};
   static int level_edges[] = {1};
   static short level[] = {0};
   ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
   static int origin_start[] = {0};
   static int origin_edges[] = {49};
   static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
   ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
   static int model_start[] = {0};
   static int model_edges[] = {43};
   static char model[] = {"LAPS - Local Analysis and Prediction System"};
   ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }
   
   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {1};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {1};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LM1 */
#ifdef __STDC__
int cre_lm1(char *fname)
#else
int cre_lm1(fname) 		/* create fname */
char *fname;
#endif
{   
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  lsm_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, lsm_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, lsm_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  lsm_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 3);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1);
   namelen_dim = ncdimdef(cdfid, "namelen", 132);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 3);
   var_len_dim = ncdimdef(cdfid, "var_len", 4);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18);

   /* define variables */
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lsm_id = ncvardef (cdfid, "lsm", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lsm_comment_id = ncvardef (cdfid, "lsm_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lsm_fcinv_id = ncvardef (cdfid, "lsm_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, lsm_id, "long_name", NC_CHAR, 13, (void *)"soil moisture");
   ncattput (cdfid, lsm_id, "units", NC_CHAR, 18, (void *)"meters**3/meter**3");
   lsm_valid_range[0] = 0;
   lsm_valid_range[1] = 0.1;
   ncattput (cdfid, lsm_id, "valid_range", NC_FLOAT, 2, (void *) lsm_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, lsm_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, lsm_id, "LAPS_var", NC_CHAR, 3, (void *)"LSM");
   ncattput (cdfid, lsm_id, "lvl_coord", NC_CHAR, 3, (void *)"   ");
   ncattput (cdfid, lsm_id, "LAPS_units", NC_CHAR, 3, (void *)"M**3/M**3");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, lsm_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
   static int fctimes_start[] = {0};
   static int fctimes_edges[] = {1};
   static short fctimes[] = {0};
   ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {15};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store level */
   static int level_start[] = {0};
   static int level_edges[] = {3};
   static short level[] = {-1, -2, -3};
   ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
   static int origin_start[] = {0};
   static int origin_edges[] = {49};
   static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
   ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
   static int model_start[] = {0};
   static int model_edges[] = {43};
   static char model[] = {"LAPS - Local Analysis and Prediction System"};
   ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {3};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {3};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LM2 */
#ifdef __STDC__
int cre_lm2(char *fname)
#else
int cre_lm2(fname) 		/* create fname */
char *fname;
#endif
{   
   int  cdfid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  civ_id, dwf_id, wx_id, evp_id, sc_id, sm_id, mwf_id, lvl_id, 
        imax_id, jmax_id, kmax_id, kdim_id, civ_comment_id, dwf_comment_id, 
        wx_comment_id, evp_comment_id, sc_comment_id, sm_comment_id, 
        mwf_comment_id, laps_domain_file_id, asctime_id, fctimes_id, 
        level_id, civ_fcinv_id, dwf_fcinv_id, wx_fcinv_id, evp_fcinv_id, 
        sc_fcinv_id, sm_fcinv_id, mwf_fcinv_id, origin_id, model_id, 
        version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  civ_valid_range[2];
   float  dwf_valid_range[2];
   float  wx_valid_range[2];
   float  evp_valid_range[2];
   float  sc_valid_range[2];
   float  sm_valid_range[2];
   float  mwf_valid_range[2];

   /* enter define mode */
   cdfid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (cdfid == -1) return cdfid;

   /* define dimensions */
   level_dim = ncdimdef(cdfid, "level", 1);
   fctimes_dim = ncdimdef(cdfid, "fctimes", 1);
   namelen_dim = ncdimdef(cdfid, "namelen", 132);
   lat_dim = ncdimdef(cdfid, "lat", NY_LONG);
   lon_dim = ncdimdef(cdfid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(cdfid, "n_2d_grids", 7);
   var_len_dim = ncdimdef(cdfid, "var_len", 4);
   coord_len_dim = ncdimdef(cdfid, "coord_len", 5);
   unit_len_dim = ncdimdef(cdfid, "unit_len", 11);
   comm_len_dim = ncdimdef(cdfid, "comm_len", 126);
   domain_len_dim = ncdimdef(cdfid, "domain_len", 12);
   asc_len_dim = ncdimdef(cdfid, "asc_len", 18);

   /* define variables */
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   civ_id = ncvardef (cdfid, "civ", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   dwf_id = ncvardef (cdfid, "dwf", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   wx_id = ncvardef (cdfid, "wx", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   evp_id = ncvardef (cdfid, "evp", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sc_id = ncvardef (cdfid, "sc", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sm_id = ncvardef (cdfid, "sm", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mwf_id = ncvardef (cdfid, "mwf", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (cdfid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (cdfid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (cdfid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (cdfid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (cdfid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   civ_comment_id = ncvardef (cdfid, "civ_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   dwf_comment_id = ncvardef (cdfid, "dwf_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   wx_comment_id = ncvardef (cdfid, "wx_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   evp_comment_id = ncvardef (cdfid, "evp_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sc_comment_id = ncvardef (cdfid, "sc_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sm_comment_id = ncvardef (cdfid, "sm_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mwf_comment_id = ncvardef (cdfid, "mwf_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (cdfid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (cdfid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (cdfid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (cdfid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   civ_fcinv_id = ncvardef (cdfid, "civ_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dwf_fcinv_id = ncvardef (cdfid, "dwf_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   wx_fcinv_id = ncvardef (cdfid, "wx_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   evp_fcinv_id = ncvardef (cdfid, "evp_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sc_fcinv_id = ncvardef (cdfid, "sc_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sm_fcinv_id = ncvardef (cdfid, "sm_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mwf_fcinv_id = ncvardef (cdfid, "mwf_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (cdfid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (cdfid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (cdfid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (cdfid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (cdfid, civ_id, "long_name", NC_CHAR, 30, (void *)"cumulative infiltration volume");
   ncattput (cdfid, civ_id, "units", NC_CHAR, 6, (void *)"meters");
   civ_valid_range[0] = -200;
   civ_valid_range[1] = 200;
   ncattput (cdfid, civ_id, "valid_range", NC_FLOAT, 2, (void *) civ_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, civ_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, civ_id, "LAPS_var", NC_CHAR, 3, (void *)"CIV");
   ncattput (cdfid, civ_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = '\115';
   ncattput (cdfid, civ_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, dwf_id, "long_name", NC_CHAR, 22, (void *)"depth to wetting front");
   ncattput (cdfid, dwf_id, "units", NC_CHAR, 6, (void *)"meters");
   dwf_valid_range[0] = -200;
   dwf_valid_range[1] = 200;
   ncattput (cdfid, dwf_id, "valid_range", NC_FLOAT, 2, (void *) dwf_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, dwf_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, dwf_id, "LAPS_var", NC_CHAR, 3, (void *)"DWF");
   ncattput (cdfid, dwf_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = '\115';
   ncattput (cdfid, dwf_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (cdfid, wx_id, "long_name", NC_CHAR, 18, (void *)"wet/dry grid point");
   ncattput (cdfid, wx_id, "units", NC_CHAR, 4, (void *)"none");
   wx_valid_range[0] = 0;
   wx_valid_range[1] = 1;
   ncattput (cdfid, wx_id, "valid_range", NC_FLOAT, 2, (void *) wx_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, wx_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, wx_id, "LAPS_var", NC_CHAR, 2, (void *)"WX");
   ncattput (cdfid, wx_id, "lvl_coord", NC_CHAR, 4, (void *)"none");
   ncattput (cdfid, wx_id, "LAPS_units", NC_CHAR, 4, (void *)"none");
   ncattput (cdfid, evp_id, "long_name", NC_CHAR, 16, (void *)"evaporation dataÀ");
   ncattput (cdfid, evp_id, "units", NC_CHAR, 13, (void *)"meters/second");
   evp_valid_range[0] = -75;
   evp_valid_range[1] = 125;
   ncattput (cdfid, evp_id, "valid_range", NC_FLOAT, 2, (void *) evp_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, evp_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, evp_id, "LAPS_var", NC_CHAR, 3, (void *)"EVP");
   ncattput (cdfid, evp_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (cdfid, evp_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (cdfid, sc_id, "long_name", NC_CHAR, 12, (void *)"snow covered");
   ncattput (cdfid, sc_id, "units", NC_CHAR, 4, (void *)"none");
   sc_valid_range[0] = 0;
   sc_valid_range[1] = 1;
   ncattput (cdfid, sc_id, "valid_range", NC_FLOAT, 2, (void *) sc_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, sc_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, sc_id, "LAPS_var", NC_CHAR, 2, (void *)"SC");
   ncattput (cdfid, sc_id, "lvl_coord", NC_CHAR, 5, (void *)"none ");
   ncattput (cdfid, sc_id, "LAPS_units", NC_CHAR, 4, (void *)"none");
   ncattput (cdfid, sm_id, "long_name", NC_CHAR, 12, (void *)"snow melting");
   ncattput (cdfid, sm_id, "units", NC_CHAR, 19, (void *)"meters**3/meters**3");
   sm_valid_range[0] = 0;
   sm_valid_range[1] = 20000;
   ncattput (cdfid, sm_id, "valid_range", NC_FLOAT, 2, (void *) sm_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, sm_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, sm_id, "LAPS_var", NC_CHAR, 2, (void *)"SM");
   ncattput (cdfid, sm_id, "lvl_coord", NC_CHAR, 4, (void *)"none");
   ncattput (cdfid, sm_id, "LAPS_units", NC_CHAR, 9, (void *)"M**3/M**3");
   ncattput (cdfid, mwf_id, "long_name", NC_CHAR, 35, (void *)"soil moisture content wetting front");
   ncattput (cdfid, mwf_id, "units", NC_CHAR, 19, (void *)"meters**3/meters**3");
   mwf_valid_range[0] = -20000;
   mwf_valid_range[1] = 20000;
   ncattput (cdfid, mwf_id, "valid_range", NC_FLOAT, 2, (void *) mwf_valid_range);
   float_val = 1.0e+37;
   ncattput (cdfid, mwf_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (cdfid, mwf_id, "LAPS_var", NC_CHAR, 3, (void *)"MWF");
   ncattput (cdfid, mwf_id, "lvl_coord", NC_CHAR, 4, (void *)"none");
   ncattput (cdfid, mwf_id, "LAPS_units", NC_CHAR, 9, (void *)"M**3/M**3");
   ncattput (cdfid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (cdfid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (cdfid, civ_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, dwf_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, wx_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, evp_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, sc_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, sm_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (cdfid, mwf_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (cdfid);

   {			/* store version */
    static long version = {2};
    ncvarput1(cdfid, version_id, (int *)0, (void *)&version);
   }

   {			/* store fctimes */
   static int fctimes_start[] = {0};
   static int fctimes_edges[] = {1};
   static short fctimes[] = {0};
   ncvarput(cdfid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {33};
    ncvarput1(cdfid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {			/* store level */
   static int level_start[] = {0};
   static int level_edges[] = {1};
   static short level[] = {0};
   ncvarput(cdfid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
   static int origin_start[] = {0};
   static int origin_edges[] = {49};
   static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
   ncvarput(cdfid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
   static int model_start[] = {0};
   static int model_edges[] = {43};
   static char model[] = {"LAPS - Local Analysis and Prediction System"};
   ncvarput(cdfid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(cdfid, imax_id, (int *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(cdfid, jmax_id, (int *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {7};
    ncvarput1(cdfid, kmax_id, (int *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {7};
    ncvarput1(cdfid, kdim_id, (int *)0, (void *)&kdim);
   }
   return cdfid;
}
/*************************************************************************/
/* file to create netCDF format file with extension VRD */
#ifdef __STDC__
int cre_vrd(char *fname)
#else
int cre_vrd(fname)              /* create fname */
char *fname;
#endif
{
   int  ncid;                   /* netCDF id */
 
   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim,
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim,
        comm_len_dim, domain_len_dim, asc_len_dim;
 
   /* variable ids */
   int  refd_id, veld_id, lvl_id, imax_id, jmax_id, kmax_id,
        kdim_id, veld_comment_id, refd_comment_id, laps_domain_file_id,
        asctime_id, fctimes_id, level_id, refd_fcinv_id,
        veld_fcinv_id, origin_id, model_id, version_id,num_variables_id;
 
   /* variable shapes */
   int dims[4];
 
   /* containers for scalar attributes */
   short  short_val;
   float  float_val;
 
   /* attribute vectors */
   float  refd_valid_range[2];
   float  veld_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;
 
   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", NZ_LONG);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 42L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);
 
   /* define variables */
 
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   refd_id = ncvardef (ncid, "refd", NC_FLOAT, 4, dims);
 
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   veld_id = ncvardef (ncid, "veld", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);
 
   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);
 
   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);
 
   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);
 
   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);
 
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   veld_comment_id = ncvardef (ncid, "veld_comment", NC_CHAR, 3, dims);
 
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   refd_comment_id = ncvardef (ncid, "refd_comment", NC_CHAR, 3, dims);
 
   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);
 
   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);
 
   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);
 
   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);
 
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   refd_fcinv_id = ncvardef (ncid, "refd_fcinv", NC_SHORT, 2, dims);
 
   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   veld_fcinv_id = ncvardef (ncid, "veld_fcinv", NC_SHORT, 2, dims);
 
   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);
 
   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);
 
   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);
 
   /* assign attributes */
   ncattput (ncid, refd_id, "long_name", NC_CHAR, 8, (void *)"3D radar");
   ncattput (ncid, refd_id, "units", NC_CHAR, 3, (void *)"DBZ");
   refd_valid_range[0] = -50;
   refd_valid_range[1] = 100;
   ncattput (ncid, refd_id, "valid_range", NC_FLOAT, 2, (void *) refd_valid_range
);
   float_val = 1e+37;
   ncattput (ncid, refd_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, refd_id, "LAPS_var", NC_CHAR, 3, (void *)"REF");
   ncattput (ncid, refd_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, refd_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (ncid, veld_id, "long_name", NC_CHAR, 8, (void *)"3D radar");
   ncattput (ncid, veld_id, "units", NC_CHAR, 13, (void *)"meters/second");
   veld_valid_range[0] = -50;
   veld_valid_range[1] = 100;
   ncattput (ncid, veld_id, "valid_range", NC_FLOAT, 2, (void *) veld_valid_range
);
   float_val = 1e+37;
   ncattput (ncid, veld_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, veld_id, "LAPS_var", NC_CHAR, 3, (void *)"VEL");
   ncattput (ncid, veld_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, veld_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times")
;
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, refd_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, veld_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
 
   /* leave define mode */
   ncendef (ncid);
 
   {                    /* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }
 
   {			/* store num_variables */
    static long num_variables = {18};
    ncvarput1(ncid, num_variables_id, (int *)0, (void *)&num_variables);
   }

   {                    /* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }
 
   {                    /* store level */
    static long level_start[] = {0};
    static long level_edges[] = {NZ};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600,
 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }
 
   {                    /* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }
 
   {                    /* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }
 
   {                    /* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }
 
   {                    /* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }
 
   {                    /* store kmax */
    static long kmax = {(NZ*2)};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }
 

   {                    /* store kdim */
    static long kdim = {(NZ*2)};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid ;
}
/*************************************************************************/
/* file to create netCDF format radar file */
#ifdef __STDC__
int cre_v_radar(char *fname)
#else
int cre_v_radar(fname)              /* create fname */
char *fname;
#endif
{

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  refd_id, veld_id, nyqd_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        veld_comment_id, refd_comment_id, nyqd_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, refd_fcinv_id, veld_fcinv_id, 
        nyqd_fcinv_id, origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  refd_valid_range[2];
   float  veld_valid_range[2];
   float  nyqd_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;


   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", NZ_LONG);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 63L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   refd_id = ncvardef (ncid, "refd", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   veld_id = ncvardef (ncid, "veld", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim; 
   dims[2] = lat_dim; 
   dims[3] = lon_dim;
   nyqd_id = ncvardef (ncid, "nyqd", NC_FLOAT, 4, dims);
 
   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   veld_comment_id = ncvardef (ncid, "veld_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   refd_comment_id = ncvardef (ncid, "refd_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   nyqd_comment_id = ncvardef (ncid, "nyqd_comment", NC_CHAR, 3, dims);
 
   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   refd_fcinv_id = ncvardef (ncid, "refd_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   veld_fcinv_id = ncvardef (ncid, "veld_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   nyqd_fcinv_id = ncvardef (ncid, "nyqd_fcinv", NC_SHORT, 2, dims);
 
   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, refd_id, "long_name", NC_CHAR, 8, (void *)"3D radar");
   ncattput (ncid, refd_id, "units", NC_CHAR, 3, (void *)"DBZ");
   refd_valid_range[0] = -50;
   refd_valid_range[1] = 100;
   ncattput (ncid, refd_id, "valid_range", NC_FLOAT, 2, (void *) refd_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, refd_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, refd_id, "LAPS_var", NC_CHAR, 3, (void *)"REF");
   ncattput (ncid, refd_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, refd_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (ncid, veld_id, "long_name", NC_CHAR, 8, (void *)"3D radar");
   ncattput (ncid, veld_id, "units", NC_CHAR, 13, (void *)"meters/second");
   veld_valid_range[0] = -50;
   veld_valid_range[1] = 100;
   ncattput (ncid, veld_id, "valid_range", NC_FLOAT, 2, (void *) veld_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, veld_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, veld_id, "LAPS_var", NC_CHAR, 3, (void *)"VEL");
   ncattput (ncid, veld_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, veld_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, nyqd_id, "long_name", NC_CHAR, 25, (void *)"3D radar nyquist velocity");
   ncattput (ncid, nyqd_id, "units", NC_CHAR, 13, (void *)"meters/second");
   nyqd_valid_range[0] = -50;
   nyqd_valid_range[1] = 100;
   ncattput (ncid, nyqd_id, "valid_range", NC_FLOAT, 2, (void *) nyqd_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, nyqd_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, nyqd_id, "LAPS_var", NC_CHAR, 3, (void *)"NYQ");
   ncattput (ncid, nyqd_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, nyqd_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, refd_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, veld_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, nyqd_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {21};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {NZ};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {(NZ*3)};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {(NZ*3)};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/*************************************************************************/
/* file to create netCDF format file with extension LGA */
#ifdef __STDC__
int cre_lga(char *fname)
#else
int cre_lga(fname)              /* create fname */
char *fname;
#endif
{

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  mgz_id, mgt_id, mgrh_id, mgu_id, mgv_id, lvl_id, imax_id, jmax_id, 
        kmax_id, kdim_id, mgz_comment_id, mgt_comment_id, mgrh_comment_id, 
        mgu_comment_id, mgv_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, mgz_fcinv_id, mgt_fcinv_id, mgrh_fcinv_id, 
        mgu_fcinv_id, mgv_fcinv_id, origin_id, model_id, version_id, 
        num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short short_val;
   float float_val;

   /* attribute vectors */
   float  mgz_valid_range[2];
   float  mgt_valid_range[2];
   float  mgrh_valid_range[2];
   float  mgu_valid_range[2];
   float  mgv_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 21L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 105L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mgz_id = ncvardef (ncid, "mgz", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mgt_id = ncvardef (ncid, "mgt", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mgrh_id = ncvardef (ncid, "mgrh", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mgu_id = ncvardef (ncid, "mgu", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mgv_id = ncvardef (ncid, "mgv", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mgz_comment_id = ncvardef (ncid, "mgz_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mgt_comment_id = ncvardef (ncid, "mgt_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mgrh_comment_id = ncvardef (ncid, "mgrh_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mgu_comment_id = ncvardef (ncid, "mgu_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mgv_comment_id = ncvardef (ncid, "mgv_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mgz_fcinv_id = ncvardef (ncid, "mgz_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mgt_fcinv_id = ncvardef (ncid, "mgt_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mgrh_fcinv_id = ncvardef (ncid, "mgrh_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mgu_fcinv_id = ncvardef (ncid, "mgu_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mgv_fcinv_id = ncvardef (ncid, "mgv_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, mgz_id, "long_name", NC_CHAR, 12, (void *)"MAPS heights");
   ncattput (ncid, mgz_id, "units", NC_CHAR, 6, (void *)"meters");
   mgz_valid_range[0] = -200;
   mgz_valid_range[1] = 200;
   ncattput (ncid, mgz_id, "valid_range", NC_FLOAT, 2, (void *) mgz_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mgz_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mgz_id, "LAPS_var", NC_CHAR, 2, (void *)"HT");
   ncattput (ncid, mgz_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mgz_id, "LAPS_units", NC_CHAR, 6, (void *)"METERS");
   ncattput (ncid, mgt_id, "long_name", NC_CHAR, 16, (void *)"MAPS temperature");
   ncattput (ncid, mgt_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   mgt_valid_range[0] = -200;
   mgt_valid_range[1] = 200;
   ncattput (ncid, mgt_id, "valid_range", NC_FLOAT, 2, (void *) mgt_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mgt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mgt_id, "LAPS_var", NC_CHAR, 2, (void *)"T3");
   ncattput (ncid, mgt_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mgt_id, "LAPS_units", NC_CHAR, 6, (void *)"KELVIN");
   ncattput (ncid, mgrh_id, "long_name", NC_CHAR, 22, (void *)"MAPS specific humidity");
   ncattput (ncid, mgrh_id, "units", NC_CHAR, 4, (void *)"none");
   mgrh_valid_range[0] = -200;
   mgrh_valid_range[1] = 200;
   ncattput (ncid, mgrh_id, "valid_range", NC_FLOAT, 2, (void *) mgrh_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mgrh_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mgrh_id, "LAPS_var", NC_CHAR, 2, (void *)"SH");
   ncattput (ncid, mgrh_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mgrh_id, "LAPS_units", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, mgu_id, "long_name", NC_CHAR, 21, (void *)"MAPS u wind component");
   ncattput (ncid, mgu_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mgu_valid_range[0] = -200;
   mgu_valid_range[1] = 200;
   ncattput (ncid, mgu_id, "valid_range", NC_FLOAT, 2, (void *) mgu_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mgu_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mgu_id, "LAPS_var", NC_CHAR, 2, (void *)"U3");
   ncattput (ncid, mgu_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mgu_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, mgv_id, "long_name", NC_CHAR, 21, (void *)"MAPS v wind component");
   ncattput (ncid, mgv_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mgv_valid_range[0] = -20000;
   mgv_valid_range[1] = 20000;
   ncattput (ncid, mgv_id, "valid_range", NC_FLOAT, 2, (void *) mgv_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mgv_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mgv_id, "LAPS_var", NC_CHAR, 2, (void *)"V3");
   ncattput (ncid, mgv_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mgv_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, mgz_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mgt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mgrh_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mgu_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mgv_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {27};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {105};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {105};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return (ncid);
}
/*************************************************************************/
/* file to create netCDF format file with extension LGF */
#ifdef __STDC__
int cre_lgf(char *fname)
#else
int cre_lgf(fname)              /* create fname */
char *fname;
#endif
{

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  mfz_id, mft_id, mfrh_id, mfu_id, mfv_id, lvl_id, imax_id, jmax_id, 
        kmax_id, kdim_id, mfz_comment_id, mft_comment_id, mfrh_comment_id, 
        mfu_comment_id, mfv_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, mfz_fcinv_id, mft_fcinv_id, mfrh_fcinv_id, 
        mfu_fcinv_id, mfv_fcinv_id, origin_id, model_id, version_id, 
        num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short short_val;
   float float_val;

   /* attribute vectors */
   float  mfz_valid_range[2];
   float  mft_valid_range[2];
   float  mfrh_valid_range[2];
   float  mfu_valid_range[2];
   float  mfv_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 21L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 105L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mfz_id = ncvardef (ncid, "mfz", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mft_id = ncvardef (ncid, "mft", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mfrh_id = ncvardef (ncid, "mfrh", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mfu_id = ncvardef (ncid, "mfu", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   mfv_id = ncvardef (ncid, "mfv", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mfz_comment_id = ncvardef (ncid, "mfz_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mft_comment_id = ncvardef (ncid, "mft_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mfrh_comment_id = ncvardef (ncid, "mfrh_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mfu_comment_id = ncvardef (ncid, "mfu_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   mfv_comment_id = ncvardef (ncid, "mfv_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mfz_fcinv_id = ncvardef (ncid, "mfz_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mft_fcinv_id = ncvardef (ncid, "mft_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mfrh_fcinv_id = ncvardef (ncid, "mfrh_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mfu_fcinv_id = ncvardef (ncid, "mfu_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   mfv_fcinv_id = ncvardef (ncid, "mfv_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, mfz_id, "long_name", NC_CHAR, 12, (void *)"MAPS heights");
   ncattput (ncid, mfz_id, "units", NC_CHAR, 6, (void *)"meters");
   mfz_valid_range[0] = -200;
   mfz_valid_range[1] = 200;
   ncattput (ncid, mfz_id, "valid_range", NC_FLOAT, 2, (void *) mfz_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mfz_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mfz_id, "LAPS_var", NC_CHAR, 2, (void *)"HT");
   ncattput (ncid, mfz_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mfz_id, "LAPS_units", NC_CHAR, 6, (void *)"METERS");
   ncattput (ncid, mft_id, "long_name", NC_CHAR, 16, (void *)"MAPS temperature");
   ncattput (ncid, mft_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   mft_valid_range[0] = -200;
   mft_valid_range[1] = 200;
   ncattput (ncid, mft_id, "valid_range", NC_FLOAT, 2, (void *) mft_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mft_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mft_id, "LAPS_var", NC_CHAR, 2, (void *)"T3");
   ncattput (ncid, mft_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mft_id, "LAPS_units", NC_CHAR, 6, (void *)"KELVIN");
   ncattput (ncid, mfrh_id, "long_name", NC_CHAR, 22, (void *)"MAPS specific humidity");
   ncattput (ncid, mfrh_id, "units", NC_CHAR, 4, (void *)"none");
   mfrh_valid_range[0] = -200;
   mfrh_valid_range[1] = 200;
   ncattput (ncid, mfrh_id, "valid_range", NC_FLOAT, 2, (void *) mfrh_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mfrh_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mfrh_id, "LAPS_var", NC_CHAR, 2, (void *)"SH");
   ncattput (ncid, mfrh_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mfrh_id, "LAPS_units", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, mfu_id, "long_name", NC_CHAR, 21, (void *)"MAPS u wind component");
   ncattput (ncid, mfu_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mfu_valid_range[0] = -200;
   mfu_valid_range[1] = 200;
   ncattput (ncid, mfu_id, "valid_range", NC_FLOAT, 2, (void *) mfu_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mfu_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mfu_id, "LAPS_var", NC_CHAR, 2, (void *)"U3");
   ncattput (ncid, mfu_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mfu_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, mfv_id, "long_name", NC_CHAR, 21, (void *)"MAPS v wind component");
   ncattput (ncid, mfv_id, "units", NC_CHAR, 13, (void *)"meters/second");
   mfv_valid_range[0] = -20000;
   mfv_valid_range[1] = 20000;
   ncattput (ncid, mfv_id, "valid_range", NC_FLOAT, 2, (void *) mfv_valid_range);
   float_val = 1e+37;
   ncattput (ncid, mfv_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, mfv_id, "LAPS_var", NC_CHAR, 2, (void *)"V3");
   ncattput (ncid, mfv_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, mfv_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, mfz_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mft_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mfrh_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mfu_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, mfv_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {27};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {43};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {105};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {105};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return (ncid);
}
/*************************************************************************/
/* file to create netCDF format file with extension LN3 */
#ifdef __STDC__
int cre_ln3(char *fname)
#else
int cre_ln3(fname)              /* create fname */
char *fname;
#endif
{

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  r04_id, r48_id, r8c_id, et_id, rco_id, vil_id, lvl_id, imax_id, 
        jmax_id, kmax_id, kdim_id, r04_comment_id, r48_comment_id, 
        r8c_comment_id, et_comment_id, rco_comment_id, vil_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, r04_fcinv_id, 
        r48_fcinv_id, r8c_fcinv_id, et_fcinv_id, rco_fcinv_id, vil_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  r04_valid_range[2];
   float  r48_valid_range[2];
   float  r8c_valid_range[2];
   float  et_valid_range[2];
   float  rco_valid_range[2];
   float  vil_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 6L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   r04_id = ncvardef (ncid, "r04", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   r48_id = ncvardef (ncid, "r48", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   r8c_id = ncvardef (ncid, "r8c", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   et_id = ncvardef (ncid, "et", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rco_id = ncvardef (ncid, "rco", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vil_id = ncvardef (ncid, "vil", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   r04_comment_id = ncvardef (ncid, "r04_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   r48_comment_id = ncvardef (ncid, "r48_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   r8c_comment_id = ncvardef (ncid, "r8c_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   et_comment_id = ncvardef (ncid, "et_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rco_comment_id = ncvardef (ncid, "rco_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vil_comment_id = ncvardef (ncid, "vil_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   r04_fcinv_id = ncvardef (ncid, "r04_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   r48_fcinv_id = ncvardef (ncid, "r48_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   r8c_fcinv_id = ncvardef (ncid, "r8c_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   et_fcinv_id = ncvardef (ncid, "et_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rco_fcinv_id = ncvardef (ncid, "rco_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vil_fcinv_id = ncvardef (ncid, "vil_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, r04_id, "long_name", NC_CHAR, 38, (void *)"low level layer composite reflectivity");
   ncattput (ncid, r04_id, "units", NC_CHAR, 3, (void *)"dBZ");
   r04_valid_range[0] = -10;
   r04_valid_range[1] = 80;
   ncattput (ncid, r04_id, "valid_range", NC_FLOAT, 2, (void *) r04_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, r04_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, r04_id, "LAPS_var", NC_CHAR, 3, (void *)"R04");
   ncattput (ncid, r04_id, "lvl_coord", NC_CHAR, 4, (void *)"    ");
   ncattput (ncid, r04_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (ncid, r48_id, "long_name", NC_CHAR, 38, (void *)"mid level layer composite reflectivity");
   ncattput (ncid, r48_id, "units", NC_CHAR, 3, (void *)"dBZ");
   r48_valid_range[0] = -10;
   r48_valid_range[1] = 80;
   ncattput (ncid, r48_id, "valid_range", NC_FLOAT, 2, (void *) r48_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, r48_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, r48_id, "LAPS_var", NC_CHAR, 3, (void *)"R48");
   ncattput (ncid, r48_id, "lvl_coord", NC_CHAR, 4, (void *)"    ");
   ncattput (ncid, r48_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (ncid, r8c_id, "long_name", NC_CHAR, 39, (void *)"high level layer composite reflectivity");
   ncattput (ncid, r8c_id, "units", NC_CHAR, 3, (void *)"dBZ");
   r8c_valid_range[0] = -10;
   r8c_valid_range[1] = 80;
   ncattput (ncid, r8c_id, "valid_range", NC_FLOAT, 2, (void *) r8c_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, r8c_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, r8c_id, "LAPS_var", NC_CHAR, 3, (void *)"R8C");
   ncattput (ncid, r8c_id, "lvl_coord", NC_CHAR, 4, (void *)"    ");
   ncattput (ncid, r8c_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (ncid, et_id, "long_name", NC_CHAR, 16, (void *)"echo tops height");
   ncattput (ncid, et_id, "units", NC_CHAR, 6, (void *)"meters");
   et_valid_range[0] = 0;
   et_valid_range[1] = 30000;
   ncattput (ncid, et_id, "valid_range", NC_FLOAT, 2, (void *) et_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, et_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, et_id, "LAPS_var", NC_CHAR, 2, (void *)"ET");
   ncattput (ncid, et_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, et_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rco_id, "long_name", NC_CHAR, 22, (void *)"composite reflectivity");
   ncattput (ncid, rco_id, "units", NC_CHAR, 3, (void *)"dBZ");
   rco_valid_range[0] = -10;
   rco_valid_range[1] = 80;
   ncattput (ncid, rco_id, "valid_range", NC_FLOAT, 2, (void *) rco_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rco_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rco_id, "LAPS_var", NC_CHAR, 3, (void *)"RCO");
   ncattput (ncid, rco_id, "lvl_coord", NC_CHAR, 4, (void *)"    ");
   ncattput (ncid, rco_id, "LAPS_units", NC_CHAR, 3, (void *)"DBZ");
   ncattput (ncid, vil_id, "long_name", NC_CHAR, 34, (void *)"radar vertically integrated liquid");
   ncattput (ncid, vil_id, "units", NC_CHAR, 14, (void *)"grams/meter**2");
   vil_valid_range[0] = 0;
   vil_valid_range[1] = 0.1;
   ncattput (ncid, vil_id, "valid_range", NC_FLOAT, 2, (void *) vil_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, vil_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, vil_id, "LAPS_var", NC_CHAR, 3, (void *)"VIL");
   ncattput (ncid, vil_id, "lvl_coord", NC_CHAR, 4, (void *)"    ");
   ncattput (ncid, vil_id, "LAPS_units", NC_CHAR, 6, (void *)"G/M**2");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, r04_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, r48_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, r8c_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, et_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rco_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, vil_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {30};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {132};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {132};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {6};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {6};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/*************************************************************************/
/* file to create netCDF format file with extension MM5 */
#ifdef __STDC__
int cre_mm5(char *fname)                /* create fname */
#else
int cre_mm5(fname)
char *fname;
#endif
{

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  rz_id, ru_id, rv_id, rw_id, rt_id, rsh_id, rh3_id, lwc_id, ice_id, 
        lvl_id, imax_id, jmax_id, kmax_id, kdim_id, rz_comment_id, 
        ru_comment_id, rv_comment_id, rw_comment_id, rt_comment_id, 
        rsh_comment_id, rh3_comment_id, lwc_comment_id, ice_comment_id, 
        laps_domain_file_id, asctime_id, fctimes_id, level_id, rz_fcinv_id, 
        ru_fcinv_id, rv_fcinv_id, rw_fcinv_id, rt_fcinv_id, rsh_fcinv_id, 
        rh3_fcinv_id, lwc_fcinv_id, ice_fcinv_id, origin_id, model_id, 
        version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  rz_valid_range[2];
   float  ru_valid_range[2];
   float  rv_valid_range[2];
   float  rw_valid_range[2];
   float  rt_valid_range[2];
   float  rsh_valid_range[2];
   float  rh3_valid_range[2];
   float  lwc_valid_range[2];
   float  ice_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 21L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 189L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rz_id = ncvardef (ncid, "rz", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ru_id = ncvardef (ncid, "ru", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rv_id = ncvardef (ncid, "rv", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rw_id = ncvardef (ncid, "rw", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rt_id = ncvardef (ncid, "rt", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rsh_id = ncvardef (ncid, "rsh", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rh3_id = ncvardef (ncid, "rh3", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lwc_id = ncvardef (ncid, "lwc", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ice_id = ncvardef (ncid, "ice", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rz_comment_id = ncvardef (ncid, "rz_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ru_comment_id = ncvardef (ncid, "ru_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rv_comment_id = ncvardef (ncid, "rv_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rw_comment_id = ncvardef (ncid, "rw_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rt_comment_id = ncvardef (ncid, "rt_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rsh_comment_id = ncvardef (ncid, "rsh_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rh3_comment_id = ncvardef (ncid, "rh3_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lwc_comment_id = ncvardef (ncid, "lwc_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ice_comment_id = ncvardef (ncid, "ice_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rz_fcinv_id = ncvardef (ncid, "rz_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ru_fcinv_id = ncvardef (ncid, "ru_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rv_fcinv_id = ncvardef (ncid, "rv_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rw_fcinv_id = ncvardef (ncid, "rw_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rt_fcinv_id = ncvardef (ncid, "rt_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rsh_fcinv_id = ncvardef (ncid, "rsh_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rh3_fcinv_id = ncvardef (ncid, "rh3_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lwc_fcinv_id = ncvardef (ncid, "lwc_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ice_fcinv_id = ncvardef (ncid, "ice_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, rz_id, "long_name", NC_CHAR, 16, (void *)"MM5 model height");
   ncattput (ncid, rz_id, "units", NC_CHAR, 6, (void *)"meters");
   rz_valid_range[0] = -200;
   rz_valid_range[1] = 200;
   ncattput (ncid, rz_id, "valid_range", NC_FLOAT, 2, (void *) rz_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rz_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rz_id, "LAPS_var", NC_CHAR, 2, (void *)"HT");
   ncattput (ncid, rz_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rz_id, "LAPS_units", NC_CHAR, 6, (void *)"METERS");
   ncattput (ncid, ru_id, "long_name", NC_CHAR, 17, (void *)"MM5 eastward wind");
   ncattput (ncid, ru_id, "units", NC_CHAR, 13, (void *)"meters/second");
   ru_valid_range[0] = -200;
   ru_valid_range[1] = 200;
   ncattput (ncid, ru_id, "valid_range", NC_FLOAT, 2, (void *) ru_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, ru_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, ru_id, "LAPS_var", NC_CHAR, 2, (void *)"U3");
   ncattput (ncid, ru_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, ru_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, rv_id, "long_name", NC_CHAR, 18, (void *)"MM5 northward wind");
   ncattput (ncid, rv_id, "units", NC_CHAR, 13, (void *)"meters/second");
   rv_valid_range[0] = -200;
   rv_valid_range[1] = 200;
   ncattput (ncid, rv_id, "valid_range", NC_FLOAT, 2, (void *) rv_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rv_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rv_id, "LAPS_var", NC_CHAR, 2, (void *)"V3");
   ncattput (ncid, rv_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rv_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, rw_id, "long_name", NC_CHAR, 9, (void *)"MM5 omega");
   ncattput (ncid, rw_id, "units", NC_CHAR, 14, (void *)"pascals/second");
   rw_valid_range[0] = -20000;
   rw_valid_range[1] = 20000;
   ncattput (ncid, rw_id, "valid_range", NC_FLOAT, 2, (void *) rw_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rw_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rw_id, "LAPS_var", NC_CHAR, 2, (void *)"OM");
   ncattput (ncid, rw_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rw_id, "LAPS_units", NC_CHAR, 4, (void *)"PA/S");
   ncattput (ncid, rt_id, "long_name", NC_CHAR, 15, (void *)"MM5 temperature");
   ncattput (ncid, rt_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   rt_valid_range[0] = 0;
   rt_valid_range[1] = 100;
   ncattput (ncid, rt_id, "valid_range", NC_FLOAT, 2, (void *) rt_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rt_id, "LAPS_var", NC_CHAR, 2, (void *)"T3");
   ncattput (ncid, rt_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   char_val = 'K';
   ncattput (ncid, rt_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rsh_id, "long_name", NC_CHAR, 21, (void *)"MM5 specific humidity");
   ncattput (ncid, rsh_id, "units", NC_CHAR, 5, (void *)"kg/kg");
   rsh_valid_range[0] = 0;
   rsh_valid_range[1] = 0.1;
   ncattput (ncid, rsh_id, "valid_range", NC_FLOAT, 2, (void *) rsh_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rsh_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rsh_id, "LAPS_var", NC_CHAR, 2, (void *)"SH");
   ncattput (ncid, rsh_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, rsh_id, "LAPS_units", NC_CHAR, 10, (void *)"          ");
   ncattput (ncid, rh3_id, "long_name", NC_CHAR, 21, (void *)"MM5 relative humidity");
   ncattput (ncid, rh3_id, "units", NC_CHAR, 7, (void *)"percent");
   rh3_valid_range[0] = 0;
   rh3_valid_range[1] = 100;
   ncattput (ncid, rh3_id, "valid_range", NC_FLOAT, 2, (void *) rh3_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rh3_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rh3_id, "LAPS_var", NC_CHAR, 3, (void *)"RH3");
   ncattput (ncid, rh3_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, rh3_id, "LAPS_units", NC_CHAR, 7, (void *)"PERCENT");
   ncattput (ncid, lwc_id, "long_name", NC_CHAR, 22, (void *)"MM5 cloud liquid water");
   ncattput (ncid, lwc_id, "units", NC_CHAR, 14, (void *)"grams/meter**3");
   lwc_valid_range[0] = 0;
   lwc_valid_range[1] = 100;
   ncattput (ncid, lwc_id, "valid_range", NC_FLOAT, 2, (void *) lwc_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lwc_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lwc_id, "LAPS_var", NC_CHAR, 3, (void *)"LWC");
   ncattput (ncid, lwc_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, lwc_id, "LAPS_units", NC_CHAR, 6, (void *)"G/M**3");
   ncattput (ncid, ice_id, "long_name", NC_CHAR, 13, (void *)"MM5 cloud ice");
   ncattput (ncid, ice_id, "units", NC_CHAR, 14, (void *)"grams/meter**3");
   ice_valid_range[0] = 0;
   ice_valid_range[1] = 100;
   ncattput (ncid, ice_id, "valid_range", NC_FLOAT, 2, (void *) ice_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, ice_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, ice_id, "LAPS_var", NC_CHAR, 3, (void *)"ICE");
   ncattput (ncid, ice_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, ice_id, "LAPS_units", NC_CHAR, 6, (void *)"G/M**3");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, rz_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, ru_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rv_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rw_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rsh_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, rh3_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, lwc_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, ice_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store num_variables */
    static long num_variables = {39};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {21};
    static short level[] = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {132};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {132};
    static char model[] = {"MM5 - Mesoscale Model - ver. 5"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {189};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {189};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/*************************************************************************/
/* file to create netCDF format file with extension MSF */
#ifdef __STDC__
int cre_msf(char *fname)                /* create fname */
#else
int cre_msf(fname)
char *fname;
#endif
{                       /* create rsf.cdf */

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  rus_id, rvs_id, rts_id, rps_id, rtd_id, rh_id, lcb_id, lct_id, msl_id, 
        lil_id, tpw_id, r01_id, rto_id, s01_id, sto_id, th_id, the_id, pbe_id, 
        nbe_id, ps_id, cce_id, vis_id, lcv_id, lmt_id, spt_id, lhe_id, li_id, 
        hi_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        rus_comment_id, rvs_comment_id, rts_comment_id, rps_comment_id, rtd_comment_id, 
        rh_comment_id, lcb_comment_id, lct_comment_id, msl_comment_id, lil_comment_id, 
        tpw_comment_id, r01_comment_id, rto_comment_id, s01_comment_id, sto_comment_id, 
        th_comment_id, the_comment_id, pbe_comment_id, nbe_comment_id, ps_comment_id,
        cce_comment_id,vis_comment_id,lcv_comment_id,lmt_comment_id,spt_comment_id,
        lhe_comment_id, li_comment_id, hi_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, rus_fcinv_id, rvs_fcinv_id, rts_fcinv_id, rps_fcinv_id, 
        rtd_fcinv_id, rh_fcinv_id, lcb_fcinv_id, lct_fcinv_id, msl_fcinv_id, lil_fcinv_id, 
        tpw_fcinv_id, r01_fcinv_id, rto_fcinv_id, s01_fcinv_id, sto_fcinv_id, th_fcinv_id, 
        the_fcinv_id, pbe_fcinv_id, nbe_fcinv_id, ps_fcinv_id, cce_fcinv_id, vis_fcinv_id,
        lcv_fcinv_id, lmt_fcinv_id, spt_fcinv_id, lhe_fcinv_id, li_fcinv_id, hi_fcinv_id, 
        origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  rus_valid_range[2];
   float  rvs_valid_range[2];
   float  rts_valid_range[2];
   float  rps_valid_range[2];
   float  rtd_valid_range[2];
   float  rh_valid_range[2];
   float  lcb_valid_range[2];
   float  lct_valid_range[2];
   float  msl_valid_range[2];
   float  lil_valid_range[2];
   float  tpw_valid_range[2];
   float  r01_valid_range[2];
   float  rto_valid_range[2];
   float  s01_valid_range[2];
   float  sto_valid_range[2];
   float  th_valid_range[2];
   float  the_valid_range[2];
   float  pbe_valid_range[2];
   float  nbe_valid_range[2];
   float  ps_valid_range[2];
   float  cce_valid_range[2];
   float  vis_valid_range[2];
   float  lcv_valid_range[2];
   float  lmt_valid_range[2];
   float  spt_valid_range[2];
   float  lhe_valid_range[2];
   float  li_valid_range[2];
   float  hi_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 28L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rus_id = ncvardef (ncid, "rus", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rvs_id = ncvardef (ncid, "rvs", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rts_id = ncvardef (ncid, "rts", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rps_id = ncvardef (ncid, "rps", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rtd_id = ncvardef (ncid, "rtd", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rh_id = ncvardef (ncid, "rh", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lcb_id = ncvardef (ncid, "lcb", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lct_id = ncvardef (ncid, "lct", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   msl_id = ncvardef (ncid, "msl", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lil_id = ncvardef (ncid, "lil", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   tpw_id = ncvardef (ncid, "tpw", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   r01_id = ncvardef (ncid, "r01", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   rto_id = ncvardef (ncid, "rto", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s01_id = ncvardef (ncid, "s01", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   sto_id = ncvardef (ncid, "sto", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   th_id = ncvardef (ncid, "th", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   the_id = ncvardef (ncid, "the", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   pbe_id = ncvardef (ncid, "pbe", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   nbe_id = ncvardef (ncid, "nbe", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   ps_id = ncvardef (ncid, "ps", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   cce_id = ncvardef (ncid, "cce", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vis_id = ncvardef (ncid, "vis", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lcv_id = ncvardef (ncid, "lcv", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lmt_id = ncvardef (ncid, "lmt", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   spt_id = ncvardef (ncid, "spt", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   lhe_id = ncvardef (ncid, "lhe", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   li_id = ncvardef (ncid, "li", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   hi_id = ncvardef (ncid, "hi", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rus_comment_id = ncvardef (ncid, "rus_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rvs_comment_id = ncvardef (ncid, "rvs_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rts_comment_id = ncvardef (ncid, "rts_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rps_comment_id = ncvardef (ncid, "rps_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rtd_comment_id = ncvardef (ncid, "rtd_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rh_comment_id = ncvardef (ncid, "rh_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lcb_comment_id = ncvardef (ncid, "lcb_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lct_comment_id = ncvardef (ncid, "lct_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   msl_comment_id = ncvardef (ncid, "msl_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lil_comment_id = ncvardef (ncid, "lil_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   tpw_comment_id = ncvardef (ncid, "tpw_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   r01_comment_id = ncvardef (ncid, "r01_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   rto_comment_id = ncvardef (ncid, "rto_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s01_comment_id = ncvardef (ncid, "s01_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   sto_comment_id = ncvardef (ncid, "sto_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   th_comment_id = ncvardef (ncid, "th_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   the_comment_id = ncvardef (ncid, "the_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   pbe_comment_id = ncvardef (ncid, "pbe_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   nbe_comment_id = ncvardef (ncid, "nbe_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   ps_comment_id = ncvardef (ncid, "ps_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   cce_comment_id = ncvardef (ncid, "cce_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vis_comment_id = ncvardef (ncid, "vis_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lcv_comment_id = ncvardef (ncid, "lcv_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lmt_comment_id = ncvardef (ncid, "lmt_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   spt_comment_id = ncvardef (ncid, "spt_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   lhe_comment_id = ncvardef (ncid, "lhe_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   li_comment_id = ncvardef (ncid, "li_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   hi_comment_id = ncvardef (ncid, "hi_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rus_fcinv_id = ncvardef (ncid, "rus_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rvs_fcinv_id = ncvardef (ncid, "rvs_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rts_fcinv_id = ncvardef (ncid, "rts_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rps_fcinv_id = ncvardef (ncid, "rps_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rtd_fcinv_id = ncvardef (ncid, "rtd_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rh_fcinv_id = ncvardef (ncid, "rh_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lcb_fcinv_id = ncvardef (ncid, "lcb_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lct_fcinv_id = ncvardef (ncid, "lct_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   msl_fcinv_id = ncvardef (ncid, "msl_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lil_fcinv_id = ncvardef (ncid, "lil_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   tpw_fcinv_id = ncvardef (ncid, "tpw_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   r01_fcinv_id = ncvardef (ncid, "r01_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   rto_fcinv_id = ncvardef (ncid, "rto_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s01_fcinv_id = ncvardef (ncid, "s01_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   sto_fcinv_id = ncvardef (ncid, "sto_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   th_fcinv_id = ncvardef (ncid, "th_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   the_fcinv_id = ncvardef (ncid, "the_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   pbe_fcinv_id = ncvardef (ncid, "pbe_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   nbe_fcinv_id = ncvardef (ncid, "nbe_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   ps_fcinv_id = ncvardef (ncid, "ps_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   cce_fcinv_id = ncvardef (ncid, "cce_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vis_fcinv_id = ncvardef (ncid, "vis_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lcv_fcinv_id = ncvardef (ncid, "lcv_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lmt_fcinv_id = ncvardef (ncid, "lmt_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   spt_fcinv_id = ncvardef (ncid, "spt_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   lhe_fcinv_id = ncvardef (ncid, "lhe_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   li_fcinv_id = ncvardef (ncid, "li_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   hi_fcinv_id = ncvardef (ncid, "hi_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, rus_id, "long_name", NC_CHAR, 26, (void *)"MM5 Fcst sfc eastward wind");
   ncattput (ncid, rus_id, "units", NC_CHAR, 13, (void *)"meters/second");
   rus_valid_range[0] = -200;
   rus_valid_range[1] = 200;
   ncattput (ncid, rus_id, "valid_range", NC_FLOAT, 2, (void *) rus_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rus_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   char_val = 'U';
   ncattput (ncid, rus_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rus_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rus_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, rvs_id, "long_name", NC_CHAR, 27, (void *)"MM5 Fcst sfc northward wind");
   ncattput (ncid, rvs_id, "units", NC_CHAR, 13, (void *)"meters/second");
   rvs_valid_range[0] = -200;
   rvs_valid_range[1] = 200;
   ncattput (ncid, rvs_id, "valid_range", NC_FLOAT, 2, (void *) rvs_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rvs_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   char_val = 'V';
   ncattput (ncid, rvs_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rvs_id, "lvl_coord", NC_CHAR, 4, (void *)"HPA ");
   ncattput (ncid, rvs_id, "LAPS_units", NC_CHAR, 3, (void *)"M/S");
   ncattput (ncid, rts_id, "long_name", NC_CHAR, 24, (void *)"MM5 Fcst sfc temperature");
   ncattput (ncid, rts_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   rts_valid_range[0] = 0;
   rts_valid_range[1] = 100;
   ncattput (ncid, rts_id, "valid_range", NC_FLOAT, 2, (void *) rts_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rts_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   char_val = 'T';
   ncattput (ncid, rts_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rts_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   char_val = 'K';
   ncattput (ncid, rts_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rps_id, "long_name", NC_CHAR, 32, (void *)"MM5 Fcst 1500 M Reduced Pressure");
   ncattput (ncid, rps_id, "units", NC_CHAR, 7, (void *)"pascals");
   rps_valid_range[0] = 0;
   rps_valid_range[1] = 0.1;
   ncattput (ncid, rps_id, "valid_range", NC_FLOAT, 2, (void *) rps_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rps_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   char_val = 'P';
   ncattput (ncid, rps_id, "LAPS_var", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rps_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, rps_id, "LAPS_units", NC_CHAR, 2, (void *)"PA");
   ncattput (ncid, rtd_id, "long_name", NC_CHAR, 37, (void *)"MM5 Fcst surface dewpoint temperature");
   ncattput (ncid, rtd_id, "units", NC_CHAR, 14, (void *)"degrees kelvin");
   rtd_valid_range[0] = -75;
   rtd_valid_range[1] = 125;
   ncattput (ncid, rtd_id, "valid_range", NC_FLOAT, 2, (void *) rtd_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rtd_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rtd_id, "LAPS_var", NC_CHAR, 2, (void *)"TD");
   ncattput (ncid, rtd_id, "lvl_coord", NC_CHAR, 4, (void *)"AGL ");
   char_val = 'K';
   ncattput (ncid, rtd_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rh_id, "long_name", NC_CHAR, 26, (void *)"MM5 Fcst relative humidity");
   ncattput (ncid, rh_id, "units", NC_CHAR, 7, (void *)"percent");
   rh_valid_range[0] = 0;
   rh_valid_range[1] = 100;
   ncattput (ncid, rh_id, "valid_range", NC_FLOAT, 2, (void *) rh_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rh_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rh_id, "LAPS_var", NC_CHAR, 2, (void *)"RH");
   ncattput (ncid, rh_id, "lvl_coord", NC_CHAR, 3, (void *)"HPA");
   ncattput (ncid, rh_id, "LAPS_units", NC_CHAR, 7, (void *)"PERCENT");
   ncattput (ncid, lcb_id, "long_name", NC_CHAR, 19, (void *)"MM5 Fcst cloud base");
   ncattput (ncid, lcb_id, "units", NC_CHAR, 6, (void *)"meters");
   lcb_valid_range[0] = 0;
   lcb_valid_range[1] = 0.1;
   ncattput (ncid, lcb_id, "valid_range", NC_FLOAT, 2, (void *) lcb_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lcb_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lcb_id, "LAPS_var", NC_CHAR, 3, (void *)"LCB");
   ncattput (ncid, lcb_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, lcb_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, lct_id, "long_name", NC_CHAR, 18, (void *)"MM5 Fcst cloud top");
   ncattput (ncid, lct_id, "units", NC_CHAR, 6, (void *)"meters");
   lct_valid_range[0] = 0;
   lct_valid_range[1] = 0.1;
   ncattput (ncid, lct_id, "valid_range", NC_FLOAT, 2, (void *) lct_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lct_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lct_id, "LAPS_var", NC_CHAR, 3, (void *)"LCT");
   ncattput (ncid, lct_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, lct_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, msl_id, "long_name", NC_CHAR, 21, (void *)"MM5 Fcst MSL pressure");
   ncattput (ncid, msl_id, "units", NC_CHAR, 7, (void *)"pascals");
   msl_valid_range[0] = -20000;
   msl_valid_range[1] = 20000;
   ncattput (ncid, msl_id, "valid_range", NC_FLOAT, 2, (void *) msl_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, msl_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, msl_id, "LAPS_var", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, msl_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, msl_id, "LAPS_units", NC_CHAR, 2, (void *)"PA");
   ncattput (ncid, lil_id, "long_name", NC_CHAR, 32, (void *)"MM5 Fcst integrated liquid water");
   ncattput (ncid, lil_id, "units", NC_CHAR, 14, (void *)"grams/meter**2");
   lil_valid_range[0] = 0;
   lil_valid_range[1] = 0.1;
   ncattput (ncid, lil_id, "valid_range", NC_FLOAT, 2, (void *) lil_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lil_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lil_id, "LAPS_var", NC_CHAR, 3, (void *)"LIL");
   ncattput (ncid, lil_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, lil_id, "LAPS_units", NC_CHAR, 6, (void *)"G/M**2");
   ncattput (ncid, tpw_id, "long_name", NC_CHAR, 44, (void *)"MM5 Fcst integrated total precipitable water");
   ncattput (ncid, tpw_id, "units", NC_CHAR, 6, (void *)"meters");
   tpw_valid_range[0] = 0;
   tpw_valid_range[1] = 0.1;
   ncattput (ncid, tpw_id, "valid_range", NC_FLOAT, 2, (void *) tpw_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, tpw_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, tpw_id, "LAPS_var", NC_CHAR, 3, (void *)"TPW");
   ncattput (ncid, tpw_id, "lvl_coord", NC_CHAR, 4, (void *)"none");
   char_val = 'M';
   ncattput (ncid, tpw_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, r01_id, "long_name", NC_CHAR, 28, (void *)"MM5 Fcst cycle precip.accum.");
   ncattput (ncid, r01_id, "units", NC_CHAR, 6, (void *)"meters");
   r01_valid_range[0] = 0;
   r01_valid_range[1] = 200;
   ncattput (ncid, r01_id, "valid_range", NC_FLOAT, 2, (void *) r01_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, r01_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, r01_id, "LAPS_var", NC_CHAR, 3, (void *)"R01");
   ncattput (ncid, r01_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, r01_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, rto_id, "long_name", NC_CHAR, 35, (void *)"MM5 Fcst storm total precip. accum.");
   ncattput (ncid, rto_id, "units", NC_CHAR, 6, (void *)"meters");
   rto_valid_range[0] = 0;
   rto_valid_range[1] = 200;
   ncattput (ncid, rto_id, "valid_range", NC_FLOAT, 2, (void *) rto_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, rto_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, rto_id, "LAPS_var", NC_CHAR, 3, (void *)"RTO");
   ncattput (ncid, rto_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, rto_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, s01_id, "long_name", NC_CHAR, 26, (void *)"MM5 Fcst cycle snow accum.");
   ncattput (ncid, s01_id, "units", NC_CHAR, 6, (void *)"meters");
   s01_valid_range[0] = 0;
   s01_valid_range[1] = 200;
   ncattput (ncid, s01_id, "valid_range", NC_FLOAT, 2, (void *) s01_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s01_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s01_id, "LAPS_var", NC_CHAR, 3, (void *)"S01");
   ncattput (ncid, s01_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, s01_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, sto_id, "long_name", NC_CHAR, 38, (void *)"MM5 Fcst storm total snow accumulation");
   ncattput (ncid, sto_id, "units", NC_CHAR, 6, (void *)"meters");
   sto_valid_range[0] = 0;
   sto_valid_range[1] = 200;
   ncattput (ncid, sto_id, "valid_range", NC_FLOAT, 2, (void *) sto_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, sto_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, sto_id, "LAPS_var", NC_CHAR, 3, (void *)"STO");
   ncattput (ncid, sto_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, sto_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, th_id, "long_name", NC_CHAR, 30, (void *)"MM5 Fcst potential temperature");
   ncattput (ncid, th_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   th_valid_range[0] = -75;
   th_valid_range[1] = 125;
   ncattput (ncid, th_id, "valid_range", NC_FLOAT, 2, (void *) th_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, th_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, th_id, "LAPS_var", NC_CHAR, 2, (void *)"TH");
   ncattput (ncid, th_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (ncid, th_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, the_id, "long_name", NC_CHAR, 41, (void *)"MM5 Fcst equivalent potential temperature");
   ncattput (ncid, the_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   the_valid_range[0] = -20000;
   the_valid_range[1] = 20000;
   ncattput (ncid, the_id, "valid_range", NC_FLOAT, 2, (void *) the_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, the_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, the_id, "LAPS_var", NC_CHAR, 3, (void *)"THE");
   ncattput (ncid, the_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (ncid, the_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, pbe_id, "long_name", NC_CHAR, 32, (void *)"MM5 Fcst positive buoyant energy");
   ncattput (ncid, pbe_id, "units", NC_CHAR, 15, (void *)"joules/kilogram");
   pbe_valid_range[0] = -20000;
   pbe_valid_range[1] = 20000;
   ncattput (ncid, pbe_id, "valid_range", NC_FLOAT, 2, (void *) pbe_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, pbe_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, pbe_id, "LAPS_var", NC_CHAR, 3, (void *)"PBE");
   ncattput (ncid, pbe_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, pbe_id, "LAPS_units", NC_CHAR, 4, (void *)"J/KG");
   ncattput (ncid, nbe_id, "long_name", NC_CHAR, 32, (void *)"MM5 Fcst negative buoyant energy");
   ncattput (ncid, nbe_id, "units", NC_CHAR, 15, (void *)"joules/kilogram");
   nbe_valid_range[0] = -20000;
   nbe_valid_range[1] = 20000;
   ncattput (ncid, nbe_id, "valid_range", NC_FLOAT, 2, (void *) nbe_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, nbe_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, nbe_id, "LAPS_var", NC_CHAR, 3, (void *)"NBE");
   ncattput (ncid, nbe_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, nbe_id, "LAPS_units", NC_CHAR, 4, (void *)"J/KG");
   ncattput (ncid, ps_id, "long_name", NC_CHAR, 25, (void *)"MM5 Fcst surface pressure");
   ncattput (ncid, ps_id, "units", NC_CHAR, 7, (void *)"pascals");
   ps_valid_range[0] = -20000;
   ps_valid_range[1] = 20000;
   ncattput (ncid, ps_id, "valid_range", NC_FLOAT, 2, (void *) ps_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, ps_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, ps_id, "LAPS_var", NC_CHAR, 2, (void *)"PS");
   ncattput (ncid, ps_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, ps_id, "LAPS_units", NC_CHAR, 2, (void *)"PA");
   ncattput (ncid, cce_id, "long_name", NC_CHAR, 22, (void *)"MM5 Fcst cloud ceiling");
   ncattput (ncid, cce_id, "units", NC_CHAR, 6, (void *)"meters");
   cce_valid_range[0] = -20000;
   cce_valid_range[1] = 20000;
   ncattput (ncid, cce_id, "valid_range", NC_FLOAT, 2, (void *) cce_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, cce_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, cce_id, "LAPS_var", NC_CHAR, 3, (void *)"CCE");
   ncattput (ncid, cce_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'M';
   ncattput (ncid, cce_id, "LAPS_units", NC_CHAR, 1, (void *)&char_val);
   ncattput (ncid, vis_id, "long_name", NC_CHAR, 19, (void *)"MM5 Fcst visibility");
   ncattput (ncid, vis_id, "units", NC_CHAR, 6, (void *)"meters");
   vis_valid_range[0] = -20000;
   vis_valid_range[1] = 20000;
   ncattput (ncid, vis_id, "valid_range", NC_FLOAT, 2, (void *) vis_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, vis_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, vis_id, "LAPS_var", NC_CHAR, 3, (void *)"VIS");
   ncattput (ncid, vis_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'M';
   ncattput (ncid, vis_id, "LAPS_units", NC_CHAR, 1, (void *)&char_val);

   ncattput (ncid, lcv_id, "long_name", NC_CHAR, 20, (void *)"MM5 Fcst cloud cover");
   ncattput (ncid, lcv_id, "units", NC_CHAR, 4, (void *)"none");
   lcv_valid_range[0] = 0.0;
   lcv_valid_range[1] = 0.100;
   ncattput (ncid, lcv_id, "valid_range", NC_FLOAT, 2, (void *) lcv_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lcv_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lcv_id, "LAPS_var", NC_CHAR, 3, (void *)"LCB");
   ncattput (ncid, lcv_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, lcv_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, lmt_id, "long_name", NC_CHAR, 22, (void *)"MM5 Fcst max echo tops");
   ncattput (ncid, lmt_id, "units", NC_CHAR, 6, (void *)"meters");
   lmt_valid_range[0] = -20000;
   lmt_valid_range[1] = 20000;
   ncattput (ncid, lmt_id, "valid_range", NC_FLOAT, 2, (void *) lmt_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lmt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lmt_id, "LAPS_var", NC_CHAR, 3, (void *)"LMT");
   ncattput (ncid, lmt_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   char_val = 'M';
   ncattput (ncid, lmt_id, "LAPS_units", NC_CHAR, 1, (void *)&char_val);

   ncattput (ncid, spt_id, "long_name", NC_CHAR, 28, (void *)"MM5 Fcst surface precip type");
   ncattput (ncid, spt_id, "units", NC_CHAR, 4, (void *)"none");
   spt_valid_range[0] = 0.0;
   spt_valid_range[1] = 0.100;
   ncattput (ncid, spt_id, "valid_range", NC_FLOAT, 2, (void *) spt_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, spt_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, spt_id, "LAPS_var", NC_CHAR, 3, (void *)"SPT");
   ncattput (ncid, spt_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, spt_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, lhe_id, "long_name", NC_CHAR, 17, (void *)"MM5 Fcst helicity");
   ncattput (ncid, lhe_id, "units", NC_CHAR, 16, (void *)"meters/second**2");
   lhe_valid_range[0] = 0;
   lhe_valid_range[1] = 0.1;
   ncattput (ncid, lhe_id, "valid_range", NC_FLOAT, 2, (void *) lhe_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, lhe_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, lhe_id, "LAPS_var", NC_CHAR, 3, (void *)"LHE");
   ncattput (ncid, lhe_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, lhe_id, "LAPS_units", NC_CHAR, 6, (void *)"M/S**2");
   ncattput (ncid, li_id, "long_name", NC_CHAR, 21, (void *)"MM5 Fcst lifted index");
   ncattput (ncid, li_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   li_valid_range[0] = -20000;
   li_valid_range[1] = 20000;
   ncattput (ncid, li_id, "valid_range", NC_FLOAT, 2, (void *) li_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, li_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, li_id, "LAPS_var", NC_CHAR, 2, (void *)"LI");
   ncattput (ncid, li_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   char_val = 'K';
   ncattput (ncid, li_id, "LAPS_units", NC_CHAR, 1,(void *) &char_val);
   ncattput (ncid, hi_id, "long_name", NC_CHAR, 19, (void *)"MM5 Fcst Heat index");
   ncattput (ncid, hi_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   hi_valid_range[0] = 210;
   hi_valid_range[1] = 366;
   ncattput (ncid, hi_id, "valid_range", NC_FLOAT, 2, (void *) hi_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, hi_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, hi_id, "LAPS_var", NC_CHAR, 2, (void *)"HI");
   ncattput (ncid, hi_id, "lvl_coord", NC_CHAR, 3, (void *)"AGL");
   ncattput (ncid, hi_id, "LAPS_units", NC_CHAR, 1, (void *)"K");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, rus_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rvs_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rts_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rps_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rtd_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rh_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lcb_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lct_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, msl_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lil_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, tpw_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, r01_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, rto_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, s01_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, sto_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, th_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, the_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, pbe_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, nbe_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, ps_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, cce_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, vis_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lcv_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lmt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, spt_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, lhe_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, li_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   ncattput (ncid, hi_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store num_variables */
    static long num_variables = {96};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {49};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {29};
    static char model[] = {"MM5 - Mesoscale Model ver. 5"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {28};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {28};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/************************************************************************/
/* file to create netCDF format file with extension SST */
#ifdef __STDC__
int cre_sst(char *fname)
#else
int cre_sst(fname)              /* create fname */
char *fname;
#endif
{
   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  st1_id, st2_id, st3_id, st4_id, st5_id, lvl_id, imax_id, jmax_id, 
        kmax_id, kdim_id, st1_comment_id, st2_comment_id, st3_comment_id, 
        st4_comment_id, st5_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, st1_fcinv_id, st2_fcinv_id, st3_fcinv_id, 
        st4_fcinv_id, st5_fcinv_id, origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  st1_valid_range[2];
   float  st2_valid_range[2];
   float  st3_valid_range[2];
   float  st4_valid_range[2];
   float  st5_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 5L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   st1_id = ncvardef (ncid, "st1", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   st2_id = ncvardef (ncid, "st2", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   st3_id = ncvardef (ncid, "st3", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   st4_id = ncvardef (ncid, "st4", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   st5_id = ncvardef (ncid, "st5", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   st1_comment_id = ncvardef (ncid, "st1_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   st2_comment_id = ncvardef (ncid, "st2_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   st3_comment_id = ncvardef (ncid, "st3_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   st4_comment_id = ncvardef (ncid, "st4_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   st5_comment_id = ncvardef (ncid, "st5_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   st1_fcinv_id = ncvardef (ncid, "st1_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   st2_fcinv_id = ncvardef (ncid, "st2_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   st3_fcinv_id = ncvardef (ncid, "st3_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   st4_fcinv_id = ncvardef (ncid, "st4_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   st5_fcinv_id = ncvardef (ncid, "st5_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, st1_id, "long_name", NC_CHAR, 18, (void *)"sea sfc temp var 1");
   ncattput (ncid, st1_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   st1_valid_range[0] = 210;
   st1_valid_range[1] = 366;
   ncattput (ncid, st1_id, "valid_range", NC_FLOAT, 2, (void *) st1_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, st1_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, st1_id, "LAPS_var", NC_CHAR, 3, (void *)"ST1");
   ncattput (ncid, st1_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, st1_id, "LAPS_units", NC_CHAR, 5, (void *)"DEG K");
   ncattput (ncid, st2_id, "long_name", NC_CHAR, 18, (void *)"sea sfc temp var 2");
   ncattput (ncid, st2_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   st2_valid_range[0] = 0;
   st2_valid_range[1] = 0.1;
   ncattput (ncid, st2_id, "valid_range", NC_FLOAT, 2, (void *) st2_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, st2_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, st2_id, "LAPS_var", NC_CHAR, 3, (void *)"ST2");
   ncattput (ncid, st2_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, st2_id, "LAPS_units", NC_CHAR, 5, (void *)"DEG K");
   ncattput (ncid, st3_id, "long_name", NC_CHAR, 18, (void *)"sea sfc temp var 3");
   ncattput (ncid, st3_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   st3_valid_range[0] = 210;
   st3_valid_range[1] = 366;
   ncattput (ncid, st3_id, "valid_range", NC_FLOAT, 2, (void *) st3_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, st3_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, st3_id, "LAPS_var", NC_CHAR, 3, (void *)"ST3");
   ncattput (ncid, st3_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, st3_id, "LAPS_units", NC_CHAR, 5, (void *)"DEG K");
   ncattput (ncid, st4_id, "long_name", NC_CHAR, 18, (void *)"sea sfc temp var 4");
   ncattput (ncid, st4_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   st4_valid_range[0] = 210;
   st4_valid_range[1] = 366;
   ncattput (ncid, st4_id, "valid_range", NC_FLOAT, 2, (void *) st4_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, st4_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, st4_id, "LAPS_var", NC_CHAR, 3, (void *)"ST4");
   ncattput (ncid, st4_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, st4_id, "LAPS_units", NC_CHAR, 5, (void *)"DEG K");
   ncattput (ncid, st5_id, "long_name", NC_CHAR, 18, (void *)"sea sfc temp var 5");
   ncattput (ncid, st5_id, "units", NC_CHAR, 14, (void *)"degrees Kelvin");
   st5_valid_range[0] = 210;
   st5_valid_range[1] = 366;
   ncattput (ncid, st5_id, "valid_range", NC_FLOAT, 2, (void *) st5_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, st5_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, st5_id, "LAPS_var", NC_CHAR, 3, (void *)"ST5");
   ncattput (ncid, st5_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, st5_id, "LAPS_units", NC_CHAR, 5, (void *)"DEG K");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, st1_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, st2_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, st3_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, st4_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, st5_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {27};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {132};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {132};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {5};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {5};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid; 
}
/************************************************************************/
/* file to create netCDF format file with extension VEG*/
#ifdef __STDC__
int cre_veg(char *fname)
#else
int cre_veg(fname)              /* create fname */
char *fname;
#endif
{
   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, n_2d_grids_dim, 
        var_len_dim, coord_len_dim, unit_len_dim, comm_len_dim, domain_len_dim, 
        asc_len_dim;

   /* variable ids */
   int  vg1_id, vg2_id, vg3_id, vg4_id, vg5_id, lvl_id, imax_id, jmax_id, 
        kmax_id, kdim_id, vg1_comment_id, vg2_comment_id, vg3_comment_id, 
        vg4_comment_id, vg5_comment_id, laps_domain_file_id, asctime_id, 
        fctimes_id, level_id, vg1_fcinv_id, vg2_fcinv_id, vg3_fcinv_id, 
        vg4_fcinv_id, vg5_fcinv_id, origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   char  char_val;
   short  short_val;
   float  float_val;

   /* attribute vectors */
   float  vg1_valid_range[2];
   float  vg2_valid_range[2];
   float  vg3_valid_range[2];
   float  vg4_valid_range[2];
   float  vg5_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 5L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vg1_id = ncvardef (ncid, "vg1", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vg2_id = ncvardef (ncid, "vg2", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vg3_id = ncvardef (ncid, "vg3", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vg4_id = ncvardef (ncid, "vg4", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   vg5_id = ncvardef (ncid, "vg5", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vg1_comment_id = ncvardef (ncid, "vg1_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vg2_comment_id = ncvardef (ncid, "vg2_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vg3_comment_id = ncvardef (ncid, "vg3_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vg4_comment_id = ncvardef (ncid, "vg4_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   vg5_comment_id = ncvardef (ncid, "vg5_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vg1_fcinv_id = ncvardef (ncid, "vg1_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vg2_fcinv_id = ncvardef (ncid, "vg2_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vg3_fcinv_id = ncvardef (ncid, "vg3_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vg4_fcinv_id = ncvardef (ncid, "vg4_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   vg5_fcinv_id = ncvardef (ncid, "vg5_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, vg1_id, "long_name", NC_CHAR, 16, (void *)"vegetation var 1");
   char_val = '\000';
   ncattput (ncid, vg1_id, "units", NC_CHAR, 0,(void *) &char_val);
   vg1_valid_range[0] = 0;
   vg1_valid_range[1] = 10000;
   ncattput (ncid, vg1_id, "valid_range", NC_FLOAT, 2, (void *) vg1_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, vg1_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, vg1_id, "LAPS_var", NC_CHAR, 3, (void *)"ST1");
   ncattput (ncid, vg1_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, vg1_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, vg2_id, "long_name", NC_CHAR, 16, (void *)"vegetation var 2");
   ncattput (ncid, vg2_id, "units", NC_CHAR, 4, (void *)"none");
   vg2_valid_range[0] = 0;
   vg2_valid_range[1] = 0.1;
   ncattput (ncid, vg2_id, "valid_range", NC_FLOAT, 2, (void *) vg2_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, vg2_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, vg2_id, "LAPS_var", NC_CHAR, 3, (void *)"VG2");
   ncattput (ncid, vg2_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, vg2_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, vg3_id, "long_name", NC_CHAR, 16, (void *)"vegetation var 3");
   char_val = '\000';
   ncattput (ncid, vg3_id, "units", NC_CHAR, 0,(void *) &char_val);
   vg3_valid_range[0] = 0;
   vg3_valid_range[1] = 10000;
   ncattput (ncid, vg3_id, "valid_range", NC_FLOAT, 2, (void *) vg3_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, vg3_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, vg3_id, "LAPS_var", NC_CHAR, 3, (void *)"VG3");
   ncattput (ncid, vg3_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, vg3_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, vg4_id, "long_name", NC_CHAR, 16, (void *)"vegetation var 4");
   char_val = '\000';
   ncattput (ncid, vg4_id, "units", NC_CHAR, 0,(void *) &char_val);
   vg4_valid_range[0] = 0;
   vg4_valid_range[1] = 10000;
   ncattput (ncid, vg4_id, "valid_range", NC_FLOAT, 2, (void *) vg4_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, vg4_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, vg4_id, "LAPS_var", NC_CHAR, 3, (void *)"VG4");
   ncattput (ncid, vg4_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, vg4_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, vg5_id, "long_name", NC_CHAR, 16, (void *)"vegetation var 5");
   char_val = '\000';
   ncattput (ncid, vg5_id, "units", NC_CHAR, 0,(void *) &char_val);
   vg5_valid_range[0] = 0;
   vg5_valid_range[1] = 10000;
   ncattput (ncid, vg5_id, "valid_range", NC_FLOAT, 2, (void *) vg5_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, vg5_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, vg5_id, "LAPS_var", NC_CHAR, 3, (void *)"VG5");
   ncattput (ncid, vg5_id, "lvl_coord", NC_CHAR, 3, (void *)"MSL");
   ncattput (ncid, vg5_id, "LAPS_units", NC_CHAR, 5, (void *)"UNDIM");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, vg1_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, vg2_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, vg3_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, vg4_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, vg5_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {27};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {132};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {132};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {5};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {5};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/************************************************************************/
/* file to create netCDF format file with extension LS8*/
#ifdef __STDC__
int cre_ls8(char *fname)
#else
int cre_ls8(fname)              /* create fname */
char *fname;
#endif
{
   int  ncid;			/* netCDF id */

   /* dimension ids */
   int  level_dim, fctimes_dim, namelen_dim, lat_dim, lon_dim, 
        n_2d_grids_dim, var_len_dim, coord_len_dim, unit_len_dim, 
        comm_len_dim, domain_len_dim, asc_len_dim;

   /* variable ids */
   int  s01_id, s02_id, s03_id, s04_id, s05_id, s06_id, s07_id, s08_id, 
        s09_id, s10_id, s11_id, s12_id, s13_id, s14_id, s15_id, s16_id, 
        s17_id, s18_id, s19_id, lvl_id, imax_id, jmax_id, kmax_id, kdim_id, 
        s01_comment_id, s02_comment_id, s03_comment_id, s04_comment_id, 
        s05_comment_id, s06_comment_id, s07_comment_id, s08_comment_id, 
        s09_comment_id, s10_comment_id, s11_comment_id, s12_comment_id, 
        s13_comment_id, s14_comment_id, s15_comment_id, s16_comment_id, 
        s17_comment_id, s18_comment_id, s19_comment_id, laps_domain_file_id, 
        asctime_id, fctimes_id, level_id, s01_fcinv_id, s02_fcinv_id, 
        s03_fcinv_id, s04_fcinv_id, s05_fcinv_id, s06_fcinv_id, s07_fcinv_id, 
        s08_fcinv_id, s09_fcinv_id, s10_fcinv_id, s11_fcinv_id, s12_fcinv_id, 
        s13_fcinv_id, s14_fcinv_id, s15_fcinv_id, s16_fcinv_id, s17_fcinv_id, 
        s18_fcinv_id, s19_fcinv_id, origin_id, model_id, version_id, num_variables_id;

   /* variable shapes */
   int dims[4];

   /* containers for scalar attributes */
   short  short_val;
   float  float_val;
   double  double_val;

   /* attribute vectors */
   float  s01_valid_range[2];
   float  s02_valid_range[2];
   float  s03_valid_range[2];
   float  s04_valid_range[2];
   float  s05_valid_range[2];
   float  s06_valid_range[2];
   float  s07_valid_range[2];
   float  s08_valid_range[2];
   float  s09_valid_range[2];
   float  s10_valid_range[2];
   float  s11_valid_range[2];
   float  s12_valid_range[2];
   float  s13_valid_range[2];
   float  s14_valid_range[2];
   float  s15_valid_range[2];
   float  s16_valid_range[2];
   float  s17_valid_range[2];
   float  s18_valid_range[2];
   float  s19_valid_range[2];

   /* enter define mode */
   ncid = nccreate(fname, NC_CLOBBER);  /* returns -1 if error */
   if (ncid == -1) return ncid;

   /* define dimensions */
   level_dim = ncdimdef(ncid, "level", 1L);
   fctimes_dim = ncdimdef(ncid, "fctimes", 1L);
   namelen_dim = ncdimdef(ncid, "namelen", 132L);
   lat_dim = ncdimdef(ncid, "lat", NY_LONG);
   lon_dim = ncdimdef(ncid, "lon", NX_LONG);
   n_2d_grids_dim = ncdimdef(ncid, "n_2d_grids", 19L);
   var_len_dim = ncdimdef(ncid, "var_len", 4L);
   coord_len_dim = ncdimdef(ncid, "coord_len", 5L);
   unit_len_dim = ncdimdef(ncid, "unit_len", 11L);
   comm_len_dim = ncdimdef(ncid, "comm_len", 126L);
   domain_len_dim = ncdimdef(ncid, "domain_len", 12L);
   asc_len_dim = ncdimdef(ncid, "asc_len", 18L);

   /* define variables */

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s01_id = ncvardef (ncid, "s01", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s02_id = ncvardef (ncid, "s02", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s03_id = ncvardef (ncid, "s03", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s04_id = ncvardef (ncid, "s04", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s05_id = ncvardef (ncid, "s05", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s06_id = ncvardef (ncid, "s06", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s07_id = ncvardef (ncid, "s07", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s08_id = ncvardef (ncid, "s08", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s09_id = ncvardef (ncid, "s09", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s10_id = ncvardef (ncid, "s10", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s11_id = ncvardef (ncid, "s11", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s12_id = ncvardef (ncid, "s12", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s13_id = ncvardef (ncid, "s13", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s14_id = ncvardef (ncid, "s14", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s15_id = ncvardef (ncid, "s15", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s16_id = ncvardef (ncid, "s16", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s17_id = ncvardef (ncid, "s17", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s18_id = ncvardef (ncid, "s18", NC_FLOAT, 4, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = lat_dim;
   dims[3] = lon_dim;
   s19_id = ncvardef (ncid, "s19", NC_FLOAT, 4, dims);

   dims[0] = n_2d_grids_dim;
   lvl_id = ncvardef (ncid, "lvl", NC_LONG, 1, dims);

   imax_id = ncvardef (ncid, "imax", NC_LONG, 0, 0);

   jmax_id = ncvardef (ncid, "jmax", NC_LONG, 0, 0);

   kmax_id = ncvardef (ncid, "kmax", NC_LONG, 0, 0);

   kdim_id = ncvardef (ncid, "kdim", NC_LONG, 0, 0);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s01_comment_id = ncvardef (ncid, "s01_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s02_comment_id = ncvardef (ncid, "s02_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s03_comment_id = ncvardef (ncid, "s03_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s04_comment_id = ncvardef (ncid, "s04_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s05_comment_id = ncvardef (ncid, "s05_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s06_comment_id = ncvardef (ncid, "s06_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s07_comment_id = ncvardef (ncid, "s07_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s08_comment_id = ncvardef (ncid, "s08_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s09_comment_id = ncvardef (ncid, "s09_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s10_comment_id = ncvardef (ncid, "s10_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s11_comment_id = ncvardef (ncid, "s11_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s12_comment_id = ncvardef (ncid, "s12_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s13_comment_id = ncvardef (ncid, "s13_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s14_comment_id = ncvardef (ncid, "s14_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s15_comment_id = ncvardef (ncid, "s15_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s16_comment_id = ncvardef (ncid, "s16_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s17_comment_id = ncvardef (ncid, "s17_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s18_comment_id = ncvardef (ncid, "s18_comment", NC_CHAR, 3, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   dims[2] = comm_len_dim;
   s19_comment_id = ncvardef (ncid, "s19_comment", NC_CHAR, 3, dims);

   dims[0] = domain_len_dim;
   laps_domain_file_id = ncvardef (ncid, "laps_domain_file", NC_CHAR, 1, dims);

   dims[0] = asc_len_dim;
   asctime_id = ncvardef (ncid, "asctime", NC_CHAR, 1, dims);

   dims[0] = fctimes_dim;
   fctimes_id = ncvardef (ncid, "fctimes", NC_SHORT, 1, dims);

   dims[0] = level_dim;
   level_id = ncvardef (ncid, "level", NC_SHORT, 1, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s01_fcinv_id = ncvardef (ncid, "s01_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s02_fcinv_id = ncvardef (ncid, "s02_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s03_fcinv_id = ncvardef (ncid, "s03_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s04_fcinv_id = ncvardef (ncid, "s04_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s05_fcinv_id = ncvardef (ncid, "s05_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s06_fcinv_id = ncvardef (ncid, "s06_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s07_fcinv_id = ncvardef (ncid, "s07_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s08_fcinv_id = ncvardef (ncid, "s08_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s09_fcinv_id = ncvardef (ncid, "s09_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s10_fcinv_id = ncvardef (ncid, "s10_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s11_fcinv_id = ncvardef (ncid, "s11_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s12_fcinv_id = ncvardef (ncid, "s12_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s13_fcinv_id = ncvardef (ncid, "s13_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s14_fcinv_id = ncvardef (ncid, "s14_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s15_fcinv_id = ncvardef (ncid, "s15_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s16_fcinv_id = ncvardef (ncid, "s16_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s17_fcinv_id = ncvardef (ncid, "s17_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s18_fcinv_id = ncvardef (ncid, "s18_fcinv", NC_SHORT, 2, dims);

   dims[0] = fctimes_dim;
   dims[1] = level_dim;
   s19_fcinv_id = ncvardef (ncid, "s19_fcinv", NC_SHORT, 2, dims);

   dims[0] = namelen_dim;
   origin_id = ncvardef (ncid, "origin", NC_CHAR, 1, dims);

   dims[0] = namelen_dim;
   model_id = ncvardef (ncid, "model", NC_CHAR, 1, dims);

   version_id = ncvardef (ncid, "version", NC_LONG, 0, 0);

   num_variables_id = ncvardef (ncid, "num_variables", NC_LONG, 0, 0);

   /* assign attributes */
   ncattput (ncid, s01_id, "long_name", NC_CHAR, 37, (void *)"GOES-8 14.71 micron sounding radiance");
   ncattput (ncid, s01_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.001471;
   ncattput (ncid, s01_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s01_valid_range[0] = 0;
   s01_valid_range[1] = 100;
   ncattput (ncid, s01_id, "valid_range", NC_FLOAT, 2, (void *) s01_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s01_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s01_id, "LAPS_var", NC_CHAR, 3, (void *)"S01");
   ncattput (ncid, s01_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s01_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s02_id, "long_name", NC_CHAR, 37, (void *)"GOES-8 14.37 micron sounding radiance");
   ncattput (ncid, s02_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.001437;
   ncattput (ncid, s02_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s02_valid_range[0] = 0;
   s02_valid_range[1] = 100;
   ncattput (ncid, s02_id, "valid_range", NC_FLOAT, 2, (void *) s02_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s02_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s02_id, "LAPS_var", NC_CHAR, 3, (void *)"S02");
   ncattput (ncid, s02_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s02_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s03_id, "long_name", NC_CHAR, 37, (void *)"GOES-8 14.06 micron sounding radiance");
   ncattput (ncid, s03_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.001406;
   ncattput (ncid, s03_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s03_valid_range[0] = 0;
   s03_valid_range[1] = 100;
   ncattput (ncid, s03_id, "valid_range", NC_FLOAT, 2, (void *) s03_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s03_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s03_id, "LAPS_var", NC_CHAR, 3, (void *)"S03");
   ncattput (ncid, s03_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s03_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s04_id, "long_name", NC_CHAR, 37, (void *)"GOES-8 13.64 micron sounding radiance");
   ncattput (ncid, s04_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.001364;
   ncattput (ncid, s04_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s04_valid_range[0] = 0;
   s04_valid_range[1] = 100;
   ncattput (ncid, s04_id, "valid_range", NC_FLOAT, 2, (void *) s04_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s04_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s04_id, "LAPS_var", NC_CHAR, 3, (void *)"S04");
   ncattput (ncid, s04_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s04_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s05_id, "long_name", NC_CHAR, 37, (void *)"GOES-8 13.37 micron sounding radiance");
   ncattput (ncid, s05_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.001337;
   ncattput (ncid, s05_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s05_valid_range[0] = 0;
   s05_valid_range[1] = 100;
   ncattput (ncid, s05_id, "valid_range", NC_FLOAT, 2, (void *) s05_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s05_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s05_id, "LAPS_var", NC_CHAR, 3, (void *)"S05");
   ncattput (ncid, s05_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s05_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s06_id, "long_name", NC_CHAR, 37, (void *)"GOES-8 12.66 micron sounding radiance");
   ncattput (ncid, s06_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.001266;
   ncattput (ncid, s06_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s06_valid_range[0] = 0;
   s06_valid_range[1] = 100;
   ncattput (ncid, s06_id, "valid_range", NC_FLOAT, 2, (void *) s06_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s06_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s06_id, "LAPS_var", NC_CHAR, 3, (void *)"S06");
   ncattput (ncid, s06_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s06_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s07_id, "long_name", NC_CHAR, 37, (void *)"GOES-8 12.02 micron sounding radiance");
   ncattput (ncid, s07_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.001202;
   ncattput (ncid, s07_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s07_valid_range[0] = 0;
   s07_valid_range[1] = 100;
   ncattput (ncid, s07_id, "valid_range", NC_FLOAT, 2, (void *) s07_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s07_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s07_id, "LAPS_var", NC_CHAR, 3, (void *)"S07");
   ncattput (ncid, s07_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s07_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s08_id, "long_name", NC_CHAR, 37, (void *)"GOES-8 11.03 micron sounding radiance");
   ncattput (ncid, s08_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.001103;
   ncattput (ncid, s08_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s08_valid_range[0] = 0;
   s08_valid_range[1] = 100;
   ncattput (ncid, s08_id, "valid_range", NC_FLOAT, 2, (void *) s08_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s08_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s08_id, "LAPS_var", NC_CHAR, 3, (void *)"S08");
   ncattput (ncid, s08_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s08_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s09_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 9.71 micron sounding radiance");
   ncattput (ncid, s09_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000971;
   ncattput (ncid, s09_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s09_valid_range[0] = 0;
   s09_valid_range[1] = 100;
   ncattput (ncid, s09_id, "valid_range", NC_FLOAT, 2, (void *) s09_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s09_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s09_id, "LAPS_var", NC_CHAR, 3, (void *)"S09");
   ncattput (ncid, s09_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s09_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s10_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 7.43 micron sounding radiance");
   ncattput (ncid, s10_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000743;
   ncattput (ncid, s10_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s10_valid_range[0] = 0;
   s10_valid_range[1] = 100;
   ncattput (ncid, s10_id, "valid_range", NC_FLOAT, 2, (void *) s10_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s10_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s10_id, "LAPS_var", NC_CHAR, 3, (void *)"S10");
   ncattput (ncid, s10_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s10_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s11_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 7.02 micron sounding radiance");
   ncattput (ncid, s11_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000702;
   ncattput (ncid, s11_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s11_valid_range[0] = 0;
   s11_valid_range[1] = 100;
   ncattput (ncid, s11_id, "valid_range", NC_FLOAT, 2, (void *) s11_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s11_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s11_id, "LAPS_var", NC_CHAR, 3, (void *)"S11");
   ncattput (ncid, s11_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s11_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s12_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 6.51 micron sounding radiance");
   ncattput (ncid, s12_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000651;
   ncattput (ncid, s12_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s12_valid_range[0] = 0;
   s12_valid_range[1] = 100;
   ncattput (ncid, s12_id, "valid_range", NC_FLOAT, 2, (void *) s12_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s12_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s12_id, "LAPS_var", NC_CHAR, 3, (void *)"S12");
   ncattput (ncid, s12_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s12_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s13_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 4.57 micron sounding radiance");
   ncattput (ncid, s13_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000457;
   ncattput (ncid, s13_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s13_valid_range[0] = 0;
   s13_valid_range[1] = 100;
   ncattput (ncid, s13_id, "valid_range", NC_FLOAT, 2, (void *) s13_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s13_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s13_id, "LAPS_var", NC_CHAR, 3, (void *)"S13");
   ncattput (ncid, s13_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s13_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s14_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 4.52 micron sounding radiance");
   ncattput (ncid, s14_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000452;
   ncattput (ncid, s14_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s14_valid_range[0] = 0;
   s14_valid_range[1] = 100;
   ncattput (ncid, s14_id, "valid_range", NC_FLOAT, 2, (void *) s14_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s14_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s14_id, "LAPS_var", NC_CHAR, 3, (void *)"S14");
   ncattput (ncid, s14_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s14_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s15_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 4.45 micron sounding radiance");
   ncattput (ncid, s15_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000445;
   ncattput (ncid, s15_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s15_valid_range[0] = 0;
   s15_valid_range[1] = 100;
   ncattput (ncid, s15_id, "valid_range", NC_FLOAT, 2, (void *) s15_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s15_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s15_id, "LAPS_var", NC_CHAR, 3, (void *)"S15");
   ncattput (ncid, s15_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s15_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s16_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 4.13 micron sounding radiance");
   ncattput (ncid, s16_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000413;
   ncattput (ncid, s16_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s16_valid_range[0] = 0;
   s16_valid_range[1] = 100;
   ncattput (ncid, s16_id, "valid_range", NC_FLOAT, 2, (void *) s16_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s16_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s16_id, "LAPS_var", NC_CHAR, 3, (void *)"S16");
   ncattput (ncid, s16_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s16_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s17_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 3.98 micron sounding radiance");
   ncattput (ncid, s17_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000398;
   ncattput (ncid, s17_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s17_valid_range[0] = 0;
   s17_valid_range[1] = 100;
   ncattput (ncid, s17_id, "valid_range", NC_FLOAT, 2, (void *) s17_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s17_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s17_id, "LAPS_var", NC_CHAR, 3, (void *)"S17");
   ncattput (ncid, s17_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s17_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s18_id, "long_name", NC_CHAR, 36, (void *)"GOES-8 3.74 micron sounding radiance");
   ncattput (ncid, s18_id, "units", NC_CHAR, 29, (void *)"W/(m*-2 sec sterradian cm*-1)");
   double_val = 0.000374;
   ncattput (ncid, s18_id, "wavelen_in_cm", NC_DOUBLE, 1,(void *) &double_val);
   s18_valid_range[0] = 0;
   s18_valid_range[1] = 100;
   ncattput (ncid, s18_id, "valid_range", NC_FLOAT, 2, (void *) s18_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s18_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s18_id, "LAPS_var", NC_CHAR, 3, (void *)"S18");
   ncattput (ncid, s18_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s18_id, "LAPS_units", NC_CHAR, 10, (void *)"units attr");
   ncattput (ncid, s19_id, "long_name", NC_CHAR, 33, (void *)"GOES-8 0.67 micron visible counts");
   ncattput (ncid, s19_id, "units", NC_CHAR, 6, (void *)"counts");
   s19_valid_range[0] = 0;
   s19_valid_range[1] = 8191;
   ncattput (ncid, s19_id, "valid_range", NC_FLOAT, 2, (void *) s19_valid_range);
   float_val = 9.9999999e+36;
   ncattput (ncid, s19_id, "_FillValue", NC_FLOAT, 1,(void *) &float_val);
   ncattput (ncid, s19_id, "LAPS_var", NC_CHAR, 3, (void *)"S19");
   ncattput (ncid, s19_id, "lvl_coord", NC_CHAR, 4, (void *)"NONE");
   ncattput (ncid, s19_id, "LAPS_units", NC_CHAR, 6, (void *)"COUNTS");
   ncattput (ncid, fctimes_id, "long_name", NC_CHAR, 14, (void *)"forecast times");
   ncattput (ncid, level_id, "long_name", NC_CHAR, 5, (void *)"level");
   short_val = 0;
   ncattput (ncid, s01_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s02_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s03_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s04_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s05_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s06_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s07_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s08_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s09_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s10_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s11_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s12_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s13_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s14_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s15_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s16_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s17_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s18_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);
   short_val = 0;
   ncattput (ncid, s19_fcinv_id, "_FillValue", NC_SHORT, 1,(void *) &short_val);

   /* leave define mode */
   ncendef (ncid);

   {			/* store version */
    static long version = {2};
    ncvarput1(ncid, version_id, (long *)0, (void *)&version);
   }

   {			/* store fctimes */
    static long fctimes_start[] = {0};
    static long fctimes_edges[] = {1};
    static short fctimes[] = {0};
    ncvarput(ncid, fctimes_id, fctimes_start, fctimes_edges, (void *)fctimes);
   }

   {			/* store num_variables */
    static long num_variables = {70};
    ncvarput1(ncid, num_variables_id, (long *)0, (void *)&num_variables);
   }

   {			/* store level */
    static long level_start[] = {0};
    static long level_edges[] = {1};
    static short level[] = {0};
    ncvarput(ncid, level_id, level_start, level_edges, (void *)level);
   }

   {			/* store origin */
    static long origin_start[] = {0};
    static long origin_edges[] = {132};
    static char origin[] = {"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO"};
    ncvarput(ncid, origin_id, origin_start, origin_edges, (void *)origin);
   }

   {			/* store model */
    static long model_start[] = {0};
    static long model_edges[] = {132};
    static char model[] = {"LAPS - Local Analysis and Prediction System"};
    ncvarput(ncid, model_id, model_start, model_edges, (void *)model);
   }

   {			/* store imax */
    static long imax = {NX};
    ncvarput1(ncid, imax_id, (long *)0, (void *)&imax);
   }

   {			/* store jmax */
    static long jmax = {NY};
    ncvarput1(ncid, jmax_id, (long *)0, (void *)&jmax);
   }

   {			/* store kmax */
    static long kmax = {19};
    ncvarput1(ncid, kmax_id, (long *)0, (void *)&kmax);
   }

   {			/* store kdim */
    static long kdim = {19};
    ncvarput1(ncid, kdim_id, (long *)0, (void *)&kdim);
   }
   return ncid;
}
/************************************************************************/
#ifdef __STDC__
void log_diag(int ii, char *cstring, char *name)
#else
void log_diag( ii, cstring, name )
int ii;
char *cstring, *name;
#endif
{
        return;
}

/*************************************************************************
* itoa
*
* The purpose of this function is to convert an integer into a null 
* terminated ascii string of maximum length n.
*	
* Programmer: James Fluke
*************************************************************************/
#ifdef __STDC__
char *itoa(int ii,register char *str,register int n)
#else
char *itoa (ii, str, n)
int ii;
register char *str;
register int n;
#endif
{
	register int div;
  
	if (n <= 0) return(NULL);
	if (n <= 1) {
		*str = '\0';
		return(NULL);
	}
	if ( ii < 0) {
		ii = -ii;
		*str = '-';
		++str;
		if (--n <= 1) {
			*str = '\0';
			return(NULL);
		}
	}
	div = 10;
	while ( div <= ii ) div *= 10;
	do {
		div /= 10;
		*str = ( ii/div ) + '0';
		++str;		ii = ii % div;
	} while ( div > 1 && --n > 1);
	*str = '\0';
	if (div > 1) return(NULL);
	return (str);
}

/*************************************************************************
*	FILL_EMPTY_GRIDS
*	Category	Product Management
*	Group		General Purpose Database
*	Module		Retrieve Grid Data
*	Purpose		Packs data array with 1e+37 at the levels where
*			data could not be retrieved
*
*	Designer/Programmer : Linda Wharton
*	Modifications : original 1/93
*
*	Input:
*		iimax      	Expected X dimension of data
*		jjmax      	Expected Y dimension of data
*		kkmax      	Expected number of data grids to retrieve
*	Output:
*		data		3D grid containing data requested
*		process_var	Array containing 1 if variable was processed,
*				  otherwise contains 0
*	Globals:
*		none
*	Returns:
*		none
*************************************************************************/
#ifdef __STDC__
#if defined(__alpha) || defined(__IP21)
void fill_empty_grids(int *iimax, int *jjmax, int *kkmax,
		      int *process_var,float *data)
#else
void fill_empty_grids(long *iimax, long *jjmax, long *kkmax,
		      long *process_var,float *data)
#endif
#else                 
#if defined(__alpha) || defined(__IP21)
void fill_empty_grids(iimax,jjmax,kkmax,process_var,data)
int *iimax; 
int *jjmax; 
int *kkmax;
int *process_var;
float *data;
#else
void fill_empty_grids(iimax,jjmax,kkmax,process_var,data)
long *iimax; 
long *jjmax; 
long *kkmax;
long *process_var;
float *data;
#endif
#endif
{
	int i, j, k;
	
	for (k = 0; k < *kkmax; k++) {
	   if (process_var[k] == 0) {
	      for (j = 0; j < *jjmax; j++) {
	         for (i = 0; i < *iimax; i++)
	           *(data + (k*(*iimax)*(*jjmax)) + (j*(*iimax)) + i) = 1e+37;
	      }
	   }
	}
	return;
}

/*************************************************************************
*	CDF_GET_COORD
*	Category	Product Management
*	Group		General Purpose Database
*	Module		Read variable from netCDF file
*	Purpose		To read the contents of a coordinate variable in an 
*				open netCDF file.
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 8/22/90
*			modified for READLAPSDATA 1/93 Linda Wharton
*
*	Input:
*		i_cdfid		The file descriptor of an open netCDF file
*		v_name		The name of the variable
*		i_size		The number of elements in the variable
*		i_cdftype	The type of variable, as specified by netCDF:
*					NC_BYTE, NC_CHAR, NC_SHORT, NC_LONG, NC_FLOAT,
*					NC_DOUBLE.
*	Output:
*		none
*	Globals:
*		none
*	Returns:
*		pointer to the values read from the netCDF file
*		-1 if there is an error
*************************************************************************/
#ifdef __STDC__
#if defined(__alpha) || defined(__IP21)
char *cdf_get_coord (int i_cdfid, char *v_name, long i_size, int i_cdftype)
#else
char *cdf_get_coord (int i_cdfid, char *v_name, int i_size, int i_cdftype)
#endif
#else
#if defined(__alpha) || defined(__IP21)
char *cdf_get_coord (i_cdfid, v_name, i_size, i_cdftype)
int i_cdfid;
char *v_name;
long i_size;
int i_cdftype;
#else
char *cdf_get_coord (i_cdfid, v_name, i_size, i_cdftype)
int i_cdfid;
char *v_name;
int i_size;
int i_cdftype;
#endif
#endif

{

	static char *v_ptr;
#ifdef CRAY
	static int nbytes[]={0,sizeof(ncbyte),sizeof(ncchar),sizeof(ncshort),sizeof(nclong),sizeof(ncfloat),sizeof(ncdouble)};
	int i_status, i_varid, i_start;
#else
#if defined(__alpha) || defined(__IP21)
	static int nbytes[]={0,1,1,2,8,4,8};
	int i_status, i_varid;
	long i_start;
#else
	static int nbytes[]={0,1,1,2,4,4,8};
	int i_status, i_varid, i_start;
#endif
#endif

	log_diag (2, "cdf_get_var: v_name = %s\n", v_name);

/* allocate memory to hold the values associated with the dimension */
	v_ptr = malloc (i_size * nbytes[i_cdftype]); 

/* get the variable id for this dimension */
	i_varid = ncvarid (i_cdfid, v_name);

	log_diag (2, "cdf_get_var: Size = %d   Var id = %d\n", i_size, i_varid);

/* read the contents of the variable into memory */
	i_start = 0;
	i_status = ncvarget (i_cdfid, i_varid, &i_start, &i_size, (void *)v_ptr);
	if (i_status == (-1)){
		printf("error reading values for array index %s\n", v_name);
		return (char *)i_status;
	}
	
	return v_ptr;
}

/************************************************************************/
#ifdef __STDC__
int cdf_get_levels (int i_cdfid, char *v_name, int i_size, short *lvls)
#else
int cdf_get_levels (i_cdfid,v_name,i_size,lvls)
int i_cdfid; 
char *v_name; 
int i_size; 
short *lvls;
#endif
{
	int i_status, i_varid, i_start;

	log_diag (2, "cdf_get_var: v_name = %s\n", v_name);

/* get the variable id for this dimension */
	i_varid = ncvarid (i_cdfid, v_name);

	log_diag (2, "cdf_get_var: Size = %d   Var id = %d\n", i_size, i_varid);

/* read the contents of the variable into memory */
	i_start = 0;
	i_status = ncvarget (i_cdfid, i_varid, &i_start, &i_size,(void *)lvls);
	if (i_status == (-1)){
		printf("error reading values for array index %s\n", v_name);
		return i_status;
	}
	else
	  return 0;
}

/*************************************************************************
*	CDF_GET_INDEX
*	Category	Product Management
*	Group		General Purpose Database
*	Module		Get indexes into variable arrays
*	Purpose		Determine if a specified value exists in coordinate array
*				in an open netCDF file.
*	Usage		Call this routine before reading grid data from the netCDF
*				file.
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 8/22/90
*
*	Input:
*		i_cdfid		The file descriptor of an open netCDF file
*		i_value		The value to search for
*		d_name		The name of the coordinate array to be searched
*	Output:
*		none
*	Globals:
*		none
*	Returns
*	   	 the location of the value in the array, if it is found
*		-1 if there is an error, or if the value is not found		
*************************************************************************/
#ifdef __STDC__
int cdf_get_index (int i_cdfid, int i_value, char *d_name)
#else
int cdf_get_index (i_cdfid, i_value, d_name)
int i_cdfid, i_value;
char *d_name;
#endif
{
	int i, i_dimid, i_status, i_dsize, i_varid, i_start;
	char *d_ptr;
	short *i_ptr;
	char *cdf_get_coord();

	log_diag (2, "cdf_get_index: d_name = %s\n", d_name);
	
/* read in the dimension id and size from the netcdf file */
	if ((i_dsize = cdf_dim_size (i_cdfid, d_name)) == (-1)){
		return -1;
	}

	log_diag (2, "cdf_get_index: dim size = %d\n", i_dsize);

/* read the contents of the coordinate variable associated with this dimension
   from the net_cdf file */
	if ((d_ptr = cdf_get_coord(i_cdfid,d_name,i_dsize,NC_SHORT))==(char *)(-1))
		return -1; 
	i_ptr = (short *)d_ptr;

/* locate the value in the array pointed to by i_ptr */
	for (i=0; ((i<i_dsize) && (i_value!=*(i_ptr+i))); i++);

/* deallocate the memory */
	free (d_ptr);

/* test to see if the value was in the array */
	if (i<i_dsize)
		return i;
	else 	
		return -1;
}

/*************************************************************************
*  CDF_DIM_SIZE
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Read dimension from NetCDF file
*	Purpose			To read the size of a dimension from an open netCDF 
*						file.
*
*	Designer/Programmer : MarySue Schultz
*	Modifications :	original 9/19/90
*
*   Input :
*		i_cdfid		the file id of an open netCDF file
*		d_name		the name of the dimension
*	Output :
*		none
*	Globals:
*		none
*	Returns:
*		dimension size
*		-1 if an error occurs
**************************************************************************/
#ifdef __STDC__
int cdf_dim_size (int i_cdfid, char *d_name)
#else
int cdf_dim_size (i_cdfid, d_name)
int i_cdfid;
char *d_name;
#endif
{
	int i, i_dimid, i_status;
#if defined(__alpha) || defined(__IP21)
        long i_dsize;
#else
        int i_dsize;
#endif

/* read in the dimension id and size */
	if ((i_dimid = ncdimid (i_cdfid, d_name)) == (-1))
		return -1;

	if ((i_status = ncdiminq(i_cdfid,i_dimid,(char *)0,&i_dsize)) == (-1))
		return -1;

	return (int)i_dsize;
}

/******************************************************************************
* CDF_CHK_LAPS_INV
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Get queue location
*	Purpose			To find the location in an i4_times queue of a specified
*						product time.
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 12/31/90
*			modified for new READLAPSDATA 1/93 Linda Wharton
*
*	Input :
*		i_cdfid		file id of an open netCDF file
*		i_lindx		level index. obtained by calling cdf_indx.
*		i_fcindx	forecast time index. obtained by calling cdf_indx.
*		g_name		A character string identifying the field
*	Output :
*		None
*	Globals :
*		None
*	Returns :
*		 1 if grid is available
*		-1 if an error occurs, or grid is not available
*******************************************************************************/
#ifdef __STDC__
short cdf_check_laps_inv (int i_cdfid, int i_lindx, int i_fcindx, 
			  char *g_name)
#else
short cdf_check_laps_inv (i_cdfid, i_lindx, i_fcindx, g_name)
int i_cdfid, i_lindx, i_fcindx;
char *g_name;
#endif
{
	int 	*i4_ptr, i, i_indx, i_status, i_varid;
        long    start[2];
	short   i_flag;
	char 	dim_name[15], var_name[15];
	void	cdf_i4times();

/* read the fctimes inventory array associated with this grid */
	sprintf (var_name, "%s%s", g_name, "_fcinv");
	log_diag (2, "cdf_chk_ninv: var_name = %s\n", var_name);

	i_varid = ncvarid (i_cdfid, var_name);

	start[0] = i_fcindx;
	start[1] = i_lindx;

	i_status = ncvarget1 (i_cdfid, i_varid, start, (void *)&i_flag);
	if (i_status == -1)
		return -1;

	log_diag (2, "cdf_chk_ninv: i_flag = %d\n", i_flag);

	if (i_flag == 1)
		return 1;
	else 
		return -1;
}

/***************************************************************************
* CDF_WRITE_GRID
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Write grid into netCDF file
*	Purpose			To write a grid (including product and data headers)
*						into an open netCDF file.
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 9/5/90
*
*	Input :
*		i_cdfid 	The file id of a netCDF file opened for writing
*		gdims		A structure containing information needed to locate the
*						data in the appropriate data array : variable id, 
*						level coordinate, forecast time coordinate, position
*						in the queue of i4times, x dimension size, y dimension
*						size.  
*		g_name		A character string identifying the field:
*						h 		height
*						li		lifted index
*						p		pressure
*						pcp		precipitation
*						rh		relative humidity
*						t		temperature
*						u		east wind component
*						v		north wind component
*						w		vertical velocity
*		gpointr		A pointer to the grid data	
*	Output :
*		None
*	Globals :
*		None
*	Returns :
*		 0 if successful
*		-1 if an error occurs
**************************************************************************/
#ifdef __STDC__
int cdf_write_grid (int i_cdfid, cdf_grid_info *gdims,char *gname, 
		    float *gptr, char *cname, char *cptr)
#else
int cdf_write_grid (i_cdfid, gdims, gname, gptr, cname, cptr)
int i_cdfid;
cdf_grid_info *gdims;
char *gname; 
float *gptr;
char *cname;
char *cptr;
#endif
{
	int i_status, i_comid,i_varid; 
#if defined(__alpha) || defined(__IP21)
        long start[4], count[4], start_c[3], count_c[3];
#else
        int start[4], count[4], start_c[3], count_c[3];
#endif

/* get the variable ids */
	log_diag (2, "Data variable name = %s\n", gname);
	if ((i_varid = ncvarid (i_cdfid, gname)) == (-1))
		return -1;
	if ((i_comid = ncvarid (i_cdfid, cname)) == (-1))
		return -1;

/* construct the arrays needed by the netcdf write routine */
	start[0] = gdims->fctime_coord;
	start[1] = gdims->level_coord;
	start[2] = 0;
	start[3] = 0;

	count[0] = 1;
	count[1] = 1;
	count[2] = gdims->y_dim;
	count[3] = gdims->x_dim;

	start_c[0] = gdims->fctime_coord;
	start_c[1] = gdims->level_coord;
	start_c[2] = 0;
	
	count_c[0] = 1;
	count_c[1] = 1;
	count_c[2] = strlen(cptr);

/* write the grid to the netcdf file */
	i_status = ncvarput (i_cdfid, i_varid, start, count, gptr);
	
	if (i_status == (-1))
		return -1;
	else {
	  i_status = ncvarput (i_cdfid, i_comid, start_c, count_c, cptr);
	  
	  if (i_status == (-1))
		return -1;
	  else
		return 0;
	}
}

/***************************************************************************
* CDF_UPDATE_LAPS_INV
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Update LAPS inventory in netCDF file
*	Purpose			To write a flag indicating that a particular 
*				  LAPS grid has been written into the netCDF 
*				  file
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 9/5/90
*			adapted for new WRITELAPSDATA 1/93 Linda Wharton
*
*	Input :
*		i_cdfid 	The file id of a netCDF file opened for writing
*		gdims		A structure containing information needed to 
*				  locate the data in the appropriate data 
*				  array : variable id, level coordinate, 
*				  forecast time coordinate, position in the 
*				  queue of i4times, x dimension size, 
*				  y dimension size.  
*		g_type		Type of grid: pr = pressure (millibars), 
*				  agl = above ground level, sfc = surface.
*		g_name		A character string identifying the field:
*				    w		wind omega
*				    u		Eastward wind 
*				    v		Northward wind
*	Output :
*		None
*	Globals :
*		None
*	Returns :
*		 0 if successful
*		-1 if an error occurs
**************************************************************************/		
#ifdef __STDC__
int cdf_update_laps_inv (int i_cdfid,cdf_grid_info *gdims, char *gname)
#else
int cdf_update_laps_inv (i_cdfid, gdims, gname)
int i_cdfid;
cdf_grid_info *gdims;
char *gname;
#endif
{
	static short i_flag=1;
	int i_status, i_varid;
	char varname[15];
#if defined(__alpha) || defined(__IP21)
        long start[2];
#else
        int start[2];
#endif

/* print out all of the dimension data */
	log_diag (2, "cdf_upd_linv: fcindx = %d   lindx = %d   x = %d   y = %d\n",
		gdims->fctime_coord, gdims->level_coord, gdims->x_dim, gdims->y_dim);

/* get the variable id of the fctimes inventory array */
	sprintf (varname, "%s%s", gname, "_fcinv");
	log_diag (2, "cdf_upd_linv: inventory variable name = %s\n", varname);
	if ((i_varid = ncvarid (i_cdfid, varname)) == (-1))
		return -1;
	else
		log_diag (2, "cdf_upd_linv: varid = %d\n", i_varid);

/* construct the arrays needed by the netcdf write routine */
	start[0] = gdims->fctime_coord;
	start[1] = gdims->level_coord;

	i_status = ncvarput1 (i_cdfid, i_varid, start, (void *)&i_flag); 
	log_diag (2, "cdf_upd_linv: ncvarput1 status = %d\n", i_status);

	if (i_status == (-1)){
		return -1;
	}
	else
		return 0;
}

/*************************************************************************
*	MAKE_C_FNAME
*	Category	Product Management
*	Group		General Purpose Database
*	Module		Retrieve/Write Grid data
*	Purpose		Create a filename for passing into C functions.
*
*	Designer/Programmer : Linda Wharton
*	Modifications : original 1/93
*
*	Input:
*		dir		Directory file is located in
*		gtime		Time of the filename
*		ext		Extension of the filename
*	Output:
*		fname		Filename to open
*	Globals:
*		none
*	Returns:
*		none
*************************************************************************/
#ifdef __STDC__
void make_c_fname(char *filname,short *s_length, char fname[])
#else
void make_c_fname(filname,s_length,fname)
char *filname;
short *s_length;
char fname[];
#endif
{
      
	int not_end, i;
      
/* receive pointer to file descriptor filname.  Convert that
   to character string fname */
   
   	strncpy(fname,"\0",1);
   	not_end = 1;  /* not end of characters is true */
   	for (i = 0; i < *s_length; i++) 
   	   if ((isgraph(*(filname + i))) && not_end)
   	      strncat(fname,(filname + i),1);
   	   else {
   	      strncat(fname,"\0",1);
   	      not_end = 0;
   	   }
      
      return;
}      

/*************************************************************************
*	GET_CDF_VAR
*	Category	Product Management
*	Group		General Purpose Database
*	Module		Retrieve Grid Data
*	Purpose		Converts LAPS variables to netCDF variables, returns
*			forecast times of the LAPS variables, returns
*			process_var with zero if var_req cannot be converted.
*
*	Designer/Programmer : Linda Wharton
*	Modifications : original 1/93
*
*	Input:
*		ext		Extension identifying data to be read
*		kkmax      	Expected number of data grids to retrieve
*		var_req		Array of LAPS variables to retrieve from the file
*	Output:
*		cdf_var		Array of netCDF variables that correspond to
*				LAPS variables in the VAR_REQ array
*		process_var	Array containing 1 if variable was processed,
*				  otherwise contains 0
*		cdf_fctime	Array containing the forecast time for
*				LAPS variables in the VAR_REQ array
*	Globals:
*		none
*	Returns:
*		count of requested variables that were not converted
*************************************************************************/
#ifdef __STDC__
#if defined(__alpha) || defined(__IP21)
int get_cdf_var(char *ext,int *kkmax,char *var_req,char *cdf_var,
		int *process_var,short cdf_fctime[])
#else
int get_cdf_var(char *ext,long *kkmax,char *var_req,char *cdf_var,
		long *process_var,short cdf_fctime[])
#endif
#else
#if defined(__alpha) || defined(__IP21)
int get_cdf_var(ext,kkmax,var_req,cdf_var,process_var,cdf_fctime)
char *ext;
int *kkmax;
char *var_req;
char *cdf_var;
int *process_var;
short cdf_fctime[];
#else
int get_cdf_var(ext,kkmax,var_req,cdf_var,process_var,cdf_fctime)
char *ext;
long *kkmax;
char *var_req;
char *cdf_var;
long *process_var;
short cdf_fctime[];
#endif
#endif
{
	enum grid {LW3, LH1, LH2, LH3, LH4, LQ3, LSX, LWM, LT1, LHE,
		   LIW, LMT, LMR, LF1, L1S, LPS, LRP, LBA, LC3, LWC,
		   LIL, LCB, LCT, LCV, LMD, LVE, LCO, LTY, LCP, LVD, 
		   LMA, LMF, Z02, RAM, RSF, RSM, VRC, LM1, LM2, VRD, 
                   V00, LGA, LGF, LN3, MM5, MSF, SST, VEG, LS8}; 
	typedef char STRING5[5];
	int count, i;
	enum grid grid_no;
        char *varptr;
	
	if (strncmp(ext,"LW3",3) == 0) grid_no = LW3;
	if (strncmp(ext,"LH1",3) == 0) grid_no = LH1;
	if (strncmp(ext,"LH2",3) == 0) grid_no = LH2;
	if (strncmp(ext,"LH3",3) == 0) grid_no = LH3;
	if (strncmp(ext,"LH4",3) == 0) grid_no = LH4;
	if (strncmp(ext,"LQ3",3) == 0) grid_no = LQ3;
	if (strncmp(ext,"LSX",3) == 0) grid_no = LSX;
	if (strncmp(ext,"LWM",3) == 0) grid_no = LWM;
	if (strncmp(ext,"LT1",3) == 0) grid_no = LT1;
	if (strncmp(ext,"LHE",3) == 0) grid_no = LHE;
	if (strncmp(ext,"LIW",3) == 0) grid_no = LIW;
	if (strncmp(ext,"LMT",3) == 0) grid_no = LMT;
	if (strncmp(ext,"LMR",3) == 0) grid_no = LMR;
	if (strncmp(ext,"LF1",3) == 0) grid_no = LF1;
	if (strncmp(ext,"L1S",3) == 0) grid_no = L1S;
	if (strncmp(ext,"LPS",3) == 0) grid_no = LPS;
	if (strncmp(ext,"LRP",3) == 0) grid_no = LRP;
	if (strncmp(ext,"LBA",3) == 0) grid_no = LBA;
	if (strncmp(ext,"LC3",3) == 0) grid_no = LC3;
	if (strncmp(ext,"LWC",3) == 0) grid_no = LWC;
	if (strncmp(ext,"LIL",3) == 0) grid_no = LIL;
	if (strncmp(ext,"LCB",3) == 0) grid_no = LCB;
	if (strncmp(ext,"LCT",3) == 0) grid_no = LCT;
	if (strncmp(ext,"LCV",3) == 0) grid_no = LCV;
	if (strncmp(ext,"LMD",3) == 0) grid_no = LMD;
	if (strncmp(ext,"LCO",3) == 0) grid_no = LCO;
	if (strncmp(ext,"LTY",3) == 0) grid_no = LTY;
	if (strncmp(ext,"LCP",3) == 0) grid_no = LCP;
	if (strncmp(ext,"LVD",3) == 0) grid_no = LVD;
	if (strncmp(ext,"LVE",3) == 0) grid_no = LVE;
	if (strncmp(ext,"LMA",3) == 0) grid_no = LMA;
	if (strncmp(ext,"LMF",3) == 0) grid_no = LMF;
	if (strncmp(ext,"Z02",3) == 0) grid_no = Z02;
	if (strncmp(ext,"RAM",3) == 0) grid_no = RAM;
	if (strncmp(ext,"RSF",3) == 0) grid_no = RSF;
	if (strncmp(ext,"RSM",3) == 0) grid_no = RSM;
	if (strncmp(ext,"VRC",3) == 0) grid_no = VRC;
	if (strncmp(ext,"LM1",3) == 0) grid_no = LM1;
	if (strncmp(ext,"LM2",3) == 0) grid_no = LM2;
	if (strncmp(ext,"VRD",3) == 0) grid_no = VRD;
	if (strncmp(ext,"V01",3) == 0) grid_no = V00;
	if (strncmp(ext,"V02",3) == 0) grid_no = V00;
	if (strncmp(ext,"V03",3) == 0) grid_no = V00;
	if (strncmp(ext,"V04",3) == 0) grid_no = V00;
	if (strncmp(ext,"V05",3) == 0) grid_no = V00;
	if (strncmp(ext,"V06",3) == 0) grid_no = V00;
	if (strncmp(ext,"V07",3) == 0) grid_no = V00;
	if (strncmp(ext,"V08",3) == 0) grid_no = V00;
	if (strncmp(ext,"V09",3) == 0) grid_no = V00;
	if (strncmp(ext,"V10",3) == 0) grid_no = V00;
	if (strncmp(ext,"V11",3) == 0) grid_no = V00;
	if (strncmp(ext,"V12",3) == 0) grid_no = V00;
	if (strncmp(ext,"V13",3) == 0) grid_no = V00;
	if (strncmp(ext,"V14",3) == 0) grid_no = V00;
	if (strncmp(ext,"V15",3) == 0) grid_no = V00;
	if (strncmp(ext,"V16",3) == 0) grid_no = V00;
	if (strncmp(ext,"V17",3) == 0) grid_no = V00;
	if (strncmp(ext,"V18",3) == 0) grid_no = V00;
	if (strncmp(ext,"V19",3) == 0) grid_no = V00;
	if (strncmp(ext,"V20",3) == 0) grid_no = V00;
	if (strncmp(ext,"LGA",3) == 0) grid_no = LGA;
	if (strncmp(ext,"LGF",3) == 0) grid_no = LGF;
	if (strncmp(ext,"LN3",3) == 0) grid_no = LN3;
	if (strncmp(ext,"MM5",3) == 0) grid_no = MM5;
	if (strncmp(ext,"MSF",3) == 0) grid_no = MSF;
	if (strncmp(ext,"SST",3) == 0) grid_no = SST;
	if (strncmp(ext,"VEG",3) == 0) grid_no = VEG;
	if (strncmp(ext,"LS8",3) == 0) grid_no = LS8;

	count = 0;	
	switch (grid_no)  {
		
	case LW3:
	   /* cdf_fctime grid left at zero, only analysis in LW3 */
           for (i = 0; i < *kkmax; i++) {
              if ((strncmp((var_req + i*4), "U",3) ==  0) ||
                  (strncmp((var_req + i*4), "U3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"u");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "V",3) ==  0) ||
                       (strncmp((var_req + i*4), "V3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"v");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "OM",3) ==  0) {
                 strcpy((cdf_var + i*5),"w");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LW3  */

	case LH1:
	   /* cdf_fctime grid left at zero, only analysis in LH1 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "PW",3) ==  0) {
                 strcpy((cdf_var + i*5),"pw");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LH1  */
	
	case LH2:
	   /* cdf_fctime grid left at zero, only analysis in LH2 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "PW1",3) ==  0) {
                 strcpy((cdf_var + i*5),"lpw");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PW2",3) ==  0) {
                 strcpy((cdf_var + i*5),"lpw");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PW3",3) ==  0) {
                 strcpy((cdf_var + i*5),"lpw");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PW",3) ==  0) {
                 strcpy((cdf_var + i*5),"lpw");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LH2  */
	
	case LH3:
	   /* cdf_fctime grid left at zero, only analysis in LH3 */
           for (i = 0; i < *kkmax; i++) {
              if ((strncmp((var_req + i*4), "RH",3) ==  0) ||
                  (strncmp((var_req + i*4), "RH3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"rh");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RHL",3) ==  0) {
                 strcpy((cdf_var + i*5),"rhl");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LH3  */
	
	case LH4:
	   /* cdf_fctime grid left at zero, only analysis in LH4 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "TPW",3) ==  0) {
                 strcpy((cdf_var + i*5),"tpw");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LH4  */
	
	case LQ3:
	   /* cdf_fctime grid left at zero, only analysis in LQ3 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "SH",3) ==  0) {
                 strcpy((cdf_var + i*5),"sh");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LQ3  */
        case LSX:
	   /* cdf_fctime grid left at zero, only analysis in LSX */
           for (i = 0; i < *kkmax; i++) {
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "U",3) ==  0)) {
                 strcpy((cdf_var + i*5),"su");
                 process_var[i] = 1;
              }                                        
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "V",3) ==  0)) {
                 strcpy((cdf_var + i*5),"sv");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "P",3) ==  0)) {
                 strcpy((cdf_var + i*5),"fp");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "T",3) ==  0)) {
                 strcpy((cdf_var + i*5),"st");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "TD",3) ==  0)) {
                 strcpy((cdf_var + i*5),"std");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "VV",3) ==  0)) {
                 strcpy((cdf_var + i*5),"vv");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "RH",3) ==  0)) {
                 strcpy((cdf_var + i*5),"srh");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "CCE",3) ==  0)) {
                 strcpy((cdf_var + i*5),"ccg");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "MSL",3) ==  0)) {
                 strcpy((cdf_var + i*5),"mp");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "TAD",3) ==  0)) {
                 strcpy((cdf_var + i*5),"ta");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "TH",3) ==  0)) {
                 strcpy((cdf_var + i*5),"pot");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "THE",3) ==  0)) {
                 strcpy((cdf_var + i*5),"ept");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "PS",3) ==  0)) {
                 strcpy((cdf_var + i*5),"sp");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "VOR",3) ==  0)) {
                 strcpy((cdf_var + i*5),"vor");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "MR",3) ==  0)) {
                 strcpy((cdf_var + i*5),"mr");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "MRC",3) ==  0)) {
                 strcpy((cdf_var + i*5),"mc");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "DIV",3) ==  0)) {
                 strcpy((cdf_var + i*5),"d");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "THA",3) ==  0)) {
                 strcpy((cdf_var + i*5),"pta");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "MRA",3) ==  0)) {
                 strcpy((cdf_var + i*5),"ma");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "LI",3) ==  0)) {
                 strcpy((cdf_var + i*5),"li");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "SPD",3) ==  0)) {
                 strcpy((cdf_var + i*5),"spd");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "CSS",3) ==  0)) {
                 strcpy((cdf_var + i*5),"cssi");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "PBE",3) ==  0)) {
                 strcpy((cdf_var + i*5),"pbe");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "NBE",3) ==  0)) {
                 strcpy((cdf_var + i*5),"nbe");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "VIS",3) ==  0)) {
                 strcpy((cdf_var + i*5),"vis");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "FWX",3) ==  0)) {
                 strcpy((cdf_var + i*5),"fwx");
                 process_var[i] = 1;
              }
              if ((process_var[i] == 0) && 
                 (strncmp((var_req + i*4), "HI",3) ==  0)) {
                 strcpy((cdf_var + i*5),"hi");
                 process_var[i] = 1;
              }
              if (process_var[i] == 0) {
                 strcpy((cdf_var + i*5),"\0");
              	 count = count + 1;
              }
           }
	   break;  /*  case LSX  */
	case LWM:
	   /* cdf_fctime grid left at zero, only analysis in LWM */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "SU",3) ==  0) {
                 strcpy((cdf_var + i*5),"u");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SV",3) ==  0) {
                 strcpy((cdf_var + i*5),"v");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LWM  */
	case LT1:
	   /* cdf_fctime grid left at zero, only analysis in LT1 */
           for (i = 0; i < *kkmax; i++) {
              if ((strncmp((var_req + i*4), "T",3) ==  0) ||
                  (strncmp((var_req + i*4), "T3",3) == 0)) {
                 strcpy((cdf_var + i*5),"t");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "HT",3) ==  0) {
                 strcpy((cdf_var + i*5),"z");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LT1  */
	case LHE:
	   /* cdf_fctime grid left at zero, only analysis in LHE */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LHE",3) ==  0) {
                 strcpy((cdf_var + i*5),"hel");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "MU",3) ==  0) {
                 strcpy((cdf_var + i*5),"mu");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "MV",3) ==  0) {
                 strcpy((cdf_var + i*5),"mv");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LHE  */
	case LIW:
	   /* cdf_fctime grid left at zero, only analysis in LIW */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LIW",3) ==  0) {
                 strcpy((cdf_var + i*5),"liw");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "W",3) ==  0) {
                 strcpy((cdf_var + i*5),"w");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LIW  */
	case LMT:
	   /* cdf_fctime grid left at zero, only analysis in LMT */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LMT",3) ==  0) {
                 strcpy((cdf_var + i*5),"etop");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LLR",3) ==  0) {
                 strcpy((cdf_var + i*5),"llr");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LMT  */
	case LMR:
	   /* cdf_fctime grid adjusted for analysis, 1 hr and 2 hr in LMR */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "R00",3) ==  0) {
                 strcpy((cdf_var + i*5),"mxrf");
                 cdf_fctime[i] = 0;
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "R06",3) ==  0) {
                 strcpy((cdf_var + i*5),"mxrf");
                 cdf_fctime[i] = 1;
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "R12",3) ==  0) {
                 strcpy((cdf_var + i*5),"mxrf");
                 cdf_fctime[i] = 2;
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LMR  */
	case LF1:
	   /* cdf_fctime grid adjusted for analysis, 1 hr and 2 hr in LF1 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "H00",3) ==  0) {
                 strcpy((cdf_var + i*5),"mxrh");
                 cdf_fctime[i] = 0;
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "H06",3) ==  0) {
                 strcpy((cdf_var + i*5),"mxrh");
                 cdf_fctime[i] = 1;
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "H12",3) ==  0) {
                 strcpy((cdf_var + i*5),"mxrh");
                 cdf_fctime[i] = 2;
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LF1  */
	case L1S:
	   /* cdf_fctime grid left at zero, only analysis in L1S */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "S01",3) ==  0) {
                 strcpy((cdf_var + i*5),"s1hr");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "STO",3) ==  0) {
                 strcpy((cdf_var + i*5),"stot");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "R01",3) ==  0) {
                 strcpy((cdf_var + i*5),"pc");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RTO",3) ==  0) {
                 strcpy((cdf_var + i*5),"pt");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case L1S  */
	case LPS:
	   /* cdf_fctime grid left at zero, only analysis in LPS */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "REF",3) ==  0) {
                 strcpy((cdf_var + i*5),"ref");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LPS  */
	case LRP:
	   /* cdf_fctime grid left at zero, only analysis in LRP */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LRP",3) ==  0) {
                 strcpy((cdf_var + i*5),"icg");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LRP  */
	case LBA:
	   /* cdf_fctime grid left at zero, only analysis in LBA */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "HT",3) ==  0) {
                 strcpy((cdf_var + i*5),"bz");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "U",3) ==  0) {
                 strcpy((cdf_var + i*5),"bu");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "V",3) ==  0) {
                 strcpy((cdf_var + i*5),"bv");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "OM",3) ==  0) {
                 strcpy((cdf_var + i*5),"bw");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LBA  */
	case LC3:
	   /* cdf_fctime grid left at zero, only analysis in LC3 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LC3",3) ==  0) {
                 strcpy((cdf_var + i*5),"camt");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LC3  */
	case LWC:
	   /* cdf_fctime grid left at zero, only analysis in LWC */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LWC",3) ==  0) {
                 strcpy((cdf_var + i*5),"lwc");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ICE",3) ==  0) {
                 strcpy((cdf_var + i*5),"ice");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PCN",3) ==  0) {
                 strcpy((cdf_var + i*5),"pcn");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LWC  */
	case LIL:
	   /* cdf_fctime grid left at zero, only analysis in LIL */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LIL",3) ==  0) {
                 strcpy((cdf_var + i*5),"ilw");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LIL  */
	case LCB:
	   /* cdf_fctime grid left at zero, only analysis in LCB */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LCB",3) ==  0) {
                 strcpy((cdf_var + i*5),"cbas");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LCT",3) ==  0) {
                 strcpy((cdf_var + i*5),"ctop");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "CCE",3) ==  0) {
                 strcpy((cdf_var + i*5),"cce");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LCB  */
	case LCT:
	   /* cdf_fctime grid left at zero, only analysis in LCT */
           for (i = 0; i < *kkmax; i++) {
              if ((strncmp((var_req + i*4), "PTY",3) ==  0) ||
                  (strncmp((var_req + i*4), "SPT",3) ==  0)) {
                 strcpy((cdf_var + i*5),"spt");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PTT",3) ==  0) {
                 strcpy((cdf_var + i*5),"ptt");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SCT",3) ==  0) {
                 strcpy((cdf_var + i*5),"sct");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LCT  */
	case LCV:
	   /* cdf_fctime grid left at zero, only analysis in LCV */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LCV",3) ==  0) {
                 strcpy((cdf_var + i*5),"ccov");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "CSC",3) ==  0) {
                 strcpy((cdf_var + i*5),"csc");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "CWT",3) ==  0) {
                 strcpy((cdf_var + i*5),"cwt");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LCV  */
	case LMD:
	   /* cdf_fctime grid left at zero, only analysis in LMD */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LMD",3) ==  0) {
                 strcpy((cdf_var + i*5),"mcd");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LMD  */
	case LCO:
	   /* cdf_fctime grid left at zero, only analysis in LCO */
           for (i = 0; i < *kkmax; i++) {
              if ((strncmp((var_req + i*4), "OM",3) ==  0) ||
                  (strncmp((var_req + i*4), "COM",3) ==  0)) {
                 strcpy((cdf_var + i*5),"cw");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LCO  */
	case LTY:
	   /* cdf_fctime grid left at zero, only analysis in LTY */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "PTY",3) ==  0) {
                 strcpy((cdf_var + i*5),"ptyp");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "CTY",3) ==  0) {
                 strcpy((cdf_var + i*5),"ctyp");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LTY  */
	case LCP:
	   /* cdf_fctime grid left at zero, only analysis in LCP */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LCP",3) ==  0) {
                 strcpy((cdf_var + i*5),"ccpc");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
	   break;  /*  case LCP  */
	case LVD:
	   /* cdf_fctime grid left at zero, only analysis in LVD */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "S8W",3) ==  0) {
                 strcpy((cdf_var + i*5),"s8w");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S8C",3) ==  0) {
                 strcpy((cdf_var + i*5),"s8c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SVS",3) ==  0) {
                 strcpy((cdf_var + i*5),"svs");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SVN",3) ==  0) {
                 strcpy((cdf_var + i*5),"svn");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ALB",3) ==  0) {
                 strcpy((cdf_var + i*5),"alb");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S3A",3) ==  0) {
                 strcpy((cdf_var + i*5),"s3a");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S3C",3) ==  0) {
                 strcpy((cdf_var + i*5),"s3c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S4A",3) ==  0) {
                 strcpy((cdf_var + i*5),"s4a");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S4C",3) ==  0) {
                 strcpy((cdf_var + i*5),"s4c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S5A",3) ==  0) {
                 strcpy((cdf_var + i*5),"s5a");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S5C",3) ==  0) {
                 strcpy((cdf_var + i*5),"s5c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S8A",3) ==  0) {
                 strcpy((cdf_var + i*5),"s8a");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SCA",3) ==  0) {
                 strcpy((cdf_var + i*5),"sca");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SCC",3) ==  0) {
                 strcpy((cdf_var + i*5),"scc");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LVD  */
	case LVE:
	   /* cdf_fctime grid left at zero, only analysis in LVD */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "D8W",3) ==  0) {
                 strcpy((cdf_var + i*5),"d8w");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "D8C",3) ==  0) {
                 strcpy((cdf_var + i*5),"d8c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SVS",3) ==  0) {
                 strcpy((cdf_var + i*5),"svs");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SVN",3) ==  0) {
                 strcpy((cdf_var + i*5),"svn");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ALB",3) ==  0) {
                 strcpy((cdf_var + i*5),"alb");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "D3A",3) ==  0) {
                 strcpy((cdf_var + i*5),"d3a");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "D3C",3) ==  0) {
                 strcpy((cdf_var + i*5),"d3c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "D4A",3) ==  0) {
                 strcpy((cdf_var + i*5),"d4a");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "D4C",3) ==  0) {
                 strcpy((cdf_var + i*5),"d4c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "D5A",3) ==  0) {
                 strcpy((cdf_var + i*5),"d5a");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "D5C",3) ==  0) {
                 strcpy((cdf_var + i*5),"d5c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "D8A",3) ==  0) {
                 strcpy((cdf_var + i*5),"d8a");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "DCA",3) ==  0) {
                 strcpy((cdf_var + i*5),"dca");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "DCC",3) ==  0) {
                 strcpy((cdf_var + i*5),"dcc");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LVE  */
	case LMA:
	   /* cdf_fctime grid left at zero, only analysis in LMA */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "HT",3) ==  0) {
                 strcpy((cdf_var + i*5),"maz");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "T",3) ==  0) {
                 strcpy((cdf_var + i*5),"mat");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RH",3) ==  0) {
                 strcpy((cdf_var + i*5),"marh");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "U",3) ==  0) {
                 strcpy((cdf_var + i*5),"mau");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "V",3) ==  0) {
                 strcpy((cdf_var + i*5),"mav");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LMA  */
	case LMF:
	   /* cdf_fctime grid left at zero, only analysis in LMF */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "HT",3) ==  0) {
                 strcpy((cdf_var + i*5),"maz");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "T",3) ==  0) {
                 strcpy((cdf_var + i*5),"mat");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RH",3) ==  0) {
                 strcpy((cdf_var + i*5),"marh");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "U",3) ==  0) {
                 strcpy((cdf_var + i*5),"mau");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "V",3) ==  0) {
                 strcpy((cdf_var + i*5),"mav");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LMF  */
	case Z02:
	   /* cdf_fctime grid left at zero, only analysis in Z02 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "S8W",3) ==  0) {
                 strcpy((cdf_var + i*5),"s8w");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S8C",3) ==  0) {
                 strcpy((cdf_var + i*5),"s8c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SVS",3) ==  0) {
                 strcpy((cdf_var + i*5),"svs");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SVN",3) ==  0) {
                 strcpy((cdf_var + i*5),"svn");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ALB",3) ==  0) {
                 strcpy((cdf_var + i*5),"alb");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case Z02  */
	case RAM:
	   /* cdf_fctime grid left at zero, only analysis in RAM */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "HT",3) ==  0) {
                 strcpy((cdf_var + i*5),"rz");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "U",3) ==  0) ||
                       (strncmp((var_req + i*4), "U3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"ru");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "V",3) ==  0) ||
                       (strncmp((var_req + i*4), "V3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"rv");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "OM",3) ==  0) {
                 strcpy((cdf_var + i*5),"rw");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "T",3) ==  0) ||
                       (strncmp((var_req + i*4), "T3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"rt");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RMR",3) ==  0) {
                 strcpy((cdf_var + i*5),"rmr");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RH3",3) ==  0) {
                 strcpy((cdf_var + i*5),"rh3");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LWC",3) ==  0) {
                 strcpy((cdf_var + i*5),"lwc");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ICE",3) ==  0) {
                 strcpy((cdf_var + i*5),"ice");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RAI",3) ==  0) {
                 strcpy((cdf_var + i*5),"rai");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SNO",3) ==  0) {
                 strcpy((cdf_var + i*5),"sno");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PIC",3) ==  0) {
                 strcpy((cdf_var + i*5),"pic");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "REF",3) ==  0) {
                 strcpy((cdf_var + i*5),"ref");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case RAM  */
	case RSF:
	   /* cdf_fctime grid left at zero, only analysis in RSF */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "U",3) ==  0) {
                 strcpy((cdf_var + i*5),"rus");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "V",3) ==  0) {
                 strcpy((cdf_var + i*5),"rvs");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "P",3) ==  0) ||
                       (strncmp((var_req + i*4), "RP",3) ==  0)) {
                 strcpy((cdf_var + i*5),"rps");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "T",3) ==  0) {
                 strcpy((cdf_var + i*5),"rts");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "TD",3) ==  0) {
                 strcpy((cdf_var + i*5),"rtd");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RH",3) ==  0) {
                 strcpy((cdf_var + i*5),"rh");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LCB",3) ==  0) {
                 strcpy((cdf_var + i*5),"lcb");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LCT",3) ==  0) {
                 strcpy((cdf_var + i*5),"lct");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "MSL",3) ==  0) {
                 strcpy((cdf_var + i*5),"msl");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LIL",3) ==  0) {
                 strcpy((cdf_var + i*5),"lil");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "TPW",3) ==  0) {
                 strcpy((cdf_var + i*5),"tpw");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "R01",3) ==  0) {
                 strcpy((cdf_var + i*5),"r01");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RTO",3) ==  0) {
                 strcpy((cdf_var + i*5),"rto");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S01",3) ==  0) {
                 strcpy((cdf_var + i*5),"s01");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "STO",3) ==  0) {
                 strcpy((cdf_var + i*5),"sto");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "TH",3) ==  0) {
                 strcpy((cdf_var + i*5),"th");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "THE",3) ==  0) {
                 strcpy((cdf_var + i*5),"the");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PBE",3) ==  0) {
                 strcpy((cdf_var + i*5),"pbe");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "NBE",3) ==  0) {
                 strcpy((cdf_var + i*5),"nbe");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PS",3) ==  0) {
                 strcpy((cdf_var + i*5),"ps");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "CCE",3) ==  0) {
                 strcpy((cdf_var + i*5),"cce");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "VIS",3) ==  0) {
                 strcpy((cdf_var + i*5),"vis");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LCV",3) ==  0) {
                 strcpy((cdf_var + i*5),"lcv");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LMT",3) ==  0) {
                 strcpy((cdf_var + i*5),"lmt");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SPT",3) ==  0) {
                 strcpy((cdf_var + i*5),"spt");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LHE",3) ==  0) {
                 strcpy((cdf_var + i*5),"lhe");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LI",3) ==  0) {
                 strcpy((cdf_var + i*5),"li");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "HI",3) ==  0) {
                 strcpy((cdf_var + i*5),"hi");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case RSF  */
	case RSM:
	   /* cdf_fctime grid left at zero, only analysis in RSM */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LSM",3) ==  0) {
                 strcpy((cdf_var + i*5),"lsm");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case RSM  */
	case VRC:
	   /* cdf_fctime grid left at zero, only analysis in VRC */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "REF",3) ==  0) {
                 strcpy((cdf_var + i*5),"now");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case VRC  */
	case LM1:
	   /* cdf_fctime grid left at zero, only analysis in LM1 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "LSM",3) ==  0) {
                 strcpy((cdf_var + i*5),"lsm");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LM1  */
	case LM2:
	   /* cdf_fctime grid left at zero, only analysis in LM2 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "CIV",3) ==  0) {
                 strcpy((cdf_var + i*5),"civ");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "DWF",3) ==  0) {
                 strcpy((cdf_var + i*5),"dwf");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "WX",3) ==  0) {
                 strcpy((cdf_var + i*5),"wx");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "EVP",3) ==  0) {
                 strcpy((cdf_var + i*5),"evp");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SC",3) ==  0) {
                 strcpy((cdf_var + i*5),"sc");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SM",3) ==  0) {
                 strcpy((cdf_var + i*5),"sm");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "MWF",3) ==  0) {
                 strcpy((cdf_var + i*5),"mwf");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LM2  */
	case VRD:
	   /* cdf_fctime grid left at zero, only analysis in VRD */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "REF",3) ==  0) {
                 strcpy((cdf_var + i*5),"refd");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "VEL",3) ==  0) {
                 strcpy((cdf_var + i*5),"veld");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case VRD  */
	case V00:
	   /* cdf_fctime grid left at zero, only analysis in V00 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "REF",3) ==  0) {
                 strcpy((cdf_var + i*5),"refd");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "VEL",3) ==  0) {
                 strcpy((cdf_var + i*5),"veld");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "NYQ",3) ==  0) {
                 strcpy((cdf_var + i*5),"nyqd");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case V00  */
	case LGA:
	   /* cdf_fctime grid left at zero, only analysis in LGA */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "HT",3) ==  0) {
                 strcpy((cdf_var + i*5),"mgz");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "T",3) ==  0) ||
                       (strncmp((var_req + i*4), "T3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"mgt");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SH",3) ==  0) {
                 strcpy((cdf_var + i*5),"mgrh");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "U",3) ==  0) ||
                       (strncmp((var_req + i*4), "U3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"mgu");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "V",3) ==  0) ||
                       (strncmp((var_req + i*4), "V3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"mgv");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LGA  */
	case LGF:
	   /* cdf_fctime grid left at zero, only analysis in LGF */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "HT",3) ==  0) {
                 strcpy((cdf_var + i*5),"mfz");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "T",3) ==  0) ||
                       (strncmp((var_req + i*4), "T3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"mft");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SH",3) ==  0) {
                 strcpy((cdf_var + i*5),"mfrh");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "U",3) ==  0) ||
                       (strncmp((var_req + i*4), "U3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"mfu");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "V",3) ==  0) ||
                       (strncmp((var_req + i*4), "V3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"mfv");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LGF  */
	case LN3:
	   /* cdf_fctime grid left at zero, only analysis in LN3 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "R04",3) ==  0) {
                 strcpy((cdf_var + i*5),"r04");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "R48",3) ==  0) {
                 strcpy((cdf_var + i*5),"r48");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "R8C",3) ==  0) {
                 strcpy((cdf_var + i*5),"r8c");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ET",3) ==  0) {
                 strcpy((cdf_var + i*5),"et");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RCO",3) ==  0) {
                 strcpy((cdf_var + i*5),"rco");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "VIL",3) ==  0) {
                 strcpy((cdf_var + i*5),"vil");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case LN3  */
	case MM5:
	   /* cdf_fctime grid left at zero, only analysis in MM5 */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "HT",3) ==  0) {
                 strcpy((cdf_var + i*5),"rz");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "U",3) ==  0) ||
                       (strncmp((var_req + i*4), "U3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"ru");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "V",3) ==  0) ||
                       (strncmp((var_req + i*4), "V3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"rv");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "OM",3) ==  0) {
                 strcpy((cdf_var + i*5),"rw");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "T",3) ==  0) ||
                       (strncmp((var_req + i*4), "T3",3) ==  0)) {
                 strcpy((cdf_var + i*5),"rt");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SH",3) ==  0) {
                 strcpy((cdf_var + i*5),"rsh");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RH3",3) ==  0) {
                 strcpy((cdf_var + i*5),"rh3");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LWC",3) ==  0) {
                 strcpy((cdf_var + i*5),"lwc");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ICE",3) ==  0) {
                 strcpy((cdf_var + i*5),"ice");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case MM5  */
	case MSF:
	   /* cdf_fctime grid left at zero, only analysis in MSF */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "U",3) ==  0) {
                 strcpy((cdf_var + i*5),"rus");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "V",3) ==  0) {
                 strcpy((cdf_var + i*5),"rvs");
                 process_var[i] = 1;
              }
              else if ((strncmp((var_req + i*4), "P",3) ==  0) ||
                       (strncmp((var_req + i*4), "RP",3) ==  0)) {
                 strcpy((cdf_var + i*5),"rps");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "T",3) ==  0) {
                 strcpy((cdf_var + i*5),"rts");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "TD",3) ==  0) {
                 strcpy((cdf_var + i*5),"rtd");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RH",3) ==  0) {
                 strcpy((cdf_var + i*5),"rh");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LCB",3) ==  0) {
                 strcpy((cdf_var + i*5),"lcb");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LCT",3) ==  0) {
                 strcpy((cdf_var + i*5),"lct");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "MSL",3) ==  0) {
                 strcpy((cdf_var + i*5),"msl");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LIL",3) ==  0) {
                 strcpy((cdf_var + i*5),"lil");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "TPW",3) ==  0) {
                 strcpy((cdf_var + i*5),"tpw");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "R01",3) ==  0) {
                 strcpy((cdf_var + i*5),"r01");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "RTO",3) ==  0) {
                 strcpy((cdf_var + i*5),"rto");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "S01",3) ==  0) {
                 strcpy((cdf_var + i*5),"s01");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "STO",3) ==  0) {
                 strcpy((cdf_var + i*5),"sto");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "TH",3) ==  0) {
                 strcpy((cdf_var + i*5),"th");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "THE",3) ==  0) {
                 strcpy((cdf_var + i*5),"the");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PBE",3) ==  0) {
                 strcpy((cdf_var + i*5),"pbe");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "NBE",3) ==  0) {
                 strcpy((cdf_var + i*5),"nbe");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "PS",3) ==  0) {
                 strcpy((cdf_var + i*5),"ps");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "CCE",3) ==  0) {
                 strcpy((cdf_var + i*5),"cce");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "VIS",3) ==  0) {
                 strcpy((cdf_var + i*5),"vis");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LCV",3) ==  0) {
                 strcpy((cdf_var + i*5),"lcv");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LMT",3) ==  0) {
                 strcpy((cdf_var + i*5),"lmt");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "SPT",3) ==  0) {
                 strcpy((cdf_var + i*5),"spt");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LHE",3) ==  0) {
                 strcpy((cdf_var + i*5),"lhe");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "LI",3) ==  0) {
                 strcpy((cdf_var + i*5),"li");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "HI",3) ==  0) {
                 strcpy((cdf_var + i*5),"hi");
                 process_var[i] = 1;
              }
              else
              	 count = count + 1;
           }
           break;  /*  case MSF  */
        case SST:
           /* cdf_fctime grid left at zero, only analysis in SST */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "ST1",3) ==  0) {
                 strcpy((cdf_var + i*5),"st1");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ST2",3) ==  0) {
                 strcpy((cdf_var + i*5),"st2");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ST3",3) ==  0) {
                 strcpy((cdf_var + i*5),"st3");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ST4",3) ==  0) {
                 strcpy((cdf_var + i*5),"st4");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "ST5",3) ==  0) {
                 strcpy((cdf_var + i*5),"st5");
                 process_var[i] = 1;
              }
              else
                 count = count + 1;
           }
           break;  /*  case SST  */
        case VEG:
           /* cdf_fctime grid left at zero, only analysis in VEG */
           for (i = 0; i < *kkmax; i++) {
              if (strncmp((var_req + i*4), "VG1",3) ==  0) {
                 strcpy((cdf_var + i*5),"vg1");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "VG2",3) ==  0) {
                 strcpy((cdf_var + i*5),"vg2");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "VG3",3) ==  0) {
                 strcpy((cdf_var + i*5),"vg3");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "VG4",3) ==  0) {
                 strcpy((cdf_var + i*5),"vg4");
                 process_var[i] = 1;
              }
              else if (strncmp((var_req + i*4), "VG5",3) ==  0) {
                 strcpy((cdf_var + i*5),"vg5");
                 process_var[i] = 1;
              }
              else
                 count = count + 1;
           }
           break;  /*  case VEG  */
        case LS8:
           /* cdf_fctime grid left at zero, only analysis in LS8 */
           for (i = 0; i < *kkmax; i++) {
             strncpy(cdf_var + i*5, var_req + i*4, 3);
             varptr = cdf_var + i*5;
             if (strncmp(varptr, "S",1) == 0) strncpy(varptr, "s", 1); 
             process_var[i] = 1;
           }
           break;  /*  case LS8  */


	}  /* end of switch */
	return count;
}
/***************************************************************************
* OPEN_CDF
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Open netCDF file
*	Purpose			To open a netCDF file for read mode
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 7/11/91
*			modified for new READLAPSDATA 1/93 Linda Wharton
*
*	Input :
*		mode		A flag indicating whether to open the file 
*				  for writing (NC_WRITE) or reading 
*				  (NC_NOWRITE). 
*		fname		Name of the netCDF file to open
*               no_laps_diag    Flag, when not 0 suppresses printf msg
*	Output :
*		None
*	Globals :
*		None
*	Returns :
*		file id of netCDF file
*		-1 if an error occurs
*****************************************************************************/


#ifdef __STDC__
#if defined(__alpha) || defined(__IP21)
int open_cdf (int mode, char *fname,int *no_laps_diag)
#else
int open_cdf (int mode, char *fname,long *no_laps_diag)
#endif
#else
#if defined(__alpha) || defined(__IP21)
int open_cdf (mode, fname, no_laps_diag)
int mode;
char *fname;
int *no_laps_diag;
#else
int open_cdf (mode, fname, no_laps_diag)
int mode;
char *fname;
long *no_laps_diag;
#endif
#endif
{
	int	 cdfid, istatus, i;

	ncopts = 0;

	log_diag (2, "open_cdf: cdf file name = %s\n", fname);

/* open the netcdf file and get the file id and the variable id */
	cdfid = ncopen (fname, mode);
	
	ncopts = NC_VERBOSE;
	
	if (cdfid == (-1)) {
	  if (*no_laps_diag == 0) { 
/*	     printf("open_cdf: cannot open file as netCDF %s.\n", fname); */
	  }
	  return -1;
	}
	else {
		log_diag (2, "open_cdf: %s file open. Cdfid = %d\n", 
			fname, cdfid);
		return cdfid;
	}
}

/*****************************************************************************
*	CDF_RETRIEVE_HDR
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Retrieve Grid Data
*	Purpose			To retrieve header info from a netcdf file
*
*	Designer/Programmer : Linda Wharton
*	Modifications : original 1/93
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
#if defined(__alpha) || defined(__IP21)
int cdf_retrieve_hdr(int i_cdfid,int *imaxn, int *jmaxn, int *kmaxn,
		     char *laps_dom_file,char *asctime, char *version,
		     char *model, char *origin, int *num_variables)
#else
int cdf_retrieve_hdr(int i_cdfid,long *imaxn, long *jmaxn, long *kmaxn,
		     char *laps_dom_file,char *asctime, char *version,
		     char *model, char *origin, long *num_variables)
#endif
#else
#if defined(__alpha) || defined(__IP21)
int cdf_retrieve_hdr(i_cdfid,imaxn,jmaxn,kmaxn,laps_dom_file,asctime, 
		     version,model,origin,num_variables)
int i_cdfid;                               
int *imaxn;
int *jmaxn;
int *kmaxn;
char *laps_dom_file;
char *asctime;
char *version;
char *model;
char *origin;
int *num_variables;
#else
int cdf_retrieve_hdr(i_cdfid,imaxn,jmaxn,kmaxn,laps_dom_file,asctime, 
		     version,model,origin,num_variables)
int i_cdfid;                               
long *imaxn;
long *jmaxn;
long *kmaxn;
char *laps_dom_file;
char *asctime;
char *version;
char *model;
char *origin;
long *num_variables;
#endif
#endif
{
	int i_status, i_varid, i, i_version, temp, str_len;
#if defined(__alpha) || defined(__IP21)
        long mindex[1], start[1], count_asc[1], count_ldom[1], count_long[1]; 
#else
        int mindex[1], start[1], count_asc[1], count_ldom[1], count_long[1]; 
#endif
	static char c_ver[5];
	char *t_ptr;

/* turn off the error handling done by the netCDF routines */
	ncopts = NC_VERBOSE;
	mindex[0] = 0;
	start[0] = 0;
	count_asc[0] = 18;
	count_ldom[0] = 12;
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
	log_diag (2, "cdf_read_grid: data variable name = %s\n", "kmax");
	if ((i_varid = ncvarid (i_cdfid, "kmax")) == (-1))
		return -1;

/* read the var from the netcdf file */
	i_status = ncvarget1 (i_cdfid, i_varid, mindex, kmaxn);
	if (i_status == (-1))
	   return -1;

/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", "version");
	if ((i_varid = ncvarid (i_cdfid, "version")) == (-1))
		return -1;
	   
/* read the version from the netcdf file */
	i_status = ncvarget1 (i_cdfid, i_varid, mindex, &i_version);
	if (i_status == (-1))
	   return -1;

/* convert integer value of version to character string */
	if (i_version < 10) {
	   strcpy(version,"  V");
	   itoa(i_version,c_ver,2);
	   strcat(version,c_ver);
	}
	else if (i_version < 100) {
	   strcpy(version," V");
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

/* get the variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", 
	          "laps_domain_file");
	if ((i_varid = ncvarid (i_cdfid, "laps_domain_file")) == (-1))
		return -1;
	   
/* read the laps domain name from the netcdf file */
	i_status = ncvarget (i_cdfid, i_varid, start, count_ldom, 
	                     laps_dom_file);
	if (i_status == (-1))
	   return -1;
	else {
	   str_len = strlen(laps_dom_file);
	   if (str_len < 11) {
	      for (i = str_len+1; i < 12; i++)
	         strcat(laps_dom_file," ");
	   }
        }
        
/* get the  variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", "model");
	if ((i_varid = ncvarid (i_cdfid, "model")) == (-1))
		return -1;
	   
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
	if ((i_varid = ncvarid (i_cdfid, "origin")) == (-1))
		return -1;
	   
/* read the origin from the netcdf file */
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
	if ((i_varid = ncvarid (i_cdfid, "num_variables")) == (-1))
		return -1;

/* read the var from the netcdf file */
	i_status = ncvarget1 (i_cdfid, i_varid, mindex, num_variables);
	if (i_status == (-1))
	   return -1;
	   	   
/* normal return */

	return 0;
}

/***************************************************************************
* WRITE_HDR_CDF
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Write header into netCDF file
*	Purpose			To write a header into an open netCDF file.
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 9/5/90
*			adapted for new WRITELAPSDATA 1/93 Linda Wharton
*
*	Input :
*		i_cdfid 	The file id of a netCDF file opened for writing
*		kdim		Variable indicating number of levels of data
*				  to be written.
*		lvl		Array containing levels of data
*		var		Array containing LAPS variables for each level
*		lvl_coord	Array containing coordinates for each level
*		units		Array containing units associated with each level
*		comment		Array containing comment for each level
*		asctime		ASCII time of data
*	Output :
*		None
*	Globals :
*		None
*	Returns :
*		 0 if successful
*		-1 if an error occurs
**************************************************************************/
#ifdef __STDC__
#if defined(__alpha) || defined(__IP21)
int write_hdr_cdf(int i_cdfid,int *kdim,int *lvl,char *laps_dom_file,
                  char *asctime)
#else
int write_hdr_cdf(int i_cdfid,long *kdim,long *lvl,char *laps_dom_file,
                  char *asctime)
#endif
#else
#if defined(__alpha) || defined(__IP21)
int write_hdr_cdf(i_cdfid,kdim,lvl,laps_dom_file,asctime)
int i_cdfid;
int *kdim;
int *lvl;
char *laps_dom_file;
char *asctime;
#else
int write_hdr_cdf(i_cdfid,kdim,lvl,laps_dom_file,asctime)
int i_cdfid;
long *kdim;
long *lvl;
char *laps_dom_file;
char *asctime;
#endif
#endif
{
	int i_status, i, i_lvlid, i_ldomfid, i_asctmid;
#if defined(__alpha) || defined(__IP21)
        long ldom_ct[1], ldom_st[1], asc_ct[1], asc_st[1],
             lvl_ct[1],lvl_st[1];
#else
        int ldom_ct[1], ldom_st[1], asc_ct[1], asc_st[1],
            lvl_ct[1],lvl_st[1];
#endif
	    
	 static long level_val;
	
	lvl_ct[0] = 1;
	ldom_ct[0] = 12;
	ldom_st[0] = 0;
	asc_ct[0] = 18;
	asc_st[0] = 0;

/* get the variable ids */

	log_diag (2, "Data variable name = %s\n", "lvl");
	if ((i_lvlid = ncvarid (i_cdfid,"lvl" )) == (-1))
		return -1;
   	if ((i_ldomfid = ncvarid (i_cdfid,"laps_domain_file" )) == (-1))
		return -1;
 	if ((i_asctmid = ncvarid (i_cdfid,"asctime" )) == (-1))
		return -1;
 
/* write header data */
          i_status = ncvarput(i_cdfid,i_ldomfid,ldom_st,ldom_ct,
        		    (void *)laps_dom_file);
	if (i_status == (-1)) return -1;

        i_status = ncvarput(i_cdfid,i_asctmid,asc_st,asc_ct,
        		    (void *)asctime);
	if (i_status == (-1)) return -1;

        for (i = 0; i < *kdim; i++) {
		lvl_st[0] = i;
		
		level_val = (long)lvl[i];
        	i_status = ncvarput1(i_cdfid,i_lvlid,lvl_st,
        			    (void *)&level_val );
		if (i_status == (-1)) return -1;				
        }
	
	if (i_status == (-1))
		return -1;
	else
		return 0;
}
/*****************************************************************************
* CDF_UPDATE_LAPS
*	Category		Product Management
*	Group			General Purpose Database
*	Module			Update the LAPS netCDF file with a new grid
*	Purpose			To coordinate the writing of an LAPS grid into
*				  the appropriate netCDF file: find the 
*				  appropriate location to write the grid, 
*				  and call the routine that does the writing.
*
*	Designer/Programmer : MarySue Schultz
*	Modifications : original 9/5/90
*			adapted for new WRITELAPSDATA 1/93 Linda Wharton
*
*	Input :
*		i_level		The level: 0,100,150,200,250,300,350,400,450,
*					   500,550,600,650,700,750,800,850,
*					   900,950,1000,1050,1100
*		i_fctime	The forecast time: 0
*		s_field		A character string identifying the field:
*				    w		wind omega
*				    u		Eastward wind 
*				    v		Northward wind
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
int cdf_update_laps (int i_cdfid, int i_level, int i_fctime, 
		     char *s_field,float *gptr, char *commnt,
		     char *comm_ptr)
#else
int cdf_update_laps (i_cdfid, i_level, i_fctime, s_field, gptr, 
		     commnt, comm_ptr)
int i_cdfid;
int i_level;
int i_fctime;
char *s_field;
float *gptr;
char *commnt;
char *comm_ptr;
#endif
{
	int i_status;
	cdf_grid_info dims;

	log_diag (2, "cdf_update_laps:grid name = %s\n", s_field);

/* get the forecast time index of the data array */
	dims.fctime_coord = cdf_get_index (i_cdfid, i_fctime, "fctimes");
	if (dims.fctime_coord == (-1)){
		ncclose (i_cdfid);
		return -1;
	}
	log_diag (2, "cdf_update_laps:fctime coord = %d\n", dims.fctime_coord);

/* get the level index of the data array */
	dims.level_coord = cdf_get_index (i_cdfid, i_level, "level");
	if (dims.level_coord == (-1)){
		ncclose (i_cdfid);
		return -1;
	}
	log_diag (2, "cdf_update_laps: level coord = %d\n", dims.level_coord);

/* get the x and y dimension sizes */
	dims.x_dim = cdf_dim_size (i_cdfid, "lon");
	dims.y_dim = cdf_dim_size (i_cdfid, "lat");

/* write the grid to the netcdf file */
	i_status = cdf_write_grid (i_cdfid, &dims, s_field, gptr,
				   commnt, comm_ptr);
	if (i_status == 0)
		log_diag (2, "cdf_update_laps: cdf write ok\n",i_status);
	else {
		log_diag (1, "cdf_update_laps: error during cdf write\n",i_status);
		return -1;		
	}

/* update the inventory in the netcdf file */
	 i_status = cdf_update_laps_inv (i_cdfid, &dims, s_field);

	return i_status;
}

/*****************************************************************************
*	CDF_RETRIEVE_LAPS_GRID
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
*		i_level		One of the following levels:
*				  0, 100, 150, 200, 250, 300, 350,400, 
*				  450, 500, 550, 600, 650, 700, 750, 800,
*				  850, 900, 950, 1000, 1050, 1100
*		i_fctime	One of the following forecast times:
*				  0
*		s_field		A character string identifying the field
*		gptr		Pointer to location in data array to write data.
*		cptr		Pointer to location in data array to write comments
*		lvlptr		Pointer to location in data array to write 
*				lvl_coordinates.
*		uptr		Pointer to location in data array to write units.
*		ext		Extension identifying data expected in file
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
int cdf_retrieve_laps_grid(int i_cdfid,int i_level,int i_fctime,
			   char *s_field, float *gptr, char *cptr,
			   char *lvlptr, char *uptr, char *ext)
#else
int cdf_retrieve_laps_grid(i_cdfid, i_level, i_fctime, s_field, gptr, 
			   cptr, lvlptr, uptr, ext)
int i_cdfid;
int i_level;
int i_fctime;
char *s_field;
float *gptr;
char *cptr;
char *lvlptr;
char *uptr;
char *ext;
#endif
{
	int i_status, i_invflag, i_varid;
        long start[4],count[4],start_c[3],count_c[3];
	char var_name[13];
	cdf_grid_info dims;

/*	printf("cdf_ret_grd: level = %d   fctime = %d   field = %s\n",
		i_level, i_fctime, s_field);  */

/* turn off the error handling done by the netCDF routines */
	ncopts = NC_VERBOSE;

/* get the forecast time index of the data array */
	dims.fctime_coord = cdf_get_index (i_cdfid, i_fctime, "fctimes");
	if (dims.fctime_coord == (-1)){
		printf("cdf_retrieve_laps: no forecast time %d\n", 
		 	  i_fctime);
		return -1;
	}
	log_diag (2, "cdf_retrieve_laps: fctime coord = %d\n", 
		  dims.fctime_coord);

/* get the level index of the data array */
	dims.level_coord = cdf_get_index (i_cdfid, i_level, "level");
	if (dims.level_coord == (-1)){
		printf("cdf_retrieve_laps: no level %d\n", i_level);
		return -1;
	}
	log_diag (2, "cdf_retrieve_laps: level coord = %d\n", dims.level_coord);

/* check the inventory variable to see if this grid is available */
	i_invflag = cdf_check_laps_inv (i_cdfid, dims.level_coord, 
					dims.fctime_coord, s_field);
	if (i_invflag == (-1)){
		printf("cdf_retrieve_laps: grid not available\n");
		return -1;
	}

/* get the x and y dimension sizes */
	dims.x_dim = cdf_dim_size (i_cdfid, "lon");
	dims.y_dim = cdf_dim_size (i_cdfid, "lat");

	log_diag (2, "cdf_retrieve_laps: x = %d   y = %d\n",dims.x_dim,
							    dims.y_dim);

/* get the data variable id */
	log_diag (2, "cdf_read_grid: data variable name = %s\n", s_field);
	if ((i_varid = ncvarid (i_cdfid, s_field)) == (-1)) {
		printf("cdf_retrieve_laps: no grid available.\n");
		return -1;
	}
			   
/* construct the arrays needed to read the grid */
	start[0] = dims.fctime_coord;
	start[1] = dims.level_coord;
	start[2] = 0;
	start[3] = 0;

	count[0] = 1;
	count[1] = 1;
	count[2] = dims.y_dim;
	count[3] = dims.x_dim;

/* read the grid from the netcdf file */
	i_status = ncvarget (i_cdfid, i_varid, start, count, gptr);
	if (i_status == (-1)) {
	   printf("cdf_retrieve_laps: error retrieving data %s grid.\n", 
		 	 *s_field);
	   return -1;
	}

/* get attributes lvl_coord and LAPS_units */
	if ((strncmp(ext,"LW3",3) == 0) && (strncmp(s_field,"w",1) != 0))
	   if (i_level == 0) {
	      i_status = ncattget (i_cdfid, i_varid, "lvl_coord_sfc", lvlptr);
	      if (i_status == (-1)) {
		 printf("cdf_retrieve_laps: error retrieving lvl_coord_sfc.\n");
	         return -1;
	      }
	   }
	   else {
	      i_status = ncattget (i_cdfid, i_varid, "lvl_coord_3D", lvlptr); 
	      if (i_status == (-1)) {
		 printf("cdf_retrieve_laps: error retrieving lvl_coord_3D.\n");
	         return -1;
	      }
	   }
	else {
	   i_status = ncattget (i_cdfid, i_varid, "lvl_coord", lvlptr);
	   if (i_status == (-1)) {
	      printf("cdf_retrieve_laps: error retrieving lvl_coord.\n");
	      return -1;
	   }
	}
	
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

/* construct the arrays needed to read the grid */
	start_c[0] = dims.fctime_coord;
	start_c[1] = dims.level_coord;
	start_c[2] = 0;

	count_c[0] = 1;
	count_c[1] = 1;
	count_c[2] = 126;

/* read the comment from the netcdf file */
	i_status = ncvarget (i_cdfid, i_varid, start_c, count_c, cptr);
	if (i_status == (-1)) {
	   printf("cdf_retrieve_laps: error retrieving comment.\n");
	   return -1;
	}
	   
/* normal return */

	return 0;
}
/*************************************************************************
*	WRITE_CDF_FILE
*	Category	Product Management
*	Group		General Purpose Database
*	Module		Write Grid data
*	Purpose		Write grid data into netCDF file.
*
*	Designer/Programmer : Linda Wharton
*	Modifications : original 1/93
*
*	Input:
*		fname		NetCDF filename to open
*		laps_dom_file	Name of the LAPS domain the grid is valid for
*		asctime		Ascii time of file
*		ext		Extension identifying data to be written
*		var		Array of LAPS variables to write out
*		lvl_coord	Level coordinates of each level to write
*		units		Units of each level to write
*		comment		Comments for each level to write
*		imax		X dimension of data in file
*		jmax		Y dimension of data in file
*		kmax		Number of grids available in file
*		kdim 		Dimension  of the var, lvl, comment,
*				  lvl_coord, units in the calling program
*		lvl		Levels of LAPS variables to write to
*		data		3D grid containing data to write
*	Output:
*		status		Returns status to calling subroutine
*	Globals:
*		none
*	Returns:
*		none
*****************************************************************************/
#ifdef __STDC__
#if defined(__alpha) || defined(__IP21)
void write_cdf_file(char *filname,short *s_length,
		    char *f_laps_dom_file,char *f_asctime,
                    char *f_ext,char *f_var,char *f_lvl_coord,char *f_units,
                    char *f_comment,int *imax,int *jmax,int *kmax,
                    int *kdim,int *lvl, float *data, int *status)
#else
void write_cdf_file(char *filname,short *s_length,
		    char *f_laps_dom_file,char *f_asctime,
                    char *f_ext,char *f_var,char *f_lvl_coord,char *f_units,
                    char *f_comment,long *imax,long *jmax,long *kmax,
                    long *kdim,long *lvl, float *data, long *status)
#endif
#else
#if defined(__alpha) || defined(__IP21)
void write_cdf_file(filname,s_length,f_laps_dom_file,f_asctime,f_ext,f_var,
		    f_lvl_coord,f_units,f_comment,imax,jmax,kmax,
                    kdim,lvl,data,status)
char *filname;
short *s_length;
char *f_laps_dom_file;
char *f_asctime;
char *f_ext;
char *f_var;
char *f_lvl_coord;
char *f_units;
char *f_comment;
int *imax;
int *jmax;
int *kmax;
int *kdim;
int *lvl;
float *data;
int *status;
#else
void write_cdf_file(filname,s_length,f_laps_dom_file,f_asctime,f_ext,f_var,
		    f_lvl_coord,f_units,f_comment,imax,jmax,kmax,
                    kdim,lvl,data,status)
char *filname;
short *s_length;
char *f_laps_dom_file;
char *f_asctime;
char *f_ext;
char *f_var;
char *f_lvl_coord;
char *f_units;
char *f_comment;
long *imax;
long *jmax;
long *kmax;
long *kdim;
long *lvl;
float *data;
long *status;
#endif
#endif
{
	enum grid {LW3, LH1, LH2, LH3, LH4, LQ3, LSX, LWM, LT1, LHE,
		   LIW, LMT, LMR, LF1, L1S, LPS, LRP, LBA, LC3, LWC,
		   LIL, LCB, LCT, LCV, LMD, LCO, LTY, LCP, LVD, LVE, 
		   LMA, LMF, Z02, RAM, RSF, RSM, VRC, LM1, LM2, VRD, 
                   V00, LGA, LGF, LN3, MM5, MSF, SST, VEG, LS8};
		   
	static char prefix[5],var[300][4],ext[32];
	static char comm_var[13],laps_dom_file[12],asctime[18];
        char lvl_coord[300][5],units[300][11],comment[300][126];
	int out_file, istat, i, i_level, fctime, process_vr, count;
        char origin[132],model[132];
        int dims[1], vid;
        long start[1],edges[1];
	enum grid grid_no;
	char fname[92], *varptr;

/****  Header info set in cre_XXX subroutine
	   1)  Record length (IMAX), 
	   2)  Number of lines (JMAX),
	   3)  Number of fields (KMAX),
	   4)  Number of 2-D layers (KDIM),
	   5)  Version number of write routine.
     
       Header info set in write_hdr_cdf subroutine
	   1)  Variable (VAR(KDIM)),
	   2)  Level (KDIM),
	   3)  Vertical coordinate (LVL_COORD(KDIM))
	   4)  Units(KDIM),
	   6)  Comment(KDIM),
	   7)  Ascii time.
****/

/* convert fortran string f_var to c string var */
 
        for (i = 0; i < *kdim; i++) {
          nstrncpy(var[i],(f_var+i*3),3);
          nstrncpy(lvl_coord[i],(f_lvl_coord+i*4),4);
          nstrncpy(units[i],(f_units+i*10),10);
          fstrncpy(comment[i],(f_comment+i*125),125);
        }
        nstrncpy(laps_dom_file,f_laps_dom_file,11);
        fstrncpy(asctime,f_asctime,17);
 
/* convert fortran file_name into C fname, and fortran f_ext into C ext  */
        nstrncpy(fname,filname,s_length);
        nstrncpy(ext,f_ext,31);
	
	if (strncmp(ext,"LW3",3) == 0) grid_no = LW3;
	if (strncmp(ext,"LH1",3) == 0) grid_no = LH1;
	if (strncmp(ext,"LH2",3) == 0) grid_no = LH2;
	if (strncmp(ext,"LH3",3) == 0) grid_no = LH3;
	if (strncmp(ext,"LH4",3) == 0) grid_no = LH4;
	if (strncmp(ext,"LQ3",3) == 0) grid_no = LQ3;
	if (strncmp(ext,"LSX",3) == 0) grid_no = LSX;
	if (strncmp(ext,"LWM",3) == 0) grid_no = LWM;
	if (strncmp(ext,"LT1",3) == 0) grid_no = LT1;
	if (strncmp(ext,"LHE",3) == 0) grid_no = LHE;
	if (strncmp(ext,"LIW",3) == 0) grid_no = LIW;
	if (strncmp(ext,"LMT",3) == 0) grid_no = LMT;
	if (strncmp(ext,"LMR",3) == 0) grid_no = LMR;
	if (strncmp(ext,"LF1",3) == 0) grid_no = LF1;
	if (strncmp(ext,"L1S",3) == 0) grid_no = L1S;
	if (strncmp(ext,"LPS",3) == 0) grid_no = LPS;
	if (strncmp(ext,"LRP",3) == 0) grid_no = LRP;
	if (strncmp(ext,"LBA",3) == 0) grid_no = LBA;
	if (strncmp(ext,"LC3",3) == 0) grid_no = LC3;
	if (strncmp(ext,"LWC",3) == 0) grid_no = LWC;
	if (strncmp(ext,"LIL",3) == 0) grid_no = LIL;
	if (strncmp(ext,"LCB",3) == 0) grid_no = LCB;
	if (strncmp(ext,"LCT",3) == 0) grid_no = LCT;
	if (strncmp(ext,"LCV",3) == 0) grid_no = LCV;
	if (strncmp(ext,"LMD",3) == 0) grid_no = LMD;
	if (strncmp(ext,"LCO",3) == 0) grid_no = LCO;
	if (strncmp(ext,"LTY",3) == 0) grid_no = LTY;
	if (strncmp(ext,"LCP",3) == 0) grid_no = LCP;
	if (strncmp(ext,"LVD",3) == 0) grid_no = LVD;
	if (strncmp(ext,"LVE",3) == 0) grid_no = LVE;
	if (strncmp(ext,"LMA",3) == 0) grid_no = LMA;
	if (strncmp(ext,"LMF",3) == 0) grid_no = LMF;
	if (strncmp(ext,"Z02",3) == 0) grid_no = Z02;
	if (strncmp(ext,"RAM",3) == 0) grid_no = RAM;
	if (strncmp(ext,"RSF",3) == 0) grid_no = RSF;
        if (strncmp(ext,"RSM",3) == 0) grid_no = RSM;
	if (strncmp(ext,"VRC",3) == 0) grid_no = VRC;
	if (strncmp(ext,"LM1",3) == 0) grid_no = LM1;
	if (strncmp(ext,"LM2",3) == 0) grid_no = LM2;
        if (strncmp(ext,"VRD",3) == 0) grid_no = VRD;
        if (strncmp(ext,"V01",3) == 0) grid_no = V00;
        if (strncmp(ext,"V02",3) == 0) grid_no = V00;
        if (strncmp(ext,"V03",3) == 0) grid_no = V00;
        if (strncmp(ext,"V04",3) == 0) grid_no = V00;
        if (strncmp(ext,"V05",3) == 0) grid_no = V00;
        if (strncmp(ext,"V06",3) == 0) grid_no = V00;
        if (strncmp(ext,"V07",3) == 0) grid_no = V00;
        if (strncmp(ext,"V08",3) == 0) grid_no = V00;
        if (strncmp(ext,"V09",3) == 0) grid_no = V00;
        if (strncmp(ext,"V10",3) == 0) grid_no = V00;
        if (strncmp(ext,"V11",3) == 0) grid_no = V00;
        if (strncmp(ext,"V12",3) == 0) grid_no = V00;
        if (strncmp(ext,"V13",3) == 0) grid_no = V00;
        if (strncmp(ext,"V14",3) == 0) grid_no = V00;
        if (strncmp(ext,"V15",3) == 0) grid_no = V00;
        if (strncmp(ext,"V16",3) == 0) grid_no = V00;
        if (strncmp(ext,"V17",3) == 0) grid_no = V00;
        if (strncmp(ext,"V18",3) == 0) grid_no = V00;
        if (strncmp(ext,"V19",3) == 0) grid_no = V00;
        if (strncmp(ext,"V20",3) == 0) grid_no = V00;
        if (strncmp(ext,"LGA",3) == 0) grid_no = LGA;
        if (strncmp(ext,"LGF",3) == 0) grid_no = LGF;
        if (strncmp(ext,"LN3",3) == 0) grid_no = LN3;
	if (strncmp(ext,"MM5",3) == 0) grid_no = MM5;
	if (strncmp(ext,"MSF",3) == 0) grid_no = MSF;
	if (strncmp(ext,"SST",3) == 0) grid_no = SST;
	if (strncmp(ext,"VEG",3) == 0) grid_no = VEG;
	if (strncmp(ext,"LS8",3) == 0) grid_no = LS8;
		
	switch (grid_no)
	{
        case LW3:
           if (*imax == NX && *jmax == NY && *kmax <= 63 && *kdim <= 63)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lw3(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if ((strncmp(var[i], "U",3) ==  0) ||
          	    (strncmp(var[i], "U3",3) ==  0)) {
          	   strcpy(prefix,"u");
          	   strcpy(comm_var,"u_comment");
          	}
          	else if ((strncmp(var[i], "V",3) ==  0) ||
          	         (strncmp(var[i], "V3",3) ==  0)) {
          	   strcpy(prefix,"v");
          	   strcpy(comm_var,"v_comment");
          	}
          	else if (strncmp(var[i], "OM",3) ==  0) {
          	   strcpy(prefix,"w");
          	   strcpy(comm_var,"w_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LW3  */
        case LH1:
           if (*imax == NX && *jmax == NY && *kmax == 1 && *kdim == 1)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lh1(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "PW",3) ==  0) {
          	   strcpy(prefix,"pw");
          	   strcpy(comm_var,"pw_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LH1  */
        case LH2:
           if (*imax == NX && *jmax == NY && *kmax <= 3 && *kdim <= 3)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lh2(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
         	if (strncmp(var[i], "PW1",3) ==  0) {
          	   strcpy(prefix,"lpw");
          	   strcpy(comm_var,"lpw_comment");
          	}
          	else if (strncmp(var[i], "PW2",3) ==  0) {
          	   strcpy(prefix,"lpw");
          	   strcpy(comm_var,"lpw_comment");
          	}
          	else if (strncmp(var[i], "PW3",3) ==  0) {
          	   strcpy(prefix,"lpw");
          	   strcpy(comm_var,"lpw_comment");
          	}
          	else if (strncmp(var[i], "PW",3) ==  0) {
          	   strcpy(prefix,"lpw");
          	   strcpy(comm_var,"lpw_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LH2  */
        case LH3:
           if (*imax == NX && *jmax == NY && *kmax <= 42 && *kdim <= 42)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lh3(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if ((strncmp(var[i], "RH",3) ==  0) ||
          	    (strncmp(var[i], "RH3",3) ==  0)) {
          	   strcpy(prefix,"rh");
          	   strcpy(comm_var,"rh_comment");
          	}
          	else if (strncmp(var[i], "RHL",3) ==  0) {
          	   strcpy(prefix,"rhl");
          	   strcpy(comm_var,"rhl_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LH3  */
        case LH4:
           if (*imax == NX && *jmax == NY && *kmax == 1 && *kdim == 1)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lh4(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "TPW",3) ==  0) {
          	   strcpy(prefix,"tpw");
          	   strcpy(comm_var,"tpw_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LH4  */
        case LQ3:
           if (*imax == NX && *jmax == NY && *kmax <= 21 && *kdim <= 21)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lq3(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "SH",3) ==  0) {
          	   strcpy(prefix,"sh");
          	   strcpy(comm_var,"sh_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LQ3  */
        case LSX:
           if (*imax == NX && *jmax == NY && *kmax <= 27 && *kdim <= 27)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lsx(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
          	process_vr = 0;
          	if (strncmp(var[i], "U",3) ==  0) {
          	   strcpy(prefix,"su");
          	   strcpy(comm_var,"su_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "V",3) ==  0)) {
          	   strcpy(prefix,"sv");
          	   strcpy(comm_var,"sv_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "P",3) ==  0)) {
          	   strcpy(prefix,"fp");
          	   strcpy(comm_var,"fp_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "T",3) ==  0)) {
          	   strcpy(prefix,"st");
          	   strcpy(comm_var,"st_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "TD",3) ==  0)) {
          	   strcpy(prefix,"std");
          	   strcpy(comm_var,"std_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "VV",3) ==  0)) {
          	   strcpy(prefix,"vv");
          	   strcpy(comm_var,"vv_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "RH",3) ==  0)) {
          	   strcpy(prefix,"srh");
          	   strcpy(comm_var,"srh_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "CCE",3) ==  0)) {
          	   strcpy(prefix,"ccg");
          	   strcpy(comm_var,"ccg_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "MSL",3) ==  0)) {
          	   strcpy(prefix,"mp");
          	   strcpy(comm_var,"mp_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "TAD",3) ==  0)) {
          	   strcpy(prefix,"ta");
          	   strcpy(comm_var,"ta_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "TH",3) ==  0)) {
          	   strcpy(prefix,"pot");
          	   strcpy(comm_var,"pot_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "THE",3) ==  0)) {
          	   strcpy(prefix,"ept");
          	   strcpy(comm_var,"ept_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "PS",3) ==  0)) {
          	   strcpy(prefix,"sp");
          	   strcpy(comm_var,"sp_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "VOR",3) ==  0)) {
          	   strcpy(prefix,"vor");
          	   strcpy(comm_var,"vor_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "MR",3) ==  0)) {
          	   strcpy(prefix,"mr");
          	   strcpy(comm_var,"mr_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "MRC",3) ==  0)) {
          	   strcpy(prefix,"mc");
          	   strcpy(comm_var,"mc_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "DIV",3) ==  0)) {
          	   strcpy(prefix,"d");
          	   strcpy(comm_var,"d_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "THA",3) ==  0)) {
          	   strcpy(prefix,"pta");
          	   strcpy(comm_var,"pta_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "MRA",3) ==  0)) {
          	   strcpy(prefix,"ma");
          	   strcpy(comm_var,"ma_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "LI",3) ==  0)) {
          	   strcpy(prefix,"li");
          	   strcpy(comm_var,"li_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "SPD",3) ==  0)) {
          	   strcpy(prefix,"spd");
          	   strcpy(comm_var,"spd_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "CSS",3) ==  0)) {
          	   strcpy(prefix,"cssi");
          	   strcpy(comm_var,"cssi_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "PBE",3) ==  0)) {
          	   strcpy(prefix,"pbe");
          	   strcpy(comm_var,"pbe_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "NBE",3) ==  0)) {
          	   strcpy(prefix,"nbe");
          	   strcpy(comm_var,"nbe_comment");
           	   process_vr = 1;
          	}
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "VIS",3) ==  0)) {
          	   strcpy(prefix,"vis");
          	   strcpy(comm_var,"vis_comment");
           	   process_vr = 1;
          	}   
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "FWX",3) ==  0)) {
          	   strcpy(prefix,"fwx");
          	   strcpy(comm_var,"fwx_comment");
           	   process_vr = 1;
          	}   
          	if ((process_vr == 0) && 
          	   (strncmp(var[i], "HI",3) ==  0)) {
          	   strcpy(prefix,"hi");
          	   strcpy(comm_var,"hi_comment");
           	   process_vr = 1;
          	}   
          	if (process_vr == 0) 
          	   count += 1;	
          	else {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LSX  */
        case LWM:
           if (*imax == NX && *jmax == NY && *kmax <= 2 && *kdim <= 2)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lwm(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "SU",3) ==  0) {
          	   strcpy(prefix,"u");
          	   strcpy(comm_var,"u_comment");
          	}
          	else if (strncmp(var[i], "SV",3) ==  0) {
          	   strcpy(prefix,"v");
          	   strcpy(comm_var,"v_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LWM  */
        case LT1:
           if (*imax == NX && *jmax == NY && *kmax <= 42 && *kdim <= 42)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lt1(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if ((strncmp(var[i], "T",3) ==  0) ||
          	    (strncmp(var[i], "T3",3) ==  0)) {
          	   strcpy(prefix,"t");
          	   strcpy(comm_var,"t_comment");
          	}
          	else if (strncmp(var[i], "HT",3) ==  0) {
          	   strcpy(prefix,"z");
          	   strcpy(comm_var,"z_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LT1  */
        case LHE:
           if (*imax == NX && *jmax == NY && *kmax <= 3 && *kdim <= 3)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lhe(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LHE",3) ==  0) {
          	   strcpy(prefix,"hel");
          	   strcpy(comm_var,"hel_comment");
          	}
          	else if (strncmp(var[i], "MU",3) ==  0) {
          	   strcpy(prefix,"mu");
          	   strcpy(comm_var,"mu_comment");
          	}
          	else if (strncmp(var[i], "MV",3) ==  0) {
          	   strcpy(prefix,"mv");
          	   strcpy(comm_var,"mv_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LHE  */
        case LIW:
           if (*imax == NX && *jmax == NY && *kmax <= 2 && *kdim <= 2)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_liw(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LIW",3) ==  0) {
          	   strcpy(prefix,"liw");
          	   strcpy(comm_var,"liw_comment");
          	}
          	else if (strncmp(var[i], "W",3) ==  0) {
          	   strcpy(prefix,"w");
          	   strcpy(comm_var,"w_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LIW  */
        case LMT:
           if (*imax == NX && *jmax == NY && *kmax <= 2 && *kdim <= 2)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lmt(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LMT",3) ==  0) {
          	   strcpy(prefix,"etop");
          	   strcpy(comm_var,"etop_comment");
          	}
          	else if (strncmp(var[i], "LLR",3) ==  0) {
          	   strcpy(prefix,"llr");
          	   strcpy(comm_var,"llr_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LMT  */
        case LMR:
           if (*imax == NX && *jmax == NY && *kmax <= 3 && *kdim <= 3)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lmr(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "R00",3) ==  0) {
          	   strcpy(prefix,"mxrf");
          	   strcpy(comm_var,"mxrf_comment");
          	   fctime = 0;
          	}
          	else if (strncmp(var[i], "R06",3) ==  0) {
          	   strcpy(prefix,"mxrf");
          	   strcpy(comm_var,"mxrf_comment");
          	   fctime = 1;
          	}
          	else if (strncmp(var[i], "R12",3) ==  0) {
          	   strcpy(prefix,"mxrf");
          	   strcpy(comm_var,"mxrf_comment");
          	   fctime = 2;
                }
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LMR  */
        case LF1:
           if (*imax == NX && *jmax == NY && *kmax <= 3 && *kdim <= 3)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lf1(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "H00",3) ==  0) {
          	   strcpy(prefix,"mxrh");
          	   strcpy(comm_var,"mxrh_comment");
          	   fctime = 0;
          	}
          	else if (strncmp(var[i], "H06",3) ==  0) {
          	   strcpy(prefix,"mxrh");
          	   strcpy(comm_var,"mxrh_comment");
          	   fctime = 1;
          	}
          	else if (strncmp(var[i], "H12",3) ==  0) {
          	   strcpy(prefix,"mxrh");
          	   strcpy(comm_var,"mxrh_comment");
          	   fctime = 2;
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LF1  */
        case L1S:
           if (*imax == NX && *jmax == NY && *kmax <= 4 && *kdim <= 4)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_l1s(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "S01",3) ==  0) {
          	   strcpy(prefix,"s1hr");
          	   strcpy(comm_var,"s1hr_comment");
          	}
          	else if (strncmp(var[i], "STO",3) ==  0) {
          	   strcpy(prefix,"stot");
          	   strcpy(comm_var,"stot_comment");
          	}
          	else if (strncmp(var[i], "R01",3) ==  0) {
          	   strcpy(prefix,"pc");
          	   strcpy(comm_var,"pc_comment");
          	}
          	else if (strncmp(var[i], "RTO",3) ==  0) {
          	   strcpy(prefix,"pt");
          	   strcpy(comm_var,"pt_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of L1S  */
        case LPS:
           if (*imax == NX && *jmax == NY && *kmax <= 21 && *kdim <= 21)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lps(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "REF",3) ==  0) {
          	   strcpy(prefix,"ref");
          	   strcpy(comm_var,"ref_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LPS  */
        case LRP:
           if (*imax == NX && *jmax == NY && *kmax <= 21 && *kdim <= 21)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lrp(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LRP",3) ==  0) {
          	   strcpy(prefix,"icg");
          	   strcpy(comm_var,"icg_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LRP  */
        case LBA:
           if (*imax == NX && *jmax == NY && *kmax <= 84 && *kdim <= 84)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lba(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "HT",3) ==  0) {
          	   strcpy(prefix,"bz");
          	   strcpy(comm_var,"bz_comment");
          	}
          	else if (strncmp(var[i], "U",3) ==  0) {
          	   strcpy(prefix,"bu");
          	   strcpy(comm_var,"bu_comment");
          	}
          	else if (strncmp(var[i], "V",3) ==  0) {
          	   strcpy(prefix,"bv");
          	   strcpy(comm_var,"bv_comment");
          	}
          	else if (strncmp(var[i], "OM",3) ==  0) {
          	   strcpy(prefix,"bw");
          	   strcpy(comm_var,"bw_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LBA  */
        case LC3:
           if (*imax == NX && *jmax == NY && *kmax <= 42 && *kdim <= 42)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lc3(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LC3",3) ==  0) {
          	   strcpy(prefix,"camt");
          	   strcpy(comm_var,"camt_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LC3  */
        case LWC:
           if (*imax == NX && *jmax == NY && *kmax <= 63 && *kdim <= 63)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lwc(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LWC",3) ==  0) {
          	   strcpy(prefix,"lwc");
          	   strcpy(comm_var,"lwc_comment");
          	}
          	else if (strncmp(var[i], "ICE",3) ==  0) {
          	   strcpy(prefix,"ice");
          	   strcpy(comm_var,"ice_comment");
          	}
          	else if (strncmp(var[i], "PCN",3) ==  0) {
          	   strcpy(prefix,"pcn");
          	   strcpy(comm_var,"pcn_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LWC  */
        case LIL:
           if (*imax == NX && *jmax == NY && *kmax == 1 && *kdim == 1)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lil(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LIL",3) ==  0) {
          	   strcpy(prefix,"ilw");
          	   strcpy(comm_var,"ilw_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LIL  */
        case LCB:
           if (*imax == NX && *jmax == NY && *kmax <= 3 && *kdim <= 3)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lcb(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LCB",3) ==  0) {
          	   strcpy(prefix,"cbas");
          	   strcpy(comm_var,"cbas_comment");
          	}
          	else if (strncmp(var[i], "LCT",3) ==  0) {
          	   strcpy(prefix,"ctop");
          	   strcpy(comm_var,"ctop_comment");
          	}
          	else if (strncmp(var[i], "CCE",3) ==  0) {
          	   strcpy(prefix,"cce");
          	   strcpy(comm_var,"cce_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LCB  */
        case LCT:
           if (*imax == NX && *jmax == NY && *kmax <= 3 && *kdim <= 3)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lct(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
                if ((strncmp(var[i], "PTY",3) ==  0) ||
                    (strncmp(var[i], "SPT",3) ==  0)) {
          	   strcpy(prefix,"spt");
          	   strcpy(comm_var,"spt_comment");
          	}
          	else if (strncmp(var[i], "PTT",3) ==  0) {
          	   strcpy(prefix,"ptt");
          	   strcpy(comm_var,"ptt_comment");
          	}
          	else if (strncmp(var[i], "SCT",3) ==  0) {
          	   strcpy(prefix,"sct");
          	   strcpy(comm_var,"sct_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LCT  */
        case LCV:
           if (*imax == NX && *jmax == NY && *kmax <= 3 && *kdim <= 3)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lcv(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LCV",3) ==  0) {
          	   strcpy(prefix,"ccov");
          	   strcpy(comm_var,"ccov_comment");
          	}
          	else if (strncmp(var[i], "CSC",3) ==  0) {
          	   strcpy(prefix,"csc");
          	   strcpy(comm_var,"csc_comment");
          	}
          	else if (strncmp(var[i], "CWT",3) ==  0) {
          	   strcpy(prefix,"cwt");
          	   strcpy(comm_var,"cwt_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LCV  */
        case LMD:
           if (*imax == NX && *jmax == NY && *kmax <= 21 && *kdim <= 21)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lmd(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var, "LMD",3) ==  0) {
          	   strcpy(prefix,"mcd");
          	   strcpy(comm_var,"mcd_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LMD  */
        case LCO:
           if (*imax == NX && *jmax == NY && *kmax <= 21 && *kdim <= 21)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lco(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if ((strncmp(var[i], "OM",3) ==  0) ||
          	    (strncmp(var[i], "COM",3) ==  0)) {
          	   strcpy(prefix,"cw");
          	   strcpy(comm_var,"cw_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LCO  */
        case LTY:
           if (*imax == NX && *jmax == NY && *kmax <= 42 && *kdim <= 42)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lty(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "PTY",3) ==  0) {
          	   strcpy(prefix,"ptyp");
          	   strcpy(comm_var,"ptyp_comment");
          	}
          	else if (strncmp(var[i], "CTY",3) ==  0) {
          	   strcpy(prefix,"ctyp");
          	   strcpy(comm_var,"ctyp_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LTY  */
        case LCP:
           if (*imax == NX && *jmax == NY && *kmax <= 21 && *kdim <= 21)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lcp(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LCP",3) ==  0) {
          	   strcpy(prefix,"ccpc");
          	   strcpy(comm_var,"ccpc_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LCP  */
        case LVD:
           if (*imax == NX && *jmax == NY && *kmax <= 14 && *kdim <= 14)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lvd(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "S8W",3) ==  0) {
          	   strcpy(prefix,"s8w");
          	   strcpy(comm_var,"s8w_comment");
          	}
          	else if (strncmp(var[i], "S8C",3) ==  0) {
          	   strcpy(prefix,"s8c");
          	   strcpy(comm_var,"s8c_comment");
          	}
          	else if (strncmp(var[i], "SVS",3) ==  0) {
          	   strcpy(prefix,"svs");
          	   strcpy(comm_var,"svs_comment");
          	}
          	else if (strncmp(var[i], "SVN",3) ==  0) {
          	   strcpy(prefix,"svn");
          	   strcpy(comm_var,"svn_comment");
          	}
          	else if (strncmp(var[i], "ALB",3) ==  0) {
          	   strcpy(prefix,"alb");
          	   strcpy(comm_var,"alb_comment");
          	}
          	else if (strncmp(var[i], "S3A",3) ==  0) {
          	   strcpy(prefix,"s3a");
          	   strcpy(comm_var,"s3a_comment");
          	}
          	else if (strncmp(var[i], "S3C",3) ==  0) {
          	   strcpy(prefix,"s3c");
          	   strcpy(comm_var,"s3c_comment");
          	}
          	else if (strncmp(var[i], "S4A",3) ==  0) {
          	   strcpy(prefix,"s4a");
          	   strcpy(comm_var,"s4a_comment");
          	}
          	else if (strncmp(var[i], "S4C",3) ==  0) {
          	   strcpy(prefix,"s4c");
          	   strcpy(comm_var,"s4c_comment");
          	}
          	else if (strncmp(var[i], "S5A",3) ==  0) {
          	   strcpy(prefix,"s5a");
          	   strcpy(comm_var,"s5a_comment");
          	}
          	else if (strncmp(var[i], "S5C",3) ==  0) {
          	   strcpy(prefix,"s5c");
          	   strcpy(comm_var,"s5c_comment");
          	}
          	else if (strncmp(var[i], "S8A",3) ==  0) {
          	   strcpy(prefix,"s8a");
          	   strcpy(comm_var,"s8a_comment");
          	}
          	else if (strncmp(var[i], "SCA",3) ==  0) {
          	   strcpy(prefix,"sca");
          	   strcpy(comm_var,"sca_comment");
          	}
          	else if (strncmp(var[i], "SCC",3) ==  0) {
          	   strcpy(prefix,"scc");
          	   strcpy(comm_var,"scc_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LVD  */
        case LVE:
           if (*imax == NX && *jmax == NY && *kmax <= 14 && *kdim <= 14)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lve(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "D8W",3) ==  0) {
          	   strcpy(prefix,"d8w");
          	   strcpy(comm_var,"d8w_comment");
          	}
          	else if (strncmp(var[i], "D8C",3) ==  0) {
          	   strcpy(prefix,"d8c");
          	   strcpy(comm_var,"d8c_comment");
          	}
          	else if (strncmp(var[i], "SVS",3) ==  0) {
          	   strcpy(prefix,"svs");
          	   strcpy(comm_var,"svs_comment");
          	}
          	else if (strncmp(var[i], "SVN",3) ==  0) {
          	   strcpy(prefix,"svn");
          	   strcpy(comm_var,"svn_comment");
          	}
          	else if (strncmp(var[i], "ALB",3) ==  0) {
          	   strcpy(prefix,"alb");
          	   strcpy(comm_var,"alb_comment");
          	}
          	else if (strncmp(var[i], "D3A",3) ==  0) {
          	   strcpy(prefix,"d3a");
          	   strcpy(comm_var,"d3a_comment");
          	}
          	else if (strncmp(var[i], "D3C",3) ==  0) {
          	   strcpy(prefix,"d3c");
          	   strcpy(comm_var,"d3c_comment");
          	}
          	else if (strncmp(var[i], "D4A",3) ==  0) {
          	   strcpy(prefix,"d4a");
          	   strcpy(comm_var,"d4a_comment");
          	}
          	else if (strncmp(var[i], "D4C",3) ==  0) {
          	   strcpy(prefix,"d4c");
          	   strcpy(comm_var,"d4c_comment");
          	}
          	else if (strncmp(var[i], "D5A",3) ==  0) {
          	   strcpy(prefix,"d5a");
          	   strcpy(comm_var,"d5a_comment");
          	}
          	else if (strncmp(var[i], "D5C",3) ==  0) {
          	   strcpy(prefix,"d5c");
          	   strcpy(comm_var,"d5c_comment");
          	}
          	else if (strncmp(var[i], "D8A",3) ==  0) {
          	   strcpy(prefix,"d8a");
          	   strcpy(comm_var,"d8a_comment");
          	}
          	else if (strncmp(var[i], "DCA",3) ==  0) {
          	   strcpy(prefix,"dca");
          	   strcpy(comm_var,"dca_comment");
          	}
          	else if (strncmp(var[i], "DCC",3) ==  0) {
          	   strcpy(prefix,"dcc");
          	   strcpy(comm_var,"dcc_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LVE  */
        case LMA:
           if (*imax == 11 && *jmax == 11 && *kmax <= 105 && *kdim <= 105)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lma(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "HT",3) ==  0) {
          	   strcpy(prefix,"maz");
          	   strcpy(comm_var,"maz_comment");
          	}
          	else if (strncmp(var[i], "T",3) ==  0) {
          	   strcpy(prefix,"mat");
          	   strcpy(comm_var,"mat_comment");
          	}
          	else if (strncmp(var[i], "RH",3) ==  0) {
          	   strcpy(prefix,"marh");
          	   strcpy(comm_var,"marh_comment");
          	}
          	else if (strncmp(var[i], "U",3) ==  0) {
          	   strcpy(prefix,"mau");
          	   strcpy(comm_var,"mau_comment");
          	}
          	else if (strncmp(var[i], "V",3) ==  0) {
          	   strcpy(prefix,"mav");
          	   strcpy(comm_var,"mav_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LMA  */
        case LMF:
           if (*imax == 11 && *jmax == 11 && *kmax <= 105 && *kdim <= 105)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lmf(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "HT",3) ==  0) {
          	   strcpy(prefix,"maz");
          	   strcpy(comm_var,"maz_comment");
          	}
          	else if (strncmp(var[i], "T",3) ==  0) {
          	   strcpy(prefix,"mat");
          	   strcpy(comm_var,"mat_comment");
          	}
          	else if (strncmp(var[i], "RH",3) ==  0) {
          	   strcpy(prefix,"marh");
          	   strcpy(comm_var,"marh_comment");
          	}
          	else if (strncmp(var[i], "U",3) ==  0) {
          	   strcpy(prefix,"mau");
          	   strcpy(comm_var,"mau_comment");
          	}
          	else if (strncmp(var[i], "V",3) ==  0) {
          	   strcpy(prefix,"mav");
          	   strcpy(comm_var,"mav_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LMF  */
        case Z02:
           if (*imax == NX && *jmax == NY && *kmax <= 5 && *kdim <= 5)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_z02(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "S8W",3) ==  0) {
          	   strcpy(prefix,"s8w");
          	   strcpy(comm_var,"s8w_comment");
          	}
          	else if (strncmp(var[i], "S8C",3) ==  0) {
          	   strcpy(prefix,"s8c");
          	   strcpy(comm_var,"s8c_comment");
          	}
          	else if (strncmp(var[i], "SVS",3) ==  0) {
          	   strcpy(prefix,"svs");
          	   strcpy(comm_var,"svs_comment");
          	}
          	else if (strncmp(var[i], "SVN",3) ==  0) {
          	   strcpy(prefix,"svn");
          	   strcpy(comm_var,"svn_comment");
          	}
          	else if (strncmp(var[i], "ALB",3) ==  0) {
          	   strcpy(prefix,"alb");
          	   strcpy(comm_var,"alb_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of Z02  */
        case RAM:
           if (*imax == NX && *jmax == NY && *kmax <= 273 && *kdim <= 273)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_ram(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "HT",3) ==  0) {
          	   strcpy(prefix,"rz");
          	   strcpy(comm_var,"rz_comment");
          	}
                else if ((strncmp(var[i], "U",3) ==  0) ||
                         (strncmp(var[i], "U3",3) ==  0)) {
          	   strcpy(prefix,"ru");
          	   strcpy(comm_var,"ru_comment");
          	}
                else if ((strncmp(var[i], "V",3) ==  0) ||
                         (strncmp(var[i], "V3",3) ==  0)) {
          	   strcpy(prefix,"rv");
          	   strcpy(comm_var,"rv_comment");
          	}
          	else if (strncmp(var[i], "OM",3) ==  0) {
          	   strcpy(prefix,"rw");
          	   strcpy(comm_var,"rw_comment");
          	}
                else if ((strncmp(var[i], "T",3) ==  0) ||
                         (strncmp(var[i], "T3",3) ==  0)) {
          	   strcpy(prefix,"rt");
          	   strcpy(comm_var,"rt_comment");
          	}
          	else if (strncmp(var[i], "RMR",3) ==  0) {
          	   strcpy(prefix,"rmr");
          	   strcpy(comm_var,"rmr_comment");
          	}
          	else if (strncmp(var[i], "RH3",3) ==  0) {
          	   strcpy(prefix,"rh3");
          	   strcpy(comm_var,"rh3_comment");
          	}
          	else if (strncmp(var[i], "LWC",3) ==  0) {
          	   strcpy(prefix,"lwc");
          	   strcpy(comm_var,"lwc_comment");
          	}
          	else if (strncmp(var[i], "ICE",3) ==  0) {
          	   strcpy(prefix,"ice");
          	   strcpy(comm_var,"ice_comment");
          	}
          	else if (strncmp(var[i], "RAI",3) ==  0) {
          	   strcpy(prefix,"rai");
          	   strcpy(comm_var,"rai_comment");
          	}
          	else if (strncmp(var[i], "SNO",3) ==  0) {
          	   strcpy(prefix,"sno");
          	   strcpy(comm_var,"sno_comment");
          	}
          	else if (strncmp(var[i], "PIC",3) ==  0) {
          	   strcpy(prefix,"pic");
          	   strcpy(comm_var,"pic_comment");
          	}
          	else if (strncmp(var[i], "REF",3) ==  0) {
          	   strcpy(prefix,"ref");
          	   strcpy(comm_var,"ref_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of RAM  */
        case RSF:
           if (*imax == NX && *jmax == NY && *kmax <= 28 && *kdim <= 28)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_rsf(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "U",3) ==  0) {
          	   strcpy(prefix,"rus");
          	   strcpy(comm_var,"rus_comment");
          	}
          	else if (strncmp(var[i], "V",3) ==  0) {
          	   strcpy(prefix,"rvs");
          	   strcpy(comm_var,"rvs_comment");
          	}
                else if ((strncmp(var[i], "P",3) ==  0) ||
                         (strncmp(var[i], "RP",3) ==  0)) {
          	   strcpy(prefix,"rps");
          	   strcpy(comm_var,"rps_comment");
          	}
          	else if (strncmp(var[i], "T",3) ==  0) {
          	   strcpy(prefix,"rts");
          	   strcpy(comm_var,"rts_comment");
          	}
          	else if (strncmp(var[i], "TD",3) ==  0) {
          	   strcpy(prefix,"rtd");
          	   strcpy(comm_var,"rtd_comment");
          	}
          	else if (strncmp(var[i], "RH",3) ==  0) {
          	   strcpy(prefix,"rh");
          	   strcpy(comm_var,"rh_comment");
          	}
          	else if (strncmp(var[i], "LCB",3) ==  0) {
          	   strcpy(prefix,"lcb");
          	   strcpy(comm_var,"lcb_comment");
          	}
          	else if (strncmp(var[i], "LCT",3) ==  0) {
          	   strcpy(prefix,"lct");
          	   strcpy(comm_var,"lct_comment");
          	}
          	else if (strncmp(var[i], "MSL",3) ==  0) {
          	   strcpy(prefix,"msl");
          	   strcpy(comm_var,"msl_comment");
          	}
          	else if (strncmp(var[i], "LIL",3) ==  0) {
          	   strcpy(prefix,"lil");
          	   strcpy(comm_var,"lil_comment");
          	}
          	else if (strncmp(var[i], "TPW",3) ==  0) {
          	   strcpy(prefix,"tpw");
          	   strcpy(comm_var,"tpw_comment");
          	}
          	else if (strncmp(var[i], "R01",3) ==  0) {
          	   strcpy(prefix,"r01");
          	   strcpy(comm_var,"r01_comment");
          	}
          	else if (strncmp(var[i], "RTO",3) ==  0) {
          	   strcpy(prefix,"rto");
          	   strcpy(comm_var,"rto_comment");
          	}
          	else if (strncmp(var[i], "S01",3) ==  0) {
          	   strcpy(prefix,"s01");
          	   strcpy(comm_var,"s01_comment");
          	}
          	else if (strncmp(var[i], "STO",3) ==  0) {
          	   strcpy(prefix,"sto");
          	   strcpy(comm_var,"sto_comment");
          	}
          	else if (strncmp(var[i], "TH",3) ==  0) {
          	   strcpy(prefix,"th");
          	   strcpy(comm_var,"th_comment");
          	}
          	else if (strncmp(var[i], "THE",3) ==  0) {
          	   strcpy(prefix,"the");
          	   strcpy(comm_var,"the_comment");
          	}
          	else if (strncmp(var[i], "PBE",3) ==  0) {
          	   strcpy(prefix,"pbe");
          	   strcpy(comm_var,"pbe_comment");
          	}
          	else if (strncmp(var[i], "NBE",3) ==  0) {
          	   strcpy(prefix,"nbe");
          	   strcpy(comm_var,"nbe_comment");
          	}
          	else if (strncmp(var[i], "PS",3) ==  0) {
          	   strcpy(prefix,"ps");
          	   strcpy(comm_var,"ps_comment");
          	}
          	else if (strncmp(var[i], "CCE",3) ==  0) {
          	   strcpy(prefix,"cce");
          	   strcpy(comm_var,"cce_comment");
          	}
          	else if (strncmp(var[i], "VIS",3) ==  0) {
          	   strcpy(prefix,"vis");
          	   strcpy(comm_var,"vis_comment");
          	}
          	else if (strncmp(var[i], "LCV",3) ==  0) {
          	   strcpy(prefix,"lcv");
          	   strcpy(comm_var,"lcv_comment");
          	}
          	else if (strncmp(var[i], "LMT",3) ==  0) {
          	   strcpy(prefix,"lmt");
          	   strcpy(comm_var,"lmt_comment");
          	}
          	else if (strncmp(var[i], "SPT",3) ==  0) {
          	   strcpy(prefix,"spt");
          	   strcpy(comm_var,"spt_comment");
          	}
          	else if (strncmp(var[i], "LHE",3) ==  0) {
          	   strcpy(prefix,"lhe");
          	   strcpy(comm_var,"lhe_comment");
          	}
          	else if (strncmp(var[i], "LI",3) ==  0) {
          	   strcpy(prefix,"li");
          	   strcpy(comm_var,"li_comment");
          	}
          	else if (strncmp(var[i], "HI",3) ==  0) {
          	   strcpy(prefix,"hi");
          	   strcpy(comm_var,"hi_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of RSF  */
        case RSM:
           if (*imax == NX && *jmax == NY && *kmax <= 11 && *kdim <= 11)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_rsm(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LSM",3) ==  0) {
          	   strcpy(prefix,"lsm");
          	   strcpy(comm_var,"lsm_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of RSM  */
        case VRC:
           if (*imax == NX && *jmax == NY && *kmax == 1 && *kdim == 1)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_vrc(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "REF",3) ==  0) {
          	   strcpy(prefix,"now");
          	   strcpy(comm_var,"now_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of VRC  */
        case LM1:
           if (*imax == NX && *jmax == NY && *kmax <= 3 && *kdim <= 3)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lm1(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "LSM",3) ==  0) {
          	   strcpy(prefix,"lsm");
          	   strcpy(comm_var,"lsm_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LM1  */
        case LM2:
           if (*imax == NX && *jmax == NY && *kmax <= 7 && *kdim <= 7)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lm2(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "CIV",3) ==  0) {
          	   strcpy(prefix,"civ");
          	   strcpy(comm_var,"civ_comment");
          	}
          	else if (strncmp(var[i], "DWF",3) ==  0) {
          	   strcpy(prefix,"dwf");
          	   strcpy(comm_var,"dwf_comment");
          	}
          	else if (strncmp(var[i], "WX",3) ==  0) {
          	   strcpy(prefix,"wx");
          	   strcpy(comm_var,"wx_comment");
          	}
          	else if (strncmp(var[i], "EVP",3) ==  0) {
          	   strcpy(prefix,"evp");
          	   strcpy(comm_var,"evp_comment");
          	}
          	else if (strncmp(var[i], "SC",3) ==  0) {
          	   strcpy(prefix,"sc");
          	   strcpy(comm_var,"sc_comment");
          	}
          	else if (strncmp(var[i], "SM",3) ==  0) {
          	   strcpy(prefix,"sm");
          	   strcpy(comm_var,"sm_comment");
          	}
          	else if (strncmp(var[i], "MWF",3) ==  0) {
          	   strcpy(prefix,"mwf");
          	   strcpy(comm_var,"mwf_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LM2  */
        case VRD:
           if (*imax == NX && *jmax == NY && *kmax <= (NZ*2) && *kdim <= (NZ*2))
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_vrd(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "REF",3) ==  0) {
          	   strcpy(prefix,"refd");
          	   strcpy(comm_var,"refd_comment");
          	}
          	else if (strncmp(var[i], "VEL",3) ==  0) {
          	   strcpy(prefix,"veld");
          	   strcpy(comm_var,"veld_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of VRD  */
        case V00:
           if (*imax == NX && *jmax == NY && *kmax <= (NZ*3) && *kdim <= (NZ*3))
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_v_radar(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "REF",3) ==  0) {
          	   strcpy(prefix,"refd");
          	   strcpy(comm_var,"refd_comment");
          	}
          	else if (strncmp(var[i], "VEL",3) ==  0) {
          	   strcpy(prefix,"veld");
          	   strcpy(comm_var,"veld_comment");
          	}
          	else if (strncmp(var[i], "NYQ",3) ==  0) {
          	   strcpy(prefix,"nyqd");
          	   strcpy(comm_var,"nyqd_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of V00  */
        case LGA:
           if (*imax == NX && *jmax == NY && *kmax <= 105 && *kdim <= 105)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lga(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "HT",3) ==  0) {
          	   strcpy(prefix,"mgz");
          	   strcpy(comm_var,"mgz_comment");
          	}
                else if ((strncmp(var[i], "T",3) ==  0) ||
                         (strncmp(var[i], "T3",3) ==  0)) {
          	   strcpy(prefix,"mgt");
          	   strcpy(comm_var,"mgt_comment");
          	}
          	else if (strncmp(var[i], "SH",3) ==  0) {
          	   strcpy(prefix,"mgrh");
          	   strcpy(comm_var,"mgrh_comment");
          	}
                else if ((strncmp(var[i], "U",3) ==  0) ||
                         (strncmp(var[i], "U3",3) ==  0)) {
          	   strcpy(prefix,"mgu");
          	   strcpy(comm_var,"mgu_comment");
          	}
                else if ((strncmp(var[i], "V",3) ==  0) ||
                         (strncmp(var[i], "V3",3) ==  0)) {
          	   strcpy(prefix,"mgv");
          	   strcpy(comm_var,"mgv_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LGA  */
        case LGF:
           if (*imax == NX && *jmax == NY && *kmax <= 105 && *kdim <= 105)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_lgf(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "HT",3) ==  0) {
          	   strcpy(prefix,"mfz");
          	   strcpy(comm_var,"mfz_comment");
          	}
                else if ((strncmp(var[i], "T",3) ==  0) ||
                         (strncmp(var[i], "T3",3) ==  0)) {
          	   strcpy(prefix,"mft");
          	   strcpy(comm_var,"mft_comment");
          	}
          	else if (strncmp(var[i], "SH",3) ==  0) {
          	   strcpy(prefix,"mfrh");
          	   strcpy(comm_var,"mfrh_comment");
          	}
                else if ((strncmp(var[i], "U",3) ==  0) ||
                         (strncmp(var[i], "U3",3) ==  0)) {
          	   strcpy(prefix,"mfu");
          	   strcpy(comm_var,"mfu_comment");
          	}
                else if ((strncmp(var[i], "V",3) ==  0) ||
                         (strncmp(var[i], "V3",3) ==  0)) {
          	   strcpy(prefix,"mfv");
          	   strcpy(comm_var,"mfv_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LGF  */
        case LN3:
           if (*imax == NX && *jmax == NY && *kmax <= 6 && *kdim <= 6)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_ln3(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "R04",3) ==  0) {
          	   strcpy(prefix,"r04");
          	   strcpy(comm_var,"r04_comment");
          	}
          	else if (strncmp(var[i], "R48",3) ==  0) {
          	   strcpy(prefix,"r48");
          	   strcpy(comm_var,"r48_comment");
          	}
          	else if (strncmp(var[i], "R8C",3) ==  0) {
          	   strcpy(prefix,"r8c");
          	   strcpy(comm_var,"r8c_comment");
          	}
          	else if (strncmp(var[i], "ET",3) ==  0) {
          	   strcpy(prefix,"et");
          	   strcpy(comm_var,"et_comment");
          	}
          	else if (strncmp(var[i], "RCO",3) ==  0) {
          	   strcpy(prefix,"rco");
          	   strcpy(comm_var,"rco_comment");
          	}
          	else if (strncmp(var[i], "VIL",3) ==  0) {
          	   strcpy(prefix,"vil");
          	   strcpy(comm_var,"vil_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LN3  */
        case MM5:
           if (*imax == NX && *jmax == NY && *kmax <= 189 && *kdim <= 189)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_mm5(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "HT",3) ==  0) {
          	   strcpy(prefix,"rz");
          	   strcpy(comm_var,"rz_comment");
          	}
                else if ((strncmp(var[i], "U",3) ==  0) ||
                         (strncmp(var[i], "U3",3) ==  0)) {
          	   strcpy(prefix,"ru");
          	   strcpy(comm_var,"ru_comment");
          	}
                else if ((strncmp(var[i], "V",3) ==  0) ||
                         (strncmp(var[i], "V3",3) ==  0)) {
          	   strcpy(prefix,"rv");
          	   strcpy(comm_var,"rv_comment");
          	}
          	else if (strncmp(var[i], "OM",3) ==  0) {
          	   strcpy(prefix,"rw");
          	   strcpy(comm_var,"rw_comment");
          	}
                else if ((strncmp(var[i], "T",3) ==  0) ||
                         (strncmp(var[i], "T3",3) ==  0)) {
          	   strcpy(prefix,"rt");
          	   strcpy(comm_var,"rt_comment");
          	}
          	else if (strncmp(var[i], "SH",3) ==  0) {
          	   strcpy(prefix,"rsh");
          	   strcpy(comm_var,"rsh_comment");
          	}
          	else if (strncmp(var[i], "RH3",3) ==  0) {
          	   strcpy(prefix,"rh3");
          	   strcpy(comm_var,"rh3_comment");
          	}
          	else if (strncmp(var[i], "LWC",3) ==  0) {
          	   strcpy(prefix,"lwc");
          	   strcpy(comm_var,"lwc_comment");
          	}
          	else if (strncmp(var[i], "ICE",3) ==  0) {
          	   strcpy(prefix,"ice");
          	   strcpy(comm_var,"ice_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of MM5  */
        case MSF:
           if (*imax == NX && *jmax == NY && *kmax <= 28 && *kdim <= 28)
               {}
           else {
           	*status = -3;	/* returns dimension error */
           	return;
           }
           out_file = cre_msf(fname);
           if (out_file == -1) {
           	*status = -2;	/* returns "error opening file" */
           	return;
           }
          
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
          	*status = -5;	/* returns error writing header info */
          	ncclose(out_file);          	
          	return;
           }

	   /* write data to grid and write out comment*/
	   count = 0;
           for (i = 0; i < *kdim; i++) {
           	process_vr = 1;
          	if (strncmp(var[i], "U",3) ==  0) {
          	   strcpy(prefix,"rus");
          	   strcpy(comm_var,"rus_comment");
          	}
          	else if (strncmp(var[i], "V",3) ==  0) {
          	   strcpy(prefix,"rvs");
          	   strcpy(comm_var,"rvs_comment");
          	}
                else if ((strncmp(var[i], "P",3) ==  0) ||
                         (strncmp(var[i], "RP",3) ==  0)) {
          	   strcpy(prefix,"rps");
          	   strcpy(comm_var,"rps_comment");
          	}
          	else if (strncmp(var[i], "T",3) ==  0) {
          	   strcpy(prefix,"rts");
          	   strcpy(comm_var,"rts_comment");
          	}
          	else if (strncmp(var[i], "TD",3) ==  0) {
          	   strcpy(prefix,"rtd");
          	   strcpy(comm_var,"rtd_comment");
          	}
          	else if (strncmp(var[i], "RH",3) ==  0) {
          	   strcpy(prefix,"rh");
          	   strcpy(comm_var,"rh_comment");
          	}
          	else if (strncmp(var[i], "LCB",3) ==  0) {
          	   strcpy(prefix,"lcb");
          	   strcpy(comm_var,"lcb_comment");
          	}
          	else if (strncmp(var[i], "LCT",3) ==  0) {
          	   strcpy(prefix,"lct");
          	   strcpy(comm_var,"lct_comment");
          	}
          	else if (strncmp(var[i], "MSL",3) ==  0) {
          	   strcpy(prefix,"msl");
          	   strcpy(comm_var,"msl_comment");
          	}
          	else if (strncmp(var[i], "LIL",3) ==  0) {
          	   strcpy(prefix,"lil");
          	   strcpy(comm_var,"lil_comment");
          	}
          	else if (strncmp(var[i], "TPW",3) ==  0) {
          	   strcpy(prefix,"tpw");
          	   strcpy(comm_var,"tpw_comment");
          	}
          	else if (strncmp(var[i], "R01",3) ==  0) {
          	   strcpy(prefix,"r01");
          	   strcpy(comm_var,"r01_comment");
          	}
          	else if (strncmp(var[i], "RTO",3) ==  0) {
          	   strcpy(prefix,"rto");
          	   strcpy(comm_var,"rto_comment");
          	}
          	else if (strncmp(var[i], "S01",3) ==  0) {
          	   strcpy(prefix,"s01");
          	   strcpy(comm_var,"s01_comment");
          	}
          	else if (strncmp(var[i], "STO",3) ==  0) {
          	   strcpy(prefix,"sto");
          	   strcpy(comm_var,"sto_comment");
          	}
          	else if (strncmp(var[i], "TH",3) ==  0) {
          	   strcpy(prefix,"th");
          	   strcpy(comm_var,"th_comment");
          	}
          	else if (strncmp(var[i], "THE",3) ==  0) {
          	   strcpy(prefix,"the");
          	   strcpy(comm_var,"the_comment");
          	}
          	else if (strncmp(var[i], "PBE",3) ==  0) {
          	   strcpy(prefix,"pbe");
          	   strcpy(comm_var,"pbe_comment");
          	}
          	else if (strncmp(var[i], "NBE",3) ==  0) {
          	   strcpy(prefix,"nbe");
          	   strcpy(comm_var,"nbe_comment");
          	}
          	else if (strncmp(var[i], "PS",3) ==  0) {
          	   strcpy(prefix,"ps");
          	   strcpy(comm_var,"ps_comment");
          	}
          	else if (strncmp(var[i], "CCE",3) ==  0) {
          	   strcpy(prefix,"cce");
          	   strcpy(comm_var,"cce_comment");
          	}
          	else if (strncmp(var[i], "VIS",3) ==  0) {
          	   strcpy(prefix,"vis");
          	   strcpy(comm_var,"vis_comment");
          	}
          	else if (strncmp(var[i], "LCV",3) ==  0) {
          	   strcpy(prefix,"lcv");
          	   strcpy(comm_var,"lcv_comment");
          	}
          	else if (strncmp(var[i], "LMT",3) ==  0) {
          	   strcpy(prefix,"lmt");
          	   strcpy(comm_var,"lmt_comment");
          	}
          	else if (strncmp(var[i], "SPT",3) ==  0) {
          	   strcpy(prefix,"spt");
          	   strcpy(comm_var,"spt_comment");
          	}
          	else if (strncmp(var[i], "LHE",3) ==  0) {
          	   strcpy(prefix,"lhe");
          	   strcpy(comm_var,"lhe_comment");
          	}
          	else if (strncmp(var[i], "LI",3) ==  0) {
          	   strcpy(prefix,"li");
          	   strcpy(comm_var,"li_comment");
          	}
          	else if (strncmp(var[i], "HI",3) ==  0) {
          	   strcpy(prefix,"hi");
          	   strcpy(comm_var,"hi_comment");
          	}
          	else  {
          	  process_vr = 0;
          	  count += 1;
          	}
          	if (process_vr == 1) {
          	   i_level = (int)lvl[i];
          	   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                	 		   (data + (i*(*imax)*(*jmax))),
                	 		   comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
          	      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of MSF  */
        
        case SST:
           if (*imax == NX && *jmax == NY && *kmax <= 5 && *kdim <= 5)
               {}
           else {
                *status = -3;   /* returns dimension error */
                return;
           }
           out_file = cre_sst(fname);
           if (out_file == -1) {
                *status = -2;   /* returns "error opening file" */
                return;
           }
         
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
                *status = -5;   /* returns error writing header info */
                ncclose(out_file);
                return;
           }

           /* write data to grid and write out comment*/
           count = 0;
           for (i = 0; i < *kdim; i++) {
                process_vr = 1;
                if (strncmp(var[i], "ST1",3) ==  0) {
                   strcpy(prefix,"st1");
                   strcpy(comm_var,"st1_comment");
                }
                else if (strncmp(var[i], "ST2",3) ==  0) {
                   strcpy(prefix,"st2");
                   strcpy(comm_var,"st2_comment");
                }
                else if (strncmp(var[i], "ST3",3) ==  0) {
                   strcpy(prefix,"st3");
                   strcpy(comm_var,"st3_comment");
                }
                else if (strncmp(var[i], "ST4",3) ==  0) {
                   strcpy(prefix,"st4");
                   strcpy(comm_var,"st4_comment");
                }
                else if (strncmp(var[i], "ST5",3) ==  0) {
                   strcpy(prefix,"st5");
                   strcpy(comm_var,"st5_comment");
                }
                else  {
                  process_vr = 0;
                  count += 1;
                }
                if (process_vr == 1) {
                   i_level = (int)lvl[i];
                   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                                           (data + (i*(*imax)*(*jmax))),
                                           comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
                      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of SST  */
        case VEG:
           if (*imax == NX && *jmax == NY && *kmax <= 5 && *kdim <= 5)
               {}
           else {
                *status = -3;   /* returns dimension error */
                return;
           }
           out_file = cre_veg(fname);
           if (out_file == -1) {
                *status = -2;   /* returns "error opening file" */
                return;
           }
         
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
                *status = -5;   /* returns error writing header info */
                ncclose(out_file);
                return;
           }

           /* write data to grid and write out comment*/
           count = 0;
           for (i = 0; i < *kdim; i++) {
                process_vr = 1;
                if (strncmp(var[i], "VG1",3) ==  0) {
                   strcpy(prefix,"vg1");
                   strcpy(comm_var,"vg1_comment");
                }
                else if (strncmp(var[i], "VG2",3) ==  0) {
                   strcpy(prefix,"vg2");
                   strcpy(comm_var,"vg2_comment");
                }
                else if (strncmp(var[i], "VG3",3) ==  0) {
                   strcpy(prefix,"vg3");
                   strcpy(comm_var,"vg3_comment");
                }
                else if (strncmp(var[i], "VG4",3) ==  0) {
                   strcpy(prefix,"vg4");
                   strcpy(comm_var,"vg4_comment");
                }
                else if (strncmp(var[i], "VG5",3) ==  0) {
                   strcpy(prefix,"vg5");
                   strcpy(comm_var,"vg5_comment");
                }
                else  {
                  process_vr = 0;
                  count += 1;
                }
                if (process_vr == 1) {
                   i_level = (int)lvl[i];
                   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                                           (data + (i*(*imax)*(*jmax))),
                                           comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
                      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of VEG  */
        case LS8:
           if (*imax == NX && *jmax == NY && *kmax <= 19 && *kdim <= 19)
               {}
           else {
                *status = -3;   /* returns dimension error */
                return;
           }
           out_file = cre_ls8(fname);
           if (out_file == -1) {
                *status = -2;   /* returns "error opening file" */
                return;
           }
         
           istat = write_hdr_cdf(out_file,kdim,lvl,laps_dom_file,asctime);
           if (istat == -1) {
                *status = -5;   /* returns error writing header info */
                ncclose(out_file);
                return;
           }

           /* write data to grid and write out comment*/
           count = 0;
           process_vr = 1;
           for (i = 0; i < *kdim; i++) {
                strncpy(prefix, var[i], 3); 
                if (isupper(prefix[1]) == 0) strncpy(prefix, "s", 1);
                strcpy(comm_var,prefix);
                strcat(comm_var,"_comment");

                if (process_vr == 1) {
                   i_level = (int)lvl[i];
                   fctime = 0;
                   istat = cdf_update_laps(out_file,i_level,fctime,prefix,
                                           (data + (i*(*imax)*(*jmax))),
                                           comm_var,comment[i]);
                   if (istat == -1) {
                      *status = -4;
                      ncclose(out_file);
                      return;
                   }
                }
           }
           break;  /*  end of LS8  */

        }	/* end of switch */
        *status = count;	/* normal return */

/*        strcpy(model,"LAPS - Local Analysis and Prediction System");
        dims[0] = 132;
        vid =  ncvarid (out_file, "model", NC_CHAR, 1, dims);
        start[0] = 0;
        edges[0] = strlen(model);
        ncvarput(out_file, vid, start, edges, (void *)model);

        strcpy(origin,"NOAA/ERL/Forecast Systems Laboratory, Boulder, CO");
        dims[0] = 132;
        vid =  ncvarid (out_file, "origin", NC_CHAR, 1, dims);
        start[0] = 0;
        edges[0] = strlen(origin);
        ncvarput(out_file, vid, start, edges, (void *)origin);
*/

        ncclose(out_file);          	
        return;
}

/*************************************************************************
*	READ_CDF_FILE
*	Category	Product Management
*	Group		General Purpose Database
*	Module		write_cdf-Read_cdf_file
*	Purpose		Read LAPS data from a netCDF format file.
*
*	Designer/Programmer : Linda Wharton
*	Modifications : original 1/93
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
#if defined(__alpha) || defined(__IP21)
void read_cdf_file (int *iimax,int *jjmax,int *kkmax,int *kdim,
		    int *imax, int *jmax, int *kmax,int *num_variables,
		    char *var_req, int *lvl_req,
		    char *filname, short *s_length,
		    char *f_ext, char *f_lvl_coord,char *f_units,char *f_comment,
		    char *f_laps_dom_file,char *f_asctime,char *version,
		    char *model, char *origin,float *data,
		    int *process_var,int *no_laps_diag,int *status)
#else
void read_cdf_file (long *iimax,long *jjmax,long *kkmax,long *kdim,
		    long *imax, long *jmax, long *kmax,long *num_variables,
		    char *var_req,long *lvl_req,
		    char *filname, short *s_length,
		    char *f_ext, char *f_lvl_coord,char *f_units,char *f_comment,
		    char *f_laps_dom_file,char *f_asctime,char *version,
		    char *model, char *origin,float *data,
		    long *process_var,long *no_laps_diag,long *status)
#endif
#else
#if defined(__alpha) || defined(__IP21)
void read_cdf_file (iimax,jjmax,kkmax,kdim,imax,jmax,kmax,num_variables,
		    var_req,lvl_req,filname,s_length,f_ext,f_lvl_coord,f_units,
		    f_comment,f_laps_dom_file,f_asctime,version,model,origin,
		    data,process_var,no_laps_diag,status)
int *iimax;
int *jjmax;
int *kkmax;
int *kdim;
int *imax;
int *jmax;
int *kmax;
int *num_variables;
char *var_req;
int *lvl_req;
char *filname;
short *s_length;
char *f_ext;
char *f_lvl_coord;
char *f_units;
char *f_comment;
char *f_laps_dom_file;
char *f_asctime;
char *version;
char *model;
char *origin;
float *data;
int *process_var;
int *no_laps_diag;
int *status;
#else
void read_cdf_file (iimax,jjmax,kkmax,kdim,imax,jmax,kmax,num_variables,
		    var_req,lvl_req,filname,s_length,f_ext,f_lvl_coord,f_units,
		    f_comment,f_laps_dom_file,f_asctime,version,model,origin,
		    data,process_var,no_laps_diag,status)
long *iimax;
long *jjmax;
long *kkmax;
long *kdim;
long *imax;
long *jmax;
long *kmax;
long *num_variables;
char *var_req;
long *lvl_req;
char *filname;
short *s_length;
char *f_ext;
char *f_lvl_coord;
char *f_units;
char *f_comment;
char *f_laps_dom_file;
char *f_asctime;
char *version;
char *model;
char *origin;
float *data;
long *process_var;
long *no_laps_diag;
long *status;
#endif
#endif
{		   
	typedef char STRING5[5];
	char cdf_var[300][5],c_var[300][4],ext[32];
        char lvl_coord[300][5],units[300][11],comment[300][126];
 	char in_var[5],laps_dom_file[12],asctime[18];
	static short cdf_fctime[300];
	int out_file, istat, i,j, i_level, i_fctime, unconv_var, cdfid;
	char fname[92],*cpt;
	
/* null out arrays for lvl_coord,units,comment,laps_dom_file,asctime,
     version,model, and origin before using   */

        for (i = 0; i < *kdim; i++) {
          strncpy(lvl_coord[i],"    ",4);
          strncpy(units[i],"          ",10);
          strncpy(comment[i],"                                                                                                                                 ",125);
        } 
        strncpy(laps_dom_file,"           ",11);
        strncpy(asctime,"                 ",17);
      
/* convert fortran string var_req to c string c_var */

        for (i = 0; i < *kdim; i++)
          nstrncpy(c_var[i],(var_req+i*3),3);

/* convert fortran file_name into C fname, and fortran f_ext into C ext  */
        nstrncpy(fname,filname,s_length);
        nstrncpy(ext,f_ext,31);


/* convert LAPS variables to corresponding netCDF variables--unconv_var
   is a count of how many variables in c_var could not be converted */

        unconv_var = get_cdf_var(ext,kkmax,c_var,cdf_var,process_var,
        			 cdf_fctime);
	
/* open netCDF file for reading */
	cdfid = open_cdf(NC_NOWRITE,fname,no_laps_diag);
	if (cdfid == -1) {
		*status = -1;	/* error opening file */
		return;
	}

	istat = cdf_retrieve_hdr(cdfid,imax,jmax,kmax,laps_dom_file,
	                         asctime,version,model,origin,
	                         num_variables);
	if (istat == -1) {
	   *status = -2;
	   ncclose(cdfid);
	   return;
	}
	if (*imax > *iimax || *jmax > *jjmax || *kkmax > *kdim) {
	   *status = -3;
	   ncclose(cdfid);
	   return;
	}
	
	for (i = 0; i < *kkmax; i++) {
	   i_level = lvl_req[i];
	   i_fctime = cdf_fctime[i];
	   strncpy(&in_var,&cdf_var[i][0],5);
	   istat = cdf_retrieve_laps_grid(cdfid,i_level,i_fctime,in_var,
	   		                  (data + i*(*iimax)*(*jjmax)),
                                          comment[i],lvl_coord[i],
					  units[i],ext);

				        /*  (comment + i*126),
				          (lvl_coord + i*5),
				          (units + i*11),ext);*/
	   if (istat == -1)  {
	      if ( process_var[i] != 0) {
	         process_var[i] = 0;
	         unconv_var += 1;
	      }
	   }
		
	}
	
	ncclose(cdfid);
	fill_empty_grids(iimax,jjmax,kkmax,process_var,data);  
	
        for (i = 0; i < *kdim; i++) {
          fstrncpy((f_comment+(i*126)),comment[i],125);
          fstrncpy((f_lvl_coord+(i*5)),lvl_coord[i],10);
          fstrncpy((f_units+(i*11)),units[i],4);
        }
        fstrncpy(f_laps_dom_file,laps_dom_file,11);
        fstrncpy(f_asctime,asctime,17);

	*status = unconv_var;
	return;
}

/*****************************************************************************/
#ifdef __STDC__
int cdf_inquire_var(int cdfid, long *n_var,char *var_avail, 
		    char *LAPS_var_avail,long *num_levels,
		    short *lvl_avail)
#else
int cdf_inquire_var(cdfid,n_var,var_avail,LAPS_var_avail,
		    num_levels,lvl_avail)
int cdfid;
long *n_var;
char *var_avail; 
char *LAPS_var_avail;
long *num_levels;
short *lvl_avail;
#endif
{
	nc_type datatype;
	int ndims, dim[6], natts, var_no, err_ct = {0}, istat, lv, i;
	unsigned char *d_ptr;
	
	for (var_no = 0; var_no < *n_var; var_no++)  {
	   istat = ncvarinq(cdfid, var_no, (var_avail + (var_no*20)),
	      		    &datatype, &ndims, dim, &natts);
	   if (natts > 3) 
	      istat = ncattget(cdfid, var_no, "LAPS_var", 
	   		       (LAPS_var_avail + (var_no*4))); 
	   if (istat == -1) err_ct += 1;
	}
	
/* read in the dimension id and size from the netcdf file */
	if ((*num_levels = cdf_dim_size (cdfid, "level")) == (-1))
	   err_ct += 1;
	else {

/* read the contents of the coordinate variable associated with this dimension
   from the net_cdf file */
	   istat = cdf_get_levels(cdfid,"level",*num_levels,lvl_avail);
	   if (istat == (-1)) {
	      err_ct += 1;
	      return err_ct;
	   }
	      
	}
	return err_ct;
}
/*****************************************************************************/
#ifdef __STDC__
void read_cdf_header(char *filname, short *s_length,
		     long *imax, long *jmax, 
                     long *kmax, long *num_variables, 
                     char *laps_dom_file, char *asctime, 
                     char *version, char *model, 
                     char *origin, char *var_avail, 
                     char *LAPS_var_avail, long *num_levels, 
                     short *lvl_avail,long *no_laps_diag,long *status)
#else
void read_cdf_header(filname,s_length,imax,jmax,kmax,num_variables,
		     laps_dom_file,asctime,version,model,
		     origin,var_avail,LAPS_var_avail,num_levels, 
                     lvl_avail,no_laps_diag,status)
char *filname;
short *s_length;
long *imax;
long *jmax;
long *kmax;
long *num_variables;
char *laps_dom_file;
char *asctime;
char *version;
char *model;
char *origin;
char *var_avail;
char *LAPS_var_avail;
long *num_levels;
short *lvl_avail;
long *no_laps_diag;
long *status;
#endif
{
	int istat, i, cdfid;
	char fname[92];

/* open netCDF file for reading */

	for (i = 0; i < 92; i++)
	   fname[i] = '\0';
	make_c_fname(filname,s_length,fname);
		
	cdfid = open_cdf(NC_NOWRITE,fname,no_laps_diag);
	if (cdfid == -1) {
		*status = -1;	/* error opening file */
		return;
	}
	   
	istat = cdf_retrieve_hdr(cdfid,imax,jmax,kmax,laps_dom_file,
	                         asctime,version,model,origin,
	                         num_variables);
	if (istat == -1) {
	   *status = -2;
	   ncclose(cdfid);
	   return;
	}
	
	istat =  cdf_inquire_var(cdfid,num_variables,var_avail,
				 LAPS_var_avail,num_levels,
				 lvl_avail);
	if (istat != 0) {
	   *status = -3;
	   ncclose(cdfid);
	   return;
	}
	
/* normal return */	
	*status = 1;
	ncclose(cdfid);
	return;
}
/*****************************************************************************/
