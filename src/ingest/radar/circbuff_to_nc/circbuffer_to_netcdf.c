/*  remap.c
*     
*     ##################################################################
*     ##################################################################
*     ######                                                      ######
*     ######                      REMAP.C                         ######
*     ######                                                      ######
*     ######    Center for Analysis and Prediction of Storms      ######
*     ######                University of Oklahoma                ######
*     ######             NOAA/ERL/Forecast Systems Lab            ######
*     ######               Copyright (c) 1994-1995                ######
*     ######                 All rights reserved                  ######
*     ######                                                      ######
*     ##################################################################
*     ##################################################################
*
*
*   Reads WSR 88D radar data from circular buffer using NSSL 
*   Realtime/Archive II radar reading routines and passes the data
*   to FSL LAPS remapping routines.
*
*   Usage:
*     remap                       reads data from live data
*     remap -f /dev/rmt0          reads data from tape drive rmt0
*
*   Advertised to work for live WSR-88D radar data as well
*   as archive II taped data.
*
*   Requires NSSL libraries:  -la2rt  for real-time
*                             -la2tp  and -ltpio for tape playback
*
*   ORIGINAL: 
*      Steve Albers, FSL and Keith Brewster, CAPS       May 1994
*   MODIFICATIONS:
*      Debug and overhaul  Keith Brewster, CAPS         June 1994
*   ADDITION
*      Network Program  Abraham Kappamammmootil, CAPS   March, 1995
*   MODIFICATIONS:
*      Mods to network program and update to new a2io release
*      Bells n' whistles added.
*                                Keith Brewster, CAPS   March 1995 
*      Call to pclose added after compression to be sure its done
*      before starting xmission process.  LAPS_XMIT added and
*      compression and transmission broken into two subroutines.
*                                Keith Brewster, CAPS   April 1995 
*      Modified logic to gaurantee call to fill_common when n_rays
*      reaches N_RAY_TILT.  Cleaned up if block indentation.
*                                Keith Brewster, CAPS   May 1995 
*/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>

/* When VERBOSE is defined, the program generates verbose output
*  that is useful for testing but annoying for operations */

/* #define VERBOSE */

/*  Files may be compressed using the program identified by LAPS_COMPRESS.
*  Currently supported LAPS_COMPRESS options are:  gzip, compress.
*
*  When the proper environment variables are defined
*  the program will try to send the remapped files 
*  to directory LAPS_DEST on machine LAPS_REMOTE owned by LAPS_USER,
*  using transmission mode LAPS_XMIT.
*
*  currently supported LAPS_XMIT methods are
*      rcp, ftp, mv
*
*  The LAPS_ variables are environment variables of the parent shell.
*
*  Sample:
*  setenv LAPS_COMPRESS gzip
*  setenv LAPS_XMIT rcp
*  setenv LAPS_USER kbrews
*  setenv LAPS_REMOTE cirrus.gcn.uoknor.edu
*  setenv LAPS_DEST /scratch/cirrus/kbrews/lapsprd/vrd
*/

/*  Note:
*   The radar name is determined from the environment variable
*   RADARNAME
*
*  Sample:
*  setenv RADARNAME KTLX
*
*  A file containing the lat,lon and elevation of the radar sites
*  should be in the directory of remap execution.
*  This file is named radarinfo.dat.
*/

/* The following defines the circular buffer to use
*  "RAW" is RAW
*  "E" is unfolded  
*/

#define RT_BUFFER "RAW"

#define NUM_REF_GATES 460
#define NUM_VEL_GATES 920
#define N_REF_TILT 174800
#define N_VEL_TILT 349600
#define N_RAY_TILT 380
#define NSIZE 920

/*  The following are for comparison to return
    values from the radar data access libraries */

#define GOOD_STATUS 1
#define BAD_STATUS 0
#define MAX_BAD_STAT 1000

#define TRUE 1
#define FALSE 0
#define MISSING -999.0

#define PLAYBACK 0
#define REALTIME 1

#ifdef FORTRANUNDERSCORE
#define lut_gen lut_gen_
#define fill_common fill_common_
#define remap_process remap_process_
#define int_to_i4time int_to_i4time_
#define make_fnam_lp make_fnam_lp_
#endif

#define allocate(s) (s= (char *) malloc(80))

/*    compr_file_sub controls 
      compression of the remapped data file             */

compr_file_sub(file_ptr,compr_ptr,i_status)
char  *file_ptr,*compr_ptr;
int i_status;
{
  char *string1;  
  FILE *comp_file;  
  int comp_stat;
  int slash;

  allocate(string1);

  if(strcmp(compr_ptr,"gzip") == 0)
  {
    strcpy(string1,"gzip -f ");
    strcat(string1,file_ptr);  
    comp_file= popen(string1,"r"); 
    if(comp_file == (FILE *) NULL)
      printf("\n Error occurred starting compression..\n");
    else
      comp_stat=pclose(comp_file); 
      if( comp_stat == 0) {
        strcat(file_ptr,".gz");  
        printf("\n File Successfully Compressed\n");
      }
      else {
        printf("\n Error occured while compressing %s\n",file_ptr);
        printf("\n Returned status = %i\n",comp_stat);
      }
  }
  else if(strcmp(compr_ptr,"compress") == 0)
  {
    strcpy(string1,"compress -f ");
    strcat(string1,file_ptr);  
    comp_file= popen(string1,"r"); 
    strcat(file_ptr,".Z");  
    fclose(comp_file);
  }
  else
  {
    printf(" %s Compression Unsupported\n",compr_ptr);
    printf(" Set LAPS_COMPRESS to gzip or compress\n");
  }

}
/*    send_file_sub controls 
      transmission of the remapped data file             */

send_file_sub(file_ptr,xmit_ptr,
              user_ptr,host_ptr,dest_ptr,i_status)
char  *file_ptr,*xmit_ptr,*user_ptr,*host_ptr,*dest_ptr;
int i_status;
{
  char *string1;  
  FILE *p1;
  int xmit_stat;
  int slash;
  char c_slash;
  char buf[200];

  allocate(string1);

/* Strip any directory info from the file name so all
*  directory info in the destination file name is from
*  dest_ptr, environment variable LAPS_DEST.
*  The source file name is unchanged */

/*  c_slash='/'; slash=atoi(c_slash);  */
  slash=47;
  string1=(char *)strrchr(file_ptr,slash);

  if (strncmp(xmit_ptr,"mv",2) == 0) {

/* Create mv command */

     if (string1 == NULL)
        sprintf(buf,"mv %s %s/%s",file_ptr,
             dest_ptr,file_ptr);
     else
        sprintf(buf,"mv %s %s%s",file_ptr,
             dest_ptr,string1);

     printf("\n Renaming file.....\n");

  }
  else if (strncmp(xmit_ptr,"rcp",3) == 0) {

/* Create the rcp command */

     if (string1 == NULL)
        sprintf(buf,"rcp %s %s@%s:%s/%s",file_ptr,user_ptr,host_ptr,
             dest_ptr,file_ptr);
     else
        sprintf(buf,"rcp %s %s@%s:%s%s",file_ptr,user_ptr,host_ptr,
             dest_ptr,string1);

     printf("\n Copying file to remote host.....\n");

   }
   else if (xmit_ptr != NULL) {

/* Create the custom user command

   For this option, a command line is created from the
   LAPS environment variables and the long and short filenames
   The LAPS_XMIT variable (represented here by xmit_ptr) is the
   name of a user-programmed script or executable that accepts
   the arguments as passed here.  These arguments should give
   the script or program enough info to work with. 

   The command line created consists of:
   LAPS_XMIT long_filename short_filename LAPS_USER LAPS_REMOTE LAPS_DEST

 */

     if (string1 == NULL)
        sprintf(buf,"%s %s %s %s %s %s",xmit_ptr,file_ptr,file_ptr,
                user_ptr,host_ptr,dest_ptr);
     else
        sprintf(buf,"%s %s %s %s %s %s",xmit_ptr,file_ptr,string1,
                user_ptr,host_ptr,dest_ptr);

     printf("\n Issuing command to user-defined data management program...\n");

   }

  fprintf(stdout,buf,"%s","\n");

  p1= popen(buf,"r");  

  if (p1 == (FILE *) NULL)
    printf("\n Error occurred starting file transmission..\n");
  else {
    xmit_stat = pclose(p1); 
    if( xmit_stat == 0 )
      printf("\n File Successfully Transmitted.\n");
    else {
      printf("\n Error occured during transmission %s\n",file_ptr);
      printf("\n Returned status = %i\n",xmit_stat);
    }
  }

}

/******************************************************
*                                                     *
*                                                     *
*              MAIN REMAPPING DRIVER                  *
*                                                     *
*                                                     *
*******************************************************/

main(argc,argv)
int argc;
char *argv[];
{ 

/* Variables used only in remap_process */

   int i_last_scan,i_first_scan;
   int i_tilt_proc;
   int i4time_vol,i_num_finished_products,i_status;

/* Variables used for data access and in fill_common */
   
   float b_ref[N_REF_TILT], *ref_ptr, *ref0_ptr;
   float b_vel[N_VEL_TILT], *vel_ptr, *vel0_ptr;

   float v_nyquist_ray_a[N_RAY_TILT], *nyq0_ptr;
   float azim[N_RAY_TILT], *azi0_ptr;
   float eleva;
   int ref_index, vel_index, io_stat;
   int n_rays, i_scan, i_tilt, n_ref_gates, n_vel_gates;
   float b_missing_data;

/* File transfer variables */

   char *compr_name;
   char *xmit_method, *remote_user, *remote_host, *dest_dir;  

/* Radar Location variables */

   int i_lat,i_lon,i_alt;
   float radar_lat, *rlat_ptr;
   float radar_lon, *rlon_ptr;
   float radar_alt, *ralt_ptr;
   char *rname_ptr;

/* Misc Local variables */

   char *sw;
   int iyr, iday, imon, ihour, imin, isec;    /* time variables */ 
   char source[80], *src_ptr;
   char string_time[9], *strtm_ptr;
   char full_fname[91], *fnm_ptr;
   int initial_ray;                     /* flag for first ray in volume */
   int alls_well, read_next, knt_bad_stat;
   int i_angle, past_angle;
   int past_scan, past_tilt;
   int len_fname;
   int compr_on, xmit_on, write_and_exit;
   int i_mode, i;
   int i_vcp;

#ifdef VERBOSE
   int ng_ref,ng_vel;          /* number of gates vel, ref returned */
   int gsp_ref,gsp_vel;          /* gate spacing, ref, vel */
#endif

/* Beginning of Executable Code */


/* Some initializations */

   n_vel_gates = NUM_VEL_GATES;
   n_ref_gates = NUM_REF_GATES;
   b_missing_data = MISSING;

   strtm_ptr = string_time;
   fnm_ptr = full_fname;
   ref0_ptr= b_ref;
   vel0_ptr= b_vel;
   azi0_ptr= azim;
   nyq0_ptr= v_nyquist_ray_a;

   rlat_ptr=&radar_lat;
   rlon_ptr=&radar_lon;
   ralt_ptr=&radar_alt;


/* Get Radar name environment variable */

   allocate(rname_ptr);
   rname_ptr=getenv("RADARNAME");

   if (rname_ptr == NULL )
   {
     printf (" Couldn't evaluate RADARNAME environment variable.\n");
     printf (" Set the 4-character radar name using:\n");
     printf (" setenv RADARNAME Kxxx\n");
     printf ("    before running remap.\n");
     exit(1);
   }

/* Get LAPS environment variables */

   allocate(compr_name);
   compr_name=getenv("LAPS_COMPRESS");
   allocate(xmit_method);
   xmit_method=getenv("LAPS_XMIT");
   allocate(remote_user);
   remote_user=getenv("LAPS_USER");
   allocate(remote_host);
   remote_host=getenv("LAPS_REMOTE");
   allocate(dest_dir);
   dest_dir=getenv("LAPS_DEST");
   if (dest_dir == NULL )
     dest_dir = ".";

   compr_on = TRUE;
   xmit_on = TRUE;
   if (remote_user == NULL || remote_host == NULL ) {
     xmit_on = FALSE;
     printf ("\n\n\n Trouble evaluating LAPS environment variables\n");
     printf ("    which specify file transfer options.\n");
     printf (" Will assume no file transferring desired.\n\n");
   }

   if (xmit_method == NULL )
     xmit_method = "rcp";

   if (compr_name == NULL )
     compr_on = FALSE;

   printf (" LAPS_COMPRESS: %s\n",compr_name);
   printf (" LAPS_XMIT:     %s\n",xmit_method);
   printf (" LAPS_USER:     %s\n",remote_user);
   printf (" LAPS_REMOTE:   %s\n",remote_host);
   printf (" LAPS_DEST:     %s\n",dest_dir);

/* Establish whether the user wants to read from a tape drive */

#ifdef VERBOSE
   printf (" Number of arguments is %i\n",argc);
#endif

   i_mode=1;
   src_ptr=source;

   if (argc < 3) {
     src_ptr = RT_BUFFER;
     if ( i_mode == REALTIME )
     {
       printf (" Reading real-time datastream from radar %s\n",rname_ptr);
       printf (" Using circular buffer %s\n",source);
     }
     else
     {
       printf (" Remapper compiled with tape playback libraries\n");
       printf (" Provide tape drive name, or use real-time version.\n");
       printf (" Usage: remaprt.exe    OR\n");
       printf (" Usage: remapa2.exe -f tape_device\n");
     }
   } else {
     if ( i_mode == PLAYBACK )
     {
       sw = argv[1];
       src_ptr = argv[2];
       printf(" The command line switch was %s and drive name %s\n",
            sw,source);
       printf (" Archive data is from radar %s\n",rname_ptr);
       printf (" Reading from tape drive %s\n",source);
     }
     else
     {
       printf (" Remapper compiled with real-time libraries\n");
       printf (" Remove tape drive name, or use real-time version.\n");
       printf (" Usage: remaprt.exe    OR\n");
       printf (" Usage: remapa2.exe -f [host:]tape_device\n");
     }
   }

/* call Archive II initialization routine */

   radar_init(src_ptr);

   i_alt=389.;
   i_alt=get_altitude();
   i_lat=get_latitude();
   i_lon=get_longitude();

   radar_alt= (float) i_alt;
   radar_lat=  0.00001 * (float) i_lat;
   radar_lon= -0.00001 * (float) i_lon;
   printf ("\n Radar altitude (m): %f\n",radar_alt);
   printf (" Radar latitude (degrees): %f\n",radar_lat);
   printf (" Radar longitude (degrees): %f\n\n",radar_lon);

/* call lut_gen FORTRAN routine */
/*
   lut_gen(rname_ptr,rlat_ptr,rlon_ptr,ralt_ptr);
*/

/* get data indices needed for other a2io library calls */

   ref_index=get_field_num("DBZ");
   printf (" Retrieved reflectivity index as %i\n",ref_index);
   vel_index=get_field_num("VEL");
   printf (" Retrieved velocity index as %i\n",vel_index);

/* Misc initializations */

   initial_ray  = TRUE;
   n_rays=0;
   i_first_scan = TRUE;
   i_tilt_proc = 0;
   i_last_scan = FALSE;
   write_and_exit = FALSE;
   read_next = TRUE;
   alls_well = TRUE;

/* Begin infinite loop to continuously read radar data */

   while (alls_well) {

/* Begin loop to fill buffer arrays with data from the circular buffer.
   Call remap routines and reset pointers at the end of a volume scan */

     if (read_next) {

#      ifdef VERBOSE
         printf( " Calling read_radial \n" );
#      endif

       io_stat=read_radial();
       if(io_stat == 1) 
       {
         printf( " Read_radial returned double eof \n");
         write_and_exit = TRUE;
       }

#      ifdef VERBOSE
         printf( " Back from read_radial \n");
#      endif

     }
     else
       read_next = TRUE;

/* Test for existence of velocity data.
   Do we also need to test for reflectivity data?  */

       if ( get_status(ref_index) == GOOD_STATUS ||
            get_status(vel_index) == GOOD_STATUS ) {
         knt_bad_stat = 0;
         i_angle = get_fixed_angle();
         i_scan = get_scan();
         i_tilt = get_tilt();

#        ifdef VERBOSE
           printf( " Good status received\n");
           printf( " i_angle = %i  i_tilt = %i \n", i_angle, i_tilt);
#        endif

         if ( initial_ray == TRUE ) {
           ref_ptr = b_ref;
           vel_ptr = b_vel;
           past_scan = i_scan;
           past_tilt = i_tilt;
           past_angle = i_angle;
           eleva = 0.01 * (float) i_angle;

           iyr = get_year();
           imon = get_month();
           iday = get_day();
           ihour = get_hour();
           imin = get_min();
           isec = get_sec();
           i4time_vol = int_to_i4time(&iyr,&imon,&iday,&ihour,&imin,&isec);
           make_fnam_lp (&i4time_vol,strtm_ptr,&i_status);

           i_vcp=get_vcp();
           printf ("  VCP number for this volume: %i\n",i_vcp);

#          ifdef VERBOSE
             printf ("  iyr = %i imon = %i iday = %i\n",iyr,imon,iday);
             printf ("i4time_vol returned %i\n",i4time_vol);
#          endif

           printf ("ihour = %i imin = %i isec = %i\n",ihour,imin,isec);
           printf ("Time is %s\n",string_time);
           initial_ray = FALSE;
         }

         if( i_tilt == past_tilt && i_scan == past_scan &&
             n_rays < N_RAY_TILT ) {

           n_rays ++;
           azim[n_rays-1] = 0.01 * (float) get_azi();
           v_nyquist_ray_a[n_rays-1] = 0.01 * (float) get_nyquist();

#          ifdef VERBOSE
             printf( "   INFO FOR n_rays = %i \n",n_rays);
             printf( "   ref_ptr = %i    vel_ptr = %i\n",
                         ref_ptr,vel_ptr);
#          endif
           if ( (n_rays-1) % 60 == 0)
             printf( " eleva = %f  azim = %f  Nyqst = %f\n",
                  eleva,azim[n_rays-1],v_nyquist_ray_a[n_rays-1]);
#          ifdef VERBOSE
             ng_ref = get_number_of_gates(ref_index);
             gsp_ref = get_first_gate(ref_index);
             ng_vel = get_number_of_gates(vel_index);
             gsp_vel = get_first_gate(vel_index);
             printf( " ref: Number of gates = %i,  first gate = %i\n",
                      ng_ref,gsp_ref);
             printf( " vel: Number of gates = %i,  first gate = %i\n",
                      ng_vel,gsp_vel);
#          endif

           io_stat = get_data_field(ref_index, ref_ptr, n_ref_gates);
           io_stat = get_data_field(vel_index, vel_ptr, n_vel_gates);

#          ifdef VERBOSE
             printf( " sample reflectivities %f  %f  %f  %f\n",
                     *(ref_ptr+20),*(ref_ptr+40),
                     *(ref_ptr+60),*(ref_ptr+80));
             printf( " sample velocities     %f  %f  %f  %f\n",
                     *(vel_ptr+20),*(vel_ptr+40),
                     *(vel_ptr+60),*(vel_ptr+80));
#          endif

           ref_ptr += NUM_REF_GATES;
           vel_ptr += NUM_VEL_GATES;

         } else {

           if( i_angle < past_angle || i_scan != past_scan )
             i_last_scan = TRUE;

/* call the FORTRAN routine to fill up the common data area  */

           printf( " Calling fill_common, i_angle = %i, past_angle = %i\n",
                     i_angle,past_angle); 
           printf( " n_rays = %i, past_tilt = %i, b_missing_data = %f\n",
                     n_rays, past_tilt, b_missing_data);

#          ifdef VERBOSE
             printf( " n_ref_gates = %i, n_vel_gates = %i\n",
                       n_ref_gates, n_vel_gates);
#          endif

/*
           fill_common(
               ref0_ptr,vel0_ptr,&n_rays,&i_tilt,
               &n_ref_gates,&n_vel_gates,
               azi0_ptr,nyq0_ptr,&eleva,&b_missing_data);
*/

/* call the FORTRAN remapper module  */

           i_tilt_proc++;

           printf(" Calling remap_process past_tilt = %i\n", past_tilt);
           printf(" Calling remap_process i_tilt_proc = %i\n", i_tilt_proc);

           printf(" i_last = %i, i_first = %i\n",
                i_last_scan,i_first_scan);
           printf(" i4time_vol = %i, i_num = %i,  i_status= %i\n",
                i4time_vol,i_num_finished_products,i_status);

/*
           remap_process(
                &i_tilt_proc,&i_last_scan,&i_first_scan,
                &i4time_vol,fnm_ptr,&len_fname,
                &i_num_finished_products,&i_status);
*/

           if (i_num_finished_products > 0 ) {
              *(fnm_ptr+len_fname)='\0';
             if(compr_on == TRUE)
               compr_file_sub(fnm_ptr,compr_name,i_status);
             if(xmit_on == TRUE)
               send_file_sub(fnm_ptr,xmit_method,
                     remote_user,remote_host,dest_dir,i_status);
           }

           i_last_scan = FALSE;
           i_first_scan = FALSE;

           if( i_angle < past_angle || i_scan != past_scan ) {
             i_first_scan = TRUE;
             i_tilt_proc = 0; 
             past_angle= i_angle;
           }
           n_rays = 0;
           initial_ray = TRUE;
           read_next = FALSE;
         }        /* end of tilt-down or n_ray overflow block */

/* For bad status, increment bad status counter and try again. */

       } else if ( knt_bad_stat < MAX_BAD_STAT && 
                   write_and_exit == FALSE ) {

#        ifdef VERBOSE
           printf( " Bad status received for data\n");
#        endif

         knt_bad_stat ++;

/* Once 1000 consecutive bad stati have been received, assume end
   of data and dump what might be in the buffer. */

       } else {
         printf ( " %i bad read status reports received \n",
                   knt_bad_stat);
         knt_bad_stat = 0;
         if ( n_rays > 0) {

           printf ( " Transferring %i available radials\n", n_rays);

           if( i_angle < past_angle )
             i_last_scan = TRUE;

/* call the FORTRAN routine to fill up the common data area  */

/*
           fill_common(
             ref0_ptr,vel0_ptr,&n_rays,&past_tilt,
             &n_ref_gates,&n_vel_gates,
             &azi0_ptr,&nyq0_ptr,&eleva,&b_missing_data);
*/

/* call the FORTRAN remapper module  */

           i_tilt_proc++;

           printf(" Calling remap_process past_tilt = %i\n", past_tilt);
           printf(" Calling remap_process i_tilt_proc = %i\n", i_tilt_proc);
           printf(" i_last = %i, i_first = %i\n",
              i_last_scan,i_first_scan);
           printf(" i4time_vol = %i, i_num = %i,  i_status= %i\n",
              i4time_vol,i_num_finished_products,i_status);

/*
           remap_process(
              &i_tilt_proc,&i_last_scan,&i_first_scan,
              &i4time_vol,fnm_ptr,&len_fname,
              &i_num_finished_products,&i_status);
*/

             if (i_num_finished_products > 0 ) {
                *(fnm_ptr+len_fname)='\0';
                if(compr_on == TRUE)
                   compr_file_sub(fnm_ptr,compr_name,i_status);
                if(xmit_on == TRUE)
                   send_file_sub(fnm_ptr,xmit_method,
                     remote_user,remote_host,dest_dir,i_status);
             }

           if(write_and_exit == TRUE) exit(0);

           i_last_scan = FALSE;
           i_first_scan = FALSE;

           if( i_angle < past_angle ) {
             i_first_scan = TRUE;
             i_tilt_proc = 0;
             past_angle = i_angle;
           }
           n_rays = 0;
           initial_ray = TRUE;

         } /* close n_rays > 0 block */
      }    /* close velocity status block */ 
   }       /* close infinite while loop */
}          /* end of main */
