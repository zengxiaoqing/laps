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
/*  remap.c

    Reads from NSSL Archive II routines and passes the data
    to FSL LAPS remapping routines.

    Usage:
      remap                       reads data from live data
      remap -f /dev/rmt0          reads data from tape drive rmt0

    Advertised to work for live WSR-88D radar data as well
    as archive II taped data.

    Requires NSSL libraries:  -la2rt  for real-time
                              -la2tp  and -ltpio for tape playback

    ORIGINAL: 
       Steve Albers, FSL and Keith Brewster, CAPS       May 1994
    MODIFICATIONS:
       Debug and overhaul  Keith Brewster, CAPS         June 1994
    ADDITION
       Network Program  Abraham Kappamammmootil, CAPS   March 1995 
    
*/
#include <stdio.h>

#include <stdlib.h>


/* When VERBOSE is defined, the program generates verbose output
*  that is useful for testing but annoying for operations */

/* #define VERBOSE */

/* When HOST1 is defined, the program will try to send the remapped
*  files to DEST1 on machine HOST1 */
#define HOST1 "stratus.gcn.uoknor.edu"
#define DEST1 "/scratch/stratus/akappam"
#define allocate(s) (s= (char *) malloc(80))

/* #define RADAR_NAME "KOUN" */
#define RADAR_NAME "KTLX"
#define RT_BUFFER "R"

#define NUM_REF_GATES 460
#define NUM_VEL_GATES 920
#define N_REF_TILT 174800
#define N_VEL_TILT 349600
#define N_RAY_TILT 380
#define NSIZE 920

#define REF_INDEX 0
#define VEL_INDEX 1
#define SPW_INDEX 2
#define SNR_INDEX 3

#define GOOD_STATUS 1
#define BAD_STATUS 0
#define MAX_BAD_STAT 1000

#define TRUE 1
#define FALSE 0
#define MISSING -999.0

#define lut_gen lut_gen_
#define fill_common fill_common_
#define remap_process remap_process_
#define int_to_i4time int_to_i4time_
#define make_fnam_lp make_fnam_lp_

/* The following are needed to use the fake a2io subroutines */
#define get_nyquist get_nyquist_
#define get_data get_data_
#define get_sec get_sec_
#define get_month get_month_
#define get_status get_status_
#define set_radar_name set_radar_name_
#define get_year get_year_
#define get_tilt get_tilt_
#define get_azi get_azi_
#define get_scan get_scan_
#define get_fixed_angle get_fixed_angle_
#define get_hour get_hour_
#define get_day get_day_
#define read_radial read_radial_
#define get_min  get_min_
#define a2_init a2_init_

send_file_sub(file_ptr,host_ptr,dest_ptr,i_status)
char  *file_ptr,*host_ptr,*dest_ptr;
int i_status;
{

  char *string1;  
  FILE *comp_file;  
  FILE *p1;
  char buffer[100];
  char buffer1[100];
  char buf[100];

   allocate(string1);

   strcpy(string1,"compress ");
   strcat(string1,file_ptr);  
   comp_file= popen(string1,"r"); 

  if(comp_file == (FILE *) NULL)
  {
   fprintf(stdout,"%s","\n Error occured while Compressing..\n");
   fprintf(stdout,"%s","\n Exiting......\n");
   exit(0);
  }
  else
   fprintf(stdout,"%s","\n File Successfully Compressed\n");

   strcat(file_ptr,".Z");  

   fclose(comp_file);


 start:  fprintf(stdout,"%s","\n Connecting to STRATUS.......\n");

 sprintf(buf,"rcp %s akappam@stratus:/scratch/stratus/akappam/%s",
                                                  file_ptr,file_ptr);

  p1= popen(buf,"r");  

   if (p1 == (FILE *) NULL)
  {
   fclose(p1);

    fprintf(stderr,"%s","Cannot RCP to the STRATUS\n");
    fprintf(stdout,"%s","\n Connecting to TORNADO.......\n");

    sprintf(buf,"rcp %s akappam@stratus:/scratch/stratus/akappam/%s",
                                                  file_ptr,file_ptr);

   p1= popen(buf,"r");  
  
   if (p1 == (FILE *) NULL)
   {
    fprintf(stderr,"%s","Cannot RCP to the STRATUS\n");
    fprintf(stdout,"%s","\n Connecting to CIRRUS.......\n");
    fclose(p1);

    sprintf(buf,"rcp %s akappam@stratus:/scratch/stratus/akappam/%s",
                                                  file_ptr,file_ptr);

    p1= popen(buf,"r");  
  
    if (p1 == (FILE *) NULL)
    {
     fprintf(stderr,"%s","Cannot RCP to the CIRRUS\n");
     fprintf(stdout,"%s","\n Connecting to CUMULUS.......\n");
     fclose(p1);

     sprintf(buf,"rcp %s akappam@stratus:/scratch/stratus/akappam/%s",
                                                  file_ptr,file_ptr);

     p1= popen(buf,"r");  
  
     if (p1 == (FILE *) NULL)
     {
      fprintf(stderr,"%s","Cannot RCP to the CUMULUS\n");
      goto start;   
     }
    }
   }
  }


 fprintf(stdout,"%s","\n File Successfully send to /scratch/stratus/akappam\n");


 fclose(p1); 
}



main(argc,argv)
int argc;
char *argv[];
{ 

/* Variables used only in remap_process */

   int i_last_scan,i_first_scan;
   int i_tilt_proc;
   int i4time_vol,i_num_finished_products,i_status;

/* Variables used in fill_common */
   
   float b_ref[N_REF_TILT], *ref_ptr, *ref0_ptr;
   float b_vel[N_VEL_TILT], *vel_ptr, *vel0_ptr;

   float v_nyquist_ray_a[N_RAY_TILT], *nyq0_ptr;
   float azim[N_RAY_TILT], *azi0_ptr;
   float eleva;
   int n_rays, i_scan, i_tilt, n_ref_gates, n_vel_gates;
   float b_missing_data;

/* Local variables */

   float dumarray[NSIZE], *dum_ptr;
   char *site = RADAR_NAME;
   char *sw;
   char *drive;
   int iyr, iday, imon, ihour, imin, isec;    /* time variables */ 
   char string_time[9], *strtm_ptr;
   char full_fname[91], *fnm_ptr;
   char *host, *dest_dir;
   int initial_ray;                     /* flag for first ray in volume */
   int alls_well, read_next, knt_bad_stat, i_angle, past_angle;
   int past_scan, past_tilt;
   int len_fname;
#ifdef VERBOSE
   int ng_ref,ng_vel;          /* number of gates vel, ref returned */
   int gsp_ref,gsp_vel;          /* gate spacing, ref, vel */
#endif

/* Beginning of Executable Code */
/* Set the radar name */
   set_radar_name(site);

/* Establish whether the user wants to read from a tapedrive */

#ifdef VERBOSE
   printf (" Number of arguments is %i\n",argc);
#endif

   if (argc < 3) {
     drive = RT_BUFFER;
     printf (" Reading real-time datastream from radar %s\n",site);
     printf (" Using circular buffer %s\n",drive);
   } else {
     sw = argv[1];
     drive = argv[2];
     printf(" The command line switch was %s and drive name %s\n",
            sw,drive);
     printf (" Reading from tape drive %s\n",drive);
   }

/* call lut_gen FORTRAN routine */
   lut_gen();
/* call Archive II initialization routine */
   a2_init(drive);

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
   dum_ptr = dumarray;

   initial_ray  = TRUE;
   n_rays=0;
   i_first_scan = TRUE;
   i_tilt_proc = 0;
   i_last_scan = FALSE;
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

       read_radial();

#      ifdef VERBOSE
         printf( " Back from read_radial \n");
#      endif

     }
     else
       read_next = TRUE;

/* Test for existence of velocity data.
   Do we also need to test for reflectivity data?  */

     if ( get_status(VEL_INDEX) == GOOD_STATUS) {
       knt_bad_stat = 0;
       i_angle = get_fixed_angle();
       i_scan = get_scan();
       i_tilt = get_tilt();

#      ifdef VERBOSE
         printf( " Good status received for velocity \n");
         printf( " i_angle = %i  i_tilt = %i \n", i_angle, i_tilt);
#      endif

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

#        ifdef VERBOSE
           printf ("  iyr = %i imon = %i iday = %i \n",iyr,imon,iday);
           printf ("i4time_vol returned %i \n",i4time_vol);
#        endif

         printf ("ihour = %i imin = %i isec = %i \n",ihour,imin,isec);
         printf ("Filename string %s \n",string_time);
         initial_ray = FALSE;
       }

       if( i_tilt == past_tilt && i_scan == past_scan ) {

         n_rays ++;
         if( n_rays < N_RAY_TILT) {
           azim[n_rays-1] = 0.01 * (float) get_azi();
           v_nyquist_ray_a[n_rays-1] = 0.01 * (float) get_nyquist();

#          ifdef VERBOSE
             printf( "   INFO FOR n_rays = %i \n",n_rays);
             printf( "   ref_ptr = %i    vel_ptr = %i \n",
                       ref_ptr,vel_ptr);
#          endif
           if ( (n_rays-1) % 60 == 0)
             printf( " eleva = %f  azim = %f  Nyqst = %f \n",
                  eleva,azim[n_rays-1],v_nyquist_ray_a[n_rays-1]);
#          ifdef VERBOSE
             ng_ref = get_number_of_gates(REF_INDEX);
             gsp_ref = get_first_gate(REF_INDEX);
             ng_vel = get_number_of_gates(VEL_INDEX);
             gsp_vel = get_first_gate(VEL_INDEX);
             printf( " ref: Number of gates = %i,  first gate = %i \n",
                    ng_ref,gsp_ref);
             printf( " vel: Number of gates = %i,  first gate = %i \n",
                    ng_vel,gsp_vel);
#          endif
           get_data(ref_ptr,vel_ptr,
                    dum_ptr,dum_ptr, NSIZE);

#          ifdef VERBOSE
             printf( " sample reflectivities %f  %f  %f  %f \n",
                   *(ref_ptr+20),*(ref_ptr+40),
                   *(ref_ptr+60),*(ref_ptr+80));
             printf( " sample velocities     %f  %f  %f  %f \n",
                   *(vel_ptr+20),*(vel_ptr+40),
                   *(vel_ptr+60),*(vel_ptr+80));
#          endif

           ref_ptr += NUM_REF_GATES;
           vel_ptr += NUM_VEL_GATES;
         }

       } else {

         if( i_angle < past_angle || i_scan > past_scan )
           i_last_scan = TRUE;

/* call the FORTRAN routine to fill up the common data area  */

           printf( " Calling fill_common, i_angle = %i, past_angle = %i \n",
                     i_angle,past_angle); 
           printf( " n_rays = %i, past_tilt = %i, b_missing_data = %f \n",
                     n_rays, past_tilt, b_missing_data);

#          ifdef VERBOSE
             printf( " n_ref_gates = %i, n_vel_gates = %i \n",
                       n_ref_gates, n_vel_gates);
#          endif

           fill_common(
             ref0_ptr,vel0_ptr,&n_rays,&i_tilt,
             &n_ref_gates,&n_vel_gates,
             azi0_ptr,nyq0_ptr,&eleva,&b_missing_data);

/* call the FORTRAN remapper module  */

           i_tilt_proc++;

           printf(" Calling remap_process past_tilt = %i \n", past_tilt);
           printf(" Calling remap_process i_tilt_proc = %i \n", i_tilt_proc);

           printf(" i_last = %i, i_first = %i, \n",
              i_last_scan,i_first_scan);
           printf(" i4time_vol = %i, i_num = %i,  i_status= %i \n",
              i4time_vol,i_num_finished_products,i_status);

           remap_process(
              &i_tilt_proc,&i_last_scan,&i_first_scan,
              &i4time_vol,fnm_ptr,&len_fname,
              &i_num_finished_products,&i_status);

#          ifdef HOST1
             if (i_num_finished_products > 0 ) {
                *(fnm_ptr+len_fname)='\0';
                host=HOST1;
                dest_dir=DEST1;
                send_file_sub(fnm_ptr,host,dest_dir,i_status);
             }
#          endif

           i_last_scan = FALSE;
           i_first_scan = FALSE;

           if( i_angle < past_angle ) {
             i_first_scan = TRUE;
             i_tilt_proc = 0; 
             past_angle= i_angle;
           }
           n_rays = 0;
           initial_ray = TRUE;
           read_next = FALSE;

         }     /* end of tilt down block */

/* For bad status, increment bad status counter and try again. */

       } else if (knt_bad_stat < MAX_BAD_STAT) {

#        ifdef VERBOSE
           printf( " Bad status received for velocity \n");
#        endif

         knt_bad_stat ++;

/* Once 1000 consecutive bad stati have been received, assume end
   of data and dump what might be in the buffer. */

       } else {
         printf ( " %i bad read status reports received \n",
                   knt_bad_stat);
         knt_bad_stat = 0;
         if ( n_rays > 0) {

           printf ( " Transferring %i available radials \n", n_rays);

           if( i_angle < past_angle )
             i_last_scan = TRUE;

/* call the FORTRAN routine to fill up the common data area  */

           fill_common(
             ref0_ptr,vel0_ptr,&n_rays,&past_tilt,
             &n_ref_gates,&n_vel_gates,
             &azi0_ptr,&nyq0_ptr,&eleva,&b_missing_data);

/* call the FORTRAN remapper module  */

           i_tilt_proc++;

           printf(" Calling remap_process past_tilt = %i \n", past_tilt);
           printf(" Calling remap_process i_tilt_proc = %i \n", i_tilt_proc);
           printf(" i_last = %i, i_first = %i, \n",
              i_last_scan,i_first_scan);
           printf(" i4time_vol = %i, i_num = %i,  i_status= %i \n",
              i4time_vol,i_num_finished_products,i_status);

           remap_process(
              &i_tilt_proc,&i_last_scan,&i_first_scan,
              &i4time_vol,fnm_ptr,&len_fname,
              &i_num_finished_products,&i_status);

#          ifdef HOST1
             if (i_num_finished_products > 0 ) {
                *(fnm_ptr+len_fname)='\0';
                host=HOST1;
                dest_dir=DEST1;
                send_file_sub(fnm_ptr,host,dest_dir,i_status);
             }
#          endif

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
