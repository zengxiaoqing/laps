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
#define allocate(s) (s= (char *) malloc(80))

main()
{
/* File transfer variables */

   char *fnm_ptr;
   char *compr_name, *remote_user, *remote_host, *dest_dir;
   int i_status;

/* Beginning of Executable Code */

   allocate(compr_name);
   compr_name=getenv("LAPS_COMPRESS");
   allocate(remote_user);
   remote_user=getenv("LAPS_USER");
   allocate(remote_host);
   remote_host=getenv("LAPS_REMOTE");
   allocate(dest_dir);
   dest_dir=getenv("LAPS_DEST");

   allocate(fnm_ptr);
   fnm_ptr="/tmp_mnt/nsslsun/nsslsun/caps/remap/test.file";

   send_file_sub(fnm_ptr,compr_name,
                 remote_user,remote_host,dest_dir,i_status);

}

send_file_sub(file_ptr,compr_ptr,
              user_ptr,host_ptr,dest_ptr,i_status)
char  *compr_ptr,*file_ptr,*user_ptr,*host_ptr,*dest_ptr;
int i_status;
{
  char *string1;
  FILE *comp_file; 
  FILE *p1;
  int loc;
  int slash;
  char c_slash;
  char buf[200];

  allocate(string1);

  if(strcmp(compr_ptr,"\0") == 0)
    fprintf(stdout,"%s","\n No File Compression\n");
  else if(strcmp(compr_ptr,"gzip") == 0)
  {
    strcpy(string1,"gzip ");
    strcat(string1,file_ptr); 
    comp_file= popen(string1,"r");
    strcat(file_ptr,".gz"); 
    if(comp_file == (FILE *) NULL)
      fprintf(stdout,"%s","\n Error occured while Compressing..\n");
    else
      fprintf(stdout,"%s","\n File Successfully Compressed\n");

    fclose(comp_file);
  }
  else if(strcmp(compr_ptr,"compress") == 0)
  {
    strcpy(string1,"compress ");
    strcat(string1,file_ptr); 
    comp_file= popen(string1,"r");
    strcat(file_ptr,".Z"); 
    if(comp_file == (FILE *) NULL)
      fprintf(stdout,"%s","\n Error occured while Compressing..\n");
    else
      fprintf(stdout,"%s","\n File Successfully Compressed\n");



    fclose(comp_file);
  }
  else
  {
    fprintf(stdout,"%s %s",compr_ptr," Compression Unsupported\n");
    fprintf(stdout,"%s"," Set LAPS_COMPRESS to gzip or compress\n");
  }

  fprintf(stdout,"%s","\n Connecting to remote host.....\n");

/* Strip any directory info from the file name so all
*  directory info in the destination file name is from
*  dest_ptr, environment variable LAPS_DEST.
*  The source file name is unchanged */

/*  c_slash='/'; slash=atoi("/\0");  */
  slash=47;
  string1=(char *)strrchr(file_ptr,slash);

/* Create the rcp command */

  if (string1 == NULL)
  sprintf(buf,"rcp %s %s@%s:%s/%s",file_ptr,user_ptr,host_ptr,
             dest_ptr,file_ptr);
  else
  sprintf(buf,"rcp %s %s@%s:%s%s",file_ptr,user_ptr,host_ptr,
             dest_ptr,string1);

  fprintf(stdout,buf,"%s","\n");

  p1= popen(buf,"r"); 
  if (p1 == (FILE *) NULL)
  {
    fclose(p1);
    fprintf(stderr,"%s","Cannot RCP to host\n");
  }

 fprintf(stdout,"%s","\n End of file transfer subroutine\n");

 fclose(p1);
}
