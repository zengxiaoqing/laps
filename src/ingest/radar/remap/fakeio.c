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
#include <strings.h>
#include <time.h>
struct tm rtime, *radtime;
struct count_struct {
     int iazi;
     int iangle;
     int itilt;
     int iscan;
} radcount;

struct loc_struct {
     char rsite[80];
     char radid[4];
     int latitude;
     int longitude;
     int altitude;
} radinfo;

void set_radar_name(site)
char site[80];
{
      extern struct loc_struct radinfo;
      strcpy(radinfo.rsite,site);
}

int get_field_num(fchar)
char fchar[3];
{
      if(strcmp(fchar,"DBZ") == 0)
        return(0);
      else if(strcmp(fchar,"VEL") == 0)
        return(1);
      else if(strcmp(fchar,"SPW") == 0)
        return(2);
      else
        return(3);
}

void radar_init(drive)
char drive[80];
{
      extern struct count_struct radcount;
      extern struct tm rtime, *radtime;
      printf("/n");
      printf(" Initializing a2io_fake for file testing only/n");
      printf("/n");
      radtime =&rtime;
      radcount.iazi=0;
      radcount.iangle=100;
      radcount.itilt=1;
      radcount.iscan=1;
}

#define MAXTILT 14
int kelev[MAXTILT] =
         {50,150,240,340,430,530,620,750,
          870,1000,1200,1400,1670,1950};

int read_radial()
{
      extern struct tm rtime, *radtime;
      extern struct count_struct radcount;
      extern int kelev[];
      int istat;
      long curtime[2], *tp;
      long timez[2], *tzp;

      tp = curtime;
      tzp = timez;

      radcount.iazi=radcount.iazi+100;
      if ( radcount.iazi > 35999 ) {
        radcount.iazi=radcount.iazi-36000;
        radcount.itilt++;
        if ( radcount.itilt > MAXTILT) {
          radcount.itilt=1;
          radcount.iscan++;
        }
      }
      radcount.iangle=kelev[radcount.itilt];
      istat=gettimeofday(tp,tzp);
      radtime=gmtime(tp);
      return(0);
}

int get_data_field(ifield,field,nsize)
int ifield;
float field[];
int nsize;
{
    int i;
    float *ptr;
    ptr=field;
 
    if( ifield == 0 )
      for (i = 0; i < nsize; i++){
        *ptr = (0.1 * (float) i);
        ptr++;
      }
    else
      for (i = 0; i < nsize; i++){
        *ptr = (0.025 * (float) i);
        ptr++;
      }
    return(0);
}

int get_year()
{
      extern struct tm *radtime;
      return(radtime->tm_year);
}

int get_month()
{
      extern struct tm *radtime;
      int imon;
      imon=radtime->tm_mon + 1;
      return(imon);
}

int get_day()
{
      extern struct tm *radtime;
      return(radtime->tm_mday);
}

int get_hour()
{
      extern struct tm *radtime;
      return(radtime->tm_hour);
}

int get_min()
{
      extern struct tm *radtime;
      return(radtime->tm_min);
}

int get_sec()
{
      extern struct tm *radtime;
      return(radtime->tm_sec);
}

int get_azi()
{
      extern struct count_struct radcount;
      return(radcount.iazi);
}

int get_fixed_angle()
{
      extern struct count_struct radcount;
      return(radcount.iangle);
}

int get_scan()
{
      extern struct count_struct radcount;
      return(radcount.iscan);
}

int get_tilt()
{
      extern struct count_struct radcount;
      return(radcount.itilt);
}

int get_status(index)
int index;
{
      return(1);
}

int get_rt_mode()
{
      return(1);
}
 
int get_nyquist()
{
      return(3500);
}

int get_vcp()
{
      return(41);
}
 
int get_latitude()
{
      return(3978666);
}
 
int get_longitude()
{
      return(10454555);
}

int get_altitude()
{
      return(1300);
}
