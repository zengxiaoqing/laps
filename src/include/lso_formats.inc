c
c Definitions for reading and writing the lso file 
c Changes:
c         J. Edwards 09-16-98  Moved from lib/read_surface_obs.f and
c                              ingest/sao/write_surface_obs.f
c         P. Stamus  12-15-98  Changed vis format (in 911) to 8.2 from 8.0
c
 900	  format(1x,a24,2x,i6,2x,i6)
 901	  format(1x,a20,2x,i8,2x,a11,2x,f8.4,1x,f9.4,2x,f8.1,2x,i8)
 903	  format(9x,a6,6x,a6,2x,a30)
c903	  format(5x,a10,2x,a10,2x,a30) !variables are char a6 in code
 905	  format(5x,f8.2,2x,f5.2,2x,f8.2,2x,f5.2,2x,f8.2,2x,f5.2)
 907	  format(5x,4(f8.1,2x),2(f6.2,2x))
 909	  format(5x,f8.2,2x,f8.2,2x,f8.2,2x,i4,2x,f8.2,2x,f8.2,2x,f8.2)
 911	  format(5x,f8.2,2x,f6.3,2x,f8.1,2x,f6.1,2x,4f8.2)
 913	  format(5x,4(f8.2,2x),f6.1,2x,f8.2,2x,f4.1)
 915	  format(5x,i4,2x,f8.2,2x,f8.2)
 917	  format(10x,a8,2x,f8.0)
 919      format(a132)
