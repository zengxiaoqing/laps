
      program genlvdlut

      implicit none

      integer  nx_l,ny_l
      integer  istatus

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
c
c ======================== START ==============================
c
      call get_grid_dim_xy(nx_l,ny_l,istatus)
      if(istatus.ne.1)then
         write(6,*)'Error getting nx_l/ny_l'
         goto 1000
      else
         write(6,*)'LAPS nx_l and ny_l obtained'
      endif
c
c---------------------------------------------------------------
c
      call config_satellite_lvd(istatus)
      if(istatus.ne.1)then
         write(*,*)'Error configuring satellite-master common'
         goto 1000
      endif
c
c------------------- MAIN PROGRAM -----------------------
c
      call genlvdlut_sub(nx_l,ny_l,istatus)
      if(istatus.lt.0)then
         write(*,*)'Error generating lut'
         goto 1000
      endif

c ---------------------------------------------------------
      call rewrite_satellite_lvd_nl(istatus)
      if(istatus.ne.1)then
         write(*,*)'Error returned from rewrite_satellite_master_nl'
         write(*,*)'Check satellite_master.nl namelist in data/static'
         goto 1000
      endif

      print*, 'Finished in genlvdlut'
1000  stop
      end
