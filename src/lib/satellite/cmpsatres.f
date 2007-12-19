
      subroutine compute_sat_res_m(rp_div,rl_div,
     &rpix,rline,start_pix,start_line,instr,r_res_m,
     &istatus)
c
c
      include 'trigd.inc'
      implicit none

      Real*8      rlat8_1,rlat8_2,rlat8_3
      Real*8      rlon8_1,rlon8_2,rlon8_3
      Real*8      elev1,elev2,scan1,scan2
      Real*8      scpx,evln
      Real*8      RL
      Real*8      RP

      Real        rlat1,rlat2,rlat3,rlon1,rlon2,rlon3
      Real        grid_spacing_scan,grid_spacing_elev
      Real        rpix,rline
      Real*8      radtodeg
      Real*8      pi
      Real        r_res_m
c     Real        resscln
c     Real        reselev
      Real*8      rl_div
      Real*8      rp_div

      Integer     instr
      Integer     istatus
      Integer     ierr
      Integer     start_line
      Integer     start_pix
c --------------------------------------------
c start
c
c   x lat3,lon3
c
c   x lat1,lon1      x  lat2,lon2
c
      istatus=0
      pi=acos(-1.)
      radtodeg=180.d0/pi
c     instr = 1  !imager=1, sounder = 2.
c pix/line absolute coordinates
c
      RP= rpix*rp_div +float(start_pix)-rp_div
      RL= rline*rl_div+float(start_line)-rl_div
      elev1 = evln(instr,RL)
      scan1 = scpx(instr,RP)
      call lpoint(elev1,scan1,rlat8_1,rlon8_1,ierr)
c
      scan2 = scpx(instr,RP+rp_div)  !e-w increment
      elev2 = evln(instr,RL+rl_div)  !n-s increment
      call lpoint(elev1,scan2,rlat8_2,rlon8_2,ierr)
      call lpoint(elev2,scan1,rlat8_3,rlon8_3,ierr)

      rlat1=rlat8_1*radtodeg
      rlat2=rlat8_2*radtodeg
      rlon1=rlon8_1*radtodeg
      rlon2=rlon8_2*radtodeg
      rlat3=rlat8_3*radtodeg
      rlon3=rlon8_3*radtodeg

c     r_res_m=sqrt((rlat3-rlat1)**2+
c    1            ((rlon2-rlon1)*cosd(rlat1))**2)
c     r_res_m=r_res_m*111.1*1000.

c     resscln=abs(rlon2-rlon1)*111.1*1000.
c     reselev=abs(rlat1-rlat3)*111.1*1000.
c     r_res_m=(resscln+reselev)/2.0                        !result in meters

      grid_spacing_scan = sqrt( (rlat1-rlat2)**2 +
     1             ((rlon1-rlon2)*cosd(rlat1))**2)
      grid_spacing_scan = grid_spacing_scan*111.1*1000.0   !result in meters

      grid_spacing_elev=sqrt( (rlat1-rlat3)**2 +
     1             ((rlon1-rlon3)*cosd(rlat1))**2)
      grid_spacing_elev=grid_spacing_elev*111.1*1000.0   !result in meters

      r_res_m=(grid_spacing_scan+grid_spacing_elev)/2.0

      write(6,*)'Image Resolution (m): ',r_res_m

      if(ierr.ne.0)istatus=-1
      return
      end
