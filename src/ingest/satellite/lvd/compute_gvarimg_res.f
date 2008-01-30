
      subroutine compute_gvarimage_resolution(rp_div,rl_div,
     &rpix,rline,start_pix,start_line,r_img_res_m,istatus)
c
c
      implicit none

      Real*8        rlat8_1,rlat8_2,rlat8_3
      Real*8        rlon8_1,rlon8_2,rlon8_3
      Real*8        elev1,elev2,scan1,scan2
      Real*8        scpx,evln
      Real*8        RL
      Real*8        RP

      Real        rlat1,rlat2,rlat3,rlon1,rlon2,rlon3
      Real        grid_spacing_scan,grid_spacing_elev
      Real        rpix,rline
      Real*8        radtodeg
      Real*8        pi
      Real        r_img_res_m
      Real*8        rl_div
      Real*8        rp_div

      Integer     instr
      Integer     istatus
      Integer     ierr
      Integer     start_line
      Integer     start_pix
      real        cosd
c --------------------------------------------
c start
c
      istatus=0

      pi=acos(-1.)
      radtodeg=180.d0/pi
      instr = 1  !imager=1, sounder = 2.
c
c compute image resolution in meters
c
      RP=float(int( ( (rpix -1.) *rp_div ) +start_pix) )
      RL=float(int( ( (rline-1.) *rl_div ) +start_line))
      
      elev1 = evln(instr,RL)
      scan1 = scpx(instr,RP)
      call lpoint(elev1,scan1,rlat8_1,rlon8_1,ierr)

      elev2 = evln(instr,RL+rl_div)
      scan2 = scpx(instr,RP+rp_div)
      call lpoint(elev1,scan2,rlat8_2,rlon8_2,ierr)
      call lpoint(elev2,scan1,rlat8_3,rlon8_3,ierr)

      rlat1=rlat8_1*radtodeg
      rlat2=rlat8_2*radtodeg
      rlon1=rlon8_1*radtodeg
      rlon2=rlon8_2*radtodeg
      rlat3=rlat8_3*radtodeg
      rlon3=rlon8_3*radtodeg

      grid_spacing_scan = sqrt( (rlat1-rlat2)**2 +
     1             ((rlon1-rlon2)*cosd(rlat1))**2)
      grid_spacing_scan = grid_spacing_scan*111.1*1000.0   !result in meters

      grid_spacing_elev=sqrt( (rlat1-rlat3)**2 +
     1             ((rlon1-rlon3)*cosd(rlat1))**2)
      grid_spacing_elev=grid_spacing_elev*111.1*1000.0   !result in meters

      r_img_res_m=(grid_spacing_scan+grid_spacing_elev)/2.0

      write(6,*)'Image Resolution (m): ',r_img_res_m

      if(ierr.ne.0)istatus=-1
      return
      end
