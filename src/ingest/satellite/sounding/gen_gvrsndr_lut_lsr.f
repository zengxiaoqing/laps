       Subroutine gen_gvrsndr_lut_lsr(c_filename_sat,nlines,nelems,
     &wavelength,resolution_x,resolution_y,rsndr_res_m,start_pix,
     &start_line,end_pix,end_line,ewCycles,ewIncs,nsCycles,nsIncs,
     &f_time,orbAt,nch,nxl,nyl,lat,lon,rpix,rline,istatus)
c
c
c
      implicit none
 
      Include       'instco.inc'

      Integer     nlines,nelems
      Integer     nxl,nyl
      Integer     nch
c
      real*8        wavelength(nch)
      real*4        scalingBias(nlines,nch)
      real*4        scalingGain(nlines,nch)
      real*4        lat(nxl,nyl)
      real*4        lon(nxl,nyl)
      real*4        rline(nxl,nyl)
      real*4        rpix(nxl,nyl)
      real*4        pi
      real*8        rl_div
      real*8        rp_div
      real*4        radtodeg
      real*4        rsndr_res_m
      real*4        r_missing_data
      real*4        rnx,rny,rnp
      real*4        rpct_out

      real*8        r8lat,r8lon
      real*8        ELEV,SCAN
      real*8        RL
      real*8        RP
      Real*8        orbAt(336)
      Real*4        time_50,time50
      Real*8        t50_8
      Real*8        t
      Real*8        f_time
      Real*8        SatSubLAT,SatSubLON
c
      Integer     ewCycles,ewIncs
      Integer     nsCycles,nsIncs
      Integer     start_line
      Integer     start_pix
      Integer     end_line
      Integer     end_pix
      Integer     resolution_x
      Integer     resolution_y


      Integer     i1,j1
      Integer     ndsize_ch
      Integer     ndsize_x(nlines)
      Integer     ndisze_y
      Integer     cstatus
      Integer     istatus
      Integer     mstatus
      Integer     nstatus
      Integer     wstatus
      Integer     IERR
      Integer     INSTR
      Integer     i,j,k,n
      Integer     n2,nn
      Integer     time_spec(2)
      Integer     imc
      Integer     npoints_out
      Integer     i4time_nearest_sat
      Integer     isat
      Integer     nrl,nrp

      Character     filename_parm*255
      Character     c_filename_sat*255
      Character*255 c_sounding_path
      Character     table_path*200
      Character*9   c_ftimesat
      Character     ctype*4
      Character     c_imc*4
c ---------------------------------------------------------
c read static information for gvarimage navigation
c
      istatus = -1 
c
c need a loop for this from i=1,nch
c
c        ctype=c_sndrch_types(i,isat)
c        filename_parm=cname(1:n2)//ctype//'.parms'
c        call read_sndr_parmfile(filename_parm,
c    &                         sBsG_OA_read,
c    &                         ewCycles,
c    &                         ewIncs,
c    &                         nsCycles,
c    &                         nsIncs,
c    &                         f_time,
c    &                         imc,
c    &                         start_line,
c    &                         start_pix,
c    &                         nlines,
c    &                         scalingBias,
c    &                         scalingGain,
c    &                         orbAt,
c    &                         SatSubLAT,
c    &                         SatSubLON,
c    &                         wavelength(i),
c    &                         rsndr_res_m,
c    &                         nstatus)

c        if(nstatus.ne.0)then
c
c Find files for GOES-8 sounding data.
c
c     call Read_sounder_db_cdf_header(c_filename_sat,nlines,
c    &                               nelems,ndsize_x,nch,
c    &                               wavelength,
c    &                               resolution_x,
c    &                               resolution_y,
c    &                               start_pix,start_line,
c    &                               end_pix,end_line,
c    &                               ewCycles,
c    &                               ewIncs,
c    &                               nsCycles,
c    &                               nsIncs,
c    &                               f_time,
c    &                               c_imc,
c    &                               orbAt,
c    &                               scalingBias,
c    &                               scalingGain,
c    &                               ndsize_ch,
c    &                               istatus)

c     if(istatus .ne.0)then
c        write(6,*)'No sounder netcdf file'
c        goto 901
c     endif

      call get_r_missing_data(r_missing_data,istatus)

      rsndr_res_m=(resolution_y+resolution_x)/2.0
      rl_div = resolution_y/8.
      rp_div = resolution_x/8.
      if(resolution_x.eq.0)then
         rp_div=1.
      endif
      if(resolution_y.eq.0)then
         rl_div=1.
      endif

c     read(c_imc(2:4),111)imc
c111   format(i3)
c     if(imc.ne.0)imc=1

      imc=1
      INSTR=2          !1=Imager, 2=Sounder
      pi=3.141592653589793
      radtodeg=180.0/pi

      call bcd_to_int(orbAt(12),time_spec)
      time_50 = time50(time_spec)
      t50_8=time_50
      t = f_time /60. + 7305. * 24. * 60.

      call SETCON(INSTR,nsCycles,nsIncs,ewCycles,ewIncs)
      call LMODEL(t,t50_8,OrbAt,imc,SatSubLAT,SatSubLON)

      write(6,*)'Sat Subpoint lat (deg) ',SatSubLAT*radtodeg
      write(6,*)'Sat Subpoint lon (deg) ',SatSubLON*radtodeg
      write(6,*)'***********************************'

      cstatus = 0
      npoints_out = 0
      do j=1,nyl
      do k=1,nxl
 
         r8lat=lat(k,j)*pi/180.d0
         r8lon=lon(k,j)*pi/180.d0

         call GPOINT(r8lat,r8lon,ELEV,SCAN,IERR)

         if(IERR.ne.0)then
c           write(6,*)'Error computing Elev/Scan in GPOINT from'
c           write(6,*)'Lat/Lon ', lat(k,j),lon(k,j)
            rline(k,j)=-2.
            rpix(k,j)=-2.
            cstatus=cstatus-1
         else

            call EVSC2L(INSTR,ELEV,SCAN,RL,RP)

            nrl=nint(rl)
            nrp=nint(rp)
            if( (nrl.gt.0.0.and.nrp.gt.0.0).and.
     &          (nrl.ge.start_line.and.nrp.ge.start_pix).and.
     &          (nrl.le.end_line.and.nrp.le.end_pix) )then

               rline(k,j)=(RL-start_line+rl_div)*rl_div
               rpix(k,j)= (RP-start_pix+rp_div )*rp_div
               i1=k
               j1=j
            else

               npoints_out = npoints_out+1
               rline(k,j)=r_missing_data
               rpix(k,j) =r_missing_data

            endif

         endif

      enddo
      enddo

      if(cstatus.lt.0)then

        write(6,*)'WARNING! Some rl/rp values not computed '
        write(6,*)'For LUT ',ctype,' status = ',cstatus

      endif

      if(npoints_out .gt. 0)then

         rnx=float(nxl)
         rny=float(nyl)
         rnp=float(npoints_out)
         write(6,*)'WARNING! Some rl/rp values out of bounds'
         rpct_out=rnp/(rnx*rny)
         if(rpct_out.gt.0.50)then
            write(6,*)'More than 50% of domain not covered ',rpct_out
            return
         else
            write(6,*)'Percent out of domain = ',rpct_out
         endif
        
      endif
c
c     do j = 1,nyl,20
c     do k = 1,nxl,20
c        write(6,*)'i,j,rpix,rline: ',k,j,rpix(k,j),rline(k,j)
c     enddo
c     enddo
c
c compute image resolution in meters
c
      if(rpix(nxl/2,nyl/2).ne.r_missing_data.and.
     &   rline(nxl/2,nyl/2).ne.r_missing_data)then
         i1=nxl/2
         j1=nyl/2
      endif
      call compute_sat_res_m(rp_div,rl_div,
     &rpix(i1,j1),rline(i1,j1),start_pix,start_line,instr,
     &rsndr_res_m,istatus)

      istatus = 0
      goto 1000

1000  return
      end
