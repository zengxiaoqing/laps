      Subroutine gen_gvarimage_lut(isat,jtype,kchl,
     &nx_l,ny_l,xlat,xlon,istatus)
c
      implicit none
c
      Include       'instco.inc'

      Integer     nx_l,ny_l
      Integer     isat
      Integer     jtype
      Integer     kchl

      Integer     idx
      parameter  (idx=2)

      Integer     nxl,nyl
      Integer     nx,ny
c
      real*4        xlat(nx_l,ny_l)
      real*4        xlon(nx_l,ny_l)
      real*4        lat(nx_l+idx,ny_l+idx)
      real*4        lon(nx_l+idx,ny_l+idx)

      real*4        rline(nx_l+idx,ny_l+idx)
      real*4        rpix(nx_l+idx,ny_l+idx)
      real*4        rel_ri(nx_l+idx,ny_l+idx)
      real*4        rel_rj(nx_l+idx,ny_l+idx)
      real*4        ri_laps(nx_l,ny_l)
      real*4        rj_laps(nx_l,ny_l)

      real*4        rl_abs(nx_l+idx,ny_l+idx)
      real*4        rp_abs(nx_l+idx,ny_l+idx)

      real*8        pi
      real*8        rl_div
      real*8        rp_div
      real*8        radtodeg
      real*4        r_img_res_m
      real*4        rls,rle,res,ree
      real*4        nwpixabs,nwlinabs
      real*4        nepixabs,nelinabs
      real*4        swpixabs,swlinabs
      real*4        sepixabs,selinabs

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
      Integer     start_line
      Integer     start_pix
      Integer     decimat
      Integer     imagedepth
      Integer     imagewidth
      Integer     i1,j1
      Integer     x_step
      Integer     y_step
      Integer     cstatus
      Integer     istatus
      Integer     ustatus
      Integer     wstatus
      Integer     IERR
      Integer     INSTR
      Integer     i,j,ii,jj
      Integer     n1,n2,nn,nc
      Integer     ils,ile
      Integer     nsCycles,nsIncs
      Integer     ewCycles,ewIncs
      Integer     time_spec(2)
      Integer     npoints_out
      Integer     ispec
      Integer     lend
      Integer     iwrite 
      data iwrite/1/

      Integer     linestart,lineend
      Integer     elemstart,elemend

      Integer     linestart_orig
      Integer     lineend_orig
      Integer     elemstart_orig
      Integer     elemend_orig

      Integer     idiff_org
      Integer     idiff_new
      Integer     jdiff_org
      Integer     jdiff_new
      integer     indx

      Character     filename*255
      Character     filename_sat*200
      Character     path*255
      Character     cname*100
      Character     cdir*100
      Character     table_path*200
      Character     cmode*10
      character     c_imc(4)*1
      character     ct*3,csattype*3
c
c The technique used here is to make a slightly larger domain
c (1 grid point larger on each side) and use this to compute
c the remapping lut. This insures that we have enough data on
c the edges when either interpolating or averaging during the
c remapping process (ie, in lvd_sat_ingest).
c
c include the static gvar navigation parameters via common.
c ---------------------------------------------------------
c
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      logical     lfirst(maxtype,maxsat)              !4 types x 2 sats (3-16-98)
      data lfirst /.false.,.false.,.false.,.false.,
     &             .false.,.false.,.false.,.false./
      save
c
c ---------------------------------------------------------
      istatus = 0
      nc=index(c_channel_types(kchl,jtype,isat),' ')-1
      if(nc.le.0)nc=3
      call lvd_file_specifier(c_channel_types(kchl,jtype,isat),
     &indx,istatus)

      csattype  = c_sat_types(jtype,isat)
      goto(61,62,63,62,62)indx

61       ct='vis'
         elemstart_orig = i_start_vis(jtype,isat)
         elemend_orig   = i_end_vis(jtype,isat)
         linestart_orig = j_start_vis(jtype,isat)
         lineend_orig   = j_end_vis(jtype,isat)
         goto 65

62       if(.not.lfirst(jtype,isat))then
            ct='ir'
            elemstart_orig = i_start_ir(jtype,isat)
            elemend_orig   = i_end_ir(jtype,isat)
            linestart_orig = j_start_ir(jtype,isat)
            lineend_orig   = j_end_ir(jtype,isat)
            lfirst(jtype,isat)=.true.
         else
            goto 1000
            istatus = 1
         endif 
         goto 65
 
63       ct='wv'
         elemstart_orig = i_start_wv(jtype,isat)
         elemend_orig   = i_end_wv(jtype,isat)
         linestart_orig = j_start_wv(jtype,isat)
         lineend_orig   = j_end_wv(jtype,isat)

65    continue
c
c get current nav parameters for file header
c
      call update_gvarimg_parms(c_sat_id(isat),
     &                         csattype,
     &             c_channel_types(kchl,jtype,isat),
     &             path_to_raw_sat(kchl,jtype,isat),
     &                         ewCycles,
     &                         ewIncs,
     &                         nsCycles,
     &                         nsIncs,
     &                         f_time,
     &                         c_imc,
     &                         start_pix,
     &                         start_line,
     &                         orbAt,
     &                         SatSubLAT,
     &                         SatSubLON,
     &                         decimat,
     &                         nx,ny,
     &                         imagewidth,imagedepth,
     &                         ustatus)

      if(ustatus.ne.0)then
         goto 901
      endif
c
c get expanded domain lats/lons
c
      nxl=nx_l+idx
      nyl=ny_l+idx

      call expand_domain(nx_l,ny_l,xlat,xlon,nxl,nyl,lat,lon,
     &istatus)
c
c if the gvar image data file ever changes, for example, data thinning is
c applied, then the following values should change as well.
c
      x_step=1.0
      y_step=1.0
      if(csattype.eq.'gwc')then
         if(ct(1:nc).eq.'vis')then
            x_step=float(decimat)*2.0   !the *2.0 is because afgwc only returns every other pixel!
            y_step=float(decimat)
         elseif(ct(1:nc).eq.'wv ')then
            x_step=float(decimat)*2.0
            y_step=float(decimat)
         else
            x_step=float(decimat)*2.0
            y_step=float(decimat)
         endif

         ewCycles=i_ewCycles(jtype,isat)
         ewIncs=i_ewIncs(jtype,isat)
         nsCycles=i_nsCycles(jtype,isat)
         nsIncs=i_nsIncs(jtype,isat)

      endif

      rp_div = 4.0*x_step
      rl_div = 4.0*y_step            !channels 2, 4, and 5 (3.9u, 11u, and 12u)
      if(ct(1:nc).eq.'wv ') then
         rl_div = 8.0*y_step         !channel 3 = water vapor
      elseif (ct(1:nc).eq.'vis') then
         rp_div = x_step             !channel 1 = visible
         rl_div = y_step
      endif

      INSTR=1          !1=Imager, 2=Sounder
      pi=3.141592653589793
      radtodeg=180.d0/pi

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
      do i=1,nxl
 
         r8lat=lat(i,j)*(pi/180.d0)
         r8lon=lon(i,j)*(pi/180.d0)

         call GPOINT(r8lat,r8lon,ELEV,SCAN,IERR)
         if(IERR.ne.0)then

c           write(6,*)'Error computing Elev/Scan in GPOINT from'
c           write(6,*)'Lat/Lon ', lat(i,j),lon(i,j)
            rline(i,j)=-2.
            rpix(i,j)=-2.
            cstatus=cstatus-1

         else

            call EVSC2L(INSTR,ELEV,SCAN,RL,RP)

c save corners
            if(i.eq.1.and.j.eq.1)then
               swpixabs=RP
               swlinabs=RL
            endif
            if(i.eq.1.and.j.eq.nyl)then
               nwpixabs=RP
               nwlinabs=RL
            endif
            if(i.eq.nxl.and.j.eq.1)then
               sepixabs=RP
               selinabs=RL
            endif
            if(i.eq.nxl.and.j.eq.nyl)then
               nepixabs=RP
               nelinabs=RL
            endif

            if( nint(rl).gt.0.0. and. nint(rl).le.25000 .and.
     &         nint(rp).gt.0.0 .and. nint(rp).le.30000) then

               rline(i,j)=(RL-start_line+rl_div)/rl_div 
               rpix(i,j)= (RP-start_pix+rp_div)/rp_div

            else

               npoints_out = npoints_out+1

            endif

         endif

      enddo
      enddo

      if(cstatus.lt.0)then
        write(6,*)'WARNING! Some rl/rp values not computed '
        write(6,*)'For LUT ',csattype,' status = ',cstatus
        goto 903
      endif

      if(npoints_out .gt. 0)then
         write(6,*)'WARNING! Some rl/rp values out of bounds'
         write(6,*)'These would be absolute sat coordinates'
      endif
c
c print corners
c
      write(6,*)'Absolute Sat Coords - Expanded Laps Domain'
      write(6,*)'------------------------------------------'
      write(6,120)nwpixabs,nwlinabs,nepixabs,nelinabs
      write(6,121)swpixabs,swlinabs,sepixabs,selinabs

120   format(5x,'NW Pix/Line ',2f10.4,10x,'NE Pix/Line ',2f10.4)
121   format(5x,'SW Pix/Line ',2f10.4,10x,'SE Pix/Line ',2f10.4)
 
      write(6,*)
c
c compute relative lut
c
      call get_sat_boundary(nxl,nyl,ny,nx,rpix,rline,
     &linestart,lineend,elemstart,elemend,
     &rls,rle,res,ree,istatus)
      if(istatus.ne.1)then
         write(6,*)'Fatal Error in get_sat_boundary - gen_gvarimg_lut'
         write(6,*)'Laps domain apparently outside satellite data!'
         goto 1000
      endif
c
      idiff_org = elemend_orig - elemstart_orig
      idiff_new = elemend-elemstart
      jdiff_org = lineend_orig - linestart_orig
      jdiff_new = lineend-linestart

      if(idiff_org.eq.idiff_new.and.jdiff_org.eq.jdiff_new)
     &then

         write(6,*)'gvar array sizes still the same'

      else

         write(6,*)'gvar array sizes have changed!'
         write(6,*)'Elem diffs: orig/new: ',idiff_org,idiff_new
         write(6,*)'Line diffs: orig/new: ',jdiff_org,jdiff_new

      endif

      if(elemstart_orig.ne.elemstart)then
         write(6,*)'Element Start .ne. Original'
      endif
      if(elemend_orig.ne.elemend)then
         write(6,*)'Element End .ne. Original'
      endif
      if(linestart_orig.ne.linestart)then
         write(6,*)'Line Start .ne. Original'
      endif
      if(lineend_orig.ne.lineend)then
         write(6,*)'Line End .ne. Original'
      endif
c
c compute ri, rj relative look up table for the block of data surrounding
c the laps domain.
c
      npoints_out = 0

      do j = 1,nyl
      do i = 1,nxl

         if(rpix(i,j).gt.0.0.and.rpix(i,j).lt.nx.and.
     &rline(i,j).gt.0.0.and.rline(i,j).lt.ny)then

            rel_ri(i,j) = rpix(i,j)  - float(elemstart)
            rel_rj(i,j) = rline(i,j) - float(linestart)

         else

            npoints_out = npoints_out + 1

         endif

      enddo
      enddo

      if(npoints_out .gt. 0)then
         write(6,*)'Relative pix/line coords not correct'
         write(6,*)'LUT not computed'
         goto  1000
      endif
c
c put the rel ri/rj lut into the laps domain size
c
      ils=idx/2+1
      ile=idx/2
      jj = 0
      do j = ils,nyl-ile
         jj = jj+1
         ii = 0
         do i = ils,nxl-ile
            ii = ii+1
            ri_laps(ii,jj) = rel_ri(i,j)
            rj_laps(ii,jj) = rel_rj(i,j)
         enddo
      enddo
c
      if(iwrite.eq.1)then

         do j = 1,ny_l,10
         do i = 1,nx_l,10
            write(6,*)'i,j,ri,rj: ',i,j,ri_laps(i,j),rj_laps(i,j)
         enddo
         enddo

      endif
c
c output
c
      call get_directory('static',path,lend)
      path=path(1:lend)//'lvd/'
      n1=index(path,' ')-1

      cname=path(1:n1)//c_sat_id(isat)//'-llij-'//ct(1:nc)
      n2=index(cname,' ')-1
      table_path = cname(1:n2)//'-'//csattype//'.lut'
      n2=index(table_path,' ')

      call write_table (table_path,nx_l,ny_l,xlat,xlon,
     &ri_laps,rj_laps,wstatus)
c
c compute image resolution in meters. This done with the original line/pix
c values since we use gimloc here.
c
      i1=nx_l/2
      j1=ny_l/2

      call compute_gvarimage_resolution(rp_div,rl_div,
     &rpix(i1,j1),rline(i1,j1),start_pix,start_line,
     &r_img_res_m,istatus)

      goto(71,72,73,72,72)indx

71       i_start_vis(jtype,isat)= elemstart
         i_end_vis(jtype,isat) = elemend
         j_start_vis(jtype,isat) = linestart
         j_end_vis(jtype,isat) = lineend
         image_depth_vis(jtype,isat) = imagedepth
         image_width_vis(jtype,isat) = imagewidth
         r_sat_sub_lat(isat) = SatSubLAT
         r_sat_sub_lon(isat) = SatSubLON
         
         goto 75

72       i_start_ir(jtype,isat) = elemstart
         i_end_ir(jtype,isat) = elemend
         j_start_ir(jtype,isat) = linestart
         j_end_ir(jtype,isat) = lineend
         image_depth_ir(jtype,isat) = imagedepth
         image_width_ir(jtype,isat) = imagewidth
         r_sat_sub_lat(isat) = SatSubLAT
         r_sat_sub_lon(isat) = SatSubLON

         goto 75

73       i_start_wv(jtype,isat) = elemstart
         i_end_wv(jtype,isat) = elemend
         j_start_wv(jtype,isat) = linestart
         j_end_wv(jtype,isat) = lineend
         image_depth_wv(jtype,isat) = imagedepth
         image_width_wv(jtype,isat) = imagewidth
         r_sat_sub_lat(isat) = SatSubLAT
         r_sat_sub_lon(isat) = SatSubLON

75    continue
 
      call write_orb_att(path,c_sat_id(isat),336,orbAt)

      istatus = 1

      goto 1000

900   write(6,*)'Error opening ascii test file'
      istatus=-1
      goto 1000

901   write(6,*)'Error in update_gvarimg_parmfile'
      istatus = -1
      goto 1000

903   write(6,*)'Error computing rline/rpix - terminating'
      istatus = -1
      goto 1000

905   write(6,*)'Error writing gvarimg parmfile'
      istatus = -1
      goto 1000

906   write(6,*)'Error writing satsector file'
      goto 1000

907   write(6,*)'Error - lvd_file_specifier'
      goto 1000

908   write(6,*)'Error in update_gvarimage_parmfile'
      istatus = -1
      goto 1000

999   write(6,*)'No update to parameter file: ',cmode

1000  return
      end
