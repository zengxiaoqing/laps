      Subroutine gen_gvarimage_lut(isat,jtype,kchl,
     &nx_l,ny_l,lat,lon,istatus)
c
      implicit none
c
      Include       'instco.inc'

      Integer     nx_l,ny_l
      Integer     isat
      Integer     jtype
      Integer     kchl

      Integer       idx
      parameter    (idx=2)

      Integer       idx2
      Integer       nxl,nyl
      Integer       nxl2,nyl2

      Integer       nx,ny
c
      real*4        lat(nx_l,ny_l)
      real*4        lon(nx_l,ny_l)

      real*4        xlat (nx_l+idx,ny_l+idx)
      real*4        xlon (nx_l+idx,ny_l+idx)

      real*4        xlat2 (nx_l+2*idx,ny_l+2*idx)
      real*4        xlon2 (nx_l+2*idx,ny_l+2*idx)
      real*4        rline (nx_l+2*idx,ny_l+2*idx)
      real*4        rpix  (nx_l+2*idx,ny_l+2*idx)
      real*4        rel_ri(nx_l+2*idx,ny_l+2*idx)
      real*4        rel_rj(nx_l+2*idx,ny_l+2*idx)
      real*4        rl_abs(nx_l+2*idx,ny_l+2*idx)
      real*4        rp_abs(nx_l+2*idx,ny_l+2*idx)

      real*4        ri_laps(nx_l,ny_l)
      real*4        rj_laps(nx_l,ny_l)

      real*8        pi
      real*8        rl_div
      real*8        rp_div
      real*8        radtodeg
      real*4        r_img_res_m
      real*4        rls,rle,res,ree
      real*4        r_thin
      real*4        nwpixabs,nwlinabs
      real*4        nepixabs,nelinabs
      real*4        swpixabs,swlinabs
      real*4        sepixabs,selinabs
      real*4        r_missing_data

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
      Integer     i1,j1
      Integer     x_step
      Integer     y_step
      Integer     xres,yres
      Integer     cstatus
      Integer     istatus
      Integer     ustatus
      Integer     wstatus
      Integer     IERR
      Integer     instr
      Integer     i,j,ii,jj
c     Integer     n1,n2,nn,nc
      Integer     n1,n2,nc
      Integer     ils,ile
      Integer     nsCycles,nsIncs
      Integer     ewCycles,ewIncs
      Integer     time_spec(2)
      Integer     npoints_out
      Integer     nijout
c     Integer     ispec
      Integer     lend

      logical     lwrite 
      data lwrite/.false./

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

c     Character     filename*255
c     Character     filename_sat*200
      Character     path*255
      Character     cname*100
c     Character     cdir*100
      Character     table_path*200
c     Character     cmode*10
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

      logical     lpoint
      logical     lfirst(maxtype,maxsat)              !4 types x 4 sats (5-26-98)
      data lfirst /.false.,.false.,.false.,.false.,
     &             .false.,.false.,.false.,.false.,
     &             .false.,.false.,.false.,.false.,
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
c
c only need to generate IR lut once
c
            return
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
     &                         csattype,l_cell_afwa,
     &             c_channel_types(kchl,jtype,isat),
     &             path_to_raw_sat(kchl,jtype,isat),
     &                         ewCycles,ewIncs,
     &                         nsCycles,nsIncs,
     &                         f_time,
     &                         c_imc,
     &                         xres,yres,
     &                         start_pix,start_line,
     &                         orbAt,
     &                         SatSubLAT,SatSubLON,
     &                         decimat,
     &                         nx,ny,
     &                         ustatus)

      if(ustatus.ne.0)then
         if(ustatus.eq.1)then
            print*,'GVAR parameters not obtained cannot proceed'
            istatus =-1
            return
         else
            goto 901
         endif
      endif
c
c get expanded domain lats/lons
c
      idx2=idx*2
      nxl=nx_l+idx
      nyl=ny_l+idx
      nxl2=nx_l+idx2
      nyl2=ny_l+idx2

      call expand_domain(nx_l,ny_l,lat,lon,nxl,nyl,xlat,xlon,
     +istatus)
      call expand_domain(nxl,nyl,xlat,xlon,nxl2,nyl2,xlat2,xlon2,
     +istatus)
c
c if the gvar image data file ever changes, for example, data thinning is
c applied, then the following values should change as well.
c
      x_step=1.0
      y_step=1.0
c
c currently AFWA goespatch is every other pixel and every scan line.
c
      r_thin=2.0

      if(csattype.eq.'gwc')then
         if(ct(1:nc).eq.'vis')then
            x_step=float(decimat)*r_thin   ! *2.0 because afwa returns every other pixel!
            y_step=float(decimat)
         elseif(ct(1:nc).eq.'wv ')then
            x_step=float(decimat)*r_thin
            y_step=float(decimat)
         else
            x_step=float(decimat)*r_thin
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

      instr=1          !1=Imager, 2=Sounder
      pi=3.141592653589793
      radtodeg=180.d0/pi

      call bcd_to_int(orbAt(12),time_spec)
      time_50 = time50(time_spec)
      t50_8=time_50
      t = f_time /60. + 7305. * 24. * 60.

      call SETCON(instr,nsCycles,nsIncs,ewCycles,ewIncs)
      call LMODEL(t,t50_8,OrbAt,imc,SatSubLAT,SatSubLON)

      write(6,*)'Sat Subpoint lat (deg) ',SatSubLAT*radtodeg
      write(6,*)'Sat Subpoint lon (deg) ',SatSubLON*radtodeg
      write(6,*)'***********************************'

      cstatus = 0
      npoints_out = 0

      do j=1,nyl2
      do i=1,nxl2
 
         r8lat=xlat2(i,j)*(pi/180.d0)
         r8lon=xlon2(i,j)*(pi/180.d0)

         call GPOINT(r8lat,r8lon,ELEV,SCAN,IERR)
         if(IERR.ne.0)then

c           write(6,*)'Error computing Elev/Scan in GPOINT from'
c           write(6,*)'Lat/Lon ', xlat2(i,j),xlon2(i,j)
            rline(i,j)=-2.
            rpix(i,j)=-2.
            cstatus=cstatus-1

         else

            call EVSC2L(instr,ELEV,SCAN,RL,RP)

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

            if( idnint(rl).gt.0 .and. idnint(rl).le.25000 .and.
     &          idnint(rp).gt.0 .and. idnint(rp).le.30000) then

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
      call get_sat_boundary(nxl2,nyl2,ny,nx,rpix,rline,
     &linestart,lineend,elemstart,elemend,
     &rls,rle,res,ree,istatus)
      if(istatus.ne.1)then
         write(6,*)'Note: Laps domain outside satellite data!'
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
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)return

      nijout = 0
      lpoint =.true.
      do j = 1,nyl2
      do i = 1,nxl2
         if(rpix(i,j).ne.r_missing_data.and.
     &rline(i,j).ne.r_missing_data)then
            rel_ri(i,j) = rpix(i,j)  - res
            rel_rj(i,j) = rline(i,j) - rls 
            if(lpoint)then
               i1=i
               j1=j
               lpoint=.false.
            endif
         else
            rel_ri(i,j) = r_missing_data
            rel_rj(i,j) = r_missing_data
            nijout=nijout+1
         endif
      enddo
      enddo
c
c put the rel ri/rj lut into the laps domain size
c
      ils=idx2/2+1
      ile=idx2/2
      jj = 0
      do j = ils,nyl2-ile
         jj = jj+1
         ii = 0
         do i = ils,nxl2-ile
            ii = ii+1
            ri_laps(ii,jj) = rel_ri(i,j)
            rj_laps(ii,jj) = rel_rj(i,j)
         enddo
      enddo
c
      if(nijout.gt.0)then
         print*,'Found ',nijout,' points outside domain'
      endif

      if(lwrite)then

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

      call write_table (table_path,nx_l,ny_l,lat,lon,
     &ri_laps,rj_laps,wstatus)

      r_sat_sub_lat(isat) = SatSubLAT
      r_sat_sub_lon(isat) = SatSubLON

      if(elemstart.le.0)elemstart=1
      if(elemend.gt.nx)elemend=nx
      if(linestart.le.0)linestart=1
      if(lineend.gt.ny)lineend=ny
c
c ------------------------------------------------------------------------
c compute image resolution in meters. This done with the original line/pix
c values since we use gimloc here.
c ------------------------------------------------------------------------
c
      call compute_sat_res_m(rp_div,rl_div,
     &rpix(i1,j1),rline(i1,j1),start_pix,start_line,
     &instr,r_img_res_m,istatus)

      goto(71,72,73,72,72)indx

71       i_start_vis(jtype,isat)= elemstart
         i_end_vis(jtype,isat) = elemend
         j_start_vis(jtype,isat) = linestart
         j_end_vis(jtype,isat) = lineend
         r_resolution_x_vis(jtype,isat) = r_img_res_m  !float(xres)
         r_resolution_y_vis(jtype,isat) = r_img_res_m  !float(yres)
         n_pixels_vis(jtype,isat) = nx
         n_lines_vis(jtype,isat)  = ny

         goto 75

72       i_start_ir(jtype,isat) = elemstart
         i_end_ir(jtype,isat) = elemend
         j_start_ir(jtype,isat) = linestart
         j_end_ir(jtype,isat) = lineend
         r_resolution_x_ir(jtype,isat) = r_img_res_m  !float(xres)
         r_resolution_y_ir(jtype,isat) = r_img_res_m  !float(yres)
         n_pixels_ir(jtype,isat) = nx
         n_lines_ir(jtype,isat)  = ny

         goto 75

73       i_start_wv(jtype,isat) = elemstart
         i_end_wv(jtype,isat) = elemend
         j_start_wv(jtype,isat) = linestart
         j_end_wv(jtype,isat) = lineend
         r_resolution_x_wv(jtype,isat) = r_img_res_m  !float(xres)
         r_resolution_y_wv(jtype,isat) = r_img_res_m  !float(yres)
         n_pixels_wv(jtype,isat) = nx
         n_lines_wv(jtype,isat)  = ny

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

1000  return
      end
