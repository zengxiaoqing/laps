      Subroutine gen_gvarimage_lut(isat,jtype,kchl,
     &nx_l,ny_l,lat,lon,ri_laps,rj_laps,istatus)
c
      implicit none
c
      Include       'instco.inc'

      Integer     nx_l,ny_l
      Integer     isat
      Integer     jtype
      Integer     kchl

      Integer       idx
      parameter    (idx=20)

      Integer       nxe,nye
      Integer       nx,ny
c
      real        lat(nx_l,ny_l)
      real        lon(nx_l,ny_l)

      real        xlat (nx_l+idx,ny_l+idx)
      real        xlon (nx_l+idx,ny_l+idx)
      real        rline (nx_l+idx,ny_l+idx)
      real        rpix  (nx_l+idx,ny_l+idx)
      real        rel_ri(nx_l+idx,ny_l+idx)
      real        rel_rj(nx_l+idx,ny_l+idx)
      real        rl_abs(nx_l+idx,ny_l+idx)
      real        rp_abs(nx_l+idx,ny_l+idx)

      real        ri_laps(nx_l,ny_l)
      real        rj_laps(nx_l,ny_l)

      real        xmn(nx_l+idx),ymn(ny_l+idx)
      real        xtn(nx_l+idx),ytn(ny_l+idx)

      real*8        pi
      real*8        rl_div, rp_div
      real*8        radtodeg, degtorad
      real        r4lat,r4lon
      real        rls,rle,res,ree
      real        jdx,jdy
      real        r_thin
      real        nwpixabs,nwlinabs
      real        nepixabs,nelinabs
      real        swpixabs,swlinabs
      real        sepixabs,selinabs
      real        rmetx,rmety
      real        golonsbp,golatsbp,goalpha
      real        grid_spacing_proj_m
      real        grid_spacing_actual_mx
      real        grid_spacing_actual_my
      real        grid_spacing_m
      real        gssumx,gssumy
      real        mdlat,mdlon
      real        sat_res_m,erad
      real        deltax,deltay
      real        factor
      real        r_ratio
      real        lat1,lat2,lon0
      real        r_missing_data
      real        time_50,time50

      real*8        r8lat,r8lon
      real*8        ELEV,SCAN
      real*8        RL, RP
      real*8        orbAt(336)
      real*8        t50_8, t, f_time
      real*8       SatSubLAT,SatSubLON
c
      Integer     start_line
      Integer     start_pix
      Integer     i1,j1
      Integer     x_step
      Integer     y_step
      Integer     xres,yres
      Integer     cstatus
      Integer     istatus
      Integer     gstatus
      Integer     ustatus
      Integer     wstatus
      Integer     IERR
      Integer     expansion,ifactor
      Integer     instr
      Integer     i,j,ii,jj
      Integer     n,n1,n2,nc,nl
      Integer     ils,ile
      Integer     nsCycles,nsIncs
      Integer     ewCycles,ewIncs
      Integer     time_spec(2)
      Integer     npoints_out
      Integer     nijout
      Integer     lend
      Integer     strpix,strscnl,stppix,stpscnl 
      Integer     reqobstm,decimat
      Integer     istrbdy1,istrbdy2,istpbdy1,istpbdy2
      Integer     bepixfc,bescnfc,iwidth,idepth,fsci

      logical     lwrite 
      data lwrite/.true./

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

      Real        rdum
      Integer     strtpix,strtline
      Integer     stoppix,stopline

      Character     cfilename*255
c     Character     filename_sat*200
      Character     path*255
      Character     cname*100
      Character     c_afwa_fname*100
      Character     cdir*200
      Character     table_path*200
c     Character     cmode*10
      character     c_imc(4)*1
      character     ct*3,csattype*3,cty*3
      Character     image_type*2
      character     c6_maproj_ret*6
      character     cchantype*3
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

c     logical     lfirst(maxtype,maxsat)              !4 types x 4 sats (5-26-98) 6 sats (1-31-03)
c     data lfirst /.false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false.,
c    &             .false.,.false.,.false.,.false./
c     save lfirst
c
c ---------------------------------------------------------
      write(6,*)' subroutine gen_gvarimage_lut:'

      istatus = 0
      cchantype = c_channel_types(kchl,jtype,isat)
      call s_len(cchantype,nc)
      if(nc.le.0)nc=3
      call lvd_file_specifier(cchantype,indx,istatus)

      csattype  = c_sat_types(jtype,isat)

      if(indx.eq.1)then

         ct='vis'
         elemstart_orig = i_start_vis(jtype,isat)
         elemend_orig   = i_end_vis(jtype,isat)
         linestart_orig = j_start_vis(jtype,isat)
         lineend_orig   = j_end_vis(jtype,isat)

      elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then

c        if(.not.lfirst(jtype,isat))then
            ct='ir'
            elemstart_orig = i_start_ir(jtype,isat)
            elemend_orig   = i_end_ir(jtype,isat)
            linestart_orig = j_start_ir(jtype,isat)
            lineend_orig   = j_end_ir(jtype,isat)
c           lfirst(jtype,isat)=.true.
c        else
c
c only need to generate IR lut once. Removed by JS: 10-5-06
c
c           return

c        endif 

      elseif(indx.eq.3)then
 
         ct='wv'
         elemstart_orig = i_start_wv(jtype,isat)
         elemend_orig   = i_end_wv(jtype,isat)
         linestart_orig = j_start_wv(jtype,isat)
         lineend_orig   = j_end_wv(jtype,isat)

      endif
c
c get current nav parameters for file header
c
      call update_gvarimg_parms(c_sat_id(isat),
     &                         csattype,l_cell_afwa,
     &                         cchantype,
     &             path_to_raw_sat(kchl,jtype,isat),
     &                         ewCycles,ewIncs,
     &                         nsCycles,nsIncs,
     &                         f_time,
     &                         c_imc,
     &                         xres,yres,
     &                         start_pix,start_line,
     &                         orbAt,
     &                         SatSubLAT,SatSubLON,
     &                         x_step,y_step,
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

      write(6,*)' Returned from update_gvarimg_parms'
      write(6,*)' xres/yres = ',xres,yres
      write(6,*)' x_step/y_step = ',x_step,y_step
c
c get expanded domain lats/lons. Some of this code follows what
c happens in gridgen_model.
c
      if(xres.eq.0.0.or.yres.eq.0.0)then
         xres=4.
         yres=4.
         if(ct(1:nc).eq.'vis')then
            xres=1.
            yres=1.
         elseif(ct(1:nc).eq.'wv ')then
            xres=8.
            yres=8.
         endif
         write(6,*)' WARNING: Reset values from update_gvarimg_parms'
         write(6,*)' xres/yres = ',xres,yres
      endif

      call get_grid_spacing(grid_spacing_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error return from get_grid_spacing'     
          return
      endif

      if(.false.)then
        gssumx=0
        gssumy=0
        do j=1,ny_l
        do i=1,nx_l
           call get_grid_spacing_actual_xy(lat(i,j),lon(i,j)
     1             ,grid_spacing_actual_mx,grid_spacing_actual_my
     1             ,istatus)
           gssumx=gssumx+grid_spacing_actual_mx
           gssumy=gssumy+grid_spacing_actual_my
        enddo
        enddo
        grid_spacing_m=(gssumx/(nx_l*ny_l)+gssumy/(nx_l*ny_l))/2.
      endif

      print*,'Grid spacing (m): ',grid_spacing_m
      if(grid_spacing_m.le.0.0)then
        print*,'Oh no, grid spacing = 0'
        istatus = -1
        return
      endif

      sat_res_m=((xres+yres)*1000.)/2.0

      r_ratio=sat_res_m/grid_spacing_m
      expansion=(nint(r_ratio)+1)*2
      if(expansion.lt.4)expansion = 4
      if(expansion.gt.20)expansion=20
      nxe=nx_l+expansion
      nye=ny_l+expansion

      call get_grid_center(mdlat,mdlon,istatus)
      if(istatus .ne. 1)then
         write(6,*)' Error calling laps routine'
         return
      endif
      write(6,*)' grid_center = ',mdlat,mdlon

      call get_c6_maproj(c6_maproj_ret,istatus)
      if(istatus .ne. 1)then
         print*,'Error return from get_c6_maproj'
         return
      endif
      call downcase(c6_maproj_ret,c6_maproj_ret) 

      call get_standard_latitudes(lat1,lat2,istatus)
      call get_standard_longitude(lon0,istatus)

      if(c6_maproj_ret .eq. 'plrstr')then
         call get_ps_parms(lat1,lat2,grid_spacing_m,lon0
     1                       ,grid_spacing_proj_m)

         deltax = grid_spacing_proj_m

      else
         deltax = grid_spacing_m
      endif

      deltay=deltax
      call POLAR_GPG(mdlat,mdlon,XMN(1),YMN(1),deltax,deltay,nxe,nye)

      do i=2,nxe
         xmn(i)=xmn(i-1)+deltax
      enddo
      xmn(nxe)=2*xmn(nxe-1)-xmn(nxe-2)
      do j=2,nye
         ymn(j)=ymn(j-1)+deltay
      enddo
      ymn(nye)=2*ymn(nye-1)-ymn(nye-2)
      do i=2,nxe
         xtn(i)=.5*(xmn(i)+xmn(i-1))
      enddo
      xtn(1)=1.5*xmn(1)-.5*xmn(2)
      do j=2,nye
         ytn(j)=.5*(ymn(j)+ymn(j-1))
      enddo
      ytn(1)=1.5*ymn(1)-.5*ymn(2)

      call get_earth_radius(erad,istatus)
      if(istatus .ne. 1)then
         write(6,*)' Error calling get_earth_radius'
         return
      endif

      do j=1,nye
         do i=1,nxe
            call xy_to_latlon(xtn(i),ytn(j),erad,xlat(i,j),xlon(i,j))
         enddo
      enddo

c if the gvar image data file ever changes, for example, data thinning is
c applied, then the following values should change as well. These now obtained
c from gvar data files.
c
c     x_step=1.0
c     y_step=1.0
c
c currently AFWA goespatch is every other pixel and every scan line.
c use r_thin to adjust appropriately
c
      print*,'Earth Location to line/elem look-up-table: '
     &,c_sat_id(isat)
      print*

      if(c_sat_id(isat)(1:4).eq.'goes')then

        r_thin = 1.0
        if(csattype.eq.'gwc')r_thin=2.0
        x_step=x_step*r_thin

c JS: Test using the 4 vars below taken directly from the netCDF file
c       ewCycles=i_ewCycles(jtype,isat)
c       ewIncs=i_ewIncs(jtype,isat)
c       nsCycles=i_nsCycles(jtype,isat)
c       nsIncs=i_nsIncs(jtype,isat)

        if(.true.)then ! new experimental code
            rp_div = xres
            rl_div = yres            !channels 2, 4, and 5 (3.9u, 11u, and 12u)
        else            ! original code
            rp_div = 4.0*x_step
            rl_div = 4.0*y_step      !channels 2, 4, and 5 (3.9u, 11u, and 12u)

            if(ct(1:nc).eq.'wv '.and.csattype.ne.'gwc') then
              rl_div = 8.0*y_step          !channel 3 = water vapor; only FSL public
            elseif (ct(1:nc).eq.'vis') then
              rp_div = x_step              !channel 1 = visible
              rl_div = y_step
            endif
        endif

        instr=1          !1=Imager, 2=Sounder
        pi=3.141592653589793
        radtodeg=180.d0/pi
        degtorad=1./radtodeg

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

        do j=1,nye
        do i=1,nxe
 
          r8lat=xlat(i,j)*degtorad
          r8lon=xlon(i,j)*degtorad

          call GPOINT(r8lat,r8lon,ELEV,SCAN,IERR)
          if(IERR.ne.0)then

c           write(6,*)'Error computing Elev/Scan in GPOINT from'
c           write(6,*)'Lat/Lon ', xlat(i,j),xlon(i,j)
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
             if(i.eq.1.and.j.eq.nye)then
               nwpixabs=RP
               nwlinabs=RL
             endif
             if(i.eq.nxe.and.j.eq.1)then
                sepixabs=RP
                selinabs=RL
             endif
             if(i.eq.nxe.and.j.eq.nye)then
                nepixabs=RP
                nelinabs=RL
             endif

             if( idnint(rl).gt.0 .and. idnint(rl).le.25000 .and.
     &           idnint(rp).gt.0 .and. idnint(rp).le.30000) then

                rline(i,j)=(RL-start_line+rl_div)/rl_div 
                rpix(i,j)= (RP-start_pix+rp_div)/rp_div

             else

                npoints_out = npoints_out+1

             endif

          endif

        enddo
        enddo

      else
c
c METEOSAT
c
        cdir=path_to_raw_sat(kchl,jtype,isat)
        cty=cchantype                             
        n=index(cdir,' ')-1
        cname=c_afwa_fname(c_sat_id(isat),cty)
        call s_len(cname,nl)
        cfilename=cdir(1:n)//cname(1:nl)

        call read_gwc_header(cfilename,l_cell_afwa,strtpix,
     &strtline,stoppix,stopline,reqobstm,image_type,golatsbp,golonsbp,
     &iwidth,idepth,goalpha,istrbdy1,istrbdy2,istpbdy1,
     &istpbdy2,bepixfc,bescnfc,fsci,decimat,gstatus)
        if(gstatus.ne.0)then
           write(6,*)'Error in read_gwc_header'
           istatus=-1
           goto 900
        endif

        pi=3.141592653589793
        radtodeg=180.d0/pi
        degtorad=pi/180.

        do j=1,nye
        do i=1,nxe

          call fc01_conv_pts_el_to_met_1(golonsbp,
     &                        stoppix,fsci,
     &                        strtline,decimat,
     &                        cty,
     &                        xlat(i,j),
     &                        xlon(i,j),
     &                        rmetx,
     &                        rmety,
     &                        istatus)

          rline(i,j)=rmety
          rpix(i,j)= rmetx

          if(istatus.ne.0)then
             print*,'Error in fc01_conv_pts_el_to_met_1'
             goto 903
          endif

        enddo
        enddo

      endif

      if(lwrite)then
         write(6,*)'Satellite points of LAPS domain corners'
         print*,   '---------------------------------------'
         print*,'lat/lon/x/y (1,1) = ',xlat(1,1),xlon(1,1),
     .rpix(1,1),rline(1,1)
         print*,'lat/lon/x/y (nx,1) = ',xlat(nxe,1),xlon(nxe,1),
     .rpix(nxe,1),rline(nxe,1)
         print*,'lat/lon/x/y (1,ny) = ',xlat(1,nye),xlon(1,nye),
     .rpix(1,nye),rline(1,nye)
         print*,'lat/lon/x/y (nx,ny) = ',xlat(nxe,nye),xlon(nxe,nye)
     .,rpix(nxe,nye),rline(nxe,nye)

!        print*,'Absolute Satellite Pix/Line Coords:'
!        print*,'-----------------------------------'
!        do j = 1,nye,100
!        do i = 1,nxe,100
!           write(6,*)'i,j,ri,rj: ',i,j,rpix(i,j),rline(i,j)
!        enddo
!        enddo

      endif

c     if(c_sat_id(isat) .eq. 'meteos')goto 777

      if(cstatus.lt.0)then
        write(6,*)'WARNING! Some rl/rp values not computed '
        write(6,*)'For Look-Up-Table ',csattype
     +,' status = ',cstatus
        write(6,*)'Earth location may not be viewable from',
     +' satellite'
        write(6,*)'Use data with caution!'

c       goto 903
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
      write(6,*)' Call get_sat_boundary...'
      call get_sat_boundary(nx_l,ny_l,nxe,nye,idx,ny,nx
     &,rpix,rline,linestart,lineend,elemstart,elemend,
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
c compute ri, rj relative look up table for expanded domain.
c
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)return

      nijout = 0
      lpoint =.true.
      do j = 1,nye
      do i = 1,nxe
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
c adjust relative ri/rj for laps domain
c
      factor=expansion/2.
      ils=int(factor)+1
      ile=int(factor+0.5)
      jj = 0
      do j = ils,nye-ile
         jj = jj+1
         ii = 0
         do i = ils,nxe-ile
            ii = ii+1
            ri_laps(ii,jj) = rel_ri(i,j)
            rj_laps(ii,jj) = rel_rj(i,j)
         enddo
      enddo
c
      if(nijout.gt.0)then
         print*,'Found ',nijout,' points outside domain'
      endif

777   if(lwrite)then

         write(6,*)' Domain corner ri/rj_laps info:'
         do j = 1,ny_l,ny_l-1
         do i = 1,nx_l,nx_l-1
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

c     call write_table (table_path,nx_l,ny_l,lat,lon,
c    &ri_laps,rj_laps,wstatus)

      r_sat_sub_lat(isat) = SatSubLAT
      r_sat_sub_lon(isat) = SatSubLON

      if(elemstart.le.0)elemstart=1
      if(elemend.gt.nx)elemend=nx
      if(linestart.le.0)linestart=1
      if(lineend.gt.ny)lineend=ny

      print*,' 5-8-08 '
      print*,'-----------------------------------'
      print*,'*** NOT USING compute_sat_res_m ***'
      print*,'-----------------------------------'
c
c ------------------------------------------------------------------------
c compute image resolution in meters. This done with the original line/pix
c values since we use gimloc here.
c ------------------------------------------------------------------------
c
c     if(c_sat_id(isat)(1:4).eq.'goes')then
c        call compute_sat_res_m(rp_div,rl_div,
c    &rpix(i1,j1),rline(i1,j1),start_pix,start_line,
c    &instr,sat_res_m,istatus)
c     else
c        sat_res_m=5000.0*decimat
c     endif

      if(indx.eq.1)then

         i_start_vis(jtype,isat)= elemstart
         i_end_vis(jtype,isat) = elemend
         j_start_vis(jtype,isat) = linestart
         j_end_vis(jtype,isat) = lineend
         r_resolution_x_vis(jtype,isat) = sat_res_m  !float(xres)
         r_resolution_y_vis(jtype,isat) = sat_res_m  !float(yres)
         n_pixels_vis(jtype,isat) = nx
         n_lines_vis(jtype,isat)  = ny

      elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then

         i_start_ir(jtype,isat) = elemstart
         i_end_ir(jtype,isat) = elemend
         j_start_ir(jtype,isat) = linestart
         j_end_ir(jtype,isat) = lineend
         r_resolution_x_ir(jtype,isat) = sat_res_m  !float(xres)
         r_resolution_y_ir(jtype,isat) = sat_res_m  !float(yres)
         n_pixels_ir(jtype,isat) = nx
         n_lines_ir(jtype,isat)  = ny

      elseif(indx.eq.3)then

         i_start_wv(jtype,isat) = elemstart
         i_end_wv(jtype,isat) = elemend
         j_start_wv(jtype,isat) = linestart
         j_end_wv(jtype,isat) = lineend
         r_resolution_x_wv(jtype,isat) = sat_res_m  !float(xres)
         r_resolution_y_wv(jtype,isat) = sat_res_m  !float(yres)
         n_pixels_wv(jtype,isat) = nx
         n_lines_wv(jtype,isat)  = ny

      endif
 
      if(.false.)call write_orb_att(path,c_sat_id(isat),336,orbAt)

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

c -------------------

       SUBROUTINE POLAR_GPG(LAT,LON,X,Y,DX,DY,NX,NY)
C
      include 'trigd.inc'
       REAL   LAT,LON,X,Y,DX,DY,
     1        ERAD,TLAT,TLON                                      ! ,PLAT,PLON,
     1        XDIF,YDIF
C
       INTEGER   NX,NY
C
       RAD=3.141592654/180.

       call get_earth_radius(erad,istatus)
       if(istatus .ne. 1)then
           write(6,*)' Error calling get_earth_radius'
           stop
       endif

       TLAT=90.0
       call get_standard_longitude(std_lon,istatus)
       if(istatus .ne. 1)then
           write(6,*)' Error calling laps routine'
           stop
       endif
       TLON=std_lon
C
C      CALL GETOPS(PLAT,PLON,LAT,LON,TLAT,TLON)
C      CALL PSTOXY(XDIF,YDIF,PLAT,PLON,ERAD)

C      call latlon_to_xy(LAT,LON,TLAT,TLON,ERAD,XDIF,YDIF)
       call latlon_to_xy(LAT,LON,ERAD,XDIF,YDIF)

C
       X=XDIF+(1.-FLOAT(NX)/2.)*DX
       Y=YDIF+(1.-FLOAT(NY)/2.)*DY
C
       RETURN
C
       END
