      subroutine readlut(csat_id,csat_typ,maxch,nch,
     &chtype,nx,ny,ri,rj,istatus)
c
      implicit none

      integer nx,ny
      integer maxch
      integer nch
      integer istatus
      integer i,n
      integer lend,lenf
      integer ispec
      integer istat

      logical   lgot_lut(maxch)

      real    ri(nx,ny,maxch)
      real    rj(nx,ny,maxch)
c     real    ri_in(nx,ny)
c     real    rj_in(nx,ny)
c     real    rdummy(nx,ny)

      character*100 cpath
      character*255 file
      character*6   csat_id
      character*3   csat_typ
      character*3   chtype(maxch)
      character*3   ct
c
c-----------------------------------------------------------------------
c
      istatus = 1

      call get_directory('static',cpath,lend)
      cpath=cpath(1:lend)//'lvd/'
      lend=index(cpath,' ')-1

c JS: 8-31-06. No longer reading a mapping look-up-table (lut).
c  
      do i=1,maxch
         lgot_lut(i)=.true.   !!!.false.
      enddo

      do i=1,nch

         call lvd_file_specifier(chtype(i),ispec,istat)
         if(istat.eq.0)then
c           if(.not.lgot_lut(ispec))then
c              lgot_lut(ispec)=.true.
c              ct=chtype(i)
               if(ispec.eq.2.or.ispec.eq.4.or.ispec.eq.5)then  !all ir channels have the same ri/rj values
                  call move(rj(1,1,2),rj(1,1,ispec),nx,ny)
                  call move(ri(1,1,2),ri(1,1,ispec),nx,ny)
c                 ct='ir'
               endif

c              n=index(ct,' ')-1
c              if(n.le.0)n=3
c              file=cpath(1:lend)//csat_id//'-llij-'
c              lenf=index(file,' ')-1
c              file=file(1:lenf)//ct(1:n)//'-'//csat_typ//'.lut'
c              n=index(file,' ')-1
c              open(12,file=file,
c    &form='unformatted',status='old',err=101)
c              write(6,*)'Reading ',file(1:n)
c              read(12,err=23,end=23) rdummy
c              read(12,err=23,end=23) rdummy
c              read(12,err=23,end=23) ri_in
c              read(12,err=23,end=23) rj_in
c              close (12)
c              call move(ri_in,ri(1,1,ispec),nx,ny)
c              call move(rj_in,rj(1,1,ispec),nx,ny)
c           else
c              print*,'Not reading mapping look-up-table'
c              print*,'Return without ri/rj values'
c           endif

         endif
      enddo

      istatus = 0
      goto 1000

23    write(6,*)'Error reading or eof ll/ij lookup table'
      close(12)
      goto 1000

101   write(6,*)'Error opening file ',file(1:n)

1000  return
      end
c
c================================================
c
      subroutine check_luts(cfname_cur,isat,jtype,
     &chtype,maxchannels,nchannels,l_lut_flag,istatus)
c
c most applications will relocalize which will delete any
c existing lut's and new ones are automatically generated.
c For WFO and AFWA, the satellite data files can possibly
c change. This is particularly true of AFWA satellite data.
c
      implicit none

      integer        maxchannels
      integer        nchannels

      character      cfname_cur*(*)
      character      chtype(maxchannels)*3
      character      cname*100
      character      c_afwa_fname*100
      character      cfname*200

      integer        i,il
      integer        isat,jtype,kch
      integer        istatus
      integer        nlin,npix
      integer        nx,ny
      integer        ispec
      integer        ilat00,ilon00
      integer        i_la1,i_lo1
      integer        idx,idy,iresx,iresy
      integer        nw_line,nw_pix
      integer        nxaf,nyaf
      integer        lenf

      integer        strpix,strscnl,stppix,stpscnl
      integer        reqobstm,imgtype
      integer        iwidth,idepth
      integer        istrbdy1,istrbdy2,istpbdy1,istpbdy2
      integer        bepixfc,bescnfc,fsci,idecimat

      integer        istrtline,istrtpix
      integer        nw_vis_line_gwc
      integer        nw_vis_pix_gwc

      logical        l_lut_flag

      real resx,resy
      real centerlat,centerlon
      real rlat00,rlon00
      real rlatnxny,rlonnxny
      real rlatdxdy,rlondxdy
      real dx,dy
      real golatsbp,golonsbp,goalpha
     
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      istatus = 1
      l_lut_flag=.false.

      do i = 1,nchannels
      call lvd_file_specifier(chtype(i),ispec,istatus)

      if(ispec.eq.1)then

         resx=r_resolution_x_vis(jtype,isat)
         resy=r_resolution_y_vis(jtype,isat)
         nlin=n_lines_vis(jtype,isat)
         npix=n_pixels_vis(jtype,isat)
         nw_line=i_nwline_vis(jtype,isat)
         nw_pix=i_nwpix_vis(jtype,isat)
      elseif(ispec.eq.2.or.ispec.eq.4.or.ispec.eq.5)then

         resx=r_resolution_x_ir(jtype,isat)
         resy=r_resolution_y_ir(jtype,isat)
         nlin=n_lines_ir(jtype,isat)
         npix=n_pixels_ir(jtype,isat)
         nw_line=i_nwline_ir(jtype,isat)
         nw_pix=i_nwpix_ir(jtype,isat)

      elseif(ispec.eq.3)then

         resx=r_resolution_x_wv(jtype,isat)
         resy=r_resolution_y_wv(jtype,isat)
         nlin=n_lines_wv(jtype,isat)
         npix=n_pixels_wv(jtype,isat)
         nw_line=i_nwline_wv(jtype,isat)
         nw_pix=i_nwpix_wv(jtype,isat)

      endif

c check if namelist parameters are current

c
c --- WFO --- 
c
      if(c_sat_types(jtype,isat).eq.'wfo')then 

         call get_wfo_nav_parms(path_to_raw_sat(ispec,jtype,isat)
     &,chtype(i),centerlat,centerlon,rlat00,rlon00,rlatnxny,rlonnxny,
     &rlatdxdy,rlondxdy,dx,dy,nx,ny,istatus)

         if(istatus.eq.0)then

           ilat00=nint(rlat00*1000.)
           ilon00=nint(rlon00*1000.)
           i_la1 =nint(r_la1(jtype,isat)*1000.)
           i_lo1 =nint(r_lo1(jtype,isat)*1000.)
           iresx =nint(resx)
           iresy =nint(resy)
           idx   =nint(dx)
           idy   =nint(dy)

           if((ilat00.ne.i_la1).or.
     &        (ilon00.ne.i_lo1) )then
             write(6,*)'rlat00/rlon00/r_la1/rlo1 ',rlat00,rlon00,
     &r_la1(jtype,isat),r_lo1(jtype,isat)
              l_lut_flag=.true.
           endif
           if( idx.ne.iresx .or. idy.ne.iresy )then
              write(6,*)'dx/dy/resx/resy/ ',dx,dy,resx,resy
              l_lut_flag=.true.
           endif
           if( nlin.ne.ny .or. npix.ne.nx )then
              write(6,*)'nx/ny/npix/nlin ',nx,ny,npix,nlin
              l_lut_flag=.true.
           endif

         endif

      elseif(c_sat_types(jtype,isat).eq.'gwc')then
c
c --- AFWA ---
c for GOES data only
         cname=c_afwa_fname(c_sat_id(isat),chtype(i))
         call s_len(cname,il)
         lenf=index(path_to_raw_sat(ispec,jtype,isat),' ')-1
         cfname=path_to_raw_sat(ispec,jtype,isat)(1:lenf)//cname(1:il)

         call read_gwc_header(cfname,l_cell_afwa,
     &strpix,strscnl,stppix,stpscnl,
     &reqobstm,imgtype,golatsbp,golonsbp,iwidth,idepth,goalpha,istrbdy1,
     &istrbdy2,istpbdy1,istpbdy2,bepixfc,bescnfc,fsci,idecimat,istatus)

         if(istatus.eq.0)then
            if(l_cell_afwa)then
               if(iwidth*256.ne.npix.or.idepth*64.ne.nlin)then
                  l_lut_flag=.true.
               endif
            else
               nxaf=stppix-strpix+1
               nyaf=stpscnl-strscnl+1
               if(nxaf.ne.npix.or.nyaf.ne.nlin)then
                  l_lut_flag=.true.
               endif
            endif
c
c NOTE: nw_line and nw_pix are hardwired by src/include/satdata_lvd.f
c       however, gen_lut_gvar.f uses the new values for these variables.
c       this will always trip the lut regeneration unless nw_line/pix match.
c
            istrtline=nw_vis_line_gwc(chtype(i),idecimat,bescnfc,fsci)
            istrtpix=nw_vis_pix_gwc(chtype(i),idecimat,bepixfc,goalpha)
            if(istrtline.ne.nw_line.or.istrtpix.ne.nw_pix)then
               l_lut_flag=.true.
            endif

         else
            print*,'gwc header not read. No lut update'
         endif

      elseif(c_sat_types(jtype,isat).eq.'gvr')then

         if(resx.eq.0.0.or.resy.eq.0.0.or.
     &      nlin.eq.0.0.or.npix.eq.0.0)then
        
            l_lut_flag = .true.
         endif

      endif


      enddo
c
c this check is for gvar type data only. once a day using
c a current O&A block.
c
      if(cfname_cur(6:8).eq.'000'.and.
     &   (c_sat_types(jtype,isat).eq.'gvr'   .or.
     &    c_sat_types(jtype,isat).eq.'gwc')       )then
         l_lut_flag=.true.
c        write(6,*)'Auto-update the gvar navigation'
      endif

c 8-31-06: JS. forces an a-ok return and no regeneration of
c mapping. All mapping now done on the fly.

      l_lut_flag=.false.

      return
      end
