      subroutine readlut(csat_id,csat_typ,maxch,nch,
     &chtype,nx,ny,ri,rj,istatus)
c
      implicit none

      integer nx,ny
      integer nx_in
      integer ny_in
      integer maxch
      integer nch
      integer istatus
      integer i,j,n,np
      integer lend,lenf
      integer ispec
      integer istat

      logical   lgot_lut(maxch)

      real*4    ri(nx,ny,maxch)
      real*4    rj(nx,ny,maxch)
      real*4    ri_in(nx,ny)
      real*4    rj_in(nx,ny)
      real*4    rdummy(nx,ny)

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

      do i=1,maxch
         lgot_lut(i)=.false.
      enddo

      do i=1,nch

         call lvd_file_specifier(chtype(i),ispec,istat)
         if(istat.eq.0)then

            if(.not.lgot_lut(ispec))then
              lgot_lut(ispec)=.true.

              ct=chtype(i)
              goto(3,2,3,2,2)ispec
2             ct='ir'

3             n=index(ct,' ')-1
              if(n.le.0)n=3
              file=cpath(1:lend)//csat_id//'-llij-'
              lenf=index(file,' ')-1
              file=file(1:lenf)//ct(1:n)//'-'//csat_typ//'.lut'

              n=index(file,' ')-1
              open(12,file=file,
     &form='unformatted',status='old',err=101)
              write(6,*)'Reading ',file(1:n)
              read(12,err=23,end=23) rdummy
              read(12,err=23,end=23) rdummy
              read(12,err=23,end=23) ri_in
              read(12,err=23,end=23) rj_in
              close (12)
              call move(ri_in,ri(1,1,ispec),nx,ny)
              call move(rj_in,rj(1,1,ispec),nx,ny)
            endif
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
      character      cname*11
      character      c_afwa_fname*11
      character      cfname*200

      integer        i
      integer        isat,jtype,kch
      integer        istatus
      integer        nlin,npix
      integer        nx,ny
      integer        ispec
      integer        ilat00,ilon00
      integer        i_la1,i_lo1
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
      real rlat00,rlon00
      real dx,dy
      real golatsbp,golonsbp,goalpha
     
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      istatus = 1
      l_lut_flag=.false.

      do i = 1,nchannels
      call lvd_file_specifier(chtype(i),ispec,istatus)

      goto(10,11,12,11,11)ispec

10       resx=r_resolution_x_vis(jtype,isat)
         resy=r_resolution_y_vis(jtype,isat)
         nlin=n_lines_vis(jtype,isat)
         npix=n_pixels_vis(jtype,isat)
         goto 13

11       resx=r_resolution_x_ir(jtype,isat)
         resy=r_resolution_y_ir(jtype,isat)
         nlin=n_lines_ir(jtype,isat)
         npix=n_pixels_ir(jtype,isat)
         goto 13

12       resx=r_resolution_x_wv(jtype,isat)
         resy=r_resolution_y_wv(jtype,isat)
         nlin=n_lines_wv(jtype,isat)
         npix=n_pixels_wv(jtype,isat)
         goto 13

13    continue

c check if namelist parameters are current

c
c --- WFO --- 
c
      if(c_sat_types(jtype,isat).eq.'wfo')then 

         call get_wfo_nav_parms(path_to_raw_sat(ispec,jtype,isat),
     &                          chtype(i),rlat00,rlon00,dx,dy,nx,ny,
     &                          istatus)

         ilat00=nint(rlat00*1000.)
         ilon00=nint(rlon00*1000.)
         i_la1 =nint(r_la1(jtype,isat)*1000.)
         i_lo1 =nint(r_lo1(jtype,isat)*1000.)

         if((ilat00.ne.i_la1).or.
     &      (ilon00.ne.i_lo1) )then
            write(6,*)'rlat00/rlon00/r_la1/rlo1 ',rlat00,rlon00,
     &r_la1(jtype,isat),r_lo1(jtype,isat)
            l_lut_flag=.true.
         endif
         if( dx.ne.resx .or. dy.ne.resy )then
            write(6,*)'dx/dy/resx/resy/ ',dx,dy,resx,resy
            l_lut_flag=.true.
         endif
         if( nlin.ne.ny .or. npix.ne.nx )then
            write(6,*)'nx/ny/npix/nlin ',nx,ny,npix,nlin
            l_lut_flag=.true.
         endif

      elseif(c_sat_types(jtype,isat).eq.'gwc')then
c
c --- AFWA ---
c
         cname=c_afwa_fname(c_sat_id(isat),chtype(i))
         lenf=index(path_to_raw_sat(ispec,jtype,isat),' ')-1
         cfname=path_to_raw_sat(ispec,jtype,isat)(1:lenf)//cname

         call read_gwc_header(cfname,strpix,strscnl,stppix,stpscnl,
     &reqobstm,imgtype,golatsbp,golonsbp,iwidth,idepth,goalpha,istrbdy1,
     &istrbdy2,istpbdy1,istpbdy2,bepixfc,bescnfc,fsci,idecimat,istatus)

         if(istatus.eq.0)then
c
c this test is only good for the old SDHS data files.
c
            if(iwidth*256.ne.npix.or.idepth*64.ne.nlin)then
               l_lut_flag=.true.
            endif
c
c this test will help for both the old and new SDHS data files.
c
            istrtline = nw_vis_line_gwc(chtype(i),bescnfc,fsci)
            istrtpix =  nw_vis_pix_gwc(chtype(i),bepixfc,goalpha)

         else
            print*,'gwc header not read. No lut update'
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
         write(6,*)'Auto-update the gvar navigation'
      endif

      return
      end
