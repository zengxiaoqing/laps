      subroutine readlut(csat_id,csat_typ,maxch,nft,
     &ntm,chtype,nx,ny,ri,rj,istatus)
c
      implicit none

      integer nx,ny
      integer nx_in
      integer ny_in
      integer maxch
      integer istatus
      integer i,j,n,np
      integer lend,lenf
      integer ispec
      integer istat
      integer nft
      integer ntm(nft)

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
      character*3   chtype(maxch,nft)
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

      do j=1,nft
      do i=1,ntm(j)

         call lvd_file_specifier(chtype(i,j),ispec,istat)
         if(istat.eq.0)then

            if(.not.lgot_lut(ispec))then
              lgot_lut(ispec)=.true.

              ct=chtype(i,j)
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
      enddo

      istatus = 0
      goto 1000

23    write(6,*)'Error reading or eof ll/ij lookup table'
      goto 1000

101   write(6,*)'Error opening file ',file(1:n)

1000  return
      end
c
c================================================
c
      subroutine check_luts(cfname_cur,isat,jtype,
     &chtype,maxchannels,nchannels,l_lut_flag,istatus)

      implicit none

      integer        maxchannels
      integer        nchannels

      character*(*)  cfname_cur
      character*3    chtype(maxchannels)

      integer        i
      integer        isat,jtype,kch
      integer        istatus
      integer        nlin,npix
      integer        nx,ny
      integer        ispec
      integer        ilat00,ilon00
      integer        i_la1,i_lo1

      logical        l_lut_flag

      real resx,resy
      real rlat00,rlon00
      real dx,dy
     
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      istatus = 0

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

      l_lut_flag=.false.
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

      endif

      enddo

      if(cfname_cur(6:8).eq.'000'.and.
     &c_sat_types(jtype,isat).eq.'gvr'.or.
     &c_sat_types(jtype,isat).eq.'gwc')then
         l_lut_flag=.true.
         write(6,*)'Auto-update the gvar navigation'
      endif

      return
      end
