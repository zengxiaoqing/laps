      subroutine get_bkgd_mdl_info(bgmodel,cmodel,fullname
     &,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     &,gproj,dlat,dlon,centrallat,centrallon,dxbg,dybg
     &,Lat0,Lat1,Lon0,sw,ne,cgrddef,istatus)
c
c JSmart 04-2001
c
c     USE laps_static

      implicit none

      include 'grid_fname.cmn'

      character*200 fullname
      character*132 cmodel
      character*30  projname
      character*13  fname13
      character*13  fname9_to_wfo_fname13
      character*4   cf
      character*2   gproj
      character*1   cgrddef
      
      integer       i,j
      integer       istatus
      integer       nxbg,nybg,nzbg
      integer       nzbg_ht
      integer       nzbg_tp
      integer       nzbg_sh
      integer       nzbg_uv
      integer       nzbg_ww
      integer       lenfn,nclen,lenn,leng
      integer       bgmodel
      integer       record
      integer       n_valtimes

      real          Lat0,Lat1
      real          Lon0
      real          La1in,La2in
      real          Lo1in,Lo2in
      real          dlat,dlon
      real          sw(2),ne(2)
      real          latdxdy,londxdy
      real          rlon00,rlat00
      real          latnxny,lonnxny
      real          centrallat,centrallon
      real          dxbg,dybg
      real          rotation

      interface

        subroutine get_eta48_dims(filename,NX,NY,NZ
     &,StdLat1,StdLat2,Lon0,La1,Lo1,La2,Lo2,istatus)
          character*200 filename
          integer NZ, NX, NY 
          integer istatus
          real    StdLat1,StdLat2
          real    Lon0
          real    La1,Lo1
          real    La2,Lo2

        end subroutine

        subroutine get_ruc2_dims(filename,NX,NY,NZ
     &,StdLat1,StdLat2,Lon0,La1,Lo1,La2,Lo2,istatus)
          character*200 filename
          integer NZ, NX, NY
          integer istatus
          real    StdLat1,StdLat2
          real    Lon0
          real    La1,Lo1
          real    La2,Lo2
        end subroutine

        subroutine get_attribute_sbn(cdfname,centralLat,centralLon,
     &rlat00,rlon00,latNxNy,lonNxNy,latdxdy,londxdy,dx,dy,nx,ny,
     &rotation,projname,istatus)
          character cdfname*200
          character projname*30
          integer   nx,ny
          real      centralLat
          real      centralLon
          real      rlat00
          real      rlon00
          real      dx,dy
          real      latNxNy
          real      lonNxNy
          real      latdxdy
          real      londxdy
          real      rotation
        end subroutine

        subroutine get_sbn_dims(cdfname,cmodel
     +,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     +,n_valtimes,istatus)
          character*132 cmodel
          character*200 cdfname
          integer nxbg,nybg
          integer nzbg_ht
          integer nzbg_tp
          integer nzbg_sh
          integer nzbg_uv
          integer nzbg_ww
          integer n_valtimes 
          integer istatus
        end subroutine

        subroutine readavnpublicdims(fname,x,y,numIsoLevel,record,
     +istatus)
          character*200 fname
          integer       numIsoLevel
          integer       record, x, y
          integer       nf_fid, nf_vid, nf_status
          integer       istatus
        end subroutine

      end interface


      print*,'Here: get_bkgd_mdl_info'
      istatus=1
      call s_len(fullname,lenfn)
      call s_len(cmodel,nclen)

      if(bgmodel.eq.0)then 
c      if(cmodel(1:nclen).eq.'LAPS_FUA')
c    &then
c        call find_domain_name(generic_data_root,grid_fnam_common,
c    &istatus)
c        call get_horiz_grid_spec(generic_data_root)
c        call s_len(grid_type,leng)
c        if(grid_type(1:5).eq. 'polar')gproj='PS'
c        if(grid_type(1:17).eq.'lambert conformal')gproj='LC'
c        if(grid_type(1:8).eq. 'mercator')gproj='MC'
         
c        nxbg=x
c        nybg=y
c        nzbg_ht=nk
c        nzbg_tp=nk
c        nzbg_sh=nk
c        nzbg_uv=nk
c        nzbg_ww=nk
c        sw(1)=la1
c        sw(2)=lo1
c        ne(1)=la2
c        ne(2)=lo2
c        Lon0=lov
c        Lat0=latin1
c        Lat1=latin2
c      elseif(cmodel(1:nclen).eq.'MODEL_FUA')then
c        call get_fua_dims()
c      endif

      endif
c ETA Public
c ----------
      if(bgmodel.eq.2.and.cmodel(1:nclen).eq.'ETA48_CONUS')
     &then
         call get_eta48_dims(fullname,nxbg,nybg,nzbg
     &         ,Lat0,Lat1,Lon0,La1in,Lo1in,La2in,Lo2in,istatus)
         if(istatus.eq.1)then
            gproj='LC'
            nzbg_ht=nzbg
            nzbg_tp=nzbg
            nzbg_sh=nzbg
            nzbg_uv=nzbg
            nzbg_ww=nzbg
            sw(1)=La1in
            sw(2)=Lo1in
            ne(1)=La2in
            ne(2)=Lo2in
         else
            print*,'Error - get_eta48_dims: ',fullname(1:lenfn)
         endif
      endif

c All SBN grids!
c ----------------
      if(bgmodel.eq.4)then

         j=lenfn
         j=j-13
         if(index(fullname(j+1:j+13),'/').eq.0 .and.
     +            fullname(j:j).eq.'/') then
            fname13=fname9_to_wfo_fname13(fullname(j+1:j+9))
c           cf=fullname(j+6:j+9)
c           fullname=fullname(1:j)//fname13//cf
            fullname=fullname(1:j)//fname13            !//cf
         else
            print*,'cannot convert fullname to WFO format'
            print*,'in get_bkgd_mdl_info ',fullname(1:lenfn)
         endif

         call s_len(fullname,lenfn)
         call get_sbn_dims(fullname,cmodel
     +,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     +,n_valtimes,istatus)

         if(istatus.ne. 1)then
            print*,'Error: get_sbn_dims'
	    print*,bgmodel,cmodel(1:nclen)
	    return
         endif

         print*,'call get_attribute_sbn'
         call get_attribute_sbn(fullname,centralLat
     +,centralLon,rlat00,rlon00,latNxNy,lonNxNy,latdxdy
     +,londxdy,dxbg,dybg,nxbg,nybg,rotation,projname,istatus)

         if(istatus.ne. 1)then
            print*,'Error: get_attribute_sbn'
            print*,bgmodel,cmodel(1:nclen)
            return
         endif

c set projection type for gridconv.f

         call s_len(projname,lenn)
         if(projname(1:lenn).eq.'LAMBERT_CONFORMAL')gproj='LC'
         if(projname(1:lenn).eq.'STEREOGRAPHIC')gproj='PS'
         if(projname(1:lenn).eq.'CYLINDRICAL_EQUIDISTANT')
     .gproj='LE'

         if(cmodel(1:nclen).eq.'RUC40_NATIVE'.or.
     .      cmodel(1:nclen).eq.'ETA48_CONUS')then

            nzbg_tp=nzbg_tp-1
            nzbg_uv=nzbg_uv-1
            nzbg_sh=nzbg_sh-1
            print*,'Retrieved SBN attributes for ',cmodel(1:nclen)

         elseif(cmodel(1:nclen).eq.'AVN_SBN_CYLEQ')then

c for global AVN, nav code expects grid 1,1 in nw corner
            print*,'set return variables'
            rlat00 =-1.*rlat00
            latNxNY=-1.*latNxNy
            nzbg_tp=nzbg_tp-2
            nzbg_uv=nzbg_uv-2
            nzbg_sh=nzbg_sh-2
            print*,'Retrieved SBN attributes for ',cmodel(1:nclen)
         elseif(cmodel(1:7).eq.'MesoEta')then
            print*,'Retrieved SBN attributes for ',cmodel(1:nclen)
         else 
            print*,'Unknown SBN model type: cmodel = ',cmodel
         endif

         Lon0=centralLon
         Lat0=centralLat
         Lat1=Lat0        !this has be the second latitude (tangent lambert) since no Lat1.
         dlat=dxbg/111.1
         dlon=dybg/111.1
         sw(1)=rlat00 
         sw(2)=rlon00
         ne(1)=latNxNy
         ne(2)=lonNxNy

      endif


c RUC Public
c ----------
      if(bgmodel.eq.5.and.cmodel(1:nclen).eq.'RUC40_NATIVE')then
         call get_ruc2_dims(fullname,nxbg,nybg,nzbg
     &              ,Lat0,Lat1,Lon0,La1in,Lo1in,La2in,Lo2in,istatus)
         if(istatus.eq.1)then
            nzbg_ht=nzbg
            nzbg_tp=nzbg
            nzbg_sh=nzbg
            nzbg_uv=nzbg
            nzbg_ww=nzbg
            Lon0=Lon0-360.
            gproj='LC'
            sw(1)=La1in
            sw(2)=-Lo1in    !-360.
            ne(1)=La2in
            ne(2)=Lo2in
         else
            print*,'Error - get_ruc2_dims: ',fullname(1:lenfn)
         endif
      endif
      
c AVN Public
c ----------
      if(bgmodel.eq.6.and.cmodel(1:nclen).eq.'AVN_FSL_NETCDF')
     &then
         call readavnpublicdims(fullname,nxbg,nybg,nzbg,record
     +,istatus)
         if(istatus.eq.1)then
            nzbg_ht=nzbg
            nzbg_tp=nzbg
            nzbg_sh=nzbg
            nzbg_uv=nzbg
            nzbg_ww=nzbg
            gproj='LL'
            Lat0=90.0   !although consistent with AVN Public, doesn't work
            Lon0=0.0    !with gridconv latlon_2_llij
            dlat=1.0 
            dlon=1.0
            cgrddef='N'
         else
            print*,'Error - readavnpublicdims '
         endif

      endif

c AVN binary (at AFWA)
c --------------------
      if(bgmodel.eq.6.and.cmodel(1:nclen).eq.'AVN_AFWA_DEGRIB')
     &then
         nxbg=360
         nybg=181
         nzbg=26
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         gproj='LL'
         Lat0=-90.0
         Lon0=0.0
         dlat=1.0
         dlon=1.0
         cgrddef='S'
      endif

c Taiwan FA and NF models
c -----------------------
      if(bgmodel.eq.3.and.cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_RE')
     &then
         nxbg  = 91
         nybg  = 91
         nzbg  = 16
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         gproj='LC'
         Lat0=10.0
         Lat1=40.0
         Lon0=+120.
         sw(1)=15.879
         sw(2)=+112.545
         ne(1)=32.384
         ne(2)=+131.172
      endif
      if(bgmodel.eq.3.and.cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_NF')
     &then
         nxbg = 145
         nybg = 139
         nzbg = 11
         nzbg_ht=nzbg
         nzbg_tp=nzbg
         nzbg_sh=nzbg
         nzbg_uv=nzbg
         nzbg_ww=nzbg
         gproj='LC'
         Lat0=10.0
         Lat1=40.0
         Lon0=+120.
         sw(1)=15.80
         sw(2)=+109.24
         ne(1)=34.987
         ne(2)=+131.60
      endif

c     if (bgmodel .eq. 7) then
c        gproj='LC'
c        nx_lc=nx
c        ny_lc=ny
c        nz_lc=nz
c        lat1=25.0
c        lat1_lc=lat1
c        lat2=25.0
c        lat2_lc=lat2
c        lon0=-95.0
c        lon0_lc=lon0
c        sw(1)=12.19
c        sw(2)=-133.459
c        ne(1)=57.29
c        ne(2)=-49.3849
c     elseif (bgmodel.eq.8)then
c        gproj='LL'
c        nx_ll=nx
c        ny_ll=ny
c        nz_ll=nz
c        lat0=-90.0
c        lon0=0.0
c        lon0_lc=lon0
c        dlat=1.0
c        dlon=1.0
c     endif


c     print*,'done in get_bkgd_mdl_info'
c     print*,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww

      return
      end
