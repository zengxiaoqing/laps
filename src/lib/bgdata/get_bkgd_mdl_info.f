      subroutine get_bkgd_mdl_info(bgmodel,cmodel,fullname
     &,mxvars,mxlvls,nx,ny,nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww
     &,gproj,dlat,dlon,centrallat,centrallon,dx,dy
     &,Lat0,Lat1,Lon0,sw,ne,istatus)
c
c JSmart 04-2001
c
      implicit none

      include 'bgdata.inc'

      integer       mxvars,mxlvls

      character*200 fullname
      character*132 cmodel
      character*2   gproj
      character*9   fname13
      character*9   fname9_to_wfo_fname13
      character*4   cf
      
      integer       i,j
      integer       nvars
      integer       istatus
      integer       nx,ny
      integer       nzbg_ht
      integer       nzbg_sh
      integer       nzbg_uv
      integer       nzbg_ww
      integer       nz
      integer       lenfn,nclen
      integer       bgmodel
      integer       record
      integer       n_valtimes
      integer       idims(mxlvls,mxvars)
      integer       levels(mxlvls,mxvars)

      real          Lat1,Lat2
      real          Lat0,Lon0
      real          La1,La2
      real          Lo1,Lo2
      real          dlat,dlon
      real          sw(2),ne(2)
      real          latdxdy,londxdy
      real          rlon00,rlat00
      real          latnxny,lonnxny
      real          centrallat,centrallon
      real          dx,dy
      real          rotation

      istatus=1
      call s_len(fullname,lenfn)
      call s_len(cmodel,nclen)

c ETA Public
c ----------
      if(bgmodel.eq.2.and.cmodel(1:nclen).eq.'ETA48_CONUS')
     &then
         call get_eta48_dims(fullname,nx,ny,nz
     &         ,Lat0,Lat1,Lon0,La1,Lo1,La2,Lo2,istatus)
         if(istatus.eq.1)then
            gproj='LC'
            nzbg_ht=nz
            nzbg_sh=nz
            nzbg_uv=nz
            nzbg_ww=nz  !double check this ... its good!
            sw(1)=La1
            sw(2)=Lo1
            ne(1)=La2
            ne(2)=Lo2
         else
            print*,'Error - get_eta48_dims: ',fullname(1:lenfn)
         endif
         goto 1000
      endif

c All SBN grids!
c ----------------
      if(bgmodel.eq.4)then

         j=lenfn
         j=j-13
         if(index(fullname(j+1:j+13),'/').eq.0 .and.
     +            fullname(j:j).eq.'/') then
            fname13=fname9_to_wfo_fname13(fullname(j+1:j+9))
            cf=fullname(j+6:j+9)
            fullname=fullname(1:j)//fname13//cf
         else
            print*,'cannot convert fullname to WFO format'
            print*,'in get_bkgd_mdl_info ',fullname(1:lenfn)
            goto 1000
         endif

         call get_sbn_dims(fullname,cmodel,mxvars,mxlvls
     +,nvars,nx,ny,nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww
     +,n_valtimes,istatus)

         if(istatus.ne. 1)then
            print*,'Error: get_sbn_dims'
	    print*,bgmodel,cmodel(1:nclen)
	    return
         endif

         call get_attribute_sbn(fullname,centralLat
     +,centralLon,rlat00,rlon00,latNxNy,lonNxNy,latdxdy,londxdy
     +,dx,dy,nx,ny,rotation,istatus)

         if(istatus.ne. 1)then
            print*,'Error: get_attribute_sbn'
            print*,bgmodel,cmodel(1:nclen)
            return
         endif

         if(cmodel(1:nclen).eq.'RUC40_NATIVE'.or.
     .      cmodel(1:nclen).eq.'ETA48_CONUS')then
            gproj='LC'

         elseif(cmodel(1:nclen).eq.'AVN_SBN_CYLEQ')then

c for global AVN, nav code expects grid 1,1 in nw corner
            rlat00 =-1.*rlat00
            latNxNY=-1.*latNxNy
            gproj='LE'
         else
            print*,'Unknown SBN model type: cmodel = ',cmodel
         endif

         Lon0=centralLon
         Lat0=centralLat
         Lat1=Lat0        !this has be the second latitude (tangent lambert) since no Lat1.
         dlat=dx/111.1
         dlon=dy/111.1
c        dx=dx*1000.
c        dy=dy*1000.
         sw(1)=rlat00 
         sw(2)=rlon00
         ne(1)=latNxNy
         ne(2)=lonNxNy

         goto 1000

      endif


c RUC Public
c ----------
      if(bgmodel.eq.5.and.cmodel(1:nclen).eq.'RUC40_NATIVE')then
         call get_ruc2_dims(fullname,nx,ny,nz
     &              ,Lat0,Lat1,Lon0,La1,Lo1,La2,Lo2,istatus)
         if(istatus.eq.1)then
            nzbg_ht=nz
            nzbg_sh=nz
            nzbg_uv=nz
            nzbg_ww=nz
            Lon0=Lon0-360.
            gproj='LC'
            sw(1)=La1
            sw(2)=-Lo1    !-360.
            ne(1)=La2
            ne(2)=Lo2
         else
            print*,'Error - get_ruc2_dims: ',fullname(1:lenfn)
         endif
         goto 1000
      endif
      
c AVN Public
c ----------
      if(bgmodel.eq.6.and.cmodel(1:nclen).eq.'AVN_FSL_NETCDF')
     &then
         call readavnpublicdims(fullname,nx,ny,nz,record
     +,istatus)
         if(istatus.eq.1)then
            nzbg_ht=nz
            nzbg_sh=nz
            nzbg_uv=nz
            nzbg_ww=nz
            gproj='LL'
            Lat0=90.0
            Lon0=0.0
            dlat=1.0
            dlon=1.0
         else
            print*,'Error - readavnpublicdims '
         endif
         goto 1000

      endif

c AVN binary (at AFWA)
c --------------------
      if(bgmodel.eq.6.and.cmodel(1:nclen).eq.'AVN_AFWA_DEGRIB')
     &then
         nx=360
         ny=181
         nz=26
         nzbg_ht=nz
         nzbg_sh=nz
         nzbg_uv=nz
         nzbg_ww=nz
         gproj='LL'
         Lat0=-90.0
         Lon0=0.0
         dlat=1.0
         dlon=1.0
         goto 1000
      endif

c Taiwan FA and NF models
c -----------------------
      if(bgmodel.eq.3.and.cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_RE')
     &then
         nx    = 91
         ny    = 91
         nz    = 16
         nzbg_ht=nz
         nzbg_sh=nz
         nzbg_uv=nz
         nzbg_ww=nz
         gproj='LC'
         Lat0=10.0
         Lat1=40.0
         Lon0=+120.
         sw(1)=15.879
         sw(2)=+112.545
         ne(1)=32.384
         ne(2)=+131.172
         goto 1000
      endif
      if(bgmodel.eq.3.and.cmodel(1:nclen).eq.'CWB_20FA_LAMBERT_NF')
     &then
         nx   = 145
         ny   = 139
         nz   = 11
         nzbg_ht=nz
         nzbg_sh=nz
         nzbg_uv=nz
         nzbg_ww=nz	 !double check this!
         gproj='LC'
         Lat0=10.0
         Lat1=40.0
         Lon0=+120.
         sw(1)=15.80
         sw(2)=+109.24
         ne(1)=34.987
         ne(2)=+131.60
         goto 1000
      endif

      print*,'get_bkgd_mdl_info not yet working for: '
     &,'   ',bgmodel,cmodel(1:nclen)

1000  return
      end
